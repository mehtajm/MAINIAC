import sys
import shutil
import tempfile
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

# Paths to external executables (bin directory is one level above 'src')
BIN_DIR = Path(__file__).resolve().parent.parent / "bin"
CONVERTER = BIN_DIR / "MafftGapConverter"
CONCATENATOR = BIN_DIR / "concatenate"


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Maine Pythonic MAINIAC",
        description="Align and number immune sequences using MAFFT and custom tools."
    )
    parser.add_argument(
        "--source_file", dest="source", type=Path, required=True,
        help="TSV/CSV source with sequence column."
    )
    parser.add_argument(
        "--seq_col_name", dest="seq_col", required=True,
        help="Name of the column with sequences."
    )
    parser.add_argument(
        "--sep", dest="sep", default="\t",
        help="Field separator (default: tab)."
    )
    parser.add_argument(
        "--ref_v", dest="ref_v", type=Path, required=True,
        help="FASTA reference for V genes."
    )
    parser.add_argument(
        "--ref_j", dest="ref_j", type=Path, required=True,
        help="FASTA reference for J genes."
    )
    parser.add_argument(
    "--output_file", dest="output_file", type=Path, required=True,
    help="Final output CSV filename (e.g. 'something.out.csv')"
    )

    parser.add_argument(
        "--keep-temp", dest="keep_temp", action="store_true",
        help="Retain intermediate temp files for debugging."
    )
    parser.add_argument(
        "--keep-align", dest="keep_align", action="store_true",
        help="Retain intermediate alignment files."
    )
    args = parser.parse_args()
    return args


def run(cmd, cwd=None):
    logging.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)


def main():
    args = parse_args()

    source = args.source.resolve()
    ref_v = args.ref_v.resolve()
    ref_j = args.ref_j.resolve()
    output_file = args.output_file.resolve()
    output_dir = output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    

    try:
        df = pd.read_csv(source, sep=args.sep, engine="python", on_bad_lines="skip")
    except Exception as e:
        logging.error(f"Failed to read {args.source}: {e}")
        sys.exit(1)

    if args.seq_col not in df.columns:
        logging.error(f"Sequence column '{args.seq_col}' not found.")
        sys.exit(1)

    sequences = df[args.seq_col].dropna().astype(str)


    # create temp dir manually
    tmpdir = tempfile.mkdtemp(prefix="mainiac_", dir="/tmp")
    tmp = Path(tmpdir)
    success = False
    try:
        # 0) write input FASTA
        input_fa = tmp / "queries.fa"
        with input_fa.open("w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">q{i}\n{seq.strip()}\n")

        # 1) Initial J alignment + split
        run(f"mafft --addfragments {input_fa} --compactmapout --thread 1 {ref_j} > /dev/null", cwd=tmp)
        j_map = tmp / f"{input_fa.name}.map"
        run(f"{CONVERTER} {input_fa} {j_map} delete.fa query j", cwd=tmp)
        (tmp / "query").rename(tmp / "query_vSeqs.fa")
        # converter creates query_jSeqs.fa

        # 2) V+CDR3 alignment + V numbering
        v_in = tmp / "query_vSeqs.fa"
        v_out = tmp / "mafft_v.out"
        run(f"mafft --addfragments {v_in} --compactmapout --thread 1 {ref_v} > {v_out}", cwd=tmp)
        v_map = tmp / f"{v_in.name}.map"
        run(f"{CONVERTER} {v_in} {v_map} {v_out} v.number.csv", cwd=tmp)

        # 3) J fragment alignment + J numbering
        j_in = tmp / "query_jSeqs.fa"
        j_out = tmp / "mafft_j2.out"
        run(f"mafft --addfragments {j_in} --compactmapout --thread 1 {ref_j} > {j_out}", cwd=tmp)
        j_map2 = tmp / f"{j_in.name}.map"
        run(f"{CONVERTER} {j_in} {j_map2} {j_out} j.number.csv", cwd=tmp)

        # 4) Concatenate results
        run(f"{CONCATENATOR} v.number.csv j.number.csv out.csv", cwd=tmp)

        # Move outputs: final CSVs
        #for fn in ["out.csv", "v.number.csv", "j.number.csv"]:
        #    shutil.move(str(tmp / fn), str(output_dir / fn))
        shutil.move(str(tmp / "out.csv"), str(output_file))

        if args.keep_align:
            # Also copy aligned FASTA & map files for inspection
            
            align_files = [
                input_fa.name,
                f"{input_fa.name}.map",
                v_in.name,
                f"{v_in.name}.map",
                j_in.name,
                f"{j_in.name}.map",
                v_out.name,
                j_out.name
            ]
            for fn in align_files:
                src = tmp / fn
                if src.exists():
                    shutil.move(str(src), str(output_dir / fn))

        success = True
        logging.info(f"Finished. Results in {args.output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)
    finally:
        if args.keep_temp or not success:
            logging.info(f"Retained temp directory: {tmpdir}")
        else:
            shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main()
