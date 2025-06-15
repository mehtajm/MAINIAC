import argparse
import subprocess
from pathlib import Path
import sys
import time
import re
import string
import shlex
import pandas as pd
from uuid import uuid4


PROGRAM_HOME = Path(__file__).resolve().parent

## go up two level then up one more to get the parent dir
## which contains both “MAINIAC” and “disease-classifier”
UTILS_HOME = PROGRAM_HOME.parent.parent / "disease-classifier"
sys.path.insert(0, str(UTILS_HOME))

from utils import merge_CDRs, parse_CDRs,create_job_id, get_ref_colnames,check_done, merge_CDRs


def csvs_to_db(csv_files, db_dir, *, len_sort=False, clonotype=False, save_msa=False):
    """
    Given a list of CSV outputs from run_mainiac.py, collate them into a single
    .npz DB plus .size (and optional TSV) under db_dir.

    Parameters
    ----------
    csv_files : List[pathlib.Path]
        The individual output_csv paths from each job.
    db_dir : str or pathlib.Path
        Directory where we should write the merged DB.
    len_sort : bool
        If True, output base name gets '-len_sort' suffix and sequences
        will be sorted by pseudo_len.
    clonotype : bool
        If True, builds a clonotype DB (includes v_call/j_call).
    save_msa : bool
        If True, also emit a `{base}.tsv` of pseudo-sequences.
    """
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)




    # 1) fix each CSV in-place so merge_CDRs can parse it
    #    enumerate so we know which subject (file) each row came from
    for subject_idx, csv_path in enumerate(csv_files):
        df = pd.read_csv(csv_path)
        # prefix the local row‐index with the subject index
        df["Id"] = f"{subject_idx}-" + df.index.astype(str)
        df["e-value"] = 0.0
        if "sequence_id" in df.columns:
            df.drop(columns=["sequence_id"], inplace=True)
        df.to_csv(csv_path, index=False)
        print(f"Fixed CSV for merge_CDRs: {csv_path.name}")

    # 2) build base name
    base = "clonotype" if clonotype else "pseudo_seq"
    if len_sort:
        base += "-len_sort"
    base_path = db_dir / base

    # 3) decide output paths
    pseudo_file = str(base_path) + ".tsv" if save_msa else None
    npz_file    = str(base_path) + ".npz"
    size_file   = str(base_path) + ".size"

    # 4) merge into one DB
    entries_no = merge_CDRs(
        csv_files,
        pseudo_file=pseudo_file,
        npz_file=npz_file,
        len_sort=len_sort,
        clonotype=clonotype,
        subjects_no=len(csv_files),
    )

    # 5) write .size
    with open(size_file, "w") as fh:
        fh.write(f"{entries_no}\n")

    print(
        f"DB files generated: {npz_file}, {size_file}"
        + (f", {pseudo_file}" if pseudo_file else "")
    )


    
def n_lsf_jobs(jobname: str, queue: str) -> int:
    """
    Get the number of LSF jobs for a given job name and queue.
    """
    job_cmd = f"bjobs -J {jobname} -q {queue} | grep -v JOBID"
    status, stdout = subprocess.getstatusoutput(job_cmd)
    if status:
        return 0
    jobs = stdout.strip().split("\n")
    return len(jobs) if jobs[0] else 0

def lsf_wait(jobname: str, max_jobs: int, queue: str):
    """
    Sleep until the number of jobs is less than max_jobs.
    """
    while n_lsf_jobs(jobname, queue) > max_jobs:
        time.sleep(10)

def parse_args():
    parser = argparse.ArgumentParser(description="Run MAINIAC on multiple input files via LSF.")
    parser.add_argument("-i","--input", type=Path, required=True,
                        help="Input file or directory containing .tsv files.")
    parser.add_argument("-o","--output_dir", type=Path, required=True,
                        help="Directory to store all output CSVs.")
    parser.add_argument("--seq_col_name", default="sequence_aa",
                        help="Name of the column with sequences.")
    parser.add_argument("--ref_v", type=Path, 
                        help="Path to V reference FASTA.")
    parser.add_argument("--ref_j", type=Path, 
                        help="Path to J reference FASTA.")
    parser.add_argument("--keep-temp", action="store_true",
                        help="Keep temp directories even on success.")
    parser.add_argument("--keep-align", action="store_true",
                        help="Keep intermediate alignment outputs.")
    default_script = Path(__file__).resolve().parent / "mainiac.py"
    parser.add_argument("--mainiac_script", type=Path, default=default_script,
                        help="Path to mainiac.py script (default: mainiac.py in same directory).")
    parser.add_argument("--no-lsf", action="store_true",
                        help="Run jobs directly instead of using bsub.")
    parser.add_argument("--max-jobs", type=int, default=100,
                        help="Maximum number of concurrent LSF jobs.")
    parser.add_argument("--queue", default="alma8-batch",
                        help="LSF queue to use (default: alma8-batch).")
    parser.add_argument("--echo", action="store_true",
                        help="Only print jobs instead of submitting.")
    parser.add_argument("--job-name", 
                        help="Base job name for LSF submissions.")
    parser.add_argument("--db",action="store_true",help="if given, collect all per-job CSVs into a single DB under this dir")
    parser.add_argument("--len_sort",action="store_true",
                        help="speeding up by sorting entries by pseudo_len")
    parser.add_argument("--clonotype", action="store_true",
                        help="use clonotype.")
    parser.add_argument("--save_msa", action="store_true",
                        help="save MSA.")
    return parser.parse_args()

def main():
    args = parse_args()
    job_id=create_job_id()
    if args.job_name:
        job_id = args.job_name
        
    input_files = []
    if args.input.is_file():
        input_files = [args.input.resolve()]
    elif args.input.is_dir():
        input_files = sorted(args.input.glob("*.tsv"))
    else:
        sys.exit(f"Invalid input: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    output_csvs=[]
    for i, f in enumerate(input_files):
        output_csv = args.output_dir / (f.stem + ".out.csv")
        output_csvs.append(output_csv)
        cmd = [
            sys.executable, str(args.mainiac_script),
            "--source_file", str(f),
            "--seq_col_name", args.seq_col_name,
            "--output_file", str(output_csv)
        ]
        if args.ref_v:
            cmd.append("--ref_v")
            cmd.append(str(args.ref_v))
        if args.ref_j:
            cmd.append("--ref_j")
            cmd.append(str(args.ref_j))

        if args.keep_temp:
            cmd.append("--keep-temp")
        if args.keep_align:
            cmd.append("--keep-align")

        
        job_cmd = " ".join(shlex.quote(arg) for arg in cmd)

        if args.no_lsf:
            launch_cmd = job_cmd
        else:
            log_path = args.output_dir / f"job_{i}.out"
            launch_cmd = f"bsub -J {job_id} -q {args.queue} -o {log_path} {job_cmd}"

        if args.echo:
        
        
            if args.no_lsf:
                print("Would submit job:", job_cmd)
            else:
                quoted_launch = f"bsub -J {job_id} -q {args.queue} -o {log_path} {job_cmd}"
                print("Would submit job:", quoted_launch)

        else:
            print(f"Submitting job: {launch_cmd}")
            subprocess.run(launch_cmd, shell=True)
            if not args.no_lsf:
                lsf_wait(job_id, args.max_jobs, args.queue)

    if not args.no_lsf:
        print("Waiting for all jobs to complete...")
        lsf_wait(job_id, 0, args.queue)

        

    if  args.db:
        check_done(job_id)
        csvs_to_db(output_csvs,args.output_dir,len_sort=args.len_sort,clonotype=args.clonotype,save_msa=args.save_msa)


    print("All jobs completed.")

if __name__ == "__main__":
    main()
