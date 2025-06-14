import argparse
import subprocess
from pathlib import Path
import sys
import time

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
    default_ref_v = Path(__file__).resolve().parent / "testData/ighv.F.X.ref"    
    parser.add_argument("--ref_v", type=Path, default=default_ref_v,
                        help="Path to V reference FASTA.")
    default_ref_j = Path(__file__).resolve().parent / "testData/ighj.ref"    
    parser.add_argument("--ref_j", type=Path, default=default_ref_j,
                        help="Path to J reference FASTA.")
    parser.add_argument("--sep", default="\t",
                        help="Input file separator (default: tab).")
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
    parser.add_argument("--job-name", default="mainiac_job",
                        help="Base job name for LSF submissions.")
    return parser.parse_args()

def main():
    args = parse_args()
    input_files = []
    if args.input.is_file():
        input_files = [args.input.resolve()]
    elif args.input.is_dir():
        input_files = sorted(args.input.glob("*.tsv"))
    else:
        sys.exit(f"Invalid input: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for i, f in enumerate(input_files):
        output_csv = args.output_dir / (f.stem + ".out.csv")
        cmd = [
            sys.executable, str(args.mainiac_script),
            "--source_file", str(f),
            "--seq_col_name", args.seq_col_name,
            "--sep", repr(args.sep).strip("'"),
            "--ref_v", str(args.ref_v),
            "--ref_j", str(args.ref_j),
            "--output_file", str(output_csv)
        ]
        if args.keep_temp:
            cmd.append("--keep-temp")
        if args.keep_align:
            cmd.append("--keep-align")

        job_cmd = " ".join(cmd)

        if args.no_lsf:
            launch_cmd = job_cmd
        else:
            log_path = args.output_dir / f"job_{i}.out"
            launch_cmd = f"bsub -J {args.job_name}_{i} -q {args.queue} -o {log_path} {job_cmd}"

        if args.echo:
            import shlex
            quoted_cmd = " ".join(shlex.quote(arg) for arg in cmd)
            if args.no_lsf:
                print("Would submit job:", quoted_cmd)
            else:
                quoted_launch = f"bsub -J {args.job_name}_{i} -q {args.queue} -o {log_path} {quoted_cmd}"
                print("Would submit job:", quoted_launch)

        else:
            print(f"Submitting job: {launch_cmd}")
            subprocess.run(launch_cmd, shell=True)
            if not args.no_lsf:
                lsf_wait(args.job_name, args.max_jobs, args.queue)

    if not args.no_lsf:
        print("Waiting for all jobs to complete...")
        lsf_wait(args.job_name, 0, args.queue)
        print("All jobs completed.")

if __name__ == "__main__":
    main()
