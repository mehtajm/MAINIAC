# Jai Mehta
# 03/26/2025
# Python wrapper for MAFFT Insertions
import os, sys, glob
import numpy as np
import pandas as pd
import getpass
import subprocess
import argparse
import importlib
import shutil
import re
import logging
from pathlib import Path
import os

bin_path = Path(__file__).resolve().parent / "backend" / ".." / ".." / "bin"
converter_exec = bin_path / "MafftGapConverter"
concat_exec = bin_path / "concatenate"

def cleanData(args):
    """Remove leading './' from all string arguments in args."""
    for key, value in vars(args).items():
        if isinstance(value, str) and value.startswith("./"):
            setattr(args, key, value[2:])
    makeFasta(args)

def makeFasta(args):
    print("Making fasta file from input file")
    try:
        data = pd.read_csv(
            args.source, delimiter=args.sep, engine="python", quotechar='"', on_bad_lines="skip"
        )
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")
    path = str(os.path.abspath(os.getcwd()))
    pypath = str(Path(__file__).parent.absolute())
    if args.seqCol not in data.columns:
        sys.exit(f"Error: Column '{args.seqCol}' not found in the input file.")

    seqs = data[args.seqCol].dropna()  # Drop NaN values to avoid writing "nan" in the output

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Change working directory
    os.chdir(args.output)

    # Write sequences to a FASTA file
    with open("query_sequences.fa", "w") as file:
        for i, seq in enumerate(seqs):
            file.write(f">q{i}\n{seq.strip()}\n")  # Correct FASTA format
    print(f'File created successfully, named "query_sequences.fa"\n')
    alignJ(args, pypath, path)

def alignJ(args, pypath, path):
    print("Aligning J sequences using MAFFT.")
    mafftAlign = []
    mafftAlign.append(f"mafft --addfragments query_sequences.fa --compactmapout --thread -1 {path}/{args.refJ} > /dev/null")
    os.system(mafftAlign[-1])
    print("Aligned sequences. Now extracting J regions and V + CDR3 regions.")
    #os.system(f"g++ -o MafftGapConverter {pypath}/backend/MafftGapConverter.cpp")
    os.system("./MafftGapConverter query_sequences.fa query_sequences.fa.map delete.fa query j")
    os.system("mv query query_vSeqs.fa && rm delete.fa")
    alignV(args, pypath, path)

def alignV(args, pypath, path):
    print("Aligning v + cdr3 sequences")
    os.system(f"mafft --addfragments query_vSeqs.fa --compactmapout --thread -1 {path}/{args.refV} > mafftoutV")
    os.system("./MafftGapConverter query_vSeqs.fa query_vSeqs.fa.map mafftoutV v.number.csv")
    print("Aligned v + cdr3 sequences to v.number.csv. Aligning j sequences to j.number.csv")
    os.system(f"mafft --addfragments query_jSeqs.fa --compactmapout --thread -1 {path}/{args.refJ} > mafftoutJ")
    os.system("./MafftGapConverter query_jSeqs.fa query_jSeqs.fa.map mafftoutJ j.number.csv")
    print("now concatenating")
    #os.system(f"g++ -o concatenate {pypath}/backend/concatenate.cpp &&"              
    os.system("./concatenate v.number.csv j.number.csv out.csv")
    #os.system("rm concatenate && rm MafftGapConverter")
    print('Completed. Concatenated csv saved to "out.csv"')
def parseArgs():
    parser = argparse.ArgumentParser(
        prog = 'Mafft-based Adaptive Immune Numbering and Immunoglobulin or Antibody Concatenation (MAINIAC)',
        description='Given amino acid sequences, the program uses MAFFT to align query sequences to reference files, and assigns a numbering scheme. ')
    parser.add_argument('--source_file', dest='source', required=True,
                        help='Source file containing sequences.')
    parser.add_argument('--seq_col_name', dest='seqCol', required=True,
                        help='The sequence column name in the source file that contains the sequences.')
    parser.add_argument('--sep', dest='sep', required=False, default='\t',
                        help="Used to specify the separator of the input file is (optional, default '\t').")
    parser.add_argument('--ref_v', dest='refV', required=True, 
                        help='Reference file of V genes, in fasta format.')
    parser.add_argument('--ref_j', dest='refJ', required=True,
                        help='Reference file of J genes, in fasta format.')
    parser.add_argument('--output_folder', dest='output', default="output",
                        help='Output folder name.')
    return parser.parse_args()



if __name__ == "__main__":
    args = parseArgs()
    cleanData(args)
