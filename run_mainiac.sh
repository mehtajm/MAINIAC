#!/usr/bin/env bash
# Quick smoke test for mainiac
python ./src/mainiac.py --source_file ./testData/test_igh.tsv --seq_col_name sequence_aa --ref_v ./testData/ighv.F.X.ref --ref_j ./testData/ighj.ref --output_file ./filtered_IGH.csv
