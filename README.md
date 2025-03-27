# MAINIAC
Mafft-based Adaptive Immune Numbering and Immunoglobulin or Antibody Concatenation (MAINIAC)

## Example usage

# First unzip test input file

```
gunzip input.gz && tar -xvf input
```
# Example program usage
```
python ./src/mainiac.py --source_file ./testData/airr_IGH-filtered.tsv --seq_col_name sequence_aa --ref_v ./testData/ighv.F.X.ref --ref_j ./testData/ighj.ref --output filtered_IGH
```
