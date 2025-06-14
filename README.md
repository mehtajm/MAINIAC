# MAINIAC

**MAFFT-based Adaptive Immune Numbering and Immunoglobulin or Antibody Concatenation (MAINIAC)**

This tool processes antibody sequences using reference alignment and gap handling, leveraging both Python and C++ components.

---

## 📦 Installation

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/mainiac.git
cd mainiac
```

### 2. Create and activate a Conda environment

You can create a new Conda environment with the required dependencies:

```bash
conda create -n mainiac python=3.9
conda activate mainiac
```

Then install dependencies (update this list as needed):

```bash
conda install pandas biopython numpy
```

> You may also use `pip` for any missing dependencies.

### 3. Compile the C++ binaries

```bash
make
```

This will compile the C++ tools (`MafftGapConverter` and `concatenate`) into the `bin/` directory.

---

## 🚀 Usage

### Prepare the input data

Unzip and extract the test input files:

```bash
gunzip input.gz && tar -xvf input
```

### Run MAINIAC

```bash
python ./src/mainiac.py \
  --source_file ./testData/airr_IGH-filtered.tsv \
  --seq_col_name sequence_aa \
  --ref_v ./testData/ighv.F.X.ref \
  --ref_j ./testData/ighj.ref \
  --output filtered_IGH
```

---

## 🛠 Development Notes

- C++ executables are expected in the `bin/` directory after running `make`.
- If you modify any `.cpp` files in `src/backend/`, rerun `make` to recompile.
- To clean compiled binaries:

```bash
make clean
```

---

## 📁 Directory Structure (simplified)

```
mainiac/
├── bin/                    # Compiled C++ binaries (after make)
├── src/
│   ├── backend/            # C++ source files
│   └── mainiac.py          # Main Python script
├── testData/               # Sample input files
├── input.gz                # Compressed test input
├── Makefile
└── README.md
```

---



## License

MIT License. See [LICENSE](LICENSE) for details.

