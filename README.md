# MAP55672: Case Studies in HPC
## Case 1: Communication-Avoiding QR Factorization

**Author:** Xinyue Niu


## Structure

`python/tsqr.ipynb` — Question 1: TSQR implementation in Python
`tsqr_c/tsqr.c` — Question 2 & 3: TSQR implementation in C using LAPACK
`tsqr_c/Makefile` — Build and run the C code

## run

**Python (Question 1):**
Open `python/tsqr.ipynb` and run all cells.

**C (Questions 2 & 3):**
```bash
cd tsqr_c
make
make run
```