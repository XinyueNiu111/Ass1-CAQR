#ifndef TSQR_H
#define TSQR_H

/*
 * TSQR: Communication-Avoiding QR Factorization
 *
 * Does TSQR on a tall-skinny matrix A (m x n)
 * which is  distributed evenly across 4 processors.
 * Uses LAPACK dgeqrf for local QR factorizations.
 *
 * Parameters:
 *   A     - input matrix, column-major (m x n) 
 *   m     - number of rows
 *   n     - number of columns
 *   R_out - output upper-triangular factor (n x n)
 */
void TSQR(double *A, int m, int n, double *R_out); 

#endif