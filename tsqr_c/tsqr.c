#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "tsqr.h"

// LAPACK: QR factorization based householder to general matrix

extern void dgeqrf_(int *m, int *n, double *A, int *lda,
                    double *tau, double *work, int *lwork,
                    int *info);

/* Helper: get upper-triangular R from dgeqrf output
 * dgeqrf stores R in the upper triangle of A 
 */
static void extract_R(double *A, int m, int n, double *R)
{
    // Zero out R
    memset(R, 0, n * n * sizeof(double));

    // Only copy upper triangle 
    for (int j = 0; j < n; j++) {
        for (int i = 0; i <= j; i++) {
            R[i + j * n] = A[i + j * m];     // column-major order
        }
    }
}

/* Helper: local QR by LAPACK dgeqrf
 * Input:  A_block (column-major)
 * Output: R_out   (n x n, upper triangular)
 */
static void local_qr(double *A_block, int rows, int n,
                     double *R_out)
{
    int lda   = rows;
    int lwork = 10 * n;
    int info  = 0;

    double *tau    = (double *)malloc(n * sizeof(double));
    double *work   = (double *)malloc(lwork * sizeof(double));
    double *A_copy = (double *)malloc(rows * n * sizeof(double));

    // Copy input matrix
    memcpy(A_copy, A_block, rows * n * sizeof(double));

    // Call LAPACK Householder QR 
    dgeqrf_(&rows, &n, A_copy, &lda, tau, work, &lwork, &info);

    if (info != 0) {
        fprintf(stderr, "dgeqrf failed with info = %d\n", info);
        exit(1);
    }

    // Extract R from upper triangle of A_copy 
    extract_R(A_copy, rows, n, R_out);

    free(tau);
    free(work);
    free(A_copy);
}

/* Helper: stack two n x n matrices vertically into 2n x n
 * Result store as column-major order
 */
static double *stack_two(double *top, double *bot, int n)
{
    double *out = (double *)malloc(2 * n * n * sizeof(double));
    for (int col = 0; col < n; col++) {
        memcpy(out + col * 2 * n,
               top + col * n,
               n * sizeof(double));
        memcpy(out + col * 2 * n + n,
               bot + col * n,
               n * sizeof(double));
    }
    return out;
}

/* TSQR: main function
 * Three steps of binary tree reduction:
 *   Step1: 4 local QRs on row blocks
 *   Step2: QR on paired R factors
 *   Step3: 1 final QR on last pair
 */
void TSQR(double *A, int m, int n, double *R_out)
{
    int bs = m / 4;  

    printf("\n Step1: Local QR on 4 row blocks \n");

    double *R1[4];
    for (int i = 0; i < 4; i++) {
        R1[i] = (double *)malloc(n * n * sizeof(double));

        // Copy block i in column-major order 
        double *block = (double *)malloc(bs * n * sizeof(double));
        for (int j = 0; j < n; j++) {
            memcpy(block + j * bs,
                A + j * m + i * bs,
                bs * sizeof(double));
            }
            local_qr(block, bs, n, R1[i]);
            free(block);
            printf("Processor %d: local QR done (rows %d to %d)\n",
                i, i * bs, (i + 1) * bs - 1);
    }

    printf("\n Step2: Pair reductions \n");

    double *R2[2];
    int pairs[2][2] = {{0, 1}, {2, 3}};

    for (int p = 0; p < 2; p++) {
        int i = pairs[p][0];
        int j = pairs[p][1];

        // Stack R1[i] on top of R1[j]
        double *stacked = stack_two(R1[i], R1[j], n);

        R2[p] = (double *)malloc(n * n * sizeof(double));
        local_qr(stacked, 2 * n, n, R2[p]);

        printf("Pair (%d,%d): QR done\n", i, j);
        free(stacked);
    }

    printf("\n Step3: Final reduction \n");

    double *stacked_final = stack_two(R2[0], R2[1], n);
    local_qr(stacked_final, 2 * n, n, R_out);

    for (int i = 0; i < 4; i++) free(R1[i]);
    for (int i = 0; i < 2; i++) free(R2[i]);
    free(stacked_final);
}

/* Verify: compute ||A^T A - R^T R||_F
 * if TSQR is correct， result should be close to zero
 */
static double verify(double *A, int m, int n, double *R)
{
    // Compute A^T A
    double *ATA = (double *)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < m; k++) {
                sum += A[k + i * m] * A[k + j * m];
            }
            ATA[i + j * n] = sum;
        }
    }

    // Compute R^T R 
    double *RTR = (double *)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k <= i && k <= j; k++) {
                sum += R[k + i * n] * R[k + j * n];
            }
            RTR[i + j * n] = sum;
        }
    }

    // Frobenius norm of difference 
    double err = 0.0;
    for (int i = 0; i < n * n; i++) {
        double d = ATA[i] - RTR[i];
        err += d * d;
    }

    free(ATA);
    free(RTR);
    return sqrt(err);
}

int main(void)
{
    // Scale with respect to m (fix n = 5)
    int n_fixed = 5;
    int m_values[] = {40, 80, 160, 320, 640, 1280, 2560};
    int num_m = 7;

    FILE *fm = fopen("scale_m.csv", "w");
    fprintf(fm, "m,n,time\n");

    for (int t = 0; t < num_m; t++) {
        int m = m_values[t];
        int n = n_fixed;

        // Build random matrix
        double *A = (double *)malloc(m * n * sizeof(double));
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                A[i + j * m] = sin((double)(i + 1) * (double)(j + 1) * 0.3);

        double *R = (double *)malloc(n * n * sizeof(double));

        // Time TSQR 
        int repeats = 10;
        clock_t start = clock();
        for (int r = 0; r < repeats; r++) {
            TSQR(A, m, n, R);
        }
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC / repeats;

        fprintf(fm, "%d,%d,%.6f\n", m, n, elapsed);
        printf("m=%d, n=%d, time=%.6f s\n", m, n, elapsed);

        free(A);
        free(R);
    }
    fclose(fm);

    // Scale with respect to n (fix m = 400) 
    int m_fixed = 400;
    int n_values[] = {5, 10, 20, 40, 80, 100};
    int num_n = 6;

    FILE *fn = fopen("scale_n.csv", "w");
    fprintf(fn, "m,n,time\n");

    for (int t = 0; t < num_n; t++) {
        int m = m_fixed;
        int n = n_values[t];

        double *A = (double *)malloc(m * n * sizeof(double));
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                A[i + j * m] = sin((double)(i + 1) * (double)(j + 1) * 0.3);

        double *R = (double *)malloc(n * n * sizeof(double));
        
        int repeats = 10;
        clock_t start = clock();
        for (int r = 0; r < repeats; r++) {
            TSQR(A, m, n, R);
        }
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC / repeats;

        fprintf(fn, "%d,%d,%.6f\n", m, n, elapsed);
        printf("m=%d, n=%d, time=%.6f s\n", m, n, elapsed);

        free(A);
        free(R);
    }
    fclose(fn);

    return 0;
}