#include <stdio.h>
#include "lapacke.h"
#include "blas.h"
#include<math.h>
#include<time.h>
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0;
    return((double)rand()/upper_bound);

}
int main()
{
    char    TRANS = 'N';
    int     n=4,i;
    int     NRHS = 1;
    int     IPIV[3];
    double gflops;
    struct timespec cstart = {0,0}, cend ={0,0};
    double *A,*B;
   // for(n=1000;n<6000;n=n+1000)
    //{
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        A=(double *) calloc(sizeof(double), n*n);
        A={9,13,5,,2,1,11,7,6,3,7,4,1,6,0,7,10};
        B=(double *) calloc(sizeof(double), n*1);
        B={1,2,3,4};
    // LU factorization
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        LAPACK_dgetrf(&n,&n,A,&LDA,IPIV,&INFO);
         clock_gettime(CLOCK_MONOTONIC, &cend);
        double cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\nCPU time for LU factorization n=%d is %f",cpu_time);
        gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
        printf("\nthe gflops used are=%f",gflops);
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;

        for(int i = 0; i < n; i++)
        {
            double tmp = B[IPIV[i]-1];
            B[IPIV[i]-1] = B[i];
            B[i] = tmp;
        }
        clock_gettime(CLOCK_MONOTONIC, &cstart);
    // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
        clock_gettime(CLOCK_MONOTONIC, &cend);
        double cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\nCPU time for LU factorization n=%d is %f",cpu_time);
        gflops=(pow(n,2))/(cpu_time*pow(10,9));
        printf("\nthe gflops used are=%f",gflops);
        UPLO = 'U';
        DIAG = 'N';
    // backward Ux = y
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
        clock_gettime(CLOCK_MONOTONIC, &cend);
        double cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\nCPU time for n=%d is %f",cpu_time);
        gflops=(pow(n,2))/(cpu_time*pow(10,9));
        printf("\nthe gflops used are=%f",gflops);

    return 0;
}

