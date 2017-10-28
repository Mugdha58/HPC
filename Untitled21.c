#include <stdio.h>
#include "lapacke.h"
#include "blas.h"
#include<math.h>
#include<time.h>
double Random_gen()
{
    double upper_bound=RAND_MAX/10.0;
    return ( (double) rand() / upper_bound );
}
int main()
{
    char    TRANS = 'N';
    int     n,i,j;
    double gflops,cpu_time;
    struct timespec cstart = {0,0}, cend ={0,0};
    double *A,*B;
    for(n=1000;n<6000;n=n+1000)
    {
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
	int 	IPIV[n];
 int     NRHS = 1;

        A=(double *) calloc(sizeof(double), n*n);
        B=(double *) calloc(sizeof(double), n*1);
        for(i=0;i<n;i++)
        for(j=0;j<n;j++)
   	 {
              A[i*n+j]=(double)Random_gen();
              B[j]=(double)Random_gen(); 
    	}	
    // LU factorization
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        LAPACK_dgetrf(&n,&n,A,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\nCPU time for LU factorization n=%d is %f",n,(cpu_time/CLOCKS_PER_SEC));
        gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
        printf("\nthe gflops used are=%f",gflops);
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;

        for(i = 0; i < n; i++)
        {
            double tmp = B[IPIV[i]-1];
            B[IPIV[i]-1] = B[i];
            B[i] = tmp;
        }
    // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
        UPLO = 'U';
        DIAG = 'N';
    // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
	free(A);
	free(B);
    }  
    return 0;
}

