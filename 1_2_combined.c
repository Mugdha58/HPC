#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "lapacke.h"
#include "blas.h"



double Random_gen ( int upper, int lower)
{
    double s;
    s = ((double)rand()/(RAND_MAX))*(upper-lower);
    return s;
}

void mydgetrf(double *A,int *pvt, double *tempv, int n){
    int i,t,j,k,maxind,temps;

    double max;
    for(i=0;i<n-1;i++){
        maxind = i;
        max=abs(A[i*n+i]);
        for(t=i+1;t<n;t++){
            if(abs(A[t*n+i])>max){
                maxind = t;
                max = abs(A[t*n+i]);
            }
        }
        if(max==0){
            printf("LU factorization failed: coefficient matrix is singular\n");
            return;
        }
        else{
            if(maxind != i){
                //Save pivoting information
                temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;
                //Swap rows
                for(k=0;k<n;k++){
                    tempv[k] = A[i*n+k];
                    A[i*n+k] = A[maxind*n+k];
                    A[maxind*n+k] = tempv[k];
                }
            }
        }
        for(j=i+1;j<n;j++){
            A[j*n+i] = A[j*n+i]/A[i*n+i];
            for(k=i+1;k<n;k++){
                A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
            }
        }
    }
}

void mydtrsm_f(int n, double *A, double *B, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,j,k,t;
    y[0] = B[pvt[0]];
    for(i=1;i<n;i++){
        sum = 0.0;
        for(k=0;k<i;k++){
            sum += y[k]*A[i*n+k];
        }
        y[i] = B[pvt[i]]-sum;
    }
}

void mydtrsm_b(int n, double *A, double *B, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,j,k,t;
    x[n-1] = y[n-1]/A[(n-1)*n+(n-1)];
    for(i=n-2;i>=0;i--){
        sum=0.0;
        for(k=i+1;k<n;k++){
            sum+= x[k]*A[i*n+k];
        }
        x[i] = (y[i]-sum)/A[i*n+i];
    }
}

void transpose(double *a, int n){
    int i,j;
    double temp;
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            temp = a[i*n+j];
            a[i*n+j] = a[j*n+i];
            a[j*n+i] = temp;
        }
    }
}

int main()
{   srand((double)time(NULL));
    double time,gflops;
    int u=10,l=1;
    //int size = (sizeof(arrayLen)/sizeof(arrayLen[0]));
    double ran = Random_gen(u,l);
    int i,j,k,n,t,temp;
    printf("Using LAPACK Library\n");
    for(n=1000;n<6000;n=n+1000)
    {
        struct timespec tstart={0,0},tend={0,0};
        char TRANS = 'N';
        int INFO = n;
        int LDA = n;
        int LDB = n;
        int N = n;
        int NRHS = 1;
        int *IPIV = (int *)calloc(sizeof(int),n);
        double  *A, *A1, *B, *B1, *x, *y, *abk, *tempv, difference, error =0.0;
        int *pvt;
        A=(double *) calloc(sizeof(double), n*n);
        B=(double *) calloc(sizeof(double), n);
        A1=(double *) calloc(sizeof(double), n*n);
        B1=(double *) calloc(sizeof(double), n);
        pvt=(int *) calloc(sizeof(int), n);
        y=(double *) calloc(sizeof(double), n);
        x=(double *) calloc(sizeof(double), n);
        tempv=(double *) calloc(sizeof(double), n);
        for(i=0;i<n;i++){
            for(j=0;j<n;j++)
            {
              A[i*n+j]=Random_gen(u,l);
              A1[i*n+j]=A[i*n+j];
            }
        }
    for(i=0;i<n;i++){
        B[i]=Random_gen(u,l);
        B1[i]=B[i];
        pvt[i]=i;
    }
    transpose(A,n);

        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
        printf("\nLAPACK LIBRARY\n");
        clock_gettime(CLOCK_MONOTONIC,&tstart);
        LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC,&tend);
        double time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);

        for(i = 0; i < N; i++)
        {
            double tmp = B[IPIV[i]-1];
        	B[IPIV[i]-1] = B[i];
        	B[i] = tmp;
        }

        // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        UPLO = 'U';
        DIAG = 'N';

        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        printf("Size N = %d\n",n);
        printf("Time Taken = %.5f seconds\n",time);
        double gflops = (2*pow(n,3))/(3*time*pow(10,9));
        printf("\nPerformance in GFLOPS = %f\n",gflops);
        printf("\n");
        clock_gettime(CLOCK_MONOTONIC,&tstart);
        mydgetrf(A1,pvt, tempv,n);
        clock_gettime(CLOCK_MONOTONIC,&tend);
        mydtrsm_f(n,A1,B1,pvt,x,y);

        mydtrsm_b(n,A1,B1,pvt,x,y);
        printf("MYDGETRF VERSION\n");
        time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        gflops = (2*pow(n,3))/(3*time*pow(10,9));
        printf("Time Taken = %.5f seconds\n",time);
        printf("\nPerformance in GFLOPS = %f\n",gflops);
        printf("\n");
        for(i=0;i<n*n;i++)
        {
           difference = abs(B[i]-x[i]);
           if(difference>error)
            error=difference;
        }
        printf("\n the error value for n=%d is %f ",n,error);
        free(A);
        free(B);
        free(A1);
        free(B1);
        free(pvt);
        free(x);
        free(y);
        free(IPIV);
        free(tempv);
    }
    return 0;
}
