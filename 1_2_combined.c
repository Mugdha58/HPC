#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
#include "lapacke.h"
#include "blas.h"
int i,j,k,n,t,temp;
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0;
    return((double)rand()/upper_bound);

}
void transpose(double *a,int n) //to implement transpose as fortran has column wise implementation and the program is implemented in row wise implementation
{
    for(i=0;i<n;i++)
        for(j=i;j<n;j++)
        {
		temp=a[i*n+j];
		a[i*n+j]=a[j*n+i];
		a[j*n+i]=temp;
	}
}
void mydgetrf(double *a,int *pvt,int n,double *tempv)
{
    int i,j,k,t,maxind,temps;
    double max;
    for(i=0;i<n-1;i++)
    {
       maxind=i;
       max=abs(a[i*n+i]);
       for(t=i+1;t<n;t++)
       {

            if(abs(a[t*n+i])>max)
             {
                 maxind=t;
                 max=abs(a[t*n+i]);

             }
       }
       if(max==0)
       {

        printf("LU factorization failed:coefficient matrix is singular");
        return;
       }
       else
       {

        if(maxind!=i)
        {   //save pivoting information
            temps=pvt[i];
            pvt[i]=pvt[maxind];
            pvt[maxind]=temps;
//            for(k=0;k<n;i++)
//            {
//                tempv[k]=a[i*n+k];
//                a[i*n+k]=a[maxind*n+k];
//                a[maxind*n+k]=tempv[k];
//            }
        }
       }
//      //factorizing
//      for(j=i+1;j<n;k++)
//      {
//          a[j*n+i]=a[j*n+i]/a[i*n+i];
//          for(k=i+1;k<n;k++)
//            a[j*n+k]=a[j*n+k]-(a[j*n+i]*a[i*n+k]);
//      }

    }
}
void mydtrsm(int n,double *a,double *B,int *pvt,double *x,double *y,int label)
{
    double sum=0.0,temp;
    if(label==0)// passing label to call forward and backward substitution separately
    {//forward substitution
    y[0]=B[pvt[0]];
    for(i=1;i<n;i++)
    {
        sum=0.0;
      for(k=0;k<i;k++)
      {
        sum+=y[k]*a[i*n+k];

      }
      y[i]=B[pvt[i]]-sum;
    }
    }
    //backward substitution
    else
    {
    x[n-1]=y[n-1]/a[(n-1)*n+(n-1)];
    for(i=n-2;i>=0;i--){

        sum=0.0;
        for(k=i+1;k<n;k++)
    {
       sum+= x[k]*a[i*n+k];

    }
    temp=y[i]-sum;
    x[i]=temp/a[i*n+i];

    }
    }

}
int main()
{
    int *pvt,n=8;
    double *a,*B,*a1,*B1,*x,*y,*tempv;
    double error=0.0,difference;
    double gflops,cpu_time;
    struct timespec cstart = {0,0}, cend ={0,0};
   // for(n=1000;n<6000;n=n+1000)
    //{
    a=(double *) calloc(sizeof(double), n*n);
    B=(double *) calloc(sizeof(double), n);
    a1=(double *) calloc(sizeof(double), n*n);
    B1=(double *) calloc(sizeof(double), n);
    pvt=(int *) calloc(sizeof(int), n);
    y=(double *) calloc(sizeof(double), n);
    x=(double *) calloc(sizeof(double), n);
    tempv=(double *) calloc(sizeof(double), n);
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
    {
              a[i*n+j]=(double)Random_gen();
              a1[i*n+j]=a[i*n+j];
    }
    for(i=0;i<n;i++){
        B[i]=(double)Random_gen();
        B1[i]=B[i];
        pvt[i]=i;
    }
    transpose(a,n);
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    mydgetrf(a,pvt,n,tempv);
    clock_gettime(CLOCK_MONOTONIC, &cend);
//    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
//    printf("\nCPU time for LU factorization n=%d is %f",n,cpu_time);
//    gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
//    printf("\nthe gflops used are=%f",gflops);
//    mydtrsm(n,a,B,pvt,x,y,0);
//    mydtrsm(n,a,B,pvt,x,y,1); // label 1 is passed so that backward substitution will be done
//     for(i=0;i<n;i++)
//    printf("\nthe result of my implementation is %f\t",x[i]);
//    char    TRANS = 'N';
//    int     NRHS = 1;
//    int *IPIV = (int *)calloc(sizeof(int),n);
//        int     INFO = n;
//        int     LDA = n;
//        int     LDB = n;
//    // LU factorization
//        clock_gettime(CLOCK_MONOTONIC, &cstart);
//        LAPACK_dgetrf(&n,&n,a1,&LDA,IPIV,&INFO);
//        clock_gettime(CLOCK_MONOTONIC, &cend);
//        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
//        printf("\nCPU time for LU factorization n=%d is %f",n,cpu_time);
//        gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
//        printf("\nthe gflops used are=%f",gflops);
//        char     SIDE = 'L';        char     UPLO = 'L';
//        char     DIAG = 'U';
//        int      M    = 1;
//        double   b    = 1.0;
//
//        for(i = 0; i < n; i++)
//        {
//            double tmp = B1[IPIV[i]-1];
//            B1[IPIV[i]-1] = B1[i];
//            B1[i] = tmp;
//        }
//    // forward  L(Ux) = B => y = Ux
//        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&b,a1, &n, B1, &n);
//        UPLO = 'U';
//        DIAG = 'N';
//    // backward Ux = y
//        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&b,a1, &n, B1, &n);
//        for(i=0;i<n;i++)
//            printf("\nthe result of library function is %f\t",B1[i]);
//     for(i=0;i<n*n;i++)
//        {
//       difference=(abs)(B1[i]-x[i]);
//       if(difference>error)
//        error=difference;
//        }
//        printf("\n the error value for n=%d is %f ",n,error);
//    free(a);
//    free(B);
//    free(x);
//    free(y);
//    free(pvt);
//    free(a1);
//    free(B1);
//    free(tempv);
  //  }
    return 0;
}

