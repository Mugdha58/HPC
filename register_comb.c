#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0,temp;
    return((double)rand()/upper_bound);

}
int main()
{
    int i,j,k,n;
    double *a,*b,*c1,*c2,cpu_time,gflops,error=0.00,difference;
    clock_t start,end;

    //implementation of dgemmo
for(n=64;n<2049;n=n*2)
{


    a=(double *) calloc(sizeof(double), n*n);
    b=(double *) calloc(sizeof(double), n*n);
    c1=(double *) calloc(sizeof(double), n*n);
    c2=(double *) calloc(sizeof(double), n*n);
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
    {
              a[i*n+j]=(double)Random_gen();
              b[i*n+j]=(double)Random_gen();
              temp=(double)Random_gen();
              c1[i*n+j]=temp;
              c2[i*n+j]=temp;


    }
    start= clock();
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            for(k=0;k<n;k++)
                c1[i*n+j]+=a[i*n+k]*b[k*n+j];

                end=clock();
    cpu_time=end-start/(CLOCKS_PER_SEC);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used by dgemmo  for n=%d are=%.16f",n,gflops);
    //implementation of dgemm1
    start= clock();
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            register double r=c2[i*n+j];
            for(k=0;k<n;k++)
                    r+=a[i*n+k]*b[k*n+j];
                    c2[i*n+j]=r;

        }
    end=clock();
    cpu_time=(double)(end-start)/(CLOCKS_PER_SEC);
    printf("\nCPU time for n=%d is %f",n,cpu_time);
    gflops=((double)(double)2*pow(n,3))/(double)(cpu_time*pow(10,9));
    printf("\nthe gflops used by dgemm1 n=%d are %.16f",n,gflops);
    for(i=0;i<n*n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
        printf("\n the error value for n=%d is %f ",n,error);
    free(a);
    free(b);
    free(c1);
    free(c2);

    }
    return 0;
}
