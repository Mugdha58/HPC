#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0;
    return((double)rand()/upper_bound);

}
int main()
{
    int i,j,k,n=2048;
    double *a,*b,*c1,*c2,difference,error=0.0,temp,gflops;
    register double sum,r;
    struct timespec cstart = {0,0}, cend ={0,0};
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
    //implementing ijk loop
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            sum = c1[i*n+j];
            for (k=0; k<n; k++)
            sum += a[i*n+k] * b[k*n+j];
            c1[i*n+j] = sum;
  }
}
	clock_gettime(CLOCK_MONOTONIC, &cend);
        double cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("CPU time for ijk loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    //implementing jik loop
     clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (j=0; j<n; j++)
    {
        for (i=0; i<n; i++)
        {
            sum =c2[i*n+j];
            for (k=0; k<n; k++)
            sum += a[i*n+k] * b[k*n+j];
            c2[i*n+j] = sum;
        }
    }
	clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for jik loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    for(i=0;i<n*n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
        printf("\n the error value is %f ",error);
     clock_gettime(CLOCK_MONOTONIC, &cstart);
    //implementing kij
    for (k=0; k<n; k++) {
        for (i=0; i<n; i++) {
            r = a[i*n+k];
            for (j=0; j<n; j++)
            c1[i*n+j] += r * b[k*n+j];
        }
    }
	clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for kij loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    //implementing ikj
    for (i=0; i<n; i++)
    {
        for (k=0; k<n; k++) {
        r = a[i*n+k];
        for (j=0; j<n; j++)
        c2[i*n+j] += r * b[k*n+j];
        }
    }
 clock_gettime(CLOCK_MONOTONIC, &cend);      
 cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for ikj loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    for(i=0;i<n*n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
        printf("\n the error value is %f ",error);
     clock_gettime(CLOCK_MONOTONIC, &cstart);
    //implementing jki
    for (j=0; j<n; j++) {
        for (k=0; k<n; k++) {
            r = b[k*n+j];
        for (i=0; i<n; i++)
        c1[i*n+j] += a[i*n+k] * r;
        }
    }
	clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for jki loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    //implemeting kji
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            r = b[k*n+j];
            for (i=0; i<n; i++)
            c2[i*n+j] += a[i*n+k] * r;
        }
    }
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for kji loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%.16f",gflops);
    for(i=0;i<n*n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
    printf("\n the error value is %f ",error);
    free(a);
    free(b);
    free(c1);
    free(c2);

    return 0;
}
