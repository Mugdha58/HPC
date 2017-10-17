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
    int i,j,k,n=2048,i1,j1,k1,block_size;
    double *a,*b,*c1,*c2,cpu_time,sum,error=0.0,difference,gflops,temp;
    register double r;
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
    for(block_size=4;block_size<17;block_size=block_size*2)
    {
            printf("BLOCK %d calculation:",block_size);
         //implementing ijk loop
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        for (i = 0; i < n; i+=block_size)
	for (j = 0; j < n; j+=block_size)
             for (k = 0; k < n; k+=block_size)
		 //block multiplication
                  for (i1 = i; i1 < i+block_size; i1++)
                      for (j1 = j; j1 < j+block_size; j1++){
                           r=c1[i1*n+j1];
                          for (k1 = k; k1 < k+block_size; k1++)
	                       r+= a[i1*n + k1]*b[k1*n + j1];
	                       c1[i1*n+j1]=r;
                      }
	 clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\n CPU time for ijk is  %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    //implementing jik loop
     clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (j = 0; j < n; j+=block_size)
    for (i = 0; i < n; i+=block_size)
        for (k = 0; k < n; k+=block_size)
		 //block multiplication
            for (j1 = j; j1 < j+block_size; j1++)
                  for (i1 = i; i1 < i+block_size; i1++){
                        r=c2[i1*n+j1];
                      for (k1 = k; k1 < k+block_size; k1++)
	                       r+= a[i1*n + k1]*b[k1*n + j1];
                            c2[i1*n+j1]=r;
                  }
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for jik loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    for(i=0;i<n*n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
        printf("\n the error value is %f ",error);


     //implementing kij
     clock_gettime(CLOCK_MONOTONIC, &cstart);
     for (k = 0; k < n; k+=block_size)
        for (i = 0; i < n; i+=block_size)
            for (j = 0; j < n; j+=block_size)
		 //block multiplication
                for (k1 = k; k1 < k+block_size; k1++)
                 for (i1 = i; i1 < i+block_size; i1++)
                        {
                        r=a[i1*n+k1];
                        for (j1 = j; j1 < j+block_size; j1++)
	                      c1[i1*n+j1] += r*b[k1*n + j1];

                        }
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for kij loop for block size is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    //implementing ikj
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (i = 0; i < n; i+=block_size)
        for (k = 0; k < n; k+=block_size)
            for (j = 0; j < n; j+=block_size)
		 //block multiplication
                 for (i1 = i; i1 < i+block_size; i1++)
                        for (k1 = k; k1 < k+block_size; k1++){
                        r=a[i1*n+k1];
                        for (j1 = j; j1 < j+block_size; j1++)
	                      c2[i1*n+j1] += r*b[k1*n + j1];
                        }
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for ikj loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    for(i=0;i<n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
    printf("\n the error value is %f ",error);
    //implementing jki
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (j = 0; j < n; j+=block_size)
        for (k = 0; k < n; k+=block_size)
            for (i = 0; i < n; i+=block_size)
            //block multiplication
            for (j1 = j; j1 < j+block_size; j1++)
                for (k1 = k; k1 < k+block_size; k1++)
                {
                r=b[k1*n+j1];
                  for (i1 = i; i1 < i+block_size; i1++)
                      c1[i1*n+j1] +=r*a[i1*n + k1];
                }
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for jki loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);

    //implementing kji
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    for (k = 0; k < n; k+=block_size)
        for (j = 0; j < n; j+=block_size)
            for (i = 0; i < n; i+=block_size)
            //block multiplication
            for (k1 = k; k1 < k+block_size; k1++)
            for (j1 = j; j1 < j+block_size; j1++)

                {
                r=b[k1*n+j1];
                  for (i1 = i; i1 < i+block_size; i1++)
                      c2[i1*n+j1] +=r*a[i1*n + k1];
                }

    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for kji loop is %f",cpu_time);
    gflops=(2*pow(n,3))/(cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    for(i=0;i<n;i++)
    {
       difference=(abs)(c1[i]-c2[i]);
       if(difference>error)
        error=difference;
    }
    printf("\n the error value is %f ",error);
    }
    free(a);
    free(b);
    free(c1);
    free(c2);

    return 0;
}

