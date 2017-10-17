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
    int i,j,k,n;
    double *a,*b,*c1,*c2,cpu_time,gflops,temp,error=0.0,difference;
    register int t,tt,ta,tta,tb,ttb;
    register double c00,c01,c10,c11,a00,a01,a10,a11,b00,b01,b10,b11;
    struct timespec cstart = {0,0}, cend ={0,0};
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    for(n=66;n<258;n=n*2)
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

for(i = 0; i < n; i += 3)
       for(j = 0; j < n; j += 3)  {
            int t = i*n+j; int tt = t+n; int ttt=tt+n;
            register double c00 = c1[t]; register double c01 = c1[t+1];  register double c02 = c1[t+2];
            register double c10=c1[tt];  register double c11=c1[tt+1];   register double c12=c1[tt+2];
            register double c20=c1[ttt];  register double c21=c1[ttt+1];   register double c22=c1[ttt+2];
            for(k = 0; k < n; k += 3) {
                /* 2 by 2 mini matrix multiplication using registers*/
                int ta = i*n+k; int tta = ta+n; int tb = k*n+j; int ttb = tb+n; int ttta=tta+n; int tttb=ttb+n;
                register double a00 = a[ta]; register double a01 = a[tta]; register double b00 = b[tb]; register double b01 = b[tb+1];
                register double a02=a[ttta]; register double b02=b[tb+2];

                c00 += a00*b00 ; c01 += a00*b01 ; c02 += a00*b02 ; c10 += a01*b00 ;
                c11+=a01*b01; c12+=a01*b02; c20+=a02*b00; c21+=a02*b01; c22+=a02*b02;

                a00 = a[ta+1]; a01 = a[tta+1]; a02=a[ttta+1]; b00 = b[ttb]; b01 = b[ttb+1]; b02=b[ttb+2];

                c00 += a00*b00 ; c01 += a00*b01 ; c02 += a00*b02 ; c10 += a01*b00 ;
                c11+=a01*b01; c12+=a01*b02; c20+=a02*b00; c21+=a02*b01; c22+=a02*b02;

                a00 = a[ta+2]; a01 = a[tta+2]; a02=a[ttta+2]; b00 = b[tttb]; b01 = b[tttb+1]; b02=b[tttb+2];

                c00 += a00*b00 ; c01 += a00*b01 ; c02 += a00*b02 ; c10 += a01*b00 ;
                c11+=a01*b01; c12+=a01*b02; c20+=a02*b00; c21+=a02*b01; c22+=a02*b02;

             }
             c1[t] = c00;
             c1[t+1] = c01;
             c1[t+2]=c02;
             c1[tt] = c10;
             c1[tt+1] = c11;
             c1[tt+2]=c12;
             c1[ttt]=c20;
             c1[ttt+1]=c21;
             c1[ttt+2]=c22;
        }
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\ntime of execution for dgemm3 n=%d is %f",n,cpu_time);
        gflops=(2*(double)pow(n,3))/(cpu_time*pow(10,9));
        printf("\n gflops used are %.16f",gflops);
        //implementing dgemmo for comparison
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                for(k=0;k<n;k++)
                    c2[i*n+j]+=a[i*n+k]*b[k*n+j];
          clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\ntime of execution for dgemm0 n=%d is %f",n,cpu_time);
        gflops=(2*(double)pow(n,3))/(cpu_time*pow(10,9));
        printf("\n gflops used are %.16f",gflops);
                    //comparison of error of dgemm3 with dgemmo
        for(i=0;i<n*n;i++)
        {
           difference=(abs)(c1[i]-c2[i]);
           if(difference>error)
            error=difference;
        }

            printf("\n the error value for n=%d is %f ",n,error);
       //implementing dgemm1 for performance comparison
       clock_gettime(CLOCK_MONOTONIC, &cstart);
       for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            register double r=c2[i*n+j];
            for(k=0;k<n;k++)
                    r+=a[i*n+k]*b[k*n+j];
                    c2[i*n+j]=r;


        }
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
         printf("\ntime of execution for dgemm1 n=%d is %f",n,cpu_time);
        gflops=(2*(double)pow(n,3))/(cpu_time*pow(10,9));
        printf("\n gflops used are %.16f",gflops);

        //implementing dgemm2 for performance comparison
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        for(i=0;i<n;i=i+2)
    {
        for(j=0;j<n;j=j+2)
        {
            t=i*n+j;
            tt=t+n;
            c00 = c1[t];
            c01 = c1[t+1];
            c10 = c1[tt];
            c11 = c1[tt+1];
            for(k=0;k<n;k=k+2)
            {
                ta = i*n+k;
                tta = ta+n;
                tb = k*n+j;
                ttb = tb+n;
                a00 = a[ta];
                a01 = a[ta+1];
                a10 = a[tta];
                a11 = a[tta+1];
                b00 = b[tb];
                b01 = b[tb+1];
                b10 = b[ttb];
                b11 = b[ttb+1];
                c00 += a00*b00 + a01*b10;
                c01 += a00*b01 + a01*b11;
                c10 += a10*b00 + a11*b10;
                c11 += a10*b01 + a11*b11;
            }
            c1[t] = c00;
            c1[t+1] = c01;
            c1[tt] = c10;
            c1[tt+1] = c11;
        }
    }
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
         printf("\ntime of execution for dgemm2 n=%d is %f",n,cpu_time);
        gflops=(2*(double)pow(n,3))/(cpu_time*pow(10,9));
        printf("\n gflops used are %.16f",gflops);

    }
    return 0;
}
