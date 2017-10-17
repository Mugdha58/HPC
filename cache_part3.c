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
    int i,j,k,n=2048,block_size=2,i1,j1,k1;
    double *a,*b,*c1,*c2;
    double cpu_time,gflops,temp,error,difference;
    register int t,tt,ta,tta,tb,ttb;
    register double c00,c01,c10,c11,a00,a01,a10,a11,b00,b01,b10,b11;
    clock_t start,end;

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

    for(block_size=4;block_size<20;block_size=block_size*2)
    {
        start= clock();
        printf("BLOCK %d calculation:",block_size);
        for(i=0;i<n;i+=block_size)
            for(j=0;j<n;j+=block_size)
                for(k=0;k<n;k+=block_size)
                {
                    for (i1 = i; i1 < i+block_size; i1+=2)
                    for (j1 = j; j1 < j+block_size; j1+=2){
                        t=i1*n+j1;
                        tt=t+n;
                        c00 = c1[t];
                        c01 = c1[t+1];
                        c10 = c1[tt];
                        c11 = c1[tt+1];
                    for (k1 = k; k1 < k+block_size; k1+=2){
                        ta = i1*n+k1;
                        tta = ta+n;
                        tb = k1*n+j1;
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


        end=clock();
        cpu_time=(double)(end-start)/(CLOCKS_PER_SEC);
        printf("\ntime of execution is %f",cpu_time);
        gflops=(2*(double)pow(n,3))/(cpu_time*pow(10,9));
        printf("\n gflops used are %f",gflops);
        //implementing dgemmo for comparison
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                for(k=0;k<n;k++)
                    c2[i*n+j]+=a[i*n+k]*b[k*n+j];
                    //comparison of dgemm2 with dgemmo
        error = 0.0;
        for(i=0;i<n*n;i++)
        {
           difference=(abs)(c1[i]-c2[i]);
           if(difference>error)
            error=difference;
        }

            printf("\n the error value for n=%d is %f ",n,error);
    }
    free(a);
    free(b);
    free(c1);
    free(c2);
    return 0;
}
