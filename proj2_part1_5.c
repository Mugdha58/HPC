#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
int i,j,k,n,t;
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0;
    return((double)rand()/upper_bound);

}
// add transpose
void mydgetrf(double *a,int *pvt,int n)
{
    int maxind,temps;
    double max,tempv;
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
            for(k=i;k<n;i++)
            {tempv=a[i*n+k];
            a[i*n+k]=a[maxind*n+k];
            a[maxind*n+k]=tempv;
            }
        }
       }
      //factorizing
      for(j=i+1;j<n;k++)
      {
          a[j*n+i]=a[j*n+i]/a[i*n+i];
          for(k=i+1;k<n;k++)
            a[j*n+k]=a[j*n+k]-(a[j*n+i]*a[i*n+k]);
      }

    }
}
void mydtrsm(int n,double *a,double *b,int *pvt,double *x,double *y,int label)
{
    double sum=0.0,temp;
    if(label==0)// passing label to call forward and backward substitution separately
    {//forward substitution
    y[0]=b[pvt[0]];
    for(i=1;i<n;i++)
    {
      for(k=0;k<i;k++)
      {
        sum+=y[k]*a[i*n+k];
        y[i]=b[pvt[i]]-sum;
      }
    }
    }
    //backward substitution
    else
    {
    x[n-1]=y[n-1]/a[(n-1)*n+(n-1)];
    for(i=n-1;i>=0;i--)
        for(k=i+1;k<n;k++)
    {
       sum+= x[k]*a[i*n+k];
       temp=y[i]-sum;
       x[i]=temp/a[i*n+i];

    }
    }

}

int main()
{
    int *pvt,n;
    double *a,*b,*abk,*x,*y;
    double gflops,cpu_time;
    struct timespec cstart = {0,0}, cend ={0,0};
    for(n=1000;n<6000;n=n+1000)
    {
    a=(double *) calloc(sizeof(double), n*n);
    b=(double *) calloc(sizeof(double), n*1);
    pvt=(int *) calloc(sizeof(int), 1*n);
    y=(double *) calloc(sizeof(double), n*n);
    x=(double *) calloc(sizeof(double), n*n);
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
    {
              a[i*n+j]=(double)Random_gen();
              b[j]=(double)Random_gen();
    }
    for(i=0;i<n;i++)
        pvt[i]=i;
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    mydgetrf(a,pvt,n);
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for LU factorization n=%d is %f",cpu_time);
    gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    mydtrsm(n,a,b,pvt,x,y,0);
    mydtrsm(n,a,b,pvt,x,y,1); // label 1 is passed so that backward substitution will be done
    free(a);
    free(b);
    free(x);
    free(y);
    free(pvt);
    }
    return 0;
}
