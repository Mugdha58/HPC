#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
int i,j,k,n,t,l,p,q,m;
double Random_gen ( )
{
    double upper_bound=RAND_MAX/10.0;
    return((double)rand()/upper_bound);

}
double neg_inverse(double *a,int n)
{
    for(i=0;i<n;i++)
    for(j=0;j<n;j++){
        if(i!=j)
            a[i*n+j]*=(-1);
    }
}
void mydgetrf(double *a,int *pvt,int n,int block)
{
    int maxind,temps,ib,end;
    double max,tempv;
    double *ll;
    for(ib=0;ib<n;ib+=block)
    {
     end=3+ib;
    for(i=ib;i<end;i++)
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
       ll = (double*)calloc(sizeof(double), block*block);
            p=0;q=0;
            for(l=ib;l<=end;l++){
                for(m=ib;m<=end;m++){
                    if(l<m){
                        ll[p*block+q] = a[p*block+q];
                    }
                    else if(l==m){
                        ll[p*block+q] = 1;
                    }
                    else{
                        ll[p*block+q] = 0;
                    }
                    q++;
                }
                p++;
                q=0;
            }
            //matInverse(ll,lln,b);
            p=0;q=0;
            for(j=ib;j<=end;j++){
                //arrA[j*n+i] = arrA[j*n+i]/arrA[i*n+i];
                for(k=end+1;k<n;k++){
                    for(m=ib;m<end;m++){
                        a[j*n+k] += ll[p*block+q] * a[m*n+k];
                        q++;
                    }
                    q=0;
                }
                p++;
                q=0;
            }
    }
      //factorizing
      for(j=end;j<n;j++)
      {
          for(k=end+1;k<n;k++){
            for(l=ib;i<=end;l++)
            a[j*n+k]=a[j*n+k]-(a[j*n+l]*a[l*n+k]);
          }
      }
    }
    //ll inverse

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
    int *pvt,n=8;
    double *a,*b,*abk,*x,*y;
    double gflops,cpu_time;
    struct timespec cstart = {0,0}, cend ={0,0};
   // for(n=1000;n<6000;n=n+1000)
    //{
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
    clock_gettime(CLOCK_MONOTONIC, &cstart);
    mygetrf(a,pvt,n,4);
    clock_gettime(CLOCK_MONOTONIC, &cend);
    cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
    printf("\nCPU time for LU factorization n=%d is %f",cpu_time);
    gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
    printf("\nthe gflops used are=%f",gflops);
    mydtrsm(n,a,b,pvt,x,y,0);
    mydtrsm(n,a,b,pvt,x,y,1); // label 1 is passed so that backward substitution will be done
     for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    printf("\nthe result of library function is %f\t",x[i*n+j]);
    char    TRANS = 'N';
    int     NRHS = 1;
    int     IPIV[3];
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
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
    // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
        UPLO = 'U';
        DIAG = 'N';
    // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&a,A, &n, B, &n);
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
            printf("\nthe result of library function is %f\t",B[i*n+j]);
    free(a);
    free(b);
    free(x);
    free(y);
    free(pvt);
    }
    return 0;
}
