
/*---------------------------------------------------------------------------*/
/**
  @file	linalgutils.h
  @brief simple functions for vector/matrix handling.
 **/
/*--------------------------------------------------------------------------*/

#ifndef _VECTUTILS_H_
#define _VECTUTILS_H_

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double *vzeros(int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=0;
  return res;
}

double *vones(int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=1;
    return res;
}

double *linspace(double a, double b, int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  double delta=(b-a)/(n-1);
  for (int i=0;i<n;i++)
    res[i]=a+i*delta;
  return res;
}

void vscale(int n, double *v, double alpha)
{
  for (int i=0;i<n;i++)
    v[i]=alpha*v[i];
}

void vprint(int n, double *v)
{
  printf("[");
  for (int i=0;i<n;i++)
    printf(" %2.2f ",v[i]);
  printf("]\n");
}

double *vinit(int n,double alpha)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=alpha;
  return res;
}

double **minit(int m, int n,double alpha)
{
  double ** res= (double **)malloc(m*sizeof(double *));
  double * data=(double *)malloc(m*n*sizeof(double));
  for (int i=0;i<m*n;i++)
    data[i]=alpha;
  for (int i=0;i<m;i++)
    res[i]=data+i*n;
  return res;
}

void mfree(double **a)
{
  free(&a[0][0]);
  free(a);
}

void mprint(int m, int n, double **a)
{
  for (int i=0;i<m;i++)
  {
  printf("[");
  for (int j=0;j<n;j++)
    printf(" %2.2f ",a[i][j]);
  printf("]\n");
  }
}

double vdot(int n , double *x,double *y)
{
  double res=0;
  for (int i=0;i<n;i++)
    res+=x[i]*y[i];
  return res;
}

void vaxpy(int n, double alpha , double *x, double *y)
{
  for (int i=0;i<n;i++)
    y[i]+=alpha*x[i];
}

void mger(int m, int n, double alpha, double **A, double *x, double * y)
{
  for (int i=0;i<m;i++)
      for (int j=0;j<n;j++)
        A[i][j]+=alpha*x[i]*y[j];
}

void mgemv(int m, int n, double alpha, double **A, double *x, double beta, double *y)
{
  for (int i=0;i<m;i++)
    {
      y[i]=beta*y[i];
      for (int j=0;j<n;j++)
        y[i]+=alpha*A[i][j]*x[j];
    }
}

void mgemm(int m, int n, int p, double alpha, double **A, double **B, double beta, double **C)
{
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
  {
    C[i][j]=beta*C[i][j];
    for (int k=0;k<p;k++)
      C[i][j]+=alpha*A[i][k]*B[k][j];
  }
}

void diagsv(int n,double **A,double *b,double *x)
{
  for (int i=0;i<n;i++)
    x[i]=b[i]/A[i][i];
}

void triinfsv(int n,double **L,double *b,double *x)
{
  x[0]=b[0]/L[0][0];
  for (int k=1;k<n;k++)
  {
    x[k]=b[k];
    for (int j=0;j<k;j++)
      x[k]-=L[k][j]*x[j];
    x[k]=x[k]/L[k][k];
  }
}

void trisupsv(int n,double **L,double *b,double *x)
{
  x[n-1]=b[n-1]/L[n-1][n-1];
  for (int k=n-2;k>=0;k--)
  {
    x[k]=b[k];
    for (int j=k+1;j<n;j++)
      x[k]-=L[k][j]*x[j];
    x[k]=x[k]/L[k][k];
  }
}


void pivot(int n, double **A, double *b)
{
  for (int k=0;k<n-1;k++)
  { // k : colonne a mettre en echelon
    for (int i=k+1;i<n;i++)
      { // i : ligne mise a jour
        for (int j=k+1;j<n;j++)
          A[i][j]-=A[k][j]*A[i][k]/A[k][k];
        // maj vec b
        b[i]-=A[i][k]*b[k]/A[k][k];
        A[i][k]=0;// pivot mis a 0
      }
  }
}


/*
  printf("pivot de A:\n");
printf("Matrice A:\n");
mprint(n,n,A);
printf("b=");vprint(n,b);
*/

double err(int n, double *v1, double *v2)
{
  double res=0;
  for (int i=0;i<n;i++)
    res+= (v1[i]-v2[i])*(v1[i]-v2[i]);
  return res;

}

double *vmap(double (*f)(double),double *x, int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=f(x[i]);
  return res;
}

void vmap2(double (*f)(double),double *x, int n, double * res)
{
  for (int i=0;i<n;i++)
    res[i]=f(x[i]);
}

double * vsub(double * a, double *b, int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=a[i]-b[i];
  return res;
}

void vsub2(double * a, double *b, int n, double *c)
{
  for (int i=0;i<n;i++)
    c[i]=a[i]-b[i];
}

double * vsum(double * a, double *b, int n)
{
  double * res= (double *)malloc(n*sizeof(double));
  for (int i=0;i<n;i++)
    res[i]=a[i]+b[i];
  return res;
}

void vsum2(double * a, double *b, int n, double *c)
{
  for (int i=0;i<n;i++)
    c[i]=a[i]+b[i];
}




#endif
