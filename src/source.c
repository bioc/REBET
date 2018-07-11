/* History: Aug 03 2017 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}

void SIGMA2_cont(double*, int*, int*, double*, int*, double*, double*, double*);
void SIGMA2_binary(double*, int*, int*, double*, int*, double*, double*, double*);
void computeXtDX(double*, double*, int*, int*, double*);

static const R_CMethodDef callMethods[] = {
  {"SIGMA2_cont", (DL_FUNC)&SIGMA2_cont, 8},
  {"SIGMA2_binary", (DL_FUNC)&SIGMA2_binary, 8},
  {"computeXtDX", (DL_FUNC)&computeXtDX, 5},
  {NULL, NULL, 0}
};


/*
static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}

static void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %g ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}
*/


/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) Calloc(n, double);
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }
  
  return(ret);

} /* END: dVec_alloc */

/* Function to allocate an integer matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) Calloc(nrow, double *);
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) Free(x[i]);
  }
  Free(x);

} /* END: matrix_free */

/* Function to put a vector (stacked by column) into a matrix */
static void vecIntoMat(vec, nrow, ncol, ret)
double *vec;  /* length nrow*ncol stacked by column */
double **ret; /* dimension nrow x ncol */
int nrow, ncol;
{
  int i, j;
  double *pvec;

  pvec = vec;
  for (i=0; i<ncol; i++) {
    for (j=0; j<nrow; j++) ret[j][i] = *pvec++;
  }

  return;

} /* END: vecIntoMat */

/* Multiply matrices (stacked by columns) */
static void mult_stVec_stVec(X, invXtX, n, m, ret)
double *X, *invXtX, **ret;
int n, m;
{
  int i, j, k, index1, index2;
  double sum, *p2;

  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      sum    = 0.0;
      index1 = i;
      index2 = m*j;
      p2     = &invXtX[index2];
      for (k=0; k<m; k++) {
        sum += X[index1]*p2[k];
        index1 += n;
      }
      ret[i][j] = sum;
    }
  }

  return;

} /* END: mult_stVec_stVec */

/* Multiply (elementwise) a column of a matrix (stacked vector) by a scalar */
static void mult_stVec_scalar(X, n, col, c, ret)
double *X, c, *ret;
int n, col;
{
  int i;
  double *p, *pret;

  p = &X[n*col];
  for (i=0, pret=ret; i<n; i++, pret++) *pret = p[i]*c;

  return;

} /* END: mult_stVec_scalar */

/* Multiply vector by the transpose of a matrix(stacked by columns)
   and store result in a vector */
static void mult_vec_TstVec(vec, TstVec, n, m, ret)
double *vec;    /* Vector of length m */
double *TstVec; /* A stacked vector (of columns) for the transpose of an nxm matrix  */
double *ret;    /* Return vector of length n */
int n, m;     
{
  int i, k, index;
  double sum, *p1, *pret;

  for (i=0, pret=ret; i<n; i++, pret++) {
      sum   = 0.0;
      index = i;
      for (k=0, p1=vec; k<m; k++, p1++) {
        sum += *p1 * TstVec[index];
        index += n;
      }
      *pret = sum;
  }

  return;

} /* END: mult_vec_TstVec */

/* Dot product of 2 vectors */
static double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double sum=0.0, *p1, *p2;

  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) sum += *p1 * *p2;

  return(sum);

} /* END: dotProd */

/* Function to compute the diagonal elements of
  SIGMA <- Gt%*% (diag(1, n) - X %*% inv.XtX %*% t(X)) %*% G * tilde.sigma2
*/
void SIGMA2_cont(G, p_nsub, p_nsnp, X, p_ncov, invXtX, p_tildeSigma, ret)
double *G;            /* Stacked vector of genotype matrix (columns stacked) */
int *p_nsub;          /* Number of subjects */
int *p_nsnp;          /* Number of snps */
double *X;            /* Stacked vector of covariate matrix (columns stacked) */
int *p_ncov;          /* Number of covariates */
double *invXtX;       /* Stacked vector of XtX^-1 matrix (columns stacked) */
double *p_tildeSigma; /* Variance */
double *ret;          /* Return vector of length nsnp */
{
  int nsnp, nsub, ncov, index, row, i, cov1Flag;
  double **XinvXtX, *Gsigma2, *tvec2, sum, tildeSigma2, temp, *p;
  double *tGrow, invXtX1;

  nsnp        = *p_nsnp;
  ncov        = *p_ncov;
  nsub        = *p_nsub;
  tildeSigma2 = *p_tildeSigma;
  if (ncov < 2) {
    invXtX1  = *invXtX;
    cov1Flag = 1;
  } else {
    invXtX1  = 0.0;
    cov1Flag = 0;
  }

  /* Vectors for scratch space */
  Gsigma2 = dVec_alloc(nsub, 0, 0.0);
  tvec2   = dVec_alloc(nsub, 0, 0.0);

  /* Matrix to hold X*invXtX */
  if (!cov1Flag) {
    XinvXtX = dMat_alloc(nsub, ncov, 0, 0.0);
  
    /* Multiply X %*% invXtX */
    mult_stVec_stVec(X, invXtX, nsub, ncov, XinvXtX);
  } else {
    XinvXtX = NULL;
  }

  /* Loop over each element of return vector */
  for (index=0; index<nsnp; index++) {
    sum = 0.0;

    /* We need the (index,index) entry of the SIGMA matrix.
       First, compute column index of G*tilde.sigma2 */
    mult_stVec_scalar(G, nsub, index, tildeSigma2, Gsigma2);

    /* Get a pointer to the row index of G-transpose. This vector is in order. */ 
    tGrow = &G[index*nsub];

    /* Compute diag(1, n) - X %*% inv.XtX %*% t(X), row by row,
       and then multiply that row by G * tilde.sigma2 = Gsigma2 and
       store it it tvec2. This will give us column index of
       (diag(1, n) - X %*% inv.XtX %*% t(X)) %*% G * tilde.sigma2
    */
    for (row=0; row<nsub; row++) {
      /* Simple case when there is only an intercept in X (ncov=1) */
      if (cov1Flag) {
        for (i=0, p=tvec2; i<nsub; i++, p++) *p = -invXtX1;
        /* tvec2 is now a row of -(X %*% inv.XtX %*% t(X)) */
      } else {
        mult_vec_TstVec(XinvXtX[row], X, nsub, ncov, tvec2);
        /* tvec2 is now a row of X %*% inv.XtX %*% t(X) */

        /* Now negate the vector*/
        for (i=0, p=tvec2; i<nsub; i++, p++) *p = -*p;
      }

      /* Now add the row of diag(1, n), since we negated the vector above */
      tvec2[row] = 1.0 + tvec2[row];
      /* tvec2 is now a row of diag(1, n) - X %*% inv.XtX %*% t(X) */

      /* Multiple this row by the column Gsigma2, then multiply this
         scalar by the rowth element of row index of G-transpose */
      temp = dotProd(tvec2, Gsigma2, nsub)*tGrow[row];
      sum += temp;
    }
    ret[index] = sum;
  }

  if (!cov1Flag) matrix_free((void **)XinvXtX, nsub);
  Free(Gsigma2);
  Free(tvec2);

  return;

} /* END: SIGMA2_cont */

/* Function to compute the diagonal elements of
  SIGMA <- Gt %*% (D - (d %o% d) * (X %*% (inv.XtDX) %*% t(X))) %*% G
*/
void SIGMA2_binary(G, p_nsub, p_nsnp, X, p_ncov, invXtX, d, ret)
double *G;            /* Stacked vector of genotype matrix (columns stacked) */
int *p_nsub;          /* Number of subjects */
int *p_nsnp;          /* Number of snps */
double *X;            /* Stacked vector of covariate matrix (columns stacked) */
int *p_ncov;          /* Number of covariates */
double *invXtX;       /* Stacked vector of XtX^-1 matrix (columns stacked) */
double *d;            /* mu*(1-mu) */
double *ret;          /* Return vector of length nsnp */
{
  int nsnp, nsub, ncov, index, row, i, cov1Flag;
  double **XinvXtX, *Gcol, *tvec2, sum, temp, *p1, *p2;
  double *tGrow, invXtX1, drow;

  nsnp        = *p_nsnp;
  ncov        = *p_ncov;
  nsub        = *p_nsub;
  if (ncov < 2) {
    invXtX1  = *invXtX;
    cov1Flag = 1;
  } else {
    invXtX1  = 0.0;
    cov1Flag = 0;
  }

  /* Vector for scratch space */
  tvec2   = dVec_alloc(nsub, 0, 0.0);

  /* Matrix to hold X*invXtX */
  if (!cov1Flag) {
    XinvXtX = dMat_alloc(nsub, ncov, 0, 0.0);
  
    /* Multiply X %*% invXtX */
    mult_stVec_stVec(X, invXtX, nsub, ncov, XinvXtX);
  } else {
    XinvXtX = NULL;
  }

  /* Loop over each element of return vector */
  for (index=0; index<nsnp; index++) {
    sum = 0.0;

    /* We need the (index,index) entry of the SIGMA matrix.
       First, get column index of G */
    Gcol = &G[nsub*index];

    /* Get a pointer to the row index of G-transpose. This vector is in order. */ 
    tGrow = Gcol;

    /* Compute D - (d %o% d) * ( X %*% inv.XtX %*% t(X)), row by row,
       and then multiply that row by G  and
       store it it tvec2. This will give us column index of
       (D - (d %o% d) * (  X %*% inv.XtX %*% t(X)) %*% G 
    */
    for (row=0; row<nsub; row++) {
      drow = d[row];      

      /* Simple case when there is only an intercept in X (ncov=1) */
      if (cov1Flag) {
        temp = -drow*invXtX1;
        for (i=0, p1=d, p2=tvec2; i<nsub; i++, p1++, p2++) *p2 = temp* *p1;
        /* tvec2 is now a row of -(X %*% inv.XtX %*% t(X)) */
      } else {
        mult_vec_TstVec(XinvXtX[row], X, nsub, ncov, tvec2);
        /* tvec2 is now a row of X %*% inv.XtX %*% t(X) */

        /* Now multiple by row of (d %o% d) and negate */
        temp = -drow;
        for (i=0, p1=d, p2=tvec2; i<nsub; i++, p1++, p2++) *p2 = temp* *p1 * *p2;
      }

      /* Now add row of D = diag(d) */
      tvec2[row] = drow + tvec2[row];
      /* tvec2 is now a row of D - (d %o% d) * (X %*% (inv.XtDX) %*% t(X)) */

      /* Multiple this row by the column G, then multiply this
         scalar by the rowth element of row index of G-transpose */
      temp = dotProd(tvec2, Gcol, nsub)*tGrow[row];
      sum += temp;
    }
    ret[index] = sum;
  }

  if (!cov1Flag) matrix_free((void **)XinvXtX, nsub);
  Free(tvec2);

  return;

} /* END: SIGMA2_binary */

/* Function to compute XtDX <- t(X)%*%D%*%X, D=diag(d) */
void computeXtDX(X, d, p_nsub, p_ncov, ret)
double *X, *d, *ret; /* ret should be ncov*ncov */
int *p_nsub, *p_ncov;
{
  int i, j, k, nsub, ncov;
  double **Xmat, sum;

  nsub = *p_nsub;
  ncov = *p_ncov;

  /* Put X into a matrix */
  Xmat = dMat_alloc(nsub, ncov, 0, 0.0);
  vecIntoMat(X, nsub, ncov, Xmat);

  for (i=0; i<ncov; i++) {
    for (j=i; j<ncov; j++) {
      sum = 0.0;
      for (k=0; k<nsub; k++) sum += Xmat[k][i]*d[k]*Xmat[k][j];

      /* (i,j) element of matrix is j*ncov+i */
      ret[j*ncov+i] = sum;
      if (i != j) ret[i*ncov+j] = sum;
    }
  }

  matrix_free((void **)Xmat, nsub);

  return;

} /* END: computeXtDX */  

void R_init_REBET(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

