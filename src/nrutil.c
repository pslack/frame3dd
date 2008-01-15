/** @file
	Memory allocation functions from Numerical Recipes in C, by Press,
	Cambridge University Press, 1988
*/

#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrutil.h"

/* print error message to stderr	*/
void nrerror(char error_text[]){

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


/* allocate storage for a vector		*/
float *vector(int nl, int nh){
	float *v;

	v=(float *)malloc((unsigned) ((nh-nl+1)*sizeof(float)));
	if (!v) {
		fprintf(stderr," nl = %d  nh = %d \n", nl, nh );
		nrerror("allocation failure in vector()");
	}
	return v-nl;
}

int *ivector(nl,nh)	/* allocate storage for a vector		*/
long nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

fcomplex *Cvector(nl,nh)
long nl,nh;		/* allocate storage for a complex vector	*/
{
	fcomplex *v;

	v=(fcomplex *)malloc((unsigned) (nh-nl+1)*sizeof(fcomplex));
	if (!v) nrerror("allocation failure in Cvector()");
	return v-nl;
}

double *dvector(nl,nh)	/* allocate storage for a vector		*/
long nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

float **matrix(nrl,nrh,ncl,nch)	/* allocate storage for a matrix	*/
int nrl,nrh,ncl,nch;
{
	int	i;
	float 	**m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

float ***D3matrix(nrl,nrh,ncl,nch,nzl,nzh)  /* storage for a 3-D matrix */
int nrl,nrh,ncl,nch,nzl,nzh;
{
	int     i,j;
	float   ***m;

	m=(float ***) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("alloc failure 1 in 3Dmatrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(float **) malloc((unsigned) (nch-ncl+1)*sizeof(float*));
		if (!m[i]) nrerror("alloc failure 2 in 3Dmatrix()");
		m[i] -= ncl;
		for (j=ncl;j<=nch;j++) {
			m[i][j]=(float *)
				malloc((unsigned) (nzh-nzl+1)*sizeof(float));
			if (!m[i][j]) nrerror("alloc failure 3 in 3Dmatrix()");
			m[i][j] -= nzl;
		}
	}
	return m;
}

fcomplex **Cmatrix(nrl,nrh,ncl,nch)	
int nrl,nrh,ncl,nch;	/* allocate storage for a Complex matrix	*/
{
	int	i;
	fcomplex **m;

	m=(fcomplex **)malloc((unsigned) (nrh-nrl+1)*sizeof(fcomplex*));
	if (!m) nrerror("allocation failure 1 in Cmatrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(fcomplex *)malloc((unsigned)(nch-ncl+1)*sizeof(fcomplex));
		if (!m[i]) nrerror("allocation failure 2 in Cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)   /* allocate storage for a matrix	*/
int nrl,nrh,ncl,nch;
{
	int	i;
	double 	**m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)	/* allocate storage for a matrix	*/
int nrl,nrh,ncl,nch;
{
	int	i;
	int 	**m;

	m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}


void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v;
long nl,nh;
{
	free((char*) (v+nl));
}

void free_Cvector(v,nl,nh)
fcomplex *v;
long nl,nh;
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int	i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int	i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int	i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_D3matrix(m,nrl,nrh,ncl,nch,nzl,nzh)
float ***m;
int nrl,nrh,ncl,nch,nzl,nzh;
{
	int     i,j;

	for(i=nrh;i>=nrl;i--) {
		for(j=nch;j>=ncl;j--) {
			free((char*) (m[i][j]+nzl));
		}
	}
}

void free_Cmatrix(m,nrl,nrh,ncl,nch)
fcomplex **m;
int nrl,nrh,ncl,nch;
{
	int	i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
	int i,j,nrow,ncol;

	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;

	m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

/*---------------------------------------------------------------------------
SHOW_VECTOR  -  display a vector of dimension [1..n]
----------------------------------------------------------------------------*/
void show_vector ( n, A )
int     n;
float   *A;
{
	FILE    *fp_m;
	int     j;

	for (j=1; j <= n; j++) {
		if (A[j] != 0) fprintf(stdout,"%14.6e", A[j] );
		else		  fprintf(stdout,"   0       ");
	}
	fprintf(stdout," ];\n\n");
	return;
}


/*---------------------------------------------------------------------------
SHOW_DVECTOR  -  display a vector of dimension [1..n]
----------------------------------------------------------------------------*/
void show_dvector ( n, A )
int     n;
double *A;
{
	FILE    *fp_m;
	int     j;

	for (j=1; j <= n; j++) {
		if (A[j] != 0)	fprintf(stdout,"%14.6e", A[j] );
		else		fprintf(stdout,"   0       ");
	}
	fprintf(stdout," ];\n\n");
	return;
}


/*---------------------------------------------------------------------------
SHOW_MATRIX  -  display a matrix of dimension [1..m][1..n]
----------------------------------------------------------------------------*/
void show_matrix ( m,n, A )
int     n,m;
float   **A;
{
	FILE    *fp_m;
	int     i,j;

	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			if (A[i][j] != 0) fprintf(stdout,"%14.6e", A[i][j] );
			else		  fprintf(stdout,"   0       ");
		}
		if (i==m)	fprintf(stdout," ];\n\n");
		else		fprintf(stdout," \n");
	}
	return;
}


/*---------------------------------------------------------------------------
SHOW_DMATRIX  - display a matrix of dimension [1..m][1..n] 
----------------------------------------------------------------------------*/
void show_dmatrix ( m,n, A )
int     n,m;
double  **A;
{
	FILE    *fp_m;
	int     i,j;

	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			if (fabs(A[i][j]) > 1.e-99) fprintf(stdout,"%11.3e", A[i][j] );
			else		  fprintf(stdout,"   0       ");
		}
		if (i==m)	fprintf(stdout," ];\n\n");
		else		fprintf(stdout," \n");
	}
	return;
}


/*---------------------------------------------------------------------------
SAVE_VECTOR  -  save a vector of dimension [1..n] to the named file 
----------------------------------------------------------------------------*/
void save_vector( n, V, file )
int     n;
float   *V;
char    file[];
{
	FILE    *fp_v;
	int     i;
	void	exit();

	if ((fp_v = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= n; i++) {
		if (V[i] != 0)	fprintf(fp_v,"%15.6e", V[i] );
		else		fprintf(fp_v,"    0         ");         
		fprintf(fp_v,"\n");
	}
	fclose(fp_v);
	return;
}

/*---------------------------------------------------------------------------
SAVE_DVECTOR  -  save an integer vector of dimension [1..n] to the named file 
----------------------------------------------------------------------------*/
void save_ivector( n, V, file )
int     n, *V;
char    file[];
{
	FILE    *fp_v;
	int     i;
	void	exit();

	if ((fp_v = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= n; i++) {
		if (V[i] != 0)	fprintf(fp_v,"%14d", V[i] );
		else		fprintf(fp_v,"   0         ");         
		fprintf(fp_v,"\n");
	}
	fclose(fp_v);
	return;
}

/*---------------------------------------------------------------------------
SAVE_MATRIX  -  save a matrix of dimension [1..m][1..n] to the named file
----------------------------------------------------------------------------*/
void save_matrix ( m,n, A, file )
int     n,m;
float   **A;
char    file[];
{
	FILE    *fp_m;
	int     i,j;
	void	exit();

	if ((fp_m = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			if (A[i][j] != 0) fprintf(fp_m,"%15.6e", A[i][j] );
			else		  fprintf(fp_m,"    0          ");
		}
		fprintf(fp_m,"\n");
	}
	fclose ( fp_m);
	return;
}


/*---------------------------------------------------------------------------
SAVE_DMATRIX  - save a matrix of dimension [1..m][1..n] to the named file
----------------------------------------------------------------------------*/
void save_dmatrix ( m,n, A, file )
int     n,m;
double  **A;
char    file[];
{
	FILE    *fp_m;
	int     i,j;
	void	exit();

	if ((fp_m = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			if (fabs(A[i][j]) > 1.e-99) fprintf(fp_m,"%21.12e", A[i][j] );
			else		            fprintf(fp_m,"    0                ");
		}
		fprintf(fp_m,"\n");
	}
	fclose ( fp_m);
	return;
}


/*---------------------------------------------------------------------------
SAVE_UT_MATRIX  - 						     23apr01 
 save a symmetric matrix of dimension [1..n][1..n] to the named file 
 use only upper-triangular part
----------------------------------------------------------------------------*/
void save_ut_matrix ( n, A, file )
int     n;
float   **A;
char    file[];
{
	FILE    *fp_m;
	int     i,j;
        void	exit();

	if ((fp_m = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= n; i++) {
	  for (j=1; j <= n; j++) {
		if ( i > j ) {
			if (fabs(A[j][i]) > 1.e-99) fprintf(fp_m,"%15.6e", A[j][i] );
			else		            fprintf(fp_m,"    0          ");
		} else {
			if (fabs(A[i][j]) > 1.e-99) fprintf(fp_m,"%15.6e", A[i][j] );
			else		            fprintf(fp_m,"    0          ");
		}
	  }
	  fprintf(fp_m,"\n");
	}
	fclose ( fp_m);
	return;
}

/*---------------------------------------------------------------------------
SAVE_UT_DMATRIX  - 						23apr01
  save a symetric matrix of dimension [1..n][1..n] to the named file 
  use only upper-triangular part
----------------------------------------------------------------------------*/
void save_ut_dmatrix ( n, A, file )
int     n;
double  **A;
char    file[];
{
	FILE    *fp_m;
	int     i,j;
        void	exit();

	if ((fp_m = fopen (file, "w")) == NULL) {
		printf (" error: cannot open file: %s \n", file );
		exit(1);
	}
	for (i=1; i <= n; i++) {
	  for (j=1; j <= n; j++) {
		if ( i > j ) {
			if (fabs(A[j][i]) > 1.e-99) fprintf(fp_m,"%21.12e", A[j][i] );
			else		            fprintf(fp_m,"    0                ");
		} else {
			if (fabs(A[i][j]) > 1.e-99) fprintf(fp_m,"%21.12e", A[i][j] );
			else		            fprintf(fp_m,"    0                ");
		}
	  }
	  fprintf(fp_m,"\n");
	}
	fclose ( fp_m);
	return;
}
