/** @file
	Memory allocation functions from Numerical Recipes in C, by Press,
	Cambridge University Press, 1988
*/
#ifndef FRAME_NRUTIL_H
#define FRAME_NRUTIL_H

typedef struct FCOMPLEX {float r,i;} fcomplex; /* */

void nrerror(char error_text[]);

float *vector(int nl, int nh);
int *ivector(long nl, long nh);
fcomplex *Cvector(long nl, long nh);
double *dvector(long nl, long nh);

float **matrix(int nrl, int nrh,int ncl, int nch);
float ***D3matrix(int nrl,int nrh, int ncl, int nch, int nzl, int nzh);
fcomplex **Cmatrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);

void free_vector(float *v,int nl, int nh);
void free_ivector(int *v,long nl, long nh);
void free_Cvector(fcomplex *v,long nl, long nh);
void free_dvector(double *v, int nl, int nh);

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl,int nrh,int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_D3matrix(float ***m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh);
void free_Cmatrix(fcomplex **m, int nrl, int nrh, int ncl, int nch);
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);

float **convert_matrix(float *a,int nrl, int nrh, int ncl, int nch);
void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);

void show_vector(int n, float *A);
void show_dvector(int n, double *A);
void show_matrix(int m, int n, float **A);
void show_dmatrix(int m, int n, double **A);

void save_vector(int n, float *V, char file[]);
void save_ivector(int n, int *V, char file[]);
void save_matrix (int m, int n, float **A, char file[]);
void save_dmatrix (int m, int n, double **A, char file[]);
void save_ut_matrix (int n, float **A, char file[]);
void save_ut_dmatrix (int n, double **A, char file[]);

#endif /* FRAME_NRUTIL_H */

