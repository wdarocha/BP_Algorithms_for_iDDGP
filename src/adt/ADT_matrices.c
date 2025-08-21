#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include "adt/ADT_math.h"
#include "adt/ADT_vectors.h"
#include "adt/ADT_matrices.h"

/* *********************************************************************************** */
/* ---------------------------------- ADT double matrices ---------------------------- */
/* *********************************************************************************** */
int **alloc_mat_d(int m, int n){

	int i;
	int **M = (int**)malloc(m*sizeof(int*));
	if (M == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	}
	for(i = 0; i < m; i++)
		M[i] = alloc_vec_d(n);
	
	return M;
}
/* *********************************************************************************** */
void dealloc_mat_d(int **M, int m){

	int i;
	for (i = 0; i < m; i++)
		free(M[i]);
	free(M);
}
/* *********************************************************************************** */
double **alloc_mat_lf(int m, int n){

	int i;
	double **M = (double**)malloc(m*sizeof(double*));
	if (M == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	}
	for(i = 0; i < m; i++)
		M[i] = alloc_vec_lf(n);
	
	return M;
}

/* *********************************************************************************** */
void dealloc_mat_lf(double **M, int m) {
	
	int i;
	
	for(i = 0; i < m; i++)
		free(M[i]);
	
	free(M);
}
/* *********************************************************************************** */
double **zeros_mat_lf(int m, int n){
 
	int i, j;
	
	double **M = alloc_mat_lf(m, n);
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M[i][j] = 0.0;
	return M;
}
/* *********************************************************************************** */
double ***alloc_3Darray_vec_lf(int n, int *vec){

	double ***T;
	int i;
	
	T = (double***)malloc(n*sizeof(double**));
	if (T == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	}
	for(i = 0; i < n; i++)
		T[i] = alloc_mat_lf(vec[i], 3);
	
	return T;
}
/* *********************************************************************************** */
void dealloc_3Darray_vec_lf(double ***T, int *vec_kv, int n){
	
	int k;
	for(k = 0; k < n; k++)
		dealloc_mat_lf(T[k], vec_kv[k]);
	free(T);
}
/* *********************************************************************************** */
double ***alloc_3Darray_lf(int l, int m, int n){
	
	double ***T;
	int i;

	T = (double***)malloc(l*sizeof(double**));
	if (T == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	}
	for(i = 0; i < l; i++)
		T[i] = alloc_mat_lf(m, n);
		
	return T;
}
/* *********************************************************************************** */
void dealloc_3Darray_lf(double ***T, int l, int m){
	
	int k;
	for(k = 0; k < l; k++)
		dealloc_mat_lf(T[k], m);
	free(T);
}
/* *********************************************************************************** */
double *get_column_mat_lf(double **M, int m, int c){
	
	double *col = alloc_vec_lf(m);
	int i;
	
	for(i = 0; i < m; i++)
		col[i] = M[i][c-1];
	
	return col;
}
/* *********************************************************************************** */
double **get_columns_mat_lf(double **M, int m, int *cols, int ncols){
	
	int k, i;
	double **A = alloc_mat_lf(m, ncols);
	
	for(k = 0; k < ncols; k++)
		for(i = 0; i < m; i++)
			A[i][k] = M[i][cols[k]-1];
		
	return A;
}
/* *********************************************************************************** */
int **mat_lf_2_mat_d(double **M, int m, int n){
	
	int i, j;
	int **A = alloc_mat_d(m, n);
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			A[i][j] = (int) M[i][j];
		
	return A;
}
/* *********************************************************************************** */
double **ones_mat_lf(int m, int n){
 
	int i, j;
	
	double **M = alloc_mat_lf(m, n);
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M[i][j] = 1.0;
	return M;
}
/* *********************************************************************************** */
void mat_e_mat_lf(double **M1, double **M0, int m, int n){

	int i, j;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M1[i][j] = M0[i][j];
	
}
/* *********************************************************************************** */
void mat_m_mat_lf(double **M, double **M1, double **M0, int m, int n){

	int i, j;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M[i][j] = M1[i][j] - M0[i][j];
	
}
/* *********************************************************************************** */
double frobenius_norm(double **A, int m, int n){

	int i, j;
	double fn = 0.0;
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			fn += A[i][j]*A[i][j];
	
	return sqrt(fn);	
}
/* *********************************************************************************** */
void l_t_mat_lf(double l, double **M, int m, int n){

	int i, j;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M[i][j] = l*M[i][j];
}
/* *********************************************************************************** */
double **trans_mat_lf(double **M, int m, int n){

	int i, j;
	double **Mt = alloc_mat_lf(n, m);
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			Mt[j][i] = M[i][j];
			
	return Mt;
}
/* *********************************************************************************** */
double **mat_times_mat_lf(double **M1, double **M2, int m, int p, int n){

	int i, j;
	double *cj;
	double **M = alloc_mat_lf(m, n);
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++){
			cj = get_column_mat_lf(M2, p, j+1);
			M[i][j] = dot_product_lf(&M1[i][0], cj, p);
			free(cj);
		}
		
	return M;
}
/* *********************************************************************************** */
void mat_2_vec_lf(double **M, int m, int n, double *v){

	int i, j, k;
	k = 0;
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			v[k++] = M[i][j];
	
}
/* *********************************************************************************** */
void vec_2_mat_lf(double **M, int m, int n, double *v){

	int i, j;
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			M[i][j] = v[i*n + j];
	
}
/* *********************************************************************************** */
void vec_2_diagmat_lf(double **M, double *v, int p){

	int i;
	for(i = 0; i < p; i++)
		M[i][i] = v[i];
	
}
/* *********************************************************************************** */
void svd_row_major_matrix(double **A, int m0, int n0, double **U, double **S, double **Vt){
	
	double *A_vec  = alloc_vec_lf(m0*n0);
	double *U_vec  = alloc_vec_lf(m0*m0);
	double *S_vec  = alloc_vec_lf(minAB_d(m0,n0));
	double *Vt_vec = alloc_vec_lf(n0*n0);
	
	mat_2_vec_lf(A, m0, n0, A_vec);
	
	// Singular value decomposition
	char jobu  = 'A'; // Compute all left singular vectors
	char jobvt = 'A'; // Compute all right singular vectors
	lapack_int m = m0; // Number of rows of matrix A
	lapack_int n = n0; // Number of columns of matrix A
	lapack_int lda  = n0; // In LAPACK_ROW_MAJOR, lda must be at least n.
	lapack_int ldu  = m0; // If jobu is 'A', ldu must be at least m.
	lapack_int ldvt = n0; // If jobvt is 'A', ldvt must be at least n.
	
	// Allocate space for the superdiagonal elements
	double *superb = alloc_vec_lf(minAB_d(m0,n0) - 1);
	
	// Compute SVD
	lapack_int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, A_vec, lda, S_vec, U_vec, ldu, Vt_vec, ldvt, superb);
	
	// Check for success
	if(info > 0){
		printf("The algorithm computing SVD failed to converge.\n");
	}
	else{
		vec_2_mat_lf(U, m0, m0, U_vec);
		vec_2_diagmat_lf(S, S_vec, minAB_d(m0, n0));
		vec_2_mat_lf(Vt, n0, n0, Vt_vec);
	}
	
	free(A_vec);
	free(U_vec);
	free(S_vec);
	free(Vt_vec);
	free(superb);
}
/* *********************************************************************************** */
/* ---------------------------------- ADT string matrices ---------------------------- */
/* *********************************************************************************** */
char **alloc_mat_s(int m, int n){

	int i;
	char **T = (char**)malloc(m*sizeof(char*));
	if (T == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	}
	for(i = 0; i < m; i++)
		T[i] = alloc_vec_s(n);
	
	return T;
}
/* *********************************************************************************** */
void dealloc_mat_s(char **T, int m){

	int i;
	for (i = 0; i < m; i++)
		free(T[i]);
	free(T);
}
