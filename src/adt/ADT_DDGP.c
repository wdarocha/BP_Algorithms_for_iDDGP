#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "adt/ADT_math.h"
#include "adt/ADT_vectors.h"
#include "adt/ADT_intervals.h"
#include "adt/ADT_matrices.h"
#include "adt/ADT_strings.h"
#include "adt/ADT_DDGP.h"
#define PI 3.14159265359
#define MYZERO 0.000001

/* *********************************************************************************** */
/**
 * @brief Computes the cardinality the set of adjacent predecessors (Ui) for all vertices.
 *
 * This function receives a list of m vertices (with values ranging from 1 to n) and counts 
 * how many times each vertex appears. It returns an array of size n containing the
 * cardinality for each vertex (|Ui|).
 *
 * @param vertices   An array of m vertex identifiers (assumed to be 1-based indices).
 * @param m		  Number of elements in the vertices array.
 * @param n		  Total number of distinct vertices.
 * @return int*	  An array of size n, where cardinality[i] holds the count for vertex (i+1).
 */
int *adjacent_predecessors_cardinality(int *vertices, int m, int n) {
	
	int i;
	int *cardinality = zeros_vec_d(n);
	
	for(i = 0; i < m; i++)
		cardinality[vertices[i] - 1]++;
	
	return cardinality;
}
/* *********************************************************************************** */
/**
 * @brief Converts a distance constraints input vector to a full instance representation.
 *
 * This function transforms the dc_vec array into:
 *   - instanceAdjlist: adjacency list containing distances (as [v_j, d_l, d_u])
 *   - protein: an array of proteinstructure, initialized with atom/residue info
 *   - kv: array with the number of predecessors for each vertex
 *   - n: total number of distinct vertices
 *
 * @param dc_vec		Input array with atomic and distance information
 * @param instanceAdjlist	Output pointer to a 3D array: [n][kv[i]][3]
 * @param protein		Output pointer to an array of proteinstructure[n]
 * @param kv			Output pointer to an array of predecessor counts per vertex
 * @param m			Number of elements in dc_vec
 * @param n			Output: number of unique vertices (updated in-place)
 */
void dcinputfile_2_instance(dcinputfile *dc_vec, double ***distanceConstraints, double ****instanceAdjlist, proteinstructure **protein, int **kv, int m, int *n) {
	
	int k, i;
	int cardUk, *adjpred = NULL, *pos_adjpred_vk = NULL;
	
	// Allocate the distance constraints matrix: distanceConstraints[m][4]
	(*distanceConstraints) = alloc_mat_lf(m, 4);
	// Allocate array to store all 'i' values from dc_vec
	int *vertices = alloc_vec_d(m);
	for (k = 0; k < m; k++) {
		vertices[k] = dc_vec[k].i;
		
		(*distanceConstraints)[k][0] = dc_vec[k].i;
		(*distanceConstraints)[k][1] = dc_vec[k].j;
		(*distanceConstraints)[k][2] = dc_vec[k].dl;
		(*distanceConstraints)[k][3] = dc_vec[k].du;
	}
	
	// Determine the cardinality of V
	(*n) = max_vec_d(vertices, m);
		
	// Compute the cardinality of the adjacent predecessors set for each vertex
	(*kv) = adjacent_predecessors_cardinality(vertices, m, *n);
			
	// Allocate 3D adjacency list: instanceAdjlist[n][kv[i]][3]
	*instanceAdjlist = alloc_3Darray_vec_lf((*n), (*kv));
	
	// Allocate the protein structure array
	(*protein) = alloc_vec_proteinstructure((*n));
	
	// Initialize the first atom directly from the first element of dc_vec
	(*protein)[0].v_k = dc_vec[0].j;
	(*protein)[0].res_seq = dc_vec[0].rj;
	strcpy((*protein)[0].atom_name, dc_vec[0].aj);
	strcpy((*protein)[0].res_name, dc_vec[0].rtj);
	(*protein)[0].x = 0.0;
	(*protein)[0].y = 0.0;
	(*protein)[0].z = 0.0;
	
	// For each vertex k > 0, find its adjacent predecessors and fill structures
	for (k = 1; k < (*n); k++) {
		// describe the set Uk
		adjpred = find_val_vec_d(vertices, m, k + 1);
		// compute |Uk|
		cardUk = sum_vec_d(adjpred, m);
		// find the indices of each element of Uk
		pos_adjpred_vk = find_ones_position_vec_d(adjpred, m, cardUk);
		// Fill metadata for the atom k
		(*protein)[k].v_k = dc_vec[pos_adjpred_vk[0]].i;
		(*protein)[k].res_seq = dc_vec[pos_adjpred_vk[0]].ri;
		strcpy((*protein)[k].atom_name, dc_vec[pos_adjpred_vk[0]].ai);
		strcpy((*protein)[k].res_name, dc_vec[pos_adjpred_vk[0]].rti);
		(*protein)[k].x = 0.0;
		(*protein)[k].y = 0.0;
		(*protein)[k].z = 0.0;
		// Fill adjacency list for current vertex k
		for (i = 0; i < cardUk; i++) {
			(*instanceAdjlist)[k][i][0] = dc_vec[pos_adjpred_vk[i]].j;   // predecessor ID
			(*instanceAdjlist)[k][i][1] = dc_vec[pos_adjpred_vk[i]].dl;  // lower bound
			(*instanceAdjlist)[k][i][2] = dc_vec[pos_adjpred_vk[i]].du;  // upper bound
		}
		
		free(adjpred);
		free(pos_adjpred_vk);
	}
	free(vertices);
}
/* *********************************************************************************** */
/**
 * @brief Returns the lower or upper bound of the distance between vertices i and j.
 *
 * The function searches the adjacency list of vertex i for a connection to j.
 * If not found, it looks in the list of vertex j for a connection to i.
 * The value returned depends on the 'l_or_u' parameter:
 *   - "l" returns the lower bound
 *   - "u" returns the upper bound
 *
 * @param adjlist	A 3D array [n][kv[i]][3] where kv[i] = |Ui| and [v_j, d_l, d_u] = [neighbor, lower_bound, upper_bound].
 * @param cardUi	Array with where each element is |Ui|.
 * @param i		Index of vertex i (1-based).
 * @param j		Index of vertex j (1-based).
 * @param l_or_u	"l" to get lower bound, anything else gets upper bound.
 * @return double	The distance bound between i and j.
 */
double dij_from_adjlist(double ***adjlist, int *cardUi, int i, int j, char l_or_u[]) {
	
	// 1 = lower bound, 2 = upper bound
	int pos = (!strcmp(l_or_u, "l")) ? 1 : 2;
	double dij;
	// get the set Ui
	double *Ui = get_column_mat_lf(adjlist[i - 1], cardUi[i - 1], 1);
	// find j in Ui
	int *vec = find_val_vec_lf(Ui, cardUi[i - 1], (double) j);
	int len = sum_vec_d(vec, cardUi[i - 1]);
	
	if (len > 0) {
		// j is a adjacente predecessor of i
		int *pos_dij = find_ones_position_vec_d(vec, cardUi[i - 1], len);
		dij = adjlist[i - 1][pos_dij[0]][pos];
		free(pos_dij);
		free(vec);
		free(Ui);
	} 
	else {
		// i is possibly a adjacente predecessor of j
		free(Ui);
		free(vec);
		// get the set Uj
		double *Uj = get_column_mat_lf(adjlist[j - 1], cardUi[j - 1], 1);
		// find i in Uj
		vec = find_val_vec_lf(Uj, cardUi[j - 1], (double) i);
		len = sum_vec_d(vec, cardUi[j - 1]);
		free(Uj);
		
		if (len == 0) {
			fprintf(stderr, "Error: no distance found between vertices %d and %d.\n", i, j);
			free(vec);
			exit(EXIT_FAILURE);
		}
		
		int *pos_dji = find_ones_position_vec_d(vec, cardUi[j - 1], len);
		dij = adjlist[j - 1][pos_dji[0]][pos];
		
		free(vec);
		free(pos_dji);
	}
	return dij;
}
/* *********************************************************************************** */
/**
 * @brief Gets information from an adjacency list to build a structured array of discretization edges.
 *
 * This function fills the array discretizationEdges_2 with the triplets (j, d^2_l, d^2_u),
 * representing the squared lower and upper distance bounds between a vertex ik+1 and its
 * three discretization clique predecessors, as specified in the cliques array.
 *
 * The format of each entry is:
 *   discretizationEdges_2[k][i] = [v_j, d^2_l, d^2_u], for i = 0, 1, 2
 *
 * @param instanceAdjlist		Adjacency list of the instance [n][kv[i]][3]
 * @param kv					  |Uv|
 * @param n					   |V|
 * @param cliques				 Array of discretization cliques [n][4], where each line contains:
 *								(v_k, v_{k1}, v_{k2}, v_{k3}) for vertex k
 * @param discretizationEdges_2  Output: structured array [n][3][3] storing (j, d^2_l, d^2_u) for each vertex
 */
void adjlist_2_discretization_edges_2(double ***instanceAdjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2) {
	
	int k, i, vk;
	double dl, du;

	// Allocate output array: [n][3][3] where each entry is [j, d^2_l, d^2_u]
	*discretizationEdges_2 = alloc_3Darray_lf(n, 3, 3);

	// Manually set edges for vertex 2 (index 1): connected to vertex 1
	(*discretizationEdges_2)[1][0][0] = cliques[1][1];
	dl = dij_from_adjlist(instanceAdjlist, kv, 2, cliques[1][1], "l");
	du = dij_from_adjlist(instanceAdjlist, kv, 2, cliques[1][1], "u");
	(*discretizationEdges_2)[1][0][1] = dl * dl;
	(*discretizationEdges_2)[1][0][2] = du * du;
	// Manually set edges for vertex 3 (index 2): connected to vertex 1
	(*discretizationEdges_2)[2][0][0] = cliques[2][1];
	dl = dij_from_adjlist(instanceAdjlist, kv, 3, cliques[2][1], "l");
	du = dij_from_adjlist(instanceAdjlist, kv, 3, cliques[2][1], "u");
	(*discretizationEdges_2)[2][0][1] = dl * dl;
	(*discretizationEdges_2)[2][0][2] = du * du;
	// Manually set edges for vertex 3 (index 2): connected to vertex 2
	(*discretizationEdges_2)[2][1][0] = cliques[2][2];
	dl = dij_from_adjlist(instanceAdjlist, kv, 3, cliques[2][2], "l");
	du = dij_from_adjlist(instanceAdjlist, kv, 3, cliques[2][2], "u");
	(*discretizationEdges_2)[2][1][1] = dl * dl;
	(*discretizationEdges_2)[2][1][2] = du * du;

	// For all k >= 4, retrieve the three predecessors from the clique and get squared bounds
	for (k = 3; k < n; k++)
		for (i = 0; i < 3; i++) {
			vk = cliques[k][i + 1]; // clique vertices at columns 1, 2, 3
			(*discretizationEdges_2)[k][i][0] = vk;
			dl = dij_from_adjlist(instanceAdjlist, kv, k + 1, vk, "l");
			du = dij_from_adjlist(instanceAdjlist, kv, k + 1, vk, "u");
			(*discretizationEdges_2)[k][i][1] = dl * dl;
			(*discretizationEdges_2)[k][i][2] = du * du;
		}
}
/* *********************************************************************************** */
void adjlist_2_prune_edges_2(double ***instanceAdjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges_2) {
	
	int i, j;
	double *Ui, dl, du;
	int *vec1i1, *vec1i2, *vec1i3, *vec1i1i2i3, numPd, numId, numPd_i, numId_i;
	
	(*pruneEdges_2) = alloc_vec_pruneedgesset(n);
	
	(*pruneEdges_2)[0].cardUkP = 0;
	(*pruneEdges_2)[0].cardUkI = 0;
	(*pruneEdges_2)[1].cardUkP = 0;
	(*pruneEdges_2)[1].cardUkI = 0;
	(*pruneEdges_2)[2].cardUkP = 0;
	(*pruneEdges_2)[2].cardUkI = 0;
	
	for(i = 3; i < n; i++)	
		if(kv[i] > 3){
			numPd = 0;
			numId = 0;
			Ui = get_column_mat_lf(instanceAdjlist[i], kv[i], 1);
			vec1i1 = find_val_vec_lf(Ui, kv[i], cliques[i][1]);
			vec1i2 = find_val_vec_lf(Ui, kv[i], cliques[i][2]);
			vec1i3 = find_val_vec_lf(Ui, kv[i], cliques[i][3]);
			vec1i1i2i3 = alloc_vec_d(kv[i]);
			vec_p_vec_d(vec1i1i2i3, vec1i1, vec1i2, kv[i]);
			vec_p_vec_d(vec1i1i2i3, vec1i1i2i3, vec1i3, kv[i]);
			free(vec1i1);
			free(vec1i2);
			free(vec1i3);
			
			for(j = 0; j < kv[i]; j++)
				if(!vec1i1i2i3[j]){
					if(fabs(instanceAdjlist[i][j][1] - instanceAdjlist[i][j][2]) < MYZERO)
						numPd++;
					else
						numId++;
				}
			
			(*pruneEdges_2)[i].cardUkP = numPd;
			(*pruneEdges_2)[i].cardUkI = numId;
			if(numPd > 0)
				(*pruneEdges_2)[i].precise = alloc_mat_lf(numPd, 3);
			if(numId > 0)
				(*pruneEdges_2)[i].interval = alloc_mat_lf(numId, 3);
			
			numPd_i = 0;
			numId_i = 0;
			
			for(j = 0; j < kv[i]; j++)
				if(!vec1i1i2i3[j]){
					dl = dij_from_adjlist(instanceAdjlist, kv, i+1, Ui[j], "l");
					du = dij_from_adjlist(instanceAdjlist, kv, i+1, Ui[j], "u");
					if(fabs(dl - du) < MYZERO){
						(*pruneEdges_2)[i].precise[numPd_i][0] = instanceAdjlist[i][j][0];
						(*pruneEdges_2)[i].precise[numPd_i][1] = instanceAdjlist[i][j][1]*instanceAdjlist[i][j][1];
						(*pruneEdges_2)[i].precise[numPd_i][2] = instanceAdjlist[i][j][2]*instanceAdjlist[i][j][2];
						numPd_i++;
					}
					else{
						(*pruneEdges_2)[i].interval[numId_i][0] = instanceAdjlist[i][j][0];
						(*pruneEdges_2)[i].interval[numId_i][1] = instanceAdjlist[i][j][1]*instanceAdjlist[i][j][1];
						(*pruneEdges_2)[i].interval[numId_i][2] = instanceAdjlist[i][j][2]*instanceAdjlist[i][j][2];
						numId_i++;
					}
				}
			free(Ui);
			free(vec1i1i2i3);
			
			double **auxMat = alloc_mat_lf(numId_i, 3);
			mat_e_mat_lf(auxMat, (*pruneEdges_2)[i].interval, numId_i, 3);
			
			double *vecW = alloc_vec_lf(numId_i);
			for(j = 0; j < numId_i; j++)
				vecW[j] = (*pruneEdges_2)[i].interval[j][2] - (*pruneEdges_2)[i].interval[j][1];
				
			int *vecIndices = alloc_vec_d(numId_i);
			for(j = 0; j < numId_i; j++)
				vecIndices[j] = j;
				
			quicksort(vecW, vecIndices, 0, numId_i - 1);
			
				for(j = 0; j < numId_i; j++){
					(*pruneEdges_2)[i].interval[j][0] = auxMat[vecIndices[j]][0];
					(*pruneEdges_2)[i].interval[j][1] = auxMat[vecIndices[j]][1];
					(*pruneEdges_2)[i].interval[j][2] = auxMat[vecIndices[j]][2];
				}
			
			dealloc_mat_lf(auxMat, numId_i);
			free(vecIndices);
			free(vecW);
		}
		else{
			(*pruneEdges_2)[i].cardUkP = 0;
			(*pruneEdges_2)[i].cardUkI = 0;
		}
}
/* *********************************************************************************** */
/**
 * @brief Constructs discretization and pruning structures from an adjacency list.
 *
 * This function initializes:
 *   - discretizationEdges_2: a structured array with squared lower/upper bounds
 *   - pruneEdges_2: a list of edges used for pruning (non-discretization constraints)
 *
 * Both structures are derived from the adjacency list and clique information.
 *
 * @param instanceAdjlist	   3D array [n][kv[i]][3] representing distance constraints.
 * @param kv					 |Uv|.
 * @param n					  |V|.
 * @param cliques				Discretization cliques [n][4] (vk, vk1, vk2, vk3).
 * @param discretizationEdges_2 Output: structured 3D array [n][3][3] for discretization edges.
 * @param pruneEdges_2		  Output: array of pruning edge structures.
 */
void adjlist_2_graph_parameters(double ***instanceAdjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2, prune_edges_set **pruneEdges_2) {
	
	adjlist_2_discretization_edges_2(instanceAdjlist, kv, n, cliques, discretizationEdges_2);
	adjlist_2_prune_edges_2(instanceAdjlist, kv, n, cliques, pruneEdges_2);
}
/* *********************************************************************************** */
/**
 * @brief Computes the parameters A[i], B[i], C[i] that define the local 3D circle for vertex i.
 *
 * Given three previously placed points xi1, xi2, xi3 (in R^3), this function constructs:
 *   - A[i] = center of the circle
 *   - B[i], C[i] = orthogonal vectors in the plane of the circle (defining its orientation)
 *
 * @param i		Index of the current vertex (0-based).
 * @param A, B, C	Output arrays of shape [n][3], modified in-place.
 * @param xi3, xi2, xi1	Points x3, x2, x1 as 3D coordinates.
 * @param lambda_210	Scalar lambda_i = projection along ê.
 * @param rho2_210	rho^2_i = squared radius of the circle.
 * @param di1i2		Euclidean distance between x1 and x2.
 */
void circunference_parameters(int i, double **A, double **B, double **C, double *xi3, double *xi2, double *xi1, double lambda_210, double rho2_210, double di1i2) {
	
	// Temporary 3D vectors
	double aux[3], v[3], ehat[3], yhat[3], zhat[3];
	
	// Compute ê = normalize(xi1 - xi2)
	vec_m_vec_lf3(aux, xi1, xi2);
	l_t_vec_lf3(ehat, 1.0 / di1i2, aux);
	
	// Compute v = xi3 - xi2
	vec_m_vec_lf3(v, xi3, xi2);
	
	// Compute ẑ = normalize(ê × v)
	cross_product_lf3(aux, ehat, v);
	l_t_vec_lf3(zhat, 1.0 / norm2_lf3(aux), aux);
	
	// Compute ŷ = normalize(ẑ × ê)
	cross_product_lf3(aux, zhat, ehat);
	l_t_vec_lf3(yhat, 1.0 / norm2_lf3(aux), aux);
	
	// Compute A[i] = xi2 + lambda_i ê
	l_t_vec_lf3(aux, lambda_210, ehat);
	vec_p_vec_lf3(&A[i][0], xi2, aux);
	
	// Compute C[i] = rho_i ẑ
	l_t_vec_lf3(&C[i][0], sqrt(rho2_210), zhat);
	
	// Compute B[i] = rho_i ŷ = (C × ê)
	// Option 1: Use orthonormal ŷ
	// l_t_vec_lf3(&(*B)[i][0], sqrt(rho2_210), yhat);
	// Option 2: Ensure orthogonality using cross(C, ê)
	cross_product_lf3(&B[i][0], &C[i][0], ehat);
}
/* *********************************************************************************** */
/**
 * @brief Embeds vertex i into R^3 using angular parametrization on a circle.
 *
 * Computes:
 *	 xi = A + B * cos(tau) + C * sin(tau)
 *
 * Where:
 *   - A is the center of the circle in R^3
 *   - B and C are orthogonal vectors spanning the plane of the circle
 *   - tau is the torsion (sampling) angle, given by its cosine and the sign of the sine.
 *
 * @param xi	Output: coordinate of the vertex in R^3 (modified in-place).
 * @param A	Vector A[3]: center of the circle.
 * @param B	Vector B[3]: first spanning vector.
 * @param C	Vector C[3]: second spanning vector.
 * @param tau	Torsion angle tau.
 */
void vertex_embedding(double *xi, double *A, double *B, double *C, double tau) {
	
	// Temporary vectors in R^3
	double vec_aux1[3], vec_aux2[3], vec_aux3[3];
	
	// vec_aux1 = B * cos(tau)
	l_t_vec_lf3(vec_aux1, cos(tau), B);
	
	// vec_aux2 = C * sin(tau)
	l_t_vec_lf3(vec_aux2, sin(tau), C);
	
	// vec_aux3 = vec_aux1 + vec_aux2
	vec_p_vec_lf3(vec_aux3, vec_aux1, vec_aux2);
	
	// xi = A + vec_aux3
	vec_p_vec_lf3(xi, A, vec_aux3);
}
/* *********************************************************************************** */
void compute_cosTauii3_parameters(int i, double **X, double ***discretizationEdges_2, int *i1, int *i2, int *i3, double *di1i2, double *lambda_210, double *rho2_210, double *lambda_213, double *rho2_213){
	
	double di1i3_2, di2i3_2, dii1_2, dii2_2;
	
	(*i1) = (int) (discretizationEdges_2[i][0][0] - 1.0);
	(*i2) = (int) (discretizationEdges_2[i][1][0] - 1.0);
	(*i3) = (int) (discretizationEdges_2[i][2][0] - 1.0);
	
	di1i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i1)][0]);
	(*di1i2) =  dxixj_lf3(&X[(*i2)][0], &X[(*i1)][0]);
	di2i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i2)][0]);
	
	dii1_2 = discretizationEdges_2[i][0][1];
	dii2_2 = discretizationEdges_2[i][1][1];
	
	(*lambda_210) = lambdai_d2(dii2_2, (*di1i2), dii1_2);
	(*rho2_210)   = rhoi2_d2(dii2_2, (*lambda_210));
	
	(*lambda_213) = lambdai_d2(di2i3_2, (*di1i2), di1i3_2);
	(*rho2_213)   = rhoi2_d2(di2i3_2, (*lambda_213));
}
/* *********************************************************************************** */
void compute_tauii3_interval(int i, double **X, double ***discretizationEdges_2, int *i1, int *i2, int *i3, double *di1i2, double *lambda_210, double *rho2_210, double *lambda_213, double *rho2_213, type_matrix_interval *tauii3) {
	
	(*i1) = (int) (discretizationEdges_2[i][0][0] - 1.0);
	(*i2) = (int) (discretizationEdges_2[i][1][0] - 1.0);
	(*i3) = (int) (discretizationEdges_2[i][2][0] - 1.0);
	
	double di1i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i1)][0]);
	       (*di1i2) =  dxixj_lf3(&X[(*i2)][0], &X[(*i1)][0]);
	double di2i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i2)][0]);
	
	double dii1_2   = discretizationEdges_2[i][0][1];
	double dii2_2   = discretizationEdges_2[i][1][1];
	double dii3_2_l = discretizationEdges_2[i][2][1];
	double dii3_2_u = discretizationEdges_2[i][2][2];
	
	(*lambda_210) = lambdai_d2(dii2_2, (*di1i2), dii1_2);
	(*rho2_210)   = rhoi2_d2(dii2_2, (*lambda_210));
	
	(*lambda_213) = lambdai_d2(di2i3_2, (*di1i2), di1i3_2);
	(*rho2_213)   = rhoi2_d2(di2i3_2, (*lambda_213));
	
	double lbdmlbdb = (*lambda_210) - (*lambda_213);
	double p30 = lbdmlbdb * lbdmlbdb + (*rho2_210) + (*rho2_213);
	double q30 = sqrt((*rho2_210) * (*rho2_213));
	
	// Case: dii3_2_l == dii3_2_u (degenerate interval)
	if (fabs(dii3_2_l - dii3_2_u) < MYZERO) {
		double tau_l = acos(cos_torsion_angle_with_constants_d2(p30, q30, dii3_2_l));
		
		build_simple_symmetric_interval(tauii3, tau_l, tau_l);
		
		if ((fabs(tau_l) < MYZERO) || (fabs(tau_l - PI) < MYZERO))
			tauii3->nneg_intervals = 0;
	}
	else { // Case: dii3_2_l < dii3_2_u (normal interval)
		double tau_l = acos(cos_torsion_angle_with_constants_d2(p30, q30, dii3_2_l));
		double tau_u = acos(cos_torsion_angle_with_constants_d2(p30, q30, dii3_2_u));
		
		build_simple_symmetric_interval(tauii3, tau_l, tau_u);
	}
}
/* *********************************************************************************** */
void compute_mde_and_lde(double **X, double **distanceConstraints, int m, double *mde, double *lde){
	
	int i, j, k;
	double dij, dijL, dijU;
	double *deVec = zeros_vec_lf(m);
	
	(*lde) = 0.0;
	(*mde) = 0.0;
	
	for (k = 0; k < m; k++) {
		i = (int) distanceConstraints[k][0] - 1;
		j = (int) distanceConstraints[k][1] - 1;
		dij = dxixj_lf3(&X[i][0], &X[j][0]);
		
		dijL = distanceConstraints[k][2];
		dijU = distanceConstraints[k][3];
				
		deVec[k] = maxAB_lf(0.0, maxAB_lf(dijL - dij, dij - dijU));
	}
	
	(*lde) = max_val_vec_lf(deVec, m);
	(*mde) = sum_vec_lf(deVec, m)/m;
	
	free(deVec);
}
/* *********************************************************************************** */
double compute_rmsd(double **X, double **Y, int n){ // X = X0, Y = Xr;
	
	double **X0 = alloc_mat_lf(n, 3);
	mat_e_mat_lf(X0, X, n, 3);
	
	double **Xr = alloc_mat_lf(n, 3);
	mat_e_mat_lf(Xr, Y, n, 3);
	
	double cm_x_Xr, cm_y_Xr, cm_z_Xr, cm_x_X0, cm_y_X0, cm_z_X0;
	double *vec_aux;
	
	// compute the CM of Xr
	vec_aux = get_column_mat_lf(Xr, n, 1);
	cm_x_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(Xr, n, 2);
	cm_y_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(Xr, n, 3);
	cm_z_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	// compute the CM of X0
	vec_aux = get_column_mat_lf(X0, n, 1);
	cm_x_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(X0, n, 2);
	cm_y_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(X0, n, 3);
	cm_z_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	// subtract the CM of the two structures
	int i;
	for(i = 0; i < n; i++){
		Xr[i][0] -= cm_x_Xr;
		Xr[i][1] -= cm_y_Xr;
		Xr[i][2] -= cm_z_Xr;
		
		X0[i][0] -= cm_x_X0;
		X0[i][1] -= cm_y_X0;
		X0[i][2] -= cm_z_X0;
	}
	
	// compute Xr^tX0
	double **Xrt = trans_mat_lf(Xr, n, 3);
	double **Xrt_x_X0 = mat_times_mat_lf(Xrt, X0, 3, n, 3);
	dealloc_mat_lf(Xrt, 3);
		
	// compute USV^t = svd(Xr^tX0)
	double **U = alloc_mat_lf(3, 3);
	double **S = alloc_mat_lf(3, 3);
	double **Vt = alloc_mat_lf(3, 3);
	
	svd_row_major_matrix(Xrt_x_X0, 3, 3, U, S, Vt);
	dealloc_mat_lf(Xrt_x_X0, 3);
	dealloc_mat_lf(S, 3);
	// compute Q = V*Ut
	double **Ut = trans_mat_lf(U, 3, 3);
	dealloc_mat_lf(U, 3);
	double **V = trans_mat_lf(Vt, 3, 3);
	dealloc_mat_lf(Vt, 3);
	double **Q = mat_times_mat_lf(V, Ut, 3, 3, 3);
	dealloc_mat_lf(Ut, 3);
	dealloc_mat_lf(V, 3);
	
	// compute XrQ^t
	double **Qt = trans_mat_lf(Q, 3, 3);
	dealloc_mat_lf(Q, 3);
	double **Xr_x_Qt = mat_times_mat_lf(Xr, Qt, n, 3, 3);
	dealloc_mat_lf(Qt, 3);
	
	// compute X0 - XrQ^t
	double **X0_m_XrQt = zeros_mat_lf(n, 3);
	mat_m_mat_lf(X0_m_XrQt, X0, Xr_x_Qt, n, 3);
	dealloc_mat_lf(Xr_x_Qt, n);
	
	double rmsd = frobenius_norm(X0_m_XrQt, n, 3)/sqrt(n);
	
	dealloc_mat_lf(X0_m_XrQt, n);
	dealloc_mat_lf(X0, n);
	dealloc_mat_lf(Xr, n);
	
	return rmsd;
}
