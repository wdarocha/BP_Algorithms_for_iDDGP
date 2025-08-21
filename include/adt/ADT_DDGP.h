/* *********************************************************************************** */
int *adjacent_predecessors_cardinality(int *vertices, int m, int n);
/* *********************************************************************************** */
void dcinputfile_2_instance(dcinputfile *dc_vec, double ***distanceConstraints, double ****instance_adjlist, proteinstructure **protein, int **kv, int m, int *n);
/* *********************************************************************************** */
double dij_from_adjlist(double ***adjlist, int *cardUi, int i, int j, char l_or_u[]);
/* *********************************************************************************** */
void adjlist_2_discretization_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2);
/* *********************************************************************************** */
void adjlist_2_prune_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges_2);
/* *********************************************************************************** */
void adjlist_2_graph_parameters(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges, prune_edges_set **pruneEdges);
/* *********************************************************************************** */
void circunference_parameters(int i, double **A, double **B, double **C, double *xi3, double *xi2, double *xi1, double lambda_210, double rho2_210, double di1i2);
/* *********************************************************************************** */
void vertex_embedding(double *xi, double *A, double *B, double *C, double tau);
/* *********************************************************************************** */
void compute_cosTauii3_parameters(int i, double **X, double ***discretizationEdges, int *i1, int *i2, int *i3, double *di1i2, double *lambda_210, double *rho2_210, double *lambda_213, double *rho2_213);
/* *********************************************************************************** */
void compute_tauii3_interval(int i, double **X, double ***discretizationEdges_2, int *i1, int *i2, int *i3, double *di1i2, double *lambda_210, double *rho2_210, double *lambda_213, double *rho2_213, type_matrix_interval *tau);
/* *********************************************************************************** */
void compute_mde_and_lde(double **X, double **distanceConstraints, int m, double *mde, double *lde);
/* *********************************************************************************** */
double compute_rmsd(double **X, double **Y, int n);
/* *********************************************************************************** */
