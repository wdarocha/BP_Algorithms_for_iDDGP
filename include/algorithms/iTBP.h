#pragma once
const char *itbp_version(void);
const char *itbp_release_date(void);
/* *********************************************************************************** */
int sampling_interval_union_uniformly_with_tolerance(double *sample, double intervals[][2], int numIntervals, int sampleSize, double tolerance);
/* *********************************************************************************** */
void sample_ddgp_interval(int i, double **branches, int *branchNum, type_matrix_interval interval, int sampleSize, double tolerance);
/* *********************************************************************************** */
void change_referential_precise_case(type_matrix_interval *tauiik, double cosTau, double phase);
/* *********************************************************************************** */
void change_referential_interval_case(type_matrix_interval *tauiik, double cosTauk_l, double cosTauk_u, double phase);
/* *********************************************************************************** */
void get_tauii3_interval_from_angle(type_matrix_interval *tau, double tauAst, double deltaTau);
/* *********************************************************************************** */
int compute_tauiik_interval(double **X, double *pruneEdges_2, double lambda_210, double rho2_210, double lambda_213, double rho2_213, double di1i2, int i3, int i2, int i1, type_matrix_interval *tauiik);
/* *********************************************************************************** */
void compute_tauiik_precise(double **X, double *pruneEdges_2, double lambda_210, double rho2_210, double lambda_213, double rho2_213, double di1i2, int i3, int i2, int i1, type_matrix_interval *tauiik);
/* *********************************************************************************** */
void iTBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double angularResolution, double timeLimit, int GivenNumOfSols, int *signTau, double *givenTau, double *givenTauDeviation, run_metrics *runMetrics, double **distanceConstraints, int num_dc, double ***allSolutions, double solutionDifferenceThreshold);
/* *********************************************************************************** */
