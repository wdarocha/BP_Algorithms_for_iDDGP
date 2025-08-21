#pragma once
const char *ibp_version(void);
const char *ibp_release_date(void);
/* *********************************************************************************** */
void sample_symmetric_ddgp_interval(int i, double **branches, int *branchNum, type_matrix_interval cosTau, int sampleSize, double tolerance);
/* *********************************************************************************** */
void iBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double tolerance, double angularResolution, double timeLimit, int GivenNumOfSols, run_metrics *runMetrics, double **distanceConstraints, int num_dc, double ***allSolutions, double solutionDifferenceThreshold, int referenceSolutionIndex);
/* *********************************************************************************** */
