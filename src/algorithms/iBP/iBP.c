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
#include "algorithms/common.h"
#include "algorithms/iBP.h"
#define MYZERO 0.000001
#define MYZERO2 0.000000001
#include "project_version.h"
#include "algorithms/iBP.h"
/* *********************************************************************************** */
const char *ibp_version(void)      { return IBP_VERSION; }
const char *ibp_release_date(void) { return IBP_RELEASE_DATE; }
/* *********************************************************************************** */
/**
 * @brief Samples an interval (positive and negative parts) and fills branches[i] with the sampled values.
 *
 * @param i		Index of the current vertex.
 * @param branches	Matrix [n][2*sampleSize] containing branch samples for each vertex.
 * @param branchNum	Array [n] storing the index where valid branches start in each row.
 * @param interval	type_matrix_interval containing positive and negative intervals.
 * @param sampleSize	Maximum number of samples per interval.
 * @param tolerance	Minimum spacing between sampled points.
 */
void sample_symmetric_ddgp_interval(int i, double **branches, int *branchNum, type_matrix_interval interval, int sampleSize, double tolerance) {
	
	int k, l;
	int numNegPoints = 0;
	int numPosPoints = 0;
	double sampleNeg[sampleSize], samplePos[sampleSize];
	
	// Positive intervals
	if (interval.npos_intervals > 0)
		numPosPoints = sampling_simple_interval_uniformly_with_tolerance(samplePos, interval.pos_interval[0], sampleSize, tolerance);
	
	// Negative intervals
	if (interval.nneg_intervals > 0)
		numNegPoints = sampling_simple_interval_uniformly_with_tolerance(sampleNeg, interval.neg_interval[0], sampleSize, tolerance);
	
	// Store the starting index of the sampled branches
	branchNum[i] = 2 * sampleSize - (numPosPoints + numNegPoints);
	
	// Copy negative samples
	l = branchNum[i];
	for (k = 0; k < numNegPoints; k++)
		branches[i][l++] = sampleNeg[k];
	
	// Copy positive samples
	for (k = 0; k < numPosPoints; k++)
		branches[i][l++] = samplePos[k];
}
/* *********************************************************************************** */
void iBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double tolerance, double angularResolution, double timeLimit, int GivenNumOfSols, run_metrics *runMetrics, double **distanceConstraints, int num_dc, double ***allSolutions, double solutionDifferenceThreshold) {
	// ---------- Variables initialization ----------
	time_t startTime, nowWall;
	clock_t cpuStart, cpuNow;
	
	int k, satisfiedPruneEdges;
	int i = 3;
	int lev = 2;
	int twoSampleSize = 2*sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	int i1, i2, i3, ik;
	int breakLoop = 0;
	int nocs = 0;
	
	long nosf = 0;
	
	double cpu_time_used = 0.0;
	double timeCheck = timeLimit - 10.0;
	double lambda_210, rho2_210, lambda_213, rho2_213, di1i2, diik_2_0, diik_2_l, diik_2_u;
	double noev = 0.000003;
	double maxMDE = -1.0;
        double maxLDE = -1.0;
        double minRMSD = 10000.0;
	double **A  = zeros_mat_lf(n, 3);
	double **B  = zeros_mat_lf(n, 3);
	double **C  = zeros_mat_lf(n, 3);
	double **Xr = zeros_mat_lf(n, 3);
	double **branches = ones_mat_lf(n, twoSampleSize);
	l_t_mat_lf(404.0, branches, n, twoSampleSize);
	
	type_matrix_interval tauii3;
	
	referential_x1_x2_x3(Xr, discretizationEdges_2);
		
	cpuStart = clock();
	startTime = time(NULL);
	while (1) {
		nowWall = time(NULL);
		if (difftime(nowWall, startTime) >= timeCheck) {
			cpuNow = clock();
			cpu_time_used = ((double) (cpuNow - cpuStart)) / CLOCKS_PER_SEC;

			if (cpu_time_used >= timeLimit) {
				printf("iBP: Time limit reached -- search terminated %s.\n",
					lev < n ? "before finding any solution" : "after finding at least one solution"
				);
				break;
			}
		}
		
		if (i == n) {
			breakLoop = handle_solution_cycle("iBP", &nosf, &maxMDE, &maxLDE, &minRMSD, &nocs, GivenNumOfSols, Xr, n, allSolutions, distanceConstraints, num_dc, solutionDifferenceThreshold, &i, exploredVertex, branches, branchNum, twoSampleSize);
			if(breakLoop)
				break;
		}
		
		if (exploredVertex[i] == 0) {
			// Compute cos(tau_{i,i3})
			compute_tauii3_interval(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &lambda_210, &rho2_210, &lambda_213, &rho2_213, &tauii3);
			
			sample_symmetric_ddgp_interval(i, branches, branchNum, tauii3, sampleSize, angularResolution);
			
			exploredVertex[i] = 1;
			circunference_parameters(i, A, B, C, &Xr[i3][0], &Xr[i2][0], &Xr[i1][0], lambda_210, rho2_210, di1i2);
		}
		else
			branches[i][branchNum[i]++] = 404.0;
		
		if (branchNum[i] < twoSampleSize) {
			vertex_embedding(&Xr[i][0], &A[i][0], &B[i][0], &C[i][0], branches[i][branchNum[i]]);
			noev += 0.000001;
		}
		else {
			tree_backtracking(&i, exploredVertex, branches, branchNum, twoSampleSize);
			if (i == 2) {
				printf("iBP: %s\n",
					nosf == 0
					? "ERROR -- No solution found after exhaustively exploring the entire search space."
					: "Solutions were found, and the entire search space has been exhaustively explored."
				);
				break;
			}
			continue;
		}
			
		satisfiedPruneEdges = 1;
		
		for (k = 0; k < pruneEdges_2[i].cardUkP; k++) {
			ik = (int) (pruneEdges_2[i].precise[k][0] - 1.0);
			
			diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
			
			diik_2_l = pruneEdges_2[i].precise[k][1];
			
			if(fabs(diik_2_l - diik_2_0) > tolerance){
				satisfiedPruneEdges = 0;
				break;
			}
		}
		
		if (satisfiedPruneEdges) {
			for (k = 0; k < pruneEdges_2[i].cardUkI; k++) {
				ik = (int) (pruneEdges_2[i].interval[k][0] - 1.0);
				
				diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
				
				diik_2_l = pruneEdges_2[i].interval[k][1];
				diik_2_u = pruneEdges_2[i].interval[k][2];
				
				if((diik_2_0 < fabs(diik_2_l - tolerance)) || (fabs(diik_2_u + tolerance) < diik_2_0)){
					satisfiedPruneEdges = 0;
					break;
				}
			}
		}
		
		if ((satisfiedPruneEdges) && (lev < ++i))
			lev = i;
	}
	
	dealloc_mat_lf(A, n);
	dealloc_mat_lf(B, n);
	dealloc_mat_lf(C, n);
	dealloc_mat_lf(Xr, n);
	dealloc_mat_lf(branches, n);
	free(exploredVertex);
	free(branchNum);
	
	if (cpu_time_used < MYZERO2) {
		cpuNow = clock();
		cpu_time_used = ((double) (cpuNow - cpuStart)) / CLOCKS_PER_SEC;
	}
	
	runMetrics->cpu_time_used = cpu_time_used;
	runMetrics->lev = lev;
	runMetrics->noev = 1000000*noev;
	runMetrics->nosf = nosf;
	runMetrics->nocs = nocs;
	runMetrics->maxMDE = maxMDE;
	runMetrics->maxLDE = maxLDE;
	runMetrics->minRMSD = minRMSD;
}
