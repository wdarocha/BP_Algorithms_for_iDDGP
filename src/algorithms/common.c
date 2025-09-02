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
#define PI 3.14159265359
#define MYZERO 0.000001

/* *********************************************************************************** */
/**
 * @brief Initializes the reference positions X[0], X[1], and X[2]. 
 *
 * Given the squared distances between the first three vertices, this function 
 * sets the coordinates of x1, x2 and x3 in a 3D referential. It assumes:
 * - x3 is at the origin: x3 = (0, 0, 0)
 * - x2 lies on the negative y-axis: x2 = (0, -d23, 0)
 * - x1 lies in the xy-plane, computed via triangle geometry.
 *
 * @param X					 X 2D array of shape [n][3] of zeros.
 * @param discretizationEdges_2 3D array containing discretization edges in the format:
 *							  [vertex][edge][data], where:
 *							  - data[0] = neighbor vertex
 *							  - data[1] = d^2_l squared lower bound
 *							  - data[2] = d^2_u squared upper bound
 */
void referential_x1_x2_x3(double **X, double ***discretizationEdges_2) {

	// d^2_{1,2}: squared distance between x1 and x2
	double d12_2 = discretizationEdges_2[1][0][1];
	// d^2_{1,3}: squared distance between x1 and x3
	double d13_2 = discretizationEdges_2[2][1][1];
	// d^2_{2,3}: squared distance between x2 and x3
	double d23_2 = discretizationEdges_2[2][0][1];
	double d23   = sqrt(d23_2);
	
	// Using the Law of Cosines to compute the projection of d12 on d23
	double d12cosTh = 0.5 * (d12_2 + d23_2 - d13_2) / d23;
	// Compute the orthogonal component (sine) using Pythagoras
	double d12sinTh = sqrt(d12_2 - d12cosTh * d12cosTh);

	// Set x1 = (x, y, 0)
	X[0][0] = -d12sinTh;		// x-component
	X[0][1] = d12cosTh - d23;	// y-component
	// Set x2 = (0, -d23, 0)
	X[1][1] = -d23;
	// Set x3 = (0, 0, 0), it is done when X is allocaded.
}
/* *********************************************************************************** */
/**
 * @brief Handles the storage and analysis of a candidate solution.
 *
 * This function saves the current solution 'Xr' if it is sufficiently distinct
 * from previously found solutions, based on an RMSD threshold. It also updates
 * MDE/LDE/RMSD statistics and performs a backtracking step.
 *
 * @param method 			Name of the method calling this routine (for logging).
 * @param nosf 				Pointer to the number of solutions found so far.
 * @param maxMDE 			Pointer to the maximum MDE found.
 * @param maxLDE 			Pointer to the maximum LDE found.
 * @param minRMSD	 		Pointer to the minimum RMSD found.
 * @param nocs 				Pointer to the number of candidate solutions stored.
 * @param GivenNumOfSols 		Target number of distinct solutions to find.
 * @param Xr 				Current solution matrix.
 * @param n 				Number of vertices (rows in Xr).
 * @param allSolutions			3D matrix storing all solutions found so far.
 * @param distanceConstraints		Distance constraints matrix.
 * @param num_dc			Number of distance constraints.
 * @param solutionDifferenceThreshold	Threshold to consider solutions distinct.
 * @param i				Pointer to the current vertex index in the exploration.
 * @param exploredVertex		Array tracking explored vertices.
 * @param branches Matrix		storing branches during the search.
 * @param branchNum Array		storing branch counts.
 * @param twoSampleSize			Maximum number of branches per vertex.
 * @return 1 if the search should stop, 0 otherwise.
 */
int handle_solution_cycle(const char *method, long int *nosf, double *maxMDE, double *maxLDE, double *minRMSD, int *nocs, int GivenNumOfSols, double **Xr, int n, double ***allSolutions, double **distanceConstraints, int num_dc, double solutionDifferenceThreshold, int *i, int *exploredVertex, double **branches, int *branchNum, int twoSampleSize) {
	
	(*nosf)++;
	
	int isDistinct = 0; 		// whether we will append to allSolutions[++(*nocs)]
	double rmsd_to_ref = 0.0; 	// RMSD w.r.t. reference (allSolutions[0])
	double mde, lde;
	
	// Case A: no stored distinct solutions yet -- accept the first by convention
	if (*nocs == 0) {
		// Initialize the reference only if still unset
		if(frobenius_norm(allSolutions[0], n, 3) < MYZERO)
			mat_e_mat_lf(allSolutions[0], Xr, n, 3); // set reference
		else
		// Reference already set: compute RMSD from reference solution
			rmsd_to_ref = compute_rmsd(allSolutions[0], Xr, n);
    		
		isDistinct = 1;	// First solution: store directly
	}
	// Case B: there are already stored distinct solutions -- test distinctness
	else {
		// Compute RMSD from reference solution
		rmsd_to_ref = compute_rmsd(allSolutions[0], Xr, n);
		
		isDistinct = 1;
		// Test if candidate is sufficiently distinct
		for (int k = 1; k <= (*nocs); k++) {
			double rmsd = compute_rmsd(allSolutions[k], Xr, n);
			if (rmsd < solutionDifferenceThreshold) {
				isDistinct = 0;
				break;  // early exit: not distinct with at least one solution
			}
		}
	}
	
	// Update global metrics (all candidates)
	compute_mde_and_lde(Xr, distanceConstraints, num_dc, &mde, &lde);
	(*maxMDE)  = maxAB_lf((*maxMDE), mde);
	(*maxLDE)  = maxAB_lf((*maxLDE), lde);
	(*minRMSD) = minAB_lf((*minRMSD), rmsd_to_ref);
		
	// If distinct enough, store new solution
	if (isDistinct)
		mat_e_mat_lf(allSolutions[++(*nocs)], Xr, n, 3);

	// Check stopping condition after storing
	if ((*nocs) == GivenNumOfSols) {
		printf("%s: The specified number of distinct solutions was found.\n", method);
		return 1;
	}
	
	// Backtracking step
	(*i)--;
	tree_backtracking(i, exploredVertex, branches, branchNum, twoSampleSize);
	if ((*i) == 2) {
		printf("%s: Solutions were found, and the entire search space has been exhaustively explored.\n", method);
		return 1;
	}
	
	return 0;
}
/* *********************************************************************************** */
/**
 * @brief Performs backtracking in a tree-based embedding process.
 *
 * This function updates the current embedding level *i by moving back in the tree.
 * It marks the current vertex as unexplored and resets the last branch tried.
 * The process continues until a vertex with remaining branches is found.
 *
 * The special marker 404.0 is assigned to the discarded branch.
 *
 * @param i			Pointer to the current embedding depth/level.
 * @param exploredVertex	Array marking whether a vertex has been explored (1) or not (0).
 * @param branches		2D array [n][2*sampleSize] storing angular samples per vertex.
 * @param branchNum		Array storing the current branch index for each vertex.
 * @param sampleSize		Total number of angular samples per vertex.
 */
void tree_backtracking(int *i, int *exploredVertex, double **branches, int *branchNum, int sampleSize) {
	
	// Mark the current vertex as unexplored and decrement the level
	exploredVertex[(*i)--] = 0;
	
	// Backtrack until a vertex with untried branches is found
	while (1) {
		branches[(*i)][branchNum[(*i)]] = 404.0;  // Mark branch as discarded
		if (branchNum[(*i)] < sampleSize)
			break;  // Found a vertex with remaining branches to try
		
		exploredVertex[(*i)--] = 0;  // Backtrack further
	}
}
