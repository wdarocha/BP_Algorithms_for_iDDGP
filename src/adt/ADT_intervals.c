#include <stdio.h>
#include <math.h>
#include "adt/ADT_math.h"
#include "adt/ADT_intervals.h"
#define MYZERO 0.000000001
#define HALFPI 1.57079632680
#define PI 3.14159265359

/* *********************************************************************************** */
/**
 * @brief Sets the bounds of an interval in a 2D array.
 *
 * @param matrix_interval	Target interval matrix [n][2].
 * @param index			Interval index (0 or 1).
 * @param a			Minimum value of the interval.
 * @param b			Maximum value of the interval.
 */
void set_interval(double matrix_interval[][2], int index, double a, double b) {
	
	matrix_interval[index][0] = a;
	matrix_interval[index][1] = b;
}
/* *********************************************************************************** */
/**
 * @brief Builds a symmetric interval structure with sign symmetry.
 *
 * This function assigns the same angular interval [a, b] to both the positive
 * and negative sine domains.
 *
 * @param interval	Pointer to the type_matrix_interval structure to be filled.
 * @param a		Lower bound of the interval.
 * @param b		Upper bound of the interval.
 */
void build_simple_symmetric_interval(type_matrix_interval *interval, double a, double b) {
	
	// One interval for the region \mathbb{R}_+
	interval->npos_intervals = 1;
	set_interval(interval->pos_interval, 0, a, b);
	
	// One interval for the region \mathbb{R}_-
	interval->nneg_intervals = 1;
	set_interval(interval->neg_interval, 0, -b, -a);
}
/* *********************************************************************************** */
/**
 * @brief Computes the intersection of two 1D intervals C = A cap B.
 *
 * If the intervals overlap (or touch within a specified numerical tolerance),
 * the intersection is stored in C and CisEmpty is set to 0.
 * Otherwise, CisEmpty is set to 1 and C is left unchanged.
 *
 * @param A		Array [2] representing interval A [a1, a2].
 * @param B		Array [2] representing interval B [b1, b2].
 * @param tolerance	Tolerance for floating-point comparison.
 * @param C		Output: Array [2] where the intersection [c1, c2] is stored.
 * @param CisEmpty	Output: 1 if the intersection is empty, 0 otherwise.
 */
void simple_intervals_intersection(const double *A, const double *B, double tolerance, double *C, int *CisEmpty) {
	
	// Compute the maximum of the lower bounds and the minimum of the upper bounds
	double a = maxAB_lf(A[0], B[0]);
	double b = minAB_lf(A[1], B[1]);
	
	// Case 1: intervals overlap normally
	if(a < b){
		(*CisEmpty) = 0;
		C[0] = a;
		C[1] = b;
	}
	// Case 2: intervals touch at a point within tolerance
	else if(fabs(a - b) < tolerance){
		(*CisEmpty) = 0;
		C[0] = minAB_lf(a, b);
		C[1] = maxAB_lf(a, b);
	}
	// Case 3: no intersection
	else
		(*CisEmpty) = 1;
}
/* *********************************************************************************** */
/**
 * @brief Computes pairwise intersections between two interval sets.
 *
 * For each interval A[i] and B[j], if they intersect (or nearly touch within tolerance),
 * the resulting intersection is stored in matrix C.
 *
 * @param A		Matrix [nA][2] of intervals.
 * @param nA		Number of intervals in A.
 * @param B		Matrix [nB][2] of intervals.
 * @param nB		Number of intervals in B.
 * @param tolerance	Tolerance for floating-point comparison.
 * @param C		Output matrix [nA * nB][2] storing resulting intersections.
 * @param num_intervals	Output: Number of valid intersections stored in C.
 */
void intersect_interval_sets(const double A[][2], int nA, const double B[][2], int nB, double tolerance, double C[][2], int *num_intervals) {
	
	int i, j;
	int isCapEmptyFlag;
	(*num_intervals) = 0;
	
	for(i = 0; i < nA; i++)
		for(j = 0; j < nB; j++){
			simple_intervals_intersection(A[i], B[j], tolerance, C[(*num_intervals)], &isCapEmptyFlag);
			if(!isCapEmptyFlag)
				(*num_intervals)++;
		}
}
/* *********************************************************************************** */
/**
 * @brief Intersects the intervals in (*A) with those in (*B) and updates (*A) in place.
 *
 * This function computes the intersection of both the positive and negative interval sets
 * between the input A and B. If any of the input sets is empty, the corresponding
 * result is cleared.
 *
 * @param A		Pointer to the type_matrix_interval to be updated.
 * @param B		Pointer to the type_matrix_interval to intersect with.
 * @param tolerance	Tolerance used for floating-point comparisons in intersections.
 */
void matrix_intervals_intersection(type_matrix_interval *A, const type_matrix_interval *B, double tolerance) {

	int capacity = maxAB_d(A->npos_intervals + B->npos_intervals, A->nneg_intervals + B->nneg_intervals);
	double result[capacity][2];  // Assuming a maximum of 4 intersections in total.
	int i, num_intervals;
	
	/* -------- Positive Intervals -------- */
	if (A->npos_intervals > 0 && B->npos_intervals > 0) {
		intersect_interval_sets(
			(const double (*)[2]) A->pos_interval, A->npos_intervals,
			(const double (*)[2]) B->pos_interval, B->npos_intervals,
			tolerance, result, &num_intervals
		);
		num_intervals = merge_intervals(result, num_intervals, tolerance);
		
		for (i = 0; i < num_intervals; i++)
			set_interval(A->pos_interval, i, result[i][0], result[i][1]);
		
		A->npos_intervals = num_intervals;
	}
	else
		A->npos_intervals = 0;  // No positive intervals in A or B
	
	/* -------- Negative Intervals -------- */
	if (A->nneg_intervals > 0 && B->nneg_intervals > 0) {
		intersect_interval_sets(
			(const double (*)[2]) A->neg_interval, A->nneg_intervals,
			(const double (*)[2]) B->neg_interval, B->nneg_intervals,
			tolerance, result, &num_intervals
		);
		num_intervals = merge_intervals(result, num_intervals, tolerance);
		
		for (i = 0; i < num_intervals; i++)
			set_interval(A->neg_interval, i, result[i][0], result[i][1]);
		
		A->nneg_intervals = num_intervals;
	}
	else
		A->nneg_intervals = 0;  // No negative intervals in A or B
}
/* *********************************************************************************** */
/**
 * @brief Merges a list of possibly overlapping or adjacent intervals within a tolerance.
 *
 * The input is a matrix of n intervals of the form [a, b]. This function merges all intervals
 * that are overlapping or within "tolerance" of each other. The result is a simplified set of
 * disjoint intervals, written in-place in "matrix_interval".
 *
 * The intervals do not need to be sorted beforehand.
 *
 * @param matrix_interval	Array [n][2] of input intervals.
 * @param num_intervals		Number of input intervals (n).
 * @param tolerance		Allowed tolerance to consider two intervals as overlapping or adjacent.
 * @return			The number of intervals after merging.
 */
int merge_intervals(double matrix_interval[][2], int num_intervals, double tolerance) {
	
	if (num_intervals <= 1)
		return num_intervals;

	int i, j;
	double key_a, key_b;
	// Step 1: Sort intervals by starting point using insertion sort (in-place)
	for (i = 1; i < num_intervals; i++) {
		key_a = matrix_interval[i][0];
		key_b = matrix_interval[i][1];
		j = i - 1;
		while ((j >= 0) && (matrix_interval[j][0] > key_a)) {
			set_interval(matrix_interval, j + 1, matrix_interval[j][0], matrix_interval[j][1]);
			j--;
		}
		set_interval(matrix_interval, j + 1, key_a, key_b);
	}

	// Step 2: Merge intervals in-place
	int merged_count = 0;
	for (i = 1; i < num_intervals; i++)
		if (matrix_interval[merged_count][1] + tolerance >= matrix_interval[i][0]) {
			// Overlapping or adjacent: extend the upper bound
			if (matrix_interval[i][1] > matrix_interval[merged_count][1])
				matrix_interval[merged_count][1] = matrix_interval[i][1];
		}
		else
			// Disjoint: move to next slot
			set_interval(matrix_interval, ++merged_count, matrix_interval[i][0], matrix_interval[i][1]);

	return merged_count + 1;
}
/* *********************************************************************************** */
/**
 * @brief Uniformly samples a 1D interval [a, b] with spacing no smaller than `tolerance`.
 *
 * If the ideal spacing between "sampleSize" points is less than "tolerance",
 * the function reduces the number of points and enforces spacing = "tolerance",
 * centering the result within the interval.
 *
 * @param sample	Output array to store sampled points (must support up to sampleSize values).
 * @param interval	Array of two doubles [a, b] representing the sampling bounds.
 * @param sampleSize	Desired number of samples (must be >= 1).
 * @param tolerance	Minimum spacing between consecutive points.
 * @return		Actual number of points written to "sample".
 */
int sampling_simple_interval_uniformly_with_tolerance(double *sample, double *interval, int sampleSize, double tolerance) {
	
	double a = interval[0];
	double b = interval[1];
	double totalWidth = b - a;
	
	// Case 1: Interval too small or only one sample requested
	if (totalWidth < tolerance || sampleSize <= 1) {
		sample[0] = (a + b) / 2.0;
		return 1;
	}
	
	// Compute ideal spacing
	double h_ideal = totalWidth / (sampleSize - 1);
	double h;
	int numPoints;
	
	if (h_ideal < tolerance) {
		// Case 2: Ideal spacing too small -- force spacing = tolerance and center
		h = tolerance;
		numPoints = (int)(totalWidth / h) + 1;
		
		double coveredWidth = (numPoints - 1) * h;
		double gap = totalWidth - coveredWidth;
		sample[0] = a + gap / 2.0;
	}
	else {
		// Case 3: Ideal spacing acceptable -- use requested number of samples
		h = h_ideal;
		numPoints = sampleSize;
		sample[0] = a;
	}
	
	// Fill sample array
	for (int i = 1; i < numPoints; i++)
		sample[i] = sample[i - 1] + h;
	
	return numPoints;
}
