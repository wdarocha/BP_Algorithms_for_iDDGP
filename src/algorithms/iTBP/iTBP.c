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
#include "algorithms/iTBP.h"
#define CF_DEG2RAD 0.01745329252
#define HALFPI 1.57079632680
#define PI 3.14159265359
#define TWOPI 6.28318530718
#define MYZERO 0.000001
#define MYZERO2 0.000000001
#include "project_version.h"
#include "algorithms/iTBP.h"
/* *********************************************************************************** */
const char *itbp_version(void)      { return ITBP_VERSION; }
const char *itbp_release_date(void) { return ITBP_RELEASE_DATE; }
/* *********************************************************************************** */
/**
 * @brief Samples a union of disjoint intervals with uniform spacing.
 * 
 * The function tries to distribute "sampleSize" points uniformly across the intervals.
 * - If the ideal spacing (totalWidth / (sampleSize - 1)) is >= tolerance, it uses sampleSize points.
 * - Otherwise, it forces spacing = tolerance and centers the points across the union.
 * 
 * Assumes:
 * - Intervals are sorted and disjoint.
 * - sampleSize and tolerance are consistent with the union total width.
 * 
 * @param sample	Output array to store sampled points.
 * @param intervals	Array of intervals [numIntervals][2].
 * @param numIntervals	Number of intervals in the union.
 * @param sampleSize	Desired number of samples.
 * @param tolerance	Minimum allowed spacing between points.
 * @return		Number of points actually sampled.
 */
 int sampling_interval_union_uniformly_with_tolerance(double *sample, double intervals[][2], int numIntervals, int sampleSize, double tolerance) {
	
	// Step 1: Compute total width of the union
	double totalWidth = 0.0;
	for (int i = 0; i < numIntervals; i++)
		totalWidth += (intervals[i][1] - intervals[i][0]);
	
	// Step 2: Degenerate case: too small or single sample requested
	if (totalWidth < tolerance || sampleSize <= 1) {
		sample[0] = (intervals[0][0] + intervals[0][1]) / 2.0;
		return 1;
	}
	
	// Step 3: Compute spacing and centering
	double h_ideal = totalWidth / (sampleSize - 1);
	double h, start;
	int numPoints;
	
	if (h_ideal < tolerance) {
		// Case 1: Force spacing = tolerance and center the sampling
		h = tolerance;
		numPoints = (int)(totalWidth / h) + 1;
		
		double coveredWidth = (numPoints - 1) * h;
		double gap = totalWidth - coveredWidth;
		start = gap / 2.0;
	}
	else {
		// Case 2: Use requested number of samples with ideal spacing
		h = h_ideal;
		numPoints = sampleSize;
		start = 0.0;
	}
	
	// Step 4: Begin sampling across the union
	int currentInterval = 0;
	double si = intervals[currentInterval][0] + start;
	int count = 1;
	sample[0] = si;
	
	while (count < numPoints) {
		si = sample[count - 1] + h;
		
		// Advance to the next interval if si passed the end
		while (si > intervals[currentInterval][1] + MYZERO) {
			double overflow = si - intervals[currentInterval][1];
			currentInterval++;
			si = intervals[currentInterval][0] + overflow;
		}
		
		sample[count++] = si;
	}
	
	return numPoints;
}
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
void sample_ddgp_interval(int i, double **branches, int *branchNum, type_matrix_interval interval, int sampleSize, double tolerance) {
	
	int k, l;
	int numNegPoints = 0;
	int numPosPoints = 0;
	double sampleNeg[sampleSize], samplePos[sampleSize];
	
	// Positive intervals
	if (interval.npos_intervals > 0) {
		if (interval.npos_intervals == 1)
			numPosPoints = sampling_simple_interval_uniformly_with_tolerance(samplePos, interval.pos_interval[0], sampleSize, tolerance);
		else
			numPosPoints = sampling_interval_union_uniformly_with_tolerance(samplePos, interval.pos_interval, interval.npos_intervals, sampleSize, tolerance);
	}
	
	// Negative intervals
	if (interval.nneg_intervals > 0) {
		if (interval.nneg_intervals == 1)
			numNegPoints = sampling_simple_interval_uniformly_with_tolerance(sampleNeg, interval.neg_interval[0], sampleSize, tolerance);
		else
			numNegPoints = sampling_interval_union_uniformly_with_tolerance(sampleNeg, interval.neg_interval, interval.nneg_intervals, sampleSize, tolerance);
	}
	
	// Store the starting index of the sampled branches
	branchNum[i] = 2 * sampleSize - (numPosPoints + numNegPoints);
	
	///*
	l = branchNum[i];
	// Copy negative samples
	for (k = 0; k < numNegPoints; k++)
		branches[i][l++] = sampleNeg[k];
	// Copy positive samples
	for (k = 0; k < numPosPoints; k++)
		branches[i][l++] = samplePos[k];
	//*/
	/*
	double neg_interval_width = 0;
	for(k = 0; k < interval.nneg_intervals; k++)
		neg_interval_width += (interval.neg_interval[k][1] - interval.neg_interval[k][0]);
	
	double pos_interval_width = 0;
	for(k = 0; k < interval.npos_intervals; k++)
		pos_interval_width += (interval.pos_interval[k][1] - interval.pos_interval[k][0]);
	
	if(neg_interval_width > pos_interval_width){
		l = branchNum[i];
		// Copy negative samples
		for (k = 0; k < numNegPoints; k++)
			branches[i][l++] = sampleNeg[k];
		// Copy positive samples
		for (k = 0; k < numPosPoints; k++)
			branches[i][l++] = samplePos[k];
	}
	else {
		l = branchNum[i];
		// Copy positive samples
		for (k = 0; k < numPosPoints; k++)
			branches[i][l++] = samplePos[k];
		// Copy negative samples
		for (k = 0; k < numNegPoints; k++)
			branches[i][l++] = sampleNeg[k];
	}
	*/
}
/* *********************************************************************************** */
/**
 * @brief Converts cosine bounds into angular intervals with respect to a phase shift.
 *
 * Given lower and upper cosine values (cosTauk_l and cosTauk_u), this function converts
 * the corresponding angular bounds into intervals in the sine-positive ([0, PI]) and
 * sine-negative ([-PI, 0]) domains, transformed by a reference phase.
 *
 * This version assumes cosTauk_l != cosTauk_u.
 *
 * @param tauiik	Output structure to store the transformed angular intervals.
 * @param cosTauk_l	Lower cosine bound (corresponds to upper angle).
 * @param cosTauk_u	Upper cosine bound (corresponds to lower angle).
 * @param phase		Reference angular shift for transformation.
 */
void change_referential_interval_case(type_matrix_interval *tauiik, double cosTauk_l, double cosTauk_u, double phase){
	
	if (fabs(cosTauk_u + 1.0) < MYZERO) {
	// absTauk_u = PI
		if (fabs(cosTauk_l - 1.0) < MYZERO) {
		// absTauk_u = PI & absTauk_l = 0
			tauiik->nneg_intervals = 1;
			set_interval(tauiik->neg_interval, 0, -PI, 0.0);
			
			tauiik->npos_intervals = 1;
			set_interval(tauiik->pos_interval, 0, 0.0, PI);
		}
		else {
		// absTauk_u = PI & absTauk_l != 0
			double absTauk_l = acos(cosTauk_l);
			double tauk_pos_l = change_referential_1_to_0(phase,  absTauk_l);
			double tauk_neg_l = change_referential_1_to_0(phase, -absTauk_l);
			int sign_tauk_pos_l = sign_lf(tauk_pos_l);
			int sign_tauk_neg_l = sign_lf(tauk_neg_l);
			
			if ((sign_tauk_pos_l > 0) && (sign_tauk_neg_l > 0)) {
				if (fabs(PI - absTauk_l) > HALFPI) { // if (|absTauk_u - absTauk_l| > HALFPI)
					tauiik->nneg_intervals = 1;
					set_interval(tauiik->neg_interval, 0, -PI, 0.0);
					
					tauiik->npos_intervals = 2;
					set_interval(tauiik->pos_interval, 0, 0.0, tauk_neg_l);
					set_interval(tauiik->pos_interval, 1, tauk_pos_l, PI);
				}
				else {
					tauiik->nneg_intervals = 0;
					
					tauiik->npos_intervals = 1;
					set_interval(tauiik->pos_interval, 0, tauk_pos_l, tauk_neg_l);
				}
			}
			else if ((sign_tauk_pos_l < 0) && (sign_tauk_neg_l < 0)) {
				if (fabs(PI - absTauk_l) > HALFPI) { // if (|absTauk_u - absTauk_l| > HALFPI)
					tauiik->nneg_intervals = 2;
					set_interval(tauiik->neg_interval, 0, -PI, tauk_neg_l);
					set_interval(tauiik->neg_interval, 1, tauk_pos_l, 0.0);
					
					tauiik->npos_intervals = 1;
					set_interval(tauiik->pos_interval, 0, 0.0, PI);
				}
				else {
					tauiik->nneg_intervals = 1;
					set_interval(tauiik->neg_interval, 0, tauk_pos_l, tauk_neg_l);
					
					tauiik->npos_intervals = 0;
				}
			}
			else if ((sign_tauk_pos_l > 0) && (sign_tauk_neg_l < 0)) {
				tauiik->nneg_intervals = 1;
				set_interval(tauiik->neg_interval, 0, -PI, tauk_neg_l);
				
				tauiik->npos_intervals = 1;
				set_interval(tauiik->pos_interval, 0, tauk_pos_l, PI);
			}
			else if ((sign_tauk_pos_l < 0) && (sign_tauk_neg_l > 0)) {
				tauiik->nneg_intervals = 1;
				set_interval(tauiik->neg_interval, 0, tauk_pos_l, 0.0);
				
				tauiik->npos_intervals = 1;
				set_interval(tauiik->pos_interval, 0, 0.0, tauk_neg_l);
			}
		}
	}
	else {
	// absTauk_u != PI
		if (fabs(cosTauk_l - 1.0) < MYZERO) {
		// absTauk_u != PI & absTauk_l = 0
			
			double absTauk_u = acos(cosTauk_u);
			double tauk_pos_u = change_referential_1_to_0(phase,  absTauk_u);
			double tauk_neg_u = change_referential_1_to_0(phase, -absTauk_u);
			int sign_tauk_pos_u = sign_lf(tauk_pos_u);
			int sign_tauk_neg_u = sign_lf(tauk_neg_u);
			
			if ((sign_tauk_pos_u > 0) && (sign_tauk_neg_u > 0)) {
				if (absTauk_u > HALFPI) { // if (|absTauk_u - absTauk_l| > HALFPI)
					tauiik->nneg_intervals = 1;
					set_interval(tauiik->neg_interval, 0, -PI, 0.0);
					
					tauiik->npos_intervals = 2;
					set_interval(tauiik->pos_interval, 0, 0.0, tauk_pos_u);
					set_interval(tauiik->pos_interval, 1, tauk_neg_u, PI);
				}
				else {
					tauiik->nneg_intervals = 0;
					
					tauiik->npos_intervals = 1;
					set_interval(tauiik->pos_interval, 0, tauk_neg_u, tauk_pos_u);
				}
			}
			else if ((sign_tauk_pos_u < 0) && (sign_tauk_neg_u < 0)) {
				if(absTauk_u > HALFPI) {  // if (|absTauk_u - absTauk_l| > HALFPI)
					tauiik->nneg_intervals = 2;
					set_interval(tauiik->neg_interval, 0, -PI, tauk_pos_u);
					set_interval(tauiik->neg_interval, 1, tauk_neg_u, 0.0);
					
					tauiik->npos_intervals = 1;
					set_interval(tauiik->pos_interval, 0, 0.0, PI);
				}
				else {
					tauiik->nneg_intervals = 1;
					set_interval(tauiik->neg_interval, 0, tauk_neg_u, tauk_pos_u);
					
					tauiik->npos_intervals = 0;
				}
			}
			else if ((sign_tauk_pos_u > 0) && (sign_tauk_neg_u < 0)) {
				tauiik->nneg_intervals = 1;
				set_interval(tauiik->neg_interval, 0, tauk_neg_u, 0.0);
				
				tauiik->npos_intervals = 1;
				set_interval(tauiik->pos_interval, 0, 0.0, tauk_pos_u);
			}
			else if ((sign_tauk_pos_u < 0) && (sign_tauk_neg_u > 0)) {
				tauiik->nneg_intervals = 1;
				set_interval(tauiik->neg_interval, 0, -PI, tauk_pos_u);
				
				tauiik->npos_intervals = 1;
				set_interval(tauiik->pos_interval, 0, tauk_neg_u, PI);
			}
		}
		else {
		// absTauk_u != PI & absTauk_l != 0
		
			double absTauk_l = acos(cosTauk_l);
			double absTauk_u = acos(cosTauk_u);
			double tauk_pos_u = change_referential_1_to_0(phase,  absTauk_u);
			double tauk_pos_l = change_referential_1_to_0(phase,  absTauk_l);
			double tauk_neg_u = change_referential_1_to_0(phase, -absTauk_u);
			double tauk_neg_l = change_referential_1_to_0(phase, -absTauk_l);
			
			int sign_tauk_pos_l = sign_lf(tauk_pos_l);
			int sign_tauk_pos_u = sign_lf(tauk_pos_u);
			int sign_tauk_neg_l = sign_lf(tauk_neg_l);
			int sign_tauk_neg_u = sign_lf(tauk_neg_u);
			/* ******************************************************** *
			label:
			sign_tauk_pos_u sign_tauk_pos_l sign_tauk_neg_u sign_tauk_neg_l
			0 means > 0
			1 means < 0
			
			0 0 0 0
			0 0 0 1
			0 0 1 0
			0 0 1 1
			
			0 1 0 0
			0 1 0 1
			0 1 1 0 XXX
			0 1 1 1
			
			1 0 0 0
			1 0 0 1 XXX
			1 0 1 0
			1 0 1 1
			
			1 1 0 0
			1 1 0 1
			1 1 1 0
			1 1 1 1
			* ******************************************************** */
			if (sign_tauk_pos_u > 0) {
				if (sign_tauk_pos_l > 0) {
					if (sign_tauk_neg_u > 0) {
						if (sign_tauk_neg_l > 0) {
							// 0 0 0 0
							tauiik->nneg_intervals = 0;
							
							tauiik->npos_intervals = 2;
							if (tauk_pos_l < tauk_neg_u) {
								set_interval(tauiik->pos_interval, 0, tauk_pos_l, tauk_pos_u);
								set_interval(tauiik->pos_interval, 1, tauk_neg_u, tauk_neg_l);
							}
							else {
								set_interval(tauiik->pos_interval, 0, tauk_neg_u, tauk_neg_l);
								set_interval(tauiik->pos_interval, 1, tauk_pos_l, tauk_pos_u);
							}
						}
						else {
							// 0 0 0 1
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, -PI, tauk_neg_l);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, tauk_pos_l, tauk_pos_u);
							set_interval(tauiik->pos_interval, 1, tauk_neg_u, PI);
						}
					}
					else {
						if (sign_tauk_neg_l > 0) {
							// 0 0 1 0
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, tauk_neg_u, 0.0);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_neg_l);
							set_interval(tauiik->pos_interval, 1, tauk_pos_l, tauk_pos_u);
						}
						else {
							// 0 0 1 1
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, tauk_neg_u, tauk_neg_l);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->pos_interval, 0, tauk_pos_l, tauk_pos_u);
						}
					}
				}
				else {
					if (sign_tauk_neg_u > 0) {
						if (sign_tauk_neg_l > 0) {
							// 0 1 0 0
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, tauk_pos_l, 0.0);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_pos_u);
							set_interval(tauiik->pos_interval, 1, tauk_neg_u, tauk_neg_l);
						}
						else {
							// 0 1 0 1
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->neg_interval, 0, -PI, tauk_neg_l);
							set_interval(tauiik->neg_interval, 1, tauk_pos_l, 0.0);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_pos_u);
							set_interval(tauiik->pos_interval, 1, tauk_neg_u, PI);
						}
					}
					else {
						if (sign_tauk_neg_l > 0)
							// 0 1 1 0
							printf("Impossible case\n");
						else {
							// 0 1 1 1
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->neg_interval, 0, tauk_neg_u, tauk_neg_l);
							set_interval(tauiik->neg_interval, 1, tauk_pos_l, 0.0);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_pos_u);
						}
					}
				}
			}
			else {
				if (sign_tauk_pos_l > 0) {
					if (sign_tauk_neg_u > 0) {
						if (sign_tauk_neg_l > 0) {
							// 1 0 0 0
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, -PI, tauk_pos_u);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, tauk_neg_u, tauk_neg_l);
							set_interval(tauiik->pos_interval, 1, tauk_pos_l, PI);
						}
						else
							// 1 0 0 1
							printf("Impossible case\n");
					}
					else {
						if (sign_tauk_neg_l > 0) {
							// 1 0 1 0
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->neg_interval, 0, -PI, tauk_pos_u);
							set_interval(tauiik->neg_interval, 1, tauk_neg_u, 0.0);
							
							tauiik->npos_intervals = 2;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_neg_l);
							set_interval(tauiik->pos_interval, 1, tauk_pos_l, PI);
						}
						else {
							// 1 0 1 1
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->pos_interval, 0, tauk_pos_l, PI);
							set_interval(tauiik->neg_interval, 0, -PI, tauk_pos_u);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->neg_interval, 1, tauk_neg_u, tauk_neg_l);
						}
					}
				}
				else {
					if (sign_tauk_neg_u > 0) {
						if (sign_tauk_neg_l > 0) {
							// 1 1 0 0
							tauiik->nneg_intervals = 1;
							set_interval(tauiik->neg_interval, 0, tauk_pos_l, tauk_pos_u);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->pos_interval, 0, tauk_neg_u, tauk_neg_l);
						}
						else {
							// 1 1 0 1
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->neg_interval, 0, -PI, tauk_neg_l);
							set_interval(tauiik->neg_interval, 1, tauk_pos_l, tauk_pos_u);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->pos_interval, 0, tauk_neg_u, PI);
						}
					}
					else {
						if (sign_tauk_neg_l > 0) {
							// 1 1 1 0
							tauiik->nneg_intervals = 2;
							set_interval(tauiik->neg_interval, 0, tauk_pos_l, tauk_pos_u);
							set_interval(tauiik->neg_interval, 1, tauk_neg_u, 0.0);
							
							tauiik->npos_intervals = 1;
							set_interval(tauiik->pos_interval, 0, 0.0, tauk_neg_l);
						}
						else {
							// 1 1 1 1
							tauiik->npos_intervals = 0;
							
							tauiik->nneg_intervals = 2;
							if (tauk_pos_l < tauk_neg_u) {
								set_interval(tauiik->neg_interval, 0, tauk_pos_l, tauk_pos_u);
								set_interval(tauiik->neg_interval, 1, tauk_neg_u, tauk_neg_l);
							}
							else {
								set_interval(tauiik->neg_interval, 0, tauk_neg_u, tauk_neg_l);
								set_interval(tauiik->neg_interval, 1, tauk_pos_l, tauk_pos_u);
							}
						}
					}
				}
			}
		}
	}
}
/* *********************************************************************************** */
/**
 * @brief Converts a degenerate cosine bound into angular intervals at a fixed point.
 *
 * This function handles the special case where cosTauk_l == cosTauk_u, thus representing
 * a single angular point. It creates a single degenerate interval in the appropriate
 * sine-positive or sine-negative domain.
 *
 * @param tauiik	Output structure to store the angular interval.
 * @param cosine	Cosine value (cosTauk_l == cosTauk_u).
 * @param phase		Reference angular shift for transformation.
 */
void change_referential_precise_case(type_matrix_interval *tauiik, double absTauk_m, double phase) {
	
	double tauk_pos = change_referential_1_to_0(phase,  absTauk_m);
	double tauk_neg = change_referential_1_to_0(phase, -absTauk_m);
	
	int sign_tauk_pos = sign_lf(tauk_pos);
	int sign_tauk_neg = sign_lf(tauk_neg);
	
	if ((sign_tauk_pos > 0) && (sign_tauk_neg > 0)) {
		double tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		double tau_k2 = maxAB_lf(tauk_pos, tauk_neg);

		tauiik->nneg_intervals = 0;
		
		if (fabs(tau_k2 - tau_k1) < MYZERO) {
			tauiik->npos_intervals = 1;
			set_interval(tauiik->pos_interval, 0, tau_k1, tau_k1);
		}
		else {
			tauiik->npos_intervals = 2;
			set_interval(tauiik->pos_interval, 0, tau_k1, tau_k1);
			set_interval(tauiik->pos_interval, 1, tau_k2, tau_k2);
		}
	}
	else if ((sign_tauk_pos < 0) && (sign_tauk_neg < 0)) {
		double tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		double tau_k2 = maxAB_lf(tauk_pos, tauk_neg);

		if (fabs(tau_k2 - tau_k1) < MYZERO) {
			tauiik->nneg_intervals = 1;
			set_interval(tauiik->neg_interval, 0, tau_k1, tau_k1);
		}
		else {
			tauiik->nneg_intervals = 2;
			set_interval(tauiik->neg_interval, 0, tau_k1, tau_k1);
			set_interval(tauiik->neg_interval, 1, tau_k2, tau_k2);
		}
		
		tauiik->npos_intervals = 0;
	}
	else if ((sign_tauk_pos > 0) && (sign_tauk_neg < 0)) {
		tauiik->nneg_intervals = 1;
		set_interval(tauiik->neg_interval, 0, tauk_neg, tauk_neg);
		
		tauiik->npos_intervals = 1;
		set_interval(tauiik->pos_interval, 0, tauk_pos, tauk_pos);
	}
	else if ((sign_tauk_pos < 0) && (sign_tauk_neg > 0)) {
		tauiik->nneg_intervals = 1;
		set_interval(tauiik->neg_interval, 0, tauk_pos, tauk_pos);
		
		tauiik->npos_intervals = 1;
		set_interval(tauiik->pos_interval, 0, tauk_neg, tauk_neg);
	}
}
/* *********************************************************************************** */
void get_tauii3_interval_from_angle(type_matrix_interval *tau, double tauAst, double deltaTau) {

	double tau_l = tauAst - deltaTau;
	double tau_u = tauAst + deltaTau;

	int sign_tau_l = sign_lf(tau_l);
	int sign_tau_u = sign_lf(tau_u);
	
	if ((sign_tau_u > 0) && (sign_tau_l > 0)) {
		if (tau_u > PI) { // Interval crosses from positive to negative region
			tau->nneg_intervals = 1;
			set_interval(tau->neg_interval, 0, -PI, tau_u - TWOPI);
			
			tau->npos_intervals = 1;
			set_interval(tau->pos_interval, 0, tau_l, PI);
		}
		else { // Interval fully in the positive region
			tau->nneg_intervals = 0;
			
			tau->npos_intervals = 1;
			set_interval(tau->pos_interval, 0, tau_l, tau_u);
		}
	}
	else if ((sign_tau_u < 0) && (sign_tau_l < 0)) {
		if (tau_l < -PI) { // Interval crosses from negative to positive region
			tau->nneg_intervals = 1;
			set_interval(tau->neg_interval, 0, -PI, tau_u);
			
			tau->npos_intervals = 1;
			set_interval(tau->pos_interval, 0, tau_l + TWOPI, PI);
		}
		else { // Interval fully in the negative region
			tau->nneg_intervals = 1;
			set_interval(tau->neg_interval, 0, tau_l, tau_u);
			
			tau->npos_intervals = 0;
		}
	}
	else if ((sign_tau_u > 0) && (sign_tau_l < 0)) {
		// Interval crosses from negative to positive region
		tau->nneg_intervals = 1;
		set_interval(tau->neg_interval, 0, tau_l, 0.0);
		
		tau->npos_intervals = 1;
		set_interval(tau->pos_interval, 0, 0.0, tau_u);
	}
	else if ((sign_tau_u < 0) && (sign_tau_l > 0)) {
		// Invalid case: upper bound negative, lower bound positive -> empty
		tau->nneg_intervals = 0;
		tau->npos_intervals = 0;
	}
}
/* *********************************************************************************** */
int compute_tauiik_interval(double **X, double *pruneEdges_2, double lambda_210, double rho2_210, double lambda_213, double rho2_213, double di1i2, int i3, int i2, int i1, type_matrix_interval *tauiik) {
	
	int ik = (int) (pruneEdges_2[0] - 1.0);
	
	double di1ik_2 = d2xixj_lf3(&X[ik][0], &X[i1][0]);
	double di2ik_2 = d2xixj_lf3(&X[ik][0], &X[i2][0]);
	double di3ik_2 = d2xixj_lf3(&X[ik][0], &X[i3][0]);
		
	double lambda_21k = lambdai_d2(di2ik_2, di1i2, di1ik_2);
	double rho2_21k   = rhoi2_d2(di2ik_2, lambda_21k);
	
	double lbdmlbdb = lambda_210 - lambda_21k;
	double pk0 = lbdmlbdb * lbdmlbdb + rho2_210 + rho2_21k;
	double qk0 = sqrt(rho2_210*rho2_21k);
	
	double diik_2_l = pruneEdges_2[1];
	double diik_2_u = pruneEdges_2[2];
	
	double cosTauk_l = cos_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	double cosTauk_u = cos_torsion_angle_with_constants_d2(pk0, qk0, diik_2_u);
	
	if (fabs(cosTauk_l - cosTauk_u) < MYZERO2)
		return 1;
	else if (fabs((cosTauk_l - cosTauk_u) - 2.0) < MYZERO2)
		return -1;
	else {
		lbdmlbdb = lambda_213 - lambda_21k;
		double pk3 = lbdmlbdb * lbdmlbdb + rho2_21k + rho2_213;
		double qk3 = sqrt(rho2_21k*rho2_213);
				
		double cosPhi = cos_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
		
		double signPhi = (double) sign_torsion_angle(&X[i3][0], &X[i2][0], &X[i1][0], &X[ik][0]);
		
		double phi = signPhi * acos(cosPhi);
		
		change_referential_interval_case(tauiik, cosTauk_l, cosTauk_u, phi);
		
		return 0;
	}
}
/* *********************************************************************************** */
void compute_tauiik_precise(double **X, double *pruneEdges_2, double lambda_210, double rho2_210, double lambda_213, double rho2_213, double di1i2, int i3, int i2, int i1, type_matrix_interval *tauiik) {
	
	int ik = (int) (pruneEdges_2[0] - 1.0);
	
	double di1ik_2 = d2xixj_lf3(&X[ik][0], &X[i1][0]);
	double di2ik_2 = d2xixj_lf3(&X[ik][0], &X[i2][0]);
	double di3ik_2 = d2xixj_lf3(&X[ik][0], &X[i3][0]);
		
	double lambda_21k = lambdai_d2(di2ik_2, di1i2, di1ik_2);
	double rho2_21k   = rhoi2_d2(di2ik_2, lambda_21k);
	
	double lbdmlbdb = lambda_210 - lambda_21k;
	double pk0 = lbdmlbdb * lbdmlbdb + rho2_210 + rho2_21k;
	double qk0 = sqrt(rho2_210*rho2_21k);
	
	double diik_2_l = pruneEdges_2[1];
	
	double cosTauiik = cos_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	
	lbdmlbdb = lambda_213 - lambda_21k;
	double pk3 = lbdmlbdb * lbdmlbdb + rho2_21k + rho2_213;
	double qk3 = sqrt(rho2_21k*rho2_213);
				
	double cosPhi = cos_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
	
	double signPhi = (double) sign_torsion_angle(&X[i3][0], &X[i2][0], &X[i1][0], &X[ik][0]);
		
	double phi = signPhi * acos(cosPhi);
		
	change_referential_precise_case(tauiik, acos(cosTauiik), phi);
}
/* *********************************************************************************** */
void iTBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double angularResolution, double timeLimit, int GivenNumOfSols, int *signTau, double *givenTau, double *givenTauDeviation, run_metrics *runMetrics, double **distanceConstraints, int num_dc, double ***allSolutions, double solutionDifferenceThreshold) {

	// ---------- Variables initialization ----------
	time_t startTime, nowWall;
	clock_t cpuStart, cpuNow;
		
	int k, isTauIntervalEmpty;
	int i = 3;
	int lev = 2;
	int twoSampleSize = 2 * sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	int i1, i2, i3;
	int breakLoop = 0;
	int nocs = 0;
	
	long nosf = 0;
	
	double cpu_time_used = 0.0;
	double timeCheck = timeLimit - 10.0;
	double lambda_210, rho2_210, lambda_213, rho2_213, di1i2;
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
	
	type_matrix_interval tauiik, tauii3;
	
	int is_iabp = (norm2_d(signTau, n) > 0) ? 0 : 1;
	char method[6];
	sprintf(method, "i%sBP", is_iabp ? "A" : "T");
	
	referential_x1_x2_x3(Xr, discretizationEdges_2);
	
	double angularResolution0 = 0.5 * CF_DEG2RAD;
		
	cpuStart = clock();
	startTime = time(NULL);
	while (1) {
		nowWall = time(NULL);
		if (difftime(nowWall, startTime) >= timeCheck) {
			cpuNow = clock();
			cpu_time_used = ((double) (cpuNow - cpuStart)) / CLOCKS_PER_SEC;

			if (cpu_time_used >= timeLimit) {
				printf("i%sBP: Time limit reached -- search terminated %s.\n",
					is_iabp ? "A" : "T",
					lev < n ? "before finding any solution" : "after finding at least one solution"
				);
				break;
			}
		}
	
		if (i == n) {
			breakLoop = handle_solution_cycle(method, &nosf, &maxMDE, &maxLDE, &minRMSD, &nocs, GivenNumOfSols, Xr, n, allSolutions, distanceConstraints, num_dc, solutionDifferenceThreshold, &i, exploredVertex, branches, branchNum, twoSampleSize);
			if(breakLoop)
				break;
		}
				
		if(exploredVertex[i] == 0){
			// Compute cos(tau_{i,i3})
			if(signTau[i] == 0)
				// iABP case
				compute_tauii3_interval(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &lambda_210, &rho2_210, &lambda_213, &rho2_213, &tauii3);
			else {
				// iTBP case
				compute_cosTauii3_parameters(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &lambda_210, &rho2_210, &lambda_213, &rho2_213);
				
				get_tauii3_interval_from_angle(&tauii3, signTau[i]*givenTau[i], givenTauDeviation[i]);
			}
			
			isTauIntervalEmpty = 0;
			
			// compute the angular interval for precise distances
			for(k = 0; k < pruneEdges_2[i].cardUkP; k++) {
				// Compute cos(tau_{i,ik}) precise case
				compute_tauiik_precise(Xr, pruneEdges_2[i].precise[k], lambda_210, rho2_210, lambda_213, rho2_213, di1i2, i3, i2, i1, &tauiik);
				
				matrix_intervals_intersection(&tauii3, &tauiik, angularResolution0);
					
				if ((tauii3.npos_intervals + tauii3.nneg_intervals) == 0) {
					isTauIntervalEmpty = 1;
					break;
				}
			}
			
			// compute the angular interval for interval distances
			if (isTauIntervalEmpty != 1) {
				for (k = 0; k < pruneEdges_2[i].cardUkI; k++) {
					// Compute cos(tau_{i,ik}) interval case
					isTauIntervalEmpty = compute_tauiik_interval(Xr, pruneEdges_2[i].interval[k], lambda_210, rho2_210, lambda_213, rho2_213, di1i2, i3, i2, i1, &tauiik);
					
					if (isTauIntervalEmpty == 1)
						break;
					else if (isTauIntervalEmpty == 0) {
						matrix_intervals_intersection(&tauii3, &tauiik, angularResolution0);
						
						if ((tauii3.npos_intervals + tauii3.nneg_intervals) == 0) {
							isTauIntervalEmpty = 1;
							break;
						}
					}
				}
			}
			// if the angular is not empty, then sample it
			if (isTauIntervalEmpty != 1) {
				sample_ddgp_interval(i, branches, branchNum, tauii3, sampleSize, angularResolution);
				exploredVertex[i] = 1;
				
				circunference_parameters(i, A, B, C, &Xr[i3][0], &Xr[i2][0], &Xr[i1][0], lambda_210, rho2_210, di1i2);
			}
			else
				branchNum[i] = twoSampleSize;
		}
		else
			branches[i][branchNum[i]++] = 404.0;
		
		if (branchNum[i] < twoSampleSize) {
			vertex_embedding(&Xr[i][0], &A[i][0], &B[i][0], &C[i][0], branches[i][branchNum[i]]);
			
			noev+= 0.000001;
			if(lev < ++i)
				lev = i;
		}
		else {
			tree_backtracking(&i, exploredVertex, branches, branchNum, twoSampleSize);
			if (i == 2) {
				printf("i%sBP: %s\n", is_iabp ? "A" : "T", 
					nosf == 0 ? 
					"ERROR -- No solution found after exhaustively exploring the entire search space." :
					"Solutions were found, and the entire search space has been exhaustively explored."
				);
				break;
			}
			continue;
		}
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
