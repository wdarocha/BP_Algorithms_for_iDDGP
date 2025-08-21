#include <math.h>
#include "adt/ADT_math.h"
#include "adt/ADT_vectors.h"
#define MYZERO 0.000001
#define PI 3.14159265359
#define TWOPI 6.28318530718

/* *********************************************************************************** */
/**
 * @brief Returns the minimum of two integer values.
 *
 * @param	a First integer value.
 * @param	b Second integer value.
 * @return	The smaller of a and b.
 */
int minAB_d(int a, int b) {

	return (a < b) ? a : b;
}
/* *********************************************************************************** */
/**
 * @brief Returns the maximum of two integer values.
 *
 * @param	a First integer value.
 * @param	b Second integer value.
 * @return	The larger of a and b.
 */
int maxAB_d(int a, int b) {

	return (a > b) ? a : b;
}
/* *********************************************************************************** */
/**
 * @brief Returns the minimum of two double values.
 *
 * @param	a First double value.
 * @param	b Second double value.
 * @return	The smaller of a and b.
 */
double minAB_lf(double a, double b) {

	return (a < b) ? a : b;
}
/* *********************************************************************************** */
/**
 * @brief Returns the maximum of two double values.
 *
 * @param	a First double value.
 * @param	b Second double value.
 * @return	The larger of a and b.
 */
double maxAB_lf(double a, double b) {

	return (a > b) ? a : b;
}
/* *********************************************************************************** */
/**
 * @brief Computes the lambda_i parameter based on squared distances.
 *
 * Geometric formula derived from the Law of Cosines for internal coordinates.
 *
 * @param dii2_2   d^2(x_i, x_{i2}).
 * @param di1i2	d^2(x_{i1}, x_{i2}).
 * @param dii1_2   d^2(x_i, x_{i2}).
 * @return double  lambda_i value.
 */
double lambdai_d2(double dii2_2, double di1i2, double dii1_2) {
	
	return (dii2_2 + di1i2 * di1i2 - dii1_2) / (2.0 * di1i2);
}
/* *********************************************************************************** */
/**
 * @brief Computes the squared radius rho^2_i.
 *
 * Used in angular calculations when placing vertex x_i in 3D space.
 *
 * @param dii2_2   d^2(x_i, x_{i2}).
 * @param lbdi	 lambda_i value.
 * @return double  rho^2_i value.
 */
double rhoi2_d2(double dii2_2, double lbdi) {
	
	return (dii2_2 - lbdi * lbdi);
}
/* *********************************************************************************** */
/**
 * @brief Computes the cosine of a torsion angle from squared distance constraints.
 *
 * This function computes the cosine of the torsion angle tau that satisfies:
 * 	cos(tau) = (p - d^2(x_i, x_{i3})) / (2q)
 * for a given squared distance between points. The result is clamped to [-1, 1]
 * according to the valid domain of the cosine function.
 *
 * This approach avoids computing square roots by working with squared distances.
 *
 * @param p		Constant term: (lambda_i - lambdabar_i)^2 + rho^2_i + rhobar^2_i.
 * @param q		Constant term: rho_i * rhobar_i.
 * @param dii3_2	Squared distance: d^2(x_i, x_{i3}).
 * @return double	Cosine of the torsion angle, clamped in [-1, 1].
 */
double cos_torsion_angle_with_constants_d2(double p, double q, double dii3_2) {
	
	double twoq = 2.0 * q;
	double dmax_2 = p + twoq;
	
	// If d^2(x_i, x_{i3}) exceeds or touches the upper bound, the cosine is -1 (tau = pi)
	if (__builtin_expect(dmax_2 < dii3_2 || fabs(dmax_2 - dii3_2) < MYZERO, 0))
		return -1.0;
	
	double dmin_2 = p - twoq;
	
	// If within the valid angular domain, compute cos(tau)
	if (dmin_2 < dii3_2)
		return ((p - dii3_2) / twoq);
	
	// If d^2(x_i, x_{i3}) is below the lower bound, the cosine is 1 (tau = 0)
	return 1.0;
}
/* *********************************************************************************** */
/**
 * @brief Determines the sign of the torsion (dihedral) angle between four points in 3D space.
 *
 * Given four points x1, x2, x3, and x4, this function computes the sign of the torsion angle
 * defined by the planes (x1, x2, x3) and (x2, x3, x4). The sign is positive if the angle is 
 * counterclockwise (right-handed), negative if clockwise, and zero if the points are coplanar.
 *
 * Geometrically, this is determined by:
 * 	sign(dot((r × v), s))  where:
 *		- v = x1 - x2
 *		- r = x3 - x2
 *		- s = x4 - x2
 *
 * @param x1 is the coordinates in R^3 of point 1
 * @param x2 is the coordinates in R^3 of point 2
 * @param x3 is the coordinates in R^3 of point 3
 * @param x4 is the coordinates in R^3 of point 4
 * @return int +1 if torsion angle is greater than zero, -1 otherwise.
 */
int sign_torsion_angle(double *x1, double *x2, double *x3, double *x4) {
	
	// Define temporary 3D vectors for vector differences
	// double v[3], r[3], s[3], sperp[3];
	double v[3], r[3];
	
	// Compute displacement vectors relative to x2
	vec_m_vec_lf3(v, x1, x2);	// v = x1 - x2
	vec_m_vec_lf3(r, x3, x2);	// r = x3 - x2
	// vec_m_vec_lf3(s, x4, x2);	// s = x4 - x2
	
	// Compute perpendicular vector: sperp = r × v
	//cross_product_lf3(sperp, r, v);
	
	// compute of the scalar triple product cosTh = dot(sperp, s)
	double cosTh = 0.0;
	cosTh += (x4[0] - x2[0]) * (r[1]*v[2] - r[2]*v[1]);
	cosTh += (x4[1] - x2[1]) * (r[2]*v[0] - r[0]*v[2]);
	cosTh += (x4[2] - x2[2]) * (r[0]*v[1] - r[1]*v[0]);
	
	//return sign_double(dot_product_lf3(sperp, s));
	return sign_lf(cosTh);
}
/* *********************************************************************************** */
/**
 * @brief Returns the sign of a double value.
 *
 * Returns:
 *   +1 if v >= 0
 *   -1 if v < 0
 *
 * This version assumes that very small values (|v| < ε) should still be treated as positive.
 *
 * @param v   Input value.
 * @return int  +1 if v ≥ 0, -1 if v < 0
 */
int sign_lf(double v) {
	
	return (v >= -MYZERO) ? 1 : -1;
}
/* *********************************************************************************** */
/**
 * @brief Converts an angle from one torsion reference frame to another.
 *
 * Given:
 *   - phase: the angular offset between reference frames 1 and 0
 *   - tau_ref1: the torsion angle expressed in reference frame 1
 *
 * Returns:
 *   - tau_ref0: the same torsion angle expressed in reference frame 0,
 *                 normalized to the interval (-pi, pi]
 *
 * @param phase       Angular offset between referentials (in radians).
 * @param tau_ref1    Angle in reference frame 1 (in radians).
 * @return double     Equivalent angle in reference frame 0, \in (-pi, pi]
 */
double change_referential_1_to_0(double phase, double tau_ref1) {
	
	// Sum angle in new frame
	double tau_ref0 = phase + tau_ref1;
	
	// Normalize to [0, 2pi)
	tau_ref0 = fmod(tau_ref0 + PI, TWOPI);
	if (tau_ref0 < 0.0)
		tau_ref0 += TWOPI;
	
	// Shift to final interval (-pi, pi]
	tau_ref0 -= PI;
	
	return tau_ref0;
}
/* *********************************************************************************** */
/**
 * @brief Computes the squared Euclidean distance between two 3D points.
 *
 * @param x1 is a double 3x1 array
 * @param x2 is a double 3x1 array
 * @return double squared distance ||x1 - x2||^2
 */
double d2xixj_lf3(double *x1, double *x2) { 
	
	double aux;
	double norm2 = 0.0;
	aux = x1[0] - x2[0];
	norm2 += aux*aux;
	aux = x1[1] - x2[1];
	norm2 += aux*aux;
	aux = x1[2] - x2[2];
	norm2 += aux*aux;
	
	return norm2;
}
/* *********************************************************************************** */
/**
 * @brief Computes the Euclidean distance between two 3D points.
 *
 * This function wraps d2xixj_lf3 to obtain the full Euclidean norm:
 *	 ||x1 - x2|| = sqrt(dot(x1 - x2, x1 - x2))
 *
 * @param x1 is a double 3x1 array
 * @param x2 is a double 3x1 array
 * @return double  Euclidean distance ||x1 - x2||
 */
double dxixj_lf3(double *x1, double *x2) { 
	
	return sqrt(d2xixj_lf3(x1, x2));
}
