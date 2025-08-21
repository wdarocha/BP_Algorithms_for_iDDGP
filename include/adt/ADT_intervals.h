/* *********************************************************************************** */
/**
 * @brief Structure to represent cosine intervals for positive and negative sine regions.
 *
 * This structure stores up to 2 intervals for both the positive and negative projections
 * of cosine, typically representing symmetric angular intervals in the cosine-sine plane.
 *
 * - npos_intervals: Number of intervals in the positive sine region (quadrants 1 and 2).
 * - nneg_intervals: Number of intervals in the negative sine region (quadrants 3 and 4).
 * - pos_interval: Array [2][2] of [min, max] intervals for cosine in the positive sine region.
 * - neg_interval: Array [2][2] of [min, max] intervals for cosine in the negative sine region.
 */
typedef struct {
	int npos_intervals;
	int nneg_intervals;
	double pos_interval[4][2];  // [interval_index][min, max]
	double neg_interval[4][2];  // [interval_index][min, max]
} type_matrix_interval;
/* *********************************************************************************** */
void set_interval(double matrix_interval[][2], int index, double a, double b);
/* *********************************************************************************** */
void build_simple_symmetric_interval(type_matrix_interval *interval, double a, double b);
/* *********************************************************************************** */
void simple_intervals_intersection(const double *A, const double *B, double tolerance, double *C, int *CisEmpty);
/* *********************************************************************************** */
void intersect_interval_sets(const double A[][2], int nA, const double B[][2], int nB, double tolerance, double C[][2], int *num_intervals);
/* *********************************************************************************** */
void matrix_intervals_intersection(type_matrix_interval *A, const type_matrix_interval *B, double tolerance);
/* *********************************************************************************** */
int merge_intervals(double matrix_interval[][2], int num_intervals, double tolerance);
/* *********************************************************************************** */
int sampling_simple_interval_uniformly_with_tolerance(double *sample, double *interval, int sampleSize, double tolerance);
/* *********************************************************************************** */
