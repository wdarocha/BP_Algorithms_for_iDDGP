typedef struct {
        double noev;
        int lev;
        long nosf;
        int nocs;
        double cpu_time_used;
        double maxMDE;
        double maxLDE;
        double minRMSD;
} run_metrics;
/* *********************************************************************************** */
void referential_x1_x2_x3(double **X, double ***discretizationEdges);
/* *********************************************************************************** */
int handle_solution_cycle(const char *method, long int *nosf, double *maxMDE, double *maxLDE, double *minRMSD, int *nocs, int GivenNumOfSols, double **Xr, int n, double ***allSolutions, double **distanceConstraints, int num_dc, int referenceSolutionIndex, double solutionDifferenceThreshold, int *i, int *exploredVertex, double **branches, int *branchNum, int twoSampleSize);
/* *********************************************************************************** */
void tree_backtracking(int *i, int *exploredVertex, double **branches, int *branchNum, int sampleSize);
/* *********************************************************************************** */
