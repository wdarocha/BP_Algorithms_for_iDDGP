#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <lapacke.h>
#include <cblas.h>
#include "adt/ADT_math.h"
#include "adt/ADT_vectors.h"
#include "adt/ADT_intervals.h"
#include "adt/ADT_matrices.h"
#include "adt/ADT_files.h"
#include "adt/ADT_strings.h"
#include "adt/ADT_DDGP.h"
#include "algorithms/common.h"
#include "algorithms/iBP.h"
#include "algorithms/iTBP.h"
#include "project_version.h"
#define CF_DEG2RAD 0.01745329252
#define MAX_LINE_LENGTH 1024

int main(int argc, char *argv[]){
	
	if(argc < 3){
		printf("Usage: %s "
		"<%%s: input_filename_path> "
		"<%%s: output_folder_path> \n", argv[0]);
		return 0;
	}
	
	char *fname = argv[1];
	char *outputFolder = argv[2];
	char *structure_id = NULL;
	char *structure_chain = NULL;
	char *method = NULL;
	char *instanceFile = NULL;
	char *cliquesFile = NULL;
	char *initialStructureFile = NULL;
	double timeLimit = 0.0;
	double distanceResolution = 0.0;
	double angularResolution = 0.0;
	int numOfSols = 0; 
	int sampleSize = 0;
	double solutionDifferenceThreshold = 0.0;
	
	read_input_file(fname, &structure_id, &structure_chain, &method, &instanceFile, &cliquesFile, &initialStructureFile, &timeLimit, &distanceResolution, &angularResolution, &numOfSols, &sampleSize, &solutionDifferenceThreshold);
		
	char aux_str[MAX_LINE_LENGTH];
	char filepath[MAX_LINE_LENGTH];
	
	int m, n, k;
	int *kv, *tauSign;
	int **cliques;
	int vec_cliques[] = {1, 2, 3, 4};
	
	double *tauSign_lf, *givenTau, *givenTauDeviation;
	double **T0, **cliques_lf, **distanceConstraints;
	double ***instanceAdjlist, ***discretizationEdges_2;
	
	dcinputfile *dc_vec;
	
	prune_edges_set *pruneEdges_2;
	
	proteinstructure *protein;
	
	run_metrics runMetrics;
	
	m = num_lines_file(instanceFile);
	
	dc_vec = dcInputFile_2_dcDataVector(instanceFile);
	
	dcinputfile_2_instance(dc_vec, &distanceConstraints, &instanceAdjlist, &protein, &kv, m, &n);
	
	free(dc_vec);
	
	if ((numOfSols == 0) || (numOfSols > 20))
		numOfSols = 20;
	
	double ***allSolutions = alloc_3Darray_lf(numOfSols + 1, n, 3); // 1-based array in the fist dimension
	
	// if the PDB structure is known it is saved as the first solution
	if (initialStructureFile != NULL) {
		double **X0 = file_2_mat_lf(initialStructureFile, ' ');

		// Save X0 in allSolutions[0]
		mat_e_mat_lf(allSolutions[0], X0, n, 3);

		// Build the output file path
		snprintf(filepath, sizeof(filepath), "%s/%s_%s.pdb", outputFolder, structure_id, structure_chain);

		// Save a pdb file with the given pdb protein chain
		save_given_protein_PDBformat(protein, n, X0, structure_id, structure_chain, filepath);
		// Free the matrix read from file
		dealloc_mat_lf(X0, n);
	}
	
	if(cliquesFile != NULL)
		T0 = file_2_mat_lf(cliquesFile, ' ');
	else{
		printf("Error: the path to the cliques file could not be found\n");
		return 0;
	}
	
	cliques_lf = get_columns_mat_lf(T0, n, vec_cliques, 4);
	cliques = mat_lf_2_mat_d(cliques_lf, n, 4);
	dealloc_mat_lf(cliques_lf, n);
	
	tauSign_lf = get_column_mat_lf(T0, n, 5);
	tauSign = vec_lf_2_vec_d(tauSign_lf, n);
	free(tauSign_lf);
	
	givenTau = get_column_mat_lf(T0, n, 6);
	l_t_vec_lf(givenTau, CF_DEG2RAD, givenTau, n);
	
	givenTauDeviation = get_column_mat_lf(T0, n, 7);
	l_t_vec_lf(givenTauDeviation, CF_DEG2RAD, givenTauDeviation, n);
	
	dealloc_mat_lf(T0, n);
	
	adjlist_2_graph_parameters(instanceAdjlist, kv, n, cliques, &discretizationEdges_2, &pruneEdges_2);
	
	dealloc_3Darray_vec_lf(instanceAdjlist, kv, n);
	free(kv);
	
	free(instanceFile);
	free(initialStructureFile);
	free(cliquesFile);
	
	if (strcmp(method, "ibp") == 0) {
		printf("****************** iBP version %s (%s) ******************\n",  ibp_version(),  ibp_release_date());
		iBP(n, discretizationEdges_2, pruneEdges_2, sampleSize, distanceResolution, angularResolution, timeLimit, numOfSols, &runMetrics, distanceConstraints, m, allSolutions, solutionDifferenceThreshold);
	}
	else if (strcmp(method, "iabp") == 0) {
		printf("****************** iABP version %s (%s) ******************\n",  itbp_version(),  itbp_release_date());
		// In iABP, all torsion signs are assumed to be zero
		free(tauSign);
		tauSign = zeros_vec_d(n);
		iTBP(n, discretizationEdges_2, pruneEdges_2, sampleSize, angularResolution, timeLimit, numOfSols, tauSign, givenTau, givenTauDeviation, &runMetrics, distanceConstraints, m, allSolutions, solutionDifferenceThreshold);
	}
	else if (strcmp(method, "itbp") == 0) {
		printf("****************** iTBP version %s (%s) ******************\n",  itbp_version(),  itbp_release_date());
		iTBP(n, discretizationEdges_2, pruneEdges_2, sampleSize, angularResolution, timeLimit, numOfSols, tauSign, givenTau, givenTauDeviation, &runMetrics, distanceConstraints, m, allSolutions, solutionDifferenceThreshold);
	}
	else {
		printf("ERROR: The selected method is invalid -- expected 'ibp', 'iabp', or 'itbp'.\n");
		return 0;
	}
	
	strcpy(aux_str, outputFolder);
	FILE *fp = fopen(strcat(aux_str, "/results.txt"), "w");
	if(!fp){
		printf("Error: file reading error.\n");
		exit(1);
	}
	
	fprintf(fp, "CPU time = %.8lf\n", runMetrics.cpu_time_used);
	fprintf(fp, "Last embedded vertex = %d/%d\n", runMetrics.lev, n);
	fprintf(fp, "Number of embedded vertices = %.0lf\n", runMetrics.noev);
	fprintf(fp, "Number of solutions found = %ld\n", runMetrics.nosf);
	fprintf(fp, "Number of considered solutions = %d\n", runMetrics.nocs);
	
	if(runMetrics.nosf > 0) {
		fprintf(fp, "maximum MDE = %.8lf\n", runMetrics.maxMDE);
		fprintf(fp, "maximum LDE = %.8lf\n", runMetrics.maxLDE);
		fprintf(fp, "minimum RMSD = %.8lf\n", runMetrics.minRMSD);
		
		snprintf(filepath, sizeof(filepath), "%s/%s.pdb", outputFolder, structure_id);
		save_protein_all_models_PDBformat(method, protein, n, allSolutions, runMetrics.nocs, structure_id, structure_chain, filepath);
	}
	else {
		fprintf(fp, "maximum MDE = ---\n");
		fprintf(fp, "maximum LDE = ---\n");
		fprintf(fp, "minimum RMSD = ---\n");
	}
	
	fclose(fp);
		
	free(structure_id);
	free(structure_chain);
	free(method);
	
	dealloc_mat_lf(distanceConstraints, m);
	free(protein);
	
	dealloc_3Darray_lf(allSolutions, numOfSols + 1, n);
	
	dealloc_mat_d(cliques, n);
	free(tauSign);
	free(givenTau);
	free(givenTauDeviation);
	
	dealloc_3Darray_lf(discretizationEdges_2, n, 3);
	for(k = 0; k < n; k++){
		if(pruneEdges_2[k].cardUkP > 0)
			dealloc_mat_lf(pruneEdges_2[k].precise, pruneEdges_2[k].cardUkP);
		if(pruneEdges_2[k].cardUkI > 0)
			dealloc_mat_lf(pruneEdges_2[k].interval, pruneEdges_2[k].cardUkI);
	}
	free(pruneEdges_2);
		
	return 0;
}
