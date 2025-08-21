/* *********************************************************************************** */
/* ---------------------------------- ADT files -------------------------------------- */
/* *********************************************************************************** */
int num_lines_file(const char filename[]);
/* *********************************************************************************** */
int num_columns_file_lines(const char filename[], char sep_char);
/* *********************************************************************************** */
double **file_2_mat_lf(const char filename[], char sep_char);
/* *********************************************************************************** */
dcinputfile *dcInputFile_2_dcDataVector(const char *filename);
/* *********************************************************************************** */
void read_input_file(const char *filename, char **structure_id, char **structure_chain, char **method, char **fname_dc, char **fname_T0, char **fname_X0, double *timeLimit, double *distanceResolution, double *angularResolution, int *numOfSols, int *sampleSize, double *differenceThreshold);
/* *********************************************************************************** */
void save_protein_model_PDBformat(const char filename[], proteinstructure *protein, int n, const char model[], const char structure_chain[]);
/* *********************************************************************************** */
void save_given_protein_PDBformat(proteinstructure *protein, int n, double **X, const char structure_id[], const char structure_chain[], const char filename[]);
/* *********************************************************************************** */
void save_protein_all_models_PDBformat(const char method[], proteinstructure *protein, int n, double ***X, int numSols, const char structure_id[], const char structure_chain[], const char filename[]);
/* *********************************************************************************** */
