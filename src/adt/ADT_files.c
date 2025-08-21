#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "adt/ADT_math.h"
#include "adt/ADT_vectors.h"
#include "adt/ADT_matrices.h"
#include "adt/ADT_files.h"
#include "adt/ADT_strings.h"
#define CF_DEG2RAD 0.01745329252
#define MAX_LINE_LENGTH 1024
/* *********************************************************************************** */
/* ---------------------------------- ADT files -------------------------------------- */
/* *********************************************************************************** */
/**
 * @brief Efficiently counts all lines in a text file, including empty and comment lines.
 *
 * This function uses buffered binary reading with fread to efficiently count newline characters.
 * It works well even for large files and handles cases where the last line may not end with a newline.
 *
 * @param filename The name of the file to read.
 * @return int The total number of lines in the file.
 */
int num_lines_file(const char filename[]) {
	
	FILE *fp = fopen(filename, "rb"); // Use binary mode for portability and speed
	if (fp == NULL) {
		fprintf(stderr, "Error: could not open file '%s' for reading.\n", filename);
		exit(EXIT_FAILURE);
	}

	const size_t BUFFER_SIZE = 1 << 16;  // 64 KB
	char buffer[BUFFER_SIZE];
	size_t bytes_read;
	int line_count = 0;
	int last_char = 0;

	while ((bytes_read = fread(buffer, 1, BUFFER_SIZE, fp)) > 0)
		for (size_t i = 0; i < bytes_read; i++) {
			if (buffer[i] == '\n')
				line_count++;
			last_char = buffer[i];
		}

	if (ferror(fp)) {
		fprintf(stderr, "Error reading file '%s'.\n", filename);
		fclose(fp);
		exit(EXIT_FAILURE);
	}

	fclose(fp);

	// If file is not empty and last character is not '\n', count the last line
	if (line_count == 0 && last_char != '\n' && last_char != 0)
		line_count = 1;
	else if (last_char != '\n')
		line_count++;
	
	return line_count;
}
/* *********************************************************************************** */
/**
 * @brief Counts columns in the *first logical line* (until first '\n') of a file.
 *        Handles very long lines by dynamically reading in chunks.
 *
 * @param filename The name of the file to read.
 * @param sep_char The separator character (e.g., ',', '\t').
 * @return int Number of columns (separators + 1), or 0 if empty.
 */
int num_columns_file_lines(const char filename[], char sep_char) {
	
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Error: could not open file '%s' for reading.\n", filename);
		exit(EXIT_FAILURE);
	}

	const size_t CHUNK = 8192;
	char buffer[CHUNK];
	int count = 0;
	int column_counted = 0;

	while (fgets(buffer, CHUNK, fp)) {
		size_t len = strlen(buffer);

		// Count separators in the chunk
		for (size_t i = 0; i < len; i++)
			if (buffer[i] == sep_char)
				count++;
			else if (buffer[i] == '\n') {
				column_counted = 1;
				break;
			}

		if (column_counted || feof(fp))
			break;
	}
	
	if (ferror(fp)) {
		fprintf(stderr, "Error reading file '%s'.\n", filename);
		fclose(fp);
		exit(EXIT_FAILURE);
	}

	fclose(fp);

	// If line was empty (no '\n' and nothing read)
	if (count == 0 && !column_counted)
		return 0;

	return count + 1;
}
/* *********************************************************************************** */
/**
 * @brief Reads a matrix of doubles from a text file using a specified column separator.
 *
 * This function determines the number of rows and columns in the input file, allocates
 * memory for the matrix, and fills it by parsing each line. It provides detailed error
 * messages if the file format is inconsistent or unreadable.
 *
 * @param filename   The name of the input file.
 * @param sep_char   The character used to separate columns (e.g., ',', ' ', '\t').
 * @return double**  A dynamically allocated matrix [m x n] containing the file data.
 */
double **file_2_mat_lf(const char filename[], char sep_char) {
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: could not open file '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    
    int m = num_lines_file(filename);
    int n = num_columns_file_lines(filename, sep_char);
    double **T = alloc_mat_lf(m, n);
    
    char buffer[MAX_LINE_LENGTH];
    char sep_str[2] = {sep_char, '\0'};
    int i = 0;
    
    while (fgets(buffer, sizeof(buffer), fp) != NULL && i < m) {
        size_t len = strlen(buffer);
        while (len > 0 && (buffer[len - 1] == '\n' || buffer[len - 1] == '\r'))
            buffer[--len] = '\0';
        if (len == 0) continue;
        
        char *token = strtok(buffer, sep_str);
        int j = 0;
        while (token && j < n) {
            if (sscanf(token, "%lf", &T[i][j]) != 1) {
                fprintf(stderr, "Error: invalid double at line %d, column %d in '%s'.\n", i + 1, j + 1, filename);
                fclose(fp);
                dealloc_mat_lf(T, m);
                exit(EXIT_FAILURE);
            }
            token = strtok(NULL, sep_str);
            j++;
        }
        if (j != n) {
            fprintf(stderr, "Error: line %d has %d columns (expected %d) in '%s'.\n", i + 1, j, n, filename);
            fclose(fp);
            dealloc_mat_lf(T, m);
            exit(EXIT_FAILURE);
        }
        i++;
    }
    
    if (i != m) {
        fprintf(stderr, "Error: expected %d lines, but read %d in '%s'.\n", m, i, filename);
        fclose(fp);
        dealloc_mat_lf(T, m);
        exit(EXIT_FAILURE);
    }
    
    fclose(fp);
    return T;
}
/* *********************************************************************************** */
/**
 * @brief Parses a file with distance constraints into a vector of dcinputfile structs.
 *
 * Each line of the file is expected to contain 10 space-separated fields:
 * i j ri rj dl du ai aj rti rtj
 * Fields may be separated by multiple spaces, which are normalized internally.
 *
 * @param filename Path to the input file.
 * @return Pointer to a dynamically allocated array of dcinputfile structs.
 */
dcinputfile *dcInputFile_2_dcDataVector(const char *filename) {
	
	// Count the number of lines in the input file
	int m = num_lines_file(filename);

	// Allocate memory for the result vector
	dcinputfile *dcvector = alloc_vec_dcInputFile(m);

	// Open the input file for reading
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: could not open file '%s' for reading.\n", filename);
		exit(EXIT_FAILURE);
	}

	char line[MAX_LINE_LENGTH];  // Buffer to hold each line of the file
	int i = 0;

	// Read the file line by line
	while (fgets(line, sizeof(line), fp) != NULL && i < m) {
		// Normalize the line: convert multiple spaces to single space
		replace_multiple_spaces(line);

		// Token array to hold the split parts of the line
		char *tokens[10];
		char *token = strtok(line, " ");
		int ntok = 0;

		// Tokenize the line by spaces and store up to 10 fields
		while (token && ntok < 10) {
			tokens[ntok++] = token;
			token = strtok(NULL, " ");
		}

		// If the line doesn't contain exactly 10 fields, report an error
		if (ntok != 10) {
			fprintf(stderr, "Parsing error on line %d of file '%s': expected 10 tokens, got %d\n", i + 1, filename, ntok);
			free(dcvector);
			fclose(fp);
			exit(EXIT_FAILURE);
		}

		// Convert and store each field into the struct
		dcvector[i].i   = atoi(tokens[0]);
		dcvector[i].j   = atoi(tokens[1]);
		dcvector[i].ri  = atoi(tokens[2]);
		dcvector[i].rj  = atoi(tokens[3]);
		dcvector[i].dl  = strtod(tokens[4], NULL);
		dcvector[i].du  = strtod(tokens[5], NULL);

		// Copy atom names and residue types into fixed-size char arrays
		strncpy(dcvector[i].ai,  tokens[6], sizeof(dcvector[i].ai) - 1);
		dcvector[i].ai[sizeof(dcvector[i].ai) - 1] = '\0';

		strncpy(dcvector[i].aj,  tokens[7], sizeof(dcvector[i].aj) - 1);
		dcvector[i].aj[sizeof(dcvector[i].aj) - 1] = '\0';

		strncpy(dcvector[i].rti, tokens[8], sizeof(dcvector[i].rti) - 1);
		dcvector[i].rti[sizeof(dcvector[i].rti) - 1] = '\0';

		strncpy(dcvector[i].rtj, tokens[9], sizeof(dcvector[i].rtj) - 1);
		dcvector[i].rtj[sizeof(dcvector[i].rtj) - 1] = '\0';

		i++;
	}

	fclose(fp);
	return dcvector;
}
/* *********************************************************************************** */
/**
 * @brief Reads and parses a 11-line structured input file with key: value format.
 *
 * Each line of the file is expected to follow the format:
 *   <field description>: <value>
 * The function extracts values from the lines and stores them in the corresponding variables.
 *
 * Special handling is used to correctly parse fields whose description contains colons,
 * such as "time limit (days-hours:minutes:seconds): 0-00:00:03".
 *
 * @param filename		Path to the input file.
 * @param structure_id		Output: structure ID string.
 * @param structure_chain	Output: chain ID string.
 * @param method		Output: method string.
 * @param fname_dc		Output: path to distance constraints file.
 * @param fname_T0		Output: path to torsion angles file.
 * @param fname_X0		Output: path to reference XYZ structure file.
 * @param timeLimit		Output: total time limit in seconds.
 * @param distanceResolution	Output: distance resolution in Ångströms.
 * @param angularResolution	Output: angular resolution in radians.
 * @param numOfSols		Output: number of solutions to find (0 = all).
 * @param sampleSize		Output: number of samples to generate per torsion.
 * @param differenceThreshold	Output: RMSD threshold (in Ångströms) used to decide whether two solutions are considered different.
 */
void read_input_file(const char *filename, char **structure_id, char **structure_chain, char **method, char **fname_dc, char **fname_T0, char **fname_X0, double *timeLimit, double *distanceResolution, double *angularResolution, int *numOfSols, int *sampleSize, double *differenceThreshold) {
	
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: could not open file '%s' for reading.\n", filename);
		exit(EXIT_FAILURE);
	}
	
	int numLines = num_lines_file(filename);
	if (numLines != 12) {
		fprintf(stderr, "Error: input file must contain exactly 11 lines (found %d).\n", numLines);
		fclose(fp);
		exit(EXIT_FAILURE);
	}
	
	char line[MAX_LINE_LENGTH];  // Buffer to hold each line
	char *content[12] = {NULL};  // Array to store extracted values
	int i;
	
	// Loop over each line and extract the value part
	for (i = 0; i < 12; i++) {
		if (fgets(line, MAX_LINE_LENGTH, fp) == NULL) {
			fprintf(stderr, "Error: unexpected end of file at line %d in input file.\n", i + 1);
			fclose(fp);
			exit(EXIT_FAILURE);
		}

		// Try to find the colon that separates field name from value.
		// If the line contains '(', we assume the field label ends at the colon after ')'
		char *colon_position = NULL;
		char *closing_paren = strrchr(line, ')');
		
		if (closing_paren)
			colon_position = strchr(closing_paren, ':');
		else
			colon_position = strchr(line, ':');  // fallback: first colon
		
		if (!colon_position) {
			fprintf(stderr, "Error: line %d is malformed (missing ':') in input file.\n", i + 1);
			fclose(fp);
			exit(EXIT_FAILURE);
		}
		
		// Skip the colon and any whitespace after it
		char *value = colon_position + 1;
		while (*value == ' ' || *value == '\t')
			value++;
		
		// Duplicate the value string and strip newline
		content[i] = custom_strdup(value);
		if (!content[i]) {
			fprintf(stderr, "Error: memory allocation failed at line %d in input file.\n", i + 1);
			fclose(fp);
			exit(EXIT_FAILURE);
		}
		content[i][strcspn(content[i], "\n")] = '\0';
	}
	
	fclose(fp);  // Done reading file
	
	// Copy strings into user-provided pointers
	*structure_id     = custom_strdup(content[0]);
	*structure_chain  = custom_strdup(content[1]);
	*method           = custom_strdup(content[2]);
	*fname_dc         = custom_strdup(content[3]);
	*fname_T0         = custom_strdup(content[4]);
	*fname_X0         = custom_strdup(content[5]);
	
	// Parse time limit in the format: days-hours:minutes:seconds
	double days, hours, minutes, seconds;
	if (sscanf(content[6], "%lf-%lf:%lf:%lf", &days, &hours, &minutes, &seconds) != 4) {
		fprintf(stderr, "Error: invalid time format on line 7 of input file (expected days-hours:minutes:seconds).\n");
		exit(EXIT_FAILURE);
	}
	*timeLimit = 86400.0 * days + 3600.0 * hours + 60.0 * minutes + seconds;
	
	// Parse remaining numeric fields
	*distanceResolution  = strtod(content[7], NULL);
	*angularResolution   = strtod(content[8], NULL) * CF_DEG2RAD;  // degrees -> radians
	*numOfSols           = (int)strtod(content[9], NULL);
	*sampleSize          = (int)strtod(content[10], NULL);
	*differenceThreshold = strtod(content[11], NULL);
	
	// Free temporary buffers
	for (i = 0; i < 12; i++)
		free(content[i]);
}
/* *********************************************************************************** */
/**
 * @brief Saves a single protein model in PDB format.
 *
 * This function appends a MODEL block to a PDB file, containing ATOM and TER entries.
 *
 * @param filename         Output file name.
 * @param protein          Array of proteinstructure with atom/residue data.
 * @param n                Number of atoms in the structure.
 * @param model            Model number (as string).
 * @param structure_chain  Chain identifier (e.g., "A").
 */
void save_protein_model_PDBformat(const char filename[], proteinstructure *protein, int n, const char model[], const char structure_chain[]) {
	
	FILE *fp = fopen(filename, "a");
	if (!fp) {
		fprintf(stderr, "Error: could not open file '%s' for writing.\n", filename);
		exit(EXIT_FAILURE);
	}
	
	fprintf(fp, "MODEL        %6s\n", model);
	
	int k;
	char name[5];        // Atom name field
	char resName_fmt[4]; // Residue name field
	char element[3];     // Element field
	const char *resName, *atomRaw;
	double x, y, z;
	int serial, resSeq;
	double occupancy = 1.00;
	double tempFactor = 1.00;
	char chainID = structure_chain[0];
	char altLoc = ' ';
	char iCode = ' ';
	char charge[3] = "  ";
	
	// Write ATOM lines
	for (k = 0; k < n; k++) {
		serial   = protein[k].v_k;
		resName  = protein[k].res_name;
		atomRaw  = protein[k].atom_name;
		resSeq   = protein[k].res_seq;
		x        = protein[k].x;
		y        = protein[k].y;
		z        = protein[k].z;

		// Prepare atom name (right-aligned, max 4 chars)
		memset(name, ' ', 4);
		name[4] = '\0';
		size_t len = strlen(atomRaw);
		if (len > 4) {
			fprintf(stderr, "Warning: atom name '%s' truncated to 4 characters (atom %d)\n", atomRaw, serial);
			strncpy(name, atomRaw, 4);
		}
		else
			strncpy(&name[4 - len], atomRaw, len);
		
		// Prepare resName and element
		snprintf(resName_fmt, sizeof(resName_fmt), "%3s", resName);
		element[0] = atomRaw[0];
		element[1] = '\0';
		
		fprintf(fp, "%-6s%5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
		        "ATOM", serial, name, altLoc, resName_fmt, chainID, resSeq, iCode,
		        x, y, z, occupancy, tempFactor, element, charge);
	}
	
	// TER line
	fprintf(fp, "%-6s%5d %-4s%c%3s %c%4d%c\n",
	        "TER", n + 1, "    ", ' ', protein[n - 1].res_name, chainID, protein[n - 1].res_seq, ' ');
	
	fprintf(fp, "ENDMDL\n");
	
	fclose(fp);
}
/* *********************************************************************************** */
/**
 * @brief Saves all protein models in PDB format.
 *
 * This function writes the PDB HEADER/TITLE/REMARK and then appends each model.
 *
 * @param method           Method name (e.g., "itbp", "iabp", or "PDB").
 * @param protein          Array of proteinstructure with atom/residue data.
 * @param n                Number of atoms in the structure.
 * @param X                3D array containing all solutions [numSols][n][3].
 * @param numSols          Number of solutions to save.
 * @param structure_id     PDB structure ID.
 * @param structure_chain  Chain identifier (e.g., "A").
 * @param filename         Output filename.
 */
void save_protein_all_models_PDBformat(const char method[], proteinstructure *protein, int n, double ***X, int numSols, const char structure_id[], const char structure_chain[], const char filename[], int is_x0) {
	
	// Write HEADER and TITLE block once
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		fprintf(stderr, "Error: could not open file '%s' for writing.\n", filename);
		exit(EXIT_FAILURE);
	}
	
	// Current date YYYY-MM-DD
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	char date[11];
	snprintf(date, sizeof(date), "%04d-%02d-%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
	
	fprintf(fp, "HEADER    %62s   %4s\n", date, structure_id);
	if (strcmp(method, "PDB") != 0)
		fprintf(fp, "TITLE     %-4s solution\n", method);
	else
		fprintf(fp, "TITLE     PDB structure\n");
	
	fprintf(fp, "REMARK   2\n");
	fprintf(fp, "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n");
	
	fclose(fp);  // Done writing header
	
	// Write each model
	char model[12];
	for (int k = (1 - is_x0); k < numSols; k++) {
		// Update coordinates in the protein structure
		for (int i = 0; i < n; i++) {
			protein[i].x = X[k][i][0];
			protein[i].y = X[k][i][1];
			protein[i].z = X[k][i][2];
		}
		
		// Model number as string
		snprintf(model, sizeof(model), "%d", maxAB_d(k, 1));
		
		// Save this model to file (append mode)
		save_protein_model_PDBformat(filename, protein, n, model, structure_chain);
	}
	
	FILE *fp2 = fopen(filename, "a");
	if (!fp2) {
		fprintf(stderr, "Error: could not open file '%s' for writing.\n", filename);
		exit(EXIT_FAILURE);
	}
	
	fprintf(fp2, "END\n");
	
	fclose(fp2);
}
