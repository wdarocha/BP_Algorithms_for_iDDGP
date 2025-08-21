#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adt/ADT_strings.h"

#define MAX_LINE_LENGTH 1024
/* *********************************************************************************** */
/* ---------------------------------- ADT strings ------------------------------------ */
/* *********************************************************************************** */
/**
 * @brief Replaces multiple consecutive spaces in a string with a single space.
 *
 * This function modifies the input string in-place. It skips over extra spaces
 * and retains only one space between words.
 *
 * @param str The input string to process.
 */
void replace_multiple_spaces(char *str) {
	
	char result[MAX_LINE_LENGTH];
	int i = 0, j = 0;
	int in_space = 0;
	
	while (str[i] != '\0' && j < MAX_LINE_LENGTH - 1) {
		if (str[i] != ' ') {
			result[j++] = str[i];
			in_space = 0;
		}
		else if (!in_space) {  // Only write first space
			result[j++] = ' ';
			in_space = 1;
		}
		i++;
	}

	// Remove trailing space, if present
	if (j > 0 && result[j - 1] == ' ')
		j--;

	result[j] = '\0';
	strcpy(str, result);
}
/* *********************************************************************************** */
/**
 * @brief Duplicates a C string by allocating a new buffer and copying its contents.
 *
 * This function behaves like POSIX strdup(): it allocates memory for a copy of the string
 * and returns a pointer to it.
 *
 * @param str   The null-terminated string to duplicate.
 * @return char* A newly allocated copy of the string, or NULL on allocation failure.
 */
char *custom_strdup(const char *str) {
	
	if (str == NULL)
		return NULL;  // Optional: handle NULL input defensively

	size_t len = strlen(str) + 1;  // Include null terminator
	char *dup = (char*)malloc(len);

	if (dup != NULL)
		memcpy(dup, str, len);

	return dup;
}
