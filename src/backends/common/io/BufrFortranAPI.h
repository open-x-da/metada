#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Open a BUFR file for processing
 *
 * Note: While the Fortran OPENBF function has 3 visible parameters, the 4th
 * parameter is an implicit string length parameter that Fortran automatically
 * adds for any character string arguments. This is required for proper
 * C/Fortran interoperability.
 *
 * @param lunit Logical unit number for the file
 * @param io_method I/O method ('IN' for input, 'OUT' for output, etc.)
 * @param lundx Table unit number
 * @param io_method_len Length of the io_method string (passed implicitly by
 * Fortran)
 */

void open_bufr_file_(const char* filename, int unit_num, int* status);
/**
 * @brief Read and process the next station report from a PREPBUFR file
 *
 * Note: The subset_len parameter is an implicit string length parameter added
 * by Fortran.
 *
 * @param lunit Logical unit number for the file
 * @param subset Subset name output
 * @param idate Date output
 * @param iret Return code (0=OK, 1=last subset, -1=EOF)
 * @param subset_len Length of the subset string buffer (passed implicitly by
 * Fortran)
 */
void readpb_(int* lunit, char* subset, int* idate, int* iret, int subset_len);

/**
 * @brief Close a BUFR file
 *
 * @param lunit Logical unit number for the file
 */
void close_bufr_file_(int unit_num);

#ifdef __cplusplus
}
#endif