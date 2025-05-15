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
 * @param hdr_out Output array containing header information
 * @param evns_out Output 4D array containing events data (MXR8PM x MXR8LV x
 * MXR8VN x MXR8VT)
 * @param nlev_out Number of levels in the report
 * @param iret Return code (0=OK, 1=last subset, -1=EOF)
 * @param subset_len Length of the subset string buffer (passed implicitly by
 * Fortran)
 */
void readpb_(int* lunit, char* subset, int* idate, double* hdr_out,
             double* evns_out, int* nlev_out, int* iret, int subset_len);

/**
 * @brief Close a BUFR file
 *
 * @param lunit Logical unit number for the file
 */
void close_bufr_file_(int unit_num);

#ifdef __cplusplus
}
#endif