#pragma once

#ifdef METADA_WITH_BUFR

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Open a BUFR file for processing
 *
 * @details Opens a BUFR file for reading or writing using the specified unit
 * number. The function handles the Fortran/C interoperability layer for BUFR
 * file operations.
 *
 * @param filename Path to the BUFR file to open
 * @param unit_num Logical unit number to associate with the file
 * @param status Output status code (0=success, non-zero=error)
 * @throws std::runtime_error If file cannot be opened
 */
void open_bufr_file_(const char* filename, int unit_num, int* status);

/**
 * @brief Read and process the next station report from a PREPBUFR file
 *
 * @details Reads a single station report from a PREPBUFR file, including all
 * associated levels and variables. The function handles the Fortran/C
 * interoperability layer for BUFR data reading.
 *
 * @param lunit Logical unit number for the file
 * @param subset Output buffer for subset name
 * @param idate Output date value
 * @param hdr_out Output array for header information (NHR8PM elements)
 * @param evns_out Output 4D array for events data (MXR8PM x MXR8LV x MXR8VN x
 * MXR8VT)
 * @param nlev_out Output number of levels in the report
 * @param iret Return code (0=OK, 1=last subset, -1=EOF)
 * @param subset_len Length of the subset string buffer
 * @throws std::runtime_error If reading fails
 */
void readpb_(int* lunit, char* subset, int* idate, double* hdr_out,
             double* evns_out, int* nlev_out, int* iret, int subset_len);

/**
 * @brief Close a BUFR file
 *
 * @details Closes a previously opened BUFR file and releases associated
 * resources. The function handles the Fortran/C interoperability layer for BUFR
 * file cleanup.
 *
 * @param unit_num Logical unit number of the file to close
 * @throws std::runtime_error If file cannot be closed
 */
void close_bufr_file_(int unit_num);

/**
 * @brief Process virtual temperature observations
 *
 * @details Handles virtual temperature observations according to PREPBUFR
 * Table 14. For VIRTMP program code 8 with reason code 3, the observation is
 * skipped. For VIRTMP program code 8 with any other reason code, steps down the
 * event stack index to find sensible temperature.
 *
 * @param lev Level index in the event stack
 * @param k Variable type index
 * @param idx Output index of the virtual temperature event (0 if not found)
 * @param flag Output flag (-1 to skip observation, 1 to process)
 */
void virtmp_(int* lev, int* k, int* idx, int* flag);

#ifdef __cplusplus
}
#endif

#endif  // METADA_WITH_BUFR