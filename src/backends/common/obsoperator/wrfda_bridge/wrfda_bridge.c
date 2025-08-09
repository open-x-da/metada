// @file wrfda_bridge.c
// @brief C side of the bridge. Only ensures the Fortran symbols are linked.

#include <stddef.h>
#include "../WRFDAObsOperator_c_api.h"

// Declarations are provided by the Fortran object with bind(C) names, so this
// TU does not need to implement anything further. Building it in the library
// ensures the symbols are part of the archive and headers are installed.


