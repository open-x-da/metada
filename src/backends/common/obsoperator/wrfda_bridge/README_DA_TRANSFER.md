# WRFDA Domain Initialization: Leveraging `da_transfer_wrftoxb`

## Overview

The WRFDA system includes a comprehensive subroutine `da_transfer_wrftoxb` that properly transfers WRF model fields to the WRFDA first guess (background) structure. This subroutine is well-tested and handles all the complex transformations needed for data assimilation.

## Why Use `da_transfer_wrftoxb`?

The `da_transfer_wrftoxb` subroutine (`var/da/da_transfer_model/da_transfer_wrftoxb.inc`) performs:

### 1. **Thermodynamic Calculations**
- Computes full pressure from base state + perturbation
- Calculates air density from inverse specific volume
- Converts potential temperature to actual temperature
- Supports two hypsometric formulations

### 2. **Derived Diagnostic Fields**
- Relative humidity and saturation vapor pressure
- Dew point temperature
- Sea level pressure (SLP)
- Total precipitable water (TPW)
- GPS refractivity and zenith total delay (ZTD)
- Brightness temperature for satellite data

### 3. **Surface Variables**
- 10-meter winds (u10, v10)
- 2-meter temperature and humidity (t2, q2)
- Surface roughness length
- Surface regime classification

### 4. **Grid and Coordinate Handling**
- Proper Arakawa-C grid staggering
- Vertical coordinate transformations
- Grid box areas and map factors
- Halo exchanges for parallel processing

### 5. **Quality Control Setup**
- Background error covariance preparation
- Vertical inner product computation
- Cloud model time steps (if enabled)
- Lifting condensation level (LCL) for radar DA

## Implementation

### Current Approach (Simplified)

```cpp
// In WRFState.hpp: initializeWRFDADomain()
// - Manually constructs grid structure
// - Calls wrfda_construct_domain_from_arrays()
// - Limited field initialization
// - Hardcoded parameters
```

**Limitations:**
- Does not compute derived diagnostics
- Missing surface variable calculations
- No GPS refractivity or radiance support
- Incomplete thermodynamic setup

### Recommended Approach (Complete)

```fortran
! In wrfda_dispatch.F90
! New function: wrfda_init_domain_from_wrf_fields()
! 
! 1. Populates full WRF grid structure from C++ arrays
! 2. Calls da_transfer_wrftoxb(xbx, grid, config_flags)
! 3. Returns fully initialized WRFDA domain
```

**Advantages:**
- ✅ All derived fields computed correctly
- ✅ Consistent with WRFDA's internal calculations  
- ✅ Proper handling of grid staggering
- ✅ Support for all observation types
- ✅ Less code to maintain
- ✅ Automatically benefits from WRFDA updates

## Migration Path

### Phase 1: Add Required WRF Fields (Current)
- [x] Create `wrfda_init_domain_from_wrf_fields()` binding
- [x] Use full `moist` array (eliminates redundant `qvapor` parameter)
- [x] Read additional 3D WRF fields from NetCDF:
  - `W` - vertical velocity
  - `MU`, `MUB` - column dry air mass
  - `P`, `PB` - pressure perturbation and base state
  - `PH`, `PHB` - geopotential perturbation and base state
  - `T_INIT` - initial temperature field
  - `HGT` - terrain height
- [x] Read moisture species (QCLOUD, QRAIN, QICE, QSNOW, QGRAUPEL)
- [x] Load 1D vertical coordinate arrays:
  - `ZNU`, `ZNW` - eta values on mass/staggered levels
  - `DN`, `DNW` - delta eta values (computed if not available)
  - `RDNW`, `RDN` - inverse delta eta (computed if not available)
  - `P_TOP` - pressure at model top
- [ ] Surface fields (TSK, TSLB, SMOIS, etc.) - for future enhancement

### Phase 2: Update WRFState::initializeWRFDADomain()
```cpp
// Replace manual construction with:
int rc = wrfda_init_domain_from_wrf_fields(
    &nx, &ny, &nz,
    u_data, v_data, w_data, t_data,
    mu_data, mub_data, p_data, pb_data,
    ph_data, phb_data, xlat_data, xlong_data, ht_data,
    znu_data, znw_data, dn_data, dnw_data,
    rdnw_data, rdn_data, &p_top,
    t_init_data, moist_data, &num_moist, psfc_data,
    &start_year, &start_month, &start_day, &start_hour);
```

**Note**: `moist_data` is a 4D array containing all moisture species where:
- `moist(:,:,:,1)` = QVAPOR (water vapor)
- `moist(:,:,:,2)` = QCLOUD (cloud water) 
- `moist(:,:,:,3)` = QRAIN (rain)
- `moist(:,:,:,4+)` = QICE, QSNOW, QGRAUPEL (optional)

### Phase 3: Remove Deprecated Code
- Remove `wrfda_construct_domain_from_arrays()`
- Remove manual unstaggering code
- Remove hardcoded parameters

## Benefits

1. **Correctness**: Uses the same calculations as WRFDA internally
2. **Maintainability**: Single source of truth
3. **Completeness**: All observation operators work correctly
4. **Consistency**: Same results as standalone WRFDA
5. **Updates**: Automatically inherit WRFDA improvements

## References

- WRFDA source: `var/da/da_transfer_model/da_transfer_wrftoxb.inc`
- WRFDA documentation: WRFDA User's Guide Chapter 5
- Current binding: `src/backends/common/obsoperator/wrfda_bridge/wrfda_dispatch.F90`

## Status

- ✅ Function binding created
- ⏳ WRF field reading in progress
- ⏳ WRFState integration pending
- ⏳ Testing and validation needed

