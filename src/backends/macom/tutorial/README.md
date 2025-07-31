# MACOM Backend Tutorial

This directory contains tutorial configuration files for the MACOM (Marine Circulation and Ocean Model) backend.

## Configuration Files

### ETKF Configuration (`etkf.yaml`)
- **Purpose**: Ensemble Transform Kalman Filter data assimilation
- **Features**:
  - 2 ensemble members with temperature variable
  - Sea surface temperature observations
  - Localization radius: 100.0 km
  - Inflation factor: 1.0
  - NetCDF output format

### LETKF Configuration (`letkf.yaml`)
- **Purpose**: Local Ensemble Transform Kalman Filter data assimilation
- **Features**:
  - 2 ensemble members with temperature variable
  - Sea surface temperature observations
  - Localization radius: 100.0 km
  - Inflation factor: 1.0
  - NetCDF output format

### ENKF Configuration (`enkf.yaml`)
- **Purpose**: Ensemble Kalman Filter data assimilation
- **Features**:
  - 2 ensemble members with temperature, salinity, and velocity variables
  - Sea surface temperature and Argo profile observations
  - Inflation factor: 1.0
  - Multiplicative inflation method
  - NetCDF output format

### LWEnKF Configuration (`lwenkf.yaml`)
- **Purpose**: Local Weighted Ensemble Kalman Filter data assimilation
- **Features**:
  - 2 ensemble members with temperature, salinity, and velocity variables
  - Sea surface temperature and Argo profile observations
  - Localization radius: 100.0 km
  - Gaussian localization function
  - Adaptive weighting scheme
  - Inflation factor: 1.1
  - NetCDF output format

### Particle Filter Configuration (`particle_filter.yaml`)
- **Purpose**: Particle Filter data assimilation
- **Features**:
  - 2 ensemble members with temperature, salinity, and velocity variables
  - Sea surface temperature and Argo profile observations
  - Systematic resampling method
  - Jittering enabled for particle diversity
  - Effective sample size threshold: 50
  - NetCDF output format

### TL/AD Checks Configuration (`tl_ad_checks.yaml`)
- **Purpose**: Tangent Linear and Adjoint consistency checks
- **Features**:
  - Single state with temperature variable for testing
  - Single observation point for gradient checks
  - Multiple finite difference epsilons for testing
  - LBFGS minimization algorithm
  - Diagonal background error covariance
  - NetCDF output format

### Variational Configuration (`variational.yaml`)
- **Purpose**: Variational data assimilation (3DVAR/4DVAR/FGAT)
- **Features**:
  - Single background state with full ocean variables
  - Single observation point for analysis
  - LBFGS minimization algorithm
  - Configurable variational type (3DVAR/4DVAR/FGAT)
  - Localization enabled with 100km radius
  - NetCDF output format with trajectory saving

### Forecast Configuration (`forecast.yaml`)
- **Purpose**: Ocean model forecast simulation
- **Features**:
  - Single state with temperature, salinity, and velocity variables
  - 24-hour forecast length
  - 360-second time step
  - NetCDF output format

## Variable Definitions

- **t**: Temperature (°C)
- **s**: Salinity (PSU)
- **u**: Zonal velocity (m/s)
- **v**: Meridional velocity (m/s)
- **w**: Vertical velocity (m/s)

## Observation Types

### Sea Surface Temperature (SST)
- **Source**: Satellite observations
- **Format**: Geographic coordinates (lat, lon, level, value, error)
- **Error**: 0.5°C standard deviation
- **Missing value**: -999.0

### Argo Profiles
- **Source**: Argo float observations
- **Format**: Geographic coordinates (lat, lon, depth, value, error)
- **Variables**: Temperature (error: 0.3°C) and Salinity (error: 0.1 PSU)
- **Missing value**: -999.0

## File Structure Requirements

### Input Files
- **Grid file**: NetCDF format containing grid information
- **State files**: NetCDF format containing model state variables
- **Observation files**: Text format with geographic coordinates

### Sample Observation Files
- **sample_obs.txt**: Single point temperature observation for testing
  - Format: lat lon level value error
  - Single observation at (26.60°N, 124.40°E, 5.0m depth)
  - Temperature: 15.2°C with 0.4°C error

### Output Files
- **Analysis results**: NetCDF format for most algorithms
- **Detailed results**: Available for particle filter

## Usage

1. Ensure all input files are available at the specified paths
2. Modify file paths in the configuration files as needed
3. Run the corresponding application with the configuration file:
   ```bash
   ./etkf config_file.yaml
   ./letkf config_file.yaml
   ./enkf config_file.yaml
   ./lwenkf config_file.yaml
   ./particle_filter config_file.yaml
   ./forecast config_file.yaml
   ```

## Configuration Parameters

### Common Parameters
- **inflation**: Ensemble spread inflation factor
- **localization_radius**: Distance for localization (km)
- **format**: Output file format (nc, txt)

### Algorithm-Specific Parameters
- **ETKF/LETKF**: Basic ensemble Kalman filter variants
- **ENKF**: Standard ensemble Kalman filter
- **LWEnKF**: Local weighted ensemble Kalman filter with adaptive weighting
- **Particle Filter**: Particle-based data assimilation with resampling 