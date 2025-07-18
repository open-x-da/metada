# Ensemble Kalman Filter (EnKF) Algorithm

## Overview

The Ensemble Kalman Filter (EnKF) is a Monte Carlo implementation of the Kalman filter that represents the state distribution by a set of ensemble members. Unlike the ETKF which uses ensemble space transformations, the EnKF uses the traditional Kalman gain formulation with observation perturbations.

## Algorithm Description

The EnKF algorithm follows these steps:

1. **Forecast Step**: Build forecast quantities (mean and perturbations)
2. **Observation Space Propagation**: Apply observation operator to each ensemble member
3. **Innovation Computation**: Compute innovation between observations and forecast mean
4. **Kalman Gain Computation**: Compute traditional Kalman gain matrix
5. **Observation Perturbations**: Generate perturbed observations for each ensemble member
6. **Analysis Update**: Update each ensemble member using Kalman gain and perturbed observations

### Mathematical Formulation

The EnKF uses the traditional Kalman gain formulation:

```
K = P^f H^T (H P^f H^T + R)^(-1)
```

Where:
- `P^f` is the background error covariance matrix
- `H` is the observation operator
- `R` is the observation error covariance matrix

The analysis update for each ensemble member i is:

```
x_i^a = x_i^f + K (y^obs + ε_i - H x_i^f)
```

Where `ε_i` are observation perturbations drawn from N(0, R).

### Key Differences from ETKF

- **Traditional Kalman Gain**: Uses the standard Kalman gain formulation
- **Observation Perturbations**: Requires perturbed observations for each ensemble member
- **Ensemble Spread**: Maintains ensemble spread through the analysis
- **Computational Cost**: More computationally intensive but conceptually simpler

## Features

- **Multiple Inflation Methods**: Multiplicative, additive, and relaxation inflation
- **Comprehensive Diagnostics**: Innovation norm, spread statistics, Kalman gain analysis
- **Numerical Stability Monitoring**: Condition number tracking and warnings
- **Backend Agnostic**: Works with any backend that satisfies the required concepts
- **Robust Error Handling**: Comprehensive validation and error checking

## Usage

### Basic Usage

```cpp
#include "EnKF.hpp"

// Create EnKF instance
EnKF<BackendTag> enkf(ensemble, obs, obs_op, config);

// Perform analysis
enkf.Analyse();

// Save results
enkf.saveEnsemble();

// Get analysis results
auto results = enkf.getAnalysisResults();
```

### Configuration

The EnKF requires the following configuration parameters:

```yaml
# Analysis settings
analysis:
  algorithm: "enkf"
  inflation: 1.0  # Inflation factor for ensemble spread
  inflation_method: "multiplicative"  # Options: multiplicative, additive, relaxation
  format: "txt"
  output_base_file: "./analysis"
```

### Analysis Results

The `getAnalysisResults()` method returns a structure containing:

- `innovation_norm`: Norm of the innovation vector
- `analysis_increment_norm`: Norm of the analysis increment
- `background_spread`: Background ensemble spread
- `analysis_spread`: Analysis ensemble spread
- `max_kalman_gain`, `min_kalman_gain`: Kalman gain statistics
- `condition_number`: Condition number of innovation covariance
- `ensemble_size`: Number of ensemble members
- `observation_count`: Number of observations
- `inflation_method`: Method used for inflation
- `inflation_factor`: Inflation factor applied

## Applications

EnKF is particularly useful for:

- **Traditional Kalman Filter Applications**: When the standard Kalman gain formulation is preferred
- **Educational Purposes**: Easier to understand than ensemble space methods
- **Legacy System Integration**: When existing systems use traditional Kalman filter formulations
- **Research Applications**: When observation perturbations are needed for specific studies

## Limitations

- **Computational Cost**: More expensive than ETKF due to matrix operations
- **Observation Perturbations**: Requires generation of observation perturbations
- **Numerical Issues**: Can suffer from numerical instability with high-dimensional systems
- **Ensemble Size**: Performance degrades with small ensemble sizes

## Inflation Methods

### Multiplicative Inflation (Default)
Scales ensemble perturbations by a factor:
```
X'^f = √α * X'^f
```

### Additive Inflation
Adds random perturbations to ensemble members:
```
x_i^f = x_i^f + ε_i
```

### Relaxation Inflation
Reduces the analysis increment to maintain spread:
```
x_i^a = x_i^f + α * K * (y^obs - H x_i^f)
```

## References

1. Evensen, G. (1994) "Sequential data assimilation with a nonlinear quasi-geostrophic model using Monte Carlo methods to forecast error statistics", Journal of Geophysical Research, 99(C5), 10143-10162.

2. Burgers, G., van Leeuwen, P. J., & Evensen, G. (1998) "Analysis scheme in the ensemble Kalman filter", Monthly Weather Review, 126(6), 1719-1724.

3. Houtekamer, P. L., & Mitchell, H. L. (1998) "Data assimilation using an ensemble Kalman filter technique", Monthly Weather Review, 126(3), 796-811.

## Examples

See the tutorial configuration file `src/backends/simple/tutorial/enkf.yaml` for a complete example setup.

## Testing

Run the EnKF tests with:

```bash
make enkf_test
./enkf_test
```

## Driver Application

The EnKF driver application is located at:
`applications/data_assimilation/ensemble/enkf.cpp`

Usage:
```bash
./enkf <config_file>
```

## Comparison with Other Algorithms

| Feature | EnKF | ETKF | LETKF | Particle Filter |
|---------|------|------|-------|-----------------|
| Kalman Gain | Traditional | Ensemble space | Local ensemble space | Weight-based |
| Observation Perturbations | Required | Not needed | Not needed | Not needed |
| Computational Cost | High | Medium | Medium | High |
| Ensemble Spread | Maintained | Maintained | Maintained | Variable |
| Numerical Stability | Moderate | Good | Good | Variable |
| Localization | Not built-in | Not built-in | Built-in | Optional |

## Debug Configuration

A debug configuration for the EnKF can be added to VS Code's launch.json:

```json
{
    "name": "Debug Simple EnKF",
    "type": "cppdbg",
    "request": "launch",
    "program": "${workspaceFolder}/build/bin/enkf",
    "args": ["${workspaceFolder}/src/backends/simple/tutorial/enkf.yaml"],
    "stopAtEntry": true,
    "cwd": "${workspaceFolder}/build",
    "environment": [
        {
            "name": "GLOG_logtostderr",
            "value": "1"
        },
        {
            "name": "GLOG_v",
            "value": "0"
        },
        {
            "name": "GLOG_colorlogtostderr",
            "value": "1"
        },
        {
            "name": "GLOG_minloglevel",
            "value": "0"
        }
    ],
    "externalConsole": false,
    "MIMode": "gdb",
    "miDebuggerPath": "C:/msys64/mingw64/bin/gdb.exe",
    "setupCommands": [
        {
            "description": "Enable pretty-printing for gdb",
            "text": "-enable-pretty-printing",
            "ignoreFailures": true
        },
        {
            "description": "Set breakpoint pending on",
            "text": "-gdb-set breakpoint pending on",
            "ignoreFailures": true
        }
    ]
}
``` 