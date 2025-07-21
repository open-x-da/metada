# Local Weighted Ensemble Kalman Filter (LWEnKF)

## Overview

The Local Weighted Ensemble Kalman Filter (LWEnKF) is an advanced ensemble-based data assimilation algorithm that extends the traditional Ensemble Kalman Filter (EnKF) by incorporating localization and adaptive weighting features. This implementation provides improved covariance estimation and better handling of spurious correlations in high-dimensional systems.

## Key Features

### Localization
- **Distance-based localization**: Reduces spurious correlations between distant state variables and observations
- **Multiple localization functions**: Gaussian, exponential, cutoff, and Gaspari-Cohn functions
- **Configurable localization radius**: Adjustable based on the problem characteristics

### Adaptive Weighting
- **Multiple weighting schemes**: Uniform, adaptive, inverse variance, and likelihood-based weighting
- **Ensemble member weighting**: Improves covariance estimation by weighting ensemble members based on their quality
- **Dynamic weight adjustment**: Weights are updated during the analysis process

### Inflation Methods
- **Multiplicative inflation**: Traditional inflation method
- **Additive inflation**: Alternative inflation approach
- **Relaxation inflation**: Relaxation-based inflation

## Algorithm Description

The LWEnKF algorithm follows these main steps:

1. **Forecast Step**: Propagate ensemble members using the model
2. **Weight Computation**: Compute adaptive weights for ensemble members
3. **Localization**: Apply distance-based localization to reduce spurious correlations
4. **Weighted Covariance**: Compute weighted ensemble covariance
5. **Kalman Gain**: Compute localized Kalman gain
6. **Analysis Update**: Update ensemble members using observations
7. **Ensemble Statistics**: Update ensemble mean and perturbations

## Mathematical Formulation

### Weighted Covariance
The weighted background covariance is computed as:
```
P^f_weighted = Σ(w_i * (x_i^f - x̄^f)(x_i^f - x̄^f)^T)
```

### Localized Kalman Gain
The localized Kalman gain is computed as:
```
K_local = P^f_local * H^T * (H * P^f_local * H^T + R)^(-1)
```

### Analysis Update
Each ensemble member is updated as:
```
x_i^a = x_i^f + K_local * (y^obs + ε_i - H * x_i^f)
```

## Configuration

### Required Parameters

```yaml
analysis:
  # Algorithm parameters
  inflation: 1.1                    # Inflation factor
  inflation_method: "multiplicative" # Inflation method
  
  # Localization parameters
  localization_radius: 0.5          # Localization radius
  localization_function: "gaussian"  # Localization function
  
  # Weighting parameters
  weighting_scheme: "adaptive"      # Weighting scheme
  
  # Output configuration
  output_base_file: "analysis"      # Output file base name
  format: "txt"                     # Output format
```

### Localization Functions

- **gaussian**: `exp(-0.5 * (d/r)^2)`
- **exponential**: `exp(-d/r)`
- **cutoff**: `1 if d ≤ r, 0 otherwise`
- **gaspari_cohn**: Gaspari-Cohn function for smooth localization

### Weighting Schemes

- **uniform**: Equal weights for all ensemble members
- **adaptive**: Weights inversely proportional to ensemble spread
- **inverse_var**: Weights inversely proportional to variance
- **likelihood**: Weights based on observation likelihood

## Usage

### Command Line
```bash
./lwenkf config.yaml
```

### Example Configuration
```yaml
backend:
  type: "simple"
  ensemble:
    size: 8
    member_files:
      - "ens_1.txt"
      - "ens_2.txt"
      - "ens_3.txt"
      - "ens_4.txt"
      - "ens_5.txt"
      - "ens_6.txt"
      - "ens_7.txt"
      - "ens_8.txt"
  
  observation:
    file: "obsC.txt"
    covariance_file: "obsB_format.txt"
  
  observation_operator:
    type: "identity"

analysis:
  inflation: 1.1
  inflation_method: "multiplicative"
  localization_radius: 0.5
  localization_function: "gaussian"
  weighting_scheme: "adaptive"
  output_base_file: "analysis"
  format: "txt"

logging:
  level: "info"
  output: "console"
```

## Output Files

The LWEnKF generates the following output files:

- `analysis_mean.txt`: Analysis mean state
- `analysis_member_0.txt` to `analysis_member_N.txt`: Individual ensemble members
- `analysis_diagnostics.txt`: Diagnostic information

## Diagnostics

The algorithm provides comprehensive diagnostics including:

- **Innovation statistics**: Innovation norm and analysis increment norm
- **Spread statistics**: Background and analysis ensemble spreads
- **Kalman gain statistics**: Maximum and minimum Kalman gain elements
- **Numerical stability**: Condition number of innovation covariance
- **Weight statistics**: Maximum, minimum, and variance of ensemble weights
- **Localization information**: Localization radius and function used

## Advantages

1. **Improved covariance estimation**: Weighted ensemble statistics provide better covariance estimates
2. **Reduced spurious correlations**: Localization eliminates unrealistic correlations
3. **Adaptive ensemble management**: Dynamic weighting improves ensemble quality
4. **Numerical stability**: Better handling of ill-conditioned problems
5. **Flexible configuration**: Multiple options for localization and weighting

## Limitations

1. **Computational cost**: Additional computations for localization and weighting
2. **Parameter tuning**: Requires careful tuning of localization radius and weighting parameters
3. **Distance computation**: Requires distance information between state variables and observations
4. **Memory usage**: Localization matrices can be memory-intensive for large problems

## References

- Chen, Y. et al. "Local Weighted Ensemble Kalman Filter"
- Hamill, T. M., Whitaker, J. S., & Snyder, C. (2001) "Distance-dependent filtering of background error covariance estimates in an ensemble Kalman filter"
- Gaspari, G., & Cohn, S. E. (1999) "Construction of correlation functions in two and three dimensions"

## Testing

Run the unit tests to verify the implementation:

```bash
make lwenkf_test
./lwenkf_test
```

## Performance Considerations

- **Localization radius**: Smaller radius reduces computational cost but may lose important correlations
- **Weighting scheme**: Adaptive weighting provides better results but increases computational cost
- **Ensemble size**: Larger ensembles provide better statistics but increase memory and computational requirements
- **Observation density**: Higher observation density benefits more from localization

## Future Improvements

1. **Advanced localization**: Implement more sophisticated localization methods
2. **Adaptive localization**: Dynamic localization radius based on ensemble statistics
3. **Hybrid methods**: Combine with variational methods for hybrid data assimilation
4. **Parallel implementation**: Optimize for parallel computing environments
5. **Observation localization**: Apply localization to observation space as well 