# Particle Filter Algorithm

## Overview

The Particle Filter (PF) is a non-parametric Bayesian data assimilation method that represents the posterior probability density function using a set of weighted particles (ensemble members). Unlike Kalman filter variants that assume Gaussian distributions, particle filters can handle non-linear and non-Gaussian systems.

## Algorithm Description

The particle filter algorithm follows these steps:

1. **Forecast Step**: Propagate particles using the model (if available)
2. **Analysis Step**: Update particle weights using observations
3. **Resampling Step**: Resample particles based on weights to prevent degeneracy
4. **Optional Jittering**: Apply small perturbations to maintain diversity

### Weight Update

The weight update follows Bayes' rule:

```
w_i^{(a)} ∝ w_i^{(f)} p(y|X_i^{(f)})
```

Where:
- `w_i^{(a)}` is the analysis weight of particle i
- `w_i^{(f)}` is the forecast weight of particle i  
- `p(y|X_i^{(f)})` is the likelihood of observation y given particle i

For Gaussian observation errors, the likelihood is:

```
p(y|X_i^{(f)}) ∝ exp(-1/2 * (y - H(X_i^{(f)}))^T R^(-1) (y - H(X_i^{(f)})))
```

### Resampling Methods

The implementation supports four resampling methods:

1. **Systematic Resampling** (default): Most efficient and commonly used
2. **Multinomial Resampling**: Simple but less efficient
3. **Stratified Resampling**: Reduces variance compared to multinomial
4. **Residual Resampling**: Combines deterministic and random sampling

### Effective Sample Size

The algorithm monitors the effective sample size (ESS) to determine when resampling is needed:

```
ESS = 1 / Σ(w_i^2)
```

Resampling is triggered when ESS falls below a configurable threshold.

## Features

- **Multiple Resampling Methods**: Systematic, multinomial, stratified, and residual
- **Adaptive Resampling**: Based on effective sample size threshold
- **Jittering**: Optional perturbation to maintain particle diversity
- **Weight Management**: External weight setting and normalization
- **Comprehensive Diagnostics**: Weight statistics, likelihood values, resampling info
- **Backend Agnostic**: Works with any backend that satisfies the required concepts

## Usage

### Basic Usage

```cpp
#include "ParticleFilter.hpp"

// Create particle filter instance
ParticleFilter<BackendTag> pf(ensemble, obs, obs_op, config);

// Perform analysis
pf.Analyse();

// Save results
pf.saveEnsemble();

// Get analysis results
auto results = pf.getAnalysisResults();
```

### Configuration

The particle filter requires the following configuration parameters:

```yaml
# Particle Filter specific parameters
resampling_threshold: 50  # Effective sample size threshold
resampling_method: "systematic"  # Resampling method
jittering_enabled: true  # Enable jittering
jittering_std: 0.01  # Jittering standard deviation

# Output configuration
output_base_file: "analysis_particle_filter"
format: "txt"
save_detailed_results: true
```

### Analysis Results

The `getAnalysisResults()` method returns a structure containing:

- `weights`: Final particle weights
- `effective_sample_size`: Effective sample size
- `max_weight`, `min_weight`: Weight statistics
- `weight_variance`: Variance of weights
- `resampling_performed`: Whether resampling occurred
- `resampling_method`: Method used for resampling
- `resampling_threshold`: Threshold for resampling
- `likelihood_values`: Likelihood values for each particle

## Applications

Particle filters are particularly useful for:

- **Non-linear systems**: Where linearization assumptions fail
- **Non-Gaussian errors**: When observation or model errors are not Gaussian
- **Multi-modal distributions**: When the posterior has multiple peaks
- **High-dimensional systems**: With appropriate localization techniques

## Limitations

- **Computational cost**: Scales with ensemble size
- **Degeneracy**: Weights can become concentrated on few particles
- **Curse of dimensionality**: Performance degrades in high dimensions
- **Tuning**: Requires careful selection of resampling threshold and method

## References

1. Arulampalam, M. S., Maskell, S., Gordon, N., & Clapp, T. (2002). A tutorial on particle filters for online nonlinear/non-Gaussian Bayesian tracking. IEEE Transactions on signal processing, 50(2), 174-188.

2. Poterjoy, J. (2016). A localized particle filter for high-dimensional nonlinear systems. Monthly Weather Review, 144(1), 59-76.

3. van Leeuwen, P. J. (2009). Particle filtering in geophysical systems. Monthly Weather Review, 137(12), 4089-4114.

## Examples

See the tutorial configuration file `src/backends/simple/tutorial/particle_filter.yaml` for a complete example setup.

## Testing

Run the particle filter tests with:

```bash
make particle_filter_test
./particle_filter_test
```

## Driver Application

The particle filter driver application is located at:
`applications/data_assimilation/ensemble/particle_filter.cpp`

Usage:
```bash
./particle_filter <config_file>
``` 