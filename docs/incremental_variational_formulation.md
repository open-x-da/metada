# Incremental Variational Formulation in METADA

## Overview

This document describes the implementation of the incremental variational formulation in METADA, which provides better numerical conditioning and follows the standard approach used in operational systems like WRFDA.

## Mathematical Formulation

### Total State Formulation (Original)
The original METADA implementation uses the total state formulation:

```
J(x) = 1/2 * (x-xb)^T B^-1 (x-xb) + 1/2 * Σᵢ (yᵢ - Hᵢ(Mᵢ(x)))^T Rᵢ^-1 (yᵢ - Hᵢ(Mᵢ(x)))
```

Where:
- `x` is the analysis state (to be determined)
- `xb` is the background state (first guess)
- `B` is the background error covariance matrix
- `yᵢ` are observations at time i
- `Hᵢ` is the observation operator at time i
- `Mᵢ` is the model propagated to time i
- `Rᵢ` is the observation error covariance matrix at time i

### Incremental Formulation (New)
The incremental formulation works with analysis increments:

```
J(δx) = 1/2 * δx^T B^-1 δx + 1/2 * Σᵢ (dᵢ - Hᵢ(Mᵢ(xb + δx)))^T Rᵢ^-1 (dᵢ - Hᵢ(Mᵢ(xb + δx)))
```

Where:
- `δx` is the analysis increment (δx = x - xb)
- `dᵢ = yᵢ - Hᵢ(Mᵢ(xb))` is the innovation vector (O-B)
- All other variables have the same meaning as above

## Key Advantages

1. **Better Numerical Conditioning**: Working with increments rather than absolute values
2. **Easier Background Error Covariance Handling**: B matrix operates directly on increments
3. **More Efficient Linear Algebra**: Smaller magnitude values in computations
4. **Standard Approach**: Matches WRFDA and other operational systems
5. **Pre-computed Innovations**: Innovation vectors are computed once at initialization

## Implementation Components

### 1. IncrementalCostFunction
- **File**: `src/framework/adapters/IncrementalCostFunction.hpp`
- **Purpose**: Implements the incremental cost function J(δx)
- **Key Features**:
  - Pre-computes innovation vectors dᵢ = yᵢ - Hᵢ(Mᵢ(xb))
  - Supports 3DVAR, FGAT, and 4DVAR formulations
  - Computes gradients using adjoint methods

### 2. IncrementalMinimization
- **File**: `src/framework/adapters/IncrementalMinimization.hpp`
- **Purpose**: Minimization adapter for incremental formulation
- **Key Features**:
  - Works with increments δx instead of total states x
  - Bridges framework types to optimization algorithms
  - Supports L-BFGS, CG, and Steepest Descent

### 3. IncrementalVariational
- **File**: `src/framework/algorithms/IncrementalVariational.hpp`
- **Purpose**: Main incremental variational algorithm
- **Key Features**:
  - Orchestrates the incremental variational process
  - Manages cost function and minimization
  - Provides analysis results with both state and increment
  - **NEW**: Comprehensive verification methods for observation operators
  - **NEW**: Individual and combined check methods for TL/AD consistency
  - **NEW**: Integration with observation operator verification framework

### 4. IncrementalGradientChecks
- **File**: `src/framework/algorithms/IncrementalGradientChecks.hpp`
- **Purpose**: Gradient and observation operator verification for incremental formulation
- **Key Features**:
  - Finite difference gradient tests for cost function
  - Multiple random direction testing
  - Comprehensive adjoint verification
  - **NEW**: Observation operator tangent linear/adjoint consistency checks
  - **NEW**: Observation operator tangent linear finite difference checks
  - **NEW**: Comprehensive observation operator verification

### 5. Verification Methods

The incremental variational implementation includes comprehensive verification methods to ensure correctness:

#### 5.1 Cost Function Gradient Verification
- **Method**: `performGradientTest()`
- **Purpose**: Verifies the analytical gradient matches finite difference approximation
- **Mathematical Check**: ∇J(δx) ≈ [J(δx + ε·d) - J(δx - ε·d)] / (2ε)

#### 5.2 Observation Operator Tangent Linear/Adjoint Consistency
- **Method**: `performObsOperatorTLADCheck()`
- **Purpose**: Verifies TL/AD consistency for observation operators
- **Mathematical Check**: ⟨H'dx, dy⟩ = ⟨dx, H'^T dy⟩

#### 5.3 Observation Operator Tangent Linear Verification
- **Method**: `performObsOperatorTangentLinearCheck()`
- **Purpose**: Verifies tangent linear implementation using finite differences
- **Mathematical Check**: H'dx ≈ [H(x_b + ε·dx) - H(x_b)] / ε

#### 5.4 Comprehensive Verification
- **Method**: `performVerificationChecks()`
- **Purpose**: Performs all verification checks in sequence
- **Includes**: Cost function gradient + observation operator checks

### 6. Incremental Variational Application
- **File**: `applications/data_assimilation/variational/incremental_variational.cpp`
- **Purpose**: Driver application for incremental variational
- **Key Features**:
  - Command-line interface
  - Configuration management
  - Results output and diagnostics

## Usage

### Building
The incremental variational components are automatically built with the main METADA build:

```bash
cd build
make incremental_variational
```

### Running
Use the incremental variational application with a configuration file:

```bash
cd build
./bin/incremental_variational ../src/backends/wrf/tutorial/incremental_var3d.yaml
```

### Configuration
The incremental variational application uses the same configuration format as the original variational application, with the same YAML structure.

## Key Differences from Original Implementation

1. **Cost Function**: Works with increments δx instead of total states x
2. **Innovation Pre-computation**: Innovation vectors are computed once at initialization
3. **Minimization**: Optimizes over increments starting from zero
4. **Results**: Provides both analysis state (xb + δx*) and analysis increment (δx*)
5. **Gradient Testing**: Specialized gradient tests for incremental formulation

## Backward Compatibility

The original total state formulation remains available and unchanged. The incremental formulation is implemented as a parallel system, allowing users to choose between:

- **Original**: `variational` application with total state formulation
- **Incremental**: `incremental_variational` application with incremental formulation

## Future Enhancements

1. **Hybrid Formulation**: Combine total and incremental approaches
2. **Incremental File Format**: Dedicated file format for increments
3. **Advanced Preconditioning**: Leverage incremental formulation for better preconditioning
4. **Ensemble Integration**: Integrate with ensemble methods

## References

1. Courtier, P., et al. (1994). "A strategy for operational implementation of 4D-Var, using an incremental approach." Quarterly Journal of the Royal Meteorological Society, 120(519), 1367-1387.

2. Rabier, F., et al. (2000). "The ECMWF operational implementation of four-dimensional variational assimilation. I: Experimental results with simplified physics." Quarterly Journal of the Royal Meteorological Society, 126(564), 1143-1170.

3. WRFDA Documentation: https://www2.mmm.ucar.edu/wrf/users/wrfda/
