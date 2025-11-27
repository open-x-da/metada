# MetaDA WRFDA CG Alignment

## Overview

This document describes the refactoring of MetaDA's inner-loop minimization workflow to align with WRFDA's `da_minimise_cg` algorithm.

## Changes Made

### 1. New WRFDA-Aligned CG Optimizer

**File**: `src/framework/base/WRFDAConjugateGradientOptimizer.hpp`

Created a new optimizer class that implements the exact workflow from WRFDA's `da_minimise_cg.inc`:

#### Key Features:
- **Gradient along search direction**: Computes `∇J(phat)` where `phat` is the search direction (matches WRFDA line 164-165)
- **WRFDA step size formula**: `step = rrmold / apdotp` where `apdotp = <∇J(phat), phat>` (matches WRFDA line 167-174)
- **Cost approximation**: Uses `J ≈ J0 + 0.5 * <ghat0, xhat>` during iterations (matches WRFDA line 218-222 when `calculate_cg_cost_fn=false`)
- **Gradient reorthonormalization**: Modified Gram-Schmidt algorithm (matches WRFDA lines 180-189 when `orthonorm_gradient=true`)

#### Workflow Alignment:

**Initialization** (matches WRFDA lines 70-81):
1. Compute initial cost `J(xhat=0)`
2. Compute initial gradient `ghat = ∇J(xhat=0)`
3. Set initial search direction `phat = -ghat`
4. Setup orthonormalization (if enabled)

**CG Iteration Loop** (matches WRFDA lines 161-235):
1. Compute gradient along search direction: `fhat = ∇J(phat)`
2. Compute step size: `step = rrmold / apdotp`
3. Update: `ghat = ghat + step * fhat`, `xhat = xhat + step * phat`
4. Reorthonormalize gradient (if enabled)
5. Update search direction: `phat = -ghat + β*phat`
6. Approximate cost (if enabled)
7. Check convergence

**Finalization** (matches WRFDA lines 258-268):
1. Compute exact final cost
2. Compute exact final gradient

### 2. Updated IncrementalMinimization

**File**: `src/framework/adapters/IncrementalMinimization.hpp`

Added support for the new WRFDA-aligned optimizer:

- Added include for `WRFDAConjugateGradientOptimizer.hpp`
- Updated `createOptimizer()` to recognize `"WRFDA-CG"` or `"WRFDA_CG"` algorithm
- Reads configuration options:
  - `use_cost_approximation` (default: `true`)
  - `orthonorm_gradient` (default: `true`)

## Configuration

To use the WRFDA-aligned CG optimizer, set in your configuration:

```yaml
minimization_algorithm: "WRFDA-CG"  # or "WRFDA_CG"
max_iterations: 100
cost_tolerance: 1e-8
gradient_tolerance: 1e-6
use_cost_approximation: true   # Use cost approximation during iterations
orthonorm_gradient: true       # Enable gradient reorthonormalization
```

## Differences from Standard CG

### Standard CG (`ConjugateGradientOptimizer`):
- Computes gradient at new solution after line search
- Uses line search to find step size
- Computes exact cost at each iteration
- No gradient reorthonormalization

### WRFDA-Aligned CG (`WRFDAConjugateGradientOptimizer`):
- Computes gradient along search direction (before update)
- Uses WRFDA's step size formula: `step = rrmold / apdotp`
- Uses cost approximation during iterations (if enabled)
- Supports gradient reorthonormalization (if enabled)

## Workflow Comparison

### WRFDA `da_minimise_cg`:
```
[1.0] Initialization:
  - da_calculate_j(xhat=0) → J0, re
  - da_calculate_gradj(xhat=0, re) → ghat
  - phat = -ghat

[2.0] CG Loop:
  - da_calculate_gradj(phat) → fhat = ∇J(phat)
  - step = rrmold / <fhat, phat>
  - ghat = ghat + step * fhat
  - xhat = xhat + step * phat
  - Reorthonormalize ghat
  - phat = -ghat + β*phat
  - J ≈ J0 + 0.5 * <ghat0, xhat>

[3.0] Finalization:
  - da_calculate_j(xhat) → J_final
  - da_calculate_gradj(xhat, re) → ghat_final
```

### MetaDA WRFDA-Aligned CG:
```
[1.0] Initialization:
  - cost_func(xhat=0) → J0
  - gradient_func(xhat=0) → ghat
  - phat = -ghat

[2.0] CG Loop:
  - gradient_func(phat) → fhat = ∇J(phat)
  - step = rrmold / <fhat, phat>
  - ghat = ghat + step * fhat
  - xhat = xhat + step * phat
  - Reorthonormalize ghat
  - phat = -ghat + β*phat
  - J ≈ J0 + 0.5 * <ghat0, xhat>

[3.0] Finalization:
  - cost_func(xhat) → J_final
  - gradient_func(xhat) → ghat_final
```

## Notes

1. **Gradient along search direction**: The key difference is computing `∇J(phat)` where `phat` is the search direction, not the current solution. This matches WRFDA's workflow exactly.

2. **Cost approximation**: During iterations, the cost is approximated using the initial gradient and current control variable. This avoids expensive cost evaluations but may be less accurate.

3. **Gradient reorthonormalization**: The modified Gram-Schmidt algorithm ensures the gradient remains orthogonal to previous search directions, preventing loss of conjugacy due to numerical errors.

4. **Step size**: WRFDA uses a direct formula based on the curvature along the search direction, rather than line search.

## Testing

To test the alignment:
1. Run MetaDA with `minimization_algorithm: "WRFDA-CG"`
2. Compare trace output with WRFDA's `da_minimise_cg` trace
3. Verify that:
   - Gradient is computed along search direction
   - Step size matches WRFDA's formula
   - Cost approximation is used (if enabled)
   - Gradient reorthonormalization is performed (if enabled)

## Future Work

- [ ] Add support for preconditioning (when `precondition_cg=true`)
- [ ] Add support for exact cost evaluation during iterations (when `calculate_cg_cost_fn=true`)
- [ ] Add detailed gradient output (when `write_detail_grad_fn=true`)
- [ ] Test against WRFDA trace files to verify exact alignment

