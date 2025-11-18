# WRFDA `da_minimise_cg` Workflow Analysis

## Settings
- `precondition_cg = false` (no preconditioning, `precon = 1.0`)
- `orthonorm_gradient = true` (gradient reorthonormalization enabled)
- `WRF_CHEM = 0` (no chemistry)
- `calculate_cg_cost_fn = false` (use cost approximation instead of full evaluation)
- `write_detail_grad_fn = false` (no detailed gradient output)

## Workflow Overview

This is the **inner loop** of the incremental variational method. It minimizes the cost function in control space using Conjugate Gradient (CG) with gradient reorthonormalization.

### Key Variables
- `cv`: Outer-loop control variable (updated between outer iterations)
- `xhat`: Inner-loop control variable (updated during CG iterations)
- `ghat`: Gradient vector `∇J(xhat+cv)`
- `phat`: Search direction (conjugate direction)
- `fhat`: Gradient along search direction `∇J(phat)`
- `qhat`: Orthonormalized gradient vectors (for reorthonormalization)

---

## Detailed Workflow

### [1.0] Initialization (lines 67-81)

1. **Compute initial cost** (line 70-71):
   ```fortran
   call da_calculate_j(it, 0, cv_size, ..., xhat, cv, re, y, j_cost, ...)
   ```
   - Computes `J(xhat, cv)` and stores in `j_cost`
   - Also computes residual `re = (O-B) - H(xhat)` in `y` structure

2. **Compute initial gradient** (line 73-74):
   ```fortran
   call da_calculate_gradj(..., xhat+cv, y, ghat, ..., re)
   ```
   - **WITH `re` argument**: Uses residual computed by `da_calculate_j`
   - Computes `ghat = ∇J(xhat+cv)` 
   - Path: `da_calculate_grady(iv, re, jo_grad_y)` → `da_transform_vtoy_adj`
   - This is the **full gradient** including all cost terms (Jb, Jo, Je, etc.)

3. **Store initial state** (lines 76-81):
   - `j0_total = j_cost%total` (initial cost)
   - `ghat0 = ghat` (initial gradient, used for cost approximation)

### [1.1] Preconditioning (lines 83-129)

- **Since `precondition_cg = false`**: This section is **skipped**
- `precon = 1.0` (identity preconditioner)

### [1.2] Initial Search Direction (line 131)

```fortran
phat = - precon * ghat  ! = -ghat (since precon = 1.0)
```

### [1.3] Initial Gradient Norm (lines 133-138)

```fortran
rrmold = da_dot_cv(cv_size, -phat, ghat, ...)  ! = da_dot_cv(ghat, ghat) = ||ghat||²
j_grad_norm_target = sqrt(rrmold)  ! = ||ghat||
```

### [1.4] Orthonormalization Setup (lines 140-144)

- **Since `orthonorm_gradient = true`**:
  - Allocates `qhat(0:ntmax(it))` array
  - Normalizes initial gradient: `qhat(0)%values = ghat / sqrt(rrmold)`

---

## [2.0] CG Iteration Loop (lines 161-235)

For each iteration `iter = 1, ..., ntmax(it)`:

### Step 1: Compute Gradient Along Search Direction (line 164-165)

```fortran
call da_calculate_gradj(..., phat, y, fhat, ...)
```

- **WITHOUT `re` argument**: Must compute residual internally
- Computes `fhat = ∇J(phat)` (gradient of cost w.r.t. search direction)
- **Path in `da_calculate_gradj` (lines 108-134)**:
  1. `da_transform_vtoy(cv_size, ..., phat, ..., y, ...)` 
     - Transforms `phat` to observation space: `y = H'(phat)`
  2. `da_calculate_grady(iv, y, jo_grad_y)`
     - Computes weighted residual: `jo_grad_y = -R⁻¹ · y`
  3. `da_transform_vtoy_adj(cv_size, ..., grad_jo, ..., jo_grad_y, ...)`
     - Adjoint transform: `grad_jo = H'^T · jo_grad_y`
  4. `grad_jo = -grad_jo` (sign compensation)
- This computes the **observation gradient** component; background gradient is `jb_factor * phat`

### Step 2: Compute Step Size (lines 167-174)

```fortran
apdotp = da_dot_cv(cv_size, fhat, phat, ...)  ! <∇J(phat), phat>
step = rrmold / apdotp  ! Optimal step size (if apdotp > 0)
```

- `apdotp` is the curvature along search direction
- Step size minimizes cost along `phat`

### Step 3: Update Control Variable and Gradient (lines 176-177)

```fortran
ghat = ghat + step * fhat  ! Update gradient
xhat = xhat + step * phat  ! Update control variable
```

- CG update: move along search direction `phat` by step `step`

### Step 4: Orthonormalize Gradient (lines 179-189)

- **Since `orthonorm_gradient = true`**:
  - Modified Gram-Schmidt reorthonormalization:
    ```fortran
    do i = iter-1, 0, -1
       gdot = da_dot_cv(cv_size, precon*ghat, qhat(i)%values, ...)
       ghat = ghat - gdot * qhat(i)%values  ! Remove component along qhat(i)
    end do
    ```
  - Ensures gradient is orthogonal to all previous search directions
  - Prevents loss of conjugacy due to numerical errors

### Step 5: Compute New Gradient Norm (lines 191-196)

```fortran
rrmnew = da_dot_cv(cv_size, precon*ghat, ghat, ...)  ! = ||ghat||² (since precon=1.0)
rrmnew_norm = sqrt(rrmnew)  ! = ||ghat||
```

### Step 6: Update Search Direction (lines 198-206)

```fortran
ratio = rrmnew / rrmold  ! Polak-Ribière β
phat = - precon * ghat + ratio * phat  ! = -ghat + β*phat
```

- Standard CG search direction update
- If `orthonorm_gradient = true`: also stores normalized gradient in `qhat(iter)`

### Step 7: Compute Cost (lines 213-223)

- **Since `calculate_cg_cost_fn = false`**:
  ```fortran
  j_total = j0_total + 0.5 * da_dot_cv(cv_size, ghat0, xhat, ...)
  ```
  - **Cost approximation** using initial gradient and current control variable
  - Avoids expensive full cost evaluation (`da_calculate_j`)
  - Based on: `J(xhat) ≈ J(0) + <∇J(0), xhat> + 0.5*<∇J(0), xhat>`

- **If `calculate_cg_cost_fn = true`**: Would call `da_calculate_j` for exact cost

### Step 8: Check Convergence (line 231)

```fortran
if (rrmnew_norm < eps(it) * j_grad_norm_target) exit
```

- Convergence when gradient norm is sufficiently small relative to initial gradient norm

---

## [3.0] Finalization (lines 237-275)

1. **Compute final cost** (line 258-259):
   ```fortran
   call da_calculate_j(it, iter, ..., xhat, cv, re, y, j_cost, ...)
   ```
   - Full cost evaluation at final `xhat`

2. **Compute final gradient** (line 261-262):
   ```fortran
   call da_calculate_gradj(..., xhat+cv, y, ghat, ..., re)
   ```
   - **WITH `re` argument**: Uses residual from final cost evaluation
   - Final gradient for diagnostics

3. **Compute final gradient norm** (line 264-268):
   ```fortran
   rrmnew_norm = sqrt(da_dot_cv(cv_size, ghat, ghat, ...))
   ```

---

## Key Observations: `da_calculate_gradj` Usage

### With `re` Argument (lines 73-74, 261-262)

**When**: Initial and final gradient computation
- `re` is the residual `(O-B) - H(xhat)` computed by `da_calculate_j`
- **Path**: `da_calculate_grady(iv, re, jo_grad_y)` → `da_transform_vtoy_adj`
- **Efficient**: Reuses residual from cost evaluation

### Without `re` Argument (line 164-165)

**When**: During CG iterations (computing gradient along search direction)
- Must compute residual internally
- **Path**: `da_transform_vtoy(cv, ..., y)` → `da_calculate_grady(iv, y, jo_grad_y)` → `da_transform_vtoy_adj`
- **More expensive**: Requires forward transform `v → y` before computing gradient

### Why Two Paths?

1. **With `re`**: When residual is already available (from cost evaluation), reuse it
2. **Without `re`**: When computing gradient of arbitrary control vector (like search direction `phat`), must compute residual from scratch

---

## Cost Function Structure

The gradient `grad = ∇J(v)` consists of:
- `grad_jb`: Background term gradient `= jb_factor * v`
- `grad_jo`: Observation term gradient `= H'^T · (-R⁻¹ · residual)`
- `grad_je`: Ensemble term gradient `= je_factor * v` (if applicable)
- `grad_jp`: VarBC term gradient (if applicable)
- `grad_js`: Satellite CV term gradient (if applicable)
- `grad_jl`: Lateral boundary term gradient (if applicable)
- `grad_jm`: Divergence constraint gradient (if applicable)
- `grad_jt`: TAMDAR VarBC term gradient (if applicable)

**Total**: `grad = grad_jo + grad_jb + grad_je + grad_jd + grad_jp + grad_js + grad_jl + grad_jm + grad_jt`

---

## Summary

The CG minimization workflow:
1. **Initialize**: Compute cost and gradient at `xhat=0`
2. **Iterate**: 
   - Compute gradient along search direction (without `re`)
   - Update control variable and gradient
   - Reorthonormalize gradient (if enabled)
   - Update search direction using CG formula
   - Approximate cost (if `calculate_cg_cost_fn=false`)
3. **Finalize**: Compute exact final cost and gradient (with `re`)

The key difference between `da_calculate_gradj` with/without `re`:
- **With `re`**: Reuses precomputed residual (efficient)
- **Without `re`**: Computes residual internally via forward transform (necessary for arbitrary control vectors)

