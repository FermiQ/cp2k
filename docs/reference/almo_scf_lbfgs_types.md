# almo_scf_lbfgs_types Module Documentation

## Overview

The `almo_scf_lbfgs_types` module in CP2K implements the Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm. L-BFGS is a quasi-Newton optimization method that is well-suited for large-scale unconstrained optimization problems, such as those encountered in electronic structure calculations (e.g., minimizing energy with respect to molecular orbital coefficients or density matrix elements). It approximates the inverse Hessian matrix using information from a limited number of previous iterations, making it memory-efficient.

This module provides the necessary data structures to store the L-BFGS history (differences in variables and gradients) and routines to create, manage, and utilize this history to compute a search direction. The implementation appears to work with variables and gradients represented as DBCSR (Distributed Block Compressed Sparse Row) matrices.

## Key Components

### Derived Type: `lbfgs_history_type`

This is the central data structure for storing the information required by the L-BFGS algorithm.

*   **`nstore` (INTEGER):** The maximum number of past iteration updates (pairs of variable differences `s_k` and gradient differences `y_k`) to store. This defines the "memory" of the L-BFGS algorithm.
*   **`istore(2)` (INTEGER, DIMENSION(2)):** An array of two integers. `istore(1)` counts the total number of variable difference vectors (`s_k`) pushed with `action=2` (see `lbfgs_history_push`). `istore(2)` counts the same for gradient difference vectors (`y_k`). These counters help manage the circular buffer for storing history.
*   **`matrix` (TYPE(`dbcsr_type`), DIMENSION(:,:,:), ALLOCATABLE):** A 3D allocatable array of DBCSR matrices.
    *   Dimension 1: Spin index (e.g., for spin-polarized calculations).
    *   Dimension 2: History slot index (from 1 to `nstore`). This is a circular buffer.
    *   Dimension 3: Type of data stored:
        *   `matrix(:,:,1)`: Stores variable differences `s_k = x_{k+1} - x_k`.
        *   `matrix(:,:,2)`: Stores gradient differences `y_k = g_{k+1} - g_k`.
    Initially, before differences are computed, these slots may temporarily hold the actual variables `x_k` and gradients `g_k`.
*   **`rho` (REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE):** A 2D allocatable array storing the scalar values `rho_k = 1.0 / (y_k^T * s_k)` for each history entry `k` and spin index.
    *   Dimension 1: Spin index.
    *   Dimension 2: History slot index.

### Public Subroutines

1.  **`lbfgs_create(history, nspins, nstore)`**:
    *   **Purpose:** Initializes the `lbfgs_history_type` object.
    *   **Arguments:**
        *   `history` (TYPE(`lbfgs_history_type`), INTENT(INOUT)): The L-BFGS history object to be initialized.
        *   `nspins` (INTEGER, INTENT(IN)): The number of spin channels.
        *   `nstore` (INTEGER, INTENT(IN)): The maximum number of history entries to store.
    *   **Functionality:** Sets `history%nstore` and allocates the `history%matrix` and `history%rho` arrays according to `nspins` and `nstore`. Initializes `history%istore` to zero.

2.  **`lbfgs_seed(history, variable, gradient)`**:
    *   **Purpose:** Seeds the L-BFGS algorithm with the initial variable vector (`x_0`) and its corresponding gradient vector (`g_0`). This is the starting point for the optimization.
    *   **Arguments:**
        *   `history` (TYPE(`lbfgs_history_type`), INTENT(INOUT)): The L-BFGS history object.
        *   `variable` (TYPE(`dbcsr_type`), DIMENSION(:), INTENT(IN)): The initial variable vector(s) (one per spin).
        *   `gradient` (TYPE(`dbcsr_type`), DIMENSION(:), INTENT(IN)): The initial gradient vector(s).
    *   **Functionality:** Calls the internal `lbfgs_history_push` subroutine with `action=1` (direct storage) for both the `variable` (as `vartype=1`) and `gradient` (as `vartype=2`) to store them in the first history slot.

3.  **`lbfgs_get_direction(history, variable, gradient, direction)`**:
    *   **Purpose:** This is the main routine called at each L-BFGS iteration. It updates the history with the latest step and computes the new search direction.
    *   **Arguments:**
        *   `history` (TYPE(`lbfgs_history_type`), INTENT(INOUT)): The L-BFGS history object.
        *   `variable` (TYPE(`dbcsr_type`), DIMENSION(:), INTENT(IN)): The current variable vector (`x_{k+1}`).
        *   `gradient` (TYPE(`dbcsr_type`), DIMENSION(:), INTENT(IN)): The current gradient vector (`g_{k+1}`).
        *   `direction` (TYPE(`dbcsr_type`), DIMENSION(:), INTENT(INOUT)): The computed L-BFGS search direction (`d_{k+1}`).
    *   **Functionality:**
        1.  Calculates `s_k = x_{k+1} - x_k` and `y_k = g_{k+1} - g_k` by calling `lbfgs_history_push` with `action=2`. The "old" values (`x_k`, `g_k`) are those currently stored in the history slot that is about to be updated.
        2.  Computes `rho_k = 1.0 / (y_k^T * s_k)` for the newly formed `s_k`, `y_k` pair using `lbfgs_history_last_rho`.
        3.  Computes the new search `direction` using the L-BFGS two-loop recursion algorithm implemented in `lbfgs_history_direction`. This uses the stored history of `s_i`, `y_i`, and `rho_i`.
        4.  Stores the current `variable` (`x_{k+1}`) and `gradient` (`g_{k+1}`) into the history for use in the next iteration by calling `lbfgs_history_push` with `action=1`.

4.  **`lbfgs_release(history)`**:
    *   **Purpose:** Deallocates all dynamically allocated memory within the `lbfgs_history_type` object.
    *   **Functionality:** Iterates through the stored DBCSR matrices in `history%matrix` and releases them using `dbcsr_release`. Then, deallocates the `history%matrix` and `history%rho` arrays.

### Internal (Private) Helper Subroutines

*   **`lbfgs_history_push(history, matrix, vartype, action)`**:
    *   Manages the storage of data into the `history%matrix` circular buffer.
    *   `vartype`: 1 for variable-related data, 2 for gradient-related data.
    *   `action`:
        *   `1`: Directly stores the input `matrix` into the history slot. This is used for `x_k` or `g_k`. The `history%istore(vartype)` is not incremented after this action (it's temporarily incremented then decremented).
        *   `2`: Computes the difference: `history_slot = input_matrix - current_history_slot_content`. This is used to calculate and store `s_k` or `y_k`. `history%istore(vartype)` is incremented.
    *   Handles the circular buffer logic using `MOD(history%istore(vartype) - 1, history%nstore) + 1` to get the current slot index.
    *   Creates DBCSR matrices on the fly during the first `nstore` calls with `action=1`.

*   **`lbfgs_history_last_rho(history)`**:
    *   Calculates `rho_k = 1.0 / (y_k^T * s_k)` for the most recent pair of `s_k` and `y_k` (which are stored as `history%matrix(ispin, istore, 1)` and `history%matrix(ispin, istore, 2)` respectively). `istore` points to this latest pair.

*   **`lbfgs_history_direction(history, gradient, direction)`**:
    *   Implements the standard L-BFGS two-loop recursion to compute the search direction `d = -H_k * g_k`.
    *   `g_k` is the input `gradient`. The output is `direction`.
    *   The algorithm:
        1.  Initialize `q = g_k`.
        2.  First loop (from `k-1` down to `k-m`, where `m` is `nterms`):
            *   `alpha_i = rho_i * s_i^T * q`
            *   `q = q - alpha_i * y_i`
        3.  Compute `gamma_k = s_{k-1}^T * y_{k-1} / (y_{k-1}^T * y_{k-1})`. The initial Hessian approximation `H_0` is taken as `gamma_k * I`.
        4.  `q = gamma_k * q` (this `q` is now `r` in Nocedal's notation, `r = H_0 * q`).
        5.  Second loop (from `k-m` up to `k-1`):
            *   `beta = rho_i * y_i^T * q`
            *   `q = q + s_i * (alpha_i - beta)`
        6.  `direction = -q`.

## Important Variables/Constants

*   `history%nstore`: Key parameter determining the memory usage and the quality of the Hessian approximation.
*   `s_k` (variable differences) and `y_k` (gradient differences): Stored in `history%matrix`. These are fundamental to L-BFGS.
*   `rho_k`: Stored in `history%rho`, essential for the two-loop recursion.

## Usage Examples

The L-BFGS routines are typically used within an iterative optimization algorithm like an SCF procedure.

```fortran
TYPE(lbfgs_history_type) :: lbfgs_opt
TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: current_x, current_g, search_d

INTEGER :: n_spins, lbfgs_mem_depth
! ... initialize n_spins, lbfgs_mem_depth ...
! ... allocate and initialize current_x, current_g, search_d with appropriate templates ...

CALL lbfgs_create(lbfgs_opt, n_spins, lbfgs_mem_depth)

! Initial point (x_0, g_0)
! ... compute initial current_x and current_g ...
CALL lbfgs_seed(lbfgs_opt, current_x, current_g)

DO iter = 1, max_iterations
  ! ... (Optional) Line search: current_x = previous_x + step_alpha * search_d ...
  ! ... Compute new current_g at current_x ...

  ! Get new search direction
  CALL lbfgs_get_direction(lbfgs_opt, current_x, current_g, search_d)

  ! ... Check for convergence based on current_g or change in current_x ...
  ! ... Update current_x using search_d (e.g., current_x = current_x + search_d if step size is 1) ...
END DO

CALL lbfgs_release(lbfgs_opt)
```

## Dependencies and Interactions

*   **`cp_dbcsr_api` and `cp_dbcsr_contrib`**: Heavily used for all operations involving DBCSR matrices, which represent the variables, gradients, and their historical differences. This includes creation, copying, scaling, addition, and dot products.
*   **`kinds`**: Provides the `dp` kind parameter for double-precision real numbers.
*   **`cp_log_handling`**: For logging (though largely commented out in the provided source).
*   **ALMO-SCF routines (e.g., in `almo_scf_optimizer.F` or `almo_scf.F`):** This L-BFGS implementation would be called by higher-level optimization routines within the ALMO-SCF framework. The "variables" could be ALMO coefficients, and "gradients" would be the corresponding energy gradients with respect to these variables.

This module provides a self-contained L-BFGS optimizer capable of handling distributed sparse matrices, making it a potentially powerful tool for optimizing variables in large quantum chemical calculations performed with CP2K.
```
