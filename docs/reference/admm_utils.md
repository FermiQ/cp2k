# admm_utils Module Documentation

## Overview

The `admm_utils` module in CP2K provides utility subroutines designed to support Auxiliary Density Matrix Method (ADMM) calculations. Specifically, these routines handle corrections to the Kohn-Sham (KS) matrix that may be necessary when certain ADMM purification methods are employed, particularly in relation to eigenvalue computations or specific steps in the Self-Consistent Field (SCF) cycle. The corrections often involve adding or subtracting pre-calculated matrix components stored within the ADMM environment (`admm_type`).

## Key Components

The module contains two public subroutines:

1.  **`admm_correct_for_eigenvalues(ispin, admm_env, ks_matrix)`**
    *   **Purpose:** Modifies the input Kohn-Sham matrix (`ks_matrix`) based on the ADMM purification method specified in `admm_env`. This is typically done to prepare the KS matrix for a process like eigenvalue calculation where certain components added by ADMM might need to be temporarily altered or replaced.
    *   **Arguments:**
        *   `ispin` (Integer, Input): The spin index for which the correction is applied.
        *   `admm_env` (TYPE(`admm_type`), POINTER, Input): The ADMM environment. This structure contains various pre-calculated matrices and settings, including `admm_env%A` (the transformation matrix `S_aux_inv * Q_mixed`), `admm_env%K(ispin)` (auxiliary KS matrix), `admm_env%H_corr(ispin)` (a corrected Hamiltonian term, often `A^T*K_aux*A`), and `admm_env%ks_to_be_merged(ispin)` (a KS component that was previously merged).
        *   `ks_matrix` (TYPE(`dbcsr_type`), POINTER, Input/Output): The Kohn-Sham matrix in the orbital basis that will be modified.
    *   **Functionality:**
        *   The routine only performs corrections if `admm_env%block_dm` (blocked density matrix ADMM) is `.FALSE.`.
        *   The behavior depends on `admm_env%purification_method`:
            *   **`do_admm_purify_cauchy_subspace`**:
                1.  A temporary sparse matrix `work` is created.
                2.  `admm_env%ks_to_be_merged(ispin)` (a previously added component) is copied to `work`.
                3.  `work` is subtracted from `ks_matrix`.
                4.  The term `H_corr_tmp = A^T * K_aux(ispin) * A` is computed (where `K_aux(ispin)` is `admm_env%K(ispin)`) and stored in `admm_env%H_corr(ispin)`.
                5.  `admm_env%H_corr(ispin)` is copied to `work`.
                6.  `work` (now `H_corr_tmp`) is added to `ks_matrix`.
            *   **`do_admm_purify_mo_diag`**:
                1.  A temporary sparse matrix `work` is created.
                2.  The term `H_corr_tmp = A^T * K_aux(ispin) * A` is computed and stored in `admm_env%H_corr(ispin)`.
                3.  `admm_env%H_corr(ispin)` is copied to `work`.
                4.  `work` (now `H_corr_tmp`) is added to `ks_matrix`.
                *(Note: The subtraction of a previously merged component like `ks_to_be_merged` seems implied for consistency but is not explicitly performed on `ks_matrix` in this path before adding `H_corr` in the provided code, unlike the `cauchy_subspace` case. This might rely on `ks_matrix` not yet having that component or it being handled elsewhere.)*
            *   **`do_admm_purify_mo_no_diag`, `do_admm_purify_none`, `do_admm_purify_cauchy`**: No action is taken.

2.  **`admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_matrix)`**
    *   **Purpose:** Reverts the modifications made by `admm_correct_for_eigenvalues` to restore the Kohn-Sham matrix to its state prior to the correction.
    *   **Arguments:** Same as `admm_correct_for_eigenvalues`.
    *   **Functionality:**
        *   The routine only performs un-corrections if `admm_env%block_dm` is `.FALSE.`.
        *   The behavior depends on `admm_env%purification_method`:
            *   **`do_admm_purify_cauchy_subspace`**:
                1.  A temporary sparse matrix `work` is created.
                2.  `admm_env%H_corr(ispin)` (which is `A^T * K_aux(ispin) * A`) is copied to `work`.
                3.  `work` is subtracted from `ks_matrix`.
                4.  `admm_env%ks_to_be_merged(ispin)` is copied to `work`.
                5.  `work` (now `ks_to_be_merged`) is added back to `ks_matrix`.
            *   **`do_admm_purify_mo_diag`**:
                1.  A temporary sparse matrix `work` is created.
                2.  `admm_env%H_corr(ispin)` is copied to `work`.
                3.  `work` is subtracted from `ks_matrix`.
            *   **`do_admm_purify_mo_no_diag`, `do_admm_purify_none`, `do_admm_purify_cauchy`**: No action is taken.

## Important Variables/Constants

*   `admm_env%purification_method`: An integer flag (from `input_constants`) determining the type of ADMM purification being used. The behavior of the utility functions is conditional on this value.
*   `admm_env%block_dm`: A logical flag. If true, these utilities generally do not apply corrections.
*   `admm_env%A`: The transformation matrix `S_aux_inv * Q_mixed` stored as a dense matrix (`cp_fm_type`).
*   `admm_env%K(ispin)`: The Kohn-Sham matrix in the auxiliary basis for the given spin, stored as a dense matrix.
*   `admm_env%H_corr(ispin)`: A dense matrix used to store `A^T * K_aux(ispin) * A`.
*   `admm_env%ks_to_be_merged(ispin)`: A dense matrix component that is part of the KS matrix in certain ADMM schemes.
*   `ks_matrix`: The primary Kohn-Sham matrix in the orbital basis, represented as a DBCSR sparse matrix.

## Usage Examples

These subroutines are not typically called directly by the user but are used internally within the SCF cycle or related computational steps (e.g., eigenvalue solvers, response property calculations) when ADMM is active.

**Conceptual Workflow:**

```fortran
! Inside an SCF iteration or a specific calculation step:
! ... ks_matrix is in a certain state ...

! If eigenvalues are needed and specific ADMM purification is on:
CALL admm_correct_for_eigenvalues(current_spin, admm_environment, orbital_ks_matrix)

! ... perform operations requiring the corrected ks_matrix (e.g., solve eigenproblem) ...

! Revert the KS matrix to its previous state:
CALL admm_uncorrect_for_eigenvalues(current_spin, admm_environment, orbital_ks_matrix)

! ... continue with the SCF or other calculations ...
```

## Dependencies and Interactions

*   **`admm_types`**: Provides the `admm_type` data structure, which is a primary input to these utilities and contains the necessary matrices and flags.
*   **`cp_dbcsr_api`**: Used for creating, copying, adding, setting to zero, and deallocating DBCSR sparse matrices (`ks_matrix`, `work`).
*   **`cp_dbcsr_operations`**: Specifically, `copy_fm_to_dbcsr` is used to convert dense matrices from `admm_env` to sparse format for operations with `ks_matrix`.
*   **`input_constants`**: Defines the integer flags for different ADMM purification methods (e.g., `do_admm_purify_cauchy_subspace`, `do_admm_purify_mo_diag`).
*   **`parallel_gemm_api`**: Used for dense matrix multiplications (`parallel_gemm`) to compute terms like `A^T * K_aux * A`.
*   **`admm_methods`**: The routines in `admm_methods.F` are responsible for the main ADMM calculations and likely set up the state of `admm_env` (including matrices like `K(ispin)`, `H_corr(ispin)`, `ks_to_be_merged(ispin)`) that these utilities depend on. The `admm_utils` routines are helpers for procedures orchestrated by `admm_methods`.

The `admm_utils` module plays a supporting role by providing mechanisms to temporarily alter and restore the Kohn-Sham matrix, facilitating compatibility between ADMM's modified Hamiltonian and standard operations like eigenvalue solvers under specific ADMM purification schemes.
```
