# almo_scf_optimizer Module Documentation

## Overview

The `almo_scf_optimizer` module in CP2K houses the core optimization algorithms used in the Absolutely Localized Molecular Orbitals Self-Consistent Field (ALMO-SCF) process and its extensions like XALMO (eXtended ALMO) and NLMO (Non-orthogonal Localized Molecular Orbitals) construction. It provides routines for iteratively refining ALMO coefficients (or related variational parameters) to minimize the energy or an appropriate objective function.

The module implements several optimization strategies:
*   DIIS-accelerated SCF for block-diagonal ALMOs.
*   Eigensolver-based methods for XALMOs (often for perturbative delocalization).
*   Iterative methods like Preconditioned Conjugate Gradient (PCG) and Trust Region for both block-diagonal ALMOs and XALMOs.
*   PCG-based optimization for constructing NLMOs by minimizing a localization functional with an orthogonality penalty.

These optimizers interact closely with other `almo_scf_*` modules, particularly `almo_scf_methods` for matrix manipulations and `almo_scf_qs` for interfacing with the broader QS (Quickstep) environment.

## Key Components

### Public Optimizer Subroutines

1.  **`almo_scf_block_diagonal(qs_env, almo_scf_env, optimizer)`**:
    *   **Purpose:** Performs the SCF optimization for strictly block-diagonal ALMOs (ALMOs localized on their respective fragments/domains).
    *   **Algorithm:** Primarily uses the Direct Inversion in the Iterative Subspace (DIIS) method for convergence acceleration, as implemented in `almo_scf_diis_types`.
    *   **Workflow:**
        1.  Calculates the block-diagonal Kohn-Sham (KS) matrix and the DIIS error vector using `almo_scf_ks_to_ks_blk`.
        2.  Pushes the current KS matrix and error vector into the DIIS history (`almo_scf_diis_push`).
        3.  Checks for convergence based on the error norm.
        4.  If not converged, extrapolates a new KS matrix using `almo_scf_diis_extrapolate`.
        5.  Obtains new ALMO coefficients (`matrix_t_blk`) by diagonalizing the blocks of the (extrapolated) KS matrix (`almo_scf_ks_blk_to_tv_blk`). Smearing is handled here if enabled.
        6.  Constructs the new density matrix from the updated ALMOs using `almo_scf_t_to_proj`.
        7.  Calculates the new total energy and KS matrix using `almo_dm_to_almo_ks`.
        8.  Repeats until convergence or max iterations.

2.  **`almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, optimizer)`**:
    *   **Purpose:** Optimizes eXtended ALMOs (XALMOs), which allow for electron delocalization between fragments, using an approach based on solving an eigenvalue problem. This is often suitable for perturbative or one-shot delocalization corrections.
    *   **Algorithm:** Typically involves DIIS for converging the projected KS matrices on overlapping domains.
    *   **Workflow:**
        1.  Calculates the KS matrix projected onto overlapping domains (`domain_ks_xx`) and the corresponding error vector using `almo_scf_ks_to_ks_xx`.
        2.  Pushes these domain KS matrices and error vectors to DIIS history.
        3.  Checks for convergence.
        4.  If not converged (and not a perturbative one-shot calculation), extrapolates new domain KS matrices.
        5.  Obtains new XALMOs (`matrix_t`) by diagonalizing these domain KS matrices (`almo_scf_ks_xx_to_tv_xx`).
        6.  Updates the global density matrix (`almo_scf_t_to_proj`).
        7.  If not a perturbative treatment, recomputes the global KS matrix and energy.
        8.  Repeats if SCF is required for XALMOs.

3.  **`almo_scf_xalmo_pcg(qs_env, almo_scf_env, optimizer, quench_t, matrix_t_in, matrix_t_out, assume_t0_q0x, perturbation_only, special_case)`**:
    *   **Purpose:** Optimizes ALMOs or XALMOs using a Preconditioned Conjugate Gradient (PCG) algorithm.
    *   **Arguments:**
        *   `optimizer`: Controls PCG settings (tolerances, max iterations, preconditioner type).
        *   `quench_t`: Quenching matrix defining sparsity/locality for XALMOs.
        *   `matrix_t_in`, `matrix_t_out`: Input and output ALMO coefficient matrices.
        *   `assume_t0_q0x`: If true, optimizes `X` in `T = T_0 + (1-R_0)*X`, where `T_0` are block-diagonal ALMOs.
        *   `perturbation_only`: If true, the KS matrix is not updated during optimization (first-order correction).
        *   `special_case`: Flags for `xalmo_case_block_diag`, `xalmo_case_fully_deloc`, or `xalmo_case_normal`.
    *   **Workflow:** Iteratively minimizes an objective function (energy + penalties).
        1.  Calculates the objective function value and its gradient using `main_var_to_xalmos_and_loss_func` and `compute_gradient`.
        2.  Computes a search direction using PCG (conjugation coefficient `beta` from `compute_cg_beta`, optional preconditioning via `compute_preconditioner` and `newton_grad_to_step` for full Hessian or `apply_domain_operators` for domain preconditioner).
        3.  Performs a line search along the direction to find an optimal step size.
        4.  Updates the variational parameters (`m_theta` which maps to `matrix_t_out`).
        5.  Repeats until convergence.

4.  **`almo_scf_xalmo_trustr(qs_env, almo_scf_env, optimizer, quench_t, matrix_t_in, matrix_t_out, perturbation_only, special_case)`**:
    *   **Purpose:** Optimizes ALMOs or XALMOs using a Trust Region method.
    *   **Arguments:** Similar to `almo_scf_xalmo_pcg`.
    *   **Workflow:**
        1.  Calculates the objective function and gradient (`main_var_to_xalmos_and_loss_func`, `compute_gradient`).
        2.  Constructs a quadratic model of the objective function within the trust radius.
        3.  Solves the trust region subproblem (e.g., using Steihaug-CG, Dogleg, or Cauchy point methods) to find a trial step (`step`). The model Hessian can be approximated or computed (e.g. via `compute_preconditioner`).
        4.  Evaluates the actual reduction in the objective function with the trial step.
        5.  Compares actual vs. predicted reduction (`rho`) to update the variational parameters (`m_theta`) and adjust the trust radius (`radius_current`).
        6.  Repeats until convergence.

5.  **`almo_scf_construct_nlmos(...)`**:
    *   **Purpose:** Constructs Non-orthogonal Localized Molecular Orbitals (NLMOs) by minimizing a localization functional subject to an orthogonality penalty.
    *   **Algorithm:** Uses a PCG-like approach where the main variational parameters (`m_theta`) define a unitary transformation applied to initial MOs (`matrix_mo_in`).
    *   **Workflow:**
        1.  Calculates the objective function (localization + penalty term) using `compute_obj_nlmos`. The penalty term typically involves `log(det(Sigma))`, where `Sigma` is the NLMO overlap.
        2.  Calculates the gradient of this objective function using `compute_gradient_nlmos`.
        3.  Uses PCG to find a search direction for `m_theta`.
        4.  Performs a line search and updates `m_theta`.
        5.  The strength of the penalty term (`penalty_amplitude`) can be adjusted iteratively to reach a target determinant value for `Sigma`.
        6.  The final NLMOs are `matrix_mo_out = matrix_mo_in * m_theta_normalized`.

### Internal Helper Routines (Selected)

*   **`main_var_to_xalmos_and_loss_func(...)`**: Transforms the primary optimization variables (e.g., elements of `m_theta` or `X`) into ALMO coefficients (`m_t_out`), computes the corresponding energy, and any penalty terms. It also calculates frequently used matrix products like `FTsiginv = F*T*Sigma_inv`.
*   **`compute_gradient(...)`**: Calculates the gradient of the objective function (energy + penalties) with respect to the primary optimization variables.
*   **`compute_preconditioner(...)`**: Constructs a preconditioner matrix (often an approximation to the Hessian inverse) for PCG or Trust Region methods. Can be domain-based or full.
*   **`apply_hessian(...)`**: Computes the Hessian-vector product `H*X`.
*   **`newton_grad_to_step(...)`**: Solves the linear system `H*delta = -grad` using PCG to get the Newton step `delta`.
*   **`compute_cg_beta(...)`**: Computes the `beta` parameter for various CG update formulas (Fletcher-Reeves, Polak-Ribiere, Hestenes-Stiefel, etc.).
*   **`compute_obj_nlmos(...)`**: Calculates the NLMO localization objective function, which includes a localization term (e.g., Berry phase, Pipek-Mezey) and an orthogonality penalty term (e.g., `log(det(Sigma))`).
*   **`compute_gradient_nlmos(...)`**: Computes the gradient of the NLMO objective function.
*   **`compute_xalmos_from_main_var(...)`**: Converts the main optimization variable (e.g., `theta` for NLMOs, or `X` for `T = T0 + (1-R0)X`) into the ALMO/NLMO coefficient matrix `T`.
*   **`xalmo_analysis(...)`**: Performs Energy Decomposition Analysis (EDA) and Charge Transfer Analysis (CTA) based on the XALMO wavefunction.
*   **`wrap_up_xalmo_scf(...)`**: Finalizes XALMO calculations, including energy updates and EDA/CTA if requested.
*   **Printing and Utility Routines**: `fixed_r_report`, `trust_r_report`, `energy_lowering_report`, `print_mathematica_matrix`, `tanh_of_elements`, `dtanh_of_elements`, `inverse_of_elements`, `print_block_sum`.

## Important Variables/Constants

*   `optimizer%optimizer_type`: Flag determining the choice of algorithm (DIIS, PCG, Trust Region).
*   `optimizer%eps_error`: Convergence threshold for the gradient norm.
*   `optimizer%max_iter`: Maximum number of SCF iterations.
*   `almo_scf_env`: The ALMO environment containing matrices like `matrix_s` (AO overlap), `matrix_ks` (KS matrix), `matrix_t_blk` (block-diagonal ALMOs), `matrix_sigma_inv` (ALMO overlap inverse), `quench_t` (XALMO quencher).
*   `m_theta`: Primary variational parameters in PCG/Trust Region methods for XALMOs or NLMOs.
*   `grad`: Gradient of the objective function.
*   `step`: Search direction in PCG/Trust Region.

## Usage Examples

These optimizer routines are called by the main ALMO driver in `almo_scf.F`.
*   `almo_scf_block_diagonal` is called if `almo_scf_env%almo_update_algorithm == almo_scf_diag`.
*   `almo_scf_xalmo_pcg` or `almo_scf_xalmo_trustr` are called by `almo_scf_main` (for block-diagonal ALMOs if `almo_update_algorithm` is PCG/TRUSTR) or by `almo_scf_delocalization` (for XALMOs if `xalmo_update_algorithm` is PCG/TRUSTR).
*   `almo_scf_xalmo_eigensolver` is called by `almo_scf_delocalization` for certain XALMO methods.
*   `almo_scf_construct_nlmos` is called by `construct_nlmos` in `almo_scf.F`.

## Dependencies and Interactions

*   **`almo_scf_types`**: Defines `almo_scf_env_type` and `optimizer_options_type`.
*   **`almo_scf_methods`**: Provides many low-level matrix manipulation and transformation routines used by the optimizers.
*   **`almo_scf_qs`**: For interfacing with the QS environment, e.g., updating the KS matrix from a new density.
*   **`almo_scf_diis_types`**: Provides the DIIS algorithm used in `almo_scf_block_diagonal` and `almo_scf_xalmo_eigensolver`.
*   **`almo_scf_lbfgs_types`**: Provides L-BFGS algorithm (though its direct use isn't immediately apparent in the public interfaces of this module, it might be used internally or planned for future use).
*   **`cp_dbcsr_api`, `cp_dbcsr_contrib`, `cp_dbcsr_cholesky`**: Essential for sparse matrix operations.
*   **`ct_methods`, `ct_types`**: Used for Cayley transformations, particularly in older versions or specific paths of delocalization methods (e.g., `harris_foulkes_correction`).
*   **`qs_loc_utils`, `qs_localization_methods`**: For localization functional components (e.g., Berry phase operator) used in NLMO construction.
*   **Numerical Libraries (LAPACK/BLAS):** Underlying dense matrix operations.

This module forms the computational core of ALMO-based SCF optimizations in CP2K, providing a flexible set of algorithms to achieve self-consistency for various types of localized orbital Ansätze.
```
