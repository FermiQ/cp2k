# admm_methods Module Documentation

## Overview

The `admm_methods` module in CP2K provides a comprehensive suite of subroutines for performing Auxiliary Density Matrix Method (ADMM) calculations, primarily when working with molecular orbitals (MOs). This is in contrast to `admm_dm_methods`, which operates directly on density matrices. The MO-based ADMM involves fitting MO coefficients to an auxiliary basis, potentially purifying these auxiliary MOs or the resulting auxiliary density matrix, calculating properties using this auxiliary representation, and then merging contributions (like the Kohn-Sham matrix or MO derivatives) back to the primary orbital basis representation.

The module includes extensive support for k-point calculations, various purification schemes, and specialized treatments for GAPW (Gaussian and Augmented Plane Wave) methods, as well as the ADMMQ, ADMMP, and ADMMS variants (Merlot et al., J. Chem. Theory Comput. 2014, 10, 12, 5590-5600).

## Key Components

### Core ADMM Workflow (MO-based)

1.  **`admm_mo_calc_rho_aux(qs_env)` / `admm_mo_calc_rho_aux_kp(qs_env)`**:
    *   **Purpose:** Calculate the auxiliary electron density (`rho_aux_fit`) using MOs. The `_kp` version handles k-point sampling.
    *   **Workflow:**
        1.  **Fit MO Coefficients (`admm_fit_mo_coeffs`)**:
            *   If geometry changed, recalculates transformation matrices `A = S_aux_fit^(-1) * S_aux_fit_vs_orb` and `B = Q^T * A` (in `fit_mo_coeffs`). `S_aux_fit` is the auxiliary basis overlap, `S_aux_fit_vs_orb` (or `Q`) is the mixed auxiliary-orbital basis overlap.
            *   Applies an MO purification/fitting scheme based on `admm_env%purification_method`:
                *   `do_admm_purify_mo_no_diag` / `do_admm_purify_cauchy_subspace`: Uses `purify_mo_cholesky` (calculates auxiliary MOs `C_aux = A * C_orb * Lambda_inv_sqrt`).
                *   `do_admm_purify_mo_diag`: Uses `purify_mo_diag` (diagonalizes Lambda, then `C_aux = A * C_orb * R * Eigvals_lambda_inv_sqrt`).
                *   Default (`do_admm_purify_none`): Uses `purify_mo_none` (`C_aux = A * C_orb`).
        2.  **Calculate Auxiliary Density Matrix**:
            *   If `admm_env%block_dm` is true, calls `blockify_density_matrix` to construct a blocked auxiliary DM.
            *   Otherwise, calls `calculate_dm_mo_no_diag` to form `P_aux = C_hat_occ * Lambda_inv * C_hat^T` (where `C_hat` are the (potentially purified) auxiliary MOs, and `Lambda_inv` comes from the MO purification step).
        3.  **Purify Auxiliary DM (Optional)**: If `admm_env%purification_method == do_admm_purify_cauchy`, calls `purify_dm_cauchy` to apply Cauchy purification to `P_aux`.
        4.  **Calculate Electron Density**: Calls `calculate_rho_elec` to compute the real-space (`rho_r_aux`) and reciprocal-space (`rho_g_aux`) electron density from `P_aux`, using the "AUX_FIT" basis (or "AUX_FIT_SOFT" for GAPW).
        5.  **GAPW Handling**: If `admm_env%do_gapw` is true:
            *   Calculates atomic density coefficients using `calculate_rho_atom_coeff` with the auxiliary basis.
            *   Prepares GAPW densities using `prepare_gapw_den`.
        6.  **ADMMQ/P/S Scaling**: If any of these methods are active, calculates scaling factors `admm_env%gsi` (ratio of electrons in orbital basis vs. auxiliary basis) and scales `P_aux` if needed (for ADMMQ/S).

2.  **`admm_mo_merge_ks_matrix(qs_env)`**:
    *   **Purpose:** Merges the Kohn-Sham (KS) matrix contributions from the auxiliary basis representation back to the primary orbital basis.
    *   **Workflow (depends on `admm_env%purification_method`):**
        *   `do_admm_purify_cauchy`: Calls `merge_ks_matrix_cauchy`.
        *   `do_admm_purify_cauchy_subspace`: Calls `merge_ks_matrix_cauchy_subspace`.
        *   `do_admm_purify_none`: Calls `merge_ks_matrix_none` (or `merge_ks_matrix_none_kp` for k-points). This typically involves `K_orb += A^T * K_aux * A`. For ADMMQ/P/S methods, additional scaling and lambda terms are applied to `K_orb` and the energy.
        *   `do_admm_purify_mo_diag` / `do_admm_purify_mo_no_diag`: Does nothing directly, as merging is often handled via MO derivative merging in OT-like schemes.

### Force and Derivative Calculations

*   **`admm_mo_merge_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, mo_derivs_aux_fit, matrix_ks_aux_fit)`**:
    *   Merges MO derivatives (e.g., `H*C`) from the auxiliary basis to the orbital basis. Used when `purification_method` is `do_admm_purify_mo_diag` (calls `merge_mo_derivs_diag`) or `do_admm_purify_mo_no_diag` (calls `merge_mo_derivs_no_diag`).
*   **`admm_update_ks_atom(qs_env, calculate_forces)`**:
    *   (Primarily for GAPW) Updates the auxiliary KS matrix with atomic contributions from DFT exchange using `update_ks_atom`.
    *   Scales forces appropriately if ADMMS or ADMMP methods are active.
    *   Separates the pure DFT exchange contribution from the HFX contribution in the auxiliary KS matrices.
*   **`calc_admm_mo_derivatives(qs_env, mo_derivs)`**:
    *   Calculates the derivatives for auxiliary MOs (e.g., `H_aux * C_aux`) based on the orbital MO derivatives (`H_orb * C_orb`) by calling `admm_mo_merge_derivs`.
*   **`calc_admm_ovlp_forces(qs_env)` / `calc_admm_ovlp_forces_kp(qs_env)`**:
    *   Calculates forces arising from the derivatives of the mixed auxiliary/orbital overlap matrices (`Q`) and auxiliary overlap matrices (`S_aux_fit`).
    *   Calls `calc_mixed_overlap_force` (non-kp) or implements the logic directly (kp version). The k-point version is complex, involving Fourier transforms of matrix products like `S_aux_inv * K_aux * A * P_orb` and `S_aux_inv * K_aux * A * P_orb * A^T`.
*   **`admm_projection_derivative(qs_env, matrix_hz, matrix_pz, fval)`**:
    *   Calculates force contributions for response calculations (e.g., TDDFT). Computes `W_Q = -2 * S_aux_inv*H_Z*A*P_Z` and `W_S = 2 * S_aux_inv*H_Z*A*P_Z*A^T`, then contracts these with overlap derivatives.
*   **`calc_mixed_overlap_force(qs_env)`**: (Non-kp version)
    *   Calculates ADMM overlap forces using various matrix products involving `A`, `S_inv`, `H'` (auxiliary MO derivatives), `Lambda_inv_sqrt`, and terms from MO purification (`R_schur_R_t`). It constructs effective "density" matrices `matrix_w_s` and `matrix_w_q` to contract with overlap derivatives `dS_aux/dR` and `dQ/dR`.
    *   Includes specific terms for ADMMQ/P/S methods.

### K-Point Specific Routines

*   **`kpoint_calc_admm_matrices(qs_env, calculate_forces)`**:
    *   Crucial for k-point calculations. It prepares the k-dependent auxiliary overlap matrix `S_aux_fit(k)` and transformation matrix `A(k) = S_aux_fit(k)^(-1) * Q(k)`.
    *   This involves:
        1.  Fourier transforming the real-space `matrix_s_aux_fit` and `matrix_s_aux_fit_vs_orb` to each k-point using `rskp_transform`.
        2.  Distributing these k-point specific matrices to the correct processor groups.
        3.  Inverting `S_aux_fit(k)` and computing `A(k)`.
        4.  Storing `A(k)` and `S_aux_fit(k)` in the `kpoint_aux_env`.

### Helper and Internal Subroutines (Selected)

*   **`fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed)`**: Calculates `S_inv` (inverse of auxiliary overlap `S_aux_fit`), `Q` (mixed overlap `S_aux_fit_vs_orb`), `A = S_inv * Q`, and `B = Q^T * A`. Handles `block_fit` option.
*   **`purify_mo_cholesky(admm_env, mos, mos_aux_fit)`**: `Lambda = C^T*B*C`, then `C_aux = A*C*Lambda_inv_chol`.
*   **`purify_mo_diag(admm_env, mos, mos_aux_fit)`**: `Lambda = C^T*B*C`, diagonalizes `Lambda = R*Eig*R^T`, then `C_aux = A*C*R*Eig_inv_sqrt`.
*   **`purify_mo_none(admm_env, mos, mos_aux_fit)`**: `C_aux = A*C`.
*   **`calculate_dm_mo_no_diag(...)`**: Computes `P_aux = (C_hat*occ) * Lambda_inv * C_hat^T`. Handles ADMMQ/P/S scaling by `gsi`.
*   **`blockify_density_matrix(...)`**: Copies blocks from orbital DM to auxiliary DM if `admm_env%block_map` is active.
*   **`purify_dm_cauchy(...)`**: Purifies auxiliary DM `P_aux` using `P_new = S_inv * R * M * R^T * S_inv`, where `R` diagonalizes `S_inv_chol * P_aux * S_inv_chol^T` and `M` contains Heaviside step functions of eigenvalues.
*   **`merge_ks_matrix_cauchy/cauchy_subspace/none/none_kp(...)`**: Implement specific KS matrix merging logic. `_none` versions typically compute `K_orb += A^T * K_aux * A` with ADMMQ/P/S adjustments. `_cauchy` involves complex terms with `R_purify` and `M_purify`.
*   **`calc_spin_dep_aux_exch_ener(...)`**: Calculates spin-specific auxiliary exchange and exact exchange energies for ADMMP/S methods at k-points.
*   **`scale_dm(qs_env, rho_ao_orb, scale_back)`**: Scales orbital DM by `gsi` or `1/gsi` for ADMMP force calculations.
*   **`admm_aux_response_density(qs_env, dm, dm_admm)`**: Computes `P_aux_resp = A * P_resp * A^T`.
*   **`delta(x)`, `Heaviside(x)`**: Elemental functions.

## Important Variables/Constants

*   `admm_env%purification_method`: Controls which purification scheme is used (e.g., `do_admm_purify_none`, `do_admm_purify_mo_diag`, `do_admm_purify_cauchy`).
*   `admm_env%A`, `admm_env%B`: Transformation matrices stored as `cp_fm_type`.
*   `admm_env%S_inv`: Inverse of auxiliary overlap matrix, stored as `cp_fm_type`.
*   `admm_env%lambda`, `admm_env%lambda_inv`, `admm_env%lambda_inv_sqrt`: Matrices related to MO purification, stored as `cp_fm_type`.
*   `admm_env%C_hat`: Auxiliary MO coefficients after fitting/purification, stored as `cp_fm_type`.
*   `admm_env%gsi`: Scaling factor for ADMMQ/P/S methods (electrons in orbital basis / electrons in aux basis).
*   `admm_env%lambda_merlot`: Lambda parameter for ADMMQ/P/S energy corrections.
*   `admm_env%do_gapw`: Logical, true if GAPW is used with ADMM.
*   `admm_env%block_dm`, `admm_env%block_fit`: Logicals for using blocked versions of DM or fitting.
*   `kpoints%kp_aux_env(ikp)%kpoint_env%amat`: Stores k-dependent `A` matrix.
*   `kpoints%kp_aux_env(ikp)%kpoint_env%smat`: Stores k-dependent `S_aux_fit` matrix.

## Usage Examples

These routines are called internally by the CP2K SCF driver when ADMM is enabled with MOs.
1.  **Initialization**: `kpoint_calc_admm_matrices` is called early to set up k-dependent `A(k)` and `S_aux(k)`.
2.  **SCF Cycle**:
    *   `admm_mo_calc_rho_aux` (or `_kp`) is called to compute the auxiliary density.
    *   Energies are calculated using this auxiliary density.
    *   The auxiliary KS matrix is formed.
    *   `admm_mo_merge_ks_matrix` is called to add the auxiliary KS contributions to the orbital KS matrix.
    *   If OT is used, `admm_mo_merge_derivs` (via `calc_admm_mo_derivatives`) is used instead of `admm_mo_merge_ks_matrix` for updating MO gradients.
3.  **Forces**: If forces are required, `calc_admm_ovlp_forces` (or `_kp`) is called to compute forces from overlap derivatives. `admm_update_ks_atom` handles atomic contributions for GAPW.

## Dependencies and Interactions

*   **`admm_types`**: Provides `admm_type` and `get_admm_env`.
*   **`qs_environment_types`, `qs_ks_types`, `qs_rho_types`, `qs_mo_types`**: For accessing core QS data like MOs, density matrices, KS matrices, overlap matrices.
*   **`kpoint_types`, `kpoint_methods`**: For handling k-point specific data and transformations (`rskp_transform`, `kpoint_density_matrices`, `kpoint_density_transform`).
*   **`cp_cfm_types`, `cp_fm_types`, `cp_dbcsr_api`, `parallel_gemm_api`**: For various dense and sparse matrix operations (ScaLAPACK-based `cp_fm`, DBCSR).
*   **`qs_collocate_density`**: For `calculate_rho_elec`.
*   **`qs_gapw_densities`, `qs_rho_atom_methods`, `qs_vxc_atom`**: For GAPW specific density and potential calculations.
*   **`input_constants`**: For purification method flags.
*   **`bibliography`**: For citing relevant papers (e.g., Merlot2014).

This module is central to the performance and accuracy of ADMM calculations in CP2K when orbital information is directly used in the ADMM procedure.
```
