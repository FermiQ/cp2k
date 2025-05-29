# admm_dm_methods Module Documentation

## Overview

The `admm_dm_methods` module implements Auxiliary Density Matrix Method (ADMM) functionalities within CP2K that specifically operate on density matrices. ADMM is a technique used to speed up electronic structure calculations, particularly density functional theory (DFT), by introducing an auxiliary basis set to represent certain quantities, thereby reducing computational bottlenecks. This module provides routines to calculate an auxiliary density matrix from the primary density matrix and to merge contributions from the auxiliary Kohn-Sham (KS) matrix back into the primary KS matrix. It also includes methods for McWeeny purification of the auxiliary density matrix and its reversal.

## Key Components

The module exposes two main public subroutines:

1.  **`admm_dm_calc_rho_aux(qs_env)`**:
    *   **Purpose:** Calculates the auxiliary density matrix (`rho_ao_aux`) from the primary density matrix (`rho_ao`).
    *   **Workflow:**
        1.  Retrieves the ADMM environment (`admm_dm`) from `qs_env`.
        2.  Selects the method for mapping the primary density matrix to the auxiliary one based on `admm_dm%method`:
            *   `do_admm_basis_projection`: Calls `map_dm_projection`.
            *   `do_admm_blocked_projection`: Calls `map_dm_blocked`.
        3.  If `admm_dm%purify` is true, it applies McWeeny purification to the auxiliary density matrix by calling `purify_mcweeny`.
        4.  Calls `update_rho_aux` to compute the electron density in real and reciprocal space (`rho_r_aux`, `rho_g_aux`) from the (potentially purified) auxiliary density matrix.

2.  **`admm_dm_merge_ks_matrix(qs_env)`**:
    *   **Purpose:** Merges the auxiliary Kohn-Sham (KS) matrix (`matrix_ks_aux_fit`) into the primary KS matrix (`matrix_ks`).
    *   **Workflow:**
        1.  Retrieves the ADMM environment.
        2.  If `admm_dm%purify` was used for `rho_aux`, it calls `revert_purify_mcweeny` to obtain the appropriate KS matrix contribution (`matrix_ks_merge`) by reversing the purification steps. Otherwise, `matrix_ks_merge` is directly `matrix_ks_aux_fit`.
        3.  Selects the method for merging based on `admm_dm%method`:
            *   `do_admm_basis_projection`: Calls `merge_dm_projection`.
            *   `do_admm_blocked_projection`: Calls `merge_dm_blocked`.
        4.  If purification was involved, deallocates the temporary `matrix_ks_merge`.

### Private Helper Subroutines

*   **`map_dm_projection(qs_env)`**:
    *   Calculates `P_aux = A * P * A^T`.
    *   `P` is the primary density matrix (`rho_ao`).
    *   `A = S_aux_fit^(-1) * S_aux_fit_vs_orb`, where `S_aux_fit` is the overlap matrix of the auxiliary fitting basis and `S_aux_fit_vs_orb` is the overlap matrix between the auxiliary fitting basis and the orbital basis.
    *   The matrix `A` is stored in `admm_dm%matrix_A` and reused if the overlap matrices haven't changed.

*   **`map_dm_blocked(qs_env)`**:
    *   Constructs `P_aux` by copying specific blocks from the primary density matrix `P`.
    *   The blocks to be copied are determined by `admm_dm%block_map`, a logical mask.

*   **`update_rho_aux(qs_env)`**:
    *   Takes the auxiliary density matrix `rho_ao_aux` and computes the corresponding electron density on the real-space grid (`rho_r_aux`) and its reciprocal-space representation (`rho_g_aux`).
    *   This is achieved by calling `calculate_rho_elec` using the "AUX_FIT" basis type.

*   **`merge_dm_projection(qs_env, matrix_ks_merge)`**:
    *   Adds the auxiliary KS contribution to the primary KS matrix: `K_orb += A^T * K_aux * A`.
    *   `K_aux` is `matrix_ks_merge`.
    *   `A` is the transformation matrix `admm_dm%matrix_A` calculated in `map_dm_projection`.

*   **`merge_dm_blocked(qs_env, matrix_ks_merge)`**:
    *   Adds `matrix_ks_merge` to the primary KS matrix (`matrix_ks`).
    *   Before addition, blocks in `matrix_ks_merge` that correspond to `admm_dm%block_map(iatom, jatom) == 0` are zeroed out.

*   **`purify_mcweeny(qs_env)`**:
    *   Applies McWeeny purification to `rho_ao_aux` (denoted as `P`) to enforce idempotency (P*S*P = P, for an orthogonal basis S=I, this is P*P=P).
    *   The iterative formula used is `P_new = 3*P_old*S*P_old - 2*P_old*S*P_old*S*P_old*S*P_old`.
    *   The history of `P` matrices at each step is stored in `admm_dm%mcweeny_history` for use in `revert_purify_mcweeny`.
    *   For spin-unpolarized cases, `P` is scaled by 0.5 before purification and 2.0 after.

*   **`revert_purify_mcweeny(qs_env, matrix_ks_merge)`**:
    *   Calculates the effective KS matrix contribution from the auxiliary system when McWeeny purification has been applied. This involves applying the chain rule (propagating derivatives) backward through the McWeeny iteration steps.
    *   It uses the stored `matrix_p` from `admm_dm%mcweeny_history` and calls `reverse_mcweeny_step` iteratively.
    *   The output `matrix_ks_merge` is the transformed KS matrix.

*   **`reverse_mcweeny_step(matrix_k, matrix_s, matrix_p)`**:
    *   Performs one step of the reverse McWeeny iteration. It updates `matrix_k` based on `matrix_s` (overlap) and `matrix_p` (density matrix from a purification step) according to the derivative of the McWeeny formula.
    *   The update involves several `dbcsr_multiply` operations:
        `K_new = 3*K*P*S + 3*S*P*K - 2*(K*P*S*P*S + S*P*K*P*S + S*P*S*K*P*S)` (conceptually, actual implementation uses intermediate matrices).

## Important Variables/Constants

*   `admm_dm%method`: Integer flag determining the ADMM projection/merging strategy.
    *   `do_admm_basis_projection`: Use projection based on overlap matrices.
    *   `do_admm_blocked_projection`: Use block-wise copying/addition.
*   `admm_dm%purify`: Logical flag indicating whether to apply McWeeny purification to the auxiliary density matrix.
*   `admm_dm%matrix_A`: DBCSR matrix storing `S_aux_fit^(-1) * S_aux_fit_vs_orb`, used in projection methods.
*   `admm_dm%block_map`: Integer matrix indicating which atom pairs' blocks are active in the blocked methods.
*   `admm_dm%mcweeny_history`: Array of linked lists (`mcweeny_history_type`) storing density matrices from purification steps, used for the reverse process.
*   `admm_dm%eps_filter`: Threshold for matrix inversion and convergence checks.
*   `rho_ao`: Primary density matrix (orbital basis).
*   `rho_ao_aux`: Auxiliary density matrix (auxiliary fitting basis).
*   `matrix_ks`: Primary Kohn-Sham matrix (orbital basis).
*   `matrix_ks_aux_fit`: Auxiliary Kohn-Sham matrix (auxiliary fitting basis).
*   `matrix_s_aux_fit`: Overlap matrix of the auxiliary fitting basis.
*   `matrix_s_mixed` (aliased as `matrix_s_aux_fit_vs_orb`): Overlap matrix between auxiliary fitting basis and orbital basis.

## Usage Examples

These subroutines are typically called from within the self-consistent field (SCF) cycle in CP2K when ADMM is enabled. They are not directly invoked by users in the input file but are part of the ADMM machinery.

1.  During SCF, after the primary density matrix `rho_ao` is formed, `admm_dm_calc_rho_aux` is called to compute `rho_ao_aux` and then `rho_r_aux`.
2.  The energy and potential contributions from the auxiliary density are calculated.
3.  The auxiliary KS matrix `matrix_ks_aux_fit` is computed.
4.  `admm_dm_merge_ks_matrix` is called to incorporate the auxiliary KS contributions back into the primary `matrix_ks`.

## Dependencies and Interactions

*   **`admm_dm_types`**: Provides `admm_dm_type` and `mcweeny_history_type`.
*   **`admm_types`**: Provides `get_admm_env` to access the ADMM environment.
*   **`cp_dbcsr_api`, `cp_dbcsr_contrib`, `cp_dbcsr_operations`**: Heavily used for all distributed sparse matrix operations (creation, multiplication, addition, norm, etc.).
*   **`qs_environment_types`**: Provides `qs_environment_type` and `get_qs_env` for accessing various parts of the QS setup (like density matrices, KS matrices, DFT control parameters).
*   **`qs_rho_types`**: Provides `qs_rho_type` for handling density information (AO matrices, real-space grids) and `qs_rho_get`/`qs_rho_set` accessors.
*   **`qs_ks_types`**: Provides `qs_ks_env_type`.
*   **`qs_collocate_density`**: Uses `calculate_rho_elec` for computing electron density on grids.
*   **`input_constants`**: Defines constants like `do_admm_basis_projection`.
*   **`iterate_matrix`**: Uses `invert_Hotelling` for matrix inversion.
*   **`cp_log_handling`**: For logging.
*   **`task_list_types`**: For task management in `calculate_rho_elec`.

The module forms a crucial part of the ADMM implementation, interfacing between the primary representation of densities and potentials and their auxiliary counterparts.
```
