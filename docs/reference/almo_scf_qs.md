# almo_scf_qs Module Documentation

## Overview

The `almo_scf_qs` module serves as a vital interface layer between the Absolutely Localized Molecular Orbitals Self-Consistent Field (ALMO-SCF) specific routines and the general Quickstep (QS) environment within CP2K. Its primary responsibilities include:

*   Creating DBCSR matrices with distributions and structures suitable for ALMO calculations.
*   Converting matrices (like density matrices and Kohn-Sham matrices) between the standard QS format/distribution and the ALMO-specific format/distribution.
*   Updating the QS environment with results obtained from ALMO calculations (e.g., new density matrix, total energy).
*   Constructing the "quencher" matrix, which defines the sparsity and interaction patterns for eXtended ALMO (XALMO) methods.
*   Calculating the energy-weighted density matrix (W matrix) required for force calculations in the ALMO framework.

This module ensures that ALMO calculations can seamlessly integrate with the rest of the CP2K QM/MM machinery and leverage standard QS functionalities for tasks like Hamiltonian construction and energy evaluation.

## Key Components

### Public Subroutines

1.  **`matrix_almo_create(matrix_new, matrix_qs, almo_scf_env, name_new, size_keys, symmetry_new, spin_key, init_domains)`**:
    *   **Purpose:** Creates a new DBCSR matrix (`matrix_new`) tailored for ALMO calculations.
    *   **Arguments:**
        *   `matrix_qs`: A template DBCSR matrix from the QS environment (often the AO overlap matrix `S`) used to derive initial distribution information.
        *   `almo_scf_env`: The ALMO SCF environment, providing context like domain definitions and distribution strategies.
        *   `name_new`: A string name for the new matrix.
        *   `size_keys`: A 2-element integer array specifying the type of dimension for rows and columns (e.g., `almo_mat_dim_aobasis` for atomic orbitals, `almo_mat_dim_occ` for occupied ALMOs).
        *   `symmetry_new`: The symmetry type of the new matrix (e.g., `dbcsr_type_symmetric`, `dbcsr_type_no_symmetry`).
        *   `spin_key`: Spin index (1 or 2) if a dimension refers to MOs; 0 for AO basis.
        *   `init_domains` (LOGICAL): If true, initializes only the block-diagonal elements according to the domain structure defined in `almo_scf_env`.
    *   **Functionality:** Sets up the distribution and block sizes of `matrix_new` based on `almo_scf_env` settings (e.g., `mat_distr_aos`, `mat_distr_mos`) and the `size_keys`.

2.  **`almo_scf_construct_quencher(qs_env, almo_scf_env)`**:
    *   **Purpose:** Constructs the quencher matrix (`almo_scf_env%quench_t`). This matrix defines the sparsity pattern for XALMO methods by specifying which inter-fragment interactions (blocks) are allowed.
    *   **Functionality:**
        *   The construction depends on `almo_scf_env%constraint_type`:
            *   `almo_constraint_block_diagonal`: Only diagonal blocks (intra-domain ALMO interactions) are created in `quench_t`.
            *   `almo_constraint_ao_overlap`: Blocks between domains are included if the maximum atomic orbital (AO) overlap between atoms of the respective domains exceeds a threshold (controlled by `almo_scf_env%quencher_s0` and `quencher_s1`).
            *   `almo_constraint_distance`: Blocks are included if the distance between the domains (defined by closest atom-atom distance or other criteria) is within a specified cutoff radius (controlled by `almo_scf_env%quencher_r0_factor`, `quencher_r1_factor`, and atomic radii).
        *   It builds a global list of interacting domain pairs and populates `almo_scf_env%domain_map` and `almo_scf_env%quench_t` accordingly. `quench_t` typically contains 1.0 for allowed blocks and is sparse otherwise.

3.  **`calculate_w_matrix_almo(matrix_w, almo_scf_env)`**:
    *   **Purpose:** Calculates the energy-weighted density matrix `W = RFR`, where `R` is the ALMO density projector (`T * Sigma_inv * T^T * S`) and `F` is the Kohn-Sham matrix. The `W` matrix is essential for computing Pulay forces in the ALMO framework.
    *   **Functionality:** Performs the matrix multiplications: `W_tmp = T * (Sigma_inv * (T^T F T) * Sigma_inv)` and then `W = W_tmp * T^T`. The resulting `W` matrix in ALMO format is then converted to the QS format and stored in the output `matrix_w`.

4.  **`init_almo_ks_matrix_via_qs(qs_env, matrix_ks, mat_distr_aos, eps_filter)`**:
    *   **Purpose:** Initializes the ALMO Kohn-Sham matrix array (`matrix_ks`) using the current KS matrix from the QS environment.
    *   **Functionality:** If the KS matrix does not exist in `qs_env`, it is allocated and initialized (typically to zero). Then, for each spin, the QS KS matrix is converted to the ALMO distribution and format using `matrix_qs_to_almo` and stored in `matrix_ks(ispin)`.

5.  **`almo_scf_update_ks_energy(qs_env, energy, energy_singles_corr)`**:
    *   **Purpose:** Updates the total energy stored in the QS environment (`qs_env%energy%total`).
    *   **Functionality:** Sets `qs_env%energy%total` to the input `energy` (if provided). If `energy_singles_corr` (a correction term, often from perturbative delocalization) is provided, it's added to `qs_env%energy%total` and stored in `qs_env%energy%singles_corr`.

6.  **`construct_qs_mos(qs_env, almo_scf_env)`**:
    *   **Purpose:** Creates and initializes `mo_set_type` objects within the `qs_env` to hold the final ALMO coefficients, making them accessible to other parts of CP2K.
    *   **Functionality:** Allocates `qs_env%mos` (an array of `mo_set_type`) and initializes each element with the appropriate dimensions (number of AOs, number of ALMOs) based on the `almo_scf_env`. This prepares `qs_env%mos` to receive the ALMO coefficients, typically in `cp_fm_type` (dense matrix) format.

7.  **`matrix_qs_to_almo(matrix_qs, matrix_almo, mat_distr_aos)`**:
    *   **Purpose:** Converts a DBCSR matrix from the standard QS distribution (`matrix_qs`) to an ALMO-specific distribution (`matrix_almo`).
    *   **Functionality:**
        *   If `mat_distr_aos` is `almo_mat_distr_atomic`, a direct copy is performed as the distributions are compatible.
        *   If `mat_distr_aos` is `almo_mat_distr_molecular`, `matrix_qs` is first desymmetrized, and then `dbcsr_complete_redistribute` is called to map the data to the molecular block distribution of `matrix_almo`.

8.  **`matrix_almo_to_qs(matrix_almo, matrix_qs, mat_distr_aos)`**:
    *   **Purpose:** Converts a DBCSR matrix from an ALMO-specific distribution (`matrix_almo`) back to the standard QS distribution (`matrix_qs`). This is the inverse operation of `matrix_qs_to_almo`.

9.  **`almo_dm_to_qs_env(qs_env, matrix_p, mat_distr_aos)`**:
    *   **Purpose:** Updates the density matrix in the main QS environment (`qs_env%rho`) with the ALMO density matrix (`matrix_p`).
    *   **Functionality:** Converts `matrix_p` (for each spin) to the QS format using `matrix_almo_to_qs` and stores it in `qs_env%rho%rho_ao`. It then calls `qs_rho_update_rho` to ensure consistency of derived density quantities (like real-space density) and flags that the density has changed.

10. **`almo_dm_to_almo_ks(qs_env, matrix_p, matrix_ks, energy_total, eps_filter, mat_distr_aos, smear, kTS_sum)`**:
    *   **Purpose:** A composite routine that updates the QS KS matrix and energy based on an ALMO density matrix, and then retrieves the updated KS matrix in ALMO format.
    *   **Functionality:**
        1.  Calls `almo_dm_to_qs_ks` to update the QS environment (density, energy, KS matrix).
        2.  Retrieves the updated KS matrix from `qs_env` and converts it to the ALMO format using `matrix_qs_to_almo`, storing it in the output `matrix_ks`.

11. **`almo_dm_to_qs_ks(qs_env, matrix_p, energy_total, mat_distr_aos, smear, kTS_sum)`**:
    *   **Purpose:** Updates the QS environment based on an ALMO density matrix and calculates the total energy.
    *   **Functionality:**
        1.  Calls `almo_dm_to_qs_env` to transfer the ALMO density matrix `matrix_p` to `qs_env`.
        2.  Calls `qs_ks_update_qs_env` which performs the main QS work: building the KS matrix from the new density and calculating the energy.
        3.  Updates the output `energy_total` with `qs_env%energy%total`, adjusting for electronic entropy (`kTS_sum`) if `smear` is active.

## Important Variables/Constants

*   `almo_scf_env_type`: The ALMO SCF environment, frequently used to get distribution schemes, domain information, and matrix templates.
*   `qs_env_type`: The main QS environment.
*   `mat_distr_aos`: Integer flag specifying the AO distribution strategy (atomic or molecular).
*   `almo_mat_dim_*`: Constants defining matrix dimension types (AO, occupied MO, virtual MO).
*   `matrix_p`: ALMO density matrix.
*   `matrix_ks`: ALMO Kohn-Sham matrix.
*   `quench_t`: Quencher matrix for XALMO.

## Usage Examples

These routines are primarily called by the higher-level ALMO SCF driver routines in `almo_scf.F` and `almo_scf_optimizer.F`.

*   `matrix_almo_create` is used during `almo_scf_init` to set up various ALMO-specific matrices.
*   `almo_scf_construct_quencher` is called during `almo_scf_init`.
*   `almo_dm_to_almo_ks` is called repeatedly within SCF loops (e.g., in `almo_scf_initial_guess`, `almo_scf_block_diagonal`) to get the energy and KS matrix corresponding to the current ALMO density.
*   `almo_dm_to_qs_env` is used at the end of an ALMO calculation in `almo_scf_post` to update the main QS density before property calculations.
*   `calculate_w_matrix_almo` is called in `almo_scf_post` if forces are needed.

## Dependencies and Interactions

*   **`almo_scf_types`**: Provides `almo_scf_env_type` and dimension constants.
*   **`qs_environment_types`, `qs_ks_methods`, `qs_rho_methods`, `qs_rho_types`, `qs_mo_types`, `qs_energy_types`, `qs_scf_types`**: For all interactions with the Quickstep environment, including getting/setting density, KS matrix, MOs, energy, and triggering updates.
*   **`cp_dbcsr_api`, `cp_dbcsr_operations`, `cp_dbcsr_cp2k_link`**: For all low-level distributed sparse matrix operations, including creation, copying, redistribution, and allocation.
*   **`input_constants`**: Defines flags for ALMO constraints, domain layouts, and distribution strategies used in quencher construction and matrix creation.
*   **System Definition Modules (`atomic_kind_types`, `cell_types`, `molecule_types`, `particle_types`):** Provide geometric and atomic information necessary for constructing the quencher matrix based on distance or overlap criteria.
*   **`qs_neighbor_list_types`**: Used to iterate through neighbor lists for distance-based quencher construction.

This module is critical for the proper functioning and integration of ALMO SCF methods within the larger CP2K framework.
```
