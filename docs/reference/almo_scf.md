# almo_scf Module Documentation

## Overview

The `almo_scf` module is the primary driver for Absolutely Localized Molecular Orbitals Self-Consistent Field (ALMO-SCF) calculations within the CP2K package. It manages the entire workflow of an ALMO-based SCF procedure, which aims to obtain molecular orbitals localized on predefined fragments (domains) of a larger system. This approach can offer computational advantages and insights into local electronic structure.

The module handles various stages:
1.  Initialization of the ALMO environment and parameters.
2.  Generation of an initial guess for the ALMOs.
3.  SCF optimization of strictly block-diagonal ALMOs (localized within their domains).
4.  Optional electron delocalization steps to account for inter-domain interactions, using methods like perturbative corrections or eXtended ALMO (XALMO) SCF.
5.  Optional construction of Non-orthogonal Localized Molecular Orbitals (NLMOs).
6.  Post-SCF processing, including storing data for extrapolation and interfacing with standard CP2K property calculations.

The module supports various SCF algorithms (DIIS, PCG, Trust Region), smearing techniques, and different strategies for defining domains and distributing matrices in parallel.

## Key Components

### Main Entry Point

*   **`almo_entry_scf(qs_env, calc_forces)`**:
    *   **Purpose:** This is the main subroutine that initiates and controls the entire ALMO SCF calculation.
    *   **Arguments:**
        *   `qs_env` (TYPE(`qs_environment_type`), POINTER): Pointer to the main Quantum ESPRESSO (QS) environment, providing access to global system information.
        *   `calc_forces` (LOGICAL, INTENT(IN)): A flag indicating whether forces need to be calculated upon completion of the SCF.
    *   **Workflow:**
        1.  Retrieves the `almo_scf_env_type` from `qs_env`.
        2.  Calls `almo_scf_init` to set up the ALMO environment, parameters, and matrices.
        3.  Calls `almo_scf_initial_guess` to generate the starting ALMOs.
        4.  Calls `almo_scf_main` to perform the SCF optimization for block-diagonal ALMOs.
        5.  Calls `almo_scf_delocalization` to optionally include inter-domain electron delocalization effects.
        6.  Calls `construct_nlmos` if requested, to generate NLMOs.
        7.  Calls `almo_scf_post` for post-SCF operations like property calculations and preparing data for the next step.
        8.  Calls `almo_scf_clean_up` to release allocated memory.

### Core Workflow Subroutines

1.  **`almo_scf_init(qs_env, almo_scf_env, calc_forces)`**:
    *   Initializes `almo_scf_env` (the ALMO specific environment). This includes setting up optimizer options, defining domains (molecular or atomic), calculating domain-specific properties (number of basis functions, electrons, virtuals), allocating numerous ALMO-specific matrices (e.g., `matrix_s_blk`, `matrix_t_blk`, `matrix_sigma_blk`), and preparing the AO overlap matrix and its functions (S<sup>-1/2</sup>, S<sup>1/2</sup>) via `almo_scf_init_ao_overlap`. It also constructs the "quencher" matrix used for XALMO methods.

2.  **`almo_scf_initial_guess(qs_env, almo_scf_env)`**:
    *   Generates the initial set of ALMO coefficients (`matrix_t_blk`).
    *   Supports several guess types:
        *   `molecular_guess`: Uses MOs from pre-calculated isolated fragment calculations.
        *   `atomic_guess`: Starts from an atomic density matrix.
        *   `restart_guess`: Reads ALMOs from a restart file.
    *   Can use Always Stable Predictor-Corrector (ASPC) extrapolation if previous step data is available.
    *   Orthogonalizes the initial ALMOs and converts them to projectors/density matrices.
    *   Calculates the initial KS matrix based on this guess using `almo_dm_to_almo_ks`.

3.  **`almo_scf_main(qs_env, almo_scf_env)`**:
    *   Performs the SCF optimization for block-diagonal ALMOs (ALMOs localized entirely on their respective fragments).
    *   Uses an optimization algorithm specified by `almo_scf_env%almo_update_algorithm` (e.g., `almo_scf_diag` for DIIS, `almo_scf_pcg` for Preconditioned Conjugate Gradient, `almo_scf_trustr` for Trust Region). These call routines from `almo_scf_optimizer` module.
    *   Stores the converged block-diagonal KS matrix and MO overlap inverse for later use.

4.  **`almo_scf_delocalization(qs_env, almo_scf_env)`**:
    *   Implements methods to allow electron delocalization between fragments.
    *   The specific method is chosen by `almo_scf_env%deloc_method`:
        *   `almo_deloc_none`: No delocalization.
        *   `almo_deloc_x`, `almo_deloc_scf`, `almo_deloc_x_then_scf`: Methods for full system delocalization (perturbative or SCF).
        *   `almo_deloc_xalmo_1diag`, `almo_deloc_xalmo_x`, `almo_deloc_xalmo_scf`: eXtended ALMO methods that allow controlled delocalization, typically using a quencher matrix to define allowed interactions. These can be perturbative or fully self-consistent.
    *   These routines typically call optimizers from `almo_scf_optimizer` (e.g., `almo_scf_xalmo_pcg`, `almo_scf_xalmo_eigensolver`).

5.  **`construct_nlmos(qs_env, almo_scf_env)`**:
    *   If `almo_scf_env%construct_nlmos` is true, this routine and its helper `construct_nlmos_wrapper` construct Non-orthogonal Localized Molecular Orbitals (NLMOs) from the (potentially delocalized) ALMOs.
    *   Involves optimizing a localization functional using `almo_scf_construct_nlmos` from `almo_scf_optimizer`.
    *   Can construct virtual NLMOs via `construct_virtuals`.
    *   May apply `nlmo_compactification` to filter small coefficients.

6.  **`almo_scf_post(qs_env, almo_scf_env)`**:
    *   Handles tasks after the ALMO SCF has converged.
    *   Stores extrapolation data for subsequent MD steps (`almo_scf_store_extrapolation_data`).
    *   Optionally orthogonalizes the final ALMOs/NLMOs.
    *   Transfers the final MO coefficients to the `qs_env%mos` structure.
    *   Calls `almo_post_scf_compute_properties` (currently a wrapper for `qs_scf_compute_properties`) for standard property calculations.
    *   If forces are calculated, computes the W matrix (Pulay force contribution) using `calculate_w_matrix_almo`.

7.  **`almo_scf_clean_up(almo_scf_env)`**:
    *   Deallocates all matrices and arrays within the `almo_scf_env` to free memory.

### Other Important Subroutines

*   **`almo_scf_store_extrapolation_data(almo_scf_env)`**: Saves relevant matrices (e.g., ALMO coefficients, density matrices) from the current SCF step for use in ASPC extrapolation in subsequent steps.
*   **`almo_scf_print_job_info(almo_scf_env, unit_nr)`**: Prints a detailed summary of the ALMO SCF setup, including optimizer settings, delocalization methods, and fragment statistics (basis functions, electrons, virtuals per fragment).
*   **`almo_scf_init_ao_overlap(matrix_s, almo_scf_env)`**: Initializes ALMO's internal copy of the AO overlap matrix (`matrix_s_blk`) and, if needed by the algorithm, its square root and inverse square root (`matrix_s_blk_sqrt`, `matrix_s_blk_sqrt_inv`).
*   **`almo_scf_env_create_matrices(almo_scf_env, matrix_s0)`**: A utility called by `almo_scf_init` to allocate and create the numerous DBCSR matrices used within the `almo_scf_env_type` based on templates and sizing information.
*   **`construct_virtuals(almo_scf_env)`**: Generates a set of virtual orbitals, orthogonal to the occupied ones, which can then be localized to form virtual NLMOs.
*   **`nlmo_compactification(qs_env, almo_scf_env, matrix)`**: Filters the NLMO coefficient matrix by setting small elements or blocks to zero to enhance sparsity and localization.
*   **`almo_post_scf_compute_properties(qs_env)`**: Currently a wrapper that calls the standard `qs_scf_compute_properties`.

## Important Variables/Constants (within `almo_scf_env_type`)

*   `almo_update_algorithm`: Flag for choosing the SCF optimization algorithm for block-diagonal ALMOs (e.g., `almo_scf_diag`, `almo_scf_pcg`).
*   `deloc_method`: Flag for choosing the electron delocalization method (e.g., `almo_deloc_none`, `almo_deloc_xalmo_scf`).
*   `xalmo_update_algorithm`: Flag for the optimization algorithm used in XALMO methods.
*   `matrix_t_blk`: DBCSR matrix storing block-diagonal ALMO coefficients.
*   `matrix_t`: DBCSR matrix storing (potentially) delocalized ALMO/NLMO coefficients.
*   `matrix_s_blk`: DBCSR matrix for the block-diagonal AO overlap.
*   `matrix_ks`: DBCSR matrix for the Kohn-Sham matrix in AO basis.
*   `matrix_sigma_blk`, `matrix_sigma_inv_0deloc`: ALMO overlap matrix (T<sup>&dagger;</sup>ST) and its inverse for block-diagonal ALMOs.
*   `quench_t`: Quencher matrix used in XALMO methods to define allowed delocalizations.
*   `domain_layout_mos`, `mat_distr_aos`, `mat_distr_mos`: Flags defining how domains are laid out and how matrices are distributed in parallel.
*   `nocc_of_domain`, `nbasis_of_domain`, `nvirt_of_domain`: Arrays describing properties of each domain.
*   `smear`, `smear_e_temp`: Parameters for Fermi-Dirac smearing.
*   Optimizer option types (e.g., `opt_block_diag_diis`, `opt_xalmo_pcg`): Store settings for the various optimizers used.

## Usage Examples

The `almo_entry_scf` subroutine is the single point of entry. It is called by the main CP2K SCF driver when an ALMO-based calculation is requested in the input file (e.g., via `SCF_TYPE ALMO_SCF` in the `QS` section). Users configure the ALMO calculation through specific keywords within the `FORCE_EVAL/DFT/SCF/ALMO_SCF` input section.

## Dependencies and Interactions

The `almo_scf` module interacts with a wide range of other CP2K modules:

*   **`almo_scf_methods`**: Provides low-level methods for ALMO transformations (e.g., density to orbitals, rescaling, orthogonalization).
*   **`almo_scf_optimizer`**: Contains the actual SCF optimization algorithms (DIIS, PCG, Trust Region) tailored for ALMO and XALMO.
*   **`almo_scf_qs`**: Handles the interface between ALMO specific quantities and standard QS quantities (e.g., converting ALMO density to QS KS matrix, constructing QS MOs from ALMOs).
*   **`almo_scf_types`**: Defines the crucial `almo_scf_env_type` and various optimizer option types.
*   **`cp_dbcsr_api` and related modules (`cp_dbcsr_contrib`, `cp_dbcsr_diag`, `cp_dbcsr_operations`):** Extensively used for all distributed sparse matrix operations.
*   **`qs_environment_types`, `qs_mo_types`, `qs_rho_types`, `qs_scf_types`**: For accessing and manipulating core QS data structures like MOs, density matrices, KS matrices, and the SCF environment.
*   **`input_constants`**: Provides named integer constants for various methods and options.
*   **`molecule_types`, `atomic_kind_types`, `particle_types`, `qs_kind_types`**: For system definition and properties.
*   **Parallel Environment (`message_passing`, `cp_blacs_env`):** For parallel execution and matrix distribution.
*   **SCF Utilities (`iterate_matrix`, `qs_initial_guess`, `qs_scf_post_scf`):** For tasks like matrix inversion, initial guess generation, and post-SCF property calculations.

The `almo_scf` module acts as a high-level coordinator, leveraging these specialized modules to perform the complex steps of an ALMO-SCF calculation.
```
