# admm_types Module Documentation

## Overview

The `admm_types` module in CP2K is responsible for defining the core data structures (Fortran derived types) that encapsulate the environment and parameters for Auxiliary Density Matrix Method (ADMM) calculations. It provides routines for creating, releasing, and accessing (getting/setting) components of these data structures. This module is fundamental for ADMM, particularly when ADMM is applied in conjunction with molecular orbitals (MOs), k-point sampling, and GAPW methods.

## Key Components

### Derived Types

1.  **`eigvals_type`**:
    *   **Purpose:** A simple container for an array of eigenvalues.
    *   **Fields:**
        *   `DATA` (REAL(dp), DIMENSION(:), POINTER): Pointer to a 1D array of double-precision eigenvalues.

2.  **`eigvals_p_type`**:
    *   **Purpose:** A pointer to an `eigvals_type` instance.
    *   **Fields:**
        *   `eigvals` (TYPE(`eigvals_type`), POINTER): Pointer to an `eigvals_type` object.

3.  **`admm_gapw_r3d_rs_type`**:
    *   **Purpose:** Holds data specific to ADMM calculations when combined with the GAPW (Gaussian and Augmented Plane Wave) method.
    *   **Fields:**
        *   `admm_kind_set` (TYPE(`qs_kind_type`), DIMENSION(:), POINTER): Pointer to an array of `qs_kind_type` objects, storing information about atomic kinds (basis sets, potentials, grids) relevant to the ADMM auxiliary basis in a GAPW context.
        *   `local_rho_set` (TYPE(`local_rho_type`), POINTER): Pointer to a `local_rho_type` object, which contains soft and hard components of atomic densities constructed from the auxiliary fitting basis.
        *   `task_list` (TYPE(`task_list_type`), POINTER): Pointer to a `task_list_type` object, used for managing parallel operations related to soft density plane wave calculations in GAPW.
        *   `oce` (TYPE(`oce_matrix_type`), POINTER): Pointer to `oce_matrix_type` storing precomputed OCE (Orthogonalized Corrected Ensemble) integrals, if applicable.

4.  **`admm_type`**:
    *   **Purpose:** This is the main derived type that encapsulates the entire ADMM environment. It holds various matrices, control parameters, and links to other data structures.
    *   **Key Fields (selected):**
        *   **Dense Matrices (`cp_fm_type` pointers):**
            *   `S_inv`, `S`: Auxiliary basis overlap matrix and its inverse.
            *   `Q`: Mixed overlap matrix (auxiliary basis vs. orbital basis).
            *   `A`: Transformation matrix, typically `S_inv * Q`.
            *   `B`: Transformation matrix, typically `Q^T * A`.
            *   Numerous work matrices (e.g., `work_orb_orb`, `work_aux_orb`, `work_aux_aux`) for intermediate calculations.
            *   Matrices related to MO purification/fitting: `lambda`, `lambda_inv`, `lambda_inv_sqrt`, `R` (eigenvectors of lambda), `C_hat` (fitted/purified auxiliary MOs).
            *   `H`, `K`: Kohn-Sham related matrices in the auxiliary basis.
        *   **Eigenvalue Structures (`eigvals_p_type` pointers):**
            *   `eigvals_lambda`: Eigenvalues of the `lambda` matrix.
            *   `eigvals_P_to_be_purified`: Eigenvalues of the density matrix undergoing purification.
        *   **Control Parameters & Flags:**
            *   `purification_method` (INTEGER): Flag indicating the ADMM purification scheme.
            *   `charge_constrain` (LOGICAL): Flag for charge constrained ADMM.
            *   `do_admmp`, `do_admmq`, `do_admms` (LOGICAL): Flags for specific ADMM variants (Merlot et al., 2014).
            *   `block_dm`, `block_fit` (LOGICAL): Flags for using blocked ADMM approaches.
            *   `block_map` (INTEGER, DIMENSION(:,:), POINTER): Map defining active blocks for blocked methods.
        *   **Dimensions & Counts:**
            *   `nao_orb`, `nao_aux_fit`: Number of atomic orbitals in orbital and auxiliary fitting bases.
            *   `nmo(2)`: Number of molecular orbitals per spin.
        *   **GAPW Specific:**
            *   `admm_gapw_env` (TYPE(`admm_gapw_r3d_rs_type`), POINTER): Pointer to the GAPW specific data.
            *   `do_gapw` (LOGICAL): True if ADMM is used with GAPW.
        *   **Density Matrix ADMM Integration:**
            *   `admm_dm` (TYPE(`admm_dm_type`), POINTER): Pointer to an `admm_dm_type` instance for density matrix-based ADMM.
        *   **K-Point Support (`kpoint_transitional_type`):**
            *   `matrix_ks_aux_fit`, `matrix_ks_aux_fit_dft`, `matrix_ks_aux_fit_hfx`: Auxiliary KS matrix and its DFT/HFX components.
            *   `matrix_s_aux_fit`, `matrix_s_aux_fit_vs_orb`: Auxiliary overlap and mixed overlap matrices. These types can hold either a single matrix (Gamma point) or distributed k-point matrices.
        *   **Other QS Components:**
            *   `mos_aux_fit` (POINTER to `mo_set_type` array): Auxiliary basis MOs.
            *   `sab_aux_fit`, `sab_aux_fit_asymm`, `sab_aux_fit_vs_orb` (POINTER to `neighbor_list_set_p_type` array): Neighbor lists for different basis combinations.
            *   `rho_aux_fit`, `rho_aux_fit_buffer` (POINTER to `qs_rho_type`): Auxiliary density and a buffer.
            *   `task_list_aux_fit` (POINTER to `task_list_type`).
            *   `xc_section_primary`, `xc_section_aux` (POINTER to `section_vals_type`): XC functional information.

### Public Subroutines

1.  **`admm_env_create(admm_env, admm_control, mos, para_env, natoms, nao_aux_fit, blacs_env_ext)`**:
    *   **Purpose:** Allocates and initializes an `admm_type` instance (`admm_env`).
    *   **Arguments:**
        *   `admm_env` (Output): Pointer to the `admm_type` to be created.
        *   `admm_control` (Input): ADMM control settings from input.
        *   `mos` (Input): Molecular orbitals of the primary basis set.
        *   `para_env` (Input): Parallel environment.
        *   `natoms` (Input): Number of atoms.
        *   `nao_aux_fit` (Input): Number of AOs in the auxiliary fitting basis.
        *   `blacs_env_ext` (Optional Input): External BLACS context.
    *   **Functionality:** Sets up the `admm_env` by allocating memory for its numerous components (especially the many `cp_fm_type` dense matrices) and initializing them based on dimensions derived from inputs and control parameters from `admm_control`.

2.  **`admm_env_release(admm_env)`**:
    *   **Purpose:** Deallocates all components of an `admm_type` instance to prevent memory leaks.
    *   **Functionality:** Systematically releases all `cp_fm_type` matrices, deallocates eigenvalue arrays, section values, and calls specialized release routines for nested types like `admm_gapw_env_release`, `admm_dm_release`, `deallocate_mo_set`, `release_neighbor_list_sets`, `kpoint_transitional_release`, `qs_rho_release`, and `deallocate_task_list`.

3.  **`admm_gapw_env_release(admm_gapw_env)`**:
    *   **Purpose:** Specifically deallocates the contents of an `admm_gapw_r3d_rs_type` instance.
    *   **Functionality:** Releases `admm_kind_set`, `local_rho_set`, `task_list`, and `oce` components.

4.  **`get_admm_env(admm_env, ...)`**:
    *   **Purpose:** Provides read-only access to the various pointers within an `admm_type` instance.
    *   **Arguments:** An input `admm_env` and a list of optional output pointer arguments corresponding to fields within `admm_type`.
    *   **Functionality:** If an optional argument is present, the corresponding pointer from `admm_env` is assigned to it. It uses helper functions from `kpoint_transitional` (`get_1d_pointer`, `get_2d_pointer`) to correctly retrieve pointers from `kpoint_transitional_type` members.

5.  **`set_admm_env(admm_env, ...)`**:
    *   **Purpose:** Allows associating external data (via pointers) with the fields within an `admm_type` instance.
    *   **Arguments:** An input/output `admm_env` and a list of optional input pointer arguments.
    *   **Functionality:** If an optional argument is present, the corresponding pointer field in `admm_env` is set to this input pointer. It uses helper functions from `kpoint_transitional` (`set_1d_pointer`, `set_2d_pointer`) for `kpoint_transitional_type` members.

## Important Variables/Constants

*   `moduleN` (Character, Private, Parameter): Name of the module, 'admm_types'.
*   Various constants from `input_constants` are used internally during `admm_env_create` to set flags like `purification_method`, `charge_constrain`, etc., based on `admm_control_type` settings.

## Usage Examples

The `admm_type` is a central data structure managed by higher-level ADMM driver routines.

```fortran
TYPE(admm_type), POINTER :: my_admm_environment
TYPE(admm_control_type), POINTER :: control_settings
TYPE(mo_set_type), DIMENSION(:), POINTER :: orbital_mos
! ... other necessary types ...

! Creation
CALL admm_env_create(my_admm_environment, control_settings, orbital_mos, &
                     current_para_env, num_atoms, num_aux_orbitals)

! Accessing components (e.g., getting the A matrix)
TYPE(cp_fm_type), POINTER :: matrix_A
CALL get_admm_env(my_admm_environment, A=matrix_A)
IF (ASSOCIATED(matrix_A)) THEN
  ! Use matrix_A
END IF

! Setting components (e.g., associating an external auxiliary MO set)
TYPE(mo_set_type), DIMENSION(:), POINTER :: external_aux_mos
! ... external_aux_mos is allocated and filled ...
CALL set_admm_env(my_admm_environment, mos_aux_fit=external_aux_mos)

! Release
CALL admm_env_release(my_admm_environment)
```

## Dependencies and Interactions

*   **`admm_dm_types`**: For the nested `admm_dm_type` and its release routine.
*   **`cp_control_types`**: Provides `admm_control_type` which dictates many settings for `admm_type`.
*   **Matrix Libraries (`cp_blacs_env`, `cp_fm_struct`, `cp_fm_types`, `cp_dbcsr_api`, `cp_dbcsr_operations`):** Essential for creating, managing, and releasing the numerous dense (`cp_fm_type`) and sparse (`dbcsr_p_type`) matrices stored within `admm_type`.
*   **`input_constants`**: Defines flags used to configure ADMM behavior.
*   **`kpoint_transitional`**: Provides mechanisms (`get_1d/2d_pointer`, `set_1d/2d_pointer`, `kpoint_transitional_release`) for handling matrix types that can be either Gamma-point or k-point distributed. Several key matrices in `admm_type` are of `kpoint_transitional_type`.
*   **Various QS Modules (`qs_kind_types`, `qs_local_rho_types`, `qs_mo_types`, `qs_neighbor_list_types`, `qs_oce_types`, `qs_rho_types`, `task_list_types`):** `admm_type` stores pointers to many objects defined in these modules (e.g., `mos_aux_fit`, `rho_aux_fit`), and their respective deallocation/release routines are called during `admm_env_release`.
*   **`admm_methods`**: This module contains the algorithmic parts of ADMM that heavily use the `admm_type` data structure created and managed here.

The `admm_types` module acts as the data backbone for ADMM calculations, providing a structured way to manage the complex set of parameters, matrices, and connections to other parts of the CP2K code.
```
