# admm_dm_types Module Documentation

## Overview

The `admm_dm_types` module in CP2K defines the data structures (Fortran derived types) and associated lifecycle management subroutines (create and release) for Auxiliary Density Matrix Method (ADMM) calculations that specifically involve density matrices. These types encapsulate settings, control parameters, and matrices required by the ADMM density matrix (`_dm`) methods.

## Key Components

### Derived Types

1.  **`mcweeny_history_type`**
    *   **Purpose:** Represents a single entry in a linked list used to store the history of density matrices during McWeeny purification.
    *   **Fields:**
        *   `m` (TYPE(`dbcsr_type`)): A DBCSR (Distributed Block Compressed Sparse Row) matrix, typically storing a density matrix from a specific McWeeny iteration.
        *   `count` (INTEGER): An integer identifier for the iteration step in the McWeeny process. Default: -1.
        *   `next` (TYPE(`mcweeny_history_type`), POINTER): Pointer to the subsequent `mcweeny_history_type` entry in the linked list. Default: `Null()`.

2.  **`mcweeny_history_p_type`**
    *   **Purpose:** A wrapper type that points to the head of a `mcweeny_history_type` linked list.
    *   **Fields:**
        *   `p` (TYPE(`mcweeny_history_type`), POINTER): Pointer to the first element of the McWeeny history list. Default: `Null()`.

3.  **`admm_dm_type`**
    *   **Purpose:** The primary data structure holding all necessary information and matrices for ADMM calculations involving density matrices.
    *   **Fields:**
        *   `purify` (LOGICAL): A flag indicating whether McWeeny purification should be applied to the auxiliary density matrix. Default: `.FALSE.`. Set based on `admm_control%purification_method`.
        *   `method` (INTEGER): An integer flag specifying the ADMM method to be employed (e.g., basis projection, blocked projection). Default: -1. Set from `admm_control%method`.
        *   `matrix_a` (TYPE(`dbcsr_type`), POINTER): Pointer to a DBCSR matrix. This typically stores the transformation matrix `A = S_aux_fit^(-1) * S_aux_fit_vs_orb` used in basis projection ADMM methods. Default: `Null()`.
        *   `eps_filter` (REAL(KIND=dp)): A small floating-point value used as a threshold in numerical operations, such as matrix inversion or filtering small matrix elements. Default: `1e-20_dp`. Set from `admm_control%eps_filter`.
        *   `mcweeny_max_steps` (INTEGER): The maximum number of iterations allowed for the McWeeny purification procedure. Default: 100.
        *   `block_map` (INTEGER, DIMENSION(:, :), POINTER): A 2D integer array that defines a mapping for block-based ADMM methods. `block_map(i, j) = 1` typically means that the block corresponding to atom `i` and atom `j` is considered "active" or part of a user-defined region. Default: `Null()`. Allocated and populated if `method` is not `do_admm_basis_projection`.
        *   `mcweeny_history` (TYPE(`mcweeny_history_p_type`), DIMENSION(:), POINTER): An allocatable array of `mcweeny_history_p_type`. Its size corresponds to the number of spin channels (`nspins`), allowing storage of McWeeny purification history for each spin component. Default: `Null()`.

### Public Subroutines

1.  **`admm_dm_create(admm_dm, admm_control, nspins, natoms)`**
    *   **Purpose:** Allocates and initializes an instance of `admm_dm_type`.
    *   **Arguments:**
        *   `admm_dm` (TYPE(`admm_dm_type`), POINTER, Output): The pointer that will hold the newly created `admm_dm_type` instance.
        *   `admm_control` (TYPE(`admm_control_type`), POINTER, Input): Pointer to the ADMM control settings (typically derived from user input).
        *   `nspins` (INTEGER, Input): The number of spin components in the calculation (e.g., 1 for restricted, 2 for unrestricted).
        *   `natoms` (INTEGER, Input): The total number of atoms in the system.
    *   **Functionality:**
        *   Allocates memory for `admm_dm`.
        *   Copies relevant settings from `admm_control` (like `purification_method`, `method`, `eps_filter`) into the corresponding fields of `admm_dm`.
        *   Allocates the `mcweeny_history` array with `nspins` elements.
        *   If the ADMM method specified is not basis projection (`do_admm_basis_projection`), it allocates `block_map` (size `natoms` x `natoms`) and populates it based on the `blocks` defined in `admm_control`. The `block_map(i,j)` is set to 1 if atoms `i` and `j` belong to the same defined block in the input.

2.  **`admm_dm_release(admm_dm)`**
    *   **Purpose:** Deallocates all memory associated with an `admm_dm_type` instance to prevent memory leaks.
    *   **Arguments:**
        *   `admm_dm` (TYPE(`admm_dm_type`), POINTER, Input/Output): The `admm_dm_type` instance to be released.
    *   **Functionality:**
        *   Checks if `admm_dm` is associated; if not, it returns immediately.
        *   If `admm_dm%matrix_a` is associated, it calls `dbcsr_release` on the matrix and then deallocates the pointer itself.
        *   If `admm_dm%block_map` is associated, it deallocates it.
        *   Deallocates the `admm_dm%mcweeny_history` array.
        *   Finally, deallocates the `admm_dm` pointer itself.
        *   Note: The routine does not explicitly traverse and deallocate the `mcweeny_history` linked lists. This is typically handled by the `revert_purify_mcweeny` subroutine in `admm_dm_methods.F` as it consumes the history.

## Important Variables/Constants

*   `moduleN` (Character, Private, Parameter): Name of the module, 'admm_dm_types'.
*   Constants from `input_constants` module:
    *   `do_admm_basis_projection`: Integer flag identifying the basis projection ADMM method.
    *   `do_admm_purify_mcweeny`: Integer flag identifying McWeeny purification.

## Usage Examples

Instances of `admm_dm_type` are created and managed by higher-level ADMM routines within CP2K.

```fortran
TYPE(admm_dm_type), POINTER :: my_admm_dm
TYPE(admm_control_type), POINTER :: my_admm_control_params
INTEGER :: num_spins, num_atoms

! ... (my_admm_control_params, num_spins, num_atoms are set) ...

! Create the admm_dm object
CALL admm_dm_create(my_admm_dm, my_admm_control_params, num_spins, num_atoms)

! ... (use my_admm_dm in ADMM calculations) ...

! Release the admm_dm object
CALL admm_dm_release(my_admm_dm)
```

## Dependencies and Interactions

*   **`cp_control_types`**: Provides `admm_control_type`, which contains the input parameters that configure the `admm_dm_type` instance.
*   **`cp_dbcsr_api`**: Provides `dbcsr_type` for sparse matrix representation and `dbcsr_release` for deallocating these matrices. `admm_dm_type` stores `matrix_a` as a `dbcsr_type`. `mcweeny_history_type` also stores a `dbcsr_type` matrix.
*   **`input_constants`**: Supplies integer constants (e.g., `do_admm_basis_projection`, `do_admm_purify_mcweeny`) that are used to set the `method` and `purify` fields of `admm_dm_type`.
*   **`kinds`**: Provides the `dp` kind parameter for double-precision real numbers (e.g., `eps_filter`).
*   **`admm_dm_methods`**: This module (documented separately) contains the subroutines that consume and operate on `admm_dm_type` objects (e.g., `admm_dm_calc_rho_aux`, `admm_dm_merge_ks_matrix`). The `revert_purify_mcweeny` subroutine in `admm_dm_methods` is responsible for deallocating the individual `mcweeny_history_type` entries.

This module is foundational for ADMM operations focusing on density matrices, providing the necessary data structures to manage state and parameters.
```
