# almo_scf_diis_types Module Documentation

## Overview

The `almo_scf_diis_types` module provides the data structures and routines for a Direct Inversion in the Iterative Subspace (DIIS) extrapolation method, specifically tailored for use within the Absolutely Localized Molecular Orbitals Self-Consistent Field (ALMO-SCF) framework in CP2K. DIIS is a widely used technique to accelerate the convergence of SCF iterations by forming an optimal linear combination of previous iteration vectors (variables) based on their corresponding error vectors.

This module supports DIIS extrapolation for variables represented either as full DBCSR (Distributed Block Compressed Sparse Row) matrices or as collections of domain-specific submatrices, making it adaptable to different aspects of ALMO calculations.

## Key Components

### Parameters

*   **`diis_error_orthogonal` (INTEGER, PARAMETER):** A constant indicating that the error vectors are orthogonal, influencing how their overlap is calculated. Value: 1.
*   **`diis_env_dbcsr` (INTEGER, PARAMETER):** A constant identifying that the DIIS environment is set up to handle full DBCSR matrices. Value: 1.
*   **`diis_env_domain` (INTEGER, PARAMETER):** A constant identifying that the DIIS environment is set up to handle domain submatrices. Value: 2.

### Derived Type: `almo_scf_diis_type`

This is the central data structure for managing the DIIS process.

*   **`diis_env_type` (INTEGER):** Specifies the type of DIIS environment, either `diis_env_dbcsr` or `diis_env_domain`.
*   **`buffer_length` (INTEGER):** The current number of variable/error vector pairs stored in the history.
*   **`max_buffer_length` (INTEGER):** The maximum number of variable/error vector pairs to be stored (i.e., the maximum size of the DIIS subspace).
*   **`m_var` (TYPE(`dbcsr_type`), DIMENSION(:), ALLOCATABLE):** An array storing the history of variable vectors (e.g., MO coefficients or density matrices) as DBCSR matrices. Used when `diis_env_type` is `diis_env_dbcsr`.
*   **`m_err` (TYPE(`dbcsr_type`), DIMENSION(:), ALLOCATABLE):** An array storing the history of error vectors as DBCSR matrices. Used when `diis_env_type` is `diis_env_dbcsr`.
*   **`d_var` (TYPE(`domain_submatrix_type`), DIMENSION(:, :), ALLOCATABLE):** A 2D array storing the history of variable vectors as domain submatrices. The first dimension is the history index, the second is the domain index. Used when `diis_env_type` is `diis_env_domain`.
*   **`d_err` (TYPE(`domain_submatrix_type`), DIMENSION(:, :), ALLOCATABLE):** A 2D array storing the history of error vectors as domain submatrices. Used when `diis_env_type` is `diis_env_domain`.
*   **`m_b` (TYPE(`domain_submatrix_type`), DIMENSION(:), ALLOCATABLE):** An array (one element per domain if `diis_env_type` is `diis_env_domain`, or a single element if `diis_env_dbcsr`) storing the DIIS B-matrix. The B-matrix elements are overlaps of error vectors, `B_ij = Tr(err_i^T * err_j)`. It's constructed to include the Lagrange multiplier constraint for normalization of coefficients.
*   **`in_point` (INTEGER):** The index indicating where the next variable/error pair will be inserted into the history buffer (circular buffer logic).
*   **`error_type` (INTEGER):** Specifies the nature of the error vector (e.g., `diis_error_orthogonal`), which can affect how overlaps are computed.

### Public Subroutines

1.  **`INTERFACE almo_scf_diis_init`**:
    *   A generic interface for initializing the DIIS environment.
    *   **`MODULE PROCEDURE almo_scf_diis_init_dbcsr`**: Initializes DIIS for DBCSR matrices.
    *   **`MODULE PROCEDURE almo_scf_diis_init_domain`**: Initializes DIIS for domain submatrices.

2.  **`almo_scf_diis_init_dbcsr(diis_env, sample_err, sample_var, error_type, max_length)`**:
    *   **Purpose:** Initializes the `diis_env` for use with full DBCSR matrices.
    *   **Arguments:**
        *   `diis_env` (Output): The DIIS environment to initialize.
        *   `sample_err`, `sample_var` (Input): Template DBCSR matrices used for creating storage arrays.
        *   `error_type` (Input): The type of error vector.
        *   `max_length` (Input): Maximum history length.
    *   **Functionality:** Sets `diis_env_type`, buffer parameters, allocates `m_var`, `m_err`, and initializes the B-matrix structure (typically a single 1x1 matrix initially, stored in `m_b(1)%mdata`).

3.  **`almo_scf_diis_init_domain(diis_env, sample_err, error_type, max_length)`**:
    *   **Purpose:** Initializes the `diis_env` for use with domain-decomposed submatrices.
    *   **Arguments:**
        *   `diis_env` (Output): The DIIS environment.
        *   `sample_err` (Input): An array of template domain submatrices.
        *   `error_type` (Input): The type of error vector.
        *   `max_length` (Input): Maximum history length.
    *   **Functionality:** Sets `diis_env_type`, buffer parameters, allocates `d_var`, `d_err`, and initializes the B-matrix structure for each domain (each initially 1x1, stored in `m_b(idomain)%mdata`).

4.  **`almo_scf_diis_push(diis_env, var, err, d_var, d_err)`**:
    *   **Purpose:** Adds a new variable and its corresponding error vector to the DIIS history.
    *   **Arguments:**
        *   `diis_env` (Input/Output): The DIIS environment.
        *   `var`, `err` (Optional Input): DBCSR variable and error matrices.
        *   `d_var`, `d_err` (Optional Input): Domain submatrix arrays for variable and error.
    *   **Functionality:**
        *   Stores the new `var`/`err` (or `d_var`/`d_err`) pair in the history buffer at `in_point`.
        *   Updates `buffer_length` (up to `max_buffer_length`).
        *   Resizes the B-matrix (`m_b(idomain)%mdata`) if the buffer was not yet full.
        *   Computes the new row and column of the B-matrix by calculating overlaps of the new error vector with all previously stored error vectors, using `almo_scf_diis_error_overlap`.
        *   Updates `in_point` for the next insertion.

5.  **`almo_scf_diis_extrapolate(diis_env, extr_var, d_extr_var)`**:
    *   **Purpose:** Computes the DIIS-extrapolated variable.
    *   **Arguments:**
        *   `diis_env` (Input/Output): The DIIS environment.
        *   `extr_var` (Optional Output): The extrapolated DBCSR variable matrix.
        *   `d_extr_var` (Optional Output): The extrapolated domain submatrix variable.
    *   **Functionality:**
        *   For each domain (or globally for `diis_env_dbcsr`):
            1.  Constructs the current B-matrix (including the -1 terms for the Lagrange multiplier constraint).
            2.  Solves the DIIS linear system `B * c = -e_1` (where `e_1 = (1,0,0...)^T`) to find the coefficients `c_i`. This is typically done by diagonalizing B using LAPACK's `dsyev`.
            3.  Computes the extrapolated variable as `extr_var = sum_i (c_i * var_i)`, where `var_i` are the stored historical variables.

6.  **`almo_scf_diis_error_overlap(diis_env, A, B, d_A, d_B) RETURN (overlap)`**:
    *   **Purpose:** A function that calculates the overlap (dot product) between two error vectors.
    *   **Arguments:**
        *   `diis_env` (Input): The DIIS environment.
        *   `A`, `B` (Optional Input/Output): DBCSR error matrices.
        *   `d_A`, `d_B` (Optional Input/Output): Domain submatrix error vectors.
    *   **Return Value:** The calculated overlap (a scalar).
    *   **Functionality:** Currently, if `diis_env%error_type` is `diis_error_orthogonal`, it computes `Tr(A^T * B)` for DBCSR matrices or `sum(d_A .* d_B)` for domain submatrices.

7.  **`almo_scf_diis_release(diis_env)`**:
    *   **Purpose:** Deallocates all dynamically allocated memory within the `diis_env`.
    *   **Functionality:** Releases DBCSR matrices in `m_var` and `m_err`, domain submatrices in `d_var`, `d_err`, and `m_b`, and then deallocates the arrays themselves.

## Important Variables/Constants

*   `diis_env_type`: Determines whether DBCSR matrices or domain submatrices are being handled.
*   `max_buffer_length`: Controls the size of the DIIS subspace and memory usage.
*   `error_type`: Defines how error vector overlaps are computed. `diis_error_orthogonal` is the primary type used.

## Usage Examples

The DIIS routines are typically used within an SCF iterative procedure:

```fortran
TYPE(almo_scf_diis_type) :: diis_optimizer
TYPE(dbcsr_type) :: current_variable, current_error, extrapolated_variable

! Initialization (e.g., for DBCSR matrices)
CALL almo_scf_diis_init_dbcsr(diis_optimizer, template_error_matrix, &
                              template_variable_matrix, diis_error_orthogonal, max_history_len)

! Inside SCF loop:
DO iter = 1, max_iter
  ! ... calculate current_variable (e.g., density matrix or MO coefficients) ...
  ! ... calculate current_error (e.g., gradient or difference vector) ...

  ! Add to DIIS history
  CALL almo_scf_diis_push(diis_optimizer, var=current_variable, err=current_error)

  ! If buffer has enough entries, extrapolate
  IF (diis_optimizer%buffer_length >= min_diis_entries) THEN
    CALL almo_scf_diis_extrapolate(diis_optimizer, extr_var=extrapolated_variable)
    ! Use extrapolated_variable as the input for the next iteration's calculation
    ! current_variable = extrapolated_variable (or a mix)
  ELSE
    ! Use current_variable as input for the next iteration
  END IF

  ! ... check for convergence ...
END DO

! Cleanup
CALL almo_scf_diis_release(diis_optimizer)
```

## Dependencies and Interactions

*   **`cp_dbcsr_api` and `cp_dbcsr_contrib`**: Used for operations on DBCSR matrices (creation, copying, addition, dot product, release).
*   **`domain_submatrix_types` and `domain_submatrix_methods`**: Used when `diis_env_type` is `diis_env_domain`, for handling collections of smaller matrices associated with different domains.
*   **`cp_log_handling`**: For logging messages.
*   **`kinds`**: For specifying double precision (`dp`).
*   **LAPACK (`dsyev`)**: Used internally by `almo_scf_diis_extrapolate` to solve the DIIS linear system via eigenvalue decomposition of the B-matrix.
*   **`almo_scf` module**: The DIIS procedures defined here are called by the main ALMO SCF driver routines in `almo_scf.F` to accelerate convergence of either the block-diagonal ALMO optimization or the XALMO optimization steps.

This module provides a crucial component for efficient SCF convergence within the ALMO framework.
```
