# almo_scf_methods Module Documentation

## Overview

The `almo_scf_methods` module in CP2K contains a collection of subroutines that provide essential functionalities for Absolutely Localized Molecular Orbitals Self-Consistent Field (ALMO-SCF) calculations. These methods are primarily concerned with matrix manipulations, transformations between different orbital representations, domain-specific operations, numerical linear algebra tasks like matrix inversion and square roots, and utility functions supporting the ALMO SCF workflow. They are used by the higher-level driver routines in `almo_scf.F` and `almo_scf_optimizer.F`.

## Key Components

This module provides a variety of public subroutines:

### Matrix Transformations and Projections

*   **`almo_scf_ks_to_ks_xx(almo_scf_env)`**:
    *   Computes projected Kohn-Sham (KS) matrices for overlapping domains (`almo_scf_env%domain_ks_xx`). This involves complex matrix products like `KS*T*Sigma_inv*T^T*S` and its variants.
    *   Also calculates the DIIS error vector (`almo_scf_env%matrix_err_xx` in AO-MO basis, `almo_scf_env%domain_err` in AO-AO domain basis) as a by-product.

*   **`almo_scf_ks_to_ks_blk(almo_scf_env)`**:
    *   Computes the block-diagonal KS matrix (`almo_scf_env%matrix_ks_blk`) from the full KS matrix. This projection ensures that interactions are only within defined ALMO domains.
    *   Also calculates the DIIS error vector for the block-diagonal ALMO SCF (`almo_scf_env%matrix_err_blk`).

*   **`almo_scf_p_blk_to_t_blk(almo_scf_env, ionic)`**:
    *   Converts a block-diagonal density matrix (`almo_scf_env%matrix_p_blk`) into block-diagonal ALMO coefficients (`almo_scf_env%matrix_t_blk`).
    *   If `ionic` is true, it diagonalizes each block of `matrix_p_blk`.
    *   Otherwise (non-ionic guess), it projects random orbitals onto `matrix_p_blk`.

*   **`almo_scf_t_to_proj(...)`**:
    *   Computes a projector `P` from ALMO coefficients `T`.
    *   If orbitals `T` are orthogonal, `P = T * T^T`.
    *   If non-orthogonal, `P = T * Sigma_inv * T^T * S`, where `Sigma = T^T * S * T` is the ALMO overlap matrix and `S` is the AO overlap matrix. This routine also computes `Sigma_inv`.
    *   Supports various algorithms for `Sigma_inv` (Hotelling, Taylor, Cholesky).

*   **`apply_projector(...)`**:
    *   Applies a projector `P = |psi_proj> Sigma_inv_proj <psi_proj|S` to a set of input orbitals `psi_in`.
    *   Can compute either `psi_out = P * psi_in` or the complementary projection `psi_out = (1-P) * psi_in`.
    *   Handles both orthogonal and non-orthogonal `psi_projector`.

### Orbital Diagonalization and Orthogonalization

*   **`almo_scf_ks_xx_to_tv_xx(almo_scf_env)`**:
    *   Computes ALMO coefficients (`almo_scf_env%matrix_t`) by diagonalizing the domain-projected KS matrices (`almo_scf_env%domain_ks_xx`).
    *   Transforms the resulting eigenvectors from the orthogonalized domain basis back to the AO basis.
    *   If smearing is enabled, stores occupied MO energies.

*   **`almo_scf_ks_blk_to_tv_blk(almo_scf_env)`**:
    *   Computes block-diagonal ALMO coefficients (`almo_scf_env%matrix_t_blk`) and virtual orbitals (`almo_scf_env%matrix_v_full_blk`).
    *   This is done by diagonalizing the blocks of the orthogonalized block-diagonal KS matrix.
    *   Also stores block-diagonal MO energies (`almo_scf_env%matrix_eoo`, `almo_scf_env%matrix_evv_full`).

*   **`orthogonalize_mos(...)`**:
    *   Orthogonalizes a given set of MOs (`ket`) using their overlap matrix `Sigma = ket^dagger * metric * ket`.
    *   The transformation is `ket_new = ket * Sigma^(-1/2)`.
    *   `Sigma^(-1/2)` is computed using `matrix_sqrt_Newton_Schulz`.
    *   Can retain locality if `retain_locality` is true (by working with a block-diagonal `Sigma`).

### Domain-Specific Operations

*   **`distribute_domains(almo_scf_env)`**:
    *   Performs load balancing by assigning domains (fragments) to different CPUs (MPI ranks).
    *   The assignment is based on an estimated computational cost per domain (proportional to `nao_domain^3`).

*   **`apply_domain_operators(...)`**:
    *   Applies domain-specific operators (e.g., preconditioner, projector) to a distributed matrix (`matrix_in`) to produce `matrix_out`.
    *   It first extracts submatrices corresponding to domains, applies local operations using `operator1` (and optionally `operator2`), and then reassembles `matrix_out`.

*   **`construct_domain_preconditioner(...)`**:
    *   Constructs domain-specific preconditioners by inverting blocks of an input matrix (`matrix_main`).
    *   Can optionally project out "bad modes" (eigenvectors with small eigenvalues) during inversion.

*   **`construct_domain_s_sqrt(...)`**:
    *   Computes domain-specific S<sup>+1/2</sup> and S<sup>-1/2</sup> matrices from the AO overlap matrix `S`. Each domain's block of `S` is processed independently.

*   **`construct_domain_s_inv(...)`**:
    *   Computes domain-specific S<sup>-1</sup> matrices from the AO overlap matrix `S`.

*   **`construct_domain_r_down(...)`**:
    *   Constructs domain sub-blocks of the covariant density matrix `R = S * T * Sigma_inv * T^dagger * S`.

### Numerical Linear Algebra Utilities

*   **`pseudo_invert_diagonal_blk(...)`**: Inverts only the diagonal blocks of a DBCSR matrix.
*   **`matrix_sqrt(A, Asqrt, Asqrtinv, N)`**: Computes the matrix square root `Asqrt` and inverse square root `Asqrtinv` of a dense symmetric matrix `A` via eigenvalue decomposition.
*   **`pseudo_invert_matrix(...)`**: Inverts a dense symmetric matrix `A`. Allows for pseudo-inversion by projecting out eigenmodes based on eigenvalue ranges or thresholds.
*   **`pseudo_matrix_power(A, Apow, power, N, ...)`**: Computes `A^power` for a dense symmetric matrix `A` via eigenvalue decomposition, with optional mode projection.
*   **`generator_to_unitary(X, U, eps_filter)`**: Computes a unitary matrix `U` from a generator matrix `X` using the Cayley transform: `U = (I - Delta) * (I + Delta)^(-1)`, where `Delta = X - X^dagger` is skew-symmetric.

### Other Utilities

*   **`fill_matrix_with_ones(matrix)`**: Sets all (allocated) blocks of a DBCSR matrix to 1.0.
*   **`almo_scf_t_rescaling(...)`**: For Fermi-Dirac smearing, rescales ALMO coefficients `T` by `sqrt(occupation_number)`. It calculates the chemical potential (`mu_of_domain`) for each domain to achieve the target number of electrons (`real_ne_of_domain`) and computes the electronic entropy contribution (`spin_kTS`).
*   **`get_overlap(...)`**: Computes the MO overlap matrix `Sigma = bra^T * metric * ket`. Corrects the diagonal for smearing if active.
*   **`set_zero_electron_blocks_in_mo_mo_matrix(...)`**: Utility to set diagonal blocks of an MO-MO matrix to a specific value if the corresponding domain has no occupied orbitals.
*   **`construct_test(...)`**: A debugging routine for testing construction and printing of domain submatrices.
*   **`xalmo_initial_guess(...)`**: Generates an initial guess for XALMO (eXtended ALMO) calculations. Can use ASPC extrapolation based on `xalmo_history`. If not extrapolating, it can start from zero delocalization amplitudes or a provided input guess.

## Important Variables/Constants

*   `almo_scf_env`: The ALMO SCF environment, frequently passed around, containing matrices and settings.
*   `matrix_t`, `matrix_t_blk`: ALMO coefficients (delocalized and block-diagonal).
*   `matrix_s`, `matrix_s_blk`: AO overlap matrix (full and block-diagonal ALMO representation).
*   `matrix_ks`, `matrix_ks_blk`: KS matrix (full and block-diagonal ALMO representation).
*   `matrix_sigma`, `matrix_sigma_inv`: ALMO overlap matrix and its inverse.
*   `domain_ks_xx`: KS matrix projected onto overlapping domains.
*   `quench_t`: Quenching matrix defining allowed interactions in XALMO.
*   `nocc_of_domain`: Array storing number of occupied orbitals per domain.

## Usage Examples

These subroutines are internal components of the ALMO SCF machinery and are called by routines in `almo_scf.F` or `almo_scf_optimizer.F`. For instance:
*   `almo_scf_ks_blk_to_tv_blk` is used within the block-diagonal ALMO SCF loop to get new orbitals from the diagonalized KS blocks.
*   `orthogonalize_mos` is used to ensure ALMOs satisfy orthogonality constraints.
*   `almo_scf_t_to_proj` is used to compute density matrices/projectors from ALMO coefficients.
*   `distribute_domains` is called during initialization for load balancing.

## Dependencies and Interactions

*   **`almo_scf_types`**: Provides `almo_scf_env_type` and `almo_scf_history_type`.
*   **`cp_dbcsr_api`, `cp_dbcsr_cholesky`, `cp_dbcsr_contrib`**: Essential for all sparse matrix (DBCSR) operations.
*   **`domain_submatrix_types`, `domain_submatrix_methods`**: For handling operations on matrices decomposed into domain-specific blocks.
*   **`fermi_utils`**: For Fermi-Dirac smearing calculations (`FermiFixed`).
*   **`input_constants`**: Provides named constants for various algorithms and options.
*   **`iterate_matrix`**: For iterative matrix functions like `matrix_sqrt_Newton_Schulz`, `invert_Hotelling`, `invert_Taylor`.
*   **`mathlib`**: Provides `invmat_symm` (dense matrix inversion) and `binomial` coefficients.
*   **LAPACK/BLAS (via `dsyev`, `dsymm`, `dgemm`):** Used for dense matrix diagonalizations and multiplications within routines like `matrix_sqrt`, `pseudo_invert_matrix`, `pseudo_matrix_power`.
*   **`almo_scf`**: The main ALMO driver module that orchestrates calls to methods in this module.
*   **`almo_scf_optimizer`**: Optimization routines may use utilities from here for specific matrix tasks or transformations.

This module provides a rich toolkit of matrix operations and transformations tailored for the specific needs of ALMO-based SCF methods, including handling of domain locality and sparse distributed matrices.
```
