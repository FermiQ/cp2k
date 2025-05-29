# ace_nlist Module Documentation

## Overview

The `ace_nlist` module provides an interface to compute Atomic Cluster Expansion (ACE) potentials. It is designed to work within the CP2K simulation package. The primary component of this module is the `ace_interface` subroutine, which orchestrates the calculation of ACE potential, forces, and virial based on atomic configurations and neighbor lists.

## Key Components

### `ace_interface` Subroutine

This is the main public subroutine within the `ace_nlist` module. It serves as the bridge between CP2K's main simulation engine and the underlying ACE library/potential calculation.

**Purpose:**
To calculate the ACE potential energy, atomic forces, and virial for a given atomic configuration.

**Inputs:**
*   `ace_natom` (Integer): The number of atoms for which the ACE potential is to be calculated.
*   `ace_atype` (Integer Array): An array of size `ace_natom` containing the types of each atom.
*   `fist_nonbond_env` (Type `fist_nonbond_env_type`, Pointer): Pointer to the non-bonded environment, which contains neighbor list information and other relevant data.
*   `cell` (Type `cell_type`, Pointer): Pointer to the simulation cell information (e.g., cell vectors).
*   `ace_data` (Type `ace_data_type`, Pointer): Pointer to a structure holding ACE-specific data, such as pre-calculated coefficients or model parameters.

**Outputs:**
*   `pot_ace` (Real, kind=dp): The calculated total ACE potential energy.
*   `ace_force` (Real Array, kind=dp): A 3x`ace_natom` array containing the x, y, and z components of the force on each atom.
*   `ace_virial` (Real Array, kind=dp): A 6-element array representing the virial tensor components (xx, yy, zz, xy, xz, yz).

**Functionality:**
1.  Retrieves neighbor list information from `fist_nonbond_env`.
2.  Updates internal data structures (`ace_data`) if the neighbor lists have changed since the last call. This involves:
    *   Allocating/reallocating arrays for neighbor indices, atom positions, types, and ghost atom information.
    *   Building the neighbor list (`nlist`) for the ACE calculation, including handling periodic boundary conditions and ghost atoms.
    *   Transforming atom indices and positions to be compatible with the ACE library (e.g., 0-based indexing).
3.  If neighbor lists have not changed, it updates only the atomic positions.
4.  Calls the `ace_model_compute` subroutine (from `ace_wrapper`) to perform the actual ACE calculation using the prepared data.
5.  Reshapes the output forces and sums the per-atom energies to get the total potential.
6.  Includes a compile-time check (`#if defined(__ACE)`) to ensure CP2K was compiled with ACE library support, aborting if not.

## Important Variables/Constants

*   `moduleN` (Character, Private, Parameter): Name of the module, 'ace_nlist'.
*   Within `ace_interface`:
    *   `natom` (Integer): Local copy of `ace_natom`.
    *   `nghost` (Integer): Number of ghost atoms.
    *   `energy` (Real Array, kind=dp): Per-atom energy, allocated and used during the call to `ace_model_compute`.
    *   `forceunroll` (Real Array, kind=dp): Unrolled (1D) array for forces before reshaping.
    *   `nonbonded` (Type `fist_neighbor_type`, Pointer): Pointer to FIST neighbor list data.
    *   `r_last_update_pbc` (Type `pos_type`, Array, Pointer): Pointer to atom positions at the last neighbor list update.

## Usage Examples

Direct usage of `ace_interface` typically occurs within the FIST module of CP2K when an ACE potential is selected in the input file. A user would not call this subroutine directly in their CP2K input.

Example conceptual call (not actual CP2K input):
```fortran
CALL ace_interface( &
    ace_natom=my_natoms, &
    ace_atype=my_atom_types, &
    pot_ace=my_ace_potential, &
    ace_force=my_ace_forces, &
    ace_virial=my_ace_virial, &
    fist_nonbond_env=my_fist_env, &
    cell=my_cell_info, &
    ace_data=my_ace_model_data &
)
```

## Dependencies and Interactions

The `ace_nlist` module relies on several other CP2K modules:

*   **`ace_wrapper`**:
    *   Uses `ace_model_compute` to perform the core ACE calculation.
*   **`cell_types`**:
    *   Uses `cell_type` for representing the simulation cell.
*   **`fist_neighbor_list_types`**:
    *   Uses `fist_neighbor_type` for handling neighbor lists.
    *   Uses `neighbor_kind_pairs_type` for information about pairs of atom kinds in neighbor lists.
*   **`fist_nonbond_env_types`**:
    *   Uses `ace_data_type` to store and manage ACE-specific data.
    *   Uses `fist_nonbond_env_get` to retrieve data from the non-bonded environment.
    *   Uses `fist_nonbond_env_type` for the non-bonded environment structure.
    *   Uses `pos_type` for storing atomic positions.
*   **`kinds`**:
    *   Uses `dp` for defining double-precision real numbers.
*   **`physcon`**:
    *   Uses `angstrom` for unit conversion (positions are converted to Angstroms for the ACE calculation).

**External Libraries:**
*   The module's functionality is critically dependent on an external ACE library. The preprocessor directive `#if defined(__ACE)` guards the core logic, and CP2K will abort if compiled without ACE support when this interface is invoked.

**Interaction Flow:**
1.  Higher-level routines in CP2K (likely within the FIST module) prepare the necessary data structures (`fist_nonbond_env`, `cell`, `ace_data`).
2.  `ace_interface` is called with this data.
3.  `ace_interface` processes the neighbor lists, potentially rebuilding them and managing ghost atoms.
4.  It calls `ace_model_compute` from the `ace_wrapper`, passing the formatted atomic data and neighbor lists.
5.  `ace_model_compute` (presumably) calls the external ACE library.
6.  Results (energy, forces, virial) are returned up the call stack.
```
