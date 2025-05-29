# ace_wrapper Module Documentation

## Overview

The `ace_wrapper` module acts as a low-level Fortran interface to a C library implementing Atomic Cluster Expansion (ACE) potentials. It provides Fortran routines that wrap C functions responsible for initializing, computing, and finalizing ACE models. This module is essential for integrating ACE potentials into the CP2K simulation package.

The core functionality is conditionally compiled based on the `__ACE` preprocessor macro. If CP2K is not compiled with ACE support, calls to these routines will result in program termination.

## Key Components

### Derived Type: `ace_model_type`

This type is used to manage the ACE model from the Fortran side.

*   `c_ptr` (Type `C_PTR`): A C pointer that holds the memory address of the ACE model object created and managed by the C library. Initialized to `C_NULL_PTR`.
*   `symbolc` (Character(LEN=2), Allocatable Array): An array storing the chemical symbols (e.g., "H", "Si") of the atom types involved in the ACE potential.

### Fortran Subroutines (Wrappers for C functions)

1.  **`ace_model_initialize(ntypec, symbolc, fname, rcutc, model)`**
    *   **Purpose:** Initializes an ACE potential model using parameters provided and a potential file.
    *   **Arguments:**
        *   `ntypec` (Integer, Input): Number of atom types.
        *   `symbolc` (Character(KIND=C_CHAR, LEN=2) Array, Input): Array of chemical symbols for each atom type.
        *   `fname` (Character(KIND=C_CHAR, LEN=*), Input): Name/path of the ACE potential file.
        *   `rcutc` (Real(KIND=8) Array, Output): A 2D array (ntypec x ntypec) to be filled with cutoff radii by the C library.
        *   `model` (Type `ace_model_type`, Output): The Fortran structure that will hold the pointer to the initialized C ACE model and the symbols.
    *   **Functionality:**
        *   Converts the Fortran `fname` string to a C-compatible null-terminated string.
        *   Calls the C function `AcePotInitialize` with the provided parameters. The C function is expected to load the potential and store a pointer to it in `model%c_ptr`.
        *   Allocates and stores `symbolc` into `model%symbolc`.

2.  **`ace_model_compute(natomc, nghostc, neic, neiatc, originc, nlistc, attypec, atposc, forcec, virialc, energyc, model)`**
    *   **Purpose:** Computes the ACE potential energy, atomic forces, and virial for a given atomic configuration and neighbor list.
    *   **Arguments:**
        *   `natomc` (Integer, Input): Number of real atoms.
        *   `nghostc` (Integer, Input): Number of ghost atoms (for periodic boundary conditions).
        *   `neic` (Integer, Input): Total number of neighbors in the neighbor list.
        *   `neiatc` (Integer Array, Input): Cumulative count of neighbors for each atom (index `0` to `natomc`). `neiatc(i)` is the starting index in `nlistc` for atom `i`.
        *   `originc` (Integer Array, Input): Array mapping ghost atoms to their original real atom indices (1-based for Fortran, but values are 0-based for C).
        *   `nlistc` (Integer Array, Input): Flat list of neighbor indices (0-based for C).
        *   `attypec` (Integer Array, Input): Atom types for all atoms (real + ghost, 0-based for C).
        *   `atposc` (Real(KIND=8) Array, Input): Atomic positions (x,y,z consecutively) for all atoms (real + ghost).
        *   `forcec` (Real(KIND=8) Array, Input/Output): Array to store computed forces on real atoms.
        *   `virialc` (Real(KIND=8) Array, Input/Output): Array to store the computed virial tensor (6 components).
        *   `energyc` (Real(KIND=8) Array, Input/Output): Array to store per-atom energies for real atoms.
        *   `model` (Type `ace_model_type`, Input): The initialized ACE model.
    *   **Functionality:**
        *   Calls the C function `AcePotCompute`, passing all the atomic and model data. The C function performs the actual calculation.

3.  **`ace_model_release(model)`**
    *   **Purpose:** Releases the memory and resources associated with an ACE model.
    *   **Arguments:**
        *   `model` (Type `ace_model_type`, Input/Output): The ACE model to be released.
    *   **Functionality:**
        *   Calls the C function `AcePotFinalize` with the `model%c_ptr`.
        *   Resets `model%c_ptr` to `C_NULL_PTR`.
        *   Deallocates `model%symbolc` if it was allocated.

### C Interface Block (`ace_interface`)

This block defines the signatures of the C functions that the Fortran wrapper subroutines call. This ensures type compatibility between Fortran and C.

*   **`AcePotInitialize(ntypec, symbolc, nlen, cname, rcutc, acedata_ptr)`**
    *   Corresponds to `ace_model_initialize`.
    *   `acedata_ptr` (Type `C_PTR`): Output, pointer to the created ACE model in C.
*   **`AcePotCompute(natomc, nghostc, neic, neiatc, originc, nlistc, attypec, atposc, forcec, virialc, energyc, acedata_ptr)`**
    *   Corresponds to `ace_model_compute`.
    *   `acedata_ptr` (Type `C_PTR`): Input, pointer to the ACE model in C.
*   **`AcePotFinalize(acedata_ptr)`**
    *   Corresponds to `ace_model_release`.
    *   `acedata_ptr` (Type `C_PTR`): Input, pointer to the ACE model in C to be finalized.

## Important Variables/Constants

*   The module uses `ISO_C_BINDING` to ensure interoperability with C, defining types like `C_PTR`, `C_INT`, `C_DOUBLE`, `C_CHAR`, `C_NULL_CHAR`, `C_NULL_PTR`.
*   Preprocessor Macro: `__ACE` determines if the ACE functionality is compiled. If not defined, the subroutines will call `CPABORT`.

## Usage Examples

This module is not intended for direct use in CP2K input files. It is a lower-level component used by other modules (like `ace_nlist` or modules within the FIST framework) that manage the setup and execution of ACE potential calculations.

Conceptual workflow:
1.  Another module calls `ace_model_initialize` to load an ACE potential file and get an `ace_model_type` instance.
2.  During simulation steps, that module prepares atomic data (positions, types, neighbor lists) in the format required by `ace_model_compute`.
3.  It then calls `ace_model_compute` to get energies, forces, and virial.
4.  At the end of the simulation or when the potential is no longer needed, `ace_model_release` is called to clean up.

## Dependencies and Interactions

*   **ISO_C_BINDING**: This standard Fortran module is crucial for defining C-compatible types and enabling the `BIND(C)` attribute for linking with C functions.
*   **`base_uses.f90`**: Included for common base definitions within CP2K, likely providing `MARK_USED` and `CPABORT`.
*   **External C Library (ACE):** The actual ACE potential calculations are performed by an external C library. The `ace_wrapper` module links to this library's functions (`AcePotInitialize`, `AcePotCompute`, `AcePotFinalize`).
*   **Higher-level CP2K Modules:** Modules like `ace_nlist.F` will use `ace_wrapper` to perform ACE calculations. For example, `ace_nlist`'s `ace_interface` calls `ace_model_compute`.

**Compilation:**
*   CP2K must be compiled with the `__ACE` preprocessor macro defined and linked against the appropriate ACE C library for this module to be functional.
```
