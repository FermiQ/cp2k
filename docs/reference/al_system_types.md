# al_system_types Module Documentation

## Overview

The `al_system_types` module in CP2K defines the Fortran derived types (data structures) necessary for implementing the Adaptive Langevin (AD_LANGEVIN) thermostat. This thermostat is a method for canonical (NVT) ensemble sampling that utilizes velocity rescaling techniques, combining aspects of Nosé-Hoover chains and Langevin dynamics. The module provides types to store thermostat variables and overall system parameters, along with subroutines for their initialization and deallocation.

The implementation references the work by Jones et al., J. Chem. Phys. 135, 084124 (2011).

## Key Components

### Derived Types

1.  **`al_thermo_type`**:
    *   **Purpose:** Represents the state and parameters of an individual thermostat, which might be associated with a specific region or set of degrees of freedom.
    *   **Fields:**
        *   `degrees_of_freedom` (INTEGER): The number of degrees of freedom controlled by this thermostat instance. Default: 0.
        *   `nkt` (REAL(KIND=dp)): The target kinetic energy for the thermostatted region (NkT, where N is `degrees_of_freedom`, k is the Boltzmann constant, and T is the target temperature). Default: 0.0_dp.
        *   `chi` (REAL(KIND=dp)): The thermostat variable, often representing a friction or momentum-like term in the extended Lagrangian. Default: 0.0_dp.
        *   `mass` (REAL(KIND=dp)): A fictitious mass associated with the thermostat variable `chi` in its equation of motion. Default: 0.0_dp.
        *   `region_kin_energy` (REAL(KIND=dp)): The actual kinetic energy of the degrees of freedom coupled to this thermostat. Default: 0.0_dp.

2.  **`al_system_type`**:
    *   **Purpose:** Encapsulates the overall Adaptive Langevin thermostat system, potentially managing multiple `al_thermo_type` instances for different regions or a global thermostat.
    *   **Fields:**
        *   `region` (INTEGER): An identifier for the region this thermostat system applies to (though individual `al_thermo_type` instances might handle sub-regions or be part of a chain). Default: 0.
        *   `glob_num_al` (INTEGER): The global number of Adaptive Langevin thermostats across all MPI tasks. Default: 0.
        *   `loc_num_al` (INTEGER): The number of Adaptive Langevin thermostats local to the current MPI task. Default: 0.
        *   `tau_nh` (REAL(KIND=dp)): The time constant associated with the Nosé-Hoover component of the thermostat. Default: 0.0_dp.
        *   `tau_langevin` (REAL(KIND=dp)): The time constant associated with the Langevin (stochastic) component of the thermostat. Default: 0.0_dp.
        *   `dt_fact` (REAL(KIND=dp)): A factor applied to the simulation timestep (`dt`) for thermostat integration, if applicable. Default: 0.0_dp.
        *   `dt` (REAL(KIND=dp)): The simulation timestep. Default: 0.0_dp.
        *   `nvt` (TYPE(`al_thermo_type`), POINTER, DIMENSION(:)): A pointer to a dynamically allocated array of `al_thermo_type` objects. Each element corresponds to a local thermostat instance. Default: `NULL()`.
        *   `map_info` (TYPE(`map_info_type`), POINTER): Pointer to a `map_info_type` (defined in `extended_system_types`). This is likely used for mapping specific atoms or degrees of freedom to their corresponding thermostat instance(s), especially if regional thermostats are used. Default: `NULL()`.

### Public Subroutines

1.  **`al_init(al, simpar, section)`**:
    *   **Purpose:** Initializes an `al_system_type` object.
    *   **Arguments:**
        *   `al` (TYPE(`al_system_type`), POINTER, Output): The Adaptive Langevin system to be initialized.
        *   `simpar` (TYPE(`simpar_type`), POINTER, Input): Global simulation parameters, primarily to obtain the timestep `dt`.
        *   `section` (TYPE(`section_vals_type`), POINTER, Input): CP2K input section values from which thermostat parameters like `TIMECON_NH` (for `tau_nh`) and `TIMECON_LANGEVIN` (for `tau_langevin`) are read.
    *   **Functionality:** Sets default values, assigns `dt` from `simpar`, reads time constants from the input section, cites the relevant scientific paper, and creates the `map_info` object.

2.  **`al_thermo_create(al)`**:
    *   **Purpose:** Allocates and performs basic initialization of the `nvt` array (of `al_thermo_type`) within an already initialized `al_system_type`.
    *   **Arguments:**
        *   `al` (TYPE(`al_system_type`), POINTER, Input/Output): The Adaptive Langevin system whose `nvt` array is to be created.
    *   **Functionality:** Allocates `al%nvt` to hold `al%loc_num_al` thermostat instances and initializes their `chi` values to 0.0. It also allocates a local `seed` array, presumably for random number generation in the Langevin dynamics, although this seed is not stored back into the `al_system_type`.

3.  **`al_dealloc(al)`**:
    *   **Purpose:** Deallocates an `al_system_type` object and all its dynamically allocated components.
    *   **Arguments:**
        *   `al` (TYPE(`al_system_type`), POINTER, Input/Output): The system to be deallocated.
    *   **Functionality:** Calls `al_thermo_dealloc` to free the `nvt` array, `release_map_info_type` to free `map_info`, and then deallocates the `al_system_type` instance itself.

### Private Subroutines (Used Internally)

*   **`al_thermo_dealloc(nvt)`**:
    *   **Purpose:** Deallocates an array of `al_thermo_type`.
    *   **Arguments:**
        *   `nvt` (TYPE(`al_thermo_type`), DIMENSION(:), POINTER, Input/Output): The array to be deallocated.

## Important Variables/Constants

*   `moduleN` (Character, Private, Parameter): Name of the module, 'al_system_types'.
*   `tau_nh`: Time constant for the Nosé-Hoover like deterministic part of the thermostat.
*   `tau_langevin`: Time constant for the stochastic Langevin part of the thermostat.

## Usage Examples

These types and subroutines are used internally by CP2K's molecular dynamics engine when the "AD_LANGEVIN" thermostat is selected in the input file.

```fortran
TYPE(al_system_type), POINTER :: ad_langevin_thermostat
TYPE(simpar_type), POINTER :: simulation_parameters
TYPE(section_vals_type), POINTER :: motion_section

! Initialization
CALL al_init(ad_langevin_thermostat, simulation_parameters, motion_section)
! ... (set al_system_type%loc_num_al and al_system_type%glob_num_al based on system setup) ...
CALL al_thermo_create(ad_langevin_thermostat)

! ... (During MD steps, other routines would use ad_langevin_thermostat to:
!      - Get thermostat parameters like tau_nh, tau_langevin.
!      - Access and update al_thermo_type variables like chi, region_kin_energy.
!      - Use map_info to apply thermostat forces/rescaling to correct atoms/DOFs.) ...

! Deallocation at the end of the simulation
CALL al_dealloc(ad_langevin_thermostat)
```

## Dependencies and Interactions

*   **`bibliography`**: Used to cite the scientific paper by Jones et al. (2011) which likely describes the Adaptive Langevin method.
*   **`extended_system_types`**: Provides `map_info_type` (and its `create_map_info_type`, `release_map_info_type` routines), which is used within `al_system_type` for mapping degrees of freedom to thermostats.
*   **`input_section_types`**: Used by `al_init` to read thermostat-specific parameters (e.g., `TIMECON_NH`, `TIMECON_LANGEVIN`) from the CP2K input file.
*   **`kinds`**: Provides the `dp` kind parameter for double-precision real numbers.
*   **`simpar_types`**: Provides `simpar_type`, from which global simulation parameters like the timestep `dt` are obtained.
*   **MD Modules**: Other modules responsible for molecular dynamics (e.g., integrators, force calculation) would interact with `al_system_type` to implement the thermostat algorithm, applying velocity modifications and evolving the thermostat variables (`chi`).

This module lays the data foundation for the Adaptive Langevin thermostat in CP2K.
```
