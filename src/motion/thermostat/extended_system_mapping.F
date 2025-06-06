!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      CJM, 20-Feb-01
!>      JGH (10-Mar-2001)
!>      CJM (10-Apr-2001)
!> \author CJM
! **************************************************************************************************
MODULE extended_system_mapping

   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                              lnhc_parameters_type,&
                                              map_info_type
   USE input_constants,                 ONLY: &
        do_thermo_communication, do_thermo_no_communication, do_thermo_only_master, &
        isokin_ensemble, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
        nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, &
        npt_ia_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble, reftraj_ensemble
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: global_constraint_type,&
                                              molecule_type
   USE simpar_types,                    ONLY: simpar_type
   USE thermostat_mapping,              ONLY: adiabatic_mapping_region,&
                                              init_baro_map_info,&
                                              thermostat_mapping_region
   USE thermostat_types,                ONLY: thermostat_info_type
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'extended_system_mapping'

   PUBLIC :: nhc_to_particle_mapping, nhc_to_barostat_mapping, &
             nhc_to_shell_mapping, nhc_to_particle_mapping_fast, &
             nhc_to_particle_mapping_slow

CONTAINS

! **************************************************************************************************
!> \brief Creates the thermostatting for the barostat
!> \param simpar ...
!> \param nhc ...
!> \par History
!>      CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!>      JGH (10-Mar-2001) : set nhc variables to zero when not in use
!> \author CJM
! **************************************************************************************************
   SUBROUTINE nhc_to_barostat_mapping(simpar, nhc)

      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(lnhc_parameters_type), POINTER                :: nhc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_barostat_mapping'

      INTEGER                                            :: handle, i, number
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      SELECT CASE (simpar%ensemble)
      CASE DEFAULT
         CPABORT('Never reach this point!')
      CASE (npt_i_ensemble, npt_f_ensemble, npt_ia_ensemble)
         map_info => nhc%map_info
         map_info%dis_type = do_thermo_only_master

         ! Counting the total number of thermostats ( 1 for NPT_I, NPT_IA, and NPT_F )
         nhc%loc_num_nhc = 1
         nhc%glob_num_nhc = 1
         IF (simpar%ensemble == npt_f_ensemble) THEN
            number = 9
         ELSE
            number = 1
         END IF

         CALL init_baro_map_info(map_info, number, nhc%loc_num_nhc)

         ALLOCATE (nhc%nvt(nhc%nhc_len, nhc%loc_num_nhc))
         ! Now that we know how many there are stick this into nhc % nkt
         ! (number of degrees of freedom times k_B T )
         DO i = 1, nhc%loc_num_nhc
            nhc%nvt(1, i)%nkt = simpar%temp_ext*number
            nhc%nvt(1, i)%degrees_of_freedom = number
            IF (debug_isotropic_limit) THEN
               nhc%nvt(1, i)%nkt = simpar%temp_ext
            END IF
         END DO

         ! getting the number of degrees of freedom times k_B T for the rest of the chain
         DO i = 2, nhc%nhc_len
            nhc%nvt(i, :)%nkt = simpar%temp_ext
         END DO

         ! Let's clean the arrays
         map_info%s_kin = 0.0_dp
         map_info%v_scale = 0.0_dp
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE nhc_to_barostat_mapping

! **************************************************************************************************
!> \brief Creates the thermostatting maps
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param nhc ...
!> \param para_env ...
!> \param gci ...
!> \par History
!>      29-Nov-00 (JGH) correct counting of DOF if constraints are off
!>      CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!>      JGH (10-Mar-2001) : set nhc variables to zero when not in use
!>      CJM(10-NOV-2001) : New parallelization with new molecule structures
!>      Teodoro Laino 09.2007 [tlaino] - University of Zurich - cleaning and updating
!> \author CJM
! **************************************************************************************************
   SUBROUTINE nhc_to_particle_mapping(thermostat_info, simpar, local_molecules, &
                                      molecule_set, molecule_kind_set, nhc, para_env, gci)

      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_constraint_type), POINTER              :: gci

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_particle_mapping'

      INTEGER                                            :: handle, i, imap, j, natoms_local, &
                                                            sum_of_thermostats
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_atom_list
      REAL(KIND=dp)                                      :: fac
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      NULLIFY (massive_atom_list, deg_of_freedom)

      SELECT CASE (simpar%ensemble)
      CASE DEFAULT
         CPABORT('Unknown ensemble!')
      CASE (nve_ensemble, isokin_ensemble, npe_f_ensemble, npe_i_ensemble, nph_uniaxial_ensemble, &
            nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble)
         CPABORT('Never reach this point!')
      CASE (nvt_ensemble, npt_i_ensemble, npt_f_ensemble, npt_ia_ensemble)

         CALL setup_nhc_thermostat(nhc, thermostat_info, deg_of_freedom, massive_atom_list, &
                                   molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
                                   simpar, sum_of_thermostats, gci)

         ! Sum up the number of degrees of freedom on each thermostat.
         ! first: initialize the target
         map_info => nhc%map_info
         map_info%s_kin = 0.0_dp
         DO i = 1, 3
            DO j = 1, natoms_local
               map_info%p_kin(i, j)%point = map_info%p_kin(i, j)%point + 1
            END DO
         END DO

         ! if thermostats are replicated but molecules distributed, we have to
         ! sum s_kin over all processors
         IF (map_info%dis_type == do_thermo_communication) CALL para_env%sum(map_info%s_kin)

         ! We know the total number of system thermostats.
         IF ((sum_of_thermostats == 1) .AND. (map_info%dis_type /= do_thermo_no_communication)) THEN
            fac = map_info%s_kin(1) - deg_of_freedom(1) - simpar%nfree_rot_transl
            IF (fac == 0.0_dp) THEN
               CPABORT('Zero degrees of freedom. Nothing to thermalize!')
            END IF
            nhc%nvt(1, 1)%nkt = simpar%temp_ext*fac
            nhc%nvt(1, 1)%degrees_of_freedom = FLOOR(fac)
         ELSE
            DO i = 1, nhc%loc_num_nhc
               imap = map_info%map_index(i)
               fac = (map_info%s_kin(imap) - deg_of_freedom(i))
               nhc%nvt(1, i)%nkt = simpar%temp_ext*fac
               nhc%nvt(1, i)%degrees_of_freedom = FLOOR(fac)
            END DO
         END IF

         ! Getting the number of degrees of freedom times k_B T for the rest
         ! of the chain
         DO i = 2, nhc%nhc_len
            nhc%nvt(i, :)%nkt = simpar%temp_ext
            nhc%nvt(i, :)%degrees_of_freedom = 1
         END DO
         DEALLOCATE (deg_of_freedom)
         DEALLOCATE (massive_atom_list)

         ! Let's clean the arrays
         map_info%s_kin = 0.0_dp
         map_info%v_scale = 0.0_dp
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE nhc_to_particle_mapping

! **************************************************************************************************
!> \brief Main general setup for Adiabatic Nose-Hoover thermostats
!> \param nhc ...
!> \param thermostat_info ...
!> \param deg_of_freedom ...
!> \param massive_atom_list ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param para_env ...
!> \param natoms_local ...
!> \param simpar ...
!> \param sum_of_thermostats ...
!> \param gci ...
!> \param shell ...
!> \author CJM -PNNL -2011
! **************************************************************************************************
   SUBROUTINE setup_adiabatic_thermostat(nhc, thermostat_info, deg_of_freedom, &
                                         massive_atom_list, molecule_kind_set, local_molecules, molecule_set, &
                                         para_env, natoms_local, simpar, sum_of_thermostats, gci, shell)

      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_atom_list
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: natoms_local
      TYPE(simpar_type), POINTER                         :: simpar
      INTEGER, INTENT(OUT)                               :: sum_of_thermostats
      TYPE(global_constraint_type), POINTER              :: gci
      LOGICAL, INTENT(IN), OPTIONAL                      :: shell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_adiabatic_thermostat'

      INTEGER                                            :: handle, nkind, number, region
      LOGICAL                                            :: do_shell
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      do_shell = .FALSE.
      IF (PRESENT(shell)) do_shell = shell
      map_info => nhc%map_info

      nkind = SIZE(molecule_kind_set)
      sum_of_thermostats = thermostat_info%sum_of_thermostats
      map_info%dis_type = thermostat_info%dis_type
      number = thermostat_info%number_of_thermostats
      region = nhc%region

      CALL adiabatic_mapping_region(map_info, deg_of_freedom, massive_atom_list, &
                                    molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
                                    simpar, number, region, gci, do_shell, thermostat_info%map_loc_thermo_gen, &
                                    sum_of_thermostats)
      ALLOCATE (nhc%nvt(nhc%nhc_len, number))

      ! Now that we know how many there are stick this into nhc%nkt
      ! (number of degrees of freedom times k_B T for the first thermostat
      !  on the chain)
      nhc%loc_num_nhc = number
      nhc%glob_num_nhc = sum_of_thermostats

      CALL timestop(handle)

   END SUBROUTINE setup_adiabatic_thermostat

! **************************************************************************************************
!> \brief Creates the thermostatting maps
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param nhc ...
!> \param para_env ...
!> \param gci ...
!> \par History
!> \author CJM
! **************************************************************************************************
   SUBROUTINE nhc_to_particle_mapping_slow(thermostat_info, simpar, local_molecules, &
                                           molecule_set, molecule_kind_set, nhc, para_env, gci)

      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_constraint_type), POINTER              :: gci

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_particle_mapping_slow'

      INTEGER                                            :: handle, i, imap, j, natoms_local, &
                                                            sum_of_thermostats
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_atom_list
      REAL(KIND=dp)                                      :: fac
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      NULLIFY (massive_atom_list, deg_of_freedom)

      SELECT CASE (simpar%ensemble)
      CASE DEFAULT
         CPABORT('Unknown ensemble!')
      CASE (nvt_adiabatic_ensemble)
         CALL setup_adiabatic_thermostat(nhc, thermostat_info, deg_of_freedom, massive_atom_list, &
                                         molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
                                         simpar, sum_of_thermostats, gci)

         ! Sum up the number of degrees of freedom on each thermostat.
         ! first: initialize the target
         map_info => nhc%map_info
         map_info%s_kin = 0.0_dp
         DO i = 1, 3
            DO j = 1, natoms_local
               IF (ASSOCIATED(map_info%p_kin(i, j)%point)) &
                  map_info%p_kin(i, j)%point = map_info%p_kin(i, j)%point + 1
            END DO
         END DO

         ! if thermostats are replicated but molecules distributed, we have to
         ! sum s_kin over all processors
         IF (map_info%dis_type == do_thermo_communication) CALL para_env%sum(map_info%s_kin)

         ! We know the total number of system thermostats.
         IF ((sum_of_thermostats == 1) .AND. (map_info%dis_type /= do_thermo_no_communication)) THEN
            fac = map_info%s_kin(1) - deg_of_freedom(1) - simpar%nfree_rot_transl
            IF (fac == 0.0_dp) THEN
               CPABORT('Zero degrees of freedom. Nothing to thermalize!')
            END IF
            nhc%nvt(1, 1)%nkt = simpar%temp_slow*fac
            nhc%nvt(1, 1)%degrees_of_freedom = FLOOR(fac)
         ELSE
            DO i = 1, nhc%loc_num_nhc
               imap = map_info%map_index(i)
               fac = (map_info%s_kin(imap) - deg_of_freedom(i))
               nhc%nvt(1, i)%nkt = simpar%temp_slow*fac
               nhc%nvt(1, i)%degrees_of_freedom = FLOOR(fac)
            END DO
         END IF

         ! Getting the number of degrees of freedom times k_B T for the rest
         ! of the chain
         DO i = 2, nhc%nhc_len
            nhc%nvt(i, :)%nkt = simpar%temp_slow
            nhc%nvt(i, :)%degrees_of_freedom = 1
         END DO
         DEALLOCATE (deg_of_freedom)
         DEALLOCATE (massive_atom_list)

         ! Let's clean the arrays
         map_info%s_kin = 0.0_dp
         map_info%v_scale = 0.0_dp
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE nhc_to_particle_mapping_slow

! **************************************************************************************************
!> \brief Creates the thermostatting maps
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param nhc ...
!> \param para_env ...
!> \param gci ...
!> \par History
!> \author CJM
! **************************************************************************************************
   SUBROUTINE nhc_to_particle_mapping_fast(thermostat_info, simpar, local_molecules, &
                                           molecule_set, molecule_kind_set, nhc, para_env, gci)

      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_constraint_type), POINTER              :: gci

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_particle_mapping_fast'

      INTEGER                                            :: handle, i, imap, j, natoms_local, &
                                                            sum_of_thermostats
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_atom_list
      REAL(KIND=dp)                                      :: fac
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      NULLIFY (massive_atom_list, deg_of_freedom)

      SELECT CASE (simpar%ensemble)
      CASE DEFAULT
         CPABORT('Unknown ensemble!')
      CASE (nvt_adiabatic_ensemble)
         CALL setup_adiabatic_thermostat(nhc, thermostat_info, deg_of_freedom, massive_atom_list, &
                                         molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
                                         simpar, sum_of_thermostats, gci)

         ! Sum up the number of degrees of freedom on each thermostat.
         ! first: initialize the target
         map_info => nhc%map_info
         map_info%s_kin = 0.0_dp
         DO i = 1, 3
            DO j = 1, natoms_local
               IF (ASSOCIATED(map_info%p_kin(i, j)%point)) &
                  map_info%p_kin(i, j)%point = map_info%p_kin(i, j)%point + 1
            END DO
         END DO

         ! if thermostats are replicated but molecules distributed, we have to
         ! sum s_kin over all processors
         IF (map_info%dis_type == do_thermo_communication) CALL para_env%sum(map_info%s_kin)

         ! We know the total number of system thermostats.
         IF ((sum_of_thermostats == 1) .AND. (map_info%dis_type /= do_thermo_no_communication)) THEN
            fac = map_info%s_kin(1) - deg_of_freedom(1) - simpar%nfree_rot_transl
            IF (fac == 0.0_dp) THEN
               CPABORT('Zero degrees of freedom. Nothing to thermalize!')
            END IF
            nhc%nvt(1, 1)%nkt = simpar%temp_fast*fac
            nhc%nvt(1, 1)%degrees_of_freedom = FLOOR(fac)
         ELSE
            DO i = 1, nhc%loc_num_nhc
               imap = map_info%map_index(i)
               fac = (map_info%s_kin(imap) - deg_of_freedom(i))
               nhc%nvt(1, i)%nkt = simpar%temp_fast*fac
               nhc%nvt(1, i)%degrees_of_freedom = FLOOR(fac)
            END DO
         END IF

         ! Getting the number of degrees of freedom times k_B T for the rest
         ! of the chain
         DO i = 2, nhc%nhc_len
            nhc%nvt(i, :)%nkt = simpar%temp_fast
            nhc%nvt(i, :)%degrees_of_freedom = 1
         END DO
         DEALLOCATE (deg_of_freedom)
         DEALLOCATE (massive_atom_list)

         ! Let's clean the arrays
         map_info%s_kin = 0.0_dp
         map_info%v_scale = 0.0_dp
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE nhc_to_particle_mapping_fast

! **************************************************************************************************
!> \brief Main general setup for Nose-Hoover thermostats
!> \param nhc ...
!> \param thermostat_info ...
!> \param deg_of_freedom ...
!> \param massive_atom_list ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param para_env ...
!> \param natoms_local ...
!> \param simpar ...
!> \param sum_of_thermostats ...
!> \param gci ...
!> \param shell ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! **************************************************************************************************
   SUBROUTINE setup_nhc_thermostat(nhc, thermostat_info, deg_of_freedom, &
                                   massive_atom_list, molecule_kind_set, local_molecules, molecule_set, &
                                   para_env, natoms_local, simpar, sum_of_thermostats, gci, shell)

      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_atom_list
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: natoms_local
      TYPE(simpar_type), POINTER                         :: simpar
      INTEGER, INTENT(OUT)                               :: sum_of_thermostats
      TYPE(global_constraint_type), POINTER              :: gci
      LOGICAL, INTENT(IN), OPTIONAL                      :: shell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_nhc_thermostat'

      INTEGER                                            :: handle, nkind, number, region
      LOGICAL                                            :: do_shell
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      do_shell = .FALSE.
      IF (PRESENT(shell)) do_shell = shell
      map_info => nhc%map_info

      nkind = SIZE(molecule_kind_set)
      sum_of_thermostats = thermostat_info%sum_of_thermostats
      map_info%dis_type = thermostat_info%dis_type
      number = thermostat_info%number_of_thermostats
      region = nhc%region

      CALL thermostat_mapping_region(map_info, deg_of_freedom, massive_atom_list, &
                                     molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
                                     simpar, number, region, gci, do_shell, thermostat_info%map_loc_thermo_gen, &
                                     sum_of_thermostats)

      ALLOCATE (nhc%nvt(nhc%nhc_len, number))

      ! Now that we know how many there are stick this into nhc%nkt
      ! (number of degrees of freedom times k_B T for the first thermostat
      !  on the chain)
      nhc%loc_num_nhc = number
      nhc%glob_num_nhc = sum_of_thermostats

      CALL timestop(handle)

   END SUBROUTINE setup_nhc_thermostat

! **************************************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param nhc ...
!> \param para_env ...
!> \param gci ...
! **************************************************************************************************
   SUBROUTINE nhc_to_shell_mapping(thermostat_info, simpar, local_molecules, &
                                   molecule_set, molecule_kind_set, nhc, para_env, gci)

      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(lnhc_parameters_type), POINTER                :: nhc
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_constraint_type), POINTER              :: gci

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_shell_mapping'

      INTEGER                                            :: handle, i, imap, j, nshell_local, &
                                                            sum_of_thermostats
      INTEGER, DIMENSION(:), POINTER                     :: deg_of_freedom, massive_shell_list
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)

      NULLIFY (massive_shell_list, deg_of_freedom)

      SELECT CASE (simpar%ensemble)
      CASE DEFAULT
         CPABORT('Unknown ensemble!')
      CASE (isokin_ensemble, nph_uniaxial_ensemble, &
            nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble)
         CPABORT('Never reach this point!')
      CASE (nve_ensemble, nvt_ensemble, npe_f_ensemble, npe_i_ensemble, npt_i_ensemble, npt_f_ensemble, &
            npt_ia_ensemble)

         CALL setup_nhc_thermostat(nhc, thermostat_info, deg_of_freedom, massive_shell_list, &
                                   molecule_kind_set, local_molecules, molecule_set, para_env, nshell_local, &
                                   simpar, sum_of_thermostats, gci, shell=.TRUE.)

         map_info => nhc%map_info
         ! Sum up the number of degrees of freedom on each thermostat.
         ! first: initialize the target, via p_kin init s_kin
         map_info%s_kin = 0.0_dp
         DO j = 1, nshell_local
            DO i = 1, 3
               map_info%p_kin(i, j)%point = map_info%p_kin(i, j)%point + 1
            END DO
         END DO

         ! If thermostats are replicated but molecules distributed, we have to
         ! sum s_kin over all processors
         IF (map_info%dis_type == do_thermo_communication) CALL para_env%sum(map_info%s_kin)

         ! Now that we know how many there are stick this into nhc%nkt
         ! (number of degrees of freedom times k_B T )
         DO i = 1, nhc%loc_num_nhc
            imap = map_info%map_index(i)
            nhc%nvt(1, i)%nkt = simpar%temp_sh_ext*map_info%s_kin(imap)
            nhc%nvt(1, i)%degrees_of_freedom = INT(map_info%s_kin(imap))
         END DO

         ! Getting the number of degrees of freedom times k_B T for the rest of the chain
         DO i = 2, nhc%nhc_len
            nhc%nvt(i, :)%nkt = simpar%temp_sh_ext
            nhc%nvt(i, :)%degrees_of_freedom = 1
         END DO
         DEALLOCATE (deg_of_freedom)
         DEALLOCATE (massive_shell_list)

         ! Let's clean the arrays
         map_info%s_kin = 0.0_dp
         map_info%v_scale = 0.0_dp
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE nhc_to_shell_mapping

END MODULE extended_system_mapping
