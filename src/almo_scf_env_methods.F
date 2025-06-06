!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief almo_scf_env methods
!> \par History
!>       2016.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_env_methods

   USE almo_scf_types,                  ONLY: almo_max_cutoff_multiplier,&
                                              almo_scf_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: &
        almo_constraint_distance, almo_deloc_none, almo_deloc_xalmo_1diag, &
        almo_domain_layout_atomic, almo_domain_layout_molecular, almo_frz_crystal, &
        almo_mat_distr_molecular, almo_scf_diag, almo_scf_skip, almo_scf_trustr, cg_hager_zhang, &
        do_bondparm_vdw, molecular_guess, tensor_orthogonal, virt_full, virt_minimal, virt_number, &
        xalmo_trial_r0_out
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_env_methods'

   PUBLIC :: almo_scf_env_create

CONTAINS

! **************************************************************************************************
!> \brief Creation and basic initialization of the almo environment
!> \param qs_env ...
!> \par History
!>       2016.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_env_create(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_env_create'

      INTEGER                                            :: handle, nallocate
      TYPE(almo_scf_env_type), POINTER                   :: almo_scf_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ALLOCATE (almo_scf_env)

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      ! parse the almo_scf section and set appropriate quantities
      CALL almo_scf_init_read_write_input(input, almo_scf_env)

      ! set up the buffer for the history of matrices
      almo_scf_env%nspins = dft_control%nspins
      almo_scf_env%almo_history%nstore = almo_scf_env%almo_extrapolation_order
      almo_scf_env%almo_history%istore = 0
      ! do not allocate zero
      nallocate = MAX(1, almo_scf_env%almo_extrapolation_order)
      ALLOCATE (almo_scf_env%almo_history%matrix_p_up_down(almo_scf_env%nspins, nallocate))
      ALLOCATE (almo_scf_env%almo_history%matrix_t(almo_scf_env%nspins))
      almo_scf_env%xalmo_history%nstore = almo_scf_env%xalmo_extrapolation_order
      almo_scf_env%xalmo_history%istore = 0
      nallocate = MAX(1, almo_scf_env%xalmo_extrapolation_order)
      ALLOCATE (almo_scf_env%xalmo_history%matrix_p_up_down(almo_scf_env%nspins, nallocate))
      !ALLOCATE (almo_scf_env%xalmo_history%matrix_x(almo_scf_env%nspins, nallocate))
      ALLOCATE (almo_scf_env%xalmo_history%matrix_t(almo_scf_env%nspins))

      ! put almo_scf_env in qs_env
      CALL set_qs_env(qs_env, almo_scf_env=almo_scf_env)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_env_create

! **************************************************************************************************
!> \brief Parses the ALMO input section
!> \param input ...
!> \param almo_scf_env ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_init_read_write_input(input, almo_scf_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init_read_write_input'

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER :: almo_analysis_section, almo_opt_diis_section, &
         almo_opt_pcg_section, almo_opt_trustr_section, almo_scf_section, matrix_iterate_section, &
         nlmo_opt_pcg_section, penalty_section, xalmo_opt_newton_pcg_section, &
         xalmo_opt_pcg_section, xalmo_opt_trustr_section

      CALL timeset(routineN, handle)

      almo_scf_section => section_vals_get_subs_vals(input, "DFT%ALMO_SCF")
      almo_opt_diis_section => section_vals_get_subs_vals(almo_scf_section, &
                                                          "ALMO_OPTIMIZER_DIIS")
      almo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
                                                         "ALMO_OPTIMIZER_PCG")
      almo_opt_trustr_section => section_vals_get_subs_vals(almo_scf_section, &
                                                            "ALMO_OPTIMIZER_TRUSTR")
      xalmo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
                                                          "XALMO_OPTIMIZER_PCG")
      xalmo_opt_trustr_section => section_vals_get_subs_vals(almo_scf_section, &
                                                             "XALMO_OPTIMIZER_TRUSTR")
      nlmo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section, &
                                                         "NLMO_OPTIMIZER_PCG")
      almo_analysis_section => section_vals_get_subs_vals(almo_scf_section, "ANALYSIS")
      xalmo_opt_newton_pcg_section => section_vals_get_subs_vals(xalmo_opt_pcg_section, &
                                                                 "XALMO_NEWTON_PCG_SOLVER")
      matrix_iterate_section => section_vals_get_subs_vals(almo_scf_section, &
                                                           "MATRIX_ITERATE")

      ! read user input
      ! common ALMO options
      CALL section_vals_val_get(almo_scf_section, "EPS_FILTER", &
                                r_val=almo_scf_env%eps_filter)
      CALL section_vals_val_get(almo_scf_section, "ALMO_SCF_GUESS", &
                                i_val=almo_scf_env%almo_scf_guess)
      CALL section_vals_val_get(almo_scf_section, "ALMO_ALGORITHM", &
                                i_val=almo_scf_env%almo_update_algorithm)
      CALL section_vals_val_get(almo_scf_section, "XALMO_ALGORITHM", &
                                i_val=almo_scf_env%xalmo_update_algorithm)
      CALL section_vals_val_get(almo_scf_section, "XALMO_TRIAL_WF", &
                                i_val=almo_scf_env%xalmo_trial_wf)
      CALL section_vals_val_get(almo_scf_section, "MO_OVERLAP_INV_ALG", &
                                i_val=almo_scf_env%sigma_inv_algorithm)
      CALL section_vals_val_get(almo_scf_section, "DELOCALIZE_METHOD", &
                                i_val=almo_scf_env%deloc_method)
      CALL section_vals_val_get(almo_scf_section, "XALMO_R_CUTOFF_FACTOR", &
                                r_val=almo_scf_env%quencher_r0_factor)
      CALL section_vals_val_get(almo_scf_section, "ALMO_EXTRAPOLATION_ORDER", &
                                i_val=almo_scf_env%almo_extrapolation_order)
      almo_scf_env%almo_extrapolation_order = MAX(0, almo_scf_env%almo_extrapolation_order)
      CALL section_vals_val_get(almo_scf_section, "XALMO_EXTRAPOLATION_ORDER", &
                                i_val=almo_scf_env%xalmo_extrapolation_order)
      almo_scf_env%xalmo_extrapolation_order = MAX(0, almo_scf_env%xalmo_extrapolation_order)
      CALL section_vals_val_get(almo_scf_section, "RETURN_ORTHOGONALIZED_MOS", &
                                l_val=almo_scf_env%return_orthogonalized_mos)
      CALL section_vals_val_get(almo_scf_section, "CONSTRUCT_NLMOS", &
                                l_val=almo_scf_env%construct_nlmos)

      CALL section_vals_val_get(matrix_iterate_section, "EPS_LANCZOS", &
                                r_val=almo_scf_env%eps_lanczos)
      CALL section_vals_val_get(matrix_iterate_section, "ORDER_LANCZOS", &
                                i_val=almo_scf_env%order_lanczos)
      CALL section_vals_val_get(matrix_iterate_section, "MAX_ITER_LANCZOS", &
                                i_val=almo_scf_env%max_iter_lanczos)
      CALL section_vals_val_get(matrix_iterate_section, "EPS_TARGET_FACTOR", &
                                r_val=almo_scf_env%matrix_iter_eps_error_factor)

      ! optimizers
      CALL section_vals_val_get(almo_opt_diis_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_block_diag_diis%eps_error)
      CALL section_vals_val_get(almo_opt_diis_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_block_diag_diis%max_iter)
      CALL section_vals_val_get(almo_opt_diis_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_block_diag_diis%eps_error_early)
      CALL section_vals_val_get(almo_opt_diis_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_block_diag_diis%max_iter_early)
      CALL section_vals_val_get(almo_opt_diis_section, "N_DIIS", &
                                i_val=almo_scf_env%opt_block_diag_diis%ndiis)

      CALL section_vals_val_get(almo_opt_pcg_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_block_diag_pcg%eps_error)
      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter)
      CALL section_vals_val_get(almo_opt_pcg_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_block_diag_pcg%eps_error_early)
      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter_early)
      CALL section_vals_val_get(almo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_block_diag_pcg%max_iter_outer_loop)
      CALL section_vals_val_get(almo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
                                r_val=almo_scf_env%opt_block_diag_pcg%lin_search_eps_error)
      CALL section_vals_val_get(almo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
                                r_val=almo_scf_env%opt_block_diag_pcg%lin_search_step_size_guess)
      CALL section_vals_val_get(almo_opt_pcg_section, "PRECOND_FILTER_THRESHOLD", &
                                r_val=almo_scf_env%opt_block_diag_pcg%neglect_threshold)
      CALL section_vals_val_get(almo_opt_pcg_section, "CONJUGATOR", &
                                i_val=almo_scf_env%opt_block_diag_pcg%conjugator)
      CALL section_vals_val_get(almo_opt_pcg_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_block_diag_pcg%preconditioner)

      CALL section_vals_val_get(almo_opt_trustr_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_block_diag_trustr%eps_error)
      CALL section_vals_val_get(almo_opt_trustr_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_block_diag_trustr%max_iter)
      CALL section_vals_val_get(almo_opt_trustr_section, "ALGORITHM", &
                                i_val=almo_scf_env%opt_block_diag_trustr%trustr_algorithm)
      CALL section_vals_val_get(almo_opt_trustr_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_block_diag_trustr%eps_error_early)
      CALL section_vals_val_get(almo_opt_trustr_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_block_diag_trustr%max_iter_early)
      CALL section_vals_val_get(almo_opt_trustr_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_block_diag_trustr%max_iter_outer_loop)
      CALL section_vals_val_get(almo_opt_trustr_section, "ETA", &
                                r_val=almo_scf_env%opt_block_diag_trustr%rho_do_not_update)
      almo_scf_env%opt_block_diag_trustr%rho_do_not_update = &
         MIN(MAX(almo_scf_env%opt_block_diag_trustr%rho_do_not_update, 0.0_dp), 0.25_dp)
      CALL section_vals_val_get(almo_opt_trustr_section, "MODEL_GRAD_NORM_RATIO", &
                                r_val=almo_scf_env%opt_block_diag_trustr%model_grad_norm_ratio)
      CALL section_vals_val_get(almo_opt_trustr_section, "INITIAL_TRUST_RADIUS", &
                                r_val=almo_scf_env%opt_block_diag_trustr%initial_trust_radius)
      CALL section_vals_val_get(almo_opt_trustr_section, "MAX_TRUST_RADIUS", &
                                r_val=almo_scf_env%opt_block_diag_trustr%max_trust_radius)
      CALL section_vals_val_get(almo_opt_trustr_section, "CONJUGATOR", &
                                i_val=almo_scf_env%opt_block_diag_trustr%conjugator)
      CALL section_vals_val_get(almo_opt_trustr_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_block_diag_trustr%preconditioner)

      CALL section_vals_val_get(xalmo_opt_trustr_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_xalmo_trustr%eps_error)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "ALGORITHM", &
                                i_val=almo_scf_env%opt_xalmo_trustr%trustr_algorithm)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_xalmo_trustr%eps_error_early)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter_early)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_xalmo_trustr%max_iter_outer_loop)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "ETA", &
                                r_val=almo_scf_env%opt_xalmo_trustr%rho_do_not_update)
      almo_scf_env%opt_xalmo_trustr%rho_do_not_update = &
         MIN(MAX(almo_scf_env%opt_xalmo_trustr%rho_do_not_update, 0.0_dp), 0.25_dp)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "MODEL_GRAD_NORM_RATIO", &
                                r_val=almo_scf_env%opt_xalmo_trustr%model_grad_norm_ratio)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "INITIAL_TRUST_RADIUS", &
                                r_val=almo_scf_env%opt_xalmo_trustr%initial_trust_radius)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "MAX_TRUST_RADIUS", &
                                r_val=almo_scf_env%opt_xalmo_trustr%max_trust_radius)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "CONJUGATOR", &
                                i_val=almo_scf_env%opt_xalmo_trustr%conjugator)
      CALL section_vals_val_get(xalmo_opt_trustr_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_xalmo_trustr%preconditioner)

      CALL section_vals_val_get(xalmo_opt_pcg_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_xalmo_pcg%eps_error)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_xalmo_pcg%eps_error_early)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter_early)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_xalmo_pcg%max_iter_outer_loop)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
                                r_val=almo_scf_env%opt_xalmo_pcg%lin_search_eps_error)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
                                r_val=almo_scf_env%opt_xalmo_pcg%lin_search_step_size_guess)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "PRECOND_FILTER_THRESHOLD", &
                                r_val=almo_scf_env%opt_xalmo_pcg%neglect_threshold)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "CONJUGATOR", &
                                i_val=almo_scf_env%opt_xalmo_pcg%conjugator)
      CALL section_vals_val_get(xalmo_opt_pcg_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_xalmo_pcg%preconditioner)

      penalty_section => section_vals_get_subs_vals(nlmo_opt_pcg_section, "PENALTY")
      CALL section_vals_val_get(nlmo_opt_pcg_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_nlmo_pcg%eps_error)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "EPS_ERROR_EARLY", &
                                r_val=almo_scf_env%opt_nlmo_pcg%eps_error_early)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER_EARLY", &
                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter_early)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_nlmo_pcg%max_iter_outer_loop)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "LIN_SEARCH_EPS_ERROR", &
                                r_val=almo_scf_env%opt_nlmo_pcg%lin_search_eps_error)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
                                r_val=almo_scf_env%opt_nlmo_pcg%lin_search_step_size_guess)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "PRECOND_FILTER_THRESHOLD", &
                                r_val=almo_scf_env%opt_nlmo_pcg%neglect_threshold)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "CONJUGATOR", &
                                i_val=almo_scf_env%opt_nlmo_pcg%conjugator)
      CALL section_vals_val_get(nlmo_opt_pcg_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_nlmo_pcg%preconditioner)
      CALL section_vals_val_get(penalty_section, &
                                "OPERATOR", &
                                i_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%operator_type)
      CALL section_vals_val_get(penalty_section, &
                                "PENALTY_STRENGTH", &
                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength)
      CALL section_vals_val_get(penalty_section, &
                                "PENALTY_STRENGTH_DECREASE_FACTOR", &
                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%penalty_strength_dec_factor)
      CALL section_vals_val_get(penalty_section, &
                                "DETERMINANT_TOLERANCE", &
                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%determinant_tolerance)
      CALL section_vals_val_get(penalty_section, &
                                "FINAL_DETERMINANT", &
                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%final_determinant)
      CALL section_vals_val_get(penalty_section, &
                                "COMPACTIFICATION_FILTER_START", &
                                r_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%compactification_filter_start)
      CALL section_vals_val_get(penalty_section, &
                                "VIRTUAL_NLMOS", &
                                l_val=almo_scf_env%opt_nlmo_pcg%opt_penalty%virtual_nlmos)

      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "EPS_ERROR", &
                                r_val=almo_scf_env%opt_xalmo_newton_pcg_solver%eps_error)
      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "MAX_ITER", &
                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%max_iter)
      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%max_iter_outer_loop)
      CALL section_vals_val_get(xalmo_opt_newton_pcg_section, "PRECONDITIONER", &
                                i_val=almo_scf_env%opt_xalmo_newton_pcg_solver%preconditioner)

      CALL section_vals_val_get(almo_analysis_section, "_SECTION_PARAMETERS_", &
                                l_val=almo_scf_env%almo_analysis%do_analysis)
      CALL section_vals_val_get(almo_analysis_section, "FROZEN_MO_ENERGY_TERM", &
                                i_val=almo_scf_env%almo_analysis%frozen_mo_energy_term)

      !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_AOS",&
      !                          i_val=almo_scf_env%domain_layout_aos)
      !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_MOS",&
      !                          i_val=almo_scf_env%domain_layout_mos)
      !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_AOS",&
      !                          i_val=almo_scf_env%mat_distr_aos)
      !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_MOS",&
      !                          i_val=almo_scf_env%mat_distr_mos)
      !CALL section_vals_val_get(almo_scf_section,"CONSTRAINT_TYPE",&
      !                          i_val=almo_scf_env%constraint_type)
      !CALL section_vals_val_get(almo_scf_section,"MU",&
      !                          r_val=almo_scf_env%mu)
      !CALL section_vals_val_get(almo_scf_section,"FIXED_MU",&
      !                          l_val=almo_scf_env%fixed_mu)
      !CALL section_vals_val_get(almo_scf_section,"EPS_USE_PREV_AS_GUESS",&
      !                          r_val=almo_scf_env%eps_prev_guess)
      !CALL section_vals_val_get(almo_scf_section,"MIXING_FRACTION",&
      !                          r_val=almo_scf_env%mixing_fraction)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_TENSOR_TYPE",&
      !                          i_val=almo_scf_env%deloc_cayley_tensor_type)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_CONJUGATOR",&
      !                          i_val=almo_scf_env%deloc_cayley_conjugator)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_MAX_ITER",&
      !                          i_val=almo_scf_env%deloc_cayley_max_iter)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_USE_OCC_ORBS",&
      !                          l_val=almo_scf_env%deloc_use_occ_orbs)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_USE_VIRT_ORBS",&
      !                          l_val=almo_scf_env%deloc_cayley_use_virt_orbs)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_LINEAR",&
      !                          l_val=almo_scf_env%deloc_cayley_linear)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_EPS_CONVERGENCE",&
      !                          r_val=almo_scf_env%deloc_cayley_eps_convergence)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_OCC_PRECOND",&
      !                          l_val=almo_scf_env%deloc_cayley_occ_precond)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_VIR_PRECOND",&
      !                          l_val=almo_scf_env%deloc_cayley_vir_precond)
      !CALL section_vals_val_get(almo_scf_section,"ALMO_UPDATE_ALGORITHM_BD",&
      !                          i_val=almo_scf_env%almo_update_algorithm)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_TRUNCATE_VIRTUALS",&
      !                          i_val=almo_scf_env%deloc_truncate_virt)
      !CALL section_vals_val_get(almo_scf_section,"DELOC_VIRT_PER_DOMAIN",&
      !                          i_val=almo_scf_env%deloc_virt_per_domain)
      !
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_EPS_CONVERGENCE",&
      !                          r_val=almo_scf_env%opt_k_eps_convergence)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_MAX_ITER",&
      !                          i_val=almo_scf_env%opt_k_max_iter)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_OUTER_MAX_ITER",&
      !                          i_val=almo_scf_env%opt_k_outer_max_iter)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE",&
      !                          r_val=almo_scf_env%opt_k_trial_step_size)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJUGATOR",&
      !                          i_val=almo_scf_env%opt_k_conjugator)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
      !                          r_val=almo_scf_env%opt_k_trial_step_size_multiplier)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_START",&
      !                          i_val=almo_scf_env%opt_k_conj_iter_start)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_START",&
      !                          i_val=almo_scf_env%opt_k_prec_iter_start)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_FREQ_RESET",&
      !                          i_val=almo_scf_env%opt_k_conj_iter_freq)
      !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_FREQ_UPDATE",&
      !                          i_val=almo_scf_env%opt_k_prec_iter_freq)
      !
      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_RADIUS_TYPE",&
      !                          i_val=almo_scf_env%quencher_radius_type)
      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_FACTOR",&
      !                          r_val=almo_scf_env%quencher_r0_factor)
      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_FACTOR",&
      !                          r_val=almo_scf_env%quencher_r1_factor)
      !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_SHIFT",&
      !!                          r_val=almo_scf_env%quencher_r0_shift)
      !!
      !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_SHIFT",&
      !!                          r_val=almo_scf_env%quencher_r1_shift)
      !!
      !!almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
      !!   almo_scf_env%quencher_r0_shift,"angstrom")
      !!almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
      !!   almo_scf_env%quencher_r1_shift,"angstrom")
      !
      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_0",&
      !                          r_val=almo_scf_env%quencher_s0)
      !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_1",&
      !                          r_val=almo_scf_env%quencher_s1)

      !CALL section_vals_val_get(almo_scf_section,"ENVELOPE_AMPLITUDE",&
      !                          r_val=almo_scf_env%envelope_amplitude)

      !! how to read lists
      !CALL section_vals_val_get(almo_scf_section,"INT_LIST01", &
      !        n_rep_val=n_rep)
      !counter_i = 0
      !DO k = 1,n_rep
      !  CALL section_vals_val_get(almo_scf_section,"INT_LIST01",&
      !          i_rep_val=k,i_vals=tmplist)
      !   DO jj = 1,SIZE(tmplist)
      !      counter_i=counter_i+1
      !      almo_scf_env%charge_of_domain(counter_i)=tmplist(jj)
      !   ENDDO
      !ENDDO

      almo_scf_env%domain_layout_aos = almo_domain_layout_molecular
      almo_scf_env%domain_layout_mos = almo_domain_layout_molecular
      almo_scf_env%mat_distr_aos = almo_mat_distr_molecular
      almo_scf_env%mat_distr_mos = almo_mat_distr_molecular

      almo_scf_env%constraint_type = almo_constraint_distance
      almo_scf_env%mu = -0.1_dp
      almo_scf_env%fixed_mu = .FALSE.
      almo_scf_env%mixing_fraction = 0.45_dp
      almo_scf_env%eps_prev_guess = almo_scf_env%eps_filter/1000.0_dp

      almo_scf_env%deloc_cayley_tensor_type = tensor_orthogonal
      almo_scf_env%deloc_cayley_conjugator = cg_hager_zhang
      almo_scf_env%deloc_cayley_max_iter = 100
      almo_scf_env%deloc_use_occ_orbs = .TRUE.
      almo_scf_env%deloc_cayley_use_virt_orbs = .FALSE.
      almo_scf_env%deloc_cayley_linear = .FALSE.
      almo_scf_env%deloc_cayley_eps_convergence = 1.0E-6_dp
      almo_scf_env%deloc_cayley_occ_precond = .TRUE.
      almo_scf_env%deloc_cayley_vir_precond = .TRUE.
      almo_scf_env%deloc_truncate_virt = virt_full
      almo_scf_env%deloc_virt_per_domain = -1

      almo_scf_env%opt_k_eps_convergence = 1.0E-5_dp
      almo_scf_env%opt_k_max_iter = 100
      almo_scf_env%opt_k_outer_max_iter = 1
      almo_scf_env%opt_k_trial_step_size = 0.05_dp
      almo_scf_env%opt_k_conjugator = cg_hager_zhang
      almo_scf_env%opt_k_trial_step_size_multiplier = 1.05_dp
      almo_scf_env%opt_k_conj_iter_start = 0
      almo_scf_env%opt_k_prec_iter_start = 0
      almo_scf_env%opt_k_conj_iter_freq = 10000000
      almo_scf_env%opt_k_prec_iter_freq = 10000000

      almo_scf_env%quencher_radius_type = do_bondparm_vdw
      almo_scf_env%quencher_r1_factor = almo_scf_env%quencher_r0_factor
      !almo_scf_env%quencher_r0_shift=0.0_dp
      !almo_scf_env%quencher_r1_shift=0.0_dp
      !almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
      !   almo_scf_env%quencher_r0_shift,"angstrom")
      !almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
      !   almo_scf_env%quencher_r1_shift,"angstrom")

      almo_scf_env%quencher_s0 = 1.0E-4_dp
      almo_scf_env%quencher_s1 = 1.0E-6_dp

      almo_scf_env%envelope_amplitude = 1.0_dp

      almo_scf_env%logical01 = .FALSE. ! md in eDOF space
      almo_scf_env%logical02 = .TRUE. ! not used
      almo_scf_env%logical03 = .TRUE. ! not used
      almo_scf_env%logical04 = .TRUE. ! use preconditioner
      almo_scf_env%logical05 = .FALSE. ! optimize theta

      almo_scf_env%real01 = almo_scf_env%eps_filter/10.0_dp ! skip gradients
      almo_scf_env%real02 = 0.0_dp ! not used
      almo_scf_env%real03 = 0.0_dp ! not used
      almo_scf_env%real04 = 0.5_dp ! mixing s-f precond

      almo_scf_env%integer01 = 10 ! start eDOF-md
      almo_scf_env%integer02 = 4 ! preconditioner type
      almo_scf_env%integer03 = 0 ! not used
      almo_scf_env%integer04 = 0 ! fixed number of line searches (no grad)
      almo_scf_env%integer05 = 0 ! not used

      ! check for conflicts between options
      IF (almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out .AND. &
          almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_trustr) THEN
         CPABORT("Trust region algorithm cannot optimize projected XALMOs")
      END IF

      CALL section_vals_val_get(almo_scf_section, "XALMO_ALGORITHM", &
                                i_val=almo_scf_env%xalmo_update_algorithm)
      CALL section_vals_val_get(almo_scf_section, "XALMO_TRIAL_WF", &
                                i_val=almo_scf_env%xalmo_trial_wf)
      IF (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag .AND. &
          almo_scf_env%xalmo_update_algorithm .NE. almo_scf_diag) THEN
         CPABORT("1-step delocalization correction requires a different algorithm")
      END IF

      IF (almo_scf_env%xalmo_trial_wf .EQ. xalmo_trial_r0_out .AND. &
          almo_scf_env%almo_update_algorithm .EQ. almo_scf_skip .AND. &
          almo_scf_env%almo_scf_guess .NE. molecular_guess) THEN
         CPABORT("R0 projector requires optimized ALMOs")
      END IF

      IF (almo_scf_env%deloc_method .EQ. almo_deloc_none .AND. &
          almo_scf_env%almo_update_algorithm .EQ. almo_scf_skip) THEN
         CPABORT("No optimization requested")
      END IF

      IF (almo_scf_env%deloc_truncate_virt .EQ. virt_number .AND. &
          almo_scf_env%deloc_virt_per_domain .LE. 0) THEN
         CPABORT("specify a positive number of virtual orbitals")
      END IF

      IF (almo_scf_env%deloc_truncate_virt .EQ. virt_minimal) THEN
         CPABORT("VIRT TRUNCATION TO MINIMAL BASIS IS NIY")
      END IF

      IF (almo_scf_env%domain_layout_mos .NE. almo_domain_layout_molecular) THEN
         CPABORT("use MOLECULAR domains")
      END IF

      IF (almo_scf_env%domain_layout_aos .NE. almo_domain_layout_molecular) THEN
         CPABORT("use MOLECULAR domains")
      END IF

      IF (almo_scf_env%mat_distr_mos .NE. almo_mat_distr_molecular) THEN
         CPABORT("use MOLECULAR distr for MOs")
      END IF

      IF (almo_scf_env%mat_distr_aos == almo_mat_distr_molecular .AND. &
          almo_scf_env%domain_layout_aos == almo_domain_layout_atomic) THEN
         CPABORT("AO blocks cannot be larger than domains")
      END IF

      IF (almo_scf_env%mat_distr_mos == almo_mat_distr_molecular .AND. &
          almo_scf_env%domain_layout_mos == almo_domain_layout_atomic) THEN
         CPABORT("MO blocks cannot be larger than domains")
      END IF

      IF (almo_scf_env%quencher_r1_factor .GT. almo_max_cutoff_multiplier) THEN
         CALL cp_abort(__LOCATION__, &
                       "XALMO_R_CUTOFF_FACTOR is larger than almo_max_cutoff_multiplier. "// &
                       "Increase the hard-coded almo_max_cutoff_multiplier")
      END IF

      ! check analysis settings
      IF (almo_scf_env%almo_analysis%do_analysis) THEN

         IF (almo_scf_env%almo_analysis%frozen_mo_energy_term == almo_frz_crystal &
             .AND. almo_scf_env%almo_scf_guess .NE. molecular_guess) THEN
            CPABORT("To compute frozen-MO energy term set ALMO_SCF_GUESS MOLECULAR")
         END IF

      END IF ! end analysis settings

      CALL timestop(handle)

   END SUBROUTINE almo_scf_init_read_write_input

END MODULE almo_scf_env_methods

