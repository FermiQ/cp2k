!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculation section for TreeMonteCarlo
!> \par History
!>      11.2012 created [Mandes Schoenherr]
!> \author Mandes
! **************************************************************************************************

MODULE tmc_calculations
   USE cell_methods,                    ONLY: init_cell
   USE cell_types,                      ONLY: cell_copy,&
                                              cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_to_string
   USE f77_interface,                   ONLY: calc_energy,&
                                              calc_force,&
                                              set_cell
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE physcon,                         ONLY: boltzmann,&
                                              joule
   USE tmc_move_types,                  ONLY: mv_type_MD
   USE tmc_stati,                       ONLY: task_type_MC,&
                                              task_type_gaussian_adaptation,&
                                              task_type_ideal_gas
   USE tmc_tree_types,                  ONLY: tree_type
   USE tmc_types,                       ONLY: tmc_atom_type,&
                                              tmc_env_type,&
                                              tmc_param_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_calculations'

   PUBLIC :: calc_potential_energy
   PUBLIC :: get_scaled_cell, get_cell_scaling
   PUBLIC :: nearest_distance
   PUBLIC :: geometrical_center, center_of_mass
   PUBLIC :: init_vel, calc_e_kin
   PUBLIC :: compute_estimated_prob
   PUBLIC :: get_subtree_efficiency
CONTAINS

! **************************************************************************************************
!> \brief start the calculation of the energy
!>        (distinguish between exact and approximate)
!> \param conf actual configurations to calculate potential energy
!> \param env_id f77_interface env id
!> \param exact_approx_pot flag if result should be stores in exact or approx
!>        energy variable
!> \param tmc_env TMC environment parameters
!> \author Mandes 01.2013
! **************************************************************************************************
   SUBROUTINE calc_potential_energy(conf, env_id, exact_approx_pot, &
                                    tmc_env)
      TYPE(tree_type), POINTER                           :: conf
      INTEGER                                            :: env_id
      LOGICAL                                            :: exact_approx_pot
      TYPE(tmc_env_type), POINTER                        :: tmc_env

      INTEGER                                            :: ierr
      LOGICAL                                            :: flag
      REAL(KIND=dp)                                      :: e_pot, rnd
      TYPE(cell_type), POINTER                           :: tmp_cell

      rnd = 0.0_dp

      CPASSERT(ASSOCIATED(conf))
      CPASSERT(env_id .GT. 0)
      CPASSERT(ASSOCIATED(tmc_env))

      SELECT CASE (tmc_env%params%task_type)
      CASE (task_type_gaussian_adaptation)
         !CALL gaussian_adaptation_energy(, )
      CASE (task_type_MC)
         IF (tmc_env%params%pressure .GE. 0.0_dp) THEN
            ALLOCATE (tmp_cell)
            CALL get_scaled_cell(cell=tmc_env%params%cell, box_scale=conf%box_scale, &
                                 scaled_cell=tmp_cell)
            CALL set_cell(env_id=env_id, new_cell=tmp_cell%hmat, ierr=ierr)
            CPASSERT(ierr .EQ. 0)
            DEALLOCATE (tmp_cell)
         END IF

         ! TODO check for minimal distances
         flag = .TRUE.
         IF (flag .EQV. .TRUE.) THEN
            IF (tmc_env%params%print_forces .OR. &
                conf%move_type .EQ. mv_type_MD) THEN
               e_pot = 0.0_dp
               conf%frc(:) = 0.0_dp
               CALL calc_force(env_id=env_id, pos=conf%pos, n_el_pos=SIZE(conf%pos), &
                               e_pot=e_pot, force=conf%frc, &
                               n_el_force=SIZE(conf%frc), ierr=ierr)
            ELSE
               e_pot = 0.0_dp
               CALL calc_energy(env_id=env_id, pos=conf%pos, n_el=SIZE(conf%pos), e_pot=e_pot, ierr=ierr)
            END IF
         ELSE
            e_pot = HUGE(e_pot)
         END IF
      CASE (task_type_ideal_gas)
         e_pot = 0.0_dp
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "worker task typ is unknown "// &
                       cp_to_string(tmc_env%params%task_type))
      END SELECT

      ! ---     wait a bit
      rnd = tmc_env%rng_stream%next()
      !rnd = 0.5
!TODO    IF(worker_random_wait.AND.exact_approx_pot)THEN
!      CALL SYSTEM_CLOCK(time0, time_rate, time_max)
!      wait_end=time0+(1.0+rnd)*worker_wait_msec*time_rate/1000.0
!      !wait_end=time0+((worker_wait_msec*time_rate+999)/1000)
!      time_wait: DO
!        CALL SYSTEM_CLOCK(time1, time_rate, time_max)
!        IF(time1<time0.OR.time1>wait_end) exit time_wait
!      END DO time_wait
!    END IF
      IF (exact_approx_pot) THEN
         conf%potential = e_pot
      ELSE
         conf%e_pot_approx = e_pot
      END IF
   END SUBROUTINE calc_potential_energy

! **************************************************************************************************
!> \brief handles properties and calculations of a scaled cell
!> \param cell original cell
!> \param box_scale scaling factors for each direction
!> \param scaled_hmat returns the scaled h matrix (matrix of cell vectors)
!> \param scaled_cell ...
!> \param vol returns the cell volume
!> \param abc ...
!> \param vec a vector, which will be folded (pbc) in the cell
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE get_scaled_cell(cell, box_scale, scaled_hmat, scaled_cell, vol, &
                              abc, vec)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(:), POINTER               :: box_scale
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: scaled_hmat
      TYPE(cell_type), OPTIONAL, POINTER                 :: scaled_cell
      REAL(KIND=dp), OPTIONAL                            :: vol
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: abc
      REAL(KIND=dp), DIMENSION(3), OPTIONAL              :: vec

      LOGICAL                                            :: new_scaled_cell
      TYPE(cell_type), POINTER                           :: tmp_cell

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(box_scale))

      new_scaled_cell = .FALSE.

      IF (.NOT. PRESENT(scaled_cell)) THEN
         ALLOCATE (tmp_cell)
         new_scaled_cell = .TRUE.
      ELSE
         tmp_cell => scaled_cell
      END IF
      CALL cell_copy(cell_in=cell, cell_out=tmp_cell)
      tmp_cell%hmat(:, 1) = tmp_cell%hmat(:, 1)*box_scale(1)
      tmp_cell%hmat(:, 2) = tmp_cell%hmat(:, 2)*box_scale(2)
      tmp_cell%hmat(:, 3) = tmp_cell%hmat(:, 3)*box_scale(3)
      CALL init_cell(cell=tmp_cell)

      IF (PRESENT(scaled_hmat)) &
         scaled_hmat(:, :) = tmp_cell%hmat

      IF (PRESENT(vec)) THEN
         vec = pbc(r=vec, cell=tmp_cell)
      END IF

      IF (PRESENT(vol)) CALL get_cell(cell=tmp_cell, deth=vol)
      IF (PRESENT(abc)) CALL get_cell(cell=tmp_cell, abc=abc)
      IF (new_scaled_cell) DEALLOCATE (tmp_cell)

   END SUBROUTINE get_scaled_cell

! **************************************************************************************************
!> \brief handles properties and calculations of a scaled cell
!> \param cell original cell
!> \param scaled_hmat returns the scaled h matrix (matrix of cell vectors)
!> \param box_scale scaling factors for each direction
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE get_cell_scaling(cell, scaled_hmat, box_scale)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: scaled_hmat
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: box_scale

      REAL(KIND=dp), DIMENSION(3)                        :: abc_new, abc_orig
      TYPE(cell_type), POINTER                           :: tmp_cell

      CPASSERT(ASSOCIATED(cell))

      ALLOCATE (tmp_cell)
      CALL cell_copy(cell_in=cell, cell_out=tmp_cell)
      tmp_cell%hmat(:, :) = scaled_hmat(:, :)
      CALL init_cell(cell=tmp_cell)
      CALL get_cell(cell=cell, abc=abc_orig)
      CALL get_cell(cell=tmp_cell, abc=abc_new)

      box_scale(:) = abc_new(:)/abc_orig(:)

      DEALLOCATE (tmp_cell)
   END SUBROUTINE get_cell_scaling

! **************************************************************************************************
!> \brief neares distance of atoms within the periodic boundary condition
!> \param x1 ...
!> \param x2 ...
!> \param cell ...
!> \param box_scale ...
!> \return ...
!> \author Mandes 11.2012
! **************************************************************************************************
   FUNCTION nearest_distance(x1, x2, cell, box_scale) RESULT(res)
      REAL(KIND=dp), DIMENSION(:)                        :: x1, x2
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: box_scale
      REAL(KIND=dp)                                      :: res

      REAL(KIND=dp), DIMENSION(3)                        :: dist_vec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_box_scale

      NULLIFY (tmp_box_scale)

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(SIZE(x1) .EQ. 3)
      CPASSERT(SIZE(x2) .EQ. 3)

      dist_vec(:) = x2(:) - x1(:) ! distance vector between atoms
      ALLOCATE (tmp_box_scale(3))
      IF (PRESENT(box_scale)) THEN
         CPASSERT(SIZE(box_scale) .EQ. 3)
         tmp_box_scale(:) = box_scale
      ELSE
         tmp_box_scale(:) = 1.0_dp
      END IF
      CALL get_scaled_cell(cell=cell, box_scale=box_scale, vec=dist_vec)
      res = SQRT(SUM(dist_vec(:)*dist_vec(:)))
      DEALLOCATE (tmp_box_scale)
   END FUNCTION nearest_distance

! **************************************************************************************************
!> \brief calculate the geometrical center of an amount of atoms
!>        array size should be multiple of dim_per_elem
!> \param pos list of atoms
!> \param center return value, the geometrical center
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE geometrical_center(pos, center)
      REAL(KIND=dp), DIMENSION(:)                        :: pos
      REAL(KIND=dp), DIMENSION(:), POINTER               :: center

      CHARACTER(LEN=*), PARAMETER :: routineN = 'geometrical_center'

      INTEGER                                            :: handle, i

      CPASSERT(ASSOCIATED(center))
      CPASSERT(SIZE(pos) .GE. SIZE(center))

      ! start the timing
      CALL timeset(routineN, handle)

      center = 0.0_dp
      DO i = 1, SIZE(pos), SIZE(center)
         center(:) = center(:) + &
                     pos(i:i + SIZE(center) - 1)/(SIZE(pos)/REAL(SIZE(center), KIND=dp))
      END DO
      ! end the timing
      CALL timestop(handle)
   END SUBROUTINE geometrical_center

! **************************************************************************************************
!> \brief calculate the center of mass of an amount of atoms
!>        array size should be multiple of dim_per_elem
!> \param pos ...
!> \param atoms ...
!> \param center ...
!> \param
!> \param
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE center_of_mass(pos, atoms, center)
      REAL(KIND=dp), DIMENSION(:)                        :: pos
      TYPE(tmc_atom_type), DIMENSION(:), OPTIONAL        :: atoms
      REAL(KIND=dp), DIMENSION(:), POINTER               :: center

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'center_of_mass'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: mass_sum, mass_tmp

      CPASSERT(ASSOCIATED(center))
      CPASSERT(SIZE(pos) .GE. SIZE(center))

      ! start the timing
      CALL timeset(routineN, handle)

      center = 0.0_dp
      mass_sum = 0.0_dp
      DO i = 1, SIZE(pos), SIZE(center)
         IF (PRESENT(atoms)) THEN
            CPASSERT(SIZE(atoms) .EQ. SIZE(pos)/SIZE(center))
            mass_tmp = atoms(INT(i/REAL(SIZE(center), KIND=dp)) + 1)%mass
            center(:) = center(:) + pos(i:i + SIZE(center) - 1)/ &
                        (SIZE(pos)/REAL(SIZE(center), KIND=dp))*mass_tmp
            mass_sum = mass_sum + mass_tmp
         ELSE
            CPWARN("try to calculate center of mass without any mass.")
            center(:) = center(:) + pos(i:i + SIZE(center) - 1)/ &
                        (SIZE(pos)/REAL(SIZE(center), KIND=dp))
            mass_sum = 1.0_dp
         END IF
      END DO
      center(:) = center(:)/mass_sum
      ! end the timing
      CALL timestop(handle)
   END SUBROUTINE center_of_mass

! **************************************************************************************************
!> \brief routine sets initial velocity, using the Box-Muller Method for Normal
!>         (Gaussian) Deviates
!> \param vel ...
!> \param atoms ...
!> \param temerature ...
!> \param rng_stream ...
!> \param rnd_seed ...
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE init_vel(vel, atoms, temerature, rng_stream, rnd_seed)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vel
      TYPE(tmc_atom_type), DIMENSION(:), POINTER         :: atoms
      REAL(KIND=dp)                                      :: temerature
      TYPE(rng_stream_type), INTENT(INOUT)               :: rng_stream
      REAL(KIND=dp), DIMENSION(3, 2, 3)                  :: rnd_seed

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: kB, mass_tmp, rnd1, rnd2

      kB = boltzmann/joule

      CPASSERT(ASSOCIATED(vel))
      CPASSERT(ASSOCIATED(atoms))

      CALL rng_stream%set(bg=rnd_seed(:, :, 1), cg=rnd_seed(:, :, 2), ig=rnd_seed(:, :, 3))
      DO i = 1, SIZE(vel)
         rnd1 = rng_stream%next()
         rnd2 = rng_stream%next()

         mass_tmp = atoms(INT(i/REAL(3, KIND=dp)) + 1)%mass

         vel(i) = SQRT(-2.0_dp*LOG(rnd1))*COS(2.0_dp*PI*rnd2)* &
                  SQRT(kB*temerature/mass_tmp)
      END DO
      CALL rng_stream%get(bg=rnd_seed(:, :, 1), cg=rnd_seed(:, :, 2), ig=rnd_seed(:, :, 3))

   END SUBROUTINE init_vel

! **************************************************************************************************
!> \brief routine calculates the kinetic energy, using the velocities
!>        and atom mass, both in atomic units
!> \param vel ...
!> \param atoms ...
!> \return ...
!> \author Mandes 11.2012
! **************************************************************************************************
   FUNCTION calc_e_kin(vel, atoms) RESULT(ekin)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vel
      TYPE(tmc_atom_type), DIMENSION(:), POINTER         :: atoms
      REAL(KIND=dp)                                      :: ekin

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mass_tmp

      CPASSERT(ASSOCIATED(vel))
      CPASSERT(ASSOCIATED(atoms))
      ekin = 0.0_dp

      DO i = 1, SIZE(vel)
         mass_tmp = atoms(INT(i/REAL(3, KIND=dp)) + 1)%mass
         ekin = ekin + 0.5_dp*mass_tmp*vel(i)*vel(i)
      END DO
   END FUNCTION calc_e_kin

! **************************************************************************************************
!> \brief assuming an (exponential) decreasing function, this function
!>        extrapolate the converged value
!> \param v1 function values
!> \param v2 function values
!> \param v3 function values
!> \param extrapolate extrapolated final value (result)
!> \param res_err error of the result
!> \author Mandes 12.2012
! **************************************************************************************************
   SUBROUTINE three_point_extrapolate(v1, v2, v3, extrapolate, res_err)
      REAL(KIND=dp)                            :: v1, v2, v3
      REAL(KIND=dp), INTENT(OUT)               :: extrapolate, res_err

      REAL(KIND=dp)                            :: e1, e2, e3
      REAL(KIND=dp)                            :: a, b, c, d12, d23, ddd

      extrapolate = HUGE(extrapolate)

      !> solve({exp(a+b)+c = e1, exp(2*a+b)+c = e2, exp(3*a+b)+c = e3}, [a, b, c])
      !> solve({a*b+c = e1, a^2*b+c = e2, a^3*b+c = e3}, [a, b, c]);
      !   [[                                   3                   2         ]]
      !   [[    -e3 + e2              (e1 - e2)                 -e2  + e1 e3 ]]
      !   [[a = --------, b = ---------------------------, c = --------------]]
      !   [[    e1 - e2       (-e3 + e2) (e3 - 2 e2 + e1)      e3 - 2 e2 + e1]]

      ! sort so that e1>=e2>=e3
      e1 = v1; e2 = v2; e3 = v3
      CALL swap(e1, e2)
      CALL swap(e1, e3)
      CALL swap(e2, e3)
      ! we need extra care if some of the difference e1-e2, e3-e2 are nearly zero,
      !  since the formulae suffer from sever loss of precision
      d12 = e1 - e2
      d23 = e2 - e3
      ddd = d12 - d23
      IF (d12 == 0 .OR. d23 == 0 .OR. ABS(ddd) == 0) THEN
         ! a degenerate case, we do no extrapolation
         extrapolate = e3
         res_err = e1 - e3
      ELSE
         a = d23/d12
         b = (d12**3/(d23*ddd))
         c = e2 - (d12*d23)/ddd
         ! extrapolation, let's only look 4 iterations ahead, more is presumably anyway not accurate
         ! fewer is maybe more stable
         extrapolate = a**7*b + c
         res_err = e3 - extrapolate
      END IF
      CPASSERT(extrapolate .NE. HUGE(extrapolate))
   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param x1 ...
!> \param x2 ...
! **************************************************************************************************
      SUBROUTINE swap(x1, x2)
      REAL(KIND=dp)                                      :: x1, x2

      REAL(KIND=dp)                                      :: tmp

         IF (x2 > x1) THEN
            tmp = x2
            x2 = x1
            x1 = tmp
         END IF
      END SUBROUTINE swap
   END SUBROUTINE three_point_extrapolate

! **************************************************************************************************
!> \brief calculates the probability of acceptance for given intervals of the
!>        exact energy
!> \param E_n_mu energy distribution of new configuration
!> \param E_n_sigma energy distribution of new configuration
!> \param E_o_mu energy distribution of old configuration
!> \param E_o_sigma energy distribution of old configuration
!> \param E_classical_diff the difference in approximated energies for the
!>        old and new configuration (E_o-E_n)
!> \param prior_mu energy distribution of the already converged
!>         energies
!> \param prior_sigma energy distribution of the already converged
!>         energies
!> \param p the random number, the criteria has to be smaller than this
!> \param beta ...
!> \return return probability of acceptance
!> \author Mandes 12.2012
! **************************************************************************************************
   FUNCTION compute_prob(E_n_mu, E_n_sigma, E_o_mu, E_o_sigma, E_classical_diff, &
                         prior_mu, prior_sigma, p, beta) RESULT(prob)
      REAL(KIND=dp)                                      :: E_n_mu, E_n_sigma, E_o_mu, E_o_sigma, &
                                                            E_classical_diff, prior_mu, &
                                                            prior_sigma, p, beta, prob

!    INTEGER       :: io,in
!    REAL(KIND=dp) :: diff,E_n,E_o,surface,lower_bound,upper_bound,delta

      prob = 0.5_dp*ERFC(-0.5_dp*SQRT(2.0_dp)*( &
                         (-prior_sigma**2 - E_o_sigma**2 - E_n_sigma**2)*LOG(p) + &
                         ((E_classical_diff - E_n_mu + E_o_mu)*prior_sigma**2 - prior_mu*(E_n_sigma**2 + E_o_sigma**2))*beta)/ &
                         (SQRT(E_o_sigma**2 + E_n_sigma**2)*SQRT(prior_sigma**2 + E_o_sigma**2 + E_n_sigma**2)*prior_sigma*beta))

      prob = MIN(1.0_dp - EPSILON(1.0_dp), MAX(EPSILON(1.0_dp), prob))

   END FUNCTION compute_prob

! **************************************************************************************************
!> \brief extimates the probability of acceptance considering the intermetiate
!>        step energies
!> \param elem_old old/parent sub tree element
!> \param elem_new new/actual sub tree element, which schould be checked
!> \param E_classical_diff difference in the classical energy of the old and
!>        new configuration
!> \param rnd_nr random number acceptance check will be done with
!> \param beta 1/(kB*T) can differ for different acceptance checks
!> \param tmc_params TMC environment parameters
!> \return estimated acceptance probability
!> \author Mandes 12.2012
! **************************************************************************************************
   FUNCTION compute_estimated_prob(elem_old, elem_new, E_classical_diff, &
                                   rnd_nr, beta, tmc_params) RESULT(prob)
      TYPE(tree_type), POINTER                           :: elem_old, elem_new
      REAL(KIND=dp)                                      :: E_classical_diff, rnd_nr, beta
      TYPE(tmc_param_type), POINTER                      :: tmc_params
      REAL(KIND=dp)                                      :: prob

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_estimated_prob'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: E_mu_tmp, E_n_mu, E_n_sigma, E_o_mu, &
                                                            E_o_sigma, E_sigma_tmp, prior_sigma

      CPASSERT(ASSOCIATED(elem_old))
      CPASSERT(ASSOCIATED(elem_new))
      CPASSERT(rnd_nr .GT. 0.0_dp)

      ! start the timing
      CALL timeset(routineN, handle)

      prob = -1.0_dp
      IF ((elem_new%scf_energies_count .GE. 3) .AND. &
          (elem_old%scf_energies_count .GE. 3) .AND. &
          tmc_params%prior_NMC_acc%counter .GE. 10) THEN
         !-- first the new element energy estimation
         ! using 3 point extrapolation of two different intervals -> more stable estimation
         ! the energies are sorted in the three_point_extrapolate routine !
         ! But with array of length 4 we have to select the 3 connected ones
         CALL three_point_extrapolate(v1=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 3, 4) + 1), &
                                      v2=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 2, 4) + 1), &
                                      v3=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 1, 4) + 1), &
                                      extrapolate=E_mu_tmp, res_err=E_sigma_tmp)
         IF ((elem_new%scf_energies_count .GT. 3)) THEN
            CALL three_point_extrapolate(v1=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 4, 4) + 1), &
                                         v2=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 3, 4) + 1), &
                                         v3=elem_new%scf_energies(MOD(elem_new%scf_energies_count - 2, 4) + 1), &
                                         extrapolate=E_n_mu, res_err=E_n_sigma)
            E_n_sigma = MAX(E_n_sigma, ABS(E_n_mu - E_mu_tmp))
         ELSE
            E_n_sigma = E_sigma_tmp
            E_n_mu = E_mu_tmp
         END IF

         !-- the old/parent element energy estimation
         CALL three_point_extrapolate(v1=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 3, 4) + 1), &
                                      v2=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 2, 4) + 1), &
                                      v3=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 1, 4) + 1), &
                                      extrapolate=E_mu_tmp, res_err=E_sigma_tmp)
         IF ((elem_old%scf_energies_count .GT. 3)) THEN
            CALL three_point_extrapolate(v1=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 4, 4) + 1), &
                                         v2=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 3, 4) + 1), &
                                         v3=elem_old%scf_energies(MOD(elem_old%scf_energies_count - 2, 4) + 1), &
                                         extrapolate=E_o_mu, res_err=E_o_sigma)
            E_o_sigma = MAX(E_o_sigma, ABS(E_o_mu - E_mu_tmp))
         ELSE
            E_o_sigma = E_sigma_tmp
            E_o_mu = E_mu_tmp
         END IF

         ! calculate the estimation for the average of the trajectory elements
         prior_sigma = SQRT(ABS(tmc_params%prior_NMC_acc%aver_2 &
                                - tmc_params%prior_NMC_acc%aver**2))

         ! calculate the probability of acceptance for those two elements with their energy
         ! swap and 2 potential moves are distinguished using the difference in classical energy and different betas
         prob = compute_prob(E_n_mu=E_n_mu, E_n_sigma=E_n_sigma, E_o_mu=E_o_mu, E_o_sigma=E_o_sigma, &
                             E_classical_diff=E_classical_diff, &
                             prior_mu=tmc_params%prior_NMC_acc%aver, prior_sigma=prior_sigma, &
                             p=rnd_nr, beta=beta)
      END IF
      ! end the timing
      CALL timestop(handle)
   END FUNCTION compute_estimated_prob

! **************************************************************************************************
!> \brief calculated the rate of used tree elements to created tree elements
!>        for every temperature
!> \param tmc_env TMC environment variables
!> \param eff result efficiency
!> \author Mandes 01.2013
! **************************************************************************************************
   SUBROUTINE get_subtree_efficiency(tmc_env, eff)
      TYPE(tmc_env_type), POINTER                        :: tmc_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eff

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(tmc_env))
      CPASSERT(ASSOCIATED(tmc_env%params))
      CPASSERT(ASSOCIATED(tmc_env%m_env))

      eff(:) = 0.0_dp

      DO i = 1, tmc_env%params%nr_temp
         IF (tmc_env%m_env%tree_node_count(i) > 0) &
            eff(i) = tmc_env%params%move_types%mv_count(0, i)/ &
                     (tmc_env%m_env%tree_node_count(i)*1.0_dp)
         eff(0) = eff(0) + tmc_env%params%move_types%mv_count(0, i)/ &
                  (SUM(tmc_env%m_env%tree_node_count(1:))*1.0_dp)
      END DO
   END SUBROUTINE get_subtree_efficiency
END MODULE tmc_calculations
