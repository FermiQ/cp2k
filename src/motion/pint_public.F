!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Public path integral routines that can be called from other modules
!> \author Lukasz Walewski
!> \date   2009-07-24
!> \note   Avoiding circular dependencies: please design new members of this
!>         module in such a way that they use pint_types module only.
! **************************************************************************************************
MODULE pint_public

   USE kinds,                           ONLY: dp
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE pint_types,                      ONLY: pint_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_public'

   PUBLIC :: pint_com_pos
   PUBLIC :: pint_levy_walk
   PUBLIC :: pint_calc_centroid

CONTAINS

! ***************************************************************************
!> \brief  Return the center of mass of the PI system
!> \param pint_env ...
!> \return ...
!> \date   2009-07-24
!> \par    History
!>           2009-11-30 fixed serious bug in pint_env%x indexing [lwalewski]
!> \author Lukasz Walewski
! **************************************************************************************************
   PURE FUNCTION pint_com_pos(pint_env) RESULT(com_r)

      TYPE(pint_env_type), INTENT(IN)                    :: pint_env
      REAL(kind=dp), DIMENSION(3)                        :: com_r

      INTEGER                                            :: ia, ib, ic
      REAL(kind=dp)                                      :: tmass

      tmass = 0.0_dp
      com_r(:) = 0.0_dp
      DO ia = 1, pint_env%ndim/3
         DO ib = 1, pint_env%p
            DO ic = 1, 3
               com_r(ic) = com_r(ic) + &
                           pint_env%x(ib, (ia - 1)*3 + ic)*pint_env%mass((ia - 1)*3 + ic)
               tmass = tmass + pint_env%mass((ia - 1)*3 + ic)
            END DO
         END DO
      END DO
      ! pint_env%mass is REAL, DIMENSION(NDIM) which means that each atom
      ! has its mass defined three times - here we hope that all three
      ! values are equal
      tmass = tmass/3.0_dp
      com_r(:) = com_r(:)/tmass
   END FUNCTION pint_com_pos

! ***************************************************************************
!> \brief  Return the center of geometry of the PI system
!> \param pint_env ...
!> \return ...
!> \date   2009-11-30
!> \author Lukasz Walewski
! **************************************************************************************************
   PURE FUNCTION pint_cog_pos(pint_env) RESULT(cntrd_r)

      TYPE(pint_env_type), INTENT(IN)                    :: pint_env
      REAL(kind=dp), DIMENSION(3)                        :: cntrd_r

      INTEGER                                            :: ia, ib, ic, natoms

      cntrd_r(:) = 0.0_dp
      natoms = pint_env%ndim/3
      DO ia = 1, natoms
         DO ib = 1, pint_env%p
            DO ic = 1, 3
               cntrd_r(ic) = cntrd_r(ic) + pint_env%x(ib, (ia - 1)*3 + ic)
            END DO
         END DO
      END DO
      cntrd_r(:) = cntrd_r(:)/REAL(pint_env%p, dp)/REAL(natoms, dp)
   END FUNCTION pint_cog_pos

! ***************************************************************************
!> \brief  Given the number of beads n and the variance t returns the
!>         positions of the beads for a non-interacting particle.
!> \param n ...
!> \param t ...
!> \param rng_gaussian ...
!> \param x ...
!> \param nout ...
!> \date   2010-12-13
!> \author Lukasz Walewski
!> \note  This routine implements Levy argorithm (see e.g. Rev. Mod. Phys.
!>         67 (1995) 279, eq. 5.35) and requires that n is a power of 2. The
!>         resulting bead positions are centered around (0,0,0).
! **************************************************************************************************
   SUBROUTINE pint_free_part_bead_x(n, t, rng_gaussian, x, nout)
!
!TODO this routine gives wrong spread of the particles, please fix before usage.
!
      INTEGER, INTENT(IN)                                :: n
      REAL(kind=dp), INTENT(IN)                          :: t
      TYPE(rng_stream_type), INTENT(INOUT)               :: rng_gaussian
      REAL(kind=dp), DIMENSION(:), POINTER               :: x
      INTEGER, INTENT(OUT)                               :: nout

      INTEGER                                            :: dl, i1, i2, ib, ic, il, ip, j, nlevels, &
                                                            np
      REAL(kind=dp)                                      :: rtmp, tcheck, vrnc, xc
      REAL(kind=dp), DIMENSION(3)                        :: cntrd_r

      nout = 0

      IF (n < 1) THEN
         RETURN
      END IF

      ! if number of beads is not a power of 2 return
      nlevels = NINT(LOG(REAL(n, KIND=dp))/LOG(2.0_dp))
      rtmp = 2**nlevels
      tcheck = ABS(REAL(n, KIND=dp) - rtmp)
      IF (tcheck > 100.0_dp*EPSILON(0.0_dp)) THEN
         RETURN
      END IF

      ! generate at least first point at (0,0,0)
      DO ic = 1, 3
         x(ic) = 0.0_dp
      END DO
      nout = 1

      ! loop over Levy levels
      vrnc = 2.0_dp*t
      DO il = 0, nlevels - 1

         np = 2**il ! number of points to be generated at this level
         dl = n/(2*np) ! interval betw points (in index numbers)
         vrnc = vrnc/2.0_dp; ! variance at this level (=t at level 0)

         ! loop over points added in this level
         DO ip = 0, np - 1

            j = (2*ip + 1)*dl ! index of currently generated point

            ! indices of two points betw which to generate a new point
            i1 = 2*dl*ip
            i2 = 2*dl*(ip + 1)
            IF (i2 .EQ. n) THEN
               i2 = 0
            END IF

            ! generate new point and save it under j
            DO ic = 1, 3
               xc = (x(3*i1 + ic) + x(3*i2 + ic))/2.0
               xc = xc + rng_gaussian%next(variance=vrnc)
               x(3*j + ic) = xc
            END DO
            nout = nout + 1

         END DO
      END DO

      ! translate the centroid to the origin
      cntrd_r(:) = 0.0_dp
      DO ib = 1, n
         DO ic = 1, 3
            cntrd_r(ic) = cntrd_r(ic) + x((ib - 1)*3 + ic)
         END DO
      END DO
      cntrd_r(:) = cntrd_r(:)/REAL(n, dp)
      DO ib = 1, n
         DO ic = 1, 3
            x((ib - 1)*3 + ic) = x((ib - 1)*3 + ic) - cntrd_r(ic)
         END DO
      END DO

   END SUBROUTINE pint_free_part_bead_x

! ***************************************************************************
!> \brief  Perform a Brownian walk of length n around x0 with the variance v.
!> \param x0 ...
!> \param n ...
!> \param v ...
!> \param x ...
!> \param rng_gaussian ...
!> \date   2011-01-06
!> \author Lukasz Walewski
!> \note  This routine implements Levy argorithm (Phys. Rev. 143 (1966) 58)
! **************************************************************************************************
   SUBROUTINE pint_levy_walk(x0, n, v, x, rng_gaussian)

      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: x0
      INTEGER, INTENT(IN)                                :: n
      REAL(kind=dp), INTENT(IN)                          :: v
      REAL(kind=dp), DIMENSION(:), POINTER               :: x
      TYPE(rng_stream_type), INTENT(INOUT)               :: rng_gaussian

      INTEGER                                            :: ib, ic
      REAL(kind=dp)                                      :: r, tau_i, tau_i1
      REAL(kind=dp), DIMENSION(3)                        :: cntrd_r

      x(1) = x0(1)
      x(2) = x0(2)
      x(3) = x0(3)
      DO ib = 1, n - 1
         DO ic = 1, 3
            r = rng_gaussian%next(variance=1.0_dp)
            tau_i = (REAL(ib, dp) - 1.0_dp)/REAL(n, dp)
            tau_i1 = (REAL(ib + 1, dp) - 1.0_dp)/REAL(n, dp)
            x(ib*3 + ic) = (x((ib - 1)*3 + ic)*(1.0_dp - tau_i1) + &
                            x(ic)*(tau_i1 - tau_i))/ &
                           (1.0_dp - tau_i) + &
                           r*v*SQRT( &
                           (tau_i1 - tau_i)* &
                           (1.0_dp - tau_i1)/ &
                           (1.0_dp - tau_i) &
                           )
         END DO
      END DO

      ! translate the centroid to the origin
      cntrd_r(:) = 0.0_dp
      DO ib = 1, n
         DO ic = 1, 3
            cntrd_r(ic) = cntrd_r(ic) + x((ib - 1)*3 + ic)
         END DO
      END DO
      cntrd_r(:) = cntrd_r(:)/REAL(n, dp)
      DO ib = 1, n
         DO ic = 1, 3
            x((ib - 1)*3 + ic) = x((ib - 1)*3 + ic) - cntrd_r(ic)
         END DO
      END DO

   END SUBROUTINE pint_levy_walk

! ***************************************************************************
!> \brief  Calculate the centroid
!> \param  pint_env path integral environment
!> \date   2013-02-10
!> \author lwalewski
! **************************************************************************************************
   PURE SUBROUTINE pint_calc_centroid(pint_env)

      TYPE(pint_env_type), INTENT(INOUT)                 :: pint_env

      INTEGER                                            :: ia, ib
      REAL(KIND=dp)                                      :: invp

      invp = 1.0_dp/pint_env%p
      pint_env%centroid(:) = 0.0_dp
      DO ia = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            pint_env%centroid(ia) = pint_env%centroid(ia) + pint_env%x(ib, ia)
         END DO
         pint_env%centroid(ia) = pint_env%centroid(ia)*invp
      END DO

   END SUBROUTINE pint_calc_centroid

END MODULE pint_public
