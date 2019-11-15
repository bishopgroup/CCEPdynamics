module global
! ------------------------------------------------------------ !
! This module contains the global variables used in the
! simulation.  This module is accessed in just about every 
! routine.
! ------------------------------------------------------------ !
implicit none
save

    ! The number of particles, 6*np, 5*np, 3*np, the number of pairs
    integer*4 :: np, n6, n5, n3, npair

    ! The id's of particles in a pair
    integer*4, dimension( :, : ), allocatable :: id

    ! The current and previous positions of the particles
    real*8, dimension( : ), allocatable :: x, x0
    ! The distnaces between particles
    real*8, dimension( :, : ), allocatable :: pd

    ! The non-Brownian velocity of the particles
    real*8, dimension( : ), allocatable :: u
    ! The velocity of the particles due to the shear field
    real*8, dimension( : ), allocatable :: udiff
    ! The non-Brownian velocity of the particles at the previous time step
    real*8, dimension( :, : ), allocatable :: uold
    ! The hydrodynamic contribution to the stresslet
    real*8, dimension( : ), allocatable :: sh

    ! The affine contribution to the particle velocities due to the shear field
    real*8, dimension( : ), allocatable :: uinf
    ! The bulk rate of strain as a vector
    real*8, dimension( : ), allocatable :: einf

    ! The external force and torque on all the particles
    real*8, dimension( 6 ) :: fext_h
    ! The external force and torque on each particle
    real*8, dimension( : ), allocatable :: fext
    ! The velocity of the particles due to the external force and torque
    real*8, dimension( : ), allocatable :: uext

    ! The resistance tensor R_FU
    real*8, dimension( :, : ), allocatable :: rfu
    ! The resistance tensor R_FE
    real*8, dimension( :, : ), allocatable :: rfe
    ! The resistance tensor R_SU
    real*8, dimension( :, : ), allocatable :: rsu
    ! The resistance tensor R_SE
    real*8, dimension( :, : ), allocatable :: rse
    ! The resistance tensor M_UF
    real*8, dimension( :, : ), allocatable :: zmuf 
    ! The resistance tensor M_US
    real*8, dimension( :, : ), allocatable :: zmus
    ! The resistance tensor M_ES
    real*8, dimension( :, : ), allocatable :: zmes

    ! The end time and current time of the simulation
    real*8 :: tf, t
    ! The specified and rescaled time steps
    real*8 :: dt_spec, dt
    ! The iteration counter for the simulation
    integer*4 :: iloop
    ! The time counter for outputting data
    integer*4 :: pos_time    
    ! The current execution clock time and clock rate
    integer*4 :: time_curr, time_init, clock_rate

    ! Distance variables specifing the minimum overlap allowed
    ! and minimum distance between particles treated in mobility tensor.
    real*8 :: drmin, drset, dset

    ! An error flagcoming from the routine cholesky 
    integer*4 :: istop
    ! If (1) use F-T-S.  Otherwise, use F-T.
    integer*4 :: sim_mode
    ! Steps the order of trajectory integrator <= 4.
    integer*4 :: istep
    ! The output mode If (1) output time code, (2) output for movie.
    integer*4 :: out_mode

    ! The self and pair contributions to the grand mobility tensor 
    real*8, dimension( 3, 3 ) :: mob_a
    real*8, dimension( 3, 3 ) :: mob_b
    real*8, dimension( 3, 3 ) :: mob_bt
    real*8, dimension( 3, 3 ) :: mob_c
    real*8, dimension( 3, 5 ) :: mob_gt
    real*8, dimension( 3, 5 ) :: mob_ht
    real*8, dimension( 5, 5 ) :: mob_m


    ! The identity tensor
    real*8, dimension( 3, 3 ) :: delta
    ! The permutation tensor
    real*8, dimension( 3, 3, 3 ) :: eps
    ! ID's for converting M_ES to a symmetric matrix
    integer*4, dimension( 2, 5 ) :: mesid

    ! Constants for use in computing the mobility tensors
    real*8 :: c1d2, c3d2, c3d4, c3d8, c1d3, c9d4, c18d5
    real*8 :: c6d5, c9d8, c9d2, c54d5, c36d5, c9d5, c9d10, c2d3
    real*8 :: c8d15, c8d3, c2d5
    parameter( c1d2 = 1.d0 / 2.d0 )
    parameter( c3d2 = 3.d0 / 2.d0 )
    parameter( c3d4 = 3.d0 / 4.d0 )
    parameter( c3d8 = 3.d0 / 8.d0 )
    parameter( c1d3 = 1.d0 / 3.d0 )
    parameter( c9d4 = 9.d0 / 4.d0 )
    parameter( c18d5 = 18.d0 / 5.d0 )
    parameter( c6d5 = 6.d0 / 5.d0 )
    parameter( c9d8 = 9.d0 / 8.d0 )
    parameter( c9d2 = 9.d0 / 2.d0 )
    parameter( c54d5 = 54.d0 / 5.d0 )
    parameter( c36d5 = 36.d0 / 5.d0 )
    parameter( c9d5 = 9.d0 / 5.d0 )
    parameter( c9d10 = 9.d0 / 10.d0 )
    parameter( c2d3 = 2.d0 / 3.d0 )
    parameter( c8d15 = 8.d0 / 15.d0 )
    parameter( c8d3 = 8.d0 / 3.d0 )
    parameter( c2d5 = 2.d0 / 5.d0 )
    

    real*8     pi, twopi, one36, four3, four15, d2p31m, d2p31m_inv
    real*8     one3, one6, one15, one12, one24, two15, c3d112, c7d120
    real*8     ug1, ug2, ug4, ug5, ug6, ug7

    parameter (pi=3.14159265358979d0)
    parameter (twopi=2.d0*pi)
    parameter (one36=1.d0/36.d0)
    parameter (four3=4.d0/3.d0)
    parameter (four15=4.d0/15.d0)
    parameter (d2p31m=2147483647.0)
    parameter (d2p31m_inv=1.0/d2p31m)

    parameter (one3=1.d0/3.d0)
    parameter (one6=1.d0/6.d0)
    parameter (one15=1.d0/15.d0)
    parameter (one12=1.d0/12.d0)
    parameter (one24=1.d0/24.d0)
    parameter (two15=2.d0/15.d0)
    parameter (c3d112=3.d0/112.d0)
    parameter (c7d120=7.d0/120.d0)
    parameter (ug1=94.d0/375.d0)
    parameter (ug2=31.d0/375.d0)
    parameter (ug4=39.d0/280.d0)
    parameter (ug5=113.d0/1500.d0)
    parameter (ug6=47.d0/840.d0)
    parameter (ug7=137.d0/1500.d0)

end module global


module lubdat
! ------------------------------------------------------------ !
! This module contains data used in some of the lubrication
! calculations.  The data is read from files in the directory
! 'lubdat'.  This module is called in 'init.f90' and 'lub.f90'.
! ------------------------------------------------------------ !
implicit none
save

    real*8, dimension( 39 ) :: rsabc, x11as, x12as, y11as, y12as, y11bs, y12bs, &
                             x11cs, x12cs, y11cs, y12cs, drsabc, xy11as, xy12as, &
                             dxy11as, dxy12as, dy11as, dy12as, dy11bs, dy12bs, &
                             xy11cs, xy12cs, dxy11cs, dxy12cs, dy11cs, dy12cs

    real*8, dimension( 47 ) :: rsgh, x11gs, x12gs, y11gs, y12gs, y11hs, y12hs, &
                             rsm, xms, yms, zms, drsgh, dx11gs, dx12gs, dy11gs, &
                             dy12gs, dy11hs, dy12hs


end module lubdat
