program sdbd
! ------------------------------------------------------------ !
! Stokesian Dynamics - swimmers
! Copyright 2010, James W. Swan and John F. Brady
! 
! This program performs Stokesian Dynamics simulations of an
! a swimming body composed of an arbitrary number of spherical
! particles.  The initial configuration of the particles is set 
! in the file 'conf.in' with format:
!     x1 y1 z1
!     x2 y2 z2
!     ...
!     xn yn zn
! The parameters for the simulation are set in the file 'pref.in'.
! The simulation depends on tabulated pair-wise lubrication
! interactions stored in several files in the directory 'lubdat'.
! This program outputs the positions of the particles in the file
! 'pos.out'.  Any other output must be added by the user.
!
! To compile this program simply execute 'make Makefile'.  This
! make file is built to compile using the Intel Fortran Compiler
! for Linux 9.1.  Modifications may be necessary to ensure
! compilation of any other system.  A few calls are made to the
! system clock in 'init.f90' which may need to be changed if 
! compiling on another system.  These should be obvious, however. 
! ------------------------------------------------------------ !
use global

implicit none


    call init ! init.f90

    call out_init ! out.f90

    x0( : ) = x( : )


    do while ( t .lt. tf )

        call check_dist 
        ! Compute distance between particle pairs and check for overlaps.
        ! located in: 'dist.f90'

        call form_mob 
        ! Generate the mobility matrices.
        ! located in: 'mob.f90'

        call add_lub 
        ! Add the lubrication forces to the mobility inverse
        ! located in: 'lub.f90'

        call new_vel 
        ! Get new particle velocities
        ! located in: 'vel.f90'
        
        call trajectory 
        ! Advance the trajectories of the particles
        ! located in: 'traj.f90'

        t = t + dt_spec 
        iloop = iloop + 1

        call out_inter 
        ! Output simulation data for this time step
        ! located in: 'out.f90'

    end do


end ! sdbd
