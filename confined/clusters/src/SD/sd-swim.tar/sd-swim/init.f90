subroutine init
! ------------------------------------------------------------ !
! This routine controls to flow of the program initialization.
!
! Inputs: none
! Output: none
! Changes: none
! ------------------------------------------------------------ !
use global

implicit none

    call read_pref
    ! Read the preferences   
 
    call set_consts 
    ! Set any constants
    call alloc_arrays 
    ! Allocate all the memory
    call read_conf
    ! Read the configuration of particles

    call set_vars 
    ! Set the different variables

    call read_tab
    ! Read the lubrication tables

    call open_files
    ! Open any output files

end subroutine init


subroutine read_pref
! ------------------------------------------------------------ !
! This routine reads in all of the preferences
!
! Inputs: none
! Outputs: none
! Changes: np, fext_h, tf, dt_spec, sim_mode, istep
!          drset, dset, pos_time, out_mode
! ------------------------------------------------------------ !
use global

implicit none

    open( unit = 1001, file = 'pref.in', status = 'unknown' )

    call read_junk
    read( 1001, * ) np

    call read_junk
    read( 1001, * ) fext_h( 1 )
    read( 1001, * ) fext_h( 2 )
    read( 1001, * ) fext_h( 3 )
    read( 1001, * ) fext_h( 4 )
    read( 1001, * ) fext_h( 5 )
    read( 1001, * ) fext_h( 6 )

    call read_junk
    read( 1001, * ) tf
    read( 1001, * ) dt_spec

    sim_mode = 1
    istep = 1
    drset = 1.98d0
    dset = 1.0d-9

    read( 1001, * ) pos_time
    read( 1001, * ) out_mode

    close( 1001 )

    
    if ( pos_time > tf / dt ) then

        write( *, * ) 'pos_time was too large.  It was set equal to tf/dt so &
            that the particle positions are output at each time step.'

        pos_time = tf / dt

    end if


end subroutine read_pref


subroutine set_consts
! ------------------------------------------------------------ !
! This routine sets constants for use in the simulation.
!
! Inputs: none
! Outputs: none
! Changes: npair, n6, n5, n3
! ------------------------------------------------------------ !
use global

implicit none

    npair = np * ( np - 1 ) / 2.

    n6 = 6 * np
    n5 = 5 * np
    n3 = 3 * np

end subroutine set_consts


subroutine alloc_arrays
! ------------------------------------------------------------ !
! This routine allocates all the memory needed for the arrays
! used in the simulations
!
! Inputs: none
! Outputs: none
! Changes: none
! ------------------------------------------------------------ !
use global

implicit none

    allocate( x( n6 ) )
    allocate( x0( n6 ) )
    allocate( pd( 8, npair ) )

    allocate( u( n6 ) )
    allocate( udiff( n6 ) )
    allocate( uold( n6, 3 ) )

    allocate( uinf( n6 ) )
    allocate( einf( n5 ) )
    allocate( sh( n5 ) )

    allocate( fext( n6 ) )
    allocate( uext( n6 ) )
    
    allocate( id( 2, npair ) )

    allocate( rfu( n6, n6 ) )
    allocate( zmuf( n6, n6 ) )

    allocate( rfe( n6, n5 ) )
    allocate( zmus( n6, n5 ) )

    allocate( rsu( n5, n6 ) )

    allocate( rse( n5, n5 ) )
    allocate( zmes( n5, n5 ) )

end subroutine alloc_arrays


subroutine set_vars
! ------------------------------------------------------------ !
! This routine sets the initial values for a number of
! variables used in the simulation.
!
! Inputs: none
! Outputs: none
! Changes: t, iloop, dt, pedt, pertdt, dtginv, time_curr, 
!          time_init, clock_rate, dset, drmin, pos_time, 
!          delta, eps, mesid, id, zmuf, zmus, zmes, rfu, rsu, 
!          rfe, rse, uinf, einf, fext
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj, kk
    integer*4, dimension( 3 ) :: time

    ! Timing
    t = 0
    iloop = 1

    dt = dt_spec

    time_curr = 0
    call SYSTEM_CLOCK( COUNT_RATE = clock_rate )
    call SYSTEM_CLOCK( COUNT = time_init )

    dset = 2.d0 + dset
    drmin = 10.d0

    pos_time = ( tf / dt_spec ) / dble( pos_time ) 
    
    ! Precomputed matrices for mobility calculations
    delta( :, : ) = 0.d0
    delta( 1, 1 ) = 1.d0
    delta( 2, 2 ) = 1.d0
    delta( 3, 3 ) = 1.d0

    eps( :, :, : ) = 0.d0
    eps( 1, 2, 3 ) = 1.d0
    eps( 2, 3, 1 ) = 1.d0
    eps( 3, 1, 2 ) = 1.d0
    eps( 1, 3, 2 ) = -1.d0
    eps( 3, 2, 1 ) = -1.d0
    eps( 2, 1, 3 ) = -1.d0

    mesid( 1, 1 ) = 1
    mesid( 2, 1 ) = 3
    mesid( 1, 2 ) = 1
    mesid( 2, 2 ) = 2
    mesid( 1, 3 ) = 1
    mesid( 2, 3 ) = 3
    mesid( 1, 4 ) = 2
    mesid( 2, 4 ) = 3
    mesid( 1, 5 ) = 2
    mesid( 2, 5 ) = 3


    ! Set up particle pairs
    kk = 1

    do ii = 1, np

        do jj = ii + 1, np

            id( 1, kk ) = ii
            id( 2, kk ) = jj 
            kk = kk + 1

        end do

    end do



    uinf( : ) = 0.d0
    einf( : ) = 0.d0

    rfu( :, : ) = 0.d0
    zmuf( :, : ) = 0.d0

    rfe( :, : ) = 0.d0
    zmus( :, : ) = 0.d0

    rse( :, : ) = 0.d0
    zmes( :, : ) = 0.d0

    istop = 1


    do ii = 1, np

        jj = 6 * ( ii - 1 )

        do kk = 1, 6

            fext( jj + kk ) = fext_h( kk )

        end do

    end do

end subroutine set_vars


subroutine read_conf
! ------------------------------------------------------------ !
! This routine reads in the initial particle configuration.
!
! Inputs: none
! Outputs: none
! Changes: x
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, i6, jj
    
    open( unit = 1001, file = 'conf.in', status = 'old' )

    do ii = 1, np

        i6 = ii * 6

        read( 1001, * ) x( i6 - 5 ), x( i6 - 4 ), x( i6 - 3 )

        x( i6 - 2 ) = 0.0
        x( i6 - 1 ) = 0.0
        x( i6 ) = 0.0

    end do

    close( 1001 )

end subroutine read_conf


subroutine read_tab
! ------------------------------------------------------------ !
! This routine reads in the tabulated lubrication interactions.
!
! Inputs: none
! Outputs: none
! Changes: rsabc, x11as, x12as, y11as, y12as, y11bs, y12bs,
!          x11cs, x12cs, y11cs, y12cs, drsabc, xy11as, xy12as,
!          dxy11as, dxy12as, dy11as, dy12as, dy11bs, dy12bs,
!          xy11cs, xy12cs, dxy11cs, dxy12cs, dy11cs, dy12cs,
!          rsgh, x11gs, x12gs, y11gs, y12gs, y11hs, y12hs, rsm
!          xms, yms, zms, drsgh, dx11gs, dx12gs, dy11gs, dy12gs
!          dy11hs, dy12hs    
! ------------------------------------------------------------ !
use lubdat 

implicit none

    integer*4 :: ii

    open( unit = 1001, file = 'lubdat/r2babc.dat', status = 'old' ) ! 42
    open( unit = 1002, file = 'lubdat/r2bgh.dat', status = 'old' ) ! 43
    open( unit = 1003, file = 'lubdat/r2bm.dat', status = 'old' ) ! 44
    open( unit = 1004, file = 'lubdat/r2babcd.dat', status = 'old' ) ! 45
    open( unit = 1005, file = 'lubdat/r2bghd.dat', status = 'old' ) ! 46


    do ii = 1, 39
        read( 1001, * ) rsabc( ii ), x11as( ii ), x12as( ii ), y11as( ii ), &
                        y12as( ii ), y11bs( ii ), y12bs( ii ), x11cs( ii ), &
                        x12cs( ii ), y11cs( ii ), y12cs( ii )


        read( 1004, * ) drsabc( ii ), xy11as( ii ), xy12as( ii ), dxy11as( ii ), &
                        dxy12as( ii ), dy11as( ii ), dy12as( ii ), dy11bs( ii ), &
                        dy12bs( ii ), xy11cs( ii ), xy12cs( ii ), dxy11cs( ii ), &
                        dxy12cs( ii ), dy11cs( ii ), dy12cs( ii )

    end do


    do ii = 1, 47
        read( 1002, * ) rsgh( ii ), x11gs( ii ), x12gs( ii ), y11gs( ii ), y12gs( ii ), &
                        y11hs( ii ), y12hs( ii )

        read( 1003, * ) rsm( ii ), xms( ii ), yms( ii ), zms( ii )

        read( 1005, * ) drsgh( ii ), dx11gs( ii ), dx12gs( ii ), dy11gs( ii ), &
                        dy12gs( ii ), dy11hs( ii ), dy12hs( ii )
    end do


    close( 1001 )
    close( 1002 )
    close( 1003 )
    close( 1004 )
    close( 1005 )

end subroutine read_tab


subroutine read_junk

    read( 1001, * )
    read( 1001, * )
    read( 1001, * )

end subroutine read_junk


subroutine open_files
! ------------------------------------------------------------ !
! This routine opens the output files used in the simulation.
!
! Inputs: none
! Outputs: none
! Changes: none
! ------------------------------------------------------------ !
implicit none

    open( unit = 101, file = 'pos.out', status = 'unknown' )
    
end subroutine open_files
