subroutine out_init
! ------------------------------------------------------------ !
! This routine outputs information about the conditions of the
! simulation.
!
! Inputs: none
! Outputs: none
! Changes: none
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj


    if ( out_mode .eq. 1 ) write( 101, * ) t

    do ii = 1, np

        jj = 6 * ( ii - 1 )

        write( 101, 1004 ) x( jj + 1 ), x( jj + 2 ), &
                           x( jj + 3 ), x( jj + 4 ), x( jj + 5 ), x( jj + 6 )

    end do

    if ( out_mode .eq. 2 ) write( 101, * ) '---'

1004 format( 1x, 6( F13.6 ) ) 

end subroutine out_init


subroutine out_inter
! ------------------------------------------------------------ !
! This routine outputs data while the simulation is running 
!
! Inputs: none
! Outputs: none
! Changes: none
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj, kk
    real*8 :: elap, per, left

    if ( mod( int( t / dt_spec ), 500 ) .eq. 0 ) then

        call SYSTEM_CLOCK( COUNT = time_curr )
        time_curr = time_curr - time_init
        elap = dble( time_curr ) / dble( clock_rate )
        per = t / tf

        write( *, 1001 ) 100.d0 * per
        write( *, 1002 ) dble( time_curr ) / dble( clock_rate )
        
        if ( t .ne. 0.d0 ) then

            left = ( 1.d0 - per ) * elap / per
            write( *, 1003 ) left

        end if

        write( *, * )

    end if


    if ( mod( iloop, pos_time ) .eq. 0 ) then

        if ( out_mode .eq. 1 ) write( 101, * ) t

        do ii = 1, np

            jj = 6 * ( ii - 1 )

            write( 101, 1004 ) x( jj + 1 ), x( jj + 2 ), &
                               x( jj + 3 ), x( jj + 4 ), x( jj + 5 ), x( jj + 6 )

        end do

        if ( out_mode .eq. 2 ) write( 101, * ) '---'

    end if

1001 format( 1x, 'Percent complete:', F24.2 )
1002 format( 1x, 'Current run time (s):', F20.2 ) 
1003 format( 1x, 'Estimate time remaining (s): ', F12.2 )
1004 format( 1x, 6( F13.6 ) ) 
1005 format( 1x, 6( F14.7 ) ) 

end subroutine out_inter


subroutine out_final
! ------------------------------------------------------------ !
! This routine outputs information after the simulation has
! finished running
!
! Inputs: none
! Outputs: none
! Changes: none
! ------------------------------------------------------------ !
use global


end subroutine out_final
