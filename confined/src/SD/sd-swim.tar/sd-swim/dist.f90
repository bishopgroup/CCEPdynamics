subroutine check_dist
! ------------------------------------------------------------ !
! This routine computes the distances between all the particle
! pair centers and the in the case of particles with centers
! closer than 4 radii, it also computes some relavant
! distances for use in the lubrication calculations.
!
! Inputs: none
! Outputs: none
! Changes: pd - contains the x, y, z, total distance and
!               the inverse distance between the particle
!               pair centers.
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj, kk
    real*8 :: dx, dy, dz, dr, dr_inv, dst, dst_inv

    do ii = 1, npair

        kk = 6 * id( 1, ii ) - 5
        jj = 6 * id( 2, ii ) - 5

        dx = x( jj ) - x( kk )
        dy = x( jj + 1 ) - x( kk + 1 ) 
        dz = x( jj + 2 ) - x( kk + 2 ) 
        dr = dsqrt( dx * dx + dy * dy + dz * dz )

        
        if ( dr .lt. drset ) then

            write( *, 100 ) ( jj + 5 ) / 6, ( kk + 5 ) / 6, dr, t
            istop = 0 
           
            stop 

        end if

        if ( dr .lt. dset ) then

            if ( dr .lt. drmin ) drmin = dr
            dr = dset

        end if

        dr_inv = 1.d0 / dr

        pd( 1, ii ) = dx * dr_inv
        pd( 2, ii ) = dy * dr_inv
        pd( 3, ii ) = dz * dr_inv
        pd( 4, ii ) = dr
        pd( 5, ii ) = dr_inv

        if ( dr .lt. 4.d0 ) then

            dst = dr - 2.d0
            dst_inv = 1.d0 / dst
            
            pd( 6, ii ) = dst
            pd( 7, ii ) = dst_inv
            pd( 8, ii ) = dlog( dst_inv )

        end if

    end do 

100 format( 1x, 'Particles overlapped: ', I4, I4, F7.4, F7.3 )

end subroutine check_dist
