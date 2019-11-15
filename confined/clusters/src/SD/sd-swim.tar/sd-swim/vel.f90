subroutine new_vel
! ------------------------------------------------------------ !
! This routine computes the contributions to the particle
! velocities due to a swimming gait and an external force/torque.
!
! Inputs: none
! Outputs: none
! Changes: u - the velocities due to the swimming gait 
!          sh - the hydrodynamic contribution to the stresslet
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj, jj2, kk, ll, mm, idostep
    real*8, dimension( 3 ) :: swimcm
    real*8, dimension( 6, n6 ) :: rbmconn
    real*8, dimension( 6, 6 ) :: rfu_swim
    real*8, dimension( n6, 6 ) :: swimhold
    real*8, dimension( 6 ) :: uswim, fswim
    real*8, dimension( n6 ) :: fconst, ugate, fgate
    real*8, dimension( 3, 3 ) :: ee, xf
    real*8 :: dx1, dx2, dx3, dissipation


    ! Find the center of mass of each swimmer
    swimcm( : ) = 0.d0

    do ii = 1, np

        jj = 6 * ii - 5

        swimcm( 1:3 ) = swimcm( 1:3 ) + x( jj:jj+2 )

    end do

    swimcm( : ) = swimcm( : ) / dble( np )
 

    ! Assemble the rigid body motion connectivity tensor (Sigma)
    rbmconn( :, : ) = 0.d0
   
    do ii = 1, np

        jj = 6 * ( ii - 1 )

        do kk = 1, 6

            rbmconn( kk, jj + kk ) = 1.d0

        end do


        dx1 = x( jj + 1 ) - swimcm( 1 )
        dx2 = x( jj + 2 ) - swimcm( 2 )
        dx3 = x( jj + 3 ) - swimcm( 3 )

        rbmconn( 4, jj + 2 ) = -dx3
        rbmconn( 4, jj + 3 ) = dx2
        rbmconn( 5, jj + 1 ) = dx3
        rbmconn( 5, jj + 3 ) = -dx1
        rbmconn( 6, jj + 1 ) = -dx2
        rbmconn( 6, jj + 2 ) = dx1

    end do    


    ! Calculate rfu_swim ( Sigma * RFU * Sigma^T)
    call dgemm( 'n', 't', n6, 6, n6, 1.d0, rfu, n6, rbmconn, &
        6, 0.d0, swimhold, n6 )

    call dgemm( 'n', 'n', 6, 6, n6, 1.d0, &
        rbmconn, 6, swimhold, n6, 0.d0, rfu_swim, 6 )


    ! Invert rfu_swim
    idostep = 3
    call cholesky( rfu_swim, t, 6, istop, idostep ) 

    if ( istop .eq. 0 ) then

        write( *, * ) 'Stopped after chol. 1 in vel.f90'

        stop

    end if


!!!!

! Beyond this point, all the resistance tensors and their inverses have been constructed
! one simply need compute the appropriate matrix-vector products (this is done with dgemv,
! from the blas package though one can hard-code this as well) to study the motion of
! colloids/rigid bodies/swimmers.  Note, the rate of strain (einf) has been set to zero.
! to study swimmers, this must be changed as suggested in the tutorial.


! At this point:
!
! rfu_swim = ( sigma * rfu * sigma^T ) ^ -1
!


    ! The implicit gate
    einf( : ) = 0.d0



    ! The explicit gate
    ugate( : ) = 0.d0



!!!!



    ! Begin calculating the velocities of the particles...


    uext( : ) = 0.d0
    fgate( : ) = 0.d0

    ! Calculate the propulsive thrust due to an explicit gate
    call dgemv( 'n', n6, n6, -1.d0, rfu, n6, ugate, 1, 0.d0, fgate, 1 )


    ! Calculate the propulsive thrust due to an implicit gate
    call dgemv( 'n', n6, n5, 1.d0, rfe, n6, einf, 1, 0.d0, uext, 1 )


    ! Add in any external forces
    uext( : ) = uext( : ) + fext( : ) + fgate( : )


    fswim( : ) = 0.d0


    ! Compute Sigma * u and store in fswim
    call dgemv( 'n', 6, n6, 1.d0, rbmconn, 6, uext, 1, 0.d0, fswim, 1 )

    uswim( : ) = 0.d0


    ! Compute rfu_swim * fswim
    call dgemv( 'n', 6, 6, 1.d0, rfu_swim, 6, fswim, 1, 0.d0, uswim, 1 )
 

    udiff( : ) = 0.d0


    ! Compute Sigma^T * uswim and store in udiff 
    call dgemv( 't', 6, n6, 1.d0, rbmconn, 6, uswim, 1, 0.d0, udiff, 1 )



    u( : ) = dt * ( udiff( : ) + ugate( : ) ) 


    ! Compute the constraining forces
    call dgemv( 'n', n6, n6, 1.d0, rfu, n6, udiff, 1, 0.d0, fconst, 1 )

    fconst( : ) = fconst( : ) + uext( : ) - 2.d0 * fext( : )

    sh( : ) = 0.d0


    ! Compute the hydrodynamic contribution to the stresslet
    call dgemv( 't', n6, n5, 1.d0, rfe, n6, udiff, 1, 0.d0, sh, 1 )

    call dgemv( 't', n6, n5, 1.d0, rfe, n6, ugate, 1, 1.d0, sh, 1 )

    call dgemv( 'n', n5, n5, -1.d0, rse, n5, einf, 1, 1.d0, sh, 1 )    


    dissipation = 0.d0

    do ii = 1, np

        jj = 5 * ii - 4
        kk = 6 * ii - 5

        ee( 1, 1 ) = 2.d0 / 3.d0 * ( einf( jj ) - 0.5d0 * einf( jj + 4 ) )
        ee( 2, 2 ) = 2.d0 / 3.d0 * ( einf( jj + 4 ) - 0.5d0 * einf( jj ) )
        ee( 3, 3 ) = -ee( 1, 1 ) - ee( 2, 2 )
        ee( 1, 2 ) = 0.5d0 * einf( jj + 1 )
        ee( 2, 1 ) = 0.5d0 * einf( jj + 1 )
        ee( 1, 3 ) = 0.5d0 * einf( jj + 2 )
        ee( 3, 1 ) = 0.5d0 * einf( jj + 2 )
        ee( 2, 3 ) = 0.5d0 * einf( jj + 3 )
        ee( 3, 2 ) = 0.5d0 * einf( jj + 3 )

        xf( 1, 1 ) = -x( kk ) * fconst( kk ) + sh( jj )
        xf( 2, 2 ) = -x( kk + 1 ) * fconst( kk + 1 ) + sh( jj + 4 )
        xf( 3, 3 ) = -x( kk + 2 ) * fconst( kk + 2 ) - sh( jj ) - sh( jj + 4 )
        xf( 1, 2 ) = -x( kk ) * fconst( kk + 1 ) + sh( jj + 1 )
        xf( 2, 1 ) = -x( kk + 1 ) * fconst( kk ) + sh( jj + 1 )
        xf( 1, 3 ) = -x( kk ) * fconst( kk + 2 ) + sh( jj + 2 )
        xf( 3, 1 ) = -x( kk + 2 ) * fconst( kk ) + sh( jj + 2 )
        xf( 2, 3 ) = -x( kk + 1 ) * fconst( kk + 2 ) + sh( jj + 3 )
        xf( 3, 2 ) = -x( kk + 2 ) * fconst( kk + 1 ) + sh( jj + 3 )

        dissipation = dissipation + dot_product( ee( 1, : ), xf( 1, : ) ) + dot_product( ee( 2, : ), xf( 2, : ) ) &
            + dot_product( ee( 3, : ), xf( 3, : ) )

    end do


    x0( : ) = x( : )

end subroutine new_vel
