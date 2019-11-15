subroutine cholesky( aa, rtime, ndim, istop, idostep )
! ------------------------------------------------------------ !
! This routine computes the cholesky decomposition and the
! inversion of a symmetric positive definite matrix.
!
! Inputs: aa - the square matrix to be inverted
!         rtime - a time index for the inversion call
!         ndim - the length or width of aa
!         istop - a variable for reporting errors
!         idostep - specifies which step of the inversion
!             process to stop after.  If (1), perform partial
!             decomposition.  If (2), perform full choleksy
!             decomposition.  If (3), perform full inversion.
! Outputs: none
! Changes: aa - contains the result of the decomp./inversion
!          istop - if there was an error, istop = 0.
! ------------------------------------------------------------ !
implicit none

    integer*4 :: ndim, istop, idostep
    integer*4 :: ii, jj, kk, km1, kp1
    real*8, dimension( ndim, ndim ) :: aa
    real*8 :: rtime, sum, test, tempo, aiinv


    do kk = 1, ndim

        kp1 = kk + 1
        km1 = kk - 1

        sum = 0.d0

        do ii = 1, km1

            sum = sum + aa( ii, kk ) * aa( ii, kk )

        end do 

        test = aa( kk, kk ) - sum

        if ( test .le. 0.d0 ) then

            write( *, 100 )  rtime, kk, test

            istop = 0

            return

        end if

        tempo = 1.d0 / dsqrt( test )

        aa( kk, kk ) = tempo

        do ii = kp1, ndim

            sum = 0.d0

            do jj = 1, km1
     
                sum = sum + aa( jj, ii ) * aa( jj, kk )

            end do

            aa( kk, ii ) = ( aa( kk, ii ) - sum ) * tempo

        end do

    end do


    if ( idostep .eq. 1 ) return


    do kk = 2, ndim

        km1 = kk - 1

        aiinv = aa( kk, kk )

        do ii = 1, km1

            sum = 0.d0

            do jj = ii, km1

                sum = sum - aa( jj, kk ) * aa( jj, ii )

            end do

            aa( kk, ii ) = aiinv * sum

        end do

    end do 


    if ( idostep .eq. 2 ) return


    do kk = 1, ndim

        do ii = 1, kk

            sum = 0.d0

            do jj = kk, ndim

                sum = sum + aa( jj, ii ) * aa( jj, kk )

            end do

            aa( ii, kk ) = sum

        end do

    end do


    do ii = 1, ndim

        do jj = 1, ii - 1

            aa( ii, jj ) = aa( jj, ii )

        end do

    end do

100 format(1x, 'Not sym. pos. def. matrix: ', F6.3, I3, F8.3 )

end subroutine cholesky
