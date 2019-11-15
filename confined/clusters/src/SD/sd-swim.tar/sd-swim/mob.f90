subroutine form_mob
! ------------------------------------------------------------ !
! This routine computes the far-field contributions to 
! the grand-mobility tensor which it then inverts using
! a cholesky decomposition.
!
! Inputs: none
! Outputs: none
! Changed: zmuf, zmus, zmes (far-field grand-resistance tensor)
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: idostep
    integer*4 :: ii, jj, kk
    integer*4 :: ph1, ph2, ph3, ph4, ph5, ph6

    zmuf( :, : ) = 0.d0
    zmus( :, : ) = 0.d0
    zmes( :, : ) = 0.d0


    ! Fill the self mobility terms
    call mobility( 0.d0, 0.d0, 0.d0, 0.d0, .true. )

    do kk = 1, np

        ph1 = 6 * ( kk - 1 )
        ph2 = ph1 + 3
        ph3 = 5 * ( kk - 1 )

        do ii = 1, 3

            zmuf( ph1 + ii, ph1 + ii ) = mob_a( ii, ii )
            zmuf( ph2 + ii, ph2 + ii ) = mob_c( ii, ii )

        end do


        do ii = 1, 5

            do jj = 1, 5

                zmes( ph3 + ii, ph3 + jj ) = mob_m( ii, jj )

            end do

        end do

    end do


    ! Fill the pair mobility terms
    do kk = 1, npair

        call mobility( pd( 5, kk ), pd( 1, kk ), pd( 2, kk ), pd( 3, kk ), .false. )

        ph1 = id( 1, kk )
        ph2 = id( 2, kk )

        ph5 = 5 * ( ph1 - 1 )
        ph6 = 5 * ( ph2 - 1 )

        ph1 = 6 * ( ph1 - 1 )
        ph2 = 6 * ( ph2 - 1 )

        ph3 = ph1 + 3
        ph4 = ph2 + 3

        do ii = 1, 3

            do jj = 1, 3

                zmuf( ph1 + ii, ph2 + jj ) = mob_a( ii, jj )
                zmuf( ph3 + ii, ph2 + jj ) = mob_b( ii, jj )
                zmuf( ph1 + ii, ph4 + jj ) = mob_bt( ii, jj )
                zmuf( ph3 + ii, ph4 + jj ) = mob_c( ii, jj )

                zmuf( ph2 + ii, ph1 + jj ) = mob_a( jj, ii )
                zmuf( ph4 + ii, ph1 + jj ) = mob_b( jj, ii ) 
                zmuf( ph2 + ii, ph3 + jj ) = mob_bt( jj, ii )
                zmuf( ph4 + ii, ph3 + jj ) = mob_c( jj, ii )

            end do

        
            do jj = 1, 5

                zmus( ph1 + ii, ph6 + jj ) = mob_gt( ii, jj )
                zmus( ph2 + ii, ph5 + jj ) = -1.d0 * mob_gt( ii, jj )

                zmus( ph3 + ii, ph6 + jj ) = mob_ht( ii, jj )
                zmus( ph4 + ii, ph5 + jj ) = mob_ht( ii, jj )

            end do

        end do


        do ii = 1, 5

            do jj = 1, 5

                zmes( ph5 + ii, ph6 + jj ) = mob_m( ii, jj )
                zmes( ph6 + ii, ph5 + jj ) = mob_m( jj, ii )

            end do

        end do

    end do

   
    ! Invert the grand-mobility tensor.  This is done in several steps
    ! which minimize the computation time.  The routine 'dgemm' is located
    ! in the file 'blas.f'.  The routine 'cholesky' is located in the file
    ! 'chol.f90'.
 
    ! Invert R1 = Muf ^ -1 => zmuf = zmuf ^ -1
    idostep = 3
    call cholesky( zmuf, t, n6, istop, idostep )
   
    if ( istop .eq. 0 ) then

        write( *, * ) 'Stopped after chol. 1 in mob.f90'

        stop

    end if


    if ( sim_mode .eq. 1 ) then

        ! Compute R2 = Mus(t) * R1 => rsu = zmus(t) * zmuf
         call dgemm( 't', 'n', n5, n6, n6, 1.d0, zmus, n6, zmuf, n6, 0.d0, rsu, n5 )

        ! Compute R3 = R2 * Mus - Mes => zmes = rsu * zmus - zmes
         call dgemm( 'n', 'n', n5, n5, n6, -1.d0, rsu, n5, zmus, n6, 1.d0, zmes, n5 )


        ! Invert  R4 = R3 ^ -1 => zmes = zmes ^ -1
         call cholesky( zmes, t, n5, istop, idostep )

        if ( istop .eq. 0 ) then

            write( *, * ) 'Stopped after chol. 2 in mob.f90'

            stop

        end if

        ! Compute R5 = -R3 * R4 => zmus = -rsu(t) * zmes
         call dgemm( 't', 'n', n6, n5, n5, -1.d0, rsu, n5, zmes, n5, 0.d0, zmus, n6 )
    
        ! Compute R6 = R1 - R5  => zmuf = zmuf - zmus * rsu
         call dgemm( 'n', 'n', n6, n6, n5, -1.d0, zmus, n6, rsu, n5, 1.d0, zmuf, n6 )

     end if

end subroutine form_mob


subroutine mobility( dr_inv, dx, dy, dz, self )
! ------------------------------------------------------------ !
! This routine computes the self and pair contributions
! to the far-field mobility tensor.
!
! Inputs: dr_inv - inverse distance between particle centers
!         dx, dy, dz - the x, y, z distance between part. cen.
!         self - if true, fill self tensor; if false fill pair.
! Outputs: none
! Changes: mob_a, mob_b, mob_c, mob_bt, mob_gt, mob_ht, mob_m
!          (These are the contributions to the grand mobility
!           tensors either self or pair)
! ------------------------------------------------------------ !
use global

implicit none

    real*8 :: dr_inv, dx, dy, dz
    real*8 :: dr_inv2, dr_inv3, dr_inv4, dr_inv5
    real*8 :: x12a, y12a, y12b, x12c, y12c, x12g
    real*8 :: y12g, y12h, x12m, y12m, z12m
    real*8, dimension( 3 ) :: e
    real*8, dimension( 3, 3 ) :: ee

    real*8, dimension( 3, 3, 3 ) :: gt, ht
    real*8, dimension( 3, 3, 3, 3 ) :: m

    integer*4 :: ii, jj, kk, ll, mm

    logical :: self


    if ( self ) then
    ! Determine the self contribution
    ! This is independent of dr_inv, dx, dy, dz

        mob_a( :, : ) = 0.d0
        mob_b( :, : ) = 0.d0
        mob_bt( :, : ) = 0.d0
        mob_c( :, : ) = 0.d0

        mob_gt( :, : ) = 0.d0
        mob_ht( :, : ) = 0.d0

        mob_m( :, : ) = 0.d0

        do ii = 1, 3

            mob_a( ii, ii ) = 1.d0

            mob_c( ii, ii ) = c3d4

        end do


        mob_m( 1, 5 ) = c9d10
        mob_m( 5, 1 ) = c9d10

        do ii = 1, 5

            mob_m( ii, ii ) = c9d5

        end do

    else
    ! Determine the pair contribution

        e( 1 ) = dx
        e( 2 ) = dy
        e( 3 ) = dz

        dr_inv2 = dr_inv * dr_inv
        dr_inv3 = dr_inv2 * dr_inv
        dr_inv4 = dr_inv3 * dr_inv
        dr_inv5 = dr_inv4 * dr_inv

        do ii = 1, 3

            ee( ii, ii ) = e( ii ) * e( ii )

            do jj = ii + 1, 3

                ee( ii, jj ) = e( ii ) * e( jj )
                ee( jj, ii ) = ee( ii, jj )

            end do

        end do

        x12a = c3d2 * dr_inv - dr_inv3
        y12a = c3d4 * dr_inv + c1d2 * dr_inv3

        y12b = -1.d0 * c3d4 * dr_inv2

        x12c = c3d4 * dr_inv3
        y12c = -1.d0 * c3d8 * dr_inv3

        x12g = c9d4 * dr_inv2 - c18d5 * dr_inv4
        y12g = c6d5 * dr_inv4

        y12h = -1.d0 * c9d8 * dr_inv3
        
        x12m = -1.d0 * c9d2 * dr_inv3 + c54d5 * dr_inv5
        y12m = c9d4 * dr_inv3 - c36d5 * dr_inv5
        z12m = c9d5 * dr_inv5

        do ii = 1, 3

            mob_a( ii, ii ) = x12a * ee( ii, ii ) + y12a * ( 1.d0 - ee( ii, ii ) )

            mob_c( ii, ii ) = x12c * ee( ii, ii ) + y12c * ( 1.d0 - ee( ii, ii ) )

            do jj = ii + 1, 3

                mob_a( ii, jj ) = x12a * ee( ii, jj ) - y12a * ee( ii, jj )
                mob_a( jj, ii ) = mob_a( ii, jj )

                mob_c( ii, jj ) = x12c * ee( ii, jj ) - y12c * ee( ii, jj )
                mob_c( jj, ii ) = mob_c( ii, jj )

            end do


            do jj = 1, 3

                kk = 6 - ii - jj

                if ( kk .eq. 0 ) kk = 3
                if ( kk .eq. 4 ) kk = 1

                mob_b( ii, jj ) = y12b * eps( ii, jj, kk ) * e( kk )
                mob_bt( jj, ii ) = -1.d0 * mob_b( ii, jj ) 

            end do

        end do


        do kk = 1, 3

            do ii = 1, 3

                do jj = 1, 3

                    ll = 6 - jj - kk
                    mm = 6 - ii - kk

                    if ( ll .eq. 0 ) ll = 3
                    if ( ll .eq. 4 ) ll = 1
                    if ( mm .eq. 0 ) mm = 3
                    if ( mm .eq. 4 ) mm = 1

                    gt( kk, ii, jj ) = -1.d0 * ( x12g * ( ee( ii, jj ) - c1d3 * delta( ii, jj ) ) &
                                     * e( kk ) + y12g * ( e( ii ) * delta( jj, kk ) + e( jj ) &
                                     * delta( ii, kk ) - 2.d0 * ee( ii, jj ) * e( kk ) ) )

                    ht( kk, ii, jj ) = y12h * ( ee( ii, ll ) * eps( jj, kk, ll ) &
                                     + ee( jj, mm ) * eps( ii, kk, mm ) ) 

                end do

            end do

        end do

       
        do ii = 1, 3

            do jj = 1, 3

                do kk = 1, 3

                    do ll = 1, 3

                m( ii, jj, kk, ll ) = c3d2 * x12m * ( ee( ii, jj ) - c1d3 * delta( ii, jj ) ) &
                        * ( ee( kk, ll ) - c1d3 * delta( kk, ll ) ) &
                        + c1d2 * y12m * ( ee( ii, kk ) * delta( jj, ll ) &
                        + ee( jj, kk ) * delta( ii, ll ) + ee( ii, ll ) * delta( jj, kk ) &
                        + ee( jj, ll ) * delta( ii, kk ) - 4.d0 * ee( ii, jj ) * ee( kk, ll ) ) &
                        + c1d2 * z12m * ( delta( ii, kk ) * delta( jj, ll ) &
                        + delta( jj, kk ) * delta( ii, ll ) - delta( ii, jj ) * delta( kk, ll ) &
                        + ee( ii, jj ) * delta( kk, ll ) + ee( kk, ll ) * delta( ii, jj ) &
                        - ee( ii, kk ) * delta( jj, ll ) - ee( jj, kk ) * delta( ii, ll ) &
                        - ee( ii, ll ) * delta( jj, kk ) - ee( jj, ll ) * delta( ii, kk ) &
                        + ee( ii, jj ) * ee( kk, ll ) )

                    end do

                end do

            end do

        end do

        ! This segment of code converts the pair contributions to the grand mobility tensor
        ! into a symmetric matrix.  Using the following conversions of the shear rate E and
        ! the stresslet S into the vectors EV and SV respectively.
        ! EV_1 = E_11 - E_33, EV_2 = 2 E_12, EV_3 = 2 E_13, EV_4 = 2 E_23, EV_5 = E_22 - E_33
        ! SV_1 = S_11, SV_2 = S_12 = S_21, SV_3 = S_13 = S_31, SV_4 = S_23 = S_32, SV_5 = S_22 
        do ii = 1, 3

            mob_gt( ii, 1 ) = gt( ii, 1, 1 ) - gt( ii, 3, 3 )
            mob_gt( ii, 2 ) = 2.d0 * gt( ii, 1, 2 )
            mob_gt( ii, 3 ) = 2.d0 * gt( ii, 1, 3 )
            mob_gt( ii, 4 ) = 2.d0 * gt( ii, 2, 3 )
            mob_gt( ii, 5 ) = gt( ii, 2, 2 ) - gt( ii, 3, 3 )

            mob_ht( ii, 1 ) = ht( ii, 1, 1 ) - ht( ii, 3, 3 )
            mob_ht( ii, 2 ) = 2.d0 * ht( ii, 1, 2 )
            mob_ht( ii, 3 ) = 2.d0 * ht( ii, 1, 3 )
            mob_ht( ii, 4 ) = 2.d0 * ht( ii, 2, 3 )
            mob_ht( ii, 5 ) = ht( ii, 2, 2 ) - ht( ii, 3, 3 )

        end do


        do ii = 1, 5

            if ( ( ii .eq. 1 ) .or. ( ii .eq. 5 ) ) then
               
                mob_m( ii, 1 ) = m( mesid( 1, ii ), mesid( 1, ii ), 1, 1 ) &
                               - m(  mesid( 1, ii ), mesid( 1, ii ), 3, 3 ) &
                               - ( m( mesid( 2, ii ), mesid( 2, ii ), 1, 1 ) &
                               - m(  mesid( 2, ii ), mesid( 2, ii ), 3, 3 ) )
                mob_m( ii, 2 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 1, 2 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 1, 2 ) ) 
                mob_m( ii, 3 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 1, 3 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 1, 3 ) )
                mob_m( ii, 4 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 2, 3 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 2, 3 ) )
                mob_m( ii, 5 ) = m( mesid( 1, ii ), mesid( 1, ii ), 2, 2 ) &
                               - m(  mesid( 1, ii ), mesid( 1, ii ), 3, 3 ) &
                               - ( m( mesid( 2, ii ), mesid( 2, ii ), 2, 2 ) &
                              - m(  mesid( 2, ii ), mesid( 2, ii ), 3, 3 ) )

            else

                mob_m( ii, 1 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 2, ii ) , 1, 1 ) &
                               - m( mesid( 1, ii ), mesid( 2, ii ), 3, 3 ) )
                mob_m( ii, 2 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 1, 2 ) 
                mob_m( ii, 3 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 1, 3 ) 
                mob_m( ii, 4 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 2, 3 ) 
                mob_m( ii, 5 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 2, ii ) , 2, 2 ) &
                               - m( mesid( 1, ii ), mesid( 2, ii ), 3, 3 ) )

            end if

        end do

    end if

end subroutine mobility
