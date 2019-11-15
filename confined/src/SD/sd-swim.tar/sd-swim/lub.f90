subroutine add_lub
! ------------------------------------------------------------ !
! This routine uses the pair-wise lubrication interactions
! to compute the many-bodied resistance matrix
!
! Input: none
! Output: none
! Changes: rfu, rfe, rse - the full resistance tensors
!              many-bodied far-field and lubrication included
! ------------------------------------------------------------ !
use global

implicit none

    real*8, dimension( 3 ) :: d
    real*8, dimension( 12, 12 ) :: tabc
    real*8, dimension( 12, 10 ) :: tght
    real*8, dimension( 10, 10 ) :: tzm

    real*8 :: dr, dst, dst_inv, dst_log

    integer*4 :: ii, jj, kk, im1, ira, irg, irm, icg, jm1, jca, jrg, jcm, jcg
    integer*4 :: il, i1, i2, jl, j1, j2, ir, jc, mm


    do ii = 1, n6

        do jj = 1, ii

            rfu( jj, ii ) = zmuf( jj, ii )

        end do

    end do


    do ii = 1, n5

        do jj = 1, n6

            rfe( jj, ii ) = zmus( jj, ii )

        end do

    end do


    do ii = 1, n5

        do jj = 1, ii

            rse( jj, ii ) = zmes( jj, ii )

        end do

    end do 


    do kk = 1, npair

        dr = pd( 4, kk )

        if ( dr .lt. 4.d0 ) then

            ii = id( 1, kk )
            jj = id( 2, kk )

            im1 = ii - 1
            ira = im1 * 6
            irg = ira
            irm = im1 * 5
            icg = irm

            jm1 = jj - 1
            jca = jm1 * 6
            jrg = jca
            jcm = jm1 * 5
            jcg = jcm

            do mm = 1, 3

                d( mm ) = pd( mm, kk )

            end do

            dst = pd( 6, kk )
            dst_inv = pd( 7, kk )
            dst_log = pd( 8, kk )

            call calc_lub( dr, d, dst, dst_inv, dst_log, tabc, tght, tzm )

 
            do jc = 1, 6
            
                jl = jc + 6
                j1 = ira + jc
                j2 = jca + jc

                do ir = 1, jc
              
                    il = ir + 6
                    i1 = ira + ir
                    i2 = jca + ir

                    rfu( i1, j1 ) = rfu( i1, j1 ) + tabc( ir, jc )
                    rfu( i2, j2 ) = rfu( i2, j2 ) + tabc( il, jl )

                end do

            end do

            
            do jc = 7, 12

                j1 = jca + jc - 6

                do ir = 1, 6

                    i1 = ira + ir

                    rfu( i1, j1 ) = rfu( i1, j1 ) + tabc( ir, jc )

                end do

            end do


            if ( sim_mode .eq. 1 ) then
    
                do jc = 1, 5

                    jl = jc + 5
                    j1 = icg + jc
                    j2 = jcg + jc

                    do ir = 1, 6

                        il = ir + 6
                        i1 = irg + ir
                        i2 = jrg + ir

                        rfe( i1, j1 ) = rfe( i1, j1 ) + tght( ir, jc )
                        rfe( i2, j2 ) = rfe( i2, j2 ) + tght( il, jl )
                        rfe( i1, j2 ) = rfe( i1, j2 ) + tght( ir, jl )
                        rfe( i2, j1 ) = rfe( i2, j1 ) + tght( il, jc )
       
                    end do

                end do 


                do jc = 1, 5

                    jl = jc + 5
                    j1 = irm + jc
                    j2 = jcm + jc

                    do ir = 1, jc

                        il = ir + 5
                        i1 = irm + ir
                        i2 = jcm + ir

                        rse( i1, j1 ) = rse( i1, j1 ) + tzm( ir, jc )
                        rse( i2, j2 ) = rse( i2, j2 ) + tzm( il, jl )

                    end do

                end do


                do jc = 6, 10
 
                    j1 = jcm + jc - 5

                    do ir = 1, 5

                        i1 = irm + ir
 
                        rse( i1, j1 ) = rse( i1, j1 ) + tzm( ir, jc )

                    end do

                end do    

            end if

        end if

    end do

    
    do ii = 1, n6

        do jj = 1, ii - 1

            rfu( ii, jj ) = rfu( jj, ii )

        end do

    end do


    do ii = 1, n5

        do jj = 1, ii - 1

            rse( ii, jj ) = rse( jj, ii )

        end do

    end do

end subroutine add_lub


subroutine calc_lub( dr, d, dst, dst_inv, dst_log, tabc, tght, tzm )
! ------------------------------------------------------------ !
! This routine computes the pair-wise lubrication interactions
! between particle pairs.
!
! Inputs: dr - the distance between centers
!         d - the center to center vector
!         dst - the distance between particle surfaces
!         dst_inv - 1 / dst
!         dst_log - log( dst_inv )
! Outputs: none
! Changes: tabc, tght, tzm contain the pair-wise lubrication
!               interactions.
! ------------------------------------------------------------ !
use global
use lubdat

implicit none

    real*8 :: dr, dst, dst_inv, dst_log
    real*8, dimension( 3 ) :: d
    real*8, dimension( 12, 12 ) :: tabc
    real*8, dimension( 12, 10 ) :: tght
    real*8, dimension( 10, 10 ) :: tzm
    
    real*8, dimension( 3, 3 ) :: ee
    real*8, dimension( 3, 3, 3, 3 ) :: m
    real*8, dimension( 3, 3, 3 ) :: gt, ht

    real*8 :: xi, xi1, dlx, xdlx, dlx1, csa1, csa2, csa3, csa4, csa5
    real*8 :: x11a, x12a, y11a, y12a, y11b, y12b, x11c, x12c, y11c, y12c
    real*8 :: csg1, csg2, x11g, x12g, y11g, y12g, y11h, y12h
    real*8 :: xm, ym, zm, cm, c1
    real*8 :: xmy11a, xmy12a, xmy11c, xmy12c

    integer*4 :: ida, ia, ib, j3, j6, j9, j5, i3, i6, i9, i5, cgh
    integer*4 :: ii, jj, kk, ll, mm

    real*8 :: cx13x11g, c2y11g, xm2y11g, comd11, comd22, comd33, c2ymx11, cun35, cun56
    real*8 :: cun712, cun1011, c13x12g, c2y12g, xm2y12g, cumd11, cumd22, cumd33, cun34
    real*8 :: cyhd12b, y12hd23, y12hd13, y12hd12, cyhd12a, y11hd23, y11hd13, y11hd12
    real*8 :: d33md11, d22md33, d11md22, con35, con56, con712, con1011, cun89, con89, con34
    real*8 :: c13x11g, c2ymx12

!***************************************************************************c
!***************************************************************************c

    if ( dr .le. 2.1d0 ) then

        xi   = dst

        xi1  = dst_inv 
        dlx  = dst_log 

        xdlx = xi * dlx
        dlx1 = dlx + xdlx

        csa1 = dlx * one6
        csa2 = xdlx * one6
        csa3 = dlx1 * one6
        csa4 = 0.25d0 * xi1 + 0.225d0 * dlx
        csa5 = dlx * one15

!*** a, btilda, and c terms for rfu.

        x11a = csa4 - 1.23041d0 + c3d112 * xdlx + 1.8918d0 * xi
        x12a = -x11a + 0.00312d0 - 0.0011d0 * xi
        y11a = csa1 - 0.39394d0 + 0.95665d0 * xi
        y12a = -y11a + 0.00463606d0 - 0.007049d0 * xi

        y11b = -csa1 + 0.408286d0 - xdlx * one12 - 0.84055d0 * xi
        y12b = -y11b + 0.00230818d0 - 0.007508d0 * xi

        x11c = 0.0479d0 - csa2 + 0.12494d0 * xi
        x12c = -0.031031d0 + csa2 - 0.174476d0 * xi
        y11c = 4.d0 * csa5 - 0.605434d0 + ug1 * xdlx + 0.939139d0 * xi
        y12c = csa5 - 0.212032d0 + ug2 * xdlx + 0.452843d0 * xi

!*** g and h terms for rsu.

        csg1 = csa4 + ug4 * xdlx
        csg2 = dlx * one12 + xdlx * one24

        x11g = csg1 - 1.16897d0 + 1.47882d0 * xi
        x12g = -csg1 + 1.178967d0 - 1.480493d0 * xi
        y11g = csg2 - 0.2041d0 + 0.442226d0 * xi
        y12g = -csg2 + 0.216365d0 - 0.469830 * xi

        y11h = 0.5d0 * csa5 - 0.143777d0 + ug7 * xdlx + 0.264207d0 * xi
        y12h = 2.d0 * csa5 - 0.298166d0 + ug5 * xdlx + 0.534123d0 * xi

!*** m term for rse.

        xm   = one3 * xi1 + 0.3d0 * dlx - 1.48163d0 + 0.335714d0 * xdlx + 1.413604d0 * xi
        ym   = csa3 - 0.423489d0 + 0.827286d0 * xi
        zm   = 0.0129151d0 - 0.042284d0 * xi

    else
           
        ida  =  idint (20.d0*dst)
        ib   = -1 + ida
        ia   =  ib + 1
           
        c1   = (dr-rsabc(ib))/(rsabc(ia)-rsabc(ib))

        x11a = (x11as(ia)-x11as(ib))*c1+x11as(ib)
        x12a = (x12as(ia)-x12as(ib))*c1+x12as(ib)
        y11a = (y11as(ia)-y11as(ib))*c1+y11as(ib)
        y12a = (y12as(ia)-y12as(ib))*c1+y12as(ib)

        y11b = (y11bs(ia)-y11bs(ib))*c1+y11bs(ib)
        y12b = (y12bs(ia)-y12bs(ib))*c1+y12bs(ib)

        y11c = (y11cs(ia)-y11cs(ib))*c1+y11cs(ib)
        y12c = (y12cs(ia)-y12cs(ib))*c1+y12cs(ib)
        x11c = (x11cs(ia)-x11cs(ib))*c1+x11cs(ib)
        x12c = (x12cs(ia)-x12cs(ib))*c1+x12cs(ib)

  
        if ( dr .lt. 2.2d0 ) then
            
            ib = -9 + idint(100.d0*dst)
            
        else

            ib = 7 + ida

        endif

        ia = ib + 1

        cgh  = (dr-rsgh(ib))/(rsgh(ia)-rsgh(ib))

        x11g = (x11gs(ia)-x11gs(ib))*cgh + x11gs(ib)
        x12g = (x12gs(ia)-x12gs(ib))*cgh + x12gs(ib)
        y11g = (y11gs(ia)-y11gs(ib))*cgh + y11gs(ib)
        y12g = (y12gs(ia)-y12gs(ib))*cgh + y12gs(ib)

        y11h = (y11hs(ia)-y11hs(ib))*cgh + y11hs(ib)
        y12h = (y12hs(ia)-y12hs(ib))*cgh + y12hs(ib)

        cm   = (dr-rsm(ib))/(rsm(ia)-rsm(ib))

        xm   = (xms(ia)-xms(ib))*cm + xms(ib)
        ym   = (yms(ia)-yms(ib))*cm + yms(ib)
        zm   = (zms(ia)-zms(ib))*cm + zms(ib)

    endif
    
!***************************************************************************c
!***************************************************************************c

    do ii = 1, 3

        ee( ii, ii ) = d( ii ) * d( ii )

        do jj = ii + 1, 3

            ee( ii, jj ) = d( ii ) * d( jj )
            ee( jj, ii ) = ee( ii, jj )

        end do

    end do

!***************************************************************************c
!***************************************************************************c
!****************************** form tabc for rfu **************************c

    xmy11a = x11a - y11a
    xmy12a = x12a - y12a
    xmy11c = x11c - y11c
    xmy12c = x12c - y12c

!*** insert upper half of a11.

    tabc( 1, 1 ) = xmy11a * ee( 1, 1 ) + y11a
    tabc( 2, 2 ) = xmy11a * ee( 2, 2 ) + y11a
    tabc( 3, 3 ) = xmy11a * ee( 3, 3 ) + y11a
    tabc( 1, 2 ) = xmy11a * ee( 1, 2 )
    tabc( 1, 3 ) = xmy11a * ee( 1, 3 )
    tabc( 2, 3 ) = xmy11a * ee( 2, 3 )

!*** insert a12.

    tabc( 1, 7 ) = xmy12a * ee( 1, 1 ) + y12a
    tabc( 2, 8 ) = xmy12a * ee( 2, 2 ) + y12a
    tabc( 3, 9 ) = xmy12a * ee( 3, 3 ) + y12a
    tabc( 1, 8 ) = xmy12a * ee( 1, 2 )
    tabc( 1, 9 ) = xmy12a * ee( 1, 3 )
    tabc( 2, 9 ) = xmy12a * ee( 2, 3 )
    tabc( 2, 7 ) = tabc( 1, 8 )
    tabc( 3, 7 ) = tabc( 1, 9 )
    tabc( 3, 8 ) = tabc( 2, 9 )

!*** insert upper half of c11.

    tabc( 4, 4 ) = xmy11c * ee( 1, 1 ) + y11c
    tabc( 5, 5 ) = xmy11c * ee( 2, 2 ) + y11c
    tabc( 6, 6 ) = xmy11c * ee( 3, 3 ) + y11c
    tabc( 4, 5 ) = xmy11c * ee( 1, 2 )
    tabc( 4, 6 ) = xmy11c * ee( 1, 3 )
    tabc( 5, 6 ) = xmy11c * ee( 2, 3 )

!*** insert c12.

    tabc( 4, 10 ) = xmy12c * ee( 1, 1 ) + y12c
    tabc( 5, 11 ) = xmy12c * ee( 2, 2 ) + y12c
    tabc( 6, 12 ) = xmy12c * ee( 3, 3 ) + y12c
    tabc( 4, 11 ) = xmy12c * ee( 1, 2 )
    tabc( 4, 12 ) = xmy12c * ee( 1, 3 )
    tabc( 5, 12 ) = xmy12c * ee( 2, 3 )
    tabc( 5, 10 ) = tabc( 4, 11 )
    tabc( 6, 10 ) = tabc( 4, 12 )
    tabc( 6, 11 ) = tabc( 5, 12 )

!*** fill in upper half of a22 (=a11).

    tabc( 7, 7 ) = tabc( 1, 1 )
    tabc( 7, 8 ) = tabc( 1, 2 )
    tabc( 7, 9 ) = tabc( 1, 3 )
    tabc( 8, 8 ) = tabc( 2, 2 )
    tabc( 8, 9 ) = tabc( 2, 3 )
    tabc( 9, 9 ) = tabc( 3, 3 )

!*** fill in upper half of c22 (=c11).

    tabc( 10, 10 ) = tabc( 4, 4 )
    tabc( 10, 11 ) = tabc( 4, 5 )
    tabc( 10, 12 ) = tabc( 4, 6 )
    tabc( 11, 11 ) = tabc( 5, 5 )
    tabc( 11, 12 ) = tabc( 5, 6 )
    tabc( 12, 12 ) = tabc( 6, 6 )

!*** insert bt11.

    tabc( 1, 4 ) = 0.d0
    tabc( 1, 5 ) = -y11b * d( 3 )
    tabc( 1, 6 ) = y11b * d( 2 )
    tabc( 2, 5 ) = 0.d0
    tabc( 2, 6 ) = -y11b * d( 1 )
    tabc( 2, 4 ) = -tabc( 1, 5 )
    tabc( 3, 4 ) = -tabc( 1, 6 )
    tabc( 3, 5 ) = -tabc( 2, 6 )
    tabc( 3, 6 ) = 0.d0

!*** insert bt12.

    tabc( 1, 10 ) = 0.d0
    tabc( 1, 11 ) = y12b * d( 3 )
    tabc( 1, 12 ) = -y12b * d( 2 )
    tabc( 2, 11 ) = 0.d0
    tabc( 2, 12 ) = y12b * d( 1 )
    tabc( 2, 10 ) = -tabc( 1, 11 )
    tabc( 3, 10 ) = -tabc( 1, 12 )
    tabc( 3, 11 ) = -tabc( 2, 12 )
    tabc( 3, 12 ) = 0.d0

!***************************************************************************c
!***************************************************************************c
!*** fill in bt22 (=-bt11) and b12 (=bt12).

    do j3 = 4, 6

       j6 = j3 + 3
       j9 = j3 + 6

       do ii = 1, 3

          i3 = ii + 3
          i6 = ii + 6

          tabc( i3, j6 ) = tabc( ii, j9 )
          tabc( i6, j9 ) = -tabc( ii, j3 )

        end do

    end do

!***************************************************************************c
!***************************************************************************c
!****************************** form tght for rfe **************************c
!*** insert gt11.

    if ( sim_mode .eq. 1 ) then

        do kk = 1, 3

            do ii = 1, 2

                do jj = ii, 3

                    ll = 6 - jj - kk
                    mm = 6 - ii - kk

                    if ( ll .eq. 0 ) ll = 3
                    if ( ll .eq. 4 ) ll = 1
                    if ( mm .eq. 0 ) mm = 3
                    if ( mm .eq. 4 ) mm = 1

                    gt( kk, ii, jj ) = x11g * ( ee( ii, jj ) - c1d3 * delta( ii, jj ) ) &
                                     * d( kk ) + y11g * ( d( ii ) * delta( jj, kk ) + d( jj ) &
                                     * delta( ii, kk ) - 2.d0 * ee( ii, jj ) * d( kk ) )

                    ht( kk, ii, jj ) = y11h * ( ee( ii, ll ) * eps( jj, kk, ll ) &
                                     + ee( jj, mm ) * eps( ii, kk, mm ) )

                end do

            end do

        end do


        do ii = 1, 3

            jj = ii + 3

            tght( ii, 1 ) = gt( ii, 1, 1 ) - gt( ii, 3, 3 )
            tght( ii, 2 ) = 2.d0 * gt( ii, 1, 2 )
            tght( ii, 3 ) = 2.d0 * gt( ii, 1, 3 )
            tght( ii, 4 ) = 2.d0 * gt( ii, 2, 3 )
            tght( ii, 5 ) = gt( ii, 2, 2 ) - gt( ii, 3, 3 )


            tght( jj, 1 ) = ht( ii, 1, 1 ) - ht( ii, 3, 3 )
            tght( jj, 2 ) = 2.d0 * ht( ii, 1, 2 )
            tght( jj, 3 ) = 2.d0 * ht( ii, 1, 3 )
            tght( jj, 4 ) = 2.d0 * ht( ii, 2, 3 )
            tght( jj, 5 ) = ht( ii, 2, 2 ) - ht( ii, 3, 3 )

        end do


        do kk = 1, 3

            do ii = 1, 2

                do jj = ii, 3

                    ll = 6 - jj - kk
                    mm = 6 - ii - kk

                    if ( ll .eq. 0 ) ll = 3
                    if ( ll .eq. 4 ) ll = 1
                    if ( mm .eq. 0 ) mm = 3
                    if ( mm .eq. 4 ) mm = 1

                    gt( kk, ii, jj ) = x12g * ( ee( ii, jj ) - c1d3 * delta( ii, jj ) ) &
                                     * d( kk ) + y12g * ( d( ii ) * delta( jj, kk ) + d( jj ) &
                                     * delta( ii, kk ) - 2.d0 * ee( ii, jj ) * d( kk ) )

                    ht( kk, ii, jj ) = y12h * ( ee( ii, ll ) * eps( jj, kk, ll ) &
                                     + ee( jj, mm ) * eps( ii, kk, mm ) )

                end do

            end do

        end do


        do ii = 1, 3

            kk = ii + 6
            jj = ii + 3

            tght( kk, 1 ) = gt( ii, 1, 1 ) - gt( ii, 3, 3 )
            tght( kk, 2 ) = 2.d0 * gt( ii, 1, 2 )
            tght( kk, 3 ) = 2.d0 * gt( ii, 1, 3 )
            tght( kk, 4 ) = 2.d0  *gt( ii, 2, 3 )
            tght( kk, 5 ) = gt( ii, 2, 2 ) - gt( ii, 3, 3 )


            tght( jj, 6 ) = ht( ii, 1, 1 ) - ht( ii, 3, 3 )
            tght( jj, 7 ) = 2.d0 * ht( ii, 1, 2 )
            tght( jj, 8 ) = 2.d0 * ht( ii, 1, 3 )
            tght( jj, 9 ) = 2.d0 * ht( ii, 2, 3 )
            tght( jj, 10 ) = ht( ii, 2, 2 ) - ht( ii, 3, 3 )

        end do

        c13x11g   = one3*x11g
        c2y11g    = 2.d0*y11g
        xm2y11g   = x11g - c2y11g
        comd11    = ee( 1, 1 )*xm2y11g
        comd22    = ee( 2, 2 )*xm2y11g
        comd33    = ee( 3, 3 )*xm2y11g
        c2ymx11   = c2y11g - c13x11g
        con34     = comd11 - c13x11g
        con56     = comd11 + y11g
        con712    = comd22 + y11g
        con89     = comd33 + y11g
        con1011   = comd22 - c13x11g

        tght(1,1) = d(1)*(comd11+c2ymx11)
        tght(1,2) = d(2)*con56
        tght(1,3) = d(3)*con56
        tght(1,4) = d(1)*ee(2,3)*xm2y11g
        tght(1,5) = d(1)*con1011
        tght(2,1) = d(2)*con34
        tght(2,2) = d(1)*con712
        tght(2,3) = tght(1,4)
        tght(2,4) = d(3)*con712
        tght(2,5) = d(2)*(comd22+c2ymx11)
        tght(3,1) = d(3)*con34
        tght(3,2) = tght(1,4)
        tght(3,3) = d(1)*con89
        tght(3,4) = d(2)*con89
        tght(3,5) = d(3)*con1011

!*** insert gt21.

        c13x12g   = one3*x12g
        c2y12g    = 2.d0*y12g
        xm2y12g   = x12g - c2y12g
        cumd11    = ee(1,1)*xm2y12g
        cumd22    = ee(2,2)*xm2y12g
        cumd33    = ee(3,3)*xm2y12g
        c2ymx12   = c2y12g - c13x12g
        cun34     = cumd11 - c13x12g
        cun56     = cumd11 + y12g
        cun712    = cumd22 + y12g
        cun89     = cumd33 + y12g
        cun1011   = cumd22 - c13x12g

        tght(7,1) = d(1)*(cumd11+c2ymx12)
        tght(7,2) = d(2)*cun56
        tght(7,3) = d(3)*cun56
        tght(7,4) = d(1)*ee(2,3)*xm2y12g
        tght(7,5) = d(1)*cun1011
        tght(8,1) = d(2)*cun34
        tght(8,2) = d(1)*cun712
        tght(8,3) = tght(7,4)
        tght(8,4) = d(3)*cun712
        tght(8,5) = d(2)*(cumd22+c2ymx12)
        tght(9,1) = d(3)*cun34
        tght(9,2) = tght(7,4)
        tght(9,3) = d(1)*cun89
        tght(9,4) = d(2)*cun89
        tght(9,5) = d(3)*cun1011

!*** insert ht11.

        d11md22   =  ee(1,1) - ee(2,2)
        d22md33   =  ee(2,2) - ee(3,3)
        d33md11   =  ee(3,3) - ee(1,1)
        y11hd12   =  y11h*ee(1,2)
        y11hd13   =  y11h*ee(1,3)
        y11hd23   =  y11h*ee(2,3)
        cyhd12a   =  2.d0*y11hd12

        tght(4,1) =  0.d0
        tght(4,2) = -y11hd13
        tght(4,3) =  y11hd12
        tght(4,4) =  y11h*d22md33
        tght(4,5) = -2.d0*y11hd23
        tght(5,1) =  2.d0*y11hd13
        tght(5,2) =  y11hd23
        tght(5,3) =  y11h*d33md11
        tght(5,4) = -y11hd12
        tght(5,5) =  0.d0
        tght(6,1) = -cyhd12a
        tght(6,2) =  y11h*d11md22
        tght(6,3) = -y11hd23
        tght(6,4) =  y11hd13
        tght(6,5) =  cyhd12a

!*** insert ht12.

        y12hd12    =  y12h*ee(1,2)
        y12hd13    =  y12h*ee(1,3)
        y12hd23    =  y12h*ee(2,3)
        cyhd12b    =  2.d0*y12hd12
        tght(4,6)  =  0.d0
        tght(4,7)  = -y12h*ee(1,3)
        tght(4,8)  =  y12h*ee(1,2)
        tght(4,9)  =  y12h*d22md33
        tght(4,10) = -2.d0*y12hd23
        tght(5,6)  =  2.d0*y12hd13
        tght(5,7)  =  y12hd23
        tght(5,8)  =  y12h*d33md11
        tght(5,9)  = -y12hd12
        tght(5,10) =  0.d0
        tght(6,6)  = -cyhd12b
        tght(6,7)  =  y12h*d11md22
        tght(6,8)  = -y12hd23
        tght(6,9)  =  y12hd13
        tght(6,10) =  cyhd12b

!***************************************************************************c
!***************************************************************************c
!*** insert gt12 (=-gt21), gt22(=-gt11), ht21 (=ht12), ht22 (=ht11).

        do ii = 1, 3
            
            i3 = ii + 3
            i6 = ii + 6
            i9 = ii + 9

            do jj = 1, 5

                j5 = jj + 5

                tght(ii,j5) = -tght(i6,jj)
                tght(i6,j5) = -tght(ii,jj)
                tght(i9,jj) = tght(i3,j5)
                tght(i9,j5) = tght(i3,jj)
 
            end do

        end do
        
!***************************************************************************c
!***************************************************************************c
!**************************** form tzm for rse *****************************c


        do ii = 1, 3

            do jj = ii, 3

                do kk = 1, 3

                    do ll = kk, 3

                m( ii, jj, kk, ll ) = c3d2 * xm * ( ee( ii, jj ) - c1d3 * delta( ii, jj ) ) &
                        * ( ee( kk, ll ) - c1d3 * delta( kk, ll ) ) &
                        + c1d2 * ym * ( ee( ii, kk ) * delta( jj, ll ) &
                        + ee( jj, kk ) * delta( ii, ll ) + ee( ii, ll ) * delta( jj, kk ) &
                        + ee( jj, ll ) * delta( ii, kk ) - 4.d0 * ee( ii, jj ) * ee( kk, ll ) ) &
                        + c1d2 * zm * ( delta( ii, kk ) * delta( jj, ll ) &
                        + delta( jj, kk ) * delta( ii, ll ) - delta( ii, jj ) * delta( kk, ll ) &
                        + ee( ii, jj ) * delta( kk, ll ) + ee( kk, ll ) * delta( ii, jj ) &
                        - ee( ii, kk ) * delta( jj, ll ) - ee( jj, kk ) * delta( ii, ll ) &
                        - ee( ii, ll ) * delta( jj, kk ) - ee( jj, ll ) * delta( ii, kk ) &
                        + ee( ii, jj ) * ee( kk, ll ) )

                    end do

                end do

            end do
 
        end do


        do ii = 1, 5

            if ( ( ii .eq. 1 ) .or. ( ii .eq. 5 ) ) then
               
                tzm( ii, 1 ) = m( mesid( 1, ii ), mesid( 1, ii ), 1, 1 ) &
                               - m(  mesid( 1, ii ), mesid( 1, ii ), 3, 3 ) &
                               - ( m( mesid( 2, ii ), mesid( 2, ii ), 1, 1 ) &
                               - m(  mesid( 2, ii ), mesid( 2, ii ), 3, 3 ) )
                tzm( ii, 2 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 1, 2 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 1, 2 ) ) 
                tzm( ii, 3 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 1, 3 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 1, 3 ) )
                tzm( ii, 4 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 1, ii ), 2, 3 ) &
                               - m( mesid( 2, ii ), mesid( 2, ii ), 2, 3 ) )
                tzm( ii, 5 ) = m( mesid( 1, ii ), mesid( 1, ii ), 2, 2 ) &
                               - m(  mesid( 1, ii ), mesid( 1, ii ), 3, 3 ) &
                               - ( m( mesid( 2, ii ), mesid( 2, ii ), 2, 2 ) &
                               - m(  mesid( 2, ii ), mesid( 2, ii ), 3, 3 ) )

            else

                tzm( ii, 1 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 2, ii ) , 1, 1 ) &
                               - m( mesid( 1, ii ), mesid( 2, ii ), 3, 3 ) )
                tzm( ii, 2 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 1, 2 ) 
                tzm( ii, 3 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 1, 3 ) 
                tzm( ii, 4 ) = 4.d0 * m( mesid( 1, ii ), mesid( 2, ii ), 2, 3 ) 
                tzm( ii, 5 ) = 2.d0 * ( m( mesid( 1, ii ), mesid( 2, ii ) , 2, 2 ) &
                               - m( mesid( 1, ii ), mesid( 2, ii ), 3, 3 ) )

            end if

        end do

!***************************************************************************c
!***************************************************************************c
!*** fill in upper half of m12 (=m11) and m22 (=m11).

        do jj = 1, 5

            j5 = jj + 5

            do ii = 1, jj

                i5 = ii + 5

                tzm( ii, j5 ) = tzm( ii, jj )
                tzm( i5, j5 ) = tzm( ii, jj )

            end do

        end do

!*** fill in the lower half of m12.

        do ii = 1, 5

            i5 = ii + 5

            do jj = ii + 1, 5

                j5 = jj + 5

                tzm( jj, i5 ) = tzm( ii, j5 )

            end do

        end do

    end if

end subroutine calc_lub
