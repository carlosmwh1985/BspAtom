      SUBROUTINE MAKE_F_ANG

!      SUBROUTINE TO CALCULATE THE ANGULAR TERMS OF THE VECTOR POTENTIAL
!        TIMES THE Y_(l+1,0)(th,ph) SPHERICAL HARMONIC, AND THEN INTEGRATED
!        OVER ph (CONSIDERED AS AN ANALYTICAL INTEGRAL)
!      f(r,th) = A(r,th) * Y_{l0+1}^m0(th,ph)

      USE MOD_TYPES
      USE MOD_BSPLINES
      USE MOD_GRID
      USE MOD_PHOTOION
      USE MOD_FAC_VARS

      IMPLICIT NONE

      INTEGER :: ithbet, ithbetmin, ithbetmax, iphbet, iphbetmin, iphbetmax, itd
      INTEGER :: igl, iglmin, iglmax, ithgl, ithglmin, ithglmax, ith, i
      INTEGER :: ibet, ibetmin, ibetmax, lf, mf, l0, m0, mi, ma, m1, m2
      REAL(DP) :: ec0, c0, f1, f2, g1, g2, fth, fkr, w, z0, z1, gz, y2
      REAL(DP) :: r, th, ph, rho, x, y, z, dx, cr, N0lp, N1lp

      COMPLEX(DPC) :: zAref, zdAref
      COMPLEX(DPC), ALLOCATABLE :: zfi(:,:), zfrth(:,:)

      WRITE(6,'(A25)') 'Calculating Angular Terms'

      l0 = lmf(0,1)
      m0 = lmf(0,2)

      mf = m0
      IF( KIND_PI == 4 ) mf = m0 + ma + mph

      lf = lmf(nlm,1)
      m1 = -lf
      m2 = lf
      ma = ABS(moam)
      nm = 1
      IF( KIND_PI >= 5 .AND. bx /= 0.D0 ) THEN
        m1 = m0
        m2 = m1
      END IF

      itd = 0
      IF( KIND_TD == 1 .OR. KIND_PI >= 7 ) THEN
        nm = nlm
        itd = 1
      END IF

!      NORMALIZATION VECTOR POTENTIAL      
      np = 0
      N0lp = 1.D0
      IF( KIND_PI == 4 ) N0lp = SQRT((2.D0*fac(np))/(PI*fac(np+ABS(moam))))

!      NUMBER OF COMPONENTS FOR INTERACTION HAMILTONIAN
      IF( KIND_PI == 3 ) THEN
        ncomp = 1
      ELSE IF( KIND_PI == 4 ) THEN
        ncomp = 3
      ELSE IF( KIND_PI == 8 .OR. KIND_PI == 9 ) THEN
        ncomp = 5                  !A0_z + A_rho + A_z + B_z + B0
      ELSE
        ncomp = 2
      END IF

      ALLOCATE( zAfth(nkp,ka,0:nfib0,nm,ncomp), zfrth(nm,ncomp) )
      zAfth = 0.D0
      zfrth = 0.D0
      
      ibetmin = 1
      ibetmax = nkp - 1

      iglmin = 1
      iglmax = ka

!      RAYLEIGH RANGE
      z0 = kph * (w0**2) / 2.D0

!      CALCULATE INTERACTION HAMILTONIANS

!$OMP PARALLEL DEFAULT(PRIVATE), SHARED(l0,m0,lmax,nlm,lmf,nm,mph,ma,moam,m1,m2, &
!$OMP&    b0,bx,w0,z0,afocus,N0lp,nfib0,ibetmin,ibetmax,iglmin,iglmax,xg,rt,tht, &
!$OMP&    pht,ncomp,zAfth,zdAfth,KIND_PI,itd,kph)

!$OMP DO

 dor:      DO ibet = ibetmin, ibetmax

        f1 = (rt(ibet+1) + rt(ibet) ) / 2.D0
        f2 = (rt(ibet+1) - rt(ibet) ) / 2.D0

 doGLr: DO igl = iglmin, iglmax

          r = f1 + xg(igl) * f2

    doth: DO ith = 0, nfib0

            th = tht(ith)
            ph = pht(ith)
            
            x = COS(th)
            z = r * x

            gz = SQRT(1.D0 + ((z / z0)**2))
            w = w0 * gz

            zfrth = 0.D0
            CALL ARTH(KIND_PI,itd,nm,nlm,lmf,mph,kph,afocus,bx,moam,b0,w0,w,r, &
     &                ith,th,ph,ncomp,zfrth)
            zfrth = N0lp * zfrth

            zAfth(ibet,igl,ith,:,1) = zfrth(:,1)                    !==4:A, >= 5:Er, 8:E_z (Lin.)
!            IF( KIND_PI >= 4 ) zAfth(ibet,igl,ith,:,2) = zfrth(:,2) 
            IF( KIND_PI >= 4 ) THEN             !==4:dA, >=5:Ez, 8:E_r,E_z, B_z, B0 (RVB)
              DO i = 2, ncomp
                zAfth(ibet,igl,ith,:,i) = zfrth(:,i)
              END DO
            END IF

          END DO doth

        END DO doGLr
      END DO dor

!$OMP END DO
!$OMP END PARALLEL

      WRITE(6,*) 'Done'
      
      DEALLOCATE( zfrth )

      END SUBROUTINE MAKE_F_ANG
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ARTH(KIND_A,itd,nm,nlm,lmf,mph,kph,afocus,bx,ma,b,w0,w,r,ith,&
     &                th,ph,ncomp,zf)
      
!      SUBROUTINE TO CALCULATE THE VALUE OF THE INTERACTION HAMILTONIAN
      
      USE MOD_TYPES
      USE MOD_FAC_VARS
      USE MOD_FTW_PH
      USE MOD_BSPLINES, ONLY: ciLP
      USE MOD_PHOTOION, ONLY: A0, A0x, A0y, A0z, B0z
      USE MOD_MATRICES, ONLY: zYlm

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: KIND_A, itd, nm, nlm, lmf(0:nlm,2), mph, ma, ith, ncomp
      REAL(DP), INTENT(IN) :: kph, afocus, r, th, ph, w0, w, b, bx
      
      COMPLEX(DPC), INTENT(OUT) :: zf(nm,ncomp)

      INTEGER :: i, ama, m1, m2, l0, m0, lf, mf, nupts, imin, imax, ish
      REAL(DP) :: x, dx, y, s, s2, z, rho, fth0, fth1, fth2, frho, cw, cr, cp, cm
      REAL(DP) :: ec0p, ec0m, c0, ePl0m0p, ePl0m0m, Pl0m0, Pl1, Pl2, qr, qz, u, du
      REAL(DP) :: T1, T2, T3, Legendre_P, Laguerre_L, THREE_J, BESSJ0, BESSJ1
      REAL(DP) :: xp, yp, zp, rhop, php, cr0, cr1, cz0, cz1, cz, ap, am, muB, fiy
      REAL(DP) :: cth, sth, c2th, cph, sph, cph2, sph2, rpsth, kz, kr, er, ez
      REAL(DP) :: rzero0, rzero1, ex, ey
      COMPLEX(DPC) :: zfkr, zf1, zf2, zTr, zTz, zTf, zTl, zT1, zT2, zTB, zc0
      COMPLEX(DPC) :: zfth0, zfth1, zfthp, zfthm, zfthB
      
      COMPLEX(DPC), ALLOCATABLE :: zfEr(:), zfEz(:), zfBf(:)

      l0 = lmf(0,1)
      m0 = lmf(0,2)

      lf = lmf(nlm,1)
      mf = lmf(nlm,2)
      IF( mf == lf ) THEN
        m1 = -lf
        m2 = lf
      ELSE
        m1 = m0
        m2 = m0
      END IF

      ama = ABS(ma)

      x = COS(th)
      y = SIN(th)
      dx = 1.D0                         !Already included in the Fibonacci Integratal Routine...
      fiy = 1.D0 !/ y
      IF( y == 0.D0 ) fiy = 1.D0

      IF( KIND_A == 3 .OR. KIND_A == 4 ) THEN
        ePl0m0p = Legendre_P(l0+1,m1,x)
        ePl0m0m = Legendre_P(l0+1,m2,x)
      END IF

      z = r * x
      rho = r * y
      kz = kph * z
      kr = kph * rho
      er = rho
      ez = z
      ex = r * SIN(th) * COS(ph)
      ey = r * SIN(th) * SIN(ph)

      IF( bx /= 0.D0 ) THEN
        x = r * SIN(th) * COS(ph)
        y = r * SIN(th) * SIN(ph)
        z = r * COS(th)
        xp = x
        yp = y * COS(bx) + z * SIN(bx)
        zp = z * COS(bx) - y * SIN(bx)
        rhop = SQRT(xp**2 + yp**2)
        kz = kph * zp
        kr = kph * rhop
        er = rhop
        ez = zp
        muB = 0.5D0
        c0 = 1.D0 / (A0 * kph)
      END IF
      
      s = SQRT(2.D0) * rho / w
      s2 = ((rho**2) + (b**2)) / (w**2)

      zf = 0.D0

!      NORMALIZATION FOR SPHERICAL HARMONIC RELATED TO INITIAL STATE
      IF( KIND_A <= 4 ) THEN
        c0 = SQRT(DBLE(2*l0+1) / (4.D0*PI)) * SQRT(fac(l0-m0) / fac(l0+m0))
        ec0p = SQRT(DBLE(2*(l0+1) + 1) / (4.D0*PI)) * SQRT(fac(l0+1-m1) / fac(l0+1+m1))
        ec0m = SQRT(DBLE(2*(l0+1) + 1) / (4.D0*PI)) * SQRT(fac(l0+1-m2) / fac(l0+1+m2))
      END IF

      cw = w0 / w
      zc0 = zi

!      TO BE USED FOR TD CALCULATION

      imin = 0
      imax = 0
      ish = 1
      IF( itd == 1 ) THEN
        imin = 1
        imax = nlm
        ish = 0
      END IF

!      EVALUATE THE VALUE OF THE INTERACTION HAMILTONIAN FOR A GIVEN r AND th

      zf = 0.D0

      IF( KIND_A == 3 ) THEN                  !GAUSSIAN BEAM

        fth0 = 0.D0
        fth1 = ec0p * ePl0m0p
        fth2 = 0.D0
        frho = EXP(-s2)
        zfkr = COS(kph*z)
        cw = 2.D0 * cw
        cr = 0.d0

      ELSE IF( KIND_A == 4 ) THEN            !LAGUERRE-GAUSSIAN BEAM

        fth0 = c0 * Pl0m0
        fth1 = ec0p * ePl0m0p
        fth2 = ec0m * ePl0m0m
        zfkr = EXP(-zi*kph*z)
        cr = cw        
        IF( b == 0.D0 ) THEN
          frho = EXP(-s2) * (s**ama) * Laguerre_L(0,ama,s)
        ELSE
          frho = EXP(-s2)
        END IF
        
      ELSE IF( KIND_A == 5 ) THEN            !RADIAL VECTOR BEAM, WITH BESSEL PROFILE

        nupts = 201
        ALLOCATE( zfEr(nupts), zfEz(nupts), zfBf(nupts) )
        zfEr = 0.D0
        zfEz = 0.D0
        zfBf = 0.D0
        du = 1.D0 / DBLE(nupts-1)

!        cp = DBLE(l0*(l0+1)) / SQRT(DBLE((2*l0+1)*(2*l0+3)))
!        cm = 0.D0
!        IF( l0 /= 0 ) cm = DBLE(l0*(l0-1)) / SQRT(DBLE((2*l0-1)*(2*l0+1)))

        fth0 = 0.D0 !Plm(l0,m0)
        fth1 = 0.D0 !cp * Pl1 - cm * Pl2

        qz = COS(afocus) * kz
        qr = SIN(afocus) * kr
        cr = kph / COS(afocus)            !For Bf

!        DO i = 1, nupts
!          u = DBLE(i-1) * du
!          zfEr(i) = 0.5D0 * zi * BESSJ1(u*qr) * EXP(zi*qz*u)
!          zfEz(i) = -0.5D0 * BESSJ0(u*qr) * EXP(zi*qz*u)
!          zfBf(i) = -0.5D0 * zi * BESSJ1(u*qr) * EXP(zi*qz*u) * u
!        END DO
!        CALL ZSIMPINT(nupts,zfEr,zfEz,zfBf,zTr,zTz,zTf)

          zTr = z !er * zTr                              !Er
        zTz = 0.D0 !TAN(afocus) * ez * zTz            !Ez
        zTB = c0 * muB * B0z * m0                        !B_z

        zc0 = 1.D0

      ELSE IF( KIND_A == 6 ) THEN            !AZIMUTHAL VECTOR BEAM, WITH BESSEL PROFILE

        nupts = 201
        ALLOCATE( zfEr(nupts), zfEz(nupts), zfBf(nupts) )
        zfEr = 0.D0
        zfEz = 0.D0
        zfBf = 0.D0
        du = 1.D0 / DBLE(nupts-1)

        qz = COS(afocus) * kz
        qr = SIN(afocus) * kr
        cr = SIN(th)
        IF( cr == 0.D0 ) cr = eps
        cr = COS(th) / cr

        IF( bx /= 0.D0 ) THEN
          cth = COS(th)
          sth = SIN(th)
          c2th = COS(2.D0*th)
          cph = COS(ph)
          sph = SIN(ph)
          cph2 = COS(ph)**2
          sph2 = SIN(ph)**2
          cr0 = -rho*(cth*(cph2 + sph2*(COS(bx)**2)) - 0.5*sth*sph*SIN(2.D0*bx)) &
     &              - z*SIN(bx)*(cth*sph*COS(bx) - sth*SIN(bx))
          cr1 = SIN(bx) * (rho*sph*cph*SIN(bx) + z*cph*COS(bx))
          rpsth = rhop * sth
          IF( rpsth == 0.D0 ) rpsth = eps
          cr0 = cr0 / rpsth
          cr1 = cr1 / rpsth
          cz0 = cth*sph*SIN(bx) + sth*COS(bx)
          cz1 = -cph * SIN(bx)
          rpsth = sth
          IF( rpsth == 0.D0 ) rpsth = eps
          cz0 = cz0 / rpsth
          cz1 = cz1 / rpsth
          ap = DBLE(l0*(l0+1)) / SQRT(DBLE((2*l0+1)*(2*l0+3)))
          am = 0.D0
          zfthp = zYlm(l0+1,m0,ith)
          IF( l0 /= 0 ) THEN
            am = -DBLE(l0*(l0-1)) / SQRT(DBLE((2*l0-1)*(2*l0+1)))
            zfthm = zYlm(l0-1,m0,ith)
          END IF
          zfth0 = zi * m0 * zYlm(l0,m0,ith)            ! = d Ylm / d ph
          zfth1 = ap * zfthp + am * zfthm            ! = SIN(th) d Ylm / d th
        END IF

!        DO i = 1, nupts
!          u = DBLE(i-1) * du
!          zfEr(i) = 0.5D0 * BESSJ1(u*qr) * EXP(zi*qz*u) * u
!          zfEz(i) = -0.5D0 * zi * BESSJ0(u*qr) * EXP(zi*qz*u) * u
!          zfBf(i) = -0.5D0 * BESSJ1(u*qr) * EXP(zi*qz*u)
!        END DO
!          CALL ZSIMPINT(nupts,zfEr,zfEz,zfBf,zT1,zT2,zT3)

!        zT1 = COS(afocus) * cr * zT1            !Br
!        zT2 = -SIN(afocus) * zT2                  !Bz
!        zT3 = 0.D0
!        TO VERIFY INDIVIDUAL TRANSITIONS
        zTr = 0.D0
        zTz = 0.D0
        zTB = c0 * muB * B0z * m0 * zYlm(l0,m0,ith)

        zc0 = m0

        IF( bx /= 0.D0 ) THEN
          zTr = -zi * COS(afocus) * zTr * (cr0 * zfth0 + cr1 * zfth1)      !Br
          zTz = -zi * SIN(afocus) * zTz * (cz0 * zfth0 + cz1 * zfth1)      !Bz
          zTB = c0 * muB * B0z * m0 * zYlm(l0,m0,ith)
          zc0 = 1.D0
        END IF

      ELSE IF( KIND_A == 7 ) THEN             !VECTOR POTENTIAL - AHARONOV--BOHM EFFECT

        rpsth = r * SIN(th)
        IF( rpsth == 0.D0 ) rpsth = eps
        
        zTr = m0 * B0z / (PI * (rpsth**2))
        zTz = 0.D0
        zTB = 0.D0

        zc0 = 1.D0
      
      ELSE IF( KIND_A == 8 ) THEN            !TOROIDAL MOMENT: LIN: POL + RVB

        rzero0 = 2.40482556D0             !1st ZERO J_0
        rzero1 = 3.83170597D0             !1st ZERO J_1

!        LINEAR POLARIZED FIELD
        
        zTl = A0z * z + A0y * ey + A0x * ex

!        RADIAL POLARIZED FIELD 

        nupts = 201
        ALLOCATE( zfEr(nupts), zfEz(nupts), zfBf(nupts) )
        zfEr = 0.D0
        zfEz = 0.D0
        zfBf = 0.D0
        du = 1.D0 / DBLE(nupts-1)

        qz = COS(afocus) * kz
        qr = SIN(afocus) * kr
        cr = kph / COS(afocus)            !For Bf

        DO i = 1, nupts
          u = DBLE(i-1) * du
          zfEr(i) = 0.5D0 * zi * BESSJ1(u*qr) * EXP(zi*qz*u)
          zfEz(i) = -0.5D0 * BESSJ0(u*qr) * EXP(zi*qz*u)
          zfBf(i) = -0.5D0 * zi * BESSJ1(u*qr) * EXP(zi*qz*u) * u
        END DO
        CALL ZSIMPINT(nupts,zfEr,zfEz,zfBf,zTr,zTz,zTf)

        IF( qr >= rzero0 ) zTz = 0.D0           !AFTER THE 1st ZERO...
        IF( qr >= rzero1 ) THEN
          zTr = 0.D0
          zTf = 0.D0
        END IF

        zTr = er * zTr                          !Er
        zTz = TAN(afocus) * ez * zTz            !Ez
        zTB = -zi * cr * zTf * fiy              !B_z

        zc0 = 1.D0

      ELSE IF( KIND_A == 9 ) THEN            !TOROIDAL MOMENT: Z + RHO

        zTl = z
        zTr = er
        zTz = 0.D0
        zTB = 0.D0

      END IF

      IF( KIND_A == 3 .OR. KIND_A == 4 ) THEN

!        FIRST COMPONENT: A.p (4), ELECTRIC FIELD (5), MAGNETIC FIELD (6)

        zf1 = dx * cw * frho * zfkr * zT1 * fth1
        zf2 = dx * cw * frho * CONJG(zfkr) * zT2 * fth2
        zf(1,1) = zc0 * (zf1 + zf2) !* zIph(:,i2))      ! q A . p

!        SECOND COMPONENT: p.A (4), MAGNETIC FIELD (5), ELECTRIC FIELD (6)

        zf1 = dx * cr * frho * zfkr * zTf * fth0
        zf2 = dx * cr * frho * CONJG(zfkr) * fth0
        zf(1,2) = zc0 * (zf1 + zf2) !* zIph(:,i4))       ! q p . A

      ELSE IF( KIND_A >= 5 ) THEN      !BESSEL VECTOR BEAMS

!       LENGTH FORM OR POINCARÉ GAUGE
        DO i = imin, imax
          lf = lmf(i,1)
          mf = lmf(i,2)
          zfth0 = zYlm(lf,mf,ith)
          IF( KIND_A == 6 .AND. bx /= 0.D0 ) zfth0 = 1.D0
          IF( KIND_A < 8 ) THEN
            zf(i+ish,1) = zc0 * (zTr + zTB) * zfth0            !Er
            zf(i+ish,2) = zc0 * (zTz + zTB) * zfth0            !Ez
          ELSE
!          IF( KIND_A >= 8 ) THEN
            cp = DBLE(lf*(lf+1)) / SQRT(DBLE((2*lf+1)*(2*lf+3)))
            zfthB = cp * zYlm(lf+1,mf,ith)
            IF( lf >= 1 ) THEN
              cm = DBLE(lf*(lf-1)) / SQRT(DBLE((2*lf-1)*(2*lf+1)))
              zfthB = zfthB - cm * zYlm(lf-1,mf,ith)
            END IF
            zf(i+ish,1) = zTl * zfth0           !Linear Polarized
            zf(i+ish,2) = zTr * zfth0           !RVB: Er
            zf(i+ish,3) = zTz * zfth0           !RVB: Ez
            zf(i+ish,4) = zTB * zfthB           !RVB: Bph
            zf(i+ish,5) = zfth0                 !B0
          END IF
        END DO

      END IF

      IF( ALLOCATED(zfEr) ) DEALLOCATE( zfEr )
      IF( ALLOCATED(zfEz) ) DEALLOCATE( zfEz )
      IF( ALLOCATED(zfBf) ) DEALLOCATE( zfBf )

      END SUBROUTINE ARTH
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ZINT_TH

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_PHOTOION
      USE MOD_MATRICES
      USE MOD_FTW_PH

      IMPLICIT NONE

      INTEGER :: ibet, ibetmin, ibetmax, igl, iglmin, iglmax, lf, mf, il, jl
      INTEGER :: mi, m1, m2, i, ithph, ni, ni_max
      REAL(DP) :: f1, f2, r
      COMPLEX(DPC) :: zIthref
      COMPLEX(DPC), ALLOCATABLE :: zfth(:), zdfth(:)
      CHARACTER(LEN=18) :: cfile

!      LAPACK VARIABLES
      CHARACTER(LEN=1) :: TRANS
      INTEGER :: M, N, LDA, INCX, INCY
      COMPLEX(DPC) :: ZALPHA, ZBETA, ZDOTU
      COMPLEX(DPC), ALLOCATABLE :: A(:,:), X(:), Y(:)

      WRITE(6,'(A29)') 'Calculating Integrals Over th'

      ibetmin = 1 !nbc1
      ibetmax = nkp - 1 !- nbc2
      iglmin = 1
      iglmax = ka

      lf = lmf(nlm,1)
      mf = lmf(nlm,2)
      m1 = -lf
      m2 = lf

      IF( KIND_PI >= 5 .AND. bx /= 0.D0 ) THEN
        m1 = mf
        m2 = m1
      END IF

      ni_max = ncomp
      IF( KIND_PI == 8 ) THEN
        ni_max = 4
      ELSE IF( KIND_PI == 9 ) THEN
        ni_max = 2
      END IF

      ALLOCATE( zfth(0:nfib0), zdfth(0:nfib0) )
      zfth = 0.D0
      zdfth = 0.D0

!      PERFORMS THE INTEGRAL ON th FOR EACH r AND FOR EACH FINAL l

      IF( ALLOCATED(zIth) ) DEALLOCATE( zIth )
      ALLOCATE( zIth(nkp,ka,nlm,nm,ncomp) )
      zIth = 0.D0

!      DO i = 0, nfib0
!        WRITE(75,'(20g20.10)') tht(i), (ciLP(i,il,m1), il=0,lf)
!      END DO

100      FORMAT('CSs/Field_l_',I2.2,'.dat')

!$OMP PARALLEL DEFAULT(PRIVATE), SHARED(nlm,nm,lmf,KIND_PI,KIND_TD,ibetmin, &
!$OMP&  ibetmax,iglmin,iglmax,rt,xg,ncomp,ni_max,zYlm,zAfth,nfib,nfib0,zIth)

!$OMP DO

dolmf:DO il = 1, nlm

        lf = lmf(il,1)
        mf = lmf(il,2)

!        IF( KIND_PI >= 5 ) THEN
!          WRITE(cfile,100) lf
!          OPEN(UNIT=50, FILE=cfile, ACTION='WRITE')
!        END IF        

  dom0: DO jl = 1, nm

          mi = mf
          IF( KIND_TD == 1 ) mi = lmf(jl,2)

    dort: DO ibet = ibetmin, ibetmax

!            f1 = ( rt(ibet+1) + rt(ibet) ) / 2.D0
!            f2 = ( rt(ibet+1) - rt(ibet) ) / 2.D0

!            iglmax = ka
!            IF( ibet == ibetmin .OR. ibet == ibetmax ) iglmax = 1

      doGL: DO igl = iglmin, iglmax

!              r = f1 + xg(igl) * f2

              zfth = 0.D0
              zdfth = 0.D0

              DO ni = 1, ni_max
              
                DO ithph = 0, nfib0
                  zfth(ithph) = CONJG(zYlm(lf,mf,ithph)) * zAfth(ibet,igl,ithph,jl,ni)
                 END DO

                CALL FIBINT(nfib,nfib0,zfth,zIthref)
              
!                IF( KIND_PI >= 5 .AND. ni == 1 .AND. jl == 1) WRITE(50,'(3G20.10)') r, zIthref
!                IF( KIND_PI == 7 .AND. ni == 1 ) WRITE(50,'(3G20.10)') r, zIthref
!                IF( KIND_PI == 8 .AND. ni == 4 .AND. lf==0 ) WRITE(50,'(3G20.10)') r, zIthref
        
                zIth(ibet,igl,il,jl,ni) = zIthref

              END DO

            END DO doGL

          END DO dort

        END DO dom0

!        CLOSE(50)

      END DO dolmf

!$OMP END DO
!$OMP END PARALLEL

      DEALLOCATE( zfth, zdfth )

      END SUBROUTINE ZINT_TH
