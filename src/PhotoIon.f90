      SUBROUTINE TRANS_AMP

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_PHOTOION
      USE MOD_MATRICES
      USE MOD_FAC_VARS
!      USE MOD_VMV_VARS

      IMPLICIT NONE

      INTEGER :: ibra, jket, n0, l0, m0, nf, lf, mf, mp0, il, jl, nmin, nmax
      INTEGER :: npts, i, iref, left, mflag, j1, j, ni, li, mi
      REAL(DP) :: An, Am, a1, a2, b1, b2, c0, c1, c2, c3, c4, c5, T3ja, T3jb
      REAL(DP) :: Gfi, Rfi, Normi, Normf, Nlp, r, bsp(k), fr, sumwf, THREE_J
      COMPLEX(DPC) :: zfvmv
      CHARACTER(LEN=26) :: cfile

      REAL(DP), ALLOCATABLE :: ciall(:)
      COMPLEX(DPC), ALLOCATABLE :: zfijall(:)

!      BLAS VARIABLES
      CHARACTER(LEN=1) :: TRANS, UPLO
      INTEGER :: M, N, LDA, INCX, INCY
      REAL(DP) :: ALPHA, BETA, DDOT
      COMPLEX(DPC) :: zALPHA, zBETA, zDOTU
      REAL(DP), ALLOCATABLE :: A(:,:), u(:), v(:), x(:)
!      COMPLEX(DPC), ALLOCATABLE :: zx(:), zy(:) !zu(:), zv(:), zw(:), zx(:)
!      COMPLEX(DPC), ALLOCATABLE :: zA(:,:)

!$    DOUBLE PRECISION :: OMP_GET_WTIME, WT_OLD, WT_NEW

!$    WT_OLD = OMP_GET_WTIME()

      WRITE(6,'(/,A33)') 'Calculating Transition Amplitudes'

      n0 = n0_ini
      l0 = lmf(0,1)
      m0 = lmf(0,2)

      WRITE(6,'(A14,3I3)') 'Initial State:', n0+l0, l0, m0

      lf = lmf(nlm,1)
      mf = lmf(nlm,2)

      WRITE(6,'(A26,2G20.10)') 'Energy limits Final State:', Enl(1,lf), Enl(n1_max,lf)

!.....LENGTH OR VELOCITY GAUGE. PLANE WAVES
 ifg:      IF( KIND_PI == 1 .OR. KIND_PI == 2 ) THEN
 
        ALLOCATE( T_fi(n0_fin:n1_fin,lf:lf) )
        T_fi = 0.D0

        ALLOCATE( A(nfun,nfun), u(nfun), v(nfun), x(nfun) )

!        VARIABLES FOR BLAS SUBROUTINE

        TRANS = 'N'
        M = nfun
        N = nfun
        ALPHA = 1.D0
        LDA = nfun
        INCX = 1
        BETA = 0.D0
        INCY = 1

        T3ja = THREE_J(lf,1,l0,-mf,mph,m0)
        IF( KIND_PI == 1 ) THEN
          T3jb = THREE_J(lf,1,l0,0,0,0)
          c1 = (-1.D0)**(lf+l0+mf) * SQRT( DBLE((2*lf+1)*(2*l0+1)) ) * T3ja * T3jb
          c0 = 1.0D0
          A(:,:) = c1 * rij(:,:,1)
        ELSE
          c0 = SQRT(DBLE(l0 + 1)) * T3ja
          c1 = 0.D0
          c2 = 0.D0
          IF( lf == l0 + 1 ) THEN
            c1 = DBLE(l0 + 1)                  !Coeff. for u(r) * (1/r) * u(r) term
            c2 = -1.D0                        !Coeff. for u(r) * d u(r)
          ELSE IF( lf == l0 - 1 ) THEN
            c1 = DBLE(l0)
            c2 = 1.D0
          END IF
          A(:,:) = c1 * rij(:,:,1) + c2 * rij(:,:,2)
        END IF
        
        T_fi = 0.D0

        x(:) = ci_ini(:) 

!        BLAS SUBROUTINE FOR MATRIX-VECTOR MULTIPLICATION
!        A IS AN m BY n MATRIX
!        y := alpha*A*x + beta*y
        CALL DGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,v,INCY)

donf12: DO ni = n0_fin, n1_fin

          An = SQRT( 2.D0 / (E_fin(ni+1) - E_fin(ni-1)))            !NORMALIZATION BY DENSITY OF STATES
          u(:) = ci_fin(:,ni)

!          BLAS FUNCTION TO PERFORM A INNER PRODUCT
          T_fi(ni,lf) = An * c0 * DDOT(N,u,INCX,v,INCY)

        END DO donf12

        DEALLOCATE( A, u, v, x )


!.....VELOCITY GAUGE, FOR GAUSSIAN OR LAGUERRE-GAUSSIAN BEAMS
      ELSE IF( KIND_PI >= 3 ) THEN ifg

        mp0 = m0 - mph
        a1 = DBLE(l0) * SQRT(DBLE(l0+1))
        a2 = (-1.D0)**(l0+m0) !SQRT(3.D0) * THREE_J(l0+1,1,l0,mp0,mph,-m0)
        c1 = a1 * a2

        b1 = -SQRT(DBLE(l0+1))
        b2 = a2
        c2 = b1 * b2

        c0 = A0                        !GLOBAL FACTOR FROM THE KINECTIC En. TERM
                                          ! x NORM. OV x Intensity Vec. Pot.

        c3 = 0.5D0 * A0
        c4 = -(c_au**2) * c3                  ! q Phi_sc, from Lorenz gauge
        c5 = 0.D0

!!!        nm = 2 * lf + 1

        IF( KIND_PI >= 5 ) THEN
          c0 = c0 !* Eph
          c1 = 1.D0                  !TERM WITH WAVE FUNCTION (Er field)
          c2 = 0.D0                  !TERM WITH WAVE FUNCTION (B field)
          c3 = 0.D0
          c4 = 1.D0                  !TERM WITH WAVE FUNCTION (Ez field)
          IF( KIND_PI == 6 ) THEN
            c0 = c0 / c_au
          ELSE IF( KIND_PI >= 8 ) THEN
            c0 = 1.D0 !SQRT(I0 / I0_au)      !AMPLITUDE LIN. POL. FIELD
            c1 = 1.D0
            c2 = 1.D0
            c3 = 1.D0 ! A01 * Eph            !PREFACTOR RVB
            c4 = 1.D0
            c5 = 0.5D0 !* B0z                !Bohr Magenton; CONSTANT MAGNETIC FIELD
!            ncomp = ncomp + 1
          END IF          
        END IF

        nmin = n0
        nmax = n0

        IF( KIND_TD == 1 .OR. bx /= 0.D0 .OR. KIND_PI == 7 ) THEN
          nmin = 1
          nmax = n1_max
        END IF

        CALL SEL_STATES(nmin,nmax)

        WRITE(6,'(A16,I5,A18,I5)') 'Num. States Bra:', nbra, ', Num. States Ket:', nket
        WRITE(6,'(A19,I2)') 'Num. of Components:', ncomp

        ALLOCATE( zT_fi(nbra,nket,ncomp), ciall(ncomp), zfijall(ncomp) )
        zT_fi = 0.D0

        ciall = (/c1, c2, c3, c4, c5/)

        iref = 1
        IF( KIND_SCP == 1 ) iref = 3
        IF( KIND_PI >= 5 ) iref = 2

!        VARIABLES FOR BLAS SUBROUTINES

        TRANS = 'N'
        UPLO = 'U'
        M = nfun
        N = nfun
        ZALPHA = 1.D0
        LDA = nfun
        INCX = 1
        ZBETA = 0.D0
        INCY = 1

        WRITE(6,'(A32)') 'Calculating Interaction Matrices'

!$      WT_OLD = OMP_GET_WTIME()

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nlm,lmf,rEki,zAij,c0,c1,c2,c3,c4,ciall,cinl, &
!$OMP&  zT_fi,iref,KIND_PI,KIND_TD,nbra,nket,nl_bra,nl_ket,nfun,ncomp,Sij,UPLO, &
!$OMP&  TRANS,M,N,zALPHA,LDA,INCX,zBETA,INCY,B0z)

      ALLOCATE( zx(N), zy(N), zv_temp(N), zA(N,N) )

!$OMP DO

 dobra: DO ibra = 1, nbra

          nf = nl_bra(ibra,1)
          lf = nl_bra(ibra,2)
          mf = nl_bra(ibra,3)
          il = nl_bra(ibra,4)

!         An = rEki(nf,lf)                        !NORMALIZATION BY DENSITY OF STATES
!         IF( KIND_TD == 1 ) An = 1.D0
          An = 1.D0

   doket: DO jket = ibra, nket

            ni = nl_bra(jket,1)
            li = nl_bra(jket,2)
            mi = nl_bra(jket,3)
            jl = nl_bra(jket,4)

!           Am = rEki(ni,li)                        !NORMALIZATION BY DENSITY OF STATES
!           IF( KIND_TD == 1 ) Am = 1.D0
            Am = 1.D0

            zA = 0.D0

!           SUBROUTINE TO PERFORM VECTOR - MATRIX -VECTOR MULTIPLICATION
!           BASED ON ZHEMV AND ZDOTU BLAS ROUTINES
            zx(1:N) = cinl(1:N,ni,li)
            zy(1:N) = cinl(1:N,nf,lf)
            
            DO i = 1, MIN(4,ncomp)
              zT_fi(ibra,jket,i) = 0.D0
              IF( ciall(i) /= 0 ) THEN
                zA(:,:) = ciall(i) * zAij(:,:,il,jl,i)
                CALL ZHVMV(UPLO,N,zx,zy,zv_temp,zA,zfvmv)
                zT_fi(ibra,jket,i) = An * zfvmv * Am
              END IF
            END DO
            
!           B_0:
            zT_fi(ibra,jket,5) = 0.D0
            IF( li == lf .AND. mi == mf .AND. mi /= 0 .AND. B0z /= 0.D0 ) THEN
              zA(:,:) = ciall(5) * mi * Sij(:,:)
              CALL ZHVMV(UPLO,N,zx,zy,zv_temp,zA,zfvmv)
              zT_fi(ibra,jket,5) = An * zfvmv * Am
            END IF

          END DO doket

        END DO dobra

!$OMP END DO

      DEALLOCATE( ZA, zx, zy, zv_temp )

!$OMP END PARALLEL

!$    WT_NEW = OMP_GET_WTIME()
!$    WRITE(6,'(/,A,1PG12.4)') 'OMP_GET_WTIME:', WT_NEW - WT_OLD

      OPEN( UNIT=60, FILE='CSs/MatElem_All.dat', ACTION='WRITE' )
      WRITE(60,*) n1_max, nbra, nket

      DO ibra = 1, nbra
        DO jket = ibra, nket
          WRITE(60,500) ibra, jket, (zT_fi(ibra,jket,i), i=1,ncomp)
        END DO
      END DO

      CLOSE(60)

500   FORMAT(2I8,X,20G20.10)

      END IF ifg

      END SUBROUTINE TRANS_AMP
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE CROSS_SECTIONS

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_PHOTOION
      USE MOD_FTW_PH, ONLY : FTFtw

      IMPLICIT NONE

      INTEGER :: n0, l0, m0, nf, lf, mf, li, mi, lj, mj, l1, l2, il, jl
      INTEGER :: ibra, jket, iref, nbi, nbj, ni, nmin, nmax, nj
      CHARACTER(LEN=27) :: csfile
      CHARACTER(LEN=22) :: crfile, czfile
      CHARACTER(LEN=23) :: cmrfile, cmzfile
      REAL(DP) :: c0, c1, d1, d2, csl, M_au, E0, Ef, FTw, wfi
      LOGICAL :: fsind

      WRITE(6,'(/,A26)') 'Calculating Cross Sections'

      M_au = (a_au**2) * 1.0D18            !a.u. TO Mb, for Cross Sections

      n0 = n0_ini
      l0 = lmf(0,1)
      m0 = lmf(0,2)

      iref = 1
!!!      IF( KIND_TD == 1 ) iref = l0*nlm + n0

      E0 = Enl(n0,l0)
      WRITE(6,*) 'E0=', E0

      lf = lmf(nlm,1)
      mf = lmf(nlm,2)

      ALLOCATE( T2fi(n1_max,nlm) )
      IF( KIND_PI >= 5 ) ALLOCATE( S2fi(n1_max,nlm) )

      T2fi = 0.D0

      IF( KIND_PI == 1 .OR. KIND_PI == 2 ) THEN            !PLANE WAVE
        c0 = 4.D0 * (PI**2) / (c_au)
      ELSE IF( KIND_PI >= 3 ) THEN      !GAUSSIAN OR LG BEAM
!        c0 = 8.D0 * PI * (PI * (rb**2)) / (c_au * (w0**2))
        c0 = 4.D0 * (PI**2) / (c_au)
      END IF
      c1 = 1.D0 / DBLE(2*l0 + 1)

      ibra = 0
dolf: DO il = 1, nlm

        lf = lmf(il,1)
        mf = lmf(il,2)

        IF( KIND_PI == 1 ) THEN

          OPEN( UNIT=30, FILE='CSs/CrossSection_Len.dat', ACTION='WRITE' )

        ELSE IF( KIND_PI == 2 ) THEN

          OPEN( UNIT=30, FILE='CSs/CrossSection_Vel.dat', ACTION='WRITE' )

        ELSE IF( KIND_PI == 3 .OR. KIND_PI == 4 ) THEN

          IF( b0 == 0.D0 .AND. mb == 0 ) THEN
            WRITE(csfile,100) lf
          ELSE
            WRITE(csfile,110) lf, ABS(mf)
            IF( mf < 0 ) WRITE(csfile,115) lf, ABS(mf)
          END IF
          OPEN( UNIT=30, FILE=csfile, ACTION='WRITE' ) 

          n0_fin = n01(lf,1)
          n1_fin = n01(lf,2)
          
        ELSE IF( KIND_PI == 5 .OR. KIND_PI == 6 ) THEN

          WRITE(csfile,110) lf, ABS(mf)        
          WRITE(crfile,200) lf, ABS(mf)
          WRITE(czfile,250) lf, ABS(mf)
          WRITE(cmrfile,300) lf, ABS(mf)
          WRITE(cmzfile,350) lf, ABS(mf)
          IF( mf < 0 ) THEN
            WRITE(csfile,115) lf, ABS(mf)
            WRITE(crfile,210) lf, ABS(mf)
            WRITE(czfile,260) lf, ABS(mf)
            WRITE(cmrfile,310) lf, ABS(mf)
            WRITE(cmzfile,360) lf, ABS(mf)
          END IF
          OPEN( UNIT=30, FILE=crfile, ACTION='WRITE' )
          OPEN( UNIT=35, FILE=czfile, ACTION='WRITE' )
          OPEN( UNIT=40, FILE=cmrfile, ACTION='WRITE' )
          OPEN( UNIT=45, FILE=cmzfile, ACTION='WRITE' )
          OPEN( UNIT=50, FILE=csfile, ACTION='WRITE' )
          IF( il == 1 ) OPEN( UNIT=37, FILE='CSs/FourierTG.dat', ACTION='WRITE' )

          n0_fin = n01(lf,1)
          n1_fin = n01(lf,2)
        
        ELSE IF( KIND_PI == 7 ) THEN
          
          WRITE(cmrfile,300) lf, ABS(mf)
          IF( mf < 0 ) WRITE(cmrfile,310) lf, ABS(mf)
          OPEN( UNIT=40, FILE=cmrfile, ACTION='WRITE' )

          nbi = n01(lf,3)
          
          n0_fin = 1
          n1_fin = nbi
          
        END IF

  doEf: DO nf = n0_fin, n1_max

          ibra = ibra + 1

          IF( KIND_PI == 1 ) THEN                  !LENGTH GAUGE
            d1 = E_fin(nf) - E_ini(n0)
            Ef = E_fin(nf)
            d2 = T_fi(nf,il)**2
          ELSE IF( KIND_PI == 2 ) THEN            !VELOCITY GAUGE
            d1 = 1.D0 / (E_fin(nf) - E_ini(n0))
            Ef = E_fin(nf)
            d2 = T_fi(nf,il)**2
          ELSE
            Ef = Enl(nf,lf)
            d1 = 1.D0 / (Ef - E0)
            d2 = ABS(zT_fi(ibra,iref,1))**2
            T2fi(nf,il) = d2
            IF( KIND_PI >= 5 ) S2fi(nf,il) = ABS(zT_fi(ibra,iref,2))**2
          END IF

          csl = M_au * c0 * c1 * d1 * d2

          wfi = Ef - E0
          FTw = 1.D0 !ABS(FTFtw(Eph,Ef-E0,ncyc))**2

          IF( ABS(wfi) < 1.0D-15 ) wfi = 1.0D-15
          IF( nf <= n1_fin .AND. KIND_PI < 5 ) WRITE(30,400) Ef, csl
          IF( KIND_PI == 5 .OR. KIND_PI == 6 ) THEN
            WRITE(30,400) Ef, wfi * FTw * T2fi(nf,il)
            WRITE(35,400) Ef, wfi * FTw * S2fi(nf,il)
            WRITE(40,410) Ef, zT_fi(ibra,iref,1)
            WRITE(45,410) Ef, zT_fi(ibra,iref,2)
            IF( il == 1 ) WRITE(37,400) Ef, FTw
            WRITE(50,410) Ef, c0*wfi*T2fi(nf,il)*M_au, c0*wfi*S2fi(nf,il)*M_au
          ELSE IF( KIND_PI == 7 ) THEN
            li = lf
            mi = mf
            DO jl = 1, nlm
              lj = lmf(jl,1)
              mj = lmf(jl,2)
              nbj = n01(lj,3)
              jket = (jl-1)*n1_max
              DO nj = 1, nbj
                jket = jket + 1
                IF( nf <= nbi ) WRITE(40,420) nf+li, li, mi, nj+lj, lj, mj, zT_fi(ibra,jket,1)
              END DO
            END DO

          END IF

!          IF( KIND_PI >= 5 ) WRITE(30,400) Ef, d2*FTw/wfi

        END DO doEf

        CLOSE(30)
        IF( KIND_PI >= 5 ) THEN
          CLOSE(35)
          CLOSE(40)
          CLOSE(45)
          CLOSE(50)
        END IF

      END DO dolf

      IF( KIND_PI >= 3 .AND. KIND_PI <= 6 ) CALL SUM_LF

100      FORMAT('CSs/CrossSection_l',I2.2,'.dat')
110      FORMAT('CSs/CrossSection_l',I2.2,'+',I2.2,'.dat')
115      FORMAT('CSs/CrossSection_l',I2.2,'-',I2.2,'.dat')
200      FORMAT('CSs/OscStr_r_',I2.2,'+',I2.2,'.dat')
210      FORMAT('CSs/OscStr_r_',I2.2,'-',I2.2,'.dat')
250      FORMAT('CSs/OscStr_z_',I2.2,'+',I2.2,'.dat')
260      FORMAT('CSs/OscStr_z_',I2.2,'-',I2.2,'.dat')
300      FORMAT('CSs/MatElem_r_',I2.2,'+',I2.2,'.dat')
310      FORMAT('CSs/MatElem_r_',I2.2,'-',I2.2,'.dat')
350      FORMAT('CSs/MatElem_z_',I2.2,'+',I2.2,'.dat')
360      FORMAT('CSs/MatElem_z_',I2.2,'-',I2.2,'.dat')

400      FORMAT(2G20.10E3)
410      FORMAT(3G20.10E3)
420      FORMAT(2(3I3,X),2G20.10)
500      FORMAT(2I8,X,10G20.10)

      END SUBROUTINE CROSS_SECTIONS
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE SUM_LF

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_PHOTOION
      USE MOD_FTW_PH

      IMPLICIT NONE

      INTEGER :: n0, l0, ni, li, nf, lf, mf, lmin, nn, il, i, j, ij, n, ntimes, KIND_ITP
      REAL(DP) :: Emin, Emax, dE, th, thi, thf, dth, v1, v2, vref, Ei, Ej, kj, hj
      REAL(DP) :: E0, Ef, dw, M_au, c0, c1, coeff_cs, coeff_dp
      COMPLEX(DPC) :: zsum, zT, zS, zfw, zgw, zsw, zph

      REAL(DP), ALLOCATABLE :: Es(:), t2i(:), t2f(:), t2_pts(:,:), sigmal(:)
      REAL(DP), ALLOCATABLE :: TA1(:), TAr2(:), TAi2(:)
      COMPLEX(DPC), ALLOCATABLE :: zYlm(:,:), zTS(:)

      vref = 1.D-10

      n0 = n0_ini
      l0 = lmf(0,1)

      li = lmf(1,1)
      lf = lmf(nlm,1)

      lmin = MIN(l0,li)

      E0 = Enl(n0,l0)

      ALLOCATE( Es(li:lf), sigmal(lmin:lmax) )

!      DETERMINE THE UPPER LIMIT FOR THE ENERGIES CLOSE TO THE THRESHOLD

      Es = 0.D0
      DO il = li, lf
        ni = n01(il,1)
        Es(il) = Enl(ni,il)
      END DO
      Emin = MAXVAL(Es)
      
!      DETERMINE THE LOWER LIMIT FOR THE HIGH-ENERGY LIMIT

      Es = 0.D0
      DO il = lmin, lmax
        nf = n01(il,2)
        Es(il) = Enl(nf,il)
      END DO
      Emax = MIN(MINVAL(Es),Emax_fin)

      WRITE(6,'(/,A20,2G15.5)') 'Global Energy Limits:', Emin, Emax

      n = ABS(nEpts)
      IF( nEpts >= 0 ) THEN
        dE = (Emax - Emin) / DBLE(n)
        ALLOCATE( Enf(0:n) )
        Enf = 0.D0
        Enf(0) = Emin
        DO i = 1, n
          Enf(i) = Emin + DBLE(i) * dE
        END DO
      ELSE
        ntimes = NINT(Emax / Eref)
        dE = (Emax - Emin) / DBLE(n)
        n = n + ntimes
        ALLOCATE( Enf(0:n) )
        Enf = 0.D0
        i = 0
        j = 0
        Ei = Eref
  doEi: DO
          Ej = Emin + DBLE(j) * dE
          IF( Ej < Ei ) THEN
            Enf(i) = Ej
          ELSE IF( Ej > Ei ) THEN
            Enf(i) = Ei
            Enf(MIN(i+1,n)) = Ej
            Ei = Ei + Eref
            i = i + 1
          ELSE
            Enf(i) = Ei
            Enf(MIN(i+1,n)) = Ej + dE
            Ei = Ei + Eref
            i = i + 1
          END IF
!           WRITE(6,'(I4,3G16.7)') i, Enf(i), Ej, Ei
          j = j + 1
           i = i + 1
          IF( i > n ) EXIT doEi
        END DO doEi
      END IF

      DEALLOCATE( Es )
      ALLOCATE( t2_pts(0:n,nlm), t2f(0:n) )
      ALLOCATE( TAr2(0:n), TAi2(0:n), zT_Ef(0:n,nlm) )
      IF( KIND_SCP == 1 ) ALLOCATE( zS_Ef(0:n,nlm), zTS(nlm) )
      t2_pts = 0.D0

      M_au = (a_au**2) * 1.0D18
      c0 = 8.D0 * (PI**2) / c_au
      c1 = 1.D0 / DBLE(2*l0 + 1)

!      OPEN( UNIT=30, FILE='Test_rTA.dat', ACTION='WRITE' )
!      OPEN( UNIT=31, FILE='Test_iTA.dat', ACTION='WRITE' )
      WRITE(6,*) 'Interpolating'
      DO il = 1, nlm
        lf = lmf(il,1)
        ni = n01(lf,1)
        nf = n1_max !n01(lf,2)
        nn = nf - ni
        ALLOCATE( t2i(0:nn), Es(0:nn), TA1(0:nn) )
        Es(0:nn) = Enl(ni:nf,lf)
!        INTERPOLATING |T|^2 TERMS
        t2i(0:nn) = T2fi(ni:nf,il)
        t2f = 0.D0
        CALL CUBSPL(nn,Es,t2i,n,Enf,t2f)
        t2_pts(:,il) = t2f(:)
!        INTERPOLATING REAL PART OF T
        TA1(0:nn) = REAL(zT_fi(ni+lf*nlm:nf+lf*nlm,n0+l0*nlm,1))
        v1 = ABS(MINVAL(TA1))
        v2 = ABS(MAXVAL(TA1))
        KIND_ITP = 1
        IF( v1 <= vref .AND. v2 <= vref ) KIND_ITP = 0
        IF( KIND_ITP == 1 ) CALL PHSGN(nn,Es,TA1)
        TAr2 = 0.D0
        CALL CUBSPL(nn,Es,TA1,n,Enf,TAr2)
!        INTERPOLATING IMAGINARY PART OF T
        TA1(0:nn) = AIMAG(zT_fi(ni+lf*nlm:nf+lf*nlm,n0+l0*nlm,1))
        v1 = ABS(MINVAL(TA1))
        v2 = ABS(MAXVAL(TA1))
        KIND_ITP = 1
        IF( v1 <= vref .AND. v2 <= vref ) KIND_ITP = 0
        IF( KIND_ITP == 1 ) CALL PHSGN(nn,Es,TA1)
        TAi2 = 0.D0
        CALL CUBSPL(nn,Es,TA1,n,Enf,TAi2)
!        PUT TOGETHER INTERPOLATED REAL AND IMAGINARY PARTS OF T
        DO i = 0, n
          zT_Ef(i,il) = CMPLX(TAr2(i),TAi2(i))
        END DO
!        INTERPOLATE MATRIX ELEMENTS OF SCALAR POTENTIAL
        IF( KIND_SCP == 1 ) THEN
          TA1(0:nn) = REAL(zT_fi(ni:nf,il,2))
          TAr2 = 0.D0
          CALL CUBSPL(nn,Es,TA1,n,Enf,TAr2)
          TA1(0:nn) = AIMAG(zT_fi(ni:nf,il,2))
          TAi2 = 0.D0
          CALL CUBSPL(nn,Es,TA1,n,Enf,TAi2)
          DO i = 0, n
            zS_Ef(i,il) = CMPLX(TAr2(i),TAi2(i))
          END DO
        END IF  
        DEALLOCATE( t2i, Es, TA1 )
      END DO

      OPEN( UNIT=35, FILE='CSs/CrossSection_l_All.dat', ACTION='WRITE' )
      OPEN( UNIT=36, FILE='CSs/IonProb_l_All.dat', ACTION='WRITE' )
!      OPEN( UNIT=37, FILE='CSs/FourierTG.dat', ACTION='WRITE' )

      DO i = 0, n
        Ef = Enf(i)
        dE = Ef - E0
        coeff_cs = M_au * c0 * c1 / dE
        IF( ncyc /= 0 ) THEN
          zfw = FTFtw(Eph,dE,ncyc)            !FT TIME PROFILE
          zgw = FTGtw(Eph,dE,ncyc)            !FT DER. TIME PROFILE
          zsw = FTStw(Eph,dE,ncyc)            !FT INT. TIME PROFILE
          coeff_dp = (ABS(zfw)**2)
!          WRITE(37,100) Ef, ABS(zfw)**2, ABS(zgw)**2, ABS(zsw)**2
        ELSE
          coeff_dp = 0.D0
          IF( dE == Eph ) coeff_dp = 1.D0
        END IF
        IF( KIND_SCP == 0 ) THEN
          WRITE(35,100) Ef, (coeff_cs*T2_pts(i,il), il=1,nlm)
          WRITE(36,100) Ef, (coeff_dp*T2_pts(i,il), il=1,nlm)
        ELSE
          zTS = 0.D0
          DO il = 1, nlm
            zTS(il) = zfw * zT_Ef(i,il) + zsw * zS_Ef(i,il)
          END DO
          coeff_cs = coeff_cs / (ABS(zfw)**2)
          coeff_dp = 1.D0
          WRITE(35,100) Ef, (coeff_cs*(ABS(zTS(il))**2), il=1,nlm)
          WRITE(36,100) Ef, (coeff_dp*(ABS(zTS(il))**2), il=1,nlm)
        END IF
      END DO

      CLOSE(35)
      CLOSE(36)
!      CLOSE(37)

100      FORMAT(1000G20.10E3) 

      DEALLOCATE( t2f, t2_pts, TAr2, TAi2 )

      ALLOCATE( zYlm(0:lmax,-lmax:lmax) )
      
      thi = 0.D0
      thf = 2.D0 * PI
      dth = (thf - thi) / DBLE(nthpts)

      Ej = Eref
      IF( nEpts < 0 ) THEN

        OPEN( UNIT=40, FILE='CSs/AngDist_All.dat', ACTION='WRITE' )

        DO i = 0, n

          Ei = Enf(i)

          IF( Ei == Ej ) THEN
            zfw = FTFtw(Eph,Ei-E0,ncyc)
            zgw = FTGtw(Eph,Ei-E0,ncyc)
            zsw = FTStw(Eph,Ef-E0,ncyc)
            kj = SQRT(2.D0 * Ej)
            hj = -1.D0 / kj
            CALL PHACOU(lmf(nlm,1),sigmal,hj)
            DO j = 0, nthpts
              th = thi + DBLE(j) * dth
              zYlm = 0.D0
              CALL Ylm_All(lmax,th,0.D0,zYlm)
              zsum = 0.D0
              DO il = 1, nlm
                lf = lmf(il,1)
                mf = lmf(il,2)
                zph = EXP(zi*(sigmal(lf) - 0.5D0*PI*lf))
                zT = zfw * zT_Ef(i,il) * zYlm(lf,mf) * zph
                zS = 0.D0
                IF( KIND_SCP == 1 ) zS = zsw * zS_Ef(i,il) * zYlm(lf,mf) * zph
!                zsum = zsum + zT_Ef(i,il) * zYlm(lf,mf)
                zsum = zsum + zT + zS
              END DO
              zsum = zsum / zfw
              WRITE(40,100) Enf(i), th, ABS(zsum)**2
            END DO
!            WRITE(30,*)
            Ej = Ej + Eref
          END IF

        END DO

      END IF      

      CLOSE(40)

      DEALLOCATE( zYlm )

      END SUBROUTINE SUM_LF
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE PHSGN(npts,x,y)

!      SUBROUTINE TO CORRECT, IF NEEDED THE SIGN OF THE PHASE FOR THE TRANSITION
!        AMPLITUDES...
!      03-05.01.2018

      USE MOD_TYPES

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: npts
      REAL(DP), INTENT(IN) :: x(0:npts)
      REAL(DP), INTENT(INOUT) :: y(0:npts)
      
      INTEGER :: i, j, jmin, jmax, dj, mref, n, nch, is1, is2
      INTEGER :: n1(0:4)
      REAL(DP) :: x1, x2, xi, dxi, h, d1, d2, yref, yj, dmin, dy, dimin

      INTEGER, ALLOCATABLE :: isgn(:)
      REAL(DP), ALLOCATABLE :: dfx(:)

      dmin = 5.D-3
      dimin = 1.D-5
      dxi = (x(npts) - x(0)) / 4.D0

      n1 = 0

!      WRITE(6,*) 'X limits:', x(0), x(npts)
!      WRITE(6,*) 'Num. Pts:', npts+1

!      COUNT FOR NUMBER OF POINTS IN EACH INTERVAL

      x1 = x(0)
      jmin = 0
      n = 0
      DO i = 1, 4
        x2 = x1 + dxi
        j = jmin
dsrchx: DO
          xi = x(j)
          IF( xi >= x1 .AND. xi <= x2 ) n = n + 1
          j = j + 1
          IF( xi > x2 .OR. j > npts ) THEN
            jmin = j - 1
            EXIT dsrchx
          END IF
        END DO dsrchx
        n1(i) = n - 1
!        WRITE(6,'(A5,I2,X,A4,G17.10,X,A4,G17.10,X,A4,I3)') 'Loop ', i, 'x1 =', x1, 'x2 =', x2, 'ni =', n1(i)-n1(i-1)
        x1 = x2
        x2 = x1 + dxi
      END DO

!      WRITE(6,*)
!      DO i = 1, 4
!        WRITE(6,*) x(n1(i))
!      END DO

!      CALCULATE FIRST DERIVATIVE

!      OPEN( UNIT=50, FILE='DerF.dat', ACTION='WRITE' )

      ALLOCATE( dfx(0:npts) )
      dfx(0) = (y(1) - y(0)) / (x(1) - x(0))
!      WRITE(50,*) x(0), y(0), dfx(0)
      DO i = 1, npts-1
        h = x(i+1) - x(i)
        d1 = (y(i+1) - y(i)) / h
        h = x(i) - x(i-1)
        d2 = (y(i) - y(i-1)) / h
        dfx(i) = (d1 + d2) / 2.D0
!        WRITE(50,*) x(i), y(i), dfx(i)
      END DO
      dfx(npts) = (y(npts) - y(npts-1)) / (x(npts) - x(npts-1))
!      WRITE(50,*) x(npts), y(npts), dfx(npts)
      
!      CLOSE(50)
      
!      ANALIZE FIRST DERIVATIVE IN EACH INTERVAL

      ALLOCATE( isgn(0:npts) )
      DO i = 0, npts
        isgn(i) = 1
        IF( dfx(i) < 0.D0 ) isgn(i) = -1
      END DO

      DO i = 2, 3
      
        jmin = n1(i-1)
        jmax = n1(i)
        mref = isgn(jmin)
        nch = 0
        DO j = jmin+1, jmax
          IF( isgn(j) /= mref ) nch = nch + 1
         END DO
      
        IF( nch > 5 ) THEN
          WRITE(6,*) 'Function oscillating a lot!', i, nch
          OPEN( UNIT=50, FILE='LookAtThis.dat', ACTION='WRITE')
          DO j = 0, npts
            WRITE(50,*) x(j), y(j)
          END DO
!          STOP
        END IF

!        LOOKING THROUGH THE MIDDLE INTERVALS
        IF( nch < 5 ) THEN

          IF( i == 2 ) THEN
            jmin = n1(i-1) - 1
            jmax = n1(i-2)
            dj = -1
          ELSE
                jmin = n1(i) + 1
                jmax = n1(i+1) - 1
                dj = 1
          END IF
          yref = y(jmin-dj)
          
!          LOOKING ONLY IN ONE POINT
          DO j = jmin, jmax, dj
            yj = y(j)
            dy = ABS(yref - yj) / ABS(yref)
            is1 = 1
            is2 = 1
            IF( i == 2 ) THEN
              d1 = yref - yj
              d2 = y(j+2) - y(j+1)
              IF( d1 < 0.D0 ) is1 = -1
              IF( d2 < 0.D0 ) is2 = -1
              IF( ABS(d1) <= dimin ) is1 = 1
              IF( ABS(d2) <= dimin ) is2 = 1
              IF( mref == 1 ) THEN            !FOR INCREASING FUNCTIONS (POSITIVE 1st DER.)
!                IF( yj > yref .AND. dy > dmin .AND. is1 /= is2 ) y(j) = -y(j)
                IF( yj > yref .AND. is1 /= is2 ) y(j) = -y(j)
              ELSE                        !FOR DECREASING FUNCTIONS
!                IF( yj < yref .AND. dy > dmin .AND. is1 /= is2 ) y(j) = -y(j)
                IF( yj > yref .AND. is1 /= is2 ) y(j) = -y(j)
              END IF
            ELSE
              d1 = yj - yref
              d2 = y(j-1) - y(j-2)
              IF( d1 < 0.D0 ) is1 = -1
              IF( d2 < 0.D0 ) is2 = -1
              IF( ABS(d1) <= dimin ) is1 = 1
              IF( ABS(d2) <= dimin ) is2 = 1
              IF( mref == 1 ) THEN            !FOR INCREASING FUNCTIONS (POSITIVE 1st DER.)
!                IF( yj < yref .AND. dy > dmin .AND. isgn(j-dj) /= isgn(j) ) y(j) = -y(j)
                IF( yj < yref .AND. is1 /= is2 ) y(j) = -y(j)
              ELSE                        !FOR DECREASING FUNCTIONS
!                IF( yj > yref .AND. dy > dmin .AND. isgn(j-dj) /= isgn(j) ) y(j) = -y(j)
                IF( yj > yref .AND. is1 /= is2 ) y(j) = -y(j)
              END IF
            END IF
            yref = yj
            mref = is2
          END DO
          
        END IF
      
      END DO

!      OPEN( UNIT=51, FILE='CorrectedFx.dat', ACTION='WRITE' )
!      DO i = 0, npts
!        WRITE(51,*) x(i), y(i)
!      END DO

      DEALLOCATE( isgn, dfx )

      END SUBROUTINE PHSGN
