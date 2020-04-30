      MODULE MOD_TYPES

      IMPLICIT NONE

      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)!4-BYTE INTEGER
      INTEGER, PARAMETER :: DP = KIND(1.0D0)            !REAL(DP) REAL
      INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))      !REAL(DP) COMPLEX

      REAL(DP), PARAMETER :: PI = DACOS(-1.0D0), eps=EPSILON(1.D0)
      COMPLEX(DPC), PARAMETER :: zi = DCMPLX(0.D0,1.D0)

      REAL(DP), PARAMETER :: c_au = 137.03599913815D0, a_au = 5.29177249D-9
      REAL(DP), PARAMETER :: I0_au = 3.50944758D16, E_au = 27.2113962D0
      REAL(DP), PARAMETER :: t_au = 2.41888433D-2, B0_au = 2.35051843D5
      REAL(DP), PARAMETER :: E_eV_J = 1.602176621D-19

      END MODULE MOD_TYPES
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_GRID

      USE MOD_TYPES

      IMPLICIT NONE

      INTEGER :: KIND_GRID, KIND_BC1, KIND_BC2, nbc1, nbc2, ndim
      INTEGER :: io, nintv_exp, nintv_lin, imax
      INTEGER :: lmax, llmax, nrpts, nthpts, nphpts
      REAL(DP) :: ra, rb, rmax, gsize

       REAL(DP), ALLOCATABLE :: xg(:), wg(:)            !GAUSS-LEGENDRE QUADRATURE

!      FOR INTEGRALS ON TH
      INTEGER :: kth, nthfun, nkthp, nointvth, kath
      REAL(DP), ALLOCATABLE :: tht(:), pht(:), xth(:), wth(:), zjt(:), rtot(:), wtot(:)
!      FOR INTEGRALS ON TH, USING THE FIBONACCI NUMERICAL INTEGRATION
      INTEGER :: nfib, nfib0, nfib1
      REAL(DP) :: dzf, dphf

!     OBSERVER POINTS
      INTEGER :: nro, ntho, npho, ithc, iphc
      REAL(DP), ALLOCATABLE :: robs(:), thobs(:), phobs(:)

      END MODULE MOD_GRID
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_BSPLINES

      USE MOD_TYPES

      IMPLICIT NONE

      INTEGER :: k, nfun, nkp, nkpgrid, nointv, nxr, ka
      REAL(DP), ALLOCATABLE :: rt(:), rtf(:), rtk(:)            !GRID FOR B-SPLINES

      REAL(DP), ALLOCATABLE :: Aind(:,:)      !COEFFs FOR CALC DERIV OF B-SPLINES
!      REAL(DP), ALLOCATABLE :: fbsp(:,:,:), dfbsp(:,:,:)

!      FOR INTEGRALS ON TH
      REAL(DP), ALLOCATABLE :: ciLP(:,:,:) !, Sthij(:,:), Afth(:,:,:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zAfth(:,:,:,:,:), zdAfth(:,:,:,:,:)

!      FOR INTEGRALS ON PH
!      REAL(DP), ALLOCATABLE :: Sphij(:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zciexp(:,:)

      CONTAINS
      
        SUBROUTINE BSPALL(r,left,bsp,dbsp,bsp1,bspp)
!        SUBROUTINE TO CALCULATE B-SPLINES AND THEIR DERIVATIVE FOR A GIVEN r

        IMPLICIT NONE

        REAL(DP), INTENT(IN) :: r
        REAL(DP), INTENT(INOUT) :: bsp1(:), bspp(:)

        INTEGER, INTENT(OUT) ::left
        REAL(DP), INTENT(OUT) :: bsp(:), dbsp(:)

        INTEGER :: j, jp, mflag
        REAL(DP) :: A1, A2, b1, b2

        bsp = 0.D0
        bsp1 = 0.D0
        CALL interv(rt,nkp,r,left,mflag)
        CALL bsplvb(nkp,rt,k,1,r,left,bsp)
        CALL bsplvb(nkp,rt,k-1,1,r,left,bsp1)

        bspp = 0.D0
        DO j = 1, k-1
          bspp(j+1) = bsp1(j)
        END DO

        dbsp = 0.D0
        DO j = 1, k
          jp = j + (left-k)
          A1 = 0.D0
          A2 = 0.D0
          IF( jp >= 1 .AND. jp <= nfun ) THEN
            A1 = Aind(jp,1)
            A2 = Aind(jp,2)
          END IF
          b1 = bspp(j)
          b2 = bspp(j+1)
          dbsp(j) = DBLE(k-1) * (A1*b1 - A2*b2)
        END DO
        
        END SUBROUTINE BSPALL
!-------------------------------------------------------------------------------
        SUBROUTINE gauleg(x1,x2,x,w,n)

!           Subroutine GAULEG, given the lower and upper limits
!           of integration x1 and x2 and given n, the routine
!           returns arrays x(1:n) and w(1:n), containing the
!           abcissas and weights for the Gaussian-Legendre n-point
!           quadrature formula.
!           (C) Copr. 1986-92 Numerical Recipes Software @1-. (P. 145)

        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: x1, x2
        double precision, dimension(n), intent(out) :: x, w
        integer :: i, j, m
        double precision :: p1, p2, p3, pp, xl, xm, z, z1, EPS1
!       double precision, parameter :: EPS=3.D-14
!         EPS CHANGED TO MACHINE-SPECIFIC VARIABLE... SEE MOD_TYPES
        EPS1=EPS*10
        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)
        do i=1,m
 	    z=cos(PI*(i-.25d0)/(n+.5d0))
 	    z1 = 0.D0
          do while(abs(z-z1) .gt. EPS1)
            p1=1.d0
	      p2=0.d0
	      do j=1,n
	        p3=p2
	        p2=p1
	        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
            end do
	      pp=n*(z*p1-p2)/(z*z-1.d0)
	      z1=z
	      z=z1-p1/pp
        end do
	  x(i)=xm-xl*z
	  x(n+1-i)=xm+xl*z
	  w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
	  w(n+1-i)=w(i)
      end do
      END SUBROUTINE gauleg

      END MODULE MOD_BSPLINES
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_MATRICES

      USE MOD_TYPES

      IMPLICIT NONE

!     TO INDICATE MATRIX SIZES
      INTEGER :: nmsize, nmsel, nvsize, nvsel

!     DIAGONALIZE HAMILTONIAN
      REAL(DP), ALLOCATABLE :: Hij(:,:), Bij(:,:), En(:), Xij(:,:)
      REAL(DP), ALLOCATABLE :: Sij(:,:), Vij(:,:), Uij(:,:,:), Tij(:,:)

!     WAVE FUNCTIONS
      REAL(DP), ALLOCATABLE :: fur(:,:,:), dfur(:,:,:), Urij(:,:)
      COMPLEX(DPC), ALLOCATABLE :: zYlm(:,:,:), zGlm(:,:,:,:,:)

!     MATRIX ELEMENTS MISCELANEOUS
      REAL(DP), ALLOCATABLE :: Ith(:,:,:), rvecij(:,:,:,:)
      REAL(DP), ALLOCATABLE :: fjrall(:,:,:,:,:), vecU(:,:,:,:), vecQ(:,:,:,:)
      REAL(DP), ALLOCATABLE :: Mijk(:,:,:), Nijk(:,:,:), Sijk(:,:,:)
      REAL(DP), ALLOCATABLE :: Rijl(:,:,:), dRijl(:,:,:)

      COMPLEX(DPC), ALLOCATABLE :: zIth(:,:,:,:,:), zrangij(:,:,:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zPIlm(:,:,:,:,:,:,:,:), zPmq(:,:,:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: psirt(:,:,:,:), zApsirt(:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zArij(:,:,:), Avecr(:,:,:,:,:), zYll(:,:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zArlm(:,:,:,:,:,:,:), zJijq(:,:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zRij(:,:), zYlm_all(:,:,:,:)

!     FORMER MOD_VMV_VARS
      COMPLEX(DPC), ALLOCATABLE :: zx(:), zy(:), zv_temp(:), zA(:,:)

!      COMPLEX(DPC), ALLOCATABLE :: zFPrij(:,:,:,:), zJangimnj(:,:,:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zJAr(:,:,:,:,:,:,:), zJA0imnj(:,:,:,:)

!      REAL(DP), ALLOCATABLE :: Mimnj(:,:,:,:,:), Nimnj(:,:,:,:,:)
!      REAL(DP), ALLOCATABLE :: Jrimnj(:,:,:,:,:), mJrimnj(:,:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zJangimnj(:,:,:,:,:,:), zJangzimnj(:,:,:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zJA0imnj(:,:,:,:), zJA0zimnj(:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zJAr(:,:,:,:,:,:,:), zJArr(:,:,:,:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: zJAimnj(:,:,:,:,:), zJAzimnj(:,:,:,:,:)
!      COMPLEX(DPC), ALLOCATABLE :: Azimnj(:,:,:,:,:,:,:,:,:), Aimnj(:,:,:,:,:,:,:,:,:)

      END MODULE MOD_MATRICES
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_PHOTOION

      USE MOD_TYPES

      IMPLICIT NONE

      INTEGER :: l_ini, l_fin, n0_ini, n1_ini, n0_fin, n1_fin, n1_max, num_i0
      INTEGER :: m_ini, m_fin, nlm, ncomp, nt0, ntf, ntf2, nbounds, nb_max, nb_num
      INTEGER :: KIND_PI, KIND_SCP, KIND_POT, KIND_ENV, KIND_TD, KIND_RK
      INTEGER :: KIND_VEC, KIND_EGR, KIND_NLM
      INTEGER :: moam, mb, mph, np, nEpts, ncyc, nm, nfields
      REAL(DP) :: Emax_fin, Eph, kph, kvec_z, A0, I0, A01, I01, w0, wb, b0, bb
      REAL(DP) :: afocus, Ephb, kphb, Eref, qvecr, qvecz, Zatom, bx, B0z
      REAL(DP) :: A0x, A0y, A0z

      INTEGER, ALLOCATABLE :: n01(:,:), lmf(:,:), ist0(:)
      REAL(DP), ALLOCATABLE :: rij(:,:,:), ci_ini(:), ci_fin(:,:)
      REAL(DP), ALLOCATABLE :: Enl(:,:), cinl(:,:,:), rEki(:,:)
      REAL(DP), ALLOCATABLE :: T_fi(:,:), T2fi(:,:), S2fi(:,:)
      REAL(DP), ALLOCATABLE :: E_ini(:), E_fin(:), Enf(:), rEkf(:)
      COMPLEX(DPC), ALLOCATABLE :: zAij(:,:,:,:,:), zT_fi(:,:,:)
      COMPLEX(DPC), ALLOCATABLE :: zT_Ef(:,:), zS_Ef(:,:), zVrij(:,:)
      COMPLEX(DPC), ALLOCATABLE :: zHint_ij(:,:,:)

!     ROGERS POT.
      INTEGER :: Numn(3), Ntot
      REAL(DP) :: aj(3,0:3), alphan(3)

!     SIMONS-FUES POT.
      REAL(DP), ALLOCATABLE :: Bl(:)

!     TIME PROPAGATION
      INTEGER :: nvec, nbra, nket, nt, ncyc2, ntsteps, nbvec, nsel
      REAL(DP) :: Tpulse, Tpulse2, t_ini, t_fin, tiau, tfau, Tpump, Ef_min, Ef_max
      REAL(DP) :: Epump, Eprobe, Tprobe, t_delay, t2iau, t2fau, td1, td2, Eph2

      INTEGER, ALLOCATABLE :: nl_bra(:,:), nl_ket(:,:), nl_sel(:,:)
      REAL(DP), ALLOCATABLE :: Eall(:), tall(:), Eball(:)
      COMPLEX(DPC), ALLOCATABLE :: zVtij(:,:), zf(:), zdfdt(:)
      COMPLEX(DPC), ALLOCATABLE :: zTor_ij(:,:,:), zcit(:,:)

      CONTAINS

        FUNCTION DELIJ(i,j)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: i, j

        REAL(DP) :: DELIJ

        DELIJ = 0.D0
        IF( i == j ) DELIJ = 1.D0

        END FUNCTION DELIJ
!-------------------------------------------------------------------------------
        FUNCTION SELPOT(r)

        IMPLICIT NONE

        REAL(DP), INTENT(IN) :: r
        REAL(DP) :: SELPOT

        INTEGER :: i, ni
        REAL(DP) :: Vr

        IF( KIND_POT == 0 ) THEN      !COULOMB POTENTIAL

          Vr = -Zatom / r

        ELSE IF( KIND_POT == 1 ) THEN      !ROGERS POTENTIAL. COEFFICIENTS FOR Ca+

          Vr = 0.D0
          ni = 0
          DO i = 1, 3
            ni = Numn(i)
            Vr = Vr + ni * EXP(-alphan(i)*r)
          END DO
          Vr = -1.D0*( Zatom - Ntot + Vr ) / r

        ELSE IF( KIND_POT == 2 ) THEN

          Vr = -Zatom / r

        END IF

        SELPOT = Vr

        END FUNCTION SELPOT
!-------------------------------------------------------------------------------
        SUBROUTINE SEL_T0_STATES
        
        IMPLICIT NONE
        
        INTEGER :: i, j, n0, l0, m0, ni, li, mi
        
        l0 = lmf(0,1)
        m0 = lmf(0,2)
        n0 = n0_ini - l0
        
        num_i0 = 1
        IF( KIND_NLM == 1 ) num_i0 = 2*l0 + 1
      
        ALLOCATE( ist0(num_i0) )

        ist0 = 0
        j = 0
        DO i = 1, nvec
          ni = nl_bra(i,1)
          li = nl_bra(i,2)
          mi = nl_bra(i,3)
          IF( ni == n0 .AND. li == l0 ) THEN
            IF( KIND_NLM == 0 .AND. mi == m0 ) THEN
              ist0(1) = i
            ELSE IF( KIND_NLM == 1 .AND. ABS(mi) <= l0 ) THEN
              j = j + 1
              ist0(j) = i
            END IF
          END IF
        END DO
        
        END SUBROUTINE SEL_T0_STATES
!-------------------------------------------------------------------------------
        FUNCTION CHAMP(ID,t,phi)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: ID                  !PUMP OR PROBE
        REAL(DP), INTENT(IN) :: t, phi

        REAL(DP) :: tp, pht, t2, t2p
        COMPLEX(DPC) :: CHAMP, zft1, zft2

        tp = t - tiau
        t2p = -10.D0
        IF( Eprobe /= 0.D0 ) t2p = t - t2iau

        zft1 = 0.D0
        zft2 = 0.D0
        CHAMP = 0.D0

        IF( KIND_ENV == 0 ) THEN            !SIN(w0 t) - CONTINUOUS RADIATION

          zft1 = 1.D0

        ELSE IF( KIND_ENV == 1 ) THEN      !COS^2()

          pht = PI * tp / Tpulse

          IF( t >= tiau .AND. t <= tfau  ) THEN
                IF( KIND_PI == 2 ) THEN
              zft1 = COS(PI*t/Tpump)**2
              zft2 = 0.D0
            ELSE
              zft1 = (2.D0*PI / (Tpump*Eph)) * COS(pht) * SIN(pht)
              zft2 = COS(pht)**2
            END IF
          END IF

        ELSE IF( KIND_ENV == 2 ) THEN      !SIN^2()
         
          IF( tp >= 0.D0 .AND. tp <= Tpump ) THEN
            IF( ID == 1 ) THEN
              zft1 = (SIN(PI*tp/Tpump)**2) * SIN(Eph*(tp-td1)) - &
     &                SIN(PI*tp/Tpump) * COS(PI*tp/Tpump) * COS(Eph*(tp-td1))
            ELSE IF( ID == 3 ) THEN
              zft1 = (SIN(PI*tp/Tpump)**2) * COS(Eph*(tp-td1)) / Eph
            END IF
          END IF
          IF( t2p >= 0.D0 .AND. t2p <= Tprobe ) THEN
            IF( ID == 2 ) THEN
              zft2 = (SIN(PI*t2p/Tprobe)**2) * SIN(Eph2*(t2p-td2) + phi) - &
     &                SIN(PI*t2p/Tprobe) * COS(PI*t2p/Tprobe) * COS(Eph2*(t2p-td2) + phi)
            ELSE IF( ID == 4 ) THEN
              zft2 = (SIN(PI*t2p/Tprobe)**2) * COS(Eph2*(t2p-td2) + phi) / Eph2
            END IF
          END IF

          IF( ID == 1 .OR. ID == 3 ) zft2 = 0.D0
          IF( ID == 2 .OR. ID == 4 ) zft1 = 0.D0
          
        ELSE IF( KIND_ENV == 3 ) THEN      !GAUSSIAN
        
          zft1 = EXP(-0.5D0*((t/Tpulse)**2))
          
        END IF

        CHAMP = Epump * zft1 + Eprobe * zft2
        
        END FUNCTION CHAMP
!-------------------------------------------------------------------------------
        SUBROUTINE ZHVMV(UPLO,N,zx,zy,zv_temp,zA,zf)
       
!        SUBROUTINE TO CALCULATE VECTOR - HERMITIAN MATRIX - VECTOR MULTIPLICATIONS

!        USE MOD_VMV_VARS
        
        IMPLICIT NONE
        
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: N !, ni, li, nf, lf, il, jl, indx
        COMPLEX(DPC), INTENT(IN) :: zx(:), zy(:), zA(:,:)
        COMPLEX(DPC), INTENT(INOUT) :: zv_temp(:)
        COMPLEX(DPC), INTENT(OUT) :: zf

        INTEGER :: INCX, INCY
        COMPLEX(DPC) :: zALPHA, zBETA, ZDOTU

        zALPHA = 1.D0
        zBETA = 0.D0
        INCX = 1
        INCY = 1

        zv_temp = 0.D0

        CALL ZHEMV(UPLO,N,zALPHA,zA,N,zx,INCX,zBETA,zv_temp,INCY)
        zf = ZDOTU(N,zy,INCX,zv_temp,INCY)
        
        END SUBROUTINE ZHVMV
!-------------------------------------------------------------------------------
        SUBROUTINE DSVMV(UPLO,N,x,y,v_temp,A,f)

!        SUBROUTINE TO CALCULATE VECTOR - SYMMETRIC MATRIX - VECTOR MULTIPLICATIONS

        IMPLICIT NONE
        
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: N
        REAL(DP), INTENT(IN) :: x(N), y(N), A(N,N)
        REAL(DP), INTENT(INOUT) :: v_temp(N)
        REAL(DP), INTENT(OUT) :: f

        INTEGER :: INCX, INCY
        REAL(DP) :: ALPHA, BETA, DDOT

        ALPHA = 1.D0
        BETA = 0.D0
        INCX = 1
        INCY = 1

        v_temp = 0.D0

        CALL DSYMV(UPLO,N,ALPHA,A,N,x,INCX,BETA,v_temp,INCY)
        f = DDOT(N,y,INCX,v_temp,INCY)
        
        END SUBROUTINE DSVMV
!-------------------------------------------------------------------------------
        SUBROUTINE ZGCVMV(TRANS,M,N,NRHS,zx,zA,zy,zf)
        
!        SUBROUTINE TO CALCULATE COMPLEX VECTOR - HERMITIAN MATRIX - VECTOR MULTIPLICATIONS
        
        IMPLICIT NONE
        
        CHARACTER(LEN=1), INTENT(IN) :: TRANS
        INTEGER, INTENT(IN) :: M, N, NRHS
        COMPLEX(DPC), INTENT(IN) :: zx(N), zA(M,N,NRHS), zy(M)
        COMPLEX(DPC), INTENT(OUT) :: zf(NRHS)

        INTEGER :: i, INCX, INCY
        COMPLEX(DPC) :: zALPHA, zBETA, ZDOTC
        COMPLEX(DPC), ALLOCATABLE :: zw(:), zB(:,:)

        ALLOCATE( zw(M), zB(M,N) )

        zALPHA = 1.D0
        zBETA = 0.D0
        INCX = 1
        INCY = 1

        DO i = 1, NRHS
          zB(:,:) = zA(:,:,i)
          CALL ZGEMV(TRANS,M,N,zALPHA,zB,N,zx,INCX,zBETA,zw,INCY)
          zf(i) = ZDOTC(N,zy,INCX,zw,INCY)
        END DO
      
        DEALLOCATE( zw, zB )
      
        END SUBROUTINE ZGCVMV
!-------------------------------------------------------------------------------
        SUBROUTINE ZVMV(TRANS,M,N,NRHS,INDX,zx,zA,zy,zf)
        
!        SUBROUTINE TO CALCULATE COMPLEX VECTOR - HERMITIAN MATRIX - VECTOR MULTIPLICATIONS
        
        IMPLICIT NONE
        
        CHARACTER(LEN=1), INTENT(IN) :: TRANS
        INTEGER, INTENT(IN) :: M, N, NRHS, INDX
        COMPLEX(DPC), INTENT(IN) :: zx(N), zA(M,N,NRHS), zy(M)
        COMPLEX(DPC), INTENT(OUT) :: zf(NRHS)

        INTEGER :: i, INCX, INCY
        COMPLEX(DPC) :: zALPHA, zBETA, ZDOTC
        COMPLEX(DPC), ALLOCATABLE :: zu(:), zv(:), zw(:), zB(:,:)

        ALLOCATE( zu(M), zv(M), zw(M), zB(M,N) )

        zALPHA = 1.D0
        zBETA = 0.D0
        INCX = 1
        INCY = 1

        IF( INDX == 1 ) THEN
          zv(:) = zx(:)
          zu(:) = CMPLX(zy(:))
        ELSE IF( INDX == 2 ) THEN
          zv(:) = CMPLX(zx(:))
          zu(:) = zy(:)
        END IF

        DO i = 1, NRHS
          zB(:,:) = zA(:,:,i)
          CALL ZGEMV(TRANS,M,N,zALPHA,zB,N,zv,INCX,zBETA,zw,INCY)
          zf(i) = ZDOTC(N,zu,INCX,zw,INCY)
        END DO
      
        DEALLOCATE( zu, zv, zw, zB )
      
        END SUBROUTINE ZVMV

      END MODULE MOD_PHOTOION
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_FAC_VARS

      USE MOD_TYPES

      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: fac(:)

      END MODULE MOD_FAC_VARS
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_TD_VARS

      USE MOD_TYPES
      
      IMPLICIT NONE

!      LAPACK/BLAS VARIABLES

      CHARACTER :: TRANS, UPLO
      INTEGER :: Ni, Nj, LDA, INCX, INCY
      COMPLEX(DPC) :: zALPHA, zBETA
      COMPLEX(DPC), ALLOCATABLE :: zxx(:), zyy(:)

      END MODULE MOD_TD_VARS
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_RK_PARAMS
      
      USE MOD_TYPES
      
      IMPLICIT NONE
      
!      DEFINES THE BUTCHER ARRAYS OF THE RUNGE_KUTTA METHOD
      INTEGER, PARAMETER :: p=4

      REAL(DP), PARAMETER :: a21=2.d0/9.d0, a31=1.d0/12.d0, a32=1.d0/4.d0, &
     &          a41=69.d0/128.d0, a42=-243.d0/128.d0, a43=135.d0/64.d0, &
     &          a51=-17.d0/12.d0, a52=27.d0/4.d0, a53=-27.d0/5.0d0, &
     &          a54=16.d0/15.d0, a61=65.d0/432.d0, a62=-5.d0/16.d0, &
     &          a63=13.d0/16.d0, a64=4.d0/27.d0, a65=5.d0/144.d0

      REAL(DP), PARAMETER :: b1=1.d0/9.d0, b2=0.d0, b3=9.d0/20.d0, b4=16.d0/45.d0, &
     &                           b5=1.d0/12.d0
     
      REAL(DP), PARAMETER :: c1=0.0d0, c2=2.d0/9.d0, c3=1.d0/3.d0, c4=3.0d0/4.d0, &
     &                       c5=1.d0, c6=5.d0/6.d0

      REAL(DP), PARAMETER :: d1=47.d0/450.d0, d2=0.d0, d3=12.d0/25.d0, &
     &                       d4=32.d0/225.d0, d5=1.d0/30.d0, d6=6.d0/25.d0

      REAL(DP), PARAMETER :: e1=-1.d0/150.d0, e2=0.d0, e3=3.d0/100.d0, &
     &          e4=-16.d0/75.d0, e5=-1.d0/20.d0, e6=6.d0/25.d0

      END MODULE MOD_RK_PARAMS
!
!-------------------------------------------------------------------------------
!
      MODULE MOD_FTW_PH
      
      USE MOD_TYPES
      
      CONTAINS
      
        FUNCTION FTFtw(w0,w,n)

!       FOURIER TRANSFORM OF A PULSE WITH SIN^2 ENVELOPE

        IMPLICIT NONE
      
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: w0, w
      
        COMPLEX(DPC) :: FTFtw
      
        REAL(DP) :: b, dw, c0, c1, c2, ph
      
        c0 = -SQRT(2.D0 / PI)
        b = w0 / (2.D0 * n)
        dw = w - w0

        c1 = 2.D0 * (b**2)
        c2 = dw * (dw**2 - 4.D0*(b**2))
        ph = PI * dw / (2.D0 * b)

        IF( dw /= 0.D0 ) THEN
          FTFtw = -c0 * (c1 / c2) * SIN(ph)
        ELSE
          FTFtw = SQRT(PI / 2.D0) / (2.D0 * b)
        END IF

        END FUNCTION FTFtw
!-------------------------------------------------------------------------------
        FUNCTION FTGtw(w0,w,n)
        
!       FOURIER TRANSFORM OF A PULSE WITH SIN^2 ENVELOPE

        IMPLICIT NONE
      
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: w0, w
      
        COMPLEX(DPC) :: FTGtw
      
        REAL(DP) :: b, dw, c0, c1, c2, ph
      
        c0 = SQRT(2.D0 / PI)
        b = w0 / (2.D0*n)
        dw = w - w0

        c1 = 2.D0 * w * (b**2)
        c2 = dw * (dw**2 - 4.D0*(b**2))
        ph = PI * dw / (2.D0 * b)

        IF( dw /= 0.D0 ) THEN
          FTGtw = zi * c0 * (c1 / c2) * SIN(ph)
        ELSE
          FTGtw = -zi * SQRT(PI / 2.D0) * w / (2.D0 * b)
        END IF

        END FUNCTION FTGtw
!-------------------------------------------------------------------------------
        FUNCTION FTStw(w0,w,n)

        IMPLICIT NONE
      
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: w0, w
      
        COMPLEX(DPC) :: FTStw
      
        REAL(DP) :: b, dw, c0, c1, c2, c3, c4, ph1, ph2
      
        b = w0 / (2.D0*n)
        dw = w - w0

        c0 = SQRT(2.D0 / PI) / (4.D0 * b * n * ((n**2)-1.D0))

        c1 = 4.D0 * (b**2) * (n**2 - 1.D0) + dw * (dw - n*b)
        c2 = dw * (dw**2 - 4.D0*(b**2))

        c3 = ((-1.D0)**n) * (2.D0*(n**2) - 1.D0)
        c4 = w

        ph1 = PI * dw / (2.D0 * b)
        ph2 = PI * w / (2.D0 * b)

        IF( dw /= 0.D0) THEN
          FTStw = zi * c0 * (-(c1 / c2) * SIN(ph1) + (c3 / c4) * SIN(ph2))
        ELSE
          c1 = PI * ((n**2) - 1.D0)
          c2 = 2.D0 * b
          FTStw = zi * c0 * ((c1 / c2) + (c3 / c4) * SIN(ph2))
        END IF

        END FUNCTION FTStw
!-------------------------------------------------------------------------------
        FUNCTION ZDFT(nt,t,zft,nw,w,zfw)

!       FOURIER TRANSFORM OF AN ARBITRARY COMPLEX FUNCTION

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nt, nw
        REAL(DP), INTENT(IN) :: t(0:nt), w(0:nw)
        COMPLEX(DPC), INTENT(INOUT) :: zft(0:nt), zfw(0:nw)
        
        INTEGER :: iw, it, ZDFT
        REAL(DP) :: wi, dw, t0, t1, t2
        COMPLEX(DPC) :: zsumfw0, zsumfw1, zsumfw2, zft0, zft1, zft2, zf

        dw = w(2) - w(1)

        DO iw = 0, nw

          wi = w(iw)

          zsumfw0 = 0.D0
          zsumfw1 = 0.D0
          zsumfw2 = 0.D0
          DO it = 1, nt-1
          
            t0 = t(it-1)
            t1 = t(it)
            t2 = t(it+1)
          
            zft0 = zft(it-1) * EXP(-zi*wi*t0)
            zft1 = zft(it) * EXP(-zi*wi*t1)
            zft2 = zft(it+1) * EXP(-zi*wi*t2)
            
            zsumfw0 = zsumfw0 + zft0
            zsumfw1 = zsumfw1 + zft1
            zsumfw2 = zsumfw2 + zft2
            
          END DO

          zf = (zsumfw0 + 4.D0*zsumfw1 + zsumfw2) / 12.D0
          
          IF( MOD(nt+1,2) == 0 ) THEN
          
            t0 = t(nt)
            t1 = t(nt-1)
            t2 = t(nt-2)
          
            zft0 = zft(nt) * EXP(-zi*wi*t0)
            zft1 = zft(nt-1) * EXP(-zi*wi*t1)
            zft2 = zft(nt-2) * EXP(-zi*wi*t2)
            
            zsumfw0 = zsumfw0 + zft0
            zsumfw1 = zsumfw1 + zft1
            zsumfw2 = zsumfw2 + zft2
        
            zf = zf + (5.D0*zsumfw0 + 8.D0*zsumfw1 - zsumfw2) * 12.D0
          END IF

          zfw(iw) = zf * dw
          
        END DO

        ZDFT = 1

        END FUNCTION ZDFT
!-------------------------------------------------------------------------------
        SUBROUTINE PHACOU(LM,PCOUL,ETA)!SUBROUTINE TO CALCULATE COULOMB PH. SHIFT

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LM
        REAL(DP), INTENT(IN) :: ETA
        REAL(DP), INTENT(OUT) :: PCOUL(0:LM)
        COMPLEX(DPC) :: Z, AL, Z2

        INTEGER :: L

        Z = CMPLX(DBLE(LM+1),ETA)
        Z2 = Z * Z
        AL = (((-1.D0/(1680.D0*Z2)+1.D0/1260.D0)/Z2-1.D0/360.D0)/Z2+1.D0/12.D0)/Z
        AL = AL + (Z-0.5)*LOG(Z) - Z + 0.5D0*LOG(2.D0*PI)
        PCOUL(LM) = AIMAG(AL)
        DO L = 1, LM
          PCOUL(LM-L) = PCOUL(LM-L+1) - ATAN2(ETA,DBLE(LM+1-L))
        END DO

        END SUBROUTINE PHACOU
!-------------------------------------------------------------------------------
        SUBROUTINE DSIMPINT(n,f,g,h,valI1,valI2,valI3)
        
!        INTEGRATE THE FUNCTION f, USING SIMPSON INTEGRAL METHOD

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: f(0:n), g(0:n), h(0:n)
        REAL(DP), INTENT(OUT) :: valI1, valI2, valI3

        INTEGER :: i
        REAL(DP) :: dx, f0, f1, f2, g0, g1, g2, h0, h1, h2, ffactor, gfactor, hfactor
        REAL(DP) :: sumf0, sumf1, sumf2, sumg0, sumg1, sumg2, sumh0, sumh1, sumh2

        dx = 1.D0 / DBLE(n)

        sumf0 = 0.D0
        sumf1 = 0.D0
        sumf2 = 0.D0
        sumg0 = 0.D0
        sumg1 = 0.D0
        sumg2 = 0.D0
        sumh0 = 0.D0
        sumh1 = 0.D0
        sumh2 = 0.D0
        DO i = 1, n-1, 2

          f0 = f(i-1)
          f1 = f(i)
          f2 = f(i+1)
          
          g0 = g(i-1)
          g1 = g(i)
          g2 = g(i+1)
          
          h0 = h(i-1)
          h1 = h(i)
          h2 = h(i+1)
          
          
          sumf0 = sumf0 + f0
          sumf1 = sumf1 + f1
          sumf2 = sumf2 + f2
          
          sumg0 = sumg0 + g0
          sumg1 = sumg1 + g1
          sumg2 = sumg2 + g2

          sumh0 = sumh0 + h0
          sumh1 = sumh1 + h1
          sumh2 = sumh2 + h2

        END DO
        
        ffactor = (sumf0 + 4.D0 * sumf1 + sumf2) / 3.D0
        gfactor = (sumg0 + 4.D0 * sumg1 + sumg2) / 3.D0
        hfactor = (sumh0 + 4.D0 * sumh1 + sumh2) / 3.D0
        
        IF( MOD(n+1,2) == 0 ) THEN
        
          f0 = f(n)
          f1 = f(n-1)
          f2 = f(n-2)
          
          g0 = g(n)
          g1 = g(n-1)
          g2 = g(n-2)

          h0 = h(n)
          h1 = h(n-1)
          h2 = h(n-2)

          ffactor = ffactor + (5.D0 * f0 + 8.D0 * f1 - f2) / 12.D0
          gfactor = gfactor + (5.D0 * g0 + 8.D0 * g1 - g2) / 12.D0
          hfactor = hfactor + (5.D0 * h0 + 8.D0 * h1 - h2) / 12.D0
          
        END IF
        
        valI1 = dx * ffactor
        valI2 = dx * gfactor
        valI3 = dx * hfactor

        END SUBROUTINE DSIMPINT
!-------------------------------------------------------------------------------
        SUBROUTINE ZSIMPINT(n,zf,zg,zh,zvalI1,zvalI2,zvalI3)
        
!        INTEGRATE THE FUNCTION f, USING SIMPSON INTEGRAL METHOD

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        COMPLEX(DPC), INTENT(IN) :: zf(n), zg(n), zh(n)
        COMPLEX(DPC), INTENT(OUT) :: zvalI1, zvalI2, zvalI3

        INTEGER :: i
        REAL(DP) :: dx
        COMPLEX(DPC) :: zf0, zf1, zf2, zg0, zg1, zg2, zh0, zh1, zh2
        COMPLEX(DPC) :: zffactor, zgfactor, zhfactor
        COMPLEX(DPC) :: zsumf0, zsumf1, zsumf2, zsumg0, zsumg1, zsumg2, zsumh0, zsumh1, zsumh2

        dx = 1.D0 / DBLE(n-1)

        zsumf0 = 0.D0
        zsumf1 = 0.D0
        zsumf2 = 0.D0
        zsumg0 = 0.D0
        zsumg1 = 0.D0
        zsumg2 = 0.D0
        zsumh0 = 0.D0
        zsumh1 = 0.D0
        zsumh2 = 0.D0
        DO i = 2, n-1, 2

          zf0 = zf(i-1)
          zf1 = zf(i)
          zf2 = zf(i+1)
          
          zg0 = zg(i-1)
          zg1 = zg(i)
          zg2 = zg(i+1)
          
          zh0 = zh(i-1)
          zh1 = zh(i)
          zh2 = zh(i+1)
          
          
          zsumf0 = zsumf0 + zf0
          zsumf1 = zsumf1 + zf1
          zsumf2 = zsumf2 + zf2
          
          zsumg0 = zsumg0 + zg0
          zsumg1 = zsumg1 + zg1
          zsumg2 = zsumg2 + zg2

          zsumh0 = zsumh0 + zh0
          zsumh1 = zsumh1 + zh1
          zsumh2 = zsumh2 + zh2

        END DO
        
        zffactor = (zsumf0 + 4.D0 * zsumf1 + zsumf2) / 3.D0
        zgfactor = (zsumg0 + 4.D0 * zsumg1 + zsumg2) / 3.D0
        zhfactor = (zsumh0 + 4.D0 * zsumh1 + zsumh2) / 3.D0
        
        IF( MOD(n,2) == 0 ) THEN
        
          zf0 = zf(n)
          zf1 = zf(n-1)
          zf2 = zf(n-2)
          
          zg0 = zg(n)
          zg1 = zg(n-1)
          zg2 = zg(n-2)

          zh0 = zh(n)
          zh1 = zh(n-1)
          zh2 = zh(n-2)

          zffactor = zffactor + (5.D0 * zf0 + 8.D0 * zf1 - zf2) / 12.D0
          zgfactor = zgfactor + (5.D0 * zg0 + 8.D0 * zg1 - zg2) / 12.D0
          zhfactor = zhfactor + (5.D0 * zh0 + 8.D0 * zh1 - zh2) / 12.D0
          
        END IF
        
        zvalI1 = dx * zffactor
        zvalI2 = dx * zgfactor
        zvalI3 = dx * zhfactor

        END SUBROUTINE ZSIMPINT
!-------------------------------------------------------------------------------
        FUNCTION FIBONACCI(indx)

        USE MOD_TYPES
      
        IMPLICIT NONE
      
        INTEGER, INTENT(IN) :: indx
        INTEGER :: FIBONACCI

        INTEGER :: i, i0, i1

        FIBONACCI = 0
        IF( indx == 0 ) THEN
          FIBONACCI = 0
        ELSE IF( indx == 1 .OR. indx == 2 ) THEN
          FIBONACCI = 1
        ELSE
          i0 = 1
          i1 = 1
          DO i = 3, indx
            FIBONACCI = i0 + i1
            i0 = i1
            i1 = FIBONACCI
          END DO
        END IF
      
        END FUNCTION FIBONACCI
!-------------------------------------------------------------------------------
        SUBROUTINE FIBINT(n,m,zf,zInt)
      
        USE MOD_TYPES
      
        IMPLICIT NONE
      
        INTEGER, INTENT(IN) :: n, m
        COMPLEX(DPC), INTENT(IN) :: zf(0:m)

        COMPLEX(DPC), INTENT(OUT) :: zInt

        INTEGER :: i, j, m0
        REAL(DP) :: z, dz, dzj
        COMPLEX(DPC) :: zsum

        m0 = m
      
        dz = 2.D0 / DBLE(m0)

        zInt = 0.D0

        zsum = 0.D0
        DO i = 0, m0
          z = -1.D0 + DBLE(i)*dz
          dzj = 1.D0 + COS(PI*z)
          zsum = zsum + dzj * zf(i)
        END DO

        zInt = 2.D0 * PI * dz * zsum

        END SUBROUTINE FIBINT
!-------------------------------------------------------------------------------
        FUNCTION DELIJ(i,j)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: i, j
        REAL(DP) :: DELIJ
        
        DELIJ = 0.D0
        IF( i == j ) DELIJ = 1.D0
        
        END FUNCTION DELIJ
!-------------------------------------------------------------------------------
        FUNCTION Plmr(l1,l2,l3,lp,l,m1,m2,m3,m)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: l1, l2, l3, lp, l, m1, m2, m3, m

        INTEGER :: mp
        REAL(DP) :: coeff, T0, T1, T2, sumTj, THREE_J, Plmr

        coeff = DBLE((2*l1+1) * (2*lp+1) * (2*l+1))
        coeff = SQRT(coeff / (4.D0*PI))

        T0 = THREE_J(l1,lp,l,0,0,0)

        sumTj = 0.D0
        IF( T0 /= 0.D0 ) THEN
          DO mp = -lp, lp
            T1 = THREE_J(lp,l2,l3,mp,m2,m3)
            T2 = THREE_J(l1,lp,l,m1,mp,-m)
            sumTj = sumTj + T1 * T2
          END DO
        END IF

        Plmr = ((-1.D0)**m) * coeff * T0 * sumTj

        END FUNCTION Plmr
!-------------------------------------------------------------------------------
        FUNCTION Glmr(l1,l2,l,m1,m2,m)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: l1, l2, l, m1, m2, m

        REAL(DP) :: coeff, T0, T1, THREE_J, Glmr

        coeff = DBLE((2*l1+1) * (2*l2+1) * (2*l+1))
        coeff = SQRT(coeff / (4.D0*PI))

        T0 = THREE_J(l1,l2,l,0,0,0)
        T1 = THREE_J(l1,l2,l,-m1,m2,-m)

        Glmr = ((-1.D0)**(m1+m)) * coeff * T0 * T1

        END FUNCTION Glmr
!-------------------------------------------------------------------------------
        FUNCTION Plmq(l1,l2,l3,lp,m1,m2,m3)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: l1, l2, l3, lp, m1, m2, m3

        INTEGER :: mp
        REAL(DP) :: coeff, T0, THREE_J, Plmq

        mp = -m1
        coeff = ((-1.D0)**m1) * DELIJ(l1,lp)

        T0 = THREE_J(lp,l2,l3,mp,m2,m3)
        Plmq = coeff * T0

        END FUNCTION Plmq

      END MODULE MOD_FTW_PH
