      SUBROUTINE READ_INPUTS

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_PHOTOION
      USE MOD_FTW_PH

      IMPLICIT NONE

      INTEGER :: i, j, il, jl
      REAL(DP) :: dx, rimax, suman, xn

!      DECLARATION NAME LISTS FOR INPUT FILE
      NAMELIST / VARS_BSP / KIND_GRID, ra, rb, rmax, k, ka, nfun, KIND_BC1, KIND_BC2, nfib
      NAMELIST / VARS_TISE / n0_ini, l_ini, m_ini, l_fin, lmax, Emax_fin, Zatom, &
     &                        KIND_POT, KIND_EGR, KIND_NLM
      NAMELIST / VARS_FIELD / KIND_PI, KIND_SCP, KIND_TD, KIND_ENV, KIND_RK,  &
     &                         KIND_VEC, A0, w0, Eph, ncyc, Eph2, ncyc2, moam, &
     &                         mph, I0, I01, b0, afocus, nEpts, nthpts, nphpts,&
     &                         Eref, bx, B0z, A01, t_delay, A0x, A0y, A0z

!      READ INPUT DATA

!      READING CONFIGURATION FOR B-SPLINE BASIS SET
!        PREPARING VARIABLES FIRST
      KIND_GRID = 0      !KNOT-PT: 0:LIN, 1:EXP, 2:EXP-LIN
      ra = 0.D0            !INITAL POINT FOR THE BOX
      rb = 0.D0            !FINAL POINT OF THE BOX
      k = 0                  !ORDER OF THE BSPLINES
      ka = 0            !ORDER OF THE GAUSS-LEGENDRE POLYNOMIAL...
      nfun = 0            !NUMBER OF B-SPLINE FUNCTIONS
      nthfun = 0            !NUM. B-SPLINES TO EXPAND FUNCTIONS ON th AND ph
      KIND_BC1 = 0      !TO INCLUDE (1) OR NOT (0) THE VERY 1st B-SPLINE
      KIND_BC2 = 0      !... THE LAST B-SPLINE
      nfib = 1            !FIBONACCI NUMERICAL INTEGRATION... NUMBER OF POINTS
      READ(5,VARS_BSP)

      IF( ka == 0 ) ka = k + 3
      IF( nthfun == 0 ) nthfun = 20

      nbc1 = k
      nbc2 = k
      IF( KIND_BC1 == 0 ) nbc1 = k - 1
      IF( KIND_BC2 == 0 ) nbc2 = k - 1

      nkp = nfun + k
      nointv = nkp - nbc1 - nbc2 + 1

      gsize = rb - ra

      IF( KIND_GRID == 2 ) THEN      !EXPONENTIAL-LINEAR GRID

        WRITE(6,'(/,A29,I5)') 'Initial Number of Functions: ', nfun

        dx = gsize / nointv
        rimax = (rmax - ra) / dx
        imax = NINT(rimax)

        nintv_exp = 3 * imax
        nintv_lin = nointv - imax

        nointv = nintv_exp + nintv_lin
        nkp = nointv + nbc1 + nbc2 - 1
        nfun = nkp - k

        WRITE(6,'(A29,I5)') 'Number of functions changed: ', nfun

      END IF

      WRITE(6,'(A40,I5)') 'Number of B-spline Functions / l: nfun =', nfun

!      READING VARIABLES TO SET THE SCHRÖDINGER EQUATION TO SOLVE

      KIND_POT = 0                  !TO SELECT COULOMB POT.
      n0_ini = 1                    !INITIAL STATE n
      l_ini = 0                     !... l
      m_ini = 0                     !... m
      l_fin = 0                     !FINAL STATE l
      lmax = 0                      !MAX. VALUE OF l ...
      Emax_fin = -1.D0              !UPPER LIMIT FOR CROSS SECTIONS
      Zatom = 1.D0                  !NUCLEAR CHARGE
      KIND_EGR = 0                  !0: EIGENSPECTRUM, 1: FROM A SELECTED NUM. OF STATES
      KIND_NLM = 0                  !POLARIZED/UNPOLARIZED INITIAL STATE
      READ(5,VARS_TISE)

      IF( l_fin > lmax ) lmax = l_fin
!      IF( n0_ini < l_ini + 1 ) THEN
!        WRITE(6,*) 'Correcting to a Physical initial Quantum Principal Number'
!        n0_ini = l_ini + 1
!      END IF

      WRITE(6,'(/,A38,I3)') 'Max. Angular Momenta Included: l_max =', lmax

      IF( KIND_POT == 1 ) THEN            !ROGERS POT., FOR Ca+
      
        Numn(1) = 2
        Numn(2) = 8
        Numn(3) = 8

        aj = 0.D0

        aj(1,0) = 0.8855D0
        aj(1,1) = 0.2549D0
        aj(1,2) = -0.0901D0

        aj(2,0) = 0.3386D0
        aj(2,1) = 1.1323D0
        aj(2,2) = -0.4904D0
      
        aj(3,0) = 0.1437D0
        aj(3,1) = 0.9129D0
        aj(3,2) = -0.6940D0
        aj(3,3) = 0.2503D0

        Ntot = 0
        alphan = 0.D0
        DO i = 1, 3
          Ntot = Ntot + Numn(i)
          xn = DBLE(Zatom - Ntot)
          IF( xn == 0.D0 ) xn = 1.D0
          suman = 0.D0
          DO j = 0, 3
            suman = suman + aj(i,j) / (xn**j)
          END DO
          alphan(i) = (xn + 1.D0) * suman
          WRITE(6,*) i, xn, alphan(i)
        END DO

      ELSE IF( KIND_POT == 2 ) THEN      !SIMONS-FUES POT., for Rb

        il = MAX(lmax,3)
        ALLOCATE( Bl(0:il) )
        Bl = 0.D0

        Bl(0) = 0.72657D0
        Bl(1) = 0.47095D0
        Bl(2) = -0.55508D0
        Bl(3) = -0.04008D0

      END IF

!      READING VARIABLES TO CONFIGURE THE FIELD, IF NEEDED...

!      KIND_PI == 0 : Only calculates electronic structure for the indicate l
!      KIND_PI == 1, 2 : Calculate PI for the indicated gauge...
!                        1 = Length gauge
!                        2 = Velocity gauge
!      KIND_PI == 3 : Calculate PI with a Gaussian beam... No dipolar approx.
!      KIND_PI == 4 : Calculate PI with a Laguerre-Gaussian beam...
!      KIND_PI == 5 : Vector Beam with Bessel profile... Electric term only
!      KIND_SCP == 0 : No Scalar-Potential considered
!                  1 : Scalar-Potential from Lorenz-gauge included

      KIND_PI = 0                       !ONLY ELECTRONIC STRUCTURE
      KIND_SCP = 0                      !NO EM SCALAR POTENTIAL
      KIND_TD = 0                       !NO TIME PROPAGATION
      KIND_RK = 6                       !ORDER RUNGE--KUTTA SOLVER
      KIND_VEC = 0                      !0: ALL STATES, 1: ONLY BOUND STATES
      A0 = 0.0D0                        !INTENSITY VECTOR POTENTIAL
      I0 = 0.0D0                        !INTENSITY OF THE FIELD
      A01 = 0.D0
      I01 = 0.D0
      w0 = 0.0D0                        !BEAM WAIST-SIZE
      Eph = 0.0D0                       !PHOTON ENERGY
      mph = 0                           !PHOTON SPIN (POLARIZATION). 0=LP (z), -1=RP, +1=LP
      moam = 0                          !TOPOLOGICAL CHARGE
      np = 0                            !RADIAL NODES
      b0 = 0.D0                         !IMPACT PARAMETER
      afocus = 0.D0
      nEpts = 10                        !NUM. En. PTS TO EXTRAPOLATE
      Eref = 0.D0                       !IF nEpts < 0, Energy reference values (i*Eref)
      nthpts = 1                        !NUM. PTS FOR ANGULAR DISTRIBUTIONS
      nphpts = 1
      ncyc = 0                          !NUMBER OF OPTICAL CYCLES. IF == 0: Dirac-Delta
      bx = 0.D0                         !beta Ang. - ROTATION OF THE FIELD AROUND x AXIS
      B0z = 0.D0                        !INTENSITY STATIC MAGNETIC FIELD
      t_delay = 0.D0                    !TIME DELAY (PUMP-PROBE)
      ncyc2 = 0
      Eph2 = 0.D0
      A0x = 0.D0                        !PLANE WAVE POL. x
      A0y = 0.D0                        !PLANE WAVE POL. y
      A0z = 1.D0                        !PLANE WAVE POL. z
      READ(5,VARS_FIELD)

      WRITE(6,'(/,A17)') 'Field Parameters:'
      IF( A0 /= 0.D0 ) WRITE(6,'(A4,G12.5)') 'A0 =', A0
      IF( I0 /= 0.D0 ) WRITE(6,'(A4,G12.5)') 'I0 =', I0
      IF( w0 /= 0.D0 ) WRITE(6,'(A4,G12.5)') 'w0 =', w0
      IF( b0 /= 0.D0 ) WRITE(6,'(A4,G12.5)') 'b0 =', b0
      IF( Eph /= 0.D0 ) WRITE(6,'(A5,G12.5)') 'Eph =', Eph
      IF( A01 /= 0.D0 ) WRITE(6,'(A5,G12.5)') 'A01 =', A01
      IF( I01 /= 0.D0 ) WRITE(6,'(A5,G12.5)') 'I01 =', I01
      IF( moam /= 0.D0 ) WRITE(6,'(A20,I3)') 'Topological Charge =', moam
      IF( afocus /= 0.D0 ) WRITE(6,'(A22,G12.5)') 'Focusing angle (Deg.):', afocus
      IF( ncyc /= 0.D0 ) WRITE(6,'(A18,I3)') 'Num. Opt. Cycles =', ncyc
      IF( KIND_SCP == 1 ) WRITE(6,'(A25)') 'Lorenz Scalar-Potential Included'
      IF( KIND_TD /= 0 ) WRITE(6,'(A19,I2)') 'Runge--Kutta Order:', KIND_RK
      IF( t_delay /= 0.D0 ) WRITE(6,'(A11,G12.5,A3)') 'Time-Delay:', t_delay, ' fs'
      IF( A0x /= 0.D0 ) WRITE(6,'(A18)') 'Laser Pulse Pol. x'
      IF( A0y /= 0.D0 ) WRITE(6,'(A18)') 'Laser Pulse Pol. y'
      IF( A0z /= 0.D0 ) WRITE(6,'(A18)') 'Laser Pulse Pol. z'

!      IF( I0 == 0.D0 ) I0 = 1.0D15
      IF( A0 == 0.D0 ) A0 = SQRT(I0 / I0_au) / Eph
!      IF( I01 == 0.D0 ) I01 = 1.0D15
      IF( A01 == 0.D0 ) A01 = SQRT(I01 / I0_au) / Eph
      kph = Eph / c_au
      IF( kph == 0.D0 .AND. (KIND_PI == 3 .OR. KIND_PI == 4) ) kph = 1.D0 / c_au
      IF( KIND_SCP == 1 .AND. Eph == 0.D0 ) Eph = 1.D0
      IF( KIND_PI == 5 .OR. KIND_PI == 6 ) THEN
        afocus = afocus * PI / 180.D0
        qvecz = kph * COS(afocus)
        qvecr = kph * SIN(afocus)
      END IF

      IF( nfib >= 0 ) THEN
        nfib0 = FIBONACCI(nfib)
        nfib1 = FIBONACCI(nfib-1)
        dzf = 2.D0 / DBLE(nfib0)
        dphf = 2.D0 * PI * DBLE(nfib1) / DBLE(nfib0)
        WRITE(6,'(A41,I5)') 'Number of Pts for Fibonacci sampling Pts:', nfib0
      END IF

      IF( bx /= 0.D0 ) bx = bx * PI / 180.D0
      IF( B0z /= 0.D0 ) B0z = B0z / B0_au

      Epump = SQRT(I0 / I0_au)
      Eprobe = 0.D0
      IF( KIND_PI >= 8 .AND. KIND_POT == 0 ) THEN
        nt0 = n0_ini
        ntf = 20
        Eph = 0.5D0 * ((1.D0/(DBLE(nt0)**2)) - 1.D0/(DBLE(ntf)**2))
        ncyc = CEILING(DBLE((ntf**2-nt0**2)) / DBLE(nt0**2-ntf**2+(nt0*ntf)**2))
        ncyc = MAX(ncyc,10)
        WRITE(6,'(A31,I5)') 'Modified Num. Opt. Cycles Pump:', ncyc
        WRITE(6,'(A28,G14.7)') 'Modified Photon Energy Pump:', Eph
        IF( I01 == 0.D0 ) I01 = I0
        Eprobe = SQRT(I01 / I0_au)
        IF( Eph2 == -1.D0 ) THEN
          Eph2 = Eph
          ncyc2 = ncyc
        ELSE
          ntf2 = ntf + 10
          Eph2 = 0.5D0 * ((1.D0/(DBLE(ntf)**2)) - 1.D0/(DBLE(ntf2)**2))
          ncyc2 = CEILING(DBLE((ntf2**2-ntf**2)) / DBLE(ntf**2-ntf2**2+(ntf*ntf2)**2))
        END IF
        ncyc2 = MAX(ncyc2,2)
        WRITE(6,'(A32,I5)') 'Modified Num. Opt. Cycles Probe:', ncyc2
        WRITE(6,'(A29,G14.7)') 'Modified Photon Energy Probe:', Eph2
        t_delay = t_delay / t_au
        afocus = afocus * PI / 180.D0
        kph = Eph2 / c_au
        qvecz = kph * COS(afocus)
        qvecr = kph * SIN(afocus)
      ELSE IF( KIND_POT /= 0 ) THEN
        Eprobe = SQRT(I01 / I0_au)
        t_delay = t_delay / t_au
        afocus = afocus * PI / 180.D0
        kph = Eph2 / c_au
        qvecz = kph * COS(afocus)
        qvecr = kph * SIN(afocus)
        WRITE(6,'(A31,I5)') 'Modified Num. Opt. Cycles Pump:', ncyc
        WRITE(6,'(A28,G14.7)') 'Modified Photon Energy Pump:', Eph
        WRITE(6,'(A32,I5)') 'Modified Num. Opt. Cycles Probe:', ncyc2
        WRITE(6,'(A29,G14.7)') 'Modified Photon Energy Probe:', Eph2
        WRITE(6,'(A28,G14.7)') 'Modified Photon Wave Number:', kph
      END IF

      WRITE(6,'(A7,G14.7)') 'Epump =', Epump
      WRITE(6,'(A8,G14.7)') 'Eprobe =', Eprobe

      END SUBROUTINE READ_INPUTS
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE READ_COUP

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES, ONLY: nfun, nointv, ka
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER :: i, ibra, jket, li, ni, nj, n0, l0, m0, nmin, nmax
      INTEGER :: IERROR
      REAL(DP) :: En

      REAL(DP), ALLOCATABLE :: f(:)

      WRITE(6,'(/,A)') 'Reading Energies'

      OPEN( UNIT=10, FILE='Enl.dat', ACTION='READ' )
      READ(10,*) nfun
      
      ALLOCATE( Enl(nfun,0:lmax), n01(0:lmax,3) )
      Enl = 0.D0
      n01 = 0

      DO li = 0, lmax
        DO ni = 1, nfun

          READ(10,*,IOSTAT=IERROR) i, En
          IF(IERROR /= 0 ) THEN
            WRITE(6,*) 'Error Reading Energies'
            STOP
          END IF
          Enl(ni,li) = En
          IF( En < 0.D0 ) n0_fin = ni
          IF( En <= Emax_fin ) n1_fin = ni
        END DO

        n01(li,1) = 1
        n01(li,2) = n1_fin
        n01(li,3) = n0_fin

      END DO

      CLOSE(10)

      WRITE(6,'(/,A17)') 'Reading Couplings'

      OPEN( UNIT=10, FILE='CSs/MatElem_All.dat', ACTION='READ' )
      READ(10,*) n1_max, nbra, nket

      n0 = n0_ini
      l0 = lmf(0,1)
      m0 = lmf(0,2)
      
      nmin = 1
      nmax = n1_max

!      nbra = n1_max * nlm
!      nket = n1_max * nlm

      nfields = 1
      IF( KIND_PI == 5 .OR. KIND_PI == 6 ) THEN
        nfields = 2
      ELSE IF( KIND_PI >= 8 ) THEN
        nfields = 5
      END IF

      WRITE(6,*) nbra, nket, nfields

      ALLOCATE( zHint_ij(nbra,nket,nfields), f(2*nfields) )
      zHint_ij = 0.D0
      f = 0.D0

      DO ni = 1, nbra
        DO nj = ni, nket
          READ(10,*,IOSTAT=IERROR) ibra, jket, (f(i), i=1,2*nfields)
          IF( IERROR /= 0 ) THEN
            WRITE(6,*) 'Error Reading Couplings File:', IERROR
            STOP
          END IF
          zHint_ij(ibra,jket,1) = DCMPLX(f(1),f(2))         !E0: Dip.
          IF( KIND_PI >= 8 ) THEN
            zHint_ij(ibra,jket,2) = DCMPLX(f(3),f(4))       !E_rho
            zHint_ij(ibra,jket,3) = DCMPLX(f(5),f(6))       !E_z
            zHint_ij(ibra,jket,4) = DCMPLX(f(7),f(8))       !B_phi
            zHint_ij(ibra,jket,5) = DCMPLX(f(9),f(10))      !B_0
          END IF
        END DO
      END DO

      CLOSE(10)

      END SUBROUTINE READ_COUP
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE READ_TDCOEFF

      USE MOD_TYPES
      USE MOD_GRID, ONLY: lmax, llmax, nfib, nfib0
      USE MOD_BSPLINES, ONLY: nfun, nointv, ka
      USE MOD_MATRICES, ONLY: rvecij, zrangij, zPmq, nmsize, nmsel, nvsize, nvsel
      USE MOD_PHOTOION
      USE MOD_FAC_VARS

      IMPLICIT NONE

      INTEGER :: i, imin, imax, il, ni, li, mi, n0, l0, m0, nmin, nmax, ibra, jket
      INTEGER :: IERROR, j, jl, nj, lj, mj, lf, mq, nr, ir, ith, nbold, nt_old
      REAL(DP) :: En, sumctf, ti, fr, fi(2)
      CHARACTER(LEN=23) :: cfile1

      REAL(DP), ALLOCATABLE :: fvec(:)

!...............................................................................

      WRITE(6,'(/,A16)') 'Reading Energies'

      OPEN( UNIT=10, FILE='Enl.dat', ACTION='READ' )
      READ(10,*) nfun
      
      ALLOCATE( Enl(nfun,0:lmax), n01(0:lmax,3) )
      Enl = 0.D0
      n01 = 0
      nb_max = 0
      DO li = 0, lmax
        DO ni = 1, nfun
          READ(10,*,IOSTAT=IERROR) i, En
          IF(IERROR /= 0 ) THEN
            WRITE(6,*) 'Error Reading Energies'
            STOP
          END IF
          Enl(ni,li) = En
          IF( En < 0.D0 ) THEN
            n0_fin = ni
            nbold = ni
          END IF
          IF( En <= Emax_fin ) n1_fin = ni
        END DO

        nb_max = MAX(nb_max,nbold)
        n01(li,1) = 1
        n01(li,2) = n1_fin
        n01(li,3) = n0_fin

      END DO

      CLOSE(10)

      CALL SEL_STATES(1,nb_max)

!...............................................................................

      WRITE(6,'(/,A20)') 'Reading Coefficients'

      OPEN( UNIT=10, FILE='CSs/MatElem_All.dat', ACTION='READ' )
      READ(10,*) n1_max, nbra, nket

      nvec = nbra

      CLOSE(10)

      n0 = n0_ini
      l0 = lmf(0,1)
      m0 = lmf(0,2)
      
      nm = nlm

      WRITE(6,'(A5,I6)') 'nvec:', nvec

      nmin = 1
      nmax = n1_max

      ALLOCATE( zf(nvec), Eall(nvec) )
      zf = 0.D0

      OPEN( UNIT=15, FILE='CSs/TDSE_COEFFs.dat', ACTION='READ' )

      imin = 0
      imax = 0
      nbounds = 0
      sumctf = 0.D0
      DO i = 1, nvec
        READ(15,*,IOSTAT=IERROR) ni, fi(1), fi(2)
        IF( IERROR /= 0 ) WRITE(6,*) 'An Error Ocurred Reading the Probs:', li, mi
        zf(i) = CMPLX(fi(1),fi(2))
        sumctf = sumctf + ABS(zf(i))**2
        ni = nl_bra(i,1)
        li = nl_bra(i,2)
        Eall(i) = Enl(ni,li)
      END DO

!      DO il = 1, nlm
!        li = lmf(il,1)
!        mi = lmf(il,2)
!        imin = imax + 1
!        imax = imax + n1_max
!        DO i = imin, imax
!          READ(15,*,IOSTAT=IERROR) ni, fi(1), fi(2) !zf(i)
!          IF( IERROR /= 0 ) WRITE(6,*) 'An Error Ocurred Reading the Probs:', li, mi
!          zf(i) = CMPLX(fi(1),fi(2))
!          sumctf = sumctf + ABS(zf(i))**2
!          Eall(i) = Enl(i-imin+1,li)
!          IF( Eall(i) < 0.D0 ) nbounds = nbounds + 1
!        END DO
!      END DO

      CLOSE(15)

      WRITE(6,'(A11,G20.10)') 'Sum c(tf) =', sumctf
      WRITE(6,'(/,A18,I4)') 'Num. Bound States:', nbounds
      WRITE(6,'(A12,I3)') 'N Bound Max:', nb_max

      nmsize = nvec
      nvsize = n1_max
      IF( KIND_VEC == 0 ) THEN      !ALL STATES
        nmsel = nvec
        nvsel = n1_max
      ELSE                          !ONLY BOUND STATES
        nmsel = nbounds
        nvsel = nb_max
      END IF

!...............................................................................

      WRITE(6,'(/,A35)') 'Reading Time-Dependent Coefficients'

      OPEN( UNIT=25, FILE='CSs/TD_Coeffs_All.dat', ACTION='READ' )
        
      ALLOCATE( fvec(2*nvec) )

      nt = 0
drct: DO
        READ(25,*,IOSTAT=IERROR) ti, (fvec(i), i=1,2*nvec)
        IF( IERROR /= 0 ) EXIT drct
        nt = nt + 1
      END DO drct

      ntsteps = 1
      IF( nt <= 1000 ) THEN
        ntsteps = 1
      ELSE IF( nt > 1000 .AND. nt <= 2000 ) THEN
        ntsteps = 2
      ELSE IF( nt > 2000 .AND. nt <= 10000 ) THEN
!        ntsteps = 10
        ntsteps = 5
!      ELSE IF( nt > 10000 .AND. nt <= 50000 ) THEN
!        ntsteps = 50
!      ELSE IF( nt > 50000 .AND. nt <= 100000 ) THEN
!        ntsteps = 100
      ELSE
!        ntsteps = 500
        ntsteps = 10
      END IF

      WRITE(6,'(A4,I8)') 'nt =', nt
      WRITE(6,'(A28,I4)') 'Num. Time Steps to Consider:', ntsteps

      nt_old = nt
      IF( ntsteps > 1 ) nt = CEILING(DBLE(nt)/DBLE(ntsteps)) + 1

      WRITE(6,'(A8,I6)') 'New nt =', nt

      ALLOCATE( zcit(nt,nvec), tall(nt) )
      zcit = 0.D0
      tall = 0.D0

      REWIND(UNIT=25)

      ni = 0
      DO i = 1, nt_old
        READ(25,*,IOSTAT=IERROR) ti, (fvec(j), j=1,2*nvec)
        IF( MOD(i-1,ntsteps) == 0 .OR. i == nt_old ) THEN
          ni = ni + 1
          tall(ni) = ti
          j = 0
          DO il = 1, 2*nvec, 2
            j = j + 1
            zcit(ni,j) = CMPLX(fvec(il),fvec(il+1))
          END DO
        END IF
      END DO

      CLOSE(25)

!      PROPAGATION TIMES

!      PUMP DURATION

      Tpulse = 2.D0 * PI / Eph

      tiau = 0.D0
      tfau = tiau + DBLE(ncyc) * Tpulse

      Tpump = tfau - tiau
      td1 = Tpump / 2.D0

!      PROBE DURATION

      Tpulse2 = 0.D0
      td2 = 0.D0
      t2iau = 0.D0
      t2fau = 0.D0

      IF( Eprobe /= 0.D0 .AND. Eph2 /= 0.D0 ) THEN
        Tpulse2 = 2.D0 * PI / Eph2
        Tprobe = ncyc2 * Tpulse2
        td2 = Tprobe / 2.D0
        t2iau = t_delay - 0.5D0*(Tprobe - Tpump)
        t2fau = t2iau + Tprobe
      END IF

!      TOTAL PROPAGATION

      IF( t2iau < tiau ) THEN
        tiau = tiau - t2iau
        tfau = tfau - t2iau
        t2fau = t2fau - t2iau
        t2iau = 0.D0
      END IF

      WRITE(6,'(A18,G12.5,A2)') 'Pump Pulse Period:', Tpulse * t_au, 'fs'
      WRITE(6,'(A20,G12.5,A2)') 'Pump Pulse Duration:', Tpump*t_au, 'fs'
      WRITE(6,'(A19,G12.5,A5,G12.5,A2)') 'Pump Pulse in: ti =', tiau * t_au, &
     &                                   ' tf =', tfau * t_au, 'fs'
      WRITE(6,'(T15,A4,G12.5,A5,G12.5,A4)') 'ti =', tiau, ' tf =', tfau, 'a.u.'
      IF( Eprobe /= 0.D0 ) THEN
        WRITE(6,*)
        WRITE(6,'(A19,G12.5,A3)') 'Probe Pulse Period:', Tpulse * t_au, ' fs'
        WRITE(6,'(A21,G12.5,A3)') 'Probe Pulse Duration:', Tpump*t_au, ' fs'
        WRITE(6,'(A20,G12.5,A5,G12.5,A3)') 'Probe Pulse in: ti =', t2iau * t_au, &
     &                                   ' tf =', t2fau * t_au, 'fs'
      WRITE(6,'(T15,A4,G12.5,A5,G12.5,A4)') 'ti =', t2iau, ' tf =', t2fau, 'a.u.'
      END IF

!...............................................................................

      WRITE(6,'(/,A29)') 'Calculating Needed Factorials'

      llmax = 10                        !TO USE LATER TO CALCULATE Avec...
      li = MAX(lmax+1,llmax+1)
      ALLOCATE( fac(0:2*li) )
      fac = 1.D0
      DO i = 2, 2*li
        fac(i) = i * fac(i-1)
      END DO

100   FORMAT('CSs/Prob_TDSE_',I2.2,'+',I2.2,'.dat')
105   FORMAT('CSs/Prob_TDSE_',I2.2,'-',I2.2,'.dat')

      END SUBROUTINE READ_TDCOEFF
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE READ_FR(INDX)
      
!      SUBROUTINE TO READ NEEDED RADIAL FUNCTIONS, INCLUDING FOR TOROIDAL MOMENT

      USE MOD_TYPES
      USE MOD_BSPLINES, ONLY: nkp, k, ka, rt, nfun, gauleg
      USE MOD_GRID, ONLY: lmax, nfib, nfib0, rtot, wtot
      USE MOD_MATRICES, ONLY: fur, dfur, rvecij, zrangij, zPmq, nmsize, nmsel, &
     &                         nvsize, nvsel, zJijq
      USE MOD_PHOTOION, ONLY: n1_max, cinl, KIND_PI, nlm, lmf, nvec, nl_bra,   &
     &                         nbounds, nl_ket, Eball, Eall, n0_ini, ist0,      &
     &                         num_i0, KIND_VEC, SEL_T0_STATES

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: INDX

      INTEGER :: i, j, igl, left, nr, lf, il, ni, li, mi, jl, nj, lj, mj, mq
      INTEGER :: n0, l0, m0, ith, ir, IERROR
      REAL(DP) :: f1, f2, fr, fi(8)

      REAL(DP), ALLOCATABLE :: xgl(:), wgl(:)

!      WRITE(6,*) 'Reading Eigenvectors...'

!...............................................................................

!      READING EIGENVECTORS AND CALCULATING REDUCED RADIAL WAVE FUNCTIONS

      nr = (nkp-1) * ka !(nkp - 2*k + 4) * ka

!     ALLOCATE( fur(nr,n1_max,0:lmax), dfur(nr,n1_max,0:lmax) )
      ALLOCATE( rtot(nr), wtot(nr), xgl(ka), wgl(ka) )

      fur = 0.D0
      dfur = 0.D0
      
      rtot = 0.D0
      wtot = 0.D0

      xgl = 0.D0
      wgl = 0.D0

      CALL gauleg(-1.D0,1.D0,xgl,wgl,ka)
      i = 0
      DO left = k-1, nkp-1-(k-2)
        f1 = ( rt(left+1) + rt(left) ) / 2.D0
        f2 = ( rt(left+1) - rt(left) ) / 2.D0
        DO igl = 1, ka
          i = i + 1
          rtot(i) = f1 + xgl(igl) * f2
          wtot(i) = f2 * wgl(igl)
        END DO
      END DO

!      i = 0
!      DO left = 1, nointv
!        CALL gauleg(rtk(left),rtk(left+1),xgl,wgl,ka,ka)
!        DO igl = 1, ka
!          i = i + 1
!          rtot(i) = xgl(igl)
!          wtot(i) = wgl(igl)
!        END DO
!      END DO

!...............................................................................

!      TOROIDAL MOMENT RELATED RADIAL FUNCTIONS

iftor:IF( KIND_PI >= 8 ) THEN

        WRITE(6,'(/,A40)') 'Reading Related Toroidal Matrix Elements'

!        nr = nointv * ka
        lf = lmax

        ALLOCATE( zrangij(0:lf,-lf:lf,0:lf,-lf:lf,3) )
        ALLOCATE( zPmq(nlm,nlm,-1:1,0:nfib0,4) )

        ALLOCATE( rvecij(n1_max,0:lf,n1_max,0:lf) )

        zrangij = 0.D0
        rvecij = 0.D0

        zPmq = 0.D0

!        TO CALCULATE MEAN VALUES OF r VECTOR

        OPEN( UNIT=20, FILE='CSs/rMatElemAng.dat', ACTION='READ' )
        OPEN( UNIT=25, FILE='CSs/rMatElemRad.dat', ACTION='READ' )

!        TO CALCULATE CURRENT DENSITY AND ITS MEAN VALUE

        OPEN( UNIT=30, FILE='CSs/CurrentMatFunAng.dat', ACTION='READ' )
!        OPEN( UNIT=35, FILE='CSs/CurrentMatFunRad.dat', ACTION='READ' )

        WRITE(6,'(A13)') 'Angular Terms'

drdrang:DO
          READ(20,*,IOSTAT=IERROR) li, mi, lj, mj, (fi(i), i=1,6)
          IF( IERROR /= 0 ) EXIT drdrang
          zrangij(li,mi,lj,mj,1) = CMPLX(fi(1),fi(2))
          zrangij(li,mi,lj,mj,2) = CMPLX(fi(3),fi(4))
          zrangij(li,mi,lj,mj,3) = CMPLX(fi(5),fi(6))
        END DO drdrang

        CLOSE(20)
        
drdfang:DO
          READ(30,*,IOSTAT=IERROR) il, jl, mq, ith, (fi(i), i=1,8)
          IF( IERROR /= 0 ) EXIT drdfang
          zPmq(il,jl,mq,ith,1) = CMPLX(fi(1),fi(2))
          zPmq(il,jl,mq,ith,2) = CMPLX(fi(3),fi(4))
          zPmq(il,jl,mq,ith,3) = CMPLX(fi(5),fi(6))
          zPmq(il,jl,mq,ith,4) = CMPLX(fi(7),fi(8))
        END DO drdfang

        CLOSE(30)

        WRITE(6,'(A12)') 'Radial Terms'

  drdr: DO
          READ(25,*,IOSTAT=IERROR) ni, li, nj, lj, fr
          IF( IERROR /= 0 ) EXIT drdr
          rvecij(ni,li,nj,lj) = fr
        END DO drdr

        CLOSE(25)

        CALL SEL_STATES(1,n1_max)

!        ALLOCATE( nl_bra(nmsize,4), nl_ket(nmsel,4), Eball(nmsel) )
        ALLOCATE( Eball(nmsel) )

        CALL SEL_T0_STATES

      END IF iftor

      IF( INDX == 1 ) THEN
      
        WRITE(6,'(/,A25)') 'Reading J Matrix Elements'

        OPEN( UNIT=40, FILE='CSs/JMat_Int.dat', ACTION='READ' )
        ALLOCATE(zJijq(nvec,nvec,-1:4,2))
        zJijq = 0.D0

 drJij: DO
          READ(40,*,IOSTAT=IERROR) il, jl, mq, (fi(i), i=1,4)
          IF( IERROR /= 0 ) EXIT drJij
          zJijq(il,jl,mq,1) = CMPLX(fi(1),fi(2))
          zJijq(il,jl,mq,2) = CMPLX(fi(3),fi(4))
        END DO drJij

        CLOSE(40)

      END IF

      END SUBROUTINE READ_FR
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE READ_EIGENVEC

      USE MOD_TYPES
      USE MOD_BSPLINES, ONLY: nfun
      USE MOD_GRID, ONLY: lmax
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER :: i, l, ni, n1, n2, n3, IERROR

      WRITE(6,'(/,A23)') 'Reading Eigenvectors...'

      OPEN( UNIT=80, FILE='Eigenvec_All.dat', ACTION='READ' )

      READ(80,*) n1, n2, n3

      IF( n1 /= nfun ) WRITE(6,*) 'Problem with number of Functions'
      IF( n2 /= n1_max ) WRITE(6,*) 'Problem with number of states'
      IF( n3 /= lmax ) WRITE(6,*) 'Problem with number of partial waves'

      ALLOCATE( cinl(nfun,n1_max,0:lmax) )
      cinl = 0.D0

      DO l = 0, lmax
        READ(80,*,IOSTAT=IERROR) n1
        IF( n1 /= l ) WRITE(6,*) 'Problem Reading L:', l, n1
        DO ni = 1, n1_max
          READ(80,*,IOSTAT=IERROR) n2, (cinl(i,ni,l), i=1,nfun)
          IF( IERROR /= 0 ) THEN
            WRITE(6,*) 'Problem Reading Eigenvectors...'
            STOP
          END IF
        END DO
      END DO
      
      CLOSE(80)

      END SUBROUTINE READ_EIGENVEC
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE SEL_STATES(nmin,nmax)

!     ROUTINE TO SELECT EIGENSTATES. CHOOSES NEW CONTINUUM STATES IF A SPECIAL
!       ENERGY GRID IS NEEDED

      USE MOD_TYPES
      USE MOD_BSPLINES
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nmin, nmax

      INTEGER :: i, j, il, jl, n0, l0, m0, ni, li, mi, nf, lf, mf, nnew, nold
      REAL(DP) :: Ethreshold, Ej, dE
      REAL(DP), ALLOCATABLE :: Egrid(:)

      n0 = n0_ini
      l0 = lmf(0,1)
      m0 = lmf(0,2)

      IF( ALLOCATED(nl_bra) ) DEALLOCATE(nl_bra)
      IF( ALLOCATED(nl_ket) ) DEALLOCATE(nl_ket)

      IF( KIND_EGR == 0 ) THEN                  !EIGENSPECTRUM

        nbra = n1_max * nlm
        nket = (nmax - nmin + 1) * nm

        ALLOCATE( nl_bra(nbra,4), nl_ket(nket,4) )

!       VECTOR TO STORE ALL INDEXES OF FINAL STATES (INTERACTION MATRIX - WISE)

        i = 1
        DO il = 1, nlm
          lf = lmf(il,1)
          mf = lmf(il,2)
          n0_fin = n01(lf,1)
          n1_fin = n1_max
          DO nf = n0_fin, n1_fin
            nl_bra(i,1) = nf
            nl_bra(i,2) = lf
            nl_bra(i,3) = mf
            nl_bra(i,4) = il
            i = i + 1
          END DO
        END DO

!       VECTOR TO STORE ALL INDEXES OF INITIAL STATES (INTERACTION MATRIX - WISE)

        i = 1
        DO jl = 1, nm
          li = l0
          mi = m0
          IF( KIND_TD == 1 ) THEN
            li = lmf(jl,1)
            mi = lmf(jl,2)
          END IF
          DO ni = nmin, nmax
            nl_ket(i,1) = ni
            nl_ket(i,2) = li
            nl_ket(i,3) = mi
            nl_ket(i,4) = jl
            i = i + 1
          END DO
        END DO

      ELSE                                      !CLOSE TO THE SELECTED ENERGY GRID

        Ethreshold = 0.D0

        ALLOCATE( Egrid(nEpts) )
        Egrid = 0.D0
        dE = (Emax_fin - Ethreshold) / (DBLE(nEpts)**2)
        DO i = 1, nEpts
          Egrid(i) = Ethreshold + dE * (i**2)
        END DO

        li = lmf(1,1)
        lf = lmf(nlm,1)

        nbounds = 0.D0
        DO il = li, lf
          DO i = 1, nfun
            IF( Enl(i,il) <= 0.D0 ) nbounds = nbounds + 1
          END DO
        END DO

        nbra = nbounds + nlm * nEpts
!        IF( m0 /= 0 ) nbra = 2 * nbounds + nlm * nEpts
        nket = nbra

        ALLOCATE( nl_bra(nbra,4), nl_ket(nket,4) )
        
        nl_bra = 0.D0
        nl_ket = 0.D0

        WRITE(6,'(/,A19,2I6)') 'Selecting States...', nbra, nket

        i = 1
        DO il = 1, nlm

          lf = lmf(il,1)
          mf = lmf(il,2)
          j = 1
          nf = 0
          nold = 0

dosel_nlm:DO

            nf = nf + 1
            Eref = Enl(nf,lf)
            Ej = Egrid(j)

            IF( i >= nbra ) WRITE(6,'(A26,I5,A,I5)') 'Achtung! Problems with i:', &
     &                      i, '/', nbra

            IF( i > nbra .OR. j > nEpts .OR. nf >= nfun ) EXIT dosel_nlm

            IF( Eref <= Ethreshold ) THEN
              WRITE(6,'(I5,3I3,X,G20.10)') i, nf, lf, mf, Eref
              nl_bra(i,1) = nf
              nl_bra(i,2) = lf
              nl_bra(i,3) = mf
              nl_bra(i,4) = il
              i = i + 1
            ELSE IF( Eref >= Ej ) THEN
              nnew = nf
              IF( nnew == nold ) nnew = nnew + 1
              WRITE(6,'(I5,"-",I6,3I4,X,G20.10)') i, j, nnew, lf, mf, Eref
              nl_bra(i,1) = nnew
              nl_bra(i,2) = lf
              nl_bra(i,3) = mf
              nl_bra(i,4) = il
              i = i + 1
              j = j + 1
              nold = nnew
            END IF

            IF( i > nbra .OR. j > nEpts .OR. nf >= nfun ) EXIT dosel_nlm
!            IF( j > nEpts ) EXIT dosel_nlm
!            IF( nf >= nfun ) EXIT dosel_nlm
              
          END DO dosel_nlm

        END DO

        WRITE(6,'(A14)') 'Selecting Done'
        
        nl_ket = nl_bra
        
        DEALLOCATE( Egrid )
        
      END IF

      END SUBROUTINE SEL_STATES
