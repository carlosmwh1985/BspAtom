      SUBROUTINE MATRIX_SVT

      USE MOD_TYPES
      USE MOD_BSPLINES
      USE MOD_GRID
      USE MOD_MATRICES
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER :: ibra, jket, ibet, ibetmin, ibetmax, igl, ifun, jfun, lf, mf
      INTEGER :: i, j, il, jl, left, iglmin, iglmax, nle
      REAL(DP) :: f1, f2, r, dr, Vpot, Vcent, sumS, sumV, sumT, sumr, Vl
      REAL(DP) :: A1, A2, b1, b2, fbra, fket, dfbra, dfket, Alm, sumc, sumd
      COMPLEX(DPC) :: zAlm, zdAlm, zfAr, zErlm, zEzlm
      
      REAL(DP), ALLOCATABLE :: bsp(:), dbsp(:), bsp1(:), bspp(:), sumU(:)
      COMPLEX(DPC), ALLOCATABLE :: zsumc(:,:), zsumd(:,:), zsume(:,:), zsumf(:,:)

      ALLOCATE( Sij(nfun,nfun), Vij(nfun,nfun), Uij(nfun,nfun,0:lmax), Tij(nfun,nfun) )
      ALLOCATE( sumU(0:lmax) )

!      KIND_PI == 0: ELECTRONIC STRUCTURE ONLY
!              == 1: LENGTH GAUGE
!              == 2: VELOCITY GAUGE
!              == 3: GAUSSIAN BEAM (VEL. - VECTOR POT.)
!              == 4: LAGUERRE--GAUSSIAN BEAM (VEL. - VECTOR POT.)
!              == 5: RVB VECTOR BEAM (LEN. - POINCARÉ GAUGE - ELECTRIC FIELD)
!              == 6: AVB VECTOR BEAM (LEN. - POINCARÉ GAUGE - ELECTRIC FIELD)
!              == 7: VECTOR POTENTIAL - AHARONOV--BOHM EFFECT

      IF( KIND_PI <= 2 ) THEN                        !VELOCITY GAUGE
        nm = 2
      ELSE IF( KIND_PI == 3 ) THEN                  !VEL. - GAUSS. BEAM
        nle = 1
      ELSE IF( KIND_PI == 4 ) THEN                  !VEL. LAGUERRE-GAUSS. BEAM
!!!        nm = 2*lmax + 1
        nle = nlm
      ELSE IF( KIND_PI >= 5 ) THEN                  !LEN. - POINCARÉ G. VECTOR BEAM
        nle = 1
        IF( KIND_PI >= 8 ) nle = nlm            !LIN. POL. + RVB
      END IF

      IF( KIND_PI <= 2) THEN
        ALLOCATE( rij(nfun,nfun,nm) )
      ELSE
        CALL ZINT_TH
        ALLOCATE( zAij(nfun,nfun,nlm,nm,ncomp), Xij(nfun,nfun) )
        ALLOCATE( zsumc(nlm,nm), zsumd(nlm,nm), zsume(nle,nm), zsumf(nle,nm) )
      END IF

      WRITE(6,'(/,A34)') 'Calculating S, V, U and T Matrices'

      Sij = 0.D0
      Vij = 0.D0
      Uij = 0.D0
      Tij = 0.D0
      IF( KIND_PI <= 2 ) rij = 0.D0
      IF( KIND_PI > 2 ) zAij = 0.D0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nfun,k,ka,nkp,rt,xg,wg,Aind,nlm,nm,zIth, &
!$OMP&  lmax,Bl,rij,Sij,Tij,Uij,Vij,Xij,zAij,KIND_POT,KIND_PI)

      ALLOCATE( bsp(k), dbsp(k), bsp1(k-1), bspp(k+1) )

!$OMP DO

dobra:DO ibra = 1, nfun
 doket: DO jket = 1, nfun

          ibetmin = MAX(ibra,jket)
          ibetmax = MIN(ibra+k-1,jket+k-1)

          sumc = 0.D0
          sumd = 0.D0
          sumr = 0.D0
          sumS = 0.D0
          sumT = 0.D0
          sumU = 0.D0
          sumV = 0.D0

          IF( KIND_PI >= 3 ) THEN
            zsumc = 0.D0
            zsumd = 0.D0
            zsume = 0.D0
            zsumf = 0.D0
          END IF

doBetween:DO ibet = ibetmin, ibetmax

            f1 = ( rt(ibet+1) + rt(ibet) ) / 2.D0
            f2 = ( rt(ibet+1) - rt(ibet) ) / 2.D0

      doGL: DO igl = 1, ka

              r = f1 + xg(igl) * f2
              dr = f2 * wg(igl)

!              CALCULATE B-SPLINE FUNCTIONS AND THEIR DERIVATIVE
              CALL BSPALL(r,left,bsp,dbsp,bsp1,bspp)

              IF( r == 0.D0 ) r = eps
              Vpot = SELPOT(r)

              ifun = ibra - (left-k)
              jfun = jket - (left-k)

              fbra = bsp(ifun)
              fket = bsp(jfun)

              dfbra = dbsp(ifun)
              dfket = dbsp(jfun)

              IF( KIND_PI >= 3 ) THEN
                DO il = 1, nlm
                  DO jl = 1, nm
                    IF( KIND_PI == 3 .OR. KIND_PI == 4 ) THEN
                      zAlm = zIth(ibet,igl,il,jl,1)
                      zfAr = zAlm / r
                      zsumc(il,jl) = zsumc(il,jl) + fbra * zfAr * fket * dr
                      zsumd(il,jl) = zsumd(il,jl) + fbra * zAlm * dfket * dr
                    ELSE
!                      zfAr = zIth(ibet,igl,il,1)                              !Vel. Form
!                      zdAlm = (zIth(ibet,igl,il,2) - zIth(ibet,igl,il,1)) / r
!                      zsumc(il) = zsumc(il) + fbra * zfAr * dfket * dr
!                      zsumd(il) = zsumd(il) + fbra * zdAlm * fket * dr
                      zErlm = zIth(ibet,igl,il,jl,1)                              !Len. Form
                      zEzlm = zIth(ibet,igl,il,jl,2)
                      zsumc(il,jl) = zsumc(il,jl) + fbra * zErlm * fket * dr
                      zsumd(il,jl) = zsumd(il,jl) + fbra * zEzlm * fket * dr
                      IF( KIND_PI >= 8 ) THEN
                        zsume(il,jl) = zsume(il,jl) + &
     &                              fbra * zIth(ibet,igl,il,jl,3) * fket * dr
                        zsumf(il,jl) = zsumf(il,jl) + &
     &                              fbra * zIth(ibet,igl,il,jl,4) * fket * dr
                      END IF
                    END IF
                  END DO
                END DO
              ELSE
                sumc = sumc + fbra * (1.D0 / r) * fket * dr
                sumd = sumd + fbra * dfket * dr
              END IF
              sumr = sumr + fbra * r * fket * dr
              sumS = sumS + fbra * fket * dr
              sumV = sumV + fbra * Vpot * fket * dr
              sumT = sumT + dfbra * 0.5D0 * dfket * dr
              DO lf = 0, lmax
                Vcent = DBLE(lf*(lf+1)) / (2.D0*(r**2))
                Vl = 0.D0
                IF( KIND_POT == 2 ) Vl = Bl(lf) / (r**2)
                sumU(lf) = sumU(lf) + fbra * (Vcent + Vl) * fket * dr
              END DO

            END DO doGL

          END DO doBetween

          IF( KIND_PI == 1 ) THEN                  !LENGTH GAUGE
            rij(ibra,jket,1) = sumr
          ELSE IF( KIND_PI == 2 ) THEN            !VELOCITY GAUGE
            rij(ibra,jket,1) = sumc
            rij(ibra,jket,2) = sumd
          ELSE IF( KIND_PI >= 3 ) THEN            !VEL. - GAUSS. OR LG BEAM
                DO il = 1, nlm
                  DO jl = 1, nm
                    zAij(ibra,jket,il,jl,1) = zsumc(il,jl)
                zAij(ibra,jket,il,jl,2) = zsumd(il,jl)
                IF( KIND_PI == 4 .OR. KIND_PI >= 8 ) THEN
                  zAij(ibra,jket,il,jl,3) = zsume(il,jl)
                  zAij(ibra,jket,il,jl,4) = zsumf(il,jl)
                END IF
              END DO
            END DO
            IF( KIND_PI >= 8 ) Xij(ibra,jket) = sumr
          END IF

!      WRITE(6,*) ibra, jket, sumU(0)

          Sij(ibra,jket) = sumS
          Vij(ibra,jket) = sumV
          Uij(ibra,jket,0:lmax) = sumU(0:lmax)
          Tij(ibra,jket) = sumT

        END DO doket
      END DO dobra

!$OMP END DO

      DEALLOCATE( bsp, dbsp, bsp1, bspp )

!$OMP END PARALLEL

      WRITE(6,'(A19,/)') 'Matrices Calculated'

      DEALLOCATE( sumU )
      
      IF( KIND_PI >= 3 ) DEALLOCATE( zsumc, zsumd, zsume, zsumf )

      END SUBROUTINE MATRIX_SVT
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE SOLVE_SYSTEM

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_MATRICES
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER :: i, j, l, ni, nlim, nE0, ntemp, i_max, nbds, nbold
      REAL(DP) :: Elim
      REAL(DP), ALLOCATABLE :: ctemp(:,:,:)

!      LAPACK VARIABLES
      INTEGER :: LWORK, INFO
      REAL(DP), ALLOCATABLE :: WORK(:)

      LWORK = 4 * nfun

      ALLOCATE( Hij(nfun,nfun), Bij(nfun,nfun), En(nfun) )
      ALLOCATE( WORK(LWORK) )

      IF( KIND_PI == 1 .OR. KIND_PI == 2 ) THEN
        ALLOCATE( E_ini(nfun), E_fin(nfun), ci_ini(nfun) )
      ELSE IF( KIND_PI >= 3 ) THEN
        ALLOCATE( Enl(nfun,0:lmax), rEki(nfun,0:lmax), n01(0:lmax,3) )
        rEki = 1.D0
      END IF

      n0_fin = -1
      n1_fin = -1
      nlim = 0
      nbds = 0

      OPEN( UNIT=75, FILE='Enl.dat', ACTION='WRITE' )
      WRITE(75,*) nfun

dol0: DO l = 0, lmax

        Hij(:,:) = Tij(:,:) + Uij(:,:,l) + Vij(:,:)      !HAMILTONIAN
        Bij = Sij                                          !OVERLAP

!        SOLVE THE GENERALIZED EIGENVALUE PROBLEM. THE OUTPUT IS THE
        CALL DSYGV(1,'V','U',nfun,Hij,nfun,Bij,nfun,En,WORK,LWORK,INFO)

        IF( INFO /= 0 ) THEN
          WRITE(6,*) 'ERROR DIAGONALIZING THE MATRIX!', INFO
          WRITE(6,100) 'l = ', l
          STOP
        END IF

        WRITE(6,100) 'l0 = ', l
        WRITE(6,110) 'HC = ESC eigenvalue solved'
        WRITE(6,120) 'n', 'Eigenvalues'
        WRITE(6,120) '-', '-----------'

        i_max = MIN(nfun,200)
        DO i = 1, nfun
          IF( i <= 20 ) WRITE(6,200) i+l, En(i)
          WRITE(75,200) i, En(i)
        END DO

        IF( l == l_ini ) CALL WRITE_WF(nfun,Hij(:,n0_ini))

   ifg: IF( KIND_PI == 1 .OR. KIND_PI == 2 ) THEN

   iflLV: IF( l == l_ini ) THEN
            E_ini(:) = En(:)
            ci_ini(:) = Hij(:,n0_ini)
          ELSE IF( l == l_fin ) THEN iflLV
            E_fin = En(:)
            IF( Emax_fin == -1.0D0 ) Emax_fin = En(nfun)
            i = 1
dsearchnLV: DO
              IF( En(i) < 0.D0 ) n0_fin = i
              IF( En(i) <= Emax_fin ) n1_fin = i
              i = i + 1
              IF( i > nfun ) EXIT dsearchnLV
            END DO dsearchnLV
            n0_fin = n0_fin + 1
            n0_fin = MIN(n0_fin,nfun-1)
            WRITE(6,'(/,A26,I2,A3,X,2I4)') 'LIMITS FOR FINAL STATE (l=', l_fin, &
     &                                      ') :', n0_fin, n1_fin
            ALLOCATE( ci_fin(nfun,n0_fin:n1_fin) )
            ci_fin(:,n0_fin:n1_fin) = Hij(:,n0_fin:n1_fin)
          END IF iflLV

        ELSE IF( KIND_PI >= 3 ) THEN ifg

          Enl(:,l) = En(:)
          IF( Emax_fin == -1.0D0 ) THEN
            Emax_fin = En(nfun)
            Elim = Emax_fin
          ELSE
            Elim = Emax_fin + 0.25D0
            IF( KIND_PI >= 8 ) Elim = Emax_fin
          END IF
          
          i = 1
          nbold = 0
dsearchnA:DO
            IF( En(i) < 0.D0 ) THEN
              n0_fin = i
              nbold = nbold + 1
            END IF
            IF( En(i) <= Emax_fin ) n1_fin = i
            IF( En(i) <= Elim ) ntemp = i
            IF( En(i) > Emax_fin .AND. En(i) > Elim ) EXIT dsearchnA
            i = i + 1
            IF( i > nfun ) EXIT dsearchnA
          END DO dsearchnA
          nbds = MAX(nbds,nbold)
          n0_fin = n0_fin + 1
          n1_fin = n1_fin + 1
          nE0 = n0_fin
          IF( KIND_PI >= 5 ) n0_fin = 1
          IF( ntemp > nlim ) nlim = ntemp
          n01(l,1) = n0_fin
          n01(l,2) = n1_fin
          n01(l,3) = nE0 - 1
          WRITE(6,'(/,A23,I3)') 'NUMBER OF BOUND STATES:', nbold
          WRITE(6,'(A14,I3,A6,2I5)') 'LIMITS FOR l =', l, ' STATE:', &
     &                                  n0_fin+l, n1_fin+l
          ntemp = MAX(n1_fin+40,nlim)
          ntemp = MIN(ntemp,nfun)
          IF( l == 0 ) THEN
            ALLOCATE( ctemp(nfun,ntemp,0:lmax) )
            ctemp = 0.D0
          END IF
          ctemp(:,1:ntemp,l) = Hij(:,1:ntemp)

    !     DENSITY OF STATES

          DO i = nE0+1, nfun-1
            rEki(i,l) = SQRT(2.D0 / (En(i+1) - En(i-1)))
          END DO
          rEki(nE0,l) = SQRT(1.D0 / (En(nE0+1) - En(nE0)))
          rEki(nfun,l) = SQRT(1.D0 / (En(nfun) - En(nfun-1)))

!         IF( l == l_ini ) ctemp(:,n0_ini,l) = Hij(:,n0_ini)

        END IF ifg

      END DO dol0

      CLOSE(75)

      n1_max = n1_fin
      IF( KIND_PI >= 3 ) THEN

        n1_max = MAX(MAXVAL(n01(:,2))+20,nlim)
        n1_max = MIN(n1_max,nfun)
        WRITE(6,'(A8,I5)') 'n1_max =', n1_max

        IF( KIND_EGR == 1 .AND. n1_max < nEpts ) THEN
          WRITE(6,'(A40)') 'WARNING! Not enough continuum states!!!'
          WRITE(6,'(A25)') 'Modifying Number of E Pts'
          nEpts = nEpts - nbds
          WRITE(6,'(A10,I4)') 'New nEpts=', nEpts
        END IF

        OPEN( UNIT=80, FILE='Eigenvec_All.dat', ACTION='WRITE' )
        WRITE(80,*) nfun, n1_max, lmax

        ALLOCATE( cinl(nfun,n1_max,0:lmax) )
        cinl = 0.D0

        DO l = 0, lmax
          cinl(1:nfun,1:n1_max,l) = ctemp(1:nfun,1:n1_max,l)
          WRITE(80,*) l
          DO ni = 1, n1_max
            WRITE(80,300) ni, (cinl(i,ni,l), i=1,nfun)
          END DO
        END DO

        DEALLOCATE( ctemp )

!        CALL CHKPHS
!        CALL WRITEWF(n01(lmax,1),n1_max,lmax)
      END IF

      DEALLOCATE( Hij, Bij, En, WORK )

100      FORMAT(/,T2,A5,I2)
110      FORMAT(T2,A26,/)
120      FORMAT(T5,A1,T9,A11)
200      FORMAT(T2,I4,T8,G22.15)
300      FORMAT(I5,5000G20.10)

      END SUBROUTINE SOLVE_SYSTEM
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE CHKPHS
      
!      Subroutine to verify that all the WFs have the same phase (positive by choice)
!      Dr. Carlos Granados, MLU Halle
!      08.01.2018

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_PHOTOION
      
      IMPLICIT NONE
      
      INTEGER :: i, j, jmin, jmax, jfun, l, n, nr, left, mflag
      REAL(DP) :: r0, r1, r, dr, sumf, bsp(k), fr(3)
      
      nr = 3
      r0 = ra
      r1 = 0.1D0
      dr = (r1 - r0) / DBLE(nr)

      fr = 0.D0

dol0: DO l = 0, lmax
  don0: DO n = 1, n1_max

    dor0: DO i = 1, nr
      
            r = r0 + DBLE(i) * dr
        
            bsp = 0.D0
            CALL interv(rt,nkp,r,left,mflag)
            CALL bsplvb(nkp,rt,k,1,r,left,bsp)

            jmin = left - nbc1 + 1
            jmax = MIN(jmin + k - 1, nfun)

            sumf = 0.D0
       doj: DO j = jmin, jmax
              jfun = j - (left-nbc1)
              sumf = sumf + cinl(j,n,l) * bsp(jfun)
            END DO doj
            fr(i) = sumf
            
          END DO dor0
          
          IF( fr(1) < 0.D0 .AND. fr(2) < 0.D0 .AND. fr(3) < 0.D0 ) cinl(:,n,l) = -cinl(:,n,l)
            
        END DO don0
      END DO dol0

      END SUBROUTINE CHKPHS
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ZVTMV(il,ni,li,mi,jl,nf,lf,mf,n,ncomp,coeffi,zsumall)

      USE MOD_TYPES
      USE MOD_MATRICES, ONLY: Sij
      USE MOD_PHOTOION, ONLY: zAij, cinl, B0z

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: il, ni, li, mi, jl, nf, lf, mf, n, ncomp
      REAL(DP), INTENT(IN) :: coeffi(ncomp)

      COMPLEX(DPC), INTENT(OUT) :: zsumall(ncomp)

      INTEGER :: i, ibra, jket
      REAL(DP) :: ci, cj
      COMPLEX(DPC) :: zTij

      zsumall = 0.D0
      DO ibra = 1, n

        ci = cinl(ibra,nf,lf)

        DO jket = 1, n
        
          cj = cinl(jket,ni,li)
        
          DO i = 1, ncomp-1
            zTij = coeffi(i) * zAij(ibra,jket,il,jl,i)
            zsumall(i) = zsumall(i) + ci * zTij * cj
          END DO
          
          IF( li == lf .AND. mi == mf .AND. mi /= 0 .AND. B0z /= 0.D0 ) THEN
            zTij = coeffi(ncomp) * mi * Sij(ibra,jket)
            zsumall(ncomp) = zsumall(ncomp) + ci * zTij * cj
          END IF
          
        END DO
        
      END DO

      END SUBROUTINE ZVTMV
