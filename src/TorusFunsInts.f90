	SUBROUTINE PIFUNS(lf,zPmq)

!	SUBROUTINE TO CALCULATE THE ANGULAR TERMS IN THE CURRENT DENSITY, USING
!	  THE SPECTRAL METHOD (NOT INTEGRATED MATRIX ELEMENTS)
!	  CALLED FROM TORMAT, IN BSP_ATOM
!	C. M. GRANADOS--CASTRO
!	10-11.2018
!	V.2.0 01-02.11.2018: INTEGRALS SIMPLIFIED

	USE MOD_TYPES
	USE MOD_GRID, ONLY: nfib, nfib0, tht, pht
	USE MOD_MATRICES, ONLY: zrangij
	USE MOD_PHOTOION, ONLY: nlm, lmf
	USE MOD_FTW_PH

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: lf
	COMPLEX(DPC), INTENT(OUT) :: zPmq(nlm,nlm,-1:1,0:nfib0,4)

	INTEGER :: l1, m1, l2, m2, l3, m3, li, mi, lj, mj, mq, i, ith, il, jl, lmin, lmax
	REAL(DP) :: th, ph, ca, cb, THREE_J
	COMPLEX(DPC) :: zf123i, zIf1, zIf2, zIf3, zYi, zYj
	COMPLEX(DPC) :: zsum11, zsum12, zsum13, zsum21, zsum22, zsum23

	COMPLEX(DPC), ALLOCATABLE :: zYlm(:,:,:) !zPivec(:,:,:,:,:,:,:,:,:), 
	COMPLEX(DPC), ALLOCATABLE :: zf1(:), zf2(:), zf3(:), zfY(:,:)

	WRITE(6,'(/,A25)') 'Calculating Angular Terms'

!	ALLOCATE( zPivec(0:lf,-lf:lf,0:lf,-lf:lf,0:lf+1,-lf-1:lf+1,-1:1,0:nfib0,-1:1) )
!	zPivec = 0.D0

!	ALLOCATE( zf1(0:nfib0), zf2(0:nfib0), zf3(0:nfib0) )

	ALLOCATE( zYlm(0:(lf+1),-lf-1:lf+1,0:nfib0), zfY(0:lf+1,-lf-1:lf+1) )
	zYlm = 0.D0
	zfY = 0.D0
	
	DO ith = 0, nfib0
	  th = tht(ith)
	  ph = pht(ith)
	  CALL Ylm_All(lf+1,th,ph,zfY)
	  zYlm(:,:,ith) = zfY(:,:)
	END DO

	DEALLOCATE( zfY )

	WRITE(6,'(A20)') 'Angular-Radial Terms'

!$OMP PARALLEL DEFAULT(PRIVATE), SHARED(lf,nfib,nfib0,tht,pht,zYlm,zrangij)

	ALLOCATE( zf1(0:nfib0), zf2(0:nfib0), zf3(0:nfib0) )

!$OMP DO

dol1:	DO l1 = 0, lf
	  DO m1 = -l1, l1
	  
    dol2: DO l2 = 0, lf
	      DO m2 = -l2, l2
	      
	 doFib: DO ith = 0, nfib0
	              
		    th = tht(ith)
	          ph = pht(ith)

                zYi = CONJG(zYlm(l1,m1,ith))
                zYj = zYlm(l2,m2,ith)

		    zf1(ith) = zYi * SIN(th) * COS(ph) * zYj
		    zf2(ith) = zYi * SIN(th) * SIN(ph) * zYj
		    zf3(ith) = zYi * COS(th) * zYj

		  END DO doFib
			
		  CALL FIBINT(nfib,nfib0,zf1,zIf1)
		  CALL FIBINT(nfib,nfib0,zf2,zIf2)
		  CALL FIBINT(nfib,nfib0,zf3,zIf3)
		  
		  zrangij(l1,m1,l2,m2,1) = zIf1
		  zrangij(l1,m1,l2,m2,2) = zIf2
		  zrangij(l1,m1,l2,m2,3) = zIf3
		  
		END DO
	    END DO dol2
	    
	  END DO
	END DO dol1

!$OMP END DO

	DEALLOCATE( zf1, zf2, zf3 )

!$OMP END PARALLEL

!	DEALLOCATE( zf1, zf2, zf3 )

	zPmq = 0.D0

	WRITE(6,'(A28)') 'Angular Sums and Total Terms'

!$OMP PARALLEL DEFAULT(PRIVATE), SHARED(nlm,lmf,nfib,nfib0,zYlm,zPmq)
!$OMP DO

doith:DO ith = 0, nfib0

  doli: DO il = 1, nlm

	    li = lmf(il,1)
	    mi = lmf(il,2)

    dolj: DO jl = 1, nlm
  
  	      lj = lmf(jl,1)
  	      mj = lmf(jl,2)

	      l1 = lj + 1
	      l2 = lj - 1

    dosumq: DO mq = -1, 1

	        zsum11 = 0.D0
	        zsum12 = 0.D0
	        zsum21 = 0.D0
	        zsum22 = 0.D0
      dosum1: DO m1 = -l1, l1			!SUM ON m
		    ca = THREE_J(l1,1,lj,m1,mq,mj)
		    cb = THREE_J(l1,1,lj,m1,mq,-mj)
        	    zsum11 = zsum11 + ca * zYlm(li,-mi,ith) * zYlm(l1,m1,ith)
		    zsum12 = zsum12 + cb * zYlm(li,mi,ith) * zYlm(l1,m1,ith)
     		  END DO dosum1
      dosum2: DO m2 = -l2, l2			!SUM ON m
     		    ca = THREE_J(l2,1,lj,m2,mq,mj)
		    cb = THREE_J(l2,1,lj,m2,mq,-mj)
        	    zsum21 = zsum21 + ca * zYlm(li,-mi,ith) * zYlm(l2,m2,ith)
	          zsum22 = zsum22 + cb * zYlm(li,mi,ith) * zYlm(l2,m2,ith)
     		  END DO dosum2

     	        zPmq(il,jl,mq,ith,1) = ((-1)**(mi+mj)) * zsum11 !PI(li,1,lj+1,-mi,q,mj)
     	        zPmq(il,jl,mq,ith,2) = zsum12 !PI(li,1,lj+1,mi,q,-mj)

     	        zPmq(il,jl,mq,ith,3) = ((-1)**(mi+mj)) * zsum21 !PI(li,1,lj-1,-mi,q,mj)
     	        zPmq(il,jl,mq,ith,4) = zsum22 !PI(li,1,lj-1,mi,q,-mj)

	      END DO dosumq
     		    
     	    
     	    END DO dolj
     	  END DO doli
     	  
     	END DO doith

!$OMP END DO
!$OMP END PARALLEL

	DEALLOCATE( zYlm )

	END SUBROUTINE PIFUNS
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE FRINT(il,ni,li,ml,nm,lm,nl,nn,ln,jl,nj,lj,zimnj,ximnj,A0imnj, &
     &                 zArximnj,zAryimnj,zAzimnj,zArxjmni,zAryjmni,zAzjmni)

!	SUBROUTINE TO CALCULATE RADIAL INTEGRALS...
!	C. M. GRANADOS--CASTRO
!	UNI-HALLE, 10-11.2018
!	V.2.0 01-02.11.2018: INTEGRALS SIMPLIFIED
!     V.3.0 21.11.2018: ANOTHER ATTEMPT TO MAKE IT MORE EFFICIENT. USING THE PROPERTY
!       OF B-SPLINES THAT FOR EACH r THERE ARE k BASIS FUNCTIONS...


	USE MOD_TYPES
	USE MOD_GRID, ONLY: xg, wg
	USE MOD_BSPLINES
	USE MOD_MATRICES, ONLY: fur, dfur!, zJAr
	USE MOD_PHOTOION, ONLY: cinl

	IMPLICIT NONE

      INTEGER, INTENT(IN) :: il, ni, li, ml, nm, lm, nl, nn, ln, nj, jl, lj
      REAL(DP), INTENT(OUT) :: zimnj, ximnj, A0imnj
      COMPLEX(DPC), INTENT(OUT) :: zArximnj, zAryimnj, zAzimnj
      COMPLEX(DPC), INTENT(OUT) :: zArxjmni, zAryjmni, zAzjmni

	INTEGER :: ibra, mbra, nket, jket, ibet, ibetmin, ibetmax, igl
	INTEGER :: i, m, n, j, left, imin, imax
	REAL(DP) :: r, dr, f1, f2, fr, gr, Bi, Bm, Bn, Bj, dBn, sum1, sum2, sum3
	REAL(DP) :: ap, am, cii, cmm, cnn, cjj
	COMPLEX(DPC) :: zA0, zAr1, zAr2, zAr3, ztAr1, ztAr2, ztAr3

      REAL(DP), ALLOCATABLE :: bsp(:), dbsp(:), bsp1(:), bspp(:)

!	nr = nointv * ka

	ap = ((-1)**ln) * SQRT(3.D0*(ln+1)/DBLE(2*ln+1))
	am = ((-1)**(ln-1)) * SQRT(3.D0*(ln)/DBLE(2*ln+1))

      ALLOCATE( bsp(k), dbsp(k), bsp1(k-1), bspp(k+1) )

      zimnj = 0.D0
      ximnj = 0.D0
      sum1 = 0.D0
      sum2 = 0.D0
      A0imnj = 0.D0
      zArximnj = 0.D0
      zAryimnj = 0.D0
      zAzimnj = 0.D0
      zArxjmni = 0.D0
      zAryjmni = 0.D0
      zAzjmni = 0.D0
dobet:DO ibet = 1, nkp-1
            
        f1 = (rt(ibet+1) + rt(ibet)) / 2.D0
        f2 = (rt(ibet+1) - rt(ibet)) / 2.D0
                
  doGL: DO igl = 1, ka
                
          r = f1 + xg(igl) * f2
          dr = f2 * wg(igl)
                  
          CALL BSPALL(r,left,bsp,dbsp,bsp1,bspp)
                  
          fr = 1.D0 / (r**2)
          gr = 1.D0 / (r**3)
          
          imin = MAX(left-k+1,1)
          imax = MIN(left,nfun)

          DO ibra = imin, imax
            cii = cinl(ibra,ni,li)
            i = ibra - imin + 1
            Bi = cii * bsp(i)
            DO mbra = imin, imax
              cmm = cinl(mbra,nm,lm)
              m = mbra - imin + 1
              Bm = cmm * bsp(m)
              DO nket = imin, imax
                cnn = cinl(nket,nn,ln)
                n = nket - imin + 1
                Bn = cnn * bsp(n)
                dBn = cnn * dbsp(n)
                DO jket = imin, imax
                  cjj = cinl(jket,nj,lj)
                  j = jket - imin + 1
                  Bj = cjj * bsp(j)

                  sum1 = sum1 + Bi * Bm * fr * dBn * Bj * dr
                  sum2 = sum2 + Bi * Bm * gr * Bn * Bj * dr

                  A0imnj = A0imnj + Bi * Bm * zA0 * Bn * Bj * dr

                  zArximnj = zArximnj + Bi * Bm * zAr1 * fr * Bn * Bj * dr
                  zAryimnj = zAryimnj + Bi * Bm * zAr2 * fr * Bn * Bj * dr
                  zAzimnj = zAzimnj + Bi * Bm * zAr3 * fr * Bn * Bj * dr

                  zArxjmni = zArxjmni + Bi * Bm * ztAr1 * fr * Bn * Bj * dr
                  zAryjmni = zAryjmni + Bi * Bm * ztAr2 * fr * Bn * Bj * dr
                  zAzjmni = zAzjmni + Bi * Bm * ztAr3 * fr * Bn * Bj * dr

                END DO
              END DO
            END DO
          END DO
          
        END DO doGL
      END DO dobet

      zimnj = ap * (sum1 - (ln+1)*sum2)
      ximnj = am * (sum1 + ln*sum2)

! $OMP PARALLEL DEFAULT(PRIVATE) SHARED(nfun,cinl,k,ka,xg,wg,rt,ap,am,il,ni,li, &
! $OMP &  ml,nm,lm,nl,nn,ln,jl,nj,lj,zJAr,zimnj,ximnj,A0imnj,zArximnj,zAryimnj, &
! $OMP &  zAzimnj,zArxjmni,zAryjmni,zAzjmni)

!      ALLOCATE( bsp(k), dbsp(k), bsp1(k-1), bspp(k+1) )

! $OMP DO REDUCTION(+:zimnj,ximnj,A0imnj,zArximnj,zAryimnj,zAzimnj,zArxjmni,     &
! $OMP &  zAryjmni,zAzjmni)

	END SUBROUTINE FRINT
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE FRMATINT(ir,r)

!	SUBROUTINE TO CALCULATE RADIAL MATRIX ELEMENTS (NOT INTEGRATED), NEEDED
!	  TO CALCULATE THE CURRENT DENSITY. INVOLVES ONLY THE GRADIENT OF THE
!	  REDUCED RADIAL WAVE FUNCTIONS. USE BSPLINES AND EXPANSION COEFFICIENTS
!	  DIRECTLY, WITHOUT USING THE FULL WAVE FUNCTIONS
!	  CALLED FROM TORMAT, IN BSP_ATOM
!	C. M. GRANADOS--CASTRO
!	UNI-HALLE, 10-11.2018
!	V.1.0 04.11.2018: INTEGRALS SIMPLIFIED
!	V.1.1 04.11.2018: MOVED TO CALCULATE FUNCTION FOR A GIVEN r.
!				THE LOOP ON r MUST BE PERFORMED EXTERNALLY...
!				THE MATRIX ELEMENTS OF r ARE NO LONGER HERE...

	USE MOD_TYPES
	USE MOD_GRID, ONLY: lmax
	USE MOD_BSPLINES, ONLY: rt, nkp, k, nfun, BSPALL
	USE MOD_PHOTOION, ONLY: n1_max, cinl
!	USE MOD_MATRICES, ONLY: fJrij

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: ir
	REAL(DP), INTENT(IN) :: r

	INTEGER :: i, j, ibra, jket, jp, ni, li, nj, lj, left, igl, nr
	INTEGER :: ifun, jfun
	REAL(DP) :: A1, A2, b1, b2, ci, cj, BiBj, BidBj, sumfr, sumgr
	REAL(DP), ALLOCATABLE :: bsp(:), dbsp(:), bsp1(:), bspp(:)

	ALLOCATE( bsp(k), dbsp(k), bsp1(k-1), bspp(k+1) )

	CALL BSPALL(r,left,bsp,dbsp,bsp1,bspp)

!	RADIAL TERMS IN GAUSS-LEGENDRE QUADRATURE

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n1_max,lmax,left,k,cinl,bsp,dbsp) !,fJrij)
!$OMP DO

doni: DO ni = 1, n1_max
  doli: DO li = 0, lmax

    donj: DO nj = 1, n1_max
	dolj: DO lj = 0, lmax

		  sumfr = 0.D0
		  sumgr = 0.D0
		  DO i = 1, nfun
		    ifun = i - (left-k)
		    DO j = 1, k
		      jfun = j - (left-k)

			BiBj = bsp(i) * bsp(j)
			BidBj = bsp(i) * dbsp(j)
			ci = 0.D0
			cj = 0.D0
			IF( (ifun >= 1 .AND. ifun <= nfun) .AND. &
     &		    (jfun >= 1 .AND. jfun <= nfun) ) THEN
			  ci = cinl(ifun,ni,li)
			  cj = cinl(jfun,nj,lj)
			END IF
			  
			sumfr = sumfr + ci * BidBj * (1.D0 / (r**2)) * cj
			sumgr = sumgr + ci * BiBj * (1.D0 / (r**3)) * cj
		    END DO
		  END DO
		    
!		  fJrij(ni,li,nj,lj,1,1) = sumfr
!		  fJrij(ni,li,nj,lj,1,2) = sumgr

		END DO dolj
	    END DO donj
		
	  END DO doli
	END DO doni

!$OMP END DO
!$OMP END PARALLEL

!	DO ni = 1, n1_max
!	  DO li = 0, lmax
!	    DO nj = 1, n1_max
!		DO lj = 0, lmax

!		  WRITE(90,150) ni, li, nj, lj, ir, fJrij(ni,li,nj,lj,1,1), &
!     &			    fJrij(ni,li,nj,lj,1,2)
!     
!		END DO
!	    END DO
!	  END DO
!	END DO

150	FORMAT(4I4,I6,2G20.10)

	DEALLOCATE( bsp, dbsp, bsp1, bspp )

	END SUBROUTINE FRMATINT
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE SPHVEC(mq,th,phi,zereqvec)

!	COEFFICIENTS NEEDED TO CALCULATE THE COMPONENTS OF er X (er X eq), IN
!	  SPHERICAL COMPONENTS

	USE MOD_TYPES
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: mq
	REAL(DP), INTENT(IN) :: th, phi
	COMPLEX(DPC), INTENT(OUT) :: zereqvec(3)

	REAL(DP) :: x, y
	COMPLEX(DPC) :: ephp, ephm, zc1, zc2, zc3

	x = COS(th)
	y = SIN(th)
	ephp = EXP(zi*phi)
	ephm = EXP(-zi*phi)
	
	IF( mq == -1 ) THEN

	  zc1 = -1.D0 + 0.5D0*(y**2)
	  zc2 = ephm * x * y / SQRT(2.D0)
	  zc3 = -((ephm*y)**2) / 2.D0
	  
	ELSE IF( mq == 0 ) THEN
	
	  zc1 = ephp * x / SQRT(2.D0)
	  zc2 = -y
	  zc3 = -ephm * x / SQRT(2.D0)
	
	ELSE IF( mq == 1 ) THEN
	
	  zc1 = -((ephp*y)**2) / 2.D0
	  zc2 = ephp * x * y / SQRT(2.D0)
	  zc3 = -1.D0 + 0.5D0*(y**2)

	ELSE
	
	  zc1 = 0.D0
	  zc2 = 0.D0
	  zc3 = 0.D0
	  
	END IF
	
	zereqvec = (/zc1, zc2, zc3/)

	END SUBROUTINE SPHVEC
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE INT3D(nr,zsumjr)

      USE MOD_TYPES
      USE MOD_GRID, ONLY: nfib, nfib0, rtot, wtot
      USE MOD_MATRICES, ONLY: psirt, zApsirt
      USE MOD_FTW_PH, ONLY: FIBINT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nr
      COMPLEX(DPC), INTENT(OUT) :: zsumjr(6)

      INTEGER :: ir, ith, mq
      REAL(DP) :: r, dr, fr
      COMPLEX(DPC) :: zpsi, cpsi, zIfr, zf2, zc0, zdfx, zdfy, zdfz
      COMPLEX(DPC) :: cdfx, cdfy, cdfz, zfAx, zfAy, zfAz, zfbAx, zfbAy, zfbAz
      COMPLEX(DPC) :: zpsib, cpsib, zdfbx, zdfby, zdfbz, cdfbx, cdfby, cdfbz
      COMPLEX(DPC), ALLOCATABLE :: zfrth(:,:), zf(:)

      zsumjr = 0.D0
      zc0 = -0.5D0 * zi

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nr,rtot,wtot,nfib,nfib0,psirt,zApsirt,zc0, &
!$OMP &  zsumjr)

      ALLOCATE( zfrth(0:nfib0,6), zf(0:nfib0) )

!$OMP DO REDUCTION(+:zsumjr)

      DO ir = 1, nr
      
        r = rtot(ir)
        IF( r == 0.D0 ) r = eps
        dr = wtot(ir)

        fr = 1.D0 !/ (r**2)

        DO ith = 0, nfib0

!         COMPLETE WAVEPACKET
        
          zpsi = psirt(ir,ith,0,1)
          cpsi = CONJG(zpsi)

          zdfx = -(psirt(ir,ith,1,2) - psirt(ir,ith,-1,2)) / SQRT(2.D0)
          zdfy = zi * (psirt(ir,ith,1,2) + psirt(ir,ith,-1,2)) / SQRT(2.D0)
          zdfz = psirt(ir,ith,0,2)

          cdfx = CONJG(zdfx)
          cdfy = CONJG(zdfy)
          cdfz = CONJG(zdfz)

          zfAx = cpsi * zApsirt(ir,ith,1) * fr * zpsi
          zfAy = cpsi * zApsirt(ir,ith,2) * fr * zpsi
          zfAz = cpsi * zApsirt(ir,ith,3) * fr * zpsi
          
          zfrth(ith,1) = zc0 * (cpsi * zdfx - zpsi * cdfx) + zfAx
          zfrth(ith,2) = zc0 * (cpsi * zdfy - zpsi * cdfy) + zfAy
          zfrth(ith,3) = zc0 * (cpsi * zdfz - zpsi * cdfz) + zfAz
          
!         BOUND STATES
          
          zpsib = psirt(ir,ith,0,3)
          cpsib = CONJG(zpsib)
          
          zdfbx = -(psirt(ir,ith,1,4) - psirt(ir,ith,-1,4)) / SQRT(2.D0)
          zdfby = zi * (psirt(ir,ith,1,4) + psirt(ir,ith,-1,4)) / SQRT(2.D0)
          zdfbz = psirt(ir,ith,0,4)

          cdfbx = CONJG(zdfbx)
          cdfby = CONJG(zdfby)
          cdfbz = CONJG(zdfbz)

          zfbAx = cpsib * zApsirt(ir,ith,1) * fr * zpsib
          zfbAy = cpsib * zApsirt(ir,ith,2) * fr * zpsib
          zfbAz = cpsib * zApsirt(ir,ith,3) * fr * zpsib
          
          zfrth(ith,4) = zc0 * (cpsib * zdfbx - zpsib * cdfbx) + zfbAx
          zfrth(ith,5) = zc0 * (cpsib * zdfby - zpsib * cdfby) + zfbAy
          zfrth(ith,6) = zc0 * (cpsib * zdfbz - zpsib * cdfbz) + zfbAz

        END DO

        DO mq = 1, 6
          zf(:) = zfrth(:,mq)
          CALL FIBINT(nfib,nfib0,zf,zIfr)
          zsumjr(mq) = zsumjr(mq) + zIfr * dr
        END DO          
        
      END DO

!$OMP END DO

      DEALLOCATE( zfrth, zf )

!$OMP END PARALLEL

      END SUBROUTINE INT3D
