	SUBROUTINE ANG_GRID

!	SUBROUTINE TO CALCULATE A LINEAR GRID ON TH...
!	THE GRID IS USED LATER TO CALCULATE B-SPLINES, AND TO PROJECT THE
!	  ASSOCIATED-LEGENDRE POLYNOMIAL ON B-SPLINES...
!
!	V 2.0 CALCULATE ANGULAR GRID TO BE USED ON THE FIBONACCI INTEGRATION METHOD

	USE MOD_TYPES
	USE MOD_BSPLINES
	USE MOD_GRID
	USE MOD_FAC_VARS
	USE MOD_FTW_PH
	USE MOD_PHOTOION
	USE MOD_MATRICES

	IMPLICIT NONE

	INTEGER :: n1, n2, i, l0, m0, lf, mf, m1, m2
	REAL(DP) :: th0, th1, dth, ph0, ph1, dph, A1, A2, z, zj
	COMPLEX(DPC), ALLOCATABLE :: zYlm_i(:,:)

	WRITE(6,'(/,A19)') 'Calculating Th-Grid'

	l0 = lmf(0,1)
	m0 = lmf(0,2)
	
	lf = lmf(nlm,1)
	mf = lmf(nlm,2)

	m1 = -lf
	m2 = lf
	IF( KIND_POT >= 5 .AND. bx /= 0.D0 ) THEN
	  m1 = m0
	  m2 = m1
	END IF

!	FACTORIALS NEEDED TO CALCULATE THE ASSOCIATED LEGENDRE FUNCTIONS
	ALLOCATE( fac(0:2*(lmax+1)) )
	fac = 1.D0
	DO i = 2, 2*(lmax+1)
	  fac(i) = DBLE(i) * fac(i-1)
	END DO

	ALLOCATE( tht(0:nfib0), pht(0:nfib0), zjt(0:nfib0) )
	ALLOCATE( zYlm(0:lf+1,-lf-1:lf+1,0:nfib0), zYlm_i(0:lf+1,-lf-1:lf+1) )
!	ALLOCATE( ciLP(0:nfib0,0:lf,m1:m2) , zciexp(0:nfib0,-2*mf:2*mf) )

!	ciLP = 0.D0
!	zciexp = 0.D0

	zYlm = 0.D0

 doi:	DO i = 0, nfib0

	  z = -1.D0 + DBLE(i)*dzf
	  ph1 = DBLE(i) * dphf

	  zj = z + (SIN(PI*z) / PI)

	  zjt(i) = zj
	  tht(i) = ACOS(zj)
	  pht(i) = ph1

	  CALL Ylm_All(lf+1,tht(i),pht(i),zYlm_i)
	  zYlm(0:lf+1,-lf-1:lf+1,i) = zYlm_i(0:lf+1,-lf-1:lf+1)

	END DO doi

	DEALLOCATE( zYlm_i )

	END SUBROUTINE ANG_GRID
