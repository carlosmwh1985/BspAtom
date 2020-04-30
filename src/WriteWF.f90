	SUBROUTINE WRITEWF(n0,n1,l0)

!	Subroutine to write down the calculated wave functions, for the indicated
!	  states
!	Dr. Carlos Granados, MLU Halle
!	08.01.2018

	USE MOD_TYPES
	USE MOD_GRID
	USE MOD_BSPLINES
	USE MOD_PHOTOION
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: n0, n1, l0
	
	INTEGER :: i, j, jmin, jmax, jfun, is, n, left, mflag
	REAL(DP) :: r, dr, sumf
	
	REAL(DP), ALLOCATABLE :: bsp(:), fr(:)
	
	OPEN( UNIT=30, FILE='WFs.dat', ACTION='WRITE' )

	nrpts = 10000
	dr = (rb - ra) / DBLE(nrpts)

	ALLOCATE( fr(0:n1-n0+1), bsp(k) )
	fr = 0.D0

	DO i = 0, nrpts

	  r = ra + DBLE(i) * dr

	  bsp = 0.D0
	  CALL interv(rt,nkp,r,left,mflag)
	  CALL bsplvb(nkp,rt,k,1,r,left,bsp)
	  
	  jmin = left - nbc1 + 1
	  jmax = MIN(jmin + k - 1, nfun)

	  DO is = 0, n1-n0+1

	    IF( is == 0 ) THEN
	      n = 1
	    ELSE
	      n = n + 1
	    END IF

	    sumf = 0.D0
	    DO j = jmin, jmax
	      jfun = j - (left-nbc1)
	      sumf = sumf + cinl(j,n,l0) * bsp(jfun)
	    END DO
	    fr(is) = sumf
	    
	    IF( is == 0 ) n = n0 - 1
	    
	  END DO

	  WRITE(30,100) r, (fr(j), j=0,n1-n0+1)

	END DO
	
	CLOSE(30)
	
100	FORMAT(100G20.10)

	END SUBROUTINE
