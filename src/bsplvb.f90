!----------------------------------------------------------------------
! Subroutine BSPLV from the Book 'A practical guide to splines'
! Carl de Boor
! Given a radial point r and the index i in the knot series, the
! routine determines the value of all non-zero splines
! MODIFIED ON 18.01.2018
! HAVING SOME PROBLEMS WHILE RUNNING ON PARALLEL CODES
!----------------------------------------------------------------------

	SUBROUTINE BSPLVB(ndim,t,jhigh,index,x,left,biatx)

	USE MOD_TYPES

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: ndim, jhigh, index, left
	REAL(DP), INTENT(IN) :: t(ndim), x

	REAL(DP), INTENT(INOUT) :: biatx(jhigh)

	INTEGER :: i, j
	REAL(DP) :: saved, term, deltal(ndim), deltar(ndim)

	IF( index == 1 ) THEN
	  j = 1
	  biatx(1) = 1.00D0
	  IF( jhigh <= j ) RETURN
	END IF

	IF( t(left+1) <= t(left) ) THEN
	  WRITE(6,*) 'FATAL ERROR - BSPLVB'
	  WRITE(6,*) left, t(left+1), t(left)
	  STOP
	END IF

	j = 1

 dok: DO
	  deltar(j) = t(left+j) - x
	  deltal(j) = x - t(left+1-j)
	  saved = 0.0D0
	  DO i = 1, j
	    term = biatx(i) / ( deltar(i) + deltal(j+1-i) )
	    biatx(i) = saved + deltar(i) * term
	    saved = deltal(j+1-i) * term
	  END DO
	  biatx(j+1) = saved
	  j = j + 1
	  IF( jhigh <= j ) EXIT dok
	END DO dok

	END SUBROUTINE BSPLVB
