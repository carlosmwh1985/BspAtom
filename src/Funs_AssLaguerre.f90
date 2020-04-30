	FUNCTION Laguerre_L(n,l,x)

!	FUNCTION TO CALCUALTED THE ASSOCIATED LAGUERRE POLYNOMIAL L_l^n(x)

!	NOTE: A MODULE WITH THE NEEDED FACTORIALS MUST BE PROVIDED

!	USE MOD_FAC_VARS

	USE MOD_TYPES

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n, l
	REAL(DP), INTENT(IN) :: x

	REAL(DP) :: Laguerre_L

	INTEGER :: i, j
	REAL(DP) :: c1, c2, c3, suml
	REAL(DP), ALLOCATABLE :: Lnl(:,:)

	IF( l < 0 ) THEN
	  WRITE(6,*) 'Not possible to calculate Associated Laguerre Polynomial'
	  Laguerre_L = 0.D0
	  RETURN
	END IF

	IF( ALLOCATED(Lnl) ) DEALLOCATE( Lnl )
	ALLOCATE( Lnl(-2:MAX(n,0),0:l) )

	Lnl = 0.D0
	Lnl(0,:) = 1.D0

 don:	DO i = 1, n
   dol: DO j = 0, l

	    c1 = 2.D0*i - 1.D0 + j - x
	    c2 = i + j - 1.D0

	    Lnl(i,j) = ( c1*Lnl(i-1,j) - c2*Lnl(i-2,j) ) / DBLE(i)

	  END DO dol
	END DO don

	Laguerre_L = Lnl(n,l)

!	suml = 0.D0
!	DO i = 0, n
!	  c1 = fac(n+l) / fac(n-i)
!	  c2 = 1.D0 / fac(l+i)
!	  c3 = 1.D0 / fac(i)
!	  suml = suml + ((-1.D0)**i) * c1 * c2 * c3 * (x**i)
!	END DO
!	Laguerre_L = suml

	END FUNCTION Laguerre_L
