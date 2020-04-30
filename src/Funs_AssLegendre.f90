	FUNCTION Legendre_P(l,m,x)

!	FUNCTION TO CALCULATE THE ASSOCIATED LEGENDRE POLYNOMIAL P_l^m(x)
!	NOT NORMALIZED FUNCTIONS

	USE MOD_TYPES
	USE MOD_FAC_VARS
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: l, m
	REAL(DP), INTENT(IN) :: x

	REAL(DP) :: Legendre_P

	INTEGER :: i, j, mp
	REAL(DP), ALLOCATABLE :: Pn(:,:)
	
	ALLOCATE( Pn(0:l,-l:l) )
	Pn = 0.D0
	
	Pn(0,0) = 1.D0
	
	IF( l >= 1 ) THEN
	  Pn(1,0) = x * Pn(0,0)
	  Pn(1,1) = DSQRT(1.D0-(x**2)) * Pn(0,0)
	  Pn(1,-1) = -1.D0/2.D0 * Pn(1,1)
	  DO i = 2, l
	    Pn(i,0) = (1.D0/DBLE(i))*(DBLE(2*i-1)*x*Pn(i-1,0)-DBLE(i-1)*Pn(i-2,0))
	    DO mp = 1, i
	      Pn(i,mp) = DBLE(2*i-1)*DSQRT(1.D0-(x**2))*Pn(i-1,mp-1) + Pn(i-2,mp)
	      Pn(i,-mp) = -(fac(i-mp)/fac(i+mp)) * Pn(i,mp)
	    END DO
	  END DO
	END IF	
	
	Legendre_P = Pn(l,m)

	DEALLOCATE( Pn )
	
	END FUNCTION Legendre_P
