	SUBROUTINE CUBSPL(n0,x0,y0,n1,x1,y1)
	
	USE MOD_TYPES
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: n0, n1
	REAL(DP), INTENT(IN) :: x0(0:n0), y0(0:n0), x1(0:n1)
	
	REAL(DP), INTENT(OUT) :: y1(0:n1)

	INTEGER :: i
	REAL(DP) :: dy0i, dy0f, d2y0(0:n0), xi, yi, xref1, xref2

!	PREPARE THE INTERPOLATION:
!	1st DERIVATIVE IN THE FIRST AND LAST POINTS
	
	dy0i = (y0(1) - y0(0)) / (x0(1) - x0(0))
	dy0f = (y0(n0) - y0(n0-1)) / (x0(n0) - x0(n0-1))

!     2nd DERIVATIVES

	d2y0 = 0.D0
      CALL SPLINE(x0,y0,n0,dy0i,dy0f,d2y0)

!	INTERPOLATION POINT TO POINT

!	OPEN( UNIT=50, FILE='TestIntNew.dat', ACTION='WRITE' )

      xref1 = x0(0)
      xref2 = x0(n0)

	DO i = 0, n1

	  xi = x1(i)

        IF( xi == xref1 ) THEN
          yi = y0(0)
        ELSE IF( xi == xref2 ) THEN
          yi = y0(n0)
        ELSE            !CUBIC SPLINE INTERPOLATION
	    CALL SPLINT(x0,y0,d2y0,n0,xi,yi)
	  END IF
	  
	  y1(i) = yi

!	  WRITE(50,*) x1(i), y1(i)
	  
	END DO

	END SUBROUTINE CUBSPL
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)

	USE MOD_TYPES

	IMPLICIT NONE

	INTEGER :: n, nmax
      REAL(DP) :: yp1, ypn, x(0:n), y(0:n), y2(0:n)

	INTEGER :: i, k
	REAL(DP) :: p, qn, sig, un, u(0:n)
	
      IF( yp1 > .99D30) THEN
        y2(0) = 0.D0
        u(0) = 0.D0
      ELSE
        y2(0) = -0.5D0
        u(0) = ( 3.D0 / (x(1)-x(0)) ) * & 
     &	   ( (y(1)-y(0)) / (x(1)-x(0)) - yp1 )
      END IF

      DO i = 1, n-1
        sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
        p = sig * y2(i-1) + 2.D0
        y2(i) = ( sig - 1.D0 ) / p
        u(i) = ( 6.D0 * ( (y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1) ) &
     &        / (x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig * u(i-1) ) / p
      END DO 

      IF( ypn > .99D30) THEN
        qn = 0.D0
        un = 0.D0
      ELSE
        qn = 0.5D0
        un = ( 3.d0 / (x(n)-x(n-1)) ) * &
     &	 ( ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)) )
      END IF

      y2(n) = ( un - qn * u(n-1) ) / ( qn * y2(n-1) + 1.D0 )

      DO k = n-1, 0, -1
        y2(k) = y2(k) * y2(k+1) + u(k)
      END DO

	END SUBROUTINE SPLINE
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE SPLINT(xa,ya,y2a,n,x,y)

	USE MOD_TYPES

	IMPLICIT NONE

	INTEGER :: n
	REAL(DP) :: x, y, xa(0:n), y2a(0:n), ya(0:n)
	INTEGER :: k, khi, klo
	REAL(DP) :: a, b, h
	
	klo = 1
      khi = n
      DO WHILE( khi-klo > 1)
        k = ( khi + klo ) / 2     
        IF( xa(k) > x ) THEN
          khi = k
        ELSE
          klo = k
        END IF
      END DO
      h = xa(khi) - xa(klo)
      IF( h == 0.D0) PAUSE 'bad xa input in splint'
      a = ( xa(khi) - x ) / h
      b = ( x - xa(klo) ) / h
      y = a * ya(klo) + b * ya(khi) + &
     &    ( (a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.D0

	END SUBROUTINE SPLINT
