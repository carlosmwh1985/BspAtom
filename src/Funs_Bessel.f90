	FUNCTION BESSJ0(x)
!	Returns the Bessel function J0(x) for any real x.
	USE MOD_TYPES

	IMPLICIT NONE

	REAL(DP) :: BESSJ0, x
	
!	Returns the Bessel function J0(x) for any real x.

	REAL(DP) :: ax, xx, y, z
     
	REAL(DP), PARAMETER :: p1=1.D0, p2=-0.1098628627D-2, p3=0.2734510407D-4
	REAL(DP), PARAMETER :: p4=-0.2073370639D-5, p5=0.2093887211D-6
	REAL(DP), PARAMETER :: q1=-0.1562499995D-1, q2=0.1430488765D-3, q3=-0.6911147651D-5
	REAL(DP), PARAMETER :: q4=0.7621095161D-6, q5=-0.934945152D-7
	REAL(DP), PARAMETER :: r1=57568490574.D0, r2=-13362590354.D0, r3=651619640.7D0
	REAL(DP), PARAMETER :: r4=-11214424.18D0, r5=77392.33017D0, r6=-184.9052456D0
	REAL(DP), PARAMETER :: s1=57568490411.D0, s2=1029532985.D0, s3=9494680.718D0
	REAL(DP), PARAMETER :: s4=59272.64853D0, s5=267.8532712D0, s6=1.D0

	IF( ABS(x) < 8.0D0) THEN	! Direct rational function fit.
	  y=x**2
	  BESSJ0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
     & /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
	ELSE					! Fitting function (6.5.9).
	  ax = ABS(x)
	  z = 8.0D0 / ax
	  y = z**2
	  xx = ax - 0.785398164D0
	  BESSJ0 = SQRT(0.636619772D0 / ax) * (COS(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
     &  *p5)))) - z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
	END IF

	END FUNCTION BESSJ0
!
!-------------------------------------------------------------------------------
!
	FUNCTION BESSJ1(x)
!	Returns the Bessel function J1(x) for any real x.
	USE MOD_TYPES

	IMPLICIT NONE
	
	REAL(DP) BESSJ1, x

	REAL(DP) ax, xx, y, z
	
	REAL(DP), PARAMETER :: p1=1.D0, p2=.183105D-2, p3=-.3516396496D-4
	REAL(DP), PARAMETER :: p4=.2457520174D-5, p5=-.240337019D-6
	REAL(DP), PARAMETER :: q1=.04687499995D0, q2=-.2002690873D-3, q3=.8449199096D-5
	REAL(DP), PARAMETER :: q4=-.88228987D-6, q5=.105787412D-6
	REAL(DP), PARAMETER :: r1=72362614232.D0, r2=-7895059235.D0, r3=242396853.1D0
	REAL(DP), PARAMETER :: r4=-2972611.439D0, r5=15704.48260D0, r6=-30.16036606D0
	REAL(DP), PARAMETER :: s1=144725228442.D0, s2=2300535178.D0, s3=18583304.74D0
	REAL(DP), PARAMETER :: s4=99447.43394D0, s5=376.9991397D0, s6=1.D0
	
	IF( ABS(x) < 8.0D0 ) THEN	! Direct rational function fit.
	  y=x**2
	  BESSJ1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
     & /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
	ELSE					! Fitting function (6.5.9).
	  ax = ABS(x)
	  z = 8.0D0 / ax
	  y = z**2
	  xx = ax - 2.356194491D0
	  BESSJ1 = SQRT(0.636619772D0/ax) * (COS(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
     &  *p5)))) - z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))) * SIGN(1.D0,x)
     END IF
     
     END FUNCTION BESSJ1
