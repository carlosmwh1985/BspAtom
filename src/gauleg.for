*****************************************************
* Subroutine GAULEG, given the lower and upper limits
* of integration x1 and x2 and given n, the routine
* returns arrays x(1:n) and w(1:n), containing the
* abcissas and weights for the Gaussian-Legendre n-point
* quadrature formula.
*****************************************************
      SUBROUTINE gauleg(x1,x2,x,w,ndim,n)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION x(ndim),w(ndim)
cc    PARAMETER (EPS=3.d-14)
C     Anpassung der absoluten Rechengenauigkeit
      CALL PTINY(EPS)
      EPS=EPS*10
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
	PI=2*DASIN(1.0D0)
	z=cos(PI*(i-.25d0)/(n+.5d0))
1       continue
        p1=1.d0
	p2=0.d0
	  do 11 j=1,n
	    p3=p2
	    p2=p1
	    p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
	  pp=n*(z*p1-p2)/(z*z-1.d0)
	  z1=z
	  z=z1-p1/pp
	if(abs(z-z1).gt.EPS)goto 1
	x(i)=xm-xl*z
	x(n+1-i)=xm+xl*z
	w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
	w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software @1-.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Determing the machine precision   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTINY(EPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      HALFU=0.5D0
   50 TEMP1=1.0D0+HALFU
      IF(TEMP1.LE. (1.d0 +1.d-32)) GO TO 100
      HALFU=0.25D0*HALFU
      GO TO 50
  100 EPS=2.0D0*HALFU
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     REELLE INTEGRATION     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DSUMM(P,NXIETA,VAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION P(1)
      VAL=0.D0
      DO 1 IK=1,NXIETA
	VAL=VAL+P(IK)
1     CONTINUE
      RETURN
      END




**************************************************************
* funtion to escalate the grid
**************************************************************
      function fungrad(size,x)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      pi=2*dasin(1.0d0)
      fungrad=size/pi*dasin(x)
      return
      end
