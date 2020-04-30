subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascENDing vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to DO with the previous
!    call.  THEN, IF
!
!      XT(ILO) <= X < XT(ILO+1),
!
!    we set LEFT = ILO and are DOne after just three comparisons.
!
!    Otherwise, we repeatedly DOuble the dIFference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!
!      XT(ILO) <= X < XT(IHI)
!
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is THEN returned.
!
!  Modified:
!
!    14 February 2007
!
!	MODIFIED AGAIN ON 18.01.2018. HAVING PROBLEMS WITH PARALLEL CODES.
!	RELATED WITH THE SAVED VARIABLE ?
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real(wp) XT(LXT), a nondecreasing sequence of values.
!
!    Input, INTEGER LXT, the dimension of XT.
!
!    Input, real(wp) X, the point whose location with
!    respect to the sequence XT is to be determined.
!
!    Output, INTEGER LEFT, the index of the bracketing value:
!      1     IF             X  <  XT(1)
!      I     IF   XT(I)  <= X  < XT(I+1)
!      LXT   IF  XT(LXT) <= X
!
!    Output, INTEGER MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <= X
!
	USE MOD_TYPES

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: LXT
	REAL(DP), INTENT(IN) :: x, xt(LXT)
	
	INTEGER, INTENT(OUT) :: left, mflag

	INTEGER :: ilo

!	FIRST CHECK THAT X LIES WITHIN THE RANGE OF DATA

	IF( x > xt(LXT) ) THEN
	  mflag = 1
	  left = 1
	  RETURN
	ELSE IF( x < xt(1) ) THEN
	  mflag = -1
	  left = 1
	  RETURN
	ELSE
	  mflag = 0
	END IF

!	NOW CHECK WHERE IS X

	IF( x == xt(LXT) ) THEN
	  left = LXT
doleft: DO
	    IF( xt(left) < xt(LXT) ) RETURN
	    LEFT = LEFT - 1
	  END DO doleft
	ELSE
	  ilo = LXT - 1
 doilo: DO
	    IF( x < xt(ilo+1) .AND. x >= xt(ilo) ) THEN
	      left = ilo
	      EXIT doilo
	    ELSE
	      ilo = ilo - 1
	      IF( ilo == 0 ) EXIT doilo
	    END IF
	  END DO doilo
	END IF

	END SUBROUTINE interv
