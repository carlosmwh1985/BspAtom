      SUBROUTINE Ylm_All(l_max,th,ph,Ylm)

!      CALCULATES SPHERICAL HARMONICS FOR L= 0, 1, ... L_MAX, FOR A GIVEN th, ph
!      SAVES IT ON ARRAY Ylm
!      WRITTEN BY CARLOS GRANADOS. UNIVERSITÉ DE LORRAINE. 2013

      USE MOD_TYPES
      USE MOD_FAC_VARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l_max
      REAL(DP), INTENT(IN) :: th, ph

      COMPLEX(DPC), INTENT(OUT) :: Ylm(0:l_max,-l_max:l_max)

      INTEGER :: i, l, m, mi
      REAL(DP) :: x, x2sq, Alm, factor, Anum, Aden
      REAL(DP) :: Plm(0:l_max,-l_max:l_max)

      Ylm = DCMPLX(0.D0,0.D0)
      Plm = 0.D0

      x = DCOS(th)

!    CALCULATES ASSOC. LEGENDRE POL. Pnm(x), WITH m=0
!    USES: (2n+1) x Pnm = (n+m)Pn-1,m + (n-m+1) Pn+1,m

      Plm(0,0) = 1.D0
iflm: IF( l_max >= 1 ) THEN
        Plm(1,0) = x
domeq0: DO i = 2, l_max
          factor = DBLE(2*i-1) * x * Plm(i-1,0) - DBLE(i-1) * Plm(i-2,0)
          Plm(i,0) = factor / DBLE(i)
        END DO domeq0

!       VALUES FOR P1,1 ^ P1,-1
        x2sq = DSIN(th)
        Plm(1,1) = x2sq
        Plm(1,-1) = -Plm(1,1) / 2.D0

!        USES: (2n+1) (1-x**2)**1/2 Pnm = Pn+1,m+1 - Pn-1,m+1
!              Pn,-m = (-1)**m (n-m)!/(n+m)! Pnm
        DO l = 2, l_max
  domgt0: DO m = 1, l
            factor = fac(l-m) / fac(l+m)
            Plm(l,m) = DBLE(2*l-1)*x2sq*Plm(l-1,m-1) + Plm(l-2,m)
              Plm(l,-m) = ((-1.D0)**m) * factor * Plm(l,m)
          END DO domgt0
        END DO

      END IF iflm

      DO l = 0, l_max
        DO m = -l, l
          Anum = DBLE(2*l+1) * fac(l-m)
          Aden = 4.D0 * PI * fac(l+m)
          Alm = ((-1.D0)**m) * DSQRT( Anum / Aden )
          Ylm(l,m) = Alm * Plm(l,m) * EXP(zi*m*ph)
!         Ylm(l,m) = DCMPLX(Alm*Plm(l,m)*DCOS(m*ph),Alm*Plm(l,m)*DSIN(m*ph))
        END DO
      END DO

      RETURN

      END SUBROUTINE Ylm_All
