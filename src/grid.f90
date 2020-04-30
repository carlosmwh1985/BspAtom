      SUBROUTINE GRID

      USE MOD_TYPES
      USE MOD_BSPLINES
      USE MOD_GRID

      IMPLICIT NONE

      INTEGER :: i, j
      REAL(DP) :: delta, hin, A1, A2, dr

!      KNOT POINTS GRID

      ALLOCATE( rt(nkp), rtk(nointv+3) )

      DO i = 1, nbc1
        rt(i) = ra
      END DO
      DO i = nkp-nbc2+1, nkp
        rt(i) = rb
      END DO

      IF( KIND_GRID == 0 ) THEN            !LINEAR GRID

        WRITE(6,'(A23)') 'Linear Knotpts Sequence'

 dolin: DO i = nbc1+1, nkp-nbc2
          rt(i) = ra + DBLE(i-nbc1) * gsize / DBLE(nointv)
        END DO dolin

      ELSE IF( KIND_GRID == 1 ) THEN      !EXPONENTIAL GRID

        WRITE(6,'(A29)') 'Exponential Knotpts Sequence'

        delta = 0.01D0
        hin = LOG(gsize/delta) / DBLE(nointv-1)
        j = 1
        rt(nbc1+1) = delta
 doexp: DO i = nbc1+2, nkp-nbc2
          rt(i) = rt(nbc1+1) * EXP(hin*j)
          j = j + 1
        END DO doexp

      ELSE IF( KIND_GRID == 2 ) THEN      !EXPONENTIAL-LINEAR

        WRITE(6,'(/,A35)') 'Exponential-Linear Knotpts Sequence'
        WRITE(6,'(A30,G12.5)') 'Limit of Exp. Sequence: rmax =', rmax

        delta = 0.01D0
        hin = LOG((rmax-ra)/delta) / DBLE(nintv_exp-1)
        j = 1
        rt(nbc1+1) = delta
        DO i = 2, nintv_exp
          rt(i+nbc1) = delta * EXP(hin*j)
          j = j + 1
        END DO

        dr = (rb - rmax) / DBLE(nintv_lin)
        DO i = nintv_exp+1, nointv
          rt(i+nbc1) = rmax + DBLE(i-nintv_exp) * dr
        END DO

      END IF

      WRITE(6,'(/,A22,I6)') 'Number of Knot Points:' , nkp
      WRITE(6,'(A27,2I3)') 'Multiplicity of END points:', nbc1, nbc2
!      WRITE(6,'(/,A5)') 'Grid:'
!      WRITE(6,100) (rt(i),i=1,nkp)

      DO i = 1, nkp-nbc1-nbc2+2
        rtk(i) = rt(i+nbc1-1)
      END  DO
!      WRITE(6,'(/,A14)') 'Physical grid:'
!      WRITE(6,100) (rtk(i), i=1,nkp-nbc1-nbc2+2)
!      WRITE(6,*)

100   FORMAT(10G12.5)

      ALLOCATE( Aind(nfun,2) )

      Aind = 0.D0            !CREATE THE COEFFICIENTS USED FOR CALCULATE THE DERIVATIVES OF B-SPLINES
      DO i = 1, nfun
        A1 = 0.D0
        A2 = 0.D0
        dr = rt(i+k-1) - rt(i)
        IF ( dr > 0.D0 ) A1 = 1.D0 / dr
        dr = rt(i+k) - rt(i+1)
        IF( dr > 0.D0 ) A2 = 1.D0 / dr
        Aind(i,1) = A1
        Aind(i,2) = A2
      END DO

!     GAUSS-LEGENDRE QUADRATURE

      ALLOCATE( xg(ka), wg(ka) )

      CALL gauleg(-1.D0,1.D0,xg,wg,ka)

      END SUBROUTINE GRID
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE SEL_LM

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_PHOTOION

      IMPLICIT NONE

      INTEGER :: il, jl, im, la, l0, m0, lf, mf, mf1, mf2, lmin

      WRITE(6,'(/,A22)') 'Selecting Final States'

      l0 = l_ini
      m0 = m_ini
      
      mf = m0 - mph
      
ifAv: IF( KIND_PI <= 2 ) THEN
      
  
        IF( KIND_PI == 0 ) THEN            !ELECTRONIC STRUCTURE ONLY
          nlm = 1
          ALLOCATE( lmf(0:nlm,2) )
          lmf(1,1) = l0
          lmf(1,2) = m0
        ELSE                              !DIPOLAR CASE: Len OR Vel GAUGES
          nlm = 2
          lf = l0 - 1
          IF( lf < 0 .OR. lf < m0 ) nlm = 1
          ALLOCATE( lmf(0:nlm,2) )
          lmf(:,2) = m0
          il = 0
          DO lf = MIN(l0-1,l0+1), l0+1, 2
            IF( lf >= 0 .AND. lf >= m0 ) THEN
              il = il + 1
              lmf(il,1) = lf
            END IF
          END DO
        END IF
        
        l_fin = lmf(nlm,1)
        
      ELSE IF( KIND_PI == 5 .OR. KIND_PI == 6 ) THEN

        nlm = 0
        DO il = 0, lmax
          IF( il >= ABS(m0) ) nlm = nlm + 1
        END DO

        ALLOCATE( lmf(0:nlm,2) )
        jl = 0
        DO il = 0, lmax
          IF( il >= ABS(m0) ) THEN
            lmf(jl+1,1) = il
            lmf(jl+1,2) = m0
            jl = jl + 1
          END IF
        END DO

      ELSE IF( KIND_PI == 8 .OR. KIND_PI == 9 ) THEN
      
        IF( KIND_NLM == 0 ) THEN          !KEEP SPECIFIED QUANTUM NUMBERS
      
          nlm = 0
          DO il = 0, lmax
!           DO im = -il, il
!             IF( il >= ABS(m0) .AND. ABS(im) == ABS(m0) ) nlm = nlm + 1
!           END DO
            IF( il >= ABS(m0) ) nlm = nlm + 1
          END DO

          ALLOCATE( lmf(0:nlm,2) )
          jl = 0
          DO il = 0, lmax
            IF( il >= ABS(m0) ) THEN
              lmf(jl+1,1) = il
              lmf(jl+1,2) = m0
              jl = jl + 1
            END IF
!           DO im = -il, il
!             IF( il >= ABS(m0) .AND. ABS(im) == ABS(m0) ) THEN
!               lmf(jl+1,1) = il
!               lmf(jl+1,2) = im
!               jl = jl + 1
!             END IF
!           END DO
          END DO

        ELSE IF( KIND_NLM == 1 ) THEN     !UNPOLARIZED INITIAL STATE (INITIAL l)

          nlm = (l0 + 1)**2               !1st PART: l <= l0 (ALL)
          nlm = nlm + (lmax - l0) * (2*l0 + 1)    !2nd PART: l > l0
          
          ALLOCATE( lmf(0:nlm,2) )
          jl = 0
          DO il = 0, lmax
            la = il
            IF( la > l0 ) la = l0
            DO im = -la, la
              jl = jl + 1
              lmf(jl,1) = il
              lmf(jl,2) = im
            END DO
          END DO
          
        END IF

      ELSE ifAv
      
          nlm = (lmax+1)**2
          ALLOCATE( lmf(0:nlm,2) )
          lmf = 0
          jl = 0
          DO il = 0, lmax
            DO im = -il, il
              jl = jl + 1
              lmf(jl,1) = il
              lmf(jl,2) = im
            END DO
          END DO        
        
!        END IF ifb0
        
      END IF ifAv
      
      lmf(0,1) = l0
      lmf(0,2) = m0
      
      WRITE(6,'(/,A22)') 'Selected final states:'
      WRITE(6,'(T3,A1,T7,A2,T12,A2)') 'i', 'lf', 'mf'
      WRITE(6,'(T2,A12)') '------------'
      DO il = 1, nlm
        WRITE(6,'(T1,I3,T6,I3,T11,I3)') il, lmf(il,1), lmf(il,2)
      END DO

      END SUBROUTINE SEL_LM
