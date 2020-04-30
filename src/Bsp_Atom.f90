      PROGRAM BSP_ATOM_PI

!      PROGRAM TO CALCULATE THE ELECTRONIC STRUCTURE OF A H ATOM
!      CONSIDERING A MOVING NUCLEI
!      BY CARLOS MARIO GRANADOS--CASTRO
!      MARTIN-LUTHER UNIVERSITAET, HALLE
!      30.01.2017

!      V.3.0.0 04.10.2017: Adapted to study PI with Gaussian Beams
!      V.3.0.1 31.10.2017: BSP_INTEGRATOR added (for integrals on theta).
!                          SUBROUTINE BSPLINES now working more as a module, calculate
!                          B-splines for any grid
!      V.3.0.2 12.12.2017: Several changes, particularly for the angular integrals
!      V.3.1.0 13.12.2017: First implementation to use Laguerre-Gaussian beams
!      V.3.2.0 14.12.2017: All numerical integrals over th and ph are calculated
!                          as complex numbers. Small improvements

!      To Compile:
!      ifort Modules.f90 Bsp_Atom.f90 ReadInputs.f90 grid.f90 Ang_Ints.f90 Ang_Ints_Aux.f90 matrices.f90 PhotoIon.f90 TorusFuns.f90 TorusFunsInts.f90 CubicSpline.f90 WriteWF.f90 Funs_AssLegendre.f90 Funs_AssLaguerre.f90 Funs_SphHarms.f90 Funs_Bessel.f90 Funs_WignerSymbols.for bsplvb.f90 interv.f90 gauleg.f -o Bsp_Atom_omp.x -qopenmp -mkl

      USE MOD_TYPES
      USE MOD_PHOTOION
      USE MOD_GRID
      USE MOD_BSPLINES

      IMPLICIT NONE

      LOGICAL :: fsind

!!!      EXTERNAL unlimit_stack

!$      INTEGER :: OMP_GET_NUM_THREADS

      WRITE(6,'(A64)') 'PROGRAM TO CALCULATE ELECTRONIC STRUCTURE AND PI CROSS SECTIONS,'
      WRITE(6,'(A17,/)') '  USING B-SPLINES'

!$OMP PARALLEL
!$OMP SINGLE
!$       WRITE(6,'(A21,I3,A8)') 'Open-MP Parallel with', OMP_GET_NUM_THREADS(), ' Threads'
!$OMP END SINGLE
!$OMP END PARALLEL


!      READ INPUT FILES
      CALL READ_INPUTS

!      DETERMINE KNOT POINTS IN ASKED SEQUENCE
      CALL GRID

!      CALCULATE BSPLINES AND THEIR FIRST DERIVATIVE
!      ALLOCATE( fbsp(nfun,5*ka,0:nfun), dfbsp(nfun,5*ka,0:nfun) )
!      fbsp = 0.D0
!      dfbsp = 0.D0
!      CALL B_SPLINES(k,nfun,5*ka,rt,Aind,fbsp,dfbsp)

!      SELECT FINAL STATES, FOLLOWING SELECTION RULES (IF APPLICABLE)
      CALL SEL_LM
      
      INQUIRE( DIRECTORY='CSs', EXIST=fsind )            !ONLY WORKS WITH ifort
      IF( .NOT. fsind ) CALL system('mkdir CSs')
      
!      CALCULATES RADIAL PROFILE OF THE VECTOR POTENTIAL IF NEEDED.
!      FOR GAUSSIAN AND LAGUERRE-GAUSSIAN BEAMS SOME ANGULAR INTEGRALS MUST BE
!      PERFORMED FIRST
      IF( KIND_PI >= 3 ) THEN
        CALL ANG_GRID                                    !CREATE ANGULAR GRID
!        CALL MAT_TH                                    !OVERLAP MATRICES FOR ANGULAR GRID
        CALL MAKE_F_ANG                                    !F(th), INTEGRAL OVER ph IF NEEDED
      END IF

!      CALCULATE OVERLAP, POTENTIAL AND KINETIC ENERGY MATRICES
      CALL MATRIX_SVT

!      SOLVE GENERALIZED EIGENVALUE PROBLEM - H ATOM FIXED NUCLEI
      CALL SOLVE_SYSTEM

!      CALCULATE DIPOLAR MATRICES AND PI CROSS SECTIONS, IF ASKED FOR
      IF( KIND_PI /= 0 ) THEN

        CALL END_PROG(1)

        CALL TRANS_AMP

        IF( KIND_PI >= 8 ) THEN
          CALL END_PROG(2)
          CALL TORMAT
        ELSE
          CALL CROSS_SECTIONS
        END IF

        CALL END_PROG(3)

      END IF

      WRITE(6,'(/,A17)') 'Program Finished!'

      END PROGRAM BSP_ATOM_PI
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE WRITE_WF(n,ci)
      
      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(DP), INTENT(IN) :: ci(n)

      INTEGER :: i, j, jmin, jmax, jfun, npts, left, mflag
      REAL(DP) :: r, dr, sumf, fr
      REAL(DP), ALLOCATABLE :: bsp(:)
      
      WRITE(6,'(/,A29,/)') 'Writing down Initial State WF'

      npts = 10000
      
      dr = (rb - ra) / DBLE(npts)
      ALLOCATE( bsp(k) )

      OPEN( UNIT=30, FILE='wf_n0.dat', ACTION='WRITE' )

      DO i = 0, npts
      
        r = ra + DBLE(i) * dr

        bsp = 0.D0  
        CALL interv(rt,nkp,r,left,mflag)
        CALL bsplvb(nkp,rt,k,1,r,left,bsp)        

        jmin = left - k + 1
        jmax = jmin + k - 1!MIN(jmin + k - 1, n)

        sumf = 0.D0
        DO j = jmin, jmax
          jfun = j - (left-k)
          fr = 0.D0
          IF( j >= 1 .AND. j <= n ) fr = ci(j)
          IF( jfun >= 1 .AND. jfun <= k ) sumf = sumf + fr * bsp(jfun)
        END DO

        WRITE(30,'(2G20.10)') r, sumf
        
      END DO

      CLOSE(30)
      
      DEALLOCATE( bsp )

      END SUBROUTINE WRITE_WF
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE END_PROG(ID_CALL)

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_MATRICES
      USE MOD_PHOTOION
      USE MOD_FAC_VARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ID_CALL

      WRITE(6,'(/,A16,I1,/)') 'Ending Call ID: ', ID_CALL

      IF( ID_CALL == 1 ) THEN                  !CALLING AFTER FINISH DIAGONALIZATION
                                          !  AND CALCULATION AF ALL MATRICES
        IF( ALLOCATED(xth) ) DEALLOCATE(xth)
        IF( ALLOCATED(wth) ) DEALLOCATE(wth)
        IF( ALLOCATED(zjt) ) DEALLOCATE(zjt)
        IF( ALLOCATED(ciLP) ) DEALLOCATE(ciLP)
!        IF( ALLOCATED(Afth) ) DEALLOCATE(Afth)
        IF( ALLOCATED(zAfth) ) DEALLOCATE(zAfth)
        IF( ALLOCATED(zdAfth) ) DEALLOCATE(zdAfth)
!        IF( ALLOCATED(zciexp) ) DEALLOCATE(zciexp)
        IF( ALLOCATED(Hij) ) DEALLOCATE(Hij)
        IF( ALLOCATED(Bij) ) DEALLOCATE(Bij)
        IF( ALLOCATED(En) ) DEALLOCATE(En)
        IF( ALLOCATED(Vij) ) DEALLOCATE(Vij)
        IF( ALLOCATED(Uij) ) DEALLOCATE(Uij)
        IF( ALLOCATED(Tij) ) DEALLOCATE(Tij)
        IF( ALLOCATED(Ith) ) DEALLOCATE(Ith)
        IF( ALLOCATED(zIth) ) DEALLOCATE(zIth)

      ELSE IF( ID_CALL == 2 ) THEN

        IF( ALLOCATED(Sij) ) DEALLOCATE(Sij)
        IF( ALLOCATED(xg) ) DEALLOCATE(xg)
        IF( ALLOCATED(wg) ) DEALLOCATE(wg)      
        IF( ALLOCATED(ci_ini) ) DEALLOCATE(ci_ini)
        IF( ALLOCATED(ci_fin) ) DEALLOCATE(ci_fin)
        IF( ALLOCATED(zT_fi) ) DEALLOCATE(zT_fi)
        
      ELSE

        IF( ALLOCATED(rt) ) DEALLOCATE(rt)
        IF( ALLOCATED(rtf) ) DEALLOCATE(rtf)
        IF( ALLOCATED(rtk) ) DEALLOCATE(rtk)
        IF( ALLOCATED(Aind) ) DEALLOCATE(Aind)
        IF( ALLOCATED(tht) ) DEALLOCATE(tht)
        IF( ALLOCATED(pht) ) DEALLOCATE(pht)
        IF( ALLOCATED(zT_Ef) ) DEALLOCATE(zT_Ef)
        IF( ALLOCATED(zS_Ef) ) DEALLOCATE(zS_Ef)
        IF( ALLOCATED(zAij) ) DEALLOCATE(zAij)
        IF( ALLOCATED(T_fi) ) DEALLOCATE(T_fi)
        IF( ALLOCATED(E_ini) ) DEALLOCATE(E_ini)
        IF( ALLOCATED(E_fin) ) DEALLOCATE(E_fin)
        IF( ALLOCATED(T2fi) ) DEALLOCATE(T2fi)
        IF( ALLOCATED(Enf) ) DEALLOCATE(Enf)
        IF( ALLOCATED(n01) ) DEALLOCATE(n01)
        IF( ALLOCATED(lmf) ) DEALLOCATE(lmf)
        IF( ALLOCATED(rij) ) DEALLOCATE(rij)
        IF( ALLOCATED(Enl) ) DEALLOCATE(Enl)
        IF( ALLOCATED(cinl) ) DEALLOCATE(cinl)
        IF( ALLOCATED(fac) ) DEALLOCATE(fac)
!        IF( ALLOCATED(fJrij) ) DEALLOCATE(fJrij)
      
      END IF
      
      END SUBROUTINE END_PROG
