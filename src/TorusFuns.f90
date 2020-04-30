      SUBROUTINE TORMAT

!     SUBROUTINE TO CALCULATE THE NEEDED FUNCTIONS TO OBTAIN THE TOROIDAL MOMENT...
!       1st: CALCULATE NEEDED ANGULAR INTEGRALS AND THE CORRESPONDING TERMS
!       2nd: CALCULATE RADIAL WAVE FUNCTIONS AND THE CORRESPONDING TERMS


      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_MATRICES, ONLY: Xij, rvecij, zrangij, zPmq !, fJrij
      USE MOD_PHOTOION, ONLY: nlm, lmf, n1_max, cinl, DSVMV

      IMPLICIT NONE
      
      INTEGER :: i, il, ni, li, mi, nj, jl, lj, mj, l1, m1, l2, m2
      INTEGER :: lf, mq, ith, ir, nr, igl, left
      REAL(DP) :: c1, c2, r, fvmv
      COMPLEX(DPC) :: zc0

      REAL(DP), ALLOCATABLE :: xgl(:), wgl(:)
      REAL(DP), ALLOCATABLE :: x(:), y(:), v_temp(:)
      
      WRITE(6,'(/,A45)') 'Calculating Related Terms for Toroidal Moment'

      lf = lmax
      nr = nointv * ka
      
!     ANGULAR TERMS

      ALLOCATE( zPmq(nlm,nlm,-1:1,0:nfib0,4) )
      ALLOCATE( zrangij(0:lf,-lf:lf,0:lf,-lf:lf,3) )
!     ALLOCATE( zj0angij(nlm,nlm,-1:1,2) )
      zPmq = 0.D0
      zrangij = 0.D0
!     zj0angij = 0.D0

      CALL PIFUNS(lf,zPmq)

!      WRITE DOWN ANGULAR TERMS... ANGULAR INTEGRALS MUST BE PERFORMED LATER,
!        ONCE THE TIME PROPAGATION IS FINISHED

      OPEN( UNIT=80, FILE='CSs/CurrentMatFunAng.dat', ACTION='WRITE' )
      OPEN( UNIT=85, FILE='CSs/rMatElemAng.dat', ACTION='WRITE' )

      WRITE(6,'(A26)') 'Writing Down Angular Terms'

dobra:DO il = 1, nlm
      
        li = lmf(il,1)
        mi = lmf(il,2)
      
 doket: DO jl = 1, nlm
      
          lj = lmf(jl,1)
          mj = lmf(jl,2)

    domq: DO mq = -1, 1
      doth: DO ith = 0, nfib0
                    
              WRITE(80,100) il, jl, mq, ith, zPmq(il,jl,mq,ith,1),      &
     &                      zPmq(il,jl,mq,ith,2), zPmq(il,jl,mq,ith,3), &
     &                      zPmq(il,jl,mq,ith,4)
     
            END DO doth

!           WRITE(85,110) il, jl, mq, zj0angij(il,jl,mq,1), zj0angij(il,jl,mq,2)

          END DO domq

          WRITE(85,110) li, mi, lj, mj, (zrangij(li,mi,lj,mj,i), i=1,3)
          
        END DO doket
      END DO dobra

100   FORMAT(3I4,I6,8G20.10)
110   FORMAT(4I4,6G20.10)

      CLOSE(80)
      CLOSE(85)

      DEALLOCATE( zPmq, zrangij )

!     RADIAL TERMS
      WRITE(6,'(/,A28)') 'Calculating Radial Functions'

      nr = nointv * ka

      ALLOCATE( rtot(nr), xgl(ka), wgl(ka) )

      rtot = 0.D0

      xgl = 0.D0
      wgl = 0.D0

!     GAUSS-LEGENDRE QUADRATURE
      i = 0
      DO left = 1, nointv
        CALL gauleg(rtk(left),rtk(left+1),xgl,wgl,ka)
        DO igl = 1, ka
          i = i + 1
          rtot(i) = xgl(igl)
        END DO
      END DO

      DEALLOCATE( xgl, wgl )

!     ALLOCATE( fJrij(n1_max,0:lmax,n1_max,0:lmax,1,2) )
!     fJrij = 0.D0

!     WRITE DOWN RADIAL TERMS... ANGULAR INTEGRALS MUST BE PERFORMED LATER,
!       ONCE THE TIME PROPAGATION IS FINISHED
!     OPEN( UNIT=90, FILE='CSs/CurrentMatFunRad.dat', ACTION='WRITE' )
      OPEN( UNIT=95, FILE='CSs/rMatElemRad.dat', ACTION='WRITE' )

!     DO ir = 1, nr
!      
!       r = rtot(ir)
!       CALL FRMATINT(ir,r)
!        
!     END DO

!     DEALLOCATE( fJrij )

!     MATRIX ELEMENTS OF r

      ALLOCATE( rvecij(n1_max,0:lmax,n1_max,0:lmax) )
      rvecij = 0.D0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n1_max,lmax,nfun,cinl,Xij,rvecij)

      ALLOCATE( x(nfun), y(nfun), v_temp(nfun) )

!$OMP DO

      DO ni = 1, n1_max
        DO li = 0, lmax

          DO nj = 1, n1_max
            DO lj = 0, lmax
      
              x(1:nfun) = cinl(1:nfun,ni,li)
              y(1:nfun) = cinl(1:nfun,nj,lj)
              CALL DSVMV('U',nfun,x,y,v_temp,Xij,fvmv)

              rvecij(ni,li,nj,lj) = fvmv
              
            END DO
          END DO
          
        END DO
      END DO

!$OMP END DO

      DEALLOCATE( x, y, v_temp )

!$OMP END PARALLEL

      DO ni = 1, n1_max
        DO li = 0, lmax

          DO nj = 1, n1_max
            DO lj = 0, lmax
            
!              DO ir = 1, nointv*ka

!                WRITE(90,150) ni, li, nj, lj, ir, &
!     &          fJrij(ni,li,nj,lj,ir,1), fJrij(ni,li,nj,lj,ir,2)

!              END DO

              WRITE(95,160) ni, li, nj, lj, rvecij(ni,li,nj,lj)

              END DO
          END DO
      
        END DO
      END DO

150   FORMAT(4I4,I6,2G20.10)
160   FORMAT(4I4,G20.10)

      CLOSE(90)
      CLOSE(95)

      DEALLOCATE( rvecij )

      END SUBROUTINE TORMAT
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE WFALL(irmin,irmax,ni_max,cinl,nmax) !,fur,dfur)

!      SUBROUTINE TO CALCULATE ALL REDUCED RADIAL WAVEFUNCTIONS AND THEIR DERIVATIVES
!      C. M. GRANADOS--CASTRO
!      10.2018

      USE MOD_TYPES
      USE MOD_GRID
      USE MOD_BSPLINES
      USE MOD_MATRICES, ONLY: fur, dfur

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: irmin, irmax, ni_max, nmax
      REAL(DP), INTENT(IN) :: cinl(nfun,ni_max,0:lmax)
!      REAL(DP), INTENT(OUT) :: fur(irmin:irmax,nmax,0:lmax)
!      REAL(DP), INTENT(OUT) :: dfur(irmin:irmax,nmax,0:lmax)

      INTEGER :: i, j, ir, jp, ifun, lj, left, mflag
      REAL(DP) :: r, A1, A2, b1, b2

      REAL(DP), ALLOCATABLE :: bsp(:), dbsp(:), bsp1(:), bspp(:)

      IF( ALLOCATED(fur) ) DEALLOCATE( fur )
      IF( ALLOCATED(dfur) ) DEALLOCATE( dfur )
      ALLOCATE( fur(irmin:irmax,nmax,0:lmax), dfur(irmin:irmax,nmax,0:lmax) )
      fur = 0.D0
      dfur = 0.D0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(rt,nkp,k,ka,nfun,irmin,irmax,lmax,ni_max, &
!$OMP&  nmax,Aind,rtot,fur,dfur,cinl)

      ALLOCATE( bsp(k), dbsp(k), bsp1(k-1), bspp(k+1) )

!$OMP DO

      DO ir = irmin, irmax

        bsp = 0.D0
        bsp1 = 0.D0

        r = rtot(ir)

        CALL BSPALL(r,left,bsp,dbsp,bsp1,bspp)

        DO j = 1, nmax
          DO i = 1, nfun
            
            ifun = i - (left-k)

            IF( ifun >= 1 .AND. ifun <= k ) THEN
              DO lj = 0, lmax
                fur(ir,j,lj) = fur(ir,j,lj) + cinl(i,j,lj) * bsp(ifun)
                dfur(ir,j,lj) = dfur(ir,j,lj) + cinl(i,j,lj) * dbsp(ifun)
              END DO
            END IF
              
          END DO
        END DO

      END DO

!$OMP END DO

      DEALLOCATE( bsp, dbsp, bsp1, bspp )

!$OMP END PARALLEL

      END SUBROUTINE WFALL
