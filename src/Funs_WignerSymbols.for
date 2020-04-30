      DOUBLE PRECISION FUNCTION THREE_J(I1,I2,I3,M1,M2,M3)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION AC(140),FAC(140)                                                                                             
      SAVE  INIT,FAC
      DATA INIT/ 0/
      
      IF(INIT .EQ. 0 ) THEN
      INIT= 1                                                       
      NMAX=140
      FAC(1)=0.D00                                                      
      DO I=2,NMAX                                                   
      FAC(I)=FAC(I-1)+DLOG(DFLOAT(I-1))                                 
      ENDDO
      ENDIF      
                                                            
      IF(INIT .EQ. 0) GO TO 11                                             
      nmax=140
      iw=6
      L4=I1+I2+I3+2                                                     
      IF(L4.GT.NMAX)GO TO 1                                             
      IF(M1+M2+M3) 2,3,2                                                
  3   IZMAX=MIN0(I1+I2-I3,I1-M1,I2+M2)+1                                
      IZMIN=MAX0(0,I2-I3-M1,I1+M2-I3)+1                                 
      IF(IZMAX-IZMIN) 2,4,4                                             
  4   L1=I1+I2-I3+1                                                     
      L2=I3+I1-I2+1                                                     
      L3=I3+I2-I1+1                                                     
      L5=I1+M1+1                                                        
      L6=I1-M1+1                                                        
      L7=I2+M2+1                                                        
      L8=I2-M2+1                                                        
      L9=I3+M3+1                                                        
      L10=I3-M3+1                                                       
      ABRA=0.5D00*(FAC(L1)+FAC(L2)+FAC(L3)-FAC(L4)+FAC(L5)+             
     1FAC(L6)+FAC(L7)+FAC(L8)+FAC(L9)+FAC(L10))                         
      K1=I3-I2+M1+1                                                     
      K2=I3-I1-M2+1                                                     
      GROS=250.D00                                                      
      DO 8II=IZMIN,IZMAX                                                
      I=II-1                                                            
      K3=L1-I                                                           
      K4=L6-I                                                           
      K5=L7-I                                                           
      K6=K1+I                                                           
      K7=K2+I                                                           
      AC(II)=FAC(I+1)+FAC(K3)+FAC(K4)+FAC(K5)+FAC(K6)+FAC(K7)           
      IF(AC(II).LT.GROS) GROS=AC(II)                                   
 8    CONTINUE                                                          
      ACCU=0.D00                                                        
      SIG=(-1.D00)**IZMIN                                               
      DO 9II=IZMIN,IZMAX                                                
      SIG=-SIG                                                          
      AC(II)=AC(II)-GROS                                                
      ACCU=ACCU+SIG*DEXP(-AC(II))                                       
 9    CONTINUE                                                          
      THREE_J=(-1.D00)**(I1-I2-M3)*DEXP(ABRA-GROS)*ACCU                       
      RETURN                                                            
 11   WRITE(3,98)                                                       
 98   FORMAT(1X,'WARNING: LOGFAC NOT CALLED BEFORE LOGCLE')             
      STOP                                                              
  2   THREE_J=0.D00                                                           
      RETURN                                                            
  1   WRITE(3,99)                                                       
 99   FORMAT(1X,'LOGCLE: VALUE OF PARAMETER TOO HIGH')                  
      RETURN                                                            
      END 
!
!-------------------------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION SJSI(J1,J2,J3,L1,L2,L3)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

 	INTEGER, PARAMETER :: MFD = 80
	DOUBLE PRECISION, DIMENSION(MFD) :: FCT

      SAVE INIT, FCT
      DATA INIT/0/

      PHASE(I)=DBLE(1-2*MOD(ABS(I),2))

      IF(INIT .EQ. 0) THEN
        INIT = 1
        FCT(1) = 1.D0
        DO 10 I = 2, MFD
          FCT(I) = FCT(I-1) * DBLE(I-1)
   10   CONTINUE
      END IF

      SJSI = 0.D0
      IF (J1+J2.LT.J3.OR.ABS(J1-J2).GT.J3
     :   .OR.J1+L2.LT.L3.OR.ABS(J1-L2).GT.L3
     :     .OR.L1+J2.LT.L3.OR.ABS(L1-J2).GT.L3
     :       .OR.L1+L2.LT.J3.OR.ABS(L1-L2).GT.J3) RETURN

      IWMIN = MAX(J1+J2+J3,J1+L2+L3,L1+L2+J3,L1+J2+L3)
      IWMAX = MIN(J1+J2+L1+L2,J2+J3+L2+L3,J1+J3+L1+L3)

      IF(IWMAX+2 .GT. MFD) THEN
        WRITE(6,'(A)') ' *** Factorials too big in SJSI'
        STOP
      END IF

      OMEGA = 0.D0
      DO 20 IW = IWMIN, IWMAX
        OMEGA = OMEGA + PHASE(IW)*FCT(IW+2) /
     :    (  FCT(1+IW-J1-J2-J3) * FCT(1+IW-J1-L2-L3)
     :     * FCT(1+IW-L1-J2-L3) * FCT(1+IW-L1-L2-J3)
     :     * FCT(J1+J2+L1+L2-IW +1) * FCT(J1+J3+L1+L3-IW +1)
     :     * FCT(J2+J3+L2+L3-IW +1) )
   20 CONTINUE

      SJSI = OMEGA * DELTA(J1,J2,J3,FCT) * DELTA(J1,L2,L3,FCT) *
     :               DELTA(L1,J2,L3,FCT) * DELTA(L1,L2,J3,FCT)

      RETURN

      END FUNCTION SJSI
!
!-------------------------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION DELTA(J1,J2,J3,FCT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

	INTEGER, PARAMETER :: MFD = 80
	DOUBLE PRECISION, DIMENSION(MFD) :: FCT

      DELTA = SQRT( FCT(J1+J2-J3+1)*FCT(J1-J2+J3+1)*FCT(-J1+J2+J3+1) /
     &              FCT(J1+J2+J3+2) )

      RETURN

      END
