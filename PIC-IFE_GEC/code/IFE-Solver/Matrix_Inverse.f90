SUBROUTINE	Matrix_Inverse(A,N)
	REAL(8), DIMENSION(N,N), INTENT(INOUT)		:: A
	INTEGER, INTENT(IN)							:: N
	INTEGER IP(N)      !记录主列号
	REAL P          !工作单元，放主元
	INTEGER  I0,R      !工作单元，放主列号
    EPS=0.01
	!write(*,*)'gjcp  ok0000000'
       DO K=1,N
          P=0
          I0=K
          IP(K)=K
          
          DO I=K,N
             IF(ABS(A(I,K)).GT.ABS(P))THEN
                P=A(I,K)
                I0=I
                IP(K)=I
             ENDIF
          ENDDO
          
!          IF(ABS(P).LE.EPS)THEN
!             WRITE(*,*)'DET=0'
!             stop !后来加的
!             GOTO 10
!          ENDIF
          
          IF(I0.NE.K)THEN
             DO J=1,N
                S=A(K,J)
                A(K,J)=A(I0,J)
                A(I0,J)=S
             ENDDO
          ENDIF
         
          A(K,K)=1./P
          
          DO I=1,N
             IF(I.NE.K)THEN
                A(I,K)=-A(I,K)*A(K,K)
                DO J=1,N
                   IF(J.NE.K)THEN
                      A(I,J)=A(I,J)+A(I,K)*A(K,J)
                     ENDIF
                ENDDO
             ENDIF
          ENDDO
         
          DO J=1,N
             IF(J.NE.K)THEN
                A(K,J)=A(K,K)*A(K,J)
             ENDIF
          ENDDO
       ENDDO
   ! write(*,*)'gjcp  ok01'    
       DO K=N-1,1,-1
          R=IP(K)
          IF(R.NE.K)THEN
             DO I=1,N
                S=A(I,R)
                A(I,R)=A(I,K)
                A(I,K)=S
             ENDDO
           ENDIF
        ENDDO

10	END SUBROUTINE
