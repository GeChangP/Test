MODULE ZSWEEP
    Contains
    SUBROUTINE ZSWEEP(NDIM,N,ZA,ZB,NEQ,EPS)
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
    
        DIMENSION ZA(NDIM,NDIM), ZB(NDIM,NEQ)
        DO 5 K=1,N
        P=0.0D0
        DO 1 I=K,N
        IF(P.GE.CDABS(ZA(I,K))) GOTO 1
        P=CDABS(ZA(I,K))
        IP=I
      1 CONTINUE
        IF(P.LE.EPS) GOTO 6
        IF(IP.EQ.K)  GOTO 7
          DO 2 J=K,N
          ZW=ZA(K,J)
          ZA(K,J)=ZA(IP,J)
      2   ZA(IP,J)=ZW
            DO 20 J=1,NEQ
            ZW=ZB(K,J)
            ZB(K,J)=ZB(IP,J)
     20     ZB(IP,J)=ZW
      7 CONTINUE
        IF(K.EQ.N) GOTO 70
        DO 3 J=K+1,N
      3 ZA(K,J)=ZA(K,J)/ZA(K,K)
     70   DO 30 J=1,NEQ
     30   ZB(K,J)=ZB(K,J)/ZA(K,K)
        DO 5 I=1,N
        IF(I.EQ.K) GOTO 5
        IF(K.EQ.N) GOTO 40
          DO 4 J=K+1,N
      4   ZA(I,J)=ZA(I,J)-ZA(I,K)*ZA(K,J)
     40   CONTINUE
          DO 45 J=1,NEQ
     45   ZB(I,J)=ZB(I,J)-ZA(I,K)*ZB(K,J)
      5 CONTINUE
        ZA(1,1)=(1.0D0,0.0D0)
        RETURN
      6 ZA(1,1)=DCMPLX(DABS(P),0.0D0)
        RETURN
        END

End MODULE