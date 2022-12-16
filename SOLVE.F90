MODULE SOLVE
    Contains

    SUBROUTINE SOLVE(NB,NT,AK)
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
    
        PARAMETER (MX=105,NP=100,NQ=101,NEQ=4,SML=1.0D-14)
        DIMENSION ZSA(MX,NP),ZSB(MX,NEQ),ZAA(NP,NP),ZBB(NP,NEQ)
        DIMENSION ZS(NP),ZD(NP),SS(NP),DD(NP)
    
        COMMON /PAI/ PI,PI05,PI2
          COMMON /ELM/ XP(MX),YP(MX),XQ(NQ),YQ(NQ)
        COMMON /VN2/ VN(3,NP)
        COMMON /FAI/ ZFI(4,NP)
    
        Z0=(0.0D0,0.0D0)
        ZI=(0.0D0,1.0D0)
        DO 10 I=1,NB
        DO 20 J=1,NB
     20 ZAA(I,J)=Z0
        DO 10 M=1,NEQ
        ZBB(I,M)=Z0
     10 CONTINUE
    
        DO 30 I=1,NT
        DO 40 J=1,NB
     40 ZSA(I,J)=Z0
        DO 50 M=1,NEQ
     50 ZSB(I,M)=Z0
        IF(I.LE.NB) ZSA(I,I)=DCMPLX(PI,0.0D0)
     30 CONTINUE
    
        DO 100 I=1,NT
        CALL SDSUB(XP(I),YP(I),NB,SS,DD)
        CALL SDCAL(XP(I),YP(I),AK,NB,ZS,ZD)
    
        DO 110 J=1,NB
        ZSA(I,J)=ZSA(I,J)+DD(J)+ZD(J)
    110 CONTINUE
    
        DO 120 M=1,3
        DO 120 J=1,NB
        ZSB(I,M)=ZSB(I,M)+(SS(J)+ZS(J))*VN(M,J)
    120 CONTINUE
        ZSB(I,4)=PI2*CDEXP(-AK*(YP(I)-ZI*XP(I)))
    100 CONTINUE
    
        DO 200 I=1,NB
        DO 210 J=1,NB
        DO 210 K=1,NT
        ZAA(I,J)=ZAA(I,J)+ZSA(K,I)*ZSA(K,J)
    210 CONTINUE
        DO 220 M=1,NEQ
        DO 220 K=1,NT
        ZBB(I,M)=ZBB(I,M)+ZSA(K,I)*ZSB(K,M)
    220 CONTINUE
    200 CONTINUE
    
        CALL ZSWEEP(NP,NB,ZAA,ZBB,NEQ,SML)
        IF(CDABS(ZAA(1,1)).LT.SML) WRITE(6,600)
    600 FORMAT(//10X,'*** ERROR: ZSWEEP IN SUBROUTINE (SOLVE)',
       &   ' WAS ABNORMALLY DONE.',/23X,'PLEASE CHECK!'///)
    
        DO 250 M=1,NEQ
        DO 250 I=1,NB
        ZFI(M,I)=ZBB(I,M)
    250 CONTINUE
        RETURN
        END


End MODULE