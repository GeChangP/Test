MODULE EZEiZ
    Contains
    SUBROUTINE EZE1Z(XX,YY,EC,ES)
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        DOUBLE PRECISION  NEW
    
        DATA PI,GAMMA/3.14159265358979D0,0.5772156649015D0/
    
        X =XX
        Y =DABS(YY)
        R =DSQRT(X*X+Y*Y)
        C =DATAN2(Y,X)
    
        IF(R.GT.25.0D0)  GO TO 30
        IF(X.GT.0.0D0.AND.R.GT.8.0D0)  GO TO 20
        IF(X.LE.0.0D0.AND.Y.GT.10.0D0) GO TO 20
    
        ER=-GAMMA-DLOG(R)+R*DCOS(C)
        EI=-C+R*DSIN(C)
        SB=-R
          DO 100 N=2,100
          FN=DFLOAT(N)
          CN=C*FN
          SB=-SB*R*(FN-1.0D0)/FN/FN
          ER=ER-SB*DCOS(CN)
          EI=EI-SB*DSIN(CN)
          IF(N.EQ.100)  GO TO 1
          IF(EI.EQ.0.0D0)  GO TO 10
          IF(DABS(SB/EI).LE.1.0D-8) GO TO 10
            GO TO 100
     10   IF(DABS(SB/ER).LE.1.0D-8) GO TO 1
    100   CONTINUE
      1 CC=DEXP(X)*DCOS(Y)
        SS=DEXP(X)*DSIN(Y)
        EC=CC*ER-SS*EI
        ES=CC*EI+SS*ER
        IF(YY.LT.0.0D0) ES=-ES
        RETURN
    
     20 Z =DCMPLX(X,Y)
        Z1=(1.0D0,0.0D0)
        ZSUB=(10.0D0,0.0D0)
        ZS  =Z+ZSUB/(Z1+ZSUB/Z)
          DO 200 J=1,9
            ZSUB=DCMPLX(DFLOAT(10-J),0.0D0)
          ZS  =Z+ZSUB/(Z1+ZSUB/ZS)
    200   CONTINUE
        ZSUB=Z1/ZS
        EC=DREAL(ZSUB)
        ES=DIMAG(ZSUB)
        IF(YY.LT.0.0D0) ES=-ES
        RETURN
    
     30 OLD=-1.0D0/R
        EXC=OLD*DCOS(C)
        EXS=OLD*DSIN(C)
          DO 300 N=2,100
            NEW=-OLD/R*DFLOAT(N-1)
            IF(EXS.EQ.0.0D0) GO TO 31
            IF(DABS(NEW/EXS).LE.1.0D-8) GO TO 31
            GO TO 32
     31   IF(EXC.EQ.0.0D0) GO TO 32
          IF(DABS(NEW/EXC).LE.1.0D-8) GO TO 33
     32   IF(DABS(OLD).LT.DABS(NEW))  GO TO 33
          OLD=NEW
            EXC=EXC+OLD*DCOS(C*DFLOAT(N))
          EXS=EXS+OLD*DSIN(C*DFLOAT(N))
    300   CONTINUE
     33 EC=-EXC
        ES=EXS
        IF(DABS(PI-DABS(C)).LT.1.0D-10) ES=-PI*DEXP(X)
        IF(YY.LT.0.0D0) ES=-ES
        RETURN
        END
    
End MODULE