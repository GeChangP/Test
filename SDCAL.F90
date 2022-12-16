MODULE SDCAL
    Contains
    SUBROUTINE SDCAL(XPI,YPI,AK,NB,ZS,ZD)
        !     Kernel Function: Wave Term
              IMPLICIT DOUBLE PRECISION (A-H,O-Y)
              IMPLICIT COMPLEX*16 (Z)
        
              PARAMETER (MX=105,NQ=101)
              DIMENSION ZS(NB),ZD(NB)
              COMMON /PAI/ PI,PI05,PI2
                COMMON /ELM/ XP(MX),YP(MX),XQ(NQ),YQ(NQ)
        
              Z0=(0.0D0,0.0D0)
              ZI=(0.0D0,1.0D0)
              DO 100 J=1,NB
              ZS(J)=Z0
              ZD(J)=Z0
          100 CONTINUE
        
              XX=XPI-XQ(1)
              YY=YPI+YQ(1)
              SGNX=DSIGN(1.0D0,XX)
              IF(DABS(XX).LT.1.0D-10) SGNX=0.0D0
                XE=-AK*YY
                  YE=-AK*DABS(XX)
                  ZETA=DCMPLX(XE,YE)
                  CALL EZE1Z(XE,YE,EC,ES)
              RFL1=0.5D0*DLOG(XX**2+YY**2)
              RFT1=DATAN2(YY,XX)
              ZFC1= EC-PI*CDEXP(ZETA)*ZI
              ZFS1=(ES-PI*CDEXP(ZETA))*SGNX
        
              DO 200 J=1,NB
              XX=XPI-XQ(J+1)
              YY=YPI+YQ(J+1)
              SGNX=DSIGN(1.0D0,XX)
              IF(DABS(XX).LT.1.0D-10) SGNX=0.0D0
                XE=-AK*YY
                  YE=-AK*DABS(XX)
                ZETA=DCMPLX(XE,YE)
                  CALL EZE1Z(XE,YE,EC,ES)
              RFL2=0.5D0*DLOG(XX**2+YY**2)
              RFT2=DATAN2(YY,XX)
              ZFC2= EC-PI*CDEXP(ZETA)*ZI
              ZFS2=(ES-PI*CDEXP(ZETA))*SGNX
        
              DX=XQ(J+1)-XQ(J)
              DY=YQ(J+1)-YQ(J)
              D =DSQRT(DX*DX+DY*DY)
              CDEL=DX/D
              SDEL=DY/D
              SUB =SDEL*(RFL2-RFL1)+CDEL*(RFT2-RFT1)
              ZSUB=SDEL*(ZFC2-ZFC1)-CDEL*(ZFS2-ZFS1)
              ZS(J)=ZS(J)+2.0D0/AK*(SUB+ZSUB)
              ZD(J)=ZD(J)-2.0D0*(ZFS2-ZFS1)
              RFL1=RFL2
              RFT1=RFT2
              ZFC1=ZFC2
              ZFS1=ZFS2
          200 CONTINUE
              RETURN
              END


END MODULE