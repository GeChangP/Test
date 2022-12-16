MODULE SDSUB
    Contains

    SUBROUTINE SDSUB(XPI,YPI,NB,SS,DD)
        !     Kernel Function: Rankine Source
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        
              PARAMETER (MX=105,NQ=101)
              DIMENSION SS(NB),DD(NB)
              COMMON /ELM/ XP(MX),YP(MX),XQ(NQ),YQ(NQ)
        
              DO 100 J=1,NB
              SWA=0.0D0
              DWA=0.0D0
              IF(DABS(YPI).LT.1.0D-8) GOTO 10
                DX=XQ(J+1)-XQ(J)
                   DY=YQ(J+1)-YQ(J)
                  D=DSQRT(DX*DX+DY*DY)
              CDEL=DX/D
              SDEL=DY/D
              XA=XPI-XQ(J  )
              XB=XPI-XQ(J+1)
        
              SL=-1.0D0
              DO 200 L=1,2
              SL=-SL
              YA=SL*YPI-YQ(J  )
              YB=SL*YPI-YQ(J+1)
              SUBA=XA*CDEL+YA*SDEL
              SUBB=XB*CDEL+YB*SDEL
              COEF=XA*SDEL-YA*CDEL
              ABSC=DABS(COEF)
                WA1=0.5D0*(SUBB*DLOG(XB*XB+YB*YB)-SUBA*DLOG(XA*XA+YA*YA))
                IF(ABSC.LT.1.0D-10) THEN
                  WA2=0.0D0
                  WA3=0.0D0
                ELSE
              WA2=ABSC*(DATAN(SUBB/ABSC)-DATAN(SUBA/ABSC))
              WA3=WA2/COEF
              ENDIF
              SWA=SWA-(WA1+WA2)*SL
              DWA=DWA+ WA3*SL
          200 CONTINUE
        
           10 SS(J)=SWA
              DD(J)=DWA
          100 CONTINUE
              RETURN
              END


End MODULE