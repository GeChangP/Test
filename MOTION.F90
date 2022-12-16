MODULE MOTION
    Contains
    SUBROUTINE MOTION(AK,IPRINT)
        IMPLICIT DOUBLE PRECISION (A-H,K,O-Y)
        IMPLICIT COMPLEX*16 (Z)
        DIMENSION ZAA(3,3),ZBB(3)
        DIMENSION AMPG(3),PHAG(3),ZMTNG(3)
        COMMON /PAI/ PI,PI05,PI2
        COMMON /FILEOUT/ IRAO,IDIF,IRAD
        COMMON /MDT/ CMAS,C22,OG,KZZ,GM
        COMMON /FCE/ ZAB(3,3),ZEXF(3)
        COMMON /MTN/ ZMTNO(3)
    
        SML=1.0D-14
    
        ZAA(1,1)=-AK*(CMAS+ZAB(1,1))
        ZAA(1,2)=-AK* ZAB(1,2)
        ZAA(1,3)=-AK*(ZAB(1,3)+OG*ZAB(1,1))
        ZBB(1  )= ZEXF(1)
    
        ZAA(2,1)=-AK*ZAB(2,1)
        ZAA(2,2)=-AK*(CMAS+ZAB(2,2))+C22
        ZAA(2,3)=-AK*(ZAB(2,3)+OG*ZAB(2,1))
        ZBB(2  )= ZEXF(2)
    
        ZAA(3,1)=-AK*(ZAB(3,1)+OG*ZAB(1,1))
        ZAA(3,2)=-AK*(ZAB(3,2)+OG*ZAB(1,2))
        ZAA(3,3)=-AK*(CMAS*KZZ**2+ZAB(3,3)+OG*ZAB(1,3)
       &         +OG*(ZAB(3,1)+OG*ZAB(1,1)))+CMAS*GM
        ZBB(3  )= ZEXF(3)+OG*ZEXF(1)
    
        CALL ZSWEEP(3,3,ZAA,ZBB,1,SML)
        IF(CDABS(ZAA(1,1)).LT.SML) WRITE(6,600)
    600 FORMAT(///10X,'+++ ERROR: ZSWEEP IN (MOTION) +++'///)
        DO 100 I=1,3
        ZMTNG(I)=ZBB(I)
    100 CONTINUE
        ZMTNO(1)=ZMTNG(1)+OG*ZMTNG(3)
        ZMTNO(2)=ZMTNG(2)
        ZMTNO(3)=ZMTNG(3)
    
        DO 200 I=1,3
        AMPG(I)=CDABS(ZMTNG(I))
        IF(I.EQ.3) AMPG(I)=AMPG(I)/AK
        PHAG(I)=DATAN2(DIMAG(ZMTNG(I)),DREAL(ZMTNG(I)))*180.0D0/PI
    200 CONTINUE
    
        write(IRAO, 620) AK, ( AMPG(I), I=1, 3), ( PHAG(I), I=1, 3)
    
        IF(IPRINT.EQ.0) RETURN
        WRITE( 6,610) AK,(AMPG(I),PHAG(I),I=1,3)
    610 FORMAT(//5X,'+++++ MOTIONS ABOUT ''G'' FOR K*B/2=',F7.3,
       &   '+++++',/21X,'AMP.',7X,'PHASE',/9X,'SWAY  ',E11.4,
       &   2X,F9.3,' (DEG)' ,/9X, 'HEAVE ',E11.4,2X,F9.3,' (OEG)',
       &   /9X, 'ROLL  ',E11.4,2X,F9.3,' (DEG)')
        RETURN
    620 FORMAT(99(1pe15.6))
        END
    


END MODULE