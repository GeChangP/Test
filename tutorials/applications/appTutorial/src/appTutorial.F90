PROGRAM MAIN
      
      
    Use OFFSET
    Use SDSUB
    Use SDCAL
    Use SOLVE
    Use FORCE
    Use MOTION
    Use ZSWEEP
    Use EZE1Z
  
  
  IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)
  COMMON /PAI/ PI,PI05,PI2
  COMMON /FILEOUT/ IRAO,IDIF,IRAD
  Integer :: nTPlus
  Double Precision :: kStart, kEnd, dk, akb
  Integer :: nK, iK

  PI     = 3.14159265358979D0
  PI05   = PI*0.5D0
  PI2    = PI*2.0D0
  IPRINT = 1
  NPRINT = 0

  write(6,*) "INPUT NB,H0,SIG1,SIG2: "

  NB = 80
  H0 = 1.D0
  SIG1 = 0.8D0
  SIG2 = 0.8D0

  WRITE(6,600) NB,H0,SIG1,SIG2

!     Center of Gravity
  OGD =0.05D0

!     Radius of Roll Gyration
  KZZB=0.35D0

  write(*,*) "NT plus value:"
  NTPlus = 3
  NT     = NB + nTPlus

  !NT=NB+3
  CALL OFFSET(NB,NT,H0,SIG1,SIG2,OGD,KZZB,NPRINT)

  WRITE(6,*) "KB Start, KB End, nK: "

  kStart = 0.01D0
  kEnd   = 5.0D0
  nK     = 300

  open(newunit=IRAO,file="MotionRAO.dat",status="replace")
  write(IRAO,'(a)') "# kB/2 |X2|/A |X3|/A |X4|/(kA)"&
  "Phs(X2) Phs(X3) Phs(X4)"

  open(newunit=IDIF,file="WaveExtForce.dat",status="replace")
  write(IRAO,'(a)') "# kB/2 |F2|/A |F3|/A |F4|/(kA)"&
  "Phs(F2) Phs(F3) Phs(F4)"

  open(newunit=IRAD,file="RadiationForce.dat",status="replace")
  write(IRAD,'(a)') "# nFreq "
  write(IRAD,'(a)') "# kB/2 "
  write(IRAD,'(a)') "# AddedMass(3, 3) "
  write(IRAD,'(a)') "# Damping(3, 3) "
  write(IRAD,'(i5)') nK

  dk = (kEnd - kStart) / ( nK - 1.D0 )
  akb = kStart
  do iK = 1, nK
      CALL SOLVE (NB,NT,AKB)
      CALL FORCE (NB,AKB,IPRINT)
      CALL MOTION(AKB,IPRINT)
      akb = akb + dk
  end do

9 STOP

600 FORMAT(//14X,48('*')&
    /19X,'2-D RADIATION AND DIFFRACTION PROBLEMS',&
    /19X,'    OF A GENERAL-SHAPED 2-D BODY',&
    /19X,'     BY INTEGRAL-EQUATION METHOD',/14X,48('*'),&
    /15X,'NUMBER OF PANELS OVER THE WHOLE BODY---(NB)=',I4,&
    /15X,'HALF-BEAM TO DRAFT RATIO---------H0(=B/2/D)=',F8.4,&
    /15X,'SECTIONAL AREA RATIO FOR RIGHT-SIG1(=S/B/D)=',F8.4,&
    /15X,'SECTIONAL AREA RATIO FOR LEFT--SIG2(=S/B/D)=',F8.4/)
  END