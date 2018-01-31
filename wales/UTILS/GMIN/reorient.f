      PROGRAM REORIENT
      IMPLICIT NONE
      INTEGER NATOMS, SEED
      DOUBLE PRECISION TX, TY, TZ, X, Y, Z, RANDOM, RAN1, COSTX, SINTX, COSTY, SINTY, COSTZ, SINTZ, ANGLE
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.283185307D0
      EXTERNAL RAN1
      LOGICAL TEST
      
      INQUIRE(FILE='randata',EXIST=TEST)
      IF (TEST) THEN
         OPEN(UNIT=8,FILE='randata',STATUS='OLD')
         READ(8,*) SEED
      ELSE
         PRINT*,'Seed?'
         READ(*,*) SEED
      ENDIF
      IF (SEED.GT.0) SEED=-SEED

      OPEN(UNIT=7,FILE='newcoords',STATUS='UNKNOWN')
      OPEN(UNIT=1,FILE='coords.ref',STATUS='OLD')

      RANDOM=RAN1(SEED)
      ANGLE=RANDOM*TWOPI ! in radians
      COSTX=COS(ANGLE)
      SINTX=SIN(ANGLE)
      RANDOM=RAN1(SEED)
      ANGLE=RANDOM*TWOPI ! in radians
      COSTY=COS(ANGLE)
      SINTY=SIN(ANGLE)
      RANDOM=RAN1(SEED)
      ANGLE=RANDOM*TWOPI ! in radians
      COSTZ=COS(ANGLE)
      SINTZ=SIN(ANGLE)

      DO 
         READ(1,*,END=20) X, Y, Z
         TY= COSTX*Y+SINTX*Z
         TZ=-SINTX*Y+COSTX*Z
         Y=TY
         Z=TZ
         TX= COSTY*X+SINTY*Z
         TZ=-SINTY*X+COSTY*Z
         X=TX
         Z=TZ
         TX= COSTZ*X+SINTZ*Y
         TY=-SINTZ*X+COSTZ*Y
         X=TX
         Y=TY
         WRITE(7,10) X, Y, Z
10       FORMAT(3F20.10)
      ENDDO

20    CLOSE(1)
      CLOSE(7)

      STOP
      END

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=2.0D0*min(AM*iy,RNMX)-1.0D0
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 1(-V%'2150)-3.
