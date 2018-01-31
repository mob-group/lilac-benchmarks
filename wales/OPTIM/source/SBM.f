      SUBROUTINE SBM(qo,NATOMS,grad,energy,GTEST,STEST)
      USE KEY
      USE SBMDATA 
      implicit NONE
      INTEGER NATOMS,i,j
      DOUBLE PRECISION qo(3*NATOMS), grad(3*NATOMS)
      DOUBLE PRECISION ENERGY

      LOGICAL :: CALLED=.FALSE.
      LOGICAL GTEST, STEST

        if(.NOT.CALLED)then

        write(*,*)
        write(*,*)  'Calculations will use a Structure-based SMOG model:'
        write(*,*)  '   For more information on SMOG models, see:'
        write(*,*)  '   Software: Noel, et al. PLoS Comput Biol 12, e1004794, 2016.'
        write(*,*)  '   Model: see reference list at smog-server.org/refs.html'
        write(*,*)

        call SBMinit(NATOMS)
        CALLED=.TRUE.
        SBMHESSATOM=-1
        endIF


! call the energy routine
      call calc_energy_SBM(qo,natoms,GRAD,energy)

      IF (STEST) THEN
         PRINT '(A)','ERROR - second derivatives not available'
         STOP
      ENDIF
      return
      end

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMinit
!***********************************************************************

      subroutine SBMINIT(NATOMS)
      USE KEY
      USE SBMDATA
      USE COMMONS, only: ATMASS
      USE GENRIGID, only: RIGIDINIT
      implicit NONE

        integer i,j,K,MaxCon,NATOMS,storage, dummy,ANr, Ib22,Ib21,IT1,JT1,
     Q KT1,IT2,JT2,KT2,IP1,JP1,KP1,LP1,IP2,JP2,KP2,LP2,nBA1,nTA1,nPA1, 
     Q nBA2,nTA2,nPA2,ind1,ind2,ANt,MDT1, MDT2, cl1,cl2,tempi,
     Q ATYPE,POS
      DOUBLE PRECISION :: Sigma, EpsC
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: TMPSTT
      DOUBLE PRECISION RDIFF,ALPHA
      INTEGER BONDTYPE,MAXBONDSPERATOM
      INTEGER M
      integer AA,BB,ANTEMP
      DOUBLE PRECISION dx,dy,dz
      double precision PI,ARG1,ARG2,ARG3,ARG4
      DOUBLE PRECISION RSig, Reps,DC
      integer TMPARG
      character TMPCHAR
      integer TMPINT,nexc,I1,I2
      INTEGER TT1,TT2
      double precision TMPREAL,concentration
      LOGICAL TARR(NATOMS)
      INTEGER EXCLUSIONS,NUMOFEXCLUSIONS
      INTEGER MAXEXCLUSIONS,MAXEXCLUSIONSELEC
      PARAMETER(MAXEXCLUSIONS=20)
      PARAMETER(MAXEXCLUSIONSELEC=5)
      DIMENSION EXCLUSIONS(NATOMS*MAXEXCLUSIONS),NUMOFEXCLUSIONS(NATOMS)
      INTEGER IOSTATUS,J1,J2,NRIGIDBODY,ATOMID
      CHARACTER(LEN=10) CHECK1
!$    INTEGER NTHREADS,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,TID 
!$OMP PARALLEL
!$    NTHREADS = OMP_GET_NUM_THREADS()
!$      TID = OMP_GET_THREAD_NUM()
!$    if(TID .eq. 0)then 
!$      write(*,*) 'OMP enabled. Number of threads:', NTHREADS
!$    endif
!$OMP END PARALLEL

      ALLOCATE(XSBM(NATOMS))
      ALLOCATE(YSBM(NATOMS))
      ALLOCATE(ZSBM(NATOMS))
      ALLOCATE(ATOMTYPES(NATOMS))

      pi = 3.14159265358979323846264338327950288419716939937510
      MaxCon=NATOMS*5
      do i=1,NATOMS
        NUMOFEXCLUSIONS(I)=0
      enddo


        MAXSEP=40
        MAXSEPSYS=0
      ALLOCATE(NNCINC(NATOMS*MAXSEP))
      DO I=1,NATOMS
        DO J=1,MAXSEP
          NNCINC((I-1)*MAXSEP+J)=0
        enddo
      enddo

! These lines read in the parameters.
        open(30, file='SBM.INP', status='old', access='sequential')
        write(*,*) 
        write(*,*) 'Reading the forcefield from SBM.INP'
        write(*,*) 
        read(30,*)
        read(30,*)
        read(30,*) PREFACTOR,DC,CONCENTRATION,DHswitch,DHcut
        read(30,*)
        read(30,*) NCswitch, NCcut
        write(*,*) 'NCswitch,NCcut'
        write(*,'(2F10.5)') NCswitch, NCcut
        write(*,*) 'Reading electrostatic parameters'
        write(*,'(5F10.5)') PREFACTOR,DC,CONCENTRATION,DHswitch,DHcut
        read(30,*) NATOMTYPES
        ALLOCATE(TMPSTT(NATOMTYPES))
        do I=1,NATOMTYPES
          read(30,*) TMPARG,rsig,reps
          if(TMPARG .ne. I)then
            write(*,*) 'ERROR: atomtypes must be sequential integers'
            STOP  
          endif
          TMPSTT(TMPARG)=reps*rsig**12
        enddo
        ! using combination rule 1

	Rdiff=NCcut-NCswitch
	alpha=12
        ALLOCATE(STT(NATOMTYPES**2)) 
        ALLOCATE(SA(NATOMTYPES**2)) 
        ALLOCATE(SB(NATOMTYPES**2)) 
        ALLOCATE(SC(NATOMTYPES**2)) 
        do I=1,NATOMTYPES
          do J=1,NATOMTYPES
            POS=(I-1)*NATOMTYPES+J
            STT(POS)=sqrt(TMPSTT(I)*TMPSTT(J))
            SB(POS)=-1.0/Rdiff**3*( 2*alpha*STT(POS)/NCcut**(alpha+1)  + 
     Q         (alpha)*(alpha+1)*STT(POS)*Rdiff/NCcut**(alpha+2))
            SA(POS)=-(alpha*(alpha+1)*STT(POS)/NCcut**(alpha+2)+
     Q         3*SB(POS)*Rdiff**2)/(2*Rdiff)
            SC(POS)=-(STT(POS)/NCcut**alpha +
     Q         SA(POS)/3.0*Rdiff**3+SB(POS)/4.0*Rdiff**4)
          enddo 
        enddo

        DEALLOCATE(TMPSTT)


        ! CONCENTRATION IS THE MONOVALENT ION CONCENTRATION kappa is in units
        ! A^-1
        KAPPA=0.32*sqrt(CONCENTRATION)
        PREFACTOR=PREFACTOR/DC
        read(30,*) ANtemp
        write(*,*) ANtemp, 'atoms'
        if(NATOMS .NE. ANTEMP)THEN
          write(*,*) 'ERROR: Number of atoms in SBM.INP and odata are not consistent'
          write(*,*) 'SBM.INP=',ANTEMP, 'odata',NATOMS
          STOP
        ENDIF
        NUMOFSBMCHARGES=0
        ALLOCATE(SBMCHARGE(ANTEMP))
        ALLOCATE(SBMCHARGEON(ANTEMP))
        do i=1, ANtemp
          read(30,*) TMPINT,ATYPE,TMPINT,TMPCHAR,TMPCHAR,SBMCHARGE(i),ATMASS(I)
            if(ATYPE .gt. NATOMTYPES)THEN
              write(*,*) 'ERROR: Unknown atomtype',ATYPE
            endif
            ATOMTYPES(I)=ATYPE
          if(SBMCHARGE(i) .ne. 0)then
             NUMOFSBMCHARGES=NUMOFSBMCHARGES+1
             SBMCHARGEON(NUMOFSBMCHARGES)=i
          endif 
        end do
        
        read(30,*) 
        read(30,*) NC
        write(*,*) NC, 'contacts'        
          if(NC .gt. MaxCon)then
             write(*,*) 'too many contacts'
             STOP
          endif

        ALLOCATE(IC(NC))
        ALLOCATE(JC(NC))
        ALLOCATE(CONTACTTYPE(NC))
        ALLOCATE(CONCOEF(NC*6))

        do i=1, NC
          read(30, *) IC(i), JC(i), CONTACTTYPE(I),ARG1, ARG2
          !! Add more contact types
          TMPARG=0
          IF(CONTACTTYPE(I) .NE. 5)THEN
            CALL INCLUDEEXCLUSIONS(NATOMS,IC(I),JC(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.FALSE.)
          ENDIF
          if(CONTACTTYPE(I) .eq. 1)then ! 6-12 interaction
            CONCOEF(I*6-5)=ARG2*ARG1**12
            CONCOEF(I*6-4)=2*ARG2*ARG1**6
          elseif(CONTACTTYPE(I) .eq. 2)then ! 10-12
            CONCOEF(I*6-5)=5*ARG2*ARG1**12
            CONCOEF(I*6-4)=6*ARG2*ARG1**10
          elseif(CONTACTTYPE(I) .eq. 5)then ! bare Gaussian
            CONCOEF(I*6-5)=ARG1
            CONCOEF(I*6-4)=ARG2
            read(30, *) ARG1 ! need one additional parameter
            CONCOEF(I*6-3)=1/(2*ARG1**2)
          elseif(CONTACTTYPE(I) .eq. 6)then ! Gaussian with wall
            CONCOEF(I*6-5)=ARG1
            CONCOEF(I*6-4)=ARG2
            read(30, *) ARG1,ARG2 ! need two additional parameters
            CONCOEF(I*6-3)=1/(2*ARG1**2)
            CONCOEF(I*6-2)=ARG2/CONCOEF(I*6-5)
          elseif(CONTACTTYPE(I) .eq. 7)then ! Dual gaussian
            CONCOEF(I*6-5)=ARG1
            CONCOEF(I*6-4)=ARG2
            read(30, *) ARG1,ARG2,ARG3,ARG4 ! need four additional parameters
            CONCOEF(I*6-3)=1/(2*ARG1**2)
            CONCOEF(I*6-2)=ARG2
            CONCOEF(I*6-1)=1/(2*ARG3**2)
            CONCOEF(I*6)=ARG4/CONCOEF(I*6-5)
          else
            write(*,*) 'ERROR: Unrecognized contacttype:',CONTACTTYPE(I)
            STOP
          endif
        end do
! make a function to check and add

          read(30,*)
          read(30,*) NBA
          write(*,*) NBA, 'bonds'
          ALLOCATE(IB1(NBA))
          ALLOCATE(IB2(NBA))
          ALLOCATE(RB(NBA))
          ALLOCATE(BK(NBA))
        do i=1, nBA
          read(30,*) Ib1(i), Ib2(i),bondtype,Rb(i), bK(i)

          !!! ADD BONDS OF TYPE 6
          IF(max(IB1(I),IB2(I))-min(IB1(I),IB2(I)) .gt. MAXSEP .AND. bondtype .EQ. 1)THEN
            WRITE(*,*) 'WARNING: DIFFERENCE IN BONDED ATOM INDICES IS GREATER THAN MAXSEP'
            WRITE(*,*) '         THIS CAN LEAD TO POOR PERFORMANCE.'
            WRITE(*,*) '         CONSIDER RECOMPILING WITH A LARGER VALUE OF MAXSEP.'
          ENDIF
          IF(bondtype .eq. 1 .or. bondtype .eq. 6)THEN
            CALL INCLUDEEXCLUSIONS(NATOMS,IB1(I),IB2(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          ELSEIF(bondtype .eq. 6)THEN
            CALL INCLUDEEXCLUSIONS(NATOMS,IB1(I),IB2(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.FALSE.)          
          ENDIF
        end do

! organizing for rapid hessian evaluations
        ALLOCATE(SBMNB(NATOMS))
        do i=1,NATOMS 
          SBMNB(I)=0
        enddo

        do i=1, nBA
          SBMNB(IB1(I))=SBMNB(IB1(I))+1
          SBMNB(IB2(I))=SBMNB(IB2(I))+1
        enddo
        MAXBONDSPERATOM=0
        do i=1,NATOMS
          IF(SBMNB(I) .GT. MAXBONDSPERATOM)THEN
                MAXBONDSPERATOM=SBMNB(I)
          ENDIF
        ENDDO
          ALLOCATE(SBMBONDLIST(NATOMS,MAXBONDSPERATOM))

        do i=1,NATOMS 
          SBMNB(I)=0
        enddo

        do i=1, nBA
          TT1=IB1(I)
          TT2=IB2(I)
          SBMNB(TT1)=SBMNB(TT1)+1
          SBMNB(TT2)=SBMNB(TT2)+1
          SBMBONDLIST(TT1,SBMNB(TT1))=I
          SBMBONDLIST(TT2,SBMNB(TT2))=I
        enddo

          read(30,*)
          read(30,*) NTA
          write(*,*) NTA, 'angles'
        ALLOCATE(IT(NTA))
        ALLOCATE(JT(NTA))
        ALLOCATE(KT(NTA))
        ALLOCATE(ANTC(NTA))
        ALLOCATE(TK(NTA))
        do i=1, nTA
          read(30,*) IT(i), JT(i), KT(i), ANTC(i), TK(i)
          CALL INCLUDEEXCLUSIONS(NATOMS,IT(I),JT(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,IT(I),KT(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,JT(I),KT(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          IF(max(IT(I),JT(I))-min(IT(I),JT(I)) .gt. MAXSEP .OR.
     Q       max(IT(I),KT(I))-min(IT(I),KT(I)) .gt. MAXSEP .OR.
     Q       max(KT(I),JT(I))-min(KT(I),JT(I)) .gt. MAXSEP )THEN
            WRITE(*,*) 'WARNING: DIFFERENCE IN BOND ANGLE ATOM INDICES IS GREATER THAN MAXSEP'
            WRITE(*,*) '         THIS CAN LEAD TO POOR PERFORMANCE.'
            WRITE(*,*) '         CONSIDER RECOMPILING WITH A LARGER VALUE OF MAXSEP.'
          ENDIF
        ENDDO

          read(30,*) 
          read(30,*) NPA
          write(*,*) NPA, 'dihedrals'
          ALLOCATE(IP(NPA))
          ALLOCATE(JP(NPA))
          ALLOCATE(KP(NPA))
          ALLOCATE(LP(NPA))
          ALLOCATE(PK(NPA))
          ALLOCATE(PHITYPE(NPA))
          ALLOCATE(PHISBM(NPA))
! this reads in the dihedral angles and calculates the cosines and sines
! in order to make the force and energy calculations easier, later.
        do i=1, npA
          read(30,*) IP(i),JP(i),KP(i),LP(i),PHITYPE(i),PHISBM(i),PK(i)
          CALL INCLUDEEXCLUSIONS(NATOMS,IP(I),JP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,IP(I),KP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,IP(I),LP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,JP(I),KP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,JP(I),LP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          
          CALL INCLUDEEXCLUSIONS(NATOMS,KP(I),LP(I),EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.TRUE.)          

          IF(max(IP(I),JP(I))-min(IP(I),JP(I)) .gt. MAXSEP .OR.
     Q       max(IP(I),KP(I))-min(IP(I),KP(I)) .gt. MAXSEP .OR.
     Q       max(IP(I),LP(I))-min(IP(I),LP(I)) .gt. MAXSEP .OR.
     Q       max(JP(I),KP(I))-min(JP(I),KP(I)) .gt. MAXSEP .OR.
     Q       max(JP(I),LP(I))-min(JP(I),LP(I)) .gt. MAXSEP .OR.
     Q       max(KP(I),LP(I))-min(KP(I),LP(I)) .gt. MAXSEP )THEN
            WRITE(*,*) 'WARNING: DIFFERENCE IN DIHEDRAL ANGLE ATOM INDICES IS GREATER THAN MAXSEP'
            WRITE(*,*) '         THIS CAN LEAD TO POOR PERFORMANCE.'
            WRITE(*,*) '         CONSIDER RECOMPILING WITH A LARGER VALUE OF MAXSEP.'
          ENDIF
        ENDDO 

        read(30,*)
        read(30,*) nexc
        write(*,*) nexc, 'exculusions'

      do i=1, nexc
        read(30,*) I1, I2
        CALL INCLUDEEXCLUSIONS(NATOMS,I1,I2,EXCLUSIONS,NUMOFEXCLUSIONS,MAXEXCLUSIONS,.FALSE.)          
      enddo

! now that we know what to exclude, let's make a static list of pairs to
! exclude
        NEXCLUSIONS=0
        DO I=1,NATOMS
           NEXCLUSIONS=NEXCLUSIONS+NUMOFEXCLUSIONS(I)
        ENDDO
        ALLOCATE(NNEXL1(NEXCLUSIONS))
        ALLOCATE(NNEXL2(NEXCLUSIONS))
        ALLOCATE(NNEXL1ELEC(NEXCLUSIONS))
        ALLOCATE(NNEXL2ELEC(NEXCLUSIONS))
        ! reset counters
        NEXCLUSIONS=0
        NEXCLUSIONSELEC=0
        DO I=1,NATOMS
          DO J=1,NUMOFEXCLUSIONS(I)
            NEXCLUSIONS=NEXCLUSIONS+1
            M=EXCLUSIONS((I-1)*MAXEXCLUSIONS+J)
            NNEXL1(NEXCLUSIONS)=I
            NNEXL2(NEXCLUSIONS)=M
            IF(SBMCHARGE(I) .NE. 0 .AND. SBMCHARGE(M) .NE. 0)THEN
              NEXCLUSIONSELEC=NEXCLUSIONSELEC+1
              NNEXL1ELEC(NEXCLUSIONSELEC)=I
              NNEXL2ELEC(NEXCLUSIONSELEC)=M
            ENDIF
          ENDDO
        ENDDO
! At this point, we will change the meaning of maxexclusions, which is
! actually NEXCLUDE
        do I=1,NATOMS
          TARR(I)=.FALSE.
        enddo

! read in position restraints
       read(30,*) 
       read(30,*) SBMPRN
       write(*,*) SBMPRN, 'position restraints'

      ALLOCATE(SBMPRI(SBMPRN))
      ALLOCATE(SBMPRK(SBMPRN*6))
      ALLOCATE(SBMPRX(SBMPRN*3))

       do I=1,SBMPRN
            J=(I-1)*6
            K=J/2
            read(30,*) SBMPRI(I),SBMPRK(J+1),SBMPRK(J+2),SBMPRK(J+3),
     Q         SBMPRK(J+4),SBMPRK(J+5),SBMPRK(J+6),
     Q         SBMPRX(K+1),SBMPRX(K+2),SBMPRX(K+3)
            if(TARR(SBMPRI(I)))then
                write(*,*) 'more than one restraint provided for atom ',SBMPRI(I)
                STOP
            endif
            TARR(SBMPRI(I))=.TRUE.
        if(SBMPRK(J+4) .ne. 0 .or. SBMPRK(J+5) .ne. 0 .or. SBMPRK(J+6) .ne.0)then  
          write(*,*) 'FATAL ERROR: Non-zero restraint cross-terms not supported'
        endif
       enddo 
       close(30)


! If RIGIDINIT is on, then we will want to keep information about rigid groups
! so that we may ignore their interactions later.

      ALLOCATE (SBMRBG(NATOMS))
      DO I=1,NATOMS
	SBMRBG(I)=-I
      ENDDO
        NRIGIDBODY=0
! approach to I/O for rbodyconfig taken from genrigid.f90
      if(RIGIDINIT)THEN
        OPEN(UNIT=222,FILE='rbodyconfig',status='old')
        DO
          READ(222,*,IOSTAT=iostatus) CHECK1
          IF (iostatus<0) THEN
            CLOSE(222)
            EXIT
          ELSE IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') then
            NRIGIDBODY=NRIGIDBODY+1
          ENDIF
        END DO
        CLOSE(222)
  
        write(*,*) NRIGIDBODY, ' rigid bodies will be used'
 
        OPEN(UNIT=222,FILE='rbodyconfig',status='old')
        DO J1 = 1, NRIGIDBODY
          READ(222,*) CHECK1,DUMMY
          DO J2 = 1, DUMMY
            READ(222,*) ATOMID
            SBMRBG(ATOMID)=J1
          ENDDO
        ENDDO
        CLOSE(222)
      ENDIF

      END 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of SBMinit^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!**************************************************************************
! INCLUDEEXCLUSION helps keep track of which exclusions to include in
! the model
!**************************************************************************

      SUBROUTINE INCLUDEEXCLUSIONS(NATOMS,ATOM1,ATOM2,EXCLIST,NEXCLUDE,
     Q MAXEXCLUDE,CHECKSEP)
      USE SBMDATA
      IMPLICIT NONE
      LOGICAL CHECKSEP
      INTEGER NATOMS,ATOM1,ATOM2,EXCLIST,NEXCLUDE,MAXEXCLUDE
      DIMENSION EXCLIST(NATOMS*MAXEXCLUDE),NEXCLUDE(NATOMS)
      INTEGER I,N,M,AA,BB,AB,POS
      LOGICAL INCLUDED
! If the atoms are within MAXSEP of each other, then keep track via
! NNCINC
      AA=min(ATOM1,ATOM2)
      BB=max(ATOM1,ATOM2)
      AB=BB-AA
      POS=(AA-1)*MAXSEP+AB
      IF(AB .le. MAXSEP)THEN
        NNCINC(POS)=1
        IF(AB .GT. MAXSEPSYS .AND. CHECKSEP)THEN
          MAXSEPSYS=AB
        ENDIF
      ELSE
        INCLUDED = .FALSE.
        N=MAX(ATOM1,ATOM2)
        M=MIN(ATOM1,ATOM2)
        DO I =1,NEXCLUDE(M)
          IF(EXCLIST((M-1)*MAXEXCLUDE+I) .eq. N)THEN
            INCLUDED = .TRUE.
          ENDIF
        ENDDO 
  
        if(.NOT. INCLUDED)THEN
          if(NEXCLUDE(M) .EQ. MAXEXCLUDE)THEN
            write(*,*) 'ERROR: TOO MANY EXCLUSIONS WITH ATOM',M
          ENDIF
          NEXCLUDE(M)=NEXCLUDE(M)+1
          EXCLIST((M-1)*MAXEXCLUDE+NEXCLUDE(M))=N
        ENDIF
      ENDIF

      END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


C
C Calculate the Forces and energies
C
      subroutine CALC_ENERGY_SBM(qo,natoms,GRAD,energy)
      USE SBMDATA

      INTEGER I,J,NATOMS

      DOUBLE PRECISION qo(3*NATOMS), grad(3*NATOMS), ENERGY

      DOUBLE PRECISION dx,dy,dz
      do i = 1, natoms
        j = (i-1)*3
        XSBM(i) = qo(j+1)
        YSBM(i) = qo(j+2)
        ZSBM(i) = qo(j+3)
        grad(j+1) = 0.0
        grad(j+2) = 0.0
        grad(j+3) = 0.0
      enddo
      energy = 0.0
      call SBMbonds(grad, energy, natoms)
      call SBMangl(grad, energy, natoms)
      call SBMDihedral(grad, energy, natoms)
      call SBMContacts(grad, energy,natoms)
      call SBMNonContacts(grad,energy,natoms)
      call SBMDHELEC(grad,energy,natoms)
      call SBMPR(grad, energy, natoms)
      end


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMBonds  computes the hookean force and energy between chosen atoms *
!***********************************************************************

      subroutine SBMBonds(grad,energy, natoms)
      USE KEY
      USE SBMDATA
      implicit NONE
      integer I2,J2,I23,J23,I,N,J,NATOMS
      DOUBLE PRECISION grad(3*NATOMS),energy
      DOUBLE PRECISION r2, f, r1
      DOUBLE PRECISION dx,dy,dz

      IF(SBMHESSATOM .EQ. -1)THEN 
!$OMP PARALLEL PRIVATE(I,I2,J2,I23,J23,DX,DY,DZ,R1,R2,F)REDUCTION(+:ENERGY,grad)
!$OMP DO
        DO I=1, nBA
           I2 = Ib1(I)
           J2 = Ib2(I)
             I23=I2*3
             J23=J2*3
             dx = XSBM(I2) - XSBM(J2)
             dy = YSBM(I2) - YSBM(J2)
             dz = ZSBM(I2) - ZSBM(J2)
             r2 = dx**2 + dy**2 + dz**2
             r1 = sqrt(r2)
             ENERGY = ENERGY + bk(I)*(r1-Rb(I))**2/2.0
             f = -bk(I)*(r1-Rb(I))/r1
             grad(I23-2) = grad(I23-2) - f * dx
             grad(I23-1) = grad(I23-1) - f * dy
             grad(I23)   = grad(I23)   - f * dz
             grad(J23-2) = grad(J23-2) + f * dx
             grad(J23-1) = grad(J23-1) + f * dy
             grad(J23)   = grad(J23)   + f * dz
       ENDDO 
!$OMP END DO
!$OMP END PARALLEL

      ELSE
        DO I=1, SBMNB(SBMHESSATOM)
           I2 = IB1(SBMBONDLIST(SBMHESSATOM,I))
           J2 = IB2(SBMBONDLIST(SBMHESSATOM,I))
             I23=I2*3
             J23=J2*3
             dx = XSBM(I2) - XSBM(J2)
             dy = YSBM(I2) - YSBM(J2)
             dz = ZSBM(I2) - ZSBM(J2)
             r2 = dx**2 + dy**2 + dz**2
             r1 = sqrt(r2)
             ENERGY = ENERGY + bk(I)*(r1-Rb(I))**2/2.0
             f = -bk(I)*(r1-Rb(I))/r1
             grad(I23-2) = grad(I23-2) - f * dx
             grad(I23-1) = grad(I23-1) - f * dy
             grad(I23)   = grad(I23)   - f * dz
             grad(J23-2) = grad(J23-2) + f * dx
             grad(J23-1) = grad(J23-1) + f * dy
             grad(J23)   = grad(J23)   + f * dz
       ENDDO 

      ENDIF

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF SBMBONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMANGL  computes the Force due to the bond angles                   *
!* This code is modeled after how AMBER performs angle forces           *
!***********************************************************************

      SUBROUTINE SBMANGL(grad, energy, NATOMS)
      USE KEY
      USE SBMDATA
      IMPLICIT NONE
      integer NATOMS,C1
      DOUBLE PRECISION grad(3*NATOMS),energy
      LOGICAL SKIP,NOCRST

      DOUBLE PRECISION RIJ,RKJ,RIK,ANT,XIJ,YIJ,
     + ZIJ,XKJ,YKJ,ZKJ, DF
      DOUBLE PRECISION CT0, CT1, CT2, DA, ST, 
     + CIK, CII, CKK, DT1, DT2, DT3, DT4, DT5, DT6, DT7, DT8, DT9,  STH
      DOUBLE PRECISION ONE,NEGONE

        INTEGER II
        INTEGER I3,J3,K3,I33,J33,K33
        ONE=1.0
        NEGONE=-1.0
!$OMP PARALLEL PRIVATE(I3,J3,K3,I33,J33,K33,RIJ,RKJ,RIK,ANT,XIJ,YIJ,ZIJ,XKJ,YKJ,
!$OMP&  ZKJ,DF,CT0,CT1,CT2,DA,ST,CIK,CII,CKK,DT1,DT2,DT3,DT4,
!$OMP&  DT5,DT6,DT7,DT8,DT9,STH) REDUCTION(+:ENERGY,grad)
!$OMP DO
          DO II = 1, nTA
            I3 = IT(II)
            J3 = JT(II)
            K3 = KT(II)
            IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. I3 .OR.
     Q        SBMHESSATOM .EQ. J3 .OR. SBMHESSATOM .EQ. K3)THEN
              I33=I3*3
              J33=J3*3
              K33=K3*3
              XIJ = XSBM(I3)-XSBM(J3)
              YIJ = YSBM(I3)-YSBM(J3)
              ZIJ = ZSBM(I3)-ZSBM(J3)
              XKJ = XSBM(K3)-XSBM(J3)
              YKJ = YSBM(K3)-YSBM(J3)
              ZKJ = ZSBM(K3)-ZSBM(J3)
              RIJ = XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
              RKJ = XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ
              RIK = SQRT(RIJ*RKJ)
              CT0 = (XIJ*XKJ+YIJ*YKJ+ZIJ*ZKJ)/RIK
              CT1 = MAX(NEGONE,CT0)
              CT2 = MIN(ONE,CT1)
              ANT = ACOS(CT2)
              DA = ANT - ANTC(II)
              DF = TK(II)*DA
              ST = -(DF)/SIN(ANT)
              STH = ST*CT2
              CIK = ST/RIK
              CII = STH/RIJ
              CKK = STH/RKJ
              DT1 = CIK*XKJ-CII*XIJ
              DT2 = CIK*YKJ-CII*YIJ
              DT3 = CIK*ZKJ-CII*ZIJ
              DT7 = CIK*XIJ-CKK*XKJ
              DT8 = CIK*YIJ-CKK*YKJ
              DT9 = CIK*ZIJ-CKK*ZKJ
              DT4 = -DT1-DT7
              DT5 = -DT2-DT8
              DT6 = -DT3-DT9
              grad(I33-2) = grad(I33-2)+ DT1
              grad(I33-1) = grad(I33-1)+ DT2
              grad(I33)   = grad(I33)  + DT3
              grad(J33-2) = grad(J33-2)+ DT4
              grad(J33-1) = grad(J33-1)+ DT5
              grad(J33)   = grad(J33)  + DT6
              grad(K33-2) = grad(K33-2)+ DT7
              grad(K33-1) = grad(K33-1)+ DT8
              grad(K33)   = grad(K33)  + DT9
              ENERGY = ENERGY + TK(II)*(ANTC(II)- ANT)**2/2.0
            ENDIF
          END DO
!$OMP END DO
!$OMP END PARALLEL

       RETURN
       END

!^^^^^^^^^^^^^^^^^^^^^^^^End of SBMANGL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMdihedral computes the dihedral angles and the forces due to them *
!**********************************************************************

      SUBROUTINE SBMdihedral(grad,energy, NATOMS)
      USE KEY
      USE SBMDATA
      implicit NONE
      integer I, N, J, NATOMS, II
      DOUBLE PRECISION grad(3*NATOMS),energy

      double precision lfac
      integer I3,J3,K3,L3,C1,I33,J33,K33,L33
      double precision  XIJ,YIJ,ZIJ,XKJ,YKJ,
     + ZKJ,XKL,YKL,ZKL,RIJ, RKJ,RKL,DX,DY,
     + DZ, GX,GY,GZ,CT,CPHI,
     + SPHI,Z1, Z2,FXI,FYI,FZI,
     + FXJ,FYJ,FZJ, FXK,FYK,FZK,
     + FXL,FYL,FZL,DF,Z10,Z20,Z12,Z11,Z22,ftem,CT0,CT1,AP0,AP1,
     + Dums,DFLIM, DF1, DF0, DR1, DR2,DR3,DR4,DR5,DR6,DRX,DRY,DRZ,
     + S,HGoverG,FGoverG,A1,A3

      double precision  TM24,TM06,tenm3,zero,one,NEGONE,two,four,six,twelve

      DOUBLE PRECISION TT1, TT2, TT3, TT4, TT1X,TT1Y,TT1Z,TT2X,TT2Y,
     + TT2Z, TT3X, TT3Y, TT3Z, TT4X, TT4Y, TT4Z

      DATA TM24,TM06,tenm3/1.0d-24,1.0d-06,1.0d-03/

      double precision pi
      pi = 3.14159265358979323846264338327950288419716939937510
      ONE=1.0
      NEGONE=-1.0

!$OMP PARALLEL PRIVATE(I,I3,J3,K3,L3,I33,J33,K33,L33,II,XIJ,YIJ,ZIJ,XKJ,YKJ, 
!$OMP& ZKJ,XKL,YKL,ZKL,RIJ,RKJ,RKL,DX,DY, 
!$OMP& DZ, GX,GY,GZ,CT,CPHI,
!$OMP& SPHI,Z1, Z2,FXI,FYI,FZI,
!$OMP& FXJ,FYJ,FZJ, FXK,FYK,FZK,
!$OMP& FXL,FYL,FZL,DF,Z10,Z20,Z12,Z11,Z22,CT0,CT1,AP0,AP1,
!$OMP& Dums,DFLIM, DF1, DF0, DR1, DR2,DR3,DR4,DR5,DR6,DRX,DRY,DRZ,
!$OMP& S,HGoverG,FGoverG,A1,A3,TT1,TT2,TT3,TT4,TT1X,TT1Y,TT1Z,TT2X,
!$OMP& TT2Y,TT2Z,TT3X,TT3Y,TT3Z,TT4X,TT4Y,TT4Z) REDUCTION(+:energy,grad)
!!
!$OMP DO
          DO II = 1,nPA

            I3 = IP(II)
            J3 = JP(II)
            K3 = KP(II)
            L3 = LP(II)
            IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. I3 .OR.
     Q        SBMHESSATOM .EQ. J3 .OR. SBMHESSATOM .EQ. K3 .OR. 
     Q        SBMHESSATOM .EQ. L3)THEN 
              I33=I3*3
              J33=J3*3
              K33=K3*3
              L33=L3*3
              XIJ = XSBM(I3)-XSBM(J3)
              YIJ = YSBM(I3)-YSBM(J3)
              ZIJ = ZSBM(I3)-ZSBM(J3)
              XKJ = XSBM(K3)-XSBM(J3)
              YKJ = YSBM(K3)-YSBM(J3)
              ZKJ = ZSBM(K3)-ZSBM(J3)
              RKJ = sqrt(XKJ**2+YKJ**2+ZKJ**2)
              XKL = XSBM(K3)-XSBM(L3)
              YKL = YSBM(K3)-YSBM(L3)
              ZKL = ZSBM(K3)-ZSBM(L3)                                  
              FGoverG=-(XIJ*XKJ+YIJ*YKJ+ZIJ*ZKJ)/RKJ
              HGoverG=(XKL*XKJ+YKL*YKJ+ ZKL*ZKJ)/RKJ
C DX is the M vector and G is the N vector
              DX = YIJ*ZKJ-ZIJ*YKJ
              DY = ZIJ*XKJ-XIJ*ZKJ
              DZ = XIJ*YKJ-YIJ*XKJ
              GX = ZKJ*YKL-YKJ*ZKL
              GY = XKJ*ZKL-ZKJ*XKL
              GZ = YKJ*XKL-XKJ*YKL
              FXI = SQRT(DX*DX+DY*DY+DZ*DZ)
              FYI = SQRT(GX*GX+GY*GY+GZ*GZ)
              CT = DX*GX+DY*GY+DZ*GZ
              z10 = 1.0/FXI
              z20 = 1.0/FYI
              Z12 = Z10*Z20
              Z1 = Z10
              Z2 = Z20
              CT0 = MIN(ONE,CT*Z12)
              CT1 = MAX(NEGONE,CT0)
              S = XKJ*(DZ*GY-DY*GZ)+YKJ*(DX*GZ-DZ*GX)+ZKJ*(DY*GX-DX*GY)
              AP0 = ACOS(CT1)
              AP1 = PI-SIGN(AP0,S)
              CT = AP1
              CPHI = COS(AP1)
              SPHI = SIN(AP1)
! Here is the energy part
            A1=CT-PHISBM(II)
            A3=A1*3
  
          if(PHITYPE(II) .eq. 1)then
            energy =  energy + PK(II)*(3.0/2.0-cos(A1)-0.5*cos(A3))
! dE/dPHI
              DF=PK(II)*(sin(A1)+1.5*sin(A3))
          elseif(PHITYPE(II) .eq. 2)then
            if(A1 .gt. PI)then
                  A1=A1-2*PI
            elseif(A1 .lt. -PI)then
                  A1=A1+2*PI
            endif
  
            energy =  energy + 0.5*PK(II)*A1**2
! dE/dPHI
              DF=PK(II)*A1
  
          else
  
          write(*,*) 'unrecognized dihedral type', PHITYPE(II)
          STOP
          endif
  
! insert the new 
  
! now, do dPhi/dX
  
! |G|/|A|**2 
            TT1 = Z1*Z1*RKJ*DF
! FG/(A**2*|G|)
            TT2 = FGoverG*Z1*Z1*DF
! HG/(B**2*|G|)
            TT3 = HGoverG*Z2*Z2*DF
! |G|/|B|**2 
            TT4 = Z2*Z2*RKJ*DF
! note: negatives are flipped from paper because A=-DX
            TT1X=TT1*DX
            TT1Y=TT1*DY
            TT1Z=TT1*DZ
            TT2X=TT2*DX
            TT2Y=TT2*DY
            TT2Z=TT2*DZ
            TT3X=TT3*GX
            TT3Y=TT3*GY
            TT3Z=TT3*GZ
            TT4X=TT4*GX
            TT4Y=TT4*GY
            TT4Z=TT4*GZ
            grad(I33-2) =  grad(I33-2)  + TT1X  
            grad(I33-1) =  grad(I33-1)  + TT1Y
            grad(I33)   =  grad(I33)    + TT1Z
            grad(J33-2) =  grad(J33-2)  - TT1X - TT2X - TT3X
            grad(J33-1) =  grad(J33-1)  - TT1Y - TT2Y - TT3Y
            grad(J33)   =  grad(J33)    - TT1Z - TT2Z - TT3Z
            grad(K33-2) =  grad(K33-2)  + TT2X + TT3X - TT4X
            grad(K33-1) =  grad(K33-1)  + TT2Y + TT3Y - TT4Y
            grad(K33)   =  grad(K33)    + TT2Z + TT3Z - TT4Z
            grad(L33-2) =  grad(L33-2)  + TT4X
            grad(L33-1) =  grad(L33-1)  + TT4Y
            grad(L33)   =  grad(L33)    + TT4Z
          ENDIF
        END DO
!$OMP END DO 
!$OMP END PARALLEL

          END


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END of SBMDihedral^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMCONTACTS: computes the force on all atoms due to contacts via a   *
!* 10-12 or 6-12 potential                                              *
!***********************************************************************

      subroutine SBMcontacts(grad,energy,NATOMS)
      USE KEY
      USE SBMDATA
      implicit NONE
      integer I,J,NATOMS

      DOUBLE PRECISION grad(3*NATOMS), energy
      DOUBLE PRECISION dx,dy,dz,CO1,CO2,CO3,CO4,CO5,CO6,E1,G1,G2,E1D,G1D,G2D,ET

      integer C1,C2,C13,C23 
      DOUBLE PRECISION  r2,rm2,rm10,rm12,rm6,f_over_r,r1 

! CONTACTTYPE
! type 1 is 6-12
! type 2 is 12-12
! type 5 as gaussian, no ex-vol
! type 6 as gaussian, w ex-vol
! type 7 as dual gaussian
!$OMP PARALLEL PRIVATE(I,CO1,CO2,CO3,CO4,CO5,CO6,C1,C2,C13,C23,DX,DY,DZ,
!$OMP& R2,RM2,RM6,RM10,RM12,R1,E1,G1,G2,E1D,G1D,G2D,F_OVER_R,ET) REDUCTION(+:ENERGY,grad)
!$OMP DO
          do i=1, NC
            C1 = IC(i)
            C2 = JC(i)
            IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. C1 .OR.
     Q        SBMHESSATOM .EQ. C2)THEN 
              dx = XSBM(C1) - XSBM(C2)
              dy = YSBM(C1) - YSBM(C2)
              dz = ZSBM(C1) - ZSBM(C2)
              r2 = dx**2 + dy**2 + dz**2
              J=I*6
              C13=C1*3
              C23=C2*3
              if(CONTACTTYPE(I) .eq. 1)then
                ! type 1 is 6-12 interaction
                rm2 = 1.0/r2
                rm6 = rm2**3
                CO1=CONCOEF(J-5)
                CO2=CONCOEF(J-4)
                ENERGY = ENERGY + rm6*(CO1*rm6-CO2)
                f_over_r = -12.0*rm6*(CO1*rm6-CO2/2.0)*rm2
              elseif(CONTACTTYPE(I) .eq. 2)then
                ! type 2 is 10-12 interaction
                rm2 = 1.0/r2
                rm10 = rm2**5
                CO1=CONCOEF(J-5)
                CO2=CONCOEF(J-4)
                ENERGY = ENERGY + rm10*(CO1*rm2-CO2)
                f_over_r = -60.0*rm10*(CO1*rm2/5.0-CO2/6.0)*rm2
              elseif(CONTACTTYPE(I) .eq. 5)then
                ! type 5 is a simple gaussian interaction (no excluded
                ! volume)
                r1 = sqrt(r2)
                CO1=CONCOEF(J-5) !A
                CO2=CONCOEF(J-4) !mu
                CO3=CONCOEF(J-3) !1/(2*sigma**2)
                E1= -CO1*exp(-CO3*(r1-CO2)**2)
                ENERGY = ENERGY + E1
                f_over_r = -E1*2*CO3*(r1-CO2)/r1
              elseif(CONTACTTYPE(I) .eq. 6)then
                ! type 6 is a gaussian interaction with excluded
                ! volume included
                r1 = sqrt(r2)
                rm2 = 1.0/r2
                rm12 = rm2**6
                CO1=CONCOEF(J-5) !A
                CO2=CONCOEF(J-4) !mu
                CO3=CONCOEF(J-3) !1/(2*sigma**2)
                CO4=CONCOEF(J-2) !a/A
                E1 =CO1*(1+CO4*rm12)
                G1 =1-exp(-CO3*(r1-CO2)**2)
                E1D=-12*CO1*CO4*rm12*rm2
                G1D=CO3*2*(1-CO2/r1)*exp(-CO3*(r1-CO2)**2)
                ENERGY = ENERGY + (E1*G1-CO1)
                f_over_r = E1*G1D + G1*E1D 
              elseif(CONTACTTYPE(I) .eq. 7)then
                ! type 7 is a double gaussian interaction with excluded
                ! volume included
                rm2 = 1.0/r2
                rm12 = rm2**6
                r1 = sqrt(r2)
                CO1=CONCOEF(J-5) !A
                CO2=CONCOEF(J-4) !mu1
                CO3=CONCOEF(J-3) !1/(2*sigma1**2)
                CO4=CONCOEF(J-2) !mu2
                CO5=CONCOEF(J-1) !1/(2*sigma2**2)
                CO6=CONCOEF(J)   ! a/A (excluded volume amplitude/gaussian amplitude)
                E1=1+CO6*rm12
                G1=1-exp(-CO3*(r1-CO2)**2)
                G2=1-exp(-CO5*(r1-CO4)**2)
                E1D= -12*CO6*rm12*rm2
                G1D=CO3*2*(1-CO2/r1)*exp(-CO3*(r1-CO2)**2)
                G2D=CO5*2*(1-CO4/r1)*exp(-CO5*(r1-CO4)**2)
                ENERGY = ENERGY + CO1*( E1 * G1 * G2-1) 
                f_over_r = CO1*(E1D*G1*G2 
     Q          +  E1 * G1D * G2
     Q          +  E1 * G1 * G2D)
              else
                WRITE(*,*) 'ERROR: Contact type not recognized:', CONTACTTYPE(I)
                STOP
              endif
  
              grad(C13-2) = grad(C13-2) + f_over_r * dx
              grad(C13-1) = grad(C13-1) + f_over_r * dy
              grad(C13)   = grad(C13)   + f_over_r * dz
              grad(C23-2) = grad(C23-2) - f_over_r * dx
              grad(C23-1) = grad(C23-1) - f_over_r * dy
              grad(C23)   = grad(C23)   - f_over_r * dz
            ENDIF
          enddo
!$OMP END DO
!$OMP END PARALLEL


      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of SBMContacts^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMNonContacts computes the forces due to non native contacts      *
!**********************************************************************

        subroutine SBMnoncontacts(grad,energy,NATOMS)
        USE KEY
        USE SBMDATA
        implicit NONE
        integer I,I3,N,J,AN,NATOMS

        DOUBLE PRECISION grad(3*NATOMS), energy
        DOUBLE PRECISION SAT,SBT,SCT,STTT

        integer C1,C2,C13,C23,DINDEX,ii,jj,kk,k,l,iii,jjj,POS
        DOUBLE PRECISION  r2, rm2, rm14, f_over_r 
        INTEGER NTHREADS
!$    INTEGER  OMP_GET_NUM_THREADS,
!$   Q  OMP_GET_THREAD_NUM,OMP_NUM_THREADS

        DOUBLE PRECISION dx,dy,dz
	integer tempN, alpha
	double precision Rdiff,Vfunc,Ffunc
	double precision Rcut2,Rswitch2
	! Ngrid is the number of atoms in that grid point
	! grid is the array of atoms in each grid
	!integer Ngrid,grid
	integer maxpergrid
	! number of atoms per grid, max
	parameter (maxpergrid=100)
	! dimensions of grid
!	parameter (maxgrid=30)
!	dimension Ngrid(maxgrid,maxgrid,maxgrid),
!     Q  grid(maxgrid,maxgrid,maxgrid,maxpergrid)
	integer MaxGridX,MaxGridY,MaxGridZ
	double precision gridsize,RD1
	double precision minX,minY,minZ,maxX,maxY,maxZ
	integer Xgrid,Ygrid,Zgrid,TID
	Rdiff=NCcut-NCswitch
	alpha=12
	GRIDSIZE=NCcut*1.01
	Rcut2=NCcut**2
	Rswitch2=NCswitch**2


        
        IF(SBMHESSATOM .EQ. -1 .OR. SBMFIRSTDD)THEN 
          SBMFIRSTDD = .FALSE.

! grid the system
	minX=10000000
	minY=10000000
	minZ=10000000
	
	maxX=-10000000
	maxY=-10000000
	maxZ=-10000000

	do i=1,NATOMS

	   if(XSBM(i) .lt. minX)then
		minX=XSBM(i)
	   elseif(XSBM(i) .gt. maxX)then
		maxX=XSBM(i)
	   endif
	   if(YSBM(i) .lt. minY)then
		minY=YSBM(i)
	   elseif(YSBM(i) .gt. maxY)then
		maxY=YSBM(i)
	   endif
	   if(ZSBM(i) .lt. minZ)then
		minZ=ZSBM(i)
	   elseif(ZSBM(i) .gt. maxZ)then
		maxZ=ZSBM(i)
	   endif

	enddo

	maxgridX=int((maxX-minX)/gridsize)+1
	maxgridY=int((maxY-minY)/gridsize)+1
	maxgridZ=int((maxZ-minZ)/gridsize)+1
        if(ALLOCATED(SBMNGRID))THEN
          DEALLOCATE(SBMNGRID)
        ENDIF
        if(ALLOCATED(SBMGRID))THEN
          DEALLOCATE(SBMGRID)
        ENDIF
        ALLOCATE(SBMNGRID(MAXGRIDX,MAXGRIDY,MAXGRIDZ))
        ALLOCATE(SBMGRID(MAXGRIDX,MAXGRIDY,MAXGRIDZ,MAXPERGRID))

!!  Add a second grid that only includes atoms that are charged

	do i=1,maxgridx
	 do j=1,maxgridy
	  do k=1,maxgridz
		SBMNGRID(i,j,k)=0
	  enddo
	 enddo
	enddo
	do i=1,NATOMS
		Xgrid=int((XSBM(i)-minX)/gridsize)+1
		Ygrid=int((YSBM(i)-minY)/gridsize)+1
		Zgrid=int((ZSBM(i)-minZ)/gridsize)+1
		SBMNGRID(Xgrid,Ygrid,Zgrid)=SBMNGRID(Xgrid,Ygrid,Zgrid)+1
		if(SBMNGRID(Xgrid,Ygrid,Zgrid) .gt. maxpergrid)then
			write(*,*) 'ERROR: too many atoms in a grid'
			write(*,*) SBMNGRID(Xgrid,Ygrid,Zgrid),Xgrid,Ygrid,Zgrid
                        STOP
	        endif
		SBMGRID(Xgrid,Ygrid,Zgrid,SBMNGRID(Xgrid,Ygrid,Zgrid))=i
	enddo

        ENDIF

	   tempN=0

! add a second loop that goes over charged atoms

        IF(SBMHESSATOM .EQ. -1)THEN 

!$OMP PARALLEL
!$OMP& PRIVATE(i,I3,ii,jj,kk,jjj,C1,C2,C13,C23,r2,rm2,rm14,f_over_r,RD1,dx,dy,dz,
!$OMP& Xgrid,Ygrid,Zgrid,TID,DINDEX,POS,SAT,SBT,SCT,STTT)  
!$OMP& REDUCTION(+:ENERGY,grad)
        TID=0
        NTHREADS=1;
!$      TID = OMP_GET_THREAD_NUM()
!$      NTHREADS = OMP_GET_NUM_THREADS()
        do i=1+TID,NATOMS,NTHREADS


          I3=I*3
	  Xgrid=int((XSBM(I)-minX)/gridsize)+1
	  Ygrid=int((YSBM(I)-minY)/gridsize)+1
	  Zgrid=int((ZSBM(I)-minZ)/gridsize)+1
	  do ii=XGRID-1,XGRID+1
	   do jj=YGRID-1,YGRID+1
	    do kk=ZGRID-1,ZGRID+1
             if(ii .ge. 1 .and. jj .ge. 1 .and. kk .ge. 1 .and.
     Q 	   ii .le. MAXGRIDX .and. jj .le. MAXGRIDY .and. kk .le.MAXGRIDZ)then
              do jjj=1,SBMNGRID(ii,jj,kk)
	       C2=SBMGRID(ii,jj,kk,jjj)
               C23=C2*3
               DINDEX=C2-I
               if(DINDEX .GT. 0)THEN
                if( DINDEX.GT.MAXSEPSYS.OR.
     Q             (DINDEX.LE.MAXSEPSYS.AND.NNCINC((I-1)*MAXSEP+DINDEX) .eq. 0))THEN
                   if(SBMRBG(I) .NE. SBMRBG(C2))THEN
                    dx = XSBM(I) - XSBM(C2)
                    dy = YSBM(I) - YSBM(C2)
                    dz = ZSBM(I) - ZSBM(C2)
                    r2 = dx**2 + dy**2 + dz**2
  	          if(r2 .le. Rcut2)then
                     POS=(ATOMTYPES(I)-1)*NATOMTYPES+ATOMTYPES(C2)
                     STTT=STT(POS)
                     SCT=SC(POS)
                     rm2 = 1/r2
                     rm14 = rm2**7
  	           energy=energy+STTT*rm2**6+SCT
  	           f_over_r=-STTT*12.0*rm14
  	           RD1=sqrt(r2)-NCswitch
  	           if(r2 .gt. Rswitch2)then
                      SAT=SA(POS)
                      SBT=SB(POS)
  	            f_over_r=f_over_r+(SAT*RD1**2+SBT*RD1**3)*sqrt(rm2)
  	            energy=energy+SAT/3.0*RD1**3+SBT/4.0*RD1**4
  	           endif
                     grad(I3-2) = grad(I3-2) + f_over_r * dx
                     grad(I3-1) = grad(I3-1) + f_over_r * dy
                     grad(I3)   = grad(I3)   + f_over_r * dz
                     grad(C23-2) = grad(C23-2) - f_over_r * dx
                     grad(C23-1) = grad(C23-1) - f_over_r * dy
                     grad(C23)   = grad(C23)   - f_over_r * dz
                  ENDIF
 	         endif
	        endif
               endif
              enddo
	     endif
	    enddo
	   enddo
	  enddo
	 enddo

!$OMP DO
       DO I=1,NEXCLUSIONS
         C1=NNEXL1(I) 
         C2=NNEXL2(I)
           C13=C1*3
           C23=C2*3
           dx = XSBM(C1) - XSBM(C2)
           dy = YSBM(C1) - YSBM(C2)
           dz = ZSBM(C1) - ZSBM(C2)
           r2 = dx**2 + dy**2 + dz**2
           if(r2 .le. Rcut2)then
             POS=(ATOMTYPES(C1)-1)*NATOMTYPES+ATOMTYPES(C2)
             STTT=STT(POS)
             SCT=SC(POS)
   
             rm2 = 1/r2
             rm14 = rm2**7
   	   energy=energy-STTT*rm2**6-SCT
   	   f_over_r=-STTT*12.0*rm14
   	   RD1=sqrt(r2)-NCswitch
   	 if(r2 .gt. Rswitch2)then
             SAT=SA(POS)
             SBT=SB(POS)
   	  f_over_r=f_over_r+(SAT*RD1**2+SBT*RD1**3)*sqrt(rm2)
   	  energy=energy-SAT/3.0*RD1**3-SBT/4.0*RD1**4
   	 endif
            grad(C13-2) = grad(C13-2) - f_over_r * dx
            grad(C13-1) = grad(C13-1) - f_over_r * dy
            grad(C13)   = grad(C13)   - f_over_r * dz
            grad(C23-2) = grad(C23-2) + f_over_r * dx
            grad(C23-1) = grad(C23-1) + f_over_r * dy
            grad(C23)   = grad(C23)   + f_over_r * dz
           endif
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
 

         ELSE  ! SBMHESSATOM IS DEFINED

          I=SBMHESSATOM

          I3=I*3
	  Xgrid=int((XSBM(I)-minX)/gridsize)+1
	  Ygrid=int((YSBM(I)-minY)/gridsize)+1
	  Zgrid=int((ZSBM(I)-minZ)/gridsize)+1
	  do ii=XGRID-1,XGRID+1
	   do jj=YGRID-1,YGRID+1
	    do kk=ZGRID-1,ZGRID+1
             if(ii .ge. 1 .and. jj .ge. 1 .and. kk .ge. 1 .and.
     Q 	   ii .le. MAXGRIDX .and. jj .le. MAXGRIDY .and. kk .le.MAXGRIDZ)then
              do jjj=1,SBMNGRID(ii,jj,kk)
	       C2=SBMGRID(ii,jj,kk,jjj)
               DINDEX=C2-I
!               if(DINDEX .GT. 0 .OR. SBMHESSATOM .EQ. I)THEN
                if( C2 .NE. I)THEN
                 C23=C2*3
                if( DINDEX.GT.MAXSEPSYS.OR.
     Q             (DINDEX.LE.MAXSEPSYS.AND.NNCINC((I-1)*MAXSEP+DINDEX) .eq. 0) .OR.
     Q             (DINDEX.LT.-MAXSEPSYS) .OR.
     Q             (DINDEX.GE.-MAXSEPSYS.AND.NNCINC((C2-1)*MAXSEP-DINDEX) .eq. 0)  )then
                   if(SBMRBG(I) .NE. SBMRBG(C2))THEN
                    dx = XSBM(I) - XSBM(C2)
                    dy = YSBM(I) - YSBM(C2)
                    dz = ZSBM(I) - ZSBM(C2)
                    r2 = dx**2 + dy**2 + dz**2
  	          if(r2 .le. Rcut2)then
                     POS=(ATOMTYPES(I)-1)*NATOMTYPES+ATOMTYPES(C2)
                     STTT=STT(POS)
                     SCT=SC(POS)
                     rm2 = 1/r2
                     rm14 = rm2**7
  	           energy=energy+STTT*rm2**6+SCT
  	           f_over_r=-STTT*12.0*rm14
  	           RD1=sqrt(r2)-NCswitch
  	           if(r2 .gt. Rswitch2)then
                      SAT=SA(POS)
                      SBT=SB(POS)
  	            f_over_r=f_over_r+(SAT*RD1**2+SBT*RD1**3)*sqrt(rm2)
  	            energy=energy+SAT/3.0*RD1**3+SBT/4.0*RD1**4
  	           elseif(r2 .lt. Rswitch2)then
  	          ! normal repulsive term
  	           else
  	          ! things should have fallen in one of the previous two...
  	            write(*,*) 'something went wrong with switching function'
  	           endif
                     grad(I3-2) = grad(I3-2) + f_over_r * dx
                     grad(I3-1) = grad(I3-1) + f_over_r * dy
                     grad(I3)   = grad(I3)   + f_over_r * dz
                     grad(C23-2) = grad(C23-2) - f_over_r * dx
                     grad(C23-1) = grad(C23-1) - f_over_r * dy
                     grad(C23)   = grad(C23)   - f_over_r * dz
                  ENDIF
 	         endif
	        endif
               endif
              enddo
	     endif
	    enddo
	   enddo
	  enddo
!!	 enddo


       DO I=1,NEXCLUSIONS
         C1=NNEXL1(I) 
         C2=NNEXL2(I)
         IF(SBMHESSATOM .EQ. C1 .OR. SBMHESSATOM .EQ. C2)THEN 
           C13=C1*3
           C23=C2*3
           dx = XSBM(C1) - XSBM(C2)
           dy = YSBM(C1) - YSBM(C2)
           dz = ZSBM(C1) - ZSBM(C2)
           r2 = dx**2 + dy**2 + dz**2
           if(r2 .le. Rcut2)then
             POS=(ATOMTYPES(C1)-1)*NATOMTYPES+ATOMTYPES(C2)
             STTT=STT(POS)
             SCT=SC(POS)
   
             rm2 = 1/r2
             rm14 = rm2**7
   	   energy=energy-STTT*rm2**6-SCT
   	   f_over_r=-STTT*12.0*rm14
   	   RD1=sqrt(r2)-NCswitch
   	 if(r2 .gt. Rswitch2)then
             SAT=SA(POS)
             SBT=SB(POS)
   	  f_over_r=f_over_r+(SAT*RD1**2+SBT*RD1**3)*sqrt(rm2)
   	  energy=energy-SAT/3.0*RD1**3-SBT/4.0*RD1**4
   	 endif
            grad(C13-2) = grad(C13-2) - f_over_r * dx
            grad(C13-1) = grad(C13-1) - f_over_r * dy
            grad(C13)   = grad(C13)   - f_over_r * dz
            grad(C23-2) = grad(C23-2) + f_over_r * dx
            grad(C23-1) = grad(C23-1) + f_over_r * dy
            grad(C23)   = grad(C23)   + f_over_r * dz
           endif
         ENDIF
       ENDDO
      ENDIF


      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^End of SBMNonContacts^^^^^^^^^^^^^^^^^^^^^^^^^

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! FUNCTIONS for DH Calculations
!*****************************************************

      DOUBLE PRECISION FUNCTION SBMDHV(kappa,r)
      USE KEY
      implicit NONE
      double precision kappa,r
      SBMDHV=exp(-kappa*r)/r
      END

      DOUBLE PRECISION FUNCTION SBMDHVP(kappa,r)
      USE KEY
      implicit NONE
      double precision kappa,r
      SBMDHVP=-kappa*exp(-kappa*r)/r-exp(-kappa*r)/(r**2)
      END
!****************************************************

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMDHELEC computes the forces due to DH electrostatics        *
!**********************************************************************

      subroutine SBMDHELEC(grad,energy,natoms)

      USE KEY
      USE SBMDATA
      implicit NONE
      integer I, N, J, AN, NATOMS

      DOUBLE PRECISION grad(3*NATOMS), energy,
     Q SBMDHVP, SBMDHV
      integer C1,C2,C13,C23,ii,jj,kk,k,l,iii,jjj
      DOUBLE PRECISION  r1,r2, rm2, rm14, f_over_r 
      DOUBLE PRECISION A,B,C, D, COEF2, COEF3
        INTEGER NTHREADS,  OMP_GET_NUM_THREADS,
     Q  OMP_GET_THREAD_NUM,OMP_NUM_THREADS
        DOUBLE PRECISION dx,dy,dz
	integer tempN, alpha
	double precision diff,Vfunc,Ffunc
	double precision DHswitch2,DHcut2
	! Ngrid is the number of atoms in that grid point
	! grid is the array of atoms in each grid
	integer Ngrid,grid,maxgrid,maxpergrid
	! number of atoms per grid, max
	parameter (maxpergrid=100)
	! dimensions of grid
	parameter (maxgrid=15)
	dimension Ngrid(maxgrid,maxgrid,maxgrid),
     Q  grid(maxgrid,maxgrid,maxgrid,maxpergrid)
	integer MaxGridX,MaxGridY,MaxGridZ
	double precision gridsize,RD1
	double precision minX,minY,minZ,maxX,maxY,maxZ
	integer Xgrid,Ygrid,Zgrid,TID
	double precision C2T

      if(NUMOFSBMCHARGES .eq. 0) RETURN

	diff=DHcut-DHswitch
	alpha=12


	GRIDSIZE=DHCUT*1.01
	DHcut2=DHcut**2
	DHswitch2=DHswitch**2

        A=-PREFACTOR*SBMDHV(kappa,DHcut)
        B=-PREFACTOR*SBMDHVP(kappa,DHcut)

        COEF3=(B-2*A/diff)/(diff**2)        
        COEF2=(B-3*diff**2*COEF3)/(2*diff)        
! add DH interactions

! grid the system
	minX=10000000
	minY=10000000
	minZ=10000000
	
	maxX=-10000000
	maxY=-10000000
	maxZ=-10000000

	do ii=1,NUMOFSBMCHARGES
           i=SBMCHARGEON(ii)
	   if(XSBM(i) .lt. minX)then
		minX=XSBM(i)
	   elseif(XSBM(i) .gt. maxX)then
		maxX=XSBM(i)
	   endif
	   if(YSBM(i) .lt. minY)then
		minY=YSBM(i)
	   elseif(YSBM(i) .gt. maxY)then
		maxY=YSBM(i)
	   endif
	   if(ZSBM(i) .lt. minZ)then
		minZ=ZSBM(i)
	   elseif(ZSBM(i) .gt. maxZ)then
		maxZ=ZSBM(i)
	   endif
	enddo

	maxgridX=int((maxX-minX)/gridsize)+1
	maxgridY=int((maxY-minY)/gridsize)+1
	maxgridZ=int((maxZ-minZ)/gridsize)+1

	if(maxgridX .ge. maxgrid .or. 
     Q  maxgridY .ge. maxgrid .or.
     Q  maxgridZ .ge. maxgrid )then
	write(*,*) 'system got too big for grid searching...'
        STOP
	endif


!!  Add a second grid that only includes atoms that are charged

	do i=1,maxgrid
	 do j=1,maxgrid
	  do k=1,maxgrid
		Ngrid(i,j,k)=0
	  enddo
	 enddo
	enddo


	do ii=1,NUMOFSBMCHARGES
                i=SBMCHARGEON(ii)

		Xgrid=int((XSBM(i)-minX)/gridsize)+1
		Ygrid=int((YSBM(i)-minY)/gridsize)+1
		Zgrid=int((ZSBM(i)-minZ)/gridsize)+1
		Ngrid(Xgrid,Ygrid,Zgrid)=Ngrid(Xgrid,Ygrid,Zgrid)+1
		if(Ngrid(Xgrid,Ygrid,Zgrid) .gt. maxpergrid)then
			write(*,*) 'ERROR: too many atoms in a grid in DH'
			write(*,*) Ngrid(Xgrid,Ygrid,Zgrid),Xgrid,Ygrid,Zgrid
                        STOP
	        endif
		grid(Xgrid,Ygrid,Zgrid,Ngrid(Xgrid,Ygrid,Zgrid))=i
	enddo

        tempN=0

!$OMP PARALLEL
!$OMP& PRIVATE(i,ii,iii,jj,kk,jjj,C1,C2,C13,C23,r2,rm2,rm14,f_over_r,RD1,dx,dy,dz,Xgrid,Ygrid,Zgrid,TID,C2T)  
!$OMP& REDUCTION(+:energy,grad)
        TID=0
        NTHREADS=1;
!$      TID = OMP_GET_THREAD_NUM()
!$      NTHREADS = OMP_GET_NUM_THREADS()
         do iii=1+TID,NUMOFSBMCHARGES,NTHREADS
          C1=SBMCHARGEON(iii)
          C13=C1*3
	  Xgrid=int((XSBM(C1)-minX)/gridsize)+1
	  Ygrid=int((YSBM(C1)-minY)/gridsize)+1
	  Zgrid=int((ZSBM(C1)-minZ)/gridsize)+1
	
	  do ii=XGRID-1,XGRID+1
	   do jj=YGRID-1,YGRID+1
	    do kk=ZGRID-1,ZGRID+1
             if(ii .ge. 1 .and. jj .ge. 1 .and. kk .ge. 1 .and.
     Q 	   ii .le. MAXGRIDX .and. jj .le. MAXGRIDY .and. kk .le.MAXGRIDZ)then
              do jjj=1,Ngrid(ii,jj,kk)
	       C2=grid(ii,jj,kk,jjj)
               IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. C1 .OR.
     Q          SBMHESSATOM .EQ. C2)THEN 
                C23=C2*3
                if(C2-C1 .gt. 0)then
                 dx = XSBM(C1) - XSBM(C2)
                 dy = YSBM(C1) - YSBM(C2)
                 dz = ZSBM(C1) - ZSBM(C2)
                 r2 = dx**2 + dy**2 + dz**2
  	         if(r2 .le. DHcut2)then
                  C2T=SBMCHARGE(C1)*SBMCHARGE(C2)
          ! add force evaluation now
                  rm2 = 1/r2
                  r1=sqrt(r2)
  		  energy=energy+PREFACTOR*C2T*SBMDHV(kappa,r1)
  		  f_over_r=PREFACTOR*C2T*SBMDHVP(kappa,r1)
   		   if(r2 .gt. DHswitch2)then
   		    RD1=r1-DHswitch
   		    f_over_r=(f_over_r+C2T*(2*COEF2*RD1+3*COEF3*RD1**2))
   		    energy=energy+COEF2*RD1**2+COEF3*RD1**3
   		   elseif(r2 .lt. DHswitch2)then
   		! normal DH term
   		   else
   		! things should have fallen in one of the previous two...
   		    write(*,*) 'something went wrong with DH switching function'
                    STOP
   		   endif
                   f_over_r=f_over_r*sqrt(rm2)
                   grad(C13-2) = grad(C13-2) + f_over_r * dx
                   grad(C13-1) = grad(C13-1) + f_over_r * dy
                   grad(C13)   = grad(C13)   + f_over_r * dz
                   grad(C23-2) = grad(C23-2) - f_over_r * dx
                   grad(C23-1) = grad(C23-1) - f_over_r * dy
                   grad(C23)   = grad(C23)   - f_over_r * dz
                  ENDIF
 	         endif
	        endif
               enddo
    	      endif
	     enddo
	    enddo
	   enddo
	  enddo

!!! Add routine to subtract exclusion interactions
!$OMP DO
       DO I=1,NEXCLUSIONS
         C1=NNEXL1(I) 
         C2=NNEXL2(I)
         IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. C1 .OR.
     Q     SBMHESSATOM .EQ. C2)THEN 
           C13=C1*3
           C23=C2*3
           dx = XSBM(C1) - XSBM(C2)
           dy = YSBM(C1) - YSBM(C2)
           dz = ZSBM(C1) - ZSBM(C2)
           r2 = dx**2 + dy**2 + dz**2
           if(r2 .le. DHcut2)then
            C2T=SBMCHARGE(C1)*SBMCHARGE(C2)
            rm2 = 1/r2
            r1=sqrt(r2)
            energy=energy-PREFACTOR*C2T*SBMDHV(kappa,r1)
            f_over_r=-PREFACTOR*C2T*SBMDHVP(kappa,r1)
            if(r2 .gt. DHswitch2)then
             RD1=r1-DHswitch
             f_over_r=f_over_r-C2T*(2*COEF2*RD1+3*COEF3*RD1**2)
             energy=energy-COEF2*RD1**2+COEF3*RD1**3
            elseif(r2 .lt. DHswitch2)then
           ! normal DH term
            else
           ! things should have fallen in one of the previous two...
            write(*,*) 'something went wrong with DH switching function'
            endif
            f_over_r=f_over_r*sqrt(rm2)
            grad(C13-2) = grad(C13-2) + f_over_r * dx
            grad(C13-1) = grad(C13-1) + f_over_r * dy
            grad(C13)   = grad(C13)   + f_over_r * dz
            grad(C23-2) = grad(C23-2) - f_over_r * dx
            grad(C23-1) = grad(C23-1) - f_over_r * dy
            grad(C23)   = grad(C23)   - f_over_r * dz
           ENDIF
         ENDIF
       ENDDO
!$OMP END DO
!$OMP END PARALLEL

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^End of SBMDHELEC^^^^^^^^^^^^^^^^^^^^^^^^^

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* SBMPR computes the forces due to position restraints               *
!**********************************************************************

      subroutine SBMPR(grad, energy, natoms)

      USE KEY
      USE SBMDATA
      implicit NONE
      integer I, J, J3,NATOMS,I1,I2
      DOUBLE PRECISION grad(3*NATOMS), energy,ET
      DOUBLE PRECISION DX,DY,DZ,K1,K2,K3,K4,K5,K6


      if(SBMPRN .eq. 0) RETURN

!$OMP PARALLEL PRIVATE(I,J,J3,DX,DY,DZ,K1,K2,K3,K4,K5,K6)REDUCTION(+:energy,grad)
!$OMP DO
      	do I=1,SBMPRN
           J=SBMPRI(I)

         IF(SBMHESSATOM .EQ. -1 .OR. SBMHESSATOM .EQ. J )THEN 
           I1=(I-1)*6
           I2=(I-1)*3
           J3=J*3
           DX=XSBM(J)-SBMPRX(I2+1)
           DY=YSBM(J)-SBMPRX(I2+2)
           DZ=ZSBM(J)-SBMPRX(I2+3)
           K1=SBMPRK(I1+1)
           K2=SBMPRK(I1+2)
           K3=SBMPRK(I1+3)
           K4=SBMPRK(I1+4)
           K5=SBMPRK(I1+5)
           K6=SBMPRK(I1+6)
           
           ENERGY = ENERGY +0.5*(K1*DX**2+K2*DY**2+K3*DZ**2) +
     Q     K4*DX*DY+K5*DX*DZ+K6*DY*DZ

           grad(J3-2) = grad(J3-2) + K1*DX + K4*DY + K5*DZ 
           grad(J3-1) = grad(J3-1) + K2*DY + K4*DX + K6*DZ
           grad(J3)   = grad(J3)   + K3*DZ + K5*DX + K6*DY
         ENDIF
        enddo
!$OMP END DO
!$OMP END PARALLEL

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^End of SBMPR^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

