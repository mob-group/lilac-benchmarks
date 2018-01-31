      PROGRAM GENCOORDS

!Computes the position and orientation of nC clusters each containing nE ellipsoids. 
!The parameters that define the arrangement are contained in a .model file.
!A specific function must be written to interpret those parameters.

!Given a set of input parameters each function must return:
!1. a Status indicating success or failure 0=good 1= fail.
!2. The position and rotation matrix for each cluster in the cluster chain defined in the lab frame. 'CMatrices'
!3. The position and rotation matrix for each ellipsoid in the 'zero cluster' i.e. the cluster at zero position and zero rotation in the lab frame. 'EMatrices'

!The program then replicates and transforms each ellipsoid in the cluster definition to each point in the cluster chain to generate the overall sequence of ellipsoids. 'AMatrices'

!The positions, size and angle axis of each ellipsoid within the zero cluster is output into a pysites.xyz file. 
!The position and angle axis of each cluster in the cluster chain is output into a coords file.
!The positions, size, rotation matrix and angle axis of each ellipsoid is output into the file ellipsoid.model.xyz 

!It is assumed that the first cluster definition is given in the same frame as the cluster chain, 
!and it is up to the user to ensure their function orients the cluster correctly in the zero position. It is used 'as is'.

      IMPLICIT NONE 

!declare variables
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARAMS
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CPOS, CAA
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ESIZE, EPOS, EAA
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ASIZE, APOS, AAA
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: ARMAT, CRMAT, ERMAT
      DOUBLE PRECISION :: FMAT(1:3,1:3), FPOS(1:3),newBasis(1:3,1:3),aMult, rMult
      INTEGER :: modelSelected=0, iErr=0, curE=1, curP=1, curC = 1, curA=1, nP=0, nC=0, nE=0, nA=0
      CHARACTER*100::modelType=''

!open the input and output files
      OPEN(4, FILE='ellipsoid.model',STATUS='old')
      OPEN(5, FILE='pysites.xyz',STATUS='replace')
      OPEN(8, FILE='ellipsoid.model.xyz',STATUS='replace')
      OPEN(7, FILE='coords',STATUS='replace')

!read the model file 
!first up is the modelType, number of parameters: nP, the number of clusters: nC and the number of ellipses per cluster: nE
      read(4,*) modelType, nP, nC, nE, aMult, rMult
      write(6,*,advance='no') 'Creating coords, pysites.xyz and ellipsoid.model.xyz'
      write(6,*,advance='no') 'for model: ', modelType,' with ',nC,' clusters and ', nE,' ellipsoids per cluster'

!compute total number ellipsoids
      nA=nC*nE
!allocate memory for input parameters and the output arrays
      allocate(PARAMS(nP))
      allocate(CPOS(nC,3))
      allocate(CAA(nC,3))
      allocate(CRMAT(nC,3,3)) 
      allocate(ESIZE(nE,6))
      allocate(EPOS(nE,3))
      allocate(EAA(nE,3))
      allocate(ERMAT(nE,3,3))
      allocate(ASIZE(nA,6))
      allocate(APOS(nA,3))
      allocate(AAA(nA,3))
      allocate(ARMAT(nA,3,3))

!Read the list of input parameters
      do curP=1,nP,1
      read(4,*) PARAMS(curP)
      write(6,*,advance='no') PARAMS(curP)
      end do

!call the function to generate the output
      modelSelected=0
  

      IF(modelType.eq.'pdb') THEN
        modelSelected=1
        iErr=0
        CALL pdb(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'coordsModel') THEN
        modelSelected=1
        iErr=0
        CALL COORDSMODEL(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 



      IF(modelType.eq.'figureGenerator') THEN
        modelSelected=1
        iErr=0
        CALL FIGGEN(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 


      IF(modelType.eq.'frenetHelix') THEN
        modelSelected=1
        iErr=0
        CALL HELIXMODEL(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'displacement') THEN
        modelSelected=1
        iErr=0
        CALL DISPLACEMENT(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'verticalHelix') THEN
        modelSelected=1
        iErr=0
        CALL VERTICALHELIX(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'distance') THEN
        modelSelected=1
        iErr=0
        CALL DISTANCE(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'angle') THEN
        modelSelected=1
        iErr=0
        CALL ANGLE(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      IF(modelType.eq.'anglerot') THEN
        modelSelected=1
        iErr=0
        CALL ANGLEROT(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 


      IF(modelType.eq.'angleNudge') THEN
        modelSelected=1
        iErr=0
        CALL ANGLENUDGE(iErr,nP,nC,nE,PARAMS,ESIZE,EPOS,ERMAT,CPOS,CRMAT)

        IF (iErr.gt.0) THEN
           write(6,*,advance='no') 'Error ', iErr, ' in model ', modelType
           STOP
        END IF
      END IF 

      if(modelSelected.eq.0) THEN
          write(6,*,advance='no') 'Unrecognised Modeltype: ', modelType
          CLOSE(4)
          CLOSE(5)
          CLOSE(7)
          CLOSE(8)
          STOP
      END IF
 
!      IF (globalRotateOn) THEN

!      CALL GENERATEBASIS(x,y,z,newBasis)

!      DO curC=1,nC,1
!           FPOS= matmul(CPOS(curC,:),newBasis(:,:))
!           CPOS(curC,:)=FPOS
!           FMAT(:,:)=matmul((newBasis(:,:)),matmul(CRMAT(curC,:,:),transpose(newBasis(:,:))))
!           CRMAT(curC,:,:)=FMAT        
!      END DO
!      END IF
 



     !Loop through the clusters and generate the full list of ellipsoids and the angle axis vectors
      curA=1
      DO curC=1,nC,1
       ! write(6,*) ''
        DO curE=1,nE,1
           !step one: construct the position of the ellipsoid using the cluster frame 
           !vectors in CRMAT and the info in EPOS which says how much of each to use.
           !data in rotation matrices is in column vectors.
           APOS(curA,:)=CPOS(curC,:)+matmul(EPOS(curE,:),transpose(CRMAT(curC,:,:)))
           !write(*,*) 'CRMAT:'
           !write(*,*) CRMAT(curC,:,:)
           write(*,*) 'EPOS:'
           write(*,*) EPOS(curE,:)
           write(*,*) 'APOS:'
           write(*,*) APOS(curA,:)
           !step two: Do the same for the orientation vectors (column vectors of ERMAT) of the ellipsoid, record them as column vectors of ARMAT
           !ARMAT(curA,:,:)=matmul(CRMAT(curC,:,:,),matmul(ERMAT(curE,:,:),(transpose(CRMAT(curC,:,:)))))
           ARMAT(curA,:,:)=matmul(CRMAT(curC,:,:),ERMAT(curE,:,:))
           write(*,*) 'ERMAT'
           write(*,*) ERMAT(curE,1,1:3)
           write(*,*) ERMAT(curE,2,1:3)
           write(*,*) ERMAT(curE,3,1:3)
           write(*,*) 'CRMAT'
           write(*,*) CRMAT(curC,1,1:3)
           write(*,*) CRMAT(curC,2,1:3)
           write(*,*) CRMAT(curC,3,1:3)
           write(*,*) 'ARMAT'
           write(*,*) ARMAT(curA,1,1:3)
           write(*,*) ARMAT(curA,2,1:3)
           write(*,*) ARMAT(curA,3,1:3)
         
           !transcribe sizes 
           ASIZE(curA,:)=ESIZE(curE,:)

           !angle axis conversions of rotation matrices
           CALL RMAT2AA(ERMAT(curE,:,:),EAA(curE,:))
           CALL RMAT2AA(ARMAT(curA,:,:),AAA(curA,:))
           CALL RMAT2AA(CRMAT(curC,:,:),CAA(curC,:))
           !write(*,*) 'EAA'
           !write(*,*) EAA(curE,:)
           !write(*,*) 'AAA'
           !write(*,*) AAA(curA,:)
           !write(*,*) 'CAA'
           !write(*,*) CAA(curC,:)
           curA=curA+1
        END DO
      END DO


!****** write the pysites file *******
      WRITE(5,*) nE
      WRITE(5,*,advance='no') modelType

!loop through nE ellipsoids
      DO curE=1,nE,1
        !output to the pysites.xyz file
        WRITE(5,100) EPOS(curE,1), EPOS(curE, 2), EPOS(curE, 3), aMult, rMult, ESIZE(curE,1)/2,ESIZE(curE,2)/2,ESIZE(curE,3)/2,ESIZE(curE,4)/2,ESIZE(curE,5)/2,ESIZE(curE,6)/2, EAA(curE,1) , EAA(curE,2), EAA(curE,3)
100   FORMAT (' O ', 3F16.9, ' ellipse ', 8F16.9, ' atom_vector ', 3F16.9)      
      END DO



!***** write the coords file *****
!Loop through nC clusters.
      DO curC=1,nC,1
!write the cluster positions to the coords file.       
      WRITE(7,102) CPOS(curC,1), CPOS(curC,2), CPOS(curC,3)
102   FORMAT (3F16.9)
      END DO

!Loop through nC clusters again, this time output orientations.
      DO curC=1,nC,1
!write the cluster orientations to the coords file.       
      WRITE(7,102) CAA(curC,1), CAA(curC,2), CAA(curC,3)
      END DO


!******** Write the ellipsoid.model.xyz  ****************
      WRITE(8,*) nA
      WRITE(8,198) 1, -999.9999999, 1
198 FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
      DO curA=1,nA,1
!output the attractive ellipsoid axes in ellipsoid.model.xyz      
      write(8,200) APOS(curA,1), APOS(curA,2), APOS(curA, 3), ASIZE(curA,1), ASIZE(curA,2),ASIZE(curA,3), &
                    ARMAT(curA,1,1), ARMAT(curA,1,2), ARMAT(curA,1,3), ARMAT(curA,2,1),ARMAT(curA,2,2),ARMAT(curA,2,3), &
                    ARMAT(curA,3,1),ARMAT(curA,3,2),ARMAT(curA,3,3), &
                    AAA(curA,1),AAA(curA,2),AAA(curA,3)
200   FORMAT ('    O',2x,3F20.10,2x,'ellipse ',12F15.8,2x,'atom_vector',3F15.8)
      END DO


!close the files
      CLOSE(4)
      CLOSE(5)
      CLOSE(7) 
      CLOSE(8)

!deallocate arrays
!allocate memory for input parameters and the output arrays
      deallocate(PARAMS)
      deallocate(CPOS)
      deallocate(CAA)
      deallocate(CRMAT)
      deallocate(ASIZE)
      deallocate(APOS)
      deallocate(AAA)
      deallocate(ARMAT)
      deallocate(EPOS)
      deallocate(EAA)
      deallocate(ERMAT)
      deallocate(ESIZE)
     END 
      INCLUDE 'aa2rmat.f90'
      INCLUDE 'rmat2aa.f90'
      INCLUDE 'coordsModel.f90'
      INCLUDE 'helixModel.f90'
      INCLUDE 'frenetHelix.f90'
      INCLUDE 'flatHelix.f90'
      INCLUDE 'displacement.f90'
      INCLUDE 'verticalHelix.f90'
      INCLUDE 'distance.f90'
      INCLUDE 'angle2Dist.f90'
      INCLUDE 'dist2Angle.f90'
      INCLUDE 'GWrap.f90'
      INCLUDE 'angle.f90'
      INCLUDE 'angleRot.f90'
      INCLUDE 'pdb.f90'
      INCLUDE 'angleNudge.f90'
      INCLUDE 'angleCluster.f90'
      INCLUDE 'rotateVec.f90'
      INCLUDE 'globalRotate.f90'
      INCLUDE 'generateBasis.f90'
      INCLUDE 'euler2Rmat.f90'
      INCLUDE 'figgen.f90'
