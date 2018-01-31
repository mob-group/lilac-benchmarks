MODULE NODES
     USE COMMONS,ONLY : SLURMT, PBST, CHECKSPT, DEBUG, GETMINFRQST, GETTSFRQST
     USE PORFUNCS
     IMPLICIT NONE
     SAVE

     INTEGER :: JPN,NNODES,NNLEN,MYSTAT,ACTUALLENPATH, SLURMVERSION
     INTEGER, PARAMETER :: LENPATH=7000
     CHARACTER(LEN=80),ALLOCATABLE,DIMENSION(:) :: NODENAME
     CHARACTER(LEN=80) :: USERNAME, RANSTRING
     CHARACTER(LEN=100) :: WORKINGDIRECTORY
     CHARACTER(LEN=7000) :: TOTALJOBSTRING
!    CHARACTER(LEN=1000) :: KILLSTRING
     CHARACTER(LEN=200) :: PATHSTRING
     LOGICAL :: SSHPARALLEL=.FALSE., YESNO

     CHARACTER(LEN=120) :: NODESTRING
     CHARACTER(LEN=120) :: NODEN
     CHARACTER(LEN=120) :: TEMPSTRING
     CHARACTER(LEN=120) :: HOSTNAME
     CHARACTER(LEN=LENPATH) :: SSHPATH
     INTEGER NSTART, NFINISH, NSTART2, NFINISH2, J1, N1, N2, LCOUNT

     CONTAINS

     SUBROUTINE GETNODES(NCPU)
          USE PORFUNCS
          IMPLICIT NONE
     
          INTEGER,INTENT(OUT) :: NCPU

          INTEGER :: I
          INTEGER,PARAMETER :: UNIT1=91
          CHARACTER(LEN=80) :: ARG

          WRITE(*,'(A,I1,A)') 'GetNodes> Will run ',JPN,' jobs per core'
          INQUIRE(FILE='nodes.info',EXIST=YESNO)
          IF (YESNO.AND.(PBST.OR.SLURMT)) THEN
             OPEN(UNIT=1,FILE='nodes.info',STATUS='OLD')
             IF (.FALSE.) THEN
                READ(1,'(A)') NODESTRING
                READ(1,*) NNODES
                IF (NNODES.EQ.1) THEN ! no list to pick apart!
                   NNLEN=JPN*NNODES
                   IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                   ALLOCATE(NODENAME(NNLEN))
                   NCPU=JPN*NNODES
                   WRITE(*,'(A,I6)') 'GetNodes> Number of simultaneous OPTIM jobs=',NCPU
                   SSHPARALLEL=.TRUE.   
                   WRITE(*,'(A,I6,A)') 'GetNodes> Following ',NNODES,' cores are available:'
                   PRINT '(A)', TRIM(ADJUSTL(NODESTRING))
                   NODENAME(1:NCPU)= TRIM(ADJUSTL(NODESTRING))
                ELSE
                   LCOUNT=0
                   DO 
                      LCOUNT=LCOUNT+1
                      IF (NODESTRING(LCOUNT:LCOUNT).EQ.'[') THEN
                         NODEN=NODESTRING(1:LCOUNT-1)
                         EXIT
                      ENDIF
                   ENDDO 
                   NNLEN=JPN*NNODES
                   IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                   ALLOCATE(NODENAME(NNLEN))
                   NCPU=JPN*NNODES
                   WRITE(*,'(A,I2)') 'GetNodes> Number of simultaneous OPTIM jobs=',NCPU
                   SSHPARALLEL=.TRUE.   
                   WRITE(*,'(A,I2,A)') 'GetNodes> Following ',NNODES,' cores are available:'
                   NSTART=LCOUNT+1
                   I=0
                   outer: DO 
                      LCOUNT=LCOUNT+1
                      IF (NODESTRING(LCOUNT:LCOUNT).EQ.'-') THEN
                         NFINISH=LCOUNT-1
                         NSTART2=LCOUNT+1
                         DO 
                            LCOUNT=LCOUNT+1
                            IF ((NODESTRING(LCOUNT:LCOUNT).EQ.',').OR.(NODESTRING(LCOUNT:LCOUNT).EQ.']')) THEN
                               NFINISH2=LCOUNT-1
   
                               READ(NODESTRING(NSTART:NFINISH),*) N1
                               READ(NODESTRING(NSTART2:NFINISH2),*) N2
   !                           PRINT '(A,6I8)','NSTART,NFINISH,NSTART2,NFINISH2,N1,N2=',NSTART,NFINISH,NSTART2,NFINISH2,N1,N2
                               
                               DO J1=N1,N2
                                  I=I+1
                                  WRITE(TEMPSTRING,'(I10)') J1
                                  NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // TRIM(ADJUSTL(TEMPSTRING))
                                  PRINT '(A)', TRIM(ADJUSTL(NODEN))  // TRIM(ADJUSTL(TEMPSTRING))
                               ENDDO
   
                               NSTART=LCOUNT+1
                               IF (NODESTRING(LCOUNT:LCOUNT).EQ.']') EXIT outer
                               EXIT 
                            ENDIF
                         ENDDO
                      ELSEIF (NODESTRING(LCOUNT:LCOUNT).EQ.',') THEN
                          NFINISH=LCOUNT-1
                          I=I+1
                          NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          PRINT '(A)', TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          NSTART=LCOUNT+1
                      ELSEIF (NODESTRING(LCOUNT:LCOUNT).EQ.']') THEN
                          NFINISH=LCOUNT-1
                          I=I+1
                          NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          PRINT '(A)', TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                         EXIT outer
                      ENDIF
                   ENDDO outer
                ENDIF
                READ(1,'(A)') USERNAME
                READ(1,'(A)') WORKINGDIRECTORY
                WRITE(*,'(2A)') 'GetNodes> Working in directory ',TRIM(ADJUSTL(WORKINGDIRECTORY))
                CLOSE(1)
                PRINT '(A,I8,A)','GetNodes> Complete list of cores to be used for ',NCPU,' jobs:'
                PRINT '(A)',NODENAME(1:NCPU)
             ELSE
                READ(1,*) NNODES
                NNLEN=JPN*NNODES
                IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                ALLOCATE(NODENAME(NNLEN))
                NCPU=JPN*NNODES
                WRITE(*,'(A,I2)') 'GetNodes> Number of simultaneous OPTIM jobs=',NCPU
                IF (.NOT.SLURMT) THEN
                   SSHPARALLEL=.TRUE.
                END IF
                WRITE(*,'(a,i2,a)') 'GetNodes> Following ',NNODES,' nodes are available:'
                DO I=1,NNODES
                     READ(1,'(a)') arg
                     PRINT '(A)', TRIM(ADJUSTL(arg))
                     NODENAME(JPN*(I-1)+1:JPN*I)=ADJUSTL(ARG)
                ENDDO
                READ(1,'(A)') USERNAME
                READ(1,'(A)') WORKINGDIRECTORY
                WRITE(*,'(2A)') 'GetNodes> Working in directory ',TRIM(ADJUSTL(WORKINGDIRECTORY))
                CLOSE(1)
             ENDIF
          ELSE
             IF (.NOT.PBST) THEN
                PRINT '(A)','getnodes> Interactive run - not checking for nodes.info file'
             ELSE
                PRINT '(A)','getnodes> No nodes.info file - assuming no OPTIM jobs required for PBS run'
             ENDIF
          ENDIF

          ! Find hostname for the node that PATHSAMPLE is running on - necessary for scp
          CALL MYSYSTEM(MYSTAT,DEBUG,'hostname > currentnode')
          INQUIRE(FILE='currentnode',EXIST=YESNO)
          IF (YESNO) THEN
             OPEN(UNIT=125,FILE='currentnode',STATUS='OLD')
             READ(125,'(A)') HOSTNAME
             CLOSE(125)
             CALL MYSYSTEM(MYSTAT,DEBUG,'rm currentnode')
          ELSE
             WRITE(*,'(A)') 'getnodes> Current node unknown - stopping'
             STOP
          END IF

          ! Read in current PATH so that this can be passed when submitting jobs using ssh
          IF (.NOT.SLURMT) THEN
              CALL MYSYSTEM(MYSTAT,DEBUG,'echo $PATH > sshpath')
              INQUIRE(FILE='sshpath',EXIST=YESNO)
              IF (YESNO) THEN
                 CALL MYSYSTEM(MYSTAT,DEBUG,'wc -c < sshpath >lenpath')
                 OPEN(UNIT=126,FILE='lenpath',STATUS='OLD')
                 READ(126,*) ACTUALLENPATH
                 CLOSE(126)
                 IF (ACTUALLENPATH > LENPATH) THEN
                    WRITE(*,'(A)') 'getnodes> Length of $PATH exceeds static allocation size'
                    WRITE(*,'(A)') ' - increase size of LENPATH and possibly TOTALJOBSTRING'
                    STOP
                 END IF
                 OPEN(UNIT=127,FILE='sshpath',STATUS='OLD')
                 READ(127,'(A)') SSHPATH
                 CLOSE(127)
                 CALL MYSYSTEM(MYSTAT,DEBUG,'rm sshpath lenpath')
              ELSE
                 WRITE(*,'(A)') 'getnodes> Unable to read in $PATH - stopping'
                 STOP
              END IF
          END IF

          IF (SLURMT) THEN
             CALL MYSYSTEM(MYSTAT,DEBUG,'sinfo --version | cut -c7-8 > slurmversion')
             IF (YESNO) THEN
                 OPEN(UNIT=128,FILE='slurmversion',STATUS='OLD')
                 READ(128,*) SLURMVERSION
                 CLOSE(128)
                 CALL MYSYSTEM(MYSTAT,DEBUG,'rm slurmversion')
             ELSE
                WRITE(*,'(A)') 'getnodes> Slurm version unknown - stopping'
                STOP
             END IF
          END IF
     END SUBROUTINE GETNODES

     SUBROUTINE SSHSUBMIT(ICPU,STAT,JOBSTRING,CONNSTR1,LDEBUG)
          USE PORFUNCS, ONLY: SYSTEM_SUBR
          USE COMMONS, ONLY: CHARMMT, ZSYM, COPYFILES, COPYOPTIMT, BHINTERPT, BISECTT, SSHT
          IMPLICIT NONE
          REAL HARVEST

          INTEGER,INTENT(IN) :: ICPU
          LOGICAL,INTENT(IN) :: LDEBUG
          INTEGER,INTENT(OUT) :: STAT
          CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1

          INTEGER :: MYSTATUS
          CHARACTER(LEN=80) :: NODE
          CHARACTER(LEN=100) :: TEMPCOPYFILESDIR

          NODE=TRIM(ADJUSTL(NODENAME(ICPU)))
!
!  There is a problem here with the return status. We cannot get the return status
!  from the OPTIM job - instead we get the return status from rsh, which will be 0 unless
!  something in the job submission actually fails. If we have a path.info.connstr1 file
!  then we can analyse it, so this rsh status is actually the important one!
!
!  Job submission now changed to use a single system call and one large job string.
!  Putting all the data transfer etc. in the original rsh, so that it runs on
!  the compute node, should avoid any other rcp, rsync or rsh !
!
!  Incredibly, we occasionally see problems with jobs interfering with each other because
!  the same process id is chosen. Add a random number to the temporary directory name...
!
          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER(HARVEST)
          WRITE(RANSTRING,'(F10.5)') HARVEST

          PATHSTRING='/scratch/' // TRIM(ADJUSTL(USERNAME)) // '/' // CONNSTR1 // '.' // TRIM(ADJUSTL(RANSTRING))

          ! Put COPYFILES into temporary directory so contents can be copied as a group to avoid chopping up COPYFILES string
          ! *.connstr1 also put into temporary directory to avoid unnecessarily long jobstring
          TEMPCOPYFILESDIR = connstr1 // '-tempcopyfiles'
          CALL MYSYSTEM(MYSTAT,LDEBUG,'mkdir ' // TRIM(ADJUSTL(TEMPCOPYFILESDIR)))
          CALL MYSYSTEM(MYSTAT,LDEBUG,'cp *.' // connstr1 // ' ' // TRIM(ADJUSTL(COPYFILES)) // ' ' // TRIM(ADJUSTL(TEMPCOPYFILESDIR)))

!         TEMPSTRING=" sed -e 's/^/kill -9 /' poo2 > poo3 ; "

! csw34> KILLSTRING created a killfile which contained the process ids of the OPTIM jobs spawned by PATHSAMPLE.
!        This was necessary as some of the clusters did not clean up stray OPTIM jobs which were running if 
!        PATHSAMPLE died. This has now been fixed. The code will remain here in case it is needed in the future.
!          KILLSTRING='echo "ssh ' // TRIM(NODE) // ' quote ps -f | grep ' // TRIM(ADJUSTL(CONNSTR1)) // &
!  &                  ' | grep -v ssh | grep -v bash | grep -v PATH | grep -v grep > poo1 ; ' // &
!  &                  ' cut -c10-15 poo1 > poo2 ; ' // TEMPSTRING // &
!  &                  ' chmod +x poo3 ; ./poo3 quote"  >> killfile'
!
! Build up the complete rsh command step by step:
! (1) make the scratch directory on the node. 
!     -p flag means no error is generated if the directory already exists.
!     -x disables X forwarding.
          IF (SSHT) THEN
             TOTALJOBSTRING= 'ssh -x ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING)) 
          ELSE
             TOTALJOBSTRING= 'rsh ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING)) 
          ENDIF
! (2) copy data from WORKINGDIRECTORY to the scratch directory on the node
!     - scp used instead of copy as cannot rely on the contents of NFS-mounted directories to look the same on both nodes 
!     Note that if any file is missing an error condition will result, and subsequent commands will fail.
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
  &          // ' ; scp -r ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY)) // '/' &
  &          // TRIM(ADJUSTL(TEMPCOPYFILESDIR)) // '/* ' // TRIM(ADJUSTL(PATHSTRING)) // ' ; rm -rf ' & 
  &          // TRIM(ADJUSTL(WORKINGDIRECTORY)) // '/' // TRIM(ADJUSTL(TEMPCOPYFILESDIR))
! (3) move to the scratch directory on the node
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(PATHSTRING))
! (3b) delete any existing path.info.* file (a very rare but not impossible condition!)
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && rm -f path.info.* '
! (4) run the OPTIM job
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; PATH=' // TRIM(ADJUSTL(SSHPATH)) // ' ' // JOBSTRING
! (5) copy results back
          IF (LDEBUG) THEN ! copy everything back 
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; scp *.' // connstr1 &
   &                      // ' ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY))
          ELSEIF (COPYOPTIMT.AND.(BHINTERPT.OR.BISECTT)) THEN ! copy path.info, OPTIM
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp OPTIM* min.data.info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':'// TRIM(ADJUSTL(WORKINGDIRECTORY))
          ELSEIF (COPYOPTIMT) THEN ! copy path.info, OPTIM, odata and finish
!            IF (CHECKSPT) THEN
!               TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
!  &            // ' ; scp OPTIM* ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY))
                IF (GETMINFRQST.OR.GETTSFRQST) THEN !need to copy frqs.dump
                    TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &                // ' ; scp frqs.* ' // TRIM(ADJUSTL(HOSTNAME)) // ':'// TRIM(ADJUSTL(WORKINGDIRECTORY))
                END IF
!            ELSE
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp OPTIM* *info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY))
!            ENDIF
          ELSE ! we only really need path.info or min.data.info
!            IF (BHINTERPT.OR.BISECTT) THEN
!               TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
!  &            // ' ; scp min.data.info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY))
!            ELSE
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp *info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':' // TRIM(ADJUSTL(WORKINGDIRECTORY))
!            ENDIF
          ENDIF
! (6) remove the scratch directory
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; rm -r ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
!     or don;t rm it for debugging
!         TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ' // ' " '
!         TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ls ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
          IF (LDEBUG) PRINT '(2A)', 'nodes> complete job string: ',TRIM(ADJUSTL(TOTALJOBSTRING)) 
! (7) submit the job for real
!         CALL SYSTEM_SUBR(TRIM(ADJUSTL(KILLSTRING)),MYSTATUS)  
          CALL SYSTEM_SUBR(TRIM(ADJUSTL(TOTALJOBSTRING)),MYSTATUS)  
          STAT=MYSTATUS

     END SUBROUTINE SSHSUBMIT

     SUBROUTINE SLURMSUBMIT(STAT,JOBSTRING,CONNSTR1,LDEBUG)
          USE PORFUNCS, ONLY: SYSTEM_SUBR
          USE COMMONS, ONLY: COPYFILES, COPYOPTIMT, BHINTERPT, BISECTT, CUDAT
          IMPLICIT NONE
          REAL HARVEST

          LOGICAL,INTENT(IN) :: LDEBUG
          INTEGER,INTENT(OUT) :: STAT
          CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1

          INTEGER :: MYSTATUS
          CHARACTER(LEN=100) :: TEMPCOPYFILESDIR

!
!  Incredibly, we occasionally see problems with jobs interfering with each other because
!  the same process id is chosen. Add a random number to the temporary directory name...
!
          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER(HARVEST)
          WRITE(RANSTRING,'(F10.5)') HARVEST

          PATHSTRING='/scratch/$USER/' // CONNSTR1 // '.' // TRIM(ADJUSTL(RANSTRING))

! Put COPYFILES into temporary directory so contents can be copied as a group to avoid chopping up COPYFILES string
! *.connstr1 also put into temporary directory to avoid unnecessarily long jobstring
          TEMPCOPYFILESDIR = connstr1 // '-tempcopyfiles'
          CALL MYSYSTEM(MYSTAT,LDEBUG,'mkdir ' // TRIM(ADJUSTL(TEMPCOPYFILESDIR)))
          CALL MYSYSTEM(MYSTAT,LDEBUG,'cp *.' // connstr1 // ' ' // TRIM(ADJUSTL(COPYFILES)) // ' ' // TRIM(ADJUSTL(TEMPCOPYFILESDIR)))

! (1) make the scratch directory on the node.
!     -p flag means no error is generated if the directory already exists.
          TOTALJOBSTRING= 'mkdir -p ' // TRIM(ADJUSTL(PATHSTRING))
! (2) copy data from the working directory to the scratch directory on the node
!     - scp used instead of copy as cannot rely on the contents of NFS-mounted directories to look the same on both nodes
!     Note that if any file is missing an error condition will result, and subsequent commands will fail.
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; scp -r ' // TRIM(ADJUSTL(HOSTNAME)) & 
  &          // ':$SLURM_SUBMIT_DIR/' // TRIM(ADJUSTL(TEMPCOPYFILESDIR)) // '/* ' // TRIM(ADJUSTL(PATHSTRING)) & 
  &          // ' ; rm -rf $SLURM_SUBMIT_DIR/' // TRIM(ADJUSTL(TEMPCOPYFILESDIR))
! (3) move to the scratch directory on the node
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(PATHSTRING))
! (3b) delete any existing path.info.* file (a very rare but not impossible condition!)
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && rm -f path.info.* '
! (4) run the OPTIM job
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ;' // JOBSTRING
! (5) copy results back
          IF (LDEBUG) THEN ! copy everything back
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; scp *.' // connstr1 &
   &                      // ' ' // TRIM(ADJUSTL(HOSTNAME)) // ':$SLURM_SUBMIT_DIR'
          ELSEIF (COPYOPTIMT.AND.(BHINTERPT.OR.BISECTT)) THEN ! copy path.info, OPTIM
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp OPTIM* min.data.info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':$SLURM_SUBMIT_DIR'
          ELSEIF (COPYOPTIMT) THEN ! copy path.info, OPTIM, odata and finish
                IF (GETMINFRQST.OR.GETTSFRQST) THEN !need to copy frqs.dump
                    TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &                // ' ; scp frqs.* ' // TRIM(ADJUSTL(HOSTNAME)) // ':$SLURM_SUBMIT_DIR'
                END IF
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp OPTIM* *info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':$SLURM_SUBMIT_DIR'
          ELSE ! we only really need path.info or min.data.info
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &            // ' ; scp *info* ' // TRIM(ADJUSTL(HOSTNAME)) // ':$SLURM_SUBMIT_DIR'
          ENDIF
! (6) remove the scratch directory
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; rm -r ' // TRIM(ADJUSTL(PATHSTRING))
          IF (LDEBUG) PRINT '(2A)', 'nodes> complete job string for srun command: ',TRIM(ADJUSTL(TOTALJOBSTRING))
! (7) submit the job for real
! Commands for srun written into this script
          OPEN(UNIT=225,FILE='submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh',STATUS='UNKNOWN')
          WRITE(225,'(A)') '#!/bin/bash'
          WRITE(225,'(A)') TRIM(ADJUSTL(TOTALJOBSTRING))
          CLOSE(225)
          CALL MYSYSTEM(MYSTAT,LDEBUG,'chmod +x submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh')
          IF (SLURMVERSION <= 14) THEN
              IF (CUDAT) THEN
                 CALL SYSTEM_SUBR('srun -n1 -N1 --exclusive --gres=gpu:1 submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh',MYSTATUS)
              ELSE
                 CALL SYSTEM_SUBR('srun -n1 -N1 --exclusive submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh',MYSTATUS)
              END IF
          ELSE IF (SLURMVERSION >= 15) THEN
              IF (CUDAT) THEN
                 CALL SYSTEM_SUBR('srun -n1 -N1 -l --gres=gpu:1 submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh',MYSTATUS)
              ELSE
                 CALL SYSTEM_SUBR('srun -n1 -N1 -l submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh',MYSTATUS)
              END IF
          ELSE
              WRITE(*,'(A)') 'slurmsubmit> Slurm version not recognised - stopping'
              STOP 
          END IF
! Remove srun script after completion
          CALL MYSYSTEM(MYSTAT,LDEBUG,'rm -rf submit_' // TRIM(ADJUSTL(RANSTRING)) // '.sh')
          STAT=MYSTATUS

     END SUBROUTINE SLURMSUBMIT

END MODULE NODES
