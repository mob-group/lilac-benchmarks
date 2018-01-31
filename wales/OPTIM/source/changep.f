C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE CHANGEP
      USE COMMONS
      USE KEY
      USE MODTWOEND
C      USE MODAMBER
      USE MODNEB
      use porfuncs
      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, NTYPEA
      COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO, CAT, YESNO
      CHARACTER WORD*16, WW*20, FNAME*19
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      INTEGER MP, LP,  XMP, XLP, ISTAT
      DOUBLE PRECISION GTOL,BSTPMIN,BSTPMAX, XGTOL,XBSTPMIN,XBSTPMAX
      COMMON /LB3/MP,LP,GTOL,BSTPMIN,BSTPMAX
      COMMON /XLB3/XMP,XLP,XGTOL,XBSTPMIN,XBSTPMAX

C     PRINT*,'In changep'
      IF (FILTH.EQ.0) THEN
         WRITE(FNAME,'(A12)') 'changeparams'
      ELSE 
         WRITE(FNAME,'(A)') 'changeparams.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

      INQUIRE(FILE=FNAME,EXIST=YESNO)
C     PRINT*,'EXIST=',YESNO
      IF (YESNO) THEN
C
C  Unit 5 seems to get left open in some paths through the program. It should
C  be closed in fetchz?
C
         CLOSE(5)
         OPEN (5,FILE=FNAME,STATUS='OLD')
      ELSE
         RETURN
      ENDIF

190   CALL INPUT(END)
C     PRINT*,'END=',END
      IF (.NOT. END) THEN
         CALL READU(WORD)
C        PRINT*,'WORD=',WORD
      ENDIF
C
C  POINTS - keyword at the end of the list of options after which
C           the Cartesian coordinates follow. Must be present unless VARIABLES, AMBER, 
C           VASP, QCHEM, CASTEP, ONETEP, CP2K, CHARMM or CPMD 
C           is present instead.
C
      IF (END) THEN
        CLOSE(5)
        CALL SYSTEM(' mv changeparams changeparams.read ')
C       PRINT*,'mv changeparams changeparams.read '
        RETURN
      ENDIF

      IF (WORD.EQ.'    ' .OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190
C
C  SEARCH specifies the value of INR, i.e. the search type.     - default n=0
C
      ELSE IF (WORD.EQ.'SEARCH') THEN
         CALL READI(INR)
         WRITE(*,'(A,I3)') '*** Changing search type to ',INR
C
C  Turn on LBFGS gradient minimization. GMAX is the convergence
C  criterion for the RMS gradient, default 0.001.
C  For BFGSTS NEVL and NEVS are the maximum iterations allowed in the searches for 
C  the largest and smallest eigenvectors, respectively and NBFGSMAX1 is the largest
C  number of BFGS steps allowed in the subsequent restricted minimization.
C  If the negative eigenvalue appears to have converged then NBFGSMAX2 steps
C  are allowed in the tangent space.
C  CONVU is used to determine convergence in such runs and BFGSCONV can be used
C  to set GMAX, the convergence criteria for the subspace optimization.
C
C  IF REOPT is true the smallest Hessian eigenvector is redetermined after the
C  EF step before the tangent space minimisation.
C
      ELSE IF (WORD.EQ.'BFGSCONV') THEN
         IF (NITEMS.GT.1) THEN
            CALL READF(GMAX)
         ENDIF
         WRITE(*,'(A,F15.5)') '*** Changing LBFGS RMS convergence criterion to ',GMAX
      ELSE IF (WORD.EQ.'BFGSTS') THEN
         BFGSTST=.TRUE.
         INR=2
         IF (NITEMS.GT.1) THEN
            CALL READI(NEVS)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READI(NBFGSMAX1)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READI(NBFGSMAX2)
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READF(CEIG)
         ENDIF
         IF (NITEMS.GT.5) THEN
            CALL READI(NEVL)
         ENDIF
         WRITE(*,'(A)') '*** Changing parameters for BFGSTS according to file changeparams'
      ELSE IF (WORD.EQ.'REOPT') THEN
         REOPT=.TRUE.
         WRITE(*,'(A)') '*** turning on eigenvector reoptimisation after EF step'
C
C  If CHECKINDEX is .TRUE. and the BFGSTS routine converges an attempt is
C  made to count the number of negative Hessian eigenvalues using projection,
C  orthogonalization and iteration. We also need the opportunity to change the
C  parameters NEVL and NEVS within BFGSTS if BFGSTS isn t true.
C  CHECKINDEX can also be used with BFGSMIN and should understand NOHESS too.
C
      ELSE IF (WORD.EQ.'CHECKINDEX') THEN
         CHECKINDEX=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NEVS)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(CEIG)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READI(NEVL)
         ENDIF
         WRITE(*,'(A)') '*** turning on Hessian index check'
C
C  If the index found by checkindex does not correspond to BFGSMIN or BFGSTS then
C  CHECKCONT causes a pushoff along the eigenvector correpsonding to the softest
C  undesired negative eigenvalue. 
C
      ELSE IF (WORD.EQ.'CHECKCONT') THEN
         CHECKCONT=.TRUE.
         WRITE(*,'(A)') '*** will continue after Hessian index check'
C
C  MAXBFGS x1 x2 x3 x4\/}: {\it x\/} specifies the maximum allowed step length in LBFGS
C  minimisations, {\it x1\/} for  normal minimisations, {\it x2\/} for Rayleigh-Ritz ratio
C  minimisation, {\it x3\/} for putting structures in closest coincidence with
C  {\bf mind} (NO LONGER USED!!), and {\it x4\/} for NEB minimisations. Default values all 0.2.
C
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)
         IF (NITEMS.GT.2) CALL READF(MAXXBFGS)
         IF (NITEMS.GT.3) CALL READF(MAXMBFGS)
         IF (NITEMS.GT.4) CALL READF(MAXNEBBFGS)
C
C  PRESSURE tells the program to perform a constant pressure optimisation
C           for SC, ME and P6 with periodic boundary conditions - default off
C
      ELSE IF (WORD.EQ.'PRESSURE') THEN
         PRESSURE=.TRUE.
         WRITE(*,'(A)') '*** turning on constant pressure optimisation'
C
C  NZERO is the number of zero eigenvalues, default 0.
C
      ELSE IF (WORD.EQ.'ZEROS') THEN
         CALL READI(NZERO)
         WRITE(*,'(A,I4)') '*** number of zero Hessian eigenvalues reset to ',NZERO
      ELSE IF (WORD.EQ.'EVCUT') THEN
         CALL READF(EVCUT)
         WRITE(*,'(A,F15.5)') '*** cutoff for zero Hessian eigenvalues reset to ',EVCUT
      ELSE IF (WORD.EQ.'PARALLEL') THEN
         PARALLEL=.TRUE.
         CALL READA(NPROC)
         WRITE(*,'(A,A)') '*** number of parallel jobs reset to ',NPROC
C
C  Double ended ts search.
C
      ELSE IF (WORD.EQ.'TWOENDS') THEN
         TWOENDS=.TRUE.
         IF (NITEMS.GT.1) CALL READF(FSTART)
         IF (NITEMS.GT.2) CALL READF(FINC)
         IF (NITEMS.GT.3) CALL READI(NTWO)
         IF (NITEMS.GT.4) CALL READF(RMSTWO)
         IF (NITEMS.GT.5) CALL READI(NTWOITER)
         IF (NITEMS.GT.6) CALL READF(TWOEVAL)
          WRITE(*,'(A,2F15.5,I4,F15.5,I4,F15.5)') '*** TWOENDS parameters reset to ',FSTART,FINC,NTWO,RMSTWO,NTWOITER,TWOEVAL
C
C  SCALE n sets the value of ISTCRT                             - default n=10
C
      ELSE IF (WORD.EQ.'SCALE') THEN
         CALL READI(ISTCRT)
         WRITE(*,'(A,I4)') '*** step scaling criterion reset to ',ISTCRT
C
C  PRINT n sets the value of IPRNT                              - default n=0
C
      ELSE IF (WORD.EQ.'PRINT') THEN
         CALL READI(IPRNT)
         WRITE(*,'(A,I4)') '*** print level reset to ',IPRNT
C
C  MODE n  specifies the eigenvector to follow                  - default n=0
C
      ELSE IF (WORD.EQ.'MODE') THEN
         CALL READI(IVEC)
         IF (NITEMS.GT.2) THEN
            CALL READI(IVEC2)
         ELSE
C           IVEC2=IVEC
         ENDIF
         WRITE(*,'(A,2I4)') '*** initial and subsequent eigenvectors followed reset to ',IVEC,IVEC2
C
C  MAXSTEP n specifies the maximum step size in real units      - default n=0.2
C
      ELSE IF (WORD.EQ.'MAXSTEP') THEN
         CALL READF(MXSTP)
         WRITE(*,'(A,F15.5)') '*** maximum step size reset to ',MXSTP
      ELSE IF (WORD.EQ.'MAXMAX') THEN
         CALL READF(MAXMAX)
         WRITE(*,'(A,F15.5)') '*** maximum value of the maximum allowed step size reset to ',MAXMAX
      ELSE IF (WORD.EQ.'MINMAX') THEN
         CALL READF(MINMAX)
         WRITE(*,'(A,F15.5)') '*** minimum value of the maximum allowed step size reset to ',MINMAX
C
C  VALUES n print the Hessian eigenvalues every n cycles        - default n=20     
C
      ELSE IF (WORD .EQ. 'VALUES') THEN
         CALL READI(NVALUES)
         WRITE(*,'(A,F15.5)') '*** Hessian eigenvalue print frequency reset to ',MINMAX
C
C  EFSTEPS n print the unscaled steps calculated for each mode
C          every n cycles                                       - default OFF
C
      ELSE IF (WORD .EQ. 'EFSTEPS') THEN
         EFSTEPST=.TRUE.
         CALL READI(EFSTEPS)
         WRITE(*,'(A,F15.5)') '*** unscaled step print frequency reset to ',EFSTEPS
C
C  STEPS n sets the number of optimisation steps to perform
C          per call to OPTIM                                    - default n=1     
C
      ELSE IF (WORD .EQ. 'STEPS') THEN
         CALL READI(NSTEPS)
         WRITE(*,'(A,I5)') '*** total number of steps reset to ',NSTEPS
C
C  DUMPVECTOR switches on dumping of eigenvectors to file 
C              vectors.dump                                     - default OFF
C  ALLSTEPS dumps the vector(s) at each step. ALLVECTORS dumps all the vectors.
C  The defaults are for only the vector corresponding to the softest non-zero
C  eigenvalue to be dumped for the last step.
C
      ELSE IF (WORD .EQ. 'DUMPVECTOR') THEN
         DUMPV=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READU(WORD)
            IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
            IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READU(WORD)
            IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
            IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
         ENDIF
         WRITE(*,'(A,F15.5)') '*** Hessian eigenvector dumping reset accoring to changeparams file '
C
C
C  GRADIENT n prints the gradients along the Hessian eigendirections
C             every n cycles                                    - default OFF
C
      ELSE IF (WORD .EQ. 'GRADIENTS') THEN
        PGRAD=.TRUE.
        CALL READI(NGRADIENTS)
         WRITE(*,'(A,F15.5)') '*** gradient printing reset according to changeparams file '
C
C  VECTORS n prints the eigenvectors every n cycles             - default OFF
C
      ELSE IF (WORD .EQ. 'VECTORS') THEN
         VECTORST=.TRUE.
         CALL READI(NVECTORS)
         WRITE(*,'(A,F15.5)') '*** Hessian eigenvector printing reset accoring to changeparams file '
C
C  SUMMARY n print a summary of the steps taken every n cycles  - default n=20   
C
      ELSE IF (WORD .EQ. 'SUMMARY') THEN
         IF (NITEMS.GT.1) CALL READI(NSUMMARY)
         WRITE(*,'(A,F15.5)') '*** summary print frequency reset accoring to changeparams file '
C
C  ADM [OFF | ON n] prints the atomic distance matrix every n 
C                   if switched on                 cycles       - default n=20      
C
      ELSE IF (WORD .EQ. 'ADM') THEN
         ADMT=.TRUE.
         CALL READI(NADM)
         WRITE(*,'(A,F15.5)') '*** atomic distance matrix print frequency reset accoring to changeparams file '
C
C  CONVERGE n m INDEX/NOINDEX sets the convergence criteria for the maximum 
C               unscaled step and RMS force                     - default n=0.0001, m=0.000001
C                                                           or m < 0.00001 .AND. n < m*100000  
C               If NOINDEX is set the Hessian index isn t checked - the default is
C               INDEX.
C
      ELSE IF (WORD .EQ. 'CONVERGE') THEN
         CALL READF(CONVU)
         IF (NITEMS.GT.2) THEN
            CALL READF(CONVR)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READU(WORD)
            IF (WORD.EQ.'NOINDEX') INDEXT=.FALSE.
         ENDIF
         WRITE(*,'(A,2F15.5)') '*** EF convergence criteria reset to ',CONVU,CONVR
C
C  SYMCUT n RMS force below which symmetry subroutine is called - default 0.001
C
      ELSE IF (WORD .EQ. 'SYMCUT') THEN
         CALL READF(SYMCUT)
         WRITE(*,'(A,2F15.5)') '*** RMS force below which symmetry is called reset to ',SYMCUT
C
C  TOLD n initial distance tolerance in symmetry subroutine     - default 0.0001
C
      ELSE IF (WORD .EQ. 'TOLD') THEN
         CALL READF(TOLD)
         WRITE(*,'(A,2F15.5)') '*** initial distance tolerance in symmetry reset to ',TOLD
C
C  TOLE n initial tolerance for the difference in principal moments 
C         of inertia divided by the sum of the principal moments 
C         in symmetry subroutine                                - default 0.0001
C
      ELSE IF (WORD .EQ. 'TOLE') THEN
         CALL READF(TOLE)
         WRITE(*,'(A,2F15.5)') '*** initial moment of inertia tolerance in symmetry reset to ',TOLE
      ELSE IF (WORD .EQ. 'AXIS') THEN
         CALL READI(NHCHECK)
         WRITE(*,'(A,I3)') '*** highest rotation axis checked in symmetry reset to ',NHCHECK
C
C  TRAD n sets the trust radius to n                            - default n=4       
C
      ELSE IF (WORD .EQ. 'TRAD') THEN
         CALL READF(TRAD)
         WRITE(*,'(A,F15.5)') '*** trust radius reset to ',TRAd
C
C  Nudged elastic band calculation using a maximum of NSTEPNEB steps with
C  NIMAGE images and RMS convergence criterion RMSNEB.
C
      ELSE IF (WORD.EQ.'NEB') THEN
         NEBT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSTEPNEB)
         IF (NITEMS.GT.2) CALL READI(NIMAGE)
         IF (NITEMS.GT.3) CALL READF(RMSNEB)
         WRITE(*,'(A,2I5,F15.5)') '*** read new parameters for old NEB ',NSTEPNEB,NIMAGE,RMSNEB

C
C  PUSHOFF x sets the magnitude of the step away from a stationary point of the
C            wrong index - default x=0.01
C
      ELSE IF (WORD .EQ. 'PUSHOFF') THEN
         CALL READF(PUSHOFF)
         WRITE(*,'(A,F15.5)') '*** pushoff reset to ',PUSHOFF
C
C  NSTEPMIN sets the minimum number of steps allowed before convergence.
C 
      ELSE IF (WORD .EQ. 'STEPMIN') THEN
         CALL READI(NSTEPMIN)
         WRITE(*,'(A,I4)') '*** minimum number of steps allowed before convergence reset to ',NSTEPMIN
C
C  PUSHCUT sets the threshold for when a PUSHOFF will be applied, i.e.
C  the RMS force must be less than PUSHCUT.
C
      ELSE IF (WORD .EQ. 'PUSHCUT') THEN
         CALL READF(PUSHCUT)
         WRITE(*,'(A,F15.5)') '*** cutoff below which pushoffs may be applied reset to ',PUSHCUT
C
C  Number of BFGS updates before resetting, default=4
C
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
         IF (NITEMS.GT.2) CALL READI(XMUPDATE)
         WRITE(*,'(A,2I4)') '*** number of LBFGS steps saved reset to ',MUPDATE,XMUPDATE
C
C  DEBUG ON/OFF sets n=1 for EFSTEPS, VALUES, SUMMARY above     - default OFF 
C
      ELSE IF (WORD .EQ. 'DEBUG') THEN
         CALL READU(WW)
         IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
           EFSTEPST=.TRUE.
           PGRAD=.TRUE.
           NGRADIENTS=1
           EFSTEPS=1
           NSUMMARY=1 
           NVALUES=1
         ENDIF
         WRITE(*,'(A,2I4)') '*** debug printing turned on '
C
C  Eigenvalue shift parameter.
C
      ELSE IF (WORD .EQ. 'SHIFT') THEN
         CALL READF(SHIFTV)
         WRITE(*,'(A,2I4)') '*** eigenvalue shift parameter changed to ',SHIFTV
C
C  Whether to put periodic images back in the primary supercell.
C
      ELSE IF (WORD .EQ. 'NORESET') THEN
         NORESET=.TRUE.
         WRITE(*,'(A,2I4)') '*** periodic images will not be put back into the primary supercell '
      ELSE IF (WORD.EQ.'NOIT') THEN
         NOIT=.TRUE.
         WRITE(*,'(A,2I4)') '*** NOIT has been set to .TRUE.'
C
C  Line minimisation of gradient along EF direction
C
      ELSE IF (WORD.EQ.'LINEMIN') THEN
         LINEMIN=.TRUE.
         WRITE(*,'(A)') '*** line minimisation for EF step turned on'
C
C  STEPS n sets the number of optimisation steps to perform
C          per call to OPTIM                                    - default n=1
C
      ELSE IF (WORD .EQ. 'BFGSSTEPS') THEN
        CALL READI(BFGSSTEPS)
         WRITE(*,'(A,I6)') '*** maximum number of BFGS steps changed to ',BFGSSTEPS
      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)
         IF (NITEMS.GT.1) CALL READF(XMAXERISE)
         WRITE(*,'(A,2G20.10)') '*** maximum energy/ev rise changed to ',MAXERISE,XMAXERISE
      ELSE IF (WORD .EQ. 'STOP') THEN
         WRITE(*,'(A)') '*** STOP directive read from changeparams'
         CALL SYSTEM(' mv changeparams changeparams.read ')
         STOP
      ELSE
        CALL REPORT('Unrecognized command '//WORD,.TRUE.)
      ENDIF

      CALL FLUSH(6)
      GOTO 190

      RETURN
      END
