
c     --------------------- intstr ----------------------

      subroutine restrt(procnt,AMHmaxsiz,nmres,maxpro,
     *                  numpro,maxcrd,numcrd,prcord,
     *                  oarchv)

c     ---------------------------------------------------

c     RESTRT reads in previously calculated structures

c     ---------------------------------------------------

        use amhglobals,  only:SO, quench,nquench,quench_crd,temtur,
     *               temtur_quench,itgrd,x_mcp,ires

        use key, ONLY : FILTH2

      implicit none

c     argument declarations:

         integer procnt,AMHmaxsiz,nmres,maxpro,
     *           maxcrd,numcrd,oarchv,numpro,i602
         CHARACTER(LEN=20) :: OTEMP 
         CHARACTER(LEN=20) :: OSTRING
         CHARACTER(LEN=2) :: SDUMMY
       

         double precision prcord(AMHmaxsiz,3,maxpro,maxcrd),x,y,z

c     internal variables:

         integer numprr,nmrss,idummy,gly_c

c        --- do loop indices ---

         integer i500,i502,i503
 
c        --- implied do loop indices ---

         integer i2,iii

c        ----  new annsch struff
         integer grid_num,inc,i501,i504

         double precision rinc,q_temp(200)

         rinc = 0.0
         numprr = 1 
c     --------------------- begin -----------------------

c     read in previous proteins 
c     q_temp records quenching temperature
c     open file containing coordinates

c     the structure file used to be t.dat 
c     I've changed so that format is compatible
c     with a segment of a movie file

       IF (FILTH2.NE.0) THEN
          WRITE(OTEMP,*) FILTH2
          WRITE(OSTRING,'(A)') 'start.' // TRIM(ADJUSTL(OTEMP))
       ELSE
          WRITE(OSTRING,'(A)') 'start'
       ENDIF
!      PRINT '(A,I8,A)','FILTH2,OSTRING=',FILTH2,OSTRING

      open(unit=80,file=OSTRING,status='old',form='formatted')

c     read in previous structures

c      read(80,67)nmrss,numcrr,numprr,nquench
      write(SO,*) 'quench, nquench=',quench,nquench
cc   67 format(3(i3,1x))
   67 format(4(i8,1x))

       if (numprr.gt.numpro) then
          write(SO,*) 'too many structures in movieseg file'
          write(SO,*) numprr,numpro
          stop
       endif

       if (nquench.gt.200 .and. quench) then
          write(SO,*) 'too many structures to quench',quench
          stop
       endif
  
        if (.not. quench) nquench = 1

c     check that the number of specified
c     residues is equal to those in the input
c     data file

c      if( nmrss.ne.nmres )then
c         write(oarchv,455)nmrss,nmres
c455      format(/'Restrt: restart -- # residues not consistent ',2(i4,1x))
c         stop
c      endif

c     number of atoms/residue consistent?

c      if( numcrr.ne.numcrd )then
c         write(oarchv,456)numcrr,numcrd
c 456    format(/'Restrt: restart -- numcrd not consistent ',2(i4,1x))
c      endif

c     read in trial structures for non-quench run

      if (.not. quench) then
              
             i602=1
            do 500 i500=1,nmres
                  READ(80,*)X,Y,Z
c                  write(6,77) X,Y,Z
77                format(3(G25.15))
                  prcord(i500,1,i602,1)=x
                  prcord(i500,2,i602,1)=y
                  prcord(i500,3,i602,1)=z

                  READ(80,*)X,Y,Z
c                  write(6,77)X,Y,Z
                  prcord(i500,1,i602,2)=x
                  prcord(i500,2,i602,2)=y
                  prcord(i500,3,i602,2)=z

                  READ(80,*)X,Y,Z
c                  WRITE(6,77) X,Y,Z
                  prcord(i500,1,i602,3)=x
                  prcord(i500,2,i602,3)=y
                  prcord(i500,3,i602,3)=z

500          continue
 
       gly_c = 0

       do iii = 1,nmres

       if (ires(iii).eq.8) then

        x_mcp(9*(iii-1)+1-(gly_c)*3) = dble(prcord(iii, 1, 1, 1))   !  CA X
        x_mcp(9*(iii-1)+2-(gly_c)*3) = dble(prcord(iii, 2, 1, 1))   !  CA Y
        x_mcp(9*(iii-1)+3-(gly_c)*3) = dble(prcord(iii, 3, 1, 1))   !  CA Z
c       x_mcp(9*(iii-1)+4) = (prcord(iii, 1, 1, 2))   !  CB X
c       x_mcp(9*(iii-1)+5) = (prcord(iii, 2, 1, 2))   !  CB Y
c       x_mcp(9*(iii-1)+6) = (prcord(iii, 3, 1, 2))   !  CB Z
        x_mcp(9*(iii-1)+4-(gly_c)*3) = dble(prcord(iii, 1, 1, 3))   !   O X
        x_mcp(9*(iii-1)+5-(gly_c)*3) = dble(prcord(iii, 2, 1, 3))   !   O Y
        x_mcp(9*(iii-1)+6-(gly_c)*3) = dble(prcord(iii, 3, 1, 3))   !   O Z
        gly_c = gly_c +1
      else
        x_mcp(9*(iii-1)+1-gly_c*3) = dble(prcord(iii, 1, 1, 1))   !  CA X
        x_mcp(9*(iii-1)+2-gly_c*3) = dble(prcord(iii, 2, 1, 1))   !  CA Y
        x_mcp(9*(iii-1)+3-gly_c*3) = dble(prcord(iii, 3, 1, 1))   !  CA Z
        x_mcp(9*(iii-1)+4-gly_c*3) = dble(prcord(iii, 1, 1, 2))   !  CB X
        x_mcp(9*(iii-1)+5-gly_c*3) = dble(prcord(iii, 2, 1, 2))   !  CB Y
        x_mcp(9*(iii-1)+6-gly_c*3) = dble(prcord(iii, 3, 1, 2))   !  CB Z
        x_mcp(9*(iii-1)+7-gly_c*3) = dble(prcord(iii, 1, 1, 3))   !   O X
        x_mcp(9*(iii-1)+8-gly_c*3) = dble(prcord(iii, 2, 1, 3))   !   O Y
        x_mcp(9*(iii-1)+9-gly_c*3) = dble(prcord(iii, 3, 1, 3))   !   O Z
      endif

      enddo

c        do i503=1, 189*3-15
c             write(6,*)' x_mcp restart ',  x_mcp(i503), i503
c        enddo

        quench_crd(:,:,:,:,1)=prcord
      
        else
          do i503=1, nquench
           do i502=1,numprr
             read(80,683)idummy,idummy,idummy,q_temp(i503),idummy
683          format(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')
c             write(6,*)'quench_t  temtur i503' , q_temp(i503)

               do i500=1,nmrss
c             write(SO,*) i500,nmrss,numprr,nquench,i503,temtur(i503)

                  read(80,*)
     *            (quench_crd(i500,i2,i502,1,i503),i2=1,3),
     *            (quench_crd(i500,i2,i502,2,i503),i2=1,3),
     *            (quench_crd(i500,i2,i502,3,i503),i2=1,3)
              enddo
            enddo
       enddo  ! i503 (loop over quench)
         do  i504=1, nquench
c        set up new annealling schedule  
         inc=0
         do 540 grid_num=1,4
            if( itgrd(grid_num).gt.0 )then
               rinc=q_temp(i504)
               rinc=rinc/float(itgrd(grid_num))
                do 501 i501=1,itgrd(grid_num)
                  temtur_quench(inc+i501,i504)=q_temp(i504) -
     *                             rinc*float(i501-1)
c                  write(6,*)'temtur_quench steps structure',
c     *                 temtur_quench(inc+i501,i504),inc+i501,i504
  501          continue
               inc=inc + itgrd(grid_num)
               endif  !  itgrd(grid_num)
540          continue  !   do 540 grid_num=1,4
             enddo   !  do i504=1, nquench
        endif ! quench

c     set number if initial structures to number of
c     structures just read in

      procnt=numprr

c     send message acknowledging that 'old' structures
c     are being used as starting structures

ctemphack
c      write(oarchv,450)numprr
c  450 format(/'restart with ',i3,' proteins')

       close(80)

c     ---------------------- done -----------------------

      return
      end
