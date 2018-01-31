      subroutine walesamh_initial()

      use amhglobals,  only:SO, nmtemp,itgrd,temgrd,temtur,ictemp,ctemp,
     *  iscltab,nmres,oarchv,nmstep,numpro,idigns,maxpro,maxcrd,prcord,ires,oconv,
     *  omovi,omoviseg,quench,nquench,quench_crd

      implicit none

c     subroutines required by main program

       external gentab,initil,intstr,scltab,zero_amh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     internal variables:

         integer jstrt,jfins,i_quench,len,ishkit,nmdifv

         character*10 save_name

      call zero_amh

c     --------------------- begin -----------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     open required files, read input parameter file, and generate header file
c      call read_input_alt() ! called BEFORE initil : Johan
      call initil
c      call read_altgamma()  ! called AFTER  initil : Johan
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up temperature-annealing schedule
c      call annsch(nmtemp,itgrd,temgrd,temtur,ictemp,ctemp)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate requisite force/potential tables
      call gentab
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     scale tables
      if (iscltab) write(SO,*) 'in scltab'
      if (iscltab) call scltab
      if (iscltab) write(SO,*) 'out scltab'

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate initial structures
      quench_crd=0.0D0
    
c      write(6,*) 'in intstr'
      call intstr
c      write(6,*) 'out intstr'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        if (.not. quench) nquench=1
        do i_quench = 1,nquench
        prcord=quench_crd(:,:,:,:,i_quench)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set indicies for the first and last residues
c     which are not fixed in crystal conformation
      jstrt=1
      jfins=nmres

c     set subsegment length
      len=jfins - jstrt + 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c++++++++++++++++++++++++++++++johan
!       call read_input_alt()
!       call read_altgamma() ! must be called after initil

c------------------------------johan
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- diagnostics ---
c      write(oarchv,121)jstrt,jfins,len
c  121 format(/' start ',i3,' end ',i3,' length ',i3)
c      write(oarchv,122)jstrt,jfins,nmstep,numpro
c  122 format('jstrt ',i3,' jfins ',i3,' mutations/T ',
c     *       i3,' numpro ',i3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate initial ensemble of proteins
c     find configuration which satisfies the constraints

      idigns=.false.

        enddo ! i_quench

      end
