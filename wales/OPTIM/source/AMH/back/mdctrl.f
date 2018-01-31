 
c     --------------------- mdctrl ----------------------
 
      subroutine mdctrl(jstrt,jfins,tempav,ishkit,
     *                  bdshak,tmpidx,
     *                  totke
     *                  )

c     ---------------------------------------------------

c     MDCTRL is the master subroutine for the molecular
c            dynamics algorithm outlined in Ryckaert and 
c            Ciccotti, Mol. Phys. 58, 1125-1136 (1986). 

c     arguments:

c        jstrt - first nonfixed site (i)
c        jfins - last nonfixed site (i)
c        tempav- if tempav, then compute averages (i)
c        maxshk- maximum number of shake iterations (i)
c        tolshk- shake tolerance (i)
c        ishkit- tracks the number of shake iterations (o)
c        bdshk - true if shake doesn't converge within
c                specified number of iterations (o)
c        tmpidx- temperature index (i)

c     ---------------------------------------------------

      use globals, only: maxtab,numcrd,numpro,qrcord,work1,
     *  velocp,prcord,zrcord,avepp,maxr,vpotnt,
     *  forse,work4,pexcld,
     *  numlng,nmres,rincinv,rincsq,oarchv,ilong,crdixn,
     *  ires,work3,timstp,eqdist,hbond,oxexcldv,
     *  i540,i1,
     *  avep,movanal,trcord,srcord,bondln,tolshk,maxshk,maxpro,
     *  maxcrd,i511,i518


      implicit none

c     argument declarations:

         logical tempav,bdshak

         integer jstrt,jfins,ishkit,tmpidx

         real totke(maxpro)
         

c     internal variables:

         logical scl_call

         integer i_axis, i_coord, i_pro, i_res

         real starhi,trgeng(maxtab,3)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- Every Step Part 1 of 3 ---
c  needed for printing out every step in a trajectory
C
c         character*10 name_a(10000)
c         character*3 res_type(maxsiz)
c         character*10 save_name
c         integer name_length, nl
c         integer pdb, pdb_num, pn
c         external get_res_name
c         DATA pdb_num /0/
ccccccccccccccccccccccccccccccccccccccccccccccccc


c     required subroutines

         external force,verlet,
     *            shkdrv


c     --------------------- begin -----------------------

!         write(SO,*) 'in mdctrl'

c     save current coordinates

      do 500 i_coord=1, numcrd             ! CA, CA, O
         do 524 i_pro=1,numpro
            do 525 i_axis=1,3
               do 526 i_res=1,jfins
                  qrcord(i_res,i_axis,i_pro,i_coord)=
     *            prcord(i_res,i_axis,i_pro,i_coord)
  526          continue
  525       continue
  524    continue
  500 continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find forces and potential energies 
c       write(SO,*) 'force    '

      scl_call=.false.

      call force(numpro,
     *            prcord,
     *            zrcord,avepp,tempav,
     *            maxr,vpotnt,
     *            forse,trgeng,
     *            pexcld,
     *            numlng,nmres,rincinv,rincsq,
     *            ilong,
     *            crdixn,
     *            ires,
     *            eqdist,
     *            hbond,oxexcldv,
     *            scl_call)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if( tempav ) then

c           record total E for T average

           do 540 i540=1,numpro
             do 100 i1=1,50
               avep(i540,1,i1)=avepp(i540,1,i1)
100             continue
  540      continue

         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update coordinates

      if (.not.movanal) then

!       write(SO,*) 'in verlet'

      call verlet(numpro,prcord,zrcord,
     *            timstp,jstrt,jfins,trcord)

c     find set of final coordinates consistent with
c     bond constraints

            call shkdrv(prcord,qrcord,srcord,
     *            zrcord,bondln,jstrt,jfins,tolshk,
     *            maxshk,bdshak,numpro,ishkit,
     *            maxpro,maxcrd,numcrd,
     *            oarchv,work1,work3,work4,ires)

      endif

c     if bdshak, then shake was ineffective; return

      if( bdshak )then
         write(oarchv,144)
  144    format(/'bad shake post verlet')
         return
      endif
c     set velocity 

      starhi=0.5/timstp

      do 514 i_coord=1,numcrd
         do 509 i_pro=1,numpro
            do 510 i_axis=1,3
                  do 511 i511=jstrt,jfins
                     velocp(i511,i_axis,i_pro,i_coord)=
     *              (prcord(i511,i_axis,i_pro,i_coord) - 
     *               trcord(i511,i_axis,i_pro,i_coord))*starhi
  511             continue
               do 518 i518=jstrt,jfins
                  trcord(i518,i_axis,i_pro,i_coord)=
     *            qrcord(i518,i_axis,i_pro,i_coord)
  518          continue
  510       continue
  509    continue
  514 continue

      if (tempav) then
        totke=0.0
        do i_coord=1,numcrd
        do i_pro=1,numpro
        do i_axis=1,3
        do i511=jstrt,jfins
          totke(i_pro)=totke(i_pro)+velocp(i511,i_axis,i_pro,i_coord)**2
        enddo
        enddo
        enddo
        enddo
        totke=0.5*totke
      endif

c     ---------------------- done -----------------------

      return
      end
