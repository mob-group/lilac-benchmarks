
c     --------------------- initv  ----------------------

      subroutine initv(prcord,qrcord,srcord,
     *                 trcord,zrcord,velocp,bondln,temtur,
     *                 jstrt,jfins,tolshk,maxshk,bdshak,
     *                 timstp,numpro,ishkit,
     *                 numcrd,oarchv,
     *                 work1,work3,work4,iseed_amh,
     *                 ires)

c     ---------------------------------------------------

c     INITV  initializes the velocities and coordinates
c            when the heat bath is turned off

c     arguments:

c        maxsiz- maximum number of protein residues (i)
c        prcord- new coordinates which satsify bond 
c                lengths (o)

c     ---------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       i_axis                  index over axis's
c       i_cat                   # of possible kinds of residues used in
c                               filtration
c       i_coord                 index over coordinate types (CA, CB...)
c       i_pro                   index over proteins in ensemble
c       i_res                   index over residues
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use globals, only: maxsiz,maxpro,maxcrd,maxcnt

      implicit none


c     passed argument declarations:

         logical bdshak

         integer numpro,jstrt,jfins,
     *           maxshk,ishkit,
     *           numcrd,oarchv,ires(maxsiz)
     
         real prcord(maxsiz,3,maxpro,maxcrd),
     *        qrcord(maxsiz,3,maxpro,maxcrd),tolshk,
     *        srcord(maxsiz,3,maxpro,maxcrd),temtur,
     *        trcord(maxsiz,3,maxpro,maxcrd),timstp,
     *        zrcord(maxsiz,3,maxpro,maxcrd),
     *        velocp(maxsiz,3,maxpro,maxcrd),
     *        bondln(maxsiz,maxcrd),work1(maxcnt),
     *        work3(maxcnt),work4(maxcnt)

         integer iseed_amh(4)



c     internal variables:

c        --- do loop indices ---

         integer i507,i512,i517

        integer i_axis, i_coord, i_pro, i_res
 
c        --- implied do loop indices ---

         integer i1

         integer numshk

         real temph,diff(2,maxsiz),dist(2)

c     required subroutines

        external gaussc,shake,shakab,shakox


        temph = 0.0
        do 600 i1=1,maxsiz
          diff(1,i1)=0.0
          diff(2,i1)=0.0
600        continue
        dist(1)=0.0
        dist(2)=0.0

c     --------------------- begin -----------------------



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find gaussianly distributed velocities 
      call gaussc(maxsiz,numpro,velocp,
     *            sqrt(temtur),jstrt,jfins,maxpro,
     *            iseed_amh,work1,maxcrd,numcrd)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set initial postions according to velocites

      do 502 i_pro=1,numpro                ! number of protiens in ensemble
         do 500 i_coord=1,numcrd         ! number of coord. types (CA, CB, ...)
            do 504 i_axis=1,3                  ! x,y,z

               do 516 i_res=1,jstrt-1
                  trcord(i_res,i_axis,i_pro,i_coord)=
     *            prcord(i_res,i_axis,i_pro,i_coord)   
  516          continue

               do 515 i_res=jstrt,jfins
                  trcord(i_res,i_axis,i_pro,i_coord)=
     *            prcord(i_res,i_axis,i_pro,i_coord) - timstp*
     *            velocp(i_res,i_axis,i_pro,i_coord)
  515          continue

               do 529 i_res=1,jfins
                  qrcord(i_res,i_axis,i_pro,i_coord)=
     *            prcord(i_res,i_axis,i_pro,i_coord)
  529          continue

  504       continue
  500    continue
  502 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





c     set bond constraints for each set of 
c     coordinates

      numshk=0

  400 continue        ! coming back from 'goto 400' later in subroutine

      numshk=numshk + 1
      if( numshk.gt.10 )then
         write(oarchv,130)numshk
  130    format(/'Shkdrv: convergence not obtained')
         stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call shake(maxsiz,prcord,trcord,srcord,
     *           bondln,numpro,jstrt+1,jfins,zrcord,
     *           maxshk,tolshk,bdshak,ishkit,maxpro,
     *           maxcrd,maxcnt,work1,work3,work4,oarchv)
      if( bdshak )return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call shakab(maxsiz,prcord,trcord,srcord,
     *            bondln,numpro,jstrt,jfins,zrcord,
     *            maxshk,tolshk,bdshak,ishkit,maxpro,
     *            maxcrd,maxcnt,work1,work3,work4,oarchv,
     *            ires)
      if( bdshak )return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call shakox(trcord,prcord,srcord,
     *              numpro,jstrt,jfins,
     *              maxshk,tolshk,ishkit,
     *              maxpro,maxcrd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     determine which if any of the constraints
c     are not satisfied

      do 512 i512=1,numpro
            do 507 i507=jstrt+1,jfins

               dist(1)=sqrt
     *       ( ( prcord(i507,1,i512,1) -
     *           prcord(i507-1,1,i512,1) )**2
     *       + ( prcord(i507,2,i512,1) -
     *           prcord(i507-1,2,i512,1) )**2
     *       + ( prcord(i507,3,i512,1) -
     *           prcord(i507-1,3,i512,1) )**2 )

               dist(2)=sqrt
     *       ( ( prcord(i507,1,i512,2) -
     *           prcord(i507,1,i512,1) )**2
     *       + ( prcord(i507,2,i512,2) -
     *           prcord(i507,2,i512,1) )**2
     *       + ( prcord(i507,3,i512,2) -
     *           prcord(i507,3,i512,1) )**2 )
            diff(1,i507)=abs(dist(1) - bondln(i507,1))
            diff(2,i507)=abs(dist(2) - bondln(i507,2))
  507       continue

      do 517 i517=jstrt+1,jfins

            if (ires(i517).eq.8) diff(2,i517)=0.0   ! don't skake CB of glycine

            if( max(diff(1,i517),diff(2,i517)).gt.tolshk )then
               go to 400
            endif
  517       continue
  512 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     reset velocities in accord w/ new 
c     coordinates

      temph=1.0/timstp
      do 522 i_pro=1,numpro
         do 503 i_coord=1,numcrd
            do 523 i_axis=1,3
               do 524 i_res=jstrt,jfins
                  velocp(i_res,i_axis,i_pro,i_coord)=temph*(
     *            prcord(i_res,i_axis,i_pro,i_coord) -
     *            trcord(i_res,i_axis,i_pro,i_coord)) 
  524          continue
  523       continue
  503    continue
  522 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     ---------------------- done -----------------------

      return
      end
