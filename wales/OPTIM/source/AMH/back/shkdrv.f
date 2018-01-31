
c     --------------------- shkdrv  ----------------------

      subroutine shkdrv(prcord,qrcord,srcord,
     *                  zrcord,bondln,jstrt,jfins,tolshk,
     *                  maxshk,bdshak,numpro,ishkit,
     *                  maxpro,maxcrd,numcrd,
     *                  oarchv,work1,work3,work4,
     *                  ires)

c     ---------------------------------------------------

c     SHKDRV driver routine for performing shake
c            algorithm on different sets of
c            coordinates

c     arguments:

c        maxsiz- maximum number of protein residues (i)
c        prcord- new coordinates which satsify bond 
c                lengths (o)

c     ---------------------------------------------------

      use globals, only:SO, maxsiz,maxcnt

      implicit none


c     argument declarations:

         logical bdshak

         integer numpro,jstrt,jfins,
     *           maxshk,ishkit,maxpro,
     *           maxcrd,numcrd,oarchv,
     *           ires(maxsiz),i1
     
         real prcord(maxsiz,3,maxpro,maxcrd),
     *        qrcord(maxsiz,3,maxpro,maxcrd),
     *        srcord(maxsiz,3,maxpro,maxcrd),
     *        zrcord(maxsiz,3,maxpro,maxcrd),
     *        bondln(maxsiz,maxcrd),work1(maxcnt),
     *        work3(maxcnt),work4(maxcnt),tolshk

c     internal variables:

c        --- do loop indices ---

         integer i507,i512,i517
 
         integer numshk

         
         real diff(2,maxsiz),dist(2)

        character*3 res_type


c        --- implied do loop indices ---

        integer iii,jjj








c     required subroutines

        external shake,shakab,shakox


          do 600 i1=1,maxsiz
            diff(1,i1)=0.0
            diff(2,i1)=0.0
600          continue
          dist(1)=0.0
          dist(2)=0.0


c     --------------------- begin -----------------------


c     attempt to satisfy bond constraints 
c     for each set of coordinates

c     set number of shake iterations to zero

      numshk=0

c     'do while' loop, i.e., while both sets of 
c     constraints are not satisfied keep shaking

  400 continue

      numshk=numshk + 1
c      write(6,*)'numshk ishkdrv ',  numshk

c     if the number of iterations have exceeded the 
c     maximum, then print message and abort

      if( numshk.gt.100 )then
         write(oarchv,130)numshk
  130    format(/'Shkdrv: You got problems bud ',
     *          '-- convergence not obtained w/ shake')
         stop
      endif

c     shake alpha-alpha coordinates

      call shake(maxsiz,prcord,qrcord,srcord,
     *           bondln,numpro,jstrt+1,jfins,zrcord,
     *           maxshk,tolshk,bdshak,ishkit,maxpro,
     *           maxcrd,maxcnt,work1,work3,work4,oarchv)

c     if convergence not obtained, hang it up;
c     also if only one coordinate type return

      if( bdshak.or.(numcrd.eq.1) )return

c     shake alpha-beta coordinates

      call shakab(maxsiz,prcord,qrcord,srcord,
     *            bondln,numpro,jstrt,jfins,zrcord,
     *            maxshk,tolshk,bdshak,ishkit,maxpro,
     *            maxcrd,maxcnt,work1,work3,work4,oarchv,
     *            ires)

c     convergence not obtained, hang it up

          if( bdshak ) then
              write(SO,*) "Yep it died in shakab"          
          endif


      if( bdshak )return

         call shakox(prcord,qrcord,srcord,
     *               numpro,jstrt,jfins,
     *               maxshk,tolshk,ishkit,
     *               maxpro,maxcrd)

c     determine which, if any, of the alpha-alpha 
c     constraints are not satisfied

c     find bond distance

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
            if (ires(i517).eq.8) diff(2,i517)=0.0
            if( max(diff(1,i517),diff(2,i517)).gt.tolshk )then
               go to 400
            endif
  517       continue
  512 continue



c     mission accomplished
           
c     ---------------------- done -----------------------

      return
      end
