
c     --------------------- shakab -----------------------

      subroutine shakab(maxsiz,prcord,qrcord,srcord,
     *                  bondln,numpro,jstrt,jfins,zrcord,
     *                  maxshk,tolshk,bdshak,ishkit,
     *                  maxpro,maxcrd,maxcnt,work1,work3,
     *                  work4,oarchv,ires)

c     ---------------------------------------------------

c     SHAKAB attempts to satisfy alpha-beta distance 
c            constraints using the iterative procedure 
c            proposed by Ryckaert et al., J. Comp. Phys. 
c            23, 327-341 (1977).

c     arguments:

c        maxsiz- maximum number of protein residues (i)
c        prcord- updated coordinates satisfying bond
c                constraints (i,o)
c        qrcord- previous coordinates (i)
c        srcord- scratch space (i)
c        bondln- specified bond lengths (i)
c        numpro- number of configurations (i)
c        jstrt - first nonfixed site (i)
c        jfins - last nonfixed site (i)
c        zrcord- work array (i)
c        maxshk- maximum number of shake iterations (i)
c        tolshk- shake tolerance (i)
c        bdshak- flag to alert driver to a ineffective
c                shake (o)
c        ishkit- flag used to track number of shake 
c                iterations (i)
c        maxpro- maximum number of proteins (i)
c        maxcrd- maximum number of coordinates (i)

c     ---------------------------------------------------

      implicit none

c     argument declarations:

         integer maxsiz,jstrt,jfins,numpro,
     *           maxshk,ishkit,maxpro,maxcrd,maxcnt,
     *           oarchv,ires(maxsiz)

         real prcord(maxsiz,3,maxpro,maxcrd),
     *        qrcord(maxsiz,3,maxpro,maxcrd),
     *        srcord(maxsiz,3,maxpro),
     *        zrcord(maxsiz,3,maxpro),
     *        bondln(maxsiz,maxcrd),tolshk,
     *        work1(maxcnt),work3(maxcnt),
     *        work4(maxcnt)
     
c     internal variables:

         logical lbad,bdshak

         character*3 res_type

 
c        --- implied do loop indices ---

        integer i_axis, i_pro, i_res, i_shk_iter,iii,jjj

         real diff

c     required subroutines

c     --------------------- begin -----------------------

c     --- diagnostics ---

c      if( idigns )then
c         
cc        echo scalar argument arguments
c

cc         check that there is room for storage of r 
cc         in srcord(1,.,.)
c
c          if( jstrt.le.0 )then
c             write(oarchv,106)jstrt
c  106        format(/'shake:jstrt too big ',i3)
c             stop
c          endif

c        print out current and previous bond lengths

c         write(oarchv,105)
c  105    format(/'shakab:prcord qrcord')
c        do 539 i_res=jstrt,jfins
c           work1(1)=(prcord(i_res,1,1,1)-
c    *                prcord(i_res,1,1,2))**2 +
c    *               (prcord(i_res,2,1,1)-
c    *                prcord(i_res,2,1,2))**2 +
c    *               (prcord(i_res,3,1,1)-
c    *                prcord(i_res,3,1,2))**2 
c           work1(1)=sqrt(work1(1))

c           work1(2)=(qrcord(i_res,1,1,1)-
c    *                qrcord(i_res,1,1,2))**2 +
c    *               (qrcord(i_res,2,1,1)-
c    *                qrcord(i_res,2,1,2))**2 +
c    *               (qrcord(i_res,3,1,1)-
c    *                qrcord(i_res,3,1,2))**2 
c           work1(2)=sqrt(work1(2))

c            write(oarchv,104)i_res,
c     *                      (prcord(i_res,i1,1,1),i1=1,3),
c     *                      (prcord(i_res,i1,1,2),i1=1,3)
c  104       format(i3,2(3(1x,1pe10.3),2x))

c            diff=abs(bondln(i_res,2)-work1(1))
c            if( (diff.gt.1.0) )
c     *         write(oarchv,191)i_res,work1(1),work1(2),
c     *                          bondln(i_res,2),diff
c  191          format('front shake ',i3,1x,4(1pe12.5,1x))

c  539    continue

cc        check if next-nearest neighbor constraints are 
cc        satisfied; if not, then print out offending sites
c
c         call dist2(maxsiz,jstrt,jfins,prcord,diff)
c
c      endif

c     --- end diagnostics ---

c     find the r of Ryckaert et al.

      do 521 i_pro=1,numpro
         do 520 i_axis=1,3

c           compute r

            do 513 i_res=jstrt,jfins
               if( ires(i_res).ne.8 )then
                 zrcord(i_res,i_axis,i_pro)=qrcord(i_res,i_axis,i_pro,2)
     *                                 - qrcord(i_res,i_axis,i_pro,1)
               else
                  zrcord(i_res,i_axis,i_pro)=0.0
               endif
  513       continue

  520    continue
  521 continue


        
c     --- diagnostics ---

c     echo r

c       write(oarchv,811)
c 811   format(/'Shakab:  protein: zrcord')
c      do 522 i_res=jstrt,jfins
c         write(30,135)i_res,(zrcord(i_res,i1,1),i1=1,3)
c  13     format(i3,1x,3(1pe10.3,1x))
c  52  continue

c     set array work4 to the (bond lengths)**2

      do 530 i_res=jstrt,jfins
         work4(i_res)=bondln(i_res,2)**2
  530 continue

c     loop over numpro protein configurations

      do 501 i_pro=1,numpro

c        set beta-glycine to alpha-glycine coordinates

         do 534 i_res=jstrt,jfins
            if( ires(i_res).eq.8 )then
               prcord(i_res,1,i_pro,2)=prcord(i_res,1,i_pro,1)
               prcord(i_res,2,i_pro,2)=prcord(i_res,2,i_pro,1)
               prcord(i_res,3,i_pro,2)=prcord(i_res,3,i_pro,1)
            endif
  534    continue

c        perform maxshk iterations in an attempt to satisfy the 
c         nearest-neighbor distance constraints
 
         do 500 i_shk_iter=1,maxshk

c           perform checkerboard breakup, i.e., analyze
c           constraints independently of one another

c           find r' of Ryckaert et al.

            do 503 i_axis=1,3
               do 502 i_res=jstrt,jfins
                  srcord(i_res,i_axis,i_pro)=
     *            prcord(i_res,i_axis,i_pro,2) -
     *            prcord(i_res,i_axis,i_pro,1)
  502          continue
  503       continue

c           find |r'|**2

            do 504 i_res=jstrt,jfins
               work1(i_res)=srcord(i_res,1,i_pro)**2 + 
     *                     srcord(i_res,2,i_pro)**2 +
     *                     srcord(i_res,3,i_pro)**2
  504       continue

c           find d**2 - r'**2

            do 505 i_res=jstrt,jfins
               work1(i_res)=work4(i_res) - work1(i_res)
  505       continue

c           check if done

            do 506 i_res=jstrt,jfins

c              if all constraints not satisfied, then continue;
c              otherwise, consider next configuration

               if( (abs(work1(i_res)).gt.tolshk) )then
                     go to 522
               endif

  506       continue

c           done

c           increment variable used to track the number
c           of required iterations

            ishkit=ishkit + i_shk_iter - 1
            go to 501

  522       continue

c           find r.r'

            do 507 i_res=jstrt,jfins

               work3(i_res)=zrcord(i_res,1,i_pro)*
     *                     srcord(i_res,1,i_pro) +
     *                     zrcord(i_res,2,i_pro)*
     *                     srcord(i_res,2,i_pro) +
     *                     zrcord(i_res,3,i_pro)*
     *                     srcord(i_res,3,i_pro)

  507       continue


cccccccccccccccccccccccccccccccccccccccccccccccccc
c           --- diagnostics ---
c           print work3 if too small
c            do 561 i_res=jstrt,jfins
c               if( abs(work3(i_res)).lt.epsiln )then
c                  diff=0.25*work1(i_res)/work3(i_res)
c                  write(oarchv,124)i_shk_iter,i_res,work3(i_res)
c                  write(oarchv,711)
c     *            (zrcord(i_res,i1,i_pro),i1=1,3),
c     *            (srcord(i_res,i1,i_pro),i1=1,3)
c  711             format('z + s',6(1x,1pe10.3))
c  124             format('shake iter ',i4,' site ',i3,
c     *                   ' work3 ',1pe10.3,' g ',1pe10.3)
c               endif
c  561       continue
c           --- end diagnostics ---
cccccccccccccccccccccccccccccccccccccccccccccccccc


c           compute g 

            do 509 i_res=jstrt,jfins
               if( ires(i_res).ne.8 )then
                  work3(i_res)=0.25*work1(i_res)/work3(i_res)
               else
                  work3(i_res)=0.0
               endif
  509       continue 



c           update coordinates

            do 511 i_axis=1,3

               do 510 i_res=jstrt,jfins
                  prcord(i_res,i_axis,i_pro,2)=
     *            prcord(i_res,i_axis,i_pro,2) +
     *            zrcord(i_res,i_axis,i_pro)*work3(i_res)
  510          continue

               do 524 i_res=jstrt,jfins
                  prcord(i_res,i_axis,i_pro,1)=
     *            prcord(i_res,i_axis,i_pro,1) -
     *            zrcord(i_res,i_axis,i_pro)*work3(i_res)
  524          continue
 
  511       continue

  500    continue

  501 continue

      continue

c     --- diagnostics ---

      lbad=.false.

c     determine which if any of the constraints
c     are not satisfied

      do 512 i_pro=1,numpro
         do 514 i_axis=1,3
            do 523 i_res=jstrt,jfins
               srcord(i_res,i_axis,i_pro)=prcord(i_res,i_axis,i_pro,2)
     *                              - prcord(i_res,i_axis,i_pro,1)
  523       continue
  514    continue
         do 515 i_res=jstrt,jfins
            work1(i_res)=srcord(i_res,1,i_pro)**2 +
     *                  srcord(i_res,2,i_pro)**2 +
     *                  srcord(i_res,3,i_pro)**2
  515    continue
         do 516 i_res=jstrt,jfins
            work1(i_res)=sqrt(work1(i_res))
  516    continue


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do 517 i_res=jstrt,jfins 
            diff=abs(work1(i_res) - bondln(i_res,2))
            if( diff.gt.tolshk )then
               lbad=.true.
               go to 519
            endif
  517    continue
  519    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           
         if( lbad )then


c           --- all constraints not satisfied 
c               after maxshk iterations ---

            write(oarchv,103)
  103       format(/'not all bonds satisfied in shakab')

            do 518 i_res=jstrt,jfins 
               diff=abs(work1(i_res) - bondln(i_res,2))
               if( diff.gt.tolshk )then
                  write(oarchv,110)i_pro,i_res,work1(i_res),
     *                             bondln(i_res,2),diff
  110             format('pro ',i2,' site ',i3,' d(cal) ',
     *                   1pe10.3,' d(exact) ',1pe10.3,
     *                   ' diff ',1pe10.3)
               endif
  518       continue
   
            bdshak=.true.
            do 533 i_axis=1,3
               do 525 i_res=jstrt,jfins
                  prcord(i_res,i_axis,i_pro,2)=
     *            qrcord(i_res,i_axis,i_pro,2)
  525             continue
  533       continue
         else

         endif
  512 continue
 
c     --- end diagnostics ---

c     ---------------------- done ----------------------

 
      return
      end
