
c     --------------------- shake -----------------------

      subroutine shake(maxsiz,prcord,qrcord,srcord,
     *                 bondln,numpro,jstrt,jfins,zrcord,
     *                 maxshk,tolshk,bdshak,ishkit,
     *                 maxpro,maxcrd,maxcnt,work1,work3,
     *                 work4,oarchv)

c     ---------------------------------------------------

c     SHAKE attempts to satisfy nearest-neighbor distance
c           constraints using iterative procedure 
c           proposed by Ryckaert et al., J. Comp. Phys. 
c           23, 327-341 (1977).

c     arguments:

c        maxsiz- maximum number of protein residues
c        prcord- updated coordinates satisfying bond
c                constraints
c        qrcord- previous coordinates
c        srcord- scratch space
c        bondln- specified bond lengths
c        numpro- number of configurations
c        jstrt - first nonfixed site
c        jfins - last nonfixed site
c        zrcord- work array
c        maxshk- maximum number of shake iterations
c        tolshk- shake tolerance
c        bdshak- flag to alert driver to a ineffective
c                shake
c        ishkit- flag used to track number of shake 
c                iterations
c        maxpro- maximum number of proteins

c     ---------------------------------------------------

      implicit none

c     argument declarations:

         integer maxsiz,jstrt,jfins,numpro,
     *           maxshk,ishkit,maxpro,maxcnt,oarchv,
     *           maxcrd

         real prcord(maxsiz,3,maxpro,maxcrd),
     *        qrcord(maxsiz,3,maxpro,maxcrd),
     *        srcord(maxsiz,3,maxpro),
     *        zrcord(maxsiz,3,maxpro),
     *        bondln(maxsiz,maxcrd),tolshk,
     *        work1(maxcnt),work3(maxcnt),
     *        work4(maxcnt)
    
         character*3 res_type
 
c     internal variables:

         logical lbad,bdshak

         integer istrt

c        --- do loop indices ---


        integer i_axis, i_max_shk_iter, i_odd_even, 
     *          i_pro, i_pro_2, i_res,i1,iii,jjj


c        --- implied do loop indices ---

         real diff

c     required subroutines

c     --------------------- begin -----------------------

c     --- diagnostics ---

c        echo scalar argument arguments

c         write(oarchv,200)maxcrd,maxcnt
c  200    format('maxcrd ',i3,' maxcnt ',i5)
c
ccc         check that there is room for storage of r 
ccc         in srcord(1,.,.)
cc
cc          if( jstrt.le.1 )then
cc             write(oarchv,106)jstrt
cc  106        format(/'shake:jstrt too big ',i3)
cc             stop
cc          endif

c        print out current and previous bond lengths

c         write(oarchv,105)
c  105    format(/'shake:prcord qrcord')
c         do 539 i_res=jstrt,jfins
c            work1(1)=(prcord(i_res,1,1,1)-
c     *                prcord(i_res-1,1,1,1))**2 +
c     *               (prcord(i_res,2,1,1)-
c     *                prcord(i_res-1,2,1,1))**2 +
c     *               (prcord(i_res,3,1,1)-
c     *                prcord(i_res-1,3,1,1))**2 
c            work1(1)=sqrt(work1(1))

c            work1(2)=(qrcord(i_res,1,1,1)-
c     *                qrcord(i_res-1,1,1,1))**2 +
c     *               (qrcord(i_res,2,1,1)-
c     *                qrcord(i_res-1,2,1,1))**2 +
c     *               (qrcord(i_res,3,1,1)-
c     *                qrcord(i_res-1,3,1,1))**2 
c            work1(2)=sqrt(work1(2))

c            write(oarchv,104)i_res,
c     *                      (prcord(i_res,i1,1,1),i1=1,3),
c     *                      (qrcord(i_res,i1,1,1),i1=1,3)
c  104       format(i3,2(3(1x,1pe10.3),2x))

c            diff=abs(bondln(i_res,1)-work1(1))
c            if( (diff.gt.1.0).or.(diff.le.5.0e-05) )
c     *         write(oarchv,191)i_res,work1(1),work1(2),
c     *                          bondln(i_res,1),diff
c  191          format('front shake ',i3,1x,4(1pe12.5,1x))

c  539    continue

c      endif
c     --- end diagnostics ---
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find the r of Ryckaert et al.
      do 521 i_pro=1,numpro
         do 520 i_axis=1,3

c           compute r

            do 513 i_res=jstrt,jfins
               zrcord(i_res,i_axis,i_pro)=
     *            qrcord(i_res,i_axis,i_pro,1) -
     *            qrcord(i_res-1,i_axis,i_pro,1)
  513       continue
  520    continue
  521 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     --- diagnostics ---

c      if( idigns )then

c        echo r

c          write(oarchv,*)i_res,' protein: srcord'
c         do 522 i_res=jstrt,jfins
c            write(30,135)i_res,(srcord(i_res,i1),i1=1,3)
c  135       format(i3,1x,3(1pe10.3,1x))
c  522    continue
c
c      endif

c     set array work4 to the (bond lengths)**2

      do 530 i_res=jstrt,jfins
         work4(i_res)=bondln(i_res,1)**2
  530 continue

c     loop over numpro protein configurations

      do 501 i_pro=1,numpro

c        perform maxshk iterations in an attempt 
c        to satisfy the nearest-neighbor distance 
c        constraints

         do 500 i_max_shk_iter=1,maxshk

c           perform checkerboard breakup, i.e., analyze
c           constraints independently of one another



            do 508 i_odd_even=1,2

               if( i_odd_even.eq.1 )then

c                 'odd' constraints

                  istrt=jstrt

               else

c                 'even' constraints

                  istrt=jstrt + 1

               endif





c              find r' of Ryckaert et al.

               do 503 i_axis=1,3
                  do 502 i_res=jstrt,jfins
                     srcord(i_res,i_axis,i_pro)=
     *               prcord(i_res,i_axis,i_pro,1) -
     *               prcord(i_res-1,i_axis,i_pro,1)
  502             continue
  503          continue

c              find |r'|**2

               do 504 i_res=jstrt,jfins
                     work1(i_res)=srcord(i_res,1,i_pro)**2 + 
     *                           srcord(i_res,2,i_pro)**2 +
     *                           srcord(i_res,3,i_pro)**2
  504          continue

c              find d**2 - r'**2

               do 505 i_res=jstrt,jfins
c                  write(30,402)i_res,work1(i_res),work4(i_res),
c     *                         abs(work4(i_res)-work1(i_res))
c  402             format(i3,3(1x,1pe10.3))
                  work1(i_res)=work4(i_res) - work1(i_res)
  505          continue

c              check if done


               do 506 i_res=jstrt,jfins

c                 if all constraints not satisfied, then continue;
c                 otherwise, consider next configuration

                  if( abs(work1(i_res)).gt.tolshk )then
                     go to 522
                  endif

  506          continue

c              done

c              increment variable used to track the number
c              of required iterations

               ishkit=ishkit + i_max_shk_iter - 1
               go to 501

  522          continue

c              find r.r'
              do 507 i_res=istrt,jfins,2

              work3(i_res)=zrcord(i_res,1,i_pro)*srcord(i_res,1,i_pro)
     *                   + zrcord(i_res,2,i_pro)*srcord(i_res,2,i_pro)
     *                   + zrcord(i_res,3,i_pro)*srcord(i_res,3,i_pro)

  507          continue

c              --- diagnostics ---

c               if( idigns )then
c 
cc                 print work3 if too small
c
c                  do 561 i_res=istrt,jfins,2
c                     if( abs(work3(i_res)).lt.epsiln )then
cc                        diff=0.25*work1(i_res)/work3(i_res)
c                        write(oarchv,124)i_max_shk_iter,i_res,work3(i_res)
c                        write(oarchv,711)
c     *                  (zrcord(i_res,i1,i_pro),i1=1,3),
c     *                  (srcord(i_res,i1,i_pro),i1=1,3)
c  711                   format('z + s',6(1x,1pe10.3))
c  124                   format('shake iter ',i4,' site ',i3,
c     *                         ' work3 ',1pe10.3,' g ',1pe10.3)
c                     endif
c  561             continue
c
c               endif

c              --- end diagnostics ---

c              compute g 

               do 509 i_res=istrt,jfins,2
                  work3(i_res)=0.25*work1(i_res)/work3(i_res)
  509          continue 

c              update coordinates

               do 511 i_axis=1,3

                  if( (istrt.gt.2).and.(i_odd_even.eq.1) )then

                     do 510 i_res=istrt+1,jfins-1,2
                        prcord(i_res,i_axis,i_pro,1)=
     *                  prcord(i_res,i_axis,i_pro,1) -
     *                  zrcord(i_res+1,i_axis,i_pro)*work3(i_res+1)
  510                continue

                  else

                     do 537 i_res=istrt-1,jfins-1,2
                        prcord(i_res,i_axis,i_pro,1)=
     *                  prcord(i_res,i_axis,i_pro,1) -
     *                  zrcord(i_res+1,i_axis,i_pro)*work3(i_res+1)
  537                continue

                  endif

                  do 524 i_res=istrt,jfins,2
                     prcord(i_res,i_axis,i_pro,1)=
     *               prcord(i_res,i_axis,i_pro,1) +
     *               zrcord(i_res,i_axis,i_pro)*work3(i_res)
  524             continue
 
                  if( (istrt.gt.2).and.(i_odd_even.eq.1) )then
                     prcord(istrt,i_axis,i_pro,1)=
     *               prcord(istrt,i_axis,i_pro,1) +
     *               zrcord(istrt,i_axis,i_pro)*work3(istrt)
                  endif

  511          continue

c               if( i_odd_even.eq.1 )then 
c                  write(30,918)i_max_shk_iter,istrt-1,
c     *            (prcord(istrt-1,i1,i_pro,1),
c     *                    work8(i1),i1=1,3)
c  918             format(2(i3,1x),6(1pe10.3,1x))
c               endif

  508       continue

  500    continue

  501 continue

      continue



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- diagnostics ---

c      if( idigns )then

         lbad=.false.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        determine which if any of the constraints
c        are not satisfied

         do 512 i_pro_2=1,min(numpro,2)
            do 514 i_axis=1,3
               do 523 i_res=jstrt,jfins
                  srcord(i_res,i_axis,i_pro_2)=
     *               prcord(i_res,i_axis,i_pro_2,1) -
     *               prcord(i_res-1,i_axis,i_pro_2,1)
  523          continue
  514       continue



            do 515 i_res=jstrt,jfins
               work1(i_res)=srcord(i_res,1,i_pro_2)**2 +
     *                     srcord(i_res,2,i_pro_2)**2 +
     *                     srcord(i_res,3,i_pro_2)**2
  515       continue

            do 516 i_res=jstrt,jfins
               work1(i_res)=sqrt(work1(i_res))
  516       continue

            do 517 i_res=jstrt,jfins 
               diff=abs(work1(i_res) - bondln(i_res,1))
               if( diff.gt.tolshk )then
                  lbad=.true.
                  go to 519
               endif
  517       continue


  519       continue
           
            if( lbad )then

c              --- all constraints not satisfied 
c                  after maxshk iterations ---

               write(oarchv,103)
  103          format(/'not all bonds satisfied in shake')

               do 518 i_res=jstrt,jfins 
                  diff=abs(work1(i_res) - bondln(i_res,1))
                  if( diff.gt.tolshk )then
                     write(oarchv,110)i_pro_2,i_res,work1(i_res),
     *                                bondln(i_res,1),diff
  110                format('pro ',i2,' site ',i3,' d(cal) ',
     *                      1pe10.3,' d(exact) ',1pe10.3,
     *                      ' diff ',1pe10.3)
                  endif
  518          continue
   
               bdshak=.true.
               do 533 i_axis=1,3
                  do 525 i_res=jstrt,jfins
                     prcord(i_res,i_axis,i_pro_2,1)=
     *               qrcord(i_res,i_axis,i_pro_2,1)
  525             continue
  533          continue
            else

c               write(30,*)'shake ok ',i_max_shk_iter-1

            endif
  512    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c      endif

c     --- end diagnostics ---

c     ---------------------- done ----------------------
     
      return
      end
