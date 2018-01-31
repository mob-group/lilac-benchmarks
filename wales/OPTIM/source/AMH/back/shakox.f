
c     --------------------- shakox -----------------------

      subroutine shakox(prcord,qrcord,srcord,
     *                 numpro,jstrt,jfins,
     *                 maxshk,tolshk,ishkit,
     *                 maxpro,maxcrd)

c     ---------------------------------------------------

c     SHAKOX attempts to satisfy nearest-neighbor distance
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
c        maxshk- maximum number of shake iterations
c        tolshk- shake tolerance
c                shake
c        ishkit- flag used to track number of shake 
c                iterations
c        maxpro- maximum number of proteins
c        i_axis        - index over axes
c        i_odd_even        - index over coordinate types CA and CB
c        i_pro        - index over proteins in ensemble
c        i_res        - index over residues
c        i_shake        - index over shake iterations
c
c
c
c     skip additions
c       eq_dist_sq        equilibrium distance squared
c       
c
c     ---------------------------------------------------

      use globals, only:SO,maxsiz,maxcnt

      implicit none

c     argument declarations:

         integer jstrt,jfins,numpro,
     *           maxshk,ishkit,maxpro,
     *           maxcrd,maxpro1

         parameter (maxpro1=25)

         real prcord(maxsiz,3,maxpro,maxcrd),
     *        qrcord(maxsiz,3,maxpro,maxcrd),
     *        srcord(maxsiz,3,maxpro),
     *        tolshk
     
c     internal variables:

         character*3 res_type


         integer i_axis, i_odd_even, i_pro, i_res, i_shake,iii,jjj

         real myst(maxcnt), myst2(maxcnt)

c        --- implied do loop indices ---

         real eqdist(2),lngth(maxsiz,3,maxpro1,2)
         real eq_dist_sq(2)

c     --------------------- begin -----------------------

CCCCCCCCCCCCCCCCCCCCC


        eqdist(1)=2.42677
        eqdist(2)=2.82146
        eq_dist_sq(1)=eqdist(1)*eqdist(1)
        eq_dist_sq(2)=eqdist(2)*eqdist(2)

c     find the r of Ryckaert et al.

      do 521 i_pro=1,numpro
         do 520 i_axis=1,3

c           compute r

          do 513 i_res=jstrt,jfins-1
            lngth(i_res,i_axis,i_pro,1)=
     *          qrcord(i_res,i_axis,i_pro,3)-qrcord(i_res,i_axis,i_pro,1)
513     continue

          do 514 i_res=jstrt,jfins-1
            lngth(i_res,i_axis,i_pro,2)=
     *          qrcord(i_res+1,i_axis,i_pro,1)-qrcord(i_res,i_axis,i_pro,3)
  514     continue

  520    continue
  521 continue


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     loop over numpro protein configurations
      do 501 i_pro=1,numpro

c        perform maxshk iterations in an attempt to satisfy the 
c         nearest-neighbor distance constraints

         do 500 i_shake=1,maxshk

c           perform checkerboard breakup, i.e., analyze
c           constraints independently of one another

            do 508 i_odd_even=1,2

c              find r' of Ryckaert et al.

            do 503 i_axis=1,3
              do 502 i_res=jstrt,jfins-1
               if (i_odd_even.eq.1) then
                 srcord(i_res,i_axis,i_pro)=
     *           prcord(i_res,i_axis,i_pro,3)-prcord(i_res,i_axis,i_pro,1)
               else
                 srcord(i_res,i_axis,i_pro)=
     *           prcord(i_res+1,i_axis,i_pro,1)-prcord(i_res,i_axis,i_pro,3) 
               endif
502           continue
503         continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

c              find |r'|**2

               do 504 i_res=jstrt,jfins-1
                     myst2(i_res)=srcord(i_res,1,i_pro)**2 + 
     *                           srcord(i_res,2,i_pro)**2 +
     *                           srcord(i_res,3,i_pro)**2

c                     write(SO,*) 'myst2 ', myst2(i_res)
504          continue

c                     write(SO,*) 'myst2 ', myst2(75)

c              find d**2 - r'**2

               do 505 i_res=jstrt,jfins-1
                  myst2(i_res)=eq_dist_sq(i_odd_even) - myst2(i_res)
505          continue

c              check if done

               do 506 i_res=jstrt,jfins-1

c                 if all constraints not satisfied, then continue;
c                 otherwise, consider next configuration

                  if( abs(myst2(i_res)).gt.tolshk )then
                     go to 522
                  endif

  506          continue

c              done

c              increment variable used to track the number
c              of required iterations

               ishkit=ishkit + i_shake - 1

               go to 501  ! go to loop for the next protien

  522          continue

c              find r.r'

               do 507 i_res=jstrt,jfins-1

                  myst(i_res)=lngth(i_res,1,i_pro,i_odd_even)
     *                            *srcord(i_res,1,i_pro) +
     *                        lngth(i_res,2,i_pro,i_odd_even)
     *                            *srcord(i_res,2,i_pro) +
     *                        lngth(i_res,3,i_pro,i_odd_even)
     *                            *srcord(i_res,3,i_pro)

507          continue

c              compute g 

               do 509 i_res=jstrt,jfins-1
                  if (myst(i_res) .ne. 0.0) then
                  myst(i_res)= (0.25*myst2(i_res)/myst(i_res))
                  end if
  509          continue 

               do 511 i_axis=1,3

                  do 524 i_res=jstrt,jfins-1
                    if (i_odd_even.eq.1) then
                      prcord(i_res,i_axis,i_pro,3)=
     *                prcord(i_res,i_axis,i_pro,3) +
     *                lngth(i_res,i_axis,i_pro,1)*myst(i_res)
 
                      prcord(i_res,i_axis,i_pro,1)=
     *                prcord(i_res,i_axis,i_pro,1) -
     *                lngth(i_res,i_axis,i_pro,1)*myst(i_res)
                    else
                           prcord(i_res+1,i_axis,i_pro,1)=
     *                prcord(i_res+1,i_axis,i_pro,1) +
     *                lngth(i_res,i_axis,i_pro,2)*myst(i_res)

                      prcord(i_res,i_axis,i_pro,3)=
     *                prcord(i_res,i_axis,i_pro,3) -
     *                lngth(i_res,i_axis,i_pro,2)*myst(i_res) 
                    endif
  524             continue

  511          continue   ! end of loop over axis

  508       continue        ! End of loop over odd, even

  500    continue        ! End of loop over max shake iterations

        write(SO,*) 'Failed to pass max iterations in shakox'
        stop 'passed max iterations in shakox'

501     continue                ! end of loop over proteins
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      continue

c     ---------------------- done ----------------------
     
      return
      end
