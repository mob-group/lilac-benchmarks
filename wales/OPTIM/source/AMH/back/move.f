
c     --------------------- move  ----------------------

      subroutine move(jstrt,jfins,
     *                ishkit,nmdifv,i_quench)
 
c     --------------------------------------------------

c     MOVE   is the driving subroutine for carrying out
c            the molecular dynamics: annealing schedule,
c            Verlet algorithm, and observable statistics

c     arguments:

c        jstrt - first modified site
c        jfins - last modifies site
c        ishkit- used to track number of shake 
c                iterations
c        nmdifv- number of intermediate structures saved
c                on movie tape
c     ---------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                               filtration
c        i_pro                        index over proteins
c        i_step                        index over time steps per temperature
c        i_temp                        index over temperatures
c        itcnt                        increment T (temperature) counter

c     nmstep is the number of time steps per temperature

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use globals, only:SO, i513,numpro,i1,avep,nmres,prcord,qrcord,
     *  srcord,trcord,zrcord,velocp,bondln,temtur,tolshk,maxshk,
     *  timstp,numcrd,oarchv,work1,work3,work4,iseed_amh,ires,nmtemp,
     *  movanal,incmov,omovi,omoviseg,nmstep,oPE_no_bias,oKE,
     *  ohdrgn,orama,ooxy,ochiral,oamh,oamhsr,orep,oPE_with_bias,
     *  occev,ooev,i526,i525,i527,i528,i514,nmdif,oamhlr,
     *  oamhmr,obias_Rg,ohdrgn_s,ohdrgn_m,ohdrgn_l,ohdrgn_seq,
     *  ononadd,ocon_P_AP,i_Etim_step,oPE_plus_KE,
     *  i_Etim,quench,nquench,n_Q_anneal_a,n_Q_anneal_b,
     *  Q0_a,Q0_b,Q0_inc_a,Q0_inc_b,Qvalue_a,Qvalue_b,
     *  maxpro,oobiassega,oobiassegb,Q0_safe_a, Q0_safe_b,
     *  itgrd,temtur_quench

      use amh_interfaces, only:E_write
      use globals_alt, only:do_send_output,count_alt,T_alt

      implicit none


c     passed argument declarations:

         integer jstrt,jfins,ishkit,i_quench,i504,i501


c     internal variables:

         logical tempav,bdshak

         integer itcnt,
     *           nmdifv,ibegn,idummy,
     *           nmovies,i
         real temph,totke(maxpro)

        integer i_pro, i_step, i_temp, i_grid_num

c     required subroutines

         external mdctrl,ynorm,OBinitv

        nmdifv = 0
        temph = 0.0

c     --------------------- begin -----------------------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize archive arrays

         do 516 i_pro=1,numpro
            do 100 i1=1,50
              avep(i_pro,1,i1)=0.0
              avep(i_pro,2,i1)=0.0
100            continue
  516    continue

c     ishkit tracks the number of shake iterations

      ishkit=0

c     if bdshak is true, then shake failed

      bdshak=.false.

c    new annealling schedule for quenching   
c    move temtur_quench into temtur from restrt.f
     
         if (quench) then
          do 504 i504=1, nquench
c           do 540 i_grid_num=1,4
            if( itgrd(1).gt.0 )then
             do 501 i501=1,itgrd(1)
              temtur(i501)=temtur_quench(i501,i_quench)
c          write(6,*)'in move temtur_quench steps structure',
c     * temtur(i501),temtur_quench(i501,i504),i_grid_num,i501,i_quench

501          continue
            endif  !  itgrd(grid_num)
c540        continue
504       continue
         endif

c     if no heatbath, then set initial velocities

c################# Start of Normal Run #################

         call initv(prcord,qrcord,srcord,
     *              trcord,zrcord,velocp,bondln,temtur(1),
     *              jstrt,jfins,tolshk,maxshk,bdshak,
     *              timstp,numpro,ishkit,
     *              numcrd,oarchv,work1,
     *              work3,work4,iseed_amh,ires)
         ibegn=1


cccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     set counter for T averages
      itcnt=1


c     write header for movie file
      if( .not.movanal )then
         nmdifv=int( float(nmtemp)/float(incmov) )

         if ( .not.quench ) then
         write(omovi,334)nmres,numcrd,numpro,nmdifv
         write(omoviseg,334)nmres,numcrd,numpro,nmdifv
         elseif (i_quench.eq.1) then
         write(omovi,334)nmres,numcrd,numpro,nquench
         endif
  334    format(4(i8,1x),' nmres nmcrd numpro nmsnap')
      endif


c     set counter for number of saved structures for
c     evaluating |R|

      nmdifv=0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c!!!!!!!!!!!!!!!  Start of loop over temperatures  !!!!!!!!!!!!!!!!!!!!!!
c     perform nmtemp timesteps for each protein


      if (movanal) then
        write(SO,*) 'opening movie file'
        open(11,file='movie_in',status='old')

          read(11,1040) nmovies
          nmtemp=nmovies
          incmov=1
1040        format(31x,i4)
      endif





      do 509 i_temp=ibegn,nmtemp
        if (i_temp.lt.n_Q_anneal_a) then
c  temp fix dividing by 0 !!!!!
         Q0_a = Q0_safe_a + 
     *           Q0_inc_a*real(n_Q_anneal_a-i_temp)/real(n_Q_anneal_a)
         Q0_b = Q0_safe_b + 
     *           Q0_inc_b*real(n_Q_anneal_b-i_temp)/real(n_Q_anneal_b)
        else
          Q0_a=Q0_safe_a
          Q0_b=Q0_safe_b
        endif

c      write(SO,474)i_quench,i_temp,temtur(i_temp),Qvalue_a,Qvalue_b
c474    format('i_quench,itemp,temp,Qa,Qb',
c     *                              i1,1x,i5,1x,f4.2,1x,f6.4,1x,f6.4)

      write(SO,474)i_temp,temtur(i_temp),Qvalue_a,Qvalue_b
474    format('itemp,temp,Qa,Qb',
     *                              1x,i5,1x,f4.2,1x,f6.4,1x,f6.4)




        if (i_temp.eq.nmtemp) then
            write(SO,*) 'simulation completed'
        endif

        if (movanal) then
          read(11,1041) idummy,idummy,idummy,temtur(i_temp),idummy
          write(SO,*) 'temperature is',temtur(i_temp)
1041      format(3(i6,1x),f8.4,1x,i5)
          do 10175 i=1,nmres
               read(11,1020) prcord(i,1,1,1),
     *          prcord(i,2,1,1),prcord(i,3,1,1),
     *          prcord(i,1,1,2),
     *          prcord(i,2,1,2),prcord(i,3,1,2),
     *          prcord(i,1,1,3),
     *          prcord(i,2,1,3),prcord(i,3,1,3)
1020            format(4x,3(f8.3,1x),4x,
     *           3(f8.3,1x),4x,3(f8.3,1x))

10175     continue
        endif

c        set RMS velocity for current T

         temph=sqrt(temtur(i_temp))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        determine whether averages should be carried
c        at this T

         if (quench) then
           if (i_temp.eq.nmtemp) then !if quench only print info at end
              tempav=.true.
           else
             tempav=.false.
           endif
         else
           if( mod(i_temp,incmov).eq.0 )then !in normal case print info every incmov steps
              tempav=.true.
           else
              tempav=.false.
           endif
         endif


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c        perform nmstep timesteps for the current temperature

         if (movanal) nmstep=1
 
         do 510 i_step=1,nmstep
            count_alt=i_temp/incmov                                   ! Start For altpot
            T_alt=temtur(i_temp)                                      !
            do_send_output=.false.                                    !
            if( (tempav .and. i_step.eq.nmstep)                       !
     *        .or. (i_Etim .and. mod(i_step,i_Etim_step).eq.0) )then  ! 
              do_send_output=.true.                                   !
            endif                                                     ! End For altpot



ccccccccccccccccccccccc
c           write(6,*) 'calling mdctrl'
c  (ohdrgn)  =  hbond potential
c  (ohdrgn_s) = hbond short range
c  (ohdrgn_m) = hbond medium range
c  (ohdrgn_l) = hbond long range
c  (orama) = rama potential
c  (ooxy) =  oxygen potential
c  (ochiral) chirality potential
c  (oPE) totalP.plot = (temp, total potential + KE,number)
c  (oKE) KE.plot  =  temp, total kinetic , number
c  (oamh) = amh.plot = amh total potential energy
c  (oamhsr) amhsr.plot= temp, amh short range potential energy, number
c  (oamhlr) amhlr.plot= temp,amh long range potential energy , number
c  (oamhmr) amhmr.plot= temp,amh medium  range potential energy, number
c  (occev) ccev.plot =  potential due to carbon (a/b) excluded volume
c  (ooev) ooev.plot  =  oxygen exculded
c  (obias) bias.plot  =  Q biasing potential
c  (obias_Rg) = radius of Gyration bias
c  (ototal) total.plot = (temp,total potential,number,without biasing)

c  avep(1,1,1,itcnt)  = hydrogen bond potential
c  avep(1,1,2,itcnt)  = rama potential
c  avep(1,1,3,itcnt)  = oxygen exculded
c  avep(1,1,4,itcnt)  = chirality potential
c  avep(1,1,5,itcnt)  = total amh short + medium + long
c  avep(1,1,6,itcnt)  = ------
c  avep(1,1,7,itcnt)  = amh short range
c  avep(1,1,8,itcnt)  = amh long range
c  avep(1,1,9,itcnt)  = potential due to carbon (a/b) excluded volume
c  avep(1,1,10,itcnt) =  Q biasing potental
c  avep(1,1,11,itcnt) = oxy potential
c  avep(1,1,12,itcnt) =  radius of Gyration bias
c  avep(1,1,13,itcnt) = amh medium range
c  avep(1,1,14,itcnt) = hydrogen bonds short range
c  avep(1,1,15,itcnt) = hydrogen bonds medium range
c  avep(1,1,16,itcnt) = hydrogen bonds long range
c  avep(1,1,17,itcnt) = nonaddative
c  avep(1,1,18,itcnt) = Q biasing potental seg




            call mdctrl(jstrt,jfins,tempav,ishkit,bdshak,itcnt,totke)
            if( (tempav .and. i_step.eq.nmstep)
     *        .or. (i_Etim .and. mod(i_step,i_Etim_step).eq.0) )then
              call E_write(avep(:,:,:),totke,temtur(i_temp),
     *                                  numpro,i_temp,incmov)
            endif

ccccccccccccccccccccccc
c           shake was unsuccessful, return and perform final
c           analysis

            if( bdshak )then
               write(oarchv,144)i_temp,temtur(i_temp),i_step
  144          format(/'Sorry bucky -- shake misfire ',
     *                i6,' T ',1pe10.3,' t ',i5)

               go to 300

            endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        save structure for later movie analysis, if the time is right
         if( mod(i_temp,incmov).eq.0 .and. (.not.movanal) 
     *        .and. (i_step.eq.nmstep-1.or.nmstep.eq.1)) then
c        if( 
c    *        mod(i_step,20).eq.0) then

c           save the curent structure

            nmdifv=nmdifv + 1

            do 526 i526=1,numpro
                if (.not. quench) then
              write(omovi,683)i526,nmdifv,i_step-1,temtur(i_temp),i_temp
 683          format(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')
              do 525 i525=1,nmres
                  write(omovi,332)
     *            (prcord(i525,i1,i526,1),i1=1,3),
     *            (prcord(i525,i1,i526,2),i1=1,3),
     *            (prcord(i525,i1,i526,3),i1=1,3)

332     format('CA: ',3(f8.3,1x),'CB: ',3(f8.3,1x),'Ox: ', 3(f8.3,1x))
  525         continue
                endif   ! quench

  526       continue


            if (i_temp.eq.nmtemp) then
            do 527 i527=1,numpro
                if (quench) then
                write(omovi,683)i527,nquench,i_quench,
     *                                   temtur(i_temp),i_temp
                else
                write(omoviseg,683)i527,nmdifv,i_step-1,
     *                                   temtur(i_temp),i_temp
                endif
                do 528 i528=1,nmres
                    if (quench) then
                     write(omovi,332)
     *              (prcord(i528,i1,i527,1),i1=1,3),
     *              (prcord(i528,i1,i527,2),i1=1,3),
     *              (prcord(i528,i1,i527,3),i1=1,3)
                    else
                    write(omoviseg,332)
     *              (prcord(i528,i1,i527,1),i1=1,3),
     *              (prcord(i528,i1,i527,2),i1=1,3),
     *              (prcord(i528,i1,i527,3),i1=1,3)
                    endif
  528           continue
  527         continue
            endif

         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




  510    continue
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c        increment counter for observables as 
c        function of T

         if( tempav )itcnt=itcnt + 1

c        perform collisions for next T

         call initv(prcord,qrcord,srcord,
     *              trcord,zrcord,velocp,bondln,
     *              temtur(i_temp),jstrt,jfins,tolshk,
     *              maxshk,bdshak,timstp,numpro,
     *              ishkit,numcrd,
     *              oarchv,work1,work3,
     *              work4,iseed_amh,ires)

  509 continue

        if (i_quench.eq.nquench) then
        close(ohdrgn)
        close(ohdrgn_s)
        close(ohdrgn_m)
        close(ohdrgn_l)
        close(ohdrgn_seq)
        close(orama)
        close(ooxy)
        close(ochiral)
        close(oPE_no_bias)
        close(oPE_with_bias)
        close(oPE_plus_KE)
        close(oKE)
        close(oamh)
        close(oamhsr)
        close(oamhlr)
        close(oamhmr)
        close(orep)
        close(ocon_P_AP)
        close(ononadd)
        close(occev)
        close(ooev)
c        close(obias)
        close(obias_Rg)
        close(oobiassega)
        close(oobiassegb)

        endif

c       write(6,*)' move '
c       write(6,*)'oarchv 30 ',oarchv
c       write(6,*)'omovi 55 ',omovi
c       write(6,*)'ohdrgn  ',ohdrgn
c       write(6,*)'ohdrgn_s 81 ',ohdrgn_s
c       write(6,*)'ohdrgn_m 82',ohdrgn_m
c       write(6,*)'ohdrgn_l 83',ohdrgn_l
c       write(6,*)'ohdrgn_seq 85',ohdrgn_seq
c       write(6,*)'orama 61 ',orama
c       write(6,*)'ooxy 62',ooxy
c       write(6,*)'ochiral 63 ',ochiral
c       write(6,*)'oamh  64',oamh
c       write(6,*)'ototal 65 ',ototal
c       write(6,*)'oamhsr 66 ',oamhsr
c       write(6,*)'oPE 67 ',oPE
c       write(6,*)'oKE  73',oKE
c       write(6,*)'oamhlr 68 ',oamhlr
c       write(6,*)'oamhmr 78 ',oamhmr
c       write(6,*)'ononadd 84 ',ononadd
c       write(6,*)'oamhlr 68 ',oamhlr
c       write(6,*)'oamhmr 78 ',oamhmr
c       write(6,*)'ononadd  84',ononadd
c       write(6,*)'ocon_P_AP 86',ocon_P_AP
c       write(6,*)'obias_Rg 77',obias_Rg
c       write(6,*)'omoviseg 75',omoviseg
c       write(6,*)'oobiassega 87 ',oobiassega
c       write(6,*)'oobiassegb 88',oobiassegb


c!!!!!!!!!!!!!!!!!!!   End of loop over temperatures !!!!!!!!!!!!!!!!!!!!!!!!

  300 continue

c     prepare 'observables' for printing as a 
c     function of T

c     set normalization factors

      itcnt=itcnt - 1

      temph=1.0/float(nmstep)



c        normalize 1st moments

      do 514 i514=1,itcnt
         call ynorm(numpro,1,numpro,avep(1,1,5),
     *              temph)
  514 continue


c     ---------------------- done -----------------------

      return
      end
