
!     --------------------- move  ----------------------

      subroutine move(jstrt,jfins,ishkit,nmdifv,i_quench)
 
!     --------------------------------------------------

!     MOVE   is the driving subroutine for carrying out
!            the molecular dynamics: annealing schedule,
!            Verlet algorithm, and observable statistics

!     arguments:

!        jstrt - first modified site
!        jfins - last modifies site
!        ishkit- used to track number of shake 
!                iterations
!        nmdifv- number of intermediate structures saved
!                on movie tape
!     ---------------------------------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                               filtration
!        i_pro                        index over proteins
!        i_step                        index over time steps per temperature
!        i_temp                        index over temperatures
!        itcnt                        increment T (temperature) counter

!     nmstep is the number of time steps per temperature

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     use globals, only:SO,weight_P_AP, i513,numpro,i1,avep,nmres,prcord,qrcord, &
       srcord,trcord,zrcord,velocp,bondln,temtur,tolshk,maxshk, &
       timstp,numcrd,oarchv,work1,work3,work4,iseed,ires,nmtemp, &
       movanal,incmov,omovi,omoviseg,nmstep,oPE_no_bias,oKE, &
       ohdrgn,orama,ooxy,ochiral,oamh,oamhsr,orep,oPE_with_bias, &
       occev,ooev,i526,i525,i527,i528,i514,nmdif,oamhlr, &
       oamhmr,obias_Rg,ohdrgn_s,ohdrgn_m,ohdrgn_l,ohdrgn_seq, &
       ononadd,ocon_P_AP,i_Etim_step,oPE_plus_KE, &
       i_Etim,quench,nquench,n_Q_anneal_a,n_Q_anneal_b, &
       Q0_a,Q0_b,Q0_inc_a,Q0_inc_b,Qvalue_a,Qvalue_b, &
       maxpro,oobiassega,oobiassegb,Q0_safe_a, Q0_safe_b, &
       itgrd,temtur_quench

!       use amh_interfaces, only:E_write,Ptr_Sub,arbitrate
       use amh_interfaces, only:E_write,Ptr_Sub
      use globals_alt, only:do_send_output,count_alt,T_alt

      implicit none

!     passed argument declarations:

         integer :: jstrt,jfins,ishkit,i_quench,i504,i501

!     internal variables:

         logical :: tempav,bdshak
 
         integer :: itcnt,nmdifv,ibegn,idummy,nmovies,i
         real :: temph,totke(maxpro)

!,arbitrate

     integer :: i_pro, i_step, i_temp, i_grid_num
     integer :: repout_move, rep_move
     double precision :: E_temp(5),T_temp(5)

     double precision :: DD(5), VARRRR(5)
     double precision:: CC(5), VARRR(5)
     INTEGER :: A(10), VARA(10) 
     INTEGER :: B(10), VARB(10)
     INTEGER :: C(10), VARC(10) 
     INTEGER :: D(10), VARD(10) 
     INTEGER :: BB(5), VARR(5)

     POINTER (ee, VARRRR)
     POINTER (tt, VARRR)
     POINTER (inet1, VARA) ! VAR is the pointee
                          ! p is the integer pointer
     POINTER (inet2, VARB) ! VAR is the pointee
     POINTER (inet3, VARC) ! VAR is the pointee
     POINTER (inet4, VARD) ! VAR is the pointee
     POINTER (port_number, VARR)

!     required subroutines
        external mdctrl,ynorm,initv
      
          ee = LOC(DD)
          tt = LOC(CC)
          inet1 = LOC (A)
          inet2 = LOC(B)
          inet3 = LOC(C)
          inet4 = LOC(D)
          port_number = LOC(BB)

          A(1) = 132
          B(1) = 239
          C(1) = 157
          D(1) = 126
          BB(1) = 123457

!         WRITE(6,*) 'Server IP info'
!         WRITE(6,*) 'inet1 = ', inet2
!         WRITE(6,*) 'inet2 = ', inet2
!         WRITE(6,*) 'inet3 = ', inet3
!         WRITE(6,*) 'inet4 = ', inet4 

        nmdifv = 0
        temph = 0.0

!     --------------------- begin -----------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initialize archive arrays

      do 513 i513=1,nmdif
         do 516 i_pro=1,numpro
            do 100 i1=1,50
              avep(i_pro,1,i1,i513)=0.0
              avep(i_pro,2,i1,i513)=0.0
100            continue
  516    continue
  513 continue

!     ishkit tracks the number of shake iterations

      ishkit=0

!     if bdshak is true, then shake failed

      bdshak=.false.

!    new annealling schedule for quenching   
!    move temtur_quench into temtur from restrt.f
     
         if (quench) then
          do 504 i504=1, nquench
!           do 540 i_grid_num=1,4
            if( itgrd(1).gt.0 )then
             do 501 i501=1,itgrd(1)
              temtur(i501)=temtur_quench(i501,i_quench)
!          write(6,*)'in move temtur_quench steps structure',
!     * temtur(i501),temtur_quench(i501,i504),i_grid_num,i501,i_quench

501          continue
            endif  !  itgrd(grid_num)
!540        continue
504       continue
         endif

!     if no heatbath, then set initial velocities

!################# Start of Normal Run #################

       call initv(prcord,qrcord,srcord, &
                  trcord,zrcord,velocp,bondln,temtur(1), &
                  jstrt,jfins,tolshk,maxshk,bdshak, &
                  timstp,numpro,ishkit, &
                  numcrd,oarchv,work1, &
                  work3,work4,iseed,ires)
         ibegn=1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!     set counter for T averages
      itcnt=1

!     write header for movie file
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


!     set counter for number of saved structures for
!     evaluating |R|

      nmdifv=0

!!!!!!!!!!!!!!!!  Start of loop over temperatures  !!!!!!!!!!!!!!!!!!!!!!
!     perform nmtemp timesteps for each protein


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
!  temp fix dividing by 0 !!!!!
         Q0_a = Q0_safe_a +  &
                Q0_inc_a*real(n_Q_anneal_a-i_temp)/real(n_Q_anneal_a)
         Q0_b = Q0_safe_b +  &
                Q0_inc_b*real(n_Q_anneal_b-i_temp)/real(n_Q_anneal_b)
        else
          Q0_a=Q0_safe_a
          Q0_b=Q0_safe_b
        endif

      write(SO,474)i_temp,temtur(i_temp),Qvalue_a,Qvalue_b
474    format('itemp,temp,Qa,Qb',  &
                         1x,i5,1x,f4.2,1x,f6.4,1x,f6.4)

        if (i_temp.eq.nmtemp) then
            write(SO,*) 'simulation completed'
        endif

        if (movanal) then
          read(11,1041) idummy,idummy,idummy,temtur(i_temp),idummy
          write(SO,*) 'temperature is',temtur(i_temp)
1041      format(3(i6,1x),f8.4,1x,i5)
          do 10175 i=1,nmres
               read(11,1020) prcord(i,1,1,1), &
               prcord(i,2,1,1),prcord(i,3,1,1), &
               prcord(i,1,1,2), &
               prcord(i,2,1,2),prcord(i,3,1,2), &
               prcord(i,1,1,3), &
               prcord(i,2,1,3),prcord(i,3,1,3)
1020       format(4x,3(f8.3,1x),4x,3(f8.3,1x),4x,3(f8.3,1x))

10175     continue
        endif

!        set RMS velocity for current T

         temph=sqrt(temtur(i_temp))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        determine whether averages should be carried
!        at this T

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        perform nmstep timesteps for the current temperature

         if (movanal) nmstep=1
 
         do 510 i_step=1,nmstep
            count_alt=i_temp/incmov                   ! Start For altpot
            T_alt=temtur(i_temp)                                   
            do_send_output=.false.                                 
            if( (tempav .and. i_step.eq.nmstep)     &              
             .or. (i_Etim .and. mod(i_step,i_Etim_step).eq.0) )then 
              do_send_output=.true.                                 
            endif                                     ! End For altpot

            call mdctrl(jstrt,jfins,tempav,ishkit,bdshak,itcnt,totke)
            if( (tempav .and. i_step.eq.nmstep) &
             .or. (i_Etim .and. mod(i_step,i_Etim_step).eq.0) )then
              call E_write(avep(:,:,:,itcnt),totke,temtur(i_temp), &
                                       numpro,i_temp,incmov)

      E_temp(1)=  &
         dble(avep(1,1,1,itcnt)+avep(1,1,2,itcnt)+avep(1,1,3,itcnt) + &
          avep(1,1,4,itcnt)+avep(1,1,5,itcnt)+avep(1,1,9,itcnt) + &
                     avep(1,1,11,itcnt)+avep(1,1,17,itcnt) + &
                       -weight_P_AP(1)*avep(1,1,40,itcnt) + &
                       -weight_P_AP(2)*avep(1,1,41,itcnt) + &
                       -weight_P_AP(3)*avep(1,1,42,itcnt))

         T_temp(1)=dble(temtur(i_temp))
          CC(1) = T_temp(1)
          DD(1) = E_temp(1)
        write(SO,*) "energy temp before arb call  ",E_temp(1),T_temp(1)

      rep_move=20
!      if ( mod(i_temp,rep_move) .eq. 0)then 
!        write(SO,*)"REPLICA MOVE i_step rep_move mod  ", &
!                         i_temp,rep_move,mod(i_temp,rep_move)
!       CALL arbitrate (ee,tt,inet1,inet2,inet3,inet4,port_number)
!      endif

      repout_move=5
!      if ( mod(i_temp,repout_move) .eq. 0) then
!        write(SO,*)"output i_temp rep_move mod  ",       &
!                    i_temp,repout_move,mod(i_temp,rep_move)
!      endif

      endif

!cccccccccccccccccccccc
!           shake was unsuccessful, return and perform final
!           analysis

            if( bdshak )then
               write(oarchv,144)i_temp,temtur(i_temp),i_step
  144          format(/'Sorry bucky -- shake misfire ', &
                     i6,' T ',1pe10.3,' t ',i5)
               go to 300

            endif


              write(oarchv,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXX'


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   save structure for later movie analysis, if the time is right
       if( mod(i_temp,incmov).eq.0 .and. (.not.movanal)  &
             .and. (i_step.eq.nmstep-1.or.nmstep.eq.1)) then
!        if( 
!    *        mod(i_step,20).eq.0) then

!           save the curent structure

            nmdifv=nmdifv + 1

            do 526 i526=1,numpro
                if (.not. quench) then
              write(omovi,683)i526,nmdifv,i_step-1,temtur(i_temp),i_temp
 683          format(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')
              do 525 i525=1,nmres
                  write(omovi,332) &
                 (prcord(i525,i1,i526,1),i1=1,3), &
                 (prcord(i525,i1,i526,2),i1=1,3), &
                 (prcord(i525,i1,i526,3),i1=1,3)

332     format('CA: ',3(f8.3,1x),'CB: ',3(f8.3,1x),'Ox: ', 3(f8.3,1x))
  525         continue
                endif   ! quench

  526       continue


            if (i_temp.eq.nmtemp) then
            do 527 i527=1,numpro
                if (quench) then
                write(omovi,683)i527,nquench,i_quench, &
                                        temtur(i_temp),i_temp
                else
                write(omoviseg,683)i527,nmdifv,i_step-1, &
                                        temtur(i_temp),i_temp
                endif
                do 528 i528=1,nmres
                    if (quench) then
                     write(omovi,332)                &
                   (prcord(i528,i1,i527,1),i1=1,3), &
                   (prcord(i528,i1,i527,2),i1=1,3), &
                   (prcord(i528,i1,i527,3),i1=1,3)
                    else
                    write(omoviseg,332)              &
                    (prcord(i528,i1,i527,1),i1=1,3), &
                    (prcord(i528,i1,i527,2),i1=1,3), &
                    (prcord(i528,i1,i527,3),i1=1,3)
                    endif
  528           continue
  527         continue
            endif

         endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




  510    continue
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!        increment counter for observables as 
!        function of T

         if( tempav )itcnt=itcnt + 1

!        perform collisions for next T

         call initv(prcord,qrcord,srcord,  &
                   trcord,zrcord,velocp,bondln, &
                   temtur(i_temp),jstrt,jfins,tolshk, &
                   maxshk,bdshak,timstp,numpro, &
                   ishkit,numcrd,  &
                   oarchv,work1,work3, &
                   work4,iseed,ires)

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
!        close(obias)
!        close(obias_Rg)
        close(oobiassega)
        close(oobiassegb)

        endif



!       write(6,*)' move '
!       write(6,*)'oarchv 30 ',oarchv
!       write(6,*)'omovi 55 ',omovi
!       write(6,*)'ohdrgn  ',ohdrgn
!       write(6,*)'ohdrgn_s 81 ',ohdrgn_s
!       write(6,*)'ohdrgn_m 82',ohdrgn_m
!       write(6,*)'ohdrgn_l 83',ohdrgn_l
!       write(6,*)'ohdrgn_seq 85',ohdrgn_seq
!       write(6,*)'orama 61 ',orama
!       write(6,*)'ooxy 62',ooxy
!       write(6,*)'ochiral 63 ',ochiral
!       write(6,*)'oamh  64',oamh
!       write(6,*)'ototal 65 ',ototal
!       write(6,*)'oamhsr 66 ',oamhsr
!       write(6,*)'oPE 67 ',oPE
!       write(6,*)'oKE  73',oKE
!       write(6,*)'oamhlr 68 ',oamhlr
!       write(6,*)'oamhmr 78 ',oamhmr
!       write(6,*)'ononadd 84 ',ononadd
!       write(6,*)'oamhlr 68 ',oamhlr
!       write(6,*)'oamhmr 78 ',oamhmr
!       write(6,*)'ononadd  84',ononadd
!       write(6,*)'ocon_P_AP 86',ocon_P_AP
!       write(6,*)'obias_Rg 77',obias_Rg
!       write(6,*)'omoviseg 75',omoviseg
!       write(6,*)'oobiassega 87 ',oobiassega
!       write(6,*)'oobiassegb 88',oobiassegb


!!!!!!!!!!!!!!!!!!!!   End of loop over temperatures !!!!!!!!!!!!!!!!!!!!!!!!

  300 continue

!     prepare 'observables' for printing as a 
!     function of T

!     set normalization factors

      itcnt=itcnt - 1

      temph=1.0/float(nmstep)

!        normalize 1st moments

      do 514 i514=1,itcnt
         call ynorm(numpro,1,numpro,avep(1,1,5,i514), &
                   temph)
514 continue


!     ---------------------- done -----------------------

      return
      end
