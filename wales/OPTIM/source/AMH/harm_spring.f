
c     --------------------- shkdrv  ----------------------

      subroutine harm_spring(prcord,jstrt,jfins,maxpro,maxcrd, 
     *                         ires,f_cord_harm,E_temp_harm)

c     ---------------------------------------------------

c     harm_spring driver routine for performing hard springs
c                 for maintain chain connectivity

c     arguments:

c        AMHmaxsiz- maximum number of protein residues (i)
c        prcord- new coordinates which satsify bond 
c                lengths (o)

c     ---------------------------------------------------

      use amhglobals,  only: AMHmaxsiz,E_harm_springs

      implicit none

c     argument declarations:

         integer jstrt,jfins,maxpro,maxcrd,ires(AMHmaxsiz)
     
         double precision prcord(AMHmaxsiz,3,maxpro,maxcrd),
     *     x_diff(AMHmaxsiz),y_diff(AMHmaxsiz),z_diff(AMHmaxsiz),
     *     f_cord_harm(AMHmaxsiz,maxcrd,maxcrd),ca_dist(AMHmaxsiz),ca_cb_dist(AMHmaxsiz),E_temp_harm,caf,
     *     r_dev,CA_distance,ca_spring_scl,CA_CB_distance,ca_cb_spring_scl,
     *     ca_ox_dist(AMHmaxsiz),CA_OX_distance,ca_ox_spring_scl,
     *     ox_caplus_dist(AMHmaxsiz),OX_CAplus_distance,ox_caplus_spring_scl

c     internal variables:

c        --- do loop indices ---

         integer i507,i512

c        --- implied do loop indices ---

c  from pttarg.f   ca-ca+1     bondln(i505,1)=3.8004
c  from pttarg.f   cb-ca       bondln(i505,2)=1.54
c  from shakox.f   ox-ca       eqdist(1)=2.42677
c  from shakox.f   ca+1-ox     eqdist(2)=2.82146

c CA-CA

        ca_spring_scl=50.0D0
        CA_distance=3.8004D0       
        E_temp_harm = 0.D0
        f_cord_harm=0.D0

           do 507 i507=jstrt,jfins-1

           caf=0.0D0

        x_diff(i507)=prcord(i507+1,1,1,1) - prcord(i507,1,1,1)
        y_diff(i507)=prcord(i507+1,2,1,1) - prcord(i507,2,1,1)
        z_diff(i507)=prcord(i507+1,3,1,1) - prcord(i507,3,1,1)

        ca_dist(i507)=dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)

        r_dev=ca_dist(i507)-CA_distance

        E_temp_harm =  E_temp_harm +  ca_spring_scl*r_dev**2
        caf=2.0D0*ca_spring_scl*r_dev/ca_dist(i507)

        f_cord_harm(i507,1,1) = f_cord_harm(i507,1,1) + x_diff(i507)*caf
        f_cord_harm(i507,2,1) = f_cord_harm(i507,2,1) + y_diff(i507)*caf
        f_cord_harm(i507,3,1) = f_cord_harm(i507,3,1) + z_diff(i507)*caf
        f_cord_harm(i507+1,1,1) = f_cord_harm(i507+1,1,1) - x_diff(i507)*caf
        f_cord_harm(i507+1,2,1) = f_cord_harm(i507+1,2,1) - y_diff(i507)*caf
        f_cord_harm(i507+1,3,1) = f_cord_harm(i507+1,3,1) - z_diff(i507)*caf

 507         continue
           
c           write(6,*)'ca ca ',E_temp_harm

c        if (memires(jtgres).ne.8)then

          ca_cb_spring_scl=50.0D0
          CA_CB_distance=1.54D0       

           do 607 i507=jstrt,jfins
            if (ires(i507).ne.8)then 

           caf=0.D0

          x_diff(i507)=prcord(i507,1,1,1) - prcord(i507,1,1,2)
          y_diff(i507)=prcord(i507,2,1,1) - prcord(i507,2,1,2)
          z_diff(i507)=prcord(i507,3,1,1) - prcord(i507,3,1,2)

            ca_cb_dist(i507)=dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)

        r_dev=ca_cb_dist(i507)-CA_CB_distance

        E_temp_harm = E_temp_harm + ca_cb_spring_scl*r_dev**2
        caf=2.0D0*ca_cb_spring_scl*r_dev/ca_cb_dist(i507)

        f_cord_harm(i507,1,1) = f_cord_harm(i507,1,1) - x_diff(i507)*caf
        f_cord_harm(i507,2,1) = f_cord_harm(i507,2,1) - y_diff(i507)*caf
        f_cord_harm(i507,3,1) = f_cord_harm(i507,3,1) - z_diff(i507)*caf
        f_cord_harm(i507,1,2) = f_cord_harm(i507,1,2) + x_diff(i507)*caf
        f_cord_harm(i507,2,2) = f_cord_harm(i507,2,2) + y_diff(i507)*caf
        f_cord_harm(i507,3,2) = f_cord_harm(i507,3,2) + z_diff(i507)*caf

         endif
c         write(6,*)'harmharm ', r_dev

 607         continue

c          write(6,*)'ca cb E ',E_temp_harm

c        write(6,*)'ca ca distance ', ca_dist(5),CA_distance 
c        write(6,*)'ca cb distance ', ca_cb_dist(5),CA_CB_distance 
c        write(6,*)'E_temp_harm ', E_temp_harm

       ca_ox_spring_scl=50.0D0
       CA_OX_distance= 2.42677D0  

c     Shake between Ca-Ox
c  from shakox.f   ca-ox       eqdist(1)=2.42677

           do 707 i507=jstrt,jfins

           caf=0.0D0

          x_diff(i507)=prcord(i507,1,1,1) - prcord(i507,1,1,3)
          y_diff(i507)=prcord(i507,2,1,1) - prcord(i507,2,1,3)
          z_diff(i507)=prcord(i507,3,1,1) - prcord(i507,3,1,3)

            ca_ox_dist(i507)=dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)

        r_dev=ca_ox_dist(i507)-CA_OX_distance

        E_temp_harm =  E_temp_harm +  ca_ox_spring_scl*r_dev**2
        caf=2.0D0*ca_ox_spring_scl*r_dev/ca_ox_dist(i507)

        f_cord_harm(i507,1,1) = f_cord_harm(i507,1,1) - x_diff(i507)*caf
        f_cord_harm(i507,2,1) = f_cord_harm(i507,2,1) - y_diff(i507)*caf
        f_cord_harm(i507,3,1) = f_cord_harm(i507,3,1) - z_diff(i507)*caf
        f_cord_harm(i507,1,3) = f_cord_harm(i507,1,3) + x_diff(i507)*caf
        f_cord_harm(i507,2,3) = f_cord_harm(i507,2,3) + y_diff(i507)*caf
        f_cord_harm(i507,3,3) = f_cord_harm(i507,3,3) + z_diff(i507)*caf

 707         continue

c         write(6,*)'ca ox E ',E_temp_harm


c     Shake between Ca-Ox

          ox_caplus_spring_scl=50.0D0
          OX_CAplus_distance= 2.82146D0
    
c  from shakox.f   ca-ox       eqdist(1)=2.42677
c  from shakox.f   ox-ca+1     eqdist(2)=2.82146

           do 807 i507=jstrt,jfins-1

           caf=0.0D0

          x_diff(i507)=prcord(i507,1,1,3) - prcord(i507+1,1,1,1)
          y_diff(i507)=prcord(i507,2,1,3) - prcord(i507+1,2,1,1)
          z_diff(i507)=prcord(i507,3,1,3) - prcord(i507+1,3,1,1)

            ox_caplus_dist(i507)=dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)

        r_dev=ox_caplus_dist(i507)-OX_CAplus_distance 

        E_temp_harm =  E_temp_harm +  ox_caplus_spring_scl*r_dev**2
        caf=2.0D0*ox_caplus_spring_scl*r_dev/ox_caplus_dist(i507)

        f_cord_harm(i507,1,3) = f_cord_harm(i507,1,3) - x_diff(i507)*caf
        f_cord_harm(i507,2,3) = f_cord_harm(i507,2,3) - y_diff(i507)*caf
        f_cord_harm(i507,3,3) = f_cord_harm(i507,3,3) - z_diff(i507)*caf
        f_cord_harm(i507+1,1,1) = f_cord_harm(i507+1,1,1) + x_diff(i507)*caf
        f_cord_harm(i507+1,2,1) = f_cord_harm(i507+1,2,1) + y_diff(i507)*caf
        f_cord_harm(i507+1,3,1) = f_cord_harm(i507+1,3,1) + z_diff(i507)*caf

 807         continue


c     There is one more harmonic potential for the terminal oxygen,
c     which is different to GMIN in order to prevent to states having 
c     the same energy.

         do 907 i507=jfins,jfins

           caf=0.0D0
           ox_caplus_distance= 2.4D0
!!   ox-cb
         x_diff(i507)=prcord(i507,1,1,3) - prcord(i507,1,1,2)
         y_diff(i507)=prcord(i507,2,1,3) - prcord(i507,2,1,2)
         z_diff(i507)=prcord(i507,3,1,3) - prcord(i507,3,1,2)

        ox_caplus_dist(i507)= dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)
      r_dev=ox_caplus_dist(i507)-OX_CAplus_distance

      E_temp_harm =  E_temp_harm +  ox_caplus_spring_scl*r_dev**2
      caf=2.0D0*ox_caplus_spring_scl*r_dev/ox_caplus_dist(i507)

      f_cord_harm(i507,1,3) = f_cord_harm(i507,1,3) - x_diff(i507)*caf
      f_cord_harm(i507,2,3) = f_cord_harm(i507,2,3) - y_diff(i507)*caf
      f_cord_harm(i507,3,3) = f_cord_harm(i507,3,3) - z_diff(i507)*caf
      f_cord_harm(i507,1,2) = f_cord_harm(i507,1,2) + x_diff(i507)*caf
      f_cord_harm(i507,2,2) = f_cord_harm(i507,2,2) + y_diff(i507)*caf
      f_cord_harm(i507,3,2) = f_cord_harm(i507,3,2) + z_diff(i507)*caf

         OX_CAplus_distance=4.40D0

       x_diff(i507)=prcord(i507,1,1,3) - prcord(i507-1,1,1,3)
       y_diff(i507)=prcord(i507,2,1,3) - prcord(i507-1,2,1,3)
       z_diff(i507)=prcord(i507,3,1,3) - prcord(i507-1,3,1,3)

        ox_caplus_dist(i507)= dsqrt(x_diff(i507)**2 + y_diff(i507)**2 + z_diff(i507)**2)

      r_dev=ox_caplus_dist(i507)-OX_CAplus_distance

      E_temp_harm =  E_temp_harm +  ox_caplus_spring_scl*r_dev**2
      caf=2.0D0*ox_caplus_spring_scl*r_dev/ox_caplus_dist(i507)

       E_harm_springs =  E_temp_harm

      f_cord_harm(i507,1,3) = f_cord_harm(i507,1,3) - x_diff(i507)*caf
      f_cord_harm(i507,2,3) = f_cord_harm(i507,2,3) - y_diff(i507)*caf
      f_cord_harm(i507,3,3) = f_cord_harm(i507,3,3) - z_diff(i507)*caf
      f_cord_harm(i507-1,1,3) = f_cord_harm(i507-1,1,3) + x_diff(i507)*caf
      f_cord_harm(i507-1,2,3) = f_cord_harm(i507-1,2,3) + y_diff(i507)*caf
      f_cord_harm(i507-1,3,3) = f_cord_harm(i507-1,3,3) + z_diff(i507)*caf

907         continue

c     ---------------------- done -----------------------

      return
      end
