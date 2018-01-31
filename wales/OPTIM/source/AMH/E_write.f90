      subroutine E_write(avep,totke,T,numpro,i_temp,incmov)

      use amhglobals,  only: weight_P_AP,oKE,oPE_no_bias,oPE_with_bias,oPE_backbone, &
            oPE_plus_KE,ohdrgn,ohdrgn_s,ohdrgn_m,&
            ohdrgn_l,ohdrgn_seq,ononadd,ocon_P_AP,orama,ooxy,ochiral,oamh,oamhsr,&
            oamhmr,oamhlr,occev,ooev,obias_Rg,maxpro,igomb,orep, &
            oobiassega,oobiassegb

      implicit none
       double precision, intent(in):: avep(:,:,:),totke(maxpro),T
      integer, intent(in):: numpro,i_temp,incmov 

       double precision, dimension(size(avep,1)):: E_P_AP,E_PE_no_bias,E_PE_with_bias,E_PE_backbone
      integer:: count,i_pro

! preliminary calculations
     
!  biasseg_a  avep(:,1,10)
!  bias seg_b avep(:,1,18)
 
!          E_P_AP(:)=-weight_P_AP(1)*avep(:,1,40)-weight_P_AP(2)*avep(:,1,41)-weight_P_AP(3)*avep(:,1,42)
!          E_PE_no_bias(:)=avep(:,1,1)+avep(:,1,2)+avep(:,1,3)+avep(:,1,4)+avep(:,1,5)+avep(:,1,9)+ &
!                avep(:,1,11)+avep(:,1,17)+E_P_AP
!          E_PE_with_bias(:)=E_PE_no_bias(:)+avep(:,1,10)+avep(:,1,12)+avep(:,1,18)
!          E_PE_backbone(:)= avep(:,1,2)+avep(:,1,3)+avep(:,1,4)+avep(:,1,9)+avep(:,1,11)
          count=i_temp/incmov

!   10 bias Qa  12 Rg bias 18 bias Qb
1000      format(i8,2x,100(f12.3,2x))

! potential/kinetic energies etc

          write(oKE,1000) count,T,totke(1:numpro)
          write(oPE_no_bias,1000)  count,T,E_PE_no_bias(1:numpro)
          write(oPE_with_bias,1000) count,T,E_PE_with_bias(1:numpro)
          write(oPE_plus_KE,1000) count,T,totke(1:numpro)+E_PE_with_bias(1:numpro)
          write(oPE_backbone,1000) count,T,E_PE_backbone(1:numpro)

! replica bias energy

          write(orep,1000) count,T,avep(1,1,43)

! hbond contributions
! hbond contributions
! ohdrgn      count,T,avep(1:numpro,1,1)  =  number , temp  ,   total hbond
! ohdrgn_s    count,T,avep(1:numpro,1,14),avep(i_pro,1,21:24)  = number, temp, pre pnas ab hbond, pnas ab short range
! ohdrgn_m    count,T,avep(1:numpro,1,15),avep(i_pro,1,25:28)  = number, temp, pre pnas ab hbond, pnas ab medium range
! ohdrgn_l    count,T,avep(1:numpro,1,16),avep(i_pro,1,29:32)  = number, temp, pre pnas ab hbond, pnas ab long range
! ohdrgn_seq  count,T,avep(1:numpro,1,16),avep(i_pro,1,33:39)  = number, temp, pre pnas ab hbond, pnas ab seq

          write(ohdrgn,1000) count,T,avep(1:numpro,1,1)
          write(ohdrgn_s,1000) count,T,(avep(i_pro,1,14),avep(i_pro,1,21:24),i_pro=1,numpro)
          write(ohdrgn_m,1000) count,T,(avep(i_pro,1,15),avep(i_pro,1,25:28),i_pro=1,numpro)
          write(ohdrgn_l,1000) count,T,(avep(i_pro,1,16),avep(i_pro,1,29:32),i_pro=1,numpro)
          write(ohdrgn_seq,1000) count,T,(avep(i_pro,1,33:39),i_pro=1,numpro)

! other contributions

          write(ononadd,1000) count,T,avep(1:numpro,1,17)
          write(ocon_P_AP,1000) count,T,(E_P_AP(i_pro),avep(i_pro,1,40), &
                                           avep(i_pro,1,41),avep(i_pro,1,42),i_pro=1,numpro)
          write(orama,1000) count,T,avep(1:numpro,1,2)
          write(ooxy,1000) count,T,avep(1:numpro,1,3)
          write(ochiral,1000) count, T,avep(1:numpro,1,4)
          write(oamh,1000) count,T,avep(1:numpro,1,5)    

!          if (.not.igomb) then
          write(oamhsr,1000) count,T,avep(1:numpro,1,7)
          write(oamhmr,1000) count,T,avep(1:numpro,1,13)
          write(oamhlr,1000) count,T,avep(1:numpro,1,8)
!          endif

          write(occev,1000) count, T,avep(1:numpro,1,9)
          write(ooev,1000) count,T,avep(1:numpro,1,11)
!        write(obias,1000) count,T,avep(1:numpro,1,10)
          write(oobiassega,1000) count,T,avep(1:numpro,1,10)
          write(oobiassegb,1000) count,T,avep(1:numpro,1,18)
!          write(obias_Rg,1000) count, T,avep(1:numpro,1,12)

      end
