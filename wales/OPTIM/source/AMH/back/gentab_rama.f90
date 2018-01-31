	subroutine  gentab_rama

        use globals,only:aminoa,maxsiz,tgsequences_amc,rama_prob, &
                 ramascl,rama_force

	implicit none

!     internal variables:

!   rama_prob(phi_index,psi_index,amino,amino-1,amino,amino+1

        integer isit1,isit2,mvm_scr
        integer iaa,ires,ipre,ipost,open_status,i1,nmres

!     required subroutines

!       external 

!  data in 10 degree wide bins from 180 to -180
 
        rama_prob(:,:,:,:,:)=0.0

!!$        do iaa = 1,20
!!$         do ipre = 1,20
!!$          do ipost = 1,20
!!$
!!$        open(mvm_scr, &
!!$          file='~/marcio/bdtrimers/'//aminoa(iaa)//'/'//aminoa(ipre)//'_'//aminoa(iaa)//'_'//aminoa(ipost)//'.scr',status='old',iostat=open_status)
!!$               if (open_status.ne.0) then
!!$                 write(6,*) 'failure to open file in gentab_rama'
!!$                 stop
!!$               endif
!!$        do isit1 = 1,36
!!$        read(mvm_scr,*)(rama_prob(isit1,isit2,ipre,iaa,ipost),isit2=1,36)
!!$        enddo
!!$        close(mvm_scr)
!!$ 
!!$          enddo ! do ipost = 1,20
!!$         enddo ! do ipre = 1,20
!!$        enddo ! do iaa = 1,20
!!$
!!$!   local filtering  
!!$
!!$        do iaa = 1,20
!!$         do ipre = 1,20
!!$          do ipost = 1,20
!!$             do isit1 = 1,36
!!$              do isit2 = 1,36
!!$
!!$               if(rama_prob(isit1,isit2,ipre,iaa,ipost).gt.0.0)then
!!$                rama_prob(isit1,isit2,ipre,iaa,ipost) = & 
!!$                -ramascl* Log(rama_prob(isit1,isit2,ipre,iaa,ipost))
!!$               endif 
!!$
!!$                if((isit1.eq.1).and.(isit2.eq.1))then
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost) 
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 4
!!$               else if((isit1.eq.1).and.(isit2.eq.36))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost)
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 4
!!$               else if((isit1.eq.36).and.(isit2.eq.36))then
!!$
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost) 
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 4
!!$
!!$               else if((isit1.eq.36).and.(isit2.eq.1))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost)  
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost)/ 4
!!$               else if((isit1.eq.1).and.((isit2.ne.1).or.(isit2.ne.36)))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost)
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 6
!!$          else if((isit1.eq.36).and.((isit2.ne.1).or.(isit2.ne.36)))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost) 
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 6
!!$          else if((isit2.eq.1) .and. ((isit1.ne.1).or.(isit1.ne.36)))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost)
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 6
!!$          else if((isit2.eq.36).and.((isit1.ne.1).or.(isit1.ne.36)))then
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost)
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 6
!!$               else 
!!$                 rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2+1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1-1,isit2,ipre,iaa,ipost) +  &
!!$                  rama_prob(isit1+1,isit2,ipre,iaa,ipost)
!!$
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) = &
!!$                  rama_prob(isit1,isit2,ipre,iaa,ipost) / 9
!!$               endif
!!$
!!$              enddo  ! do isit1 = 1,36
!!$            enddo ! do isit2 = 1,36
!!$          enddo ! do ipost = 1,20
!!$         enddo ! do ipre = 1,20
!!$        enddo ! do iaa = 1,20
!!$
!!$! calc rama_force
!!$
!!$          do iaa = 1,20
!!$           do ipre = 1,20
!!$            do ipost = 1,20
!!$             do isit1 = 1,36
!!$              do isit2 = 1,36
!!$
!!$             if ( isit2 .eq. 1) then
!!$               rama_force(isit1,isit2,ipre,iaa,ipost) =  & 
!!$                  (rama_prob(isit1,isit2,ipre,iaa,ipost) - & 
!!$                   rama_prob(isit1,36,ipre,iaa,ipost) )/2
!!$             else if ( isit2 .eq. 1) then 
!!$               rama_force(isit1,isit2,ipre,iaa,ipost) = &
!!$                  (rama_prob(isit1,isit2,ipre,iaa,ipost) - &
!!$                   rama_prob(36,isit2,ipre,iaa,ipost) )/2
!!$             else 
!!$               rama_force(isit1,isit2,ipre,iaa,ipost) = &
!!$                 (rama_prob(isit1,isit2,ipre,iaa,ipost) - & 
!!$                  rama_prob(isit1,isit2-1,ipre,iaa,ipost) )/2
!!$             endif
!!$
!!$              enddo ! do isit2 = 1,36
!!$             enddo ! do isit1 = 1,36
!!$            enddo ! do ipost = 1,20
!!$           enddo ! do ipre = 1,20
!!$          enddo ! do iaa = 1,20

        return
	end
