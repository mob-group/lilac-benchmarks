!  vr274
! this modules defines a wrapper interface for internal coordinates function
! to reduce/simplify the if blocks in the core routines, e.g. mylbfgs
module internals_wrapper
! implicit none
    
    private :: intwrap_kd, intwrap_nnz
    integer :: intwrap_kd = -1
    integer :: intwrap_nnz = -1

contains

    logical function intwrap_useinternals()
        use modcharmm, only : chrmmt, INTMINT
        USE INTCOMMONS, ONLY : NATINT
        use key, only : AMBERT, NABT
        implicit none
        intwrap_useinternals = (AMBERT.OR.NABT.OR.CHRMMT).AND.(INTMINT)
        return
    end function
!   
!   \brief check initalization for wrapper
!
    subroutine intwrap_checkinitialized()
        use modcharmm, only : chrmmt
        USE INTCOMMONS, ONLY : NATINT
        use key, only : AMBERT, NABT
        implicit none
        if((AMBERT.OR.NABT).AND.(.NOT.NATINT)) then
            print *,"Only natural internals supported in AMBER, pleas use keyword NATINT"
            STOP
        END IF
        if(intwrap_kd == -1) then
            IF (CHRMMT) THEN 
                call GETKD(intwrap_kd) ! get width of sparse band in G matrix KD
                call GETNNZ(intwrap_nnz) ! get number of non-zero elements in B-matrix
            else if ((AMBERT.OR.NABT).AND.NATINT) then
                call AMB_GETKDNAT(intwrap_kd) ! get width of sparse band in G matrix KD
                call AMBGETNNZNAT(intwrap_nnz) ! get number of non-zero elements in B-matrix
            else
                print *,"ERROR: Internal coordinate wrapper interface not implemented for non-charmm interfaces"
                stop
             endif
        endif
    end subroutine

!    
!   \param xcart coordinates in cartesians
!   \param gcart gradient in cartesians
!   \param xint  coorinates in internals
!   \param gint  gradient in internals
!   \param nocoor false if xint should be filled with internal coordinates
!   \param nocoor false to transform the gradient
!
    subroutine intwrap_transform(xcart,gcart,xint,gint,nocoor,noderv)
        use commons, only : NATOMS, NINTS
        use modcharmm, only : chrmmt
        USE KEY, ONLY : INTEPSILON, AMBERT, NABT
        implicit none
        double precision xcart(3*NATOMS), gcart(3*NATOMS), xint(NINTS), gint(NINTS)
        logical, intent(in) :: nocoor, noderv

        call intwrap_checkinitialized()

        if (CHRMMT) then
            call TRANSFORM(xcart,gcart,xint,gint,NINTS,3*NATOMS,intwrap_nnz,NOCOOR,NODERV,intwrap_kd,INTEPSILON)
        else if(AMBERT.OR.NABT) then
            call AMBTRANSFORM(xcart,gcart,xint,gint,NINTS,3*NATOMS,intwrap_nnz,NOCOOR,NODERV,intwrap_kd,INTEPSILON)
        else
            print *,"ERROR: Internal coordinate wrapper interface not implemented for non-charmm interfaces"
            stop
        endif
    end subroutine

    subroutine intwrap_transback(xint_new,xint_old,xcart)
        use commons, only : NATOMS, NINTS
        use modcharmm, only : chrmmt
        USE KEY, ONLY : INTEPSILON, AMBERT
        implicit none
        double precision xcart(3*NATOMS), xint_new(NINTS), xint_old(NINTS)

        call intwrap_checkinitialized()

        if (CHRMMT) then
            call TRANSBACK(xint_new, xint_old, xcart, NINTS, 3*NATOMS,intwrap_nnz,intwrap_kd,INTEPSILON)
        else
            print *,"ERROR: Internal coordinate wrapper interface not implemented for non-charmm interfaces"
            stop
        endif
    end subroutine

!     subroutine intwrap_transdelta(delta_cart,delta_int,xcart)
!         use commons, only : NATOMS, NINTS
!         use modcharmm, only : chrmmt
!         USE KEY, ONLY : INTEPSILON
!         implicit none
!         double precision, xcart(3*NATOMS), gcart(3*NATOMS), xint(NINTS), gint(NINTS)
! 
!         call intwrap_checkinitialized()
! 
!         if (CHRMMT) then
!             call TRANSFORM(xcart,gcart,xint,gint,NINTS,3*NATOMS,intwrap_nnz,NOCOOR,NODERV,intwrap_kd,INTEPSILON)
!         else
!             print *,"ERROR: Internal coordinate wrapper interface not implemented for non-charmm interfaces"
!             stop
!         endif
!     end subroutine


!   \brief transform a step from internals to cartesians
!   \param delta_int      stepo in internals, input
!   \param delta_cart   step in cartesians, output
!   \param coords       current coordinates in cartesians
!   \param failed       true if iterative transformation failed
!   \param ptest2       I don't have a clue what this is for
    subroutine intwrap_transbackdelta(delta_int, delta_cart, coords, failed, ptest2)
        use commons, only : NATOMS, NINTS
        use modcharmm, only : chrmmt
        USE KEY, ONLY : INTEPSILON, AMBERT, NABT
        implicit none
        double precision delta_cart(3*NATOMS), coords(3*NATOMS), delta_int(NINTS)
        logical failed, ptest2
        
        call intwrap_checkinitialized()

        if (CHRMMT) then
            call TRANSBACKDELTA(delta_int,delta_cart,coords,NINTS,3*NATOMS,intwrap_nnz,intwrap_kd,failed,ptest2,INTEPSILON)
        else if (AMBERT.OR.NABT) then
            call AMB_TRANSBACKDELTA(delta_int,delta_cart,coords,NINTS,3*NATOMS,intwrap_nnz,intwrap_kd,failed,ptest2,INTEPSILON)
        else
            print *,"ERROR: Internal coordinate wrapper interface not implemented for non-charmm interfaces"
            stop
        endif
    end subroutine
end module
