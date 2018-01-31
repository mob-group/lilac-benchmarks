MODULE PORFUNCS
     implicit none
     contains
          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be
               implicit none ! connected for formatted sequential output; ISTAT is ignored
               integer,intent(in) :: UNIT
               integer,intent(out),optional :: ISTAT
!              call flush(UNIT)
          end subroutine flush
 
          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated
               implicit none
 
               integer,intent(in) :: position
               character(len=*),intent(out) :: value
 
 
               call getarg(position,value)
          end subroutine getarg_subr
 
          subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated
               implicit none
               integer,intent(out) :: n
 
               integer iargc
 
               n = iargc()
          end subroutine iargc_subr
 
          subroutine fork_subr(pid) ! returns zero in the child process, PID of child in parent process
               implicit none
               integer, intent(inout) :: pid
 
 
          end subroutine fork_subr
 
          subroutine system_subr(JobString,ExitStatus)
               implicit none
 
               character(len=*),intent(in) :: JobString
               integer,intent(out) :: ExitStatus
 
          end subroutine system_subr
 
          subroutine wait_subr(pid,ExitStatus)
               implicit none
 
               integer,intent(inout) :: pid,ExitStatus
 
          end subroutine wait_subr
 
          subroutine getpid_subr(pid)
               implicit none
 
               integer,intent(out) :: pid
 
               pid=getpid()
          end subroutine getpid_subr
END MODULE PORFUNCS
