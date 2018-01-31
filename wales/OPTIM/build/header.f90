MODULE OPTIMHEADER
     implicit none
     contains
          subroutine PrintHeader
               implicit none

               write(*,'(/a)') ' OPTIM version Unversioned directory, Copyright (C) David J. Wales'
               write(*,'(a)') ' OPTIM comes with ABSOLUTELY NO WARRANTY; for details supply WARRANTY as an input keyword.'
               write(*,'(a)') ' This is free software, and you are welcome to redistribute it'
               write(*,'(a/)') ' under certain conditions; provide keyword COPYRIGHT to see the details.'
          end subroutine PrintHeader
END MODULE OPTIMHEADER
