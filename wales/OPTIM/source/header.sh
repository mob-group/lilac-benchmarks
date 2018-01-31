#!/bin/bash
#
#   OptimHeader module generator for easier version control.
#   Copyright (C) 2003-2005 Semen A. Trygubenko and David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
sourcedir=`dirname $0`

echo "MODULE OPTIMHEADER"
echo "     implicit none"
echo "     contains"
echo "          subroutine PrintHeader"
echo "               implicit none"
echo
echo "               write(*,'(/a)') ' OPTIM version `bash $sourcedir/version.sh`, Copyright (C) David J. Wales'"
echo "               write(*,'(a)') ' OPTIM comes with ABSOLUTELY NO WARRANTY; for details supply WARRANTY as an input keyword.'"
echo "               write(*,'(a)') ' This is free software, and you are welcome to redistribute it'"
echo "               write(*,'(a/)') ' under certain conditions; provide keyword COPYRIGHT to see the details.'"
echo "          end subroutine PrintHeader"
echo "END MODULE OPTIMHEADER"
