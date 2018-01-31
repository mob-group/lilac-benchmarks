!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! Checked trigonometric functions

      function checked_acos(x)

         use msevb_common

         implicit none

         DOUBLE PRECISION CHECKED_ACOS, X

         if (x<minusTrigUnity) then 
            checked_acos = pi
         elseif (x>trigUnity) then 
            checked_acos = 0.0d0
         else
            checked_acos = DACOS(x)
         endif      
         
         return 

      end function
