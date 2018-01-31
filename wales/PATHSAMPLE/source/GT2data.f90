!
!  Copyright (C) 2005-2006 Semen Trygubenko
!  This file is part of GT.
!
!  GT is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  GT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with GT; see the file COPYING. If not, write to the
!  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
!  MA 02111-1307 USA
!
module DataModule
     implicit none
     save
     integer :: Nnodes,Nsources,Nsinks,Ninodes,NunprocessedNodes
     double precision,dimension(:),allocatable :: tau
     double precision,dimension(:,:),allocatable :: P
     character(len=20) :: KeywordsFile="GTkeywords",InputFile="GTinput",s1,s2
     type Node
          integer :: degree,d,index,indexorig ! degree is for FH internal use; d is the node degree
          type(NodeList),pointer :: A
          type(RealList),pointer :: P
          type(Node),pointer :: left,right,parent,child
          double precision :: tau,prob
          logical :: FH,mark
          character(len=1) :: t ! node type: t - sink, s - source, i - intermediate, d - disconnected, p - processed
     end type Node
     type IntegerList
          type(IntegerList),pointer :: f,b
          integer :: value
     end type IntegerList
     type RealList
          type(RealList),pointer :: f,b
          double precision :: value
     end type RealList
     type NodeList
          type(NodeList),pointer :: f,b
          type(Node),pointer :: value
     end type NodeList
     type NodeArray
          type(Node),pointer :: p
     end type NodeArray
     type(NodeArray),allocatable :: G(:)
end module DataModule
