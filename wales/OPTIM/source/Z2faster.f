c     This is an untested wrapper to my z2-routine that
c     uses linked cells to efficiently find which atoms
c     are within interaction range of a given atom.

      subroutine z2faster(np,xx,Lx,Ly,Lz,etot,gg,STEST)

c     Fairly efficient Z2 evaluator for large systems
c     INPUT
c       np -- number of particles
c       xx -- coordinate vector xx({1,2,3},i) = {x,y,z}(i)
c              Ok to call with xx of dimension xx(3*np),
c              in which case {x,y,z}(i) = xx(3*(i-1) + {1,2,3})
c       Lx -- Box length in x-direction
c       Ly -- Box length in y-direction
c       Lz -- Box length in z-direction
c     OUTPUT
c       etot -- Energy of configuration
c       gg   -- Force (same dimensions as xx)
c
      implicit none
      integer np
      real*8 xx(3,np),Lx,Ly,Lz,etot,gg(3,np)
      LOGICAL STEST

      real*8 RCUT,e,de,boxl(2,3)
      integer i,nx,ny,nz
      parameter ( RCUT =  2.644877840738571D0   )

      IF (STEST) THEN
         PRINT '(A)','z2faster> second derivatives not coded'
         STOP
      ENDIF

      nx = max(3,int(Lx/(RCUT+0.01d0)))
      ny = max(3,int(Ly/(RCUT+0.01d0)))
      nz = max(3,int(Lz/(RCUT+0.01d0)))
      boxl(1,1) = Lx
      boxl(1,2) = Ly
      boxl(1,3) = Lz
      do i=1,3
         boxl(2,i) = 1.0d0 / boxl(1,i)
      enddo

      call fgcalc(np,nx,ny,nz,xx,boxl,e,gg,de)
      etot = e+de
      end


c     Divide space in boxes, and compute which particle belong
c     to which box. Boxes are supposedly have sides .ge.
c     potential reange, so only neighboring boxes need be considered
c     for potential evaluation.
      subroutine build_boxes(np,nx,ny,nz,xx,boxl,first,next)
c     Divide the unit cube simulation box into nx*ny*nz cells,
c     and put all particles into their respective cell.
      implicit none
      integer np,nx, ny,nz,i,j,k,nb(3),idx(3)
      integer first(nx,ny,nz), next(np)
      real*8 xx(3,np),boxl(2,3)

      nb(1) = nx
      nb(2) = ny
      nb(3) = nz

      do i=1,np
         next(i) = 0
      enddo

      do i=1,nb(3)
         do j=1,nb(2)
            do k=1,nb(1)
               first(k,j,i) = 0
            enddo
         enddo
      enddo
      do i=1,np
         do j=1,3
            idx(j) = int(xx(j,i)*boxl(2,j)*nb(j))
            if(xx(j,i) .lt. 0d0) then
               idx(j) = 1 + mod(idx(j)-1-nb(j)*(idx(j)-1),nb(j))
            else
               idx(j) = 1 + mod(idx(j),nb(j))
            endif
            if(idx(j).lt.1 .or. idx(j).gt.nb(j)) then
               write(*,*) i,j,xx(j,i),idx(j),nb(j),boxl(1,j),boxl(2,j)
               stop
            endif
         enddo
         next(i) = first(idx(1),idx(2),idx(3))
         first(idx(1),idx(2),idx(3)) = i
      enddo
      end

c     Distance routine, takes periodic boundary conditions into account
      subroutine z2dist(xx,yy,boxl,dr)
      implicit none
      real*8 xx(3),yy(3),boxl(2,3),dr(3)

      integer i
      real*8 t

c     Compute distance between particles p and q
      do i=1,3
         t = yy(i) - xx(i)
         dr(i) = t - boxl(1,i)*nint(t*boxl(2,i))
      enddo
      end      

c     g77 does not have ceiling implemented; this can be substituted
c      real*8 function ceiling(x)
c      implicit none
c      real*8 x,y
c      y = int(x)
c      if(y.ne.x .and. x.ge.0d0) y = y + 1d0
c      ceiling = y
c      end

      real*8 function epair(x,y,boxl)
      implicit none
      real*8 x(3),y(3),boxl(2,3),dr(3),r
      real*8 rinv,s1,s2,s3,A,B,KF,ALPH,SIG,POW,ST,RCUT,rcut2,twoK,Bsig
      parameter ( A    =  1.04D0      )
      parameter ( B    =  4200000.0D0 )
      parameter ( KF   =  4.139D0     )
      parameter ( ALPH =  0.33D0      )
      parameter ( SIG  =  0.348D0     )
      parameter ( POW  = -14.5D0      )
c     parameter ( ST=0.133915D0 ) ! Values below optimized by David
c     parameter ( RCUT=2.645D0  ) ! Wales for exact zero at RCUT
      parameter ( ST   =  0.1339154253770228D0  )
      parameter ( RCUT =  2.644877840738571D0   )
      parameter ( rcut2 =  RCUT**2              )
      parameter ( twoK  =  2.0d0 * KF           )
c     parameter ( Bsig = B/SIG**POW             )
      parameter ( Bsig  =  0.94656039153728358d0)

      call z2dist(x,y,boxl,dr)
      r = dr(1)**2 + dr(2)**2 + dr(3)**2
      if(r .lt. rcut2) then
c     Analytical potential
         r = sqrt(r)
         rinv = 1.0d0 / r
         s1 = A * exp(ALPH * r) * rinv**3
         s2 = cos(twoK * r)
         s3 = r**POW
         epair = s1*s2 + Bsig*s3 + ST
      else
         epair = 0d0
      endif
      end


c     Compute energy and gradient. Potential energy is
c     e+de on output. Typically, de/e <~ 1e-16, so this
c     sum holds energy with more resolution than a single
c     64-but floating point number.
c     Algorithm:
c     - Divide space in boxes, and distribute particles.
c     - Go through boxes, and use current box and neighbors
c       to find all particles interacting with particles in
c       current box.
c     - Sum energy contributions with Kahan summation, producing
c       higher precision. For this to work, the code need be
c       compiled without aggressive optimizatino (no changes of
c       semantics, reordering of calculations allowed). Also,
c       floating point result must be stored back in memory after
c       each expression, so that excess precision potentially
c       available in the cpu is gotten rid of before next expression
c       is evaluated. For g77, this is enabled by -ffloat-store,
c       and for the Intel fortran compiler with -fltconsistency
c       Usually optimization level -O2 is safe.
      subroutine fgcalc(np,nx,ny,nz,xx,boxl,e,gg,de)
      implicit none
      integer np,nx,ny,nz
      real*8 xx(3,np),boxl(2,3),e,gg(3,np),de
      
      integer first(nx,ny,nz),next(np),nb(3),divec(3,13)
      integer i,j,k,p,q,ii,i2,j2,k2,l
      real*8 v,r,gij,dr(3),dv,v1,t,h

c     Uncomment if g77, also uncomment ceiling implementation above
c      external ceiling
c      real*8 ceiling

c     For use when computing potential by table lookup
c      real*8 f
c      integer idx
c      include 'tz2c.h'


      real*8 rinv,s1,s2,s3,A,B,KF,ALPH,SIG,POW,ST,RCUT,rcut2,twoK,Bsig
      parameter ( A    =  1.04D0      )
      parameter ( B    =  4200000.0D0 )
      parameter ( KF   =  4.139D0     )
      parameter ( ALPH =  0.33D0      )
      parameter ( SIG  =  0.348D0     )
      parameter ( POW  = -14.5D0      )
c     parameter ( ST=0.133915D0 ) ! Values below optimized by David
c     parameter ( RCUT=2.645D0  ) ! Wales for exact zero at RCUT
      parameter ( ST   =  0.1339154253770228D0  )
      parameter ( RCUT =  2.644877840738571D0   )

      parameter ( rcut2 =  RCUT**2              )
      parameter ( twoK  =  2.0d0 * KF           )
c     parameter ( Bsig = B/SIG**POW             )
      parameter ( Bsig  =  0.94656039153728358d0)


c     Offsets to 13 uppleft boxes surrounding each cell.
      data divec / -1,-1,-1 , -1,-1, 0 , -1,-1, 1 , -1, 0,-1
     :     ,       -1, 0, 0 , -1, 0, 1 , -1, 1,-1 , -1, 1, 0
     :     ,       -1, 1, 1 ,  0,-1,-1 ,  0,-1, 0 ,  0,-1, 1
     :     ,        0, 0,-1 /

      nb(1) = nx
      nb(2) = ny
      nb(3) = nz
c     To ensure there are no negative results in modulus
c     computations later on, and to accomodate to the fact
c     that fortran indices start with one. divec is restored
c     to its original state upon exit...
      do i=1,13
         do j=1,3
            divec(j,i) = divec(j,i) + nb(j) - 1
         enddo
      enddo

      do i=1,np
         do j=1,3
            gg(j,i) = 0d0
         enddo
      enddo
      v = 0d0
      dv = 0d0

      call build_boxes(np,nb(1),nb(2),nb(3),xx,boxl,first,next)
      
      do i=1,nb(3)
         do j=1,nb(2)
            do k=1,nb(1)
               if(first(k,j,i) .gt. 0) then
c     BEGIN Pairs of particles in different boxes

c     Loop over upper left 13 of the 26 surrounding boxes
                  do ii=1,13
                     i2 = mod(i+divec(3,ii),nb(3)) + 1
                     j2 = mod(j+divec(2,ii),nb(2)) + 1
                     k2 = mod(k+divec(1,ii),nb(1)) + 1

                     if(i.eq.i2.and.j.eq.j2.and.k.eq.k2) then
                        write(*,*)
     :                       'newstweb.f: Box double counting: ',
     :                       i,j,k,ii
                        stop
                     endif

                     if(first(k2,j2,i2) .gt. 0) then
                        p = first(k,j,i)
                        do while(p .gt. 0)
                           q = first(k2,j2,i2)
                           do while(q .gt. 0)
                              call z2dist(xx(1,p),xx(1,q),boxl,dr)
                              r = dr(1)**2 + dr(2)**2 + dr(3)**2
                              if(r .lt. rcut2) then
c                             Analytical potential
                                 r = sqrt(r)
                                 rinv = 1.0d0 / r
                                 s1 = A * exp(ALPH * r) * rinv**3
                                 s2 = cos(twoK * r)
                                 s3 = r**POW
               
c                             Contribution to potential energy phi(r)
c     v = v + s1*s2 + Bsig*s3 + ST
                                 h = s1*s2 + Bsig*s3 + ST
c     Kahan summation, broken up so -ffloat-store can do its job
c     v1 = v+h+dv;
c     d  = -(v1-v)+h+dv;
c     v = v1;
                                 t = v+h
                                 v1 = t+dv
                                 t = -(v1-v)
                                 t = t+h
                                 dv = t+dv
                                 v = v1
c                             Contribution to gradient 1/r * Dphi(r)/Dr
                                 gij = (s1*((ALPH-3.0d0*rinv)*s2 -
     :                                twoK*sin(twoK*r))
     :                                + POW*Bsig*s3*rinv)*rinv

c!c                             Potential evaluated from table 
c!c                             Contribution to potential energy phi(r)
c!                                 idx = ceiling(r*idrtab);
c!                                 f = r*idrtab - idx;
c!                                 v = v + pot(idx) + f*pot1(idx)
c!
c!c                             Contribution to gradient 1/r * Dphi(r)/Dr
c!                                 gij = dpot(idx) + f*dpot1(idx)
                                 do l=1,3
                                    gg(l,p) = gg(l,p) - gij*dr(l)
                                    gg(l,q) = gg(l,q) + gij*dr(l)
                                 enddo
                              endif
                              q = next(q)
                           enddo
                           p = next(p)
                        enddo
                     endif
                  enddo
c     END

c     BEGIN Particle pairs within box k,j,i
                  p = first(k,j,i)
                  do while(p .gt. 0)
                     q = next(p)
                     do while(q .gt. 0)
                        call z2dist(xx(1,p),xx(1,q),boxl,dr)
                        r = dr(1)**2 + dr(2)**2 + dr(3)**2
                        if(r .lt. rcut2) then
c                          Analytical potential
                           r = sqrt(r)
                           rinv = 1.0d0 / r
                           s1 = A * exp(ALPH * r) * rinv**3
                           s2 = cos(twoK * r)
                           s3 = r**POW
               
c                         Contribution to potential energy phi(r)
c     v = v + s1*s2 + Bsig*s3 + ST
                           h = s1*s2 + Bsig*s3 + ST
c     Kahan summation, broken up so -ffloat-store can do its job
c     v1 = v+h+dv;
c     d  = -(v1-v)+h+dv;
c     v = v1;
                           t = v+h
                           v1 = t+dv
                           t = -(v1-v)
                           t = t+h
                           dv = t+dv
                           v = v1

c                         Contribution to gradient 1/r * Dphi(r)/Dr
                                 gij = (s1*((ALPH-3.0d0*rinv)*s2 -
     :                                twoK*sin(twoK*r))
     :                                + POW*Bsig*s3*rinv)*rinv
c!c                       Potential evaluated from table 
c!c                       Contribution to potential energy phi(r)
c!                           idx = ceiling(r*idrtab);
c!                           f = r*idrtab - idx;
c!                           v = v + pot(idx) + f*pot1(idx)
c!                           
c!c                          Contribution to gradient 1/r * Dphi(r)/Dr
c!                           gij = dpot(idx) + f*dpot1(idx)
                           do l=1,3
                              gg(l,p) = gg(l,p) - gij*dr(l)
                              gg(l,q) = gg(l,q) + gij*dr(l)
                           enddo
                        endif
                        q = next(q)
                     enddo
                     p = next(p)
                  enddo
c     END
               endif
            enddo
         enddo
      enddo

c     v1 = v+dv
c     dv  = -(v1-v)+dv
c     v  = v1
      v1 = v+dv
      t = v1-v
      dv = -t+dv
      v = v1

      e = v
      de = dv

c     Restore to orginal state, since divec might be saved
c     between calls and not reinitialized on the next entry...
      do i=1,13
         do j=1,3
            divec(j,i) = divec(j,i) - nb(j) + 1
         enddo
      enddo
      end
