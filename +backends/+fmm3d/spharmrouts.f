C     Stuff from Leslie/Zydrunas FMMLIB3D, and localexp3d.
C     pulled by Barnett 7/31/15
C     modified by Greengard 8/1/15


C***********************************************************************
      subroutine projloc3d(nterms,ldl,nquad,nquadm,xnodes,wts,
     1           phival,local)
C***********************************************************************
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl    = dimension parameter for local expansion
C           nquad  = number of quadrature nodes in theta direction.
C           nquadm = number of quadrature nodes in phi direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(j,i) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           NOTE:    We assume phi_i = i*2*pi/nquadm, as do the 
C                    routines in projection.f. However, we permit
C                    different numbers of nodes in theta and phi.
C***********************************************************************
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C     NOTE:
C
C     yrecursion.f produces Ynm with a nonstandard scaling:
C     (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C     is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C     In the first loop below, you see
C
Cccc	    marray(jj,m) = sum*2*pi/nquadm
C	    marray(jj,m) = sum/(2*nquadm)
C
C     The latter has incorporated the 1/(4*pi) normalization factor
C     into the azimuthal quadrature weight (2*pi/nquadm).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 zk,phival(nquadm,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(marray(nquad,-nterms:nterms))
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
            local(l,m) = 0.0d0
         enddo
      enddo
c
c     create marray (intermediate array)
c
      do m=-nterms,nterms
         emul = cdexp(imag*m*2*pi/nquadm)
         do jj=1,nquad
            sum = 0
            ephi = 1.0d0
            do kk = 1,nquadm
ccc               ephi = cdexp(imag*m*(kk-1)*2*pi/nquadm)
               sum = sum + phival(kk,jj)*dconjg(ephi)
               ephi = ephi*emul
            enddo
ccc	    marray(jj,m) = sum*2*pi/nquadm
            marray(jj,m) = sum/(2*nquadm)
         enddo
         emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquad
         cthetaj = xnodes(jj)
         call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
            zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end



C     from localexp3d/h3dwrappers.f =====================================

c**********************************************************************
      subroutine flattenlocexpz(nterms,in,out,stride)
c     unpacks a complex *16 2d local expansion array (in) into a 1d list (out).
c     Output array is written to using given stride.
c     Alex Barnett 3/28/12
      implicit none
      integer i,n,m,nterms,stride
      complex *16 in(0:nterms,-nterms:nterms), out(*)
      
      i = 1
      do n=0,nterms
         do m=-n,n
            out(i) = in(n,m)
            i = i+stride
         enddo
      enddo
      end

c**********************************************************************
      subroutine stacklocexpz(nterms,in,out,stride)
c     packs a complex *16 1d local expansion list (in) into a 2d array (out).
c     Input array is read using given stride.
c     Alex Barnett 2/4/13
      implicit none
      integer i,n,m,nterms,stride
      complex *16 out(0:nterms,-nterms:nterms), in((nterms+1)**2)
      
      i = 1
      do n=0,nterms
         do m=-n,n
            out(n,m) = in(i)
            i = i+stride
         enddo
      enddo
      end
C
C
C
C***********************************************************************
      subroutine shevalsphere(local,phival,
     1           nterms,lmp,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion on an 
C     nquad x nquadm grid on the unit sphere. 
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local    : coefficients of spherical harmonic exp.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for local
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi,ephi1,ephik
C
      integer l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquad))
c
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(m,jj) = 0.0d0
      enddo
      enddo
c
      do jj=1,nquad
         ctheta = xnodes(jj)
         stheta = dsqrt(1.0d0 - ctheta**2)
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,jj) = phitemp(m,jj) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
c
      do jj = 1,nquad
      do kk = 1,nquadm
         phival(kk,jj) = 0.0d0
         ephik = cdexp(2*pi*(kk-1)*imag/nquadm)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(kk,jj) = phival(kk,jj) + phitemp(m,jj)*ephi
            ephi = ephi*ephik
         enddo
      enddo
      enddo
      return
      end
C
C
C
C
C***********************************************************************
      subroutine shevalspherecircs(local,phival,
     1           nterms,lmp,nquad,nquadm,xnodes,nquadc)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion at
C     the equatorial circles for z-axis rotation to all standard spherical
C     nodes (on an nquad x nquadm grid). nquadc is the number of points
C     desired on each great circle.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local    : coefficients of spherical harmonic exp.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for local
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C     nquadc   : number of points on equator for each rotation
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : (l,j,i) corresponds to lth node on equatorial circle
C                        for normal corresponding to (theta_i,phi_j).
C     NOTE the ordering of output grids.
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadc,nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:,:)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi,ephi1,ephik
C
      integer l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquadc,nquad))
c
      do ii=1,nquad
      do ll=1,nquadc
      do m=-nterms,nterms
         phitemp(m,ll,ii) = 0.0d0
      enddo
      enddo
      enddo
ccc      call prin2(' initalized phitemp *',pi,0)
c
      do ii=1,nquad
      do ll=1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         phil = 2*pi*ll/nquadc
         phil = 2*pi*(ll-1)/nquadc
         cphil = dcos(phil)
         ctheta = -salphai*cphil
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,ll,ii) = phitemp(m,ll,ii) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
      enddo
ccc      call prin2(' computed phitemp *',pi,0)
c
      do ii = 1,nquad
      do jj = 1,nquadm
      do ll = 1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         betaj = 2*pi*jj/nquadm
         betaj = 2*pi*(jj-1)/nquadm
         cbetaj = dcos(betaj)
         sbetaj = dsin(betaj)
ccc         phil = 2*pi*ll/nquadc
         phil = 2*pi*(ll-1)/nquadc
         cphil = dcos(phil)
         sphil = dsin(phil)
         x = cbetaj*calphai*cphil-sbetaj*sphil
         y = (sbetaj*calphai*cphil+cbetaj*sphil)
ccc         zg = -salphai*cphil
ccc         if ( ((ii.eq.1).and.(jj.eq.1)).and.(ll.eq.1)) then
ccc         call prin2('111 x is *', x,1)
ccc         call prin2('111 y is *', y,1)
ccc         call prin2('111 zg is *', zg,1)
ccc         endif
         phi = datan2(y,x)
c         
         phival(ll,jj,ii) = 0.0d0
         ephik = cdexp(imag*phi)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(ll,jj,ii) = phival(ll,jj,ii) + phitemp(m,ll,ii)*ephi
            ephi = ephi*ephik
         enddo
ccc         if ( ((ii.eq.1).and.(jj.eq.1)).and.(ll.eq.1)) then
ccc         call prin2('111 phival is *', phival(ll,jj,ii),2)
ccc         endif
      enddo
      enddo
      enddo
ccc      call prin2(' computed phival *',pi,0)
      return
      end
C
c
c
      subroutine mkgreatcircles(rk,ngridt,xnodesth,ngridp,xnodesph,
     1           ngridc,xg,yg,zg)
c
c     create Cartesian coordinates for equatorial circles on sphere of
c     radius rk with z-axis rotated to location 
c     (xnodesth(i),xnodesph(j)).
c
c     INPUT:
c
c     rk         radius of sphere
c     ngridt     number of discretization nodes in theta in lab frame
c     xnodesth   discretization nodes in theta in lab frame
c     ngridp     number of discretization nodes in phi in lab frame
c     xnodesph   discretization nodes in phi in lab frame
c     ngridc     number of nodes on great circle
c
c     OUTPUT:
c
c     xg(l,j,i)  xcoord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     yg(l,j,i)  y-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     zg(l,j,i)  z-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c
c     NOTE the ordering of output grids.
c
      implicit none
      integer i,j,l,ngridc,ngridt,ngridp
      real *8 pi,rk,alphai,calphai,salphai
      real *8 betaj,cbetaj,sbetaj,phil,cphil,sphil
      real *8 xnodesth(ngridt)
      real *8 xnodesph(ngridp)
      real *8 xg(ngridc,ngridp,ngridt)
      real *8 yg(ngridc,ngridp,ngridt)
      real *8 zg(ngridc,ngridp,ngridt)
c
      pi = 4.0d0*datan(1.0d0)
c
      do i = 1,ngridt
         calphai = xnodesth(i)
         salphai = dsqrt(1 - calphai**2)
         do j = 1,ngridp
            betaj = xnodesph(j)
            cbetaj = dcos(betaj)
            sbetaj = dsin(betaj)
            do l = 1,ngridc
ccc               phil = 2*pi*l/ngridc
               phil = 2*pi*(l-1)/ngridc
               cphil = dcos(phil)
               sphil = dsin(phil)
               xg(l,j,i) = rk*(cbetaj*calphai*cphil-sbetaj*sphil)
               yg(l,j,i) = rk*(sbetaj*calphai*cphil+cbetaj*sphil)
               zg(l,j,i) = rk*(-salphai*cphil)
            enddo
         enddo
      enddo
      return
      end
