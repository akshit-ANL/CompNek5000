c---------------------------------------------------------------------
      subroutine compress_cfl
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real cfl,acfl,dcfl
      integer icalld
      save    icalld
      data    icalld /0/
      real    mdr
      save    mdr

      if (icalld.eq.0) call mindr(mdr)
      icalld = 1
      
      call cmp_acfl(acfl,aspd,dt,mdr) ! acoustic CFL
      call cmp_dcfl(dcfl,enmu,dt,mdr) ! diffusive CFL
      call compute_cfl(cfl,abs(vx)+aspd,abs(vy)+aspd,abs(vz)+aspd,dt) 

      if(nid.eq.0) then
        write(6,77) cfl,acfl,dcfl
 77     format('                                      C=',0pF7.3
     &        ,', aC=',0pF7.3,', dC=',0pF7.3)
      endif

      return
      end
c--------------------------------------------------------------------
      subroutine cmp_acfl(cfl,va,dt,mdr)  ! compute acoustic cfl

c     Given speed of sound, compute cfl

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'WZ'
      include 'SOLN'
      include 'COMPRESS'

      real dt, cfl, mdr, imd
      real va(lx1,ly1,lz1,lelt)
      integer  e

      if(mdr.gt.1.e-8) then
        imd = 1./mdr
      else
        write(6,*) 'min distance too small, mdr= ', mdr
        call exitti('min dist too small$',ndim)
      endif

      cfl = 1.e-30

      do e=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
        ctmp = (va(ix,iy,iz,e)+spd(ix,iy,iz,e))*dt*real(lx1-2)
     $         /enhe(ix,iy,iz,e)
        if(cfl.le.ctmp) cfl=ctmp
      enddo
      enddo
      enddo
      enddo

      cfl = glmax(cfl,1)

      return
      end
c--------------------------------------------------------------------
      subroutine cmp_dcfl(cfl,va,dt,mdr)  ! compute diffusive cfl

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'WZ'
      include 'SOLN'
      include 'COMPRESS'

      real dt, cfl, mdr, imd
      real va(lx1,ly1,lz1,lelt)
      integer  e

      if(mdr.gt.1.e-8) then
        imd = 1./mdr
      else
        write(6,*) 'min distance too small, mdr= ', mdr
        call exitti('min dist too small$',ndim)
      endif

      cfl = 1.e-30

      do e=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
        ctmp = (va(ix,iy,iz,e))*dt*(real(lx1-2)*real(lx1-2))
     $         /(enhe(ix,iy,iz,e)*enhe(ix,iy,iz,e))
        if(cfl.le.ctmp) cfl=ctmp
      enddo
      enddo
      enddo
      enddo
      
      cfl = glmax(cfl,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine mindr(mdr) 

c     Find minimum distance between grid points

      include 'SIZE'
      include 'TOTAL'

      real    mdr, dr1, dx, dy
      real    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4
      real    x5,x6,y5,y6,z0,z1,z2,z3,z4,z5,z6
      integer e

      mdr=1.e5
      if(if3d) then
        do e=1,nelt
          do iz=2,nz1-1
          do iy=2,ny1-1
          do ix=2,nx1-1
            x0 = xm1(ix  ,iy  ,iz  ,e)
            x1 = xm1(ix  ,iy-1,iz  ,e)
            x2 = xm1(ix+1,iy  ,iz  ,e)
            x3 = xm1(ix  ,iy+1,iz  ,e)
            x4 = xm1(ix-1,iy  ,iz  ,e)
            x5 = xm1(ix  ,iy  ,iz-1,e)
            x6 = xm1(ix  ,iy  ,iz+1,e)
            y0 = ym1(ix  ,iy  ,iz  ,e)
            y1 = ym1(ix  ,iy-1,iz  ,e)
            y2 = ym1(ix+1,iy  ,iz  ,e)
            y3 = ym1(ix  ,iy+1,iz  ,e)
            y4 = ym1(ix-1,iy  ,iz  ,e)
            y5 = ym1(ix  ,iy  ,iz-1,e)
            y6 = ym1(ix  ,iy  ,iz+1,e)
            z0 = zm1(ix  ,iy  ,iz  ,e)
            z1 = zm1(ix  ,iy-1,iz  ,e)
            z2 = zm1(ix+1,iy  ,iz  ,e)
            z3 = zm1(ix  ,iy+1,iz  ,e)
            z4 = zm1(ix-1,iy  ,iz  ,e)
            z5 = zm1(ix  ,iy  ,iz-1,e)
            z6 = zm1(ix  ,iy  ,iz+1,e)
            dr1 = dist3(x0,y0,z0,x1,y1,z1)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist3(x0,y0,z0,x2,y2,z2)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist3(x0,y0,z0,x3,y3,z3)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist3(x0,y0,z0,x4,y4,z4)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist3(x0,y0,z0,x5,y5,z5)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist3(x0,y0,z0,x6,y6,z6)
            if(dr1.lt.mdr) mdr=dr1
          enddo
          enddo
          enddo
        enddo
      else
        do e=1,nelt
          do iy=2,ny1-1
          do ix=2,nx1-1
            x0 = xm1(ix  ,iy  ,1,e)
            x1 = xm1(ix  ,iy-1,1,e)
            x2 = xm1(ix+1,iy  ,1,e)
            x3 = xm1(ix  ,iy+1,1,e)
            x4 = xm1(ix-1,iy  ,1,e)
            y0 = ym1(ix  ,iy  ,1,e)
            y1 = ym1(ix  ,iy-1,1,e)
            y2 = ym1(ix+1,iy  ,1,e)
            y3 = ym1(ix  ,iy+1,1,e)
            y4 = ym1(ix-1,iy  ,1,e)
            dr1 = dist2(x0,y0,x1,y1)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist2(x0,y0,x2,y2)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist2(x0,y0,x3,y3)
            if(dr1.lt.mdr) mdr=dr1
            dr1 = dist2(x0,y0,x4,y4)
            if(dr1.lt.mdr) mdr=dr1
          enddo
          enddo
        enddo
      endif
      if(lx1.eq.2) then
        do e=1,nelt
          dx = xm1(2,1,1,e) - xm1(1,1,1,e)
          dy = ym1(1,2,1,e) - ym1(1,1,1,e)
          mdr = min(dy,dx)
        enddo
      endif
      mdr = glmin(mdr,1)
      if(mdr.ge.1.e5.or.mdr.le.1.e-8) then
        write(6,*) 'wrong mdr',mdr,' ,nid',nid
        call exitti('min dist not right$',nid)
      endif

      return
      end
c-----------------------------------------------------------------------
      real function dist3(x1,y1,z1,x2,y2,z2)
      real x1,y1,z1,x2,y2,z2,dx,dy,dz

      dx = x1-x2
      dy = y1-y2
      dz = z1-z2
      dist3 = sqrt(dx*dx+dy*dy+dz*dz)

      return
      end
c-----------------------------------------------------------------------
      real function dist2(x1,y1,x2,y2)
      real x1,y1,x2,y2,dx,dy

      dx = x1-x2
      dy = y1-y2
      dist2 = sqrt(dx*dx+dy*dy)

      return
      end
c---------------------------------------------------------------------      
