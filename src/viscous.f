c---------------------------------------------------------------------
      subroutine viscous(rkRHS,rkU)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real rkU  (lx1,ly1,lz1,lelt,N_EQ)
      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      
      call eval_visc_U(rkRHS)
      call eval_visc_EY(rkRHS)

      return
      end
c---------------------------------------------------------------------
      subroutine eval_dudx
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      if(ifavisc.or.ifviscous) then
        call liftGrad(dudx,dudy,dudz,vx,uxpi,uxf)
        call liftGrad(dvdx,dvdy,dvdz,vy,uypi,uyf)
        if(if3d) call liftGrad(dwdx,dwdy,dwdz,vz,uzpi,uzf)
      endif

      return
      end
c---------------------------------------------------------------------  
      subroutine eval_visc_U(rkRHS)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      integer ntot,i,j,e,iface,k,lnf
      real dgFF,dgGG,dgHH

      real tau11(lx1,ly1,lz1,lelt)    ! tau_xx
      real tau12(lx1,ly1,lz1,lelt)    ! tau_xy
      real tau13(lx1,ly1,lz1,lelt)    ! tau_xz
      real tau22(lx1,ly1,lz1,lelt)    ! tau_yy
      real tau23(lx1,ly1,lz1,lelt)    ! tau_yz
      real tau33(lx1,ly1,lz1,lelt)    ! tau_zz
      
      real tauE1(lx1,ly1,lz1,lelt) 
      real tauE2(lx1,ly1,lz1,lelt)
      real tauE3(lx1,ly1,lz1,lelt) 
      
      ntot = nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt

ccc tau
      if(if3d) then
        do i=1,ntot
          tau11(i,1,1,1) = ( 4./3.*dudx(i,1,1,1)
     $                     - 2./3.*dvdy(i,1,1,1) 
     $                     - 2./3.*dwdz(i,1,1,1) ) * vis1(i,1,1,1)
          tau12(i,1,1,1) = (dudy(i,1,1,1)+dvdx(i,1,1,1))*vis1(i,1,1,1)
          tau13(i,1,1,1) = (dudz(i,1,1,1)+dwdx(i,1,1,1))*vis1(i,1,1,1)
          tau22(i,1,1,1) = ( 4./3.*dvdy(i,1,1,1)
     $                     - 2./3.*dudx(i,1,1,1) 
     $                     - 2./3.*dwdz(i,1,1,1) ) * vis1(i,1,1,1)
          tau23(i,1,1,1) = (dvdz(i,1,1,1)+dwdy(i,1,1,1))*vis1(i,1,1,1)
          tau33(i,1,1,1) = ( 4./3.*dwdz(i,1,1,1)
     $                     - 2./3.*dudx(i,1,1,1) 
     $                     - 2./3.*dvdy(i,1,1,1) ) * vis1(i,1,1,1)
        enddo
      else
        do i=1,ntot
          tau11(i,1,1,1) = ( 4./3.*dudx(i,1,1,1)
     $                     - 2./3.*dvdy(i,1,1,1) ) * vis1(i,1,1,1)
          tau12(i,1,1,1) = (dudy(i,1,1,1)+dvdx(i,1,1,1))*vis1(i,1,1,1)
          tau22(i,1,1,1) = ( 4./3.*dvdy(i,1,1,1)
     $                     - 2./3.*dudx(i,1,1,1) ) * vis1(i,1,1,1)
        enddo
      endif

ccc rhoUVW: div(tau)
      call divQ_weak(tau11,tau12,tau13,rkRHS(1,1,1,1,U_EQ),'U1')
      call divQ_weak(tau12,tau22,tau23,rkRHS(1,1,1,1,V_EQ),'U1')
      if(if3d) then
        call divQ_weak(tau13,tau23,tau33,rkRHS(1,1,1,1,W_EQ),'U1')
      endif

ccc v.tau
      if(if3d) then
        do i=1,ntot
          tauE1(i,1,1,1) = vx(i,1,1,1)*tau11(i,1,1,1)
     $                   + vy(i,1,1,1)*tau12(i,1,1,1)
     $                   + vz(i,1,1,1)*tau13(i,1,1,1)
          tauE2(i,1,1,1) = vx(i,1,1,1)*tau12(i,1,1,1)
     $                   + vy(i,1,1,1)*tau22(i,1,1,1)
     $                   + vz(i,1,1,1)*tau23(i,1,1,1)
          tauE3(i,1,1,1) = vx(i,1,1,1)*tau13(i,1,1,1)
     $                   + vy(i,1,1,1)*tau23(i,1,1,1)
     $                   + vz(i,1,1,1)*tau33(i,1,1,1)
        enddo
      else  
        do i=1,ntot
          tauE1(i,1,1,1) = vx(i,1,1,1)*tau11(i,1,1,1)
     $                   + vy(i,1,1,1)*tau12(i,1,1,1)
          tauE2(i,1,1,1) = vx(i,1,1,1)*tau12(i,1,1,1)
     $                   + vy(i,1,1,1)*tau22(i,1,1,1)
        enddo
      endif

ccc rhoE: div(v.tau)
      call divQ_weak(tauE1,tauE2,tauE3,rkRHS(1,1,1,1,E_EQ),'E1')

      return
      end
c---------------------------------------------------------------------
      subroutine eval_visc_EY(rkRHS)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)

      real qx(lx1,ly1,lz1,lelt)
      real qy(lx1,ly1,lz1,lelt)
      real qz(lx1,ly1,lz1,lelt)

      real Jx(lx1,ly1,lz1,lelt,nspec)
      real Jy(lx1,ly1,lz1,lelt,nspec)
      real Jz(lx1,ly1,lz1,lelt,nspec)

      real dTdx(lx1,ly1,lz1,lelt)
      real dTdy(lx1,ly1,lz1,lelt)
      real dTdz(lx1,ly1,lz1,lelt)

      real dXdx(lx1,ly1,lz1,lelt)
      real dXdy(lx1,ly1,lz1,lelt)
      real dXdz(lx1,ly1,lz1,lelt)

      integer ntot,i,j,e,iface,k,lnf
      
      real dgFF,dgGG,dgHH

      ntot = nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt

ccc heat diffusion
      call liftGrad(dTdx,dTdy,dTdz,temp1,temp1pi,temp1f)
      do i=1,ntot
        qx(i,1,1,1) = cond1(i,1,1,1)*dTdx(i,1,1,1)
        qy(i,1,1,1) = cond1(i,1,1,1)*dTdy(i,1,1,1)
        if(if3d) qz(i,1,1,1) = cond1(i,1,1,1)*dTdz(i,1,1,1)
      enddo

ccc Ji
      do b=1,nspec
        call liftGrad(dXdx,dXdy,dXdz,Xi(1,1,1,1,b),Xipi(1,b),Xif(1,b))
        do i=1,ntot
          Jx(i,1,1,1,b) = diff1(i,1,1,1,b)*mw(b)
     $                         *dXdx(i,1,1,1)/wmean(i,1,1,1)
          Jy(i,1,1,1,b) = diff1(i,1,1,1,b)*mw(b)
     $                         *dXdy(i,1,1,1)/wmean(i,1,1,1)
          if(if3d) then
            Jz(i,1,1,1,b) = diff1(i,1,1,1,b)*mw(b)
     $                           *dXdz(i,1,1,1)/wmean(i,1,1,1)
          endif
        enddo
      enddo

ccc Ji correction for conservation of mass
      do i=1,ntot
        sumJx = 0.
        sumJy = 0.
        if(if3d) sumJz = 0.
        do b=1,nspec
          sumJx = sumJx + Jx(i,1,1,1,b)
          sumJy = sumJy + Jy(i,1,1,1,b)
          if(if3d) sumJz = sumJz + Jz(i,1,1,1,b)
        enddo
        do b=1,nspec
          Jx(i,1,1,1,b) = Jx(i,1,1,1,b) - Yi(i,1,1,1,b)*sumJx
          Jy(i,1,1,1,b) = Jy(i,1,1,1,b) - Yi(i,1,1,1,b)*sumJy
          if(if3d) then
            Jz(i,1,1,1,b) = Jz(i,1,1,1,b) - Yi(i,1,1,1,b)*sumJz  
          endif
        enddo
      enddo

ccc add two diffusion terms q = q + sum(hi*Ji)
      do i=1,ntot
        do b=1,nspec
          qx(i,1,1,1) = qx(i,1,1,1)+hh(i,1,1,1,b)*Jx(i,1,1,1,b)
          qy(i,1,1,1) = qy(i,1,1,1)+hh(i,1,1,1,b)*Jy(i,1,1,1,b)
          if(if3d) then
            qz(i,1,1,1) = qz(i,1,1,1)+hh(i,1,1,1,b)*Jz(i,1,1,1,b)
          endif
        enddo
      enddo

ccc rhoE: div(q)
      call divQ_weak(qx,qy,qz,rkRHS(1,1,1,1,E_EQ),'E1')

ccc rhoYi: div(Ji)
      do b=1,nspec
        call divQ_weak(Jx(1,1,1,1,b),Jy(1,1,1,1,b),Jz(1,1,1,1,b)
     &                ,rkRHS(1,1,1,1,E_EQ+b),'Y1')
      enddo

      return
      end
c---------------------------------------------------------------------