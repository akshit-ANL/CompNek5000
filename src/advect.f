c---------------------------------------------------------------------      
      subroutine advect(rkRHS)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)

      call advect_volInt(rkRHS)
      call advect_interface(rkRHS)

      return
      end
c---------------------------------------------------------------------  
      subroutine advect_interface(rkRHS)

ccc Interface contribution to divergence of convective fluxes
ccc Includes a vol terms which is the result of the SBP reverse integration by parts
ccc See equation 28 of Zhang 2024

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      real dgFF,dgGG,dgHH
      real llf_d,llf_nd_x,llf_nd_y,llf_nd_z
      real volt_x,volt_y,volt_z

      integer i,j,e,iface,k,lnf

      lnf = nx1*nz1*2*ndim*nelt

      call rzero(uf,lnf)

      call eval_lmf

ccc rhoU
      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k = k+1
         volt_x = (-rhof(k)*uxf(k)*uxf(k)-prf(k))*unx(i,1,iface,e)
         volt_y = (-rhof(k)*uxf(k)*uyf(k))*uny(i,1,iface,e)
         if(if3d) then
           volt_z = (-rhof(k)*uxf(k)*uzf(k))*unz(i,1,iface,e)
         else
           volt_z = 0.
         endif
         if(ifIntLLFfull) then
            llf_nd_x = 0.5*flux1pi(k,U_EQ)*unx(i,1,iface,e)
            llf_nd_y = 0.5*flux2pi(k,U_EQ)*uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = 0.5*flux3pi(k,U_EQ)*unz(i,1,iface,e)
            else
              llf_nd_z = 0.
            endif  
         else
            llf_nd_x = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*uxpi(k)
     &                 +0.5*prpi(k))*unx(i,1,iface,e)
            llf_nd_y = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*uypi(k))
     &                 *uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*uzpi(k))
     &                 *unz(i,1,iface,e)
            else 
              llf_nd_z = 0.
            endif
         endif
         llf_d = -0.5*lmf(k)*(conspi(k,U_EQ)-consf(k,U_EQ)*2.)
         uf(k)=(volt_x+volt_y+volt_z+llf_d+llf_nd_x+llf_nd_y+llf_nd_z)
     &          *area(i,1,iface,e)
        enddo
       enddo
      enddo
      do j=1,ndg_facex
        i=dg_face(j)
        rkRHS(i,1,1,1,U_EQ) = rkRHS(i,1,1,1,U_EQ)+uf(j)/bm1(i,1,1,1)
      enddo

ccc rhoV
      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k = k+1
         volt_x = (-rhof(k)*uxf(k)*uyf(k))*unx(i,1,iface,e)
         volt_y = (-rhof(k)*uyf(k)*uyf(k)-prf(k))*uny(i,1,iface,e)
         if(if3d) then
           volt_z = (-rhof(k)*uyf(k)*uzf(k))*unz(i,1,iface,e)
         else
           volt_z = 0.
         endif
         if(ifIntLLFfull) then
            llf_nd_x = 0.5*flux1pi(k,V_EQ)*unx(i,1,iface,e)
            llf_nd_y = 0.5*flux2pi(k,V_EQ)*uny(i,1,iface,e)
            if(if3d) then 
              llf_nd_z = 0.5*flux3pi(k,V_EQ)*unz(i,1,iface,e)
            else
              llf_nd_z = 0.
            endif  
         else 
            llf_nd_x = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*uypi(k))
     &                 *unx(i,1,iface,e)
            llf_nd_y = (0.5*0.5*0.5*rhopi(k)*uypi(k)*uypi(k)
     &                 +0.5*prpi(k))*uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = (0.5*0.5*0.5*rhopi(k)*uypi(k)*uzpi(k))
     &                 *unz(i,1,iface,e)
            else 
              llf_nd_z = 0.
            endif
         endif
         llf_d = -0.5*lmf(k)*(conspi(k,V_EQ)-consf(k,V_EQ)*2.)
         uf(k)=(volt_x+volt_y+volt_z+llf_d+llf_nd_x+llf_nd_y+llf_nd_z)
     &          *area(i,1,iface,e)
        enddo
       enddo
      enddo
      do j=1,ndg_facex
        i=dg_face(j)
        rkRHS(i,1,1,1,V_EQ) = rkRHS(i,1,1,1,V_EQ)+uf(j)/bm1(i,1,1,1)
      enddo

ccc rhoW
      if(if3d) then
      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k = k+1
         volt_x = (-rhof(k)*uxf(k)*uzf(k))*unx(i,1,iface,e)
         volt_y = (-rhof(k)*uyf(k)*uzf(k))*uny(i,1,iface,e)
         volt_z = (-rhof(k)*uzf(k)*uzf(k)-prf(k))*unz(i,1,iface,e)
         if(ifIntLLFfull) then
            llf_nd_x = 0.5*flux1pi(k,W_EQ)*unx(i,1,iface,e)
            llf_nd_y = 0.5*flux2pi(k,W_EQ)*uny(i,1,iface,e)
            llf_nd_z = 0.5*flux3pi(k,W_EQ)*unz(i,1,iface,e)
         else 
            llf_nd_x = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*uzpi(k))
     &                 *unx(i,1,iface,e)
            llf_nd_y = (0.5*0.5*0.5*rhopi(k)*uypi(k)*uzpi(k))
     &                 *uny(i,1,iface,e)
            llf_nd_z = (0.5*0.5*0.5*rhopi(k)*uzpi(k)*uzpi(k)
     &                 +0.5*prpi(k))*unz(i,1,iface,e)
         endif
         llf_d = -0.5*lmf(k)*(conspi(k,W_EQ)-consf(k,W_EQ)*2.)
         uf(k)=(volt_x+volt_y+volt_z+llf_d+llf_nd_x+llf_nd_y+llf_nd_z)
     &          *area(i,1,iface,e)
        enddo
       enddo
      enddo
      do j=1,ndg_facex
        i=dg_face(j)
        rkRHS(i,1,1,1,W_EQ) = rkRHS(i,1,1,1,W_EQ)+uf(j)/bm1(i,1,1,1)
      enddo
      endif

ccc rhoE
      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k = k+1
         volt_x = (-rhof(k)*uxf(k)*dgehf(k))*unx(i,1,iface,e)
         volt_y = (-rhof(k)*uyf(k)*dgehf(k))*uny(i,1,iface,e)
         if(if3d) then
           volt_z = (-rhof(k)*uzf(k)*dgehf(k))*unz(i,1,iface,e)
         else
           volt_z = 0.
         endif
         if(ifIntLLFfull) then
            llf_nd_x = 0.5*flux1pi(k,E_EQ)*unx(i,1,iface,e)
            llf_nd_y = 0.5*flux2pi(k,E_EQ)*uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = 0.5*flux3pi(k,E_EQ)*unz(i,1,iface,e)
            else
              llf_nd_z = 0.
            endif  
         else 
            llf_nd_x = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*dgehpi(k))
     &                 *unx(i,1,iface,e)
            llf_nd_y = (0.5*0.5*0.5*rhopi(k)*uypi(k)*dgehpi(k))
     &                 *uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = (0.5*0.5*0.5*rhopi(k)*uzpi(k)*dgehpi(k))
     &                 *unz(i,1,iface,e)
            else 
              llf_nd_z = 0.
            endif
         endif
         llf_d = -0.5*lmf(k)*(conspi(k,E_EQ)-consf(k,E_EQ)*2.)
         uf(k)=(volt_x+volt_y+volt_z+llf_d+llf_nd_x+llf_nd_y+llf_nd_z)
     &          *area(i,1,iface,e)
        enddo
       enddo
      enddo
      do j=1,ndg_facex
        i=dg_face(j)
        rkRHS(i,1,1,1,E_EQ) = rkRHS(i,1,1,1,E_EQ)+uf(j)/bm1(i,1,1,1)
      enddo

ccc rhoYi
      do b=1,nspec
      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k = k+1
         volt_x = (-rhof(k)*Yif(k,b)*uxf(k))*unx(i,1,iface,e)
         volt_y = (-rhof(k)*Yif(k,b)*uyf(k))*uny(i,1,iface,e)
         if(if3d) then
           volt_z = (-rhof(k)*Yif(k,b)*uzf(k))*unz(i,1,iface,e)
         else
           volt_z = 0.
         endif
         if(ifIntLLFfull) then
            llf_nd_x = 0.5*flux1pi(k,E_EQ+b)*unx(i,1,iface,e)
            llf_nd_y = 0.5*flux2pi(k,E_EQ+b)*uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = 0.5*flux3pi(k,E_EQ+b)*unz(i,1,iface,e)
            else
              llf_nd_z = 0.
            endif  
         else
            llf_nd_x = (0.5*0.5*0.5*rhopi(k)*uxpi(k)*Yipi(k,b))
     &                 *unx(i,1,iface,e)
            llf_nd_y = (0.5*0.5*0.5*rhopi(k)*uypi(k)*Yipi(k,b))
     &                 *uny(i,1,iface,e)
            if(if3d) then
              llf_nd_z = (0.5*0.5*0.5*rhopi(k)*uzpi(k)*Yipi(k,b))
     &                 *unz(i,1,iface,e)
            else 
              llf_nd_z = 0.
            endif
         endif
         llf_d = -0.5*lmf(k)*(conspi(k,E_EQ+b)-consf(k,E_EQ+b)*2.)
         uf(k)=(volt_x+volt_y+volt_z+llf_d+llf_nd_x+llf_nd_y+llf_nd_z)
     &          *area(i,1,iface,e)
        enddo
       enddo
      enddo
      do j=1,ndg_facex
        i=dg_face(j)
        rkRHS(i,1,1,1,E_EQ+b) = rkRHS(i,1,1,1,E_EQ+b)+uf(j)/bm1(i,1,1,1)
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine advect_volInt(rkRHS)

ccc Weak volume derivatives with two-point flux
ccc this routine might need to modify rxm1 to include average of 2 points i-k
ccc See equation 25 of Zhang 2024

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      
      integer e,i,j,k,l,b

      real rho1,rho2,vx1,vx2,vy1,vy2,vz1,vz2
      real pr1,pr2,h1,h2,iT1,iT2
      real ke1,ke2
      real y1(nspec),y2(nspec)
      real fx(N_EQ),fy(N_EQ),fz(N_EQ)

      do e=1,nelt

ccc r direction
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1 ! loop through all GLL points
           do b=1,N_EQ
             rkRHS(i,j,k,e,b) = 0.
           enddo
           rho1 = vtrans(i,j,k,e,1) ! set (i,j,k) point 1 values
           vx1  = vx(i,j,k,e)
           vy1  = vy(i,j,k,e)
           vz1  = vz(i,j,k,e)
           pr1  = pr(i,j,k,e)
           h1   = dgeh(i,j,k,e) 
           iT1  = 1./temp1(i,j,k,e)
           ke1  = vx1*vx1+vy1*vy1+vz1*vz1
           do b=1,nspec
             y1(b) = Yi(i,j,k,e,b)
           enddo
           do l=1,lx1 ! loop through the r direction
             rho2 = vtrans(l,j,k,e,1) ! set (l,j,k) point 2 values
             vx2  = vx(l,j,k,e)   
             vy2  = vy(l,j,k,e)
             vz2  = vz(l,j,k,e)  
             pr2  = pr(l,j,k,e)
             h2   = dgeh(l,j,k,e)
             iT2  = 1./temp1(l,j,k,e)
             ke2  = vx2*vx2+vy2*vy2+vz2*vz2
             do b=1,nspec
               y2(b) = Yi(l,j,k,e,b)
             enddo
             if(ifVolEC) then
               call getFlux_EC(fx,fy,fz
     &                        ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                        ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             else 
               call getFlux_KEP(fx,fy,fz
     &                         ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                         ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             endif
             do b=1,N_EQ           
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFx/dr*dr/dx
     &          +2.*dxm1(i,l)*fx(b)*rxm1(i,j,k,e)/jacm1(i,j,k,e)    
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFy/dr*dr/dy
     &          +2.*dxm1(i,l)*fy(b)*rym1(i,j,k,e)/jacm1(i,j,k,e)
               if(if3d) then
                 rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFz/dr*dr/dz
     &            +2.*dxm1(i,l)*fz(b)*rzm1(i,j,k,e)/jacm1(i,j,k,e)
               endif   
             enddo   
           enddo ! close l loop
         enddo ! close i loop
         enddo ! close j loop
         enddo ! close k loop

ccc s direction
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1 ! loop through all GLL points
           rho1 = vtrans(i,j,k,e,1) ! set (i,j,k) point 1 values
           vx1  = vx(i,j,k,e)
           vy1  = vy(i,j,k,e)
           vz1  = vz(i,j,k,e)
           pr1  = pr(i,j,k,e)
           h1   = dgeh(i,j,k,e) 
           iT1  = 1./temp1(i,j,k,e)
           ke1  = vx1*vx1+vy1*vy1+vz1*vz1
           do b=1,nspec
             y1(b) = Yi(i,j,k,e,b)
           enddo  
           do l=1,lx1 ! loop through the s direction
             rho2 = vtrans(i,l,k,e,1) ! set (i,l,k) point 2 values 
             vx2  = vx(i,l,k,e)   
             vy2  = vy(i,l,k,e) 
             vz2  = vz(i,l,k,e)  
             pr2  = pr(i,l,k,e)
             h2   = dgeh(i,l,k,e)
             iT2  = 1./temp1(i,l,k,e)
             ke2  = vx2*vx2+vy2*vy2+vz2*vz2
             do b=1,nspec
               y2(b) = Yi(i,l,k,e,b)
             enddo
             if(ifVolEC) then
               call getFlux_EC(fx,fy,fz
     &                        ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                        ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             else 
               call getFlux_KEP(fx,fy,fz
     &                         ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                         ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             endif           
             do b=1,N_EQ
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFx/ds*ds/dx
     &          +2.*fx(b)*dym1(j,l)*sxm1(i,j,k,e)/jacm1(i,j,k,e)
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFy/ds*ds/dy
     &          +2.*fy(b)*dym1(j,l)*sym1(i,j,k,e)/jacm1(i,j,k,e)
               if(if3d) then
                 rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFz/ds*ds/dz
     &            +2.*fz(b)*dym1(j,l)*szm1(i,j,k,e)/jacm1(i,j,k,e)
               endif
             enddo      
           enddo ! close l loop
         enddo ! close i loop
         enddo ! close j loop
         enddo ! close k loop

ccc t direction
         if(if3d) then
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1 ! loop through all GLL points
           rho1 = vtrans(i,j,k,e,1) ! set (i,j,k) point 1 values
           vx1  = vx(i,j,k,e)
           vy1  = vy(i,j,k,e)
           vz1  = vz(i,j,k,e)
           pr1  = pr(i,j,k,e)
           h1   = dgeh(i,j,k,e) 
           iT1  = 1./temp1(i,j,k,e)
           ke1  = vx1*vx1+vy1*vy1+vz1*vz1
           do b=1,nspec
             y1(b) = Yi(i,j,k,e,b)
           enddo  
           do l=1,lx1 ! loop through the t direction
             rho2 = vtrans(i,j,l,e,1) ! set (i,j,l) point 2 values 
             vx2  = vx(i,j,l,e)   
             vy2  = vy(i,j,l,e) 
             vz2  = vz(i,j,l,e)  
             pr2  = pr(i,j,l,e)
             h2   = dgeh(i,j,l,e)
             iT2  = 1./temp1(i,j,l,e)
             ke2  = vx2*vx2+vy2*vy2+vz2*vz2
             do b=1,nspec
               y2(b) = Yi(i,j,l,e,b)
             enddo
             if(ifVolEC) then
               call getFlux_EC(fx,fy,fz
     &                        ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                        ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             else 
               call getFlux_KEP(fx,fy,fz
     &                         ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                         ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
             endif           
             do b=1,N_EQ
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFx/dt*dt/dx
     &          +2.*fx(b)*dzm1(k,l)*txm1(i,j,k,e)/jacm1(i,j,k,e)
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFy/dt*dt/dy
     &          +2.*fy(b)*dzm1(k,l)*tym1(i,j,k,e)/jacm1(i,j,k,e)
               rkRHS(i,j,k,e,b) = rkRHS(i,j,k,e,b) ! dFz/dt*dt/dz
     &          +2.*fz(b)*dzm1(k,l)*tzm1(i,j,k,e)/jacm1(i,j,k,e)
             enddo
           enddo ! close l loop
         enddo ! close i loop
         enddo ! close j loop
         enddo ! close k loop
         endif ! close if3d conditional

      enddo ! close e element loop
      
      return
      end
c---------------------------------------------------------------------
      subroutine getFlux_KEP(fx,fy,fz
     &                      ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                      ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
      
ccc Pirozzoli kinetic energy preserving (KEP) flux

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rho1,rho2,vx1,vx2,vy1,vy2,vz1,vz2
      real pr1,pr2,h1,h2,iT1,iT2
      real ke1,ke2
      real y1(nspec),y2(nspec)
      real fx(N_EQ),fy(N_EQ),fz(N_EQ)
      real ecAvg,ecLog
      integer b

      ! rho*Y_i
      do b=1,nspec
        fx(E_EQ+b) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(y1(b),y2(b))
        fy(E_EQ+b) = ecAvg(rho1,rho2)*ecAvg(vy1,vy2)*ecAvg(y1(b),y2(b))
        if(if3d) then
        fz(E_EQ+b) = ecAvg(rho1,rho2)*ecAvg(vz1,vz2)*ecAvg(y1(b),y2(b))
        endif
      enddo

      ! rho*u
      fx(U_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(vx1,vx2)
     &         + ecAvg(pr1,pr2)
      fy(U_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(vy1,vy2)
      if(if3d) then
        fz(U_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(vz1,vz2)
      endif

      ! rho*v   
      fx(V_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(vy1,vy2)
      fy(V_EQ) = ecAvg(rho1,rho2)*ecAvg(vy1,vy2)*ecAvg(vy1,vy2)
     &         + ecAvg(pr1,pr2)
      if(if3d) then
        fz(V_EQ) = ecAvg(rho1,rho2)*ecAvg(vy1,vy2)*ecAvg(vz1,vz2)  
      endif

      ! rho*w - need to fill in for 3D
      if(if3d) then
        fx(W_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(vz1,vz2)
        fy(W_EQ) = ecAvg(rho1,rho2)*ecAvg(vy1,vy2)*ecAvg(vz1,vz2)
        fz(W_EQ) = ecAvg(rho1,rho2)*ecAvg(vz1,vz2)*ecAvg(vz1,vz2)
     &           + ecAvg(pr1,pr2)
      endif

      ! rho*e  
      fx(E_EQ) = ecAvg(rho1,rho2)*ecAvg(vx1,vx2)*ecAvg(h1,h2)
      fy(E_EQ) = ecAvg(rho1,rho2)*ecAvg(vy1,vy2)*ecAvg(h1,h2)
      if(if3d) then
        fz(E_EQ) = ecAvg(rho1,rho2)*ecAvg(vz1,vz2)*ecAvg(h1,h2)
      endif

	    return 
	    end
c---------------------------------------------------------------------
      subroutine getFlux_EC(fx,fy,fz
     &                     ,rho1,vx1,vy1,vz1,pr1,h1,iT1,ke1,y1
     &                     ,rho2,vx2,vy2,vz2,pr2,h2,iT2,ke2,y2 )
      
ccc Gouasmi et al. (2020) multi-component entropy-conservative flux

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rho1,rho2,vx1,vx2,vy1,vy2,vz1,vz2
      real pr1,pr2,h1,h2,iT1,iT2
      real ke1,ke2
      real y1(nspec),y2(nspec)
      real fx(N_EQ),fy(N_EQ),fz(N_EQ)
      real ecAvg,ecLog
      integer b

      real fjT(npa),Tavg,Tavg2,Tprod
      real sum1,qq
      integer c,g

      ! rho*Y_i
      do b=1,nspec
        fx(E_EQ+b) = ecLog(rho1*y1(b),rho2*y2(b))*ecAvg(vx1,vx2)
        fy(E_EQ+b) = ecLog(rho1*y1(b),rho2*y2(b))*ecAvg(vy1,vy2)
        if(if3d) then
          fz(E_EQ+b) = ecLog(rho1*y1(b),rho2*y2(b))*ecAvg(vz1,vz2)
        endif
      enddo

      ! rho*U
      fx(U_EQ) = 0.
      fy(U_EQ) = 0.
      if(if3d) fz(U_EQ) = 0.
      do b=1,nspec
        fx(U_EQ) = fx(U_EQ) + ecAvg(vx1,vx2)*fx(E_EQ+b)  
     &           + (1./ecAvg(iT1,iT2))*(R_u/mw(b)) ! p_hat
     &           * ecAvg(rho1*y1(b),rho2*y2(b))   
        fy(U_EQ) = fy(U_EQ) + ecAvg(vx1,vx2)*fy(E_EQ+b)
        if(if3d) fz(U_EQ) = fz(U_EQ) + ecAvg(vx1,vx2)*fz(E_EQ+b)
      enddo

      ! rho*V   
      fx(V_EQ) = 0.
      fy(V_EQ) = 0.
      if(if3d) fz(V_EQ) = 0.
      do b=1,nspec
        fx(V_EQ) = fx(V_EQ) + ecAvg(vy1,vy2)*fx(E_EQ+b)   
        fy(V_EQ) = fy(V_EQ) + ecAvg(vy1,vy2)*fy(E_EQ+b)  
     &           + (1./ecAvg(iT1,iT2))*(R_u/mw(b)) ! p_hat
     &           * ecAvg(rho1*y1(b),rho2*y2(b)) 
        if(if3d) fz(V_EQ) = fz(V_EQ) + ecAvg(vy1,vy2)*fz(E_EQ+b) 
      enddo

      ! rho*W
      if(if3d) then
      fx(W_EQ) = 0.
      fy(W_EQ) = 0.
      fz(W_EQ) = 0.
      do b=1,nspec
        fx(W_EQ) = fx(W_EQ) + ecAvg(vz1,vz2)*fx(E_EQ+b)   
        fy(W_EQ) = fy(W_EQ) + ecAvg(vz1,vz2)*fy(E_EQ+b)
        fz(W_EQ) = fz(W_EQ) + ecAvg(vz1,vz2)*fz(E_EQ+b) 
     &           + (1./ecAvg(iT1,iT2))*(R_u/mw(b)) ! p_hat
     &           * ecAvg(rho1*y1(b),rho2*y2(b)) 
      enddo
      endif

      ! rho*E 
      Tavg = ecAvg(1./iT1,1./iT2)
      Tavg2 = ecAvg(1./iT1**2.,1./iT2**2.)
      Tprod = (1./iT1)*(1./iT2)
      fjT(1) = 1.
      fjT(2) = Tavg
      fjT(3) = (2./3.)*Tavg*Tavg+(1./3.)*Tavg2
      fjT(4) = Tavg*Tavg2
      fx(E_EQ) = 0.
      fy(E_EQ) = 0.
      if(if3d) fz(E_EQ) = 0.
      do b=1,nspec
        g = (b-1)*(npa+1)
        sum1 = 0.
        do c=2,npa
          sum1 = sum1+bi(g+c+1)*fjT(c-1)*Tprod
        enddo
        qq = bi(g+1) + bi(g+2)/ecLog(iT1,iT2)
     &     + sum1 - 0.5*ecAvg(ke1,ke2)
        fx(E_EQ) = fx(E_EQ) + qq * fx(E_EQ+b)
        fy(E_EQ) = fy(E_EQ) + qq * fy(E_EQ+b)
        if(if3d) fz(E_EQ) = fz(E_EQ) + qq * fz(E_EQ+b)
      enddo         
      if(if3d) then
        fx(E_EQ) = fx(E_EQ) + ecAvg(vx1,vx2)*fx(U_EQ)
     &                      + ecAvg(vy1,vy2)*fx(V_EQ) 
     &                      + ecAvg(vz1,vz2)*fx(W_EQ) 
        fy(E_EQ) = fy(E_EQ) + ecAvg(vx1,vx2)*fy(U_EQ)
     &                      + ecAvg(vy1,vy2)*fy(V_EQ)
     &                      + ecAvg(vz1,vz2)*fy(W_EQ)
        fz(E_EQ) = fz(E_EQ) + ecAvg(vx1,vx2)*fz(U_EQ)
     &                      + ecAvg(vy1,vy2)*fz(V_EQ)
     &                      + ecAvg(vz1,vz2)*fz(W_EQ)
      else
        fx(E_EQ) = fx(E_EQ) + ecAvg(vx1,vx2)*fx(U_EQ)
     &                      + ecAvg(vy1,vy2)*fx(V_EQ) 
        fy(E_EQ) = fy(E_EQ) + ecAvg(vx1,vx2)*fy(U_EQ)
     &                      + ecAvg(vy1,vy2)*fy(V_EQ)
      endif

	    return 
	    end
c---------------------------------------------------------------------
      function ecAvg(al,ar)
	    real al,ar 
	    ecAvg = 0.5*(al+ar)
	    return 
	    end
c---------------------------------------------------------------------
      function ecLog(al,ar)
      ! log average from Ismail and Roe (2009)
      ! expansion to remove singularity at al=ar
      real al,ar
      real xi,f2,f,eps2
      eps2 = 1.e-4
      f2 = (al*(al-2.*ar)+ar*ar)/(al*(al+2.*ar)+ar*ar) ! f^2
      !if(al.eq.0..and.ar.eq.0.) then
      if(al.le.0..or.ar.le.0.) then
        ecLog = 0.
      else if(f2.ge.eps2) then
        ecLog = (ar-al)/log(ar/al)
      else
        xi = ar/al
        f = (xi-1.)/(xi+1.)
        ecLog = (ar+al)
     &      / (2.+(2./3.)*f**2.+(2./5.)*f**4.+(2./7.)*f**6.)
      endif
      return 
      end
c---------------------------------------------------------------------
