c---------------------------------------------------------------------     
      subroutine dgGrad_strong(dgw1,dgu1,dirc)

      include 'SIZE'
      include 'TOTAL'

      real dgw1(lx1,ly1,lz1,lelt)
      real dgu1(lx1,ly1,lz1,lelt)

      integer dirc

c ccc using standard Nek5000 strong form derivative
c       if(dirc.eq.1) call dudxyz(dgw1,dgu1,rxm1,sxm1,txm1,jacm1,2,dirc)
c       if(dirc.eq.2) call dudxyz(dgw1,dgu1,rym1,sym1,tym1,jacm1,2,dirc)
c       if(dirc.eq.3) call dudxyz(dgw1,dgu1,rzm1,szm1,tzm1,jacm1,2,dirc)

ccc updated routines - produces same results as code above
      if(dirc.eq.1) call dudxyz_strong(dgw1,dgu1,rxm1,sxm1,txm1,dirc)
      if(dirc.eq.2) call dudxyz_strong(dgw1,dgu1,rym1,sym1,tym1,dirc)
      if(dirc.eq.3) call dudxyz_strong(dgw1,dgu1,rzm1,szm1,tzm1,dirc)

      return
      end
c---------------------------------------------------------------------   
      subroutine dudxyz_strong(dgw1,dgu1,dgrm,dgsm,dgtm,dirc)

ccc this subroutine produces the same result as dudxyz in Nek5000/core

      include 'SIZE'
      include 'TOTAL'

      real dgw1(lx1,ly1,lz1,lelt)
      real dgu1(lx1,ly1,lz1,lelt)

      real dgrm(lx1,ly1,lz1,lelt)
      real dgsm(lx1,ly1,lz1,lelt)
      real dgtm(lx1,ly1,lz1,lelt)

      integer dirc,ntot
      integer e,i,j,k,l
      integer nxyz

      ntot = nx1*ny1*nz1*nelt
      nxyz = nx1*ny1*nz1

      do e=1,nelt

         do k=1,lz1
         do j=1,ly1 
         do i=1,lx1 ! loop through all GLL points
          
           dgw1(i,j,k,e) = 0.
           
           do l=1,lx1 ! loop through the r,s,t direction
             
             dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/dr*dr/dx
     &         + dxm1(i,l)*dgu1(l,j,k,e)*dgrm(i,j,k,e)  
     &         * wxm1(i)*wym1(j)*wzm1(k)

             dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dx
     &         + dym1(j,l)*dgu1(i,l,k,e)*dgsm(i,j,k,e)
     &         * wxm1(i)*wym1(j)*wzm1(k)

             if(if3d) then
               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/dt*dt/dx
     &           + dzm1(k,l)*dgu1(i,j,l,e)*dgtm(i,j,k,e)
     &           * wxm1(i)*wym1(j)*wzm1(k)
             endif

           enddo ! close l loop
         
         enddo ! close i loop
         enddo ! close j loop
         enddo ! close k loop
         
         call invcol2(dgw1(1,1,1,e),w3m1,nxyz) ! divide by w_i*w_j*w_k
      
      enddo ! close e loop
      
      call col2(dgw1,jacmi,ntot) ! divide by Jacobian

      return
      end
c---------------------------------------------------------------------      
      subroutine dgGrad_weak(dgw1,dgu1,dirc)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS' ! remove this later

      real dgw1(lx1,ly1,lz1,lelt)
      real dgu1(lx1,ly1,lz1,lelt)

      real dgrm(lx1,ly1,lz1,lelt)
      real dgsm(lx1,ly1,lz1,lelt)
      real dgtm(lx1,ly1,lz1,lelt)

      real varf(lf),varpi(lf)

      real intX(lf)
      real intY(lf)
      real intZ(lf)

      integer dirc,ntot
      integer e,i,j,k,l
      integer nxyz

      ntot = nx1*ny1*nz1*nelt
      nxyz = nx1*ny1*nz1

      do e=1,nelt

         do k=1,lz1
         do j=1,ly1 
         do i=1,lx1 ! loop through all GLL points
          
           dgw1(i,j,k,e) = 0.
           
           do l=1,lx1 ! loop through the r direction
             
             if(dirc.eq.1) then

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/dr*dr/dx
     &           - dxm1(l,i)*dgu1(l,j,k,e)*rxm1(i,j,k,e)  
     &           * wxm1(l)*wym1(j)*wzm1(k)

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dx
     &           - dym1(l,j)*dgu1(i,l,k,e)*sxm1(i,j,k,e)
     &           * wxm1(l)*wym1(j)*wzm1(k)

               if(if3d) then
                 dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dx
     &             - dzm1(l,k)*dgu1(i,j,l,e)*txm1(i,j,k,e)
     &             * wxm1(l)*wym1(j)*wzm1(k)
               endif

             else if(dirc.eq.2) then

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/dr*dr/dy
     &           - dxm1(l,i)*dgu1(l,j,k,e)*rym1(i,j,k,e)  
     &           * wxm1(i)*wym1(l)*wzm1(k)

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dy
     &           - dym1(l,j)*dgu1(i,l,k,e)*sym1(i,j,k,e)
     &           * wxm1(i)*wym1(l)*wzm1(k)

               if(if3d) then
                 dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dy
     &             - dzm1(l,k)*dgu1(i,j,l,e)*tym1(i,j,k,e)
     &             * wxm1(i)*wym1(l)*wzm1(k)
               endif

             else if(dirc.eq.3) then

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/dr*dr/dz
     &           - dxm1(l,i)*dgu1(l,j,k,e)*rzm1(i,j,k,e)  
     &           * wxm1(i)*wym1(j)*wzm1(l)

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dz
     &           - dym1(l,j)*dgu1(i,l,k,e)*szm1(i,j,k,e)
     &           * wxm1(i)*wym1(j)*wzm1(l)

               dgw1(i,j,k,e) = dgw1(i,j,k,e) ! dF/ds*ds/dz
     &           - dzm1(l,k)*dgu1(i,j,l,e)*tzm1(i,j,k,e)
     &           * wxm1(i)*wym1(j)*wzm1(l)

             endif 

           enddo ! close l loop
         
         enddo ! close i loop
         enddo ! close j loop
         enddo ! close k loop
         
         call invcol2(dgw1(1,1,1,e),w3m1,nxyz) ! divide by w_i*w_j*w_k
      
      enddo ! close e loop
      
      call col2(dgw1,jacmi,ntot) ! divide by Jacobian

      return
      end
c---------------------------------------------------------------------