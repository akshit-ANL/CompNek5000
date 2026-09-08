c---------------------------------------------------------------------                 
      subroutine eval_lmf

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real asf(lf)
      integer lnf,i,e,iface,k
      real lambda_max,asf_max

      lnf = nx1*nz1*2*ldim*nelt

      call full2face(asf,aspd)

c      do i=1,lnf
c       lmf(i) = sqrt(uxf(i)*uxf(i)+uyf(i)*uyf(i)+uzf(i)*uzf(i))
c       lmf(i) = max(lmf(i),asf(i))
c      enddo
      
      do e=1,nelt
      do iface=1,2*ldim
      do i =1,lx1*lz1
        k = k+1
        lmf(k) = asf(k)+abs(unx(i,1,iface,e)*uxf(k)
     $                     +uny(i,1,iface,e)*uyf(k)
     $                     +unz(i,1,iface,e)*uzf(k) )
      enddo          
      enddo
      enddo
      
      call fgslib_gs_op(dg_hndlx,lmf,1,4,0) ! max

      return
      end 
c---------------------------------------------------------------------
      subroutine sum_neighbor_vars(rkU) ! evaluate F+F_neighbor
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real rkU(lx1,ly1,lz1,lelt,N_EQ)

      real f1(lx1,ly1,lz1,lelt,N_EQ)
      real f2(lx1,ly1,lz1,lelt,N_EQ)
      real f3(lx1,ly1,lz1,lelt,N_EQ)

      integer i,lnf

      lnf = nx1*nz1*2*ndim*nelt

ccc primative variables
      call sum_neighbor(vx,uxpi,uxf)
      call sum_neighbor(vy,uypi,uyf)
      if(if3d) call sum_neighbor(vz,uzpi,uzf)
      call sum_neighbor(vtrans,rhopi,rhof)
      call sum_neighbor(pr,prpi,prf)
      call sum_neighbor(dgeh,dgehpi,dgehf)
      call sum_neighbor(temp1,temp1pi,temp1f)
      do b=1,nspec
        call sum_neighbor(Yi(1,1,1,1,b),Yipi(1,b),Yif(1,b))
        call sum_neighbor(Xi(1,1,1,1,b),Xipi(1,b),Xif(1,b))
      enddo
      
ccc conserved variables
      do i=1,N_EQ
	      call sum_neighbor(rkU(1,1,1,1,i),conspi(1,i),consf(1,i))
      enddo

ccc face values
      if(istep.eq.1) then
        call full2face(xm1f,xm1)
        call full2face(ym1f,ym1)
        if(if3d) call full2face(zm1f,zm1)
      endif

ccc face values for characteristic BCs
      call full2face(ssf,ss)

ccc full LF flux
      if(ifIntLLFfull) then
        call compute_flux(f1,f2,f3)
        do i=1,N_EQ
          call sum_neighbor(f1(1,1,1,1,i),flux1pi(1,i),flux1f(1,i))
          call sum_neighbor(f2(1,1,1,1,i),flux2pi(1,i),flux2f(1,i))
          if(if3d) then
            call sum_neighbor(f3(1,1,1,1,i),flux3pi(1,i),flux3f(1,i))
          endif
        enddo
      endif

      return
      end
c---------------------------------------------------------------------      
      subroutine sum_neighbor(var,varpi,varf) ! evaluate F+F_neighbor for single variable
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real var(lx1,ly1,lz1,lelt)
      real varf(lf),varpi(lf)

      integer lnf

      lnf = nx1*nz1*2*ndim*nelt
	  
      call full2face(varpi,var)
      call copy(varf,varpi,lnf)
      call fgslib_gs_op(dg_hndlx,varpi,1,1,0) ! U+U_neighbor

      return
      end
c---------------------------------------------------------------------
      subroutine compute_flux(f1,f2,f3) ! compute fluxes
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      integer n

      real f1(lx1,ly1,lz1,lelt,N_EQ)
      real f2(lx1,ly1,lz1,lelt,N_EQ)
      real f3(lx1,ly1,lz1,lelt,N_EQ)

      n = nx1*ny1*nz1*nelt

      do i=1,n

ccc rhoU
        f1(i,1,1,1,U_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*vx(i,1,1,1)
     &                   + pr(i,1,1,1)
        f2(i,1,1,1,U_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*vy(i,1,1,1)
        if(if3d) then
          f3(i,1,1,1,U_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*vz(i,1,1,1)
        endif  

ccc rhoV
        f1(i,1,1,1,V_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*vy(i,1,1,1)
        f2(i,1,1,1,V_EQ) = vtrans(i,1,1,1,1)*vy(i,1,1,1)*vy(i,1,1,1)
     &                   + pr(i,1,1,1)
        if(if3d) then
          f3(i,1,1,1,V_EQ) = vtrans(i,1,1,1,1)*vy(i,1,1,1)*vz(i,1,1,1)
        endif

ccc rhoW
        if(if3d) then
          f1(i,1,1,1,W_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*vz(i,1,1,1)
          f2(i,1,1,1,W_EQ) = vtrans(i,1,1,1,1)*vy(i,1,1,1)*vz(i,1,1,1)
          f3(i,1,1,1,W_EQ) = vtrans(i,1,1,1,1)*vz(i,1,1,1)*vz(i,1,1,1) 
     &                     + pr(i,1,1,1)
        endif

ccc rhoE
        f1(i,1,1,1,E_EQ) = vtrans(i,1,1,1,1)*vx(i,1,1,1)*dgeh(i,1,1,1)
        f2(i,1,1,1,E_EQ) = vtrans(i,1,1,1,1)*vy(i,1,1,1)*dgeh(i,1,1,1)
        if(if3d) then
          f3(i,1,1,1,E_EQ)=vtrans(i,1,1,1,1)*vz(i,1,1,1)*dgeh(i,1,1,1)
        endif

ccc rhoYi
        do b=1,nspec
          f1(i,1,1,1,E_EQ+b) = vtrans(i,1,1,1,1)*vx(i,1,1,1)
     &                            *Yi(i,1,1,1,b)
          f2(i,1,1,1,E_EQ+b) = vtrans(i,1,1,1,1)*vy(i,1,1,1)
     &                            *Yi(i,1,1,1,b)
          if(if3d) then
            f3(i,1,1,1,E_EQ+b) = vtrans(i,1,1,1,1)*vz(i,1,1,1)
     &                            *Yi(i,1,1,1,b)
          endif
        enddo    

      enddo


      return 
      end
c---------------------------------------------------------------------