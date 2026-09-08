c---------------------------------------------------------------------      
      subroutine entrovisc(rkRHS,rkU)
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real rkU  (lx1,ly1,lz1,lelt,N_EQ)
      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      real enq1 (lx1,ly1,lz1,lelt)
      real enq2 (lx1,ly1,lz1,lelt)
      real enq3 (lx1,ly1,lz1,lelt)

      integer n,b

      n = nx1*ny1*nz1*nelt

      call eval_enmu

      do b=1,N_EQ
        call liftGrad(enq1,enq2,enq3,rkU(1,1,1,1,b)
     &                ,conspi(1,b),consf(1,b))
        call opcolv(enq1,enq2,enq3,enmu)
c        call divQ_strong(enq1,enq2,enq3,rkRHS(1,1,1,1,b),'AV')
        call divQ_weak(enq1,enq2,enq3,rkRHS(1,1,1,1,b),'AV')
      enddo
      
      return
      end
c---------------------------------------------------------------------
      subroutine liftGrad(enq1,enq2,enq3,var,varpi,varf)
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real  enq1(lx1,ly1,lz1,lelt)
      real  enq2(lx1,ly1,lz1,lelt)
      real  enq3(lx1,ly1,lz1,lelt)

      real var(lx1,ly1,lz1,lelt)
      real varf(lf),varpi(lf)

      real intX(lf)
      real intY(lf)
      real intZ(lf)
      
      integer i,j,e,iface,k,lnf
      
      lnf = nx1*nz1*2*ndim*nelt
      
      call dgGrad_strong(enq1,var,1)
      call dgGrad_strong(enq2,var,2)
      if(if3d) call dgGrad_strong(enq3,var,3)

      ! see equation 50 of Zhang and Peet 2023
      k = 0
      do e=1,nelt
      do iface=1,2*ldim
      do i =1,lx1*lz1
        k = k+1
        intX(k) = unx(i,1,iface,e)*(0.5*varpi(k)-varf(k)) ! U-U_neighbor (weak BR1)
     $          * area(i,1,iface,e)
        intY(k) = uny(i,1,iface,e)*(0.5*varpi(k)-varf(k)) ! U-U_neighbor (weak BR1)
     $          * area(i,1,iface,e)
        if(if3d) then
          intZ(k) = unz(i,1,iface,e)*(0.5*varpi(k)-varf(k)) ! U-U_neighbor (weak BR1)
     $            * area(i,1,iface,e)
        endif
      enddo          
      enddo
      enddo

      do j=1,ndg_facex
        i=dg_face(j)
        enq1(i,1,1,1) = enq1(i,1,1,1) + intX(j)/bm1(i,1,1,1) 
        enq2(i,1,1,1) = enq2(i,1,1,1) + intY(j)/bm1(i,1,1,1)
        if(if3d) enq3(i,1,1,1) = enq3(i,1,1,1) + intZ(j)/bm1(i,1,1,1)
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine divQ_strong(enq1,enq2,enq3,dgu,bctype)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real enq1(lx1,ly1,lz1,lelt)
      real enq2(lx1,ly1,lz1,lelt)
      real enq3(lx1,ly1,lz1,lelt)

      real dgu(lx1,ly1,lz1,lelt)

      real dgu1(lx1,ly1,lz1,lelt)
      real dgu2(lx1,ly1,lz1,lelt)
      real dgu3(lx1,ly1,lz1,lelt)

      real qxpi(lf),qypi(lf),qzpi(lf)
      real qxf(lf),qyf(lf),qzf(lf)

      integer ntot,i,j,e,iface,k,lnf
      real dgFF,dgGG,dgHH

      character bctype*2

      ntot = nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt
      
      call sum_neighbor(enq1,qxpi,qxf)
      call sum_neighbor(enq2,qypi,qyf)
      if(if3d) call sum_neighbor(enq3,qzpi,qzf)

      call compbc_visc(qxpi,bctype)
      call compbc_visc(qypi,bctype)
      if(if3d) call compbc_visc(qzpi,bctype)
      
      k = 0
      do e=1,nelt
      do iface=1,2*ldim
      do i =1,lx1*lz1
        k = k+1
        dgFF = unx(i,1,iface,e)*(0.5*qxpi(k)-qxf(k)) ! U-U_neighbor
        dgGG = uny(i,1,iface,e)*(0.5*qypi(k)-qyf(k)) ! U-U_neighbor
        if(if3d) then 
          dgHH = unz(i,1,iface,e)*(0.5*qzpi(k)-qzf(k)) ! U-U_neighbor
        else
          dgHH = 0.
        endif 
        uf(k) = (dgFF+dgGG+dgHH)*area(i,1,iface,e)
      enddo          
      enddo
      enddo

      do j=1,ndg_facex
        i=dg_face(j)
        dgu(i,1,1,1) = dgu(i,1,1,1) - uf(j)/bm1(i,1,1,1)
      enddo
      
      call dgGrad_strong(dgu1,enq1,1)
      call dgGrad_strong(dgu2,enq2,2)
      if(if3d) call dgGrad_strong(dgu3,enq3,3)
      
      do i=1,ntot
        dgu(i,1,1,1) = dgu(i,1,1,1) - dgu1(i,1,1,1) - dgu2(i,1,1,1)
        if(if3d) dgu(i,1,1,1) = dgu(i,1,1,1) - dgu3(i,1,1,1) 
      enddo 

      return
      end
c---------------------------------------------------------------------
      subroutine divQ_weak(enq1,enq2,enq3,dgu,bctype)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real enq1(lx1,ly1,lz1,lelt)
      real enq2(lx1,ly1,lz1,lelt)
      real enq3(lx1,ly1,lz1,lelt)

      real dgu(lx1,ly1,lz1,lelt)

      real dgu1(lx1,ly1,lz1,lelt)
      real dgu2(lx1,ly1,lz1,lelt)
      real dgu3(lx1,ly1,lz1,lelt)

      real qxpi(lf),qypi(lf),qzpi(lf)
      real qxf(lf),qyf(lf),qzf(lf)

      integer ntot,i,j,e,iface,k,lnf
      real dgFF,dgGG,dgHH

      character bctype*2

      ntot = nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt
      
      call sum_neighbor(enq1,qxpi,qxf)
      call sum_neighbor(enq2,qypi,qyf)
      if(if3d) call sum_neighbor(enq3,qzpi,qzf)

      call compbc_visc(qxpi,bctype)
      call compbc_visc(qypi,bctype)
      if(if3d) call compbc_visc(qzpi,bctype)
      
      ! see equation 53 of Zhang and Peet 2023
      k = 0
      do e=1,nelt
      do iface=1,2*ldim
      do i =1,lx1*lz1
        k = k+1
        dgFF = 0.5*qxpi(k)*unx(i,1,iface,e) ! 0.5*(U+U_neighbor)
        dgGG = 0.5*qypi(k)*uny(i,1,iface,e) ! 0.5*(U+U_neighbor)
        if(if3d) then 
          dgHH = 0.5*qzpi(k)*unz(i,1,iface,e) ! 0.5*(U+U_neighbor)
        else
          dgHH = 0.
        endif
        uf(k) = (dgFF+dgGG+dgHH)*area(i,1,iface,e)
      enddo          
      enddo
      enddo

      do j=1,ndg_facex
        i=dg_face(j)
        dgu(i,1,1,1) = dgu(i,1,1,1) - uf(j)/bm1(i,1,1,1)
      enddo
      
      call dgGrad_weak(dgu1,enq1,1)
      call dgGrad_weak(dgu2,enq2,2)
      if(if3d) call dgGrad_weak(dgu3,enq3,3)
      
      do i=1,ntot
        dgu(i,1,1,1) = dgu(i,1,1,1) - dgu1(i,1,1,1) - dgu2(i,1,1,1)
        if(if3d) dgu(i,1,1,1) = dgu(i,1,1,1) - dgu3(i,1,1,1) 
      enddo 

      return
      end
c---------------------------------------------------------------------
      subroutine eval_enmu

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      integer ntot,nxyz,i,e,ix,iy
      integer tpx1,tpy1,tpx2,tpy2
      
      real dgtmp1(lx1,ly1,lz1,lelt)
      real dgtmp2(lx1,ly1,lz1,lelt)
      real dgtmp3(lx1,ly1,lz1,lelt)
      real dguvw(lx1,ly1,lz1,lelt)
      real enRs,enCs,cmax,enPs,drhoSdt
      real avgSt,maxSt,tmpSt(lx1,ly1,lz1)
      real enlamda(lx1,ly1,lz1,lelt)
      real enmutmp(lx1,ly1,lz1,lelt)
      real enmumax(lelt)
      real entp1,entp2,entp3,entp4,entp5,entp6,entp
      real rhoStp(lx1,ly1,lz1,lelt)
 
      real eps_i,Beta_i,dtt,delu,s_L
      ntot = nx1*ny1*nz1*nelt
      nxyz = nx1*ny1*nz1

      Beta_i = 0.01   ! Decrease -> smoother
      eps_i = 1e-15

      enCs = uparam(2)
      cmax = uparam(3)
      enPs = uparam(4)

      call col3(dgtmp1,rhoS,vx,ntot)
      call col3(dgtmp2,rhoS,vy,ntot)
      if (if3d) call col3(dgtmp3,rhoS,vz,ntot)
      call opdiv(dguvw,dgtmp1,dgtmp2,dgtmp3)

      do e=1,nelt
      do i=1,nxyz
        delu = dudx(i,1,1,e) + dvdy(i,1,1,e)
        if(if3d) delu = delu + dwdz(i,1,1,e)
        s_L = spd(i,1,1,e)/enhe(i,1,1,e)
        thta(i,1,1,e) = (delu*delu)/(delu*delu+Beta_i*(s_L*s_L)+eps_i)
      enddo
      enddo

      do e=1,nelt
        avgSt = vlsc2(rhoS(1,1,1,e),bm1(1,1,1,e),nxyz)
        avgSt = avgSt/vlsum(bm1(1,1,1,e),nxyz)
        do i=1,nxyz
          rhoStp(i,1,1,e) = rhoS(i,1,1,e)-avgSt
        enddo
      enddo
      
      maxSt = glamax(rhoStp,ntot)

      do e=1,nelt
        do i=1,nxyz
          if(istep.eq.1) then
            drhoSdt = (rhoS(i,1,1,e)-rhoS1(i,1,1,e))/dt
          else
            drhoSdt = (1.5*rhoS(i,1,1,e)-2.*rhoS1(i,1,1,e)
     &                        +0.5*rhoS2(i,1,1,e))/dt
          endif
          enRs = drhoSdt + dguvw(i,1,1,e)/bm1(i,1,1,e)
          enmutmp(i,1,1,e) = enCs*(enhe(i,1,1,e)**2.)*abs(enRs)/maxSt
          enlamda(i,1,1,e) = aspd(i,1,1,e)+spd(i,1,1,e)
        enddo
        enmumax(e) = vlmax(enlamda(1,1,1,e),nxyz)
        enmumax(e) = cmax*enhemax(e)*enmumax(e)
      enddo
      
      
      if(if3d) then

        do e=1,nelt
        do ix=1,nx1
        do iy=1,ny1
        do iz=1,nz1
          tpx1=ix+1
          tpx2=ix-1
          tpy1=iy+1
          tpy2=iy-1
          tpz1=iz+1
          tpz2=iz-1
          if (tpx1.gt.nx1) tpx1=ix-1
          if (tpx2.eq.0)   tpx2=ix+1
          if (tpy1.gt.ny1) tpy1=iy-1
          if (tpy2.eq.0)   tpy2=iy+1
          if (tpz1.gt.nz1) tpz1=iz-1
          if (tpz2.eq.0)   tpz2=iz+1
          entp1 = min(enmumax(e),enmutmp(tpx1,iy,iz,e))
          entp2 = min(enmumax(e),enmutmp(tpx2,iy,iz,e))
          entp3 = min(enmumax(e),enmutmp(ix,tpy1,iz,e))
          entp4 = min(enmumax(e),enmutmp(ix,tpy2,iz,e))
          entp5 = min(enmumax(e),enmutmp(ix,iy,tpz1,e))
          entp6 = min(enmumax(e),enmutmp(ix,iy,tpz1,e))
          entp  = min(enmumax(e),enmutmp(ix,iy,iz,e))
          enmu(ix,iy,iz,e) = enPs*(6.*entp+entp1+entp2+entp3+
     $                entp4+entp5+entp6)/12.*thta(ix,iy,iz,e)
        enddo
        enddo
        enddo
        enddo

      else

        do e=1,nelt
        do ix=1,nx1
        do iy=1,ny1
          tpx1=ix+1
          tpx2=ix-1
          tpy1=iy+1
          tpy2=iy-1
          if (tpx1.gt.nx1) tpx1=ix-1
          if (tpx2.eq.0)   tpx2=ix+1
          if (tpy1.gt.ny1) tpy1=iy-1
          if (tpy2.eq.0)   tpy2=iy+1
          entp1 = min(enmumax(e),enmutmp(tpx1,iy,1,e))
          entp2 = min(enmumax(e),enmutmp(tpx2,iy,1,e))
          entp3 = min(enmumax(e),enmutmp(ix,tpy1,1,e))
          entp4 = min(enmumax(e),enmutmp(ix,tpy2,1,e))
          entp  = min(enmumax(e),enmutmp(ix,iy,1,e))
          enmu(ix,iy,1,e) = enPs*(4.*entp+entp1+entp2+entp3+entp4)/8.
     $                    * thta(ix,iy,1,e)
        enddo
        enddo
        enddo

      endif


      return
      end
c---------------------------------------------------------------------
      subroutine eval_enhe

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      integer ntot,nxyz,i,e,ix,iy,iz,iface,ia
      integer tpx1,tpy1,tpx2,tpy2
      real entp1,entp2,entp3,entp4
      real entp5,entp6,entp7,entp8
      real entp9,entp10,entp11,entp12

      ntot = nx1*ny1*nz1*nelt
      nxyz = nx1*ny1*nz1

      if(if3d) then
        
        do e=1,nelt
          entp1 =0.
          entp2 =0.
          entp3 =0.
          entp4 =0.
          entp5 =0.
          entp6 =0.
          entp7 =0.
          entp8 =0.
          entp9 =0.
          entp10=0.
          entp11=0.
          entp12=0.
          entp1=sqrt((xm1(1,1,1,e)-xm1(lx1,1,1,e))**2.
     $              +(ym1(1,1,1,e)-ym1(lx1,1,1,e))**2.
     $              +(zm1(1,1,1,e)-zm1(lx1,1,1,e))**2.)

          entp2=sqrt((xm1(1,1,1,e)-xm1(1,ly1,1,e))**2.
     $              +(ym1(1,1,1,e)-ym1(1,ly1,1,e))**2.
     $              +(zm1(1,1,1,e)-zm1(1,ly1,1,e))**2.)

          entp3=sqrt((xm1(lx1,1,1,e)-xm1(lx1,ly1,1,e))**2.
     $              +(ym1(lx1,1,1,e)-ym1(lx1,ly1,1,e))**2.
     $              +(zm1(lx1,1,1,e)-zm1(lx1,ly1,1,e))**2.)
     
          entp4=sqrt((xm1(1,ly1,1,e)-xm1(lx1,ly1,1,e))**2.
     $              +(ym1(1,ly1,1,e)-ym1(lx1,ly1,1,e))**2.
     $              +(zm1(1,ly1,1,e)-zm1(lx1,ly1,1,e))**2.)

          entp5=sqrt((xm1(1,1,lz1,e)-xm1(lx1,1,lz1,e))**2.
     $              +(ym1(1,1,lz1,e)-ym1(lx1,1,lz1,e))**2.
     $              +(zm1(1,1,lz1,e)-zm1(lx1,1,lz1,e))**2.)

          entp6=sqrt((xm1(1,1,lz1,e)-xm1(1,ly1,lz1,e))**2.
     $              +(ym1(1,1,lz1,e)-ym1(1,ly1,lz1,e))**2.
     $              +(zm1(1,1,lz1,e)-zm1(1,ly1,lz1,e))**2.)

          entp7=sqrt((xm1(lx1,1,lz1,e)-xm1(lx1,ly1,lz1,e))**2.
     $              +(ym1(lx1,1,lz1,e)-ym1(lx1,ly1,lz1,e))**2.
     $              +(zm1(lx1,1,lz1,e)-zm1(lx1,ly1,lz1,e))**2.)

          entp8=sqrt((xm1(1,ly1,lz1,e)-xm1(lx1,ly1,lz1,e))**2.
     $              +(ym1(1,ly1,lz1,e)-ym1(lx1,ly1,lz1,e))**2.
     $              +(zm1(1,ly1,lz1,e)-zm1(lx1,ly1,lz1,e))**2.)

          entp9=sqrt((xm1(lx1,ly1,1,e)-xm1(lx1,ly1,lz1,e))**2.
     $              +(ym1(lx1,ly1,1,e)-ym1(lx1,ly1,lz1,e))**2.
     $              +(zm1(lx1,ly1,1,e)-zm1(lx1,ly1,lz1,e))**2.)

          entp10=sqrt((xm1(1,ly1,1,e)-xm1(1,ly1,lz1,e))**2.
     $               +(ym1(1,ly1,1,e)-ym1(1,ly1,lz1,e))**2.
     $               +(zm1(1,ly1,1,e)-zm1(1,ly1,lz1,e))**2.)

          entp11=sqrt((xm1(1,1,1,e)-xm1(1,1,lz1,e))**2.
     $               +(ym1(1,1,1,e)-ym1(1,1,lz1,e))**2.
     $               +(zm1(1,1,1,e)-zm1(1,1,lz1,e))**2.)

          entp12=sqrt((xm1(lx1,1,1,e)-xm1(lx1,1,lz1,e))**2.
     $               +(ym1(lx1,1,1,e)-ym1(lx1,1,lz1,e))**2.
     $               +(zm1(lx1,1,1,e)-zm1(lx1,1,lz1,e))**2.)
          
          enhemax(e) = min(entp1,entp2,entp3,entp4,
     $                     entp5,entp6,entp7,entp8,
     $                     entp9,entp10,entp11,entp12)
          enhemax(e) = enhemax(e)/real(lx1-2)
          
          do i=1,nxyz
            enhe(i,1,1,e) = enhemax(e)
          enddo
        
        enddo


      else ! 2D

        do e=1,nelt
          entp1=0.
          entp2=0.
          entp3=0.
          entp4=0.
          do i=1,nx1*nz1
            entp1=entp1+area(i,1,1,e)
            entp2=entp2+area(i,1,2,e)
            entp3=entp3+area(i,1,3,e)
            entp4=entp4+area(i,1,4,e)
          enddo
          enhemax(e)=min(entp1,entp2,entp3,entp4)
          enhemax(e)=enhemax(e)/real(lx1-2)
          do i=1,nxyz
            enhe(i,1,1,e)=enhemax(e)
          enddo
        enddo

      endif

      return
      end
c---------------------------------------------------------------------
      subroutine eval_lmf_visc_reg
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real enhef(lf)
      real alpha_i
      integer lnf,i,e,iface,k
      
      lnf = nx1*nz1*2*ldim*nelt
      alpha_i = 0.0          ! 1 or 0
      call full2face(enhef,enhe)

      do e=1,nelt
      do iface=1,2*ldim
      do i=1,lx1*lz1
        k = k+1
        lmfvis(k) = alpha_i*real(lx1-1)**2.0/(2.0*enhef(k))
      enddo          
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------
