c---------------------------------------------------------------------      
      subroutine compbc

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      call compbc_setup
      call compbc_flux

      return
      end
c---------------------------------------------------------------------      
      subroutine compbc_setup

c     This routine computes ghost node values bf which share
c     a boundary with interior faces. The pi variables fed in
c     currently represent the interior values since (U+U_neighbor)=U_interior.

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      integer i,e,iface,k,lnf
      integer ntot,nxz,id_face
      integer j,c,b
      character cb*3
      real x,y,z
      real ux,uy,uz
      real yy(nspec),xx(nspec)
      real wmeanbf
      real cpbari

      nxz = nx1*nz1
      ntot = nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt

      if(istep.eq.1) then
        call rzero(rhobf,lnf)
        call rzero(temp1bf,lnf)
        call rzero(uxbf,lnf)
        call rzero(uybf,lnf)
        call rzero(uzbf,lnf)
        call rzero(prbf,lnf)
        call rzero(dgehbf,lnf)
        do b=1,nspec
          call rzero(Yibf(1,b),lnf)
        enddo
      endif

      do e=1,nelt
      do iface=1,2*ldim
        
        cb = cbc(iface,e,1)
        id_face = bc(5,iface,e,1)
        k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz

ccc adiabatic no-slip wall
        if(cb.eq.'W  ') then
           do i=1,nxz
             uxbf(k) = -uxpi(k)
             uybf(k) = -uypi(k) 
             uzbf(k) = -uzpi(k)
             
             ! testing delete for lid-driven cavity only
             !if(uny(i,1,iface,e).gt.0.75) then ! lid-driven cavity delete
             !  uxbf(k) = 2.*27.83-uxpi(k)
             !endif

             prbf(k) = prpi(k)
             temp1bf(k) = temp1pi(k)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               wmeanbf = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             rhobf(k) = rhopi(k)
             dgehbf(k) = dgehpi(k)
             k = k+1
           enddo
        endif

ccc slip wall
        if(cb.eq.'SYM') then
           do i=1,nxz         
             unInt   = uxpi(k)*unx(i,1,iface,e)
     $               + uypi(k)*uny(i,1,iface,e)
     $               + uzpi(k)*unz(i,1,iface,e)     
             uxbf(k) = uxpi(k)-2.*unInt*unx(i,1,iface,e)
             uybf(k) = uypi(k)-2.*unInt*uny(i,1,iface,e)
             uzbf(k) = uzpi(k)-2.*unInt*unz(i,1,iface,e)
             prbf(k) = prpi(k)
             temp1bf(k) = temp1pi(k)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               wmeanbf = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             rhobf(k) = rhopi(k)
             dgehbf(k) = dgehpi(k)
             k = k+1
            enddo
        endif

ccc subsonic inflow (Dirichlet pressure)
        if(cb.eq.'V  ') then
           do i=1,nxz         
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxbf(k),uybf(k),uzbf(k),prbf(k),temp1bf(k),yy)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = yy(b)
               wmeanbf = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             call eval_EOS_D(prbf(k),rhobf(k),temp1bf(k),wmeanbf)
             temp0 = temp1bf(k)
             ubf = 0. ! internal energy
             do b=1,nspec
               j = (b-1)*(npa+1)
               ui1 = 0.
               do c=1,npa+1
                 ui1 = ui1 + bi(j+c)*temp0**(c-1)
               enddo
               ubf = ubf + yy(b) * ui1
             enddo
             spd1 = sqrt( uxbf(k)*uxbf(k)
     &                  + uybf(k)*uybf(k)
     &                  + uzbf(k)*uzbf(k) )
             kinE = 0.5*spd1**2.  ! kinetic energy
             totE = ubf + kinE    ! e=u+ke
             rhoE = rhobf(k)*totE ! rhoE=rho*e            
             call eval_dgeh(rhoE,prbf(k),rhobf(k),dgehbf(k))
             k = k+1
            enddo
        endif

ccc subsonic inflow (pressure from interior) - not working for premixed flame
        if(cb.eq.'VVV  ') then
           do i=1,nxz      
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxbf(k),uybf(k),uzbf(k),prbf(k),temp1bf(k),yy)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = yy(b)
               wmeanbf = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             prbf(k) = prpi(k)
             call eval_EOS_D(prbf(k),rhobf(k),temp1bf(k),wmeanbf) 
             temp0 = temp1bf(k)
             ubf = 0. ! internal energy
             do b=1,nspec
               j = (b-1)*(npa+1)
               ui1 = 0.
               do c=1,npa+1
                 ui1 = ui1 + bi(j+c)*temp0**(c-1)
               enddo
               ubf = ubf + yy(b) * ui1
             enddo
             spd1 = sqrt( uxbf(k)*uxbf(k)
     &                  + uybf(k)*uybf(k)
     &                  + uzbf(k)*uzbf(k) )
             kinE = 0.5*spd1**2.  ! kinetic energy
             totE = ubf + kinE    ! e=u+ke
             rhoE = rhobf(k)*totE ! rhoE=rho*e            
             call eval_dgeh(rhoE,prbf(k),rhobf(k),dgehbf(k))
             k = k+1
            enddo
        endif

ccc subsonic outflow - fixed pressure (modifies p and rho)
        if(cb.eq.'O  ') then
           do i=1,nxz
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxbf(k),uybf(k),uzbf(k),prbf(k),temp1bf(k),yy)
             uxbf(k)  = uxpi(k)          
             uybf(k)  = uypi(k)
             uzbf(k)  = uzpi(k)
             temp1bf(k) = temp1pi(k)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               wmeanbf  = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             call eval_EOS_D(prbf(k),rhobf(k),temp1bf(k),wmeanbf)
             dgehbf(k) = dgehpi(k)     
             k = k+1
           enddo
        endif

ccc subsonic outflow - Neumann for entire state vector
        if(cb.eq.'ON ') then
           do i=1,nxz
             prbf(k) = prpi(k)
             uxbf(k) = uxpi(k)          
             uybf(k) = uypi(k)
             uzbf(k) = uzpi(k)
             temp1bf(k) = temp1pi(k)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               wmeanbf  = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             rhobf(k) = rhopi(k)
             dgehbf(k) = dgehpi(k)
             k = k+1
           enddo
        endif


ccc characteristic subsonic inflow
        if(cb.eq.'v  ') then     
           do i=1,nxz

             ! exterior state (upstream inflow)
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxExt,uyExt,uzExt,pExt,temp0,yy)
             wmeanExt = 0.
             do b=1,nspec
               Yibf(k,b) = yy(b) 
               wmeanExt = wmeanExt + Yibf(k,b)/mw(b)
             enddo
             wmeanExt = 1./wmeanExt
             do b=1,nspec
               Xibf(k,b) = wmeanExt*Yibf(k,b)/mw(b)
             enddo
             call eval_EOS_D(pExt,rhoExt,temp0,wmeanExt) ! density
             cpExt = 0. ! cp, cv, gamma
             do b=1,nspec
               j = (b-1)*npa
               cpi = 0.
               do c=1,npa
                cpi = cpi + aa(j+c)*temp0**(c-1)
               enddo
               cpExt = cpExt + yy(b)*cpi
             enddo
             gammaExt = cpExt/(cpExt-R_u/wmeanExt)
             cvExt = cpExt - R_u/wmeanExt
             cpbarExt = 0. ! cpbar, gammabar
             do b=1,nspec
               j = (b-1)*npa
               cpbari = 0.
               do c=1,npa
                 cpbari = cpbari + aa(j+c)*(temp0**c-298.15**c)/c
               enddo
               cpbari = cpbari/(temp0-298.15)
               cpbarExt = cpbarExt + yy(b)*cpbari  
             enddo
             gammabarExt = cpbarExt/(cpbarExt-R_u/wmeanExt)
             cvbarExt = cpbarExt - R_u/wmeanExt
             !cExt = sqrt(gammaInt*temp0*R_u/wmeanExt)
             cExt = sqrt(gammaExt*temp0*R_u/wmeanExt)

             ! interior state
             wmeanInt = 0.
             do b=1,nspec
               yy(b) = Yipi(k,b)
               wmeanInt = wmeanInt + Yibf(k,b)/mw(b)
             enddo
             wmeanInt = 1./wmeanInt
             temp0 = temp1pi(k)
             cpInt = 0. ! cp, cv, gamma
             do b=1,nspec
               j = (b-1)*npa
               cpi = 0.
               do c=1,npa
                cpi = cpi + aa(j+c)*temp0**(c-1)
               enddo
               cpInt = cpInt + yy(b)*cpi
             enddo
             gammaInt = cpInt/(cpInt-R_u/wmeanInt)
             cpbarInt = 0. ! cpbar, cvbar, gammabar
             do b=1,nspec
               j = (b-1)*npa
               cpbari = 0.
               do c=1,npa
                 cpbari = cpbari + aa(j+c)*(temp0**c-298.15**c)/c
               enddo
               cpbari = cpbari/(temp0-298.15)
               cpbarInt = cpbarInt + yy(b)*cpbari  
             enddo
             gammabarInt = cpbarInt/(cpbarInt-R_u/wmeanInt)
             cvbarInt = cpbarInt - R_u/wmeanInt
             cInt = sqrt(gammaInt*temp1pi(k)*R_u/wmeanInt)
             unInt = uxpi(k)*unx(i,1,iface,e) ! vn
     &             + uypi(k)*uny(i,1,iface,e)  

             ! Riemann invariants
             w1 = unInt - 2.*cExt/(gammabarExt-1.)
             w2 = cvbarExt*log(pExt/(rhoExt**gammabarExt))
             w3 = 0.
             w4 = 0.
             w5 = unInt + 2.*cInt/(gammabarInt-1.)
             
             ! set thermodynamic state
             rhobf(k) = rhoExt
             cstar    = 0.25*(gammabarExt-1.)*(w5-w1)
             prbf(k)  = rhobf(k)*cstar**2./gammabarExt
             temp1bf(k) = prbf(k)*wmeanExt/(rhobf(k)*R_u)
             uxbf(k)  = 0.5*(w1+w5)*unx(i,1,iface,e)+uxExt-uxpi(k)         
             uybf(k)  = 0.
             uzbf(k)  = 0.

             ! set energy
             temp0 = temp1bf(k)
             ubf = 0. ! internal energy
             do b=1,nspec
               j = (b-1)*(npa+1)
               ui1 = 0.
               do c=1,npa+1
                 ui1 = ui1 + bi(j+c)*temp0**(c-1)
               enddo
               ubf = ubf + yy(b) * ui1
             enddo
             spd1 = sqrt( uxbf(k)*uxbf(k)
     &                  + uybf(k)*uybf(k)
     &                  + uzbf(k)*uzbf(k) )
             kinE = 0.5*spd1**2.  ! kinetic energy
             totE = ubf + kinE    ! e=u+ke
             rhoE = rhobf(k)*totE ! rhoE = rho * e             
             call eval_dgeh(rhoE,prbf(k),rhobf(k),dgehbf(k))

             ifWrite = 0
             if(k.eq.317.and.ifWrite.eq.1) then
               write(*,*)'Starting inflow...'
               write(*,*)'cpInt=',cpInt
               write(*,*)'gammaInt=',gammaInt
               write(*,*)'cpbarInt=',cpbarInt
               write(*,*)'gammabarInt=',gammabarInt
               write(*,*)'cvbarInt=',cvbarInt
               write(*,*)'cInt=',cInt
               write(*,*)'unInt=',unInt
               write(*,*)'rhoExt=',rhoExt
               write(*,*)'cpbarExt=',cpbarExt
               write(*,*)'gammabarExt=',gammabarExt
               write(*,*)'cExt=',cExt
               write(*,*)'w1=',w1
               write(*,*)'w2=',w2
               write(*,*)'w5=',w5
               write(*,*)'rhopi(k)=',rhopi(k)
               write(*,*)'rhobf(k)=',rhobf(k)
               write(*,*)'cstar=',cstar
               write(*,*)'prpi(k)=',prpi(k)
               write(*,*)'prbf(k)=',prbf(k)
               write(*,*)'temp1pi(k)=',temp1pi(k)
               write(*,*)'temp1bf(k)=',temp1bf(k)
               write(*,*)'uxbf(k)=',uxbf(k)
               write(*,*)'dgehpi(k)=',dgehpi(k)
               write(*,*)'dgehbf(k)=',dgehbf(k)
               write(*,*)''
             endif

             ! testing old BC - DELETE
             ifOldBC = 0
             if(ifOldBC.eq.1) then
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxbf(k),uybf(k),uzbf(k),prbf(k),temp1bf(k),yy)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = yy(b)
               wmeanbf = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             call eval_EOS_D(prbf(k),rhobf(k),temp1bf(k),wmeanbf)
             call eval_dgeh(consf(k,E_EQ),prbf(k),rhobf(k),dgehbf(k))
             endif            
             
             ifWrite = 1
             if(k.eq.317.and.ifWrite.eq.1) then
               write(*,*)'HERE old:'
               write(*,*)'uxbf(k)=',uxbf(k)
               write(*,*)''
             endif

             k = k+1
            enddo
        endif





ccc characteristic subsonic outflow
        if(cb.eq.'o  ') then
           do i=1,nxz
             
             ! exterior state (downstream outflow)
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxExt,uyExt,uzExt,pExt,temp0,yy)
             wmeanExt = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               yy(b) = Yibf(k,b)
               wmeanExt = wmeanExt + Yibf(k,b)/mw(b)
             enddo
             wmeanExt = 1./wmeanExt
             do b=1,nspec
               Xibf(k,b) = wmeanExt*Yibf(k,b)/mw(b)
               xx(b) = Xibf(k,b)
             enddo            
             
             ! interior state
             temp0 = temp1pi(k)
             cpInt = 0. ! cp, cv, gamma
             do b=1,nspec
               j = (b-1)*npa
               cpi = 0.
               do c=1,npa
                cpi = cpi + aa(j+c)*temp0**(c-1)
               enddo
               cpInt = cpInt + yy(b)*cpi
             enddo
             gammaInt = cpInt/(cpInt-R_u/wmeanExt)
             cpbarInt = 0. ! cpbar, cvbar, gammabar
             do b=1,nspec
               j = (b-1)*npa
               cpbari = 0.
               do c=1,npa
                 cpbari = cpbari + aa(j+c)*(temp0**c-298.15**c)/c
               enddo
               cpbari = cpbari/(temp0-298.15)
               cpbarInt = cpbarInt + yy(b)*cpbari  
             enddo
             gammabarInt = cpbarInt/(cpbarInt-R_u/wmeanExt)
             cvbarInt = cpbarInt - R_u/wmeanExt
             !cInt = sqrt(gammaInt*temp1pi(k)*R_u/wmeanExt)
             cInt = sqrt(gammabarInt*temp1pi(k)*R_u/wmeanExt) ! testing delete
             unInt = uxpi(k)*unx(i,1,iface,e) ! vn
     &             + uypi(k)*uny(i,1,iface,e)             
             
             ! exterior state (outflow)
             new_max = 1000 ! Newton's method for T=T(s)s
             new_tol = 1.e-8
             dTT = 1.
             f = 0
             temp0 = 500.
             do while(abs(dTT).gt.new_tol.and.f.le.new_max)
               f = f + 1
               ss0 = -ssf(k)
               dsdT = 0.
               do b=1,nspec
                 j = (b-1)*npa
                 ss0 = ss0 + yy(b)*s0(b)+yy(b)*aa(j+1)*log(temp0/298.15)
     &                     - yy(b)*(R_u/mw(b))*log(pExt/101325.)
     &                     - yy(b)*(R_u/mw(b))*log(max(xx(b),1.e-16))
                 dsdT = dsdT + yy(b)*aa(j+1)/temp0
                 do c=2,npa
                   ss0 = ss0 + yy(b)*aa(j+c)*
     &                       (temp0**(c-1)-298.15**(c-1))/(c-1)
                   dsdT = dsdT + yy(b)*aa(j+c)*temp0**(c-2)
                 enddo
               enddo
               dTT = -ss0/dsdT
               temp0 = temp0 + dTT
             enddo
             call eval_EOS_D(pExt,rhoExt,temp0,wmeanExt) ! density
             cpbarExt = 0. ! cpbar, gammabar
             do b=1,nspec
               j = (b-1)*npa
               cpbari = 0.
               do c=1,npa
                 cpbari = cpbari + aa(j+c)*(temp0**c-298.15**c)/c
               enddo
               cpbari = cpbari/(temp0-298.15)
               cpbarExt = cpbarExt + yy(b)*cpbari  
             enddo
             gammabarExt = cpbarExt/(cpbarExt-R_u/wmeanExt)
             !cExt = sqrt(gammaInt*temp0*R_u/wmeanExt)
             cExt = sqrt(gammabarInt*temp0*R_u/wmeanExt)
             
             ! Riemann invariants
             w2 = cvbarInt*log(prpi(k)/(rhopi(k)**gammabarInt))
             w3 = 0.
             w4 = 0.
             w5 = unInt + 2.*cInt/(gammabarInt-1.)
             w1=w5-(4./(gammabarExt-1.))*sqrt(gammabarExt*pExt/rhoExt)                      
             !w1=w5-(4./(gammabarExt-1.))*sqrt(gammaInt*pExt/rhoExt) ! testing delete

             ! set thermodynamic state
             rhobf(k) = ((cExt**2./gammabarInt)*exp(-w2/cvbarInt))
     &                  **(1./(gammabarInt-1.))
             cstar    = 0.25*(gammabarInt-1.)*(w5-w1)
             prbf(k)  = rhobf(k)*cstar**2./gammabarInt
             !prbf(k)  = rhobf(k)*cstar**2./gammaInt ! testing delete
             temp1bf(k) = prbf(k)*wmeanExt/(rhobf(k)*R_u)
             uxbf(k)  = 0.5*(w1+w5)*unx(i,1,iface,e)          
             uybf(k)  = 0.
             uzbf(k)  = 0.            
             
             ! set energy
             temp0 = temp1bf(k)
             ubf = 0. ! internal energy
             do b=1,nspec
               j = (b-1)*(npa+1)
               ui1 = 0.
               do c=1,npa+1
                 ui1 = ui1 + bi(j+c)*temp0**(c-1)
               enddo
               ubf = ubf + yy(b) * ui1
             enddo
             spd1 = sqrt( uxbf(k)*uxbf(k)
     &                  + uybf(k)*uybf(k)
     &                  + uzbf(k)*uzbf(k) )
             kinE = 0.5*spd1**2.  ! kinetic energy
             totE = ubf + kinE    ! e=u+ke
             rhoE = rhobf(k)*totE ! rhoE = rho * e             
             call eval_dgeh(rhoE,prbf(k),rhobf(k),dgehbf(k))
             
             ifWrite = 1
             if(k.eq.80.and.ifWrite.eq.1) then
               write(*,*)'Starting outflow...'
               write(*,*)'cpInt=',cpInt
               write(*,*)'gammaInt=',gammaInt
               write(*,*)'cpbarInt=',cpbarInt
               write(*,*)'gammabarInt=',gammabarInt
               write(*,*)'cvbarInt=',cvbarInt
               write(*,*)'cInt=',cInt
               write(*,*)'unInt=',unInt
               write(*,*)'rhoExt=',rhoExt
               write(*,*)'cpbarExt=',cpbarExt
               write(*,*)'gammabarExt=',gammabarExt
               write(*,*)'cExt=',cExt
               write(*,*)'w1=',w1
               write(*,*)'w2=',w2
               write(*,*)'w5=',w5
               write(*,*)'rhopi(k)=',rhopi(k)
               write(*,*)'rhobf(k)=',rhobf(k)
               write(*,*)'cstar=',cstar
               write(*,*)'prpi(k)=',prpi(k)
               write(*,*)'prbf(k)=',prbf(k)
               write(*,*)'temp1pi(k)=',temp1pi(k)
               write(*,*)'temp1bf(k)=',temp1bf(k)
               write(*,*)'uxbf(k)=',uxbf(k)
               write(*,*)'uxpi(k)=',uxpi(k)
               write(*,*)'dgehpi(k)=',dgehpi(k)
               write(*,*)'dgehbf(k)=',dgehbf(k)
               write(*,*)''
             endif

             ! testing old BC - DELETE
             ifOldBC = 0
             if(ifOldBC.eq.1) then
             call userbc_comp(xm1f(k),ym1f(k),zm1f(k),id_face,cb
     $               ,uxbf(k),uybf(k),uzbf(k),prbf(k),temp1bf(k),yy)
             uxbf(k)  = uxpi(k)          
             uybf(k)  = uypi(k)
             uzbf(k)  = uzpi(k)
             temp1bf(k) = temp1pi(k)
             wmeanbf = 0.
             do b=1,nspec
               Yibf(k,b) = Yipi(k,b)
               wmeanbf  = wmeanbf + Yibf(k,b)/mw(b)
             enddo
             wmeanbf = 1./wmeanbf
             do b=1,nspec
               Xibf(k,b) = wmeanbf*Yibf(k,b)/mw(b)
             enddo
             call eval_EOS_D(prbf(k),rhobf(k),temp1bf(k),wmeanbf)
             ! call eval_dgeh(consf(k,E_EQ),prbf(k),rhobf(k),dgehbf(k)) ! old
             dgehbf(k) = dgehpi(k)
             endif

             k = k+1
           enddo
        endif








      enddo
      enddo

      return
      end
c---------------------------------------------------------------------      
      subroutine compbc_flux

c     This routine adds the ghost node values bf to the pi (U+U_neighbor) values  

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      integer i,j,e,iface,k,lnf
      integer ntot,nxz
      character cb*3
      
      nxz = nx1*nz1
      ntot= nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt

      do e=1,nelt
      do iface=1,2*ldim
      cb = cbc(iface,e,1)
      k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz
      if((cb.eq.'SYM').or.(cb.eq.'W  ').or.(cb.eq.'v  ').or.
     $   (cb.eq.'O  ').or.(cb.eq.'o  ').or.(cb.eq.'V  ').or.
     $   (cb.eq.'ON ') ) then
         do i=1,nxz
           rhopi(k) = rhopi(k)+rhobf(k)
           uxpi(k) = uxpi(k)+uxbf(k)
           uypi(k) = uypi(k)+uybf(k)
           if(if3d) uzpi(k) = uzpi(k)+uzbf(k) 
           prpi(k) = prpi(k)+prbf(k)
           dgehpi(k) = dgehpi(k)+dgehbf(k)
           temp1pi(k) = temp1pi(k)+temp1bf(k)
           do b=1,nspec
             Yipi(k,b) = Yipi(k,b)+Yibf(k,b)
             Xipi(k,b) = Xipi(k,b)+Xibf(k,b)
           enddo
           conspi(k,U_EQ) = conspi(k,U_EQ)+rhobf(k)*uxbf(k)
           conspi(k,V_EQ) = conspi(k,V_EQ)+rhobf(k)*uybf(k)
           if(if3d) conspi(k,W_EQ) = conspi(k,W_EQ)+rhobf(k)*uzbf(k)
           conspi(k,E_EQ) = conspi(k,E_EQ)+dgehbf(k)*rhobf(k)-prbf(k)
           do b=1,nspec
             conspi(k,E_EQ+b) = conspi(k,E_EQ+b)+rhobf(k)*Yibf(k,b)
           enddo          
           if(ifIntLLFfull) then
             flux1pi(k,U_EQ)=flux1pi(k,U_EQ) + rhobf(k)*uxbf(k)*uxbf(k)
     &                                       + prbf(k)
             flux2pi(k,U_EQ)=flux2pi(k,U_EQ) + rhobf(k)*uxbf(k)*uybf(k)
             if(if3d) then
               flux3pi(k,U_EQ)=flux3pi(k,U_EQ)+rhobf(k)*uxbf(k)*uzbf(k)
             endif
             flux1pi(k,V_EQ)=flux1pi(k,V_EQ) + rhobf(k)*uxbf(k)*uybf(k)
             flux2pi(k,V_EQ)=flux2pi(k,V_EQ) + rhobf(k)*uybf(k)*uybf(k)
     &                                       + prbf(k) 
             if(if3d) then
               flux3pi(k,V_EQ)=flux3pi(k,V_EQ)+rhobf(k)*uybf(k)*uzbf(k)
             endif
             flux1pi(k,E_EQ)=flux1pi(k,E_EQ)+rhobf(k)*uxbf(k)*dgehbf(k)
             flux2pi(k,E_EQ)=flux2pi(k,E_EQ)+rhobf(k)*uybf(k)*dgehbf(k)
             if(if3d) then
             flux3pi(k,E_EQ)=flux3pi(k,E_EQ)+rhobf(k)*uzbf(k)*dgehbf(k)
             endif
             if(if3d) then
               flux1pi(k,W_EQ)=flux1pi(k,W_EQ)+rhobf(k)*uxbf(k)*uzbf(k)
               flux2pi(k,W_EQ)=flux2pi(k,W_EQ)+rhobf(k)*uybf(k)*uzbf(k)
               flux3pi(k,W_EQ)=flux3pi(k,W_EQ)+rhobf(k)*uzbf(k)*uzbf(k)
     &                                       + prbf(k)
             endif
             do b=1,nspec
               flux1pi(k,E_EQ+b) = flux1pi(k,E_EQ+b)
     &                           + rhobf(k)*uxbf(k)*Yibf(k,b)
               flux2pi(k,E_EQ+b) = flux2pi(k,E_EQ+b)
     &                           + rhobf(k)*uybf(k)*Yibf(k,b)
               if(if3d) then
                 flux3pi(k,E_EQ+b) = flux3pi(k,E_EQ+b)
     &                             + rhobf(k)*uzbf(k)*Yibf(k,b)
               endif
             enddo 
           endif
           k = k+1
         enddo
      endif
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------   



c---------------------------------------------------------------------
c                BCs for viscous fluxes
c---------------------------------------------------------------------  
      subroutine compbc_visc(dgface,bctype)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real dgface(lf)
      integer i,j,e,iface,k,lnf
      integer ntot,nxz,dir
      character cb*3

      character bctype*2
      
      nxz = nx1*nz1
      ntot= nx1*ny1*nz1*nelt
      lnf = nx1*nz1*2*ndim*nelt

      
      if(bctype.eq.'AV') then

        do e=1,nelt
        do iface=1,2*ldim
          cb = cbc(iface,e,1)
          k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz
          if((cb.eq.'SYM').or.(cb.eq.'W  ').or.(cb.eq.'v  ').or.
     $       (cb.eq.'O  ').or.(cb.eq.'o  ').or.(cb.eq.'V  ').or.
     $       (cb.eq.'ON ') ) then
             do i=1,nxz
               dgface(k) = 0.
               k = k+1
             enddo
          endif
        enddo
        enddo

      else if(bctype.eq.'E1') then
      
        do e=1,nelt
        do iface=1,2*ldim
          cb = cbc(iface,e,1)
          k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz
          if((cb.eq.'v  ').or.(cb.eq.'V  ')) then
            do i=1,nxz
              dgface(k) = 2.*dgface(k)
              k = k+1
            enddo
          else if((cb.eq.'o  ').or.(cb.eq.'O  ').or.(cb.eq.'ON ')) then
            do i=1,nxz
              dgface(k) = 2.*dgface(k)
              k = k+1
            enddo   
          else if((cb.eq.'SYM').or.(cb.eq.'W  ')) then
            do i=1,nxz
              dgface(k) = 0.
              k = k+1
            enddo          
          endif
        enddo
        enddo

      else if(bctype.eq.'Y1') then

        do e=1,nelt
        do iface=1,2*ldim
          cb = cbc(iface,e,1)
          k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz
          if((cb.eq.'v  ').or.(cb.eq.'V  ')) then
            do i=1,nxz
              dgface(k) = 2.*dgface(k)
              k = k+1
            enddo
          else if((cb.eq.'o  ').or.(cb.eq.'O  ').or.(cb.eq.'ON ')) then
            do i=1,nxz
              dgface(k) = 2.*dgface(k)
              k = k+1
            enddo   
          else if((cb.eq.'SYM').or.(cb.eq.'W  ')) then
            do i=1,nxz
              dgface(k) = 0.
              k = k+1
            enddo          
          endif
        enddo
        enddo

      else if(bctype.eq.'U1') then

        do e=1,nelt
        do iface=1,2*ldim
          cb = cbc(iface,e,1)
          k = 1 + nxz*(iface-1) + (e-1)*2*ldim*nxz
          if((cb.eq.'SYM').or.(cb.eq.'W  ').or.(cb.eq.'v  ').or.
     $       (cb.eq.'O  ').or.(cb.eq.'o  ').or.(cb.eq.'V  ').or.
     $       (cb.eq.'ON ') ) then
             do i=1,nxz
               dgface(k) = 2.*dgface(k) ! 0.*dgface(k)
               k = k+1
             enddo
          endif
        enddo
        enddo

      endif

      return
      end
c---------------------------------------------------------------------  
