c---------------------------------------------------------------------
c     HLL / HLLC interface fluxes, and a blend of the two.
c
c     Selected with userparam05 (default 0 = original LLF path):
c        1 = HLL      2 = HLLC      3 = HLLC blended into HLL
c     Blend (3): alpha = clamp((s_p-sw0)/(sw1-sw0),0,1), F = alpha*F_HLL
c     + (1-alpha)*F_HLLC, with the face pressure-jump sensor
c     s_p = |p_r-p_l|/(p_r+p_l). sw0, sw1 are userparam06, userparam07
c     (defaults 0.01 and 0.1). alpha is further gated by two guards:
c       (a) Mach gate: alpha=0 if either side is subsonic (preserves HLLC
c           at flame fronts, contacts, and pre-shock stratified layers)
c       (b) Density contact sensor: alpha *= max(0,1-s_rho/sw_rho) where
c           s_rho = |rho_r-rho_l|/(rho_r+rho_l), sw_rho = userparam08
c           (default 0.1). Suppresses HLL at sharp density interfaces so
c           HLLC can preserve stratified fuel/air contacts even when a
c           shock drives a large pressure jump at that face.
c     alpha is a function of the two face states only, so both elements
c     see the same flux (conservative).
c
c     State vector order is the solver's: (rhoU,rhoV,rhoW,rhoE,rhoYi).
c     Mass is carried by the species equations (rho = sum rhoYi).
c---------------------------------------------------------------------
      subroutine advect_interface_riemann(rkRHS)

ccc Interface contribution (F* - F^-).n to div(F); same sign and
ccc scaling as advect_interface, with F* from an approximate Riemann
ccc solver. Left state = own face values, right = (U+U_nb) - own.

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)

      integer e,iface,i,k,b,ieq,ig,ifl
      real sw0,sw1,sw_rho,xn,yn,zn
      real rho_l,u_l,v_l,w_l,p_l,e_l,g_l,y_l(nspec)
      real rho_r,u_r,v_r,w_r,p_r,e_r,g_r,y_r(nspec)
      real fsn(N_EQ),fln(N_EQ)
      real un_l,un_r,hfac

      ifl    = nint(uparam(5))
      sw0    = uparam(6)
      sw1    = uparam(7)
      sw_rho = uparam(8)
      if(sw0.le.0.)    sw0    = 0.01
      if(sw1.le.sw0)   sw1    = 10.*sw0
      if(sw_rho.le.0.) sw_rho = 0.1

      k = 0
      do e=1,nelt
       do iface=1,2*ldim
        do i=1,lx1*lz1
         k  = k+1
         ig = dg_face(k)

         xn = unx(i,1,iface,e)
         yn = uny(i,1,iface,e)
         zn = 0.
         if(if3d) zn = unz(i,1,iface,e)

         rho_l = rhof(k)
         u_l   = uxf(k)
         v_l   = uyf(k)
         w_l   = 0.
         if(if3d) w_l = uzf(k)
         p_l   = prf(k)
         e_l   = consf(k,E_EQ)
         g_l   = gamf(k)

         rho_r = rhopi(k)-rhof(k)
         u_r   = uxpi(k)-uxf(k)
         v_r   = uypi(k)-uyf(k)
         w_r   = 0.
         if(if3d) w_r = uzpi(k)-uzf(k)
         p_r   = prpi(k)-prf(k)
         e_r   = conspi(k,E_EQ)-consf(k,E_EQ)
         g_r   = gampi(k)-gamf(k)

         do b=1,nspec
           y_l(b) = Yif(k,b)
           y_r(b) = Yipi(k,b)-Yif(k,b)
         enddo

         call riemann_flux(fsn,ifl,sw0,sw1,sw_rho,xn,yn,zn
     &        ,rho_l,u_l,v_l,w_l,p_l,e_l,g_l,y_l
     &        ,rho_r,u_r,v_r,w_r,p_r,e_r,g_r,y_r)

ccc physical normal flux of the interior state
         un_l = u_l*xn+v_l*yn+w_l*zn
         fln(U_EQ) = rho_l*u_l*un_l+p_l*xn
         fln(V_EQ) = rho_l*v_l*un_l+p_l*yn
         fln(W_EQ) = rho_l*w_l*un_l+p_l*zn
         fln(E_EQ) = (e_l+p_l)*un_l
         do b=1,nspec
           fln(E_EQ+b) = rho_l*y_l(b)*un_l
         enddo

         hfac = area(i,1,iface,e)/bm1(ig,1,1,1)
         do ieq=1,N_EQ
           rkRHS(ig,1,1,1,ieq) = rkRHS(ig,1,1,1,ieq)
     &                         + (fsn(ieq)-fln(ieq))*hfac
         enddo

        enddo
       enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine riemann_flux(fsn,ifl,sw0,sw1,sw_rho,xn,yn,zn
     &     ,rho_l,u_l,v_l,w_l,p_l,e_l,g_l,y_l
     &     ,rho_r,u_r,v_r,w_r,p_r,e_r,g_r,y_r)

ccc Numerical flux through a face with unit normal (xn,yn,zn), from the
ccc left state L to the right state R. e* is total energy per volume.

      include 'SIZE'
      include 'COMPRESS'

      real fsn(N_EQ)
      integer ifl
      real sw0,sw1,sw_rho,xn,yn,zn
      real rho_l,u_l,v_l,w_l,p_l,e_l,g_l,y_l(nspec)
      real rho_r,u_r,v_r,w_r,p_r,e_r,g_r,y_r(nspec)

      real fhll(N_EQ),fhllc(N_EQ)
      real un_l,un_r,a_l,a_r,s_l,s_r,s_p,s_rho,alpha
      integer b

      un_l = u_l*xn+v_l*yn+w_l*zn
      un_r = u_r*xn+v_r*yn+w_r*zn
      a_l  = sqrt(g_l*p_l/rho_l)
      a_r  = sqrt(g_r*p_r/rho_r)

ccc Davis wave speed estimates
      s_l = min(un_l-a_l,un_r-a_r)
      s_r = max(un_l+a_l,un_r+a_r)

      alpha = 0.
      if(ifl.eq.1) then
         alpha = 1.
      elseif(ifl.eq.3) then
ccc pressure-jump sensor
         s_p   = abs(p_r-p_l)/(p_r+p_l)
         alpha = min(1.,max(0.,(s_p-sw0)/(sw1-sw0)))
ccc Mach gate: keep HLLC (alpha=0) if either side is subsonic
         if(abs(un_l)/a_l.lt.1..or.abs(un_r)/a_r.lt.1.) alpha = 0.
ccc density contact sensor: suppress HLL at sharp density interfaces
         if(alpha.gt.0.) then
            s_rho = abs(rho_r-rho_l)/(rho_r+rho_l)
            alpha = alpha*max(0.,1.-s_rho/sw_rho)
         endif
      endif

      if(alpha.gt.0.) then
         call flux_hll(fhll,s_l,s_r,xn,yn,zn
     &        ,rho_l,u_l,v_l,w_l,p_l,e_l,y_l
     &        ,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)
      endif
      if(alpha.lt.1.) then
         call flux_hllc(fhllc,s_l,s_r,xn,yn,zn
     &        ,rho_l,u_l,v_l,w_l,p_l,e_l,y_l
     &        ,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)
      endif

      if(alpha.le.0.) then
         call copy(fsn,fhllc,N_EQ)
      elseif(alpha.ge.1.) then
         call copy(fsn,fhll,N_EQ)
      else
         do b=1,N_EQ
           fsn(b) = alpha*fhll(b)+(1.-alpha)*fhllc(b)
         enddo
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine phys_flux(f,xn,yn,zn,rho,u,v,w,p,e,y)

ccc F(U).n for one state

      include 'SIZE'
      include 'COMPRESS'

      real f(N_EQ),xn,yn,zn,rho,u,v,w,p,e,y(nspec)
      real un
      integer b

      un = u*xn+v*yn+w*zn
      f(U_EQ) = rho*u*un+p*xn
      f(V_EQ) = rho*v*un+p*yn
      f(W_EQ) = rho*w*un+p*zn
      f(E_EQ) = (e+p)*un
      do b=1,nspec
        f(E_EQ+b) = rho*y(b)*un
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine cons_state(q,rho,u,v,w,e,y)

      include 'SIZE'
      include 'COMPRESS'

      real q(N_EQ),rho,u,v,w,e,y(nspec)
      integer b

      q(U_EQ) = rho*u
      q(V_EQ) = rho*v
      q(W_EQ) = rho*w
      q(E_EQ) = e
      do b=1,nspec
        q(E_EQ+b) = rho*y(b)
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine flux_hll(f,s_l,s_r,xn,yn,zn
     &     ,rho_l,u_l,v_l,w_l,p_l,e_l,y_l,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)

      include 'SIZE'
      include 'COMPRESS'

      real f(N_EQ),s_l,s_r,xn,yn,zn
      real rho_l,u_l,v_l,w_l,p_l,e_l,y_l(nspec)
      real rho_r,u_r,v_r,w_r,p_r,e_r,y_r(nspec)

      real f_l(N_EQ),f_r(N_EQ),q_l(N_EQ),q_r(N_EQ)
      integer b

      call phys_flux(f_l,xn,yn,zn,rho_l,u_l,v_l,w_l,p_l,e_l,y_l)
      call phys_flux(f_r,xn,yn,zn,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)

      if(s_l.ge.0.) then
         call copy(f,f_l,N_EQ)
      elseif(s_r.le.0.) then
         call copy(f,f_r,N_EQ)
      else
         call cons_state(q_l,rho_l,u_l,v_l,w_l,e_l,y_l)
         call cons_state(q_r,rho_r,u_r,v_r,w_r,e_r,y_r)
         do b=1,N_EQ
           f(b) = (s_r*f_l(b)-s_l*f_r(b)
     &             +s_l*s_r*(q_r(b)-q_l(b)))/(s_r-s_l)
         enddo
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine flux_hllc(f,s_l,s_r,xn,yn,zn
     &     ,rho_l,u_l,v_l,w_l,p_l,e_l,y_l,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)

ccc Toro's HLLC: contact speed sM, star states in the two fans.

      include 'SIZE'
      include 'COMPRESS'

      real f(N_EQ),s_l,s_r,xn,yn,zn
      real rho_l,u_l,v_l,w_l,p_l,e_l,y_l(nspec)
      real rho_r,u_r,v_r,w_r,p_r,e_r,y_r(nspec)

      real fK(N_EQ),qK(N_EQ),qs(N_EQ)
      real un_l,un_r,sM,sK,fac,rhoK,uK,vK,wK,pK,eK,unK,yK(nspec)
      integer b

      un_l = u_l*xn+v_l*yn+w_l*zn
      un_r = u_r*xn+v_r*yn+w_r*zn

      if(s_l.ge.0.) then
         call phys_flux(f,xn,yn,zn,rho_l,u_l,v_l,w_l,p_l,e_l,y_l)
         return
      elseif(s_r.le.0.) then
         call phys_flux(f,xn,yn,zn,rho_r,u_r,v_r,w_r,p_r,e_r,y_r)
         return
      endif

      sM = (p_r-p_l+rho_l*un_l*(s_l-un_l)-rho_r*un_r*(s_r-un_r))
     &   / (rho_l*(s_l-un_l)-rho_r*(s_r-un_r))

      if(sM.ge.0.) then
         sK=s_l
         rhoK=rho_l
         uK=u_l
         vK=v_l
         wK=w_l
         pK=p_l
         eK=e_l
         unK=un_l
         do b=1,nspec
           yK(b)=y_l(b)
         enddo
      else
         sK=s_r
         rhoK=rho_r
         uK=u_r
         vK=v_r
         wK=w_r
         pK=p_r
         eK=e_r
         unK=un_r
         do b=1,nspec
           yK(b)=y_r(b)
         enddo
      endif

      call phys_flux(fK,xn,yn,zn,rhoK,uK,vK,wK,pK,eK,yK)
      call cons_state(qK,rhoK,uK,vK,wK,eK,yK)

      fac = rhoK*(sK-unK)/(sK-sM)
      qs(U_EQ) = fac*(uK+(sM-unK)*xn)
      qs(V_EQ) = fac*(vK+(sM-unK)*yn)
      qs(W_EQ) = fac*(wK+(sM-unK)*zn)
      qs(E_EQ) = fac*(eK/rhoK+(sM-unK)*(sM+pK/(rhoK*(sK-unK))))
      do b=1,nspec
        qs(E_EQ+b) = fac*yK(b)
      enddo

      do b=1,N_EQ
        f(b) = fK(b)+sK*(qs(b)-qK(b))
      enddo

      return
      end
c---------------------------------------------------------------------
