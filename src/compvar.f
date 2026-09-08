c---------------------------------------------------------------------
      subroutine eval_dependent_vars(rkU)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'   

      real rkU(lx1,ly1,lz1,lelt,N_EQ)
      integer n,i,j,f

      real yy(nspec)
      real xx(nspec)
      real totE,kinE
      real wdot1(nspec)
      real temp0
      integer new_max
      real new_tol

      real viscosity,thermalconductivity,diff0(nspec)

      real ppi

      integer icalld
      save icalld
      data icalld /0/

      n = nx1*ny1*nz1*nelt

      do i=1,n

ccc density
        vtrans(i,1,1,1,1) = 0.
        do b=1,nspec
           vtrans(i,1,1,1,1) = vtrans(i,1,1,1,1) + rkU(i,1,1,1,E_EQ+b)
        enddo

ccc clip and re-scale species
        if(ifclip) then
           ysum = 0.
           do b=1,nspec
             Yi(i,1,1,1,b) = 
     &            max(rkU(i,1,1,1,E_EQ+b),0.)/vtrans(i,1,1,1,1)
             ysum = ysum + Yi(i,1,1,1,b)
           enddo
           do b=1,nspec
             Yi(i,1,1,1,b) = Yi(i,1,1,1,b)/ysum
             rkU(i,1,1,1,E_EQ+b) = vtrans(i,1,1,1,1)*Yi(i,1,1,1,b)
          enddo
        endif

ccc species
        wminv = 0.
        do b=1,nspec
           Yi(i,1,1,1,b)=rkU(i,1,1,1,E_EQ+b)/vtrans(i,1,1,1,1)
           yy(b) = Yi(i,1,1,1,b)
           wminv = wminv + yy(b)/mw(b)
        enddo
        wmean(i,1,1,1) = 1./wminv
        do b=1,nspec
          Xi(i,1,1,1,b) = wmean(i,1,1,1)*Yi(i,1,1,1,b)/mw(b)
          xx(b) = Xi(i,1,1,1,b)
        enddo
        
ccc velocity
        vx(i,1,1,1)=rkU(i,1,1,1,U_EQ)/vtrans(i,1,1,1,1)
        if(if1D) rkU(i,1,1,1,V_EQ) = 0.
        vy(i,1,1,1)=rkU(i,1,1,1,V_EQ)/vtrans(i,1,1,1,1)
        if(ldim.ne.3) rkU(i,1,1,1,W_EQ) = 0. 
        vz(i,1,1,1)=rkU(i,1,1,1,W_EQ)/vtrans(i,1,1,1,1)
        spd(i,1,1,1)=sqrt(vx(i,1,1,1)*vx(i,1,1,1)+
     $                    vy(i,1,1,1)*vy(i,1,1,1)+
     $                    vz(i,1,1,1)*vz(i,1,1,1))

ccc energy
        kinE = 0.5*spd(i,1,1,1)**2.                 ! [J/kg]
        totE = rkU(i,1,1,1,E_EQ)/vtrans(i,1,1,1,1)  ! [J/kg]
        uu(i,1,1,1) = totE-kinE                     ! [J/kg]

ccc Newton's method for T
        new_max = 1000
        new_tol = 1.e-8
        dTT = 1.
        f = 0
        temp0 = 500.
        starting_temp = temp0
        do while(abs(dTT).gt.new_tol.and.f.le.new_max)
          f = f + 1
          u0 = -uu(i,1,1,1)
          dudT = 0.
          do b=1,nspec
            j = (b-1)*(npa+1)
            do c=1,npa+1
              u0 = u0 + yy(b)*bi(j+c)*temp0**(c-1)
              dudT = dudT + yy(b)*(c-1)*bi(j+c)*temp0**(c-2)
            enddo
          enddo
          dTT = -u0/dudT
          temp0 = temp0 + dTT
        enddo

ccc temperature clipping
        if(ifclip) then
           if(temp0.lt.298) temp0=298. ! clipping T
        endif

ccc temp/pressure/dgeh
        temp1(i,1,1,1) = temp0
        call eval_EOS_P(pr(i,1,1,1),vtrans(i,1,1,1,1)
     $              ,temp1(i,1,1,1),wmean(i,1,1,1))
        call eval_dgeh(rkU(i,1,1,1,E_EQ),pr(i,1,1,1)
     $              ,vtrans(i,1,1,1,1),dgeh(i,1,1,1))

ccc cp/enthalpy
        cp1 = 0.
        do b=1,nspec
          j = (b-1)*npa
          cpi = 0.
          hh(i,1,1,1,b) = h0(b)
          do c=1,npa
            cpi = cpi + aa(j+c)*temp0**(c-1)
            hh(i,1,1,1,b) = hh(i,1,1,1,b)
     &                    + aa(j+c)*(temp0**c-298.15**c)/c
          enddo
          cp1 = cp1 + yy(b)*cpi
        enddo

ccc internal energy
        ushift(i,1,1,1) = uu(i,1,1,1)
        do b=1,nspec
          j = (b-1)*(npa+1)
          ui(i,1,1,1,b) = 0.
          do c=1,npa+1
            ui(i,1,1,1,b) = ui(i,1,1,1,b)+bi(j+c)*temp0**(c-1)
          enddo
          ushift(i,1,1,1) = ushift(i,1,1,1) - yy(b)*bi(j+1)
        enddo

ccc gamma/c/mach
        gamma(i,1,1,1) = cp1/(cp1-R_u/wmean(i,1,1,1))
        aspd(i,1,1,1)=sqrt(pr(i,1,1,1)*
     $                   gamma(i,1,1,1)/vtrans(i,1,1,1,1))       
        mach(i,1,1,1) = spd(i,1,1,1)/aspd(i,1,1,1)

ccc reaction rates
        if(ifforces.and.(.not.ifsplitchem)) then
           call setmassfractions(yy)
           call setstate_uv(uu(i,1,1,1), 1./vtrans(i,1,1,1,1))
           call getnetproductionrates(wdot1)
           wdot_sum = 0.
           do b=1,nspec
             wdot(i,1,1,1,b) = mw(b)*wdot1(b)*1.e3
             wdot_sum = wdot_sum + wdot(i,1,1,1,b)
           enddo
           wdot(i,1,1,1,nspec) = wdot(i,1,1,1,nspec) - wdot_sum 
        endif

ccc viscosity/thermal conductivity/mass diffusivity
        if(ifviscous) then
           call setmassfractions(yy)
           call setstate_uv(uu(i,1,1,1), 1./vtrans(i,1,1,1,1))
           vis1(i,1,1,1)  = viscosity()
           cond1(i,1,1,1) = thermalconductivity()
           call getmixdiffcoeffs(diff0)
           do b=1,nspec ! b cannot be declared as integer
             diff1(i,1,1,1,b) = vtrans(i,1,1,1,1)*diff0(b)
           enddo
        endif

ccc entropy (multi-component calorically imperfect gas)
        ss(i,1,1,1) = 0.
        do b=1,nspec
          ppi = vtrans(i,1,1,1,1)*max(1.e-16,yy(b))*R_u*temp0/mw(b) ! partial pressure
          j = (b-1)*npa
          si(i,1,1,1,b) = s0(b)+aa(j+1)*log(temp0/298.15)
     &                  - (R_u/mw(b))*log(ppi/101325.) 
          do c=2,npa
            si(i,1,1,1,b) = si(i,1,1,1,b) 
     &                    + aa(j+c)*(temp0**(c-1)-298.15**(c-1))/(c-1)
          enddo       
          ss(i,1,1,1) = ss(i,1,1,1) + yy(b)*si(i,1,1,1,b)
        enddo
        rhoS(i,1,1,1) = vtrans(i,1,1,1,1)*ss(i,1,1,1)

ccc entropy (calorically perfect gas)
c        rhoS(i,1,1,1) = vtrans(i,1,1,1,1)/(gamma(i,1,1,1)-1.)*
c     $      log(pr(i,1,1,1)/(vtrans(i,1,1,1,1)**gamma(i,1,1,1)))

ccc testing cpbar and gammabar: u=p/(rho*(gammabar-1))
c        cpbar = 0.
c        do b=1,nspec
c          j = (b-1)*npa
c          hi = h0(b)
c          do c=1,npa
c            hi = hi + aa(j+c)*(temp0**c-298.15**c)/c
c          enddo
c          cpbari = hi/temp0
c          cpbar = cpbar + yy(b)*cpbari       
c        enddo
c        gammabar = cpbar/(cpbar-R_u/wmean(i,1,1,1))
cc        t(i,1,1,1,18) = pr(i,1,1,1)/(vtrans(i,1,1,1,1)*(gammabar-1.))
cc        t(i,1,1,1,19) = uu(i,1,1,1)
c        t(i,1,1,1,18) = gammabar
c        t(i,1,1,1,19) = gamma(i,1,1,1)

ccc testing cpbar and gammabar: u=p/(rho*(gammabar-1))+sum(Yi*[hi0-cpbari*T0])
c        t(i,1,1,1,18) = 0.
c        cpbar = 0.
c        do b=1,nspec
c          j = (b-1)*npa
c          cpbari = 0.
c          do c=1,npa
c            cpbari = cpbari + aa(j+c)*(temp0**c-298.15**c)/c
c          enddo
c          cpbari = cpbari/(temp0-298.15)
c          cpbar = cpbar + yy(b)*cpbari  
c          t(i,1,1,1,18) = t(i,1,1,1,18) + yy(b)*(h0(b)-cpbari*298.15)        
c        enddo
c        gammabar = cpbar/(cpbar-R_u/wmean(i,1,1,1))
c        t(i,1,1,1,18) = t(i,1,1,1,18) 
c     &                + pr(i,1,1,1)/(vtrans(i,1,1,1,1)*(gammabar-1.))
c        t(i,1,1,1,19) = uu(i,1,1,1)
c        t(i,1,1,1,18) = gammabar
c        t(i,1,1,1,19) = gamma(i,1,1,1)


      enddo


ccc minimum entropy in domain at t=0
      if(icalld.eq.0) then
        s0min = glmin(ss,n)
        icalld = 1
      endif


      return
      end
c---------------------------------------------------------------------
      subroutine eval_EOS_D(PP, DD, TT, W)
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      real PP, DD, TT, W
      DD=PP*W/(R_u*TT)
      return
      end
c---------------------------------------------------------------------
      subroutine eval_EOS_T(PP, DD, TT, W)
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      real PP, DD, TT, W
      TT=PP*W/(R_u*DD)
      return
      end
c---------------------------------------------------------------------
      subroutine eval_EOS_P(PP, DD, TT, W)
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      real PP, DD, TT, W
      PP=DD*TT*R_u/W
      return
      end
c---------------------------------------------------------------------
      subroutine eval_dgeh(rhoE, PP, DD, dgeh1)
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      real rhoE, PP, DD, dgeh1
      dgeh1=(rhoE+PP)/DD
      return
      end
c---------------------------------------------------------------------
      subroutine get_bi()

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      integer b,j,g,c  

      do b=1,nspec

        sum1 = 0
        j = (b-1)*npa        ! index for aa
        do c=1,npa
          sum1 = sum1 + aa(j+c)*(298.15**c)/c
        enddo
        
        g = (b-1)*(npa+1)    ! index for bi
        do c=1,npa+1
          if(c.eq.1) then
             bi(g+c) = h0(b) - sum1
          else if(c.eq.2) then
             bi(g+c) = aa(j+c-1)-R_u/mw(b)
          else
             bi(g+c) = aa(j+c-1)/(c-1)
          endif
        enddo  
        
      enddo

      return
      end   
c---------------------------------------------------------------------

