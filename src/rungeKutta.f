c---------------------------------------------------------------------      
      subroutine SSPRK3(dtf)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      real rktp1(lx1,ly1,lz1,lelt,N_EQ)
      
      real dtf
      
      integer i,ntot
      ntot  = nx1*ny1*nz1*nelt

      if(istep.eq.1) call eval_dependent_vars(t)
      call copy(rhoS2,rhoS1,ntot)   
      do i=1,ntot
        rhoS1(i,1,1,1) = rhoS(i,1,1,1)
      enddo

ccc stage 1
      call getdgrhs(rkRHS,t)
      
      do i=1,ntot
        rktp1(i,1,1,1,U_EQ) = t(i,1,1,1,U_EQ)-dtf*rkRHS(i,1,1,1,U_EQ)
        rktp1(i,1,1,1,V_EQ) = t(i,1,1,1,V_EQ)-dtf*rkRHS(i,1,1,1,V_EQ)
        if(if3d) then
          rktp1(i,1,1,1,W_EQ) = t(i,1,1,1,W_EQ)-dtf*rkRHS(i,1,1,1,W_EQ)
        endif
        rktp1(i,1,1,1,E_EQ) = t(i,1,1,1,E_EQ)-dtf*rkRHS(i,1,1,1,E_EQ)
        do b=1,nspec
          rktp1(i,1,1,1,E_EQ+b) = t(i,1,1,1,E_EQ+b)
     $                          -dtf*rkRHS(i,1,1,1,E_EQ+b)
        enddo
      enddo

      if(ifpplim) call pp_limiter(rktp1)

ccc stage 2
      call getdgrhs(rkRHS,rktp1)
      
      do i=1,ntot
        rktp1(i,1,1,1,U_EQ) = 0.75*t(i,1,1,1,U_EQ)
     $                      + 0.25*rktp1(i,1,1,1,U_EQ)
     $                      - 0.25*dtf*rkRHS(i,1,1,1,U_EQ)
        rktp1(i,1,1,1,V_EQ) = 0.75*t(i,1,1,1,V_EQ)
     $                      + 0.25*rktp1(i,1,1,1,V_EQ)
     $                      - 0.25*dtf*rkRHS(i,1,1,1,V_EQ)
        if(if3d) then
          rktp1(i,1,1,1,W_EQ) = 0.75*t(i,1,1,1,W_EQ)
     $                        + 0.25*rktp1(i,1,1,1,W_EQ)
     $                        - 0.25*dtf*rkRHS(i,1,1,1,W_EQ)
        endif
        rktp1(i,1,1,1,E_EQ) = 0.75*t(i,1,1,1,E_EQ)
     $                      + 0.25*rktp1(i,1,1,1,E_EQ)
     $                      - 0.25*dtf*rkRHS(i,1,1,1,E_EQ)
        do b=1,nspec 
          rktp1(i,1,1,1,E_EQ+b) = 0.75*t(i,1,1,1,E_EQ+b)
     $                          + 0.25*rktp1(i,1,1,1,E_EQ+b)
     $                          - 0.25*dtf*rkRHS(i,1,1,1,E_EQ+b)
        enddo
      enddo

      if(ifpplim) call pp_limiter(rktp1)

ccc stage 3
      call getdgrhs(rkRHS,rktp1)
      
      do i=1,ntot
       t(i,1,1,1,U_EQ) = t(i,1,1,1,U_EQ)/3.+rktp1(i,1,1,1,U_EQ)*2./3.
     $                 - dtf*rkRHS(i,1,1,1,U_EQ)*2./3.
       t(i,1,1,1,V_EQ) = t(i,1,1,1,V_EQ)/3.+rktp1(i,1,1,1,V_EQ)*2./3.
     $                 - dtf*rkRHS(i,1,1,1,V_EQ)*2./3.
       if(if3d) then
         t(i,1,1,1,W_EQ) = t(i,1,1,1,W_EQ)/3.+rktp1(i,1,1,1,W_EQ)*2./3.
     $                   - dtf*rkRHS(i,1,1,1,W_EQ)*2./3.
       endif
       t(i,1,1,1,E_EQ) = t(i,1,1,1,E_EQ)/3.+rktp1(i,1,1,1,E_EQ)*2./3.
     $                -dtf*rkRHS(i,1,1,1,E_EQ)*2./3. 
       do b=1,nspec
         t(i,1,1,1,E_EQ+b) = t(i,1,1,1,E_EQ+b)/3.
     $                     + rktp1(i,1,1,1,E_EQ+b)*2./3.
     $                     - dtf*rkRHS(i,1,1,1,E_EQ+b)*2./3.
       enddo
      enddo

      if(ifpplim) call pp_limiter(t)

ccc update variables
      call eval_dependent_vars(t)

      return
      end
c---------------------------------------------------------------------
      subroutine integrate_chem(tstart,tend)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real k1  (lx1,ly1,lz1,lelt,nspec)
      real k2  (lx1,ly1,lz1,lelt,nspec)
      
      real ysum(lx1,ly1,lz1,lelt)

      real tstart
      real tend
      real time_chem
      real dt_chem
      integer c
      integer n

      n  = nx1*ny1*nz1*nelt

      time_chem = tstart
      dt_chem = uparam(1)

      c = 0
      do while (time_chem.lt.tend)
        c = c + 1
        if(tend-time_chem.lt.dt_chem.or.
     &     tend-time_chem-dt_chem.lt.1.e-12) then
           dt_chem = tend-time_chem
        endif
        time_chem = time_chem + dt_chem

        call chemRHS(k1)
        do i=1,n
          do b=1,nspec
            Yi(i,1,1,1,b) = Yi(i,1,1,1,b)
     &                 +0.5*dt_chem*k1(i,1,1,1,b)/vtrans(i,1,1,1,1)
          enddo
        enddo

        call chemRHS(k2)        

        do i=1,n
          do b=1,nspec
            t(i,1,1,1,E_EQ+b)=t(i,1,1,1,E_EQ+b)
     &            +0.5*dt_chem*(k1(i,1,1,1,b)+k2(i,1,1,1,b))
          enddo
        enddo

      enddo

      call chemRHS(wdot)

      return
      end
c---------------------------------------------------------------------
      subroutine chemRHS(omg)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'   

      real rkU(lx1,ly1,lz1,lelt,N_EQ)
      real omg(lx1,ly1,lz1,lelt,nspec)
      real yy(nspec)
      real kinE,totE,intE
      real wdot1(nspec)
      real wdot2(nspec)
      real wdot_sum

      integer i,n,b

      n = nx1*ny1*nz1*nelt

      do i=1,n
        do b=1,nspec
          yy(b) = Yi(i,1,1,1,b)
        enddo
        call setmassfractions(yy)
        kinE = 0.5*spd(i,1,1,1)**2                  ! [J/kg]
        totE = t(i,1,1,1,E_EQ)/vtrans(i,1,1,1,1)    ! [J/kg]
        intE = totE-kinE                            ! [J/kg]
        call setstate_uv(intE, 1./vtrans(i,1,1,1,1))
        call getnetproductionrates(wdot1)
        wdot_sum = 0.
        do b=1,nspec
          wdot2(b) = mw(b)*wdot1(b)*1.e3
          omg(i,1,1,1,b) = wdot2(b)
          wdot_sum = wdot_sum + omg(i,1,1,1,b)
        enddo
        omg(i,1,1,1,nspec) = omg(i,1,1,1,nspec) - wdot_sum
      enddo

      return
      end
c---------------------------------------------------------------------    
      subroutine pp_limiter(QQ)

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
      
      real QQ(lx1,ly1,lz1,lelt,N_EQ)
      real rho1,rho1_min,rho1_avg
      real rhoY_min,rhoY(lx1*ly1*lz1,nspec)  
      real q(N_EQ),qavg(N_EQ)
      real rU_min,rU_shift
      real rS_min,rS_shift
      real eps2,thetak
      real vlsc2
      integer n,nxyz
      integer i,e,b

      nxyz = lx1*ly1*lz1
      n = nx1*ny1*nz1*nelt 

      eps2 = 1.e-10


ccc clip values below machine precision
      do e=1,nelv
        do b=1,nspec
          do i=1,nxyz
            if(abs(QQ(i,1,1,e,E_EQ+b)).lt.1.e-20) then
              QQ(i,1,1,e,E_EQ+b) = 0.
            endif  
          enddo
        enddo
      enddo



ccc rho
      do e=1,nelv

        ! rho, rho_min
        rho1_min = 1.
        do i=1,nxyz
          rho1 = 0.
          do b=1,nspec
            rho1 = rho1 + QQ(i,1,1,e,E_EQ+b)
          enddo
          if(rho1.lt.rho1_min) rho1_min = rho1 ! evaluate minimum density in element
        enddo 

        ! modify species
        if(rho1_min.lt.eps2) then
           rho1_avg = 0.
           do b=1,nspec
             qavg(E_EQ+b) = 
     &                 vlsc2(QQ(1,1,1,e,b),bm1(1,1,1,e),nxyz)/volel(e)
             rho1_avg = rho1_avg + qavg(E_EQ+b)
           enddo
           thetak = (rho1_avg-eps2)/(rho1_avg-rho1_min)
           do b=1,nspec
             do i=1,nxyz
               QQ(i,1,1,e,E_EQ+b) = qavg(E_EQ+b)
     &                   + thetak*(QQ(i,1,1,e,E_EQ+b)-qavg(E_EQ+b))
             enddo 
           enddo
        endif  

      enddo



ccc rho*Y_i
      do e=1,nelv    

        do b=1,nspec

          ! (rho*Y)_min
          rhoY_min = 1.
          do i=1,nxyz
            if(QQ(i,1,1,e,E_EQ+b).lt.rhoY_min) then
              rhoY_min = QQ(i,1,1,e,E_EQ+b)
            endif    
          enddo
          
          ! modify conservative variables
          if(rhoY_min.lt.0.) then
             qavg(E_EQ+b) =
     &            vlsc2(QQ(1,1,1,e,E_EQ+b),bm1(1,1,1,e),nxyz)/volel(e)        
             thetak = qavg(E_EQ+b)/(qavg(E_EQ+b)-rhoY_min)
             do i=1,nxyz
               QQ(i,1,1,e,E_EQ+b) = qavg(E_EQ+b) 
     &                   + thetak*(QQ(i,1,1,e,E_EQ+b)-qavg(E_EQ+b))
             enddo
          endif

        enddo

      enddo



ccc rho*u_shift
      do e=1,nelv 
          
        ! rho*u_shift for each GLL point
        rU_min = 1.
        do i=1,nxyz
          do b=1,N_EQ
            q(b) = QQ(i,1,1,e,b)
          enddo
          call eval_rU_shift(q,rU_shift)
          if(rU_shift.lt.rU_min) rU_min = rU_shift ! min rU_shift in the element
        enddo

        ! modify conservative variables
        if(rU_min.lt.eps2) then
           do b=1,N_EQ
             qavg(b) = vlsc2(QQ(1,1,1,e,b),bm1(1,1,1,e),nxyz)/volel(e) ! average of conserved variables in the element
           enddo
           call eval_rU_shift(qavg,rU_shift) ! rU_shift based on avg conserved variables
           thetak = (rU_shift-eps2)/(rU_shift-rU_min)          
           do b=1,N_EQ
             do i=1,nxyz
               QQ(i,1,1,e,b) = qavg(b)+thetak*(QQ(i,1,1,e,b)-qavg(b)) ! rescale solution
             enddo
           enddo
        endif

      enddo



ccc rho*s_shift
      do e=1,nelv

        ! rho*s_shift for each GLL point
        rS_min = 1.
        do i=1,nxyz
          do b=1,N_EQ
            q(b) = QQ(i,1,1,e,b)
          enddo
          call eval_rS_shift(q,rS_shift)
          if(rS_shift.lt.rS_min) rS_min = rS_shift ! min rS_shift in the element
        enddo

        ! modify conservative variables
        if(rS_min.lt.0.) then
           do b=1,N_EQ
             qavg(b) = vlsc2(QQ(1,1,1,e,b),bm1(1,1,1,e),nxyz)/volel(e) ! average of conserved variables in the element
           enddo
           call eval_rS_shift(qavg,rS_shift) ! rS_shift based on avg conserved variables
           thetak = rS_shift/(rS_shift-rS_min)
           thetak = min(max(thetak,0.),1.)
           do b=1,N_EQ
             do i=1,nxyz
               QQ(i,1,1,e,b) = qavg(b)+thetak*(QQ(i,1,1,e,b)-qavg(b)) ! rescale solution
             enddo
           enddo
        endif

      enddo


      return
      end
c---------------------------------------------------------------------
      subroutine eval_rU_shift(q,rU_shift)

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real q(N_EQ)
      real yy(nspec),rU_shift,ushift1
      real spd1,kinE,totE,intE,rho1
      integer b,j

ccc density
      rho1 = 0.
      do b=1,nspec
        rho1 = rho1 + q(E_EQ+b)
      enddo

ccc rhoY, yy
      do b=1,nspec
        yy(b) = q(E_EQ+b)/rho1
      enddo

ccc velocity/speed
      spd1 = 0.
      do b=1,ldim
        spd1 = spd1 + (q(U_EQ-1+b)*q(U_EQ-1+b))/rho1**2.
      enddo
      spd1 = sqrt(spd1)

ccc internal energy
      kinE = 0.5*spd1**2.
      totE = q(E_EQ)/rho1
      intE = totE-kinE

ccc u shift
      ushift1 = intE
      do b=1,nspec
        j = (b-1)*(npa+1)
        ushift1 = ushift1 - yy(b)*bi(j+1)
      enddo
      
ccc rho*ushift
      rU_shift = ushift1*rho1

      return
      end   
c---------------------------------------------------------------------
      subroutine eval_rS_shift(q,rS_shift)

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real new_max,new_tol,dTT,temp0,starting_temp
      real u0,ustart,dudT
      integer f,c

      real q(N_EQ)
      real yy(nspec),rS_shift,ssi(nspec)
      real spd1,kinE,totE,intE,rho1
      real ppi
      integer b,j

ccc density
      rho1 = 0.
      do b=1,nspec
        rho1 = rho1 + q(E_EQ+b)
      enddo

ccc rhoY, yy
      do b=1,nspec
        yy(b) = q(E_EQ+b)/rho1
      enddo

ccc velocity/speed
      spd1 = 0.
      do b=1,ldim
        spd1 = spd1 + (q(U_EQ-1+b)*q(U_EQ-1+b))/rho1**2.
      enddo
      spd1 = sqrt(spd1)

ccc internal energy
      kinE = 0.5*spd1**2.
      totE = q(E_EQ)/rho1
      intE = totE-kinE

ccc Newton's method for T
      new_max = 1000
      new_tol = 1.e-8
      dTT = 1.
      f = 0
      temp0 = 500.
      starting_temp = temp0
      do while(abs(dTT).gt.new_tol.and.f.le.new_max)
        f = f + 1
        u0 = -intE
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

ccc rS_shift - same as compvar and then subtracting rho*s_0
      rS_shift = 0.
      do b=1,nspec
        ppi = rho1*max(1.e-16,yy(b))*R_u*temp0/mw(b) ! partial pressure
        j = (b-1)*npa
        ssi(b) = s0(b)+aa(j+1)*log(temp0/298.15)
     &         - (R_u/mw(b))*log(ppi/101325.)
        do c=2,npa
          ssi(b) = ssi(b) 
     &     + aa(j+c)*(temp0**(c-1)-298.15**(c-1))/(c-1)
        enddo
        rS_shift = rS_shift + yy(b)*ssi(b) ! s = sum(Y_i*s_i)
      enddo
      rS_shift = rho1*(rS_shift-s0min) ! xi = rho*(s - sigma)

      return
      end   
c---------------------------------------------------------------------

