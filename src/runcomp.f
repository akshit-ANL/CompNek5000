c---------------------------------------------------------------------
      subroutine runcompdg
      
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'COMPRESS'
      
      real mass0,mass1,mchange

      integer icalld
      save icalld
      data icalld /0/
      
ccc error handling
      if(npa.gt.5.and.ifVolEC) then
        call exitti('ABORT: recompile with npa<=5. npa=$',nzq)
      endif  
       
ccc setup
      if(icalld.eq.0) then
         ifflow = .false.
         ifheat = .false.
         do ifld=1,N_EQ
            ifdgfld(ifld) = .true.
         enddo
         call dg_setup
         do ifld=1,N_EQ
            if(ifdgfld(ifld)) then
               call dg_setup2(tmask(1,1,1,1,ifld))
            endif
         enddo
         icalld = 1
         call eval_enhe
         call printFlag
         return
      endif

      if(ifforces.and.ifsplitchem) then

ccc flow substep 1
c         if(nid.eq.0) write(*,*)'Starting flow step 1'
c         call check_mass_min
c         call check_mass_int(mass0)   
         call compress_cfl
         call SSPRK3(dt/2.)
c         call check_mass_min
c         call check_mass_int(mass1)
c         mchange = (mass1-mass0)/mass0
c         if(nid.eq.0) write(*,*)'Flow step 1 % mass change =',mchange
c         if(nid.eq.0) write(*,*)'Completed flow step 1'

ccc chemistry substep
c         if(nid.eq.0) write(*,*)'Starting chem step'
c         call check_mass_min
c         call check_mass_int(mass0)
         call integrate_chem(time,time+dt)
c         call check_mass_min
c         call check_mass_int(mass1)
c         mchange = (mass1-mass0)/mass0
c         if(nid.eq.0) write(*,*)'Chem step % mass change =',mchange
c         if(nid.eq.0) write(*,*)'Completed chem step'

ccc flow substep 2
c         if(nid.eq.0) write(*,*)'Starting flow step 2'
c         call check_mass_min
c         call check_mass_int(mass0)
         call compress_cfl
         call SSPRK3(dt/2.)
c         call check_mass_min
c         call check_mass_int(mass1)
c         mchange = (mass1-mass0)/mass0
c         if(nid.eq.0) write(*,*)'Flow step 2 % mass change =',mchange
c         if(nid.eq.0) write(*,*)'Completed flow step 2'

      else

ccc flow/chem step
         call compress_cfl
         call SSPRK3(dt)

      endif


      return
      end
c---------------------------------------------------------------------
      subroutine check_mass_int(m0)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real m0
      integer n

      n = nx1*ny1*nz1*nelt

      m0 = glsc2(vtrans(1,1,1,1,1),bm1,n)

      return
      end
c----------------------------------------------------------------------- 
      subroutine check_mass_min

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      integer n
      integer b
      real S_min,S_max,glmin,glmax

      integer i

      n = nx1*ny1*nz1*nelt

ccc Y_i
      do b=1,nspec
        S_min = glmin(Yi(1,1,1,1,b),n)
        S_max = glmax(Yi(1,1,1,1,b),n)
        if(nid.eq.0) then
           if(b.eq.7) then
           write(*,3) 'S',b,' min/max = [',S_min,',',S_max,']'
  3        format(A1,I2.2,A12,E14.8,A1,E14.8,A1)
           endif
        endif
        if(b.eq.7) then
          do i=1,n
            if(Yi(i,1,1,1,b).eq.S_min.and.
     &         Yi(i,1,1,1,b).lt.0.) then
              write(*,*)'Yi=',Yi(i,1,1,1,b)
     &                 ,'x=',xm1(i,1,1,1)
     &                 ,'y=',ym1(i,1,1,1)
     &                 ,'z=',zm1(i,1,1,1)     
            endif
          enddo
        endif
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine printFlag()

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      if(nid.eq.0) then
        write(*,*)'ifviscous    =',ifviscous
        write(*,*)'ifforces     =',ifforces
        write(*,*)'ifavisc      =',ifavisc
        write(*,*)'ifsplitchem  =',ifsplitchem
        write(*,*)'ifclip       =',ifclip
        write(*,*)'if1d         =',if1d
        write(*,*)'ifpplim      =',ifpplim
        write(*,*)'ifIntLLFfull =',ifIntLLFfull
        write(*,*)'ifVolEC      =',ifVolEC
        if(ifsplitchem) then
          write(*,*)'dt_chem      =',uparam(1)
        endif
        if(ifavisc) then
          write(*,*)'cstar        =',uparam(2)
          write(*,*)'cmax         =',uparam(3)
          write(*,*)'P            =',uparam(4)
        endif
      endif

      return
      end
c----------------------------------------------------------------------- 

