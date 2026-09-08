c---------------------------------------------------------------------      
      subroutine compf(rkRHS)
      
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'

      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)
      real force(N_EQ)
      integer ntot,i
      
      ntot = nx1*ny1*nz1*nelt

      do i=1,ntot
        do b=1,nspec
          rkRHS(i,1,1,1,E_EQ+b) = rkRHS(i,1,1,1,E_EQ+b)
     $                            -wdot(i,1,1,1,b)
        enddo
      enddo
      
      return
      end
c---------------------------------------------------------------------   
