c---------------------------------------------------------------------      
      subroutine getdgrhs(rkRHS,rkU)

      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'   

      real rkU  (lx1,ly1,lz1,lelt,N_EQ)
      real rkRHS(lx1,ly1,lz1,lelt,N_EQ)

      call eval_dependent_vars(rkU)
      call sum_neighbor_vars(rkU)
      call compbc
      call advect(rkRHS)
      if(ifavisc.or.ifviscous) call eval_dudx
      if(ifavisc) call entrovisc(rkRHS,rkU)
      if(ifviscous) call viscous(rkRHS,rkU)
      if(ifforces.and.(.not.ifsplitchem)) call compf(rkRHS)

      return
      end
c---------------------------------------------------------------------    
