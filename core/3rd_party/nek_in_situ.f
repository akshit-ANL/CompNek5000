c-----------------------------------------------------------------------
      subroutine in_situ_init()
#ifdef VISIT
      call visit_init()
#elif CATALYST2
      include 'PARALLEL'
      call nek_catalyst_initialize(nid)
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_check()
#ifdef VISIT
      call visit_check()
#elif CATALYST2
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      integer ndim,has_temperature

      ndim = 2
      if (if3d) ndim = 3
      has_temperature = 0
      if (ifheat) has_temperature = 1
      call nek_catalyst_execute(xm1,ym1,zm1,vx,vy,vz,t,
     &     lx1,ly1,lz1,nelt,ndim,has_temperature,istep,time)
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_end()
#ifdef VISIT
      call visit_end()
#elif CATALYST2
      call nek_catalyst_finalize()
#endif
      end
c-----------------------------------------------------------------------

