c---------------------------------------------------------------------
c     Link-time stand-ins for the Cantera Fortran interface (core/ftnlib.cpp).
c     Used by Cantera-free cases (single-species perfect gas, no chemistry,
c     no transport). Any call aborts, so a case that needs Cantera fails
c     loudly instead of silently computing garbage.
c---------------------------------------------------------------------
      subroutine cantera_stub_abort(name)
      character*(*) name
      write(6,*) 'ABORT: Cantera call in Cantera-free build: ',name
      call exitt
      return
      end
c---------------------------------------------------------------------
      subroutine setmassfractions(yy)
      real yy(1)
      call cantera_stub_abort('setmassfractions')
      return
      end
c---------------------------------------------------------------------
      subroutine setstate_uv(u,v)
      real u,v
      call cantera_stub_abort('setstate_uv')
      return
      end
c---------------------------------------------------------------------
      subroutine getnetproductionrates(w)
      real w(1)
      call cantera_stub_abort('getnetproductionrates')
      return
      end
c---------------------------------------------------------------------
      subroutine getmixdiffcoeffs(d)
      real d(1)
      call cantera_stub_abort('getmixdiffcoeffs')
      return
      end
c---------------------------------------------------------------------
      real function viscosity()
      viscosity = 0.
      call cantera_stub_abort('viscosity')
      return
      end
c---------------------------------------------------------------------
      real function thermalconductivity()
      thermalconductivity = 0.
      call cantera_stub_abort('thermalconductivity')
      return
      end
c---------------------------------------------------------------------
