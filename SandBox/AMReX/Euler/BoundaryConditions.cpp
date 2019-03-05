#include <AMReX_BCUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

extern "C"
{
  // function name MUST be all lowercase!!!
  void mf_euler_applyboundaryconditions
    ( MultiFab* MF_inout, const Geometry* GEOM_in, const Vector<BCRec>& BC )
  {
          MultiFab& MF  = *MF_inout;
    const Geometry GEOM = *GEOM_in;

    //std::cout << "Call FillDomainBoundary" << std::endl;
    FillDomainBoundary( MF, GEOM, BC );
  }
}
