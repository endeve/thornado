#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_makefinemask_thornado
      ( iMultiFab*& Mask,
        const BoxArray& CrseBA,
        const DistributionMapping& CrseDM,
        const BoxArray& FineBA,
        int iCoarse, int iFine )
    {
        Mask
          = new iMultiFab( makeFineMask( CrseBA, CrseDM, FineBA,
                                         amrex::IntVect(2), iCoarse, iFine ) );
    }
}
