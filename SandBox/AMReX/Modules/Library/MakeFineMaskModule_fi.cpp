#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_makefinemask_thornado
      ( iMultiFab*& Mask,
        const BoxArray& CrseBA,
        const DistributionMapping& CrseDM,
        const BoxArray& FineBA,
        const int iCoarse,
        const int iFine,
        const int swX,
        const Geometry * geom )
    {

        Mask
          = new iMultiFab
                  ( makeFineMask
                    ( CrseBA, CrseDM, amrex::IntVect(swX), FineBA,
                      amrex::IntVect(2), geom -> periodicity(),
                      iCoarse, iFine ) );
    }

    void amrex_fi_destroyfinemask_thornado
      ( iMultiFab*& Mask )
    {
        delete Mask;
        Mask = NULL;
    }
}
