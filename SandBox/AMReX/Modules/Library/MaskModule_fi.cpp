#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_createfinemask_thornado
      ( iMultiFab*& FineMask,
        const BoxArray& CrseBA,
        const DistributionMapping& CrseDM,
        const BoxArray& FineBA,
        const int iCoarse,
        const int iFine,
        const int swX,
        const Geometry * geom )
    {
        if( &CrseBA == NULL )
        {
            std::cout << "CrseBA IS NULL!\n";
            amrex::Abort();
        }
        if( &FineBA == NULL )
        {
            std::cout << "FineBA IS NULL!\n";
            amrex::Abort();
        }
        FineMask
          = new iMultiFab
                  ( makeFineMask
                    ( CrseBA, CrseDM, amrex::IntVect(swX), FineBA,
                      amrex::IntVect(2), geom -> periodicity(),
                      iCoarse, iFine ) );
    }

    void amrex_fi_createpointmask_thornado
      ( iMultiFab*& PointMask,
        const Geometry * geom,
        const int iCovered,
        const int iNotCovered,
        const int iPhysicalBoundary,
        const int iInterior )

    {
        PointMask -> BuildMask( geom -> Domain(),
                                geom -> periodicity(),
                                iCovered,
                                iNotCovered,
                                iPhysicalBoundary,
                                iInterior );
    }
}
