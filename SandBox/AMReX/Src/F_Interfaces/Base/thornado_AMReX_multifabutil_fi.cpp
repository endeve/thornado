#include <thornado_AMReX_MultiFabUtil.H>

using namespace amrex;

extern "C"
{
  void amrex_fi_average_down_dg_conservative
    ( const MultiFab * FineMF  ,       MultiFab * CrseMF,
      const MultiFab * FineMF_G, const MultiFab * CrseMF_G,
      int nComp, int RefRatio, int nDOFX, int nFine,
      void * vpFineToCoarseProjectionMatrix )
  {
    auto * pFineToCoarseProjectionMatrix
           = reinterpret_cast<Real*>(vpFineToCoarseProjectionMatrix);
    Array4<Real const> FineToCoarseProjectionMatrix
                         ( pFineToCoarseProjectionMatrix,
                           {0,0,0}, {1,nFine,nDOFX}, nDOFX );

    amrex::average_down_dg_conservative
      ( *FineMF, *CrseMF, *FineMF_G, *CrseMF_G,
        nComp, RefRatio, nDOFX, FineToCoarseProjectionMatrix );
  }
  void amrex_fi_average_down_dg_pointwise
    ( const MultiFab * FineMF, MultiFab * CrseMF,
      int nComp, int RefRatio, int nDOFX, int nFine,
      void * vpFineToCoarseProjectionMatrix )
  {
    auto * pFineToCoarseProjectionMatrix
           = reinterpret_cast<Real*>(vpFineToCoarseProjectionMatrix);
    Array4<Real const> FineToCoarseProjectionMatrix
                         ( pFineToCoarseProjectionMatrix,
                           {0,0,0}, {1,nFine,nDOFX}, nDOFX );

    amrex::average_down_dg_pointwise
      ( *FineMF, *CrseMF, nComp, RefRatio, nDOFX,
        FineToCoarseProjectionMatrix );
  }
  void amrex_fi_average_down_cg
    ( const MultiFab * FineMF, MultiFab * CrseMF,
      int nComp, int RefRatio, int nDOFX, int nFine,
      void * vpG2L, void * vpL2G, void * vpF2C )
  {
    auto * pG2L
           = reinterpret_cast<Real*>(vpG2L);
    Array4<Real> G2L
                   ( pG2L, {0,0,0}, {1,1,nDOFX}, nDOFX );

    auto * pL2G
           = reinterpret_cast<Real*>(vpL2G);
    Array4<Real> L2G
                   ( pL2G, {0,0,0}, {1,1,nDOFX}, nDOFX );

    auto * pF2C
           = reinterpret_cast<Real*>(vpF2C);
    Array4<Real> F2C
                   ( pF2C, {0,0,0}, {1,nDOFX,nFine}, nDOFX );

    amrex::average_down_cg
      ( *FineMF, *CrseMF, nComp, RefRatio, nDOFX,
        nFine, G2L, L2G, F2C );
  }
}