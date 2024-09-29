#ifndef THORNADO_AMREX_MultiFabUtil_H_
#define THORNADO_AMREX_MultiFabUtil_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFabUtil_C.H>

#include <AMReX_MultiFabUtilI.H>

namespace amrex
{
  //! Average fine DG-based MultiFab onto crse DG-based MultiFab,
  //! maintaining conservation of multifabs S_fine and S_crse.
  //! This routine DOES NOT assume that the crse BoxArray is
  //! a coarsened version of the fine BoxArray.
  void average_down_dg_conservative
         ( const MultiFab & FineMF  ,       MultiFab & CrseMF,
           const MultiFab & FineMF_G, const MultiFab & CrseMF_G,
           int nComp, const IntVect & RefRatio, int nDOFX,
           Array4<Real const> FineToCoarseProjectionMatrix );
  void average_down_dg_conservative
         ( const MultiFab & FineMF  ,       MultiFab & CrseMF,
           const MultiFab & FineMF_G, const MultiFab & CrseMF_G,
           int nComp, int RefRatio, int nDOFX,
           Array4<Real const> FineToCoarseProjectionMatrix );

  //! Average fine DG-based MultiFab onto crse DG-based MultiFab
  //! using an L2-projection.
  //! This routine DOES NOT assume that the crse BoxArray is
  //! a coarsened version of the fine BoxArray.
  void average_down_dg_pointwise
         ( const MultiFab & FineMF, MultiFab & CrseMF,
           int nComp, const IntVect & RefRatio, int nDOFX,
           Array4<Real const> FineToCoarseProjectionMatrix );
  void average_down_dg_pointwise
         ( const MultiFab & FineMF, MultiFab & CrseMF,
           int nComp, int RefRatio, int nDOFX,
           Array4<Real const> FineToCoarseProjectionMatrix );

  //! Average fine DG-based MultiFab onto crse DG-based MultiFab,
  //! enforcing continuity across element interfaces.
  //! This routine DOES NOT assume that the crse BoxArray is
  //! a coarsened version of the fine BoxArray.
  void average_down_cg
         ( const MultiFab & FineMF, MultiFab & CrseMF,
           int nComp, const IntVect& RefRatio, int nDOFX, int nFine,
           Array4<Real> G2L, Array4<Real> L2G, Array4<Real> F2C );
  void average_down_cg
         ( const MultiFab & FineMF, MultiFab & CrseMF,
           int nComp, int RefRatio, int nDOFX, int nFine,
           Array4<Real> G2L, Array4<Real> L2G, Array4<Real> F2C );
}
#endif
