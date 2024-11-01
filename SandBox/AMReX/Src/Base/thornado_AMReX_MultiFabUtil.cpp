#include <thornado_AMReX_MultiFabUtil.H>
#if (AMREX_SPACEDIM==1)
#include <thornado_AMReX_MultiFabUtil_1D_C.H>
#elif (AMREX_SPACEDIM==2)
#include <thornado_AMReX_MultiFabUtil_2D_C.H>
#else
#include <thornado_AMReX_MultiFabUtil_3D_C.H>
#endif
#include <AMReX_Random.H>
#include <sstream>
#include <iostream>

namespace amrex
{

  // Average fine nodal DG-based MultiFab onto crse nodal DG-based MultiFab.
  // We do NOT assume that the coarse layout is a coarsened version
  // of the fine layout.
  void average_down_dg_conservative
    ( const MultiFab & FineMF,         MultiFab & CrseMF,
      const MultiFab & FineMF_G, const MultiFab & CrseMF_G,
      int nComp, int RefRatio, int nDOFX,
      Array4<Real const> FineToCoarseProjectionMatrix )
  {
    average_down_dg_conservative
      ( FineMF, CrseMF, FineMF_G, CrseMF_G, nComp,
        RefRatio * IntVect::TheUnitVector(), nDOFX,
        FineToCoarseProjectionMatrix );
  }

  void average_down_dg_conservative
    ( const MultiFab & FineMF  ,       MultiFab & CrseMF,
      const MultiFab & FineMF_G, const MultiFab & CrseMF_G,
      int nComp, const IntVect & RefRatio, int nDOFX,
      Array4<Real const> FineToCoarseProjectionMatrix )
  {

    BL_PROFILE("amrex::average_down_dg_conservative");

    if ( FineMF.is_nodal() || CrseMF.is_nodal() )
    {
      amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
    }

    AMREX_ASSERT( CrseMF.nComp() == FineMF.nComp() );

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    const BoxArray            & FineBA         = FineMF.boxArray();
    const DistributionMapping & FineDM         = FineMF.DistributionMap();
    BoxArray                    crse_S_fine_BA = FineBA;
    crse_S_fine_BA.coarsen( RefRatio );

    MultiFab crse_S_fine
      ( crse_S_fine_BA, FineDM, nComp, 0, MFInfo(), FArrayBoxFactory() );

    // Create MultiFab for SqrtGm on coarse level that uses same DM
    // as dst MultiFab
    MultiFab crse_G_fine
      ( crse_S_fine_BA, FineDM, nDOFX, 0, MFInfo(), FArrayBoxFactory() );
    crse_G_fine.ParallelCopy( CrseMF_G, 0, 0, nDOFX, 0, 0 );

#ifdef AMREX_USE_GPU

    if ( Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate() ) {
      auto const & crsema  = crse_S_fine.arrays();
      auto const & finema  = FineMF.const_arrays();
      auto const & finemaG = FineMF_G.const_arrays();
      auto const & crsemaG = CrseMF_G.const_arrays();
      ParallelFor(crse_S_fine, IntVect(0), nComp,
      [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
      {
          amrex_avgdown_dg_conservative
            ( i, j, k, n,
              crsema [box_no], finema [box_no],
              crsemaG[box_no], finemaG[box_no], RefRatio,
              nDOFX, FineToCoarseProjectionMatrix );
      });
      Gpu::streamSynchronize();
    } else

#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for( MFIter mfi( crse_S_fine, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
      {
        //  NOTE: The tilebox is defined at the coarse level.
        const Box& bx = mfi.tilebox();
        Array4<Real>       const & crsearr  = crse_S_fine.array (mfi);
        Array4<Real const> const & finearr  = FineMF.const_array(mfi);
        Array4<Real const> const & finearrG = FineMF_G.const_array(mfi);
        Array4<Real const> const & crsearrG = crse_G_fine.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
          amrex_avgdown_dg_conservative
            ( i, j, k, nComp,
              crsearr, finearr, crsearrG, finearrG,
              RefRatio, nDOFX, FineToCoarseProjectionMatrix );
        });
      }
    }

    CrseMF.ParallelCopy( crse_S_fine, 0, 0, nComp );
  } // end void average_down_dg_conservative

  // Average fine nodal DG-based MultiFab onto crse nodal DG-based MultiFab.
  // We do NOT assume that the coarse layout is a coarsened version
  // of the fine layout.
  void average_down_dg_pointwise
    ( const MultiFab & FineMF, MultiFab & CrseMF,
      int nComp, int RefRatio, int nDOFX,
      Array4<Real const> FineToCoarseProjectionMatrix )
  {
    average_down_dg_pointwise
      ( FineMF, CrseMF, nComp,
        RefRatio*IntVect::TheUnitVector(), nDOFX,
        FineToCoarseProjectionMatrix );
  }

  void average_down_dg_pointwise
         ( const MultiFab & FineMF, MultiFab & CrseMF,
           int nComp, const IntVect & RefRatio, int nDOFX,
           Array4<Real const> FineToCoarseProjectionMatrix )
  {

    BL_PROFILE("amrex::average_down_dg_pointwise");

    if ( FineMF.is_nodal() || CrseMF.is_nodal() )
    {
      amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
    }

    AMREX_ASSERT( CrseMF.nComp() == FineMF.nComp() );

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    const BoxArray            & FineBA         = FineMF.boxArray();
    const DistributionMapping & FineDM         = FineMF.DistributionMap();
    BoxArray                    crse_S_fine_BA = FineBA;
    crse_S_fine_BA.coarsen( RefRatio );

    MultiFab crse_S_fine
      ( crse_S_fine_BA, FineDM, nComp, 0, MFInfo(), FArrayBoxFactory() );

#ifdef AMREX_USE_GPU
    if ( Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate() ) {
        auto const& crsema = crse_S_fine.arrays();
        auto const& finema = FineMF.const_arrays();
        ParallelFor(crse_S_fine, IntVect(0), nComp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            amrex_avgdown_dg_pointwise
              ( i, j, k, n, crsema[box_no], finema[box_no], RefRatio,
                nDOFX, FineToCoarseProjectionMatrix );
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for ( MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
          //  NOTE: The tilebox is defined at the coarse level.
          const Box& bx = mfi.tilebox();
          Array4<Real> const& crsearr = crse_S_fine.array(mfi);
          Array4<Real const> const& finearr = FineMF.const_array(mfi);
          AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
          {
              amrex_avgdown_dg_pointwise
                ( i, j, k, nComp, crsearr, finearr, RefRatio,
                  nDOFX, FineToCoarseProjectionMatrix );
          });
      }
    }

    CrseMF.ParallelCopy( crse_S_fine, 0, 0, nComp );
  } // end void average_down_dg_pointwise

  // Average fine nodal DG-based MultiFab onto crse nodal DG-based MultiFab.
  // Enforce continuity across interfaces
  // We do NOT assume that the coarse layout is a coarsened version
  // of the fine layout.
  void average_down_cg
    ( const MultiFab & FineMF, MultiFab & CrseMF,
      int nComp, int RefRatio, int nDOFX, int nFine,
      Array4<Real> G2L, Array4<Real> L2G, Array4<Real> F2C )
  {
    average_down_cg
      ( FineMF, CrseMF, nComp,
        RefRatio * IntVect::TheUnitVector(), nDOFX,
        nFine, G2L, L2G, F2C );
  }

  void average_down_cg
    ( const MultiFab & FineMF, MultiFab & CrseMF,
      int nComp, const IntVect & RefRatio, int nDOFX,
      int nFine, Array4<Real> G2L, Array4<Real> L2G, Array4<Real> F2C )
  {
    BL_PROFILE("amrex::average_down_cg");

    if ( FineMF.is_nodal() || CrseMF.is_nodal() )
    {
        amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
    }

    AMREX_ASSERT( CrseMF.nComp() == FineMF.nComp() );

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    const BoxArray            & FineBA         = FineMF.boxArray();
    const DistributionMapping & FineDM         = FineMF.DistributionMap();
    BoxArray                    crse_S_fine_BA = FineBA;
    crse_S_fine_BA.coarsen( RefRatio );

    MultiFab crse_S_fine
      ( crse_S_fine_BA, FineDM, nComp, 0, MFInfo(), FArrayBoxFactory() );

#ifdef AMREX_USE_GPU

    if ( Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate() ) {
        auto const& crsema = crse_S_fine.arrays();
        auto const& finema = FineMF.const_arrays();
        ParallelFor(crse_S_fine, IntVect(0), nComp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            amrex_avgdown_cg
              ( i, j, k, n, crsema[box_no], finema[box_no], RefRatio,
                nDOFX, nFine, G2L, L2G, F2C );
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for ( MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
        //  NOTE: The tilebox is defined at the coarse level.
        const Box& bx = mfi.tilebox();
        Array4<Real> const& crsearr = crse_S_fine.array(mfi);
        Array4<Real const> const& finearr = FineMF.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            amrex_avgdown_cg
              ( i, j, k, nComp, crsearr, finearr, RefRatio,
                nDOFX, nFine, G2L, L2G, F2C );
        });
      }
    }

        CrseMF.ParallelCopy( crse_S_fine, 0, 0, nComp );
  } // end void average_down_cg

}
