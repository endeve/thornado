#ifndef THORNADO_AMREX_FillPatchUtil_I_H_
#define THORNADO_AMREX_FillPatchUtil_I_H_
#include <AMReX_Config.H>
#include <AMReX_FillPatchUtil_I.H>

namespace amrex {

  template <typename MF, typename Interp>
  std::enable_if_t<IsFabArray<MF>::value && !std::is_same_v<Interp,MFInterpolater>>
  FillPatchInterp( MF       & FineMF, MF const & FineMF_G,
                   MF const & CrseMF, MF const & CrseMF_G,
                   int nComp, IntVect const & swX, const Geometry & CrseGeom,
                   Box const & FineBox, const IntVect & RefRatio,
                   Interp * mapper, const Vector<BCRec>& bcs, int bcscomp,
                   int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix )
  {
      BL_PROFILE("FillPatchInterp(Fab)");

      Box const& cdomain = amrex::convert( CrseGeom.Domain(), FineMF.ixType() );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
          Vector<BCRec> bcr(nComp);
          for ( MFIter mfi( FineMF ); mfi.isValid(); ++mfi )
          {
              auto      & CrseFab   = CrseMF  [mfi];
              auto      & CrseFab_G = CrseMF_G[mfi];
              const Box & CrseBox   = CrseFab.box();

              auto      & FineFab   = FineMF  [mfi];
              auto      & FineFab_G = FineMF_G[mfi];
              Box const & dbx       = amrex::grow( mfi.validbox(), swX )
                                        & FineBox;

              amrex::setBC( CrseBox, cdomain, bcscomp, 0, nComp, bcs, bcr );

              mapper -> interpConservative
                          ( CrseFab, CrseFab_G, FineFab, FineFab_G,
                            nComp, dbx, RefRatio,
                            nDOFX, CoarseToFineProjectionMatrix, RunOn::Gpu );
          }
      }
  }

  template <typename MF, typename Interp>
  std::enable_if_t<IsFabArray<MF>::value && !std::is_same_v<Interp,MFInterpolater>>
  FillPatchInterp( MF       & FineMF,
                   MF const & CrseMF,
                   int nComp, IntVect const & swX, const Geometry & CrseGeom,
                   Box const & FineBox, const IntVect & RefRatio,
                   Interp * mapper, const Vector<BCRec>& bcs, int bcscomp,
                   int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix )
  {
      BL_PROFILE("FillPatchInterp(Fab)");

      Box const& cdomain = amrex::convert( CrseGeom.Domain(), FineMF.ixType() );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
          Vector<BCRec> bcr(nComp);
          for ( MFIter mfi( FineMF ); mfi.isValid(); ++mfi )
          {
              auto      & CrseFab = CrseMF  [mfi];
              const Box & CrseBox = CrseFab.box();

              auto      & FineFab = FineMF  [mfi];
              Box const & dbx     = amrex::grow( mfi.validbox(), swX )
                                      & FineBox;

              amrex::setBC( CrseBox, cdomain, bcscomp, 0, nComp, bcs, bcr );

              mapper -> interpPointWise
                          ( CrseFab, FineFab,
                            nComp, dbx, RefRatio,
                            nDOFX, CoarseToFineProjectionMatrix, RunOn::Gpu );
          }
      }
  }

namespace detail {

    template <typename MF, typename BC, typename Interp>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels_doit
      ( MF & mf, MF & mf_G, IntVect const & swX, Real Time,
        const Vector<MF*> & CrseMF, const Vector<MF*> & CrseMF_G,
        const Vector<Real> & CrseTime,
        const Vector<MF*> & FineMF, const Vector<MF*>& FineMF_G,
        const Vector<Real>& FineTime,
        int sComp, int dComp, int nComp,
        const Geometry & CrseGeom, const Geometry & FineGeom,
        BC & CrseBC, int CrseBCcomp,
        BC & FineBC, int FineBCcomp,
        const IntVect & RefRatio,
        Interp * mapper,
        const Vector<BCRec> & bcs, int bcscomp,
        int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix,
        EB2::IndexSpace const* index_space )
    {
        BL_PROFILE("FillPatchTwoLevels");

        if ( swX.max() > 0 || mf.getBDKey() != FineMF[0]->getBDKey() )
        {
            const InterpolaterBoxCoarsener & coarsener
                    = mapper->BoxCoarsener( RefRatio );

            // Test for Face-centered data
            if ( AMREX_D_TERM(  mf.ixType().nodeCentered(0),
                              + mf.ixType().nodeCentered(1),
                              + mf.ixType().nodeCentered(2) ) == 1 )
            {
                if ( !dynamic_cast<Interpolater*>(mapper) ){
                    amrex::Abort("This interpolater has not yet implemented a version for face-based data");
                }

                // Convert to cell-centered MF meta-data for FPInfo.
                MF mf_cc_dummy
                     ( amrex::convert( mf.boxArray(),
                                       IntVect::TheZeroVector() ),
                       mf.DistributionMap(), nComp, swX,
                       MFInfo().SetAlloc( false ) );

                MF fmf_cc_dummy
                     ( amrex::convert( FineMF[0]->boxArray(),
                                       IntVect::TheZeroVector() ),
                       FineMF[0]->DistributionMap(), nComp, swX,
                       MFInfo().SetAlloc( false ) );

                const FabArrayBase::FPinfo & fpc
                         = FabArrayBase::TheFPinfo
                             ( fmf_cc_dummy, mf_cc_dummy, swX, coarsener,
                               FineGeom, CrseGeom, index_space );

                if ( ! fpc.ba_crse_patch.empty() )
                {
                    MF mf_crse_patch
                         = make_mf_crse_patch<MF>
                             ( fpc, nComp, mf.boxArray().ixType() );
                    // Must make sure fine exists under needed coarse faces.
                    // It stores values for the final (interior) interpolation,
                    // which is done from this fine MF that's been partially filled
                    // (with only faces overlying coarse having valid data).
                    MF mf_refined_patch
                         = make_mf_refined_patch<MF>
                             ( fpc, nComp, mf.boxArray().ixType(), RefRatio );
                    auto solve_mask
                           = make_mf_crse_mask<iMultiFab>
                               ( fpc, nComp, mf.boxArray().ixType(), RefRatio );

                    mf_set_domain_bndry( mf_crse_patch, CrseGeom );
                    FillPatchSingleLevel
                      ( mf_crse_patch, Time, CrseMF, CrseTime,
                        sComp, 0, nComp, CrseGeom, CrseBC, CrseBCcomp );

                    mf_set_domain_bndry( mf_refined_patch, FineGeom );
                    FillPatchSingleLevel
                      ( mf_refined_patch, Time, FineMF, FineTime,
                        sComp, 0, nComp, FineGeom, FineBC, FineBCcomp);

                    // Aliased MFs, used to allow CPC caching.
                    MF mf_known
                         ( amrex::coarsen
                             ( FineMF[0]->boxArray(), RefRatio ),
                               FineMF[0]->DistributionMap(),
                               nComp, swX, MFInfo().SetAlloc( false ) );
                    MF mf_solution
                         ( amrex::coarsen
                             ( mf_refined_patch.boxArray(), RefRatio ),
                               mf_refined_patch.DistributionMap(),
                               nComp, 0, MFInfo().SetAlloc( false ) );

                    const FabArrayBase::CPC
                            mask_cpc( mf_solution, IntVect::TheZeroVector(),
                                      mf_known, IntVect::TheZeroVector(),
                                      FineGeom.periodicity() );

                    solve_mask.setVal( 1 );                   // Values to solve.
                    solve_mask.setVal( 0, mask_cpc, 0, 1 );   // Known values.

                    InterpFace
                      ( mapper, mf_crse_patch, 0, mf_refined_patch, 0, nComp,
                        RefRatio, solve_mask, CrseGeom, FineGeom,
                        bcscomp, RunOn::Gpu, bcs);

                    bool aliasing = false;
                    for ( auto const& fmf_a : FineMF ) {
                        aliasing = aliasing || (&mf == fmf_a);
                    }
                    if ( aliasing ) {
                        mf.ParallelCopyToGhost
                          ( mf_refined_patch, 0, dComp, nComp, IntVect{0}, swX );
                    } else {
                        mf.ParallelCopy
                          ( mf_refined_patch, 0, dComp, nComp, IntVect{0}, swX );
                    }
                }
            }
            else
            {
                const FabArrayBase::FPinfo &fpc
                         = FabArrayBase::TheFPinfo
                             ( *FineMF[0], mf, swX, coarsener,
                               FineGeom, CrseGeom, index_space );
                const FabArrayBase::FPinfo &fpcG
                         = FabArrayBase::TheFPinfo
                             ( *FineMF_G[0], mf_G, swX, coarsener,
                               FineGeom, CrseGeom, index_space );

                if ( ! fpc.ba_crse_patch.empty() )
                {
                    MF mfG_crse_patch = make_mf_crse_patch<MF>( fpcG, nDOFX );
                    MF mf_crse_patch  = make_mf_crse_patch<MF>( fpc , nComp );
                    MF mf_fine_patch  = make_mf_fine_patch<MF>( fpc , nComp );

                    mf_set_domain_bndry( mfG_crse_patch, CrseGeom );
                    mf_set_domain_bndry( mf_crse_patch , CrseGeom );

                    FillPatchSingleLevel
                      ( mfG_crse_patch, Time, CrseMF_G, CrseTime, sComp, 0,
                        nDOFX, CrseGeom, CrseBC, CrseBCcomp );
                    FillPatchSingleLevel
                      ( mf_crse_patch, Time, CrseMF, CrseTime, sComp, 0,
                        nComp, CrseGeom, CrseBC, CrseBCcomp );

                    MF mfG_fine_patch
                         ( mf_fine_patch.boxArray(),
                           mf_fine_patch.DistributionMap(), nDOFX, 1 );
                    mfG_fine_patch.setVal( (Real)0.0 );
                    mfG_fine_patch.ParallelCopy( mf_G, 0, 0, nDOFX, swX, swX );

                    FillPatchInterp
                      ( mf_fine_patch, mfG_fine_patch,
                        mf_crse_patch, mfG_crse_patch,
                        nComp, IntVect(0), CrseGeom,
                        amrex::grow( amrex::convert
                                       ( FineGeom.Domain(), mf.ixType() ), swX ),
                        RefRatio, mapper, bcs, bcscomp,
                        nDOFX, CoarseToFineProjectionMatrix );

                    mf.ParallelCopy( mf_fine_patch, 0, dComp, nComp,
                                     IntVect{0}, swX );
                }
            }
        }

        FillPatchSingleLevel
          ( mf, swX, Time, FineMF, FineTime, sComp, dComp, nComp,
            FineGeom, FineBC, FineBCcomp );
    }

    template <typename MF, typename BC, typename Interp>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels_doit
      ( MF & mf, IntVect const & swX, Real Time,
        const Vector<MF*> & CrseMF,
        const Vector<Real> & CrseTime,
        const Vector<MF*> & FineMF,
        const Vector<Real>& FineTime,
        int sComp, int dComp, int nComp,
        const Geometry & CrseGeom, const Geometry & FineGeom,
        BC & CrseBC, int CrseBCcomp,
        BC & FineBC, int FineBCcomp,
        const IntVect & RefRatio,
        Interp * mapper,
        const Vector<BCRec> & bcs, int bcscomp,
        int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix,
        EB2::IndexSpace const* index_space )
    {
        BL_PROFILE("FillPatchTwoLevels");

        if ( swX.max() > 0 || mf.getBDKey() != FineMF[0]->getBDKey() )
        {
            const InterpolaterBoxCoarsener & coarsener
                    = mapper->BoxCoarsener( RefRatio );

            // Test for Face-centered data
            if ( AMREX_D_TERM(  mf.ixType().nodeCentered(0),
                              + mf.ixType().nodeCentered(1),
                              + mf.ixType().nodeCentered(2) ) == 1 )
            {
                if ( !dynamic_cast<Interpolater*>(mapper) ){
                    amrex::Abort("This interpolater has not yet implemented a version for face-based data");
                }

                // Convert to cell-centered MF meta-data for FPInfo.
                MF mf_cc_dummy
                     ( amrex::convert( mf.boxArray(),
                                       IntVect::TheZeroVector() ),
                       mf.DistributionMap(), nComp, swX,
                       MFInfo().SetAlloc( false ) );

                MF fmf_cc_dummy
                     ( amrex::convert( FineMF[0]->boxArray(),
                                       IntVect::TheZeroVector() ),
                       FineMF[0]->DistributionMap(), nComp, swX,
                       MFInfo().SetAlloc( false ) );

                const FabArrayBase::FPinfo & fpc
                         = FabArrayBase::TheFPinfo
                             ( fmf_cc_dummy, mf_cc_dummy, swX, coarsener,
                               FineGeom, CrseGeom, index_space );

                if ( ! fpc.ba_crse_patch.empty() )
                {
                    MF mf_crse_patch
                         = make_mf_crse_patch<MF>
                             ( fpc, nComp, mf.boxArray().ixType() );
                    // Must make sure fine exists under needed coarse faces.
                    // It stores values for the final (interior) interpolation,
                    // which is done from this fine MF that's been partially filled
                    // (with only faces overlying coarse having valid data).
                    MF mf_refined_patch
                         = make_mf_refined_patch<MF>
                             ( fpc, nComp, mf.boxArray().ixType(), RefRatio );
                    auto solve_mask
                           = make_mf_crse_mask<iMultiFab>
                               ( fpc, nComp, mf.boxArray().ixType(), RefRatio );

                    mf_set_domain_bndry( mf_crse_patch, CrseGeom );
                    FillPatchSingleLevel
                      ( mf_crse_patch, Time, CrseMF, CrseTime,
                        sComp, 0, nComp, CrseGeom, CrseBC, CrseBCcomp );

                    mf_set_domain_bndry( mf_refined_patch, FineGeom );
                    FillPatchSingleLevel
                      ( mf_refined_patch, Time, FineMF, FineTime,
                        sComp, 0, nComp, FineGeom, FineBC, FineBCcomp);

                    // Aliased MFs, used to allow CPC caching.
                    MF mf_known
                         ( amrex::coarsen
                             ( FineMF[0]->boxArray(), RefRatio ),
                               FineMF[0]->DistributionMap(),
                               nComp, swX, MFInfo().SetAlloc( false ) );
                    MF mf_solution
                         ( amrex::coarsen
                             ( mf_refined_patch.boxArray(), RefRatio ),
                               mf_refined_patch.DistributionMap(),
                               nComp, 0, MFInfo().SetAlloc( false ) );

                    const FabArrayBase::CPC
                            mask_cpc( mf_solution, IntVect::TheZeroVector(),
                                      mf_known, IntVect::TheZeroVector(),
                                      FineGeom.periodicity() );

                    solve_mask.setVal( 1 );                   // Values to solve.
                    solve_mask.setVal( 0, mask_cpc, 0, 1 );   // Known values.

                    InterpFace
                      ( mapper, mf_crse_patch, 0, mf_refined_patch, 0, nComp,
                        RefRatio, solve_mask, CrseGeom, FineGeom,
                        bcscomp, RunOn::Gpu, bcs);

                    bool aliasing = false;
                    for ( auto const& fmf_a : FineMF ) {
                        aliasing = aliasing || (&mf == fmf_a);
                    }
                    if ( aliasing ) {
                        mf.ParallelCopyToGhost
                          ( mf_refined_patch, 0, dComp, nComp, IntVect{0}, swX );
                    } else {
                        mf.ParallelCopy
                          ( mf_refined_patch, 0, dComp, nComp, IntVect{0}, swX );
                    }
                }
            }
            else
            {
                const FabArrayBase::FPinfo &fpc
                         = FabArrayBase::TheFPinfo
                             ( *FineMF[0], mf, swX, coarsener,
                               FineGeom, CrseGeom, index_space );

                if ( ! fpc.ba_crse_patch.empty() )
                {
                    MF mf_crse_patch = make_mf_crse_patch<MF>( fpc , nComp );
                    MF mf_fine_patch = make_mf_fine_patch<MF>( fpc , nComp );

                    mf_set_domain_bndry( mf_crse_patch , CrseGeom );

                    FillPatchSingleLevel
                      ( mf_crse_patch, Time, CrseMF, CrseTime, sComp, 0,
                        nComp, CrseGeom, CrseBC, CrseBCcomp );

                    FillPatchInterp
                      ( mf_fine_patch,
                        mf_crse_patch,
                        nComp, IntVect(0), CrseGeom,
                        amrex::grow( amrex::convert
                                       ( FineGeom.Domain(), mf.ixType() ), swX ),
                        RefRatio, mapper, bcs, bcscomp,
                        nDOFX, CoarseToFineProjectionMatrix );

                    mf.ParallelCopy( mf_fine_patch, 0, dComp, nComp,
                                     IntVect{0}, swX );
                }
            }
        }

        FillPatchSingleLevel
          ( mf, swX, Time, FineMF, FineTime, sComp, dComp, nComp,
            FineGeom, FineBC, FineBCcomp );
    }

} // namespace detail

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
FillPatchTwoLevels
  ( MF & mf, MF & mf_G, Real Time,
    const Vector<MF*> & CrseMF, const Vector<MF*> & CrseMF_G,
    const Vector<Real> & CrseTime,
    const Vector<MF*> & FineMF, const Vector<MF*> & FineMF_G,
    const Vector<Real> & FineTime,
    int sComp, int dComp, int nComp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix )
{
    EB2::IndexSpace const* index_space = nullptr;

    detail::FillPatchTwoLevels_doit
      ( mf, mf_G, mf.nGrowVect(), Time,
        CrseMF, CrseMF_G, CrseTime,
        FineMF, FineMF_G, FineTime,
        sComp, dComp, nComp, CrseGeom, FineGeom,
        CrseBC, CrseBCcomp, FineBC, FineBCcomp, RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix,
        index_space );
}

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
FillPatchTwoLevels
  ( MF & mf, Real Time,
    const Vector<MF*> & CrseMF,
    const Vector<Real> & CrseTime,
    const Vector<MF*> & FineMF,
    const Vector<Real> & FineTime,
    int sComp, int dComp, int nComp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix )
{
    EB2::IndexSpace const* index_space = nullptr;

    detail::FillPatchTwoLevels_doit
      ( mf, mf.nGrowVect(), Time,
        CrseMF, CrseTime,
        FineMF, FineTime,
        sComp, dComp, nComp, CrseGeom, FineGeom,
        CrseBC, CrseBCcomp, FineBC, FineBCcomp, RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix,
        index_space );
}

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel
  ( MF & mf, MF & mf_G, Real Time,
    const MF & CrseMF, const MF & CrseMF_G,
    int scomp, int dcomp, int ncomp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX,
    Array4<Real const> CoarseToFineProjectionMatrix )
{
    InterpFromCoarseLevel
      ( mf, mf_G, mf.nGrowVect(), Time,
        CrseMF, CrseMF_G, scomp, dcomp, ncomp,
        CrseGeom, FineGeom,
        CrseBC, CrseBCcomp,
        FineBC, FineBCcomp,
        RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix );
}

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel
  ( MF & mf, const MF & mf_G, IntVect const & swX, Real Time,
    const MF & CrseMF, const MF & CrseMF_G,
    int sComp, int dComp, int nComp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX,
    Array4<Real const> CoarseToFineProjectionMatrix )
{
    using FAB = typename MF::FABType::value_type;

    const InterpolaterBoxCoarsener & coarsener
            = mapper->BoxCoarsener( RefRatio );

    const BoxArray            & ba = mf.boxArray();
    const DistributionMapping & dm = mf.DistributionMap();

    const IndexType & typ = ba.ixType();

    BL_ASSERT( typ == CrseMF.boxArray().ixType() );

    Box fdomain_g( amrex::convert( FineGeom.Domain(), mf.ixType() ) );
    for ( int iDimX = 0; iDimX < AMREX_SPACEDIM; ++iDimX ) {
        if ( FineGeom.isPeriodic( iDimX ) ) {
            fdomain_g.grow( iDimX, swX[iDimX] );
        }
    }

    BoxArray ba_crse_patch( ba.size() );
    {  // TODO: later we might want to cache this
        for ( int i = 0, N = ba.size(); i < N; ++i )
        {
            Box bx = amrex::convert(amrex::grow(ba[i],swX), typ);
            bx &= fdomain_g;
            ba_crse_patch.set(i, coarsener.doit(bx));
        }
    }

    MF mf_crse_patch;
    MF mfG_crse_patch;
    mf_crse_patch.define ( ba_crse_patch, dm, nComp, 0 );
    mfG_crse_patch.define( ba_crse_patch, dm, nDOFX, 0 );

    detail::mf_set_domain_bndry( mf_crse_patch , CrseGeom );
    detail::mf_set_domain_bndry( mfG_crse_patch, CrseGeom );

    mf_crse_patch.ParallelCopy
                    ( CrseMF, sComp, 0, nComp, CrseGeom.periodicity() );
    mfG_crse_patch.ParallelCopy
                    ( CrseMF_G, sComp, 0, nDOFX, CrseGeom.periodicity() );

    CrseBC( mf_crse_patch, 0, nComp,
            mf_crse_patch.nGrowVect(), Time, CrseBCcomp);

    FillPatchInterp
      ( mf           , mf_G,
        mf_crse_patch, mfG_crse_patch,
        nComp, swX, CrseGeom,
        fdomain_g,
        RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi( mf ); mfi.isValid(); ++mfi )
    {
        FAB & dfab   = mf[mfi];
        Box dfab_bx = dfab.box();
        dfab_bx.grow(swX-mf.nGrowVect());
        const Box & dbx = dfab_bx & fdomain_g;
    }

    FineBC( mf, dComp, nComp, swX, Time, FineBCcomp );
}

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel
  ( MF & mf, Real Time,
    const MF & CrseMF,
    int scomp, int dcomp, int ncomp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX,
    Array4<Real const> CoarseToFineProjectionMatrix )
{
    InterpFromCoarseLevel
      ( mf, mf.nGrowVect(), Time,
        CrseMF, scomp, dcomp, ncomp,
        CrseGeom, FineGeom,
        CrseBC, CrseBCcomp,
        FineBC, FineBCcomp,
        RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix );
}

template <typename MF, typename BC, typename Interp>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel
  ( MF & mf, IntVect const & swX, Real Time,
    const MF & CrseMF,
    int sComp, int dComp, int nComp,
    const Geometry & CrseGeom, const Geometry & FineGeom,
    BC & CrseBC, int CrseBCcomp,
    BC & FineBC, int FineBCcomp,
    const IntVect & RefRatio,
    Interp * mapper,
    const Vector<BCRec> & bcs, int bcscomp,
    int nDOFX,
    Array4<Real const> CoarseToFineProjectionMatrix )
{
    using FAB = typename MF::FABType::value_type;

    const InterpolaterBoxCoarsener & coarsener
            = mapper->BoxCoarsener( RefRatio );

    const BoxArray            & ba = mf.boxArray();
    const DistributionMapping & dm = mf.DistributionMap();

    const IndexType & typ = ba.ixType();

    BL_ASSERT( typ == CrseMF.boxArray().ixType() );

    Box fdomain_g( amrex::convert( FineGeom.Domain(), mf.ixType() ) );
    for ( int iDimX = 0; iDimX < AMREX_SPACEDIM; ++iDimX ) {
        if ( FineGeom.isPeriodic( iDimX ) ) {
            fdomain_g.grow( iDimX, swX[iDimX] );
        }
    }

    BoxArray ba_crse_patch( ba.size() );
    {  // TODO: later we might want to cache this
        for ( int i = 0, N = ba.size(); i < N; ++i )
        {
            Box bx = amrex::convert(amrex::grow(ba[i],swX), typ);
            bx &= fdomain_g;
            ba_crse_patch.set(i, coarsener.doit(bx));
        }
    }

    MF mf_crse_patch;
    mf_crse_patch.define ( ba_crse_patch, dm, nComp, 0 );

    detail::mf_set_domain_bndry( mf_crse_patch , CrseGeom );

    mf_crse_patch.ParallelCopy
                    ( CrseMF, sComp, 0, nComp, CrseGeom.periodicity() );

    CrseBC( mf_crse_patch, 0, nComp,
            mf_crse_patch.nGrowVect(), Time, CrseBCcomp);

    FillPatchInterp
      ( mf           ,
        mf_crse_patch,
        nComp, swX, CrseGeom,
        fdomain_g,
        RefRatio, mapper, bcs, bcscomp,
        nDOFX, CoarseToFineProjectionMatrix );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi( mf ); mfi.isValid(); ++mfi )
    {
        FAB & dfab   = mf[mfi];
        Box dfab_bx = dfab.box();
        dfab_bx.grow(swX-mf.nGrowVect());
        const Box & dbx = dfab_bx & fdomain_g;
    }

    FineBC( mf, dComp, nComp, swX, Time, FineBCcomp );
}

}

#endif
