#include <AMReX_FPhysBC.H>
#include <thornado_AMReX_FillPatchUtil.H>
#include <thornado_AMReX_FillPatchUtil_I.H>
#include <thornado_AMReX_Interpolater.H>

using namespace amrex;

namespace
{
  // THIS MUST BE CONSISTENT WITH amrex_interpolater_module in AMReX_interpolater_mod.F90!!!
  Vector<thornadoInterpolater*> interp = {
      &amrex::dg_interp, // 0
      &amrex::cg_interp  // 1
  };
}

namespace {

  extern "C"
  {

      void amrex_fi_fillpatch_dgconservative_two
             ( MultiFab * MF, MultiFab * MF_G, Real Time,
               MultiFab * pCrseMF[], MultiFab * pCrseMF_G[],
               Real CrseTime[], int nCrse,
               MultiFab * pFineMF[], MultiFab * pFineMF_G[],
               Real FineTime[], int nFine,
               int sComp, int dComp, int nComp,
               const Geometry* pCrseGeom, const Geometry * pFineGeom,
               FPhysBC::fill_physbc_funptr_t fpCrseFillPhysBC,
               FPhysBC::fill_physbc_funptr_t fpFineFillPhysBC,
               int RefRatio, int interp_id,
               int * pLoBC[], int * pHiBC[],
               int nFineV, int nDOFX, void * vpCoarseToFineProjectionMatrix )
      {
          Vector<BCRec> bcs;
          for ( int iComp = 0; iComp < nComp; ++iComp ) {
              bcs.emplace_back( pLoBC[iComp+sComp], pHiBC[iComp+sComp] );
          }

          FPhysBC CrseBC( fpCrseFillPhysBC, pCrseGeom );
          FPhysBC FineBC( fpFineFillPhysBC, pFineGeom );

          auto *pCoarseToFineProjectionMatrix
                 = reinterpret_cast<Real*>(vpCoarseToFineProjectionMatrix);
          Array4<Real const> CoarseToFineProjectionMatrix
                               ( pCoarseToFineProjectionMatrix,
                                 {0,0,0}, {1,nFineV,nDOFX}, nDOFX );

          amrex::FillPatchTwoLevels
                   ( *MF, *MF_G, Time,
                     Vector<MultiFab*>{pCrseMF  , pCrseMF  +nCrse},
                     Vector<MultiFab*>{pCrseMF_G, pCrseMF_G+nCrse},
                     Vector<Real>     {CrseTime , CrseTime +nCrse},
                     Vector<MultiFab*>{pFineMF  , pFineMF  +nFine},
                     Vector<MultiFab*>{pFineMF_G, pFineMF_G+nFine},
                     Vector<Real>     {FineTime , FineTime +nFine},
                     sComp, dComp, nComp,
                     *pCrseGeom, *pFineGeom,
                     CrseBC, 0, FineBC, 0,
                     IntVect{AMREX_D_DECL(RefRatio,RefRatio,RefRatio)},
                     interp[interp_id], bcs, 0,
                     nDOFX, CoarseToFineProjectionMatrix );
      }

      void amrex_fi_fillpatch_dgpointwise_two
             ( MultiFab * MF, Real Time,
               MultiFab * pCrseMF[],
               Real CrseTime[], int nCrse,
               MultiFab * pFineMF[],
               Real FineTime[], int nFine,
               int sComp, int dComp, int nComp,
               const Geometry* pCrseGeom, const Geometry * pFineGeom,
               FPhysBC::fill_physbc_funptr_t fpCrseFillPhysBC,
               FPhysBC::fill_physbc_funptr_t fpFineFillPhysBC,
               int RefRatio, int interp_id,
               int * pLoBC[], int * pHiBC[],
               int nFineV, int nDOFX, void * vpCoarseToFineProjectionMatrix )
      {
          Vector<BCRec> bcs;
          for ( int iComp = 0; iComp < nComp; ++iComp ) {
              bcs.emplace_back( pLoBC[iComp+sComp], pHiBC[iComp+sComp] );
          }

          FPhysBC CrseBC( fpCrseFillPhysBC, pCrseGeom );
          FPhysBC FineBC( fpFineFillPhysBC, pFineGeom );

          auto *pCoarseToFineProjectionMatrix
                 = reinterpret_cast<Real*>(vpCoarseToFineProjectionMatrix);
          Array4<Real const> CoarseToFineProjectionMatrix
                               ( pCoarseToFineProjectionMatrix,
                                 {0,0,0}, {1,nFineV,nDOFX}, nDOFX );

          amrex::FillPatchTwoLevels
                   ( *MF, Time,
                     Vector<MultiFab*>{pCrseMF  , pCrseMF  +nCrse},
                     Vector<Real>     {CrseTime , CrseTime +nCrse},
                     Vector<MultiFab*>{pFineMF  , pFineMF  +nFine},
                     Vector<Real>     {FineTime , FineTime +nFine},
                     sComp, dComp, nComp,
                     *pCrseGeom, *pFineGeom,
                     CrseBC, 0, FineBC, 0,
                     IntVect{AMREX_D_DECL(RefRatio,RefRatio,RefRatio)},
                     interp[interp_id], bcs, 0,
                     nDOFX, CoarseToFineProjectionMatrix );
      }

      void amrex_fi_fillcoarsepatch_dgconservative
             ( MultiFab * MF, MultiFab * MF_G, Real Time,
               const MultiFab * pCrseMF, const MultiFab * pCrseMF_G,
               int sComp, int dComp, int nComp,
               const Geometry * pCrseGeom, const Geometry * pFineGeom,
               FPhysBC::fill_physbc_funptr_t fpFillPhysBCCrse,
               FPhysBC::fill_physbc_funptr_t fpFillPhysBCFine,
               int RefRatio, int interp_id,
               int * pLoBC[], int * pHiBC[],
               int nFineV, int nDOFX, void * vpCoarseToFineProjectionMatrix )
      {
          Vector<BCRec> bcs;
          for ( int iComp = 0; iComp < nComp; ++iComp) {
              bcs.emplace_back( pLoBC[iComp+sComp], pHiBC[iComp+sComp] );
          }

          FPhysBC CrseBC( fpFillPhysBCCrse, pCrseGeom );
          FPhysBC FineBC( fpFillPhysBCFine, pFineGeom );

          auto *pCoarseToFineProjectionMatrix
                 = reinterpret_cast<Real*>(vpCoarseToFineProjectionMatrix);
          Array4<Real const> CoarseToFineProjectionMatrix
                               ( pCoarseToFineProjectionMatrix,
                                 {0,0,0}, {1,nFineV,nDOFX}, nDOFX );

          amrex::InterpFromCoarseLevel
                   ( *MF, *MF_G, Time, *pCrseMF, *pCrseMF_G,
                     sComp, dComp, nComp,
                     *pCrseGeom, *pFineGeom,
                     CrseBC, 0, FineBC, 0,
                     IntVect{AMREX_D_DECL(RefRatio,RefRatio,RefRatio)},
                     interp[interp_id], bcs, 0,
                     nDOFX, CoarseToFineProjectionMatrix );
      }

      void amrex_fi_fillcoarsepatch_dgpointwise
             ( MultiFab * MF, Real Time,
               const MultiFab * pCrseMF,
               int sComp, int dComp, int nComp,
               const Geometry * pCrseGeom, const Geometry * pFineGeom,
               FPhysBC::fill_physbc_funptr_t fpFillPhysBCCrse,
               FPhysBC::fill_physbc_funptr_t fpFillPhysBCFine,
               int RefRatio, int interp_id,
               int * pLoBC[], int * pHiBC[],
               int nFineV, int nDOFX, void * vpCoarseToFineProjectionMatrix )
      {
          Vector<BCRec> bcs;
          for ( int iComp = 0; iComp < nComp; ++iComp) {
              bcs.emplace_back( pLoBC[iComp+sComp], pHiBC[iComp+sComp] );
          }

          FPhysBC CrseBC( fpFillPhysBCCrse, pCrseGeom );
          FPhysBC FineBC( fpFillPhysBCFine, pFineGeom );

          auto *pCoarseToFineProjectionMatrix
                 = reinterpret_cast<Real*>(vpCoarseToFineProjectionMatrix);
          Array4<Real const> CoarseToFineProjectionMatrix
                               ( pCoarseToFineProjectionMatrix,
                                 {0,0,0}, {1,nFineV,nDOFX}, nDOFX );

          amrex::InterpFromCoarseLevel
                   ( *MF, Time, *pCrseMF,
                     sComp, dComp, nComp,
                     *pCrseGeom, *pFineGeom,
                     CrseBC, 0, FineBC, 0,
                     IntVect{AMREX_D_DECL(RefRatio,RefRatio,RefRatio)},
                     interp[interp_id], bcs, 0,
                     nDOFX, CoarseToFineProjectionMatrix );
      }

  } // end extern C
} // end namespace
