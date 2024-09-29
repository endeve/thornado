#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <thornado_AMReX_Interpolater.H>
#include <thornado_AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

#include <climits>

namespace amrex {

/*
 * DGInterp only works with ref ratio of 2. Not tested for GPU
 *
 * CGInterp only works with ref ratio of 2. Not tested for GPU
 */

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
DGInterp dg_interp;
CGInterp cg_interp;

Box
DGInterp::CoarseBox (const Box& fine,
                           int        ratio)
{
    return amrex::coarsen(fine,ratio);
}

Box
DGInterp::CoarseBox (const Box&     fine,
                           const IntVect& ratio)
{
    return amrex::coarsen(fine,ratio);
}

void
DGInterp::interpConservative
  ( const FArrayBox &  CrseFab                     ,
    const FArrayBox &  CrseFab_G                   ,
    FArrayBox       &  FineFab                     ,
    const FArrayBox &  FineFab_G                   ,
    int                nComp                       ,
    const Box       &  fine_region                 ,
    const IntVect   &  RefRatio                    ,
    int                nDOFX                       ,
    Array4<Real const> CoarseToFineProjectionMatrix,
    RunOn              runon                        )
{
    BL_PROFILE("DGInterp::interpConservative()");

    Array4<Real>       const & FineArr   = FineFab  .      array();
    Array4<Real const> const & FineArr_G = FineFab_G.const_array();
    Array4<Real const> const & CrseArr   = CrseFab  .const_array();
    Array4<Real const> const & CrseArr_G = CrseFab_G.const_array();

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
    {
        amrex::dginterpConservative_interp
          ( tbx, FineArr, FineArr_G, nComp, CrseArr, CrseArr_G, RefRatio,
            nDOFX, CoarseToFineProjectionMatrix );
    });
}

void
DGInterp::interpPointWise
  ( const FArrayBox &  CrseFab                     ,
    FArrayBox       &  FineFab                     ,
    int                nComp                       ,
    const Box       &  fine_region                 ,
    const IntVect   &  RefRatio                    ,
    int                nDOFX                       ,
    Array4<Real const> CoarseToFineProjectionMatrix,
    RunOn              runon                        )
{
    BL_PROFILE("DGInterp::interpPointWise()");

    Array4<Real>       const & FineArr = FineFab.      array();
    Array4<Real const> const & CrseArr = CrseFab.const_array();

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
    {
        amrex::dginterpPointWise_interp
          ( tbx, FineArr, nComp, CrseArr, RefRatio,
            nDOFX, CoarseToFineProjectionMatrix );
    });
}

Box
CGInterp::CoarseBox (const Box& fine,
                           int        ratio)
{
    return amrex::coarsen(fine,ratio);
}

Box
CGInterp::CoarseBox (const Box&     fine,
                           const IntVect& ratio)
{
    return amrex::coarsen(fine,ratio);
}

//void
//CGInterp::interp
//  ( const FArrayBox     & /*crse        */,
//    int                   /*crse_comp   */,
//    FArrayBox           & /*fine        */,
//    int                   /*fine_comp   */,
//    int                   /*ncomp       */,
//    const Box           & /*fine_region */,
//    const IntVect       & /*ratio       */,
//    const Geometry      & /*crse_geom   */,
//    const Geometry      & /*fine_geom   */,
//    Vector<BCRec> const & /*bcr         */,
//    int                   /*actual_comp */,
//    int                   /*actual_state*/,
//    RunOn                 /*runon       */ )
//{
//    BL_PROFILE("CGInterp::interp()");
///*
//    int nDOFX = amrex::DG::nDOFX;
//    int nFine = amrex::DG::nFineV;
//
//    auto *pProjectionMatrix
//           = reinterpret_cast<Real*>(amrex::DG::ProjectionMatrix1D);
//    Array4<Real> ProjectionMatrix
//                   ( pProjectionMatrix, {0,0,0}, {nFine,nDOFX,nDOFX}, 1 );
//
//    for( int iFine = 0; iFine < nFine; iFine++ ) {
//    for( int iNX   = 0; iNX   < nDOFX; iNX++   ) {
//    for( int jNX   = 0; jNX   < nDOFX; jNX++   ) {
//        ProjectionMatrix(iFine,iNX,jNX,0)
//          = amrex::DG::ProjectionMatrix[iFine][iNX][jNX];
//    }}}
//
//    auto *pWeightsX_q
//           = reinterpret_cast<Real*>(amrex::DG::WeightsX_q);
//    Array4<Real> WeightsX_q( pWeightsX_q, {0,0,0}, {nDOFX,1,1}, 1 );
//    for( int iNX = 0; iNX < amrex::DG::nDOFX; iNX++ ) {
//      WeightsX_q(iNX,0,0,0) = amrex::DG::WeightsX_q[iNX];
//    }
//
//    Array4<Real const> const& crsearr = crse.const_array();
//    Array4<Real> const& finearr = fine.array();;
//
//    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
//    {
//        amrex::cginterp_interp
//          ( tbx, finearr, fine_comp, ncomp, crsearr, crse_comp, ratio,
//            nDOFX, ProjectionMatrix, WeightsX_q );
//    });
//*/
//
//}
void
CGInterp::interpConservative
  ( const FArrayBox &  CrseFab                     ,
    const FArrayBox &  CrseFab_G                   ,
    FArrayBox       &  FineFab                     ,
    const FArrayBox &  FineFab_G                   ,
    int                nComp                       ,
    const Box       &  fine_region                 ,
    const IntVect   &  RefRatio                    ,
    int                nDOFX                       ,
    Array4<Real const> CoarseToFineProjectionMatrix,
    RunOn              runon                        )
{
//    BL_PROFILE("CGInterp::interpConservative()");
//
//    Array4<Real>       const & FineArr   = FineFab  .      array();
//    Array4<Real const> const & FineArr_G = FineFab_G.const_array();
//    Array4<Real const> const & CrseArr   = CrseFab  .const_array();
//    Array4<Real const> const & CrseArr_G = CrseFab_G.const_array();
//
//    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
//    {
//        amrex::cginterpConservative_interp
//          ( tbx, FineArr, FineArr_G, nComp, CrseArr, CrseArr_G, RefRatio,
//            nDOFX, CoarseToFineProjectionMatrix );
//    });
}

void
CGInterp::interpPointWise
  ( const FArrayBox &  CrseFab                     ,
    FArrayBox       &  FineFab                     ,
    int                nComp                       ,
    const Box       &  fine_region                 ,
    const IntVect   &  RefRatio                    ,
    int                nDOFX                       ,
    Array4<Real const> CoarseToFineProjectionMatrix,
    RunOn              runon                        )
{
//    BL_PROFILE("CGInterp::interpPointWise()");
//
//    Array4<Real>       const & FineArr = FineFab.      array();
//    Array4<Real const> const & CrseArr = CrseFab.const_array();
//
//    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
//    {
//        amrex::cginterpPointWise_interp
//          ( tbx, FineArr, nComp, CrseArr, RefRatio,
//            nDOFX, CoarseToFineProjectionMatrix );
//    });
}

}
