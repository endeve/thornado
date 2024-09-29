#ifndef THORNADO_AMREX_INTERP_2D_C_H_
#define THORNADO_AMREX_INTERP_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterpConservative_interp
  ( Box                const & FineBox,
    Array4<Real>       const & FineArr,
    Array4<Real const> const & FineArrG,
    int                const   nComp,
    Array4<Real const> const & CrseArr,
    Array4<Real const> const & CrseArrG,
    IntVect            const & RefRatio,
    int                        nDOFX,
    Array4<Real const>         CoarseToFineProjectionMatrix ) noexcept
{
  const Real Zero = 0.0;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound( FineBox );
  const auto hiFine = amrex::ubound( FineBox );

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {
    // Loop over fine cells in y-direction
    for( int jFine = loFine.y; jFine <= hiFine.y; jFine++ )
    {
      // Get coarse element corresponding to fine element jFine
      const int jCrse = amrex::coarsen( jFine, RefRatio[1] );

      // Loop over fine cells in x-direction
      for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
      {
        // Get coarse element corresponding to fine element iFine
        const int iCrse = amrex::coarsen( iFine, RefRatio[0] );

        if     ( ( iFine % 2 == 0 ) &
                 ( jFine % 2 == 0 ) ) { iProj = 0; }
        else if( ( iFine % 2 != 0 ) &
                 ( jFine % 2 == 0 ) ) { iProj = 1; }
        else if( ( iFine % 2 == 0 ) &
                 ( jFine % 2 != 0 ) ) { iProj = 2; }
        else                          { iProj = 3; }

        // Loop over DOF of fine element
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {
          // Initialize fine data to zero
          FineArr(iFine,jFine,0,nDOFX*iField+iNX) = Zero;

          // Project coarse data onto fine data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            FineArr(iFine,jFine,0,nDOFX*iField+iNX)
              += CoarseToFineProjectionMatrix(0,iProj,iNX,jNX)
                   * CrseArr (iCrse,jCrse,0,nDOFX*iField+jNX)
                   * CrseArrG(iCrse,jCrse,0,jNX);
          } // jNX

          FineArr(iFine,jFine,0,nDOFX*iField+iNX)
            /= FineArrG(iFine,jFine,0,iNX);

        } // iNX
      } // iFine
    } // jFine
  } // iField
} // end void dginterpConservative_interp

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dginterpPointWise_interp
  ( Box                const & FineBox,
    Array4<Real>       const & FineArr,
    int                const   nComp,
    Array4<Real const> const & CrseArr,
    IntVect            const & RefRatio,
    int                        nDOFX,
    Array4<Real const>         CoarseToFineProjectionMatrix ) noexcept
{
  const Real Zero = 0.0;

  // Get bounds on fine grid
  const auto loFine = amrex::lbound( FineBox );
  const auto hiFine = amrex::ubound( FineBox );

  // Number of fields in MultiFab
  int nFields = nComp / nDOFX;

  // Index for projection matrix
  int iProj;

  // Loop over fields
  for( int iField = 0; iField < nFields; iField++ )
  {
    // Loop over fine cells in y-direction
    for( int jFine = loFine.y; jFine <= hiFine.y; jFine++ )
    {
      // Get coarse element corresponding to fine element jFine
      const int jCrse = amrex::coarsen( jFine, RefRatio[1] );

      // Loop over fine cells in x-direction
      for( int iFine = loFine.x; iFine <= hiFine.x; iFine++ )
      {
        // Get coarse element corresponding to fine element iFine
        const int iCrse = amrex::coarsen( iFine, RefRatio[0] );

        if     ( ( iFine % 2 == 0 ) &
                 ( jFine % 2 == 0 ) ) { iProj = 0; }
        else if( ( iFine % 2 != 0 ) &
                 ( jFine % 2 == 0 ) ) { iProj = 1; }
        else if( ( iFine % 2 == 0 ) &
                 ( jFine % 2 != 0 ) ) { iProj = 2; }
        else                          { iProj = 3; }

        // Loop over DOF of fine element
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {
          // Initialize fine data to zero
          FineArr(iFine,jFine,0,nDOFX*iField+iNX) = Zero;

          // Project coarse data onto fine data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            FineArr(iFine,jFine,0,nDOFX*iField+iNX)
              += CoarseToFineProjectionMatrix(0,iProj,iNX,jNX)
                   * CrseArr(iCrse,jCrse,0,nDOFX*iField+jNX);
          } // jNX
        } // iNX
      } // iFine
    } // jFine
  } // iField
} // end void dginterpPointWise_interp

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cginterpConservative_interp
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio,
    int nDOFX, Array4<Real> ProjectionMatrix,
    Array4<Real> WeightsX_q ) noexcept
{
  amrex::ignore_unused
  ( bx,
    fine, fcomp, ncomp,
    crse, ccomp,
    ratio,
    nDOFX, ProjectionMatrix,
    WeightsX_q );
  amrex::Abort( "Src/AmrCore/AMReX_Interp_2D_C.H: cginterp_interp not implemented. Calling amrex::Abort()" );
} // end void cginterpConservative_interp

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cginterpPointWise_interp
  ( Box const& bx,
    Array4<Real> const& fine, const int fcomp, const int ncomp,
    Array4<Real const> const& crse, const int ccomp,
    IntVect const& ratio,
    int nDOFX, Array4<Real> ProjectionMatrix,
    Array4<Real> WeightsX_q ) noexcept
{
  amrex::ignore_unused
  ( bx,
    fine, fcomp, ncomp,
    crse, ccomp,
    ratio,
    nDOFX, ProjectionMatrix,
    WeightsX_q );
  amrex::Abort( "Src/AmrCore/AMReX_Interp_2D_C.H: cginterp_interp not implemented. Calling amrex::Abort()" );
} // end void cginterpPointWise_interp

}  // namespace amrex

#endif
