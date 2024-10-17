#ifndef THORNADO_AMREX_MULTIFAB_UTIL_3D_C_H_
#define THORNADO_AMREX_MULTIFAB_UTIL_3D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_dg_conservative
    ( int iCrse, int jCrse, int kCrse, int nComp,
      Array4<Real>       const & CrseArr,
      Array4<Real const> const & FineArr,
      Array4<Real const> const & CrseArrG,
      Array4<Real const> const & FineArrG,
      IntVect const & RefRatio,
      int nDOFX,
      Array4<Real const> FineToCoarseProjectionMatrix ) noexcept
  {
    const Real Zero = 0.0;

    // Number of fields in MultiFab
    int nFields = nComp / nDOFX;

    int iRefRatio = RefRatio[0];
    int jRefRatio = RefRatio[1];
    int kRefRatio = RefRatio[2];

    // Corresponding fine element
    const int iFine = iCrse * iRefRatio;
    const int jFine = jCrse * jRefRatio;
    const int kFine = kCrse * kRefRatio;

    // Index for projection matrix
    int iProj;

    // Loop over fields
    for( int iField = 0; iField < nFields; iField++ )
    {
      // Loop over DOF of coarse element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize coarse data to zero
        CrseArr(iCrse,jCrse,kCrse,nDOFX*iField+iNX) = Zero;

        iProj = -1;

        // Loop over fine cells
        for( int kFn = 0; kFn < kRefRatio; kFn++ ) {
        for( int jFn = 0; jFn < jRefRatio; jFn++ ) {
        for( int iFn = 0; iFn < iRefRatio; iFn++ ) {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            CrseArr(iCrse,jCrse,kCrse,nDOFX*iField+iNX)
              += FineToCoarseProjectionMatrix(0,iProj,iNX,jNX)
                   * FineArr (iFine+iFn,jFine+jFn,kFine+kFn,nDOFX*iField+jNX)
                   * FineArrG(iFine+iFn,jFine+jFn,kFine+kFn,jNX);
          } // jNX
        }}} // Fine cells

        CrseArr(iCrse,jCrse,kCrse,nDOFX*iField+iNX)
          /= CrseArrG(iCrse,jCrse,kCrse,iNX);

      } // iNX
    } // iField
  } // END void amrex_avgdown_dg_conservative

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_dg_pointwise
    ( int iCrse, int jCrse, int kCrse, int nComp,
      Array4<Real>       const & CrseArr,
      Array4<Real const> const & FineArr,
      IntVect const & RefRatio,
      int nDOFX,
      Array4<Real const> FineToCoarseProjectionMatrix ) noexcept
  {
    const Real Zero = 0.0;

    // Number of fields in MultiFab
    int nFields = nComp / nDOFX;

    int iRefRatio = RefRatio[0];
    int jRefRatio = RefRatio[1];
    int kRefRatio = RefRatio[2];

    // Corresponding fine element
    const int iFine = iCrse * iRefRatio;
    const int jFine = jCrse * jRefRatio;
    const int kFine = kCrse * kRefRatio;

    // Index for projection matrix
    int iProj;

    // Loop over fields
    for( int iField = 0; iField < nFields; iField++ )
    {
      // Loop over DOF of coarse element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize coarse data to zero
        CrseArr(iCrse,jCrse,kCrse,nDOFX*iField+iNX) = Zero;

        iProj = -1;

        // Loop over fine cells
        for( int kFn = 0; kFn < kRefRatio; kFn++ ) {
        for( int jFn = 0; jFn < jRefRatio; jFn++ ) {
        for( int iFn = 0; iFn < iRefRatio; iFn++ ) {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            CrseArr(iCrse,jCrse,kCrse,nDOFX*iField+iNX)
              += FineToCoarseProjectionMatrix(0,iProj,iNX,jNX)
                   * FineArr(iFine+iFn,jFine+jFn,kFine+kFn,nDOFX*iField+jNX);
          } // jNX
        }}} // Fine cells
      } // iNX
    } // iField
  } // END void amrex_avgdown_dg_pointwise

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_cg
    ( int iCrse, int jCrse, int kCrse, int nComp,
      Array4<Real> const& crse,
      Array4<Real const> const& fine,
      IntVect const& ratio,
      int nDOFX,
      int nFine,
      Array4<Real> G2L,
      Array4<Real> L2G,
      Array4<Real> F2C ) noexcept
  {
    amrex::ignore_unused
    ( iCrse, jCrse, kCrse, nComp,
      crse,
      fine,
      ratio,
      nDOFX,
      nFine,
      G2L,
      L2G,
      F2C );
    amrex::Abort( "Src/Base/AMReX_MultiFabUtil_3D_C.H: amrex_avgdown_cg not implemented. Calling amrex::Abort()" );
  } // END void amrex_avgdown_cg

} // namespace amrex

#endif
