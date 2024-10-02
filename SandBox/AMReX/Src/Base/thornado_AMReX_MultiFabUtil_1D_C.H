#ifndef THORNADO_AMREX_MULTIFAB_UTIL_1D_C_H_
#define THORNADO_AMREX_MULTIFAB_UTIL_1D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_dg_conservative
    ( int iCrse, int, int, int nComp,
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

    // Corresponding fine element
    const int iFine = iCrse * iRefRatio;

    // Index for projection matrix
    int iProj;

    // Loop over fields
    for( int iField = 0; iField < nFields; iField++ )
    {
      // Loop over DOF of coarse element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize coarse data to zero
        CrseArr(iCrse,0,0,nDOFX*iField+iNX) = Zero;

        iProj = -1;

        // Loop over fine cells
        for( int iFn = 0; iFn < iRefRatio; iFn++ ) {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            CrseArr(iCrse,0,0,nDOFX*iField+iNX)
              += FineToCoarseProjectionMatrix(0,iProj,iNX,jNX)
                   * FineArr (iFine+iFn,0,0,nDOFX*iField+jNX)
                   * FineArrG(iFine+iFn,0,0,jNX);
          } // jNX
        } // Fine cells

        CrseArr(iCrse,0,0,nDOFX*iField+iNX)
          /= CrseArrG(iCrse,0,0,iNX);

      } // iNX
    } // iField
  } // END void amrex_avgdown_dg_conservative

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_dg_pointwise
    ( int iCrse, int, int, int nComp,
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

    // Corresponding fine element
    const int iFine = iCrse * iRefRatio;

    // Index for projection matrix
    int iProj;

    // Loop over fields
    for( int iField = 0; iField < nFields; iField++ )
    {
      // Loop over DOF of coarse element
      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        // Initialize coarse data to zero
        CrseArr(iCrse,0,0,nDOFX*iField+iNX) = Zero;

        iProj = -1;

        // Loop over fine cells
        for( int iFn = 0; iFn < iRefRatio; iFn++ ) {

          iProj += 1;

          // Project fine data onto coarse data
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            CrseArr(iCrse,0,0,nDOFX*iField+iNX)
              += FineToCoarseProjectionMatrix(0,iProj,iNX,jNX)
                   * FineArr(iFine+iFn,0,0,nDOFX*iField+jNX);
          } // jNX
        } // Fine cells
      } // iNX
    } // iField
  } // END void amrex_avgdown_dg_pointwise

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void amrex_avgdown_cg
    ( int iCrse, int, int, int nComp,
      Array4<Real>       const & CrseArr,
      Array4<Real const> const & FineArr,
      IntVect const & RefRatio,
      int nDOFX,
      int nFine,
      Array4<Real> G2L,
      Array4<Real> L2G,
      Array4<Real> F2C ) noexcept
  {
    amrex::ignore_unused
      ( iCrse, nComp, CrseArr, FineArr, RefRatio, nDOFX, nFine, G2L, L2G, F2C );
  /*
    const Real Zero = 0.0;

    // Corresponding fine element
    int iRefRatio = RefRatio[0];
    const int iFine = iCrse * iRefRatio;

    // Number of fields in MultiFab
    const int nGF = nComp / nDOFX;

    Real * GhG = new(std::nothrow) Real[nDOFX];
    Real * GhL = new(std::nothrow) Real[nDOFX];

    Real ** ghG = new(std::nothrow) Real*[nFine];
    Real ** ghL = new(std::nothrow) Real*[nFine];
    for( int iFn = 0; iFn < nFine; iFn++ )
    {
      ghG[iFn] = new(std::nothrow) Real[nDOFX];
      ghL[iFn] = new(std::nothrow) Real[nDOFX];
    }

    for( int iGF = 0; iGF < nGF; iGF++ )
    {

      for( int iFn = 0; iFn < nFine; iFn++ )
      {
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {
          ghG[iFn][iNX] = FineArr(iFine+iFn,0,0,nDOFX*iGF+iNX);
        }
      }

      for( int iFn = 0; iFn < nFine; iFn++ )
      {
        for( int iNX = 0; iNX < nDOFX; iNX++ )
        {
          ghL[iFn][iNX] = Zero;
          for( int jNX = 0; jNX < nDOFX; jNX++ )
          {
            ghL[iFn][iNX] += G2L(0,0,iNX,jNX) * ghG[iFn][jNX];
          }
        }
      }

      for( int i = 0; i < nDOFX; i++ )
      {
        GhL[i] = Zero;
        for( int iFn = 0; iFn < nFine; iFn++ )
        {
          for( int k = 0; k < nDOFX; k++ )
          {
            GhL[i] += F2C(0,i,iFn,k) * ghL[iFn][k];
          }
        }
      }

      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        GhG[iNX] = Zero;
        for( int jNX = 0; jNX < nDOFX; jNX++ )
        {
          GhG[iNX] += L2G(0,0,iNX,jNX) * GhL[jNX];
        }
      }

      for( int iNX = 0; iNX < nDOFX; iNX++ )
      {
        CrseArr(iCrse,0,0,nDOFX*iGF+iNX) = GhG[iNX];
      }

    } // iGF

    for( int iFn = 0; iFn < nFine; iFn++ )
    {
      delete[] ghL[iFn];
      delete[] ghG[iFn];
    }
    delete[] ghL;
    delete[] ghG;
    delete[] GhL;
    delete[] GhG;
  */
  } // END void amrex_avgdown_cg

} // end namespace amrex

#endif
