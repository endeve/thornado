#ifndef THORNADO_AMREX_FillPatchUtil_H_
#define THORNADO_AMREX_FillPatchUtil_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <thornado_AMReX_Interpolater.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cmath>
#include <limits>

#include <AMReX_FillPatchUtil.H>

namespace amrex
{
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
      int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix );

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
      int nDOFX, Array4<Real const> CoarseToFineProjectionMatrix );
}

#endif
