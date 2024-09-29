#ifndef THORNADO_AMREX_FLUXREG_C_H_
#define THORNADO_AMREX_FLUXREG_C_H_
#include <AMReX_Config.H>

#if   ( AMREX_SPACEDIM == 1 )
#include <thornado_AMReX_FluxReg_1D_C.H>
#elif ( AMREX_SPACEDIM == 2 )
#include <thornado_AMReX_FluxReg_2D_C.H>
#else
#include <thornado_AMReX_FluxReg_3D_C.H>
#endif

#endif
