
#include <AMReX_BArena.H>
#include <thornado_AMReX_FluxRegister.H>
#include <thornado_AMReX_FluxReg_C.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

int
thornadoFluxRegister::nComp () const noexcept
{
    return ncomp;
}

void
thornadoFluxRegister::CrseInit_DG
  ( const MultiFab & SurfaceFlux,
    int              iDimX,
    int              nFields,
    int              nDOFX_X1,
    int              nDOFX_X2,
    int              nDOFX_X3,
    Real           * WeightsX_X1,
    Real           * WeightsX_X2,
    Real           * WeightsX_X3,
    FrOp             op )
{

    Real * WeightsX_X = nullptr;
    int nDOFX_X       = -1;

    if( iDimX == 0 )
    {
        WeightsX_X = WeightsX_X1;
        nDOFX_X    = nDOFX_X1;
    }
    else if( iDimX == 1 )
    {
        WeightsX_X = WeightsX_X2;
        nDOFX_X    = nDOFX_X2;
    }
    else
    {
        WeightsX_X = WeightsX_X3;
        nDOFX_X    = nDOFX_X3;
    }

    int nComp = nDOFX_X * nFields;

    /* Define MultiFab for FluxRegister */
    MultiFab mf_reg( SurfaceFlux.boxArray(), SurfaceFlux.DistributionMap(),
                     nComp, 0, MFInfo(), SurfaceFlux.Factory() );
    mf_reg.setVal( (Real)0.0 );

    int iX_B0[3];
    int iX_E0[3];

    /* Populate destination MultiFab */
#ifdef AMREX_USE_OMP
#pragma omp parallel private( iX_B0, iX_E0 ) if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi( SurfaceFlux, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
    {
        const Box& bx   = mfi.tilebox();
        auto       reg  = mf_reg.array( mfi );
        auto const sf_C = SurfaceFlux.const_array( mfi );

//        AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, nComp, i, j, k, n,
//        {
//          reg(i,j,k,n) = sf_C(i,j,k,n+SrcComp) * mult;
//        });

        /* bx.loVect() has only amrex_spacedim components */

        iX_B0[0] = bx.loVect()[0];
        iX_E0[0] = bx.hiVect()[0];

        iX_B0[1] = 0;
        iX_E0[1] = 0;
        if( AMREX_SPACEDIM > 1 ) {
        iX_B0[1] = bx.loVect()[1];
        iX_E0[1] = bx.hiVect()[1];
        }

        iX_B0[2] = 0;
        iX_E0[2] = 0;
        if( AMREX_SPACEDIM > 2 ) {
        iX_B0[2] = bx.loVect()[2];
        iX_E0[2] = bx.hiVect()[2];
        }

        for( int iCrse = iX_B0[0]; iCrse <= iX_E0[0]; iCrse++ ) {
        for( int jCrse = iX_B0[1]; jCrse <= iX_E0[1]; jCrse++ ) {
        for( int kCrse = iX_B0[2]; kCrse <= iX_E0[2]; kCrse++ ) {

            for( int iField = 0; iField < nFields; iField++ )
            {
                for( int iNX_X = 0; iNX_X < nDOFX_X; iNX_X++ )
                {
                  reg(iCrse,jCrse,kCrse,iNX_X+iField*nDOFX_X)
                    = -WeightsX_X[iNX_X]
                         * sf_C(iCrse,jCrse,kCrse,iNX_X+iField*nDOFX_X);
                } /* iNX_X */
            } /* iField */

        }}} /* iCrse, jCrse, kCrse */

    } /* END for MFIter  */

    /* face_lo = (0), face_hi = (1) */
    const Orientation face_lo( iDimX, Orientation::low  );
    const Orientation face_hi( iDimX, Orientation::high );

    /* pass <==> which side of face */
    for ( int pass = 0; pass < 2; pass++ )
    {
        /* if pass == 0 then face = face_lo; else face = face_hi */
        const Orientation face = ( ( pass == 0 ) ? face_lo : face_hi );

        if ( op == thornadoFluxRegister::COPY )
        {
            bndry[face].copyFrom( mf_reg, 0, 0, 0, nComp );
        }
// This `else` never happens because `op` is always set to the default
// (i.e., FluxRegister::COPY) in the Fortran interface, so it has been
// commented out
        else
        {
std::cout<<"THIS SHOULD NEVER PRINT! Src/AmrCore/AMReX_FluxRegister.cpp\n";
/*
            FabSet fs( bndry[face].boxArray(), bndry[face].DistributionMap(),
                       nComp );

            fs.setVal( 0 );

            fs.copyFrom( mf, 0, 0, 0, nComp );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                auto const sf_C = fs.const_array(mfi);
                auto       reg  = bndry[face].array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (bx, nComp, i, j, k, n,
                {
                    reg(i,j,k,n) += sf_C(i,j,k,n);
                });
            }
*/
        }
    } /* END for ( int pass = 0; pass < 2; pass++ ) */

} /* END void FluxRegister::CrseInit_DG */

void
thornadoFluxRegister::FineAdd_DG
  ( const MultiFab & SurfaceFluxes,
    int              iDimX,
    int              nFields,
    Real             FaceRatio,
    int              nDOFX_X1,
    int              nDOFX_X2,
    int              nDOFX_X3,
    int              nFineX_X1,
    int              nFineX_X2,
    int              nFineX_X3,
    Real           * WeightsX_X1,
    Real           * WeightsX_X2,
    Real           * WeightsX_X3,
    void           * vpLX_X1_Refined,
    void           * vpLX_X2_Refined,
    void           * vpLX_X3_Refined )
{
    auto *pLX_X1_Refined
           = reinterpret_cast<Real*>(vpLX_X1_Refined);
    Array4<Real> LX_X1_Refined
                   ( pLX_X1_Refined,
                     {0,0,0}, {1,nDOFX_X1,nFineX_X1}, nDOFX_X1 );
    auto *pLX_X2_Refined
           = reinterpret_cast<Real*>(vpLX_X2_Refined);
    Array4<Real> LX_X2_Refined
                   ( pLX_X2_Refined,
                     {0,0,0}, {1,nDOFX_X2,nFineX_X2}, nDOFX_X2 );
    auto *pLX_X3_Refined
           = reinterpret_cast<Real*>(vpLX_X3_Refined);
    Array4<Real> LX_X3_Refined
                   ( pLX_X3_Refined,
                     {0,0,0}, {1,nDOFX_X3,nFineX_X3}, nDOFX_X3 );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(SurfaceFluxes); mfi.isValid(); ++mfi )
    {
        const int k = mfi.index();
        FineAdd_DG( SurfaceFluxes[mfi], iDimX, nFields,
                    FaceRatio, nDOFX_X1, nDOFX_X2, nDOFX_X3,
                    WeightsX_X1, WeightsX_X2, WeightsX_X3,
                    LX_X1_Refined, LX_X2_Refined, LX_X3_Refined,
                    k, RunOn::Gpu );
    }
} /* END void FluxRegister::FineAdd_DG */

void
thornadoFluxRegister::FineAdd_DG
  ( const FArrayBox    & SurfaceFluxes,
    int                  iDimX,
    int                  nFields,
    Real                 FaceRatio,
    int                  nDOFX_X1,
    int                  nDOFX_X2,
    int                  nDOFX_X3,
    Real               * WeightsX_X1,
    Real               * WeightsX_X2,
    Real               * WeightsX_X3,
    Array4<Real const>   LX_X1_Refined,
    Array4<Real const>   LX_X2_Refined,
    Array4<Real const>   LX_X3_Refined,
    int                  BoxNumber,
    RunOn                runon) noexcept
{
    Array4<Real const> LX_X;
    Real * WeightsX_X = nullptr;
    int nDOFX_X       = -1;


    if( iDimX == 0 )
    {
        WeightsX_X = WeightsX_X1;
        LX_X       = LX_X1_Refined;
        nDOFX_X    = nDOFX_X1;
    }
    else if( iDimX == 1 )
    {
        WeightsX_X = WeightsX_X2;
        LX_X       = LX_X2_Refined;
        nDOFX_X    = nDOFX_X2;
    }
    else if( iDimX == 2 )
    {
        WeightsX_X = WeightsX_X3;
        LX_X       = LX_X3_Refined;
        nDOFX_X    = nDOFX_X3;
    }
    else
    {
        amrex::Abort( "Invalid value for iDimX" );
    }

    FArrayBox& loreg = bndry[Orientation(iDimX,Orientation::low)][BoxNumber];
    FArrayBox& hireg = bndry[Orientation(iDimX,Orientation::high)][BoxNumber];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = SurfaceFluxes.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineadd_dg( tlobx, loarr, farr,
                                  iDimX, nFields, nDOFX_X, WeightsX_X, LX_X,
                                  FaceRatio,
                                  local_ratio );
          },
          hibox, thibx,
          {
              fluxreg_fineadd_dg( thibx, hiarr, farr,
                                  iDimX, nFields, nDOFX_X, WeightsX_X, LX_X,
                                  FaceRatio,
                                  local_ratio );
          }
        );
    }
    else
    {
        fluxreg_fineadd_dg( lobox, loarr, farr,
                            iDimX, nFields, nDOFX_X, WeightsX_X, LX_X,
                            FaceRatio,
                            local_ratio );
        fluxreg_fineadd_dg( hibox, hiarr, farr,
                            iDimX, nFields, nDOFX_X, WeightsX_X, LX_X,
                            FaceRatio,
                            local_ratio );
    }

} /* END void FluxRegister::FineAdd_DG */

void
thornadoFluxRegister::Reflux_DG ( MultiFab&       MF_G,
                          MultiFab&       MF_dU,
                          const Geometry& geom,
                          int             nDOFX,
                          int             nDOFX_X1,
                          int             nDOFX_X2,
                          int             nDOFX_X3,
                          int             nFields,
                          int             iGF_SqrtGm,
                          Array4<int>     NodeNumberTableX_X1,
                          Array4<int>     NodeNumberTableX_X2,
                          Array4<int>     NodeNumberTableX_X3,
                          Array4<Real>    WeightsX_q,
                          Array4<Real>    LX_X1_Up,
                          Array4<Real>    LX_X1_Dn,
                          Array4<Real>    LX_X2_Up,
                          Array4<Real>    LX_X2_Dn,
                          Array4<Real>    LX_X3_Up,
                          Array4<Real>    LX_X3_Dn,
                          Real            dX1,
                          Real            dX2,
                          Real            dX3 )
{
    for( OrientationIter fi; fi; ++fi )
    {
        const Orientation& face = fi();
        Reflux_DG( MF_G, MF_dU, geom,
                   nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, nFields, iGF_SqrtGm,
                   NodeNumberTableX_X1, NodeNumberTableX_X2,
                   NodeNumberTableX_X3, WeightsX_q,
                   LX_X1_Up, LX_X1_Dn,
                   LX_X2_Up, LX_X2_Dn,
                   LX_X3_Up, LX_X3_Dn,
                   dX1, dX2, dX3,
                   face );
    }
} /* END void FluxRegister::Reflux_DG */

void
thornadoFluxRegister::Reflux_DG ( MultiFab&       MF_G,
                          MultiFab&       MF_dU,
                          const Geometry& geom,
                          int             nDOFX,
                          int             nDOFX_X1,
                          int             nDOFX_X2,
                          int             nDOFX_X3,
                          int             nFields,
                          int             iGF_SqrtGm,
                          Array4<int>     NodeNumberTableX_X1,
                          Array4<int>     NodeNumberTableX_X2,
                          Array4<int>     NodeNumberTableX_X3,
                          Array4<Real>    WeightsX_q,
                          Array4<Real>    LX_X1_Up,
                          Array4<Real>    LX_X1_Dn,
                          Array4<Real>    LX_X2_Up,
                          Array4<Real>    LX_X2_Dn,
                          Array4<Real>    LX_X3_Up,
                          Array4<Real>    LX_X3_Dn,
                          Real            dX1,
                          Real            dX2,
                          Real            dX3,
                          Orientation     face )
{
    BL_PROFILE("thornadoFluxRegister::Reflux_DG()");

    int nComp = nDOFX_X1 * nFields;
    int iDimX = face.coordDir();

    MultiFab MF_dF( amrex::convert( MF_dU.boxArray(),
                                   IntVect::TheDimensionVector(iDimX) ),
                   MF_dU.DistributionMap(), nComp, 0,
                   MFInfo(), MF_dU.Factory() );
    MF_dF.setVal( (Real)0.0 );
    bndry[face].copyTo( MF_dF, 0, 0, 0, nComp, geom.periodicity() );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for( MFIter mfi( MF_dU, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
    {
        const Box& bx                = mfi.tilebox();
        Array4<Real>       const& G  = MF_G.array(mfi);
        Array4<Real>       const& dU = MF_dU.array(mfi);
        Array4<Real const> const& dF = MF_dF.const_array(mfi);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
        {
            fluxreg_reflux_dg
              ( tbx, G, dU, dF, nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, nFields,
                iGF_SqrtGm,
                NodeNumberTableX_X1, NodeNumberTableX_X2, NodeNumberTableX_X3,
                WeightsX_q,
                LX_X1_Up, LX_X1_Dn, LX_X2_Up, LX_X2_Dn,
                LX_X3_Up, LX_X3_Dn, dX1, dX2, dX3, face );
        });
    }
} /* END void thornadoFluxRegister::Reflux_DG */

}
