#include <AMReX_TimeIntegrator.H>
#include <assert.h>

using namespace amrex;

using Thornado_MultiFabs_t = Vector<MultiFab>;

using Thornado_rhs_funptr_t = void(*)(
    const int*, const int*, const int*, const int*, // iZ_B0, iZ_E0, iZ_B1, iZ_E1
    Real*, Real*, Real*, Real*, Real*,              // GE, GX, U_F, U_R, dU_R
    bool*,                                          // SuppressBC_Option  (match C_BOOL type!)
    Real*, Real*, Real*,                            // SurfaceFlux_X1/X2/X3_Option
    int*);                                          // FineMask_Option

using Thornado_applybc_funptr_t = void(*)(
    const int*, const int*, const int*, const int*, // iZ_B0, iZ_E0, iZ_B1, iZ_E1
    Real*,                                          // U
    int*);                                          // iApplyBC_Option

void amrex2thornado_X
       ( int nFields,
         const amrex::Box& bx,
         amrex::Array4<Real const> const& Data_amrex,
         Real* const Data_thornado /* flat, Fortran-ordered */ )
{
    // Extents of the thornado spatial box (same in both arrays, since
    // iX_B1/iX_E1 == iX_B/iX_E in the Fortran routine).
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);
    const int nDOFX = Data_amrex.nComp() / nFields; // nC = nDOFX * nFields
    const int nX[3] = {hi.x - lo.x + 1,
                       hi.y - lo.y + 1,
                       hi.z - lo.z + 1};

    // Strides for Data_thornado(iDOF, iX1, iX2, iX3, iFd), 1-based dims,
    // column-major: iDOF fastest, iFd slowest.
    const int strideX1 = nDOFX;
    const int strideX2 = strideX1 * nX[0];
    const int strideX3 = strideX2 * nX[1];
    const int strideFd = strideX3 * nX[2];

    amrex::ParallelFor( bx, Data_amrex.nComp(),
      [=] AMREX_GPU_DEVICE (int iX1, int iX2, int iX3, int iC) noexcept
    {
        // Base flat offset (0-based) for Data_thornado(1, iX1, iX2, iX3, iFd).
        const int iFd   = iC / nDOFX;
        const int iDOFX = iC % nDOFX;
        const int base_t =
              (iX1 - lo.x) * strideX1
            + (iX2 - lo.y) * strideX2
            + (iX3 - lo.z) * strideX3
            + (iFd - 0)    * strideFd;

        const int offset_t = iDOFX;
        AMREX_ASSERT(base_t + offset_t < strideFd * nFields );
        Data_thornado[base_t + offset_t] = Data_amrex(iX1, iX2, iX3, iC);
    });
}

void amrex2thornado_Z
       ( int nDOFX, int nDOFE, int nE0, int nFields, int nSpecies, int swE,
         const amrex::Box& bx,
         amrex::Array4<Real const> const& Data_amrex,
         Real* const Data_thornado /* flat, Fortran-ordered */ )
{
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    // Size of Data_amrex:    (                       nZ2, nZ3, nZ4, nComp            )
    // Size of Data_thornado: (nDOFE*nDOFE, nE+swE*2, nZ2, nZ3, nZ4, nFields, nSpecies)
    // Note: Data_thornado includes ghosts in the energy dimension, unlike Data_amrex

    AMREX_ASSERT(Data_amrex.nComp() == nDOFX * nDOFE * nE0 * nFields * nSpecies);

    const int nZ2 = hi.x - lo.x + 1;
    const int nZ3 = hi.y - lo.y + 1;
    const int nZ4 = hi.z - lo.z + 1;

    // Strides for Data_thornado(iNodeZ, iZ1, iZ2, iZ3, iZ4, iFd, iS),
    // column-major: iNodeZ fastest, iS slowest.
    const int nE1 = nE0 + 2*swE;
    const int strideZ1 = nDOFE * nDOFX;
    const int strideZ2 = strideZ1 * nE1;
    const int strideZ3 = strideZ2 * nZ2;
    const int strideZ4 = strideZ3 * nZ3;
    // const int strideFd = strideZ4 * nZ4;
    // const int strideS  = strideFd * nFields;

    const int ncomp = nE0 * nDOFE * nDOFX; // flattened (iZ1,iNodeE,iNodeX)

    for (int iS = 0; iS < nSpecies; ++iS) {
    for (int iFd = 0; iFd < nFields; ++iFd) {
        // Base flat offset (0-based) for Data_thornado(1,iZ_B1(1),iZ_B1(2),
        // iZ_B1(3),iZ_B1(4),iFd,iS); position within (iZ2,iZ3,iZ4,iZ1,iNodeZ)
        // is added inside the kernel.
        const int base_siFd = (iFd + iS * nFields) * strideZ1 * nE1 * nZ2 * nZ3 * nZ4;

        // Component offset (0-based iD) contributed by (iS,iFd), matching
        // the Fortran formula's leading two terms.
        const int iD_siFd = iS  * nFields * ncomp
                          + iFd * ncomp;

        amrex::ParallelFor( bx, ncomp,
          [=] AMREX_GPU_DEVICE (int iZ2, int iZ3, int iZ4, int iC) noexcept
        {
            // Decompose iC (0-based) -> (iZ1, iNodeE, iNodeX), inverse of
            // the Fortran Horner-rule packing (iNodeX fastest).
            int tmp          = iC;
            const int iNodeX = tmp % nDOFX; tmp /= nDOFX;
            const int iNodeE = tmp % nDOFE; tmp /= nDOFE;
            const int iZ1    = tmp;

            const int iNodeZ = iNodeX * nDOFE + iNodeE;

            // 1-based component index into Data_amrex, matching Fortran iD.
            const int iD = iD_siFd + iC;

            const int offset_t = iNodeZ
                               + (iZ1 + swE ) * strideZ1
                               + (iZ2 - lo.x) * strideZ2
                               + (iZ3 - lo.y) * strideZ3
                               + (iZ4 - lo.z) * strideZ4;

            Data_thornado[base_siFd + offset_t] = Data_amrex(iZ2, iZ3, iZ4, iD);
            AMREX_ASSERT(base_siFd + offset_t < nDOFX*nDOFE * nE1 * nZ2 * nZ3 * nZ4 * nFields * nSpecies );
        });

//      DO iZ1    = iE_B0, iE_E0 ! always want iZ1 to not include ghost cells
//      DO iNodeE = 1    , nDOFE
//      DO iNodeX = 1    , nDOFX
//
//        iD = ( iS - 1 ) * nFields * ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
//           + ( iFd - 1 )          * ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
//           + ( iZ1 - iE_B0 )                              * nDOFE * nDOFX &
//           + ( iNodeE - 1 )                                       * nDOFX &
//           +   iNodeX
//        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE !( iNodeX - 1 ) * nDOFE + iNodeE
//
//        Data_thornado(iNodeZ,iZ1,iZ2,iZ3,iZ4,iFd,iS) &
//          = Data_amrex(iZ2,iZ3,iZ4,iD)
    }
    }

}



extern "C" {
  void amrex_fi_integrate_twomoment(
      Geometry *GEOM,
      Real *GE,         // Geometry fields [energy]
      MultiFab* MF_uGF,  // Geometry fields [space]
      MultiFab* MF_uCF,  // Conserved fluid fields
      MultiFab* MF_uCR,  // Conserved radiation fields
      Real t,
      Real dt,
      Thornado_rhs_funptr_t compute_rhs_collision,
      Thornado_rhs_funptr_t compute_rhs_streaming,
      Thornado_applybc_funptr_t applybc_twomoment,
      Thornado_applybc_funptr_t applybc_euler,
      int swX[3], int nGF, int nCF, int nCR, int nDOFX, int nDOFE, int nSpecies,
      int iE_E0, int iE_B0, int iE_E1, int iE_B1) {

#if 0
      if (MF_uGF) {
        amrex::Print() << "nComp: " << MF_uGF->nComp() << "\n";
      } else {
        amrex::Print() << "The pointer to MF_uGF is null!\n";
      }
      std::printf("C++ side ptr: %p\n", MF_uGF);
      MultiFab* mf = static_cast<MultiFab*>(MF_uGF);
      std::printf("vtable/first 8 bytes look like: %p\n", *reinterpret_cast<void**>(mf));
#endif

      // Encapsulate data and set each component to an alias of existing MultiFab
      Thornado_MultiFabs_t fields(3);
      fields[0] = MultiFab(*MF_uGF, make_alias, 0, MF_uGF->nComp());
      fields[1] = MultiFab(*MF_uCF, make_alias, 0, MF_uCF->nComp());
      fields[2] = MultiFab(*MF_uCR, make_alias, 0, MF_uCR->nComp());


#if 0
      if (fields[0].is_finite()) {
        amrex::Print() << "nComp: " << MF_uGF->nComp() << "\n";
      } else {
        amrex::Print() << "The pointer to MF_uGF is null!\n";
      }
      amrex::Print() << "nComp: " << fields[0].nComp() << "\n";
#endif

      /* === Begin RHS function ======================================================= */
        auto rhs_function = [&](
            Thornado_MultiFabs_t& u_rhs,
            Thornado_MultiFabs_t& u_data,
            const Real /* time */) {

          // TODO: Apply boundary conditions
          // TODO:   FillPatch
          // TODO:   ApplyPositivityLimiter
          // TODO:   CreateMesh_MF
          // TODO:   CreateFineMask
          // u_data.FillBoundary(geom.periodicity());

          amrex::Print() << "Starting RHS evaluation\n";
          // Unpack data container to get each MultiFab
          auto& G = u_data[0];
          auto& F = u_data[1];
          auto& R = u_data[2];
          auto& dR = u_rhs[2];
          u_rhs[0].setVal(0);
          u_rhs[1].setVal(0);
          u_rhs[2].setVal(0);
          G.FillBoundary(GEOM->periodicity());
          F.FillBoundary(GEOM->periodicity());
          R.FillBoundary(GEOM->periodicity());

          amrex::Print() << "Iterating over multifabs\n";
          for ( MFIter mfi(R); mfi.isValid(); ++mfi )
          {
            // FIXME: grow tilebox by max(swX)
            const Box& bx = mfi.growntilebox(swX[0]);
            const Box& tb = mfi.tilebox();

            amrex::Print() << "Setting bounds from bx\n";
            const int iX_B0[3] = {lbound(tb).x,
                                  lbound(tb).y,
                                  lbound(tb).z};
            const int iX_E0[3] = {ubound(tb).x,
                                  ubound(tb).y,
                                  ubound(tb).z};
            const int iX_B1[3] = {lbound(tb).x - swX[0],
                                  lbound(tb).y - swX[1],
                                  lbound(tb).z - swX[2]};
            const int iX_E1[3] = {ubound(tb).x + swX[0],
                                  ubound(tb).y + swX[1],
                                  ubound(tb).z + swX[2]};
            const int nE = (iE_E0    - iE_B0    + 1);
            const int nX = (iX_E1[0] - iX_B1[0] + 1)
                         * (iX_E1[1] - iX_B1[1] + 1)
                         * (iX_E1[2] - iX_B1[2] + 1);
            amrex::Print() << "Set bounds from bx: " << nX << " elements\n";

            const int swE = iE_E1 - iE_E0;
            AMREX_ASSERT(swE == iE_B0 - iE_B1);


            // To use the native thornado kernels, the data in/out are passed as
            // multi-dimensional arrays of Reals
            // We allocate the arrays (GX, U_F, U_R, dU_R) for this purpose
            // They should be indexed as the Fortran kernels expect
            Real GX  [nDOFX * nX * nGF                                ];
            Real U_F [nDOFX * nX * nCF                                ];
            Real U_R [nDOFX * nX * nCR * nSpecies * nDOFE * (nE+2*swE)];
            Real dU_R[nDOFX * nX * nCR * nSpecies * nDOFE * (nE+2*swE)];
            amrex::Print() << "Allocated Real arrays in RHS\n";

#if 0
            Real SurfaceFlux_X1[nDOFE * NDOFX_X1, iZ_E0[0]-iZ_B0[0],
                                                  iZ_E0[1]-iZ_B0[1]+1,
                                                  iZ_E0[2]-iZ_B0[2],
                                                  iZ_E0[3]-iZ_B0[3],
                                                  nCR, nSpecies];

            Real SurfaceFlux_X1[nDOFE * NDOFX_X1, iZ_E0[0]-iZ_B0[0],
                                                  iZ_E0[1]-iZ_B0[1],
                                                  iZ_E0[2]-iZ_B0[2]+1,
                                                  iZ_E0[3]-iZ_B0[3],
                                                  nCR, nSpecies];

            Real SurfaceFlux_X1[nDOFE * NDOFX_X1, iZ_E0[0]-iZ_B0[0],
                                                  iZ_E0[1]-iZ_B0[1]
                                                  iZ_E0[2]-iZ_B0[2],
                                                  iZ_E0[3]-iZ_B0[3]+1,,
                                                  nCR, nSpecies];
#endif
            // FIXME: Energy index starts from 0 while spatial index starts from 1
            const int iZ_B0[4] = {iE_B0, iX_B0[0]+1, iX_B0[1]+1, iX_B0[2]+1};
            const int iZ_E0[4] = {iE_E0, iX_E0[0]+1, iX_E0[1]+1, iX_E0[2]+1};
            const int iZ_B1[4] = {iE_B1, iX_B1[0]+1, iX_B1[1]+1, iX_B1[2]+1};
            const int iZ_E1[4] = {iE_E1, iX_E1[0]+1, iX_E1[1]+1, iX_E1[2]+1};


            amrex::Print() << "Converting MultiFab data to Fortran array\n";
            amrex::Print() << "G[mfi] begins " << G[mfi].array().begin << "\n";
            amrex::Print() << "G[mfi] ends " << G[mfi].array().end<< "\n";
            int iLo_MF[4] = {G[mfi].array().begin[0],
                             G[mfi].array().begin[1],
                             G[mfi].array().begin[2],
                             G[mfi].array().begin[3]};
            // int *iLo_MF = G[mfi].array().begin.getVect();
            amrex2thornado_X(nGF, bx, G[mfi].array(), GX);
            amrex::Print() << "Converted G[mfi] to Fortran array\n";

            amrex2thornado_X(nCF, bx, F[mfi].array(), U_F);
            amrex::Print() << "Converted F[mfi] to Fortran array\n";

            // amrex2thornado_Z(nCR, nSpecies, nE, iE_B0, iE_E0, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, R[mfi].dataPtr(), U_R);
            amrex::Print() << "R[mfi] has " << R[mfi].nComp() << " components \n";
            amrex2thornado_Z(nDOFX, nDOFE, nE, nCR, nSpecies, swE, bx, R[mfi].array(), U_R);
            amrex::Print() << "Converted R[mfi] to Fortran array\n";

            /* TODO Apply boundary conditions
            CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

            CALL ApplyBoundaryConditions_TwoMoment_MF &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

            CALL ApplyBoundaryConditions_Euler_MF &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, C, Edge_Map )

             */

            //FIXME: Handle inner BC through edge map
            //FIXME: Handle iApplyBC_Option
            applybc_twomoment(iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, nullptr);
            applybc_euler(iX_B0, iX_E0, iX_B1, iX_E1, U_F, nullptr);

            amrex::Print() << "Computing increment\n";
            compute_rhs_streaming(
                iZ_B0, iZ_E0, iZ_B1, iZ_E1,
                GE, GX, U_F, U_R, dU_R,
                nullptr, nullptr, nullptr, nullptr, nullptr
                );
            // Pass nullptr to optional arguments:
            // SuppressBC_Option,
            // SurfaceFlux_X1/X2/X3_Option,
            // FineMask_Option
            amrex::Print() << "Done computing increment\n";

            /* TODO: Work out remaining procedures in ComputeIncrement_MF
            CALL ComputeOffGridFlux_MF &
                   ( iZ_B0, iZ_E0, &
                     FineMask(iX_B0(1):iX_E0(1), &
                              iX_B0(2):iX_E0(2), &
                              iX_B0(3):iX_E0(3), 1:1), &
                     SurfaceFlux_X1, SurfaceFlux_X2 )

            CALL IncrementOffGridTally_TwoMoment( iLevel, iX_B0, iX_E0 )
            !DO i=1,nCR
            !       OffGridFlux_TwoMoment_MF(i,iLevel) &
            !         = OffGridFlux_TwoMoment_MF(i,iLevel) &
            !         + OffGridFlux_TwoMoment(i)
            !       OffGridFlux_TwoMoment_MF(i+nCR,iLevel) &
            !         = OffGridFlux_TwoMoment_MF(i+nCR,iLevel) &
            !         + OffGridFlux_TwoMoment(i+nCR)
            !END DO
             */

            //FIXME
#if 0
            // Convert data back to AMReX
            // AMReX: dR
            // thornado: dU_R
            thornado2amrex_Z(nCR, nSpecies, nE,
                iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF,
                iZ_B0, iZ_E0,
                dR[mfi].dataPtr(), dU_R);
#endif

            // TODO: Add SurfaceFlux corrections
            // Check with Nil
        } // End for MFIter mfi(R)
      };
      /* === End RHS function ======================================================= */
      amrex::Print() << "Defined RHS function\n";

      // Set up AMReX's time integrator
      TimeIntegrator<Thornado_MultiFabs_t> integrator(fields, t);
      amrex::Print() << "Initialized integrator\n";
      integrator.set_rhs(rhs_function);
      amrex::Print() << "Set RHS function\n";
      integrator.set_time_step(dt);
      amrex::Print() << "Set time step dt\n";

      Real evolution_start_time = ParallelDescriptor::second();

      int nsteps = 1;
      for (int step = 1; step <= nsteps; ++step)
      {
          // Set time to evolve to
          t += dt;

          Real step_start_time = ParallelDescriptor::second();

          // Advance to output time
          amrex::Print() << "Evolving for 1 step\n";
          integrator.evolve(fields, t);

          Real step_stop_time = ParallelDescriptor::second() - step_start_time;
          ParallelDescriptor::ReduceRealMax(step_stop_time);

          // Tell the I/O Processor to write out which step we're doing
          amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds; dt = " << dt << " time = " << t << "\n";

#if 0
          // Write a plotfile of the current data (plot_int was defined in the inputs file)
          if (plot_int > 0 && step%plot_int == 0)
          {
              const std::string& pltfile = amrex::Concatenate("plt",step,5);
              WriteSingleLevelPlotfile(pltfile, u, {"u"}, geom, t, step);
          }
#endif
      }

      Real evolution_stop_time = ParallelDescriptor::second() - evolution_start_time;
      ParallelDescriptor::ReduceRealMax(evolution_stop_time);
      amrex::Print() << "Total evolution time = " << evolution_stop_time << " seconds\n";

  }; // void amrex_fi_integrate

}
