#include <AMReX_TimeIntegrator.H>
#include <assert.h>

using namespace amrex;

using Thornado_MultiFabs_t = Vector<MultiFab>;

using Thornado_rhs_funptr_t = void(*)(
    const int*, const int*, const int*, const int*, // iZ_B0, iZ_E0, iZ_B1, iZ_E1
    bool*, bool*,                                   // lo, hi
    Real*, Real*, Real*, Real*,                     // GE, GX, U_F, U_R, dU_R
    bool*);                                         // Verbose_Option

using Thornado_pl_funptr_t = void(*)(
    const int*, const int*, const int*, const int*, // iZ_B0, iZ_E0, iZ_B1, iZ_E1
    bool*, bool*,                                   // lo, hi
    Real*, Real*, Real*,                            // GE, GX, U_F, U_R
    int*);                                          // swX_Option

void have_same_bounds(const Box& b1, const Box& b2, bool lo[3], bool hi[3]){
  for (int i=0; i<3; ++i) {
    lo[i] = b1.loVect3d()[i] == b2.loVect3d()[i];
    hi[i] = b1.hiVect3d()[i] == b2.hiVect3d()[i];
  }
}

extern "C" {
  void amrex_fi_integrate_twomoment(
      Geometry *GEOM,
      Real *GE,          // Geometry fields [energy]
      MultiFab* MF_uGF,  // Geometry fields [space]
      MultiFab* MF_uCF,  // Conserved fluid fields
      MultiFab* MF_uCR,  // Conserved radiation fields
      Real t,
      Real dt,
      Thornado_rhs_funptr_t compute_rhs_collision,
      Thornado_rhs_funptr_t compute_rhs_streaming,
      Thornado_pl_funptr_t applypositivitylimiter,
      int swX[3], int iE_E0, int iE_B0, int iE_E1, int iE_B1) {

      // Encapsulate data and set each component to an alias of existing MultiFab
      Thornado_MultiFabs_t fields(3);
      fields[0] = MultiFab(*MF_uGF, make_alias, 0, MF_uGF->nComp());
      fields[1] = MultiFab(*MF_uCF, make_alias, 0, MF_uCF->nComp());
      fields[2] = MultiFab(*MF_uCR, make_alias, 0, MF_uCR->nComp());

      /* === Begin RHS function ======================================================= */
      auto rhs_function = [&](
          Thornado_MultiFabs_t& u_rhs,
          Thornado_MultiFabs_t& u_data,
          const Real /* time */) {

          // Unpack data container to get each MultiFab
          auto& G = u_data[0];
          auto& F = u_data[1];
          auto& R = u_data[2];
          auto& dR = u_rhs[2];
          u_rhs[0].setVal(0);
          u_rhs[1].setVal(0);
          u_rhs[2].setVal(0);

          // Periodic boundary conditions on lower end
          G.FillBoundary(GEOM->periodicity());
          F.FillBoundary(GEOM->periodicity());
          R.FillBoundary(GEOM->periodicity());

          for ( MFIter mfi(R); mfi.isValid(); ++mfi )
          {
            const Box& tb = mfi.tilebox();
            const Box& bx = mfi.validbox();

            // FIXME: Reconcile Fortran and C++ array indexing start (0 vs 1)
            const int iZ_B0[4] = {iE_B0,
                                  1 + lbound(tb).x,
                                  1 + lbound(tb).y,
                                  1 + lbound(tb).z};
            const int iZ_E0[4] = {iE_E0,
                                  1 + ubound(tb).x,
                                  1 + ubound(tb).y,
                                  1 + ubound(tb).z};
            const int iZ_B1[4] = {iE_B1,
                                  1 + lbound(tb).x - swX[0],
                                  1 + lbound(tb).y - swX[1],
                                  1 + lbound(tb).z - swX[2]};
            const int iZ_E1[4] = {iE_E1,
                                  1 + ubound(tb).x + swX[0],
                                  1 + ubound(tb).y + swX[1],
                                  1 + ubound(tb).z + swX[2]};

            bool lo[3], hi[3];
            have_same_bounds(tb, GEOM->Domain(), lo, hi);

            amrex::Print() << "Computing increment for streaming\n";
            compute_rhs_streaming(
                iZ_B0, iZ_E0, iZ_B1, iZ_E1, lo, hi,
                G[mfi].dataPtr(), F[mfi].dataPtr(),
                R[mfi].dataPtr(), dR[mfi].dataPtr(), nullptr
                );
        } // End for MFIter mfi(R)
      };
      /* === End RHS function ======================================================= */

      /* === Begin post-stage ======================================================= */
      auto post_stage_function = [&](
        Thornado_MultiFabs_t& u_data,
        const Real /* time */) {
          auto& G = u_data[0];
          auto& F = u_data[1];
          auto& R = u_data[2];

          for ( MFIter mfi(R); mfi.isValid(); ++mfi )
          {
            const Box& tb = mfi.tilebox();
            const Box& bx = mfi.validbox();

            // FIXME: Reconcile Fortran and C++ array indexing start (0 vs 1)
            const int iZ_B0[4] = {iE_B0,
                                  1 + lbound(tb).x,
                                  1 + lbound(tb).y,
                                  1 + lbound(tb).z};
            const int iZ_E0[4] = {iE_E0,
                                  1 + ubound(tb).x,
                                  1 + ubound(tb).y,
                                  1 + ubound(tb).z};
            const int iZ_B1[4] = {iE_B1,
                                  1 + lbound(tb).x - swX[0],
                                  1 + lbound(tb).y - swX[1],
                                  1 + lbound(tb).z - swX[2]};
            const int iZ_E1[4] = {iE_E1,
                                  1 + ubound(tb).x + swX[0],
                                  1 + ubound(tb).y + swX[1],
                                  1 + ubound(tb).z + swX[2]};

            bool lo[3], hi[3];
            have_same_bounds(tb, GEOM->Domain(), lo, hi);

            amrex::Print() << "Applying positivity limiter\n";
            applypositivitylimiter(
                iZ_B0, iZ_E0, iZ_B1, iZ_E1, lo, hi,
                G[mfi].dataPtr(),
                F[mfi].dataPtr(),
                R[mfi].dataPtr(),
                nullptr);
          } // End for MFIter mfi(R)
      };
      /* === End post-stage function ======================================================= */

      // Set up AMReX's time integrator
      TimeIntegrator<Thornado_MultiFabs_t> integrator(fields, t);
      integrator.set_rhs(rhs_function);
      integrator.set_post_stage_action(post_stage_function);
      integrator.set_post_step_action(post_stage_function);
      integrator.set_time_step(dt);

      Real evolution_start_time = ParallelDescriptor::second();

      int nsteps = 1;
      for (int step = 1; step <= nsteps; ++step)
      {
          // Set time to evolve to
          t += dt;

          Real step_start_time = ParallelDescriptor::second();

          // Advance to output time
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
