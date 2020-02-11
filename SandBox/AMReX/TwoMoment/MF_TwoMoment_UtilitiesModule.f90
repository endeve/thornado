MODULE MF_TwoMoment_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor


  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeTimeStep


CONTAINS




  SUBROUTINE MF_ComputeTimeStep( TimeStepMin )


    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels-1)
    REAL(amrex_real) :: TimeStep(0:nLevels-1)
    INTEGER :: iLevel   


    TimeStepMin = HUGE( 1.0e0_amrex_real )    

    DO iLevel = 0, nLevels-1

      TimeStep( iLevel ) = 0.1_amrex_real
 
      TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

    END DO




  END SUBROUTINE MF_ComputeTimeStep



END MODULE MF_TwoMoment_UtilitiesModule
