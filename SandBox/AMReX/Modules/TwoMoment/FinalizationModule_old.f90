MODULE FinalizationModule

  ! --- AMReX Modules ---
  USE amrex_amr_module,     ONLY: &
    amrex_geometry, &
    amrex_geometry_destroy, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---
  USE RadiationFieldsModule,                ONLY: &
    DestroyRadiationFields
  USE TwoMoment_TimersModule, ONLY: &
    FinalizeTimers

  ! --- Local Modules ---
  USE InputParsingModule,                 ONLY: &
    nLevels, MyAmrFinalize
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_TwoMoment_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_TimeSteppingModule_Relativistic,    ONLY: &
    Finalize_IMEX_RK_MF


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram



 CONTAINS



  SUBROUTINE FinalizeProgram(GEOM)


    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)

    INTEGER :: iLevel

    CALL DestroyRadiationFields

    CALL Finalize_IMEX_RK_MF

    DO iLevel = 0, nLevels-1
      CALL amrex_geometry_destroy( GEOM(iLevel) )
    END DO

    CALL FinalizeSlopeLimiter_TwoMoment_MF

    CALL FinalizePositivityLimiter_TwoMoment_MF

    CALL MyAmrFinalize

    CALL FinalizeTimers

    CALL amrex_amrcore_finalize

    CALL amrex_finalize


  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
