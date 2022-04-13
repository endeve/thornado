MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module,                ONLY: &
    amrex_finalize
  USE amrex_geometry_module,            ONLY: &
    amrex_geometry, &
    amrex_geometry_destroy
  USE amrex_amrcore_module,             ONLY: &
    amrex_amrcore_finalize
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ReferenceElementModuleX,          ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshX,    &
    DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    FinalizeEquationOfState
  USE Euler_SlopeLimiterModule,         ONLY: &
    FinalizeSlopeLimiter_Euler
  USE Euler_PositivityLimiterModule,    ONLY: &
    FinalizePositivityLimiter_Euler

  ! --- Local Modules ---

  USE InputParsingModule,               ONLY: &
    nLevels, &
    FinalizeParameters
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_FinalizeFluid_SSPRK
  USE MF_Euler_TallyModule,             ONLY: &
    MF_FinalizeTally_Euler
  USE TimersModule_AMReX_Euler,         ONLY: &
    FinalizeTimers_AMReX_Euler,    &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler,  &
    Timer_AMReX_Euler_Finalize

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram( GEOM )

    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)

    INTEGER :: iLevel, iDim

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL MF_FinalizeTally_Euler

    CALL MF_FinalizeFluid_SSPRK

    CALL FinalizePositivityLimiter_Euler

    CALL FinalizeSlopeLimiter_Euler

    CALL FinalizeEquationOfState

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_geometry_destroy( GEOM(iLevel) )

    END DO

    CALL FinalizeParameters

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeTimers_AMReX_Euler

    CALL amrex_amrcore_finalize

    CALL amrex_finalize

 END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
