MODULE FinalizationModule

  ! --- AMReX Modules ---
  USE amrex_amr_module,     ONLY: &
    amrex_geometry, &
    amrex_geometry_destroy, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---
  USE ReferenceElementModuleX,          ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshType, DestroyMesh
  USE FluidFieldsModule,                ONLY: &
    DestroyFluidFields
  USE EquationOfStateModule,            ONLY: &
    FinalizeEquationOfState
  USE Euler_SlopeLimiterModule,         ONLY: &
    FinalizeSlopeLimiter_Euler
  USE Euler_PositivityLimiterModule,    ONLY: &
    FinalizePositivityLimiter_Euler

  ! --- Local Modules ---
  USE MyAmrModule,                 ONLY: &
    nLevels, MyAmrFinalize
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    MF_FinalizeFluid_SSPRK
  USE TimersModule_AMReX_Euler,    ONLY: &
    TimersStart_AMReX_Euler, TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Finalize

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram( GEOM, MeshX )

    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)
    TYPE(MeshType),        INTENT(inout) :: MeshX(1:3)

    INTEGER :: iLevel, iDim

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL MF_FinalizeFluid_SSPRK

    CALL DestroyFluidFields

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

    CALL MyAmrFinalize

    CALL amrex_amrcore_finalize

    CALL amrex_finalize

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Finalize )

 END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
