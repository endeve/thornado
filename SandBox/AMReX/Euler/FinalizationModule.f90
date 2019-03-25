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
   USE MeshModule,                      ONLY: &
    MeshType, DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    FinalizeEquationOfState
  USE Euler_SlopeLimiterModule,         ONLY: &
    Euler_FinalizeSlopeLimiter
  USE Euler_PositivityLimiterModule,    ONLY: &
    Euler_FinalizePositivityLimiter

  ! --- Local Modules ---
  USE MyAmrModule,                 ONLY: &
    nLevels, MyAmrFinalize
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    MF_FinalizeFluid_SSPRK

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram( GEOM, MeshX )

    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels)
    TYPE(MeshType),        INTENT(inout) :: MeshX(1:3)

    INTEGER :: iLevel, iDim

    CALL Euler_FinalizePositivityLimiter

    CALL Euler_FinalizeSlopeLimiter

    CALL FinalizeEquationOfState

    CALL MF_FinalizeFluid_SSPRK

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    DO iDim = 1, 3
      CALL DestroyMesh( MeshX(iDim) )
    END DO

    DO iLevel = 0, nLevels
      CALL amrex_geometry_destroy( GEOM(iLevel) )
    END DO

    CALL MyAmrFinalize

    CALL amrex_amrcore_finalize

    CALL amrex_finalize

 END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
