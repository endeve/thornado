MODULE FinalizationModule

  ! --- thornado Modules ---
  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshType, DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    FinalizeEquationOfState
  USE SlopeLimiterModule_Euler, ONLY: &
    FinalizeSlopeLimiter_Euler

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    MyAmrFinalize

  ! --- AMReX Modules ---
  USE amrex_amr_module
  USE amrex_amrcore_module

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram( nLevels, GEOM, MeshX )

    INTEGER,               INTENT(in) :: nLevels
    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels)
    TYPE(MeshType),        INTENT(inout) :: MeshX(1:3)

    INTEGER :: iLevel, iDim

    DO iLevel = 0, nLevels
      CALL amrex_geometry_destroy( GEOM(iLevel) )
    END DO

    DO iDim = 1, 3
      CALL DestroyMesh( MeshX(iDim) )
    END DO

    CALL FinalizeSlopeLimiter_Euler

    CALL FinalizeEquationOfState

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementX

    CALL MyAmrFinalize

    CALL amrex_amrcore_finalize

    CALL amrex_finalize

 END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
