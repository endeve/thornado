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

  ! --- Local Modules ---

  USE InputParsingModule,               ONLY: &
    nLevels, &
    FinalizeParameters
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_FinalizeMagnetofluid_SSPRK

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram( GEOM )

    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)

    INTEGER :: iLevel, iDim

    CALL MF_FinalizeMagnetofluid_SSPRK

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

    CALL amrex_amrcore_finalize

    CALL amrex_finalize

 END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
