MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE Euler_MeshRefinementModule, ONLY: &
    FinalizeMeshRefinement_Euler

  ! --- Local Modules ---

  USE MF_FieldsModule_Geometry, ONLY: &
    DestroyFields_Geometry_MF
  USE MF_FieldsModule_Euler, ONLY: &
    DestroyFields_Euler_MF
  USE InputParsingModule, ONLY: &
    StepNo, &
    t_old, &
    t_new

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE FinalizeProgram

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( StepNo )

    CALL FinalizeMeshRefinement_Euler

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
