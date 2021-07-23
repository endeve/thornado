MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX

  ! --- Local Modules ---

  USE MF_FieldsModule, ONLY: &
    DestroyFields_MF
  USE InputParsingModule, ONLY: &
    StepNo, &
    dt, &
    t_old, &
    t_new

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram

    CALL FinalizeReferenceElementX

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL DestroyFields_MF

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
