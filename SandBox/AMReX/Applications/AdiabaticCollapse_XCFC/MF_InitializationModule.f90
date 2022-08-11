MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE MF_InitializationModule_AdiabaticCollapse_XCFC, ONLY: &
    InitializeFields_AdiabaticCollapse_XCFC_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF, MF_uCF

    CALL InitializeFields_AdiabaticCollapse_XCFC_MF &
           ( iLevel, MF_uGF, MF_uCF )

  END SUBROUTINE InitializeFields_MF


END MODULE MF_InitializationModule
