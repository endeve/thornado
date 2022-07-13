MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE MF_InitializationModule_Relativistic_IDEAL, ONLY: &
    InitializeFields_Euler_Relativistic_IDEAL_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Euler_MF

CONTAINS


  SUBROUTINE InitializeFields_Euler_MF( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

#ifdef HYDRO_RELATIVISTIC

    CALL InitializeFields_Euler_Relativistic_IDEAL_MF &
           ( iLevel, MF_uGF, MF_uCF )

#endif

  END SUBROUTINE InitializeFields_Euler_MF

END MODULE MF_InitializationModule
