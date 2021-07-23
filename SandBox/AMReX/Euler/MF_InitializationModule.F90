MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

  USE MF_InitializationModule_Relativistic_TABLE

#else

  USE MF_InitializationModule_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

  USE MF_InitializationModule_Relativistic_IDEAL

#else

  USE MF_InitializationModule_NonRelativistic_IDEAL

#endif

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF


CONTAINS


  SUBROUTINE InitializeFields_MF( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    CALL InitializeFields_Euler_MF( iLevel, MF_uGF, MF_uCF )

  END SUBROUTINE InitializeFields_MF


END MODULE MF_InitializationModule
