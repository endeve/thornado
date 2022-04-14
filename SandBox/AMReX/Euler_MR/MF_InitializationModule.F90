MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    ProgramName
  USE MF_InitializationModule_Relativistic_IDEAL, ONLY: &
    InitializeFields_Euler_Relativistic_IDEAL_MF
  USE MF_InitializationModule_AdiabaticCollapse_XCFC, ONLY: &
    InitializeFields_Euler_AdiabaticCollapse_XCFC_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Euler_MF

CONTAINS


  SUBROUTINE InitializeFields_Euler_MF( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

#ifdef HYDRO_RELATIVISTIC

    IF( TRIM( ProgramName ) .EQ. 'AdiabaticCollapse_XCFC' )THEN

      CALL InitializeFields_Euler_AdiabaticCollapse_XCFC_MF &
             ( MF_uGF, MF_uCF )

    ELSE

      CALL InitializeFields_Euler_Relativistic_IDEAL_MF &
             ( iLevel, MF_uGF, MF_uCF )

    END IF

#endif

  END SUBROUTINE InitializeFields_Euler_MF

END MODULE MF_InitializationModule
