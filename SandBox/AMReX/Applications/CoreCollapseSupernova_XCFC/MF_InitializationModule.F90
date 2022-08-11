MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    ProgramName
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE MF_InitializationModule_CoreCollapseSupernova_XCFC, ONLY: &
    InitializeFields_CoreCollapseSupernova_XCFC_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF( iLevel, MF_uGF, MF_uCF, MF_uCR )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF, MF_uCF, MF_uCR

#ifdef HYDRO_RELATIVISTIC

    IF( TRIM( ProgramName ) .EQ. 'CoreCollapseSupernova_XCFC' )THEN

      CALL InitializeFields_CoreCollapseSupernova_XCFC_MF &
             ( iLevel, MF_uGF, MF_uCF, MF_uCR )

    ELSE

      CALL DescribeError_Euler_MF( 99, 'ProgramName not found' )

    END IF

#endif

  END SUBROUTINE InitializeFields_MF

END MODULE MF_InitializationModule
