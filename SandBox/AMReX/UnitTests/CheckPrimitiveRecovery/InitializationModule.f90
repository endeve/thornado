MODULE InitializationModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_init
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init

  ! --- Local Modules ---

  USE MF_EquationOfStateModule, ONLY: &
    InitializeEquationOfState_MF
  USE InputParsingModule, ONLY: &
    InitializeParameters

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeParameters

    CALL InitializeEquationOfState_MF

  END SUBROUTINE InitializeProgram


END MODULE InitializationModule
