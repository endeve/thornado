PROGRAM ApplicationDriver

  ! --- Local Modules ---

  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CALL InitializeProgram

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
