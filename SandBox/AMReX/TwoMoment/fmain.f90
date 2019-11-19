PROGRAM main

  USE InitializationModule,             ONLY: &
    InitializeProgram
  USE FinalizationModule,               ONLY: &
    FinalizeProgram



  IMPLICIT NONE

  CALL InitializeProgram

  CALL FinalizeProgram


END PROGRAM main
