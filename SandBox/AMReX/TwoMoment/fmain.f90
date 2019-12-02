PROGRAM main

  USE InitializationModule,             ONLY: &
    InitializeProgram
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MyAmrModule,                      ONLY: &
    nLevels,   &
    StepNo,    &
    t,         &
    dt,        &
    t_end,     &
    CFL,       &
    t_wrt,     &
    dt_wrt,    &
    t_chk,     &
    dt_chk,    &
    iCycleD,   &
    iCycleW,   &
    iCycleChk, &
    GEOM


  IMPLICIT NONE

  CALL InitializeProgram
  print*, 'Yay'
  CALL FinalizeProgram( GEOM )


END PROGRAM main
