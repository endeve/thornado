PROGRAM main


  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule,     ONLY: & 
    MF_ComputeTimeStep
  USE MyAmrDataModule,                  ONLY: &
    MF_uCR, &
    MF_uPR
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
  
  DO WHILE( ALL( t .LT. t_end ) )
     
    CALL MF_ComputeTimeStep( dt )
    
    IF( ALL( t + dt .LE. t_end ) )THEN
      t = t + dt
    ELSE
      dt = t_end - [t]
      t  = [t_end]
    END IF

  END DO
  
  CALL FinalizeProgram( GEOM )
  

END PROGRAM main
