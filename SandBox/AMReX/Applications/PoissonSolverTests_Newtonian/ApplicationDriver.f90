PROGRAM main

  ! --- AMReX Modules ---

  ! --- thornado Modules ---

  ! --- Local Modules ---

  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler
  USE MF_TimersModule, ONLY: &
    TimeIt_AMReX

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TimeIt_AMReX       = .TRUE.
  TimeIt_AMReX_Euler = .TRUE.

  CALL InitializeProgram

  CALL FinalizeProgram

END PROGRAM main
