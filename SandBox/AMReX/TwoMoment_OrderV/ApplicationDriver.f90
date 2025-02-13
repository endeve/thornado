PROGRAM main

  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule,     ONLY: &
    ComputeFromConserved_TwoMoment_MF
  USE MF_FieldsModule_Geometry,                  ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler,                  ONLY: &
    MF_uCF, &
    MF_uPF
  USE MF_FieldsModule_TwoMoment,                  ONLY: &
    MF_uPR, &
    MF_uCR
  USE InitializationModule,             ONLY: &
    InitializeProgram
  USE FinalizationModule,               ONLY: &
    FinalizeProgram


  IMPLICIT NONE

  CALL InitializeProgram

  CALL ComputeFromConserved_TwoMoment_MF( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

  CALL FinalizeProgram


END PROGRAM main
