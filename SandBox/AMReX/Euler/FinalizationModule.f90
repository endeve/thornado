MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE EquationOfStateModule, ONLY: &
    FinalizeEquationOfState

  ! --- Local Modules ---

  USE MF_FieldsModule, ONLY: &
    MF_uGF_new, &
    MF_uCF_new, &
    MF_uPF_new, &
    MF_uAF_new, &
    MF_uDF_new, &
    DestroyFields_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    FinalizeFluid_SSPRK_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile
  USE InputParsingModule, ONLY: &
    StepNo, &
    dt, &
    t_old, &
    t_new
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Finalize, &
    Timer_AMReX_Euler_InputOutput, &
    FinalizeTimers_AMReX_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF_new, MF_uCF_new, MF_uPF_new, MF_uAF_new )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF_new, &
             MF_uGF_Option = MF_uGF_new, &
             MF_uCF_Option = MF_uCF_new, &
             MF_uPF_Option = MF_uPF_new, &
             MF_uAF_Option = MF_uAF_new, &
             MF_uDF_Option = MF_uDF_new )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeEquationOfState

    CALL FinalizeReferenceElementX

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL DestroyFields_MF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeTimers_AMReX_Euler

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
