MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE Euler_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_Euler
  USE Euler_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_Euler
  USE EquationOfStateModule, ONLY: &
    FinalizeEquationOfState

  ! --- Local Modules ---

  USE MF_FieldsModule, ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    DestroyFields_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    FinalizeFluid_SSPRK_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF, &
    FinalizeTally_Euler_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    dt, &
    t_old, &
    t_new, &
    lo_bc, &
    hi_bc
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
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCF_Option = MF_uCF, &
             MF_uPF_Option = MF_uPF, &
             MF_uAF_Option = MF_uAF, &
             MF_uDF_Option = MF_uDF )

    CALL WriteFieldsAMReX_Checkpoint &
           ( StepNo, nLevels, dt, t_new, &
             MF_uGF % BA % P, &
             MF_uGF % P, &
             MF_uCF % P )

    CALL ComputeTally_Euler_MF( t_new, MF_uGF, MF_uCF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeTally_Euler_MF

    CALL FinalizeSlopeLimiter_Euler

    CALL FinalizePositivityLimiter_Euler

    CALL FinalizeEquationOfState

    CALL FinalizeReferenceElementX

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL DestroyFields_MF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeTimers_AMReX_Euler( WriteAtIntermediateTime_Option = .TRUE. )

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
