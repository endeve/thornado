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
    WriteFieldsAMReX_PlotFile
  USE InputParsingModule, ONLY: &
    StepNo, &
    dt, &
    t_old, &
    t_new, &
    lo_bc, &
    hi_bc, &
    lo_bc_uCF, &
    hi_bc_uCF
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Finalize, &
    Timer_AMReX_Euler_InputOutput, &
    FinalizeTimers_AMReX_Euler

use mf_utilitiesmodule,only:showvariablefrommultifab,filename

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram


CONTAINS


  SUBROUTINE FinalizeProgram

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

write(filename,'(A)') 'nN03_nX032_UniGrid_F.dat'
!write(filename,'(A)') 'nN03_nX032_MultiGrid_F_RefluxOff.dat'
!write(filename,'(A)') 'nN03_nX032_MultiGrid_F_RefluxOn.dat'
open(100,file=trim(filename))
close(100)
call showvariablefrommultifab(mf_ucf,1,writetofile_option=.true.)

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCF_Option = MF_uCF, &
             MF_uPF_Option = MF_uPF, &
             MF_uAF_Option = MF_uAF, &
             MF_uDF_Option = MF_uDF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeSlopeLimiter_Euler

    CALL FinalizePositivityLimiter_Euler

    CALL FinalizeEquationOfState

    CALL FinalizeReferenceElementX

    DEALLOCATE( hi_bc_uCF )
    DEALLOCATE( lo_bc_uCF )
    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

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
