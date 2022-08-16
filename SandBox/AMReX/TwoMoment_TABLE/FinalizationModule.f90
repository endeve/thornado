MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    FinalizeReferenceElementE
  USE ReferenceElementModule_Lagrange, ONLY: &
    FinalizeReferenceElement_Lagrange
  USE ReferenceElementModule, ONLY: &
    FinalizeReferenceElement
  USE GeometryFieldsModuleE, ONLY: &
    DestroyGeometryFieldsE
  USE EquationOfStateModule, ONLY: &
    FinalizeEquationOfState
  USE TwoMoment_SlopeLimiterModule_Relativistic, ONLY: &
    FinalizeSlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_Relativistic, ONLY: &
    FinalizePositivityLimiter_TwoMoment
  USE TwoMoment_OpacityModule_Relativistic, ONLY: &
    DestroyOpacities
  USE TwoMoment_TimersModule_Relativistic, ONLY: &
    FinalizeTimers

  ! --- Local Modules ---

  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF, &
    DestroyFields_Geometry_MF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    DestroyFields_Euler_MF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR, &
    MF_uPR, &
    DestroyFields_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_TwoMoment_TimeSteppingModule_Relativistic, ONLY: &
    MF_FinalizeField_IMEX_RK
  USE MF_TwoMoment_UtilitiesModule, ONLY: &
    MF_ComputeFromConserved
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    FinalizeTally_Euler_MF
  USE MF_TwoMoment_TallyModule, ONLY: &
    MF_FinalizeTally_TwoMoment
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    dt, &
    t_old, &
    t_new, &
    lo_bc, &
    hi_bc

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE FinalizeProgram

    CALL MF_ComputeFromConserved &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

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
             iWriteFields_uGF = 1, &
             iWriteFields_uCF = 1, &
             iWriteFields_uCR = 0, &
             pMF_uGF_Option = MF_uGF % P, &
             pMF_uCF_Option = MF_uCF % P )

    CALL MF_FinalizeField_IMEX_RK

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL MF_FinalizeTally_TwoMoment

    CALL FinalizeTally_Euler_MF

    CALL FinalizeSlopeLimiter_TwoMoment

    CALL FinalizePositivityLimiter_TwoMoment

    CALL DestroyOpacities

    CALL FinalizeEquationOfState

    CALL DestroyGeometryFieldsE

    CALL FinalizeReferenceElement_Lagrange
    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElementE_Lagrange
    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

    CALL DestroyFields_TwoMoment_MF
    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL FinalizeTimers

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
