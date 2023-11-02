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
  USE TwoMoment_OpacityModule, ONLY: &
    DestroyOpacities
  USE TwoMoment_TimersModule, ONLY: &
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
    MF_uGR, &
    DestroyFields_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_TwoMoment_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_TimeSteppingModule_Relativistic, ONLY: &
    Finalize_IMEX_RK_MF
  USE MF_TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment_MF, &
    ComputeGray_TwoMoment_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    FinalizeTally_Euler_MF
  USE MF_TwoMoment_TallyModule, ONLY: &
    FinalizeTally_TwoMoment_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    dt, &
    t_old, &
    t_new

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE FinalizeProgram



    CALL Finalize_IMEX_RK_MF

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL FinalizeTally_TwoMoment_MF

    CALL FinalizeTally_Euler_MF

    CALL FinalizeSlopeLimiter_TwoMoment_MF

    CALL FinalizePositivityLimiter_TwoMoment_MF

    !CALL DestroyOpacities

    CALL FinalizeEquationOfState

    CALL DestroyGeometryFieldsE

    CALL FinalizeReferenceElement_Lagrange
    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElementE_Lagrange
    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    CALL DestroyFields_TwoMoment_MF
    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL FinalizeTimers

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
