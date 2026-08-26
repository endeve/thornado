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
  USE TwoMoment_TimersModule, ONLY: &
    FinalizeTimers

  USE Euler_MeshRefinementModule, ONLY: &
    FinalizeMeshRefinement_Euler_Aniso

  !USE AnisotropicRefinementModule, ONLY: &
  !  UseAnisotropicRefinement, &
  !  FinalizeAnisotropicRefinement

  ! --- Local Modules ---

  USE MF_EquationOfStateModule, ONLY: &
    FinalizeEquationOfState_MF
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
    MF_uAR, &
    MF_uGR, &
    DestroyFields_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  !USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    !FinalizeSlopeLimiter_TwoMoment_MF
  !USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    !FinalizePositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeFromConserved_TwoMoment_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    FinalizeTally_Euler_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    EulerMomentumX1_Initial, &
    EulerMomentumX1_OffGrid, &
    EulerMomentumX2_Initial, &
    EulerMomentumX2_OffGrid, &
    EulerMomentumX3_Initial, &
    EulerMomentumX3_OffGrid, &
    EulerEnergy_Initial, &
    EulerEnergy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    dt, &
    t_old, &
    t_new
  USE MF_TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Finalize_IMEX_RK_MF
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab
  USE MF_TwoMoment_TallyModule

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE FinalizeProgram


    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL ComputeFromConserved_TwoMoment_MF &
           (  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )
    
    StepNo = StepNo + 1

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCF_Option = MF_uCF, &
             MF_uPF_Option = MF_uPF, &
             MF_uAF_Option = MF_uAF, &
             MF_uDF_Option = MF_uDF, &
             MF_uPR_Option = MF_uPR, &
             MF_uCR_Option = MF_uCR, &
             MF_uGR_Option = MF_uGR )

    CALL WriteFieldsAMReX_Checkpoint

    CALL ShowVariableFromMultiFab( MF_uPR, 1, &
                             WriteToFile_Option = .TRUE., &
                             FileNameBase_Option = 'MF_uPR' )

    CALL ShowVariableFromMultiFab( MF_uCR, 1, &
                             WriteToFile_Option = .TRUE., &
                             FileNameBase_Option = 'MF_uCR' )
    CALL FinalizeTally_TwoMoment_MF

    CALL Finalize_IMEX_RK_MF

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    !CALL FinalizeSlopeLimiter_TwoMoment_MF

    !CALL FinalizePositivityLimiter_TwoMoment_MF

    CALL FinalizeEquationOfState_MF

    CALL DestroyGeometryFieldsE

    CALL FinalizeReferenceElement_Lagrange
    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElementE_Lagrange
    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    !IF( UseAnisotropicRefinement ) &
    !  CALL FinalizeMeshRefinement_Euler_Aniso

    !CALL FinalizeAnisotropicRefinement

    CALL DestroyFields_TwoMoment_MF
    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL FinalizeTimers

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
