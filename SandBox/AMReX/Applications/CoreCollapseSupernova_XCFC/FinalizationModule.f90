MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleZ, ONLY: &
    FinalizeReferenceElementZ
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE MeshModule, ONLY: &
    MeshE
  USE EquationOfStateModule, ONLY: &
    FinalizeEquationOfState
  USE Euler_MeshRefinementModule, ONLY: &
    FinalizeMeshRefinement_Euler

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
    DestroyFields_TwoMoment_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_Euler_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_TwoMoment_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_TwoMoment_MF
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
  USE MF_GravitySolutionModule_XCFC_Poseidon, ONLY: &
    FinalizeGravitySolver_XCFC_Poseidon_MF

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
             iWriteFields_uGF = 1, &
             iWriteFields_uCF = 1, &
             iWriteFields_uCR = 0, &
             pMF_uGF_Option = MF_uGF % P, &
             pMF_uCF_Option = MF_uCF % P )

    CALL ComputeTally_Euler_MF( t_new, MF_uGF, MF_uCF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL Finalize_IMEX_RK_MF

    CALL FinalizeGravitySolver_XCFC_Poseidon_MF

    CALL FinalizeTally_TwoMoment_MF

    CALL FinalizeTally_Euler_MF

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL FinalizeSlopeLimiter_TwoMoment_MF

    CALL FinalizePositivityLimiter_TwoMoment_MF

    CALL FinalizeSlopeLimiter_Euler_MF

    CALL FinalizePositivityLimiter_Euler_MF

    CALL FinalizeEquationOfState

    CALL DestroyGeometryFieldsE

    CALL FinalizeMeshRefinement_Euler

    CALL FinalizeReferenceElement_Lagrange
    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElementE_Lagrange
    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    CALL DestroyMesh( MeshE )

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

    CALL DestroyFields_TwoMoment_MF
    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Finalize )

    CALL FinalizeTimers

    CALL FinalizeTimers_AMReX_Euler( WriteAtIntermediateTime_Option = .TRUE. )

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
