MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
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
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_Euler_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    FinalizeFluid_SSPRK_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF, &
    FinalizeTally_Euler_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    Energy_Initial, &
    Energy_OffGrid, &
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
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    FinalizeGravitySolver_XCFC_MF
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_Finalize, &
    FinalizeTimers_AMReX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE FinalizeProgram

    CALL TimersStart_AMReX( Timer_AMReX_Finalize )

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
             [ BaryonicMass_Initial  , BaryonicMass_OffGrid   ], &
             [ Energy_Initial        , Energy_OffGrid         ], &
             [ ElectronNumber_Initial, ElectronNumber_OffGrid ], &
             [ ADMMass_Initial       , ADMMass_OffGrid        ], &
             MF_uGF % BA % P, &
             iWriteFields_uGF = 1, &
             iWriteFields_uCF = 1, &
             iWriteFields_uCR = 0, &
             pMF_uGF_Option = MF_uGF % P, &
             pMF_uCF_Option = MF_uCF % P )

    CALL ComputeTally_Euler_MF( t_new, MF_uGF, MF_uCF )

    CALL FinalizeFluid_SSPRK_MF

    CALL FinalizeGravitySolver_XCFC_MF

    CALL FinalizeTally_Euler_MF

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL FinalizeSlopeLimiter_Euler_MF

    CALL FinalizePositivityLimiter_Euler_MF

    CALL FinalizeEquationOfState

    CALL FinalizeMeshRefinement_Euler

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    CALL DestroyFields_Euler_MF
    CALL DestroyFields_Geometry_MF

    CALL TimersStop_AMReX( Timer_AMReX_Finalize )

    CALL FinalizeTimers_AMReX

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
