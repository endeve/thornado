MODULE FinalizationModule

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ReferenceElementModuleX, ONLY: &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    FinalizeReferenceElementX_Lagrange
  USE MHD_MeshRefinementModule, ONLY: &
    FinalizeMeshRefinement_MHD

  ! --- Local Modules ---

  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF, &
    DestroyFields_Geometry_MF
  USE MF_FieldsModule_MHD, ONLY: &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM, &
    DestroyFields_MHD_MF
  USE MF_EquationOfStateModule_MHD, ONLY: &
    FinalizeEquationOfState_MF
  USE MF_MHD_SlopeLimiterModule, ONLY: &
    FinalizeSlopeLimiter_MHD_MF
  USE MF_MHD_PositivityLimiterModule, ONLY: &
    FinalizePositivityLimiter_MHD_MF
  USE MF_TimeSteppingModule_SSPRK_MHD, ONLY: &
    FinalizeFluid_SSPRK_MF
  USE MF_MHD_UtilitiesModule, ONLY: &
    ComputeFromConserved_MHD_MF
  USE InputOutputModuleAMReX_MHD, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_MHD_TallyModule, ONLY: &
    ComputeTally_MHD_MF, &
    FinalizeTally_MHD_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    MHDMomentumX1_Initial, &
    MHDMomentumX1_OffGrid, &
    MHDMomentumX2_Initial, &
    MHDMomentumX2_OffGrid, &
    MHDMomentumX3_Initial, &
    MHDMomentumX3_OffGrid, &
    MHDEnergy_Initial, &
    MHDEnergy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    MHDMagFieldX1_Initial, &
    MHDMagFieldX1_OffGrid, &
    MHDMagFieldX2_Initial, &
    MHDMagFieldX2_OffGrid, &
    MHDMagFieldX3_Initial, &
    MHDMagFieldX3_OffGrid, &
    MHDCleaningField_Initial, &
    MHDCleaningField_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    dt, &
    t_old, &
    t_new
  USE MF_TimersModule_MHD, ONLY: &
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

    CALL ComputeFromConserved_MHD_MF &
           ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCM_Option = MF_uCM, &
             MF_uPM_Option = MF_uPM, &
             MF_uAM_Option = MF_uAM, &
             MF_uDM_Option = MF_uDM )

    CALL WriteFieldsAMReX_Checkpoint &
           ( StepNo, nLevels, dt, t_new, &
             [ BaryonicMass_Initial   , BaryonicMass_OffGrid    ], &
             [ MHDMomentumX1_Initial, MHDMomentumX1_OffGrid ], &
             [ MHDMomentumX2_Initial, MHDMomentumX2_OffGrid ], &
             [ MHDMomentumX3_Initial, MHDMomentumX3_OffGrid ], &
             [ MHDEnergy_Initial    , MHDEnergy_OffGrid     ], &
             [ ElectronNumber_Initial , ElectronNumber_OffGrid  ], &
             [ MHDMagFieldX1_Initial, MHDMagFieldX1_OffGrid ], &
             [ MHDMagFieldX2_Initial, MHDMagFieldX2_OffGrid ], &
             [ MHDMagFieldX3_Initial, MHDMagFieldX3_OffGrid ], &
             [ MHDCleaningField_Initial, MHDCleaningField_OffGrid ], &
             [ ADMMass_Initial        , ADMMass_OffGrid         ], &
             MF_uGF % BA % P, &
             iWriteFields_uGF = 1, &
             iWriteFields_uCM = 1, &
             iWriteFields_uCR = 0, &
             pMF_uGF_Option = MF_uGF % P, &
             pMF_uCM_Option = MF_uCM % P )

    CALL ComputeTally_MHD_MF( t_new, MF_uGF, MF_uCM, MF_uAM )

    CALL FinalizeFluid_SSPRK_MF

    CALL FinalizeTally_MHD_MF

    DEALLOCATE( t_new )
    DEALLOCATE( t_old )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

    CALL FinalizeSlopeLimiter_MHD_MF

    CALL FinalizePositivityLimiter_MHD_MF

    CALL FinalizeEquationOfState_MF

    CALL FinalizeMeshRefinement_MHD

    CALL FinalizeReferenceElementX_Lagrange
    CALL FinalizeReferenceElementX

    CALL DestroyFields_MHD_MF
    CALL DestroyFields_Geometry_MF

    CALL TimersStop_AMReX( Timer_AMReX_Finalize )

    CALL FinalizeTimers_AMReX( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL amrex_amrcore_finalize()

    CALL amrex_finalize()

  END SUBROUTINE FinalizeProgram

END MODULE FinalizationModule
