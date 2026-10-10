MODULE TwoMoment_TimeSteppingModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    iE_B0, iE_B1, iE_E0, iE_E1
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, uDF
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    OffGridFlux_Euler_X1_Inner, &
    OffGridFlux_Euler_X1_Outer, &
    OffGridFlux_Euler_X2_Inner, &
    OffGridFlux_Euler_X2_Outer, &
    OffGridFlux_Euler_X3_Inner, &
    OffGridFlux_Euler_X3_Outer
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies, &
    uDR
#if defined( TWOMOMENT_ORDER_V )
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    nQ, nC, nLeptonFamilies, nLeptonFlavorCheck
#endif
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_IMEX, &
    Timer_TimeStepper, &
    Timer_Euler, &
    Timer_Poisson
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment, &
    dEnergyMomentum_PL_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    SourceAwareImplicitOpacityAvailable, &
    OffGridFlux_TwoMoment, OffGridLeptonFlux_TwoMoment
  USE TwoMoment_TallyModule, ONLY: &
    IncrementOffGridTally_Euler, &
    IncrementOffGridTally_TwoMoment, &
    IncrementPositivityLimiterTally_TwoMoment

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP)              :: OffGridFlux_U(nCF)
    REAL(DP)              :: OffGridFlux_M(2*nCR)
    REAL(DP)              :: OffGridLeptonFlux_M(3)
    REAL(DP), ALLOCATABLE :: dU_IM(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_EX(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_EX(:,:,:,:,:,:,:)
  END TYPE StageDataType

  LOGICAL                          :: EvolveEuler
  LOGICAL                          :: EvolveTwoMoment
  INTEGER                          :: nStages
  REAL(DP),            ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP),            ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP),            ALLOCATABLE :: U0(:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Ui(:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: M0(:,:,:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Mi(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)

#if defined( TWOMOMENT_ORDER_V )
  INTEGER, PARAMETER, PUBLIC :: iIMEX_Assembled             = 1
  INTEGER, PARAMETER, PUBLIC :: iIMEX_PreImplicit           = 2
  INTEGER, PARAMETER, PUBLIC :: iIMEX_PostImplicit          = 3
  INTEGER, PARAMETER, PUBLIC :: iIMEX_PostImplicitLimiter   = 4
  INTEGER, PARAMETER, PUBLIC :: iIMEX_CompletedStep          = 5
  INTEGER, PARAMETER, PUBLIC :: nIMEXStageStates            = 5
#endif

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK
#if defined( TWOMOMENT_ORDER_V )
  PUBLIC :: GetNumberOfStages_IMEX_RK
  PUBLIC :: GetImplicitStageCoefficient_IMEX_RK
  PUBLIC :: GetExplicitStageCoefficient_IMEX_RK
  PUBLIC :: GetExplicitStageWeight_IMEX_RK
  PUBLIC :: GetImplicitStageWeight_IMEX_RK
  PUBLIC :: ExplicitRHSUsed_IMEX_RK
  PUBLIC :: HasFinalAssembly_IMEX_RK
#endif

  INTERFACE
    SUBROUTINE ImplicitIncrement &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R, &
        AccumulateMatterCoupling_Mask, uDR_Option, MatterCoupling_Mask_Option &
#if defined( TWOMOMENT_ORDER_V )
        , MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
          IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
          LeptonFlavorCheck_Rate_Option &
#endif
        )
      USE KindModule           , ONLY: DP
      USE ProgramHeaderModule  , ONLY: nDOFX, nDOFE, nDOFZ
      USE GeometryFieldsModuleE, ONLY: nGE
      USE GeometryFieldsModule , ONLY: nGF
      USE FluidFieldsModule    , ONLY: nCF
      USE RadiationFieldsModule, ONLY: nCR, nSpecies, nDR
#if defined( TWOMOMENT_ORDER_V )
      USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
        nQ, nC, nLeptonFamilies, nLeptonFlavorCheck
#endif
      INTEGER,  INTENT(in)    :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        dt
      REAL(DP), INTENT(in)    :: &
        GE  (1:nDOFE, &
             iZ_B1(1):iZ_E1(1), &
             1:nGE)
      REAL(DP), INTENT(in)    :: &
        GX  (1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nGF)
      REAL(DP), INTENT(in) :: &
        U_F (1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCF)
      REAL(DP), INTENT(out) :: &
        dU_F(1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCF)
      REAL(DP), INTENT(in) :: &
        U_R (1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCR, &
             1:nSpecies)
      REAL(DP), INTENT(out) :: &
        dU_R(1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCR, &
             1:nSpecies)
      LOGICAL, INTENT(in) :: &
        AccumulateMatterCoupling_Mask
      REAL(DP), INTENT(inout), OPTIONAL :: &
        uDR_Option &
            (iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nDR)
      LOGICAL, INTENT(inout), OPTIONAL :: &
        MatterCoupling_Mask_Option &
            (iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4))
#if defined( TWOMOMENT_ORDER_V )
      REAL(DP), INTENT(inout), OPTIONAL :: &
        MatterCollisionExchange_Rate_Option &
            (1:nQ, &
             1:1+nSpecies+nC, &
             1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4))
      REAL(DP), INTENT(inout), OPTIONAL :: &
        MatterUpdate_Rate_Option &
            (1:nQ, &
             1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4))
      REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
      REAL(DP), INTENT(out), OPTIONAL :: &
        MatterCollisionRHS_Option &
            (1:nQ, &
             1:1+nSpecies+nC, &
             1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4))
      REAL(DP), INTENT(out), OPTIONAL :: &
        LeptonFlavorCheck_Rate_Option &
            (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
             iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4))
#endif
    END SUBROUTINE ImplicitIncrement
  END INTERFACE

  INTERFACE
    SUBROUTINE GravitySolver( iX_B0, iX_E0, iX_B1, iX_E1, GX, D )
      USE KindModule          , ONLY: DP
      USE ProgramHeaderModule , ONLY: nDOFX
      USE GeometryFieldsModule, ONLY: nGF
      INTEGER,  INTENT(in)    :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(inout) :: &
        GX(1:, &
           iX_B1(1):, &
           iX_B1(2):, &
           iX_B1(3):, &
           1:)
      REAL(DP), INTENT(in)    :: &
        D (1:, &
           iX_B1(1):, &
           iX_B1(2):, &
           iX_B1(3):)
    END SUBROUTINE GravitySolver
  END INTERFACE

  INTERFACE AllocateArray
    MODULE PROCEDURE AllocateArray5D
    MODULE PROCEDURE AllocateArray7D
  END INTERFACE AllocateArray

  INTERFACE DeallocateArray
    MODULE PROCEDURE DeallocateArray5D
    MODULE PROCEDURE DeallocateArray7D
  END INTERFACE DeallocateArray

  INTERFACE CopyArray
    MODULE PROCEDURE CopyArray5D
    MODULE PROCEDURE CopyArray7D
  END INTERFACE CopyArray

  INTERFACE AddToArray
    MODULE PROCEDURE AddToArray5D
    MODULE PROCEDURE AddToArray7D
  END INTERFACE AddToArray

CONTAINS


  SUBROUTINE Update_IMEX_RK &
    ( dt, GE, GX, U, M, ComputeIncrement_TwoMoment_Implicit, SolveGravity &
#if defined( TWOMOMENT_ORDER_V )
      , MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
        MatterCollisionRHS_FinalStage_Option, &
        RadiationState_BeforeFinalLimiter_Option, &
        IMEX_Stage_FluidState_Option, IMEX_Stage_RadiationState_Option, &
        IMEX_Stage_FluidImplicitRate_Option, &
        IMEX_Stage_FluidExplicitRate_Option, &
        IMEX_Stage_RadiationImplicitRate_Option, &
        IMEX_Stage_RadiationExplicitRate_Option, &
        IMEX_Stage_LeptonFlavorCheck_Rate_Option, &
        UpdateTallies_Option, OffGridFluidIncrement_Option, &
        OffGridRadiationIncrement_Option, &
        OffGridRadiationLeptonFlavorIncrement_Option, &
        ComputeIncrement_TwoMoment_Explicit_Option &
#endif
      )

    REAL(DP), INTENT(in) :: &
      dt
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iE_B1:iE_E1, &
         1:nGE)
    REAL(DP), INTENT(inout) :: &
      GX(1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3), &
         1:nGF)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(inout) :: &
      M(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    PROCEDURE(ImplicitIncrement) :: &
      ComputeIncrement_TwoMoment_Implicit
    PROCEDURE(GravitySolver), OPTIONAL :: &
      SolveGravity
#if defined( TWOMOMENT_ORDER_V )
    PROCEDURE(ComputeIncrement_TwoMoment_Explicit), OPTIONAL :: &
      ComputeIncrement_TwoMoment_Explicit_Option
    ! Stage-weighted collision exchange over the completed IMEX timestep.
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
        (1:nQ, &
         1:1+nSpecies+nC, &
         1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3))
    ! Stage-weighted complete matter update returned by the implicit calls.
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
        (1:nQ, &
         1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3))
    ! Unweighted collision RHS evaluated in the final implicit stage.
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterCollisionRHS_FinalStage_Option &
        (1:nQ, &
         1:1+nSpecies+nC, &
         1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(out), OPTIONAL :: &
      RadiationState_BeforeFinalLimiter_Option &
        (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)
    REAL(DP), INTENT(out), OPTIONAL :: &
      IMEX_Stage_FluidState_Option(:,:,:,:,:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: &
      IMEX_Stage_RadiationState_Option(:,:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      IMEX_Stage_FluidImplicitRate_Option(:,:,:,:,:,:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      IMEX_Stage_FluidExplicitRate_Option(:,:,:,:,:,:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      IMEX_Stage_RadiationImplicitRate_Option(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      IMEX_Stage_RadiationExplicitRate_Option(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: &
      IMEX_Stage_LeptonFlavorCheck_Rate_Option(:,:,:,:,:,:,:)
    LOGICAL, INTENT(in), OPTIONAL :: UpdateTallies_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      OffGridFluidIncrement_Option(1:nCF)
    REAL(DP), INTENT(out), OPTIONAL :: &
      OffGridRadiationIncrement_Option(1:2*nCR)
    REAL(DP), INTENT(out), OPTIONAL :: &
      OffGridRadiationLeptonFlavorIncrement_Option(3)
#endif

    INTEGER  :: iS, jS
    REAL(DP) :: dU_OffGrid(nCF)
    REAL(DP) :: dM_OffGrid(2*nCR)
    REAL(DP) :: dM_LeptonFlavor_OffGrid(3)
    REAL(DP) :: dM_PL(nCR)
    LOGICAL  :: FinalAssembly
#if defined( TWOMOMENT_ORDER_V )
    LOGICAL :: AccumulateMatterDiagnostics, CaptureIMEXStages, UpdateTallies
#endif

    CALL TimersStart( Timer_IMEX )

    dU_OffGrid = Zero
    dM_OffGrid = Zero
    dM_LeptonFlavor_OffGrid = Zero
    dM_PL = Zero
#if defined( TWOMOMENT_ORDER_V )
    FinalAssembly = HasFinalAssembly_IMEX_RK()
    AccumulateMatterDiagnostics = .FALSE.
    CaptureIMEXStages &
      =      PRESENT( IMEX_Stage_FluidState_Option ) &
        .OR. PRESENT( IMEX_Stage_RadiationState_Option )
    UpdateTallies = .TRUE.
    IF( PRESENT( UpdateTallies_Option ) ) &
      UpdateTallies = UpdateTallies_Option
#else
    FinalAssembly &
      =      ANY( a_IM(nStages,:) .NE. w_IM(:) ) &
        .OR. ANY( a_EX(nStages,:) .NE. w_EX(:) )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GE, GX, U, M )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( GE, GX, U, M )
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: MatterCollisionExchange_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( MatterCollisionExchange_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterUpdate_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: MatterUpdate_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( MatterUpdate_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterCollisionRHS_FinalStage_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: MatterCollisionRHS_FinalStage_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( MatterCollisionRHS_FinalStage_Option )
#endif
    END IF

    IF( PRESENT( RadiationState_BeforeFinalLimiter_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: RadiationState_BeforeFinalLimiter_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( RadiationState_BeforeFinalLimiter_Option )
#endif
      CALL CopyArray( RadiationState_BeforeFinalLimiter_Option, One, M )
    END IF

    IF( PRESENT( IMEX_Stage_FluidState_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: IMEX_Stage_FluidState_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( IMEX_Stage_FluidState_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationState_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( alloc: IMEX_Stage_RadiationState_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( IMEX_Stage_RadiationState_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_FluidImplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( to: IMEX_Stage_FluidImplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA COPYIN( IMEX_Stage_FluidImplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_FluidExplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( to: IMEX_Stage_FluidExplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA COPYIN( IMEX_Stage_FluidExplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationImplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( to: IMEX_Stage_RadiationImplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA COPYIN( IMEX_Stage_RadiationImplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationExplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA MAP( to: IMEX_Stage_RadiationExplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA COPYIN( IMEX_Stage_RadiationExplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_LeptonFlavorCheck_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: IMEX_Stage_LeptonFlavorCheck_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA CREATE( IMEX_Stage_LeptonFlavorCheck_Rate_Option )
#endif
    END IF
#endif

    CALL CopyArray( U0, One, U )
    CALL CopyArray( M0, One, M )

    DO iS = 1, nStages

      SourceAwareImplicitOpacityAvailable = .FALSE.

      CALL CopyArray( Ui, One, U0 )
      CALL CopyArray( Mi, One, M0 )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * a_IM(iS,jS), StageData(jS) % dU_IM )
          CALL AddToArray( One, Mi, dt * a_IM(iS,jS), StageData(jS) % dM_IM )

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * a_EX(iS,jS), StageData(jS) % dU_EX )
          CALL AddToArray( One, Mi, dt * a_EX(iS,jS), StageData(jS) % dM_EX )

        END IF

        IF( jS == iS - 1 )THEN

#if defined( TWOMOMENT_ORDER_V )
          IF( CaptureIMEXStages ) &
            CALL CaptureIMEXStageState &
                   ( Ui, Mi, iS, iIMEX_Assembled, &
                     IMEX_Stage_FluidState_Option, &
                     IMEX_Stage_RadiationState_Option )
#endif

          ! --- Apply Limiters ---

          IF( EvolveEuler )THEN

            CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

            CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

            IF( PRESENT( SolveGravity ) )THEN

              CALL TimersStart( Timer_Poisson )

#if   defined( THORNADO_OMP_OL )
              !$OMP TARGET UPDATE FROM( Ui )
#elif defined( THORNADO_OACC   )
              !$ACC UPDATE HOST       ( Ui )
#endif

              CALL SolveGravity &
                     ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

#if   defined( THORNADO_OMP_OL )
              !$OMP TARGET UPDATE TO( GX )
#elif defined( THORNADO_OACC   )
              !$ACC UPDATE DEVICE   ( GX )
#endif

              CALL TimersStop( Timer_Poisson )

            END IF

          END IF

          IF( EvolveTwoMoment )THEN

            CALL ApplySlopeLimiter_TwoMoment &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

            CALL ApplyPositivityLimiter_TwoMoment &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi, &
                     uDR_Option = uDR )

            dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

          END IF

#if defined( TWOMOMENT_ORDER_V )
          IF( CaptureIMEXStages ) &
            CALL CaptureIMEXStageState &
                   ( Ui, Mi, iS, iIMEX_PreImplicit, &
                     IMEX_Stage_FluidState_Option, &
                     IMEX_Stage_RadiationState_Option )
#endif

        END IF

      END DO ! jS = 1, iS - 1

#if defined( TWOMOMENT_ORDER_V )
      IF( iS == 1 .AND. CaptureIMEXStages )THEN
        CALL CaptureIMEXStageState &
               ( Ui, Mi, iS, iIMEX_Assembled, &
                 IMEX_Stage_FluidState_Option, &
                 IMEX_Stage_RadiationState_Option )
        CALL CaptureIMEXStageState &
               ( Ui, Mi, iS, iIMEX_PreImplicit, &
                 IMEX_Stage_FluidState_Option, &
                 IMEX_Stage_RadiationState_Option )
      END IF
#endif

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

#if defined( TWOMOMENT_ORDER_V )
        IF(      w_IM(iS) .NE. Zero &
            .AND. (      PRESENT( MatterCollisionExchange_Rate_Option ) &
                    .OR. PRESENT( MatterUpdate_Rate_Option ) &
                    .OR. PRESENT( MatterCollisionRHS_FinalStage_Option ) &
                    .OR. PRESENT( &
                          IMEX_Stage_LeptonFlavorCheck_Rate_Option ) ) )THEN

          IF( PRESENT( IMEX_Stage_LeptonFlavorCheck_Rate_Option ) )THEN

            IF(      iS == nStages &
                .AND. PRESENT( MatterCollisionRHS_FinalStage_Option ) )THEN

              CALL ComputeIncrement_TwoMoment_Implicit &
                     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                       Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM, &
                       AccumulateMatterCoupling_Mask = AccumulateMatterDiagnostics, &
                       uDR_Option = uDR, &
                       MatterCollisionExchange_Rate_Option = MatterCollisionExchange_Rate_Option, &
                       MatterUpdate_Rate_Option = MatterUpdate_Rate_Option, &
                       IMEX_ImplicitWeight_Option = w_IM(iS), &
                       MatterCollisionRHS_Option &
                         = MatterCollisionRHS_FinalStage_Option, &
                       LeptonFlavorCheck_Rate_Option &
                         = IMEX_Stage_LeptonFlavorCheck_Rate_Option &
                             (:,:,:,:,:,:,iS) )

            ELSE

              CALL ComputeIncrement_TwoMoment_Implicit &
                     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                       Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM, &
                       AccumulateMatterCoupling_Mask = AccumulateMatterDiagnostics, &
                       uDR_Option = uDR, &
                       MatterCollisionExchange_Rate_Option = MatterCollisionExchange_Rate_Option, &
                       MatterUpdate_Rate_Option = MatterUpdate_Rate_Option, &
                       IMEX_ImplicitWeight_Option = w_IM(iS), &
                       LeptonFlavorCheck_Rate_Option &
                         = IMEX_Stage_LeptonFlavorCheck_Rate_Option &
                             (:,:,:,:,:,:,iS) )

            END IF

          ELSE

            IF(      iS == nStages &
                .AND. PRESENT( MatterCollisionRHS_FinalStage_Option ) )THEN

              CALL ComputeIncrement_TwoMoment_Implicit &
                     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                       Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM, &
                       AccumulateMatterCoupling_Mask = AccumulateMatterDiagnostics, &
                       uDR_Option = uDR, &
                       MatterCollisionExchange_Rate_Option = MatterCollisionExchange_Rate_Option, &
                       MatterUpdate_Rate_Option = MatterUpdate_Rate_Option, &
                       IMEX_ImplicitWeight_Option = w_IM(iS), &
                       MatterCollisionRHS_Option &
                         = MatterCollisionRHS_FinalStage_Option )

            ELSE

              CALL ComputeIncrement_TwoMoment_Implicit &
                     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                       Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM, &
                       AccumulateMatterCoupling_Mask = AccumulateMatterDiagnostics, &
                       uDR_Option = uDR, &
                       MatterCollisionExchange_Rate_Option = MatterCollisionExchange_Rate_Option, &
                       MatterUpdate_Rate_Option = MatterUpdate_Rate_Option, &
                       IMEX_ImplicitWeight_Option = w_IM(iS) )

            END IF

          END IF

          AccumulateMatterDiagnostics = .TRUE.

        ELSE
#endif

          CALL ComputeIncrement_TwoMoment_Implicit &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                   Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM, &
                   AccumulateMatterCoupling_Mask = .FALSE., &
                   uDR_Option = uDR )

#if defined( TWOMOMENT_ORDER_V )
        END IF
#endif

        CALL AddToArray( One, Ui, dt * a_IM(iS,iS), StageData(iS) % dU_IM )
        CALL AddToArray( One, Mi, dt * a_IM(iS,iS), StageData(iS) % dM_IM )

#if defined( TWOMOMENT_ORDER_V )
        IF( PRESENT( IMEX_Stage_FluidImplicitRate_Option ) ) &
          CALL CopyArray &
                 ( IMEX_Stage_FluidImplicitRate_Option(:,:,:,:,:,iS), &
                   One, StageData(iS) % dU_IM )
        IF( PRESENT( IMEX_Stage_RadiationImplicitRate_Option ) ) &
          CALL CopyArray &
                 ( IMEX_Stage_RadiationImplicitRate_Option(:,:,:,:,:,:,:,iS), &
                   One, StageData(iS) % dM_IM )
#endif

#if defined( TWOMOMENT_ORDER_V )
        IF( CaptureIMEXStages ) &
          CALL CaptureIMEXStageState &
                 ( Ui, Mi, iS, iIMEX_PostImplicit, &
                   IMEX_Stage_FluidState_Option, &
                   IMEX_Stage_RadiationState_Option )
#endif

        IF( EvolveEuler )THEN

          CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

        END IF

        IF( EvolveTwoMoment )THEN

#if defined( TWOMOMENT_ORDER_V )
          IF(      iS == nStages .AND. .NOT. FinalAssembly &
              .AND. PRESENT( RadiationState_BeforeFinalLimiter_Option ) ) &
            CALL CopyArray( RadiationState_BeforeFinalLimiter_Option, One, Mi )
#endif

          CALL ApplyPositivityLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi, &
                   uDR_Option = uDR )

          dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

        END IF

#if defined( TWOMOMENT_ORDER_V )
        IF( CaptureIMEXStages ) &
          CALL CaptureIMEXStageState &
                 ( Ui, Mi, iS, iIMEX_PostImplicitLimiter, &
                   IMEX_Stage_FluidState_Option, &
                   IMEX_Stage_RadiationState_Option )
#endif

      END IF

#if defined( TWOMOMENT_ORDER_V )
      IF( CaptureIMEXStages &
          .AND. .NOT. ( ANY( a_IM(:,iS) .NE. Zero ) &
                  .OR. ( w_IM(iS) .NE. Zero ) ) )THEN
        CALL CaptureIMEXStageState &
               ( Ui, Mi, iS, iIMEX_PostImplicit, &
                 IMEX_Stage_FluidState_Option, &
                 IMEX_Stage_RadiationState_Option )
        CALL CaptureIMEXStageState &
               ( Ui, Mi, iS, iIMEX_PostImplicitLimiter, &
                 IMEX_Stage_FluidState_Option, &
                 IMEX_Stage_RadiationState_Option )
      END IF
#endif

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        IF( EvolveEuler )THEN

          CALL TimersStart( Timer_Euler )

          CALL ComputeIncrement_Euler_DG_Explicit &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, GX, &
                   Ui, uDF, StageData(iS) % dU_EX )

          StageData(iS) % OffGridFlux_U &
            =   OffGridFlux_Euler_X1_Outer &
              - OffGridFlux_Euler_X1_Inner &
              + OffGridFlux_Euler_X2_Outer &
              - OffGridFlux_Euler_X2_Inner &
              + OffGridFlux_Euler_X3_Outer &
              - OffGridFlux_Euler_X3_Inner

          CALL TimersStop( Timer_Euler )

        END IF

        IF( EvolveTwoMoment )THEN

          IF( PRESENT( ComputeIncrement_TwoMoment_Explicit_Option ) )THEN

            CALL ComputeIncrement_TwoMoment_Explicit_Option &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                     Ui, Mi, StageData(iS) % dM_EX )

          ELSE

            CALL ComputeIncrement_TwoMoment_Explicit &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                     Ui, Mi, StageData(iS) % dM_EX )

          END IF

          StageData(iS) % OffGridFlux_M = OffGridFlux_TwoMoment
          StageData(iS) % OffGridLeptonFlux_M &
            = OffGridLeptonFlux_TwoMoment

        END IF

#if defined( TWOMOMENT_ORDER_V )
        IF( PRESENT( IMEX_Stage_FluidExplicitRate_Option ) ) &
          CALL CopyArray &
                 ( IMEX_Stage_FluidExplicitRate_Option(:,:,:,:,:,iS), &
                   One, StageData(iS) % dU_EX )
        IF( PRESENT( IMEX_Stage_RadiationExplicitRate_Option ) ) &
          CALL CopyArray &
                 ( IMEX_Stage_RadiationExplicitRate_Option(:,:,:,:,:,:,:,iS), &
                   One, StageData(iS) % dM_EX )
#endif

      END IF

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    IF( FinalAssembly )THEN

      CALL CopyArray( Ui, One, U0 )
      CALL CopyArray( Mi, One, M0 )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * w_IM(iS), StageData(iS) % dU_IM )
          CALL AddToArray( One, Mi, dt * w_IM(iS), StageData(iS) % dM_IM )

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * w_EX(iS), StageData(iS) % dU_EX )
          CALL AddToArray( One, Mi, dt * w_EX(iS), StageData(iS) % dM_EX )

        END IF

      END DO

      IF( EvolveEuler )THEN

        CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

        CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

      END IF

      IF( EvolveTwoMoment )THEN

        CALL ApplySlopeLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

#if defined( TWOMOMENT_ORDER_V )
        IF( PRESENT( RadiationState_BeforeFinalLimiter_Option ) ) &
          CALL CopyArray( RadiationState_BeforeFinalLimiter_Option, One, Mi )
#endif

        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi, &
                 uDR_Option = uDR )

        dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

      END IF

    END IF

    IF( EvolveEuler .AND. PRESENT( SolveGravity ) )THEN

      CALL TimersStart( Timer_Poisson )

#if   defined( THORNADO_OMP_OL )
              !$OMP TARGET UPDATE FROM( Ui )
#elif defined( THORNADO_OACC   )
              !$ACC UPDATE HOST       ( Ui )
#endif

      CALL SolveGravity &
             ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

#if   defined( THORNADO_OMP_OL )
              !$OMP TARGET UPDATE TO( GX )
#elif defined( THORNADO_OACC   )
              !$ACC UPDATE DEVICE   ( GX )
#endif

      CALL TimersStop( Timer_Poisson )

    END IF

#if defined( TWOMOMENT_ORDER_V )
    IF( CaptureIMEXStages .AND. FinalAssembly ) &
      CALL CaptureIMEXStageState &
             ( Ui, Mi, nStages, iIMEX_CompletedStep, &
               IMEX_Stage_FluidState_Option, &
               IMEX_Stage_RadiationState_Option )
#endif

    CALL CopyArray( U, One, Ui )
    CALL CopyArray( M, One, Mi )

    IF( EvolveEuler )THEN

      DO iS = 1, nStages

        dU_OffGrid &
          = dU_OffGrid + dt * w_EX(iS) * StageData(iS) % OffGridFlux_U

      END DO

      IF( UpdateTallies ) &
        CALL IncrementOffGridTally_Euler( dU_OffGrid )

    END IF

    IF( EvolveTwoMoment )THEN

      DO iS = 1, nStages

        dM_OffGrid &
          = dM_OffGrid + dt * w_EX(iS) * StageData(iS) % OffGridFlux_M
        dM_LeptonFlavor_OffGrid &
          = dM_LeptonFlavor_OffGrid &
              + dt * w_EX(iS) * StageData(iS) % OffGridLeptonFlux_M

      END DO

      IF( UpdateTallies )THEN
        CALL IncrementOffGridTally_TwoMoment( dM_OffGrid )
        CALL IncrementPositivityLimiterTally_TwoMoment( dM_PL )
      END IF

    END IF

    IF( PRESENT( OffGridFluidIncrement_Option ) ) &
      OffGridFluidIncrement_Option = dU_OffGrid
    IF( PRESENT( OffGridRadiationIncrement_Option ) ) &
      OffGridRadiationIncrement_Option = dM_OffGrid
    IF( PRESENT( OffGridRadiationLeptonFlavorIncrement_Option ) ) &
      OffGridRadiationLeptonFlavorIncrement_Option &
        = dM_LeptonFlavor_OffGrid

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: GX, U, M ) &
    !$OMP MAP( release: GE )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( GX, U, M ) &
    !$ACC DELETE( GE )
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: MatterCollisionExchange_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( MatterCollisionExchange_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterUpdate_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: MatterUpdate_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( MatterUpdate_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterCollisionRHS_FinalStage_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: MatterCollisionRHS_FinalStage_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( MatterCollisionRHS_FinalStage_Option )
#endif
    END IF

    IF( PRESENT( RadiationState_BeforeFinalLimiter_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: RadiationState_BeforeFinalLimiter_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( RadiationState_BeforeFinalLimiter_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_FluidState_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_FluidState_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_FluidState_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationState_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_RadiationState_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_RadiationState_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_FluidImplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_FluidImplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_FluidImplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_FluidExplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_FluidExplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_FluidExplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationImplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_RadiationImplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_RadiationImplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_RadiationExplicitRate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA MAP( from: IMEX_Stage_RadiationExplicitRate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_RadiationExplicitRate_Option )
#endif
    END IF

    IF( PRESENT( IMEX_Stage_LeptonFlavorCheck_Rate_Option ) )THEN
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: IMEX_Stage_LeptonFlavorCheck_Rate_Option )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA COPYOUT( IMEX_Stage_LeptonFlavorCheck_Rate_Option )
#endif
    END IF
#endif

    CALL TimersStop( Timer_IMEX )

  END SUBROUTINE Update_IMEX_RK


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE CaptureIMEXStageState &
    ( U_State, M_State, iStage, iState, &
      FluidStates_Option, RadiationStates_Option )

    REAL(DP), INTENT(in) :: U_State(:,:,:,:,:)
    REAL(DP), INTENT(in) :: M_State(:,:,:,:,:,:,:)
    INTEGER,  INTENT(in) :: iStage, iState
    REAL(DP), INTENT(out), OPTIONAL :: FluidStates_Option(:,:,:,:,:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: RadiationStates_Option(:,:,:,:,:,:,:,:,:)

    IF( PRESENT( FluidStates_Option ) ) &
      CALL CopyArray &
             ( FluidStates_Option(:,:,:,:,:,iStage,iState), One, U_State )

    IF( PRESENT( RadiationStates_Option ) ) &
      CALL CopyArray &
             ( RadiationStates_Option(:,:,:,:,:,:,:,iStage,iState), &
               One, M_State )

  END SUBROUTINE CaptureIMEXStageState


  INTEGER FUNCTION GetNumberOfStages_IMEX_RK()

    GetNumberOfStages_IMEX_RK = nStages

  END FUNCTION GetNumberOfStages_IMEX_RK


  REAL(DP) FUNCTION GetImplicitStageCoefficient_IMEX_RK( iStage, jStage )

    INTEGER, INTENT(in) :: iStage
    INTEGER, INTENT(in), OPTIONAL :: jStage

    GetImplicitStageCoefficient_IMEX_RK = a_IM(iStage,iStage)
    IF( PRESENT( jStage ) ) &
      GetImplicitStageCoefficient_IMEX_RK = a_IM(iStage,jStage)

  END FUNCTION GetImplicitStageCoefficient_IMEX_RK


  REAL(DP) FUNCTION GetExplicitStageCoefficient_IMEX_RK( iStage, jStage )

    INTEGER, INTENT(in) :: iStage, jStage

    GetExplicitStageCoefficient_IMEX_RK = a_EX(iStage,jStage)

  END FUNCTION GetExplicitStageCoefficient_IMEX_RK


  REAL(DP) FUNCTION GetExplicitStageWeight_IMEX_RK( iStage )

    INTEGER, INTENT(in) :: iStage

    GetExplicitStageWeight_IMEX_RK = w_EX(iStage)

  END FUNCTION GetExplicitStageWeight_IMEX_RK


  REAL(DP) FUNCTION GetImplicitStageWeight_IMEX_RK( iStage )

    INTEGER, INTENT(in) :: iStage

    GetImplicitStageWeight_IMEX_RK = w_IM(iStage)

  END FUNCTION GetImplicitStageWeight_IMEX_RK


  LOGICAL FUNCTION ExplicitRHSUsed_IMEX_RK( iStage )

    INTEGER, INTENT(in) :: iStage

    ExplicitRHSUsed_IMEX_RK &
      = ANY( a_EX(:,iStage) .NE. Zero ) .OR. ( w_EX(iStage) .NE. Zero )

  END FUNCTION ExplicitRHSUsed_IMEX_RK


  LOGICAL FUNCTION HasFinalAssembly_IMEX_RK()

    HasFinalAssembly_IMEX_RK &
      =      ANY( a_IM(nStages,:) .NE. w_IM(:) ) &
        .OR. ANY( a_EX(nStages,:) .NE. w_EX(:) )

  END FUNCTION HasFinalAssembly_IMEX_RK
#endif


  SUBROUTINE Initialize_IMEX_RK &
    ( Scheme, EvolveEuler_Option, EvolveTwoMoment_Option )

    CHARACTER(LEN=*), INTENT(in)           :: Scheme
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveEuler_Option
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveTwoMoment_Option

    INTEGER :: i

    EvolveEuler = .FALSE.
    IF( PRESENT( EvolveEuler_Option ) )THEN
      EvolveEuler = EvolveEuler_Option
    END IF

    EvolveTwoMoment = .TRUE.
    IF( PRESENT( EvolveTwoMoment_Option ) )THEN
      EvolveTwoMoment = EvolveTwoMoment_Option
    END IF

    WRITE(*,*)
    WRITE(*,'(A6,A19,A)' ) '', 'IMEX-RK Scheme: '   , TRIM( Scheme )
    WRITE(*,'(A6,A19,L1)') '', 'Evolve Euler: '     , EvolveEuler
    WRITE(*,'(A6,A19,L1)') '', 'Evolve Two-Moment: ', EvolveTwoMoment

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'SSPRK1' )

        nStages = 1

        CALL AllocateButcherTables

        a_EX(1,1) = 0.0_DP
        w_EX(1)   = 1.0_DP

      CASE ( 'SSPRK2' )

        nStages = 2

        CALL AllocateButcherTables

        a_EX(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_EX(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 'SSPRK3' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_EX(1:3)   = [ 1.00_DP, 1.00_DP, 4.00_DP ] / 6.0_DP

      CASE ( 'BackwardEuler' )

        nStages = 1

        CALL AllocateButcherTables

        a_IM(1,1) = 1.0_DP
        w_IM(1)   = 1.0_DP

      CASE ( 'IMEX_ARS_111' )

        nStages = 2

        CALL AllocateButcherTables

        ! --- Coefficients from Ascher et al. (1997) ---

        a_EX(2,1) = 1.0_DP
        w_EX(1)   = a_EX(2,1)

        a_IM(2,2) = 1.0_DP
        w_IM(2)   = a_IM(2,2)

      CASE ( 'IMEX_PDARS' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 0.5_DP
        a_IM(3,3) = 0.5_DP

        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE ( 'IMEX_RKCB2' )

        ! --- Coefficients from Cavaglieri & Bewley (2015) ---
	! --- JCP, 286, 172-193 ------------------------------

        nStages = 3

        CALL AllocateButcherTables

        a_EX(2,1) = 0.4_DP
        a_EX(3,2) = 1.0_DP

        w_EX(2)   = 5.0_DP / 6.0_DP
        w_EX(3)   = 1.0_DP / 6.0_DP

        a_IM(2,2) = 0.4_DP
        a_IM(3,2) = 5.0_DP / 6.0_DP
        a_IM(3,3) = 1.0_DP / 6.0_DP

        w_IM(2)   = 5.0_DP / 6.0_DP
        w_IM(3)   = 1.0_DP / 6.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A6,A,A)') &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'BackwardEuler'
        WRITE(*,'(A6,A)') '', 'IMEX_ARS_111'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS'
        WRITE(*,'(A6,A)') '', 'IMEX_RKCB2'
        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Implicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_IM(1:nStages)

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)

    CALL AllocateArray( U0 )
    CALL AllocateArray( Ui )
    CALL AllocateArray( M0 )
    CALL AllocateArray( Mi )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages

      CALL AllocateArray( StageData(i) % dU_IM )
      CALL AllocateArray( StageData(i) % dU_EX )
      CALL AllocateArray( StageData(i) % dM_IM )
      CALL AllocateArray( StageData(i) % dM_EX )

      StageData(i) % OffGridFlux_U = Zero
      StageData(i) % OffGridFlux_M = Zero
      StageData(i) % OffGridLeptonFlux_M = Zero

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#endif

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    INTEGER :: i

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#endif

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

    CALL DeallocateArray( U0 )
    CALL DeallocateArray( Ui )
    CALL DeallocateArray( M0 )
    CALL DeallocateArray( Mi )

    DO i = 1, nStages

      CALL DeallocateArray( StageData(i) % dU_IM )
      CALL DeallocateArray( StageData(i) % dU_EX )
      CALL DeallocateArray( StageData(i) % dM_IM )
      CALL DeallocateArray( StageData(i) % dM_EX )

    END DO

    DEALLOCATE( StageData )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables

    ! --- Implicit Coefficients ---

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = Zero
    w_IM = Zero
    a_IM = Zero

    ! --- Explicit Coefficients ---

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = Zero
    w_EX = Zero
    a_EX = Zero

  END SUBROUTINE AllocateButcherTables


  SUBROUTINE AllocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

    ALLOCATE &
      ( Array5D(nDOFX, &
                iX_B1(1):iX_E1(1), &
                iX_B1(2):iX_E1(2), &
                iX_B1(3):iX_E1(3), &
                nCF) )

    Array5D = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Array5D )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( Array5D )
#endif

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AllocateArray5D


  SUBROUTINE DeallocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Array5D )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( Array5D )
#endif

    DEALLOCATE( Array5D )

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE DeallocateArray5D


  SUBROUTINE AllocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1), &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3), &
                iZ_B1(4):iZ_E1(4), &
                nCR,nSpecies) )

    Array7D = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Array7D )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( Array7D )
#endif

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AllocateArray7D


  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Array7D )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( Array7D )
#endif

    DEALLOCATE( Array7D )

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE DeallocateArray7D


  SUBROUTINE CopyArray5D( U, alpha, V )

    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, V, iX_B1, iX_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        U(iNodeX,iX1,iX2,iX3,iCF) &
          = alpha * V(iNodeX,iX1,iX2,iX3,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE CopyArray5D


  SUBROUTINE CopyArray7D( U, alpha, V )

    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( U, V, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = alpha * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE CopyArray7D


  SUBROUTINE AddToArray5D( alpha, U, beta, V )

    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(in)    :: &
      beta
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, V, iX_B1, iX_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        U(iNodeX,iX1,iX2,iX3,iCF) &
          = alpha * U(iNodeX,iX1,iX2,iX3,iCF) &
            + beta * V(iNodeX,iX1,iX2,iX3,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AddToArray5D


  SUBROUTINE AddToArray7D( alpha, U, beta, V )

    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      beta
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( U, V, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = alpha * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + beta * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AddToArray7D


END MODULE TwoMoment_TimeSteppingModule
