MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
  USE mpi, ONLY: MPI_COMM_WORLD, MPI_ABORT
  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    FourPi
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions, &
    Timer_Collisions_PrimitiveFluid, &
    Timer_Collisions_PrimitiveTwoMoment, &
    Timer_Collisions_Solve
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, iCF_Nm, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, iPF_Nm, &
    nAF, iAF_T, iAF_E , iAF_Ye, iAF_Ym, iAF_P
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    nDR, iDR_iter_outer, iDR_iter_inner
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE MeshModule, ONLY: &
    MeshE
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    SolveNeutrinoMatterCoupling, &
    InitializeNeutrinoMatterSolver, &
    FinalizeNeutrinoMatterSolver, &
    InitializeNeutrinoMatterSolverParameters
#if defined( TWOMOMENT_ORDER_V )
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    ComputeMatterRHS_Collisions_Flattened_OrderV => ComputeMatterRHS_Collisions_OrderV, &
    ComputeMeanInverseMeanFreePaths, &
    ComputeTransportOpacity, &
    AssembleTransportOpacity, &
    ComputeOpacities_Packed, &
    ComputeRates_Packed, &
    UseNewtonSolve, NewtonRadiationIncrement, NewtonFluidIncrement, &
    nQ, nC, &
    nLeptonFamilies, nLeptonFlavorCheck, &
    iLFC_InputConversion, iLFC_EntryProjection, &
    iLFC_CollisionSource, iLFC_RawMap, &
    iLFC_MapProjection, iLFC_FixedPointResidual, &
    iLFC_EquationResidual, &
    iQ_Le, iQ_Lm, iQ_E, iQ_H, iQ_M1, iQ_M2, iQ_M3
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    AllocateSourceAwareTransportOpacity, &
    StoreSourceAwareTransportOpacity, &
    CacheSourceAwareTransportOpacity_Flattened, &
    RestoreSourceAwareTransportOpacity, &
    SetSourceAwareLLFActive, &
    MarkStreamingBoundaryConditionsApplied, &
    UseSourceAwareLLF
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    SetOpacityMask, &
    nOp
  USE NeutrinoLeptonInteractionsModule, ONLY: &
    ComputeNeutrinoLeptonPairingResiduals
#endif
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, &
    ComputeConserved_TwoMoment
#if   defined( TWOMOMENT_ORDER_1 )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeConserved_Euler_NonRelativistic

#elif defined( TWOMOMENT_ORDER_V )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeConserved_Euler_NonRelativistic

#elif defined( TWOMOMENT_RELATIVISTIC )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputePressureFromPrimitive_TABLE
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    ComputeConserved_Euler_Relativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit
  PUBLIC :: GetNeutrinoMatterSolverConvergenceDiagnostics
  PUBLIC :: SetNeutrinoMatterSolverConvergenceDiagnostics
#if defined( TWOMOMENT_ORDER_V )
  PUBLIC :: ComputeMatterRHS_Collisions_OrderV
  PUBLIC :: ComputeMeanInverseMeanFreePaths_OrderV
  PUBLIC :: ComputeTransportOpacityFromConserved_OrderV
  PUBLIC :: ComputeIncrement_TwoMoment_Explicit_AP
  PUBLIC :: nQ, nC
  PUBLIC :: nLeptonFamilies, nLeptonFlavorCheck
  PUBLIC :: iLFC_InputConversion, iLFC_EntryProjection
  PUBLIC :: iLFC_CollisionSource, iLFC_RawMap
  PUBLIC :: iLFC_MapProjection, iLFC_FixedPointResidual
  PUBLIC :: iLFC_EquationResidual
  PUBLIC :: nConservationHandoffCheck
  PUBLIC :: iCHC_InputConversion_L, iCHC_PrimitivePair_L
  PUBLIC :: iCHC_FinalConversion_L, iCHC_MappedPair_L, iCHC_ReturnedPair_L
  PUBLIC :: iCHC_InputConversion_E, iCHC_PrimitivePair_E
  PUBLIC :: iCHC_FinalConversion_E, iCHC_MappedPair_E, iCHC_ReturnedPair_E
#endif

#if defined( TWOMOMENT_ORDER_V )
  INTEGER, PARAMETER :: iCHC_InputConversion_L = 1
  INTEGER, PARAMETER :: iCHC_PrimitivePair_L   = 2
  INTEGER, PARAMETER :: iCHC_FinalConversion_L = 3
  INTEGER, PARAMETER :: iCHC_MappedPair_L      = 4
  INTEGER, PARAMETER :: iCHC_ReturnedPair_L    = 5
  INTEGER, PARAMETER :: iCHC_InputConversion_E = 6
  INTEGER, PARAMETER :: iCHC_PrimitivePair_E   = 7
  INTEGER, PARAMETER :: iCHC_FinalConversion_E = 8
  INTEGER, PARAMETER :: iCHC_MappedPair_E      = 9
  INTEGER, PARAMETER :: iCHC_ReturnedPair_E    = 10
  INTEGER, PARAMETER :: nConservationHandoffCheck  = 10
#endif

  INTEGER               :: nE_G, nX_G
  INTEGER               :: nZ(4), nX(3), nE
  INTEGER               :: iE_B0, iE_E0, iX_B0(3), iX_E0(3)
  INTEGER               :: iE_B1, iE_E1, iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: GE_N(:,:)
  REAL(DP), ALLOCATABLE :: EnergyQuadratureWeight_N(:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)

  REAL(DP), ALLOCATABLE, TARGET :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: PR_N(:,:,:,:)

  REAL(DP), DIMENSION(:), CONTIGUOUS, POINTER :: N_P, G1_P, G2_P, G3_P
  REAL(DP), DIMENSION(:), CONTIGUOUS, POINTER :: J_P, H1_P, H2_P, H3_P

  INTEGER,  DIMENSION(:), ALLOCATABLE :: PositionIndexZ

  INTEGER, ALLOCATABLE :: nIterations_Inner(:)
  INTEGER, ALLOCATABLE :: nIterations_Outer(:)
  INTEGER, ALLOCATABLE :: nIterations_Prim(:)
  INTEGER, ALLOCATABLE :: nIterations_Inner_Final(:)
  INTEGER, ALLOCATABLE :: nIterations_Outer_Final(:)
  REAL(DP), ALLOCATABLE :: RelativeResidual_Inner_First(:)
#if defined( TWOMOMENT_ORDER_V )
  REAL(DP), ALLOCATABLE :: MatterEquationResidual_Rate_Final(:,:)
  REAL(DP), ALLOCATABLE :: InitialRadiationLimiter_Rate_Final(:,:)
  REAL(DP), ALLOCATABLE :: RadiationResidual_Relative_Final(:)
  REAL(DP), ALLOCATABLE :: &
    CollisionRecoveryRadiationResidual_Relative_Final(:)
  REAL(DP), ALLOCATABLE :: CollisionRecoveryMatterHeating_Rate_Final(:)
  INTEGER,  ALLOCATABLE :: CollisionRecoveryStatus_Final(:)
  REAL(DP), ALLOCATABLE :: ConservationHandoffCheck_Rate_Final(:,:)
  REAL(DP), ALLOCATABLE :: LeptonFlavorCheck_Rate_Final(:,:,:)
#endif

  LOGICAL :: RetainSolverConvergenceDiagnostics = .FALSE.
  INTEGER, ALLOCATABLE :: nIterations_Inner_Final_Last(:)
  INTEGER, ALLOCATABLE :: nIterations_Outer_Final_Last(:)
  REAL(DP), ALLOCATABLE :: RelativeResidual_Inner_First_Last(:)
#if defined( TWOMOMENT_ORDER_V )
  REAL(DP), ALLOCATABLE :: MatterEquationResidual_Rate_Final_Last(:,:)
  REAL(DP), ALLOCATABLE :: InitialRadiationLimiter_Rate_Final_Last(:,:)
  REAL(DP), ALLOCATABLE :: RadiationResidual_Relative_Final_Last(:)
  REAL(DP), ALLOCATABLE :: &
    CollisionRecoveryRadiationResidual_Relative_Final_Last(:)
  REAL(DP), ALLOCATABLE :: &
    CollisionRecoveryMatterHeating_Rate_Final_Last(:)
  INTEGER,  ALLOCATABLE :: CollisionRecoveryStatus_Final_Last(:)
  REAL(DP), ALLOCATABLE :: ConservationHandoffCheck_Rate_Final_Last(:,:)
  REAL(DP), ALLOCATABLE :: LeptonFlavorCheck_Rate_Final_Last(:,:,:)
#endif

  LOGICAL, PARAMETER :: ReportConvergenceData = .FALSE.

CONTAINS


  ! --- Public Subroutines ---


  SUBROUTINE SetNeutrinoMatterSolverConvergenceDiagnostics( Retain )

    LOGICAL, INTENT(in) :: Retain

    IF( Retain .AND. ALLOCATED( nIterations_Inner_Final_Last ) )THEN
      DEALLOCATE( nIterations_Inner_Final_Last )
      DEALLOCATE( nIterations_Outer_Final_Last )
      DEALLOCATE( RelativeResidual_Inner_First_Last )
#if defined( TWOMOMENT_ORDER_V )
      DEALLOCATE( MatterEquationResidual_Rate_Final_Last )
      DEALLOCATE( InitialRadiationLimiter_Rate_Final_Last )
      DEALLOCATE( RadiationResidual_Relative_Final_Last )
      DEALLOCATE( CollisionRecoveryRadiationResidual_Relative_Final_Last )
      DEALLOCATE( CollisionRecoveryMatterHeating_Rate_Final_Last )
      DEALLOCATE( CollisionRecoveryStatus_Final_Last )
      DEALLOCATE( ConservationHandoffCheck_Rate_Final_Last )
      DEALLOCATE( LeptonFlavorCheck_Rate_Final_Last )
#endif
    END IF

    RetainSolverConvergenceDiagnostics = Retain

  END SUBROUTINE SetNeutrinoMatterSolverConvergenceDiagnostics


  SUBROUTINE GetNeutrinoMatterSolverConvergenceDiagnostics &
    ( InnerIterationsFinal, OuterIterationsFinal, &
      FirstInnerRelativeResidual &
#if defined( TWOMOMENT_ORDER_V )
      , MatterEquationResidualRateFinal, &
        InitialRadiationLimiterRateFinal, &
        RadiationResidualRelativeFinal, &
        CollisionRecoveryRadiationResidualRelativeFinal, &
        CollisionRecoveryMatterHeatingRateFinal, &
        CollisionRecoveryStatusFinal, &
        ConservationHandoffCheckRateFinal, LeptonFlavorCheckRateFinal &
#endif
      )

    INTEGER, ALLOCATABLE, INTENT(out) :: InnerIterationsFinal(:)
    INTEGER, ALLOCATABLE, INTENT(out) :: OuterIterationsFinal(:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: FirstInnerRelativeResidual(:)
#if defined( TWOMOMENT_ORDER_V )
    REAL(DP), ALLOCATABLE, INTENT(out) :: MatterEquationResidualRateFinal(:,:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: InitialRadiationLimiterRateFinal(:,:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: RadiationResidualRelativeFinal(:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: &
      CollisionRecoveryRadiationResidualRelativeFinal(:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: &
      CollisionRecoveryMatterHeatingRateFinal(:)
    INTEGER,  ALLOCATABLE, INTENT(out) :: CollisionRecoveryStatusFinal(:)
    REAL(DP), ALLOCATABLE, INTENT(out) :: ConservationHandoffCheckRateFinal(:,:)
    REAL(DP), ALLOCATABLE, INTENT(out), OPTIONAL :: &
      LeptonFlavorCheckRateFinal(:,:,:)
#endif

    IF( .NOT. ALLOCATED( nIterations_Inner_Final_Last ) ) &
      CALL AbortCollisionDiscretization &
             ( 'No retained neutrino-matter solver convergence diagnostics' )

    ALLOCATE( InnerIterationsFinal &
                (SIZE(nIterations_Inner_Final_Last)) )
    ALLOCATE( OuterIterationsFinal &
                (SIZE(nIterations_Outer_Final_Last)) )
    ALLOCATE( FirstInnerRelativeResidual &
                (SIZE(RelativeResidual_Inner_First_Last)) )
#if defined( TWOMOMENT_ORDER_V )
    ALLOCATE( MatterEquationResidualRateFinal, &
              MOLD = MatterEquationResidual_Rate_Final_Last )
    ALLOCATE( InitialRadiationLimiterRateFinal, &
              MOLD = InitialRadiationLimiter_Rate_Final_Last )
    ALLOCATE( RadiationResidualRelativeFinal, &
              MOLD = RadiationResidual_Relative_Final_Last )
    ALLOCATE( CollisionRecoveryRadiationResidualRelativeFinal, &
              MOLD = CollisionRecoveryRadiationResidual_Relative_Final_Last )
    ALLOCATE( CollisionRecoveryMatterHeatingRateFinal, &
              MOLD = CollisionRecoveryMatterHeating_Rate_Final_Last )
    ALLOCATE( CollisionRecoveryStatusFinal, &
              MOLD = CollisionRecoveryStatus_Final_Last )
    ALLOCATE( ConservationHandoffCheckRateFinal, &
              MOLD = ConservationHandoffCheck_Rate_Final_Last )
    IF( PRESENT( LeptonFlavorCheckRateFinal ) ) &
      ALLOCATE( LeptonFlavorCheckRateFinal, &
                MOLD = LeptonFlavorCheck_Rate_Final_Last )
#endif

    InnerIterationsFinal       = nIterations_Inner_Final_Last
    OuterIterationsFinal       = nIterations_Outer_Final_Last
    FirstInnerRelativeResidual = RelativeResidual_Inner_First_Last
#if defined( TWOMOMENT_ORDER_V )
    MatterEquationResidualRateFinal = MatterEquationResidual_Rate_Final_Last
    InitialRadiationLimiterRateFinal = InitialRadiationLimiter_Rate_Final_Last
    RadiationResidualRelativeFinal = RadiationResidual_Relative_Final_Last
    CollisionRecoveryRadiationResidualRelativeFinal &
      = CollisionRecoveryRadiationResidual_Relative_Final_Last
    CollisionRecoveryMatterHeatingRateFinal &
      = CollisionRecoveryMatterHeating_Rate_Final_Last
    CollisionRecoveryStatusFinal = CollisionRecoveryStatus_Final_Last
    ConservationHandoffCheckRateFinal = ConservationHandoffCheck_Rate_Final_Last
    IF( PRESENT( LeptonFlavorCheckRateFinal ) ) &
      LeptonFlavorCheckRateFinal = LeptonFlavorCheck_Rate_Final_Last
#endif

  END SUBROUTINE GetNeutrinoMatterSolverConvergenceDiagnostics


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE ComputeIncrement_TwoMoment_Explicit_AP &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL,  INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    INTEGER :: iZ_B_Opacity(4), iZ_E_Opacity(4)
    INTEGER :: iX_B0_Local(3), iX_E0_Local(3)
    INTEGER :: iX_B1_Local(3), iX_E1_Local(3)
    IF( UseSourceAwareLLF )THEN

      iX_B0_Local = iZ_B0_Input(2:4)
      iX_E0_Local = iZ_E0_Input(2:4)
      iX_B1_Local = iZ_B1_Input(2:4)
      iX_E1_Local = iZ_E1_Input(2:4)

      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0_Local, iX_E0_Local, &
               iX_B1_Local, iX_E1_Local, U_F )

      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0_Input, iZ_E0_Input, &
               iZ_B1_Input, iZ_E1_Input, U_R )

      iZ_B_Opacity(1)   = iZ_B0_Input(1)
      iZ_E_Opacity(1)   = iZ_E0_Input(1)
      iZ_B_Opacity(2:4) = iZ_B1_Input(2:4)
      iZ_E_Opacity(2:4) = iZ_E1_Input(2:4)

      CALL ComputeTransportOpacityFromConserved_OrderV &
             ( iZ_B_Opacity, iZ_E_Opacity, iZ_B1_Input, iZ_E1_Input, &
               GE, GX, U_F, U_R )

      CALL RestoreSourceAwareTransportOpacity

      CALL MarkStreamingBoundaryConditionsApplied

    END IF

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R, dU_R, &
             LLFPenaltyFactor_Option = LLFPenaltyFactor_Option, &
             SurfaceTermFactor_Option = SurfaceTermFactor_Option, &
             VolumeTermFactor_Option = VolumeTermFactor_Option, &
             IncludeObserverCorrections_Option &
               = IncludeObserverCorrections_Option )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit_AP


  SUBROUTINE ComputeTransportOpacityFromConserved_OrderV &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, TransportOpacity_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out), OPTIONAL :: TransportOpacity_Option(:,:,:)
    INTEGER :: iN_X
    LOGICAL :: PatchHasCollisions
    REAL(DP), ALLOCATABLE :: TransportOpacity_N(:,:,:)

    CALL InitializeCollisions &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input )

    ALLOCATE( TransportOpacity_N(nE_G,nSpecies,nX_G) )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GE, GX, U_F, U_R ) &
    !$OMP MAP( alloc: TransportOpacity_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GE, GX, U_F, U_R ) &
    !$ACC CREATE( TransportOpacity_N )
#endif

    CALL MapDataForCollisions &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), CF_N(iN_X,iCF_S2), CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), CF_N(iN_X,iCF_Ne), CF_N(iN_X,iCF_Nm), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), PF_N(iN_X,iPF_V2), PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), PF_N(iN_X,iPF_Ne), PF_N(iN_X,iPF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, J_P, H1_P, H2_P, H3_P, &
             PF_N(:,iPF_V1), PF_N(:,iPF_V2), PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), PositionIndexZ, nIterations_Prim )

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym) )

    CALL ComputeTransportOpacity &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym), &
             GX_N(:,iGF_SqrtGm), &
             PR_N(:,:,:,iPR_D), &
             PR_N(:,:,:,iPR_I1), PR_N(:,:,:,iPR_I2), PR_N(:,:,:,iPR_I3), &
             TransportOpacity_N, &
             PatchHasCollisions_Option = PatchHasCollisions )

    IF( PRESENT( TransportOpacity_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( TransportOpacity_N )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( TransportOpacity_N )
    !$ACC WAIT
#endif
      TransportOpacity_Option = TransportOpacity_N
    END IF

    CALL SetSourceAwareLLFActive( PatchHasCollisions )

    IF( PatchHasCollisions )THEN

      CALL AllocateSourceAwareTransportOpacity &
             ( iZ_B0_Input, iZ_E0_Input )

      CALL StoreSourceAwareTransportOpacity &
             ( iZ_B0_Input, iZ_E0_Input, nX, TransportOpacity_N )

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: GE, GX, U_F, U_R, TransportOpacity_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( GE, GX, U_F, U_R, TransportOpacity_N )
    !$ACC WAIT
#endif

    DEALLOCATE( TransportOpacity_N )
    CALL FinalizeCollisions

  END SUBROUTINE ComputeTransportOpacityFromConserved_OrderV


  SUBROUTINE ComputeMeanInverseMeanFreePaths_OrderV &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, PF, PR, MeanInvMFP )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      PR(1:nDOFZ, &
         iZ_B1_Input(1):iZ_E1_Input(1), &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(out) :: MeanInvMFP(:,:,:)

    CALL InitializeCollisions &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             ComputeFirstMomentRates_Option = .TRUE. )

    CALL MapDataForDiagnostics( GE, GX, PF, PR )

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym) )

    CALL ComputeMeanInverseMeanFreePaths &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym), &
             GX_N(:,iGF_SqrtGm), GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), GX_N(:,iGF_Gm_dd_33), &
             PR_N(:,:,:,iPR_D), &
             PR_N(:,:,:,iPR_I1), PR_N(:,:,:,iPR_I2), PR_N(:,:,:,iPR_I3), &
             MeanInvMFP )

    CALL FinalizeCollisions

  END SUBROUTINE ComputeMeanInverseMeanFreePaths_OrderV


  SUBROUTINE ComputeMatterRHS_Collisions_OrderV &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, PF, PR, Q_Matter, NNS_RateDecomposition_Option, &
      QH_ByProcessSpecies_Option, EmAb_RateDecomposition_Option, &
      EmAb_Spectrum_Option, RadiationCollisionRHS_Option, &
      GeneralLeptonPairingResidual_Option, &
      GeneralLeptonPairingScale_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      PR(1:nDOFZ, &
         iZ_B1_Input(1):iZ_E1_Input(1), &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(out) :: Q_Matter(:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: NNS_RateDecomposition_Option(:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: QH_ByProcessSpecies_Option(:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: EmAb_RateDecomposition_Option(:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: EmAb_Spectrum_Option(:,:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: RadiationCollisionRHS_Option(:,:,:,:)
    REAL(DP), INTENT(out), OPTIONAL :: &
      GeneralLeptonPairingResidual_Option(:,:), &
      GeneralLeptonPairingScale_Option(:,:)

    LOGICAL, ALLOCATABLE :: MatterCoupling_Mask_Diagnostic(:)
    LOGICAL, ALLOCATABLE :: Opacity_Mask_Diagnostic(:,:)
    REAL(DP), ALLOCATABLE :: GeneralLeptonWeight2_Diagnostic(:)

    CALL InitializeCollisions &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             ComputeFirstMomentRates_Option = .TRUE. )

    CALL MapDataForDiagnostics( GE, GX, PF, PR )

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym) )

    ALLOCATE( MatterCoupling_Mask_Diagnostic(nX_G) )
    ALLOCATE( Opacity_Mask_Diagnostic(nOp,nX_G) )

    CALL SetOpacityMask &
           ( nX_G, nDOFX, PF_N(:,iPF_D), AF_N(:,iAF_T), &
             AF_N(:,iAF_Ye), AF_N(:,iAF_Ym), &
             Opacity_Mask_Diagnostic )

    CALL ComputeOpacities_Packed &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym), &
             GX_N(:,iGF_SqrtGm), nX_P0 = nX_G, &
             T0 = AF_N(:,iAF_T), Ye0 = AF_N(:,iAF_Ye), &
             Ym0 = AF_N(:,iAF_Ym), &
             Opacity_Mask_Input = Opacity_Mask_Diagnostic )

    CALL ComputeRates_Packed &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), PR_N(:,:,:,iPR_D), &
             PR_N(:,:,:,iPR_I1), PR_N(:,:,:,iPR_I2), PR_N(:,:,:,iPR_I3), &
             nX_P0 = nX_G, Opacity_Mask_Input = Opacity_Mask_Diagnostic )

    IF( PRESENT(GeneralLeptonPairingResidual_Option) .OR. &
        PRESENT(GeneralLeptonPairingScale_Option) )THEN
      IF( .NOT. PRESENT(GeneralLeptonPairingResidual_Option) .OR. &
          .NOT. PRESENT(GeneralLeptonPairingScale_Option) ) &
        ERROR STOP 'Both general-lepton pairing outputs are required.'
      ALLOCATE( GeneralLeptonWeight2_Diagnostic(nE_G) )
      GeneralLeptonWeight2_Diagnostic &
        = EnergyQuadratureWeight_N * GE_N(:,iGE_Ep2)
      CALL ComputeNeutrinoLeptonPairingResiduals &
             ( GeneralLeptonWeight2_Diagnostic, &
               GeneralLeptonPairingResidual_Option, &
               GeneralLeptonPairingScale_Option )
      DEALLOCATE( GeneralLeptonWeight2_Diagnostic )
    END IF

    MatterCoupling_Mask_Diagnostic = ANY( Opacity_Mask_Diagnostic, DIM = 1 )

    CALL ComputeMatterRHS_Collisions_Flattened_OrderV &
           ( PR_N(:,:,:,iPR_D), &
             PR_N(:,:,:,iPR_I1), PR_N(:,:,:,iPR_I2), PR_N(:,:,:,iPR_I3), &
             PF_N(:,iPF_V1), PF_N(:,iPF_V2), PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), MatterCoupling_Mask_Diagnostic, &
             Q_Matter, &
             NNS_RateDecomposition_Option = NNS_RateDecomposition_Option, &
             QH_ByProcessSpecies_Option = QH_ByProcessSpecies_Option, &
             EmAb_RateDecomposition_Option = EmAb_RateDecomposition_Option, &
             EmAb_Spectrum_Option = EmAb_Spectrum_Option, &
             RadiationCollisionRHS_Option = RadiationCollisionRHS_Option )

    DEALLOCATE( MatterCoupling_Mask_Diagnostic )
    DEALLOCATE( Opacity_Mask_Diagnostic )

    CALL FinalizeCollisions

  END SUBROUTINE ComputeMatterRHS_Collisions_OrderV
#endif


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, &
      uDR_Option, MatterCoupling_Mask_Option &
#if defined( TWOMOMENT_ORDER_V )
      , MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
        IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
        LeptonFlavorCheck_Rate_Option &
#endif
      )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      dt
    REAL(DP), INTENT(in) :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in) :: &
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
    ! MatterCollisionExchange_Rate_Option accumulates the collision RHS with
    ! the supplied IMEX implicit weight.  Over a completed timestep it is the
    ! stage-weighted collision exchange expressed as a rate; this is the
    ! quantity exported through the Qxxx plot variables.
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
          (1:nQ, &
           1:1+nSpecies+nC, &
           1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4))
    ! MatterUpdate_Rate_Option accumulates the complete fluid update returned
    ! by the implicit call, including changes not represented by the collision
    ! RHS decomposition.
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
          (1:nQ, &
           1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
    ! The unweighted collision RHS for this implicit stage.  This is retained
    ! for correctness checks and is not the completed-timestep Qxxx quantity.
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

    INTEGER :: iN_X, iN_E, iS
    INTEGER :: iNodeX, iX1, iX2, iX3
    LOGICAL, ALLOCATABLE :: MatterCoupling_Mask_N(:)
#if defined( TWOMOMENT_ORDER_V )
    INTEGER :: iG, iQ
    REAL(DP) :: IMEX_ImplicitWeight
    REAL(DP), ALLOCATABLE :: MatterCollisionRHS_N(:,:,:)
    REAL(DP), ALLOCATABLE :: TransportOpacity_N(:,:,:)
#endif

    CALL TimersStart( Timer_Collisions )

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    ALLOCATE( MatterCoupling_Mask_N(nX_G) )
#if defined( TWOMOMENT_ORDER_V )
    IF(      PRESENT( MatterCollisionExchange_Rate_Option ) &
        .OR. PRESENT( MatterCollisionRHS_Option ) ) &
      ALLOCATE( MatterCollisionRHS_N(nQ,1+nSpecies+nC,nX_G) )

    IMEX_ImplicitWeight = One
    IF( PRESENT( IMEX_ImplicitWeight_Option ) ) &
      IMEX_ImplicitWeight = IMEX_ImplicitWeight_Option

    IF( UseSourceAwareLLF ) &
      ALLOCATE( TransportOpacity_N(nE_G,nSpecies,nX_G) )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: dU_F, dU_R, MatterCoupling_Mask_N )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE(     dU_F, dU_R, MatterCoupling_Mask_N )
#endif
#if defined( TWOMOMENT_ORDER_V )
    IF(      PRESENT( MatterCollisionExchange_Rate_Option ) &
        .OR. PRESENT( MatterCollisionRHS_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: MatterCollisionRHS_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( MatterCollisionRHS_N )
#endif
    END IF

    IF( PRESENT( MatterCollisionExchange_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: MatterCollisionExchange_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( MatterCollisionExchange_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterCollisionRHS_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: MatterCollisionRHS_Option )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( MatterCollisionRHS_Option )
#endif
    END IF

    IF( PRESENT( MatterUpdate_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: MatterUpdate_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( MatterUpdate_Rate_Option )
#endif
    END IF

    IF( UseSourceAwareLLF )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: TransportOpacity_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( TransportOpacity_N )
#endif
    END IF

    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: LeptonFlavorCheck_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( LeptonFlavorCheck_Rate_Option )
#endif
    END IF

#endif

    CALL MapDataForCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    ! --- Compute Primitive Fluid ---

    CALL TimersStart( Timer_Collisions_PrimitiveFluid )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif



#if   defined( TWOMOMENT_ORDER_1 )

    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

#elif defined( TWOMOMENT_ORDER_V )

    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

#elif defined( TWOMOMENT_RELATIVISTIC )


    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_Relativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

#endif


    CALL TimersStop( Timer_Collisions_PrimitiveFluid )

    ! --- Compute Primitive Moments ---

    CALL TimersStart( Timer_Collisions_PrimitiveTwoMoment )

#if   defined( TWOMOMENT_ORDER_1 )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_ORDER_V )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             PositionIndexZ, &
             nIterations_Prim )

#elif defined( TWOMOMENT_RELATIVISTIC )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             GX_N(:,iGF_Alpha   ), &
             GX_N(:,iGF_Beta_1  ), &
             GX_N(:,iGF_Beta_2  ), &
             GX_N(:,iGF_Beta_3  ), &
             PositionIndexZ, &
             nIterations_Prim )

#endif

    CALL TimersStop( Timer_Collisions_PrimitiveTwoMoment )

    ! --- EOS Table Lookup ---

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym) )

#if defined( TWOMOMENT_ORDER_V )
    IF( RetainSolverConvergenceDiagnostics ) &
      CALL InitializeConservationHandoffCheck()
#endif

    ! --- Neutrino-Matter Coupling Solve ---

    CALL TimersStart( Timer_Collisions_Solve )

    CALL SolveNeutrinoMatterCoupling &
           ( dt, &
             PR_N(:,:,:,iCR_N ), &
             PR_N(:,:,:,iCR_G1), &
             PR_N(:,:,:,iCR_G2), &
             PR_N(:,:,:,iCR_G3), &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             PF_N(:,iPF_D ), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_Ye), &
             AF_N(:,iAF_Ym), &
             AF_N(:,iAF_E ), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             GX_N(:,iGF_Alpha   ), &
             GX_N(:,iGF_Beta_1  ), &
             GX_N(:,iGF_Beta_2  ), &
             GX_N(:,iGF_Beta_3  ), &
             CR_N(:,:,:,iCR_N ), &
             CR_N(:,:,:,iCR_G1), &
             CR_N(:,:,:,iCR_G2), &
             CR_N(:,:,:,iCR_G3), &
             nIterations_Inner, &
             nIterations_Outer, &
             nIterations_Inner_Final, &
             nIterations_Outer_Final, &
             RelativeResidual_Inner_First, &
             MatterCoupling_Mask_N &
#if defined( TWOMOMENT_ORDER_V )
             , RetainSolverConvergenceDiagnostics, &
               MatterEquationResidual_Rate_Final, &
               InitialRadiationLimiter_Rate_Final, &
               RadiationResidual_Relative_Final, &
               CollisionRecoveryRadiationResidual_Relative_Final, &
               CollisionRecoveryMatterHeating_Rate_Final, &
               CollisionRecoveryStatus_Final, &
               LeptonFlavorCheck_Rate_Final &
#endif
             )

    CALL TimersStop( Timer_Collisions_Solve )

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( LeptonFlavorCheck_Rate_Option, &
    !$ACC          LeptonFlavorCheck_Rate_Final )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) PRIVATE( iN_X )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX
      DO iG      = 1, nLeptonFlavorCheck
      DO iQ      = 1, nLeptonFamilies

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        LeptonFlavorCheck_Rate_Option(iQ,iG,iNodeX,iX1,iX2,iX3) &
          = LeptonFlavorCheck_Rate_Final(iQ,iG,iN_X)

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    IF( UseSourceAwareLLF )THEN

      CALL AssembleTransportOpacity( TransportOpacity_N )

      CALL CacheSourceAwareTransportOpacity_Flattened &
             ( iZ_B0, iZ_E0, nX, TransportOpacity_N )

    END IF

    IF(      PRESENT( MatterCollisionExchange_Rate_Option ) &
        .OR. PRESENT( MatterCollisionRHS_Option ) )THEN

      CALL ComputeMatterRHS_Collisions_Flattened_OrderV &
             ( PR_N(:,:,:,iCR_N ), &
               PR_N(:,:,:,iCR_G1), PR_N(:,:,:,iCR_G2), PR_N(:,:,:,iCR_G3), &
               PF_N(:,iPF_V1), PF_N(:,iPF_V2), PF_N(:,iPF_V3), &
               GX_N(:,iGF_Gm_dd_11), GX_N(:,iGF_Gm_dd_22), &
               GX_N(:,iGF_Gm_dd_33), MatterCoupling_Mask_N, &
               MatterCollisionRHS_N, &
               UseAcceptedEmAbRate_Option = UseNewtonSolve )

    END IF
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iN_X ) FIRSTPRIVATE( IMEX_ImplicitWeight )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC FIRSTPRIVATE( IMEX_ImplicitWeight ) &
    !$ACC PRESENT( MatterCollisionExchange_Rate_Option, MatterCollisionRHS_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iN_X ) FIRSTPRIVATE( IMEX_ImplicitWeight )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX
      DO iG      = 1, 1+nSpecies+nC
      DO iQ      = 1, nQ

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        IF( AccumulateMatterCoupling_Mask )THEN
          MatterCollisionExchange_Rate_Option(iQ,iG,iNodeX,iX1,iX2,iX3) &
            = MatterCollisionExchange_Rate_Option(iQ,iG,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * MatterCollisionRHS_N(iQ,iG,iN_X)
        ELSE
          MatterCollisionExchange_Rate_Option(iQ,iG,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * MatterCollisionRHS_N(iQ,iG,iN_X)
        END IF

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    END IF
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterCollisionRHS_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( MatterCollisionRHS_Option, MatterCollisionRHS_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iN_X )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX
      DO iG      = 1, 1+nSpecies+nC
      DO iQ      = 1, nQ

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        MatterCollisionRHS_Option(iQ,iG,iNodeX,iX1,iX2,iX3) &
          = MatterCollisionRHS_N(iQ,iG,iN_X)

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    END IF
#endif

    IF( PRESENT( MatterCoupling_Mask_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iNodeX, iN_X )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iNodeX, iN_X ) &
    !$ACC PRESENT( MatterCoupling_Mask_Option, &
    !$ACC          MatterCoupling_Mask_N, iX_B0, iX_E0, nX )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iNodeX, iN_X )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        iNodeX = 1

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        IF( AccumulateMatterCoupling_Mask )THEN
          MatterCoupling_Mask_Option(iX1,iX2,iX3) &
            = MatterCoupling_Mask_Option(iX1,iX2,iX3) .OR. &
              MatterCoupling_Mask_N(iN_X)
        ELSE
          MatterCoupling_Mask_Option(iX1,iX2,iX3) &
            = MatterCoupling_Mask_N(iN_X)
        END IF

      END DO
      END DO
      END DO

    END IF

    ! --- Store Iteration Counts (If Requested) ---

    IF( PRESENT( uDR_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( uDR_Option, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        uDR_Option(iX1,iX2,iX3,iDR_iter_outer) = 0.0_DP
        uDR_Option(iX1,iX2,iX3,iDR_iter_inner) = 0.0_DP

      END DO
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( uDR_Option, iX_B0, iX_E0, nX, &
    !$ACC          nIterations_Outer, nIterations_Inner )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iN_X )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        IF ( nIterations_Outer(iN_X) == 0 ) THEN
          nIterations_Inner(iN_X) = 0
        ELSE
          nIterations_Inner(iN_X) &
            = FLOOR(   REAL( nIterations_Inner(iN_X), DP ) &
                     / REAL( nIterations_Outer(iN_X), DP ) )

        ENDIF

        uDR_Option(iX1,iX2,iX3,iDR_iter_outer) &
          = MAX( uDR_Option(iX1,iX2,iX3,iDR_iter_outer), &
                 REAL( nIterations_Outer(iN_X), DP ) )

        uDR_Option(iX1,iX2,iX3,iDR_iter_inner) &
          = MAX( uDR_Option(iX1,iX2,iX3,iDR_iter_inner), &
                 REAL( nIterations_Inner(iN_X), DP ) )

      END DO
      END DO
      END DO
      END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( PR_N, CR_N, PF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

#if   defined( TWOMOMENT_ORDER_1 )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_ORDER_V )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_RELATIVISTIC )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33), &
               Zero, Zero, Zero, &
               GX_N(iN_X,iGF_Alpha   ), &
               GX_N(iN_X,iGF_Beta_1  ), &
               GX_N(iN_X,iGF_Beta_2  ), &
               GX_N(iN_X,iGF_Beta_3  ) )

#endif

    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( PF_N, CF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif



#if   defined( TWOMOMENT_ORDER_1 )

    DO iN_X = 1, nX_G

      PF_N(iN_X,iPF_E ) = AF_N(iN_X,iAF_E ) * PF_N(iN_X,iPF_D)
      PF_N(iN_X,iPF_Ne) = AF_N(iN_X,iAF_Ye) * PF_N(iN_X,iPF_D) / AtomicMassUnit

      CALL ComputeConserved_Euler_NonRelativistic &
             ( PF_N(iN_X,iPF_D),  &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO


#elif defined( TWOMOMENT_ORDER_V )

    DO iN_X = 1, nX_G

      PF_N(iN_X,iPF_E ) = AF_N(iN_X,iAF_E ) * PF_N(iN_X,iPF_D)
      PF_N(iN_X,iPF_Ne) = AF_N(iN_X,iAF_Ye) * PF_N(iN_X,iPF_D) / AtomicMassUnit
      PF_N(iN_X,iPF_Nm) = AF_N(iN_X,iAF_Ym) * PF_N(iN_X,iPF_D) / AtomicMassUnit

      CALL ComputeConserved_Euler_NonRelativistic &
             ( PF_N(iN_X,iPF_D),  &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO


#elif defined( TWOMOMENT_RELATIVISTIC )


    CALL ComputeThermodynamicStates_Primitive_TABLE &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), AF_N(:,iAF_Ym), &
             PF_N(:,iPF_E), AF_N(:,iAF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm) )

    CALL ComputePressureFromPrimitive_TABLE & 
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), PF_N(:,iPF_Nm), &
             AF_N(:,iAF_P) )


    DO iN_X = 1, nX_G

      CALL ComputeConserved_Euler_Relativistic &
             ( PF_N(iN_X,iPF_D),  &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               PF_N(iN_X,iPF_Nm), &
               CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               CF_N(iN_X,iCF_Nm), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33), &
               AF_N(iN_X,iAF_P) )

    END DO


#endif


    CALL ComputeAndMapIncrement &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, U_F, U_R, &
             MatterCoupling_Mask_N, dU_F, dU_R )

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterUpdate_Rate_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iN_X ) FIRSTPRIVATE( IMEX_ImplicitWeight )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC FIRSTPRIVATE( IMEX_ImplicitWeight ) &
    !$ACC PRESENT( MatterUpdate_Rate_Option, dU_F, PF_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iN_X ) FIRSTPRIVATE( IMEX_ImplicitWeight )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        IF( AccumulateMatterCoupling_Mask )THEN
          MatterUpdate_Rate_Option(iQ_Le,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_Le,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_Ne)
          MatterUpdate_Rate_Option(iQ_Lm,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_Lm,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_Nm)
          MatterUpdate_Rate_Option(iQ_E,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_E,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_E)
          MatterUpdate_Rate_Option(iQ_M1,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_M1,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S1)
          MatterUpdate_Rate_Option(iQ_M2,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_M2,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S2)
          MatterUpdate_Rate_Option(iQ_M3,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_M3,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S3)
          MatterUpdate_Rate_Option(iQ_H,iNodeX,iX1,iX2,iX3) &
            = MatterUpdate_Rate_Option(iQ_H,iNodeX,iX1,iX2,iX3) &
              + IMEX_ImplicitWeight &
                * ( dU_F(iNodeX,iX1,iX2,iX3,iCF_E) &
                    - PF_N(iN_X,iPF_V1) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S1) &
                    - PF_N(iN_X,iPF_V2) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S2) &
                    - PF_N(iN_X,iPF_V3) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S3) )
        ELSE
          MatterUpdate_Rate_Option(iQ_Le,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_Ne)
          MatterUpdate_Rate_Option(iQ_Lm,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_Nm)
          MatterUpdate_Rate_Option(iQ_E,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_E)
          MatterUpdate_Rate_Option(iQ_M1,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S1)
          MatterUpdate_Rate_Option(iQ_M2,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S2)
          MatterUpdate_Rate_Option(iQ_M3,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight * dU_F(iNodeX,iX1,iX2,iX3,iCF_S3)
          MatterUpdate_Rate_Option(iQ_H,iNodeX,iX1,iX2,iX3) &
            = IMEX_ImplicitWeight &
              * ( dU_F(iNodeX,iX1,iX2,iX3,iCF_E) &
                  - PF_N(iN_X,iPF_V1) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S1) &
                  - PF_N(iN_X,iPF_V2) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S2) &
                  - PF_N(iN_X,iPF_V3) * dU_F(iNodeX,iX1,iX2,iX3,iCF_S3) )
        END IF

      END DO
      END DO
      END DO
      END DO

    END IF


    IF( UseSourceAwareLLF )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: TransportOpacity_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( TransportOpacity_N )
#endif
      DEALLOCATE( TransportOpacity_N )
    END IF
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_F, dU_R ) &
    !$OMP MAP( release: GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               MatterCoupling_Mask_N )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_F, dU_R ) &
    !$ACC DELETE(       GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC               MatterCoupling_Mask_N )
#endif

    CALL FinalizeCollisions

#if defined( TWOMOMENT_ORDER_V )
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: MatterCollisionExchange_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( MatterCollisionExchange_Rate_Option )
#endif
    END IF

    IF( PRESENT( MatterCollisionRHS_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: MatterCollisionRHS_Option )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( MatterCollisionRHS_Option )
#endif
    END IF

    IF(      PRESENT( MatterCollisionExchange_Rate_Option ) &
        .OR. PRESENT( MatterCollisionRHS_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: MatterCollisionRHS_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( MatterCollisionRHS_N )
#endif
      DEALLOCATE( MatterCollisionRHS_N )
    END IF

    IF( PRESENT( MatterUpdate_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: MatterUpdate_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( MatterUpdate_Rate_Option )
#endif
    END IF

    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: LeptonFlavorCheck_Rate_Option )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( LeptonFlavorCheck_Rate_Option )
#endif
    END IF
#endif
    DEALLOCATE( MatterCoupling_Mask_N )

    CALL TimersStop( Timer_Collisions )

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  ! --- Private Subroutines ---


  SUBROUTINE InitializeCollisions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, ComputeFirstMomentRates_Option )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)
    LOGICAL, INTENT(in), OPTIONAL :: ComputeFirstMomentRates_Option

    INTEGER :: iE, iN_X, iN_E, iNodeE, iS, iZ, nZ_G

    iE_B0 = iZ_B0(1)  ; iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1)  ; iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1
    nZ = iZ_E0 - iZ_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )
    nZ_G = nE_G * nX_G * nSpecies

    ALLOCATE( GE_N(nE_G,nGE) )
    ALLOCATE( EnergyQuadratureWeight_N(nE_G) )
    ALLOCATE( GX_N(nX_G,nGF) )

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )

    ALLOCATE( CR_N(nE_G,nSpecies,nX_G,nCR) )
    ALLOCATE( PR_N(nE_G,nSpecies,nX_G,nCR) )

    ALLOCATE( PositionIndexZ(nZ_G) )

    DO iN_E = 1, nE_G
      iE     = ( iN_E - 1 ) / nDOFE + iE_B0
      iNodeE = MOD( iN_E - 1, nDOFE ) + 1
      EnergyQuadratureWeight_N(iN_E) &
        = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) &
            / ( PlanckConstant * SpeedOfLight )**3
    END DO

    iZ = 0
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      iZ = iZ + 1

      PositionIndexZ(iZ) = iN_X

    END DO
    END DO
    END DO

    ALLOCATE( nIterations_Inner(nX_G) )
    ALLOCATE( nIterations_Outer(nX_G) )
    ALLOCATE( nIterations_Inner_Final(nX_G) )
    ALLOCATE( nIterations_Outer_Final(nX_G) )
    ALLOCATE( RelativeResidual_Inner_First(nX_G) )
#if defined( TWOMOMENT_ORDER_V )
    ALLOCATE( MatterEquationResidual_Rate_Final(nQ,nX_G) )
    ALLOCATE( InitialRadiationLimiter_Rate_Final(nQ,nX_G) )
    ALLOCATE( RadiationResidual_Relative_Final(nX_G) )
    ALLOCATE( CollisionRecoveryRadiationResidual_Relative_Final(nX_G) )
    ALLOCATE( CollisionRecoveryMatterHeating_Rate_Final(nX_G) )
    ALLOCATE( CollisionRecoveryStatus_Final(nX_G) )
    ALLOCATE( LeptonFlavorCheck_Rate_Final &
                (nLeptonFamilies,nLeptonFlavorCheck,nX_G) )
    IF( RetainSolverConvergenceDiagnostics ) &
      ALLOCATE &
        ( ConservationHandoffCheck_Rate_Final &
            (nConservationHandoffCheck,nX_G) )
#endif
    ALLOCATE( nIterations_Prim (nZ_G) )

    nIterations_Inner(:) = 0
    nIterations_Outer(:) = 0
    nIterations_Inner_Final(:) = 0
    nIterations_Outer_Final(:) = 0
    RelativeResidual_Inner_First(:) = Zero
#if defined( TWOMOMENT_ORDER_V )
    MatterEquationResidual_Rate_Final(:,:) = Zero
    InitialRadiationLimiter_Rate_Final(:,:) = Zero
    RadiationResidual_Relative_Final(:) = Zero
    CollisionRecoveryRadiationResidual_Relative_Final(:) = Zero
    CollisionRecoveryMatterHeating_Rate_Final(:) = Zero
    CollisionRecoveryStatus_Final(:) = -2
    LeptonFlavorCheck_Rate_Final(:,:,:) = Zero
    IF( RetainSolverConvergenceDiagnostics ) &
      ConservationHandoffCheck_Rate_Final(:,:) = Zero
#endif
    nIterations_Prim (:) = 0

    N_P (1:nZ_G) => CR_N(:,:,:,iCR_N )
    G1_P(1:nZ_G) => CR_N(:,:,:,iCR_G1)
    G2_P(1:nZ_G) => CR_N(:,:,:,iCR_G2)
    G3_P(1:nZ_G) => CR_N(:,:,:,iCR_G3)

    J_P (1:nZ_G) => PR_N(:,:,:,iCR_N )
    H1_P(1:nZ_G) => PR_N(:,:,:,iCR_G1)
    H2_P(1:nZ_G) => PR_N(:,:,:,iCR_G2)
    H3_P(1:nZ_G) => PR_N(:,:,:,iCR_G3)

    ! --- Neutrino-Matter Solver Parameter Initialization ---
    ! --- ( can be moved to the program init ) --------------

    CALL InitializeNeutrinoMatterSolverParameters()

    CALL InitializeNeutrinoMatterSolver &
           ( iZ_B0, iZ_E0, ComputeFirstMomentRates_Option )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$OMP          nIterations_Inner, nIterations_Outer, nIterations_Prim, &
    !$OMP          nIterations_Inner_Final, nIterations_Outer_Final, &
    !$OMP          RelativeResidual_Inner_First ) &
    !$OMP MAP( to: EnergyQuadratureWeight_N ) &
    !$OMP MAP( alloc: GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$ACC         nIterations_Inner, nIterations_Outer, nIterations_Prim, &
    !$ACC         nIterations_Inner_Final, nIterations_Outer_Final, &
    !$ACC         RelativeResidual_Inner_First ) &
    !$ACC COPYIN( EnergyQuadratureWeight_N ) &
    !$ACC CREATE( GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#endif

#if defined( TWOMOMENT_ORDER_V )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: MatterEquationResidual_Rate_Final, &
    !$OMP             InitialRadiationLimiter_Rate_Final, &
    !$OMP             RadiationResidual_Relative_Final, &
    !$OMP             CollisionRecoveryRadiationResidual_Relative_Final, &
    !$OMP             CollisionRecoveryMatterHeating_Rate_Final, &
    !$OMP             CollisionRecoveryStatus_Final, &
    !$OMP             LeptonFlavorCheck_Rate_Final )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC CREATE( MatterEquationResidual_Rate_Final, &
    !$ACC         InitialRadiationLimiter_Rate_Final, &
    !$ACC         RadiationResidual_Relative_Final, &
    !$ACC         CollisionRecoveryRadiationResidual_Relative_Final, &
    !$ACC         CollisionRecoveryMatterHeating_Rate_Final, &
    !$ACC         CollisionRecoveryStatus_Final, &
    !$ACC         LeptonFlavorCheck_Rate_Final )
#endif

    IF( RetainSolverConvergenceDiagnostics )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: ConservationHandoffCheck_Rate_Final )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA CREATE( ConservationHandoffCheck_Rate_Final )
#endif
    END IF
#endif

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

#if defined( TWOMOMENT_ORDER_V )
    IF( RetainSolverConvergenceDiagnostics )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MatterEquationResidual_Rate_Final, &
    !$OMP                     InitialRadiationLimiter_Rate_Final, &
    !$OMP                     RadiationResidual_Relative_Final, &
    !$OMP                     CollisionRecoveryRadiationResidual_Relative_Final, &
    !$OMP                     CollisionRecoveryMatterHeating_Rate_Final, &
    !$OMP                     CollisionRecoveryStatus_Final, &
    !$OMP                     ConservationHandoffCheck_Rate_Final, &
    !$OMP                     LeptonFlavorCheck_Rate_Final )
#elif defined(THORNADO_OACC  )
    !$ACC UPDATE HOST( MatterEquationResidual_Rate_Final, &
    !$ACC              InitialRadiationLimiter_Rate_Final, &
    !$ACC              RadiationResidual_Relative_Final, &
    !$ACC              CollisionRecoveryRadiationResidual_Relative_Final, &
    !$ACC              CollisionRecoveryMatterHeating_Rate_Final, &
    !$ACC              CollisionRecoveryStatus_Final, &
    !$ACC              ConservationHandoffCheck_Rate_Final, &
    !$ACC              LeptonFlavorCheck_Rate_Final )
#endif
    END IF
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: nIterations_Inner, nIterations_Outer, nIterations_Prim, &
    !$OMP            nIterations_Inner_Final, nIterations_Outer_Final, &
    !$OMP            RelativeResidual_Inner_First ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$OMP               GE_N, EnergyQuadratureWeight_N, GX_N, &
    !$OMP               CF_N, PF_N, AF_N, CR_N, PR_N )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( nIterations_Inner, nIterations_Outer, nIterations_Prim, &
    !$ACC          nIterations_Inner_Final, nIterations_Outer_Final, &
    !$ACC          RelativeResidual_Inner_First ) &
    !$ACC DELETE( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$ACC         GE_N, EnergyQuadratureWeight_N, GX_N, &
    !$ACC         CF_N, PF_N, AF_N, CR_N, PR_N )
#endif

#if defined( TWOMOMENT_ORDER_V )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: MatterEquationResidual_Rate_Final, &
    !$OMP               InitialRadiationLimiter_Rate_Final, &
    !$OMP               RadiationResidual_Relative_Final, &
    !$OMP               CollisionRecoveryRadiationResidual_Relative_Final, &
    !$OMP               CollisionRecoveryMatterHeating_Rate_Final, &
    !$OMP               CollisionRecoveryStatus_Final, &
    !$OMP               LeptonFlavorCheck_Rate_Final )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( MatterEquationResidual_Rate_Final, &
    !$ACC         InitialRadiationLimiter_Rate_Final, &
    !$ACC         RadiationResidual_Relative_Final, &
    !$ACC         CollisionRecoveryRadiationResidual_Relative_Final, &
    !$ACC         CollisionRecoveryMatterHeating_Rate_Final, &
    !$ACC         CollisionRecoveryStatus_Final, &
    !$ACC         LeptonFlavorCheck_Rate_Final )
#endif

    IF( RetainSolverConvergenceDiagnostics )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ConservationHandoffCheck_Rate_Final )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA DELETE( ConservationHandoffCheck_Rate_Final )
#endif
    END IF
#endif

    IF( RetainSolverConvergenceDiagnostics )THEN

      IF( ALLOCATED( nIterations_Inner_Final_Last ) )THEN
        IF( SIZE(nIterations_Inner_Final_Last) /= nX_G )THEN
          DEALLOCATE( nIterations_Inner_Final_Last )
          DEALLOCATE( nIterations_Outer_Final_Last )
          DEALLOCATE( RelativeResidual_Inner_First_Last )
#if defined( TWOMOMENT_ORDER_V )
          DEALLOCATE( MatterEquationResidual_Rate_Final_Last )
          DEALLOCATE( InitialRadiationLimiter_Rate_Final_Last )
          DEALLOCATE( RadiationResidual_Relative_Final_Last )
          DEALLOCATE( CollisionRecoveryRadiationResidual_Relative_Final_Last )
          DEALLOCATE( CollisionRecoveryMatterHeating_Rate_Final_Last )
          DEALLOCATE( CollisionRecoveryStatus_Final_Last )
          DEALLOCATE( ConservationHandoffCheck_Rate_Final_Last )
          DEALLOCATE( LeptonFlavorCheck_Rate_Final_Last )
#endif
        END IF
      END IF

      IF( .NOT. ALLOCATED( nIterations_Inner_Final_Last ) )THEN
        ALLOCATE( nIterations_Inner_Final_Last(nX_G) )
        ALLOCATE( nIterations_Outer_Final_Last(nX_G) )
        ALLOCATE( RelativeResidual_Inner_First_Last(nX_G) )
#if defined( TWOMOMENT_ORDER_V )
        ALLOCATE( MatterEquationResidual_Rate_Final_Last(nQ,nX_G) )
        ALLOCATE( InitialRadiationLimiter_Rate_Final_Last(nQ,nX_G) )
        ALLOCATE( RadiationResidual_Relative_Final_Last(nX_G) )
        ALLOCATE( CollisionRecoveryRadiationResidual_Relative_Final_Last(nX_G) )
        ALLOCATE( CollisionRecoveryMatterHeating_Rate_Final_Last(nX_G) )
        ALLOCATE( CollisionRecoveryStatus_Final_Last(nX_G) )
        ALLOCATE( ConservationHandoffCheck_Rate_Final_Last &
                    (nConservationHandoffCheck,nX_G) )
        ALLOCATE( LeptonFlavorCheck_Rate_Final_Last &
                    (nLeptonFamilies,nLeptonFlavorCheck,nX_G) )
#endif
      END IF

      nIterations_Inner_Final_Last = nIterations_Inner_Final
      nIterations_Outer_Final_Last = nIterations_Outer_Final
      RelativeResidual_Inner_First_Last = RelativeResidual_Inner_First
#if defined( TWOMOMENT_ORDER_V )
      MatterEquationResidual_Rate_Final_Last &
        = MatterEquationResidual_Rate_Final
      InitialRadiationLimiter_Rate_Final_Last &
        = InitialRadiationLimiter_Rate_Final
      RadiationResidual_Relative_Final_Last &
        = RadiationResidual_Relative_Final
      CollisionRecoveryRadiationResidual_Relative_Final_Last &
        = CollisionRecoveryRadiationResidual_Relative_Final
      CollisionRecoveryMatterHeating_Rate_Final_Last &
        = CollisionRecoveryMatterHeating_Rate_Final
      CollisionRecoveryStatus_Final_Last &
        = CollisionRecoveryStatus_Final
      ConservationHandoffCheck_Rate_Final_Last = ConservationHandoffCheck_Rate_Final
      LeptonFlavorCheck_Rate_Final_Last = LeptonFlavorCheck_Rate_Final
#endif

    END IF

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', MINVAL( nIterations_Outer(:) )
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', MAXVAL( nIterations_Outer(:) )
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', DBLE( SUM( nIterations_Outer(:) ) ) / DBLE( nX_G )
      WRITE(*,*)

    END IF

    DEALLOCATE( GE_N, EnergyQuadratureWeight_N, GX_N )
    DEALLOCATE( CF_N, PF_N, AF_N )
    DEALLOCATE( CR_N, PR_N )
    DEALLOCATE( PositionIndexZ )
    DEALLOCATE( nIterations_Inner )
    DEALLOCATE( nIterations_Outer )
    DEALLOCATE( nIterations_Inner_Final )
    DEALLOCATE( nIterations_Outer_Final )
    DEALLOCATE( RelativeResidual_Inner_First )
#if defined( TWOMOMENT_ORDER_V )
    DEALLOCATE( MatterEquationResidual_Rate_Final )
    DEALLOCATE( InitialRadiationLimiter_Rate_Final )
    DEALLOCATE( RadiationResidual_Relative_Final )
    DEALLOCATE( CollisionRecoveryRadiationResidual_Relative_Final )
    DEALLOCATE( CollisionRecoveryMatterHeating_Rate_Final )
    DEALLOCATE( CollisionRecoveryStatus_Final )
    IF( ALLOCATED( ConservationHandoffCheck_Rate_Final ) ) &
      DEALLOCATE( ConservationHandoffCheck_Rate_Final )
    DEALLOCATE( LeptonFlavorCheck_Rate_Final )
#endif
    DEALLOCATE( nIterations_Prim )

    NULLIFY( N_P, G1_P, G2_P, G3_P )
    NULLIFY( J_P, H1_P, H2_P, H3_P )

    CALL FinalizeNeutrinoMatterSolver

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE MapDataForCollisions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in) :: &
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
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iGE, iGF, iCF, iCR, iS
    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeZ
    INTEGER :: iN_E, iN_X

    ! --- Momentum Space Geometry ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iE, iNodeE )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iE, iNodeE ) &
    !$ACC PRESENT( GE_N, GE )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE, iNodeE )
#endif
    DO iGE  = 1, nGE
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      GE_N(iN_E,iGE) = GE(iNodeE,iE,iGE)

    END DO
    END DO

    ! --- Position Space Geometry ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, GX_N, GX )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iGF  = 1, nGF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

    ! --- Fluid Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, U_F )
#elif defined(THORNADO_OMP  )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iCF  = 1, nCF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Radiation Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#endif
    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

      CR_N(iN_E,iS,iN_X,iCR) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapDataForCollisions


  SUBROUTINE MapDataForDiagnostics( GE, GX, PF, PR )

    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iE_B1:iE_E1,1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3),1:nPF)
    REAL(DP), INTENT(in) :: &
      PR(1:nDOFZ, iE_B1:iE_E1, iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3),1:nPR,1:nSpecies)

    INTEGER :: iGF, iN_E, iN_X, iNodeE, iNodeX, iNodeZ
    INTEGER :: iPF, iPR, iS, iE, iX1, iX2, iX3

    ! GE is deliberately accepted to keep this routine's interface aligned
    ! with the collision wrappers.  The two diagnostics do not consume it.

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( GX, GX_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iGF  = 1, nGF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( PF, PF_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iPF  = 1, nPF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      PF_N(iN_X,iPF) = PF(iNodeX,iX1,iX2,iX3,iPF)

    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iE, iNodeE, iNodeX, iNodeZ, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iE, iNodeE, iNodeX, iNodeZ, iX1, iX2, iX3 ) &
    !$ACC PRESENT( PR, PR_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iE, iNodeE, iNodeX, iNodeZ, iX1, iX2, iX3 )
#endif
    DO iS   = 1, nSpecies
    DO iPR  = 1, nPR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE ) + iE_B0
      iNodeE = MOD(  iN_E-1,          nDOFE ) + 1
      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1
      iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

      PR_N(iN_E,iS,iN_X,iPR) = PR(iNodeZ,iE,iX1,iX2,iX3,iPR,iS)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapDataForDiagnostics


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE InitializeConservationHandoffCheck

    INTEGER :: iE, iN_E, iN_X, iNodeE, iS
    REAL(DP) :: PhaseWeight2, PhaseWeight3
    REAL(DP) :: V_d_1, V_d_2, V_d_3, vDotI

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iE, iNodeE, iS, PhaseWeight2, PhaseWeight3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotI )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iE, iNodeE, iS, PhaseWeight2, PhaseWeight3, &
    !$ACC          V_d_1, V_d_2, V_d_3, vDotI ) &
    !$ACC PRESENT( GE_N, EnergyQuadratureWeight_N, GX_N, PF_N, PR_N, &
    !$ACC          ConservationHandoffCheck_Rate_Final )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iE, iNodeE, iS, PhaseWeight2, PhaseWeight3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotI )
#endif
    DO iN_X = 1, nX_G

      ConservationHandoffCheck_Rate_Final(:,iN_X) = Zero

      V_d_1 = GX_N(iN_X,iGF_Gm_dd_11) * PF_N(iN_X,iPF_V1) &
                / SpeedOfLight
      V_d_2 = GX_N(iN_X,iGF_Gm_dd_22) * PF_N(iN_X,iPF_V2) &
                / SpeedOfLight
      V_d_3 = GX_N(iN_X,iGF_Gm_dd_33) * PF_N(iN_X,iPF_V3) &
                / SpeedOfLight

      DO iS = 1, nSpecies
      DO iN_E = 1, nE_G

        iE     = ( iN_E - 1 ) / nDOFE + iE_B0
        iNodeE = MOD( iN_E - 1, nDOFE ) + 1
        PhaseWeight2 = EnergyQuadratureWeight_N(iN_E) &
                         * GE_N(iN_E,iGE_Ep2)
        PhaseWeight3 = EnergyQuadratureWeight_N(iN_E) &
                         * GE_N(iN_E,iGE_Ep3)
        vDotI = V_d_1 * PR_N(iN_E,iS,iN_X,iPR_I1) &
              + V_d_2 * PR_N(iN_E,iS,iN_X,iPR_I2) &
              + V_d_3 * PR_N(iN_E,iS,iN_X,iPR_I3)

        ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_L,iN_X) &
          = ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_L,iN_X) &
              + PhaseWeight2 * LeptonNumber(iS) &
                  * ( PR_N(iN_E,iS,iN_X,iPR_D) + vDotI )
        ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_E,iN_X) &
          = ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_E,iN_X) &
              + PhaseWeight3 &
                  * ( PR_N(iN_E,iS,iN_X,iPR_D) + 2.0_DP * vDotI )

      END DO
      END DO

    END DO

  END SUBROUTINE InitializeConservationHandoffCheck
#endif


  SUBROUTINE ComputeAndMapIncrement &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, U_F, U_R, &
      MatterCoupling_Mask_N, dU_F, dU_R )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      dt
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
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
    LOGICAL, INTENT(in) :: MatterCoupling_Mask_N(:)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iCF, iCR, iS
    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeZ
    INTEGER :: iN_E, iN_X

    ! --- Fluid Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( nX, iX_B0, iX_E0, dU_F, CF_N, U_F, &
    !$ACC          MatterCoupling_Mask_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      iN_X = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

      IF ( MatterCoupling_Mask_N(iN_X) ) THEN
#if defined(TWOMOMENT_ORDER_V)
        IF( UseNewtonSolve )THEN
          dU_F(iNodeX,iX1,iX2,iX3,iCF) = NewtonFluidIncrement(iN_X,iCF)/dt
        ELSE
#endif
          dU_F(iNodeX,iX1,iX2,iX3,iCF) &
            = ( CF_N(iN_X,iCF) - U_F(iNodeX,iX1,iX2,iX3,iCF) ) / dt
#if defined(TWOMOMENT_ORDER_V)
        END IF
#endif
      ELSE
        dU_F(iNodeX,iX1,iX2,iX3,iCF) = Zero
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, dU_R, CR_N, U_R, &
    !$ACC          MatterCoupling_Mask_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#endif
    DO iS = 1, nSpecies
    DO iCR = 1, nCR
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE = iE_B0, iE_E0
    DO iNodeZ = 1, nDOFZ

      iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX   ) + 1
      iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

      iN_X = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
      iN_E = iNodeE &
             + ( iE  - iE_B0    ) * nDOFE

      IF ( MatterCoupling_Mask_N(iN_X) ) THEN
#if defined(TWOMOMENT_ORDER_V)
        IF( UseNewtonSolve )THEN
          dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = NewtonRadiationIncrement(iN_E,iS,iN_X,iCR)/dt
        ELSE
#endif
          dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
            = ( CR_N(iN_E,iS,iN_X,iCR) - U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) ) / dt
#if defined(TWOMOMENT_ORDER_V)
        END IF
#endif
      ELSE
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = Zero
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if defined( TWOMOMENT_ORDER_V )
    IF( RetainSolverConvergenceDiagnostics ) &
      CALL FinalizeConservationHandoffCheck &
             ( iZ_B1, iZ_E1, dt, U_F, U_R, dU_F, dU_R, &
               MatterCoupling_Mask_N )
#endif

  END SUBROUTINE ComputeAndMapIncrement


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE FinalizeConservationHandoffCheck &
    ( iZ_B1, iZ_E1, dt, U_F, U_R, dU_F, dU_R, MatterCoupling_Mask_N )

    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: dt
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      dU_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      dU_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: MatterCoupling_Mask_N(:)

    INTEGER :: iE, iN_E, iN_X, iNodeE, iNodeX, iNodeZ, iS
    INTEGER :: iSpatial, iX1, iX2, iX3
    REAL(DP) :: PhaseWeight2, PhaseWeight3
    REAL(DP) :: PrimitiveInitialL, PrimitiveInitialE
    REAL(DP) :: ConservedInitialL, ConservedInitialE
    REAL(DP) :: PrimitiveFinalL, PrimitiveFinalE
    REAL(DP) :: ConservedFinalL, ConservedFinalE
    REAL(DP) :: ReturnedFinalL, ReturnedFinalE
    REAL(DP) :: MatterChangeL, MatterChangeE
    REAL(DP) :: ReturnedMatterChangeL, ReturnedMatterChangeE
    REAL(DP) :: V_u_1_Old, V_u_2_Old, V_u_3_Old
    REAL(DP) :: V_d_1_Final, V_d_2_Final, V_d_3_Final
    REAL(DP) :: V_u_1_Final, V_u_2_Final, V_u_3_Final
    REAL(DP) :: V_u_1_Returned, V_u_2_Returned, V_u_3_Returned
    REAL(DP) :: DensityReturned, vDotI

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iSpatial, iNodeX, iX1, iX2, iX3, iE, iNodeE, &
    !$OMP          iNodeZ, iS, PhaseWeight2, PhaseWeight3, &
    !$OMP          PrimitiveInitialL, PrimitiveInitialE, &
    !$OMP          ConservedInitialL, ConservedInitialE, &
    !$OMP          PrimitiveFinalL, PrimitiveFinalE, &
    !$OMP          ConservedFinalL, ConservedFinalE, &
    !$OMP          ReturnedFinalL, ReturnedFinalE, &
    !$OMP          MatterChangeL, MatterChangeE, &
    !$OMP          ReturnedMatterChangeL, ReturnedMatterChangeE, &
    !$OMP          V_u_1_Old, V_u_2_Old, V_u_3_Old, &
    !$OMP          V_d_1_Final, V_d_2_Final, V_d_3_Final, &
    !$OMP          V_u_1_Final, V_u_2_Final, V_u_3_Final, &
    !$OMP          V_u_1_Returned, V_u_2_Returned, V_u_3_Returned, &
    !$OMP          DensityReturned, vDotI )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iSpatial, iNodeX, iX1, iX2, iX3, iE, iNodeE, &
    !$ACC          iNodeZ, iS, PhaseWeight2, PhaseWeight3, &
    !$ACC          PrimitiveInitialL, PrimitiveInitialE, &
    !$ACC          ConservedInitialL, ConservedInitialE, &
    !$ACC          PrimitiveFinalL, PrimitiveFinalE, &
    !$ACC          ConservedFinalL, ConservedFinalE, &
    !$ACC          ReturnedFinalL, ReturnedFinalE, &
    !$ACC          MatterChangeL, MatterChangeE, &
    !$ACC          ReturnedMatterChangeL, ReturnedMatterChangeE, &
    !$ACC          V_u_1_Old, V_u_2_Old, V_u_3_Old, &
    !$ACC          V_d_1_Final, V_d_2_Final, V_d_3_Final, &
    !$ACC          V_u_1_Final, V_u_2_Final, V_u_3_Final, &
    !$ACC          V_u_1_Returned, V_u_2_Returned, V_u_3_Returned, &
    !$ACC          DensityReturned, vDotI ) &
    !$ACC PRESENT( nX, GE_N, EnergyQuadratureWeight_N, GX_N, &
    !$ACC          CF_N, PF_N, CR_N, PR_N, U_F, U_R, &
    !$ACC          dU_F, dU_R, MatterCoupling_Mask_N, &
    !$ACC          ConservationHandoffCheck_Rate_Final )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iSpatial, iNodeX, iX1, iX2, iX3, iE, iNodeE, &
    !$OMP          iNodeZ, iS, PhaseWeight2, PhaseWeight3, &
    !$OMP          PrimitiveInitialL, PrimitiveInitialE, &
    !$OMP          ConservedInitialL, ConservedInitialE, &
    !$OMP          PrimitiveFinalL, PrimitiveFinalE, &
    !$OMP          ConservedFinalL, ConservedFinalE, &
    !$OMP          ReturnedFinalL, ReturnedFinalE, &
    !$OMP          MatterChangeL, MatterChangeE, &
    !$OMP          ReturnedMatterChangeL, ReturnedMatterChangeE, &
    !$OMP          V_u_1_Old, V_u_2_Old, V_u_3_Old, &
    !$OMP          V_d_1_Final, V_d_2_Final, V_d_3_Final, &
    !$OMP          V_u_1_Final, V_u_2_Final, V_u_3_Final, &
    !$OMP          V_u_1_Returned, V_u_2_Returned, V_u_3_Returned, &
    !$OMP          DensityReturned, vDotI )
#endif
    DO iN_X = 1, nX_G

      IF( .NOT. MatterCoupling_Mask_N(iN_X) )THEN
        ConservationHandoffCheck_Rate_Final(:,iN_X) = Zero
        CYCLE
      END IF

      iNodeX  = MOD( iN_X - 1, nDOFX ) + 1
      iSpatial = ( iN_X - 1 ) / nDOFX
      iX1 = iX_B0(1) + MOD( iSpatial, nX(1) )
      iSpatial = iSpatial / nX(1)
      iX2 = iX_B0(2) + MOD( iSpatial, nX(2) )
      iX3 = iX_B0(3) + iSpatial / nX(2)

      PrimitiveInitialL = &
        ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_L,iN_X)
      PrimitiveInitialE = &
        ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_E,iN_X)
      ConservedInitialL = Zero
      ConservedInitialE = Zero
      PrimitiveFinalL = Zero
      PrimitiveFinalE = Zero
      ConservedFinalL = Zero
      ConservedFinalE = Zero
      ReturnedFinalL = Zero
      ReturnedFinalE = Zero

      V_u_1_Old = U_F(iNodeX,iX1,iX2,iX3,iCF_S1) &
        / ( GX_N(iN_X,iGF_Gm_dd_11) &
            * U_F(iNodeX,iX1,iX2,iX3,iCF_D) * SpeedOfLight )
      V_u_2_Old = U_F(iNodeX,iX1,iX2,iX3,iCF_S2) &
        / ( GX_N(iN_X,iGF_Gm_dd_22) &
            * U_F(iNodeX,iX1,iX2,iX3,iCF_D) * SpeedOfLight )
      V_u_3_Old = U_F(iNodeX,iX1,iX2,iX3,iCF_S3) &
        / ( GX_N(iN_X,iGF_Gm_dd_33) &
            * U_F(iNodeX,iX1,iX2,iX3,iCF_D) * SpeedOfLight )
      V_u_1_Final = PF_N(iN_X,iPF_V1) / SpeedOfLight
      V_u_2_Final = PF_N(iN_X,iPF_V2) / SpeedOfLight
      V_u_3_Final = PF_N(iN_X,iPF_V3) / SpeedOfLight
      V_d_1_Final = GX_N(iN_X,iGF_Gm_dd_11) * PF_N(iN_X,iPF_V1) &
        / SpeedOfLight
      V_d_2_Final = GX_N(iN_X,iGF_Gm_dd_22) * PF_N(iN_X,iPF_V2) &
        / SpeedOfLight
      V_d_3_Final = GX_N(iN_X,iGF_Gm_dd_33) * PF_N(iN_X,iPF_V3) &
        / SpeedOfLight
      DensityReturned = U_F(iNodeX,iX1,iX2,iX3,iCF_D) &
        + dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_D)
      V_u_1_Returned = &
        ( U_F(iNodeX,iX1,iX2,iX3,iCF_S1) &
            + dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_S1) ) &
          / ( GX_N(iN_X,iGF_Gm_dd_11) * DensityReturned * SpeedOfLight )
      V_u_2_Returned = &
        ( U_F(iNodeX,iX1,iX2,iX3,iCF_S2) &
            + dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_S2) ) &
          / ( GX_N(iN_X,iGF_Gm_dd_22) * DensityReturned * SpeedOfLight )
      V_u_3_Returned = &
        ( U_F(iNodeX,iX1,iX2,iX3,iCF_S3) &
            + dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_S3) ) &
          / ( GX_N(iN_X,iGF_Gm_dd_33) * DensityReturned * SpeedOfLight )

      DO iS = 1, nSpecies
      DO iN_E = 1, nE_G

        iE     = ( iN_E - 1 ) / nDOFE + iE_B0
        iNodeE = MOD( iN_E - 1, nDOFE ) + 1
        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
        PhaseWeight2 = EnergyQuadratureWeight_N(iN_E) &
                         * GE_N(iN_E,iGE_Ep2)
        PhaseWeight3 = EnergyQuadratureWeight_N(iN_E) &
                         * GE_N(iN_E,iGE_Ep3)

        ConservedInitialL = ConservedInitialL &
          + PhaseWeight2 * LeptonNumber(iS) &
              * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)
        ConservedInitialE = ConservedInitialE + PhaseWeight3 &
          * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
              + V_u_1_Old * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) &
              + V_u_2_Old * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) &
              + V_u_3_Old * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS) )

        vDotI = V_d_1_Final * PR_N(iN_E,iS,iN_X,iPR_I1) &
              + V_d_2_Final * PR_N(iN_E,iS,iN_X,iPR_I2) &
              + V_d_3_Final * PR_N(iN_E,iS,iN_X,iPR_I3)
        PrimitiveFinalL = PrimitiveFinalL &
          + PhaseWeight2 * LeptonNumber(iS) &
              * ( PR_N(iN_E,iS,iN_X,iPR_D) + vDotI )
        PrimitiveFinalE = PrimitiveFinalE + PhaseWeight3 &
          * ( PR_N(iN_E,iS,iN_X,iPR_D) + 2.0_DP * vDotI )
        ConservedFinalL = ConservedFinalL &
          + PhaseWeight2 * LeptonNumber(iS) * CR_N(iN_E,iS,iN_X,iCR_N)
        ConservedFinalE = ConservedFinalE + PhaseWeight3 &
          * ( CR_N(iN_E,iS,iN_X,iCR_N) &
              + V_u_1_Final * CR_N(iN_E,iS,iN_X,iCR_G1) &
              + V_u_2_Final * CR_N(iN_E,iS,iN_X,iCR_G2) &
              + V_u_3_Final * CR_N(iN_E,iS,iN_X,iCR_G3) )
        ReturnedFinalL = ReturnedFinalL + PhaseWeight2 * LeptonNumber(iS) &
          * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
              + dt * dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) )
        ReturnedFinalE = ReturnedFinalE + PhaseWeight3 &
          * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
              + dt * dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
              + V_u_1_Returned &
                  * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) &
                      + dt * dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) ) &
              + V_u_2_Returned &
                  * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) &
                      + dt * dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) ) &
              + V_u_3_Returned &
                  * ( U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS) &
                      + dt * dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS) ) )

      END DO
      END DO

      MatterChangeL = CF_N(iN_X,iCF_Ne) &
        - U_F(iNodeX,iX1,iX2,iX3,iCF_Ne)
      MatterChangeE = CF_N(iN_X,iCF_E) &
        - U_F(iNodeX,iX1,iX2,iX3,iCF_E)
      ReturnedMatterChangeL = dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_Ne)
      ReturnedMatterChangeE = dt * dU_F(iNodeX,iX1,iX2,iX3,iCF_E)

      ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_L,iN_X) &
        = ( PrimitiveInitialL - ConservedInitialL ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_PrimitivePair_L,iN_X) &
        = ( MatterChangeL + PrimitiveFinalL - PrimitiveInitialL ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_FinalConversion_L,iN_X) &
        = ( ConservedFinalL - PrimitiveFinalL ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_MappedPair_L,iN_X) &
        = ( MatterChangeL + ConservedFinalL - ConservedInitialL ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_ReturnedPair_L,iN_X) &
        = ( ReturnedMatterChangeL + ReturnedFinalL - ConservedInitialL ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_InputConversion_E,iN_X) &
        = ( PrimitiveInitialE - ConservedInitialE ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_PrimitivePair_E,iN_X) &
        = ( MatterChangeE + PrimitiveFinalE - PrimitiveInitialE ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_FinalConversion_E,iN_X) &
        = ( ConservedFinalE - PrimitiveFinalE ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_MappedPair_E,iN_X) &
        = ( MatterChangeE + ConservedFinalE - ConservedInitialE ) / dt
      ConservationHandoffCheck_Rate_Final(iCHC_ReturnedPair_E,iN_X) &
        = ( ReturnedMatterChangeE + ReturnedFinalE - ConservedInitialE ) / dt

    END DO

  END SUBROUTINE FinalizeConservationHandoffCheck
#endif


  SUBROUTINE AbortCollisionDiscretization( Message )

    CHARACTER(LEN=*), INTENT(in) :: Message

    INTEGER :: ErrorCode

    WRITE(ERROR_UNIT,'(A)') &
      '[TwoMoment collision discretization] ERROR: ' // TRIM( Message )
    FLUSH(ERROR_UNIT)
    CALL MPI_ABORT( MPI_COMM_WORLD, -1, ErrorCode )

  END SUBROUTINE AbortCollisionDiscretization


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
