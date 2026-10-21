MODULE TwoMoment_NewtonSolverModule_OrderV

  USE KindModule, ONLY: DP
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE LinearAlgebraModule, ONLY: LinearFactorBatched, LinearApplyBatched
  USE TwoMoment_RetainedArithmeticModule
  USE TwoMoment_CollisionOperatorModule_OrderV
  USE TwoMoment_RadiationStateRecoveryModule_OrderV
  USE TwoMoment_TimersModule, ONLY: TimersStart, TimersStop, &
    Timer_Newton_Residual,Timer_Newton_Jacobian,Timer_Newton_LinearSolve,Timer_Newton_Search

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SolveNeutrinoMatterCoupling_Newton_OrderV_Batch

  ABSTRACT INTERFACE
    SUBROUTINE QueryCollisionCoefficients_OrderV_Batch &
      ( Mode, Indices, Requested, Transfer, Process, Species, QA, QB, &
        KT, Mu, Susceptibility, Chi, Sigma, Gain, Loss, Error )
      IMPORT :: DP, RetainedReal
      INTEGER, INTENT(in) :: Mode, Indices(:), Process(:), Species(:), QA(:), QB(:)
      LOGICAL, INTENT(in) :: Requested(:)
      TYPE(RetainedReal), INTENT(in) :: Transfer(:,:)
      REAL(DP), INTENT(out) :: KT(:), Mu(:,:), Susceptibility(:), Chi(:,:), Sigma(:,:)
      REAL(DP), INTENT(out) :: Gain(:,:), Loss(:,:)
      INTEGER, INTENT(out) :: Error(:)
    END SUBROUTINE QueryCollisionCoefficients_OrderV_Batch
  END INTERFACE

CONTAINS

  SUBROUTINE SolveNeutrinoMatterCoupling_Newton_OrderV_Batch &
    ( NE, NS, ND, Indices, Active, Enabled, Energy, WeightE, Interp, W, U, D0, I0, V0, Fluid0, &
      DT, Cap, Tolerance, MaxIterations, Query, X, NativeIncrement, Iterations, &
      RadiationError, ConservationError, CoefficientQueries, Failure, TraceUnit_Option, EmAbRate_Option, Check_Option, &
      FluidSourceTolerance_Option )

    INTEGER, INTENT(in) :: NE, NS, ND, Indices(:), Active(:), MaxIterations
    LOGICAL, INTENT(in) :: Enabled(9)
    REAL(DP), INTENT(in) :: Energy(NE), WeightE(ND), Interp(ND+2,ND)
    REAL(DP), INTENT(in) :: W(:,:), U(:,:), D0(:,:), I0(:,:,:), V0(:,:), Fluid0(:,:)
    REAL(DP), INTENT(in) :: DT, Cap, Tolerance
    PROCEDURE(QueryCollisionCoefficients_OrderV_Batch) :: Query
    TYPE(RetainedReal), INTENT(out), TARGET :: X(:,:), NativeIncrement(:,:,:)
    INTEGER, INTENT(out) :: Iterations(:), CoefficientQueries(:), Failure(:)
    REAL(DP), INTENT(out) :: RadiationError(:), ConservationError(:)
    INTEGER, INTENT(in), OPTIONAL :: TraceUnit_Option
    REAL(DP), INTENT(out), OPTIONAL, TARGET :: EmAbRate_Option(:,:)
    REAL(DP), INTENT(out), OPTIONAL, TARGET :: Check_Option(:,:,:)
    REAL(DP), INTENT(in), OPTIONAL :: FluidSourceTolerance_Option
    REAL(DP) :: FluidSourceTolerance
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: ProcessIncrement(:,:,:), ProcessMoment(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: TrialIncrement(:,:,:), TrialMoment(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: DirectionIncrement(:,:,:), DirectionMoment(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: ProcessNumber(:,:), TrialNumber(:,:), DirectionNumber(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: ProcessEnergy(:), TrialEnergy(:), DirectionEnergy(:)
    LOGICAL, ALLOCATABLE, TARGET :: SharedEndpointValid(:)
    INTEGER, ALLOCATABLE, TARGET :: SharedProjectionInfo(:)
    TYPE(IncomingRadiationCache), ALLOCATABLE, TARGET :: Incoming(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: IncomingEnergy(:),RecoveryResponse(:,:,:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: EnergyRatio(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: RealEnergyRatio(:,:)
    INTEGER, ALLOCATABLE, TARGET :: Components(:,:)
    LOGICAL, ALLOCATABLE, TARGET :: ComponentsValid(:),ContinueSearch(:)
    INTEGER, ALLOCATABLE, TARGET :: SpectralRows(:,:),BlockMembers(:),BlockFirst(:),BlockLength(:)
    INTEGER, ALLOCATABLE, TARGET :: BlockInfo(:)
    INTEGER, ALLOCATABLE :: GroupFirst(:),GroupCount(:),GroupSize(:),GroupLU(:),GroupRHS(:),GroupVector(:),GroupInfo(:)
    REAL(DP), ALLOCATABLE, TARGET :: SpectralLU(:),BlockRHS(:),BlockCorrection(:)
    INTEGER :: SpectralSize,NGroups,Group,Block,LUSize,RHSSize,VectorSize,InfoSize

    INTEGER :: M, N, NF, K, NA, NC, NX, NB, NV, B, J, Q, F, Loop, Cut, Pass, Ierr
    INTEGER :: Status, VendorStatus, Mode, Therm
    INTEGER :: Choice
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: NativeF(:,:), TransferG(:,:), TrialF(:,:), TrialG(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: RetainedStep(:,:), Velocity(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: Independent(:,:), TrialIndependent(:,:), NewtonStep(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: DirectionScale(:,:), Gradient(:,:), ConstraintJacobian(:,:,:)
    REAL(DP) :: ColumnWork(SIZE(X,1)), ScaleWork(SIZE(X,1)), VCWork(3)
    REAL(DP) :: ObjectiveScale, GradientLength, GradientNorm, ActionNorm
    REAL(DP), PARAMETER :: Multipliers(6)=[0.0_DP,1.0_DP,0.1_DP,0.01_DP,0.001_DP,0.0001_DP]
    INTEGER, ALLOCATABLE :: A(:), EB(:), Kind(:), Process(:), Species(:), QA(:), QB(:), FM(:,:)
    INTEGER, ALLOCATABLE, TARGET :: LiftRows(:,:)
    INTEGER, ALLOCATABLE, TARGET :: Info(:), Errors(:), Age(:), Accepted(:), FailedTrials(:)
    INTEGER, ALLOCATABLE, TARGET :: SpectralPivot(:), BorderPivot(:), SearchCounts(:,:)
    LOGICAL, ALLOCATABLE, TARGET :: Running(:), Requested(:), Refresh(:), ThermalNeeded(:)
    LOGICAL, ALLOCATABLE, TARGET :: Searching(:), AcceptedTrial(:), Fresh(:), NeedsConstraint(:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: Trial(:,:), Transfer(:,:), ProbeTransfer(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: Target(:,:), ProbeTarget(:,:), DTarget(:,:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: C0(:,:), C1(:,:,:), Source(:,:), KappaPair(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: NewC0(:,:), NewC1(:,:,:), NewSource(:,:), NewKappa(:,:)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: DRadiation(:,:), IRadiation(:,:,:), ScratchIncrement(:,:,:)
    REAL(DP), ALLOCATABLE, TARGET :: KT(:), Mu(:,:), Susceptibility(:), ProbeKT(:), ProbeMu(:,:), ProbeChiMu(:)
    REAL(DP), ALLOCATABLE, TARGET :: Chi(:,:), Sigma(:,:), Gain(:,:), Loss(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: ProbeChi(:,:), ProbeSigma(:,:), ProbeGain(:,:), ProbeLoss(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: DChi(:,:,:), DSigma(:,:,:), DGain(:,:,:), DLoss(:,:,:)
    REAL(DP), ALLOCATABLE, TARGET :: Residual(:,:), Physical(:,:), NewResidual(:,:), NewPhysical(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: HScale(:,:), GScale(:,:), NewHScale(:,:), NewGScale(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: AbsoluteSource(:,:), NewAbsolute(:,:), DensityLoss(:,:), NewDensityLoss(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: Kappa(:,:), CN(:,:,:), KD(:,:,:), SD(:,:,:), SI(:,:,:)
    REAL(DP), ALLOCATABLE, TARGET :: ThermC(:,:,:,:), ThermS(:,:,:), Z0(:,:,:), BInverse(:,:,:,:)
    REAL(DP), ALLOCATABLE, TARGET :: SourceValue(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: Matrix(:,:,:), RHS(:,:), OriginalRHS(:,:), Correction(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: RowScale(:,:), ColumnScale(:,:), Step(:,:), Action(:,:), Weights(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: InitialNorm(:,:),  Merit(:), Slope(:), Fraction(:), FDStep(:)
    REAL(DP), ALLOCATABLE, TARGET :: FluidSourceError(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: SpectralRHS(:,:,:), BorderLU(:,:,:), BorderRHS(:,:)
    REAL(DP) :: IValues(SIZE(D0,1),3), TempSource(NS/2+4)
    REAL(DP) :: Predicted, Ratio, TrialPhi
    TYPE(RetainedReal) :: Term, MuonTransfer

    M = NE*NS
    N = 4*M
    NF = NS/2
    K = NF+4
    NA = SIZE(Active)
    NC = M+NA+K
    NX = N+NA+K
    NB = SIZE(Indices)
    FluidSourceTolerance=Tolerance
    IF( PRESENT(FluidSourceTolerance_Option) ) FluidSourceTolerance=FluidSourceTolerance_Option
    IF( .NOT. IEEE_IS_FINITE(FluidSourceTolerance) .OR. FluidSourceTolerance<=0.0_DP )THEN
      Failure=9
      RETURN
    END IF
    IF( Cap/=1.0_DP )THEN
      Failure=9
      RETURN
    END IF
    CALL BuildCollisionEvents_OrderV(NE,NS,Enabled,A,EB,Kind,Process,Species,QA,QB,FM)
    NV = SIZE(A)
    CALL BuildSpectralBlocks(M,A,EB,SpectralRows,BlockMembers,BlockFirst,BlockLength)
    SpectralSize=SUM(BlockLength**2)
    ALLOCATE(GroupFirst(SIZE(BlockLength)),GroupCount(SIZE(BlockLength)),GroupSize(SIZE(BlockLength)), &
      GroupLU(SIZE(BlockLength)),GroupRHS(SIZE(BlockLength)),GroupVector(SIZE(BlockLength)),GroupInfo(SIZE(BlockLength)))
    NGroups=0
    Block=1
    LUSize=0
    RHSSize=0
    VectorSize=0
    InfoSize=0
    DO WHILE( Block<=SIZE(BlockLength) )
      NGroups=NGroups+1
      GroupFirst(NGroups)=Block
      GroupSize(NGroups)=BlockLength(Block)
      GroupCount(NGroups)=COUNT(BlockLength==BlockLength(Block))
      GroupLU(NGroups)=LUSize+1
      GroupRHS(NGroups)=RHSSize+1
      GroupVector(NGroups)=VectorSize+1
      GroupInfo(NGroups)=InfoSize+1
      LUSize=LUSize+GroupSize(NGroups)**2*GroupCount(NGroups)*NB
      RHSSize=RHSSize+GroupSize(NGroups)*(NA+K+1)*GroupCount(NGroups)*NB
      VectorSize=VectorSize+GroupSize(NGroups)*GroupCount(NGroups)*NB
      InfoSize=InfoSize+GroupCount(NGroups)*NB
      Block=Block+GroupCount(NGroups)
    END DO
    ALLOCATE(SpectralLU(LUSize),BlockRHS(RHSSize),BlockCorrection(VectorSize),BlockInfo(InfoSize))
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to:SpectralRows,BlockMembers,BlockFirst,BlockLength) &
    !$OMP MAP(alloc:SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(SpectralRows,BlockMembers,BlockFirst,BlockLength) &
    !$ACC CREATE(SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#endif
    ALLOCATE(EnergyRatio(M,NB),RealEnergyRatio(M,NB),Components(M,NB),ComponentsValid(NB),ContinueSearch(NB))
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(alloc:EnergyRatio,RealEnergyRatio,Components,ComponentsValid,ContinueSearch)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE(EnergyRatio,RealEnergyRatio,Components,ComponentsValid,ContinueSearch)
#endif
      ALLOCATE(ProcessIncrement(M,4,NB),ProcessMoment(K,NB),TrialIncrement(M,4,NB),TrialMoment(K,NB))
      ALLOCATE(DirectionIncrement(M,4,NB),DirectionMoment(K,NB))
      ALLOCATE(ProcessNumber(M,NB),TrialNumber(M,NB),DirectionNumber(M,NB))
      ALLOCATE(ProcessEnergy(NB),TrialEnergy(NB),DirectionEnergy(NB))
      ALLOCATE(SharedEndpointValid(NB))
      ALLOCATE(SharedProjectionInfo(NB))
      ALLOCATE(Incoming(M,NB),IncomingEnergy(NB),RecoveryResponse(4,4,M,NB))
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA MAP(alloc:ProcessIncrement,ProcessMoment,TrialIncrement,TrialMoment, &
      !$OMP   DirectionIncrement,DirectionMoment,ProcessNumber,TrialNumber,DirectionNumber, &
      !$OMP   ProcessEnergy,TrialEnergy,DirectionEnergy,SharedEndpointValid,SharedProjectionInfo, &
      !$OMP   Incoming,IncomingEnergy,RecoveryResponse)
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA CREATE(ProcessIncrement,ProcessMoment,TrialIncrement,TrialMoment, &
      !$ACC   DirectionIncrement,DirectionMoment,ProcessNumber,TrialNumber,DirectionNumber, &
      !$ACC   ProcessEnergy,TrialEnergy,DirectionEnergy,SharedEndpointValid,SharedProjectionInfo, &
      !$ACC   Incoming,IncomingEnergy,RecoveryResponse)
#endif
    ALLOCATE(NativeF(N,NB),TransferG(K,NB),TrialF(N,NB),TrialG(K,NB),RetainedStep(NX,NB),Velocity(3,NB))
    ALLOCATE(Independent(N,NB),TrialIndependent(N,NB),NewtonStep(NX,NB))
    ALLOCATE(Info(NB),Errors(NB),Age(NB),Accepted(NB),FailedTrials(NB))
    ALLOCATE(Running(NB),Requested(NB),Refresh(NB),Fresh(NB),Searching(NB),AcceptedTrial(NB))
    ALLOCATE(NeedsConstraint(NB),ThermalNeeded(3))
    ALLOCATE(Trial(NX,NB),Transfer(K,NB),ProbeTransfer(K,NB),Target(M,NB),ProbeTarget(M,NB),DTarget(M,3,NB))
    ALLOCATE(C0(M,NB),C1(M,3,NB),Source(K,NB),KappaPair(M,NB))
    ALLOCATE(NewC0(M,NB),NewC1(M,3,NB),NewSource(K,NB),NewKappa(M,NB))
    ALLOCATE(DRadiation(M,NB),IRadiation(M,3,NB),ScratchIncrement(M,4,NB))
    ALLOCATE(KT(NB),Mu(NS,NB),Susceptibility(NB),ProbeKT(NB),ProbeMu(NS,NB),ProbeChiMu(NB))
    ALLOCATE(Chi(M,NB),Sigma(M,NB),Gain(NV,NB),Loss(NV,NB))
    ALLOCATE(ProbeChi(M,NB),ProbeSigma(M,NB),ProbeGain(NV,NB),ProbeLoss(NV,NB))
    ALLOCATE(DChi(M,3,NB),DSigma(M,3,NB),DGain(NV,3,NB),DLoss(NV,3,NB))
    ALLOCATE(Residual(NX,NB),Physical(NX,NB),NewResidual(NX,NB),NewPhysical(NX,NB))
    ALLOCATE(HScale(K,NB),GScale(K,NB),NewHScale(K,NB),NewGScale(K,NB))
    ALLOCATE(AbsoluteSource(K,NB),NewAbsolute(K,NB),DensityLoss(M,NB),NewDensityLoss(M,NB))
    ALLOCATE(Kappa(M,NB),CN(SpectralSize,1,NB),KD(SpectralSize,1,NB),SD(K,M,NB),SI(1,M,NB))
    ALLOCATE(ThermC(M,4,3,NB),ThermS(K,3,NB),Z0(K,K,NB),LiftRows(K,NB),BInverse(3,3,M,NB))
    ALLOCATE(SourceValue(K,NB))
    ALLOCATE(Matrix(NC,NC,NB),RHS(NC,NB),OriginalRHS(NC,NB),Correction(NC,NB))
    ALLOCATE(RowScale(NC,NB),ColumnScale(NC,NB),Step(NX,NB),Action(NX,NB),Weights(NX,NB))
    ALLOCATE(InitialNorm(NS,NB),Merit(NB),Slope(NB),Fraction(NB),FDStep(NB))
    ALLOCATE(FluidSourceError(K,NB))
    ALLOCATE(SpectralRHS(M,NA+K+1,NB),BorderLU(NA+K,NA+K,NB),BorderRHS(NA+K,NB))
    ALLOCATE(SpectralPivot(M*NB),BorderPivot((NA+K)*NB))
    ALLOCATE(SearchCounts(6,NB))
    ThermalNeeded = .FALSE.
    DO J = 1, NA
      IF( Active(J) <= MIN(NF,2) ) ThermalNeeded(Active(J)) = .TRUE.
    END DO
    ThermalNeeded(3) = .TRUE.
    Failure = 0
    IF( PRESENT(Check_Option) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA MAP(alloc:Check_Option)
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA CREATE(Check_Option)
#endif
    END IF
    IF( PRESENT(EmAbRate_Option) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA MAP(alloc:EmAbRate_Option)
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA CREATE(EmAbRate_Option)
#endif
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to: A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices, &
    !$OMP   Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0,ThermalNeeded) &
    !$OMP MAP(alloc: Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh, &
    !$OMP   Searching,AcceptedTrial,NeedsConstraint,Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget, &
    !$OMP   C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa,DRadiation,IRadiation,ScratchIncrement, &
    !$OMP   KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss, &
    !$OMP   DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual,NewPhysical, &
    !$OMP   HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa, &
    !$OMP   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction, &
    !$OMP   RowScale,ColumnScale,Step,Action,Weights,InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$OMP   X,NativeIncrement, &
    !$OMP   Iterations,RadiationError,ConservationError,FluidSourceError,CoefficientQueries,Failure, &
    !$OMP   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$OMP   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep,SearchCounts)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices, &
    !$ACC   Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0,ThermalNeeded)
    !$ACC ENTER DATA CREATE(Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh, &
    !$ACC   Searching,AcceptedTrial,NeedsConstraint,Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget, &
    !$ACC   C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa,DRadiation,IRadiation,ScratchIncrement, &
    !$ACC   KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss, &
    !$ACC   DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual,NewPhysical, &
    !$ACC   HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa, &
    !$ACC   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction, &
    !$ACC   RowScale,ColumnScale,Step,Action,Weights,InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$ACC   X,NativeIncrement, &
    !$ACC   Iterations,RadiationError,ConservationError,FluidSourceError,CoefficientQueries,Failure, &
    !$ACC   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$ACC   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep,SearchCounts)
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(Q)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(Q)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO PRIVATE(Q)
#endif
    DO B = 1, NB
      ComponentsValid(B)=.FALSE.
      DO Q=1,M
        EnergyRatio(Q,B)=Retain(U(Q,B))/W(Q,B)
        RealEnergyRatio(Q,B)=U(Q,B)/W(Q,B)
      END DO
      X(:,B) = 0.0_DP
      Transfer(:,B) = 0.0_DP
      Failure(B) = 0
      Iterations(B) = 0
      CoefficientQueries(B) = 0
      SearchCounts(:,B) = 0
      InitialNorm(:,B) = 0.0_DP
      FluidSourceError(:,B) = HUGE(1.0_DP)
        ProcessIncrement(:,:,B)=0.0_DP
        ProcessMoment(:,B)=0.0_DP
        ProcessNumber(:,B)=0.0_DP
        ProcessEnergy(B)=0.0_DP
        SharedProjectionInfo(B)=0
        CALL PrepareIncomingRadiationCache(M,U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Incoming(:,B),IncomingEnergy(B))
      CALL SelectConservativeSupport(M,NF,W(:,B),U(:,B),D0(:,B),I0(:,:,B),Cap,Z0(:,:,B),Info(B), &
            LiftRows_Option=LiftRows(:,B))
    END DO
    CALL ReadInfo
    IF( ANY(Info /= 0) )THEN
      Failure = 1
      GOTO 900
    END IF

    ! Normalize to the incoming physical residual once. Trials retain their
    ! own current coefficient packets; no opacity package is lagged.
    Requested = .TRUE.
    CALL SendRequested
    CALL Query(0,Indices,Requested,Transfer,Process,Species,QA,QB,KT,Mu,Susceptibility,Chi,Sigma,Gain,Loss,Errors)
    CALL ReadErrors
    IF( ANY(Errors /= 0) )THEN
      Failure = 2
      GOTO 900
    END IF
    CALL MakeTargets(.FALSE.)
    CALL EvaluateCurrent(.TRUE.)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO B = 1, NB
      CoefficientQueries(B) = CoefficientQueries(B)+1
    END DO

    ! No presolve and no predictor. The incoming zero increment is already
    ! conservative. Shared increments use the coupled fixed-reference search;
    ! the original modes retain their independent-radiation search.
    Running=.TRUE.
    Refresh=.TRUE.
    Age=3
    Accepted=0
    CALL SendControl
    DO Loop=1,2*MaxIterations+2
      CALL CheckConvergence
      IF( .NOT. ANY(Running) ) EXIT
      IF( ANY(Accepted>=MaxIterations .AND. Running) ) EXIT
      CALL SetTransfers
      CALL RefreshDerivatives
      IF( ANY(Failure/=0) ) EXIT
      CALL AssembleAndSolve
      IF( ANY(Failure/=0) ) EXIT
      Searching=Running
      AcceptedTrial=.FALSE.
      NeedsConstraint=.FALSE.
      FailedTrials=0
      CALL SendSearch
      CALL TimersStart(Timer_Newton_Search)
      CALL SearchDirection(11,.FALSE.)
      ! At a reused Jacobian, two expensive trial packets trigger a refresh
      ! at the SAME iterate before any further fallback directions.
      DO B=1,NB
        NeedsConstraint(B)=Searching(B) .AND. (Fresh(B) .OR. FailedTrials(B)<2)
      END DO
      CALL SendSearch
      CALL ProjectedSearch
      CALL TimersStop(Timer_Newton_Search)
      DO B=1,NB
        IF( Searching(B) .AND. .NOT. Fresh(B) ) Refresh(B)=.TRUE.
      END DO
      CALL SendSearch
      IF( PRESENT(TraceUnit_Option) )THEN
        IF( TraceUnit_Option/=0 )THEN
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM(RadiationError,ConservationError,FluidSourceError,Accepted,CoefficientQueries,Fraction)
#elif defined(THORNADO_OACC)
          !$ACC UPDATE SELF(RadiationError,ConservationError,FluidSourceError,Accepted,CoefficientQueries,Fraction)
#endif
          DO B=1,NB
            WRITE(TraceUnit_Option,'(A,3I8,3ES23.14E3,I9)') 'cpn ',Loop,Indices(B),Accepted(B), &
              RadiationError(B),ConservationError(B),Fraction(B),CoefficientQueries(B)
            WRITE(TraceUnit_Option,'(A,2I8,*(ES23.14E3))') 'cpn_source ',Indices(B),Accepted(B), &
              MAXVAL(FluidSourceError(:,B)),FluidSourceError(:,B)
          END DO
          FLUSH(TraceUnit_Option)
        END IF
      END IF
      IF( ANY(Failure/=0) ) EXIT
    END DO
    CALL CheckConvergence
    ! Do not misreport an unfinished earlier point as the cause when another
    ! point has already reported a concrete failure.
    IF( ALL(Failure==0) )THEN
      WHERE( Running ) Failure=8
    END IF

900 CONTINUE
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM(Failure)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE SELF(Failure)
#endif
    CALL WriteFailedLinearization
    IF( PRESENT(Check_Option) )THEN
      IF( ALL(Failure == 0) )THEN
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(J)
#endif
        DO B = 1, NB
          Check_Option(:,1,B) = Physical(N+1:N+K,B)
          Check_Option(:,2,B) = 0.0_DP
          DO J = 1, NA
            Check_Option(Active(J),2,B) = Physical(N+K+J,B)
          END DO
          Check_Option(:,3,B) = RetainedValue(Source(:,B))
        END DO
      END IF
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(from:Check_Option)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA COPYOUT(Check_Option)
#endif
    END IF
    IF( PRESENT(EmAbRate_Option) )THEN
      IF( ALL(Failure == 0) )THEN
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(Q,Term)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(Q,Term)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(Q,Term)
#endif
        DO B = 1, NB
          DO Q = 1, M
            Term = Chi(Q,B)*(Target(Q,B)-(D0(Q,B)+X(Q,B)))
            EmAbRate_Option(Q,B) = RetainedValue(Term)
          END DO
        END DO
      END IF
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(from:EmAbRate_Option)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA COPYOUT(EmAbRate_Option)
#endif
    END IF
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP(from: X,NativeIncrement,Iterations,RadiationError,ConservationError,FluidSourceError, &
    !$OMP   CoefficientQueries,SearchCounts) &
    !$OMP MAP(release: A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices,Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0,ThermalNeeded, &
    !$OMP   Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh,Searching,AcceptedTrial,NeedsConstraint, &
    !$OMP   Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget,C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa, &
    !$OMP   DRadiation,IRadiation,ScratchIncrement,KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss, &
    !$OMP   ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual, &
    !$OMP   NewPhysical,HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa,CN,KD,SD,SI, &
    !$OMP   ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction,RowScale,ColumnScale,Step,Action,Weights, &
    !$OMP   InitialNorm,Merit,Slope,Fraction,FDStep,Failure, &
    !$OMP   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$OMP   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT(X,NativeIncrement,Iterations,RadiationError,ConservationError,FluidSourceError, &
    !$ACC   CoefficientQueries,SearchCounts)
    !$ACC EXIT DATA DELETE(A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices,Energy,WeightE,Interp, &
    !$ACC   W,U,D0,I0,V0,Fluid0,ThermalNeeded, &
    !$ACC   Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh,Searching,AcceptedTrial,NeedsConstraint, &
    !$ACC   Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget,C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa, &
    !$ACC   DRadiation,IRadiation,ScratchIncrement,KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss, &
    !$ACC   ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual, &
    !$ACC   NewPhysical,HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa,CN,KD,SD,SI, &
    !$ACC   ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction,RowScale,ColumnScale,Step,Action,Weights, &
    !$ACC   InitialNorm,Merit,Slope,Fraction,FDStep,Failure, &
    !$ACC   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$ACC   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep)
#endif

    IF( ALLOCATED(DirectionScale) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(release:DirectionScale,Gradient,ConstraintJacobian)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE(DirectionScale,Gradient,ConstraintJacobian)
#endif
    END IF

      IF( PRESENT(TraceUnit_Option) )THEN
        IF( TraceUnit_Option/=0 )THEN
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM(ProcessIncrement,ProcessMoment)
#elif defined(THORNADO_OACC)
          !$ACC UPDATE SELF(ProcessIncrement,ProcessMoment)
#endif
          IF( ALL(Failure==0) ) CALL WriteSharedTrace
        END IF
      END IF
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(release:ProcessIncrement,ProcessMoment,TrialIncrement,TrialMoment, &
      !$OMP   DirectionIncrement,DirectionMoment,ProcessNumber,TrialNumber,DirectionNumber, &
      !$OMP   ProcessEnergy,TrialEnergy,DirectionEnergy,SharedEndpointValid,SharedProjectionInfo, &
      !$OMP   Incoming,IncomingEnergy,RecoveryResponse)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE(ProcessIncrement,ProcessMoment,TrialIncrement,TrialMoment, &
      !$ACC   DirectionIncrement,DirectionMoment,ProcessNumber,TrialNumber,DirectionNumber, &
      !$ACC   ProcessEnergy,TrialEnergy,DirectionEnergy,SharedEndpointValid,SharedProjectionInfo, &
      !$ACC   Incoming,IncomingEnergy,RecoveryResponse)
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP(release:EnergyRatio,RealEnergyRatio,Components,ComponentsValid,ContinueSearch)
    !$OMP TARGET EXIT DATA MAP(release:SpectralRows,BlockMembers,BlockFirst,BlockLength, &
    !$OMP   SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE(EnergyRatio,RealEnergyRatio,Components,ComponentsValid,ContinueSearch)
    !$ACC EXIT DATA DELETE(SpectralRows,BlockMembers,BlockFirst,BlockLength, &
    !$ACC   SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#endif

    IF( PRESENT(TraceUnit_Option) )THEN
      IF( TraceUnit_Option /= 0 )THEN
        WRITE(TraceUnit_Option,'(A)') &
          '# cpn_counts point failure accepted_iterations ordinary_trial_queries projected_fallbacks' // &
          ' projected_trial_queries derivative_queries projection_attempts total_coefficient_queries gradient_builds'
        WRITE(TraceUnit_Option,'(A)') &
          '# total_coefficient_queries = 1 initial + ordinary_trial_queries' // &
          ' + projected_trial_queries + derivative_queries'
        DO B = 1, NB
          WRITE(TraceUnit_Option,'(A,10I10)') 'cpn_counts ',Indices(B),Failure(B),Iterations(B), &
            SearchCounts(1:5,B),CoefficientQueries(B),SearchCounts(6,B)
        END DO
        IF( NF >= 2 .AND. PRESENT(Check_Option) .AND. ALL(Failure == 0) )THEN
          WRITE(TraceUnit_Option,'(A)') &
            '# cpn_muon point collision_step_internal ym_in delta_ym_applied delta_ym_direct_source' // &
            ' transfer_residual conservation_residual radiation_relative_residual'
          WRITE(TraceUnit_Option,'(A)') &
            '# delta_ym_applied - delta_ym_direct_source = -transfer_residual;' // &
            ' conservation_residual = radiation_family_increment - transfer'
          DO B = 1, NB
            MuonTransfer = 0.0_DP
            DO J = 1, NA
              IF( Active(J) == 2 ) MuonTransfer = X(N+J,B)
            END DO
            WRITE(TraceUnit_Option,'(A,I10,7ES23.14E3)') 'cpn_muon ',Indices(B),DT,Fluid0(3,B), &
              -RetainedValue(MuonTransfer),-DT*Check_Option(2,3,B), &
              Check_Option(2,2,B),Check_Option(2,1,B),RadiationError(B)
          END DO
        END IF
        FLUSH(TraceUnit_Option)
      END IF
    END IF

  CONTAINS

    SUBROUTINE WriteSharedTrace
      TYPE(RetainedReal) :: Raw(K),Term,ZeroC0(M),ZeroC1(M,3),ZeroS(K),Mapped(M,4)
      REAL(DP) :: RR(NX),PP(NX),HS(K),GS(K),AA(K),ZeroZ(N,K),MappedError,MappedScale
      INTEGER :: B,Q,J,F,S
      WRITE(TraceUnit_Option,'(A)') &
        '# cpn_shared point moment paired_increment materialized_moment_minus_paired'
      WRITE(TraceUnit_Option,'(A)') &
        '# cpn_shared_kinematic point max_increment_discrepancy max_increment_magnitude'
      ZeroC0=0.0_DP
      ZeroC1=0.0_DP
      ZeroS=0.0_DP
      AA=0.0_DP
      ZeroZ=0.0_DP
      DO B=1,NB
        Raw=0.0_DP
        DO Q=1,M
          S=(Q-1)/NE+1
          F=(S+1)/2
          Raw(F)=Raw(F)+MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*(W(Q,B)*ProcessIncrement(Q,1,B))
          Raw(NF+1:)=Raw(NF+1:)+U(Q,B)*ProcessIncrement(Q,:,B)
        END DO
        DO J=1,K
          WRITE(TraceUnit_Option,'(A,2I10,2ES23.14E3)') 'cpn_shared ',Indices(B),J, &
            RetainedValue(ProcessMoment(J,B)),RetainedValue(Raw(J)-ProcessMoment(J,B))
        END DO
        CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B),ZeroZ, &
          0.0_DP,ZeroC0,ZeroC1,ZeroS,RR,PP,HS,GS,AA,Mapped, &
            LiftRows_Option=LiftRows(:,B))
        MappedError=MAXVAL(ABS(RetainedValue(Mapped-ProcessIncrement(:,:,B))))
        MappedScale=MAXVAL(ABS(RetainedValue(ProcessIncrement(:,:,B))))
        WRITE(TraceUnit_Option,'(A,I10,2ES23.14E3)') 'cpn_shared_kinematic ',Indices(B),MappedError,MappedScale
      END DO
    END SUBROUTINE WriteSharedTrace

    SUBROUTINE ReadInfo
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Info)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Info)
#endif
    END SUBROUTINE ReadInfo

    SUBROUTINE ReadErrors
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Errors)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Errors)
#endif
    END SUBROUTINE ReadErrors

    SUBROUTINE SendRequested
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(Requested)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(Requested)
#endif
    END SUBROUTINE SendRequested

    SUBROUTINE SendControl
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(Running,Refresh,Age,Accepted)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(Running,Refresh,Age,Accepted)
#endif
    END SUBROUTINE SendControl

    SUBROUTINE SendSearch
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(Searching,AcceptedTrial,NeedsConstraint,FailedTrials,Refresh)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(Searching,AcceptedTrial,NeedsConstraint,FailedTrials,Refresh)
#endif
    END SUBROUTINE SendSearch

    SUBROUTINE SetTransfers
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(J)
#endif
      DO B = 1, NB
        Transfer(:,B) = 0.0_DP
        DO J = 1, NA
          Transfer(Active(J),B) = X(N+J,B)
        END DO
      END DO
    END SUBROUTINE SetTransfers

    SUBROUTINE MakeTargets( Probing )
      LOGICAL, INTENT(in) :: Probing
      IF( Probing )THEN
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          IF( .NOT. Requested(B) .OR. Errors(B) /= 0 ) CYCLE
          CALL ComputeEquilibriumTarget_OrderV_Point(NE,NS,ND,Energy,WeightE,Interp,ProbeKT(B),ProbeMu(:,B), &
            ProbeTarget(:,B))
        END DO
      ELSE
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          CALL ComputeEquilibriumTarget_OrderV_Point(NE,NS,ND,Energy,WeightE,Interp,KT(B),Mu(:,B),Target(:,B))
        END DO
      END IF
    END SUBROUTINE MakeTargets

    SUBROUTINE EvaluateCurrent( SetInitial )
      LOGICAL, INTENT(in) :: SetInitial
      CALL TimersStart(Timer_Newton_Residual)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(J)
#endif
      DO B = 1, NB
        DRadiation(:,B) = D0(:,B)+X(1:M,B)
        DO J = 1, 3
          IRadiation(:,J,B) = I0(:,J,B)+X(J*M+1:(J+1)*M,B)
        END DO
        CALL ComputeCollisionRates_OrderV_Point(M,NF,A,EB,Kind,FM,W(:,B),U(:,B),Gain(:,B),Loss(:,B),Chi(:,B),Sigma(:,B), &
          Target(:,B),DRadiation(:,B),IRadiation(:,:,B),C0(:,B),C1(:,:,B),Source(:,B), &
          AbsoluteSource(:,B),KappaPair(:,B),DensityLoss(:,B),EnergyRatio_Option=EnergyRatio(:,B))
        Kappa(:,B) = RetainedValue(KappaPair(:,B))
        SourceValue(:,B) = RetainedValue(Source(:,B))
          CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B),Z0(:,:,B),DT, &
            C0(:,B),C1(:,:,B),Source(:,B),Residual(:,B),Physical(:,B),HScale(:,B),GScale(:,B), &
            AbsoluteSource(:,B),NativeIncrement(:,:,B),NativeF(:,B),TransferG(:,B), &
            CollisionIncrement_Option=ProcessIncrement(:,:,B),CollisionMoment_Option=ProcessMoment(:,B), &
            Incoming_Option=Incoming(:,B), &
            LiftRows_Option=LiftRows(:,B))
        Velocity(:,B)=Retain(V0(:,B))
        DO J=1,NA
          IF( Active(J)>NF+1 ) Velocity(Active(J)-NF-1,B)=Velocity(Active(J)-NF-1,B)-X(N+J,B)
        END DO
        CALL ComputeIndependentRadiationResidual(M,NF,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          NativeF(:,B),TransferG(:,B),Independent(:,B), &
            LiftRows_Option=LiftRows(:,B))
        CALL ComputeIndependentMetrics(NE,NS,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),X(:,B),DT,Physical(:,B),Independent(:,B), &
          Kappa(:,B),DensityLoss(:,B),HScale(:,B),GScale(:,B),InitialNorm(:,B),SetInitial, &
          RadiationError(B),ConservationError(B),Weights(:,B),Merit(B))
          CALL ComputeCoupledSourceMetrics(M,NF,Active,Fluid0(:,B),X(:,B),Physical(:,B),GScale(:,B), &
            Tolerance,FluidSourceTolerance,Independent(:,B),FluidSourceError(:,B),Weights(:,B),Merit(B))
      END DO
      CALL TimersStop(Timer_Newton_Residual)
    END SUBROUTINE EvaluateCurrent

    SUBROUTINE CheckConvergence
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO B=1,NB
        IF( RadiationError(B)<Tolerance .AND. MAXVAL(FluidSourceError(:,B))<FluidSourceTolerance &
              .AND. ConservationError(B)<1.0_DP ) Running(B)=.FALSE.
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Running,Accepted)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Running,Accepted)
#endif
    END SUBROUTINE CheckConvergence

    SUBROUTINE RefreshDerivatives
      INTEGER :: Attempt, Field
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO B = 1, NB
        Fresh(B) = Running(B) .AND. (Refresh(B) .OR. Age(B) >= 3)
        IF( Fresh(B) )THEN
          DChi(:,:,B) = 0.0_DP
          DSigma(:,:,B) = 0.0_DP
          DGain(:,:,B) = 0.0_DP
          DLoss(:,:,B) = 0.0_DP
          DTarget(:,:,B) = 0.0_DP
          Age(B) = 0
        END IF
        Refresh(B) = .FALSE.
      END DO
      DO Therm = 1, 3
        IF( .NOT. ThermalNeeded(Therm) ) CYCLE
        Field = Therm
        IF( Therm == 3 ) Field = NF+1
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          Requested(B) = Fresh(B)
          ProbeTransfer(:,B) = Transfer(:,B)
          IF( Therm == 1 ) FDStep(B) = 1.0e-5_DP*MAX(ABS(Fluid0(2,B)),0.01_DP)
          IF( Therm == 2 ) FDStep(B) = 1.0e-4_DP &
            *MAX(ABS(RetainedValue(Fluid0(3,B)-Transfer(2,B))),ABS(Susceptibility(B)),TINY(1.0_DP))
          IF( Therm == 3 ) FDStep(B) = 1.0e-5_DP*MAX(ABS(Fluid0(1,B)-1.0_DP),1.0e-4_DP)
          ProbeTransfer(Field,B) = ProbeTransfer(Field,B)+FDStep(B)
        END DO
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(Requested)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE SELF(Requested)
#endif
        IF( .NOT. ANY(Requested) ) CYCLE
        Mode = 0
        CALL Query(Mode,Indices,Requested,ProbeTransfer,Process,Species,QA,QB, &
          ProbeKT,ProbeMu,ProbeChiMu,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,Errors)
        CALL ReadErrors
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          IF( .NOT. Requested(B) ) CYCLE
          CoefficientQueries(B) = CoefficientQueries(B)+1
          SearchCounts(4,B) = SearchCounts(4,B)+1
        END DO
        IF( ANY(Errors /= 0 .AND. Requested) )THEN
          ! One-sided EOS-domain probe, as in the checked reference.
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
          !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
          !$OMP PARALLEL DO
#endif
          DO B = 1, NB
            IF( .NOT. Requested(B) .OR. Errors(B) == 0 ) CYCLE
            FDStep(B) = -FDStep(B)
            ProbeTransfer(Field,B) = Transfer(Field,B)+FDStep(B)
          END DO
          CALL Query(Mode,Indices,Requested,ProbeTransfer,Process,Species,QA,QB, &
            ProbeKT,ProbeMu,ProbeChiMu,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,Errors)
          CALL ReadErrors
          ! The retry requests a second complete coefficient packet for every
          ! point in Requested, including points whose first query succeeded.
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
          !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
          !$OMP PARALLEL DO
#endif
          DO B = 1, NB
            IF( .NOT. Requested(B) ) CYCLE
            CoefficientQueries(B) = CoefficientQueries(B)+1
            SearchCounts(4,B) = SearchCounts(4,B)+1
          END DO
          IF( ANY(Errors /= 0 .AND. Requested) )THEN
            Failure = 5
            RETURN
          END IF
        END IF
        CALL MakeTargets(.TRUE.)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          IF( .NOT. Requested(B) ) CYCLE
          DTarget(:,Therm,B) = (ProbeTarget(:,B)-Target(:,B))/FDStep(B)
          DChi(:,Therm,B) = (ProbeChi(:,B)-Chi(:,B))/FDStep(B)
          DSigma(:,Therm,B) = (ProbeSigma(:,B)-Sigma(:,B))/FDStep(B)
          DGain(:,Therm,B) = (ProbeGain(:,B)-Gain(:,B))/FDStep(B)
          DLoss(:,Therm,B) = (ProbeLoss(:,B)-Loss(:,B))/FDStep(B)
        END DO
      END DO
    END SUBROUTINE RefreshDerivatives

    SUBROUTINE AssembleAndSolve
      INTEGER :: Refine
      CALL TimersStart(Timer_Newton_Jacobian)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !$OMP PRIVATE(IValues,TempSource)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE(IValues,TempSource)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(IValues,TempSource)
#endif
      DO B = 1, NB
        Matrix(:,:,B) = 0.0_DP
        RHS(:,B) = 0.0_DP
        Info(B) = 0
        IF( .NOT. Running(B) )THEN
          DO J = 1, NC
            Matrix(J,J,B) = 1.0_DP
          END DO
          OriginalRHS(:,B) = 0.0_DP
          RowScale(:,B) = 1.0_DP
          ColumnScale(:,B) = 1.0_DP
          CYCLE
        END IF
        DRadiation(:,B) = D0(:,B)+X(1:M,B)
        DO J = 1, 3
          IRadiation(:,J,B) = I0(:,J,B)+X(J*M+1:(J+1)*M,B)
        END DO
        IValues = RetainedValue(IRadiation(:,:,B))
        CALL ComputeCollisionRatesAndDerivatives_OrderV_Point(M,NF,A,EB,Kind,FM,W(:,B),U(:,B),Gain(:,B),Loss(:,B),Chi(:,B), &
          DRadiation(:,B),IRadiation(:,:,B),IValues,Kappa(:,B),DGain(:,:,B),DLoss(:,:,B),DChi(:,:,B),DSigma(:,:,B), &
          Target(:,B),DTarget(:,:,B),ThermalNeeded,CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
          ThermC(:,:,:,B),ThermS(:,:,B),ScratchIncrement(:,1:3,B),NewC1(:,:,B), &
          EnergyRatio_Option=EnergyRatio(:,B),RealEnergyRatio_Option=RealEnergyRatio(:,B), &
          SpectralRows_Option=SpectralRows)
        TempSource = RetainedValue(Source(:,B))
        CALL AssembleCollisionCondensedSystem(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
          Z0(:,:,B),DT,TempSource,Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B),ThermC(:,:,:,B), &
          ThermS(:,:,B),Residual(:,B),Matrix(:,:,B),RHS(:,B),BInverse(:,:,:,B),RowScale(:,B),ColumnScale(:,B),Info(B), &
          RetainedValue(TransferG(:,B)), &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
        CALL PrepareSharedCollisionLinearSystem(M,NF,NA,Matrix(:,:,B),RHS(:,B))
        OriginalRHS(:,B) = RHS(:,B)
      END DO
      CALL ReadInfo
      CALL TimersStop(Timer_Newton_Jacobian)
      IF( ANY(Info /= 0) )THEN
        Failure = 6
        RETURN
      END IF
      CALL TimersStart(Timer_Newton_LinearSolve)
      CALL SolveBorderedSystem
      CALL TimersStop(Timer_Newton_LinearSolve)
      IF( ANY(Failure /= 0) ) RETURN
      CALL TimersStart(Timer_Newton_Jacobian)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(TempSource,VCWork,J)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(TempSource,VCWork,J)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(TempSource,VCWork,J)
#endif
      DO B = 1, NB
        Step(:,B) = 0.0_DP
        IF( .NOT. Running(B) ) CYCLE
        TempSource = RetainedValue(Source(:,B))
        CALL ExpandCollisionCondensedStep_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B),Z0(:,:,B), &
          DT,TempSource,Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B),ThermC(:,:,:,B),ThermS(:,:,B), &
          Residual(:,B),BInverse(:,:,:,B),RHS(:,B)/ColumnScale(:,B),Step(:,B),Action(:,B), &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
        Step(N+NA+1:NX,B)=0.0_DP
        NewtonStep(:,B)=Step(:,B)
        RetainedStep(:,B)=Retain(Step(:,B))
        ! Expansion already returned this action with lift directions disabled.
        VCWork=0.0_DP
        DO J=1,NA
          IF( Active(J)>NF+1 ) VCWork(Active(J)-NF-1)=-Step(N+J,B)
        END DO
        CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          TransferG(:,B),VCWork,Action(:,B), &
            LiftRows_Option=LiftRows(:,B))
        CALL PrepareCoupledSearchWeights(M,NF,Active,Fluid0(:,B),X(:,B), Physical(:,B),Independent(:,B),Step(:,B),Action(:,B),Tolerance,FluidSourceTolerance,Weights(:,B),Merit(B))
        Slope(B)=ComputeCoupledSlope(N,Weights(:,B),Independent(:,B),Physical(:,B),Action(:,B))
      END DO
      CALL TimersStop(Timer_Newton_Jacobian)
      CALL WriteLinearTrace('newton')
    END SUBROUTINE AssembleAndSolve

    SUBROUTINE SolveBorderedSystem
      INTEGER :: Refine, BorderSize
      BorderSize = NA+K
      CALL FactorSpectralBlocks
      IF( ANY(Failure/=0) ) RETURN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO B = 1, NB
        BorderLU(:,:,B) = Matrix(M+1:NC,M+1:NC,B) &
          -MATMUL(Matrix(M+1:NC,1:M,B),SpectralRHS(:,2:BorderSize+1,B))
        BorderRHS(:,B) = RHS(M+1:NC,B)-MATMUL(Matrix(M+1:NC,1:M,B),SpectralRHS(:,1,B))
      END DO
      CALL LinearFactorBatched(BorderSize,BorderLU(1,1,1),BorderSize,BorderPivot,Info,NB,VendorStatus)
      CALL ReadInfo
      IF( VendorStatus /= 0 .OR. ANY(Info /= 0) )THEN
        Failure = 6
        RETURN
      END IF
      CALL LinearApplyBatched('N',BorderSize,1,BorderLU(1,1,1),BorderSize,BorderPivot,BorderRHS(1,1),BorderSize,NB,VendorStatus)
      IF( VendorStatus /= 0 )THEN
        Failure = 6
        RETURN
      END IF
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO B = 1, NB
        RHS(M+1:NC,B) = BorderRHS(:,B)
        RHS(1:M,B) = SpectralRHS(:,1,B)-MATMUL(SpectralRHS(:,2:BorderSize+1,B),BorderRHS(:,B))
      END DO
      DO Refine = 1, 2
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          Correction(:,B) = OriginalRHS(:,B)-MATMUL(Matrix(:,:,B),RHS(:,B))
        END DO
        CALL ApplySpectralCorrection
        IF( ANY(Failure/=0) ) RETURN
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          BorderRHS(:,B) = Correction(M+1:NC,B)-MATMUL(Matrix(M+1:NC,1:M,B),Correction(1:M,B))
        END DO
        CALL LinearApplyBatched('N',BorderSize,1,BorderLU(1,1,1),BorderSize,BorderPivot,BorderRHS(1,1),BorderSize,NB,VendorStatus)
        IF( VendorStatus /= 0 )THEN
          Failure = 6
          RETURN
        END IF
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          Correction(1:M,B) = Correction(1:M,B)-MATMUL(SpectralRHS(:,2:BorderSize+1,B),BorderRHS(:,B))
          Correction(M+1:NC,B) = BorderRHS(:,B)
          RHS(:,B) = RHS(:,B)+Correction(:,B)
        END DO
      END DO
      ! CN and KD remain unmodified in their compact block-row storage.
      ! Both refinement passes still use the original full scaled matrix.
    END SUBROUTINE SolveBorderedSystem

    SUBROUTINE FactorSpectralBlocks
      INTEGER :: G,BS,GF,GN,LB,RB,PB,IB,TT,BB,Block,II,JJ,Row,Col,Offset
      DO G=1,NGroups
        BS=GroupSize(G)
        GF=GroupFirst(G)
        GN=GroupCount(G)*NB
        LB=GroupLU(G)
        RB=GroupRHS(G)
        PB=GroupVector(G)
        IB=GroupInfo(G)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(BB,Block,II,JJ,Row,Col,Offset)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(BB,Block,II,JJ,Row,Col,Offset)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(BB,Block,II,JJ,Row,Col,Offset)
#endif
        DO TT=1,GN
          BB=MOD(TT-1,NB)+1
          Block=GF+(TT-1)/NB
          Offset=BlockFirst(Block)-1
          DO JJ=1,BS
            Col=BlockMembers(Offset+JJ)
            DO II=1,BS
              Row=BlockMembers(Offset+II)
              SpectralLU(LB+(TT-1)*BS*BS+(JJ-1)*BS+II-1)=Matrix(Row,Col,BB)
            END DO
          END DO
          DO II=1,BS
            Row=BlockMembers(Offset+II)
            BlockRHS(RB+(TT-1)*BS*(NA+K+1)+II-1)=RHS(Row,BB)
            DO JJ=1,NA+K
              BlockRHS(RB+(TT-1)*BS*(NA+K+1)+JJ*BS+II-1)=Matrix(Row,M+JJ,BB)
            END DO
          END DO
        END DO
        CALL LinearFactorBatched(BS,SpectralLU(LB),BS,SpectralPivot(PB),BlockInfo(IB),GN,VendorStatus)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(BlockInfo)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE SELF(BlockInfo)
#endif
        IF( VendorStatus/=0 .OR. ANY(BlockInfo(IB:IB+GN-1)/=0) )THEN
          Failure=6
          RETURN
        END IF
        CALL LinearApplyBatched('N',BS,NA+K+1,SpectralLU(LB),BS,SpectralPivot(PB), &
          BlockRHS(RB),BS,GN,VendorStatus)
        IF( VendorStatus/=0 )THEN
          Failure=6
          RETURN
        END IF
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(BB,Block,II,JJ,Row,Offset)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(BB,Block,II,JJ,Row,Offset)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(BB,Block,II,JJ,Row,Offset)
#endif
        DO TT=1,GN
          BB=MOD(TT-1,NB)+1
          Block=GF+(TT-1)/NB
          Offset=BlockFirst(Block)-1
          DO JJ=1,NA+K+1
          DO II=1,BS
            Row=BlockMembers(Offset+II)
            SpectralRHS(Row,JJ,BB)=BlockRHS(RB+(TT-1)*BS*(NA+K+1)+(JJ-1)*BS+II-1)
          END DO
          END DO
        END DO
      END DO
    END SUBROUTINE FactorSpectralBlocks

    SUBROUTINE ApplySpectralCorrection
      INTEGER :: G,BS,GF,GN,LB,PB,TT,BB,Block,II,Row,Offset
      DO G=1,NGroups
        BS=GroupSize(G)
        GF=GroupFirst(G)
        GN=GroupCount(G)*NB
        LB=GroupLU(G)
        PB=GroupVector(G)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(BB,Block,II,Row,Offset)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(BB,Block,II,Row,Offset)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(BB,Block,II,Row,Offset)
#endif
        DO TT=1,GN
          BB=MOD(TT-1,NB)+1
          Block=GF+(TT-1)/NB
          Offset=BlockFirst(Block)-1
          DO II=1,BS
            Row=BlockMembers(Offset+II)
            BlockCorrection(PB+(TT-1)*BS+II-1)=Correction(Row,BB)
          END DO
        END DO
        CALL LinearApplyBatched('N',BS,1,SpectralLU(LB),BS,SpectralPivot(PB),BlockCorrection(PB),BS,GN,VendorStatus)
        IF( VendorStatus/=0 )THEN
          Failure=6
          RETURN
        END IF
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(BB,Block,II,Row,Offset)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(BB,Block,II,Row,Offset)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(BB,Block,II,Row,Offset)
#endif
        DO TT=1,GN
          BB=MOD(TT-1,NB)+1
          Block=GF+(TT-1)/NB
          Offset=BlockFirst(Block)-1
          DO II=1,BS
            Row=BlockMembers(Offset+II)
            Correction(Row,BB)=BlockCorrection(PB+(TT-1)*BS+II-1)
          END DO
        END DO
      END DO
    END SUBROUTINE ApplySpectralCorrection

    SUBROUTINE MakeProcessDirection(Projected)
      LOGICAL, INTENT(in) :: Projected
      TYPE(RetainedReal) :: DD(M),DI(M,3),DC(M),DS(M),DTar(M),TD(3),TC(K),Factor
      INTEGER :: F,J
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(DD,DI,DC,DS,DTar,TD,TC,Factor,F,J,VCWork)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(DD,DI,DC,DS,DTar,TD,TC,Factor,F,J,VCWork)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(DD,DI,DC,DS,DTar,TD,TC,Factor,F,J,VCWork)
#endif
      DO B=1,NB
        IF( .NOT. Searching(B) .OR. AcceptedTrial(B) ) CYCLE
        TC=0.0_DP
        DO J=1,NA
          TC(Active(J))=RetainedStep(N+J,B)
        END DO
        TD=0.0_DP
        TD(1:MIN(NF,2))=TC(1:MIN(NF,2))
        TD(3)=TC(NF+1)
        DO J=1,3
          TD(3)=TD(3)-Velocity(J,B)*TC(NF+1+J)
        END DO
        DD=RetainedStep(1:M,B)
        DRadiation(:,B)=D0(:,B)+X(1:M,B)
        DO J=1,3
          DI(:,J)=RetainedStep(J*M+1:(J+1)*M,B)
          IRadiation(:,J,B)=I0(:,J,B)+X(J*M+1:(J+1)*M,B)
        END DO
        DC=0.0_DP
        DS=0.0_DP
        DTar=0.0_DP
        DO F=1,3
          DC=DC+TD(F)*DChi(:,F,B)
          DS=DS+TD(F)*DSigma(:,F,B)
          DTar=DTar+TD(F)*DTarget(:,F,B)
        END DO
        CALL ComputeCollisionRateDirection_OrderV_Point(M,NF,A,EB,Kind,FM,W(:,B),U(:,B),Gain(:,B),Loss(:,B), &
          Chi(:,B),Sigma(:,B),Target(:,B),DRadiation(:,B),IRadiation(:,:,B),DD,DI, &
          DGain(:,:,B),DLoss(:,:,B),TD,DC,DS,DTar,DirectionIncrement(:,1,B), &
          DirectionIncrement(:,2:4,B),DirectionMoment(:,B),NumberSource_Option=DirectionNumber(:,B), &
          NonScatteringEnergy_Option=DirectionEnergy(B),EnergyRatio_Option=EnergyRatio(:,B))
        DirectionIncrement(:,:,B)=DT*DirectionIncrement(:,:,B)
        DirectionMoment(:,B)=DT*DirectionMoment(:,B)
        DirectionNumber(:,B)=DT*DirectionNumber(:,B)
        DirectionEnergy(B)=DT*DirectionEnergy(B)
        ! TrialEnergy is overwritten by every line-search trial. Until then,
        ! use it to report the muon-family change introduced by projection,
        ! separately from the linearized physical source.
        TrialEnergy(B)=0.0_DP
        IF( NF>=2 ) TrialEnergy(B)=DirectionMoment(2,B)
        IF( .NOT. ComponentsValid(B) )THEN
          CALL BuildScatteringComponents(M,NF,A,EB,Kind,Gain(:,B),Loss(:,B),Components(:,B))
          ComponentsValid(B)=.TRUE.
        END IF
        CALL ProjectSharedProcessEndpoint(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Cap, &
          A,EB,Kind,Gain(:,B),Loss(:,B),Chi(:,B),Sigma(:,B),X(:,B),DirectionIncrement(:,:,B),DirectionMoment(:,B), &
          Trial(:,B),Info(B),NumberIncrement_Option=DirectionNumber(:,B),NonScatteringEnergy_Option=DirectionEnergy(B), &
          Incoming_Option=Incoming(:,B),IncomingEnergy_Option=IncomingEnergy(B), &
            Response_Option=RecoveryResponse(:,:,:,B),EnergyRatio_Option=EnergyRatio(:,B), &
            Components_Option=Components(:,B))
        IF( NF>=2 ) TrialEnergy(B)=DirectionMoment(2,B)-TrialEnergy(B)
        SharedEndpointValid(B)=Info(B)==0
        SharedProjectionInfo(B)=Info(B)
        ScratchIncrement(:,:,B)=DirectionIncrement(:,:,B)-ProcessIncrement(:,:,B)
        TrialMoment(:,B)=DirectionMoment(:,B)-ProcessMoment(:,B)
        CALL ComputeSharedProcessTangent(M,NF,Active,U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
          ScratchIncrement(:,:,B),TrialMoment(:,B),RetainedStep(:,B),Info(B), &
            Response_Option=RecoveryResponse(:,:,:,B))
        IF( Info(B)/=0 )THEN
          Slope(B)=0.0_DP
          Fraction(B)=0.0_DP
          CYCLE
        END IF
        ! Globalize the actual paired candidate curve, including a projected
        ! fallback. Its tangent need not equal the provisional linear direction.
        Step(:,B)=RetainedValue(RetainedStep(:,B))
        CALL ComputeCollisionLinearAction_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B),Z0(:,:,B), &
          DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B),ThermC(:,:,:,B),ThermS(:,:,B), &
          Step(:,B),.FALSE.,Action(:,B), &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
        VCWork=0.0_DP
        DO J=1,NA
          IF( Active(J)>NF+1 ) VCWork(Active(J)-NF-1)=-Step(N+J,B)
        END DO
        CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          TransferG(:,B),VCWork,Action(:,B), &
            LiftRows_Option=LiftRows(:,B))
        Slope(B)=ComputeCoupledSlope(N,Weights(:,B),Independent(:,B),Physical(:,B),Action(:,B))
        ! This is a curved kinematic candidate, not an affine primitive step.
        ! A tangent cone bound can discard a valid full endpoint in an empty
        ! bin. Check the recovered endpoint before any thermodynamic query.
        Fraction(B)=1.0_DP
      END DO
      CALL WriteLinearTrace('process-curve')
      IF( PRESENT(TraceUnit_Option) )THEN
        IF( TraceUnit_Option/=0 )THEN
          CALL ReadInfo
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM(SharedProjectionInfo,TrialEnergy,DirectionMoment,SourceValue)
#elif defined(THORNADO_OACC)
          !$ACC UPDATE SELF(SharedProjectionInfo,TrialEnergy,DirectionMoment,SourceValue)
#endif
          DO B=1,NB
            IF( Searching(B) .AND. .NOT. AcceptedTrial(B) .AND. NF>=2 ) &
              WRITE(TraceUnit_Option,'(A,2I10,3ES23.14E3)') 'cpn_shared_muon_projection ',Indices(B),Accepted(B), &
                RetainedValue(TrialEnergy(B)),RetainedValue(DirectionMoment(2,B)),DT*SourceValue(2,B)
            IF( SharedProjectionInfo(B)==0 ) CYCLE
            WRITE(TraceUnit_Option,'(A,2I10)') 'cpn_shared_projection ',Indices(B),SharedProjectionInfo(B)
          END DO
          IF( .NOT. Projected .AND. Indices(1)==1 .AND. Accepted(1)==0 )THEN
#if defined(THORNADO_OMP_OL)
            !$OMP TARGET UPDATE FROM(DirectionIncrement,DirectionMoment,Trial,D0,I0,W,U,V0,Chi,Sigma)
#elif defined(THORNADO_OACC)
            !$ACC UPDATE SELF(DirectionIncrement,DirectionMoment,Trial,D0,I0,W,U,V0,Chi,Sigma)
#endif
            WRITE(TraceUnit_Option,'(A,*(ES23.14E3))') 'cpn_shared_packet_velocity ',V0(:,1)
            WRITE(TraceUnit_Option,'(A,*(ES23.14E3))') 'cpn_shared_packet_moments ',RetainedValue(DirectionMoment(:,1))
            DO J=1,M
              WRITE(TraceUnit_Option,'(A,I6,14ES23.14E3)') 'cpn_shared_packet_bin ',J,W(J,1),U(J,1),D0(J,1),I0(J,:,1), &
                RetainedValue(DirectionIncrement(J,:,1)),Chi(J,1),Sigma(J,1), &
                RetainedValue(D0(J,1)+Trial(J,1)), &
                SQRT(SUM(RetainedValue(I0(J,:,1)+Trial([M+J,2*M+J,3*M+J],1))**2))
            END DO
          END IF
        END IF
      END IF
    END SUBROUTINE MakeProcessDirection

    SUBROUTINE SearchDirection( Cuts, Projected )
      INTEGER, INTENT(in) :: Cuts
      LOGICAL, INTENT(in) :: Projected
      INTEGER :: L
      CALL MakeProcessDirection(Projected)
      DO L = 1, Cuts
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(J)
#endif
        DO B = 1, NB
          ContinueSearch(B)=Searching(B) .AND. .NOT. AcceptedTrial(B)
          Requested(B) = Searching(B) .AND. .NOT. AcceptedTrial(B) .AND. Fraction(B) >= 2.0_DP**(-18) &
            .AND. IEEE_IS_FINITE(Slope(B)) &
            .AND. (Slope(B)<0.0_DP .OR. (.TRUE. .AND. L==1))
          Info(B)=0
          IF( .NOT. Requested(B) )THEN
            IF( Searching(B) ) Fraction(B)=0.5_DP*Fraction(B)
            CYCLE
          END IF
          IF( Requested(B) )THEN
              TrialIncrement(:,:,B)=(1.0_DP-Fraction(B))*ProcessIncrement(:,:,B) &
                +Fraction(B)*DirectionIncrement(:,:,B)
              TrialMoment(:,B)=(1.0_DP-Fraction(B))*ProcessMoment(:,B)+Fraction(B)*DirectionMoment(:,B)
              TrialNumber(:,B)=(1.0_DP-Fraction(B))*ProcessNumber(:,B)+Fraction(B)*DirectionNumber(:,B)
              TrialEnergy(B)=(1.0_DP-Fraction(B))*ProcessEnergy(B)+Fraction(B)*DirectionEnergy(B)
              IF( L/=1 .OR. Fraction(B)/=1.0_DP .OR. .NOT. SharedEndpointValid(B) )THEN
                Trial(:,B)=X(:,B)
                CALL RecoverSharedProcessState(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Cap, &
                  TrialIncrement(:,:,B),TrialMoment(:,B),Trial(:,B),Info(B),Incoming(:,B),IncomingEnergy(B), &
            Response_Option=RecoveryResponse(:,:,:,B))
              END IF
            IF( Info(B)/=0 ) Requested(B)=.FALSE.
          END IF
          IF( .NOT. CollisionStateAdmissible_OrderV_Point(M,D0(:,B),I0(:,:,B),Trial(:,B),Cap) ) Requested(B) = .FALSE.
          IF( Searching(B) .AND. .NOT. Requested(B) ) Fraction(B)=0.5_DP*Fraction(B)
          ProbeTransfer(:,B) = 0.0_DP
          DO J = 1, NA
            ProbeTransfer(Active(J),B) = Trial(N+J,B)
          END DO
        END DO
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(Requested,ContinueSearch)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE SELF(Requested,ContinueSearch)
#endif
        IF( .NOT. ANY(ContinueSearch) ) EXIT
        IF( .TRUE. .AND. L==1 .AND. PRESENT(TraceUnit_Option) )THEN
          IF( TraceUnit_Option/=0 )THEN
#if defined(THORNADO_OMP_OL)
            !$OMP TARGET UPDATE FROM(Info,Slope,Fraction)
#elif defined(THORNADO_OACC)
            !$ACC UPDATE SELF(Info,Slope,Fraction)
#endif
            DO B=1,NB
              WRITE(TraceUnit_Option,'(A,2I10,2ES23.14E3)') &
                'cpn_shared_trial ',Indices(B),Info(B),Slope(B),Fraction(B)
            END DO
          END IF
        END IF
        IF( .NOT. ANY(Requested) ) CYCLE
        Mode = 0
        CALL Query(Mode,Indices,Requested,ProbeTransfer,Process,Species,QA,QB, &
          ProbeKT,ProbeMu,ProbeChiMu,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,Errors)
        CALL ReadErrors
        CALL MakeTargets(.TRUE.)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J,Predicted,Ratio,TrialPhi)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J,Predicted,Ratio,TrialPhi)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(J,Predicted,Ratio,TrialPhi)
#endif
        DO B = 1, NB
          IF( .NOT. Requested(B) ) CYCLE
          CoefficientQueries(B) = CoefficientQueries(B)+1
          IF( Projected )THEN
            SearchCounts(3,B) = SearchCounts(3,B)+1
          ELSE
            SearchCounts(1,B) = SearchCounts(1,B)+1
          END IF
          FailedTrials(B) = FailedTrials(B)+1
          IF( Errors(B) == 0 )THEN
            DRadiation(:,B) = D0(:,B)+Trial(1:M,B)
            DO J = 1, 3
              IRadiation(:,J,B) = I0(:,J,B)+Trial(J*M+1:(J+1)*M,B)
            END DO
            CALL ComputeCollisionRates_OrderV_Point(M,NF,A,EB,Kind,FM,W(:,B),U(:,B),ProbeGain(:,B),ProbeLoss(:,B), &
              ProbeChi(:,B),ProbeSigma(:,B),ProbeTarget(:,B),DRadiation(:,B),IRadiation(:,:,B), &
              NewC0(:,B),NewC1(:,:,B),NewSource(:,B),NewAbsolute(:,B),NewKappa(:,B),NewDensityLoss(:,B), &
              EnergyRatio_Option=EnergyRatio(:,B))
              CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Trial(:,B),Z0(:,:,B),DT, &
                NewC0(:,B),NewC1(:,:,B),NewSource(:,B),NewResidual(:,B),NewPhysical(:,B),NewHScale(:,B),NewGScale(:,B), &
                NewAbsolute(:,B),ScratchIncrement(:,:,B),TrialF(:,B),TrialG(:,B), &
                CollisionIncrement_Option=TrialIncrement(:,:,B),CollisionMoment_Option=TrialMoment(:,B), &
                Incoming_Option=Incoming(:,B), &
            LiftRows_Option=LiftRows(:,B))
            DO J=1,3
              Velocity(J,B)=Retain(V0(J,B))
            END DO
            DO J=1,NA
              IF( Active(J)>NF+1 ) Velocity(Active(J)-NF-1,B)=Velocity(Active(J)-NF-1,B)-Trial(N+J,B)
            END DO
            CALL ComputeIndependentRadiationResidual(M,NF,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
              TrialF(:,B),TrialG(:,B),TrialIndependent(:,B), &
            LiftRows_Option=LiftRows(:,B))
            TrialPhi=ComputeCoupledMerit(N,Weights(:,B),TrialIndependent(:,B),NewPhysical(:,B))
            Predicted=Fraction(B)*Slope(B)
            ! A bounded normalized composition equation can have an enormous
            ! slope near a vanishing incoming fraction. Do not demand a merit
            ! decrease larger than its nonnegative lower bound permits. This
            ! is ordinary Armijo again once the step is sufficiently small.
            Predicted=MAX(Predicted,-Merit(B))
            ! A normalized source can be flat before a sign crossing even
            ! for the exact root step. Permit ONE finite endpoint test when
            ! the local slope is nonnegative, requiring actual coupled merit
            ! reduction. A rejected non-descent endpoint is not backtracked.
            IF( .TRUE. .AND. Slope(B)>=0.0_DP ) Predicted=-Merit(B)
            IF( IEEE_IS_FINITE(TrialPhi) .AND. TrialPhi <= Merit(B)+1.0e-4_DP*Predicted )THEN
              Predicted=Merit(B)-ComputeCoupledMerit(N,Weights(:,B), &
                Independent(:,B)+Fraction(B)*Action(1:N,B),Physical(:,B)+Fraction(B)*Action(:,B))
              Ratio = 0.0_DP
              IF( Predicted > 0.0_DP ) Ratio = (Merit(B)-TrialPhi)/Predicted
              X(:,B) = Trial(:,B)
                ProcessIncrement(:,:,B)=TrialIncrement(:,:,B)
                ProcessMoment(:,B)=TrialMoment(:,B)
                ProcessNumber(:,B)=TrialNumber(:,B)
                ProcessEnergy(B)=TrialEnergy(B)
              Target(:,B) = ProbeTarget(:,B)
              Chi(:,B) = ProbeChi(:,B)
              Sigma(:,B) = ProbeSigma(:,B)
              Gain(:,B) = ProbeGain(:,B)
              Loss(:,B) = ProbeLoss(:,B)
              ComponentsValid(B)=.FALSE.
              KT(B) = ProbeKT(B)
              Mu(:,B) = ProbeMu(:,B)
              Susceptibility(B) = ProbeChiMu(B)
              ! Adopt the accepted result and its OWN scales, without a
              ! second EOS, opacity or rate evaluation at the next iterate.
              Residual(:,B) = NewResidual(:,B)
              Physical(:,B) = NewPhysical(:,B)
              Independent(:,B)=TrialIndependent(:,B)
              NativeF(:,B)=TrialF(:,B)
              TransferG(:,B)=TrialG(:,B)
              HScale(:,B) = NewHScale(:,B)
              GScale(:,B) = NewGScale(:,B)
              AbsoluteSource(:,B) = NewAbsolute(:,B)
              DensityLoss(:,B) = NewDensityLoss(:,B)
              KappaPair(:,B) = NewKappa(:,B)
              Kappa(:,B) = RetainedValue(NewKappa(:,B))
              Source(:,B) = NewSource(:,B)
              SourceValue(:,B) = RetainedValue(NewSource(:,B))
              C0(:,B) = NewC0(:,B)
              C1(:,:,B) = NewC1(:,:,B)
              NativeIncrement(:,:,B) = ScratchIncrement(:,:,B)
              CALL ComputeIndependentMetrics(NE,NS,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),X(:,B),DT,Physical(:,B),Independent(:,B), &
                Kappa(:,B),DensityLoss(:,B),HScale(:,B),GScale(:,B),InitialNorm(:,B),.FALSE., &
                RadiationError(B),ConservationError(B),Weights(:,B),Merit(B))
                CALL ComputeCoupledSourceMetrics(M,NF,Active,Fluid0(:,B),X(:,B),Physical(:,B),GScale(:,B), &
                  Tolerance,FluidSourceTolerance,Independent(:,B),FluidSourceError(:,B),Weights(:,B),Merit(B))
              AcceptedTrial(B) = .TRUE.
              Searching(B) = .FALSE.
              Accepted(B) = Accepted(B)+1
              Iterations(B) = Accepted(B)
              Age(B) = Age(B)+1
              Refresh(B) = Ratio < 0.1_DP .OR. Fraction(B) < 0.25_DP
              CYCLE
            END IF
          END IF
          Velocity(:,B)=Retain(V0(:,B))
          DO J=1,NA
            IF( Active(J)>NF+1 ) Velocity(Active(J)-NF-1,B)=Velocity(Active(J)-NF-1,B)-X(N+J,B)
          END DO
          Fraction(B) = 0.5_DP*Fraction(B)
          IF( .TRUE. .AND. Slope(B)>=0.0_DP ) Fraction(B)=0.0_DP
          IF( .NOT. Fresh(B) .AND. FailedTrials(B) >= 2 ) Fraction(B) = 0.0_DP
        END DO
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Searching,AcceptedTrial,Fresh,Refresh,FailedTrials)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Searching,AcceptedTrial,Fresh,Refresh,FailedTrials)
#endif
    END SUBROUTINE SearchDirection

    SUBROUTINE ProjectedSearch
      ! All active spatial points execute together. Only the small dual
      ! system is point-local; there is no serial loop over point solves.
      IF( .NOT. ANY(NeedsConstraint) ) RETURN
      IF( .NOT. ALLOCATED(DirectionScale) )THEN
        ALLOCATE(DirectionScale(N+NA,NB),Gradient(N+NA,NB),ConstraintJacobian(K,N+NA,NB))
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET ENTER DATA MAP(alloc:DirectionScale,Gradient,ConstraintJacobian)
#elif defined(THORNADO_OACC)
        !$ACC ENTER DATA CREATE(DirectionScale,Gradient,ConstraintJacobian)
#endif
      END IF
      ! Prepare only the constraint Jacobian and projection scales here.
      ! The backup gradient is needed only if projected Newton fails.
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !$OMP PRIVATE(J,Q,ColumnWork,VCWork)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE(J,Q,ColumnWork,VCWork)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(J,Q,ColumnWork,VCWork)
#endif
      DO B=1,NB
        IF( .NOT. NeedsConstraint(B) ) CYCLE
        SearchCounts(2,B) = SearchCounts(2,B)+1
        DO Q=1,M
          DirectionScale(Q,B)=MAX(RetainedValue(D0(Q,B)+X(Q,B)),D0(Q,B), &
            ABS(Physical(Q,B))/(1.0_DP+DT*DensityLoss(Q,B)))
          DO J=1,3
            DirectionScale(J*M+Q,B)=DirectionScale(Q,B)
          END DO
        END DO
        DO J=1,N+NA
          CALL ComputeCollisionLinearColumn_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
            Z0(:,:,B),DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
            ThermC(:,:,:,B),ThermS(:,:,B),J,.FALSE.,ColumnWork, &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
          ConstraintJacobian(:,J,B)=ColumnWork(N+1:N+K)
          IF( J<=N ) CYCLE
          VCWork=0.0_DP
          IF( Active(J-N)>NF+1 ) VCWork(Active(J-N)-NF-1)=-1.0_DP
          CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
            TransferG(:,B),VCWork,ColumnWork, &
            LiftRows_Option=LiftRows(:,B))
          DirectionScale(J,B)=1.0_DP/MAX(MAXVAL(ABS(Weights(:,B)*ColumnWork)),1.0e-100_DP)
        END DO
      END DO
      DO Choice=1,SIZE(Multipliers)
        ! SearchDirection synchronizes these masks after each candidate.
        ! Avoid both the gradient kernel and empty search launches once
        ! every remaining point has accepted or requires a fresh Jacobian.
        IF( .NOT. ANY(NeedsConstraint .AND. .NOT. AcceptedTrial &
                      .AND. (Fresh .OR. FailedTrials<2)) ) EXIT
        IF( Choice==2 ) CALL BuildProjectedGradient
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(J,VCWork)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(J,VCWork)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(J,VCWork)
#endif
        DO B=1,NB
          IF( .NOT. NeedsConstraint(B) .OR. AcceptedTrial(B) ) CYCLE
          IF( .NOT. Fresh(B) .AND. FailedTrials(B)>=2 ) CYCLE
          DRadiation(:,B)=D0(:,B)+X(1:M,B)
          DO J=1,3
            IRadiation(:,J,B)=I0(:,J,B)+X(J*M+1:(J+1)*M,B)
          END DO
          Step(:,B)=0.0_DP
          IF( Choice==1 )THEN
            Step(:,B)=NewtonStep(:,B)
          ELSE
            Step(1:N+NA,B)=Multipliers(Choice)*Gradient(:,B)
          END IF
          RetainedStep(:,B)=0.0_DP
          SearchCounts(5,B) = SearchCounts(5,B)+1
          CALL ProjectConservativeDirection(M,K,NA,DRadiation(:,B),IRadiation(:,:,B),DirectionScale(:,B), &
            ConstraintJacobian(:,:,B),Physical(N+1:N+K,B),Step(1:N+NA,B),RetainedStep(1:N+NA,B),Info(B))
          Fraction(B)=0.0_DP
          IF( Info(B)/=0 ) CYCLE
          Step(:,B)=RetainedValue(RetainedStep(:,B))
          CALL ComputeCollisionLinearAction_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
            Z0(:,:,B),DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
            ThermC(:,:,:,B),ThermS(:,:,B),Step(:,B),.FALSE.,Action(:,B), &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
          VCWork=0.0_DP
          DO J=1,NA
            IF( Active(J)>NF+1 ) VCWork(Active(J)-NF-1)=-Step(N+J,B)
          END DO
          CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
            TransferG(:,B),VCWork,Action(:,B), &
            LiftRows_Option=LiftRows(:,B))
          Slope(B)=ComputeCoupledSlope(N,Weights(:,B),Independent(:,B),Physical(:,B),Action(:,B))
        END DO
        CALL SearchDirection(18,.TRUE.)
      END DO
      DO B=1,NB
        IF( Searching(B) .AND. Fresh(B) ) Failure(B)=7
      END DO
    END SUBROUTINE ProjectedSearch

    SUBROUTINE BuildProjectedGradient
      ! Called once, only for points still searching after projected Newton.
      ! Reuse the projection scales at the unchanged iterate; retain the
      ! original gradient and step-length formulas and candidate ordering.
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !$OMP PRIVATE(J,ColumnWork,VCWork,ObjectiveScale,GradientNorm,ActionNorm,GradientLength)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE(J,ColumnWork,VCWork,ObjectiveScale,GradientNorm,ActionNorm,GradientLength)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(J,ColumnWork,VCWork,ObjectiveScale,GradientNorm,ActionNorm,GradientLength)
#endif
      DO B=1,NB
        IF( .NOT. NeedsConstraint(B) .OR. AcceptedTrial(B) ) CYCLE
        IF( .NOT. Fresh(B) .AND. FailedTrials(B)>=2 ) CYCLE
        SearchCounts(6,B)=SearchCounts(6,B)+1
        ObjectiveScale=MAX(SQRT(2.0_DP*Merit(B)),1.0e-100_DP)
        DO J=1,N+NA
          CALL ComputeCollisionLinearColumn_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
            Z0(:,:,B),DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
            ThermC(:,:,:,B),ThermS(:,:,B),J,.FALSE.,ColumnWork, &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
          VCWork=0.0_DP
          IF( J>N )THEN
            IF( Active(J-N)>NF+1 ) VCWork(Active(J-N)-NF-1)=-1.0_DP
          END IF
          CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
            TransferG(:,B),VCWork,ColumnWork, &
            LiftRows_Option=LiftRows(:,B))
          Gradient(J,B)=ComputeCoupledSlope(N,Weights(:,B),Independent(:,B),Physical(:,B),ColumnWork) &
            /ObjectiveScale*(DirectionScale(J,B)/ObjectiveScale)
        END DO
        Step(:,B)=0.0_DP
        Step(1:N+NA,B)=DirectionScale(:,B)*Gradient(:,B)
        CALL ComputeCollisionLinearAction_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
          Z0(:,:,B),DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
          ThermC(:,:,:,B),ThermS(:,:,B),Step(:,B),.FALSE.,ColumnWork, &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows)
        VCWork=0.0_DP
        DO J=1,NA
          IF( Active(J)>NF+1 ) VCWork(Active(J)-NF-1)=-Step(N+J,B)
        END DO
        CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          TransferG(:,B),VCWork,ColumnWork, &
            LiftRows_Option=LiftRows(:,B))
        GradientNorm=SUM(Gradient(:,B)**2)
        ActionNorm=2.0_DP*ComputeCoupledMerit(N,Weights(:,B),ColumnWork(1:N),ColumnWork)/ObjectiveScale**2
        GradientLength=GradientNorm/MAX(ActionNorm,TINY(1.0_DP))
        Gradient(:,B)=-GradientLength*DirectionScale(:,B)*Gradient(:,B)
      END DO
    END SUBROUTINE BuildProjectedGradient

    FUNCTION DenseCaptureLift(Point) RESULT(Dense)
      INTEGER, INTENT(in) :: Point
      REAL(DP) :: Dense(N,K)
      Dense=0.0_DP
      Dense(LiftRows(:,Point),:)=Z0(:,:,Point)
    END FUNCTION DenseCaptureLift

    FUNCTION DenseCaptureSpectral(Packed) RESULT(Dense)
      REAL(DP), INTENT(in) :: Packed(:,:)
      REAL(DP) :: Dense(M,M)
      INTEGER :: I,J
      Dense=0.0_DP
      DO J=1,M
      DO I=1,M
        IF( SpectralRows(2,I)/=SpectralRows(2,J) ) CYCLE
        Dense(I,J)=Packed(SpectralRows(3,J)+SpectralRows(1,I),1)
      END DO
      END DO
    END FUNCTION DenseCaptureSpectral

    FUNCTION DenseCaptureAngular(Point) RESULT(Dense)
      INTEGER, INTENT(in) :: Point
      REAL(DP) :: Dense(K,3*M)
      INTEGER :: J,Q
      Dense=0.0_DP
      DO J=1,3
      DO Q=1,M
        Dense(NF+1+J,(J-1)*M+Q)=SI(1,Q,Point)
      END DO
      END DO
    END FUNCTION DenseCaptureAngular

    SUBROUTINE WriteFailedLinearization
      CHARACTER(1024) :: TraceName,SnapshotName
      CHARACTER(24) :: PointName
      INTEGER :: Point,Unit,Status
      IF( .NOT. PRESENT(TraceUnit_Option) ) RETURN
      IF( TraceUnit_Option==0 .OR. TraceUnit_Option==6 ) RETURN
      IF( .NOT. ANY(Failure==7) ) RETURN
      INQUIRE(UNIT=TraceUnit_Option,NAME=TraceName,IOSTAT=Status)
      IF( Status/=0 ) RETURN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(X,Target,Chi,Sigma,Gain,Loss,Source,Physical,Residual, &
      !$OMP   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,NewtonStep,DChi,DSigma,DGain,DLoss,DTarget,Kappa, &
      !$OMP   ProcessIncrement,ProcessMoment,KT,Susceptibility)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(X,Target,Chi,Sigma,Gain,Loss,Source,Physical,Residual, &
      !$ACC   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,NewtonStep,DChi,DSigma,DGain,DLoss,DTarget,Kappa, &
      !$ACC   ProcessIncrement,ProcessMoment,KT,Susceptibility)
#endif
      DO Point=1,NB
        IF( Failure(Point)/=7 ) CYCLE
        WRITE(PointName,'(I0)') Indices(Point)
        SnapshotName=TRIM(TraceName)//'.newton_point'//TRIM(PointName)//'.bin'
        OPEN(NEWUNIT=Unit,FILE=TRIM(SnapshotName),ACCESS='STREAM',FORM='UNFORMATTED', &
          STATUS='REPLACE',IOSTAT=Status)
        IF( Status/=0 ) CYCLE
        WRITE(Unit) M,NF,NA,NV,Indices(Point)
        WRITE(Unit) DT,Cap,Tolerance,FluidSourceTolerance,KT(Point),Susceptibility(Point)
        WRITE(Unit) W(:,Point),U(:,Point),D0(:,Point),I0(:,:,Point),V0(:,Point),Fluid0(:,Point)
        WRITE(Unit) Active,A,EB,Kind,FM
        WRITE(Unit) X(:,Point)%High,X(:,Point)%Low,Target(:,Point)%High,Target(:,Point)%Low
        WRITE(Unit) Chi(:,Point),Sigma(:,Point),Gain(:,Point),Loss(:,Point),Kappa(:,Point)
        WRITE(Unit) Source(:,Point)%High,Source(:,Point)%Low,Residual(:,Point),Physical(:,Point)
        WRITE(Unit) DenseCaptureSpectral(CN(:,:,Point)),DenseCaptureSpectral(KD(:,:,Point)), &
          SD(:,:,Point),DenseCaptureAngular(Point),ThermC(:,:,:,Point),ThermS(:,:,Point)
        WRITE(Unit) DenseCaptureLift(Point),NewtonStep(:,Point),DChi(:,:,Point),DSigma(:,:,Point), &
          DGain(:,:,Point),DLoss(:,:,Point),DTarget(:,:,Point)%High,DTarget(:,:,Point)%Low
        WRITE(Unit) ProcessIncrement(:,:,Point)%High,ProcessIncrement(:,:,Point)%Low, &
          ProcessMoment(:,Point)%High,ProcessMoment(:,Point)%Low
        CLOSE(Unit)
      END DO
    END SUBROUTINE WriteFailedLinearization

    SUBROUTINE WriteLinearTrace(Label)
      CHARACTER(*), INTENT(in) :: Label
      REAL(DP) :: RadiationDefect,SourceDefect,RelativeDefect
      INTEGER :: Point
      IF( .NOT. PRESENT(TraceUnit_Option) ) RETURN
      IF( TraceUnit_Option==0 ) RETURN
      ! No transfers or extra reductions when tracing is disabled. These are
      ! linear-model diagnostics, not another acceptance or stopping threshold.
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Independent,Physical,Action,Weights,Merit,Slope,Running,Searching,Accepted)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Independent,Physical,Action,Weights,Merit,Slope,Running,Searching,Accepted)
#endif
      DO Point=1,NB
        IF( .NOT. Running(Point) ) CYCLE
        IF( Label/='newton' .AND. .NOT. Searching(Point) ) CYCLE
        RadiationDefect=SQRT(SUM((Weights(1:N,Point)*(Independent(:,Point)+Action(1:N,Point)))**2))
        SourceDefect=SQRT(SUM((Weights(N+K+1:,Point) &
          *(Physical(N+K+1:,Point)+Action(N+K+1:,Point)))**2))
        RelativeDefect=SQRT(RadiationDefect**2+SourceDefect**2)/MAX(SQRT(2.0_DP*Merit(Point)),SQRT(TINY(1.0_DP)))
        WRITE(TraceUnit_Option,'(A,2I8,1X,A,4ES23.14E3)') 'cpn_linear ',Indices(Point),Accepted(Point),Label, &
          RadiationDefect,SourceDefect,RelativeDefect,Slope(Point)
      END DO
    END SUBROUTINE WriteLinearTrace

  END SUBROUTINE SolveNeutrinoMatterCoupling_Newton_OrderV_Batch

END MODULE TwoMoment_NewtonSolverModule_OrderV

