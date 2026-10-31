MODULE TwoMoment_NewtonSolverModule_OrderV

  USE KindModule, ONLY: DP
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE LinearAlgebraModule, ONLY: LinearFactorBatched, LinearApplyBatched
  USE TwoMoment_RetainedArithmeticModule, ONLY: RetainedReal
  USE TwoMoment_CollisionOperatorModule_OrderV
  USE TwoMoment_RadiationStateRecoveryModule_OrderV
  USE TwoMoment_TimersModule, ONLY: TimersStart, TimersStop, &
    Timer_Newton_Residual,Timer_Newton_Jacobian,Timer_Newton_LinearSolve,Timer_Newton_Search

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_FMA
  IMPLICIT NONE
  PRIVATE
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE AssignReal
  END INTERFACE
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE AddPair, AddReal, RealAdd
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE SubtractPair, SubtractReal, RealSubtract, NegatePair
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE MultiplyPair, MultiplyReal, RealMultiply
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE DividePair, DivideReal, RealDivide
  END INTERFACE


  PUBLIC :: SolveNeutrinoMatterCoupling_Newton_OrderV_Batch
  PUBLIC :: FinalizeNewtonSolverWorkspace

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

  TYPE :: NewtonWorkspaceType
    INTEGER :: Key(9)=0
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessIncrement(:,:,:)=>NULL(),ProcessMoment(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: TrialIncrement(:,:,:)=>NULL(),TrialMoment(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DirectionIncrement(:,:,:)=>NULL(),DirectionMoment(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessNumber(:,:)=>NULL(),TrialNumber(:,:)=>NULL(),DirectionNumber(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessEnergy(:)=>NULL(),TrialEnergy(:)=>NULL(),DirectionEnergy(:)=>NULL()
    LOGICAL, POINTER, CONTIGUOUS :: SharedEndpointValid(:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: SharedProjectionInfo(:)=>NULL()
    TYPE(IncomingRadiationCache), POINTER, CONTIGUOUS :: Incoming(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: IncomingEnergy(:)=>NULL(),RecoveryResponse(:,:,:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: EnergyRatio(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: RealEnergyRatio(:,:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: Components(:,:)=>NULL()
    LOGICAL, POINTER, CONTIGUOUS :: ComponentsValid(:)=>NULL(),ContinueSearch(:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: BlockInfo(:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: SpectralLU(:)=>NULL(),BlockRHS(:)=>NULL(),BlockCorrection(:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NativeF(:,:)=>NULL(),TransferG(:,:)=>NULL(),TrialF(:,:)=>NULL(),TrialG(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: RetainedStep(:,:)=>NULL(),Velocity(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: Independent(:,:)=>NULL(),TrialIndependent(:,:)=>NULL(),NewtonStep(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: DirectionScale(:,:)=>NULL(),Gradient(:,:)=>NULL(),ConstraintJacobian(:,:,:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: LiftRows(:,:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: Info(:)=>NULL(),Errors(:)=>NULL(),Age(:)=>NULL(),Accepted(:)=>NULL(),FailedTrials(:)=>NULL()
    INTEGER, POINTER, CONTIGUOUS :: SpectralPivot(:)=>NULL(),BorderPivot(:)=>NULL(),SearchCounts(:,:)=>NULL()
    LOGICAL, POINTER, CONTIGUOUS :: Running(:)=>NULL(),Requested(:)=>NULL(),Refresh(:)=>NULL(),ThermalNeeded(:)=>NULL()
    LOGICAL, POINTER, CONTIGUOUS :: Searching(:)=>NULL(),AcceptedTrial(:)=>NULL(),Fresh(:)=>NULL(),NeedsConstraint(:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Trial(:,:)=>NULL(),Transfer(:,:)=>NULL(),ProbeTransfer(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Target(:,:)=>NULL(),ProbeTarget(:,:)=>NULL(),DTarget(:,:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: C0(:,:)=>NULL(),C1(:,:,:)=>NULL(),Source(:,:)=>NULL(),KappaPair(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NewC0(:,:)=>NULL(),NewC1(:,:,:)=>NULL(),NewSource(:,:)=>NULL(),NewKappa(:,:)=>NULL()
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DRadiation(:,:)=>NULL(),IRadiation(:,:,:)=>NULL(),ScratchIncrement(:,:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: KT(:)=>NULL(),Mu(:,:)=>NULL(),Susceptibility(:)=>NULL(),ProbeKT(:)=>NULL(),ProbeMu(:,:)=>NULL(),ProbeChiMu(:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: Chi(:,:)=>NULL(),Sigma(:,:)=>NULL(),Gain(:,:)=>NULL(),Loss(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: ProbeChi(:,:)=>NULL(),ProbeSigma(:,:)=>NULL(),ProbeGain(:,:)=>NULL(),ProbeLoss(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: DChi(:,:,:)=>NULL(),DSigma(:,:,:)=>NULL(),DGain(:,:,:)=>NULL(),DLoss(:,:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: Residual(:,:)=>NULL(),Physical(:,:)=>NULL(),NewResidual(:,:)=>NULL(),NewPhysical(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: HScale(:,:)=>NULL(),GScale(:,:)=>NULL(),NewHScale(:,:)=>NULL(),NewGScale(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: AbsoluteSource(:,:)=>NULL(),NewAbsolute(:,:)=>NULL(),DensityLoss(:,:)=>NULL(),NewDensityLoss(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: Kappa(:,:)=>NULL(),CN(:,:,:)=>NULL(),KD(:,:,:)=>NULL(),SD(:,:,:)=>NULL(),SI(:,:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: ThermC(:,:,:,:)=>NULL(),ThermS(:,:,:)=>NULL(),Z0(:,:,:)=>NULL(),BInverse(:,:,:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: SourceValue(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: Matrix(:,:,:)=>NULL(),RHS(:,:)=>NULL(),OriginalRHS(:,:)=>NULL(),Correction(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: RowScale(:,:)=>NULL(),ColumnScale(:,:)=>NULL(),Step(:,:)=>NULL(),Action(:,:)=>NULL(),Weights(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: InitialNorm(:,:)=>NULL(),Merit(:)=>NULL(),Slope(:)=>NULL(),Fraction(:)=>NULL(),FDStep(:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: FluidSourceError(:,:)=>NULL()
    REAL(DP), POINTER, CONTIGUOUS :: SpectralRHS(:,:,:)=>NULL(),BorderLU(:,:,:)=>NULL(),BorderRHS(:,:)=>NULL()
  END TYPE NewtonWorkspaceType
  TYPE(NewtonWorkspaceType), SAVE :: Workspace

CONTAINS


  ELEMENTAL SUBROUTINE AssignReal( A, B )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(out) :: A
    REAL(DP), INTENT(in) :: B
    A % High = B
    A % Low = 0.0_DP
  END SUBROUTINE AssignReal

  ELEMENTAL FUNCTION Retain( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    REAL(DP), INTENT(in), OPTIONAL :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: X, Y, Z

    IF( .NOT. PRESENT(B) )THEN
      C % High = A+0.0_DP
      C % Low = 0.0_DP
      RETURN
    END IF
    Y = B
    X = A + Y
    Z = X - A
    C % High = X
    C % Low = ( A - ( X - Z ) ) + ( Y - Z )
  END FUNCTION Retain

  ELEMENTAL FUNCTION RetainedValue( A ) RESULT(B)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP) :: B
    B = A % High + A % Low
  END FUNCTION RetainedValue

  ELEMENTAL FUNCTION AddPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A % High + B % High
    V = S - A % High
    E = ( A % High - ( S - V ) ) + ( B % High - V )
    C = Retain( S, ( E + A % Low ) + B % Low )
  END FUNCTION AddPair

  ELEMENTAL FUNCTION AddReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A % High+B
    V = S-A % High
    E = (A % High-(S-V))+(B-V)
    C = Retain(S,E+A % Low)
  END FUNCTION AddReal

  ELEMENTAL FUNCTION RealAdd( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A+B % High
    V = S-A
    E = (A-(S-V))+(B % High-V)
    C = Retain(S,E+B % Low)
  END FUNCTION RealAdd

  ELEMENTAL FUNCTION NegatePair( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    C % High = -A % High
    C % Low = -A % Low
  END FUNCTION NegatePair

  ELEMENTAL FUNCTION SubtractPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A % High-B % High
    V = S-A % High
    E = (A % High-(S-V))-(B % High+V)
    C = Retain(S,(E+A % Low)-B % Low)
  END FUNCTION SubtractPair

  ELEMENTAL FUNCTION SubtractReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A % High-B
    V = S-A % High
    E = (A % High-(S-V))-(B+V)
    C = Retain(S,E+A % Low)
  END FUNCTION SubtractReal

  ELEMENTAL FUNCTION RealSubtract( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A-B % High
    V = S-A
    E = (A-(S-V))-(B % High+V)
    C = Retain(S,E-B % Low)
  END FUNCTION RealSubtract

  ELEMENTAL FUNCTION MultiplyPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: P,E
    ! Match the checked prototype on CPUs. Explicit association is required;
    ! compile this arithmetic without reassociation/unsafe fast-math.
    P = A % High * B % High
    E = IEEE_FMA(A % High,B % High,-P)
    C = Retain( P, E + A % High * B % Low &
                      + A % Low * B % High + A % Low * B % Low )
  END FUNCTION MultiplyPair

  ELEMENTAL FUNCTION MultiplyReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: P,E
    P = A % High*B
    E = IEEE_FMA(A % High,B,-P)
    C = Retain(P,E+A % Low*B)
  END FUNCTION MultiplyReal

  ELEMENTAL FUNCTION RealMultiply( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: P,E
    P = A*B % High
    E = IEEE_FMA(A,B % High,-P)
    C = Retain(P,E+A*B % Low)
  END FUNCTION RealMultiply

  ELEMENTAL FUNCTION DividePair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C, R
    REAL(DP) :: Q
    Q = A % High / B % High
    R = A - B * Q
    C = Retain( Q, RetainedValue(R) / B % High )
  END FUNCTION DividePair

  ELEMENTAL FUNCTION DivideReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C, R, P
    REAL(DP) :: Q
    Q = A % High/B
    P = RetainedProduct(B,Q)
    R = A-P
    C = Retain(Q,RetainedValue(R)/B)
  END FUNCTION DivideReal

  ELEMENTAL FUNCTION RealDivide( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = Retain(A) / B
  END FUNCTION RealDivide

  ELEMENTAL FUNCTION RetainedAbs( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    C = A
    IF( A % High < 0.0_DP ) C = -A
  END FUNCTION RetainedAbs

  ELEMENTAL FUNCTION RetainedSqrt( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    INTEGER :: I
    C = SQRT(MAX(0.0_DP,A % High))
    IF( C % High == 0.0_DP ) RETURN
    DO I = 1, 2
      C = C + ( A - C * C ) / ( 2.0_DP * C )
    END DO
  END FUNCTION RetainedSqrt

  ELEMENTAL FUNCTION RetainedExp( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C, R
    TYPE(RetainedReal), PARAMETER :: LnStep = RetainedReal(4.33216987849965804E-02_DP,1.44940425865393722E-18_DP)
    TYPE(RetainedReal), PARAMETER :: Coefficient(0:13) = [ &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(5.00000000000000000E-01_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.66666666666666657E-01_DP,9.25185853854297066E-18_DP), &
      RetainedReal(4.16666666666666644E-02_DP,2.31296463463574266E-18_DP), &
      RetainedReal(8.33333333333333322E-03_DP,1.15648231731787138E-19_DP), &
      RetainedReal(1.38888888888888894E-03_DP,-5.30054395437357706E-20_DP), &
      RetainedReal(1.98412698412698413E-04_DP,1.72095582934207053E-22_DP), &
      RetainedReal(2.48015873015873016E-05_DP,2.15119478667758816E-23_DP), &
      RetainedReal(2.75573192239858925E-06_DP,-1.85839327404647208E-22_DP), &
      RetainedReal(2.75573192239858883E-07_DP,2.37677146222502973E-23_DP), &
      RetainedReal(2.50521083854417202E-08_DP,-1.44881407093591197E-24_DP), &
      RetainedReal(2.08767569878681002E-09_DP,-1.20734505911325997E-25_DP), &
      RetainedReal(1.60590438368216133E-10_DP,1.25852945887520981E-26_DP) ]
    TYPE(RetainedReal), PARAMETER :: FractionalPower(0:15) = [ &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.04427378242741375E+00_DP,8.55188970553796489E-17_DP), &
      RetainedReal(1.09050773266525769E+00_DP,-3.04678207981247115E-17_DP), &
      RetainedReal(1.13878863475669156E+00_DP,8.91281267602540778E-17_DP), &
      RetainedReal(1.18920711500272103E+00_DP,3.98201523146564611E-17_DP), &
      RetainedReal(1.24185781207348400E+00_DP,4.65802759183693679E-17_DP), &
      RetainedReal(1.29683955465100964E+00_DP,2.53825027948883150E-17_DP), &
      RetainedReal(1.35425554693689265E+00_DP,7.70094837980298946E-17_DP), &
      RetainedReal(1.41421356237309515E+00_DP,-9.66729331345291345E-17_DP), &
      RetainedReal(1.47682614593949935E+00_DP,-3.48399455689279580E-17_DP), &
      RetainedReal(1.54221082540794074E+00_DP,7.94983480969762086E-17_DP), &
      RetainedReal(1.61049033194925428E+00_DP,2.47071925697978879E-17_DP), &
      RetainedReal(1.68179283050742900E+00_DP,8.19901002058149652E-17_DP), &
      RetainedReal(1.75625216037329945E+00_DP,2.96014069544887331E-17_DP), &
      RetainedReal(1.83400808640934243E+00_DP,3.28310722424562720E-17_DP), &
      RetainedReal(1.91520656139714740E+00_DP,-1.06199460561959626E-16_DP) ]
    INTEGER :: I, N, J, Power

    IF( A % High < -745.0_DP )THEN
      C = 0.0_DP
      RETURN
    END IF
    ! exp(A) = 2**Power * 2**(J/16) * exp(R), |R| <= log(2)/32.
    ! The degree-13 Taylor remainder is below 5.87e-35.
    N = NINT(A % High/LnStep % High)
    J = MODULO(N,16)
    Power = (N-J)/16
    R = A-LnStep*REAL(N,DP)
    C = Coefficient(13)
    DO I = 12, 0, -1
      C = Coefficient(I)+R*C
    END DO
    C = C*FractionalPower(J)
    C % High = SCALE(C % High,Power)
    C % Low = SCALE(C % Low,Power)
  END FUNCTION RetainedExp

  ELEMENTAL FUNCTION RetainedFermi( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C, Tail
    Tail = RetainedExp(-RetainedAbs(A))
    IF( A % High >= 0.0_DP )THEN
      C = Tail / ( 1.0_DP + Tail )
    ELSE
      C = 1.0_DP / ( 1.0_DP + Tail )
    END IF
  END FUNCTION RetainedFermi

  PURE SUBROUTINE RetainedAccumulate( Sum, Term )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(inout) :: Sum
    TYPE(RetainedReal), INTENT(in) :: Term
    Sum = Sum + Term
  END SUBROUTINE RetainedAccumulate

  ELEMENTAL FUNCTION RetainedProduct( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    C % High = A*B
    C % Low = IEEE_FMA(A,B,-C % High)
  END FUNCTION RetainedProduct


  SUBROUTINE FinalizeNewtonSolverWorkspace
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessIncrement(:,:,:),ProcessMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: TrialIncrement(:,:,:),TrialMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DirectionIncrement(:,:,:),DirectionMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessNumber(:,:),TrialNumber(:,:),DirectionNumber(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessEnergy(:),TrialEnergy(:),DirectionEnergy(:)
    LOGICAL, POINTER, CONTIGUOUS :: SharedEndpointValid(:)
    INTEGER, POINTER, CONTIGUOUS :: SharedProjectionInfo(:)
    TYPE(IncomingRadiationCache), POINTER, CONTIGUOUS :: Incoming(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: IncomingEnergy(:),RecoveryResponse(:,:,:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: EnergyRatio(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: RealEnergyRatio(:,:)
    INTEGER, POINTER, CONTIGUOUS :: Components(:,:)
    LOGICAL, POINTER, CONTIGUOUS :: ComponentsValid(:),ContinueSearch(:)
    INTEGER, POINTER, CONTIGUOUS :: BlockInfo(:)
    REAL(DP), POINTER, CONTIGUOUS :: SpectralLU(:),BlockRHS(:),BlockCorrection(:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NativeF(:,:),TransferG(:,:),TrialF(:,:),TrialG(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: RetainedStep(:,:),Velocity(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Independent(:,:),TrialIndependent(:,:),NewtonStep(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: DirectionScale(:,:),Gradient(:,:),ConstraintJacobian(:,:,:)
    INTEGER, POINTER, CONTIGUOUS :: LiftRows(:,:)
    INTEGER, POINTER, CONTIGUOUS :: Info(:),Errors(:),Age(:),Accepted(:),FailedTrials(:)
    INTEGER, POINTER, CONTIGUOUS :: SpectralPivot(:),BorderPivot(:),SearchCounts(:,:)
    LOGICAL, POINTER, CONTIGUOUS :: Running(:),Requested(:),Refresh(:),ThermalNeeded(:)
    LOGICAL, POINTER, CONTIGUOUS :: Searching(:),AcceptedTrial(:),Fresh(:),NeedsConstraint(:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Trial(:,:),Transfer(:,:),ProbeTransfer(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Target(:,:),ProbeTarget(:,:),DTarget(:,:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: C0(:,:),C1(:,:,:),Source(:,:),KappaPair(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NewC0(:,:),NewC1(:,:,:),NewSource(:,:),NewKappa(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DRadiation(:,:),IRadiation(:,:,:),ScratchIncrement(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: KT(:),Mu(:,:),Susceptibility(:),ProbeKT(:),ProbeMu(:,:),ProbeChiMu(:)
    REAL(DP), POINTER, CONTIGUOUS :: Chi(:,:),Sigma(:,:),Gain(:,:),Loss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: ProbeChi(:,:),ProbeSigma(:,:),ProbeGain(:,:),ProbeLoss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: DChi(:,:,:),DSigma(:,:,:),DGain(:,:,:),DLoss(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Residual(:,:),Physical(:,:),NewResidual(:,:),NewPhysical(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: HScale(:,:),GScale(:,:),NewHScale(:,:),NewGScale(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: AbsoluteSource(:,:),NewAbsolute(:,:),DensityLoss(:,:),NewDensityLoss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Kappa(:,:),CN(:,:,:),KD(:,:,:),SD(:,:,:),SI(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: ThermC(:,:,:,:),ThermS(:,:,:),Z0(:,:,:),BInverse(:,:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: SourceValue(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Matrix(:,:,:),RHS(:,:),OriginalRHS(:,:),Correction(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: RowScale(:,:),ColumnScale(:,:),Step(:,:),Action(:,:),Weights(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: InitialNorm(:,:),Merit(:),Slope(:),Fraction(:),FDStep(:)
    REAL(DP), POINTER, CONTIGUOUS :: FluidSourceError(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: SpectralRHS(:,:,:),BorderLU(:,:,:),BorderRHS(:,:)
    IF( .NOT. ASSOCIATED(Workspace%Matrix) ) RETURN
    ProcessIncrement=>Workspace%ProcessIncrement
    ProcessMoment=>Workspace%ProcessMoment
    TrialIncrement=>Workspace%TrialIncrement
    TrialMoment=>Workspace%TrialMoment
    DirectionIncrement=>Workspace%DirectionIncrement
    DirectionMoment=>Workspace%DirectionMoment
    ProcessNumber=>Workspace%ProcessNumber
    TrialNumber=>Workspace%TrialNumber
    DirectionNumber=>Workspace%DirectionNumber
    ProcessEnergy=>Workspace%ProcessEnergy
    TrialEnergy=>Workspace%TrialEnergy
    DirectionEnergy=>Workspace%DirectionEnergy
    SharedEndpointValid=>Workspace%SharedEndpointValid
    SharedProjectionInfo=>Workspace%SharedProjectionInfo
    Incoming=>Workspace%Incoming
    IncomingEnergy=>Workspace%IncomingEnergy
    RecoveryResponse=>Workspace%RecoveryResponse
    EnergyRatio=>Workspace%EnergyRatio
    RealEnergyRatio=>Workspace%RealEnergyRatio
    Components=>Workspace%Components
    ComponentsValid=>Workspace%ComponentsValid
    ContinueSearch=>Workspace%ContinueSearch
    BlockInfo=>Workspace%BlockInfo
    SpectralLU=>Workspace%SpectralLU
    BlockRHS=>Workspace%BlockRHS
    BlockCorrection=>Workspace%BlockCorrection
    NativeF=>Workspace%NativeF
    TransferG=>Workspace%TransferG
    TrialF=>Workspace%TrialF
    TrialG=>Workspace%TrialG
    RetainedStep=>Workspace%RetainedStep
    Velocity=>Workspace%Velocity
    Independent=>Workspace%Independent
    TrialIndependent=>Workspace%TrialIndependent
    NewtonStep=>Workspace%NewtonStep
    DirectionScale=>Workspace%DirectionScale
    Gradient=>Workspace%Gradient
    ConstraintJacobian=>Workspace%ConstraintJacobian
    LiftRows=>Workspace%LiftRows
    Info=>Workspace%Info
    Errors=>Workspace%Errors
    Age=>Workspace%Age
    Accepted=>Workspace%Accepted
    FailedTrials=>Workspace%FailedTrials
    SpectralPivot=>Workspace%SpectralPivot
    BorderPivot=>Workspace%BorderPivot
    SearchCounts=>Workspace%SearchCounts
    Running=>Workspace%Running
    Requested=>Workspace%Requested
    Refresh=>Workspace%Refresh
    ThermalNeeded=>Workspace%ThermalNeeded
    Searching=>Workspace%Searching
    AcceptedTrial=>Workspace%AcceptedTrial
    Fresh=>Workspace%Fresh
    NeedsConstraint=>Workspace%NeedsConstraint
    Trial=>Workspace%Trial
    Transfer=>Workspace%Transfer
    ProbeTransfer=>Workspace%ProbeTransfer
    Target=>Workspace%Target
    ProbeTarget=>Workspace%ProbeTarget
    DTarget=>Workspace%DTarget
    C0=>Workspace%C0
    C1=>Workspace%C1
    Source=>Workspace%Source
    KappaPair=>Workspace%KappaPair
    NewC0=>Workspace%NewC0
    NewC1=>Workspace%NewC1
    NewSource=>Workspace%NewSource
    NewKappa=>Workspace%NewKappa
    DRadiation=>Workspace%DRadiation
    IRadiation=>Workspace%IRadiation
    ScratchIncrement=>Workspace%ScratchIncrement
    KT=>Workspace%KT
    Mu=>Workspace%Mu
    Susceptibility=>Workspace%Susceptibility
    ProbeKT=>Workspace%ProbeKT
    ProbeMu=>Workspace%ProbeMu
    ProbeChiMu=>Workspace%ProbeChiMu
    Chi=>Workspace%Chi
    Sigma=>Workspace%Sigma
    Gain=>Workspace%Gain
    Loss=>Workspace%Loss
    ProbeChi=>Workspace%ProbeChi
    ProbeSigma=>Workspace%ProbeSigma
    ProbeGain=>Workspace%ProbeGain
    ProbeLoss=>Workspace%ProbeLoss
    DChi=>Workspace%DChi
    DSigma=>Workspace%DSigma
    DGain=>Workspace%DGain
    DLoss=>Workspace%DLoss
    Residual=>Workspace%Residual
    Physical=>Workspace%Physical
    NewResidual=>Workspace%NewResidual
    NewPhysical=>Workspace%NewPhysical
    HScale=>Workspace%HScale
    GScale=>Workspace%GScale
    NewHScale=>Workspace%NewHScale
    NewGScale=>Workspace%NewGScale
    AbsoluteSource=>Workspace%AbsoluteSource
    NewAbsolute=>Workspace%NewAbsolute
    DensityLoss=>Workspace%DensityLoss
    NewDensityLoss=>Workspace%NewDensityLoss
    Kappa=>Workspace%Kappa
    CN=>Workspace%CN
    KD=>Workspace%KD
    SD=>Workspace%SD
    SI=>Workspace%SI
    ThermC=>Workspace%ThermC
    ThermS=>Workspace%ThermS
    Z0=>Workspace%Z0
    BInverse=>Workspace%BInverse
    SourceValue=>Workspace%SourceValue
    Matrix=>Workspace%Matrix
    RHS=>Workspace%RHS
    OriginalRHS=>Workspace%OriginalRHS
    Correction=>Workspace%Correction
    RowScale=>Workspace%RowScale
    ColumnScale=>Workspace%ColumnScale
    Step=>Workspace%Step
    Action=>Workspace%Action
    Weights=>Workspace%Weights
    InitialNorm=>Workspace%InitialNorm
    Merit=>Workspace%Merit
    Slope=>Workspace%Slope
    Fraction=>Workspace%Fraction
    FDStep=>Workspace%FDStep
    FluidSourceError=>Workspace%FluidSourceError
    SpectralRHS=>Workspace%SpectralRHS
    BorderLU=>Workspace%BorderLU
    BorderRHS=>Workspace%BorderRHS
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP(release: ThermalNeeded,FluidSourceError,SearchCounts, &
    !$OMP   Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh,Searching,AcceptedTrial,NeedsConstraint, &
    !$OMP   Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget,C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa, &
    !$OMP   DRadiation,IRadiation,ScratchIncrement,KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss, &
    !$OMP   ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual, &
    !$OMP   NewPhysical,HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa,CN,KD,SD,SI, &
    !$OMP   ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction,RowScale,ColumnScale,Step,Action,Weights, &
    !$OMP   InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$OMP   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$OMP   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE(ThermalNeeded,FluidSourceError,SearchCounts, &
    !$ACC   Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh,Searching,AcceptedTrial,NeedsConstraint, &
    !$ACC   Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget,C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa, &
    !$ACC   DRadiation,IRadiation,ScratchIncrement,KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss, &
    !$ACC   ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual, &
    !$ACC   NewPhysical,HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa,CN,KD,SD,SI, &
    !$ACC   ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction,RowScale,ColumnScale,Step,Action,Weights, &
    !$ACC   InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$ACC   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$ACC   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep)
#endif

    IF( ASSOCIATED(DirectionScale) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(release:DirectionScale,Gradient,ConstraintJacobian)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE(DirectionScale,Gradient,ConstraintJacobian)
#endif
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
    !$OMP TARGET EXIT DATA MAP(release:SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE(EnergyRatio,RealEnergyRatio,Components,ComponentsValid,ContinueSearch)
    !$ACC EXIT DATA DELETE(SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#endif

    IF( ASSOCIATED(DirectionScale) )THEN
      DEALLOCATE(DirectionScale,Gradient,ConstraintJacobian)
    END IF
    DEALLOCATE(ProcessIncrement,ProcessMoment,TrialIncrement,TrialMoment,DirectionIncrement)
    DEALLOCATE(DirectionMoment,ProcessNumber,TrialNumber,DirectionNumber,ProcessEnergy)
    DEALLOCATE(TrialEnergy,DirectionEnergy,SharedEndpointValid,SharedProjectionInfo,Incoming)
    DEALLOCATE(IncomingEnergy,RecoveryResponse,EnergyRatio,RealEnergyRatio,Components)
    DEALLOCATE(ComponentsValid,ContinueSearch,BlockInfo,SpectralLU,BlockRHS)
    DEALLOCATE(BlockCorrection,NativeF,TransferG,TrialF,TrialG)
    DEALLOCATE(RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep)
    DEALLOCATE(LiftRows,Info,Errors,Age,Accepted)
    DEALLOCATE(FailedTrials,SpectralPivot,BorderPivot,SearchCounts,Running)
    DEALLOCATE(Requested,Refresh,ThermalNeeded,Searching,AcceptedTrial)
    DEALLOCATE(Fresh,NeedsConstraint,Trial,Transfer,ProbeTransfer)
    DEALLOCATE(Target,ProbeTarget,DTarget,C0,C1)
    DEALLOCATE(Source,KappaPair,NewC0,NewC1,NewSource)
    DEALLOCATE(NewKappa,DRadiation,IRadiation,ScratchIncrement,KT)
    DEALLOCATE(Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu)
    DEALLOCATE(Chi,Sigma,Gain,Loss,ProbeChi)
    DEALLOCATE(ProbeSigma,ProbeGain,ProbeLoss,DChi,DSigma)
    DEALLOCATE(DGain,DLoss,Residual,Physical,NewResidual)
    DEALLOCATE(NewPhysical,HScale,GScale,NewHScale,NewGScale)
    DEALLOCATE(AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa)
    DEALLOCATE(CN,KD,SD,SI,ThermC)
    DEALLOCATE(ThermS,Z0,BInverse,SourceValue,Matrix)
    DEALLOCATE(RHS,OriginalRHS,Correction,RowScale,ColumnScale)
    DEALLOCATE(Step,Action,Weights,InitialNorm,Merit)
    DEALLOCATE(Slope,Fraction,FDStep,FluidSourceError,SpectralRHS)
    DEALLOCATE(BorderLU,BorderRHS)
    NULLIFY(Workspace%ProcessIncrement,Workspace%ProcessMoment,Workspace%TrialIncrement,Workspace%TrialMoment)
    NULLIFY(Workspace%DirectionIncrement,Workspace%DirectionMoment,Workspace%ProcessNumber,Workspace%TrialNumber)
    NULLIFY(Workspace%DirectionNumber,Workspace%ProcessEnergy,Workspace%TrialEnergy,Workspace%DirectionEnergy)
    NULLIFY(Workspace%SharedEndpointValid,Workspace%SharedProjectionInfo,Workspace%Incoming,Workspace%IncomingEnergy)
    NULLIFY(Workspace%RecoveryResponse,Workspace%EnergyRatio,Workspace%RealEnergyRatio,Workspace%Components)
    NULLIFY(Workspace%ComponentsValid,Workspace%ContinueSearch,Workspace%BlockInfo,Workspace%SpectralLU)
    NULLIFY(Workspace%BlockRHS,Workspace%BlockCorrection,Workspace%NativeF,Workspace%TransferG)
    NULLIFY(Workspace%TrialF,Workspace%TrialG,Workspace%RetainedStep,Workspace%Velocity)
    NULLIFY(Workspace%Independent,Workspace%TrialIndependent,Workspace%NewtonStep,Workspace%DirectionScale)
    NULLIFY(Workspace%Gradient,Workspace%ConstraintJacobian,Workspace%LiftRows,Workspace%Info)
    NULLIFY(Workspace%Errors,Workspace%Age,Workspace%Accepted,Workspace%FailedTrials)
    NULLIFY(Workspace%SpectralPivot,Workspace%BorderPivot,Workspace%SearchCounts,Workspace%Running)
    NULLIFY(Workspace%Requested,Workspace%Refresh,Workspace%ThermalNeeded,Workspace%Searching)
    NULLIFY(Workspace%AcceptedTrial,Workspace%Fresh,Workspace%NeedsConstraint,Workspace%Trial)
    NULLIFY(Workspace%Transfer,Workspace%ProbeTransfer,Workspace%Target,Workspace%ProbeTarget)
    NULLIFY(Workspace%DTarget,Workspace%C0,Workspace%C1,Workspace%Source)
    NULLIFY(Workspace%KappaPair,Workspace%NewC0,Workspace%NewC1,Workspace%NewSource)
    NULLIFY(Workspace%NewKappa,Workspace%DRadiation,Workspace%IRadiation,Workspace%ScratchIncrement)
    NULLIFY(Workspace%KT,Workspace%Mu,Workspace%Susceptibility,Workspace%ProbeKT)
    NULLIFY(Workspace%ProbeMu,Workspace%ProbeChiMu,Workspace%Chi,Workspace%Sigma)
    NULLIFY(Workspace%Gain,Workspace%Loss,Workspace%ProbeChi,Workspace%ProbeSigma)
    NULLIFY(Workspace%ProbeGain,Workspace%ProbeLoss,Workspace%DChi,Workspace%DSigma)
    NULLIFY(Workspace%DGain,Workspace%DLoss,Workspace%Residual,Workspace%Physical)
    NULLIFY(Workspace%NewResidual,Workspace%NewPhysical,Workspace%HScale,Workspace%GScale)
    NULLIFY(Workspace%NewHScale,Workspace%NewGScale,Workspace%AbsoluteSource,Workspace%NewAbsolute)
    NULLIFY(Workspace%DensityLoss,Workspace%NewDensityLoss,Workspace%Kappa,Workspace%CN)
    NULLIFY(Workspace%KD,Workspace%SD,Workspace%SI,Workspace%ThermC)
    NULLIFY(Workspace%ThermS,Workspace%Z0,Workspace%BInverse,Workspace%SourceValue)
    NULLIFY(Workspace%Matrix,Workspace%RHS,Workspace%OriginalRHS,Workspace%Correction)
    NULLIFY(Workspace%RowScale,Workspace%ColumnScale,Workspace%Step,Workspace%Action)
    NULLIFY(Workspace%Weights,Workspace%InitialNorm,Workspace%Merit,Workspace%Slope)
    NULLIFY(Workspace%Fraction,Workspace%FDStep,Workspace%FluidSourceError,Workspace%SpectralRHS)
    NULLIFY(Workspace%BorderLU,Workspace%BorderRHS)
    Workspace%Key=0
  END SUBROUTINE FinalizeNewtonSolverWorkspace


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
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessIncrement(:,:,:),ProcessMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: TrialIncrement(:,:,:),TrialMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DirectionIncrement(:,:,:),DirectionMoment(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessNumber(:,:),TrialNumber(:,:),DirectionNumber(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: ProcessEnergy(:),TrialEnergy(:),DirectionEnergy(:)
    LOGICAL, POINTER, CONTIGUOUS :: SharedEndpointValid(:)
    INTEGER, POINTER, CONTIGUOUS :: SharedProjectionInfo(:)
    TYPE(IncomingRadiationCache), POINTER, CONTIGUOUS :: Incoming(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: IncomingEnergy(:),RecoveryResponse(:,:,:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: EnergyRatio(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: RealEnergyRatio(:,:)
    INTEGER, POINTER, CONTIGUOUS :: Components(:,:)
    LOGICAL, POINTER, CONTIGUOUS :: ComponentsValid(:),ContinueSearch(:)
    INTEGER, ALLOCATABLE, TARGET :: SpectralRows(:,:),BlockMembers(:),BlockFirst(:),BlockLength(:)
    INTEGER, POINTER, CONTIGUOUS :: BlockInfo(:)
    INTEGER, ALLOCATABLE :: GroupFirst(:),GroupCount(:),GroupSize(:),GroupLU(:),GroupRHS(:),GroupVector(:),GroupInfo(:)
    REAL(DP), POINTER, CONTIGUOUS :: SpectralLU(:),BlockRHS(:),BlockCorrection(:)
    INTEGER :: SpectralSize,NGroups,Group,Block,LUSize,RHSSize,VectorSize,InfoSize

    INTEGER :: M, N, NF, K, NA, NC, NX, NB, NV, B, J, Q, F, Loop, Cut, Pass, Ierr
    INTEGER :: Status, VendorStatus, Mode, Therm
    INTEGER :: Choice,WorkspaceKey(9)
    LOGICAL :: ReuseWorkspace,FullDirectionAction
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NativeF(:,:),TransferG(:,:),TrialF(:,:),TrialG(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: RetainedStep(:,:),Velocity(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Independent(:,:),TrialIndependent(:,:),NewtonStep(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: DirectionScale(:,:),Gradient(:,:),ConstraintJacobian(:,:,:)
    REAL(DP) :: ColumnWork(SIZE(X,1)), ScaleWork(SIZE(X,1)), VCWork(3)
    REAL(DP) :: ObjectiveScale, GradientLength, GradientNorm, ActionNorm
    REAL(DP), PARAMETER :: Multipliers(6)=[0.0_DP,1.0_DP,0.1_DP,0.01_DP,0.001_DP,0.0001_DP]
    INTEGER, ALLOCATABLE :: A(:), EB(:), Kind(:), Process(:), Species(:), QA(:), QB(:), FM(:,:)
    INTEGER, POINTER, CONTIGUOUS :: LiftRows(:,:)
    INTEGER, POINTER, CONTIGUOUS :: Info(:),Errors(:),Age(:),Accepted(:),FailedTrials(:)
    INTEGER, POINTER, CONTIGUOUS :: SpectralPivot(:),BorderPivot(:),SearchCounts(:,:)
    LOGICAL, POINTER, CONTIGUOUS :: Running(:),Requested(:),Refresh(:),ThermalNeeded(:)
    LOGICAL, POINTER, CONTIGUOUS :: Searching(:),AcceptedTrial(:),Fresh(:),NeedsConstraint(:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Trial(:,:),Transfer(:,:),ProbeTransfer(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: Target(:,:),ProbeTarget(:,:),DTarget(:,:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: C0(:,:),C1(:,:,:),Source(:,:),KappaPair(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: NewC0(:,:),NewC1(:,:,:),NewSource(:,:),NewKappa(:,:)
    TYPE(RetainedReal), POINTER, CONTIGUOUS :: DRadiation(:,:),IRadiation(:,:,:),ScratchIncrement(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: KT(:),Mu(:,:),Susceptibility(:),ProbeKT(:),ProbeMu(:,:),ProbeChiMu(:)
    REAL(DP), POINTER, CONTIGUOUS :: Chi(:,:),Sigma(:,:),Gain(:,:),Loss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: ProbeChi(:,:),ProbeSigma(:,:),ProbeGain(:,:),ProbeLoss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: DChi(:,:,:),DSigma(:,:,:),DGain(:,:,:),DLoss(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Residual(:,:),Physical(:,:),NewResidual(:,:),NewPhysical(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: HScale(:,:),GScale(:,:),NewHScale(:,:),NewGScale(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: AbsoluteSource(:,:),NewAbsolute(:,:),DensityLoss(:,:),NewDensityLoss(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Kappa(:,:),CN(:,:,:),KD(:,:,:),SD(:,:,:),SI(:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: ThermC(:,:,:,:),ThermS(:,:,:),Z0(:,:,:),BInverse(:,:,:,:)
    REAL(DP), POINTER, CONTIGUOUS :: SourceValue(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: Matrix(:,:,:),RHS(:,:),OriginalRHS(:,:),Correction(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: RowScale(:,:),ColumnScale(:,:),Step(:,:),Action(:,:),Weights(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: InitialNorm(:,:),Merit(:),Slope(:),Fraction(:),FDStep(:)
    REAL(DP), POINTER, CONTIGUOUS :: FluidSourceError(:,:)
    REAL(DP), POINTER, CONTIGUOUS :: SpectralRHS(:,:,:),BorderLU(:,:,:),BorderRHS(:,:)
    REAL(DP) :: IValues(SIZE(D0,1),3), TempSource(NS/2+4)
    REAL(DP) :: Predicted, Ratio, TrialPhi
    TYPE(RetainedReal) :: Term, MuonTransfer
    TYPE(CollisionGeometryCache) :: CurrentGeometry(SIZE(Indices)),TrialGeometry(SIZE(Indices))

    TYPE(RetainedReal) :: CurrentNumberSource(NE*NS,SIZE(Indices)),CurrentEnergySource(SIZE(Indices))

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
      RHSSize=RHSSize+GroupSize(NGroups)*(NA+1)*GroupCount(NGroups)*NB
      VectorSize=VectorSize+GroupSize(NGroups)*GroupCount(NGroups)*NB
      InfoSize=InfoSize+GroupCount(NGroups)*NB
      Block=Block+GroupCount(NGroups)
    END DO
    WorkspaceKey=[M,NF,NA,NB,NV,SpectralSize,LUSize,RHSSize,InfoSize]
    ReuseWorkspace=ASSOCIATED(Workspace%Matrix)
    IF( ReuseWorkspace ) ReuseWorkspace=ALL(Workspace%Key==WorkspaceKey)
    IF( .NOT. ReuseWorkspace ) CALL FinalizeNewtonSolverWorkspace
    ProcessIncrement=>Workspace%ProcessIncrement
    ProcessMoment=>Workspace%ProcessMoment
    TrialIncrement=>Workspace%TrialIncrement
    TrialMoment=>Workspace%TrialMoment
    DirectionIncrement=>Workspace%DirectionIncrement
    DirectionMoment=>Workspace%DirectionMoment
    ProcessNumber=>Workspace%ProcessNumber
    TrialNumber=>Workspace%TrialNumber
    DirectionNumber=>Workspace%DirectionNumber
    ProcessEnergy=>Workspace%ProcessEnergy
    TrialEnergy=>Workspace%TrialEnergy
    DirectionEnergy=>Workspace%DirectionEnergy
    SharedEndpointValid=>Workspace%SharedEndpointValid
    SharedProjectionInfo=>Workspace%SharedProjectionInfo
    Incoming=>Workspace%Incoming
    IncomingEnergy=>Workspace%IncomingEnergy
    RecoveryResponse=>Workspace%RecoveryResponse
    EnergyRatio=>Workspace%EnergyRatio
    RealEnergyRatio=>Workspace%RealEnergyRatio
    Components=>Workspace%Components
    ComponentsValid=>Workspace%ComponentsValid
    ContinueSearch=>Workspace%ContinueSearch
    BlockInfo=>Workspace%BlockInfo
    SpectralLU=>Workspace%SpectralLU
    BlockRHS=>Workspace%BlockRHS
    BlockCorrection=>Workspace%BlockCorrection
    NativeF=>Workspace%NativeF
    TransferG=>Workspace%TransferG
    TrialF=>Workspace%TrialF
    TrialG=>Workspace%TrialG
    RetainedStep=>Workspace%RetainedStep
    Velocity=>Workspace%Velocity
    Independent=>Workspace%Independent
    TrialIndependent=>Workspace%TrialIndependent
    NewtonStep=>Workspace%NewtonStep
    DirectionScale=>Workspace%DirectionScale
    Gradient=>Workspace%Gradient
    ConstraintJacobian=>Workspace%ConstraintJacobian
    LiftRows=>Workspace%LiftRows
    Info=>Workspace%Info
    Errors=>Workspace%Errors
    Age=>Workspace%Age
    Accepted=>Workspace%Accepted
    FailedTrials=>Workspace%FailedTrials
    SpectralPivot=>Workspace%SpectralPivot
    BorderPivot=>Workspace%BorderPivot
    SearchCounts=>Workspace%SearchCounts
    Running=>Workspace%Running
    Requested=>Workspace%Requested
    Refresh=>Workspace%Refresh
    ThermalNeeded=>Workspace%ThermalNeeded
    Searching=>Workspace%Searching
    AcceptedTrial=>Workspace%AcceptedTrial
    Fresh=>Workspace%Fresh
    NeedsConstraint=>Workspace%NeedsConstraint
    Trial=>Workspace%Trial
    Transfer=>Workspace%Transfer
    ProbeTransfer=>Workspace%ProbeTransfer
    Target=>Workspace%Target
    ProbeTarget=>Workspace%ProbeTarget
    DTarget=>Workspace%DTarget
    C0=>Workspace%C0
    C1=>Workspace%C1
    Source=>Workspace%Source
    KappaPair=>Workspace%KappaPair
    NewC0=>Workspace%NewC0
    NewC1=>Workspace%NewC1
    NewSource=>Workspace%NewSource
    NewKappa=>Workspace%NewKappa
    DRadiation=>Workspace%DRadiation
    IRadiation=>Workspace%IRadiation
    ScratchIncrement=>Workspace%ScratchIncrement
    KT=>Workspace%KT
    Mu=>Workspace%Mu
    Susceptibility=>Workspace%Susceptibility
    ProbeKT=>Workspace%ProbeKT
    ProbeMu=>Workspace%ProbeMu
    ProbeChiMu=>Workspace%ProbeChiMu
    Chi=>Workspace%Chi
    Sigma=>Workspace%Sigma
    Gain=>Workspace%Gain
    Loss=>Workspace%Loss
    ProbeChi=>Workspace%ProbeChi
    ProbeSigma=>Workspace%ProbeSigma
    ProbeGain=>Workspace%ProbeGain
    ProbeLoss=>Workspace%ProbeLoss
    DChi=>Workspace%DChi
    DSigma=>Workspace%DSigma
    DGain=>Workspace%DGain
    DLoss=>Workspace%DLoss
    Residual=>Workspace%Residual
    Physical=>Workspace%Physical
    NewResidual=>Workspace%NewResidual
    NewPhysical=>Workspace%NewPhysical
    HScale=>Workspace%HScale
    GScale=>Workspace%GScale
    NewHScale=>Workspace%NewHScale
    NewGScale=>Workspace%NewGScale
    AbsoluteSource=>Workspace%AbsoluteSource
    NewAbsolute=>Workspace%NewAbsolute
    DensityLoss=>Workspace%DensityLoss
    NewDensityLoss=>Workspace%NewDensityLoss
    Kappa=>Workspace%Kappa
    CN=>Workspace%CN
    KD=>Workspace%KD
    SD=>Workspace%SD
    SI=>Workspace%SI
    ThermC=>Workspace%ThermC
    ThermS=>Workspace%ThermS
    Z0=>Workspace%Z0
    BInverse=>Workspace%BInverse
    SourceValue=>Workspace%SourceValue
    Matrix=>Workspace%Matrix
    RHS=>Workspace%RHS
    OriginalRHS=>Workspace%OriginalRHS
    Correction=>Workspace%Correction
    RowScale=>Workspace%RowScale
    ColumnScale=>Workspace%ColumnScale
    Step=>Workspace%Step
    Action=>Workspace%Action
    Weights=>Workspace%Weights
    InitialNorm=>Workspace%InitialNorm
    Merit=>Workspace%Merit
    Slope=>Workspace%Slope
    Fraction=>Workspace%Fraction
    FDStep=>Workspace%FDStep
    FluidSourceError=>Workspace%FluidSourceError
    SpectralRHS=>Workspace%SpectralRHS
    BorderLU=>Workspace%BorderLU
    BorderRHS=>Workspace%BorderRHS
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to:SpectralRows,BlockMembers,BlockFirst,BlockLength)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(SpectralRows,BlockMembers,BlockFirst,BlockLength)
#endif
    IF( .NOT. ReuseWorkspace )THEN
    ALLOCATE(SpectralLU(LUSize),BlockRHS(RHSSize),BlockCorrection(VectorSize),BlockInfo(InfoSize))
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(alloc:SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE(SpectralLU,BlockRHS,BlockCorrection,BlockInfo)
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
    ALLOCATE(SpectralRHS(M,NA+1,NB),BorderLU(NA,NA,NB),BorderRHS(NA,NB))
    ALLOCATE(SpectralPivot(M*NB),BorderPivot(NA*NB))
    ALLOCATE(SearchCounts(6,NB))
    Workspace%ProcessIncrement=>ProcessIncrement
    Workspace%ProcessMoment=>ProcessMoment
    Workspace%TrialIncrement=>TrialIncrement
    Workspace%TrialMoment=>TrialMoment
    Workspace%DirectionIncrement=>DirectionIncrement
    Workspace%DirectionMoment=>DirectionMoment
    Workspace%ProcessNumber=>ProcessNumber
    Workspace%TrialNumber=>TrialNumber
    Workspace%DirectionNumber=>DirectionNumber
    Workspace%ProcessEnergy=>ProcessEnergy
    Workspace%TrialEnergy=>TrialEnergy
    Workspace%DirectionEnergy=>DirectionEnergy
    Workspace%SharedEndpointValid=>SharedEndpointValid
    Workspace%SharedProjectionInfo=>SharedProjectionInfo
    Workspace%Incoming=>Incoming
    Workspace%IncomingEnergy=>IncomingEnergy
    Workspace%RecoveryResponse=>RecoveryResponse
    Workspace%EnergyRatio=>EnergyRatio
    Workspace%RealEnergyRatio=>RealEnergyRatio
    Workspace%Components=>Components
    Workspace%ComponentsValid=>ComponentsValid
    Workspace%ContinueSearch=>ContinueSearch
    Workspace%BlockInfo=>BlockInfo
    Workspace%SpectralLU=>SpectralLU
    Workspace%BlockRHS=>BlockRHS
    Workspace%BlockCorrection=>BlockCorrection
    Workspace%NativeF=>NativeF
    Workspace%TransferG=>TransferG
    Workspace%TrialF=>TrialF
    Workspace%TrialG=>TrialG
    Workspace%RetainedStep=>RetainedStep
    Workspace%Velocity=>Velocity
    Workspace%Independent=>Independent
    Workspace%TrialIndependent=>TrialIndependent
    Workspace%NewtonStep=>NewtonStep
    Workspace%DirectionScale=>DirectionScale
    Workspace%Gradient=>Gradient
    Workspace%ConstraintJacobian=>ConstraintJacobian
    Workspace%LiftRows=>LiftRows
    Workspace%Info=>Info
    Workspace%Errors=>Errors
    Workspace%Age=>Age
    Workspace%Accepted=>Accepted
    Workspace%FailedTrials=>FailedTrials
    Workspace%SpectralPivot=>SpectralPivot
    Workspace%BorderPivot=>BorderPivot
    Workspace%SearchCounts=>SearchCounts
    Workspace%Running=>Running
    Workspace%Requested=>Requested
    Workspace%Refresh=>Refresh
    Workspace%ThermalNeeded=>ThermalNeeded
    Workspace%Searching=>Searching
    Workspace%AcceptedTrial=>AcceptedTrial
    Workspace%Fresh=>Fresh
    Workspace%NeedsConstraint=>NeedsConstraint
    Workspace%Trial=>Trial
    Workspace%Transfer=>Transfer
    Workspace%ProbeTransfer=>ProbeTransfer
    Workspace%Target=>Target
    Workspace%ProbeTarget=>ProbeTarget
    Workspace%DTarget=>DTarget
    Workspace%C0=>C0
    Workspace%C1=>C1
    Workspace%Source=>Source
    Workspace%KappaPair=>KappaPair
    Workspace%NewC0=>NewC0
    Workspace%NewC1=>NewC1
    Workspace%NewSource=>NewSource
    Workspace%NewKappa=>NewKappa
    Workspace%DRadiation=>DRadiation
    Workspace%IRadiation=>IRadiation
    Workspace%ScratchIncrement=>ScratchIncrement
    Workspace%KT=>KT
    Workspace%Mu=>Mu
    Workspace%Susceptibility=>Susceptibility
    Workspace%ProbeKT=>ProbeKT
    Workspace%ProbeMu=>ProbeMu
    Workspace%ProbeChiMu=>ProbeChiMu
    Workspace%Chi=>Chi
    Workspace%Sigma=>Sigma
    Workspace%Gain=>Gain
    Workspace%Loss=>Loss
    Workspace%ProbeChi=>ProbeChi
    Workspace%ProbeSigma=>ProbeSigma
    Workspace%ProbeGain=>ProbeGain
    Workspace%ProbeLoss=>ProbeLoss
    Workspace%DChi=>DChi
    Workspace%DSigma=>DSigma
    Workspace%DGain=>DGain
    Workspace%DLoss=>DLoss
    Workspace%Residual=>Residual
    Workspace%Physical=>Physical
    Workspace%NewResidual=>NewResidual
    Workspace%NewPhysical=>NewPhysical
    Workspace%HScale=>HScale
    Workspace%GScale=>GScale
    Workspace%NewHScale=>NewHScale
    Workspace%NewGScale=>NewGScale
    Workspace%AbsoluteSource=>AbsoluteSource
    Workspace%NewAbsolute=>NewAbsolute
    Workspace%DensityLoss=>DensityLoss
    Workspace%NewDensityLoss=>NewDensityLoss
    Workspace%Kappa=>Kappa
    Workspace%CN=>CN
    Workspace%KD=>KD
    Workspace%SD=>SD
    Workspace%SI=>SI
    Workspace%ThermC=>ThermC
    Workspace%ThermS=>ThermS
    Workspace%Z0=>Z0
    Workspace%BInverse=>BInverse
    Workspace%SourceValue=>SourceValue
    Workspace%Matrix=>Matrix
    Workspace%RHS=>RHS
    Workspace%OriginalRHS=>OriginalRHS
    Workspace%Correction=>Correction
    Workspace%RowScale=>RowScale
    Workspace%ColumnScale=>ColumnScale
    Workspace%Step=>Step
    Workspace%Action=>Action
    Workspace%Weights=>Weights
    Workspace%InitialNorm=>InitialNorm
    Workspace%Merit=>Merit
    Workspace%Slope=>Slope
    Workspace%Fraction=>Fraction
    Workspace%FDStep=>FDStep
    Workspace%FluidSourceError=>FluidSourceError
    Workspace%SpectralRHS=>SpectralRHS
    Workspace%BorderLU=>BorderLU
    Workspace%BorderRHS=>BorderRHS
    Workspace%Key=WorkspaceKey
    END IF
    FullDirectionAction=.FALSE.
    IF( PRESENT(TraceUnit_Option) ) FullDirectionAction=TraceUnit_Option/=0
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

    IF( .NOT. ReuseWorkspace )THEN
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to:ThermalNeeded) &
    !$OMP MAP(alloc: Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh, &
    !$OMP   Searching,AcceptedTrial,NeedsConstraint,Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget, &
    !$OMP   C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa,DRadiation,IRadiation,ScratchIncrement, &
    !$OMP   KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss, &
    !$OMP   DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual,NewPhysical, &
    !$OMP   HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa, &
    !$OMP   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction, &
    !$OMP   RowScale,ColumnScale,Step,Action,Weights,InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$OMP   FluidSourceError, &
    !$OMP   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$OMP   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep,SearchCounts)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(ThermalNeeded)
    !$ACC ENTER DATA CREATE(Info,Errors,Age,Accepted,FailedTrials,Running,Requested,Refresh,Fresh, &
    !$ACC   Searching,AcceptedTrial,NeedsConstraint,Trial,Transfer,ProbeTransfer,Target,ProbeTarget,DTarget, &
    !$ACC   C0,C1,Source,KappaPair,NewC0,NewC1,NewSource,NewKappa,DRadiation,IRadiation,ScratchIncrement, &
    !$ACC   KT,Mu,Susceptibility,ProbeKT,ProbeMu,ProbeChiMu,Chi,Sigma,Gain,Loss,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss, &
    !$ACC   DChi,DSigma,DGain,DLoss,Residual,Physical,NewResidual,NewPhysical, &
    !$ACC   HScale,GScale,NewHScale,NewGScale,AbsoluteSource,NewAbsolute,DensityLoss,NewDensityLoss,Kappa, &
    !$ACC   CN,KD,SD,SI,ThermC,ThermS,Z0,LiftRows,BInverse,SourceValue,Matrix,RHS,OriginalRHS,Correction, &
    !$ACC   RowScale,ColumnScale,Step,Action,Weights,InitialNorm,Merit,Slope,Fraction,FDStep, &
    !$ACC   FluidSourceError, &
    !$ACC   SpectralRHS,BorderLU,BorderRHS,SpectralPivot,BorderPivot, &
    !$ACC   NativeF,TransferG,TrialF,TrialG,RetainedStep,Velocity,Independent,TrialIndependent,NewtonStep,SearchCounts)
#endif

    END IF
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to:A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices, &
    !$OMP   Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0) &
    !$OMP MAP(alloc:X,NativeIncrement,Iterations,RadiationError,ConservationError,CoefficientQueries,Failure)
    !$OMP TARGET UPDATE TO(ThermalNeeded)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices, &
    !$ACC   Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0)
    !$ACC ENTER DATA CREATE(X,NativeIncrement,Iterations,RadiationError,ConservationError,CoefficientQueries,Failure)
    !$ACC UPDATE DEVICE(ThermalNeeded)
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
        CurrentGeometry(B)%FrameEnergy=0.0_DP
        CurrentGeometry(B)%MomentumScale=0.0_DP
        CurrentGeometry(B)%Valid=.TRUE.
        TrialGeometry(B)%Valid=.FALSE.
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

    ! Start from the conservative zero increment. Try bounded coupled
    ! fixed-point iterations before the Newton search.
    Running=.TRUE.
    Refresh=.TRUE.
    Age=3
    Accepted=0
    CALL SendControl
    CALL TrySharedPicard
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
    !$OMP TARGET UPDATE FROM(FluidSourceError,SearchCounts)
    !$OMP TARGET EXIT DATA MAP(from:X,NativeIncrement,Iterations,RadiationError,ConservationError,CoefficientQueries) &
    !$OMP MAP(release:Failure,A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices,Energy,WeightE,Interp,W,U,D0,I0,V0,Fluid0, &
    !$OMP   SpectralRows,BlockMembers,BlockFirst,BlockLength)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE SELF(FluidSourceError,SearchCounts)
    !$ACC EXIT DATA COPYOUT(X,NativeIncrement,Iterations,RadiationError,ConservationError,CoefficientQueries)
    !$ACC EXIT DATA DELETE(Failure,A,EB,Kind,Process,Species,QA,QB,FM,Active,Indices,Energy,WeightE,Interp, &
    !$ACC   W,U,D0,I0,V0,Fluid0,SpectralRows,BlockMembers,BlockFirst,BlockLength)
#endif
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
          AbsoluteSource(:,B),KappaPair(:,B),DensityLoss(:,B),EnergyRatio_Option=EnergyRatio(:,B), &
          NumberSource_Option=CurrentNumberSource(:,B),NonScatteringEnergy_Option=CurrentEnergySource(B))
        Kappa(:,B) = RetainedValue(KappaPair(:,B))
        SourceValue(:,B) = RetainedValue(Source(:,B))
          CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B),Z0(:,:,B),DT, &
            C0(:,B),C1(:,:,B),Source(:,B),Residual(:,B),Physical(:,B),HScale(:,B),GScale(:,B), &
            AbsoluteSource(:,B),NativeIncrement(:,:,B),NativeF(:,B),TransferG(:,B), &
            CollisionIncrement_Option=ProcessIncrement(:,:,B),CollisionMoment_Option=ProcessMoment(:,B), &
            Incoming_Option=Incoming(:,B), &
            LiftRows_Option=LiftRows(:,B),Geometry_Option=CurrentGeometry(B))
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

    SUBROUTINE TrySharedPicard
      INTEGER :: J,Q,PreflightCount,Passed,CheapLoop
      INTEGER :: AttemptCounts(SIZE(Indices)),AcceptedCounts(SIZE(Indices)),SolvedCounts(4)
      LOGICAL :: Eligible(SIZE(Indices)),IsRoot,HadUnresolved
      REAL(DP) :: TrialPhi,OldMerit,OldRadiationError,OldSourceError(K),Contraction
      INTEGER :: StoppedCount,SourceField
      CALL CheckConvergence
      Eligible=Running
      AttemptCounts=0
      AcceptedCounts=0
      SolvedCounts=0
      StoppedCount=0
      DO CheapLoop=1,MIN(4,MaxIterations)
      Requested=Eligible .AND. Running
      IF( .NOT. ANY(Requested) ) EXIT
      Info=0
      DO B=1,NB
        IF( .NOT. Requested(B) ) CYCLE
        Trial(:,B)=X(:,B)
        DirectionIncrement(:,1,B)=DT*C0(:,B)
        DirectionIncrement(:,2:4,B)=DT*C1(:,:,B)
        DirectionMoment(:,B)=DT*Source(:,B)
        DirectionNumber(:,B)=DT*CurrentNumberSource(:,B)
        DirectionEnergy(B)=DT*CurrentEnergySource(B)
        CALL ProjectSharedProcessEndpoint(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Cap, &
          A,EB,Kind,Gain(:,B),Loss(:,B),Chi(:,B),Sigma(:,B),X(:,B),DirectionIncrement(:,:,B),DirectionMoment(:,B), &
          Trial(:,B),Info(B),NumberIncrement_Option=DirectionNumber(:,B),NonScatteringEnergy_Option=DirectionEnergy(B), &
          Incoming_Option=Incoming(:,B),IncomingEnergy_Option=IncomingEnergy(B), &
          Response_Option=RecoveryResponse(:,:,:,B),EnergyRatio_Option=EnergyRatio(:,B),Geometry_Option=TrialGeometry(B))
        IF( Info(B)/=0 )THEN
          Requested(B)=.FALSE.
          Eligible(B)=.FALSE.
          CYCLE
        END IF
        ProbeTransfer(:,B)=0.0_DP
        DO J=1,NA
          ProbeTransfer(Active(J),B)=Trial(N+J,B)
        END DO
      END DO
      PreflightCount=COUNT(Requested)
      IF( PreflightCount==0 ) EXIT
      CALL Query(0,Indices,Requested,ProbeTransfer,Process,Species,QA,QB, &
        ProbeKT,ProbeMu,ProbeChiMu,ProbeChi,ProbeSigma,ProbeGain,ProbeLoss,Errors)
      CALL MakeTargets(.TRUE.)
      Passed=0
      DO B=1,NB
        IF( .NOT. Requested(B) ) CYCLE
        CoefficientQueries(B)=CoefficientQueries(B)+1
        AttemptCounts(B)=AttemptCounts(B)+1
        IF( Errors(B)/=0 )THEN
          Eligible(B)=.FALSE.
          CYCLE
        END IF
        DRadiation(:,B)=D0(:,B)+Trial(1:M,B)
        DO J=1,3
          IRadiation(:,J,B)=I0(:,J,B)+Trial(J*M+1:(J+1)*M,B)
        END DO
        CALL ComputeCollisionRates_OrderV_Point(M,NF,A,EB,Kind,FM,W(:,B),U(:,B),ProbeGain(:,B),ProbeLoss(:,B), &
          ProbeChi(:,B),ProbeSigma(:,B),ProbeTarget(:,B),DRadiation(:,B),IRadiation(:,:,B), &
          NewC0(:,B),NewC1(:,:,B),NewSource(:,B),NewAbsolute(:,B),NewKappa(:,B),NewDensityLoss(:,B), &
          EnergyRatio_Option=EnergyRatio(:,B),NumberSource_Option=TrialNumber(:,B), &
          NonScatteringEnergy_Option=TrialEnergy(B))
        CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),Trial(:,B), &
          Z0(:,:,B),DT,NewC0(:,B),NewC1(:,:,B),NewSource(:,B),NewResidual(:,B),NewPhysical(:,B),NewHScale(:,B), &
          NewGScale(:,B),NewAbsolute(:,B),ScratchIncrement(:,:,B),TrialF(:,B),TrialG(:,B), &
          CollisionIncrement_Option=DirectionIncrement(:,:,B),CollisionMoment_Option=DirectionMoment(:,B), &
          Incoming_Option=Incoming(:,B),LiftRows_Option=LiftRows(:,B),Geometry_Option=TrialGeometry(B))
        Velocity(:,B)=Retain(V0(:,B))
        DO J=1,NA
          IF( Active(J)>NF+1 ) Velocity(Active(J)-NF-1,B)=Velocity(Active(J)-NF-1,B)-Trial(N+J,B)
        END DO
        CALL ComputeIndependentRadiationResidual(M,NF,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          TrialF(:,B),TrialG(:,B),TrialIndependent(:,B),LiftRows_Option=LiftRows(:,B))
        OldMerit=Merit(B)
        OldRadiationError=RadiationError(B)
        OldSourceError=FluidSourceError(:,B)
        TrialPhi=ComputeCoupledMerit(N,Weights(:,B),TrialIndependent(:,B),NewPhysical(:,B))
        CALL ComputeIndependentMetrics(NE,NS,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),Trial(:,B),DT,NewPhysical(:,B), &
          TrialIndependent(:,B),RetainedValue(NewKappa(:,B)),NewDensityLoss(:,B),NewHScale(:,B),NewGScale(:,B), &
          InitialNorm(:,B),.FALSE.,RadiationError(B),ConservationError(B),Weights(:,B),Merit(B))
        CALL ComputeCoupledSourceMetrics(M,NF,Active,Fluid0(:,B),Trial(:,B),NewPhysical(:,B),NewGScale(:,B), &
          Tolerance,FluidSourceTolerance,TrialIndependent(:,B),FluidSourceError(:,B),Weights(:,B),Merit(B))
        IsRoot=RadiationError(B)<Tolerance .AND. MAXVAL(FluidSourceError(:,B))<FluidSourceTolerance &
            .AND. ConservationError(B)<1.0_DP
        IF( IsRoot .OR. (IEEE_IS_FINITE(TrialPhi) .AND. TrialPhi<OldMerit) )THEN
          AcceptedCounts(B)=AcceptedCounts(B)+1
          IF( IsRoot )THEN
            Passed=Passed+1
            SolvedCounts(CheapLoop)=SolvedCounts(CheapLoop)+1
          END IF
          X(:,B)=Trial(:,B)
          CurrentGeometry(B)=TrialGeometry(B)
          ProcessIncrement(:,:,B)=DirectionIncrement(:,:,B)
          ProcessMoment(:,B)=DirectionMoment(:,B)
          ProcessNumber(:,B)=DirectionNumber(:,B)
          ProcessEnergy(B)=DirectionEnergy(B)
          Target(:,B)=ProbeTarget(:,B)
          Chi(:,B)=ProbeChi(:,B)
          Sigma(:,B)=ProbeSigma(:,B)
          Gain(:,B)=ProbeGain(:,B)
          Loss(:,B)=ProbeLoss(:,B)
          KT(B)=ProbeKT(B)
          Mu(:,B)=ProbeMu(:,B)
          Susceptibility(B)=ProbeChiMu(B)
          Residual(:,B)=NewResidual(:,B)
          Physical(:,B)=NewPhysical(:,B)
          Independent(:,B)=TrialIndependent(:,B)
          NativeF(:,B)=TrialF(:,B)
          TransferG(:,B)=TrialG(:,B)
          HScale(:,B)=NewHScale(:,B)
          GScale(:,B)=NewGScale(:,B)
          AbsoluteSource(:,B)=NewAbsolute(:,B)
          DensityLoss(:,B)=NewDensityLoss(:,B)
          KappaPair(:,B)=NewKappa(:,B)
          Kappa(:,B)=RetainedValue(NewKappa(:,B))
          Source(:,B)=NewSource(:,B)
          SourceValue(:,B)=RetainedValue(NewSource(:,B))
          C0(:,B)=NewC0(:,B)
          C1(:,:,B)=NewC1(:,:,B)
          CurrentNumberSource(:,B)=TrialNumber(:,B)
          CurrentEnergySource(B)=TrialEnergy(B)
          NativeIncrement(:,:,B)=ScratchIncrement(:,:,B)
          Accepted(B)=Accepted(B)+1
          Iterations(B)=Accepted(B)
          IF( IsRoot )THEN
            Running(B)=.FALSE.
          ELSE
            Contraction=0.0_DP
            HadUnresolved=.FALSE.
            IF( OldRadiationError>Tolerance )THEN
              Contraction=RadiationError(B)/MAX(OldRadiationError,Tolerance)
              HadUnresolved=.TRUE.
            END IF
            DO J=1,NA
              SourceField=Active(J)
              IF( OldSourceError(SourceField)<=FluidSourceTolerance ) CYCLE
              Contraction=MAX(Contraction,FluidSourceError(SourceField,B) &
                /MAX(OldSourceError(SourceField),FluidSourceTolerance))
              HadUnresolved=.TRUE.
            END DO
            Eligible(B)=HadUnresolved .AND. Contraction<0.1_DP
            IF( .NOT. Eligible(B) ) StoppedCount=StoppedCount+1
          END IF
          ComponentsValid(B)=.FALSE.
        ELSE
          Eligible(B)=.FALSE.
          CALL ComputeIndependentMetrics(NE,NS,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),X(:,B),DT,Physical(:,B), &
            Independent(:,B),Kappa(:,B),DensityLoss(:,B),HScale(:,B),GScale(:,B),InitialNorm(:,B),.FALSE., &
            RadiationError(B),ConservationError(B),Weights(:,B),Merit(B))
          CALL ComputeCoupledSourceMetrics(M,NF,Active,Fluid0(:,B),X(:,B),Physical(:,B),GScale(:,B), &
            Tolerance,FluidSourceTolerance,Independent(:,B),FluidSourceError(:,B),Weights(:,B),Merit(B))
          Velocity(:,B)=Retain(V0(:,B))
          DO J=1,NA
            IF( Active(J)>NF+1 ) Velocity(Active(J)-NF-1,B)=Velocity(Active(J)-NF-1,B)-X(N+J,B)
          END DO
        END IF
      END DO
      END DO
      CALL SendControl
    END SUBROUTINE TrySharedPicard

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
        CALL MakeTargets(.TRUE.)
        CALL StoreProbeDerivatives
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
            Requested(B)=Requested(B) .AND. Errors(B)/=0
            IF( .NOT. Requested(B) ) CYCLE
            FDStep(B) = -FDStep(B)
            ProbeTransfer(Field,B) = Transfer(Field,B)+FDStep(B)
          END DO
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM(Requested)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE SELF(Requested)
#endif
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
            Failure = 5
            RETURN
          END IF
          CALL MakeTargets(.TRUE.)
          CALL StoreProbeDerivatives
        END IF
      END DO
    END SUBROUTINE RefreshDerivatives

    SUBROUTINE StoreProbeDerivatives
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO
#endif
        DO B = 1, NB
          IF( .NOT. Requested(B) .OR. Errors(B)/=0 ) CYCLE
          DTarget(:,Therm,B) = (ProbeTarget(:,B)-Target(:,B))/FDStep(B)
          DChi(:,Therm,B) = (ProbeChi(:,B)-Chi(:,B))/FDStep(B)
          DSigma(:,Therm,B) = (ProbeSigma(:,B)-Sigma(:,B))/FDStep(B)
          DGain(:,Therm,B) = (ProbeGain(:,B)-Gain(:,B))/FDStep(B)
          DLoss(:,Therm,B) = (ProbeLoss(:,B)-Loss(:,B))/FDStep(B)
        END DO
    END SUBROUTINE StoreProbeDerivatives

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
        RHS(:,B) = 0.0_DP
        Info(B) = 0
        IF( .NOT. Running(B) )THEN
          Matrix(:,:,B)=0.0_DP
          DO J = 1, M
            Matrix(J,J,B) = 1.0_DP
          END DO
          DO J=1,K
            Matrix(M+J,M+NA+J,B)=1.0_DP
          END DO
          DO J=1,NA
            Matrix(M+K+J,M+J,B)=1.0_DP
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
          LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B))
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
          LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B), &
          SearchWeightsOnly_Option=.NOT. FullDirectionAction)
        Step(N+NA+1:NX,B)=0.0_DP
        NewtonStep(:,B)=Step(:,B)
        RetainedStep(:,B)=Retain(Step(:,B))
        IF( FullDirectionAction )THEN
        ! Expansion already returned this action with lift directions disabled.
        VCWork=0.0_DP
        DO J=1,NA
          IF( Active(J)>NF+1 ) VCWork(Active(J)-NF-1)=-Step(N+J,B)
        END DO
        CALL ComputeIndependentRadiationColumn(M,NF,Active,W(:,B),U(:,B),Velocity(:,B),Z0(:,:,B), &
          TransferG(:,B),VCWork,Action(:,B), &
            LiftRows_Option=LiftRows(:,B))
        END IF
        CALL PrepareCoupledSearchWeights(M,NF,Active,Fluid0(:,B),X(:,B), Physical(:,B),Independent(:,B),Step(:,B),Action(:,B),Tolerance,FluidSourceTolerance,Weights(:,B),Merit(B))
        Slope(B)=0.0_DP
        IF( FullDirectionAction ) Slope(B)=ComputeCoupledSlope(N,Weights(:,B),Independent(:,B),Physical(:,B),Action(:,B))
      END DO
      CALL TimersStop(Timer_Newton_Jacobian)
      CALL WriteLinearTrace('newton')
    END SUBROUTINE AssembleAndSolve

    SUBROUTINE SolveBorderedSystem
      INTEGER :: Refine, BorderSize
      REAL(DP), POINTER, CONTIGUOUS :: BorderMatrices(:,:)
      BorderSize = NA
      BorderMatrices(1:NA,1:NA*NB)=>BorderLU
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
        BorderLU(:,:,B) = Matrix(M+K+1:NC,M+1:M+NA,B) &
          -MATMUL(Matrix(M+K+1:NC,1:M,B),SpectralRHS(:,2:BorderSize+1,B))
        BorderRHS(:,B) = RHS(M+K+1:NC,B)-MATMUL(Matrix(M+K+1:NC,1:M,B),SpectralRHS(:,1,B))
      END DO
      CALL LinearFactorBatched(BorderSize,BorderMatrices,BorderSize,BorderPivot,Info,NB,VendorStatus)
      CALL ReadInfo
      IF( VendorStatus /= 0 .OR. ANY(Info /= 0) )THEN
        Failure = 6
        RETURN
      END IF
      CALL LinearApplyBatched('N',BorderSize,1,BorderMatrices,BorderSize,BorderPivot,BorderRHS,BorderSize,NB,VendorStatus)
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
        RHS(M+1:M+NA,B) = BorderRHS(:,B)
        RHS(M+NA+1:NC,B) = 0.0_DP
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
          BorderRHS(:,B) = Correction(M+K+1:NC,B)-MATMUL(Matrix(M+K+1:NC,1:M,B),Correction(1:M,B))
        END DO
        CALL LinearApplyBatched('N',BorderSize,1,BorderMatrices,BorderSize,BorderPivot,BorderRHS,BorderSize,NB,VendorStatus)
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
          Correction(M+1:M+NA,B) = BorderRHS(:,B)
          Correction(M+NA+1:NC,B) = 0.0_DP
          RHS(:,B) = RHS(:,B)+Correction(:,B)
        END DO
      END DO
      ! CN and KD remain unmodified in their compact block-row storage.
      ! Both refinement passes still use the original full scaled matrix.
    END SUBROUTINE SolveBorderedSystem

    SUBROUTINE FactorSpectralBlocks
      INTEGER :: G,BS,GF,GN,LB,RB,PB,IB,TT,BB,Block,II,JJ,Row,Col,Offset
      REAL(DP), POINTER, CONTIGUOUS :: Matrices(:,:),Vectors(:,:)
      DO G=1,NGroups
        BS=GroupSize(G)
        GF=GroupFirst(G)
        GN=GroupCount(G)*NB
        LB=GroupLU(G)
        RB=GroupRHS(G)
        PB=GroupVector(G)
        IB=GroupInfo(G)
        Matrices(1:BS,1:BS*GN)=>SpectralLU(LB:LB+BS*BS*GN-1)
        Vectors(1:BS,1:(NA+1)*GN)=>BlockRHS(RB:RB+BS*(NA+1)*GN-1)
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
            BlockRHS(RB+(TT-1)*BS*(NA+1)+II-1)=RHS(Row,BB)
            DO JJ=1,NA
              BlockRHS(RB+(TT-1)*BS*(NA+1)+JJ*BS+II-1)=Matrix(Row,M+JJ,BB)
            END DO
          END DO
        END DO
        CALL LinearFactorBatched(BS,Matrices,BS,SpectralPivot(PB:PB+BS*GN-1),BlockInfo(IB:IB+GN-1),GN,VendorStatus)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(BlockInfo)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE SELF(BlockInfo)
#endif
        IF( VendorStatus/=0 .OR. ANY(BlockInfo(IB:IB+GN-1)/=0) )THEN
          Failure=6
          RETURN
        END IF
        CALL LinearApplyBatched('N',BS,NA+1,Matrices,BS,SpectralPivot(PB:PB+BS*GN-1), &
          Vectors,BS,GN,VendorStatus)
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
          DO JJ=1,NA+1
          DO II=1,BS
            Row=BlockMembers(Offset+II)
            SpectralRHS(Row,JJ,BB)=BlockRHS(RB+(TT-1)*BS*(NA+1)+(JJ-1)*BS+II-1)
          END DO
          END DO
        END DO
      END DO
    END SUBROUTINE FactorSpectralBlocks

    SUBROUTINE ApplySpectralCorrection
      INTEGER :: G,BS,GF,GN,LB,PB,TT,BB,Block,II,Row,Offset
      REAL(DP), POINTER, CONTIGUOUS :: Matrices(:,:),Vectors(:,:)
      DO G=1,NGroups
        BS=GroupSize(G)
        GF=GroupFirst(G)
        GN=GroupCount(G)*NB
        LB=GroupLU(G)
        PB=GroupVector(G)
        Matrices(1:BS,1:BS*GN)=>SpectralLU(LB:LB+BS*BS*GN-1)
        Vectors(1:BS,1:GN)=>BlockCorrection(PB:PB+BS*GN-1)
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
        CALL LinearApplyBatched('N',BS,1,Matrices,BS,SpectralPivot(PB:PB+BS*GN-1),Vectors,BS,GN,VendorStatus)
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
            Components_Option=Components(:,B),Geometry_Option=TrialGeometry(B))
        IF( NF>=2 ) TrialEnergy(B)=DirectionMoment(2,B)-TrialEnergy(B)
        SharedEndpointValid(B)=Info(B)==0
        SharedProjectionInfo(B)=Info(B)
        ScratchIncrement(:,:,B)=DirectionIncrement(:,:,B)-ProcessIncrement(:,:,B)
        TrialMoment(:,B)=DirectionMoment(:,B)-ProcessMoment(:,B)
        CALL ComputeSharedProcessTangent(M,NF,Active,U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
          ScratchIncrement(:,:,B),TrialMoment(:,B),RetainedStep(:,B),Info(B), &
            Response_Option=RecoveryResponse(:,:,:,B),Incoming_Option=Incoming(:,B))
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
            Response_Option=RecoveryResponse(:,:,:,B),Geometry_Option=TrialGeometry(B))
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
            LiftRows_Option=LiftRows(:,B),Geometry_Option=TrialGeometry(B))
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
              CurrentGeometry(B)=TrialGeometry(B)
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
      IF( .NOT. ASSOCIATED(DirectionScale) )THEN
        ALLOCATE(DirectionScale(N+NA,NB),Gradient(N+NA,NB),ConstraintJacobian(K,N+NA,NB))
        Workspace%DirectionScale=>DirectionScale
        Workspace%Gradient=>Gradient
        Workspace%ConstraintJacobian=>ConstraintJacobian
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
        CALL ComputeCollisionConservationJacobian_OrderV_Point(M,NF,Active,W(:,B),U(:,B), &
          D0(:,B),I0(:,:,B),V0(:,B),X(:,B),ConstraintJacobian(:,:,B))
        DO J=N+1,N+NA
          CALL ComputeCollisionLinearColumn_OrderV_Point(M,NF,Active,W(:,B),U(:,B),D0(:,B),I0(:,:,B),V0(:,B),X(:,B), &
            Z0(:,:,B),DT,SourceValue(:,B),Kappa(:,B),CN(:,:,B),KD(:,:,B),SD(:,:,B),SI(:,:,B), &
            ThermC(:,:,:,B),ThermS(:,:,B),J,.FALSE.,ColumnWork, &
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B))
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
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B))
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
            LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B))
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
          LiftRows_Option=LiftRows(:,B),SpectralRows_Option=SpectralRows,Incoming_Option=Incoming(:,B))
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

