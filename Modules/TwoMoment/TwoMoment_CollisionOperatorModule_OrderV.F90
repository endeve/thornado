MODULE TwoMoment_CollisionOperatorModule_OrderV

  USE KindModule, ONLY: DP
  USE TwoMoment_RetainedArithmeticModule

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: BuildCollisionEvents_OrderV, ComputeCollisionPressure_OrderV_Point, ComputeCollisionPressureIncrement_OrderV_Point
  PUBLIC :: ComputeCollisionClosureDerivatives_OrderV_Point, ComputeEquilibriumTarget_OrderV_Point
  PUBLIC :: ComputeCollisionRates_OrderV_Point, ComputeCollisionRateDerivatives_OrderV_Point
  PUBLIC :: ComputeCollisionRatesAndDerivatives_OrderV_Point, ComputeCollisionDensityJacobian_OrderV_Point
  PUBLIC :: ComputeCollisionResidual_OrderV_Point, ComputeCollisionColumn_OrderV_Point
  PUBLIC :: SelectCollisionMomentSupport_OrderV_Point, CollisionStateAdmissible_OrderV_Point
  PUBLIC :: ComputeCollisionLinearColumn_OrderV_Point, ComputeCollisionLinearAction_OrderV_Point, AssembleCollisionCondensedSystem_OrderV_Point
  PUBLIC :: AssembleCollisionCondensedSystem
  PUBLIC :: ExpandCollisionCondensedStep_OrderV_Point, ComputeCollisionMetrics_OrderV_Point, ComputeCollisionStepBound_OrderV_Point
  PUBLIC :: RefineCollisionConservation_OrderV_Point
  PUBLIC :: ComputeCollisionDensityIncrement_OrderV_Point
  PUBLIC :: ComputeCollisionRateDirection_OrderV_Point
  PUBLIC :: IncomingRadiationCache, PrepareIncomingRadiationCache
  PUBLIC :: BuildSpectralBlocks
  PUBLIC :: PrepareSharedCollisionLinearSystem

  ! No pointers or allocatable components: one device-resident record per
  ! incoming spectral bin. Keep both words of every retained quantity.
  TYPE :: IncomingRadiationCache
    TYPE(RetainedReal) :: FluxFactor(3), FluxFactorSquared, FluxFactorNorm, Closure
    TYPE(RetainedReal) :: Pressure(3,3), Number, NumberFrame(3), FluxFrame(3,3)
  END TYPE IncomingRadiationCache

CONTAINS

  PURE SUBROUTINE PrepareSharedCollisionLinearSystem(M,NF,NA,Matrix,RHS)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,NA
    REAL(DP), INTENT(inout) :: Matrix(:,:),RHS(:)
    INTEGER :: J,K
    K=NF+4
    ! The shared-process state already constructs conservation from paired
    ! increments. Aggregate moment constraints and their auxiliary radiation
    ! lifts are redundant here; their cancellation error can corrupt a tiny
    ! allowed source. Solve the collision and fluid-source linearization with
    ! zero lift, retaining the existing condensed/batched storage layout.
    Matrix(:,M+NA+1:M+NA+K)=0.0_DP
    DO J=1,K
      Matrix(M+J,:)=0.0_DP
      Matrix(M+J,M+NA+J)=1.0_DP
      RHS(M+J)=0.0_DP
    END DO
  END SUBROUTINE PrepareSharedCollisionLinearSystem

  SUBROUTINE BuildSpectralBlocks(M,A,B,Rows,Members,First,Length)
    INTEGER, INTENT(in) :: M,A(:),B(:)
    INTEGER, ALLOCATABLE, INTENT(out) :: Rows(:,:),Members(:),First(:),Length(:)
    INTEGER :: Label(M),Roots(M),Counts(M),Q,E,L,R,J,N,Offset,Swap,MatrixOffset
    DO Q=1,M
      Label(Q)=Q
    END DO
    ! Use every structurally enabled edge, including zero-valued kernels.
    ! No coefficient threshold or independent thermodynamic approximation.
    DO E=1,SIZE(A)
      L=A(E)
      DO WHILE( Label(L)/=L )
        L=Label(L)
      END DO
      R=B(E)
      DO WHILE( Label(R)/=R )
        R=Label(R)
      END DO
      Label(MAX(L,R))=MIN(L,R)
    END DO
    DO Q=1,M
      R=Q
      DO WHILE( Label(R)/=R )
        R=Label(R)
      END DO
      Label(Q)=R
    END DO
    N=0
    DO Q=1,M
      IF( Label(Q)/=Q ) CYCLE
      N=N+1
      Roots(N)=Q
      Counts(N)=COUNT(Label==Q)
    END DO
    ! Equal-sized blocks share one vendor batch. Members retain their
    ! original spectral ordering, including non-contiguous flavor blocks.
    DO L=1,N
    DO R=L+1,N
      IF( Counts(R)<=Counts(L) ) CYCLE
      Swap=Counts(L); Counts(L)=Counts(R); Counts(R)=Swap
      Swap=Roots(L); Roots(L)=Roots(R); Roots(R)=Swap
    END DO
    END DO
    ALLOCATE(Rows(3,M),Members(M),First(N),Length(N))
    Offset=0
    MatrixOffset=0
    DO L=1,N
      First(L)=Offset+1
      Length(L)=Counts(L)
      J=0
      DO Q=1,M
        IF( Label(Q)/=Roots(L) ) CYCLE
        J=J+1
        Members(Offset+J)=Q
        Rows(:,Q)=[J,L,MatrixOffset+(J-1)*Counts(L)]
      END DO
      Offset=Offset+J
      MatrixOffset=MatrixOffset+J*J
    END DO
  END SUBROUTINE BuildSpectralBlocks

  PURE FUNCTION SpectralEntry(Matrix,I,J,Rows_Option) RESULT(Value)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: Matrix(:,:)
    INTEGER, INTENT(in) :: I,J
    INTEGER, INTENT(in), OPTIONAL :: Rows_Option(:,:)
    REAL(DP) :: Value
    Value=0.0_DP
    IF( PRESENT(Rows_Option) )THEN
      IF( Rows_Option(2,I)/=Rows_Option(2,J) ) RETURN
      Value=Matrix(Rows_Option(3,J)+Rows_Option(1,I),1)
    ELSE
      Value=Matrix(I,J)
    END IF
  END FUNCTION SpectralEntry

  PURE SUBROUTINE PrepareIncomingRadiationCache(M,U,D0,I0,V0,Incoming,IncomingEnergy)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M
    REAL(DP), INTENT(in) :: U(M),D0(M),I0(M,3),V0(3)
    TYPE(IncomingRadiationCache), INTENT(out) :: Incoming(M)
    TYPE(RetainedReal), INTENT(out) :: IncomingEnergy
    INTEGER :: Q,J,L
    IncomingEnergy=0.0_DP
    DO Q=1,M
      Incoming(Q)%FluxFactorSquared=0.0_DP
      DO J=1,3
        Incoming(Q)%FluxFactor(J)=Retain(I0(Q,J))/D0(Q)
        Incoming(Q)%FluxFactorSquared=Incoming(Q)%FluxFactorSquared &
          +Incoming(Q)%FluxFactor(J)*Incoming(Q)%FluxFactor(J)
      END DO
      Incoming(Q)%FluxFactorNorm=RetainedSqrt(Incoming(Q)%FluxFactorSquared)
      Incoming(Q)%Closure=(3.0_DP-Incoming(Q)%FluxFactorNorm &
        +3.0_DP*Incoming(Q)%FluxFactorNorm*Incoming(Q)%FluxFactorNorm)/5.0_DP
      CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),Incoming(Q)%Pressure)
      Incoming(Q)%Number=Retain(D0(Q))
      DO J=1,3
        Incoming(Q)%NumberFrame(J)=Retain(V0(J))*I0(Q,J)
        Incoming(Q)%Number=Incoming(Q)%Number+Incoming(Q)%NumberFrame(J)
        DO L=1,3
          Incoming(Q)%FluxFrame(J,L)=Retain(V0(L))*Incoming(Q)%Pressure(J,L)
        END DO
      END DO
      IncomingEnergy=IncomingEnergy+U(Q)*Retain(D0(Q))
    END DO
  END SUBROUTINE PrepareIncomingRadiationCache

  SUBROUTINE BuildCollisionEvents_OrderV &
    ( NE, NS, Enabled, A, B, Kind, Process, Species, QA, QB, FamilyMoment )

    INTEGER, INTENT(in) :: NE, NS
    LOGICAL, INTENT(in) :: Enabled(9)
    INTEGER, ALLOCATABLE, INTENT(out) :: A(:), B(:), Kind(:), Process(:)
    INTEGER, ALLOCATABLE, INTENT(out) :: Species(:), QA(:), QB(:), FamilyMoment(:,:)
    INTEGER :: P, S, R, Q, N, J, Partner, SignB, F

    N = 0
    DO P = 1, 9
      IF( .NOT. Enabled(P) ) CYCLE
      IF( P <= 2 .OR. P == 6 )THEN
        N = N + NS * NE * ( NE - 1 ) / 2
      ELSE IF( P <= 5 )THEN
        N = N + NS / 2 * NE * NE
      ELSE IF( NS >= 4 )THEN
        N = N + 2 * NE * NE
      END IF
    END DO
    ALLOCATE(A(N),B(N),Kind(N),Process(N),Species(N),QA(N),QB(N))
    ALLOCATE(FamilyMoment(NS/2,N))
    FamilyMoment = 0
    J = 0
    DO P = 1, 9
      IF( .NOT. Enabled(P) ) CYCLE
      DO S = 1, NS
        IF( P >= 3 .AND. P <= 5 .AND. MOD(S,2) == 0 ) CYCLE
        IF( P >= 7 .AND. ( S > 2 .OR. NS < 4 ) ) CYCLE
        Partner = S
        IF( P >= 3 .AND. P <= 5 ) Partner = S + 1
        IF( P == 7 .OR. P == 8 ) Partner = S + 2
        IF( P == 9 ) Partner = 5 - S
        DO R = 1, NE
        DO Q = 1, NE
          IF( ( P <= 2 .OR. P == 6 ) .AND. Q >= R ) CYCLE
          J = J + 1
          A(J) = ( S - 1 ) * NE + R
          B(J) = ( Partner - 1 ) * NE + Q
          Kind(J) = 0
          IF( (P >= 3 .AND. P <= 5) .OR. P == 9 ) Kind(J) = 1
          Process(J) = P
          Species(J) = S
          QA(J) = R
          QB(J) = Q
          SignB = -1
          IF( Kind(J) == 1 ) SignB = 1
          F = ( S + 1 ) / 2
          ! Integer stoichiometry, including both partners, before any sum.
          FamilyMoment(F,J) = MERGE(1,-1,MOD(S,2)==1)
          F = ( Partner + 1 ) / 2
          FamilyMoment(F,J) = FamilyMoment(F,J) &
            + SignB * MERGE(1,-1,MOD(Partner,2)==1)
        END DO
        END DO
      END DO
    END DO

  END SUBROUTINE BuildCollisionEvents_OrderV

  PURE SUBROUTINE ComputeCollisionPressure_OrderV_Point( D, I, K )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: D, I(:)
    TYPE(RetainedReal), INTENT(out) :: K(3,3)
    TYPE(RetainedReal) :: R2, H, C, P, Z(3)
    INTEGER :: J, L
    R2 = 0.0_DP
    DO J = 1, 3
      Z(J) = I(J)/D
      R2 = R2 + Z(J)*Z(J)
    END DO
    H = RetainedSqrt(R2)
    C = ( 3.0_DP - H + 3.0_DP*H*H ) / 5.0_DP
    DO L = 1, 3
    DO J = 1, 3
      P = Z(J)*Z(L)
      IF( J == L ) P = P - R2/3.0_DP
      K(J,L) = C*P
      IF( J == L ) K(J,L) = K(J,L) + Retain(1.0_DP)/3.0_DP
      K(J,L) = D*K(J,L)
    END DO
    END DO
  END SUBROUTINE ComputeCollisionPressure_OrderV_Point

  PURE SUBROUTINE ComputeCollisionPressureIncrement_OrderV_Point( D0, I0, DD, DI, DK, Incoming_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: D0, I0(:)
    TYPE(RetainedReal), INTENT(in) :: DD, DI(3)
    TYPE(RetainedReal), INTENT(out) :: DK(3,3)
    TYPE(IncomingRadiationCache), INTENT(in), OPTIONAL :: Incoming_Option
    TYPE(RetainedReal) :: D, R02, R2, DR2, H0, DH, H, C0, C, DC, P0, DP0, Z0(3), Z(3), DZ(3), Tensor0
    INTEGER :: J, L
    D = D0 + DD
    IF( PRESENT(Incoming_Option) )THEN
      Z0=Incoming_Option%FluxFactor
      R02=Incoming_Option%FluxFactorSquared
      H0=Incoming_Option%FluxFactorNorm
      C0=Incoming_Option%Closure
    ELSE
      R02 = 0.0_DP
      DO J=1,3
        Z0(J)=Retain(I0(J))/D0
        R02=R02+Z0(J)*Z0(J)
      END DO
      H0=RetainedSqrt(R02)
      C0=(3.0_DP-H0+3.0_DP*H0*H0)/5.0_DP
    END IF
    DR2 = 0.0_DP
    DO J = 1, 3
      DZ(J) = (DI(J)-Z0(J)*DD)/D
      Z(J) = (I0(J)+DI(J))/D
      DR2 = DR2+2.0_DP*Z0(J)*DZ(J)+DZ(J)*DZ(J)
    END DO
    R2 = 0.0_DP
    DO J = 1, 3
      R2 = R2+Z(J)*Z(J)
    END DO
    H = RetainedSqrt(R2)
    DH = 0.0_DP
    IF( RetainedValue(H+H0) > 0.0_DP ) DH = DR2/(H+H0)
    C = (3.0_DP-H+3.0_DP*H*H)/5.0_DP
    DC = DH*(-1.0_DP+3.0_DP*(2.0_DP*H0+DH))/5.0_DP
    DO L = 1, 3
    DO J = 1, 3
      P0 = Z0(J)*Z0(L)
      DP0 = Z0(J)*DZ(L)+DZ(J)*Z0(L)+DZ(J)*DZ(L)
      IF( J == L )THEN
        P0 = P0 - R02/3.0_DP
        DP0 = DP0 - DR2/3.0_DP
      END IF
      Tensor0 = C0*P0
      IF( J == L ) Tensor0 = Tensor0+Retain(1.0_DP)/3.0_DP
      DK(J,L) = Tensor0*DD+D*(C*DP0+DC*P0)
    END DO
    END DO
  END SUBROUTINE ComputeCollisionPressureIncrement_OrderV_Point

  PURE SUBROUTINE ComputeCollisionClosureDerivatives_OrderV_Point( D, I, KD, KI )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: D, I(3)
    REAL(DP), INTENT(out) :: KD(3,3), KI(3,3,3)
    REAL(DP) :: R2, R, H, C, CP, P, Normal(3), Z(3)
    INTEGER :: J, L, A
    Z = I/D
    R2 = SUM(Z*Z)
    R = SQRT(R2)
    H = R
    C = ( 3.0_DP - H + 3.0_DP*H*H )/5.0_DP
    CP = ( -1.0_DP + 6.0_DP*H )/5.0_DP
    Normal = 0.0_DP
    IF( R > 0.0_DP ) Normal = Z/R
    DO L = 1, 3
    DO J = 1, 3
      P = Z(J)*Z(L)
      IF( J == L ) P = P - R2/3.0_DP
      KD(J,L) = -( C + H*CP )*P
      IF( J == L ) KD(J,L) = KD(J,L) + 1.0_DP/3.0_DP
      DO A = 1, 3
        KI(J,L,A) = CP*Normal(A)*P
        IF( J == A ) KI(J,L,A) = KI(J,L,A) + C*Z(L)
        IF( L == A ) KI(J,L,A) = KI(J,L,A) + C*Z(J)
        IF( J == L ) KI(J,L,A) = KI(J,L,A) - 2.0_DP*C*Z(A)/3.0_DP
      END DO
    END DO
    END DO
  END SUBROUTINE ComputeCollisionClosureDerivatives_OrderV_Point

  PURE SUBROUTINE ComputeEquilibriumTarget_OrderV_Point( NE, NS, ND, Energy, Weight, Interp, KT, Mu, Target )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: NE, NS, ND
    REAL(DP), INTENT(in) :: Energy(NE), Weight(ND), Interp(ND+2,ND), KT, Mu(NS)
    TYPE(RetainedReal), INTENT(out) :: Target(NE,NS)
    INTEGER :: S, Q
    ! The active target is the thermodynamic Fermi--Dirac distribution.
    ! Do not restore a DG cap here during a naming/organization cleanup.
    DO S = 1, NS
      DO Q = 1, NE
        Target(Q,S) = RetainedFermi((Retain(Energy(Q))-Mu(S))/KT)
      END DO
    END DO
  END SUBROUTINE ComputeEquilibriumTarget_OrderV_Point

  PURE LOGICAL FUNCTION CollisionStateAdmissible_OrderV_Point( M, D0, I0, X, Cap ) RESULT(Valid)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M
    REAL(DP), INTENT(in) :: D0(M), I0(M,3), Cap
    TYPE(RetainedReal), INTENT(in) :: X(:)
    TYPE(RetainedReal) :: D, R2, I, Slack
    INTEGER :: Q, J
    Valid = .FALSE.
    DO Q = 1, M
      D = D0(Q)+X(Q)
      IF( D % High <= 0.0_DP ) RETURN
      Slack = Cap-D
      IF( Slack % High < -8.0_DP*EPSILON(1.0_DP) ) RETURN
      R2 = 0.0_DP
      DO J = 1, 3
        I = I0(Q,J)+X(J*M+Q)
        ! Reject an obviously invalid component before dividing by a tiny
        ! density. An invalid trial must not overflow its own cone test.
        Slack = (1.0_DP+64.0_DP*EPSILON(1.0_DP))*D-RetainedAbs(I)
        IF( Slack % High < 0.0_DP ) RETURN
        I = I/D
        R2 = R2 + I*I
      END DO
      Slack = 1.0_DP-RetainedSqrt(R2)
      IF( Slack % High < -64.0_DP*EPSILON(1.0_DP) ) RETURN
    END DO
    Valid = .TRUE.
  END FUNCTION CollisionStateAdmissible_OrderV_Point

  PURE SUBROUTINE ComputeCollisionRates_OrderV_Point &
    ( M, NF, A, B, Kind, FamilyMoment, W, U, Gain, Loss, Chi, Sigma, Target, D, I, &
      C0, C1, Source, AbsoluteSource, Kappa, DensityLoss, EnergyRatio_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M, NF, A(:), B(:), Kind(:), FamilyMoment(:,:)
    REAL(DP), INTENT(in) :: W(M), U(M), Gain(:), Loss(:), Chi(M), Sigma(M)
    TYPE(RetainedReal), INTENT(in) :: Target(M), D(M), I(M,3)
    TYPE(RetainedReal), INTENT(out) :: C0(M), C1(M,3), Source(NF+4), Kappa(M)
    REAL(DP), INTENT(out) :: AbsoluteSource(NF+4), DensityLoss(M)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: EnergyRatio_Option(M)
    TYPE(RetainedReal) :: CC, Event, CA, CB, X, Y, Term, GPair, LPair, EM
    INTEGER :: Q, S, F, J, L, IA, IB, NE
    REAL(DP) :: SignB
    NE = M/(2*NF)
    Source = 0.0_DP
    AbsoluteSource = 0.0_DP
    DO Q = 1, M
      CC = Chi(Q)*(Target(Q)-D(Q))
      C0(Q) = CC
      Kappa(Q) = Chi(Q)+Sigma(Q)
      DensityLoss(Q) = Chi(Q)
      S = (Q-1)/NE+1
      F = (S+1)/2
      Term = W(Q)*MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*CC
      Source(F) = Source(F)+Term
      AbsoluteSource(F) = AbsoluteSource(F)+ABS(RetainedValue(Term))
      Term = U(Q)*CC
      Source(NF+1) = Source(NF+1)+Term
      AbsoluteSource(NF+1) = AbsoluteSource(NF+1)+ABS(RetainedValue(Term))
    END DO
    DO L = 1, SIZE(A)
      IA = A(L)
      IB = B(L)
      X = D(IA)
      Y = D(IB)
      GPair = Retain(Gain(L))*W(IA)
      LPair = Retain(Loss(L))*W(IA)
      IF( Kind(L) == 0 )THEN
        Event = GPair*Y*(1.0_DP-X)-LPair*X*(1.0_DP-Y)
        CA = GPair*Y+LPair*(1.0_DP-Y)
        CB = GPair*(1.0_DP-X)+LPair*X
        SignB = -1.0_DP
      ELSE
        Event = GPair*(1.0_DP-X)*(1.0_DP-Y)-LPair*X*Y
        CA = GPair*(1.0_DP-Y)+LPair*Y
        CB = GPair*(1.0_DP-X)+LPair*X
        SignB = 1.0_DP
      END IF
      ! One retained event supplies BOTH radiation legs and allowed sources.
      C0(IA) = C0(IA)+Event/W(IA)
      C0(IB) = C0(IB)+SignB*Event/W(IB)
      Kappa(IA) = Kappa(IA)+CA/W(IA)
      Kappa(IB) = Kappa(IB)+CB/W(IB)
      DensityLoss(IA) = DensityLoss(IA)+RetainedValue(CA)/W(IA)
      DensityLoss(IB) = DensityLoss(IB)+RetainedValue(CB)/W(IB)
      DO F = 1, NF
        IF( FamilyMoment(F,L) == 0 ) CYCLE
        Term = REAL(FamilyMoment(F,L),DP)*Event
        Source(F) = Source(F)+Term
        AbsoluteSource(F) = AbsoluteSource(F)+ABS(RetainedValue(Term))
      END DO
      IF( PRESENT(EnergyRatio_Option) )THEN
        EM = EnergyRatio_Option(IA)+SignB*EnergyRatio_Option(IB)
      ELSE
        EM = Retain(U(IA))/W(IA)+SignB*(Retain(U(IB))/W(IB))
      END IF
      IF( EM%High /= 0.0_DP .OR. EM%Low /= 0.0_DP )THEN
        Term = EM*Event
        Source(NF+1) = Source(NF+1)+Term
        AbsoluteSource(NF+1) = AbsoluteSource(NF+1)+ABS(RetainedValue(Term))
      END IF
    END DO
    DO J = 1, 3
    DO Q = 1, M
      C1(Q,J) = -Kappa(Q)*I(Q,J)
      Term = U(Q)*C1(Q,J)
      Source(NF+1+J) = Source(NF+1+J)+Term
      AbsoluteSource(NF+1+J) = AbsoluteSource(NF+1+J)+ABS(RetainedValue(Term))
    END DO
    END DO
  END SUBROUTINE ComputeCollisionRates_OrderV_Point

  PURE SUBROUTINE ComputeCollisionRateDirection_OrderV_Point &
    ( M, NF, A, B, Kind, FamilyMoment, W, U, Gain, Loss, Chi, Sigma, Target, D, I, &
      DD, DI, DGain, DLoss, ThermDirection, DChi, DSigma, DTarget, C0, C1, Source, &
      NumberSource_Option, NonScatteringEnergy_Option, EnergyRatio_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Evaluate C + DC[d] by reaction.  The very same retained event supplies
    ! both radiation partners and its allowed fluid moments.  In particular,
    ! do not obtain Source by taking moments of the assembled C0/C1 arrays.
    INTEGER, INTENT(in) :: M, NF, A(:), B(:), Kind(:), FamilyMoment(:,:)
    REAL(DP), INTENT(in) :: W(M), U(M), Gain(:), Loss(:), Chi(M), Sigma(M)
    TYPE(RetainedReal), INTENT(in) :: Target(M), D(M), I(M,3), DD(M), DI(M,3)
    REAL(DP), INTENT(in) :: DGain(:,:), DLoss(:,:)
    TYPE(RetainedReal), INTENT(in) :: ThermDirection(3), DChi(M), DSigma(M), DTarget(M)
    TYPE(RetainedReal), INTENT(out) :: C0(M), C1(M,3), Source(NF+4)
    TYPE(RetainedReal), INTENT(out), OPTIONAL :: NumberSource_Option(M)
    TYPE(RetainedReal), INTENT(out), OPTIONAL :: NonScatteringEnergy_Option
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: EnergyRatio_Option(M)
    TYPE(RetainedReal) :: Event, X, Y, DX, DY, G, L, DG, DL, Term, Kappa(M), DKappa(M), EM
    REAL(DP) :: SignB
    INTEGER :: Q, J, F, S, NE, E, IA, IB
    NE = M/(2*NF)
    Source = 0.0_DP
    IF( PRESENT(NonScatteringEnergy_Option) ) NonScatteringEnergy_Option=0.0_DP
    DO Q = 1, M
      Term = Chi(Q)*(Target(Q)-D(Q)) &
               + DChi(Q)*(Target(Q)-D(Q)) + Chi(Q)*(DTarget(Q)-DD(Q))
      C0(Q) = Term
      IF( PRESENT(NumberSource_Option) ) NumberSource_Option(Q)=Term
      Kappa(Q) = Chi(Q)+Retain(Sigma(Q))
      DKappa(Q) = DChi(Q)+DSigma(Q)
      S = (Q-1)/NE+1
      F = (S+1)/2
      Source(F) = Source(F)+MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*(W(Q)*Term)
      Source(NF+1) = Source(NF+1)+U(Q)*Term
      IF( PRESENT(NonScatteringEnergy_Option) ) &
        NonScatteringEnergy_Option=NonScatteringEnergy_Option+U(Q)*Term
    END DO
    DO E = 1, SIZE(A)
      IA = A(E)
      IB = B(E)
      X = D(IA)
      Y = D(IB)
      DX = DD(IA)
      DY = DD(IB)
      G = Retain(Gain(E))*W(IA)
      L = Retain(Loss(E))*W(IA)
      DG = 0.0_DP
      DL = 0.0_DP
      DO J=1,3
        DG = DG+ThermDirection(J)*DGain(E,J)
        DL = DL+ThermDirection(J)*DLoss(E,J)
      END DO
      DG = DG*W(IA)
      DL = DL*W(IA)
      IF( Kind(E) == 0 )THEN
        Event = G*Y*(1.0_DP-X)-L*X*(1.0_DP-Y) &
          + DG*Y*(1.0_DP-X)-DL*X*(1.0_DP-Y) &
          + (G*(1.0_DP-X)+L*X)*DY-(G*Y+L*(1.0_DP-Y))*DX
        SignB = -1.0_DP
        Kappa(IA) = Kappa(IA)+(G*Y+L*(1.0_DP-Y))/W(IA)
        DKappa(IA) = DKappa(IA)+(DG*Y+DL*(1.0_DP-Y)+(G-L)*DY)/W(IA)
      ELSE
        Event = G*(1.0_DP-X)*(1.0_DP-Y)-L*X*Y &
          + DG*(1.0_DP-X)*(1.0_DP-Y)-DL*X*Y &
          - (G*(1.0_DP-Y)+L*Y)*DX-(G*(1.0_DP-X)+L*X)*DY
        SignB = 1.0_DP
        Kappa(IA) = Kappa(IA)+(G*(1.0_DP-Y)+L*Y)/W(IA)
        DKappa(IA) = DKappa(IA)+(DG*(1.0_DP-Y)+DL*Y+(L-G)*DY)/W(IA)
      END IF
      Kappa(IB) = Kappa(IB)+(G*(1.0_DP-X)+L*X)/W(IB)
      DKappa(IB) = DKappa(IB)+(DG*(1.0_DP-X)+DL*X+(L-G)*DX)/W(IB)
      C0(IA) = C0(IA)+Event/W(IA)
      C0(IB) = C0(IB)+SignB*Event/W(IB)
      IF( PRESENT(NumberSource_Option) )THEN
        ! Same-species scattering has zero number source in every connected
        ! spectral component. Preserve this before materializing its two legs.
        IF( (IA-1)/NE/=(IB-1)/NE .OR. Kind(E)/=0 )THEN
          NumberSource_Option(IA)=NumberSource_Option(IA)+Event/W(IA)
          NumberSource_Option(IB)=NumberSource_Option(IB)+SignB*Event/W(IB)
        END IF
      END IF
      DO F = 1, NF
        IF( FamilyMoment(F,E) == 0 ) CYCLE
        Source(F) = Source(F)+REAL(FamilyMoment(F,E),DP)*Event
      END DO
      IF( PRESENT(EnergyRatio_Option) )THEN
        EM = EnergyRatio_Option(IA)+SignB*EnergyRatio_Option(IB)
      ELSE
        EM = Retain(U(IA))/W(IA)+SignB*(Retain(U(IB))/W(IB))
      END IF
      IF( EM%High /= 0.0_DP .OR. EM%Low /= 0.0_DP ) Source(NF+1) = Source(NF+1)+EM*Event
      IF( PRESENT(NonScatteringEnergy_Option) )THEN
        IF( (IA-1)/NE/=(IB-1)/NE .OR. Kind(E)/=0 ) &
          NonScatteringEnergy_Option=NonScatteringEnergy_Option+EM*Event
      END IF
    END DO
    DO Q = 1, M
      DO J = 1, 3
        Term = -Kappa(Q)*I(Q,J)-DKappa(Q)*I(Q,J)-Kappa(Q)*DI(Q,J)
        C1(Q,J) = Term
        Source(NF+1+J) = Source(NF+1+J)+U(Q)*Term
      END DO
    END DO
  END SUBROUTINE ComputeCollisionRateDirection_OrderV_Point

  PURE SUBROUTINE ComputeCollisionRateDerivatives_OrderV_Point &
    ( M, NF, A, B, Kind, FamilyMoment, W, U, Gain, Loss, Chi, D, I, Kappa, CN, KD, SD, SI )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M, NF, A(:), B(:), Kind(:), FamilyMoment(:,:)
    REAL(DP), INTENT(in) :: W(M), U(M), Gain(:), Loss(:), Chi(M), I(M,3), Kappa(M)
    TYPE(RetainedReal), INTENT(in) :: D(M)
    REAL(DP), INTENT(out) :: CN(M,M), KD(M,M), SD(NF+4,M), SI(NF+4,3*M)
    REAL(DP) :: DA, DB, SignB, EM, G, L0
    INTEGER :: Q, S, F, NE, L, IA, IB, J
    CN = 0.0_DP
    KD = 0.0_DP
    SD = 0.0_DP
    SI = 0.0_DP
    NE = M/(2*NF)
    DO Q = 1, M
      CN(Q,Q) = -Chi(Q)
      S = (Q-1)/NE+1
      F = (S+1)/2
      SD(F,Q) = -W(Q)*MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*Chi(Q)
      SD(NF+1,Q) = -U(Q)*Chi(Q)
    END DO
    DO L = 1, SIZE(A)
      IA = A(L)
      IB = B(L)
      G = Gain(L)*W(IA)
      L0 = Loss(L)*W(IA)
      IF( Kind(L) == 0 )THEN
        DA = RetainedValue(-G*D(IB)-L0*(1.0_DP-D(IB)))
        DB = RetainedValue(G*(1.0_DP-D(IA))+L0*D(IA))
        SignB = -1.0_DP
        KD(IA,IB) = KD(IA,IB)+(G-L0)/W(IA)
        KD(IB,IA) = KD(IB,IA)+(L0-G)/W(IB)
      ELSE
        DA = RetainedValue(-G*(1.0_DP-D(IB))-L0*D(IB))
        DB = RetainedValue(-G*(1.0_DP-D(IA))-L0*D(IA))
        SignB = 1.0_DP
        KD(IA,IB) = KD(IA,IB)+(L0-G)/W(IA)
        KD(IB,IA) = KD(IB,IA)+(L0-G)/W(IB)
      END IF
      CN(IA,IA) = CN(IA,IA)+DA/W(IA)
      CN(IA,IB) = CN(IA,IB)+DB/W(IA)
      CN(IB,IA) = CN(IB,IA)+SignB*DA/W(IB)
      CN(IB,IB) = CN(IB,IB)+SignB*DB/W(IB)
      DO F = 1, NF
        IF( FamilyMoment(F,L) == 0 ) CYCLE
        SD(F,IA) = SD(F,IA)+REAL(FamilyMoment(F,L),DP)*DA
        SD(F,IB) = SD(F,IB)+REAL(FamilyMoment(F,L),DP)*DB
      END DO
      EM = U(IA)/W(IA)+SignB*U(IB)/W(IB)
      SD(NF+1,IA) = SD(NF+1,IA)+EM*DA
      SD(NF+1,IB) = SD(NF+1,IB)+EM*DB
    END DO
    DO J = 1, 3
      DO Q = 1, M
        SI(NF+1+J,(J-1)*M+Q) = -U(Q)*Kappa(Q)
        DO L = 1, M
          SD(NF+1+J,L) = SD(NF+1+J,L)-U(Q)*I(Q,J)*KD(Q,L)
        END DO
      END DO
    END DO
  END SUBROUTINE ComputeCollisionRateDerivatives_OrderV_Point

  PURE SUBROUTINE ComputeCollisionDensityJacobian_OrderV_Point(M,A,B,Kind,W,Gain,Loss,Chi,D,CN,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Restore the density Jacobian after its storage was used by vendor LU.
    ! Use exactly the event order and arithmetic of the derivative kernel.
    INTEGER, INTENT(in) :: M,A(:),B(:),Kind(:)
    REAL(DP), INTENT(in) :: W(M),Gain(:),Loss(:),Chi(M)
    TYPE(RetainedReal), INTENT(in) :: D(M)
    REAL(DP), INTENT(out) :: CN(:,:)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    REAL(DP) :: G,L0,DA,DB,SignB
    INTEGER :: Q,L,IA,IB,Row(M),Offset(M),Col(M)
    Row=[(Q,Q=1,M)]
    Offset=0
    Col=Row
    IF( PRESENT(SpectralRows_Option) )THEN
      Row=SpectralRows_Option(1,:)
      Offset=SpectralRows_Option(3,:)
      Col=1
    END IF
    CN=0.0_DP
    DO Q=1,M
      CN(Offset(Q)+Row(Q),Col(Q))=-Chi(Q)
    END DO
    DO L=1,SIZE(A)
      IA=A(L)
      IB=B(L)
      G=Gain(L)*W(IA)
      L0=Loss(L)*W(IA)
      IF( Kind(L)==0 )THEN
        DA=RetainedValue(-G*D(IB)-L0*(1.0_DP-D(IB)))
        DB=RetainedValue(G*(1.0_DP-D(IA))+L0*D(IA))
        SignB=-1.0_DP
      ELSE
        DA=RetainedValue(-G*(1.0_DP-D(IB))-L0*D(IB))
        DB=RetainedValue(-G*(1.0_DP-D(IA))-L0*D(IA))
        SignB=1.0_DP
      END IF
      CN(Offset(IA)+Row(IA),Col(IA))=CN(Offset(IA)+Row(IA),Col(IA))+DA/W(IA)
      CN(Offset(IB)+Row(IA),Col(IB))=CN(Offset(IB)+Row(IA),Col(IB))+DB/W(IA)
      CN(Offset(IA)+Row(IB),Col(IA))=CN(Offset(IA)+Row(IB),Col(IA))+SignB*DA/W(IB)
      CN(Offset(IB)+Row(IB),Col(IB))=CN(Offset(IB)+Row(IB),Col(IB))+SignB*DB/W(IB)
    END DO
  END SUBROUTINE ComputeCollisionDensityJacobian_OrderV_Point

  PURE SUBROUTINE ComputeCollisionRatesAndDerivatives_OrderV_Point &
    (M,NF,A,B,Kind,FamilyMoment,W,U,Gain,Loss,Chi,D,I,IDP,Kappa, &
     DGain,DLoss,DChi,DSigma,Target,DTarget,ThermalNeeded, &
     CN,KD,SD,SI,ThermC,ThermS,C0Work,KappaWork,EnergyRatio_Option,RealEnergyRatio_Option,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,A(:),B(:),Kind(:),FamilyMoment(:,:)
    REAL(DP), INTENT(in) :: W(M),U(M),Gain(:),Loss(:),Chi(M),IDP(M,3),Kappa(M)
    REAL(DP), INTENT(in) :: DGain(:,:),DLoss(:,:),DChi(M,3),DSigma(M,3)
    TYPE(RetainedReal), INTENT(in) :: D(M),I(M,3),Target(M),DTarget(M,3)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: EnergyRatio_Option(M)
    REAL(DP), INTENT(in), OPTIONAL :: RealEnergyRatio_Option(M)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    LOGICAL, INTENT(in) :: ThermalNeeded(3)
    REAL(DP), INTENT(out) :: CN(:,:),KD(:,:),SD(NF+4,M),SI(:,:)
    REAL(DP), INTENT(out) :: ThermC(M,4,3),ThermS(NF+4,3)
    ! Caller-owned device scratch, reused later by candidate evaluation.
    TYPE(RetainedReal), INTENT(out) :: C0Work(M,3),KappaWork(M,3)
    TYPE(RetainedReal) :: TS(NF+4,3),X,Y,HX,HY,GP,LP,Event,CA,CB,Term,EMR,C1
    REAL(DP) :: DA,DB,SignB,EM,G,L0
    INTEGER :: Q,S,F,NE,L,IA,IB,J,H,Row(M),Offset(M),Col(M)
    Row=[(Q,Q=1,M)]
    Offset=0
    Col=Row
    IF( PRESENT(SpectralRows_Option) )THEN
      Row=SpectralRows_Option(1,:)
      Offset=SpectralRows_Option(3,:)
      Col=1
    END IF
    CN=0.0_DP
    KD=0.0_DP
    SD=0.0_DP
    SI=0.0_DP
    ThermC=0.0_DP
    ThermS=0.0_DP
    TS=0.0_DP
    C0Work=0.0_DP
    KappaWork=0.0_DP
    NE=M/(2*NF)
    DO Q=1,M
      CN(Offset(Q)+Row(Q),Col(Q))=-Chi(Q)
      S=(Q-1)/NE+1
      F=(S+1)/2
      SD(F,Q)=-W(Q)*MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*Chi(Q)
      SD(NF+1,Q)=-U(Q)*Chi(Q)
      DO H=1,3
        IF( .NOT. ThermalNeeded(H) ) CYCLE
        Term=DChi(Q,H)*(Target(Q)-D(Q))
        C0Work(Q,H)=Term
        KappaWork(Q,H)=DChi(Q,H)+DSigma(Q,H)
        TS(F,H)=TS(F,H)+W(Q)*MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*Term
        TS(NF+1,H)=TS(NF+1,H)+U(Q)*Term
      END DO
    END DO
    DO L=1,SIZE(A)
      IA=A(L)
      IB=B(L)
      X=D(IA)
      Y=D(IB)
      HX=1.0_DP-X
      HY=1.0_DP-Y
      G=Gain(L)*W(IA)
      L0=Loss(L)*W(IA)
      IF( Kind(L)==0 )THEN
        DA=RetainedValue(-G*Y-L0*HY)
        DB=RetainedValue(G*HX+L0*X)
        SignB=-1.0_DP
        KD(Offset(IB)+Row(IA),Col(IB))=KD(Offset(IB)+Row(IA),Col(IB))+(G-L0)/W(IA)
        KD(Offset(IA)+Row(IB),Col(IA))=KD(Offset(IA)+Row(IB),Col(IA))+(L0-G)/W(IB)
      ELSE
        DA=RetainedValue(-G*HY-L0*Y)
        DB=RetainedValue(-G*HX-L0*X)
        SignB=1.0_DP
        KD(Offset(IB)+Row(IA),Col(IB))=KD(Offset(IB)+Row(IA),Col(IB))+(L0-G)/W(IA)
        KD(Offset(IA)+Row(IB),Col(IA))=KD(Offset(IA)+Row(IB),Col(IA))+(L0-G)/W(IB)
      END IF
      CN(Offset(IA)+Row(IA),Col(IA))=CN(Offset(IA)+Row(IA),Col(IA))+DA/W(IA)
      CN(Offset(IB)+Row(IA),Col(IB))=CN(Offset(IB)+Row(IA),Col(IB))+DB/W(IA)
      CN(Offset(IA)+Row(IB),Col(IA))=CN(Offset(IA)+Row(IB),Col(IA))+SignB*DA/W(IB)
      CN(Offset(IB)+Row(IB),Col(IB))=CN(Offset(IB)+Row(IB),Col(IB))+SignB*DB/W(IB)
      DO F=1,NF
        IF( FamilyMoment(F,L)==0 ) CYCLE
        SD(F,IA)=SD(F,IA)+REAL(FamilyMoment(F,L),DP)*DA
        SD(F,IB)=SD(F,IB)+REAL(FamilyMoment(F,L),DP)*DB
      END DO
      IF( PRESENT(RealEnergyRatio_Option) )THEN
        EM=RealEnergyRatio_Option(IA)+SignB*RealEnergyRatio_Option(IB)
      ELSE
        EM=U(IA)/W(IA)+SignB*U(IB)/W(IB)
      END IF
      SD(NF+1,IA)=SD(NF+1,IA)+EM*DA
      SD(NF+1,IB)=SD(NF+1,IB)+EM*DB
      IF( PRESENT(EnergyRatio_Option) )THEN
        EMR=EnergyRatio_Option(IA)+SignB*EnergyRatio_Option(IB)
      ELSE
        EMR=Retain(U(IA))/W(IA)+SignB*(Retain(U(IB))/W(IB))
      END IF
      DO H=1,3
        IF( .NOT. ThermalNeeded(H) ) CYCLE
        GP=Retain(DGain(L,H))*W(IA)
        LP=Retain(DLoss(L,H))*W(IA)
        IF( Kind(L)==0 )THEN
          Event=GP*Y*HX-LP*X*HY
          CA=GP*Y+LP*HY
          CB=GP*HX+LP*X
        ELSE
          Event=GP*HX*HY-LP*X*Y
          CA=GP*HY+LP*Y
          CB=GP*HX+LP*X
        END IF
        C0Work(IA,H)=C0Work(IA,H)+Event/W(IA)
        C0Work(IB,H)=C0Work(IB,H)+SignB*Event/W(IB)
        KappaWork(IA,H)=KappaWork(IA,H)+CA/W(IA)
        KappaWork(IB,H)=KappaWork(IB,H)+CB/W(IB)
        DO F=1,NF
          IF( FamilyMoment(F,L)==0 ) CYCLE
          TS(F,H)=TS(F,H)+REAL(FamilyMoment(F,L),DP)*Event
        END DO
        IF( EMR%High/=0.0_DP .OR. EMR%Low/=0.0_DP ) &
          TS(NF+1,H)=TS(NF+1,H)+EMR*Event
      END DO
    END DO
    DO J=1,3
      DO Q=1,M
        IF( SIZE(SI,1)==1 )THEN
          SI(1,Q)=-U(Q)*Kappa(Q)
        ELSE
          SI(NF+1+J,(J-1)*M+Q)=-U(Q)*Kappa(Q)
        END IF
        DO L=1,M
          SD(NF+1+J,L)=SD(NF+1+J,L)-U(Q)*IDP(Q,J)*SpectralEntry(KD,Q,L,SpectralRows_Option)
        END DO
      END DO
    END DO
    DO H=1,3
      IF( .NOT. ThermalNeeded(H) ) CYCLE
      DO J=1,3
      DO Q=1,M
        C1=-KappaWork(Q,H)*I(Q,J)
        ThermC(Q,J+1,H)=RetainedValue(C1)
        TS(NF+1+J,H)=TS(NF+1+J,H)+U(Q)*C1
      END DO
      END DO
      ! Match the former post-rate target derivative and summation order.
      DO Q=1,M
        Term=Chi(Q)*DTarget(Q,H)
        C0Work(Q,H)=C0Work(Q,H)+Term
        S=(Q-1)/NE+1
        F=(S+1)/2
        TS(F,H)=TS(F,H)+MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*W(Q)*Term
        TS(NF+1,H)=TS(NF+1,H)+U(Q)*Term
      END DO
      ThermC(:,1,H)=RetainedValue(C0Work(:,H))
      ThermS(:,H)=RetainedValue(TS(:,H))
    END DO
  END SUBROUTINE ComputeCollisionRatesAndDerivatives_OrderV_Point

  PURE SUBROUTINE ComputeCollisionDensityIncrement_OrderV_Point(D,DI,I0,V,DV,DN,DeltaVI)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Shared arithmetic for native radiation and fluid family/energy moments.
    TYPE(RetainedReal), INTENT(in) :: D,DI(3),V(3),DV(3)
    REAL(DP), INTENT(in) :: I0(3)
    TYPE(RetainedReal), INTENT(out) :: DN,DeltaVI
    TYPE(RetainedReal) :: DeltaI
    INTEGER :: J
    DN=D
    DeltaVI=0.0_DP
    DO J=1,3
      DeltaI=V(J)*DI(J)+DV(J)*I0(J)
      DN=DN+DeltaI
      DeltaVI=DeltaVI+DeltaI
    END DO
  END SUBROUTINE ComputeCollisionDensityIncrement_OrderV_Point

  PURE SUBROUTINE ComputeCollisionResidual_OrderV_Point &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, C0, C1, Source0, &
      R, Physical, HScale, GScale, AbsoluteSource, NativeIncrement, &
      NativeResidual_Option, TransferResidual_Option, Exchange_Option, &
      CollisionIncrement_Option, CollisionMoment_Option, Incoming_Option, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M, NF, Active(:)
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:), C0(M), C1(M,3), Source0(NF+4)
    REAL(DP), INTENT(in) :: AbsoluteSource(NF+4)
    REAL(DP), INTENT(out) :: R(:), Physical(:), HScale(NF+4), GScale(NF+4)
    TYPE(RetainedReal), INTENT(out) :: NativeIncrement(M,4)
    TYPE(RetainedReal), INTENT(out), OPTIONAL :: NativeResidual_Option(4*M)
    TYPE(RetainedReal), INTENT(out), OPTIONAL :: TransferResidual_Option(NF+4)
    TYPE(RetainedReal), INTENT(out), OPTIONAL :: Exchange_Option(NF+4)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: CollisionIncrement_Option(M,4), CollisionMoment_Option(NF+4)
    TYPE(IncomingRadiationCache), INTENT(in), OPTIONAL :: Incoming_Option(M)
    TYPE(RetainedReal) :: T(NF+4), DV(3), V(3), DI(3), DK(3,3), K0(3,3)
    TYPE(RetainedReal) :: Moment(NF+4), Frame(NF+4), Source(NF+4), Lambda(NF+4)
    TYPE(RetainedReal) :: SharedMoment(NF+4)
    TYPE(RetainedReal) :: DeltaI, DeltaKV, DeltaVI, DN, Term, Sum, Lift, VDotL, V2
    INTEGER :: K, NA, Q, S, F, J, L, N
    REAL(DP) :: G
    K = NF+4
    NA = SIZE(Active)
    N = 4*M
    T = 0.0_DP
    DO J = 1, NA
      T(Active(J)) = X(N+J)
    END DO
    DV = -T(NF+2:K)
    V = V0+DV
    Moment = 0.0_DP
    SharedMoment = 0.0_DP
    IF( PRESENT(CollisionMoment_Option) ) SharedMoment = CollisionMoment_Option
    Frame = 0.0_DP
    HScale = 0.0_DP
    DO Q = 1, M
      DO J = 1, 3
        DI(J) = X(J*M+Q)
      END DO
      IF( PRESENT(Incoming_Option) )THEN
        K0=Incoming_Option(Q)%Pressure
        CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DI,DK,Incoming_Option(Q))
      ELSE
        CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),K0)
        CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DI,DK)
      END IF
      CALL ComputeCollisionDensityIncrement_OrderV_Point(X(Q),DI,I0(Q,:),V,DV,DN,DeltaVI)
      NativeIncrement(Q,1) = DN
      S = (Q-1)/(M/(2*NF))+1
      F = (S+1)/2
      G = MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      Moment(F) = Moment(F)+G*W(Q)*DN
      HScale(F) = HScale(F)+W(Q)*(ABS(RetainedValue(X(Q))))
      DO J = 1, 3
        HScale(F) = HScale(F)+W(Q)*(ABS(RetainedValue(V(J)*DI(J)))+ABS(RetainedValue(DV(J)*I0(Q,J))))
      END DO
      Moment(NF+1) = Moment(NF+1)+U(Q)*(X(Q)+2.0_DP*DeltaVI)
      IF( PRESENT(CollisionMoment_Option) ) SharedMoment(NF+1) = SharedMoment(NF+1)+U(Q)*DeltaVI
      HScale(NF+1) = HScale(NF+1)+U(Q)*ABS(RetainedValue(X(Q)))
      DO J = 1, 3
        HScale(NF+1) = HScale(NF+1)+2.0_DP*U(Q) &
          *(ABS(RetainedValue(V(J)*DI(J)))+ABS(RetainedValue(DV(J)*I0(Q,J))))
      END DO
      DO J = 1, 3
        DeltaKV = 0.0_DP
        DO L = 1, 3
          DeltaKV = DeltaKV+V(L)*DK(J,L)+DV(L)*K0(J,L)
        END DO
        NativeIncrement(Q,J+1) = DI(J)+DeltaKV
        Term = DI(J)+V(J)*X(Q)+DV(J)*D0(Q)+DeltaKV
        Moment(NF+1+J) = Moment(NF+1+J)+U(Q)*Term
        IF( PRESENT(CollisionMoment_Option) ) &
          SharedMoment(NF+1+J) = SharedMoment(NF+1+J)+U(Q)*(V(J)*X(Q)+DV(J)*D0(Q))
        HScale(NF+1+J) = HScale(NF+1+J)+U(Q) &
          *(ABS(RetainedValue(DI(J)))+ABS(RetainedValue(V(J)*X(Q)))+ABS(RetainedValue(DV(J)*D0(Q))))
        DO L = 1, 3
          HScale(NF+1+J) = HScale(NF+1+J)+U(Q) &
            *(ABS(RetainedValue(V(L)*DK(J,L)))+ABS(RetainedValue(DV(L)*K0(J,L))))
        END DO
        Frame(NF+1) = Frame(NF+1)+U(Q)*(DV(J)*I0(Q,J)-V(J)*DeltaKV)
        Frame(NF+1+J) = Frame(NF+1+J)+U(Q)*(DV(J)*D0(Q)-V(J)*DeltaVI)
      END DO
      IF( PRESENT(CollisionIncrement_Option) )THEN
        NativeIncrement(Q,:) = CollisionIncrement_Option(Q,:)
        DN = NativeIncrement(Q,1)
      END IF
      R(Q) = RetainedValue(DN-DT*C0(Q))
      IF( PRESENT(NativeResidual_Option) ) &
        NativeResidual_Option(Q) = DN-DT*C0(Q)
      DO J = 1, 3
        R(J*M+Q) = RetainedValue(NativeIncrement(Q,J+1)-DT*C1(Q,J))
        IF( PRESENT(NativeResidual_Option) ) &
          NativeResidual_Option(J*M+Q) = NativeIncrement(Q,J+1)-DT*C1(Q,J)
      END DO
    END DO
    IF( PRESENT(CollisionMoment_Option) ) Moment = SharedMoment
    IF( PRESENT(Exchange_Option) ) Exchange_Option = Moment
    Physical(1:N) = R(1:N)
    Source = Source0
    DO J = 1, 3
      Source(NF+1) = Source(NF+1)+V(J)*Source0(NF+1+J)
      Source(NF+1+J) = Source0(NF+1+J)+V(J)*Source0(NF+1)
    END DO
    DO J = 1, K
      R(N+J) = RetainedValue(Moment(J)-T(J))
      HScale(J) = HScale(J)+ABS(RetainedValue(T(J)))
      GScale(J) = ABS(RetainedValue(T(J)))+DT*AbsoluteSource(J)+ABS(RetainedValue(Frame(J)))
    END DO
    DO J = 1, 3
      GScale(NF+1) = GScale(NF+1)+DT*ABS(RetainedValue(V(J)))*AbsoluteSource(NF+1+J)
      GScale(NF+1+J) = GScale(NF+1+J)+DT*ABS(RetainedValue(V(J)))*AbsoluteSource(NF+1)
    END DO
    IF( PRESENT(TransferResidual_Option) ) TransferResidual_Option = 0.0_DP
    DO J = 1, NA
      L = Active(J)
      R(N+K+J) = RetainedValue(T(L)-DT*Source(L)-Frame(L))
      IF( PRESENT(TransferResidual_Option) ) &
        TransferResidual_Option(L) = T(L)-DT*Source(L)-Frame(L)
    END DO
    Physical(N+1:N+K+NA) = R(N+1:N+K+NA)
    Lambda = X(N+NA+1:N+NA+K)
    V2 = 0.0_DP
    VDotL = 0.0_DP
    DO J = 1, 3
      V2 = V2+V(J)*V(J)
      VDotL = VDotL+V(J)*Lambda(NF+1+J)
    END DO
    Sum = (Lambda(NF+1)-VDotL)/(1.0_DP-V2)
    DO J = 1, 3
      Lambda(NF+1+J) = Lambda(NF+1+J)-V(J)*Sum
    END DO
    Lambda(NF+1) = Sum
    IF( PRESENT(LiftRows_Option) )THEN
      DO Q = 1, SIZE(LiftRows_Option)
        Lift = 0.0_DP
        DO J = 1, K
          Lift = Lift+Z0(Q,J)*Lambda(J)
        END DO
        L = LiftRows_Option(Q)
        R(L) = R(L)+RetainedValue(Lift)
      END DO
    ELSE
      DO Q = 1, N
        Lift = 0.0_DP
        DO J = 1, K
          Lift = Lift+Z0(Q,J)*Lambda(J)
        END DO
        R(Q) = R(Q)+RetainedValue(Lift)
      END DO
    END IF
  END SUBROUTINE ComputeCollisionResidual_OrderV_Point

  PURE FUNCTION ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X) RESULT(V)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: V0(3)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP) :: V(3)
    INTEGER :: J,F
    V=V0
    DO J=1,SIZE(Active)
      F=Active(J)-NF-1
      IF( F>=1 .AND. F<=3 ) V(F)=RetainedValue(V0(F)-X(4*M+J))
    END DO
  END FUNCTION ComputeCollisionVelocity_OrderV_Point

  PURE SUBROUTINE ComputeCollisionColumn_OrderV_Point &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, Source0, &
      DD, DI, TC, LC, RC0, RC1, SC0, IncludeLiftDerivative, Column, Closure_Option, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M, NF, Active(:)
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source0(NF+4), DD(M), DI(M,3), TC(NF+4), LC(NF+4)
    REAL(DP), INTENT(in) :: RC0(M), RC1(M,3), SC0(NF+4)
    LOGICAL, INTENT(in) :: IncludeLiftDerivative
    REAL(DP), INTENT(out) :: Column(:)
    REAL(DP), INTENT(in), OPTIONAL :: Closure_Option(3,3,4,M)
    REAL(DP) :: T(NF+4), V(3), DV(3), VC(3), D, I(3), K0(3,3), KNow(3,3)
    REAL(DP) :: KD(3,3), KI(3,3,3), KC(3,3), KV(3), KVC(3), VI, VIC, DN
    REAL(DP) :: Moment(NF+4), Frame(NF+4), SC(NF+4), Lambda(NF+4), L0(NF+4)
    REAL(DP) :: LDot, LE, LEC, V2, Moved, G, DN0
    TYPE(RetainedReal) :: KP(3,3), DKP(3,3), DIC(3)
    INTEGER :: N, K, NA, Q, J, L, A, S, F
    LOGICAL :: MovingFrame
    N = 4*M
    K = NF+4
    NA = SIZE(Active)
    T = 0.0_DP
    DO J = 1, NA
      T(Active(J)) = RetainedValue(X(N+J))
    END DO
    DV = -T(NF+2:K)
    V = ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X)
    VC = -TC(NF+2:K)
    MovingFrame = ANY(VC /= 0.0_DP)
    Moment = 0.0_DP
    Frame = 0.0_DP
    DO Q = 1, M
      D = RetainedValue(D0(Q)+X(Q))
      DO J = 1, 3
        I(J) = RetainedValue(I0(Q,J)+X(J*M+Q))
        DIC(J) = X(J*M+Q)
      END DO
      KC = 0.0_DP
      IF( DD(Q) /= 0.0_DP .OR. ANY(DI(Q,:) /= 0.0_DP) )THEN
        IF( PRESENT(Closure_Option) )THEN
          KD = Closure_Option(:,:,1,Q)
          KI = Closure_Option(:,:,2:4,Q)
        ELSE
          CALL ComputeCollisionClosureDerivatives_OrderV_Point(D,I,KD,KI)
        END IF
        KC = KD*DD(Q)
        DO A = 1, 3
          KC = KC+KI(:,:,A)*DI(Q,A)
        END DO
      END IF
      DN = DD(Q)+SUM(V*DI(Q,:))+SUM(VC*I)
      KVC = MATMUL(KC,V)
      KV = 0.0_DP
      VI = 0.0_DP
      VIC = SUM(V*DI(Q,:))+SUM(VC*I)
      ! Retained pressure increments are needed only for a frame-changing
      ! column. Radiation/thermal columns otherwise use closure derivatives.
      IF( MovingFrame )THEN
        CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),KP)
        K0 = RetainedValue(KP)
        CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DIC,DKP)
        CALL ComputeCollisionPressure_OrderV_Point(Retain(D),Retain(I),KP)
        KNow = RetainedValue(KP)
        KV = MATMUL(RetainedValue(DKP),V)+MATMUL(K0,DV)
        KVC = KVC+MATMUL(KNow,VC)
        DO J = 1, 3
          VI = VI+V(J)*RetainedValue(DIC(J))+DV(J)*I0(Q,J)
        END DO
      END IF
      S = (Q-1)/(M/(2*NF))+1
      F = (S+1)/2
      G = MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      Moment(F) = Moment(F)+G*W(Q)*DN
      Moment(NF+1) = Moment(NF+1)+U(Q)*(DD(Q)+2.0_DP*VIC)
      Moment(NF+2:K) = Moment(NF+2:K)+U(Q)*(DI(Q,:)+V*DD(Q)+VC*D+KVC)
      Frame(NF+1) = Frame(NF+1)+U(Q)*(SUM(VC*I0(Q,:))-SUM(VC*KV)-SUM(V*KVC))
      Frame(NF+2:K) = Frame(NF+2:K)+U(Q)*(VC*D0(Q)-VC*VI-V*VIC)
      Column(Q) = DN-DT*RC0(Q)
      DO J = 1, 3
        Column(J*M+Q) = DI(Q,J)+KVC(J)-DT*RC1(Q,J)
      END DO
    END DO
    SC = SC0
    SC(NF+1) = SC0(NF+1)+SUM(V*SC0(NF+2:K))+SUM(VC*Source0(NF+2:K))
    SC(NF+2:K) = SC0(NF+2:K)+V*SC0(NF+1)+VC*Source0(NF+1)
    Column(N+1:N+K) = Moment-TC
    DO J = 1, NA
      L = Active(J)
      Column(N+K+J) = TC(L)-DT*SC(L)-Frame(L)
    END DO
    Lambda = LC
    V2 = SUM(V*V)
    IF( IncludeLiftDerivative )THEN
      L0 = RetainedValue(X(N+NA+1:N+NA+K))
      LE = (L0(NF+1)-SUM(V*L0(NF+2:K)))/(1.0_DP-V2)
      LEC = (-SUM(VC*L0(NF+2:K))+2.0_DP*SUM(V*VC)*LE)/(1.0_DP-V2)
      Lambda(NF+1) = Lambda(NF+1)+LEC
      Lambda(NF+2:K) = Lambda(NF+2:K)-VC*LE-V*LEC
    END IF
    LE = (LC(NF+1)-SUM(V*LC(NF+2:K)))/(1.0_DP-V2)
    Lambda(NF+1) = Lambda(NF+1)+LE-LC(NF+1)
    Lambda(NF+2:K) = Lambda(NF+2:K)-V*LE
    IF( PRESENT(LiftRows_Option) )THEN
      Column(LiftRows_Option) = Column(LiftRows_Option)+MATMUL(Z0,Lambda)
    ELSE
      Column(1:N) = Column(1:N)+MATMUL(Z0,Lambda)
    END IF
  END SUBROUTINE ComputeCollisionColumn_OrderV_Point

  PURE SUBROUTINE ComputeCollisionLinearColumn_OrderV_Point &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, Source, Kappa, CN, KD, SD, SI, &
      ThermC, ThermS, Field, IncludeLift, Column, Closure_Option, &
      LiftRows_Option,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    INTEGER, INTENT(in) :: M, NF, Active(:), Field
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4), Kappa(M), CN(:,:), KD(:,:)
    REAL(DP), INTENT(in) :: SD(NF+4,M), SI(:,:), ThermC(M,4,3), ThermS(NF+4,3)
    LOGICAL, INTENT(in) :: IncludeLift
    REAL(DP), INTENT(out) :: Column(:)
    REAL(DP), INTENT(in), OPTIONAL :: Closure_Option(3,3,4,M)
    REAL(DP) :: DD(M), DI(M,3), TC(NF+4), LC(NF+4), RC0(M), RC1(M,3), SC(NF+4)
    REAL(DP) :: T(NF+4), V(3), Factor
    INTEGER :: K, NA, Q, J, F, Therm
    K = NF+4
    NA = SIZE(Active)
    DD = 0.0_DP
    DI = 0.0_DP
    TC = 0.0_DP
    LC = 0.0_DP
    RC0 = 0.0_DP
    RC1 = 0.0_DP
    SC = 0.0_DP
    IF( Field <= M )THEN
      DD(Field) = 1.0_DP
      DO Q=1,M
        RC0(Q)=SpectralEntry(CN,Q,Field,SpectralRows_Option)
      END DO
      DO J = 1, 3
      DO Q = 1, M
        RC1(Q,J) = -RetainedValue(I0(Q,J)+X(J*M+Q))*SpectralEntry(KD,Q,Field,SpectralRows_Option)
      END DO
      END DO
      SC = SD(:,Field)
    ELSE IF( Field <= 4*M )THEN
      J = (Field-1)/M
      Q = MOD(Field-1,M)+1
      DI(Q,J) = 1.0_DP
      RC1(Q,J) = -Kappa(Q)
      IF( SIZE(SI,1)==1 )THEN
        SC(NF+1+J)=SI(1,Q)
      ELSE
        SC = SI(:,Field-M)
      END IF
    ELSE IF( Field <= 4*M+NA )THEN
      F = Active(Field-4*M)
      TC(F) = 1.0_DP
      Therm = 0
      Factor = 1.0_DP
      IF( F <= 2 .AND. F <= NF ) Therm = F
      IF( F == NF+1 ) Therm = 3
      IF( F > NF+1 )THEN
        T = 0.0_DP
        DO J = 1, NA
          T(Active(J)) = RetainedValue(X(4*M+J))
        END DO
        V = ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X)
        Therm = 3
        Factor = -V(F-NF-1)
      END IF
      IF( Therm > 0 )THEN
        RC0 = Factor*ThermC(:,1,Therm)
        RC1 = Factor*ThermC(:,2:4,Therm)
        SC = Factor*ThermS(:,Therm)
      END IF
    ELSE
      LC(Field-4*M-NA) = 1.0_DP
    END IF
    CALL ComputeCollisionColumn_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source, &
      DD,DI,TC,LC,RC0,RC1,SC,IncludeLift,Column,Closure_Option, &
      LiftRows_Option=LiftRows_Option)
  END SUBROUTINE ComputeCollisionLinearColumn_OrderV_Point

  PURE SUBROUTINE ComputeCollisionLinearAction_OrderV_Point &
    ( M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,CN,KD,SD,SI,ThermC,ThermS,Direction,IncludeLift,Column, &
      LiftRows_Option,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Z0(:,:),DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4),Kappa(M),CN(:,:),KD(:,:),SD(NF+4,M),SI(:,:)
    REAL(DP), INTENT(in) :: ThermC(M,4,3),ThermS(NF+4,3),Direction(:)
    LOGICAL, INTENT(in) :: IncludeLift
    REAL(DP), INTENT(out) :: Column(:)
    REAL(DP) :: DD(M),DI(M,3),TC(NF+4),LC(NF+4),RC0(M),RC1(M,3),SC(NF+4)
    REAL(DP) :: T(NF+4),V(3),Factors(3),Damping(M),AngularSource(NF+4)
    INTEGER :: J,F,NA,N,K,Q,L
    N = 4*M
    K = NF+4
    NA = SIZE(Active)
    DD = Direction(1:M)
    DO J = 1,3
      DI(:,J) = Direction(J*M+1:(J+1)*M)
    END DO
    TC = 0.0_DP
    T = 0.0_DP
    DO J = 1,NA
      TC(Active(J)) = Direction(N+J)
      T(Active(J)) = RetainedValue(X(N+J))
    END DO
    V = ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X)
    Factors = 0.0_DP
    Factors(1) = TC(1)
    IF( NF >= 2 ) Factors(2) = TC(2)
    Factors(3) = TC(NF+1)-SUM(V*TC(NF+2:K))
    LC = Direction(N+NA+1:N+NA+K)
    IF( .NOT. IncludeLift ) LC = 0.0_DP
    IF( PRESENT(SpectralRows_Option) )THEN
      RC0=0.0_DP
      Damping=0.0_DP
      DO L=1,M
      DO Q=1,M
        IF( SpectralRows_Option(2,Q)/=SpectralRows_Option(2,L) ) CYCLE
        RC0(Q)=RC0(Q)+CN(SpectralRows_Option(3,L)+SpectralRows_Option(1,Q),1)*DD(L)
        Damping(Q)=Damping(Q)+KD(SpectralRows_Option(3,L)+SpectralRows_Option(1,Q),1)*DD(L)
      END DO
      END DO
    ELSE
      RC0 = MATMUL(CN,DD)
      Damping = MATMUL(KD,DD)
    END IF
    DO J = 1,3
      RC1(:,J) = -Kappa*DI(:,J)-RetainedValue(I0(:,J)+X(J*M+1:(J+1)*M))*Damping
    END DO
    IF( SIZE(SI,1)==1 )THEN
      AngularSource=0.0_DP
      DO J=1,3
      DO Q=1,M
        AngularSource(NF+1+J)=AngularSource(NF+1+J)+SI(1,Q)*Direction(J*M+Q)
      END DO
      END DO
      SC=MATMUL(SD,DD)+AngularSource
    ELSE
      SC = MATMUL(SD,DD)+MATMUL(SI,Direction(M+1:N))
    END IF
    DO F = 1,3
      RC0 = RC0+Factors(F)*ThermC(:,1,F)
      RC1 = RC1+Factors(F)*ThermC(:,2:4,F)
      SC = SC+Factors(F)*ThermS(:,F)
    END DO
    CALL ComputeCollisionColumn_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source, &
      DD,DI,TC,LC,RC0,RC1,SC,IncludeLift,Column, &
      LiftRows_Option=LiftRows_Option)
  END SUBROUTINE ComputeCollisionLinearAction_OrderV_Point

  PURE SUBROUTINE FluxInverse( D, I, V, Damping, DT, Inverse, Info, Closure_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: D, I(3), V(3), Damping, DT
    REAL(DP), INTENT(out) :: Inverse(3,3)
    INTEGER, INTENT(out) :: Info
    REAL(DP), INTENT(in), OPTIONAL :: Closure_Option(3,3,4)
    REAL(DP) :: KD(3,3), KI(3,3,3), B(3,3), Swap(3), Factor
    INTEGER :: J, L, P
    IF( PRESENT(Closure_Option) )THEN
      KD=Closure_Option(:,:,1)
      KI=Closure_Option(:,:,2:4)
    ELSE
      CALL ComputeCollisionClosureDerivatives_OrderV_Point(D,I,KD,KI)
    END IF
    DO L = 1, 3
      B(:,L) = MATMUL(KI(:,:,L),V)
      B(L,L) = B(L,L)+1.0_DP+DT*Damping
    END DO
    Inverse = 0.0_DP
    DO J = 1, 3
      Inverse(J,J) = 1.0_DP
    END DO
    Info = 1
    DO J = 1, 3
      P = J-1+MAXLOC(ABS(B(J:3,J)),DIM=1)
      IF( ABS(B(P,J)) <= TINY(1.0_DP) ) RETURN
      Swap = B(J,:)
      B(J,:) = B(P,:)
      B(P,:) = Swap
      Swap = Inverse(J,:)
      Inverse(J,:) = Inverse(P,:)
      Inverse(P,:) = Swap
      Factor = B(J,J)
      B(J,:) = B(J,:)/Factor
      Inverse(J,:) = Inverse(J,:)/Factor
      DO L = 1, 3
        IF( L == J ) CYCLE
        Factor = B(L,J)
        B(L,:) = B(L,:)-Factor*B(J,:)
        Inverse(L,:) = Inverse(L,:)-Factor*Inverse(J,:)
      END DO
    END DO
    Info = 0
  END SUBROUTINE FluxInverse

  PURE SUBROUTINE FluxColumnAction &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, Source, Kappa, DI, Column, Closure_Option, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M, NF, Active(:)
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4), Kappa(M), DI(M,3)
    REAL(DP), INTENT(out) :: Column(:)
    REAL(DP), INTENT(in), OPTIONAL :: Closure_Option(3,3,4,M)
    REAL(DP) :: RC1(M,3), SC(NF+4), ZeroD(M), ZeroT(NF+4)
    INTEGER :: J
    ZeroD = 0.0_DP
    ZeroT = 0.0_DP
    SC = 0.0_DP
    DO J = 1, 3
      RC1(:,J) = -Kappa*DI(:,J)
      SC(NF+1+J) = SUM(U*RC1(:,J))
    END DO
    CALL ComputeCollisionColumn_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source, &
      ZeroD,DI,ZeroT,ZeroT,ZeroD,RC1,SC,.FALSE.,Column,Closure_Option, &
      LiftRows_Option=LiftRows_Option)
  END SUBROUTINE FluxColumnAction

  PURE SUBROUTINE AssembleCollisionCondensedSystem_OrderV_Point &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, Source, Kappa, CN, KD, SD, SI, &
      ThermC, ThermS, R, Matrix, RHS, BInverse, RowScale, ColumnScale, Info, FrameResidual_Option, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M, NF, Active(:)
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4), Kappa(M), CN(M,M), KD(M,M), SD(NF+4,M), SI(NF+4,3*M)
    REAL(DP), INTENT(in) :: ThermC(M,4,3), ThermS(NF+4,3), R(:)
    REAL(DP), INTENT(out) :: Matrix(:,:), RHS(:), BInverse(3,3,M), RowScale(:), ColumnScale(:)
    INTEGER, INTENT(out) :: Info
    REAL(DP), INTENT(in), OPTIONAL :: FrameResidual_Option(NF+4)
    REAL(DP) :: Column(SIZE(R)), FluxAction(SIZE(R)), DI(M,3), V(3), T(NF+4), I(3)
    REAL(DP) :: InverseG(NF+4), FrameColumn(NF+4)
    REAL(DP) :: Closure(3,3,4,M)
    INTEGER :: J, Q, K, NA, NC, N, Field, L
    K = NF+4
    NA = SIZE(Active)
    NC = M+NA+K
    N = 4*M
    T = 0.0_DP
    DO J = 1, NA
      T(Active(J)) = RetainedValue(X(N+J))
    END DO
    V = ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X)
    IF( PRESENT(FrameResidual_Option) )THEN
      InverseG=FrameResidual_Option
      InverseG(NF+1)=(FrameResidual_Option(NF+1)-SUM(V*FrameResidual_Option(NF+2:K)))/(1.0_DP-SUM(V*V))
      InverseG(NF+2:K)=FrameResidual_Option(NF+2:K)-V*InverseG(NF+1)
    END IF
    DO Q = 1, M
      DO J = 1, 3
        I(J) = RetainedValue(I0(Q,J)+X(J*M+Q))
      END DO
      CALL FluxInverse(RetainedValue(D0(Q)+X(Q)),I,V,Kappa(Q),DT,BInverse(:,:,Q),Info)
      IF( Info /= 0 ) RETURN
      ! These derivatives depend on this iterate, not on the matrix column.
      ! Keep the original tensor contractions and arithmetic order below.
      CALL ComputeCollisionClosureDerivatives_OrderV_Point(RetainedValue(D0(Q)+X(Q)),I, &
        Closure(:,:,1,Q),Closure(:,:,2:4,Q))
    END DO
    DO J = 1, NC
      Field = J
      IF( J > M ) Field = J+3*M
      CALL ComputeCollisionLinearColumn_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,CN,KD,SD,SI, &
        ThermC,ThermS,Field,.TRUE.,Column,Closure, &
      LiftRows_Option=LiftRows_Option)
      IF( PRESENT(FrameResidual_Option) .AND. J>M .AND. J<=M+NA )THEN
        L=Active(J-M)-NF-1
        IF( L>=1 .AND. L<=3 )THEN
          ! -d(T_v)/dq_P T_v^{-1} g in the transfer border. This term
          ! is essential away from the solution of the independent system.
          FrameColumn=0.0_DP
          FrameColumn(NF+1)=InverseG(NF+1+L)
          FrameColumn(NF+1+L)=InverseG(NF+1)
          Column(N+K+1:N+K+NA)=Column(N+K+1:N+K+NA)+FrameColumn(Active)
        END IF
      END IF
      DO Q = 1, M
        DI(Q,:) = -MATMUL(BInverse(:,:,Q),Column([M+Q,2*M+Q,3*M+Q]))
      END DO
      CALL FluxColumnAction(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,DI,FluxAction,Closure, &
      LiftRows_Option=LiftRows_Option)
      Matrix(1:M,J) = Column(1:M)+FluxAction(1:M)
      Matrix(M+1:NC,J) = Column(N+1:)+FluxAction(N+1:)
    END DO
    DO Q = 1, M
      DI(Q,:) = -MATMUL(BInverse(:,:,Q),R([M+Q,2*M+Q,3*M+Q]))
    END DO
    CALL FluxColumnAction(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,DI,FluxAction,Closure, &
      LiftRows_Option=LiftRows_Option)
    RHS(1:M) = -R(1:M)-FluxAction(1:M)
    RHS(M+1:NC) = -R(N+1:)-FluxAction(N+1:)
    RowScale = MAX(MAXVAL(ABS(Matrix),DIM=2),TINY(1.0_DP))
    DO J = 1, NC
      Matrix(:,J) = Matrix(:,J)/RowScale
    END DO
    RHS = RHS/RowScale
    ColumnScale = MAX(MAXVAL(ABS(Matrix),DIM=1),TINY(1.0_DP))
    DO J = 1, NC
      Matrix(:,J) = Matrix(:,J)/ColumnScale(J)
    END DO
    Info = 0
  END SUBROUTINE AssembleCollisionCondensedSystem_OrderV_Point

  PURE SUBROUTINE AssembleCollisionCondensedSystem &
    (M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,CN,KD,SD,SI,ThermC,ThermS, &
     R,Matrix,RHS,BInverse,RowScale,ColumnScale,Info,FrameResidual_Option,LiftRows_Option,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Z0(:,:),DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4),Kappa(M),CN(:,:),KD(:,:),SD(NF+4,M),SI(:,:)
    REAL(DP), INTENT(in) :: ThermC(M,4,3),ThermS(NF+4,3),R(:)
    REAL(DP), INTENT(out) :: Matrix(:,:),RHS(:),BInverse(3,3,M),RowScale(:),ColumnScale(:)
    INTEGER, INTENT(out) :: Info
    REAL(DP), INTENT(in), OPTIONAL :: FrameResidual_Option(NF+4)
    ! Row-wise construction of precisely the same angular Schur complement.
    ! Only one bin's closure tensor is live, instead of 36*M private doubles.
    REAL(DP) :: FM(NF+4,SIZE(Matrix,2)),FF(4,SIZE(Matrix,2)),FS(3,SIZE(Matrix,2))
    REAL(DP) :: Lambda(NF+4,SIZE(Active)+NF+4),VC(3,SIZE(Active)+NF+4),VCol(3)
    REAL(DP) :: Closure(3,3,4),KC(3,3),K0(3,3),KN(3,3),KV(3),KVC(3),FKVC(3)
    REAL(DP) :: V(3),DV(3),T(NF+4),SC(NF+4),TC(NF+4),LC(NF+4),L0(NF+4)
    REAL(DP) :: D,I(3),DI(3),RC1(3),DD,RC0,DN,VIC,VI,FDN,Factor,LE,LEC,V2,G
    REAL(DP) :: Lift(4),InverseG(NF+4),FrameColumn(NF+4),FluxAction(SIZE(R)),RDI(M,3)
    TYPE(RetainedReal) :: KP(3,3),DKP(3,3),DIC(3)
    INTEGER :: N,K,NA,NC,Q,J,L,A,F,S,H,Slots(4)
    LOGICAL :: NeedFrame,UseLift(SIZE(Active)+NF+4)
    N=4*M
    K=NF+4
    NA=SIZE(Active)
    NC=M+NA+K
    Matrix=0.0_DP
    FM=0.0_DP
    FF=0.0_DP
    FS=0.0_DP
    T=0.0_DP
    DO J=1,NA
      T(Active(J))=RetainedValue(X(N+J))
    END DO
    DV=-T(NF+2:K)
    V=ComputeCollisionVelocity_OrderV_Point(M,NF,Active,V0,X)
    V2=SUM(V*V)
    L0=RetainedValue(X(N+NA+1:N+NA+K))
    VC=0.0_DP
    Lambda=0.0_DP
    DO J=1,NA+K
      IF( J<=NA )THEN
        L=Active(J)-NF-1
        IF( L>=1 .AND. L<=3 ) VC(L,J)=-1.0_DP
      END IF
      LC=0.0_DP
      IF( J>NA ) LC(J-NA)=1.0_DP
      Lambda(:,J)=LC
      LE=(L0(NF+1)-SUM(V*L0(NF+2:K)))/(1.0_DP-V2)
      LEC=(-SUM(VC(:,J)*L0(NF+2:K))+2.0_DP*SUM(V*VC(:,J))*LE)/(1.0_DP-V2)
      Lambda(NF+1,J)=Lambda(NF+1,J)+LEC
      Lambda(NF+2:K,J)=Lambda(NF+2:K,J)-VC(:,J)*LE-V*LEC
      LE=(LC(NF+1)-SUM(V*LC(NF+2:K)))/(1.0_DP-V2)
      Lambda(NF+1,J)=Lambda(NF+1,J)+LE-LC(NF+1)
      Lambda(NF+2:K,J)=Lambda(NF+2:K,J)-V*LE
      UseLift(J)=ANY(Lambda(:,J)/=0.0_DP)
    END DO
    NeedFrame=ANY(VC/=0.0_DP)
    DO Q=1,M
      D=RetainedValue(D0(Q)+X(Q))
      DO A=1,3
        DIC(A)=X(A*M+Q)
        I(A)=RetainedValue(I0(Q,A)+DIC(A))
      END DO
      CALL ComputeCollisionClosureDerivatives_OrderV_Point(D,I,Closure(:,:,1),Closure(:,:,2:4))
      CALL FluxInverse(D,I,V,Kappa(Q),DT,BInverse(:,:,Q),Info,Closure)
      IF( Info/=0 ) RETURN
      KV=0.0_DP
      VI=0.0_DP
      IF( NeedFrame )THEN
        CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),KP)
        K0=RetainedValue(KP)
        CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DIC,DKP)
        CALL ComputeCollisionPressure_OrderV_Point(Retain(D),Retain(I),KP)
        KN=RetainedValue(KP)
        KV=MATMUL(RetainedValue(DKP),V)+MATMUL(K0,DV)
        DO A=1,3
          VI=VI+V(A)*RetainedValue(DIC(A))+DV(A)*I0(Q,A)
        END DO
      END IF
      Slots=0
      IF( PRESENT(LiftRows_Option) )THEN
        DO L=1,SIZE(LiftRows_Option)
          DO A=1,4
            IF( LiftRows_Option(L)==(A-1)*M+Q ) Slots(A)=L
          END DO
        END DO
      END IF
      S=(Q-1)/(M/(2*NF))+1
      F=(S+1)/2
      G=MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      DO J=1,NC
        IF( PRESENT(SpectralRows_Option) )THEN
          IF( J<=M )THEN
            IF( SpectralRows_Option(2,Q)/=SpectralRows_Option(2,J) ) CYCLE
          END IF
        END IF
        VCol=0.0_DP
        IF( J>M ) VCol=VC(:,J-M)
        DD=0.0_DP
        RC0=0.0_DP
        RC1=0.0_DP
        IF( J<=M )THEN
          IF( Q==J ) DD=1.0_DP
          RC0=SpectralEntry(CN,Q,J,SpectralRows_Option)
          RC1=-I*SpectralEntry(KD,Q,J,SpectralRows_Option)
        ELSE IF( J<=M+NA )THEN
          L=Active(J-M)
          H=0
          Factor=1.0_DP
          IF( L<=MIN(2,NF) ) H=L
          IF( L==NF+1 ) H=3
          IF( L>NF+1 )THEN
            H=3
            Factor=-V(L-NF-1)
          END IF
          IF( H>0 )THEN
            RC0=Factor*ThermC(Q,1,H)
            RC1=Factor*ThermC(Q,2:4,H)
          END IF
        END IF
        KC=0.0_DP
        IF( DD/=0.0_DP )THEN
          KC=Closure(:,:,1)*DD
          DO A=1,3
            KC=KC+Closure(:,:,A+1)*0.0_DP
          END DO
        END IF
        KVC=MATMUL(KC,V)
        IF( ANY(VCol/=0.0_DP) ) KVC=KVC+MATMUL(KN,VCol)
        VIC=SUM(VCol*I)
        DN=DD+VIC
        Lift=0.0_DP
        IF( J>M )THEN
        IF( UseLift(J-M) )THEN
          DO A=1,4
            IF( PRESENT(LiftRows_Option) )THEN
              IF( Slots(A)==0 ) CYCLE
              Lift(A)=SUM(Z0(Slots(A),:)*Lambda(:,J-M))
            ELSE
              Lift(A)=SUM(Z0((A-1)*M+Q,:)*Lambda(:,J-M))
            END IF
          END DO
        END IF
        END IF
        DI=-MATMUL(BInverse(:,:,Q),KVC-DT*RC1+Lift(2:4))
        Matrix(Q,J)=DN-DT*RC0+Lift(1)+SUM(V*DI)
        Matrix(M+F,J)=Matrix(M+F,J)+G*W(Q)*DN
        Matrix(M+NF+1,J)=Matrix(M+NF+1,J)+U(Q)*(DD+2.0_DP*VIC)
        Matrix(M+NF+2:M+K,J)=Matrix(M+NF+2:M+K,J)+U(Q)*(V*DD+VCol*D+KVC)
        FrameColumn=0.0_DP
        FrameColumn(NF+1)=U(Q)*(SUM(VCol*I0(Q,:))-SUM(VCol*KV)-SUM(V*KVC))
        FrameColumn(NF+2:K)=U(Q)*(VCol*D0(Q)-VCol*VI-V*VIC)
        Matrix(M+K+1:NC,J)=Matrix(M+K+1:NC,J)+FrameColumn(Active)
        KC=0.0_DP
        IF( ANY(DI/=0.0_DP) )THEN
          KC=Closure(:,:,1)*0.0_DP
          DO A=1,3
            KC=KC+Closure(:,:,A+1)*DI(A)
          END DO
        END IF
        FKVC=MATMUL(KC,V)
        FDN=SUM(V*DI)
        FM(F,J)=FM(F,J)+G*W(Q)*FDN
        FM(NF+1,J)=FM(NF+1,J)+U(Q)*(2.0_DP*FDN)
        FM(NF+2:K,J)=FM(NF+2:K,J)+U(Q)*(DI+FKVC)
        FF(1,J)=FF(1,J)+U(Q)*(-SUM(V*FKVC))
        FF(2:4,J)=FF(2:4,J)+U(Q)*(-V*FDN)
        FS(:,J)=FS(:,J)+U(Q)*(-Kappa(Q)*DI)
      END DO
    END DO
    IF( PRESENT(FrameResidual_Option) )THEN
      InverseG=FrameResidual_Option
      InverseG(NF+1)=(FrameResidual_Option(NF+1)-SUM(V*FrameResidual_Option(NF+2:K)))/(1.0_DP-SUM(V*V))
      InverseG(NF+2:K)=FrameResidual_Option(NF+2:K)-V*InverseG(NF+1)
    END IF
    DO J=1,NC
      VCol=0.0_DP
      IF( J>M ) VCol=VC(:,J-M)
      TC=0.0_DP
      SC=0.0_DP
      IF( J<=M )THEN
        SC=SD(:,J)
      ELSE IF( J<=M+NA )THEN
        L=Active(J-M)
        TC(L)=1.0_DP
        H=0
        Factor=1.0_DP
        IF( L<=MIN(2,NF) ) H=L
        IF( L==NF+1 ) H=3
        IF( L>NF+1 )THEN
          H=3
          Factor=-V(L-NF-1)
        END IF
        IF( H>0 ) SC=Factor*ThermS(:,H)
      END IF
      FrameColumn=SC
      FrameColumn(NF+1)=SC(NF+1)+SUM(V*SC(NF+2:K))+SUM(VCol*Source(NF+2:K))
      FrameColumn(NF+2:K)=SC(NF+2:K)+V*SC(NF+1)+VCol*Source(NF+1)
      Matrix(M+1:M+K,J)=Matrix(M+1:M+K,J)-TC+FM(:,J)
      Matrix(M+K+1:NC,J)=TC(Active)-DT*FrameColumn(Active)-Matrix(M+K+1:NC,J)
      IF( PRESENT(FrameResidual_Option) .AND. J>M .AND. J<=M+NA )THEN
        L=Active(J-M)-NF-1
        IF( L>=1 .AND. L<=3 )THEN
          FrameColumn=0.0_DP
          FrameColumn(NF+1)=InverseG(NF+1+L)
          FrameColumn(NF+1+L)=InverseG(NF+1)
          Matrix(M+K+1:NC,J)=Matrix(M+K+1:NC,J)+FrameColumn(Active)
        END IF
      END IF
      FrameColumn=0.0_DP
      FrameColumn(NF+1)=SUM(V*FS(:,J))
      FrameColumn(NF+2:K)=FS(:,J)
      SC=0.0_DP
      SC(NF+1:K)=FF(:,J)
      Matrix(M+K+1:NC,J)=Matrix(M+K+1:NC,J)+(-DT*FrameColumn(Active)-SC(Active))
    END DO
    DO Q=1,M
      RDI(Q,:)=-MATMUL(BInverse(:,:,Q),R([M+Q,2*M+Q,3*M+Q]))
    END DO
    CALL FluxColumnAction(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,RDI,FluxAction, &
      LiftRows_Option=LiftRows_Option)
    RHS(1:M)=-R(1:M)-FluxAction(1:M)
    RHS(M+1:NC)=-R(N+1:)-FluxAction(N+1:)
    RowScale=MAX(MAXVAL(ABS(Matrix),DIM=2),TINY(1.0_DP))
    DO J=1,NC
      Matrix(:,J)=Matrix(:,J)/RowScale
    END DO
    RHS=RHS/RowScale
    ColumnScale=MAX(MAXVAL(ABS(Matrix),DIM=1),TINY(1.0_DP))
    DO J=1,NC
      Matrix(:,J)=Matrix(:,J)/ColumnScale(J)
    END DO
    Info=0
  END SUBROUTINE AssembleCollisionCondensedSystem

  PURE SUBROUTINE ExpandCollisionCondensedStep_OrderV_Point &
    ( M, NF, Active, W, U, D0, I0, V0, X, Z0, DT, Source, Kappa, CN, KD, SD, SI, &
      ThermC, ThermS, R, BInverse, Condensed, Step, PhysicalAction, &
      LiftRows_Option,SpectralRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in), OPTIONAL :: SpectralRows_Option(:,:)
    INTEGER, INTENT(in) :: M, NF, Active(:)
    REAL(DP), INTENT(in) :: W(M), U(M), D0(M), I0(M,3), V0(3), Z0(:,:), DT
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Source(NF+4), Kappa(M), CN(:,:), KD(:,:), SD(NF+4,M), SI(:,:)
    REAL(DP), INTENT(in) :: ThermC(M,4,3), ThermS(NF+4,3), R(:), BInverse(3,3,M), Condensed(:)
    REAL(DP), INTENT(out) :: Step(:), PhysicalAction(:)
    REAL(DP) :: Column(SIZE(R)), Current(SIZE(R)), DI(M,3), Action(SIZE(R))
    INTEGER :: J, Q, N, Field
    N = 4*M
    Step(1:M) = Condensed(1:M)
    Step(N+1:) = Condensed(M+1:)
    Step(M+1:N) = 0.0_DP
    Current = 0.0_DP
    PhysicalAction = 0.0_DP
    CALL ComputeCollisionLinearAction_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,CN,KD,SD,SI, &
      ThermC,ThermS,Step,.TRUE.,Current, &
      LiftRows_Option=LiftRows_Option,SpectralRows_Option=SpectralRows_Option)
    DO Q = 1, M
      DI(Q,:) = -MATMUL(BInverse(:,:,Q),R([M+Q,2*M+Q,3*M+Q])+Current([M+Q,2*M+Q,3*M+Q]))
    END DO
    DO J = 1, 3
      Step(J*M+1:(J+1)*M) = DI(:,J)
    END DO
    CALL ComputeCollisionLinearAction_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z0,DT,Source,Kappa,CN,KD,SD,SI, &
      ThermC,ThermS,Step,.FALSE.,PhysicalAction, &
      LiftRows_Option=LiftRows_Option,SpectralRows_Option=SpectralRows_Option)
  END SUBROUTINE ExpandCollisionCondensedStep_OrderV_Point

  PURE SUBROUTINE ComputeCollisionMetrics_OrderV_Point &
    ( NE, NS, Active, W, U, D0, I0, X, DT, Physical, Kappa, DensityLoss, HScale, GScale, &
      InitialNorm, SetInitial, RadiationError, ConservationError, Weights, Merit )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: NE, NS, Active(:)
    REAL(DP), INTENT(in) :: W(NE*NS), U(NE*NS), D0(NE*NS), I0(NE*NS,3), DT, Physical(:)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Kappa(NE*NS), DensityLoss(NE*NS), HScale(NS/2+4), GScale(NS/2+4)
    REAL(DP), INTENT(inout) :: InitialNorm(NS)
    LOGICAL, INTENT(in) :: SetInitial
    REAL(DP), INTENT(out) :: RadiationError, ConservationError, Weights(:), Merit
    REAL(DP) :: Norm, Scale, Population, Maximum, SumD, SumR, Denom, Tau, MomentumScale
    REAL(DP) :: Residual(NE), Density(NE), Norms(4), NormD, Relative, DenomH
    REAL(DP) :: Totals(NS/2+4), Expected(NS/2+4)
    INTEGER :: M, N, K, Q, S, J, L, F
    M = NE*NS
    N = 4*M
    K = NS/2+4
    RadiationError = 0.0_DP
    Weights = 0.0_DP
    DO S = 1, NS
      DO J = 0, 3
        DO Q = 1, NE
          L = (S-1)*NE+Q
          Denom = 1.0_DP+DT*Kappa(L)
          IF( J == 0 ) Denom = 1.0_DP+DT*DensityLoss(L)
          Residual(Q) = Physical(J*M+L)/Denom
          Density(Q) = RetainedValue(D0(L)+X(L))
        END DO
        Maximum = MAXVAL(ABS(Residual))
        Norms(J+1) = 0.0_DP
        IF( Maximum > 0.0_DP ) Norms(J+1) = Maximum*SQRT(SUM(W(1:NE)*(Residual/Maximum)**2))
      END DO
      Norm = MAXVAL(Norms)
      IF( SetInitial ) InitialNorm(S) = Norm
      Maximum = MAXVAL(ABS(Density))
      NormD = 0.0_DP
      IF( Maximum > 0.0_DP ) NormD = Maximum*SQRT(SUM(W(1:NE)*(Density/Maximum)**2))
      Scale = MAX(NormD,0.01_DP*InitialNorm(S))
      IF( Scale > 0.0_DP ) RadiationError = MAX(RadiationError,Norm/Scale)
      DO J = 0, 3
      DO Q = 1, NE
        L = (S-1)*NE+Q
        Denom = 1.0_DP+DT*Kappa(L)
        IF( J == 0 ) Denom = 1.0_DP+DT*DensityLoss(L)
        Weights(J*M+L) = SQRT(W(Q))/(MAX(Scale,TINY(1.0_DP))*Denom)
      END DO
      END DO
    END DO
    ConservationError = 0.0_DP
    Totals = 0.0_DP
    Expected = 0.0_DP
    DO L = 1, M
      S = (L-1)/NE+1
      F = (S+1)/2
      Totals(F) = Totals(F)+W(L)*ABS(D0(L))
      Expected(F) = Expected(F)+W(L)*ABS(Physical(L))
      Totals(NS/2+1) = Totals(NS/2+1)+ABS(U(L)*D0(L))
      Expected(NS/2+1) = Expected(NS/2+1)+ABS(U(L)*Physical(L))
      DO J = 1, 3
        Totals(NS/2+1+J) = Totals(NS/2+1+J)+ABS(U(L)*D0(L))
        Expected(NS/2+1+J) = Expected(NS/2+1+J)+ABS(U(L)*Physical(J*M+L))
      END DO
    END DO
    DO F = 1, K
      ConservationError = MAX(ConservationError,ABS(Physical(N+F)) &
        /MAX(128.0_DP*EPSILON(1.0_DP)*HScale(F),1.0e-300_DP))
      ! Same increment-based scales as the residual. Family transfer rows
      ! are never normalized to the full electron/radiation population.
      DenomH = MAX(Totals(F)+HScale(F)+Expected(F),1.0e-300_DP)
      Weights(N+F) = 1.0_DP/DenomH
    END DO
    ! All orthonormal momentum components have the same units. Normalize
    ! them as a vector: a zero transverse component is not an infinitely
    ! accurate independent target. Its derivative need not be zero.
    Maximum = MAXVAL(GScale(NS/2+2:K))
    MomentumScale = Totals(NS/2+1)
    IF( Maximum > 0.0_DP ) MomentumScale = Maximum*SQRT(SUM((GScale(NS/2+2:K)/Maximum)**2))
    DO J = 1, SIZE(Active)
      F = Active(J)
      Denom = GScale(F)
      IF( F > NS/2+1 ) Denom = MomentumScale
      IF( Denom <= 0.0_DP ) Denom = Totals(F)
      Weights(N+K+J) = 1.0_DP/MAX(Denom,1.0e-300_DP)
    END DO
    Merit = 0.5_DP*SUM((Weights*Physical)**2)
  END SUBROUTINE ComputeCollisionMetrics_OrderV_Point

  PURE SUBROUTINE RefineCollisionConservation_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,Cap,X,Corrections)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Cap
    TYPE(RetainedReal), INTENT(inout) :: X(:)
    INTEGER, INTENT(out) :: Corrections
    REAL(DP) :: Support(4*M,NF+4),Matrix(NF+4,NF+4),RHS(NF+4),Swap(NF+4)
    REAL(DP) :: D(M),I(M,3),Column(SIZE(X)),ZeroD(M),ZeroI(M,3),ZeroS(NF+4)
    REAL(DP) :: R(SIZE(X)),Physical(SIZE(X)),HScale(NF+4),GScale(NF+4)
    REAL(DP) :: Previous,Current,Pivot,Factor,Row(NF+4)
    TYPE(RetainedReal) :: Trial(SIZE(X)),C0(M),C1(M,3),Source(NF+4),Increment(M,4),Term
    INTEGER :: N,K,Pass,J,L,P,Q,Status
    N=4*M
    K=NF+4
    Corrections=0
    ZeroD=0.0_DP
    ZeroI=0.0_DP
    ZeroS=0.0_DP
    C0=0.0_DP
    C1=0.0_DP
    Source=0.0_DP
    Support=0.0_DP
    ! Refine only the moment constraints at FIXED fluid transfer. No EOS,
    ! opacity, collision-source or fluid-state modification is performed.
    CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Support,0.0_DP,C0,C1,Source, &
      R,Physical,HScale,GScale,ZeroS,Increment)
    Previous=MAXVAL(ABS(Physical(N+1:N+K))/MAX(128.0_DP*EPSILON(1.0_DP)*HScale,1.0e-300_DP))
    IF( Previous<1.0_DP ) RETURN
    DO Pass=1,4
      D=RetainedValue(D0+X(1:M))
      DO J=1,3
        I(:,J)=RetainedValue(I0(:,J)+X(J*M+1:(J+1)*M))
      END DO
      CALL SelectCollisionMomentSupport_OrderV_Point(M,NF,W,U,D,I,Cap,Support,Status)
      IF( Status/=0 ) RETURN
      DO J=1,K
        CALL ComputeCollisionColumn_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Support,0.0_DP,ZeroS, &
          Support(1:M,J),RESHAPE(Support(M+1:N,J),[M,3]),ZeroS,ZeroS,ZeroD,ZeroI,ZeroS,.FALSE.,Column)
        Matrix(:,J)=Column(N+1:N+K)
      END DO
      RHS=Physical(N+1:N+K)
      Row=MAX(MAXVAL(ABS(Matrix),DIM=2),TINY(1.0_DP))
      DO J=1,K
        Matrix(J,:)=Matrix(J,:)/Row(J)
      END DO
      RHS=RHS/Row
      ! At most seven rows, as in the device-local support construction.
      DO J=1,K
        P=J-1+MAXLOC(ABS(Matrix(J:K,J)),DIM=1)
        IF( ABS(Matrix(P,J))<=128.0_DP*EPSILON(1.0_DP) ) RETURN
        Swap=Matrix(J,:)
        Matrix(J,:)=Matrix(P,:)
        Matrix(P,:)=Swap
        Pivot=RHS(J)
        RHS(J)=RHS(P)
        RHS(P)=Pivot
        Pivot=Matrix(J,J)
        Matrix(J,:)=Matrix(J,:)/Pivot
        RHS(J)=RHS(J)/Pivot
        DO L=1,K
          IF( L==J ) CYCLE
          Factor=Matrix(L,J)
          Matrix(L,:)=Matrix(L,:)-Factor*Matrix(J,:)
          RHS(L)=RHS(L)-Factor*RHS(J)
        END DO
      END DO
      Trial=X
      DO Q=1,N
        Term=0.0_DP
        DO J=1,K
          Term=Term+Retain(Support(Q,J))*RHS(J)
        END DO
        Trial(Q)=Trial(Q)-Term
      END DO
      IF( .NOT. CollisionStateAdmissible_OrderV_Point(M,D0,I0,Trial,Cap) ) RETURN
      CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,Trial,Support,0.0_DP,C0,C1,Source, &
        R,Physical,HScale,GScale,ZeroS,Increment)
      Current=MAXVAL(ABS(Physical(N+1:N+K))/MAX(128.0_DP*EPSILON(1.0_DP)*HScale,1.0e-300_DP))
      IF( .NOT.(Current<Previous) ) RETURN
      X=Trial
      Corrections=Corrections+1
      IF( Current<1.0_DP ) RETURN
      Previous=Current
    END DO
  END SUBROUTINE RefineCollisionConservation_OrderV_Point

  PURE FUNCTION ComputeCollisionStepBound_OrderV_Point( M, D0, I0, X, Step, Cap, RetainedStep_Option ) RESULT(Alpha)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M
    REAL(DP), INTENT(in) :: D0(M), I0(M,3), Step(:), Cap
    TYPE(RetainedReal), INTENT(in) :: X(:)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: RetainedStep_Option(:)
    REAL(DP) :: Alpha, Left, Right, Mid
    TYPE(RetainedReal) :: Trial(SIZE(X))
    INTEGER :: J
    Trial = X+Step
    IF( PRESENT(RetainedStep_Option) ) Trial=X+RetainedStep_Option
    Alpha = 1.0_DP
    IF( CollisionStateAdmissible_OrderV_Point(M,D0,I0,Trial,Cap) ) RETURN
    ! Convex line/cone intersection; bisection also covers all vector flux
    ! directions and avoids cancellation in a nearly tangent quadratic root.
    Left = 0.0_DP
    Right = 1.0_DP
    DO J = 1, 54
      Mid = 0.5_DP*(Left+Right)
      Trial = X+Mid*Step
      IF( PRESENT(RetainedStep_Option) ) Trial=X+Mid*RetainedStep_Option
      IF( CollisionStateAdmissible_OrderV_Point(M,D0,I0,Trial,Cap) )THEN
        Left = Mid
      ELSE
        Right = Mid
      END IF
    END DO
    Alpha = 0.99_DP*Left
  END FUNCTION ComputeCollisionStepBound_OrderV_Point

  PURE SUBROUTINE SelectCollisionMomentSupport_OrderV_Point( M, NF, W, U, D, I, Cap, Z, Info )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: M, NF
    REAL(DP), INTENT(in) :: W(M), U(M), D(M), I(M,3), Cap
    REAL(DP), INTENT(out) :: Z(4*M,NF+4)
    INTEGER, INTENT(out) :: Info
    REAL(DP) :: Basis(NF+4,4*M), Residual(NF+4,4*M), Q(NF+4,NF+4)
    REAL(DP) :: Selected(NF+4,NF+4), RHS(NF+4,NF+4), Row(NF+4), Capacity, Norm, Best
    REAL(DP) :: Swap(NF+4), Pivot, Factor
    INTEGER :: K, N, J, L, S, F, Column, Permutation(NF+4), P
    ! Setup only. Row-balanced pivoted modified Gram-Schmidt chooses the
    ! same well-scaled kind of support as the reference's pivoted QR.
    K = NF+4
    N = 4*M
    Info = 1
    IF( ANY(D <= 0.0_DP) ) RETURN
    Basis = 0.0_DP
    DO J = 1, M
      S = (J-1)/(M/(2*NF))+1
      F = (S+1)/2
      Basis(F,J) = MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*W(J)
      Basis(NF+1,J) = U(J)
      DO L = 1, 3
        Basis(NF+1+L,L*M+J) = U(J)
      END DO
    END DO
    Residual = Basis
    DO J = 1, M
      Capacity = MIN(D(J),Cap-D(J))
      Residual(:,J) = Residual(:,J)*MAX(0.0_DP,Capacity)
      Capacity = D(J)*MAX(0.0_DP,1.0_DP-SQRT(SUM((I(J,:)/D(J))**2)))
      DO L = 1, 3
        Residual(:,L*M+J) = Residual(:,L*M+J)*Capacity
      END DO
    END DO
    Row = MAXVAL(ABS(Residual),DIM=2)
    IF( ANY(Row <= TINY(1.0_DP)) ) RETURN
    DO J = 1, N
      Residual(:,J) = Residual(:,J)/Row
    END DO
    DO L = 1, K
      Best = -1.0_DP
      Column = 1
      DO J = 1, N
        Norm = SUM(Residual(:,J)**2)
        IF( Norm > Best )THEN
          Best = Norm
          Column = J
        END IF
      END DO
      IF( Best <= 128.0_DP*EPSILON(1.0_DP)**2 ) RETURN
      Permutation(L) = Column
      Q(:,L) = Residual(:,Column)/SQRT(Best)
      DO J = 1, N
        Residual(:,J) = Residual(:,J)-Q(:,L)*SUM(Q(:,L)*Residual(:,J))
      END DO
    END DO
    DO L = 1, K
      Selected(:,L) = Basis(:,Permutation(L))/Row
    END DO
    RHS = 0.0_DP
    DO L = 1, K
      RHS(L,L) = 1.0_DP/Row(L)
    END DO
    ! At most seven rows; this setup kernel runs on the device. The large
    ! spectral systems use the vendor batched factor/apply interfaces.
    DO J = 1, K
      P = J-1+MAXLOC(ABS(Selected(J:K,J)),DIM=1)
      IF( ABS(Selected(P,J)) <= 128.0_DP*EPSILON(1.0_DP) ) RETURN
      Swap = Selected(J,:)
      Selected(J,:) = Selected(P,:)
      Selected(P,:) = Swap
      Swap = RHS(J,:)
      RHS(J,:) = RHS(P,:)
      RHS(P,:) = Swap
      Pivot = Selected(J,J)
      Selected(J,:) = Selected(J,:)/Pivot
      RHS(J,:) = RHS(J,:)/Pivot
      DO L = 1, K
        IF( L == J ) CYCLE
        Factor = Selected(L,J)
        Selected(L,:) = Selected(L,:)-Factor*Selected(J,:)
        RHS(L,:) = RHS(L,:)-Factor*RHS(J,:)
      END DO
    END DO
    Info = 0
    Z = 0.0_DP
    IF( Info /= 0 ) RETURN
    DO L = 1, K
      Z(Permutation(L),:) = RHS(L,:)
    END DO
  END SUBROUTINE SelectCollisionMomentSupport_OrderV_Point

END MODULE TwoMoment_CollisionOperatorModule_OrderV
