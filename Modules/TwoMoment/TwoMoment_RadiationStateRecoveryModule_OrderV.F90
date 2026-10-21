MODULE TwoMoment_RadiationStateRecoveryModule_OrderV

  USE KindModule, ONLY: DP
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE TwoMoment_RetainedArithmeticModule
  USE TwoMoment_CollisionOperatorModule_OrderV, ONLY: ComputeCollisionResidual_OrderV_Point, &
    ComputeCollisionClosureDerivatives_OrderV_Point, CollisionStateAdmissible_OrderV_Point, ComputeCollisionMetrics_OrderV_Point, &
    ComputeCollisionPressure_OrderV_Point, ComputeCollisionPressureIncrement_OrderV_Point, ComputeCollisionDensityIncrement_OrderV_Point, IncomingRadiationCache

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ComputeIndependentRadiationResidual, ComputeIndependentRadiationColumn
  PUBLIC :: ComputeConservationGeometry, ReconstructConservativeRadiation
  PUBLIC :: ProjectRadiationCone, ProjectConservativeDirection
  PUBLIC :: SolveSmallSymmetric, SelectConservationPivots
  PUBLIC :: ComputeIndependentMetrics
  PUBLIC :: ComputeFluidSourceMetrics
  PUBLIC :: ComputeCoupledSourceMetrics, ComputeCoupledMerit, ComputeCoupledSlope
  PUBLIC :: PrepareCoupledSearchWeights
  PUBLIC :: SelectConservativeSupport
  PUBLIC :: RecoverRadiationExchange, ReconstructRadiationInvariants
  PUBLIC :: RecoverSharedProcessState
  PUBLIC :: ComputeSharedProcessTangent
  PUBLIC :: ProjectSharedProcessEndpoint
  PUBLIC :: BuildScatteringComponents

CONTAINS

  PURE SUBROUTINE SolveRetainedSmall(N,A,B,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: N
    TYPE(RetainedReal), INTENT(inout) :: A(N,N), B(:,:)
    INTEGER, INTENT(out) :: Info
    TYPE(RetainedReal) :: Row(N), RHS(SIZE(B,2)), Factor
    INTEGER :: J,L,P
    Info=1
    DO J=1,N
      P=J-1+MAXLOC(ABS(RetainedValue(A(J:N,J))),DIM=1)
      IF( ABS(RetainedValue(A(P,J)))<=TINY(1.0_DP) ) RETURN
      Row=A(J,:)
      A(J,:)=A(P,:)
      A(P,:)=Row
      RHS=B(J,:)
      B(J,:)=B(P,:)
      B(P,:)=RHS
      DO L=J+1,N
        Factor=A(L,J)/A(J,J)
        A(L,J+1:N)=A(L,J+1:N)-Factor*A(J,J+1:N)
        B(L,:)=B(L,:)-Factor*B(J,:)
      END DO
    END DO
    DO J=N,1,-1
      DO L=J+1,N
        B(J,:)=B(J,:)-A(J,L)*B(L,:)
      END DO
      B(J,:)=B(J,:)/A(J,J)
    END DO
    Info=0
  END SUBROUTINE SolveRetainedSmall

  PURE SUBROUTINE ComputeSharedProcessTangent(M,NF,Active,U,D0,I0,V0,X,DP0,DQM,DX,Info, &
      Response_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    TYPE(RetainedReal), INTENT(inout), OPTIONAL :: Response_Option(4,4,M)
    REAL(DP), INTENT(in) :: U(M),D0(M),I0(M,3),V0(3)
    TYPE(RetainedReal), INTENT(in) :: X(:),DP0(M,4),DQM(NF+4)
    TYPE(RetainedReal), INTENT(out) :: DX(:)
    INTEGER, INTENT(out) :: Info
    TYPE(RetainedReal) :: V(3),D,I(3),KN(3,3),A4(4,4),B4(4,4),A3(3,3),B3(3,1)
    TYPE(RetainedReal) :: ED,SumR,SumB(3),DV(3),Term,DE
    REAL(DP) :: KD(3,3),KI(3,3,3),VC(3)
    INTEGER :: Q,J,L,Pass,N
    LOGICAL :: Move(3)
    N=4*M
    V=Retain(V0)
    Move=.FALSE.
    DO J=1,SIZE(Active)
      L=Active(J)-NF-1
      IF( L<1 .OR. L>3 ) CYCLE
      V(L)=V(L)-X(N+J)
      Move(L)=.TRUE.
    END DO
    VC=RetainedValue(V)
    ED=0.0_DP
    SumR=0.0_DP
    SumB=0.0_DP
    DX=0.0_DP
    DO Pass=1,2
      DO Q=1,M
        IF( Pass==2 .AND. PRESENT(Response_Option) )THEN
          B4=Response_Option(:,:,Q)
        ELSE
          D=D0(Q)+X(Q)
          DO J=1,3
            I(J)=I0(Q,J)+X(J*M+Q)
          END DO
          CALL ComputeCollisionPressure_OrderV_Point(D,I,KN)
          CALL ComputeCollisionClosureDerivatives_OrderV_Point(RetainedValue(D),RetainedValue(I),KD,KI)
          A4=0.0_DP
          A4(1,1)=1.0_DP
          A4(1,2:4)=V
          DO J=1,3
            A4(J+1,1)=Retain(SUM(KD(J,:)*VC))
            A4(J+1,J+1)=1.0_DP
            DO L=1,3
              A4(J+1,L+1)=A4(J+1,L+1)+SUM(KI(J,:,L)*VC)
            END DO
          END DO
          B4(:,1)=DP0(Q,:)
          B4(1,2:4)=-I
          B4(2:4,2:4)=-KN
          CALL SolveRetainedSmall(4,A4,B4,Info)
          IF( Info/=0 ) RETURN
          IF( PRESENT(Response_Option) ) Response_Option(:,:,Q)=B4
        END IF
        IF( Pass==1 )THEN
          ED=ED+U(Q)*D
          SumR=SumR+U(Q)*B4(1,1)
          SumB=SumB+U(Q)*B4(1,2:4)
        ELSE
          DO L=1,4
            Term=B4(L,1)
            DO J=1,3
              Term=Term+B4(L,J+1)*DV(J)
            END DO
            DX((L-1)*M+Q)=Term
          END DO
        END IF
      END DO
      IF( Pass==2 ) EXIT
      DO J=1,3
        A3(J,:)=V(J)*SumB
        A3(J,J)=A3(J,J)+1.0_DP+ED
        B3(J,1)=-DQM(NF+1+J)-V(J)*SumR
        IF( .NOT. Move(J) )THEN
          A3(J,:)=0.0_DP
          A3(J,J)=1.0_DP
          B3(J,1)=0.0_DP
        END IF
      END DO
      CALL SolveRetainedSmall(3,A3,B3,Info)
      IF( Info/=0 ) RETURN
      DV=B3(:,1)
    END DO
    DE=DQM(NF+1)
    DO Q=1,M
    DO J=1,3
      DE=DE+U(Q)*(V(J)*DX(J*M+Q)+DV(J)*(I0(Q,J)+X(J*M+Q)))
    END DO
    END DO
    DO J=1,SIZE(Active)
      L=Active(J)
      IF( L<=NF )THEN
        DX(N+J)=DQM(L)
      ELSE IF( L==NF+1 )THEN
        DX(N+J)=DE
      ELSE
        DX(N+J)=-DV(L-NF-1)
      END IF
    END DO
    Info=0
  END SUBROUTINE ComputeSharedProcessTangent

  PURE SUBROUTINE BuildScatteringComponents(M,NF,A,B,Kind,Gain,Loss,Label)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,A(:),B(:),Kind(:)
    REAL(DP), INTENT(in) :: Gain(:),Loss(:)
    INTEGER, INTENT(out) :: Label(M)
    INTEGER :: Q,E,L,R,Root,NE
    NE=M/(2*NF)
    DO Q=1,M
      Label(Q)=Q
    END DO
    DO E=1,SIZE(A)
      IF( Kind(E)/=0 .OR. MAX(Gain(E),Loss(E))<=0.0_DP ) CYCLE
      IF( (A(E)-1)/NE/=(B(E)-1)/NE ) CYCLE
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
      Root=Q
      DO WHILE( Label(Root)/=Root )
        Root=Label(Root)
      END DO
      Label(Q)=Root
    END DO
  END SUBROUTINE BuildScatteringComponents

  PURE SUBROUTINE NormalizeScatteringComponents &
    ( M,NF,W,U,D0,I0,V0,A,B,Kind,Gain,Loss,NumberIncrement,P,QM,NonScatteringEnergy_Option,Incoming_Option, &
      EnergyRatio_Option,Components_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,A(:),B(:),Kind(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Gain(:),Loss(:)
    TYPE(RetainedReal), INTENT(in) :: NumberIncrement(M)
    TYPE(RetainedReal), INTENT(inout) :: P(M,4),QM(NF+4)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: NonScatteringEnergy_Option
    TYPE(IncomingRadiationCache), INTENT(in), OPTIONAL :: Incoming_Option(M)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: EnergyRatio_Option(M)
    INTEGER, INTENT(in), OPTIONAL :: Components_Option(M)
    INTEGER :: Label(M),E,Q,J,L,R,Root,NE,Count,Pivot
    TYPE(RetainedReal) :: OldN(M),Positive(M),Target,SumPositive,Correction,Capacity,CountQ,SumOther
    TYPE(RetainedReal) :: Floor,SumFloor,SumExcess,Available
    TYPE(RetainedReal) :: OldP(M),Event,EnergyCoefficient,NewP,Scale
    REAL(DP) :: Best
    LOGICAL :: HasNegative
    NE=M/(2*NF)
    DO Q=1,M
      IF( PRESENT(Incoming_Option) )THEN
        OldN(Q)=Incoming_Option(Q)%Number
      ELSE
        OldN(Q)=Retain(D0(Q))
        DO J=1,3
          OldN(Q)=OldN(Q)+Retain(V0(J))*I0(Q,J)
        END DO
      END IF
    END DO
    ! Connected components retain number separately; no interpolation or
    ! redistribution across a zero-coefficient edge is invented.
    IF( PRESENT(Components_Option) )THEN
      Label=Components_Option
    ELSE
      CALL BuildScatteringComponents(M,NF,A,B,Kind,Gain,Loss,Label)
    END IF
    DO Root=1,M
      IF( Label(Root)/=Root ) CYCLE
      Count=0
      Target=0.0_DP
      SumPositive=0.0_DP
      SumFloor=0.0_DP
      HasNegative=.FALSE.
      Pivot=Root
      Best=-1.0_DP
      DO Q=1,M
        IF( Label(Q)/=Root ) CYCLE
        Count=Count+1
        Target=Target+W(Q)*(OldN(Q)+NumberIncrement(Q))
        Positive(Q)=W(Q)*(OldN(Q)+P(Q,1))
        Floor=W(Q)*MAX(128.0_DP*TINY(1.0_DP), &
          128.0_DP*EPSILON(1.0_DP)**2*MAX(RetainedValue(OldN(Q)),0.0_DP))
        SumFloor=SumFloor+Floor
        OldP(Q)=P(Q,1)
        IF( RetainedValue(Positive(Q)-Floor)<0.0_DP )THEN
          HasNegative=.TRUE.
          ! The closure recovery requires a positive density. Its numerical
          ! guard is local to this bin, never to the component population.
          ! In particular, repairing one bin must not fill every empty tail.
          Positive(Q)=Floor
        END IF
        SumPositive=SumPositive+Positive(Q)
        IF( RetainedValue(Positive(Q))<=Best ) CYCLE
        Best=RetainedValue(Positive(Q))
        Pivot=Q
      END DO
      IF( Count==1 )THEN
        ! No number-redistributing process reaches this bin at all.
        P(Root,1)=NumberIncrement(Root)
        CYCLE
      END IF
      IF( .NOT. HasNegative )THEN
        IF( ABS(RetainedValue(SumPositive-Target)) &
            <=512.0_DP*EPSILON(1.0_DP)**2*ABS(RetainedValue(SumPositive)) ) CYCLE
      END IF
      IF( RetainedValue(Target)<=0.0_DP ) CYCLE
      IF( RetainedValue(SumPositive)<=0.0_DP ) CYCLE
      Correction=Target-SumPositive
      Floor=W(Pivot)*MAX(128.0_DP*TINY(1.0_DP), &
        128.0_DP*EPSILON(1.0_DP)**2*MAX(RetainedValue(OldN(Pivot)),0.0_DP))
      Capacity=Positive(Pivot)-Floor
      ! Concentrate a small number-balance correction in the largest bin.
      ! Global rescaling needlessly changes all otherwise admissible tails.
      ! A genuinely large deficit instead removes the same relative fraction
      ! from the positive bins; no additive component-wide floor is used.
      IF( RetainedValue(Capacity+Correction)>0.0_DP )THEN
        Positive(Pivot)=Positive(Pivot)+Correction
      ELSE
        Available=Target-SumFloor
        SumExcess=SumPositive-SumFloor
        IF( RetainedValue(Available)<=0.0_DP .OR. RetainedValue(SumExcess)<=0.0_DP ) CYCLE
        DO Q=1,M
          IF( Label(Q)/=Root ) CYCLE
          Floor=W(Q)*MAX(128.0_DP*TINY(1.0_DP), &
            128.0_DP*EPSILON(1.0_DP)**2*MAX(RetainedValue(OldN(Q)),0.0_DP))
          Positive(Q)=Floor+Available*((Positive(Q)-Floor)/SumExcess)
        END DO
      END IF
      SumOther=0.0_DP
      DO Q=1,M
        IF( Label(Q)/=Root ) CYCLE
        IF( Q==Pivot ) CYCLE
        CountQ=Positive(Q)
        SumOther=SumOther+CountQ
        P(Q,1)=CountQ/W(Q)-OldN(Q)
      END DO
      P(Pivot,1)=(Target-SumOther)/W(Pivot)-OldN(Pivot)
      ! These are number-conserving scattering corrections. Accumulate the
      ! energy exchange by shared transfers to a pivot, never by subtracting
      ! two large completed energy moments.
      DO Q=1,M
        IF( Label(Q)/=Root .OR. Q==Pivot ) CYCLE
        Event=W(Q)*(P(Q,1)-OldP(Q))
        IF( PRESENT(EnergyRatio_Option) )THEN
          EnergyCoefficient=EnergyRatio_Option(Q)-EnergyRatio_Option(Pivot)
        ELSE
          EnergyCoefficient=Retain(U(Q))/W(Q)-Retain(U(Pivot))/W(Pivot)
        END IF
        QM(NF+1)=QM(NF+1)+EnergyCoefficient*Event
      END DO
    END DO
    IF( PRESENT(NonScatteringEnergy_Option) )THEN
      QM(NF+1)=NonScatteringEnergy_Option
      DO Q=1,M
        QM(NF+1)=QM(NF+1)+U(Q)*(P(Q,1)-NumberIncrement(Q))
      END DO
    END IF
  END SUBROUTINE NormalizeScatteringComponents

  PURE SUBROUTINE ProjectSharedProcessEndpoint &
    ( M,NF,Active,W,U,D0,I0,V0,Cap,A,B,Kind,Gain,Loss,Chi,Sigma,X0,P,QM,X,Info, &
      NumberIncrement_Option,NonScatteringEnergy_Option,Incoming_Option,IncomingEnergy_Option, &
      Response_Option,EnergyRatio_Option,Components_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:),A(:),B(:),Kind(:)
    TYPE(RetainedReal), INTENT(inout), OPTIONAL :: Response_Option(4,4,M)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: EnergyRatio_Option(M)
    INTEGER, INTENT(in), OPTIONAL :: Components_Option(M)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Cap,Gain(:),Loss(:),Chi(M),Sigma(M)
    TYPE(RetainedReal), INTENT(in) :: X0(:)
    TYPE(RetainedReal), INTENT(inout) :: P(M,4),QM(NF+4)
    TYPE(RetainedReal), INTENT(out) :: X(:)
    INTEGER, INTENT(out) :: Info
    TYPE(RetainedReal), INTENT(inout), OPTIONAL :: NumberIncrement_Option(M)
    TYPE(RetainedReal), INTENT(inout), OPTIONAL :: NonScatteringEnergy_Option
    TYPE(IncomingRadiationCache), INTENT(in), OPTIONAL :: Incoming_Option(M)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: IncomingEnergy_Option
    REAL(DP) :: Estimate(M),Need,Best,Capacity,Padding
    TYPE(RetainedReal) :: Event,EnergyCoefficient,V(3),D,I(3),INorm,DI(3),DD,OldK(3,3),NewK(3,3),DeltaG,Needed
    TYPE(RetainedReal) :: NativeN,NativeG(3),OldG(3),NewP1,Factor,NativeSlack,FluxScale
    INTEGER :: Pass,Q,S,NE,E,Partner,Donor,J,L,F
    LOGICAL :: Adjusted
    NE=M/(2*NF)
    DO Pass=1,12
      IF( PRESENT(NumberIncrement_Option) ) &
        CALL NormalizeScatteringComponents(M,NF,W,U,D0,I0,V0,A,B,Kind,Gain,Loss,NumberIncrement_Option,P,QM, &
          NonScatteringEnergy_Option,Incoming_Option,EnergyRatio_Option,Components_Option)
      ! Project the NATIVE flux before C2P. An unrealizable native flux can
      ! produce a negative intermediate primitive D even with positive N;
      ! treating that as a density defect would create needless number moves.
      DO Q=1,M
        IF( Sigma(Q)<=0.0_DP ) CYCLE
        NativeN=Retain(D0(Q))+P(Q,1)
        OldG=Retain(I0(Q,:))
        IF( PRESENT(Incoming_Option) )THEN
          ! Preserve the original addition order with the retained increment;
          ! (N_in + P) need not equal ((D_in + P) + v_in.I_in) bitwise.
          DO J=1,3
            NativeN=NativeN+Incoming_Option(Q)%NumberFrame(J)
            DO L=1,3
              OldG(J)=OldG(J)+Incoming_Option(Q)%FluxFrame(J,L)
            END DO
          END DO
        ELSE
          CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),OldK)
          DO J=1,3
            NativeN=NativeN+Retain(V0(J))*I0(Q,J)
            DO L=1,3
              OldG(J)=OldG(J)+Retain(V0(L))*OldK(J,L)
            END DO
          END DO
        END IF
        IF( RetainedValue(NativeN)<=0.0_DP ) CYCLE
        NativeG=OldG+P(Q,2:4)
        INorm=0.0_DP
        DO J=1,3
          INorm=INorm+NativeG(J)*NativeG(J)
        END DO
        INorm=RetainedSqrt(INorm)
        FluxScale=NativeN
        DO J=1,3
          FluxScale=FluxScale+RetainedAbs(OldG(J))+RetainedAbs(P(Q,J+1))
        END DO
        NativeSlack=(128.0_DP*EPSILON(1.0_DP)**2)*FluxScale
        IF( RetainedValue(INorm-NativeN+NativeSlack)<=0.0_DP ) CYCLE
        Factor=0.0_DP
        IF( RetainedValue(NativeN-NativeSlack)>0.0_DP .AND. RetainedValue(INorm)>0.0_DP ) &
          Factor=(NativeN-NativeSlack)/INorm
        DO J=1,3
          NewP1=NativeG(J)*Factor-OldG(J)
          DeltaG=NewP1-P(Q,J+1)
          P(Q,J+1)=NewP1
          QM(NF+1+J)=QM(NF+1+J)+U(Q)*DeltaG
        END DO
      END DO
      ! Every native first-moment collision term in this model exchanges
      ! momentum with the fluid. Re-sum these retained process increments;
      ! do not cancel a large previous angular source against its repair.
      QM(NF+2:NF+4)=0.0_DP
      DO Q=1,M
        QM(NF+2:NF+4)=QM(NF+2:NF+4)+U(Q)*P(Q,2:4)
      END DO
      X=X0
      ! When an almost empty bin is populated, starting C2P at its old tiny
      ! D gives a badly scaled closure linearization. The unboosted native
      ! moments are an admissible O(v) initial guess. Keep the existing
      ! incremental guess for bins whose population has not changed greatly.
      DO Q=1,M
        NativeN=Retain(D0(Q))+P(Q,1)
        NativeG=Retain(I0(Q,:))+P(Q,2:4)
        IF( PRESENT(Incoming_Option) )THEN
          DO J=1,3
            NativeN=NativeN+Incoming_Option(Q)%NumberFrame(J)
            DO L=1,3
              NativeG(J)=NativeG(J)+Incoming_Option(Q)%FluxFrame(J,L)
            END DO
          END DO
        ELSE
          CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),OldK)
          DO J=1,3
            NativeN=NativeN+Retain(V0(J))*I0(Q,J)
            DO L=1,3
              NativeG(J)=NativeG(J)+Retain(V0(L))*OldK(J,L)
            END DO
          END DO
        END IF
        D=D0(Q)+X(Q)
        IF( RetainedValue(NativeN)<=0.0_DP ) CYCLE
        IF( RetainedValue(NativeN)<=2.0_DP*RetainedValue(D) &
              .AND. RetainedValue(NativeN)>=0.5_DP*RetainedValue(D) ) CYCLE
        X(Q)=NativeN-D0(Q)
        DO J=1,3
          X(J*M+Q)=NativeG(J)-I0(Q,J)
        END DO
      END DO
      CALL RecoverSharedProcessState(M,NF,Active,W,U,D0,I0,V0,Cap,P,QM,X,Info,Incoming_Option,IncomingEnergy_Option,Response_Option)
      IF( Info==0 ) RETURN
      IF( Info<1001 .OR. Info>3000 ) RETURN
      Estimate=RetainedValue(D0+X(1:M))
      Adjusted=.FALSE.
      ! Project density violations by a SAME-SPECIES scattering transfer.
      ! Its two number legs share one event; its energy goes to the opposite
      ! fluid source. No family transfer is changed by this operation.
      DO Q=1,M
        S=(Q-1)/NE+1
        Padding=MAX(128.0_DP*EPSILON(1.0_DP)**2*D0(Q),128.0_DP*TINY(1.0_DP))
        Need=0.0_DP
        Needed=0.0_DP
        IF( Estimate(Q)<=0.0_DP ) Needed=-(D0(Q)+X(Q))+Padding
        IF( Estimate(Q)>Cap+8.0_DP*EPSILON(1.0_DP) ) Needed=Cap-(D0(Q)+X(Q))-Padding
        Need=RetainedValue(Needed)
        IF( Need==0.0_DP ) CYCLE
        Best=0.0_DP
        Donor=0
        DO E=1,SIZE(A)
          IF( Kind(E)/=0 .OR. MAX(Gain(E),Loss(E))<=0.0_DP ) CYCLE
          Partner=0
          IF( A(E)==Q ) Partner=B(E)
          IF( B(E)==Q ) Partner=A(E)
          IF( Partner==0 ) CYCLE
          IF( (Partner-1)/NE+1/=S ) CYCLE
          IF( Need>0.0_DP )THEN
            Capacity=W(Partner)*MAX(Estimate(Partner),0.0_DP)
          ELSE
            Capacity=W(Partner)*MAX(Cap-Estimate(Partner),0.0_DP)
          END IF
          IF( Capacity<=Best ) CYCLE
          Best=Capacity
          Donor=Partner
        END DO
        Event=W(Q)*Needed
        IF( Donor==0 .OR. W(Q)*ABS(Need)>=Best )THEN
          ! A pair reaction may supply both number legs when same-species
          ! redistribution has insufficient capacity. Its family source is
          ! identically zero, including this projection correction.
          Donor=0
          Best=0.0_DP
          DO E=1,SIZE(A)
            IF( Kind(E)/=1 .OR. MAX(Gain(E),Loss(E))<=0.0_DP ) CYCLE
            Partner=0
            IF( A(E)==Q ) Partner=B(E)
            IF( B(E)==Q ) Partner=A(E)
            IF( Partner==0 ) CYCLE
            IF( ((Partner-1)/NE+2)/2/=(S+1)/2 ) CYCLE
            IF( Need>0.0_DP )THEN
              Capacity=W(Partner)*MAX(Cap-Estimate(Partner),0.0_DP)
            ELSE
              Capacity=W(Partner)*MAX(Estimate(Partner),0.0_DP)
            END IF
            IF( Capacity<=Best ) CYCLE
            Best=Capacity
            Donor=Partner
          END DO
          IF( Donor/=0 .AND. W(Q)*ABS(Need)<Best )THEN
            P(Q,1)=P(Q,1)+Event/W(Q)
            P(Donor,1)=P(Donor,1)+Event/W(Donor)
            IF( PRESENT(NumberIncrement_Option) )THEN
              NumberIncrement_Option(Q)=NumberIncrement_Option(Q)+Event/W(Q)
              NumberIncrement_Option(Donor)=NumberIncrement_Option(Donor)+Event/W(Donor)
            END IF
            IF( PRESENT(EnergyRatio_Option) )THEN
              EnergyCoefficient=EnergyRatio_Option(Q)+EnergyRatio_Option(Donor)
            ELSE
              EnergyCoefficient=Retain(U(Q))/W(Q)+Retain(U(Donor))/W(Donor)
            END IF
            QM(NF+1)=QM(NF+1)+EnergyCoefficient*Event
            IF( PRESENT(NonScatteringEnergy_Option) ) &
              NonScatteringEnergy_Option=NonScatteringEnergy_Option+EnergyCoefficient*Event
            Estimate(Q)=Estimate(Q)+Need
            Estimate(Donor)=Estimate(Donor)+RetainedValue(Event/W(Donor))
            Adjusted=.TRUE.
            CYCLE
          END IF
          ! Charged-current repair is permitted only with a nonzero actual
          ! absorption coefficient. Its family/energy change is paired, not
          ! discarded or attributed to a number-conserving reaction.
          IF( Chi(Q)<=0.0_DP ) RETURN
          P(Q,1)=P(Q,1)+Event/W(Q)
          IF( PRESENT(NumberIncrement_Option) ) NumberIncrement_Option(Q)=NumberIncrement_Option(Q)+Event/W(Q)
          F=(S+1)/2
          QM(F)=QM(F)+MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*Event
          IF( PRESENT(EnergyRatio_Option) )THEN
            EnergyCoefficient=EnergyRatio_Option(Q)
          ELSE
            EnergyCoefficient=Retain(U(Q))/W(Q)
          END IF
          QM(NF+1)=QM(NF+1)+EnergyCoefficient*Event
          IF( PRESENT(NonScatteringEnergy_Option) ) &
            NonScatteringEnergy_Option=NonScatteringEnergy_Option+EnergyCoefficient*Event
          Estimate(Q)=Estimate(Q)+Need
          Adjusted=.TRUE.
          CYCLE
        END IF
        P(Q,1)=P(Q,1)+Event/W(Q)
        P(Donor,1)=P(Donor,1)-Event/W(Donor)
        IF( PRESENT(EnergyRatio_Option) )THEN
          EnergyCoefficient=EnergyRatio_Option(Q)-EnergyRatio_Option(Donor)
        ELSE
          EnergyCoefficient=Retain(U(Q))/W(Q)-Retain(U(Donor))/W(Donor)
        END IF
        QM(NF+1)=QM(NF+1)+EnergyCoefficient*Event
        Estimate(Q)=Estimate(Q)+Need
        Estimate(Donor)=Estimate(Donor)-RetainedValue(Event/W(Donor))
        Adjusted=.TRUE.
      END DO
      IF( Adjusted ) CYCLE
      ! An angular repair is an isoenergetic momentum transfer. Keep the
      ! native number increment fixed and pair its native flux change with
      ! the fluid momentum. The moving-frame work is recovered afterward.
      V=Retain(V0)
      DO J=1,SIZE(Active)
        L=Active(J)-NF-1
        IF( L>=1 .AND. L<=3 ) V(L)=V(L)-X0(4*M+J)
      END DO
      DO Q=1,M
        IF( Sigma(Q)<=0.0_DP ) CYCLE
        D=D0(Q)+X(Q)
        IF( RetainedValue(D)<=0.0_DP ) CYCLE
        INorm=0.0_DP
        DO J=1,3
          I(J)=I0(Q,J)+X(J*M+Q)
          INorm=INorm+I(J)*I(J)
        END DO
        INorm=RetainedSqrt(INorm)
        IF( RetainedValue(INorm)<=RetainedValue(D) ) CYCLE
        DI=I*(D/INorm-1.0_DP)
        DD=0.0_DP
        DO J=1,3
          DD=DD-V(J)*DI(J)
        END DO
        IF( RetainedValue(D+DD)<=0.0_DP ) RETURN
        CALL ComputeCollisionPressure_OrderV_Point(D,I,OldK)
        CALL ComputeCollisionPressure_OrderV_Point(D+DD,I+DI,NewK)
        DO J=1,3
          DeltaG=DI(J)
          DO L=1,3
            DeltaG=DeltaG+V(L)*(NewK(J,L)-OldK(J,L))
          END DO
          P(Q,J+1)=P(Q,J+1)+DeltaG
          QM(NF+1+J)=QM(NF+1+J)+U(Q)*DeltaG
        END DO
        Adjusted=.TRUE.
      END DO
      IF( .NOT. Adjusted ) RETURN
    END DO
    Info=4000
  END SUBROUTINE ProjectSharedProcessEndpoint

  PURE SUBROUTINE RecoverSharedProcessState &
    (M,NF,Active,W,U,D0,I0,V0,Cap,P,QM,X,Info,Incoming_Option,IncomingEnergy_Option, &
      Response_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! P contains the integrated native collision terms, QM their separately
    ! paired, process-resolved comoving moments.  Neither is reconstructed from
    ! the rounded primitive endpoint. No EOS or opacity evaluation occurs here.
    INTEGER, INTENT(in) :: M,NF,Active(:)
    TYPE(RetainedReal), INTENT(inout), OPTIONAL :: Response_Option(4,4,M)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Cap
    TYPE(RetainedReal), INTENT(in) :: P(M,4),QM(NF+4)
    TYPE(RetainedReal), INTENT(inout) :: X(:)
    INTEGER, INTENT(out) :: Info
    TYPE(IncomingRadiationCache), INTENT(in), OPTIONAL :: Incoming_Option(M)
    TYPE(RetainedReal), INTENT(in), OPTIONAL :: IncomingEnergy_Option
    TYPE(RetainedReal) :: V(3),DV(3),DI(3),DK(3,3),K0(3,3),KN(3,3)
    TYPE(RetainedReal) :: A4(4,4),B4(4,4),A3(3,3),B3(3,1),DVStep(3)
    TYPE(RetainedReal) :: ED,ED0,SumR,SumB(3),F(4),DN,DVI,Term
    TYPE(RetainedReal) :: VelocityResidual(3),EnergyChange
    REAL(DP) :: KD(3,3),KI(3,3,3),VC(3),D,IV(3),Bound,Scale
    LOGICAL :: Converged,Finished,Move(3)
    INTEGER :: Q,J,L,Sweep,Pass,N
    N=4*M
    DV=0.0_DP
    Move=.FALSE.
    DO J=1,SIZE(Active)
      L=Active(J)-NF-1
      IF( L<1 .OR. L>3 ) CYCLE
      DV(L)=-X(N+J)
      Move(L)=.TRUE.
    END DO
    IF( PRESENT(IncomingEnergy_Option) )THEN
      ED0=IncomingEnergy_Option
    ELSE
      ED0=0.0_DP
      DO Q=1,M
        ED0=ED0+U(Q)*Retain(D0(Q))
      END DO
    END IF
    Info=1
    DO Sweep=1,20
      V=V0+DV
      VC=RetainedValue(V)
      ED=0.0_DP
      SumR=0.0_DP
      SumB=0.0_DP
      Finished=.TRUE.
      ! Two passes avoid an O(M) array of local 4x4 inverse/response matrices
      ! in each GPU thread. Only the small Schur system is retained.
      DO Pass=1,2
        DO Q=1,M
          IF( Pass==2 .AND. PRESENT(Response_Option) )THEN
            B4=Response_Option(:,:,Q)
          ELSE
            D=RetainedValue(D0(Q)+X(Q))
            IF( D<=TINY(1.0_DP) .OR. .NOT. IEEE_IS_FINITE(D) )THEN
              Info=1000+Q
              RETURN
            END IF
            DO J=1,3
              DI(J)=X(J*M+Q)
              IV(J)=RetainedValue(I0(Q,J)+DI(J))
            END DO
            IF( PRESENT(Incoming_Option) )THEN
              K0=Incoming_Option(Q)%Pressure
              CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DI,DK,Incoming_Option(Q))
            ELSE
              CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),K0)
              CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DI,DK)
            END IF
            KN=K0+DK
            CALL ComputeCollisionDensityIncrement_OrderV_Point(X(Q),DI,I0(Q,:),V,DV,DN,DVI)
            F(1)=DN-P(Q,1)
            DO J=1,3
              F(J+1)=DI(J)-P(Q,J+1)
              DO L=1,3
                F(J+1)=F(J+1)+V(L)*DK(J,L)+DV(L)*K0(J,L)
              END DO
            END DO
            IF( Pass==1 )THEN
              ED=ED+U(Q)*X(Q)
              Scale=ABS(RetainedValue(X(Q)))+ABS(RetainedValue(P(Q,1)))
              DO J=1,3
                Scale=Scale+ABS(RetainedValue(V(J)*DI(J)))+ABS(RetainedValue(DV(J)*I0(Q,J)))
              END DO
              Bound=MAX(128.0_DP*EPSILON(1.0_DP)**2*Scale,1.0e-300_DP)
              Finished=Finished .AND. ABS(RetainedValue(F(1)))<=Bound
              DO J=1,3
                Scale=ABS(RetainedValue(DI(J)))+ABS(RetainedValue(P(Q,J+1)))
                DO L=1,3
                  Scale=Scale+ABS(RetainedValue(V(L)*DK(J,L)))+ABS(RetainedValue(DV(L)*K0(J,L)))
                END DO
                ! A nearly free-streaming pressure component can be tiny by
                ! cancellation of O(D) tensor terms. Its rounding bound must
                ! include the parent increment scale, not demand relative
                ! retained accuracy of that vanishing tensor component.
                Scale=Scale+SUM(ABS(VC))*(ABS(RetainedValue(X(Q)))+SUM(ABS(RetainedValue(DI)))) &
                  +SUM(ABS(RetainedValue(DV)))*D0(Q)
                Bound=MAX(128.0_DP*EPSILON(1.0_DP)**2*Scale,1.0e-300_DP)
                Finished=Finished .AND. ABS(RetainedValue(F(J+1)))<=Bound
              END DO
            END IF
            CALL ComputeCollisionClosureDerivatives_OrderV_Point(D,IV,KD,KI)
            A4=0.0_DP
            A4(1,1)=1.0_DP
            A4(1,2:4)=V
            DO J=1,3
              A4(J+1,1)=Retain(SUM(KD(J,:)*VC))
              A4(J+1,J+1)=1.0_DP
              DO L=1,3
                A4(J+1,L+1)=A4(J+1,L+1)+SUM(KI(J,:,L)*VC)
              END DO
            END DO
            B4(:,1)=-F
            B4(1,2:4)=-(I0(Q,:)+DI)
            B4(2:4,2:4)=-KN
            CALL SolveRetainedSmall(4,A4,B4,Info)
            IF( Info/=0 ) RETURN
            IF( PRESENT(Response_Option) ) Response_Option(:,:,Q)=B4
          END IF
          IF( Pass==1 )THEN
            SumR=SumR+U(Q)*B4(1,1)
            SumB=SumB+U(Q)*B4(1,2:4)
          ELSE
            Term=B4(1,1)
            DO J=1,3
              Term=Term+B4(1,J+1)*DVStep(J)
            END DO
            X(Q)=X(Q)+Term
            DO L=1,3
              Term=B4(L+1,1)
              DO J=1,3
                Term=Term+B4(L+1,J+1)*DVStep(J)
              END DO
              X(L*M+Q)=X(L*M+Q)+Term
            END DO
          END IF
        END DO
        IF( Pass==2 ) EXIT
        DO J=1,3
          VelocityResidual(J)=DV(J)+QM(NF+1+J)+V(J)*ED+DV(J)*ED0
          Scale=ABS(RetainedValue(DV(J)))+ABS(RetainedValue(QM(NF+1+J))) &
                  +ABS(RetainedValue(V(J)*ED))+ABS(RetainedValue(DV(J)*ED0))
          Bound=MAX(128.0_DP*EPSILON(1.0_DP)**2*Scale,1.0e-300_DP)
          IF( Move(J) ) Finished=Finished .AND. ABS(RetainedValue(VelocityResidual(J)))<=Bound
          B3(J,1)=-VelocityResidual(J)-V(J)*SumR
          A3(J,:)=V(J)*SumB
          A3(J,J)=A3(J,J)+1.0_DP+ED0+ED
          IF( .NOT. Move(J) )THEN
            A3(J,:)=0.0_DP
            A3(J,J)=1.0_DP
            B3(J,1)=0.0_DP
          END IF
        END DO
        IF( Finished ) EXIT
        CALL SolveRetainedSmall(3,A3,B3,Info)
        IF( Info/=0 ) RETURN
        DVStep=B3(:,1)
      END DO
      IF( Finished ) EXIT
      DV=DV+DVStep
    END DO
    Info=1
    IF( .NOT. Finished )THEN
      Info=2000
      RETURN
    END IF
    IF( .NOT. CollisionStateAdmissible_OrderV_Point(M,D0,I0,X,Cap) )THEN
      Info=3000
      RETURN
    END IF
    V=V0+DV
    EnergyChange=QM(NF+1)
    DO Q=1,M
      DO J=1,3
        EnergyChange=EnergyChange+U(Q)*(V(J)*X(J*M+Q)+DV(J)*I0(Q,J))
      END DO
    END DO
    DO J=1,SIZE(Active)
      L=Active(J)
      IF( L<=NF )THEN
        X(N+J)=QM(L)
      ELSE IF( L==NF+1 )THEN
        X(N+J)=EnergyChange
      ELSE
        X(N+J)=-DV(L-NF-1)
      END IF
    END DO
    X(N+SIZE(Active)+1:)=0.0_DP
    Info=0
  END SUBROUTINE RecoverSharedProcessState

  PURE SUBROUTINE RecoverRadiationExchange(M,NF,Active,W,U,D0,I0,V0,X,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! The trial radiation determines the fluid exchange.  At fixed primitive
    ! radiation the Order-v momentum balance is a linear 3 by 3 system.
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3)
    TYPE(RetainedReal), INTENT(inout) :: X(:)
    INTEGER, INTENT(out) :: Info
    TYPE(RetainedReal) :: A(3,3),B(3),Row(3),Temp,Factor,DI(3),DK(3,3),K0(3,3),KN(3,3)
    TYPE(RetainedReal) :: Exchange(NF+4),V(3),DV(3),DN,DeltaVI
    REAL(DP) :: G
    INTEGER :: Q,J,L,P,N,K,S,F
    N=4*M
    K=NF+4
    A=0.0_DP
    B=0.0_DP
    DO J=1,3
      A(J,J)=1.0_DP
    END DO
    DO Q=1,M
      DO J=1,3
        DI(J)=X(J*M+Q)
      END DO
      CALL ComputeCollisionPressure_OrderV_Point(Retain(D0(Q)),Retain(I0(Q,:)),K0)
      CALL ComputeCollisionPressureIncrement_OrderV_Point(D0(Q),I0(Q,:),X(Q),DI,DK)
      KN=K0+DK
      DO J=1,3
        B(J)=B(J)+U(Q)*(DI(J)+V0(J)*X(Q))
        A(J,J)=A(J,J)+U(Q)*(D0(Q)+X(Q))
        DO L=1,3
          B(J)=B(J)+U(Q)*(V0(L)*DK(J,L))
          A(J,L)=A(J,L)+U(Q)*KN(J,L)
        END DO
      END DO
    END DO
    ! A frozen fluid momentum component has zero exchange, not an extra
    ! equation constraining the radiation pressure tensor.
    DO J=1,3
      IF( ANY(Active==NF+1+J) ) CYCLE
      A(J,:)=0.0_DP
      A(J,J)=1.0_DP
      B(J)=0.0_DP
    END DO
    Info=1
    DO J=1,3
      P=J-1+MAXLOC(ABS(RetainedValue(A(J:3,J))),DIM=1)
      IF( ABS(RetainedValue(A(P,J)))<=TINY(1.0_DP) ) RETURN
      Row=A(J,:)
      A(J,:)=A(P,:)
      A(P,:)=Row
      Temp=B(J)
      B(J)=B(P)
      B(P)=Temp
      DO L=J+1,3
        Factor=A(L,J)/A(J,J)
        A(L,J+1:3)=A(L,J+1:3)-Factor*A(J,J+1:3)
        B(L)=B(L)-Factor*B(J)
      END DO
    END DO
    DO J=3,1,-1
      DO L=J+1,3
        B(J)=B(J)-A(J,L)*B(L)
      END DO
      B(J)=B(J)/A(J,J)
    END DO
    DO J=1,SIZE(Active)
      IF( Active(J)>NF+1 ) X(N+J)=B(Active(J)-NF-1)
    END DO
    DV=-B
    V=V0+DV
    Exchange=0.0_DP
    DO Q=1,M
      DO J=1,3
        DI(J)=X(J*M+Q)
      END DO
      CALL ComputeCollisionDensityIncrement_OrderV_Point(X(Q),DI,I0(Q,:),V,DV,DN,DeltaVI)
      S=(Q-1)/(M/(2*NF))+1
      F=(S+1)/2
      G=MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      Exchange(F)=Exchange(F)+G*W(Q)*DN
      Exchange(NF+1)=Exchange(NF+1)+U(Q)*(X(Q)+2.0_DP*DeltaVI)
    END DO
    DO J=1,SIZE(Active)
      ! Keep the velocity from the analytic momentum solve. Updating it
      ! again from a rounded moment would change the frame of the family
      ! and energy sums just computed, reintroducing an avoidable mismatch.
      IF( Active(J)<=NF+1 ) X(N+J)=Exchange(Active(J))
    END DO
    Info=0
  END SUBROUTINE RecoverRadiationExchange

  PURE SUBROUTINE ReconstructRadiationInvariants(M,NF,Active,W,U,D0,I0,V0,Cap,X,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Only inactive invariants constrain radiation here. Active exchanges are
    ! recovered from radiation, never matched to an independent source value.
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Cap
    TYPE(RetainedReal), INTENT(inout) :: X(:)
    INTEGER, INTENT(out) :: Info
    REAL(DP) :: H(NF+4),Scale(NF+4),JH(NF+4,4*M)
    LOGICAL :: Inactive(NF+4)
    INTEGER :: Pass
    Inactive=.TRUE.
    Inactive(Active)=.FALSE.
    DO Pass=0,10
      Info=1
      IF( .NOT. CollisionStateAdmissible_OrderV_Point(M,D0,I0,X,Cap) ) RETURN
      CALL RecoverRadiationExchange(M,NF,Active,W,U,D0,I0,V0,X,Info)
      IF( Info/=0 ) RETURN
      CALL ComputeConservationGeometry(M,NF,Active,W,U,D0,I0,V0,X,H,Scale,JH)
      IF( ALL(ABS(H)<=MAX(32.0_DP*EPSILON(1.0_DP)**2*Scale,1.0e-300_DP) .OR. .NOT.Inactive) ) RETURN
      IF( Pass==10 ) EXIT
      ! Keep the just-recovered active exchange while removing the inactive
      ! moments. Re-recover all active transfers after this radiation change.
      CALL ReconstructConservativeRadiation(M,NF,Active,W,U,D0,I0,V0,Cap,X,Info,RetainedPrecision_Option=.TRUE.)
      IF( Info/=0 ) RETURN
    END DO
    Info=1
  END SUBROUTINE ReconstructRadiationInvariants

  PURE SUBROUTINE SelectConservativeSupport(M,NF,W,U,D,I,Cap,Z,Info, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(out), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M,NF
    REAL(DP), INTENT(in) :: W(M),U(M),D(M),I(M,3),Cap
    REAL(DP), INTENT(out) :: Z(:,:)
    INTEGER, INTENT(out) :: Info
    REAL(DP) :: Basis(NF+4,4*M),Capacity(4*M),A(NF+4,NF+4),B(NF+4),Rows(NF+4)
    INTEGER :: Pivots(NF+4),Q,S,F,J,K
    K=NF+4
    Basis=0.0_DP
    Z=0.0_DP
    DO Q=1,M
      S=(Q-1)/(M/(2*NF))+1
      F=(S+1)/2
      Basis(F,Q)=MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)*W(Q)
      Basis(NF+1,Q)=U(Q)
      Capacity(Q)=MIN(D(Q),Cap-D(Q))
      DO J=1,3
        Basis(NF+1+J,J*M+Q)=U(Q)
        Capacity(J*M+Q)=D(Q)-Norm3(I(Q,:))
      END DO
    END DO
    CALL SelectConservationPivots(K,4*M,Basis,Capacity,Pivots,Info)
    IF( Info/=0 ) RETURN
    IF( PRESENT(LiftRows_Option) ) LiftRows_Option=Pivots
    DO J=1,K
      A=Basis(:,Pivots)
      Rows=MAXVAL(ABS(A),DIM=2)
      DO Q=1,K
        A(Q,:)=A(Q,:)/Rows(Q)
      END DO
      B=0.0_DP
      B(J)=1.0_DP/Rows(J)
      CALL SolveSmall(K,A,B,Info)
      IF( Info/=0 ) RETURN
      IF( PRESENT(LiftRows_Option) )THEN
        Z(:,J)=B
      ELSE
        Z(Pivots,J)=B
      END IF
    END DO
  END SUBROUTINE SelectConservativeSupport

  PURE SUBROUTINE ComputeIndependentMetrics(NE,NS,Active,W,U,D0,I0,X,DT,Physical,Independent, &
    Kappa,DensityLoss,HScale,GScale,InitialNorm,SetInitial,RadiationError,ConservationError,Weights,Merit)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: NE,NS,Active(:)
    REAL(DP), INTENT(in) :: W(NE*NS),U(NE*NS),D0(NE*NS),I0(NE*NS,3),DT,Physical(:),Independent(4*NE*NS)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(in) :: Kappa(NE*NS),DensityLoss(NE*NS),HScale(NS/2+4),GScale(NS/2+4)
    REAL(DP), INTENT(inout) :: InitialNorm(NS)
    LOGICAL, INTENT(in) :: SetInitial
    REAL(DP), INTENT(out) :: RadiationError,ConservationError,Weights(:),Merit
    INTEGER :: S,J,First,Last,N
    N=4*NE*NS
    CALL ComputeCollisionMetrics_OrderV_Point(NE,NS,Active,W,U,D0,I0,X,DT,Physical,Kappa,DensityLoss,HScale,GScale, &
      InitialNorm,SetInitial,RadiationError,ConservationError,Weights,Merit)
    Weights(N+1:)=0.0_DP
    RadiationError=0.0_DP
    DO J=0,3
    DO S=1,NS
      First=J*NE*NS+(S-1)*NE+1
      Last=First+NE-1
      RadiationError=MAX(RadiationError,SQRT(SUM((Weights(First:Last)*Independent(First:Last))**2)))
    END DO
    END DO
    Merit=0.5_DP*SUM((Weights(1:N)*Independent)**2)
  END SUBROUTINE ComputeIndependentMetrics

  PURE SUBROUTINE ComputeFluidSourceMetrics(M,NF,Active,Fluid0,X,Physical,GScale,Error,Weights_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: Fluid0(3),Physical(:),GScale(NF+4)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(out) :: Error(NF+4)
    REAL(DP), INTENT(out), OPTIONAL :: Weights_Option(:)
    TYPE(RetainedReal) :: Current,Mapped
    REAL(DP) :: Denom,Maximum,MomentumScale
    INTEGER :: J,F,N,K
    N=4*M
    K=NF+4
    Error=0.0_DP
    IF( PRESENT(Weights_Option) ) Weights_Option=0.0_DP
    ! Momentum is scaled as a vector so an exactly zero transverse source
    ! does not become an independent relative-accuracy requirement.
    Maximum=MAXVAL(GScale(NF+2:K))
    MomentumScale=0.0_DP
    IF( Maximum>0.0_DP ) MomentumScale=Maximum*SQRT(SUM((GScale(NF+2:K)/Maximum)**2))
    DO J=1,SIZE(Active)
      F=Active(J)
      Denom=GScale(F)
      IF( F<=MIN(NF,2) )THEN
        ! Use the composition scale of the nested iterate-change check,
        ! applied here to the accepted-state SOURCE residual. T is a
        ! radiation gain: Y_new=Y_in-T, Y_map=Y_new+g, g=T-h*S.
        Current=Retain(Fluid0(F+1))-X(N+J)
        Mapped=Current+Physical(N+K+J)
        Denom=ABS(RetainedValue(Current))+ABS(RetainedValue(Mapped))+ABS(Fluid0(F+1))
      ELSE IF( F>NF+1 )THEN
        Denom=MomentumScale
      END IF
      Error(F)=ABS(Physical(N+K+J))/MAX(Denom,SQRT(TINY(1.0_DP)))
      IF( PRESENT(Weights_Option) ) Weights_Option(J)=1.0_DP/MAX(Denom,SQRT(TINY(1.0_DP)))
    END DO
  END SUBROUTINE ComputeFluidSourceMetrics

  PURE SUBROUTINE ComputeCoupledSourceMetrics(M,NF,Active,Fluid0,X,Physical,GScale, &
    RadiationTolerance,SourceTolerance,Independent,Error,Weights,Merit)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: Fluid0(3),Physical(:),GScale(NF+4),RadiationTolerance,SourceTolerance
    REAL(DP), INTENT(in) :: Independent(4*M)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(out) :: Error(NF+4),Merit
    REAL(DP), INTENT(inout) :: Weights(:)
    INTEGER :: N,K
    N=4*M
    K=NF+4
    ! Keep the established radiation weights. Scale the source equations by
    ! their OWN stopping tolerance, with a common radiation-tolerance factor.
    ! Conservation rows remain constraints, not terms to trade against a source.
    Weights(N+1:)=0.0_DP
    CALL ComputeFluidSourceMetrics(M,NF,Active,Fluid0,X,Physical,GScale,Error,Weights(N+K+1:))
    Weights(N+K+1:)=Weights(N+K+1:)*(RadiationTolerance/SourceTolerance)
    Merit=ComputeCoupledMerit(N,Weights,Independent,Physical)
  END SUBROUTINE ComputeCoupledSourceMetrics

  PURE SUBROUTINE PrepareCoupledSearchWeights(M,NF,Active,Fluid0,X,Physical,Independent, &
    Direction,Action,RadiationTolerance,SourceTolerance,Weights,Merit)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: Fluid0(3),Physical(:),Independent(4*M),Direction(:),Action(:)
    REAL(DP), INTENT(in) :: RadiationTolerance,SourceTolerance
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(inout) :: Weights(:)
    REAL(DP), INTENT(out) :: Merit
    TYPE(RetainedReal) :: Current,Mapped
    REAL(DP) :: Denom,DY,DMap
    INTEGER :: J,F,N,K,Row
    N=4*M
    K=NF+4
    DO J=1,SIZE(Active)
      F=Active(J)
      IF( F>MIN(NF,2) ) CYCLE
      Row=N+K+J
      Current=Retain(Fluid0(F+1))-X(N+J)
      Mapped=Current+Physical(Row)
      DY=-Direction(N+J)
      DMap=DY+Action(Row)
      ! A fixed reference for this entire search, including fallbacks. The
      ! triangle inequality bounds the linearized composition scale along
      ! the unit Newton step. This avoids a vanishing initial reference but
      ! does not change the actual-state composition stopping criterion.
      Denom=ABS(RetainedValue(Current))+ABS(RetainedValue(Mapped))+ABS(Fluid0(F+1)) &
        +ABS(DY)+ABS(DMap)
      Weights(Row)=(RadiationTolerance/SourceTolerance)/MAX(Denom,SQRT(TINY(1.0_DP)))
    END DO
    Merit=ComputeCoupledMerit(N,Weights,Independent,Physical)
  END SUBROUTINE PrepareCoupledSearchWeights

  PURE FUNCTION ComputeCoupledMerit(N,Weights,Independent,Physical) RESULT(Merit)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: N
    REAL(DP), INTENT(in) :: Weights(:),Independent(N),Physical(:)
    REAL(DP) :: Merit
    Merit=0.5_DP*(SUM((Weights(1:N)*Independent)**2) &
      +SUM((Weights(N+1:)*Physical(N+1:))**2))
  END FUNCTION ComputeCoupledMerit

  PURE FUNCTION ComputeCoupledSlope(N,Weights,Independent,Physical,Action) RESULT(Slope)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: N
    REAL(DP), INTENT(in) :: Weights(:),Independent(N),Physical(:),Action(:)
    REAL(DP) :: Slope
    Slope=SUM((Weights(1:N)*Independent)*(Weights(1:N)*Action(1:N))) &
      +SUM((Weights(N+1:)*Physical(N+1:))*(Weights(N+1:)*Action(N+1:)))
  END FUNCTION ComputeCoupledSlope


  PURE SUBROUTINE ApplyMomentMap(M,NF,W,U,R,Moment)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF
    REAL(DP), INTENT(in) :: W(M),U(M)
    TYPE(RetainedReal), INTENT(in) :: R(4*M)
    TYPE(RetainedReal), INTENT(out) :: Moment(NF+4)
    INTEGER :: Q,S,F,J
    REAL(DP) :: G
    Moment=0.0_DP
    DO Q=1,M
      S=(Q-1)/(M/(2*NF))+1
      F=(S+1)/2
      G=MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      Moment(F)=Moment(F)+G*W(Q)*R(Q)
      Moment(NF+1)=Moment(NF+1)+U(Q)*R(Q)
      DO J=1,3
        Moment(NF+1+J)=Moment(NF+1+J)+U(Q)*R(J*M+Q)
      END DO
    END DO
  END SUBROUTINE ApplyMomentMap

  PURE SUBROUTINE InvertExchangeFrame(NF,V,G,Y)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: NF
    TYPE(RetainedReal), INTENT(in) :: V(3),G(NF+4)
    TYPE(RetainedReal), INTENT(out) :: Y(NF+4)
    TYPE(RetainedReal) :: V2,Dot
    INTEGER :: J
    V2=0.0_DP
    Dot=0.0_DP
    DO J=1,3
      V2=V2+V(J)*V(J)
      Dot=Dot+V(J)*G(NF+1+J)
    END DO
    Y=G
    Y(NF+1)=(G(NF+1)-Dot)/(1.0_DP-V2)
    DO J=1,3
      Y(NF+1+J)=G(NF+1+J)-V(J)*Y(NF+1)
    END DO
  END SUBROUTINE InvertExchangeFrame

  PURE SUBROUTINE ComputeIndependentRadiationResidual(M,NF,W,U,V,Z0,F,G,Fhat, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M,NF
    REAL(DP), INTENT(in) :: W(M),U(M),Z0(:,:)
    TYPE(RetainedReal), INTENT(in) :: V(3),F(4*M),G(NF+4)
    REAL(DP), INTENT(out) :: Fhat(4*M)
    TYPE(RetainedReal) :: Moment(NF+4),Inverse(NF+4),Term
    INTEGER :: Q,J
    CALL ApplyMomentMap(M,NF,W,U,F,Moment)
    CALL InvertExchangeFrame(NF,V,G,Inverse)
    Inverse=Inverse-Moment
    IF( PRESENT(LiftRows_Option) )THEN
      Fhat=RetainedValue(F)
      DO Q=1,SIZE(LiftRows_Option)
        Term=F(LiftRows_Option(Q))
        DO J=1,NF+4
          Term=Term+Z0(Q,J)*Inverse(J)
        END DO
        Fhat(LiftRows_Option(Q))=RetainedValue(Term)
      END DO
      RETURN
    END IF
    DO Q=1,4*M
      Term=F(Q)
      DO J=1,NF+4
        Term=Term+Z0(Q,J)*Inverse(J)
      END DO
      Fhat(Q)=RetainedValue(Term)
    END DO
  END SUBROUTINE ComputeIndependentRadiationResidual

  PURE SUBROUTINE ComputeIndependentRadiationColumn(M,NF,Active,W,U,V,Z0,G,VC,Column, &
      LiftRows_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in), OPTIONAL :: LiftRows_Option(:)
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),Z0(:,:),VC(3)
    TYPE(RetainedReal), INTENT(in) :: V(3),G(NF+4)
    REAL(DP), INTENT(inout) :: Column(:)
    TYPE(RetainedReal) :: GC(NF+4),Inverse(NF+4),F(4*M)
    INTEGER :: J,N,K
    N=4*M
    K=NF+4
    GC=0.0_DP
    DO J=1,SIZE(Active)
      GC(Active(J))=Column(N+K+J)
    END DO
    CALL InvertExchangeFrame(NF,V,G,Inverse)
    DO J=1,3
      GC(NF+1)=GC(NF+1)-VC(J)*Inverse(NF+1+J)
      GC(NF+1+J)=GC(NF+1+J)-VC(J)*Inverse(NF+1)
    END DO
    F=Retain(Column(1:N))
    CALL ComputeIndependentRadiationResidual(M,NF,W,U,V,Z0,F,GC,Column(1:N), &
      LiftRows_Option=LiftRows_Option)
  END SUBROUTINE ComputeIndependentRadiationColumn

  PURE SUBROUTINE ComputeConservationGeometry(M,NF,Active,W,U,D0,I0,V0,X,H,Scale,JH)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3)
    TYPE(RetainedReal), INTENT(in) :: X(:)
    REAL(DP), INTENT(out) :: H(NF+4),Scale(NF+4),JH(NF+4,4*M)
    TYPE(RetainedReal) :: C0(M),C1(M,3),Source(NF+4),Increment(M,4),V(3)
    REAL(DP) :: R(SIZE(X)),Physical(SIZE(X)),GS(NF+4),Z(4*M,NF+4),Absolute(NF+4)
    REAL(DP) :: D,I(3),KD(3,3),KI(3,3,3),Velocity(3),G
    INTEGER :: J,L,Q,S,F,K
    K=NF+4
    C0=0.0_DP
    C1=0.0_DP
    Source=0.0_DP
    Absolute=0.0_DP
    Z=0.0_DP
    CALL ComputeCollisionResidual_OrderV_Point(M,NF,Active,W,U,D0,I0,V0,X,Z,0.0_DP,C0,C1,Source, &
      R,Physical,Scale,GS,Absolute,Increment)
    H=Physical(4*M+1:4*M+K)
    V=Retain(V0)
    DO J=1,SIZE(Active)
      L=Active(J)-NF-1
      IF( L>=1 .AND. L<=3 ) V(L)=V(L)-X(4*M+J)
    END DO
    Velocity=RetainedValue(V)
    JH=0.0_DP
    DO Q=1,M
      D=RetainedValue(D0(Q)+X(Q))
      DO J=1,3
        I(J)=RetainedValue(I0(Q,J)+X(J*M+Q))
      END DO
      CALL ComputeCollisionClosureDerivatives_OrderV_Point(D,I,KD,KI)
      S=(Q-1)/(M/(2*NF))+1
      F=(S+1)/2
      G=MERGE(1.0_DP,-1.0_DP,MOD(S,2)==1)
      JH(F,Q)=G*W(Q)
      JH(NF+1,Q)=U(Q)
      JH(NF+2:K,Q)=U(Q)*(Velocity+MATMUL(KD,Velocity))
      DO J=1,3
        JH(F,J*M+Q)=G*W(Q)*Velocity(J)
        JH(NF+1,J*M+Q)=2.0_DP*U(Q)*Velocity(J)
        JH(NF+2:K,J*M+Q)=U(Q)*MATMUL(KI(:,:,J),Velocity)
        JH(NF+1+J,J*M+Q)=JH(NF+1+J,J*M+Q)+U(Q)
      END DO
    END DO
  END SUBROUTINE ComputeConservationGeometry

  PURE SUBROUTINE SelectConservationPivots(K,N,J,Capacity,Pivots,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Column-pivoted Householder QR; only the K pivot indices are required.
    INTEGER, INTENT(in) :: K,N
    REAL(DP), INTENT(in) :: J(K,N),Capacity(N)
    INTEGER, INTENT(out) :: Pivots(K),Info
    REAL(DP) :: A(K,N),Row(K),Norms(N),Vector(K),Swap(K),Alpha,Beta,Norm
    INTEGER :: Order(N),L,Q,P,Temp
    DO Q=1,N
      A(:,Q)=J(:,Q)*MAX(Capacity(Q),TINY(1.0_DP))
      Order(Q)=Q
    END DO
    Row=MAXVAL(ABS(A),DIM=2)
    Info=1
    IF( ANY(Row<=0.0_DP) ) RETURN
    DO Q=1,N
      A(:,Q)=A(:,Q)/Row
    END DO
    DO L=1,K
      DO Q=L,N
        Norms(Q)=SUM(A(L:K,Q)**2)
      END DO
      P=L-1+MAXLOC(Norms(L:N),DIM=1)
      Swap=A(:,L)
      A(:,L)=A(:,P)
      A(:,P)=Swap
      Temp=Order(L)
      Order(L)=Order(P)
      Order(P)=Temp
      Norm=SQRT(SUM(A(L:K,L)**2))
      IF( Norm<=TINY(1.0_DP) ) RETURN
      Alpha=-SIGN(Norm,A(L,L))
      Vector=0.0_DP
      Vector(L:K)=A(L:K,L)
      Vector(L)=Vector(L)-Alpha
      Beta=2.0_DP/SUM(Vector(L:K)**2)
      DO Q=L+1,N
        A(L:K,Q)=A(L:K,Q)-Beta*SUM(Vector(L:K)*A(L:K,Q))*Vector(L:K)
      END DO
      A(L,L)=Alpha
      A(L+1:K,L)=0.0_DP
      Pivots(L)=Order(L)
    END DO
    Info=0
  END SUBROUTINE SelectConservationPivots

  PURE SUBROUTINE SolveSmall(K,A,B,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: K
    REAL(DP), INTENT(inout) :: A(K,K),B(K)
    INTEGER, INTENT(out) :: Info
    REAL(DP) :: Row(K),Temp,Factor
    INTEGER :: J,L,P
    Info=1
    DO J=1,K
      P=J-1+MAXLOC(ABS(A(J:K,J)),DIM=1)
      IF( ABS(A(P,J))<=TINY(1.0_DP) ) RETURN
      Row=A(J,:)
      A(J,:)=A(P,:)
      A(P,:)=Row
      Temp=B(J)
      B(J)=B(P)
      B(P)=Temp
      DO L=J+1,K
        Factor=A(L,J)/A(J,J)
        A(L,J+1:K)=A(L,J+1:K)-Factor*A(J,J+1:K)
        B(L)=B(L)-Factor*B(J)
      END DO
    END DO
    DO J=K,1,-1
      B(J)=(B(J)-SUM(A(J,J+1:K)*B(J+1:K)))/A(J,J)
    END DO
    Info=0
  END SUBROUTINE SolveSmall

  PURE SUBROUTINE ReconstructConservativeRadiation(M,NF,Active,W,U,D0,I0,V0,Cap,X,Info,RetainedPrecision_Option)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,NF,Active(:)
    REAL(DP), INTENT(in) :: W(M),U(M),D0(M),I0(M,3),V0(3),Cap
    TYPE(RetainedReal), INTENT(inout) :: X(:)
    INTEGER, INTENT(out) :: Info
    LOGICAL, INTENT(in), OPTIONAL :: RetainedPrecision_Option
    REAL(DP) :: Roundoff
    REAL(DP) :: H(NF+4),Scale(NF+4),JH(NF+4,4*M),Capacity(4*M)
    REAL(DP) :: A(NF+4,NF+4),B(NF+4),Row(NF+4),D,I(3),Slack
    INTEGER :: Pivots(NF+4),Pass,Q,J,K
    K=NF+4
    Roundoff=32.0_DP*EPSILON(1.0_DP)
    IF( PRESENT(RetainedPrecision_Option) )THEN
      IF( RetainedPrecision_Option ) Roundoff=Roundoff*EPSILON(1.0_DP)
    END IF
    Info=1
    DO Pass=0,10
      IF( .NOT. CollisionStateAdmissible_OrderV_Point(M,D0,I0,X,Cap) ) RETURN
      CALL ComputeConservationGeometry(M,NF,Active,W,U,D0,I0,V0,X,H,Scale,JH)
      IF( ALL(ABS(H)<=MAX(Roundoff*Scale,1.0e-300_DP)) )THEN
        Info=0
        RETURN
      END IF
      IF( Pass==10 ) EXIT
      DO Q=1,M
        D=RetainedValue(D0(Q)+X(Q))
        DO J=1,3
          I(J)=RetainedValue(I0(Q,J)+X(J*M+Q))
        END DO
        Slack=MAX(D-Norm3(I),0.0_DP)
        Capacity(Q)=MIN(D,Cap-D,Slack)
        DO J=1,3
          Capacity(J*M+Q)=Slack
        END DO
      END DO
      CALL SelectConservationPivots(K,4*M,JH,Capacity,Pivots,Info)
      IF( Info/=0 ) RETURN
      A=JH(:,Pivots)
      Row=MAXVAL(ABS(A),DIM=2)
      IF( ANY(Row<=0.0_DP) )THEN
        Info=1
        RETURN
      END IF
      DO J=1,K
        A(:,J)=A(:,J)/Row
      END DO
      B=-H/Row
      CALL SolveSmall(K,A,B,Info)
      IF( Info/=0 ) RETURN
      DO J=1,K
        X(Pivots(J))=X(Pivots(J))+B(J)
      END DO
    END DO
    Info=1
  END SUBROUTINE ReconstructConservativeRadiation

  PURE REAL(DP) FUNCTION Norm3(X) RESULT(Norm)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: X(3)
    REAL(DP) :: Scale
    Scale=MAXVAL(ABS(X))
    Norm=0.0_DP
    IF( Scale>0.0_DP ) Norm=Scale*SQRT(SUM((X/Scale)**2))
  END FUNCTION Norm3

  PURE SUBROUTINE ProjectRadiationCone(D,I,Cap,PD,PI,Jac)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Euclidean projection on the capped Lorentz cone. For radial flux
    ! this reduces to triangle_project in the Python reference.
    REAL(DP), INTENT(in) :: D,I(3),Cap
    REAL(DP), INTENT(out) :: PD,PI(3),Jac(4,4)
    REAL(DP) :: R,T,Unit(3),Radial,Tangent
    INTEGER :: J,L
    R=Norm3(I)
    Jac=0.0_DP
    IF( D>=0.0_DP .AND. D<=Cap .AND. R<=D )THEN
      PD=D
      PI=I
      DO J=1,4
        Jac(J,J)=1.0_DP
      END DO
      RETURN
    END IF
    IF( R<=Cap .AND. D>=Cap )THEN
      PD=Cap
      PI=I
      DO J=2,4
        Jac(J,J)=1.0_DP
      END DO
      RETURN
    END IF
    T=MAX(0.0_DP,MIN(Cap,0.5_DP*(D+R)))
    PD=T
    PI=0.0_DP
    IF( R==0.0_DP ) RETURN
    Unit=I/R
    PI=T*Unit
    Radial=0.0_DP
    IF( T>0.0_DP .AND. T<Cap )THEN
      Jac(1,1)=0.5_DP
      Jac(1,2:4)=0.5_DP*Unit
      Jac(2:4,1)=0.5_DP*Unit
      Radial=0.5_DP
    END IF
    Tangent=T/R
    DO J=1,3
      DO L=1,3
        Jac(J+1,L+1)=(Radial-Tangent)*Unit(J)*Unit(L)
      END DO
      Jac(J+1,J+1)=Jac(J+1,J+1)+Tangent
    END DO
  END SUBROUTINE ProjectRadiationCone

  PURE SUBROUTINE SolveSmallSymmetric(K,A,B,X)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    ! Pseudoinverse of the positive semidefinite dual Hessian. A rank
    ! deficient active face is legitimate, not a failed radiation solve.
    ! Jacobi rotations avoid normal equations and are device callable.
    INTEGER, INTENT(in) :: K
    REAL(DP), INTENT(in) :: A(K,K),B(K)
    REAL(DP), INTENT(out) :: X(K)
    REAL(DP) :: D(K,K),V(K,K),Scale,Off,Tau,T,C,S,DPQ,DLeft,DRight,VP,VQ,Cutoff
    INTEGER :: Sweep,P,Q,J
    Scale=MAXVAL(ABS(A))
    X=0.0_DP
    IF( Scale<=TINY(1.0_DP) ) RETURN
    D=A/Scale
    V=0.0_DP
    DO J=1,K
      V(J,J)=1.0_DP
    END DO
    DO Sweep=1,50
      Off=0.0_DP
      DO P=1,K-1
      DO Q=P+1,K
        DPQ=D(P,Q)
        Off=MAX(Off,ABS(DPQ))
        IF( ABS(DPQ)<=EPSILON(1.0_DP)*MAX(ABS(D(P,P)),ABS(D(Q,Q)))) CYCLE
        Tau=(D(Q,Q)-D(P,P))/(2.0_DP*DPQ)
        T=SIGN(1.0_DP,Tau)/(ABS(Tau)+SQRT(1.0_DP+Tau*Tau))
        C=1.0_DP/SQRT(1.0_DP+T*T)
        S=T*C
        D(P,P)=D(P,P)-T*DPQ
        D(Q,Q)=D(Q,Q)+T*DPQ
        D(P,Q)=0.0_DP
        D(Q,P)=0.0_DP
        DO J=1,K
          IF( J/=P .AND. J/=Q )THEN
            DLeft=D(J,P)
            DRight=D(J,Q)
            D(J,P)=C*DLeft-S*DRight
            D(P,J)=D(J,P)
            D(J,Q)=S*DLeft+C*DRight
            D(Q,J)=D(J,Q)
          END IF
          VP=V(J,P)
          VQ=V(J,Q)
          V(J,P)=C*VP-S*VQ
          V(J,Q)=S*VP+C*VQ
        END DO
      END DO
      END DO
      IF( Off<=EPSILON(1.0_DP) ) EXIT
    END DO
    Cutoff=EPSILON(1.0_DP)*MAXVAL(ABS(D))
    DO J=1,K
      IF( D(J,J)>Cutoff ) X=X+V(:,J)*(SUM(V(:,J)*B)/Scale/D(J,J))
    END DO
  END SUBROUTINE SolveSmallSymmetric

  PURE SUBROUTINE EvaluateProjection(M,K,NA,Origin,Target,Scale,E,B,Lambda,Y,Defect,Q,Obj,Endpoint)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,K,NA
    REAL(DP), INTENT(in) :: Origin(4*M),Target(4*M+NA),Scale(4*M+NA),E(K,4*M+NA),B(K),Lambda(K)
    REAL(DP), INTENT(out) :: Y(4*M+NA),Defect(K),Q(K,K),Obj,Endpoint(4*M)
    REAL(DP) :: Z(4*M+NA),D,I(3),Jac(4,4),LocalE(K,4)
    INTEGER :: P,L,Index(4),N
    N=4*M
    Z=Target-MATMUL(TRANSPOSE(E),Lambda)
    Y=Z
    Q=MATMUL(E(:,N+1:),TRANSPOSE(E(:,N+1:)))
    DO P=1,M
      Index=[P,M+P,2*M+P,3*M+P]
      CALL ProjectRadiationCone(Origin(P)+Z(P),Origin(Index(2:4))+Z(Index(2:4)), &
        1.0_DP/Scale(P),D,I,Jac)
      Endpoint(P)=D*Scale(P)
      Endpoint(Index(2:4))=I*Scale(Index(2:4))
      Y(P)=D-Origin(P)
      Y(Index(2:4))=I-Origin(Index(2:4))
      LocalE=E(:,Index)
      Q=Q+MATMUL(MATMUL(LocalE,Jac),TRANSPOSE(LocalE))
    END DO
    Defect=MATMUL(E,Y)-B
    Obj=-(0.5_DP*SUM((Y-Target)**2)+SUM(Lambda*Defect))
  END SUBROUTINE EvaluateProjection

  PURE SUBROUTINE ProjectConservativeDirection(M,K,NA,D,I,Scale,JH,H,Proposal,Step,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,K,NA
    TYPE(RetainedReal), INTENT(in) :: D(M),I(M,3)
    REAL(DP), INTENT(in) :: Scale(4*M+NA),JH(K,4*M+NA),H(K),Proposal(4*M+NA)
    TYPE(RetainedReal), INTENT(out) :: Step(4*M+NA)
    INTEGER, INTENT(out) :: Info
    REAL(DP) :: E(K,4*M+NA),Row(K),B(K),Lambda(K),Delta(K),Target(4*M+NA),Origin(4*M)
    REAL(DP) :: Y(4*M+NA),Defect(K),Q(K,K),Obj,Endpoint(4*M),TrialY(4*M+NA)
    REAL(DP) :: TrialDefect(K),TrialQ(K,K),TrialObj,TrialEnd(4*M),Fraction,Descent
    INTEGER :: J,Iteration,Cut
    LOGICAL :: Accepted
    Info=1
    Step=0.0_DP
    IF( ANY(Scale<=0.0_DP) ) RETURN
    Origin(1:M)=RetainedValue(D)/Scale(1:M)
    DO J=1,3
      Origin(J*M+1:(J+1)*M)=RetainedValue(I(:,J))/Scale(J*M+1:(J+1)*M)
    END DO
    Target=Proposal/Scale
    DO J=1,4*M+NA
      E(:,J)=JH(:,J)*Scale(J)
    END DO
    Row=MAX(MAXVAL(ABS(E),DIM=2),TINY(1.0_DP))
    DO J=1,4*M+NA
      E(:,J)=E(:,J)/Row
    END DO
    B=-H/Row
    Lambda=0.0_DP
    DO Iteration=0,80
      CALL EvaluateProjection(M,K,NA,Origin,Target,Scale,E,B,Lambda,Y,Defect,Q,Obj,Endpoint)
      IF( MAXVAL(ABS(Defect)/MAX(1.0_DP,ABS(B)))<=2.0e-13_DP )THEN
        Step(1:M)=Retain(Endpoint(1:M))-D
        DO J=1,3
          Step(J*M+1:(J+1)*M)=Retain(Endpoint(J*M+1:(J+1)*M))-I(:,J)
        END DO
        Step(4*M+1:)=Retain(Scale(4*M+1:)*Y(4*M+1:))
        CALL RepairConservativeDirection(M,K,NA,D,I,Scale,JH,H,Step,Info)
        RETURN
      END IF
      IF( Iteration==80 ) EXIT
      CALL SolveSmallSymmetric(K,Q,Defect,Delta)
      Descent=SUM(Defect*Delta)
      IF( .NOT. IEEE_IS_FINITE(Descent) .OR. Descent<=0.0_DP )THEN
        Delta=Defect/MAX(MAXVAL(SUM(ABS(Q),DIM=1)),1.0_DP)
        Descent=SUM(Defect*Delta)
      END IF
      Accepted=.FALSE.
      DO Cut=0,23
        Fraction=2.0_DP**(-Cut)
        CALL EvaluateProjection(M,K,NA,Origin,Target,Scale,E,B,Lambda+Fraction*Delta, &
          TrialY,TrialDefect,TrialQ,TrialObj,TrialEnd)
        IF( TrialObj<=Obj-1.0e-4_DP*Fraction*Descent .OR. &
            SUM(TrialDefect**2)<(1.0_DP-1.0e-4_DP*Fraction)**2*SUM(Defect**2) )THEN
          Lambda=Lambda+Fraction*Delta
          Accepted=.TRUE.
          EXIT
        END IF
      END DO
      IF( .NOT. Accepted ) RETURN
    END DO
  END SUBROUTINE ProjectConservativeDirection

  PURE SUBROUTINE RepairConservativeDirection(M,K,NA,D,I,Scale,JH,H,Step,Info)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: M,K,NA
    TYPE(RetainedReal), INTENT(in) :: D(M),I(M,3)
    REAL(DP), INTENT(in) :: Scale(4*M+NA),JH(K,4*M+NA),H(K)
    TYPE(RetainedReal), INTENT(inout) :: Step(4*M+NA)
    INTEGER, INTENT(out) :: Info
    TYPE(RetainedReal) :: EndD,EndI(3),TargetD,TargetI(3),Defect
    REAL(DP) :: Capacity(4*M),A(K,K),B(K),Row(K),INorm,DValue,Slack,Change
    INTEGER :: Pass,P,J,L,Pivots(K),Index(3)
    LOGICAL :: Feasible,Zero
    Info=1
    DO Pass=1,4
      DO P=1,M
        Index=[M+P,2*M+P,3*M+P]
        EndD=D(P)+Step(P)
        EndI=I(P,:)+Step(Index)
        TargetD=EndD
        IF( RetainedValue(EndD)<=0.0_DP ) TargetD=0.0_DP
        IF( RetainedValue(1.0_DP-EndD)<0.0_DP ) TargetD=1.0_DP
        INorm=Norm3(RetainedValue(EndI))
        TargetI=EndI
        IF( INorm>RetainedValue(TargetD) ) TargetI=EndI*(TargetD/INorm)
        Change=MAX(ABS(RetainedValue(TargetD-EndD)),MAXVAL(ABS(RetainedValue(TargetI-EndI))))/Scale(P)
        IF( Change>1.0e-11_DP ) RETURN
        Step(P)=TargetD-D(P)
        Step(Index)=TargetI-I(P,:)
        DValue=RetainedValue(TargetD)
        Slack=MAX(MIN(DValue,1.0_DP-DValue,DValue-Norm3(RetainedValue(TargetI))),0.0_DP)
        Capacity(P)=Slack
        Capacity(Index)=Slack
      END DO
      CALL SelectConservationPivots(K,4*M,JH(:,1:4*M),Capacity,Pivots,Info)
      IF( Info/=0 ) RETURN
      A=JH(:,Pivots)
      Row=MAXVAL(ABS(A),DIM=2)
      DO J=1,K
        Defect=Retain(H(J))
        DO L=1,4*M+NA
          Defect=Defect+JH(J,L)*Step(L)
        END DO
        B(J)=-RetainedValue(Defect)/Row(J)
        A(J,:)=A(J,:)/Row(J)
      END DO
      CALL SolveSmall(K,A,B,Info)
      IF( Info/=0 ) RETURN
      DO J=1,K
        Step(Pivots(J))=Step(Pivots(J))+B(J)
      END DO
      Feasible=.TRUE.
      Zero=.FALSE.
      DO P=1,M
        Index=[M+P,2*M+P,3*M+P]
        EndD=D(P)+Step(P)
        EndI=I(P,:)+Step(Index)
        DValue=RetainedValue(EndD)
        Feasible=Feasible .AND. DValue>=0.0_DP .AND. RetainedValue(1.0_DP-EndD)>=0.0_DP &
          .AND. DValue>=Norm3(RetainedValue(EndI))
        Zero=Zero .OR. DValue==0.0_DP
      END DO
      IF( Feasible )THEN
        IF( Zero ) Step=Step*0.99_DP
        Info=0
        RETURN
      END IF
    END DO
    Info=1
  END SUBROUTINE RepairConservativeDirection

END MODULE TwoMoment_RadiationStateRecoveryModule_OrderV
