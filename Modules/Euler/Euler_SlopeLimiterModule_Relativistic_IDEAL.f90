MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_TroubledCellIndicator
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX, nNodes, nNodesX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q
  USE UtilitiesModule, ONLY: &
    MinModB, NodeNumberX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, P_X2, P_X3, IndPX_Q
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_E, &
    Shock
  USE Euler_BoundaryConditionsModule, ONLY: &
    Euler_ApplyBoundaryConditions
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL, ONLY: &
    Euler_ComputeCharacteristicDecomposition_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializeSlopeLimiter_Relativistic
  PUBLIC :: Euler_FinalizeSlopeLimiter_Relativistic
  PUBLIC :: Euler_ApplySlopeLimiter_Relativistic

  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseCharacteristicLimiting
  LOGICAL  :: UseConservativeCorrection
  LOGICAL  :: UseTroubledCellIndicator
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: LimiterThreshold
  REAL(DP) :: LimiterThresholdParameter
  REAL(DP) :: I_6x6(1:6,1:6)
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)


CONTAINS


  SUBROUTINE Euler_InitializeSlopeLimiter_Relativistic &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option,  &
      UseConservativeCorrection_Option, &
      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option

    INTEGER :: i
    LOGICAL :: Verbose

    BetaTVD = One
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    BetaTVB = Zero
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    SlopeTolerance = 1.0d-3
    IF( PRESENT( SlopeTolerance_Option ) ) &
      SlopeTolerance = SlopeTolerance_Option

    UseSlopeLimiter = .TRUE.
    IF( PRESENT( UseSlopeLimiter_Option ) ) &
      UseSlopeLimiter = UseSlopeLimiter_Option

    UseCharacteristicLimiting = .FALSE.
    IF( PRESENT( UseCharacteristicLimiting_Option ) ) &
      UseCharacteristicLimiting = UseCharacteristicLimiting_Option

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    LimiterThresholdParameter = 0.03_DP
    IF( PRESENT( LimiterThresholdParameter_Option ) ) &
      LimiterThresholdParameter = LimiterThresholdParameter_Option
    LimiterThreshold = LimiterThresholdParameter * 2.0_DP**( nNodes - 2 )

    UseConservativeCorrection = .TRUE.
    IF( PRESENT( UseConservativeCorrection_Option ) ) &
      UseConservativeCorrection = UseConservativeCorrection_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(A)') '    INFO: Euler_InitializeSlopeLimiter_Relativistic'
      WRITE(*,'(A)') '    -----------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVD: ' , &
        BetaTVD
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVB: ' , &
        BetaTVB
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'SlopeTolerance: ' , &
        SlopeTolerance
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'LimiterThreshold: ' , &
        LimiterThreshold
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseConservativeCorrection: ' , &
        UseConservativeCorrection
    END IF

    IF( UseTroubledCellIndicator ) &
      CALL InitializeTroubledCellIndicator

    I_6x6 = Zero
    DO i = 1, 6
      I_6x6(i,i) = One
    END DO

  END SUBROUTINE Euler_InitializeSlopeLimiter_Relativistic


  SUBROUTINE Euler_FinalizeSlopeLimiter_Relativistic

    IF( UseTroubledCellIndicator ) &
      CALL FinalizeTroubledCellIndicator

  END SUBROUTINE Euler_FinalizeSlopeLimiter_Relativistic


  SUBROUTINE Euler_ApplySlopeLimiter_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL  :: LimitedCell(nCF,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    LOGICAL  :: SuppressBC
    INTEGER  :: iX1, iX2, iX3, iGF, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: SlopeDifference(nCF)
    REAL(DP) :: G_K(nGF)
    REAL(DP) :: dU (nCF,nDimsX)
    REAL(DP) :: U_M(nCF,0:2*nDimsX,nDOFX)
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)
    REAL(DP) :: R_X2(nCF,nCF), invR_X2(nCF,nCF)
    REAL(DP) :: R_X3(nCF,nCF), invR_X3(nCF,nCF)
    REAL(DP) :: V_K(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL Euler_ApplyBoundaryConditions &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )
    CALL DetectTroubledCells &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )
    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( Shock(iX1,iX2,iX3) .LT. LimiterThreshold ) CYCLE

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      ! --- Cell Volume ---

      V_K(iX1,iX2,iX3) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_SqrtGm) )

      ! --- Cell Average of Conserved Fluid ---

      DO iCF = 1, nCF

        U_K(iCF,iX1,iX2,iX3) &
          = SUM( WeightsX_q * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                   * U(:,iX1,iX2,iX3,iCF) ) / V_K(iX1,iX2,iX3)

      END DO

      ! --- Cell Average of Geometry (Spatial Metric, Lapse Function,
      !     and Shift Vector) ---

      DO iGF = iGF_Gm_dd_11, iGF_Beta_3
        G_K(iGF) = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF) )
      END DO

      ! --- Map to Modal Representation ---

      DO iCF = 1, nCF

        CALL MapNodalToModal_Fluid( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

        ! --- Cell Average of Neighbors in X1 Direction ---

        U_M(iCF,1,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1-1,iX2,iX3,iCF) )

        U_M(iCF,2,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1+1,iX2,iX3,iCF) )

        IF( nDimsX .GT. 1 )THEN

          ! --- Cell Average of Neighbors in X2 Direction ---

          U_M(iCF,3,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2-1,iX3,iCF) )

          U_M(iCF,4,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2+1,iX3,iCF) )

        END IF

        IF( nDimsX .GT. 2 )THEN

          ! --- Cell Average of Neighbors in X3 Direction ---

          U_M(iCF,5,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2,iX3-1,iCF) )

          U_M(iCF,6,1) = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2,iX3+1,iCF) )

        END IF

      END DO

      IF( UseCharacteristicLimiting )THEN

        ! --- Compute Eigenvectors ---

        CALL Euler_ComputeCharacteristicDecomposition_Relativistic &
               ( 1, G_K, U_M(:,0,1), R_X1, invR_X1 )

        IF( nDimsX .GT. 1 )THEN

          CALL Euler_ComputeCharacteristicDecomposition_Relativistic &
                 ( 2, G_K, U_M(:,0,1), R_X2, invR_X2 )

        END IF

        IF( nDimsX .GT. 2 )THEN

          CALL Euler_ComputeCharacteristicDecomposition_Relativistic &
                 ( 3, G_K, U_M(:,0,1), R_X3, invR_X3 )

        END IF

      ELSE

        ! --- Componentwise Limiting ---

        R_X1 = I_6x6; invR_X1 = I_6x6
        R_X2 = I_6x6; invR_X2 = I_6x6
        R_X3 = I_6x6; invR_X3 = I_6x6

      END IF

      ! --- Compute Limited Slopes ---

      dU(:,1) &
        = MinModB &
            ( MATMUL( invR_X1, U_M(:,0,2) ), &
              BetaTVD * MATMUL( invR_X1, (U_M(:,0,1)-U_M(:,1,1)) ), &
              BetaTVD * MATMUL( invR_X1, (U_M(:,2,1)-U_M(:,0,1)) ), &
              dX1, BetaTVB )

      IF( nDimsX .GT. 1 )THEN

        dU(:,2) &
          = MinModB &
              ( MATMUL( invR_X2, U_M(:,0,3) ), &
                BetaTVD * MATMUL( invR_X2, (U_M(:,0,1)-U_M(:,3,1)) ), &
                BetaTVD * MATMUL( invR_X2, (U_M(:,4,1)-U_M(:,0,1)) ), &
                dX2, BetaTVB )

      END IF

      IF( nDimsX .GT. 2 )THEN

        dU(:,3) &
          = MinModB &
              ( MATMUL( invR_X3, U_M(:,0,4) ), &
                BetaTVD * MATMUL( invR_X3, (U_M(:,0,1)-U_M(:,5,1)) ), &
                BetaTVD * MATMUL( invR_X3, (U_M(:,6,1)-U_M(:,0,1)) ), &
                dX3, BetaTVB )

      END IF

      IF( UseCharacteristicLimiting )THEN

        ! --- Transform Back from Characteristic Variables ---

        dU(:,1) = MATMUL( R_X1, dU(:,1) )

        IF( nDimsX .GT. 1 )THEN

          dU(:,2) = MATMUL( R_X2, dU(:,2) )

        END IF

        IF( nDimsX .GT. 2 )THEN

          dU(:,3) = MATMUL( R_X3, dU(:,3) )

        END IF

      END IF

      ! --- Compare Limited Slopes to Original Slopes ---

      DO iCF = 1, nCF

        SlopeDifference(iCF) = ABS( U_M(iCF,0,2) - dU(iCF,1) )

        IF( nDimsX .GT. 1 )THEN

          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(iCF,0,3) - dU(iCF,2) ) )

        END IF

        IF( nDimsX .GT. 2 )THEN

          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(iCF,0,4) - dU(iCF,3) ) )

        END IF

      END DO

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      DO iCF = 1, nCF

        IF( SlopeDifference(iCF) .GT. SlopeTolerance * ABS( U_M(iCF,0,1) ) )THEN

          U_M(iCF,0,2:nDOFX) = Zero

          U_M(:,0,2) = dU(:,1)

          IF( nDimsX .GT. 1 ) U_M(iCF,0,3) = dU(iCF,2)

          IF( nDimsX .GT. 2 ) U_M(iCF,0,4) = dU(iCF,3)

          CALL MapModalToNodal_Fluid( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

          LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

        END IF

      END DO

    END DO
    END DO
    END DO

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, V_K, U, U_K, LimitedCell )

  END SUBROUTINE Euler_ApplySlopeLimiter_Relativistic


  SUBROUTINE InitializeTroubledCellIndicator

    INTEGER  :: iNode, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: jNode, jNodeX1, jNodeX2, jNodeX3
    REAL(DP) :: WeightX

    ALLOCATE( WeightsX_X1_P(nDOFX), WeightsX_X1_N(nDOFX) )
    ALLOCATE( WeightsX_X2_P(nDOFX), WeightsX_X2_N(nDOFX) )
    ALLOCATE( WeightsX_X3_P(nDOFX), WeightsX_X3_N(nDOFX) )

    ! --- Compute Weights for Extrapolating Neighbors into Target Cell ---

    DO jNode = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNode)
      jNodeX2 = NodeNumberTableX(2,jNode)
      jNodeX3 = NodeNumberTableX(3,jNode)

      WeightsX_X1_P(jNode) = Zero
      WeightsX_X1_N(jNode) = Zero
      WeightsX_X2_P(jNode) = Zero
      WeightsX_X2_N(jNode) = Zero
      WeightsX_X3_P(jNode) = Zero
      WeightsX_X3_N(jNode) = Zero

      DO iNode = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNode)
        iNodeX2 = NodeNumberTableX(2,iNode)
        iNodeX3 = NodeNumberTableX(3,iNode)

        WeightX = WeightsX1  (iNodeX1) &
                  * WeightsX2(iNodeX2) &
                  * WeightsX3(iNodeX3)

        WeightsX_X1_P(jNode) &
          = WeightsX_X1_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) + One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X1_N(jNode) &
          = WeightsX_X1_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) - One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_P(jNode) &
          = WeightsX_X2_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) + One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_N(jNode) &
          = WeightsX_X2_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) - One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X3_P(jNode) &
          = WeightsX_X3_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) + One ) )

        WeightsX_X3_N(jNode) &
          = WeightsX_X3_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) - One ) )

      END DO

    END DO

  END SUBROUTINE InitializeTroubledCellIndicator


  SUBROUTINE FinalizeTroubledCellIndicator

    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeTroubledCellIndicator


  SUBROUTINE DetectTroubledCells( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: V_K (0:2*nDimsX)
    REAL(DP) :: U_K (0:2*nDimsX,nCF)
    REAL(DP) :: U_K0(0:2*nDimsX,nCF)

    Shock = Zero

    IF( .NOT. UseTroubledCellIndicator )THEN

      Shock = 1.1_DP * LimiterThreshold
      RETURN

    END IF

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Compute Cell Volumes and Cell Averages ---------
      ! --- in Target Cell and Neighbors in X1 Direction ---

      V_K(0) = DOT_PRODUCT &
                 ( WeightsX_q, G(:,iX1,  iX2,iX3,iGF_SqrtGm) )

      V_K(1) = DOT_PRODUCT &
                 ( WeightsX_q, G(:,iX1-1,iX2,iX3,iGF_SqrtGm) )

      V_K(2) = DOT_PRODUCT &
                 ( WeightsX_q, G(:,iX1+1,iX2,iX3,iGF_SqrtGm) )

      DO iCF = 1, nCF

        U_K(0,iCF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1,iX2,iX3,iCF) ) / V_K(0)

        U_K(1,iCF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(1)

        U_K0(1,iCF) &
          = DOT_PRODUCT &
              ( WeightsX_X1_P, &
                G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(0)

        U_K(2,iCF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(2)

        U_K0(2,iCF) &
          = DOT_PRODUCT &
              ( WeightsX_X1_N, &
                G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(0)

      END DO

      ! --- Compute Cell Volumes and Cell Averages ---
      ! --- in Neighbors in X2 Direction -------------

      IF( nDimsX .GT. 1 )THEN

        V_K(3) = DOT_PRODUCT &
                   ( WeightsX_q, G(:,iX1,iX2-1,iX3,iGF_SqrtGm) )

        V_K(4) = DOT_PRODUCT &
                   ( WeightsX_q, G(:,iX1,iX2+1,iX3,iGF_SqrtGm) )

        DO iCF = 1, nCF

          U_K(3,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(3)

          U_K0(3,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_X2_P, &
                  G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(0)

          U_K(4,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(4)

          U_K0(4,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_X2_N, &
                  G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(0)

        END DO

      END IF

      ! --- Compute Cell Volumes and Cell Averages ---
      ! --- in Neighbors in X3 Direction -------------

      IF( nDimsX .GT. 2 )THEN

        V_K(5) = DOT_PRODUCT &
                   ( WeightsX_q, G(:,iX1,iX2,iX3-1,iGF_SqrtGm) )

        V_K(6) = DOT_PRODUCT &
                   ( WeightsX_q, G(:,iX1,iX2,iX3+1,iGF_SqrtGm) )

        DO iCF = 1, nCF

          U_K(5,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(5)

          U_K0(5,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_X3_P, &
                  G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(0)

          U_K(6,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3+1,iCF) ) / V_K(6)

          U_K0(6,iCF) &
            = DOT_PRODUCT &
                ( WeightsX_X3_N, &
                  G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3+1,iCF) ) / V_K(0)

        END DO

      END IF

      ! --- Use Conserved Density to Detect Troubled Cell ---

      Shock(iX1,iX2,iX3) &
        = SUM( ABS( U_K(0,iCF_D) - U_K0(1:2*nDimsX,iCF_D) ) ) &
            / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_D) ) )

      ! --- Use Conserved Energy  to Detect Troubled Cell ---

      Shock(iX1,iX2,iX3) &
        = MAX( Shock(iX1,iX2,iX3), &
               SUM( ABS( U_K(0,iCF_E) - U_K0(1:2*nDimsX,iCF_E) ) ) &
                 / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_E) ) ) )

    END DO
    END DO
    END DO

  END SUBROUTINE DetectTroubledCells


  SUBROUTINE ApplyConservativeCorrection &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, V_K, U, U_K, LimitedCell )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      V_K(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in)    :: &
      U_K(1:nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL, INTENT(in)     :: &
      LimitedCell(nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    INTEGER  :: iX1, iX2, iX3, iCF, iPol, iDimX
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: Correction
    REAL(DP) :: LegendreX(nDOFX,nDOFX)
    REAL(DP) :: U_M(nCF,0:2**nDimsX,nDOFX)

    IF( .NOT. UseConservativeCorrection ) RETURN

    DO iPol = 1, nDOFX ! Only need for iPol = 2,3,4 (FIXME)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

        LegendreX(iNodeX,iPol) &
          = P_X1  (IndPX_Q(1,iPol)) % P( MeshX(1) % Nodes(iNodeX1) ) &
            * P_X2(IndPX_Q(2,iPol)) % P( MeshX(2) % Nodes(iNodeX2) ) &
            * P_X3(IndPX_Q(3,iPol)) % P( MeshX(3) % Nodes(iNodeX3) )

      END DO
      END DO
      END DO

    END DO

    ! --- Applies a correction to the 0-th order ---
    ! --- mode to maintain the cell average.     ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iCF = 1, nCF

        IF( LimitedCell(iCF,iX1,iX2,iX3) )THEN

          CALL MapNodalToModal_Fluid( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

          Correction = Zero
          DO iDimX = 1, nDimsX

            Correction &
              = Correction &
                  + ( U_M(iCF,0,iDimX+1) &
                      * SUM( WeightsX_q(:) * LegendreX(:,iDimX+1) &
                             * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) ) &
                    / V_K(iX1,iX2,iX3)

          END DO

          U_M(iCF,0,1) = U_K(iCF,iX1,iX2,iX3) - Correction

          CALL MapModalToNodal_Fluid( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

        END IF

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyConservativeCorrection


END MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL
