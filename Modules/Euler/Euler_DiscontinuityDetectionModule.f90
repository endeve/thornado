MODULE Euler_DiscontinuityDetectionModule

  USE KindModule,                     ONLY: &
    DP,   &
    Zero, &
    One,  &
    Two,  &
    SqrtTiny
  USE ProgramHeaderModule,            ONLY: &
    nDOFX,  &
    nDimsX, &
    nNodes
  USE ReferenceElementModuleX,        ONLY: &
    NodeNumberTableX, &
    WeightsX_q,       &
    NodesX1,          &
    NodesX2,          &
    NodesX3
  USE MeshModule,                     ONLY: &
    MeshX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, &
    L_X2, &
    L_X3
  USE GeometryFieldsModule,           ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_SqrtGm
  USE FluidFieldsModule,              ONLY: &
    nCF,     &
    iCF_D,   &
    iCF_S1,  &
    iCF_S2,  &
    iCF_S3,  &
    iCF_E,   &
    iCF_Ne,  &
    nPF,     &
    iPF_D,   &
    iPF_V1,  &
    iPF_V2,  &
    iPF_V3,  &
    iPF_E,   &
    iPF_Ne,  &
    iAF_P,   &
    iDF_TCI, &
    iDF_Sh
  USE Euler_UtilitiesModule,          ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule,          ONLY: &
    ComputePressureFromPrimitive
  USE TimersModule_Euler,             ONLY: &
    TimersStart_Euler,                 &
    TimersStop_Euler,                  &
    Timer_Euler_TroubledCellIndicator, &
    Timer_Euler_ShockDetector

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTroubledCellIndicator_Euler
  PUBLIC :: FinalizeTroubledCellIndicator_Euler
  PUBLIC :: DetectTroubledCells_Euler
  PUBLIC :: DetectShocks_Euler

  ! --- For troubled-cell indicator ---
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)

  LOGICAL  :: UseTroubledCellIndicator
  REAL(DP) :: LimiterThreshold


CONTAINS


  SUBROUTINE InitializeTroubledCellIndicator_Euler &
    ( UseTroubledCellIndicator_Option, LimiterThreshold_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UseTroubledCellIndicator_Option
    REAL(DP), INTENT(in), OPTIONAL :: LimiterThreshold_Option

    INTEGER  :: iNode, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: jNode, jNodeX1, jNodeX2, jNodeX3
    REAL(DP) :: WeightX

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    LimiterThreshold = 0.03_DP * 2.0_DP**( nNodes - 2 )
    IF( PRESENT( LimiterThreshold_Option ) ) &
      LimiterThreshold = LimiterThreshold_Option

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

        WeightX = WeightsX_q(iNode)

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

  END SUBROUTINE InitializeTroubledCellIndicator_Euler


  SUBROUTINE FinalizeTroubledCellIndicator_Euler

    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeTroubledCellIndicator_Euler


  SUBROUTINE DetectTroubledCells_Euler( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: U_K (0:2*nDimsX,nCF)
    REAL(DP) :: U_K0(0:2*nDimsX,nCF)
    REAL(DP) :: Y_K (0:2*nDimsX)
    REAL(DP) :: Y_K0(0:2*nDimsX)
    REAL(DP) :: Y(1:nDOFX, &
                  iX_B1(1):iX_E1(1), &
                  iX_B1(2):iX_E1(2), &
                  iX_B1(3):iX_E1(3))

    IF( .NOT. UseTroubledCellIndicator )THEN

      D(:,:,:,:,iDF_TCI) = 1.1_DP * LimiterThreshold
      RETURN

    END IF

    D(:,:,:,:,iDF_TCI) = Zero

    Y = U(:,:,:,:,iCF_Ne) / MAX( U(:,:,:,:,iCF_D), SqrtTiny )

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Compute Cell Averages ---
      ! --- in Target Cell and Neighbors in X1 Direction -------------

      DO iCF = 1, nCF

        U_K (0,iCF) &
          = DOT_PRODUCT( WeightsX_q,    U(:,iX1,iX2,iX3,iCF) )

        U_K (1,iCF) &
          = DOT_PRODUCT( WeightsX_q,    U(:,iX1-1,iX2,iX3,iCF) )

        U_K0(1,iCF) &
          = DOT_PRODUCT( WeightsX_X1_P, U(:,iX1-1,iX2,iX3,iCF) )

        U_K (2,iCF) &
          = DOT_PRODUCT( WeightsX_q,    U(:,iX1+1,iX2,iX3,iCF) )

        U_K0(2,iCF) &
          = DOT_PRODUCT( WeightsX_X1_N, U(:,iX1+1,iX2,iX3,iCF) )

      END DO

      Y_K (0) &
        = DOT_PRODUCT( WeightsX_q,    Y(:,iX1  ,iX2,iX3) )
      Y_K (1) &
        = DOT_PRODUCT( WeightsX_q,    Y(:,iX1-1,iX2,iX3) )
      Y_K0(1) &
        = DOT_PRODUCT( WeightsX_X1_P, Y(:,iX1-1,iX2,iX3) )
      Y_K (2) &
        = DOT_PRODUCT( WeightsX_q,    Y(:,iX1+1,iX2,iX3) )
      Y_K0(2) &
        = DOT_PRODUCT( WeightsX_X1_N, Y(:,iX1+1,iX2,iX3) )

      ! --- Compute Cell Averages ---
      ! --- in Neighbors in X2 Direction -------------

      IF( nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          U_K (3,iCF) &
            = DOT_PRODUCT( WeightsX_q,    U(:,iX1,iX2-1,iX3,iCF) )

          U_K0(3,iCF) &
            = DOT_PRODUCT( WeightsX_X2_P, U(:,iX1,iX2-1,iX3,iCF) )

          U_K (4,iCF) &
            = DOT_PRODUCT( WeightsX_q,    U(:,iX1,iX2+1,iX3,iCF) )

          U_K0(4,iCF) &
            = DOT_PRODUCT( WeightsX_X2_N, U(:,iX1,iX2+1,iX3,iCF) )

        END DO

        Y_K (3) &
          = DOT_PRODUCT( WeightsX_q,    Y(:,iX1,iX2-1,iX3) )
        Y_K0(3) &
          = DOT_PRODUCT( WeightsX_X2_P, Y(:,iX1,iX2-1,iX3) )
        Y_K (4) &
          = DOT_PRODUCT( WeightsX_q,    Y(:,iX1,iX2+1,iX3) )
        Y_K0(4) &
          = DOT_PRODUCT( WeightsX_X2_N, Y(:,iX1,iX2+1,iX3) )

      END IF

      ! --- Compute Cell Volumes and Cell Averages ---
      ! --- in Neighbors in X3 Direction -------------

      IF( nDimsX .GT. 2 )THEN

        DO iCF = 1, nCF

          U_K (5,iCF) &
            = DOT_PRODUCT( WeightsX_q,    U(:,iX1,iX2,iX3-1,iCF) )

          U_K0(5,iCF) &
            = DOT_PRODUCT( WeightsX_X3_P, U(:,iX1,iX2,iX3-1,iCF) )

          U_K (6,iCF) &
            = DOT_PRODUCT( WeightsX_q,    U(:,iX1,iX2,iX3+1,iCF) )

          U_K0(6,iCF) &
            = DOT_PRODUCT( WeightsX_X3_N, U(:,iX1,iX2,iX3+1,iCF) )

        END DO

        Y_K (5) &
          = DOT_PRODUCT( WeightsX_q,    Y(:,iX1,iX2,iX3-1) )
        Y_K0(5) &
          = DOT_PRODUCT( WeightsX_X3_P, Y(:,iX1,iX2,iX3-1) )
        Y_K (6) &
          = DOT_PRODUCT( WeightsX_q,    Y(:,iX1,iX2,iX3+1) )
        Y_K0(6) &
          = DOT_PRODUCT( WeightsX_X2_N, Y(:,iX1,iX2,iX3+1) )

      END IF

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(:,iX1,iX2,iX3,iDF_TCI) &
        = SUM( ABS( U_K(0,iCF_D) - U_K0(1:2*nDimsX,iCF_D) ) ) &
            / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_D) ) )

      ! --- Use Conserved Energy  to Detect Troubled Cell ---

      D(:,iX1,iX2,iX3,iDF_TCI) &
        = MAX( MAXVAL(D(:,iX1,iX2,iX3,iDF_TCI) ), &
               SUM( ABS( U_K(0,iCF_E) - U_K0(1:2*nDimsX,iCF_E) ) ) &
                 / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_E) ) ) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(:,iX1,iX2,iX3,iDF_TCI) &
        = MAX( MAXVAL( D(:,iX1,iX2,iX3,iDF_TCI) ), &
               SUM( 1.0d2 * ABS( Y_K(0) - Y_K0(1:2*nDimsX) ) ) &
                 / MAX( MAXVAL( ABS( Y_K(0:2*nDimsX) ) ), SqrtTiny ) )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

  END SUBROUTINE DetectTroubledCells_Euler


  SUBROUTINE DetectShocks_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: V_K
    REAL(DP) :: uPF_K(nPF)
    REAL(DP) :: uCF_K(nCF)
    REAL(DP) :: uGF_K(nGF)
    REAL(DP) :: P_K(2), VX_K(2)
    REAL(DP) :: GradP(nDimsX), DivV
    REAL(DP) :: dX1, dX2, dX3

    ! --- Shock detector, adapted from
    !     Fryxell et al., (2000), ApJS, 131, 273 ---

    CALL TimersStart_Euler( Timer_Euler_ShockDetector )

    D(:,:,:,:,iDF_Sh) = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      GradP = Zero

      ! --- Lower neighbor in X1 direction ---

      V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1-1,iX2,iX3,iGF_SqrtGm) )

      DO iCF = 1, nCF

        uCF_K(iCF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1-1,iX2,iX3,iCF) ) / V_K

      END DO

      DO iGF = 1, nGF

        uGF_K(iGF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                  * G(:,iX1-1,iX2,iX3,iGF) ) / V_K

      END DO

      CALL ComputePrimitive_Euler &
           ( uCF_K(iCF_D ), &
             uCF_K(iCF_S1), &
             uCF_K(iCF_S2), &
             uCF_K(iCF_S3), &
             uCF_K(iCF_E ), &
             uCF_K(iCF_Ne), &
             uPF_K(iPF_D ), &
             uPF_K(iPF_V1), &
             uPF_K(iPF_V2), &
             uPF_K(iPF_V3), &
             uPF_K(iPF_E ), &
             uPF_K(iPF_Ne), &
             uGF_K(iGF_Gm_dd_11), &
             uGF_K(iGF_Gm_dd_22), &
             uGF_K(iGF_Gm_dd_33) )

      VX_K(1) = uGF_K(iGF_Alpha) * uPF_K(iPF_V1) + uGF_K(iGF_Beta_1)

      CALL ComputePressureFromPrimitive &
             ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(1) )

      ! --- Upper neighbor in X1 direction ---

      V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1+1,iX2,iX3,iGF_SqrtGm) )

      DO iCF = 1, nCF

        uCF_K(iCF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                  * U(:,iX1+1,iX2,iX3,iCF) ) / V_K

      END DO

      DO iGF = 1, nGF

        uGF_K(iGF) &
          = DOT_PRODUCT &
              ( WeightsX_q, &
                G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                  * G(:,iX1+1,iX2,iX3,iGF) ) / V_K

      END DO

      CALL ComputePrimitive_Euler &
           ( uCF_K(iCF_D ), &
             uCF_K(iCF_S1), &
             uCF_K(iCF_S2), &
             uCF_K(iCF_S3), &
             uCF_K(iCF_E ), &
             uCF_K(iCF_Ne), &
             uPF_K(iPF_D ), &
             uPF_K(iPF_V1), &
             uPF_K(iPF_V2), &
             uPF_K(iPF_V3), &
             uPF_K(iPF_E ), &
             uPF_K(iPF_Ne), &
             uGF_K(iGF_Gm_dd_11), &
             uGF_K(iGF_Gm_dd_22), &
             uGF_K(iGF_Gm_dd_33) )

      VX_K(2) = uGF_K(iGF_Alpha) * uPF_K(iPF_V1) + uGF_K(iGF_Beta_1)

      CALL ComputePressureFromPrimitive &
             ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(2) )

      ! --- Compute pressure gradient and divergence of velocity (X1) ---

      GradP(1) = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

      DivV     = ( VX_K(2) - VX_K(1) ) / ( Two * dX1 )

      IF( nDimsX .GT. 1 )THEN

        ! --- Lower neighbor in X2 direction ---

        V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2-1,iX3,iGF_SqrtGm) )

        DO iCF = 1, nCF

          uCF_K(iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2-1,iX3,iCF) ) / V_K

        END DO

        DO iGF = 1, nGF

          uGF_K(iGF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                    * G(:,iX1,iX2-1,iX3,iGF) ) / V_K

        END DO

        CALL ComputePrimitive_Euler &
             ( uCF_K(iCF_D ), &
               uCF_K(iCF_S1), &
               uCF_K(iCF_S2), &
               uCF_K(iCF_S3), &
               uCF_K(iCF_E ), &
               uCF_K(iCF_Ne), &
               uPF_K(iPF_D ), &
               uPF_K(iPF_V1), &
               uPF_K(iPF_V2), &
               uPF_K(iPF_V3), &
               uPF_K(iPF_E ), &
               uPF_K(iPF_Ne), &
               uGF_K(iGF_Gm_dd_11), &
               uGF_K(iGF_Gm_dd_22), &
               uGF_K(iGF_Gm_dd_33) )

        VX_K(1) = uGF_K(iGF_Alpha) * uPF_K(iPF_V2) + uGF_K(iGF_Beta_2)

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(1) )

        ! --- Upper neighbor in X2 direction ---

        V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2+1,iX3,iGF_SqrtGm) )

        DO iCF = 1, nCF

          uCF_K(iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                    * U(:,iX1,iX2+1,iX3,iCF) ) / V_K

        END DO

        DO iGF = 1, nGF

          uGF_K(iGF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                    * G(:,iX1,iX2+1,iX3,iGF) ) / V_K

        END DO

        CALL ComputePrimitive_Euler &
             ( uCF_K(iCF_D ), &
               uCF_K(iCF_S1), &
               uCF_K(iCF_S2), &
               uCF_K(iCF_S3), &
               uCF_K(iCF_E ), &
               uCF_K(iCF_Ne), &
               uPF_K(iPF_D ), &
               uPF_K(iPF_V1), &
               uPF_K(iPF_V2), &
               uPF_K(iPF_V3), &
               uPF_K(iPF_E ), &
               uPF_K(iPF_Ne), &
               uGF_K(iGF_Gm_dd_11), &
               uGF_K(iGF_Gm_dd_22), &
               uGF_K(iGF_Gm_dd_33) )

        VX_K(2) = uGF_K(iGF_Alpha) * uPF_K(iPF_V2) + uGF_K(iGF_Beta_2)

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(2) )

        ! --- Compute pressure gradient and divergence of velocity (X2) ---

        GradP(2) = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

        DivV     = DivV + ( VX_K(2) - VX_K(1) ) / ( Two * dX2 )

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Lower neighbor in X3 direction ---

        V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3-1,iGF_SqrtGm) )

        DO iCF = 1, nCF

          uCF_K(iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3-1,iCF) ) / V_K

        END DO

        DO iGF = 1, nGF

          uGF_K(iGF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                    * G(:,iX1,iX2,iX3-1,iGF) ) / V_K

        END DO

        CALL ComputePrimitive_Euler &
             ( uCF_K(iCF_D ), &
               uCF_K(iCF_S1), &
               uCF_K(iCF_S2), &
               uCF_K(iCF_S3), &
               uCF_K(iCF_E ), &
               uCF_K(iCF_Ne), &
               uPF_K(iPF_D ), &
               uPF_K(iPF_V1), &
               uPF_K(iPF_V2), &
               uPF_K(iPF_V3), &
               uPF_K(iPF_E ), &
               uPF_K(iPF_Ne), &
               uGF_K(iGF_Gm_dd_11), &
               uGF_K(iGF_Gm_dd_22), &
               uGF_K(iGF_Gm_dd_33) )

        VX_K(1) = uGF_K(iGF_Alpha) * uPF_K(iPF_V3) + uGF_K(iGF_Beta_3)

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(1) )

        ! --- Upper neighbor in X3 direction ---

        V_K = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3+1,iGF_SqrtGm) )

        DO iCF = 1, nCF

          uCF_K(iCF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                    * U(:,iX1,iX2,iX3+1,iCF) ) / V_K

        END DO

        DO iGF = 1, nGF

          uGF_K(iGF) &
            = DOT_PRODUCT &
                ( WeightsX_q, &
                  G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                    * G(:,iX1,iX2,iX3+1,iGF) ) / V_K

        END DO

        CALL ComputePrimitive_Euler &
             ( uCF_K(iCF_D ), &
               uCF_K(iCF_S1), &
               uCF_K(iCF_S2), &
               uCF_K(iCF_S3), &
               uCF_K(iCF_E ), &
               uCF_K(iCF_Ne), &
               uPF_K(iPF_D ), &
               uPF_K(iPF_V1), &
               uPF_K(iPF_V2), &
               uPF_K(iPF_V3), &
               uPF_K(iPF_E ), &
               uPF_K(iPF_Ne), &
               uGF_K(iGF_Gm_dd_11), &
               uGF_K(iGF_Gm_dd_22), &
               uGF_K(iGF_Gm_dd_33) )

        VX_K(2) = uGF_K(iGF_Alpha) * uPF_K(iPF_V3) + uGF_K(iGF_Beta_3)

        CALL ComputePressureFromPrimitive &
               ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), P_K(2) )

        ! --- Compute pressure gradient and divergence of velocity (X3) ---

        GradP(3) = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

        DivV     = DivV + ( VX_K(2) - VX_K(1) ) / ( Two * dX3 )

      END IF

      IF( SQRT( SUM( GradP**2 ) ) .GT. 1.0_DP / 3.0_DP .AND. DivV .LT. Zero ) &
        D(:,iX1,iX2,iX3,iDF_Sh) = One

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ShockDetector )

  END SUBROUTINE DetectShocks_Euler


END MODULE Euler_DiscontinuityDetectionModule

