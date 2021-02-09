MODULE Euler_DiscontinuityDetectionModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Third, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX, &
    nNodes
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX_q, &
    NodesX1, &
    NodesX2, &
    NodesX3
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, &
    L_X2, &
    L_X3
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    iDF_TCI, &
    iDF_Sh_X1, &
    iDF_Sh_X2, &
    iDF_Sh_X3
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE UtilitiesModule, ONLY: &
    IsCornerCell
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_DD_TCI, &
    Timer_Euler_DD_ShockDetector, &
    Timer_Euler_CopyIn, &
    Timer_Euler_CopyOut, &
    Timer_Euler_Permute
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

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

  LOGICAL,  PUBLIC :: UseTroubledCellIndicator
  REAL(DP), PUBLIC :: LimiterThreshold

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET( UseTroubledCellIndicator, LimiterThreshold )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE( UseTroubledCellIndicator, LimiterThreshold )
#endif


CONTAINS


  SUBROUTINE InitializeTroubledCellIndicator_Euler &
    ( UseTroubledCellIndicator_Option, &
      LimiterThresholdParameter_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UseTroubledCellIndicator_Option
    REAL(DP), INTENT(in), OPTIONAL :: LimiterThresholdParameter_Option

    INTEGER  :: iNode, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: jNode, jNodeX1, jNodeX2, jNodeX3
    REAL(DP) :: WeightX
    REAL(DP) :: LimiterThresholdParameter

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    LimiterThresholdParameter = 0.03_DP
    IF( PRESENT( LimiterThresholdParameter_Option ) ) &
      LimiterThresholdParameter = LimiterThresholdParameter_Option
    LimiterThreshold = LimiterThresholdParameter * Two**( nNodes - 2 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( UseTroubledCellIndicator, LimiterThreshold )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE   ( UseTroubledCellIndicator, LimiterThreshold )
#endif

    IF( .NOT. UseTroubledCellIndicator ) RETURN

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

    IF( .NOT. UseTroubledCellIndicator ) RETURN

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

    INTEGER  :: iNX, iX1, iX2, iX3, iCF
    REAL(DP) :: U_K (0:2*nDimsX,nCF)
    REAL(DP) :: U_K0(0:2*nDimsX,nCF)
    REAL(DP) :: Y_K (0:2*nDimsX)
    REAL(DP) :: Y_K0(0:2*nDimsX)
    REAL(DP) :: Y(1:nDOFX, &
                  iX_B1(1):iX_E1(1), &
                  iX_B1(2):iX_E1(2), &
                  iX_B1(3):iX_E1(3))

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B1, iX_E1, D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B1, iX_E1, D )
#endif

    IF( .NOT. UseTroubledCellIndicator )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B1, iX_E1, D )
#endif
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        D(iNX,iX1,iX2,iX3,iDF_TCI) = 1.1_DP * LimiterThreshold

      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    D ) &
    !$OMP MAP( release: iX_B1, iX_E1 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      D ) &
    !$ACC DELETE(       iX_B1, iX_E1 )
#endif

      RETURN

    END IF

    D(:,:,:,:,iDF_TCI) = Zero

    Y = U(:,:,:,:,iCF_Ne) / MAX( U(:,:,:,:,iCF_D), SqrtTiny )

    CALL TimersStart_Euler( Timer_Euler_DD_TCI )

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

    CALL TimersStop_Euler( Timer_Euler_DD_TCI )

  END SUBROUTINE DetectTroubledCells_Euler


  SUBROUTINE DetectShocks_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3
    REAL(DP) :: GradP, DivV
    INTEGER  :: nK_X, nCF_X, nGF_X

    REAL(DP) :: SqrtGm(1:nDOFX,   iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    ! 1: Gm11, 2: Gm22, 3: Gm33, 4: Alpha, 5: Beta1, 6: Beta2, 7: Beta3
    REAL(DP) :: G_X(1:nDOFX,1:7,  iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: GK(1:7,           iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: U_X(1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: UK(1:nCF,         iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: PK(1:nPF,         iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: VK(3,             iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: PrK(              iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: Vol(              iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    INTEGER :: iErr(              iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    INTEGER :: iX1arr(            iX_B1(1):iX_E1(1))
    INTEGER :: iX2arr(            iX_B1(2):iX_E1(2))
    INTEGER :: iX3arr(            iX_B1(3):iX_E1(3))

    CALL TimersStart_Euler( Timer_Euler_DD_ShockDetector )

    nK_X = PRODUCT( [ iX_E1(1) - iX_B1(1) + 1, &
                      iX_E1(2) - iX_B1(2) + 1, &
                      iX_E1(3) - iX_B1(3) + 1 ] )

    nCF_X = nCF * nK_X
    nGF_X = 7   * nK_X

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D ) &
    !$OMP MAP( alloc: SqrtGm, G_X, GK, U_X, UK, PK, VK, PrK, Vol, &
    !$OMP             iX1arr, iX2arr, iX3arr, iErr )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D ) &
    !$ACC CREATE(     SqrtGm, G_X, GK, U_X, UK, PK, VK, PrK, Vol, &
    !$ACC             iX1arr, iX2arr, iX3arr, iErr )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( iX_B1, iX_E1, iX1arr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
     DO iX1 = iX_B1(1), iX_E1(1)

       iX1arr(iX1) = iX1

     END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( iX_B1, iX_E1, iX2arr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
     DO iX2 = iX_B1(2), iX_E1(2)

       iX2arr(iX2) = iX2

     END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( iX_B1, iX_E1, iX3arr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
     DO iX3 = iX_B1(3), iX_E1(3)

       iX3arr(iX3) = iX3

     END DO

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B1, iX_E1, SqrtGm, G_X, U_X, G, U, D, &
    !$ACC          iX1arr, iX2arr, iX3arr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      D(iNX,iX1,iX2,iX3,iDF_Sh_X1) = Zero
      D(iNX,iX1,iX2,iX3,iDF_Sh_X2) = Zero
      D(iNX,iX1,iX2,iX3,iDF_Sh_X3) = Zero

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      G_X(iNX,1,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,2,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,3,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,4,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Alpha   ) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,5,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Beta_1  ) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,6,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Beta_2  ) * SqrtGm(iNX,iX1,iX2,iX3)
      G_X(iNX,7,iX1,iX2,iX3) &
        = G(iNX,iX1,iX2,iX3,iGF_Beta_3  ) * SqrtGm(iNX,iX1,iX2,iX3)

      U_X(iNX,iCF_D ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_D ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S1,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_S1) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S2,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_S2) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S3,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_S3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_E ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_E ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_Ne,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCF_Ne) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Compute integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nK_X, One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_X, One, U_X, nDOFX, &
             WeightsX_q, 1, Zero, UK, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nGF_X, One, G_X, nDOFX, &
             WeightsX_q, 1, Zero, GK, 1 )

    ! --- Form cell averages ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B1, iX_E1, GK, UK, Vol, PK, PrK, VK, &
    !$ACC          iX1arr, iX2arr, iX3arr, iErr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      IF( IsCornerCell &
            ( iX_B1, iX_E1, iX1arr(iX1), iX2arr(iX2), iX3arr(iX3) ) ) CYCLE

      iErr(iX1,iX2,iX3) = 0

      GK(1,iX1,iX2,iX3) = GK(1,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(2,iX1,iX2,iX3) = GK(2,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(3,iX1,iX2,iX3) = GK(3,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(4,iX1,iX2,iX3) = GK(4,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(5,iX1,iX2,iX3) = GK(5,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(6,iX1,iX2,iX3) = GK(6,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(7,iX1,iX2,iX3) = GK(7,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

      UK(iCF_D ,iX1,iX2,iX3) = UK(iCF_D ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCF_S1,iX1,iX2,iX3) = UK(iCF_S1,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCF_S2,iX1,iX2,iX3) = UK(iCF_S2,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCF_S3,iX1,iX2,iX3) = UK(iCF_S3,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCF_E ,iX1,iX2,iX3) = UK(iCF_E ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCF_Ne,iX1,iX2,iX3) = UK(iCF_Ne,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

      CALL ComputePrimitive_Euler &
           ( UK(iCF_D ,iX1,iX2,iX3), &
             UK(iCF_S1,iX1,iX2,iX3), &
             UK(iCF_S2,iX1,iX2,iX3), &
             UK(iCF_S3,iX1,iX2,iX3), &
             UK(iCF_E ,iX1,iX2,iX3), &
             UK(iCF_Ne,iX1,iX2,iX3), &
             PK(iPF_D ,iX1,iX2,iX3), &
             PK(iPF_V1,iX1,iX2,iX3), &
             PK(iPF_V2,iX1,iX2,iX3), &
             PK(iPF_V3,iX1,iX2,iX3), &
             PK(iPF_E ,iX1,iX2,iX3), &
             PK(iPF_Ne,iX1,iX2,iX3), &
             GK(1     ,iX1,iX2,iX3), &
             GK(2     ,iX1,iX2,iX3), &
             GK(3     ,iX1,iX2,iX3), &
             iErr(     iX1,iX2,iX3) )

      CALL ComputePressureFromPrimitive &
             ( PK(iPF_D ,iX1,iX2,iX3), &
               PK(iPF_E ,iX1,iX2,iX3), &
               PK(iPF_Ne,iX1,iX2,iX3), &
               PrK(      iX1,iX2,iX3) )

      VK(1,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPF_V1,iX1,iX2,iX3) + GK(5,iX1,iX2,iX3)

      VK(2,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPF_V2,iX1,iX2,iX3) + GK(6,iX1,iX2,iX3)

      VK(3,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPF_V3,iX1,iX2,iX3) + GK(7,iX1,iX2,iX3)

    END DO
    END DO
    END DO

    ! --- Shock detector, adapted from
    !     Fryxell et al., (2000), ApJS, 131, 273 ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, PrK, VK, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Compute pressure gradient and divergence of velocity (X1) ---

      GradP = ABS( PrK(iX1+1,iX2,iX3) - PrK(iX1-1,iX2,iX3) ) &
                / MIN( PrK(iX1+1,iX2,iX3), PrK(iX1-1,iX2,iX3) )

      DivV  = VK(1,iX1+1,iX2,iX3) - VK(1,iX1-1,iX2,iX3)

      IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

        DO iNX = 1, nDOFX

          D(iNX,iX1,iX2,iX3,iDF_Sh_X1) = One

        END DO

      END IF

      IF( nDimsX .GT. 1 )THEN

        ! --- Compute pressure gradient and divergence of velocity (X2) ---

        GradP = ABS( PrK(iX1,iX2+1,iX3) - PrK(iX1,iX2-1,iX3) ) &
                  / MIN( PrK(iX1,iX2+1,iX3), PrK(iX1,iX2-1,iX3) )

        DivV  = VK(2,iX1,iX2+1,iX3) - VK(2,iX1,iX2-1,iX3)

        IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

          DO iNX = 1, nDOFX

            D(iNX,iX1,iX2,iX3,iDF_Sh_X2) = One

          END DO

        END IF

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Compute pressure gradient and divergence of velocity (X3) ---

        GradP = ABS( PrK(iX1,iX2,iX3+1) - PrK(iX1,iX2,iX3-1) ) &
                  / MIN( PrK(iX1,iX2,iX3+1), PrK(iX1,iX2,iX3-1) )

        DivV  = VK(3,iX1,iX2,iX3+1) - VK(3,iX1,iX2,iX3-1)

        IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

          DO iNX = 1, nDOFX

            D(iNX,iX1,iX2,iX3,iDF_Sh_X3) = One

          END DO

        END IF

      END IF

    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    D, iErr ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
    !$OMP               SqrtGm, G_X, GK, U_X, UK, PK, VK, PrK, Vol, &
    !$OMP               iX1arr, iX2arr, iX3arr )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      D, iErr ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
    !$ACC               SqrtGm, G_X, GK, U_X, UK, PK, VK, PrK, Vol, &
    !$ACC               iX1arr, iX2arr, iX3arr )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    IF( ANY( iErr .NE. 0 ) )THEN

      PRINT*, 'Shock Detector'

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        IF( IsCornerCell &
              ( iX_B1, iX_E1, iX1arr(iX1), iX2arr(iX2), iX3arr(iX3) ) ) CYCLE

        IF( iErr(iX1,iX2,iX3) .NE. 0 )THEN

          PRINT*, 'iX1, iX2, iX3, iErr = ', iX1, iX2, iX3, iErr
          CALL DescribeError_Euler( iErr(iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DD_ShockDetector )

  END SUBROUTINE DetectShocks_Euler


END MODULE Euler_DiscontinuityDetectionModule
