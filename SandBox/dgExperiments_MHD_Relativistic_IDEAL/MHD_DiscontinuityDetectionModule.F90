MODULE MHD_DiscontinuityDetectionModule

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
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
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
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nPM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nDM, &
    iDM_TCI, &
    iDM_Sh_X1, &
    iDM_Sh_X2, &
    iDM_Sh_X3
  USE MHD_UtilitiesModule, ONLY: &
    ComputePrimitive_MHD
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE UtilitiesModule, ONLY: &
    IsCornerCell
  USE MHD_ErrorModule, ONLY: &
    DescribeError_MHD
  USE MeshModule, ONLY: &
    MeshX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTroubledCellIndicator_MHD
  PUBLIC :: FinalizeTroubledCellIndicator_MHD
  PUBLIC :: DetectTroubledCells_MHD
  PUBLIC :: DetectShocks_MHD

  ! --- For troubled-cell indicator ---
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)

  LOGICAL,  PUBLIC :: UseTroubledCellIndicator
  REAL(DP), PUBLIC :: LimiterThreshold

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_MHD_NOGPU)
  !$OMP DECLARE TARGET( UseTroubledCellIndicator, LimiterThreshold )
#elif defined(THORNADO_OACC) && !defined(THORNADO_MHD_NOGPU)
  !$ACC DECLARE CREATE( UseTroubledCellIndicator, LimiterThreshold )
#endif


CONTAINS


  SUBROUTINE InitializeTroubledCellIndicator_MHD &
    ( UseTroubledCellIndicator_Option, &
      LimiterThresholdParameter_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UseTroubledCellIndicator_Option
    REAL(DP), INTENT(in), OPTIONAL :: LimiterThresholdParameter_Option

    INTEGER  :: iNode, iNX1, iNX2, iNX3
    INTEGER  :: jNode, jNX1, jNX2, jNX3
    REAL(DP) :: WeightX
    REAL(DP) :: LimiterThresholdParameter

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    LimiterThresholdParameter = 0.03_DP
    IF( PRESENT( LimiterThresholdParameter_Option ) ) &
      LimiterThresholdParameter = LimiterThresholdParameter_Option
    LimiterThreshold = LimiterThresholdParameter * Two**( nNodes - 2 )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_MHD_NOGPU)
    !$OMP TARGET UPDATE TO( UseTroubledCellIndicator, LimiterThreshold )
#elif defined(THORNADO_OACC) && !defined(THORNADO_MHD_NOGPU)
    !$ACC UPDATE DEVICE   ( UseTroubledCellIndicator, LimiterThreshold )
#endif

    IF( .NOT. UseTroubledCellIndicator ) RETURN

    ALLOCATE( WeightsX_X1_P(nDOFX), WeightsX_X1_N(nDOFX) )
    ALLOCATE( WeightsX_X2_P(nDOFX), WeightsX_X2_N(nDOFX) )
    ALLOCATE( WeightsX_X3_P(nDOFX), WeightsX_X3_N(nDOFX) )

    ! --- Compute Weights for Extrapolating Neighbors into Target Cell ---

    DO jNode = 1, nDOFX

      jNX1 = NodeNumberTableX(1,jNode)
      jNX2 = NodeNumberTableX(2,jNode)
      jNX3 = NodeNumberTableX(3,jNode)

      WeightsX_X1_P(jNode) = Zero
      WeightsX_X1_N(jNode) = Zero
      WeightsX_X2_P(jNode) = Zero
      WeightsX_X2_N(jNode) = Zero
      WeightsX_X3_P(jNode) = Zero
      WeightsX_X3_N(jNode) = Zero

      DO iNode = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNode)
        iNX2 = NodeNumberTableX(2,iNode)
        iNX3 = NodeNumberTableX(3,iNode)

        WeightX = WeightsX_q(iNode)

        WeightsX_X1_P(jNode) &
          = WeightsX_X1_P(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) + One ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) ) )

        WeightsX_X1_N(jNode) &
          = WeightsX_X1_N(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) - One ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) ) )

        WeightsX_X2_P(jNode) &
          = WeightsX_X2_P(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) + One ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) ) )

        WeightsX_X2_N(jNode) &
          = WeightsX_X2_N(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) - One ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) ) )

        WeightsX_X3_P(jNode) &
          = WeightsX_X3_P(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) + One ) )

        WeightsX_X3_N(jNode) &
          = WeightsX_X3_N(jNode) &
              + WeightX &
                * ( L_X1  (jNX1) % P( NodesX1(iNX1) ) &
                    * L_X2(jNX2) % P( NodesX2(iNX2) ) &
                    * L_X3(jNX3) % P( NodesX3(iNX3) - One ) )

      END DO

    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_MHD_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP          WeightsX_X2_P, WeightsX_X2_N, &
    !$OMP          WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC) && !defined(THORNADO_MHD_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC          WeightsX_X2_P, WeightsX_X2_N, &
    !$ACC          WeightsX_X3_P, WeightsX_X3_N )
#endif

  END SUBROUTINE InitializeTroubledCellIndicator_MHD


  SUBROUTINE FinalizeTroubledCellIndicator_MHD

    IF( .NOT. UseTroubledCellIndicator ) RETURN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_MHD_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP               WeightsX_X2_P, WeightsX_X2_N, &
    !$OMP               WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC) && !defined(THORNADO_MHD_NOGPU)
    !$ACC EXIT DATA &
    !$ACC DELETE(     WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC             WeightsX_X2_P, WeightsX_X2_N, &
    !$ACC             WeightsX_X3_P, WeightsX_X3_N )
#endif

    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeTroubledCellIndicator_MHD


  SUBROUTINE DetectTroubledCells_MHD( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3

    IF( .NOT. UseTroubledCellIndicator )THEN

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        D(iNX,iX1,iX2,iX3,iDM_TCI) = 1.1_DP * LimiterThreshold

      END DO
      END DO
      END DO
      END DO

    ELSE

      IF( nDimsX .EQ. 1 ) &
        CALL DetectTroubledCells_MHD_nDimsX_1 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

      IF( nDimsX .EQ. 2 ) &
        CALL DetectTroubledCells_MHD_nDimsX_2 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

      IF( nDimsX .EQ. 3 ) &
        CALL DetectTroubledCells_MHD_nDimsX_3 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    END IF

  END SUBROUTINE DetectTroubledCells_MHD


  SUBROUTINE DetectTroubledCells_MHD_nDimsX_1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iFd
    INTEGER  :: nF_K
    INTEGER, PARAMETER :: nF = 3 ! [ D, E, Ye ]
    REAL(DP) :: YeBoost

    REAL(DP) :: U_X (1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: U_K(         1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: Max_UK(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))

    REAL(DP) :: U_X1(1:nDOFX,1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: U_K_X1(      1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)
    REAL(DP) :: U_K0_X1(     1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nF_K = nF * PRODUCT( iX_E0 - iX_B0 + 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X(iNX,1,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X(iNX,2,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X(iNX,3,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                 / U(iNX,iX1,iX2,iX3,iCM_D)

    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      U_X1(iNX,1,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X1(iNX,2,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X1(iNX,3,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Compute cell-averages (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X, &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    ! --- Compute cell-averages of neighbors (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_X1_P, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_X1_N, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iFd = 1, nF

      Max_UK(iFd,iX1,iX2,iX3) = MAX( ABS( U_K   (iFd,iX1,iX2,iX3) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,1) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,2) ) )

    END DO
    END DO
    END DO
    END DO

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            YeBoost &
              * (   ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,2) ) ) &
            / MAX( Max_UK   (3,iX1,iX2,iX3), SqrtTiny ) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE DetectTroubledCells_MHD_nDimsX_1


  SUBROUTINE DetectTroubledCells_MHD_nDimsX_2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iFd
    INTEGER  :: nF_K
    INTEGER, PARAMETER :: nF = 3 ! [ D, E, Ye ]
    REAL(DP) :: YeBoost

    REAL(DP) :: U_X (1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: U_K(         1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: Max_UK(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))

    REAL(DP) :: U_X1(1:nDOFX,1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: U_K_X1(      1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)
    REAL(DP) :: U_K0_X1(     1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)

    REAL(DP) :: U_X2(1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: U_K_X2(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)  :iX_E0(2),1:2)
    REAL(DP) :: U_K0_X2(     1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)  :iX_E0(2),1:2)

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nF_K = nF * PRODUCT( iX_E0 - iX_B0 + 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X(iNX,1,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X(iNX,2,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X(iNX,3,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                 / U(iNX,iX1,iX2,iX3,iCM_D)

    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      U_X1(iNX,1,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X1(iNX,2,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X1(iNX,3,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X2(iNX,1,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X2(iNX,2,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X2(iNX,3,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Compute cell-averages (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X, &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    ! --- Compute cell-averages of neighbors (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_X1_P, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_X1_N, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    ! --- Compute cell-averages (X2) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2  (1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2  (1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),2), 1 )

    ! --- Compute cell-averages of neighbors (X2) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2   (1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), &
             nDOFX, WeightsX_X2_P, 1, Zero, &
             U_K0_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2   (1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), &
             nDOFX, WeightsX_X2_N, 1, Zero, &
             U_K0_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),2), 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iFd = 1, nF

      Max_UK(iFd,iX1,iX2,iX3) = MAX( ABS( U_K   (iFd,iX1,iX2,iX3) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,1) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,2) ), &
                                     ABS( U_K_X2(iFd,iX1,iX3,iX2,1) ), &
                                     ABS( U_K_X2(iFd,iX1,iX3,iX2,2) ) )

    END DO
    END DO
    END DO
    END DO

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            YeBoost &
              * (   ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,2) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X2(3,iX1,iX3,iX2,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X2(3,iX1,iX3,iX2,2) ) ) &
            / MAX( Max_UK   (3,iX1,iX2,iX3), SqrtTiny ) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE DetectTroubledCells_MHD_nDimsX_2


  SUBROUTINE DetectTroubledCells_MHD_nDimsX_3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iFd
    INTEGER  :: nF_K
    INTEGER, PARAMETER :: nF = 3 ! [ D, E, Ye ]
    REAL(DP) :: YeBoost

    REAL(DP) :: U_X (1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: U_K(         1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))
    REAL(DP) :: Max_UK(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3))

    REAL(DP) :: U_X1(1:nDOFX,1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: U_K_X1(      1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)
    REAL(DP) :: U_K0_X1(     1:nF,iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(1)  :iX_E0(1),1:2)

    REAL(DP) :: U_X2(1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: U_K_X2(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)  :iX_E0(2),1:2)
    REAL(DP) :: U_K0_X2(     1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(3)  :iX_E0(3), &
                                  iX_B0(2)  :iX_E0(2),1:2)

    REAL(DP) :: U_X3(1:nDOFX,1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: U_K_X3(      1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3),1:2)
    REAL(DP) :: U_K0_X3(     1:nF,iX_B0(1)  :iX_E0(1), &
                                  iX_B0(2)  :iX_E0(2), &
                                  iX_B0(3)  :iX_E0(3),1:2)

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nF_K = nF * PRODUCT( iX_E0 - iX_B0 + 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X(iNX,1,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X(iNX,2,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X(iNX,3,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                 / U(iNX,iX1,iX2,iX3,iCM_D)

    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      U_X1(iNX,1,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X1(iNX,2,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X1(iNX,3,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X2(iNX,1,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X2(iNX,2,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X2(iNX,3,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      U_X3(iNX,1,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_D)
      U_X3(iNX,2,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_E)
      U_X3(iNX,3,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM_Ne) &
                                  / MAX( U(iNX,iX1,iX2,iX3,iCM_D), &
                                         SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Compute cell-averages (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X, &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1  (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    ! --- Compute cell-averages of neighbors (X1) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), &
             nDOFX, WeightsX_X1_P, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X1   (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
             nDOFX, WeightsX_X1_N, 1, Zero, &
             U_K0_X1(  1,iX_B0(2),iX_B0(3),iX_B0(1),2), 1 )

    ! --- Compute cell-averages (X2) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2  (1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2  (1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),2), 1 )

    ! --- Compute cell-averages of neighbors (X2) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2   (1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), &
             nDOFX, WeightsX_X2_P, 1, Zero, &
             U_K0_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X2   (1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), &
             nDOFX, WeightsX_X2_N, 1, Zero, &
             U_K0_X2(  1,iX_B0(1),iX_B0(3),iX_B0(2),2), 1 )

    ! --- Compute cell-averages (X3) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X3  (1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X3(  1,iX_B0(1),iX_B0(2),iX_B0(3),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X3  (1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K_X3(  1,iX_B0(1),iX_B0(2),iX_B0(3),2), 1 )

    ! --- Compute cell-averages of neighbors (X3) ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X3   (1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), &
             nDOFX, WeightsX_X3_P, 1, Zero, &
             U_K0_X3(  1,iX_B0(1),iX_B0(2),iX_B0(3),1), 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X3   (1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), &
             nDOFX, WeightsX_X3_N, 1, Zero, &
             U_K0_X3(  1,iX_B0(1),iX_B0(2),iX_B0(3),2), 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iFd = 1, nF

      Max_UK(iFd,iX1,iX2,iX3) = MAX( ABS( U_K   (iFd,iX1,iX2,iX3) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,1) ), &
                                     ABS( U_K_X1(iFd,iX2,iX3,iX1,2) ), &
                                     ABS( U_K_X2(iFd,iX1,iX3,iX2,1) ), &
                                     ABS( U_K_X2(iFd,iX1,iX3,iX2,2) ), &
                                     ABS( U_K_X3(iFd,iX1,iX2,iX3,1) ), &
                                     ABS( U_K_X3(iFd,iX1,iX2,iX3,2) ) )

    END DO
    END DO
    END DO
    END DO

    ! --- Troubled-Cell Indicator from Fu & Shu (2017) ---
    ! --- JCP, 347, 305 - 327 ----------------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X3(1,iX1,iX2,iX3,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X3(1,iX1,iX2,iX3,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X3(2,iX1,iX2,iX3,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X3(2,iX1,iX2,iX3,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNX,iX1,iX2,iX3,iDM_TCI) &
        = MAX( D(iNX,iX1,iX2,iX3,iDM_TCI), &
            YeBoost &
              * (   ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,2) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X2(3,iX1,iX3,iX2,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X2(3,iX1,iX3,iX2,2) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X3(3,iX1,iX2,iX3,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X3(3,iX1,iX2,iX3,2) ) ) &
            / MAX( Max_UK   (3,iX1,iX2,iX3), SqrtTiny ) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE DetectTroubledCells_MHD_nDimsX_3


  SUBROUTINE DetectShocks_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, EvolveOnlyMagnetic, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic
    REAL(DP), INTENT(inout) :: &
      D(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)

    INTEGER  :: iNX, iX1, iX2, iX3
    REAL(DP) :: GradP, DivV
    INTEGER  :: nK_X, nCM_X, nGF_X

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

    REAL(DP) :: U_X(1:nDOFX,1:nCM,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: UK(1:nCM,         iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: PK(1:nPM,         iX_B1(1):iX_E1(1), &
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

    INTEGER :: ITERATION(         iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))
    INTEGER :: iErr(              iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    INTEGER :: iX1arr(            iX_B1(1):iX_E1(1))
    INTEGER :: iX2arr(            iX_B1(2):iX_E1(2))
    INTEGER :: iX3arr(            iX_B1(3):iX_E1(3))

    nK_X = PRODUCT( [ iX_E1(1) - iX_B1(1) + 1, &
                      iX_E1(2) - iX_B1(2) + 1, &
                      iX_E1(3) - iX_B1(3) + 1 ] )

    nCM_X = nCM * nK_X
    nGF_X = 7   * nK_X

     DO iX1 = iX_B1(1), iX_E1(1)

       iX1arr(iX1) = iX1

     END DO

     DO iX2 = iX_B1(2), iX_E1(2)

       iX2arr(iX2) = iX2

     END DO

     DO iX3 = iX_B1(3), iX_E1(3)

       iX3arr(iX3) = iX3

     END DO

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1       , nDOFX

      D(iNX,iX1,iX2,iX3,iDM_Sh_X1) = Zero
      D(iNX,iX1,iX2,iX3,iDM_Sh_X2) = Zero
      D(iNX,iX1,iX2,iX3,iDM_Sh_X3) = Zero

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

      U_X(iNX,iCM_D  ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_D  ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_S1 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_S1 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_S2 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_S2 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_S3 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_S3 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_E  ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_E  ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_Ne ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_Ne ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_B1 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_B1 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_B2 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_B2 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_B3 ,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_B3 ) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCM_Chi,iX1,iX2,iX3) &
        = U(iNX,iX1,iX2,iX3,iCM_Chi) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    ! --- Compute integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nK_X, One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCM_X, One, U_X, nDOFX, &
             WeightsX_q, 1, Zero, UK, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nGF_X, One, G_X, nDOFX, &
             WeightsX_q, 1, Zero, GK, 1 )

    ! --- Form cell averages ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      ITERATION(iX1,iX2,iX3) = 0
      iErr     (iX1,iX2,iX3) = 0

      IF( IsCornerCell &
            ( iX_B1, iX_E1, iX1arr(iX1), iX2arr(iX2), iX3arr(iX3) ) ) CYCLE

      GK(1,iX1,iX2,iX3) = GK(1,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(2,iX1,iX2,iX3) = GK(2,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(3,iX1,iX2,iX3) = GK(3,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(4,iX1,iX2,iX3) = GK(4,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(5,iX1,iX2,iX3) = GK(5,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(6,iX1,iX2,iX3) = GK(6,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      GK(7,iX1,iX2,iX3) = GK(7,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

      UK(iCM_D  ,iX1,iX2,iX3) = UK(iCM_D  ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_S1 ,iX1,iX2,iX3) = UK(iCM_S1 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_S2 ,iX1,iX2,iX3) = UK(iCM_S2 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_S3 ,iX1,iX2,iX3) = UK(iCM_S3 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_E  ,iX1,iX2,iX3) = UK(iCM_E  ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_Ne ,iX1,iX2,iX3) = UK(iCM_Ne ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_B1 ,iX1,iX2,iX3) = UK(iCM_B1 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_B2 ,iX1,iX2,iX3) = UK(iCM_B2 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_B3 ,iX1,iX2,iX3) = UK(iCM_B3 ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      UK(iCM_Chi,iX1,iX2,iX3) = UK(iCM_Chi,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

      CALL ComputePrimitive_MHD &
           ( UK(iCM_D  ,iX1,iX2,iX3), &
             UK(iCM_S1 ,iX1,iX2,iX3), &
             UK(iCM_S2 ,iX1,iX2,iX3), &
             UK(iCM_S3 ,iX1,iX2,iX3), &
             UK(iCM_E  ,iX1,iX2,iX3), &
             UK(iCM_Ne ,iX1,iX2,iX3), &
             UK(iCM_B1 ,iX1,iX2,iX3), &
             UK(iCM_B2 ,iX1,iX2,iX3), &
             UK(iCM_B3 ,iX1,iX2,iX3), &
             UK(iCM_Chi,iX1,iX2,iX3), &
             PK(iPM_D  ,iX1,iX2,iX3), &
             PK(iPM_V1 ,iX1,iX2,iX3), &
             PK(iPM_V2 ,iX1,iX2,iX3), &
             PK(iPM_V3 ,iX1,iX2,iX3), &
             PK(iPM_E  ,iX1,iX2,iX3), &
             PK(iPM_Ne ,iX1,iX2,iX3), &
             PK(iPM_B1 ,iX1,iX2,iX3), &
             PK(iPM_B2 ,iX1,iX2,iX3), &
             PK(iPM_B3 ,iX1,iX2,iX3), &
             PK(iPM_Chi,iX1,iX2,iX3), &
             GK(1     ,iX1,iX2,iX3), &
             GK(2     ,iX1,iX2,iX3), &
             GK(3     ,iX1,iX2,iX3), &
             GK(4     ,iX1,iX2,iX3), &
             GK(5     ,iX1,iX2,iX3), &
             GK(6     ,iX1,iX2,iX3), &
             GK(7     ,iX1,iX2,iX3), &
             EvolveOnlyMagnetic )

      CALL ComputePressureFromPrimitive &
             ( PK(iPM_D ,iX1,iX2,iX3), &
               PK(iPM_E ,iX1,iX2,iX3), &
               PK(iPM_Ne,iX1,iX2,iX3), &
               PrK(      iX1,iX2,iX3) )

      VK(1,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPM_V1,iX1,iX2,iX3) + GK(5,iX1,iX2,iX3)

      VK(2,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPM_V2,iX1,iX2,iX3) + GK(6,iX1,iX2,iX3)

      VK(3,iX1,iX2,iX3) &
        = GK(4,iX1,iX2,iX3) * PK(iPM_V3,iX1,iX2,iX3) + GK(7,iX1,iX2,iX3)

    END DO
    END DO
    END DO

    ! --- Shock detector, adapted from
    !     Fryxell et al., (2000), ApJS, 131, 273 ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Compute pressure gradient and divergence of velocity (X1) ---

      GradP = ABS( PrK(iX1+1,iX2,iX3) - PrK(iX1-1,iX2,iX3) ) &
                / MIN( PrK(iX1+1,iX2,iX3), PrK(iX1-1,iX2,iX3) )

      DivV = VK(1,iX1+1,iX2,iX3) - VK(1,iX1-1,iX2,iX3)

      IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

        DO iNX = 1, nDOFX

          D(iNX,iX1,iX2,iX3,iDM_Sh_X1) = One

        END DO

      END IF

      IF( nDimsX .GT. 1 )THEN

        ! --- Compute pressure gradient and divergence of velocity (X2) ---

        GradP = ABS( PrK(iX1,iX2+1,iX3) - PrK(iX1,iX2-1,iX3) ) &
                  / MIN( PrK(iX1,iX2+1,iX3), PrK(iX1,iX2-1,iX3) )

        DivV = VK(2,iX1,iX2+1,iX3) - VK(2,iX1,iX2-1,iX3)

        IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

          DO iNX = 1, nDOFX

            D(iNX,iX1,iX2,iX3,iDM_Sh_X2) = One

          END DO

        END IF

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Compute pressure gradient and divergence of velocity (X3) ---

        GradP = ABS( PrK(iX1,iX2,iX3+1) - PrK(iX1,iX2,iX3-1) ) &
                  / MIN( PrK(iX1,iX2,iX3+1), PrK(iX1,iX2,iX3-1) )

        DivV = VK(3,iX1,iX2,iX3+1) - VK(3,iX1,iX2,iX3-1)

        IF( GradP .GT. Third .AND. DivV .LT. Zero )THEN

          DO iNX = 1, nDOFX

            D(iNX,iX1,iX2,iX3,iDM_Sh_X3) = One

          END DO

        END IF

      END IF

    END DO
    END DO
    END DO

#ifdef HYDRO_RELATIVISTIC

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) '    MODULE: MHD_DiscontinuityDetectionModule'
      WRITE(*,*) 'SUBROUTINE: DetectShocks_MHD'

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      END DO
      END DO
      END DO

    END IF

#endif

  END SUBROUTINE DetectShocks_MHD


END MODULE MHD_DiscontinuityDetectionModule
