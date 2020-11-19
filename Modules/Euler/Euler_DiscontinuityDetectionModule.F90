MODULE Euler_DiscontinuityDetectionModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
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
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_TroubledCellIndicator, &
    Timer_Euler_ShockDetector, &
    Timer_Euler_CopyIn, &
    Timer_Euler_CopyOut, &
    Timer_Euler_Permute

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP          WeightsX_X2_P, WeightsX_X2_N, &
    !$OMP          WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC          WeightsX_X2_P, WeightsX_X2_N, &
    !$ACC          WeightsX_X3_P, WeightsX_X3_N )
#endif

  END SUBROUTINE InitializeTroubledCellIndicator_Euler


  SUBROUTINE FinalizeTroubledCellIndicator_Euler

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: WeightsX_X1_P, WeightsX_X1_N, &
    !$OMP               WeightsX_X2_P, WeightsX_X2_N,  &
    !$OMP               WeightsX_X3_P, WeightsX_X3_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(     WeightsX_X1_P, WeightsX_X1_N, &
    !$ACC             WeightsX_X2_P, WeightsX_X2_N,  &
    !$ACC             WeightsX_X3_P, WeightsX_X3_N )
#endif

    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeTroubledCellIndicator_Euler


  SUBROUTINE DetectTroubledCells_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    IF     ( nDimsX .EQ. 1 )THEN

      CALL DetectTroubledCells_Euler_nDimsX1 &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    ELSE IF( nDimsX .EQ. 2 )THEN

      CALL DetectTroubledCells_Euler_nDimsX2 &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    ELSE

      CALL DetectTroubledCells_Euler_nDimsX3 &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    END IF

  END SUBROUTINE DetectTroubledCells_Euler


  SUBROUTINE DetectTroubledCells_Euler_nDimsX1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iFd
    INTEGER  :: nX(3), nF_K
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

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, D )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    IF( .NOT. UseTroubledCellIndicator )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = 1.1_DP * LimiterThreshold

      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(   D ) &
    !$ACC DELETE(       iX_B0, iX_E0 )
#endif

      CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

      RETURN

    ELSE

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B1, iX_E1, U ), &
    !$OMP MAP( alloc: U_X1  , U_X, &
    !$OMP             U_K_X1, U_K, Max_UK, &
    !$OMP             U_K0_X1 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B1, iX_E1, U ) &
    !$ACC CREATE(     U_X1  , U_X, &
    !$ACC             U_K_X1, U_K, Max_UK, &
    !$ACC             U_K0_X1 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = Zero

      END DO
      END DO
      END DO
      END DO

    END IF

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nX   = iX_E0 - iX_B0 + 1
    nF_K = nF * PRODUCT( nX )

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, U_X, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X(iNodeX,1,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X(iNodeX,2,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X(iNodeX,3,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                    / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X1, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      U_X1(iNodeX,1,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X1(iNodeX,2,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X1(iNodeX,3,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                     / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Compute cell-averages  ---

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, Max_UK, U_K, U_K_X1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_K, U_K0_X1, Max_UK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
            YeBoost &
              * (   ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,1) ) &
                  + ABS( U_K(3,iX1,iX2,iX3) - U_K0_X1(3,iX2,iX3,iX1,2) ) ) &
            / MAX( Max_UK   (3,iX1,iX2,iX3), SqrtTiny ) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$OMP               U_X1  , U_X, &
    !$OMP               U_K_X1, U_K, Max_UK, &
    !$OMP               U_K0_X1 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$ACC               U_X1   , U_X, &
    !$ACC               U_K_X1 , U_K, Max_UK, &
    !$ACC               U_K0_X1 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

  END SUBROUTINE DetectTroubledCells_Euler_nDimsX1


  SUBROUTINE DetectTroubledCells_Euler_nDimsX2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iFd
    INTEGER  :: nX(3), nF_K
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

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, D )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    IF( .NOT. UseTroubledCellIndicator )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = 1.1_DP * LimiterThreshold

      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(   D ) &
    !$ACC DELETE(       iX_B0, iX_E0 )
#endif

      CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

      RETURN

    ELSE

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B1, iX_E1, U ), &
    !$OMP MAP( alloc: U_X1   , U_X2   , U_X, &
    !$OMP             U_K_X1 , U_K_X2 , U_K, Max_UK, &
    !$OMP             U_K0_X1, U_K0_X2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B1, iX_E1, U ) &
    !$ACC CREATE(     U_X1   , U_X2   , U_X, &
    !$ACC             U_K_X1 , U_K_X2 , U_K, Max_UK, &
    !$ACC             U_K0_X1, U_K0_X2 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = Zero

      END DO
      END DO
      END DO
      END DO

    END IF

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nX   = iX_E0 - iX_B0 + 1
    nF_K = nF * PRODUCT( nX )

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, U_X, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X(iNodeX,1,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X(iNodeX,2,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X(iNodeX,3,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                    / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X1, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      U_X1(iNodeX,1,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X1(iNodeX,2,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X1(iNodeX,3,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                    / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X2, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X2(iNodeX,1,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X2(iNodeX,2,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X2(iNodeX,3,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                     / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Compute cell-averages  ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X, &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K, 1 )

    ! --- X1 ---

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

    ! --- X2 ---

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, Max_UK, U_K, U_K_X1, U_K_X2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_K, U_K0_X1, U_K0_X2, Max_UK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
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

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$OMP               U_X1   , U_X2   , U_X, &
    !$OMP               U_K_X1 , U_K_X2 , U_K, Max_UK, &
    !$OMP               U_K0_X1, U_K0_X2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$ACC               U_X1   , U_X2   , U_X, &
    !$ACC               U_K_X1 , U_K_X2 , U_K, Max_UK, &
    !$ACC               U_K0_X1, U_K0_X2 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

  END SUBROUTINE DetectTroubledCells_Euler_nDimsX2


  SUBROUTINE DetectTroubledCells_Euler_nDimsX3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iFd
    INTEGER  :: nX(3), nF_K
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

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, D )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

    IF( .NOT. UseTroubledCellIndicator )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = 1.1_DP * LimiterThreshold

      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(   D ) &
    !$ACC DELETE(       iX_B0, iX_E0 )
#endif

      CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

      RETURN

    ELSE

    CALL TimersStart_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B1, iX_E1, U ), &
    !$OMP MAP( alloc: U_X1   , U_X2   , U_X3  , U_X, &
    !$OMP             U_K_X1 , U_K_X2 , U_K_X3, U_K, Max_UK, &
    !$OMP             U_K0_X1, U_K0_X2, U_K0_X3 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B1, iX_E1, U ) &
    !$ACC CREATE(     U_X1   , U_X2   , U_X3  , U_X, &
    !$ACC             U_K_X1 , U_K_X2 , U_K_X3, U_K, Max_UK, &
    !$ACC             U_K0_X1, U_K0_X2, U_K0_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        D(iNodeX,iX1,iX2,iX3,iDF_TCI) = Zero

      END DO
      END DO
      END DO
      END DO

    END IF

#ifdef MICROPHYSICS_WEAKLIB

    YeBoost = 1.0e2_DP

#else

    YeBoost = Zero

#endif

    nX   = iX_E0 - iX_B0 + 1
    nF_K = nF * PRODUCT( nX )

    CALL TimersStart_Euler( Timer_Euler_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRESENT( iX_B0, iX_E0, U_X, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X(iNodeX,1,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X(iNodeX,2,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X(iNodeX,3,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                    / U(iNodeX,iX1,iX2,iX3,iCF_D)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X1, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNodeX = 1, nDOFX

      U_X1(iNodeX,1,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X1(iNodeX,2,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X1(iNodeX,3,iX2,iX3,iX1) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                     / MAX( U(iNodeX,iX1,iX2,iX3,iCF_D), &
                                            SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X2, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X2(iNodeX,1,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X2(iNodeX,2,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X2(iNodeX,3,iX1,iX3,iX2) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                     / MAX( U(iNodeX,iX1,iX2,iX3,iCF_D), &
                                            SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, iX_B1, iX_E1, U_X3, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      U_X3(iNodeX,1,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_D)
      U_X3(iNodeX,2,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_E)
      U_X3(iNodeX,3,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_Ne) &
                                     / MAX( U(iNodeX,iX1,iX2,iX3,iCF_D), &
                                            SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Permute )

    ! --- Compute cell-averages  ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nF_K, One, &
             U_X, &
             nDOFX, WeightsX_q, 1, Zero, &
             U_K, 1 )

    ! --- X1 ---

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

    ! --- X2 ---

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

    ! --- X3 ---

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, Max_UK, U_K, U_K_X1, U_K_X2, U_K_X3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_K, U_K0_X1, U_K0_X2, U_K0_X3, Max_UK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      ! --- Use Conserved Density to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = (   ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X1(1,iX2,iX3,iX1,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X2(1,iX1,iX3,iX2,2) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X3(1,iX1,iX2,iX3,1) ) &
            + ABS( U_K(1,iX1,iX2,iX3) - U_K0_X3(1,iX1,iX2,iX3,2) ) ) &
          / Max_UK    (1,iX1,iX2,iX3)

      ! --- Use Conserved Energy to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
            (   ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X1(2,iX2,iX3,iX1,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X2(2,iX1,iX3,iX2,2) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X3(2,iX1,iX2,iX3,1) ) &
              + ABS( U_K(2,iX1,iX2,iX3) - U_K0_X3(2,iX1,iX2,iX3,2) ) ) &
            / Max_UK    (2,iX1,iX2,iX3) )

      ! --- Use Electron Fraction to Detect Troubled Cell ---

      D(iNodeX,iX1,iX2,iX3,iDF_TCI) &
        = MAX( D(iNodeX,iX1,iX2,iX3,iDF_TCI), &
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

    CALL TimersStart_Euler( Timer_Euler_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$OMP               U_X1   , U_X2   , U_X3   , U_X, &
    !$OMP               U_K_X1 , U_K_X2 , U_K_X3 , U_K, Max_UK, &
    !$OMP               U_K0_X1, U_K0_X2, U_K0_X3 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, U, &
    !$ACC               U_X1   , U_X2   , U_X3  , U_X, &
    !$ACC               U_K_X1 , U_K_X2 , U_K_X3, U_K, Max_UK, &
    !$ACC               U_K0_X1, U_K0_X2, U_K0_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

  END SUBROUTINE DetectTroubledCells_Euler_nDimsX3


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
    REAL(DP) :: GradP, DivV

    ! --- Shock detector, adapted from
    !     Fryxell et al., (2000), ApJS, 131, 273 ---

    CALL TimersStart_Euler( Timer_Euler_ShockDetector )

    D(:,:,:,:,iDF_Sh_X1) = Zero
    D(:,:,:,:,iDF_Sh_X2) = Zero
    D(:,:,:,:,iDF_Sh_X3) = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

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

      GradP = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

      DivV  = VX_K(2) - VX_K(1)

      IF( GradP .GT. Third .AND. DivV .LT. Zero ) &
        D(:,iX1,iX2,iX3,iDF_Sh_X1) = One

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

        GradP = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

        DivV  = VX_K(2) - VX_K(1)

        IF( GradP .GT. Third .AND. DivV .LT. Zero ) &
          D(:,iX1,iX2,iX3,iDF_Sh_X2) = One

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

        GradP = ABS( P_K(2) - P_K(1) ) / MIN( P_K(2), P_K(1) )

        DivV  = VX_K(2) - VX_K(1)

        IF( GradP .GT. Third .AND. DivV .LT. Zero ) &
          D(:,iX1,iX2,iX3,iDF_Sh_X3) = One

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ShockDetector )

  END SUBROUTINE DetectShocks_Euler


END MODULE Euler_DiscontinuityDetectionModule

