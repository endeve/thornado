MODULE SlopeLimiterModule_Euler

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX, nNodes
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q
  USE UtilitiesModule, ONLY: &
    MinModB
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_E, &
    Shock
  USE BoundaryConditionsModule_Beta, ONLY: &
    ApplyBoundaryConditions_Fluid
  USE CharacteristicDecompositionModule, ONLY: &
    ComputeCharacteristicDecomposition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler
  PUBLIC :: FinalizeSlopeLimiter_Euler
  PUBLIC :: ApplySlopeLimiter_Euler

  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseTroubledCellIndicator
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: LimiterThreshold
  REAL(DP) :: LimiterThresholdParameter
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseTroubledCellIndicator_Option, &
      LimiterThresholdParameter_Option )

    REAL(DP), INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseTroubledCellIndicator_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option

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

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator &
        = UseTroubledCellIndicator_Option

    LimiterThresholdParameter = 0.03_DP
    IF( PRESENT( LimiterThresholdParameter_Option ) ) &
      LimiterThresholdParameter &
        = LimiterThresholdParameter_Option

    LimiterThreshold = LimiterThresholdParameter * 2.0_DP**( nNodes - 2 )

    WRITE(*,*)
    WRITE(*,'(A)') '  INFO: InitializeSlopeLimiter_Euler:'
    WRITE(*,'(A)') '  -----------------------------------'
    WRITE(*,*)
    WRITE(*,'(A4,A26,L1)'   ) '', 'UseSlopeLimiter: ' , &
      UseSlopeLimiter
    WRITE(*,*)
    WRITE(*,'(A4,A26,ES9.3E2)' ) '', 'BetaTVD: ' , &
      BetaTVD
    WRITE(*,'(A4,A26,ES9.3E2)' ) '', 'BetaTVB: ' , &
      BetaTVB
    WRITE(*,'(A4,A26,ES9.3E2)' ) '', 'SlopeTolerance: ' , &
      SlopeTolerance
    WRITE(*,*)
    WRITE(*,'(A4,A26,L1)'   ) '', 'UseTroubledCellIndicator: ' , &
      UseTroubledCellIndicator
    WRITE(*,*)
    WRITE(*,'(A4,A26,ES9.3E2)' ) '', 'LimiterThreshold: ' , &
      LimiterThreshold

    IF( UseTroubledCellIndicator )THEN

      CALL InitializeTroubledCellIndicator

    END IF

  END SUBROUTINE InitializeSlopeLimiter_Euler


  SUBROUTINE FinalizeSlopeLimiter_Euler

    IF( UseTroubledCellIndicator )THEN

      CALL FinalizeTroubledCellIndicator

    END IF

  END SUBROUTINE FinalizeSlopeLimiter_Euler


  SUBROUTINE ApplySlopeLimiter_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: SlopeDifference(nCF)
    REAL(DP) :: dU (nCF,nDimsX)
    REAL(DP) :: U_M(nCF,0:2*nDimsX,nDOFX)
    REAL(DP) :: U_K(nCF,0:1)
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL ApplyBoundaryConditions_Fluid &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          IF( Shock(iX1,iX2,iX3) < LimiterThreshold ) CYCLE

          dX1 = MeshX(1) % Width(iX1)
          dX2 = MeshX(2) % Width(iX2)
          dX3 = MeshX(3) % Width(iX3)

          ! --- Map to Modal Representation ---

          DO iCF = 1, nCF

            CALL MapNodalToModal_Fluid &
                   ( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

            ! --- Cell Average of Neighbors in X1 Direction ---

            U_M(iCF,1,1) &
              = DOT_PRODUCT( WeightsX_q(:), U(:,iX1-1,iX2,iX3,iCF) )

            U_M(iCF,2,1) &
              = DOT_PRODUCT( WeightsX_q(:), U(:,iX1+1,iX2,iX3,iCF) )

            IF( nDimsX > 1 )THEN

              ! --- Cell Average of Neighbors in X2 Direction ---

              U_M(iCF,3,1) &
                = DOT_PRODUCT( WeightsX_q(:), U(:,iX1,iX2-1,iX3,iCF) )

              U_M(iCF,4,1) &
                = DOT_PRODUCT( WeightsX_q(:), U(:,iX1,iX2+1,iX3,iCF) )

            END IF

            IF( nDimsX > 2 )THEN

              ! --- Cell Average of Neighbors in X3 Direction ---

              U_M(iCF,5,1) &
                = DOT_PRODUCT( WeightsX_q(:), U(:,iX1,iX2,iX3-1,iCF) )

              U_M(iCF,6,1) &
                = DOT_PRODUCT( WeightsX_q(:), U(:,iX1,iX2,iX3+1,iCF) )

            END IF

          END DO

          ! --- Compute Eigenvectors ---

          CALL ComputeCharacteristicDecomposition &
                 ( 1, U_M(:,0,1), R_X1, invR_X1 )

          ! --- Compute Limited Slopes ---

          dU(:,1) = MinModB &
                      ( U_M(:,0,2), &
                        BetaTVD * ( U_M(:,0,1) - U_M(:,1,1) ), &
                        BetaTVD * ( U_M(:,2,1) - U_M(:,0,1) ), &
                        dX1, BetaTVB )

          IF( nDimsX > 1 )THEN

            dU(:,2) = MinModB &
                        ( U_M(:,0,3), &
                          BetaTVD * ( U_M(:,0,1) - U_M(:,3,1) ), &
                          BetaTVD * ( U_M(:,4,1) - U_M(:,0,1) ), &
                          dX2, BetaTVB )

          END IF

          IF( nDimsX > 2 )THEN

            dU(:,3) = MinModB &
                        ( U_M(:,0,4), &
                          BetaTVD * ( U_M(:,0,1) - U_M(:,5,1) ), &
                          BetaTVD * ( U_M(:,6,1) - U_M(:,0,1) ), &
                          dX3, BetaTVB )

          END IF

          ! --- Compare Limited Slopes to Original Slopes ---

          DO iCF = 1, nCF

            SlopeDifference(iCF) &
              = ABS( U_M(iCF,0,2) - dU(iCF,1) )

            IF( nDimsX > 1 )THEN

              SlopeDifference(iCF) &
                = MAX( SlopeDifference(iCF), &
                       ABS( U_M(iCF,0,3) - dU(iCF,2) ) )

            END IF

            IF( nDimsX > 2 )THEN

              SlopeDifference(iCF) &
                = MAX( SlopeDifference(iCF), &
                       ABS( U_M(iCF,0,4) - dU(iCF,3) ) )

            END IF

          END DO

          ! --- Replace Slopes and Discard High-Order Components ---
          ! --- if Limited Slopes Deviate too Much from Original ---

          DO iCF = 1, nCF

            IF( SlopeDifference(iCF) &
                  > SlopeTolerance * ABS( U_M(iCF,0,1) ) )THEN

              U_M(iCF,0,2:nDOFX) = Zero

              U_M(iCF,0,2) = dU(iCF,1)

              IF( nDimsX > 1 ) &
                U_M(iCF,0,3) = dU(iCF,2)

              IF( nDimsX > 2 ) &
                U_M(iCF,0,4) = dU(iCF,3)

              CALL MapModalToNodal_Fluid &
                     ( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

            END IF

          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ApplySlopeLimiter_Euler


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
                     ( WeightsX_q(:), G(:,iX1,  iX2,iX3,iGF_SqrtGm) )

          V_K(1) = DOT_PRODUCT &
                     ( WeightsX_q(:), G(:,iX1-1,iX2,iX3,iGF_SqrtGm) )

          V_K(2) = DOT_PRODUCT &
                     ( WeightsX_q(:), G(:,iX1+1,iX2,iX3,iGF_SqrtGm) )

          DO iCF = 1, nCF

            U_K(0,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1,iX2,iX3,iCF) ) / V_K(0)

            U_K(1,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(1)

            U_K0(1,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_X1_P(:), &
                    G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(0)

            U_K(2,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(2)

            U_K0(2,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_X1_N(:), &
                    G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(0)

          END DO

          ! --- Compute Cell Volumes and Cell Averages ---
          ! --- in Neighbors in X2 Direction -------------

          IF( nDimsX > 1 )THEN

            V_K(3) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2-1,iX3,iGF_SqrtGm) )

            V_K(4) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2+1,iX3,iGF_SqrtGm) )

            DO iCF = 1, nCF

              U_K(3,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(3)

              U_K0(3,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X2_P(:), &
                      G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(0)

              U_K(4,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(4)

              U_K0(4,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X2_N(:), &
                      G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(0)

            END DO

          END IF

          ! --- Compute Cell Volumes and Cell Averages ---
          ! --- in Neighbors in X3 Direction -------------

          IF( nDimsX > 2 )THEN

            V_K(5) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2,iX3-1,iGF_SqrtGm) )

            V_K(6) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2,iX3+1,iGF_SqrtGm) )

            DO iCF = 1, nCF

              U_K(5,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(5)

              U_K0(5,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X3_P(:), &
                      G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(0)

              U_K(6,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3+1,iCF) ) / V_K(6)

              U_K0(6,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X3_N(:), &
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


END MODULE SlopeLimiterModule_Euler
