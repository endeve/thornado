MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Three, SqrtTiny
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
    MeshX, &
    NodeCoordinate
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
    iDF_Sh
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL, ONLY: &
    ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_SlopeLimiter, &
    Timer_Euler_TroubledCellIndicator

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplySlopeLimiter_Euler_Relativistic_IDEAL

  LOGICAL      :: UseSlopeLimiter
  LOGICAL      :: UseCharacteristicLimiting
  LOGICAL      :: UseConservativeCorrection
  LOGICAL      :: UseTroubledCellIndicator
  CHARACTER(4) :: SlopeLimiterMethod

  ! --- WENO Limiter ---
  REAL(DP), ALLOCATABLE :: OrthonormalBasis(:,:,:)
  REAL(DP), ALLOCATABLE :: VandermondeMatrix(:,:)

  ! --- TVD Limiter ---
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: LimiterThreshold
  REAL(DP) :: LimiterThresholdParameter
  REAL(DP) :: I_6x6(1:6,1:6)

  ! --- For troubled-cell indicator ---
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)

  LOGICAL :: DEBUG = .FALSE.

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, SlopeLimiterMethod_Option, &
      LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, Verbose_Option )

    REAL(DP), INTENT(in),     OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in),     OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in),     OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option,  &
      UseConservativeCorrection_Option, &
      Verbose_Option
    REAL(DP), INTENT(in),     OPTIONAL :: &
      LimiterThresholdParameter_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: &
      SlopeLimiterMethod_Option

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

    SlopeLimiterMethod = 'TVD'
    IF( PRESENT( SlopeLimiterMethod_Option ) )&
      SlopeLimiterMethod = SlopeLimiterMethod_Option

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
      WRITE(*,'(A)') &
        '    INFO: InitializeSlopeLimiter_Euler_Relativistic_IDEAL'
      WRITE(*,'(A)') &
        '    -----------------------------------------------------'
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
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,*)
      WRITE(*,'(A4,A27,A)')         '', 'SlopeLimiterMethod: ', &
        TRIM( SlopeLimiterMethod )
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

    IF( TRIM( SlopeLimiterMethod ).EQ. 'WENO' )THEN

      CALL InitializeSlopeLimiter_Euler_WENO

    END IF

  END SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)          :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    SELECT CASE( TRIM( SlopeLimiterMethod ) )

      CASE( 'TVD' )

write(*,*) 'TVD'
        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

      CASE( 'WENO' )

write(*,*) 'WENO'
        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

      CASE DEFAULT

        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    END SELECT

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)          :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
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
    REAL(DP) :: U_M(nCF,0:2*nDimsX,nDOFX,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)
    REAL(DP) :: R_X2(nCF,nCF), invR_X2(nCF,nCF)
    REAL(DP) :: R_X3(nCF,nCF), invR_X3(nCF,nCF)
    REAL(DP) :: V_K(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_Sh) .LT. LimiterThreshold ) ) CYCLE

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

      ! --- Map to Modal Representation ---

      DO iCF = 1, nCF

        CALL MapNodalToModal_Fluid( U(:,iX1,iX2,iX3,iCF), &
                                    U_M(iCF,0,:,iX1,iX2,iX3) )

        ! --- Cell Average of Neighbors in X1 Direction ---

        U_M(iCF,1,1,iX1,iX2,iX3) &
          = DOT_PRODUCT( WeightsX_q, U(:,iX1-1,iX2,iX3,iCF) )

        U_M(iCF,2,1,iX1,iX2,iX3) &
          = DOT_PRODUCT( WeightsX_q, U(:,iX1+1,iX2,iX3,iCF) )

        IF( nDimsX .GT. 1 )THEN

          ! --- Cell Average of Neighbors in X2 Direction ---

          U_M(iCF,3,1,iX1,iX2,iX3) &
            = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2-1,iX3,iCF) )

          U_M(iCF,4,1,iX1,iX2,iX3) &
            = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2+1,iX3,iCF) )

        END IF

        IF( nDimsX .GT. 2 )THEN

          ! --- Cell Average of Neighbors in X3 Direction ---

          U_M(iCF,5,1,iX1,iX2,iX3) &
            = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2,iX3-1,iCF) )

          U_M(iCF,6,1,iX1,iX2,iX3) &
            = DOT_PRODUCT( WeightsX_q, U(:,iX1,iX2,iX3+1,iCF) )

        END IF

      END DO

      IF( UseCharacteristicLimiting )THEN

        ! --- Cell Average of Geometry (Spatial Metric, Lapse Function,
        !     and Shift Vector) ---

        DO iGF = iGF_Gm_dd_11, iGF_Beta_3
          G_K(iGF) = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF) )
        END DO

        ! --- Compute Eigenvectors ---

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 1, G_K, U_M(:,0,1,iX1,iX2,iX3), R_X1, invR_X1 )

        IF( nDimsX .GT. 1 )THEN

          CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
                 ( 2, G_K, U_M(:,0,1,iX1,iX2,iX3), R_X2, invR_X2 )

        END IF

        IF( nDimsX .GT. 2 )THEN

          CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
                 ( 3, G_K, U_M(:,0,1,iX1,iX2,iX3), R_X3, invR_X3 )

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
            ( MATMUL( invR_X1, U_M(:,0,2,iX1,iX2,iX3) ), &
              BetaTVD * MATMUL( invR_X1, &
                          (U_M(:,0,1,iX1,iX2,iX3)-U_M(:,1,1,iX1,iX2,iX3)) ), &
              BetaTVD * MATMUL( invR_X1, &
                          (U_M(:,2,1,iX1,iX2,iX3)-U_M(:,0,1,iX1,iX2,iX3)) ), &
              dX1, BetaTVB )

      IF( nDimsX .GT. 1 )THEN

        dU(:,2) &
          = MinModB &
              ( MATMUL( invR_X2, U_M(:,0,3,iX1,iX2,iX3) ), &
                BetaTVD * MATMUL( invR_X2, &
                            (U_M(:,0,1,iX1,iX2,iX3)-U_M(:,3,1,iX1,iX2,iX3)) ), &
                BetaTVD * MATMUL( invR_X2, &
                            (U_M(:,4,1,iX1,iX2,iX3)-U_M(:,0,1,iX1,iX2,iX3)) ), &
                dX2, BetaTVB )

      END IF

      IF( nDimsX .GT. 2 )THEN

        dU(:,3) &
          = MinModB &
              ( MATMUL( invR_X3, U_M(:,0,4,iX1,iX2,iX3) ), &
                BetaTVD * MATMUL( invR_X3, &
                            (U_M(:,0,1,iX1,iX2,iX3)-U_M(:,5,1,iX1,iX2,iX3)) ), &
                BetaTVD * MATMUL( invR_X3, &
                            (U_M(:,6,1,iX1,iX2,iX3)-U_M(:,0,1,iX1,iX2,iX3)) ), &
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

        SlopeDifference(iCF) = ABS( U_M(iCF,0,2,iX1,iX2,iX3) - dU(iCF,1) )

        IF( nDimsX .GT. 1 )THEN

          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(iCF,0,3,iX1,iX2,iX3) - dU(iCF,2) ) )

        END IF

        IF( nDimsX .GT. 2 )THEN

          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(iCF,0,4,iX1,iX2,iX3) - dU(iCF,3) ) )

        END IF

      END DO

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      DO iCF = 1, nCF

        IF( SlopeDifference(iCF) &
              .GT. SlopeTolerance * ABS( U_M(iCF,0,1,iX1,iX2,iX3) ) )THEN

          U_M(iCF,0,2:nDOFX,iX1,iX2,iX3) = Zero

          U_M(iCF,0,2,iX1,iX2,iX3) = dU(iCF,1)

          IF( nDimsX .GT. 1 ) U_M(iCF,0,3,iX1,iX2,iX3) = dU(iCF,2)

          IF( nDimsX .GT. 2 ) U_M(iCF,0,4,iX1,iX2,iX3) = dU(iCF,3)

          LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

        END IF

      END DO

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      IF( LimitedCell(iCF,iX1,iX2,iX3) )THEN
        CALL MapModalToNodal_Fluid( U(:,iX1,iX2,iX3,iCF), &
                                    U_M(iCF,0,:,iX1,iX2,iX3) )
      END IF
    END DO
    END DO
    END DO
    END DO

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, V_K, U, U_K, LimitedCell )

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD


  SUBROUTINE InitializeSlopeLimiter_Euler_WENO

    INTEGER  :: iNodeX, jNodeX
    REAL(DP) :: eta

    IF( nDimsX .GT. 1 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nDimsX > 1. Stopping...'
      STOP

    END  IF

    IF( nNodesX(1) .GT. 3 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nNodesX(1) > 3. Stopping...'
      STOP

    END  IF

    ALLOCATE( OrthonormalBasis(nDOFX,nDOFX,nDOFX) )
    ALLOCATE( VandermondeMatrix(nDOFX,nDOFX) )

    IF( nNodesX(1) .EQ. 2 )THEN
      DO iNodeX = 1, nDOFX

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis(iNodeX,1,1) &
          = One
        OrthonormalBasis(iNodeX,2,1) &
          = SQRT(12.0_DP) * eta

        ! --- 1st derivative ---
        OrthonormalBasis(iNodeX,1,2) &
          = Zero
        OrthonormalBasis(iNodeX,2,2) &
          = SQRT(12.0_DP)

      END DO
    END IF

    IF( nNodesX(1) .EQ. 3 )THEN
      DO iNodeX = 1, nDOFX

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis(iNodeX,1,1) &
          = One
        OrthonormalBasis(iNodeX,2,1) &
          = SQRT(12.0_DP) * eta
        OrthonormalBasis(iNodeX,3,1) &
          = SQRT(180.0_DP) * ( eta**2 - One / 12.0_DP )

        ! --- 1st derivative ---
        OrthonormalBasis(iNodeX,1,2) &
          = Zero
        OrthonormalBasis(iNodeX,2,2) &
          = SQRT(12.0_DP)
        OrthonormalBasis(iNodeX,3,2) &
          = SQRT(180.0_DP) * Two * eta

        ! --- 2nd derivative ---
        OrthonormalBasis(iNodeX,1,3) &
          = Zero
        OrthonormalBasis(iNodeX,2,3) &
          = Zero
        OrthonormalBasis(iNodeX,3,3) &
          = SQRT(180.0_DP) * Two

      END DO
    END IF

    IF( nNodesX(1) .EQ. 4 )THEN
      DO iNodeX = 1, nDOFX

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis(iNodeX,1,1) &
          = One
        OrthonormalBasis(iNodeX,2,1) &
          = SQRT(12.0_DP) * eta
        OrthonormalBasis(iNodeX,3,1) &
          = SQRT(180.0_DP) * ( eta**2 - One / 12.0_DP )
        OrthonormalBasis(iNodeX,4,1) &
          = SQRT(2800.0_DP) * ( eta**3 - 15.0_DP / 100.0_DP * eta )

        ! --- 1st derivative ---
        OrthonormalBasis(iNodeX,1,2) &
          = One
        OrthonormalBasis(iNodeX,2,2) &
          = SQRT(12.0_DP)
        OrthonormalBasis(iNodeX,3,2) &
          = SQRT(180.0_DP) * Two * eta
        OrthonormalBasis(iNodeX,4,2) &
          = SQRT(2800.0_DP) * ( Three * eta**2 - 15.0_DP / 100.0_DP )

        ! --- 2nd derivative ---
        OrthonormalBasis(iNodeX,1,3) &
          = Zero
        OrthonormalBasis(iNodeX,2,3) &
          = Zero
        OrthonormalBasis(iNodeX,3,3) &
          = SQRT(180.0_DP) * Two
        OrthonormalBasis(iNodeX,4,3) &
          = SQRT(2800.0_DP) * 6.0_DP * eta

        ! --- 3rd derivative ---
        OrthonormalBasis(iNodeX,1,4) &
          = Zero
        OrthonormalBasis(iNodeX,2,4) &
          = Zero
        OrthonormalBasis(iNodeX,3,4) &
          = Zero
        OrthonormalBasis(iNodeX,4,4) &
          = SQRT(2800.0_DP) * 6.0_DP

      END DO
    END IF

    DO iNodeX = 1, nDOFX ! Polynomial

      DO jNodeX = 1, nDOFX ! Grid points

        VandermondeMatrix(iNodeX,jNodeX) &
          = WeightsX_q(jNodeX) * OrthonormalBasis(jNodeX,iNodeX,1)

      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter_Euler_WENO


  SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

    IF( UseTroubledCellIndicator ) &
      CALL FinalizeTroubledCellIndicator

    IF( TRIM( SlopeLimiterMethod ) .EQ. 'WENO' )THEN

      DEALLOCATE( VandermondeMatrix )
      DEALLOCATE( OrthonormalBasis )

    END IF

  END SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)          :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3, &
                iNodeX, jNodeX, iCF
    REAL(DP) :: dX1, dX2, dX3
    LOGICAL  :: SuppressBC
    REAL(DP) :: U_M(nDOFX,nCF,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))
    REAL(DP) :: UU(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),nCF)
    LOGICAL  :: LimitedCell(iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    ! --- WENO Limiter ---
    REAL(DP) :: q(nDOFX,nDOFX,nCF)
    REAL(DP) :: LinearWeights(2*(nNodesX(1)-1))
    REAL(DP) :: NonLinearWeights(2*(nNodesX(1)-1),nCF)
    REAL(DP) :: SmoothnessIndicators(2*(nNodesX(1)-1),nCF)
    REAL(DP) :: p (nDOFX,2*(nNodesX(1)-1),nCF)
    REAL(DP) :: pC(nDOFX,2*(nNodesX(1)-1),nCF)
    REAL(DP) :: pCoeffs(2*(nNodesX(1)-1),nNodesX(1),nCF)
    REAL(DP) :: pNew(nDOFX,nCF)

    ! --- Characteristic limiting ---
    INTEGER  :: iGF
    REAL(DP) :: R_X1(nCF,nCF,0:2), invR_X1(nCF,nCF,0:2)
    REAL(DP) :: G_K(nGF)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    ! --- Get linear weights ---
    DO iNodeX = 1, 2*(nNodesX(1)-1)

      IF( MOD( iNodeX, 2 ) .NE. 0 )THEN

        LinearWeights(iNodeX) = 0.01_DP

      ELSE

        LinearWeights(iNodeX) = 0.99_DP

      END IF

    END DO

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_Sh) .LT. LimiterThreshold ) ) CYCLE

      LimitedCell(iX1,iX2,iX3) = .TRUE.

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      DO iCF = 1, nCF

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1,iX2,iX3,iCF), &
                 U_M(:,iCF,iX1,iX2,iX3) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1-1,iX2,iX3,iCF), &
                 U_M(:,iCF,iX1-1,iX2,iX3) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1+1,iX2,iX3,iCF), &
                 U_M(:,iCF,iX1+1,iX2,iX3) )

      END DO

      ! --- Step 2.1: Compute left and right eigenvectors ---

      IF( UseCharacteristicLimiting )THEN

        DO iGF = iGF_Gm_dd_11, iGF_Beta_3

          G_K(iGF) = SUM( WeightsX_q * G(:,iX1,iX2,iX3,iGF) )

        END DO

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 1, G_K, U_M(1,:,iX1  ,iX2,iX3), &
                 R_X1(:,:,0), invR_X1(:,:,0) )

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 1, G_K, U_M(1,:,iX1-1,iX2,iX3), &
                 R_X1(:,:,1), invR_X1(:,:,1) )

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 1, G_K, U_M(1,:,iX1+1,iX2,iX3), &
                 R_X1(:,:,2), invR_X1(:,:,2) )

        ! --- Step 2.2: Project into characteristic fields ---

        DO iNodeX = 1, nDOFX ! Loop over modes

          U_M(iNodeX,:,iX1  ,iX2,iX3) &
            = MATMUL( invR_X1(:,:,0), U_M(iNodeX,:,iX1  ,iX2,iX3) )

          U_M(iNodeX,:,iX1-1,iX2,iX3) &
            = MATMUL( invR_X1(:,:,1), U_M(iNodeX,:,iX1-1,iX2,iX3) )

          U_M(iNodeX,:,iX1+1,iX2,iX3) &
            = MATMUL( invR_X1(:,:,2), U_M(iNodeX,:,iX1+1,iX2,iX3) )

        END DO

      END IF

      ! --- Step 1.1: Compute 'q' polynomials ---

      DO iCF = 1, nCF

        DO iNodeX = 1, nDOFX ! Loop over modes

          DO jNodeX = 1, nDOFX ! Loop over grid points

            q(jNodeX,iNodeX,iCF) &
              = SUM( U_M(1:iNodeX,iCF,iX1,iX2,iX3) &
                       * OrthonormalBasis(jNodeX,1:iNodeX,1) )

          END DO

        END DO

      END DO

      pCoeffs = Zero

      ! --- 2nd order ---

      ! --- Step 1.2.1 ---

      DO iCF = 1, nCF

        CALL ComputeCoefficients_Order2 &
               ( LinearWeights(1:2), pCoeffs(1:2,1:2,iCF) )

        DO iNodeX = 1, nDOFX ! Loop over grid points

          p(iNodeX,1,iCF) &
            = DOT_PRODUCT( pCoeffs(1,1:2,iCF), q(iNodeX,1:2,iCF) )

          p(iNodeX,2,iCF) &
            = DOT_PRODUCT( pCoeffs(2,1:2,iCF), q(iNodeX,1:2,iCF) )

        END DO

        ! --- Step 1.3: Compute smoothness indicators ---

        ! --- Hard-code smoothness indicator beta_01 ---

        SmoothnessIndicators(1,iCF) &
          = 12.0_DP &
              * MIN( U_M(2,iCF,iX1-1,iX2,iX3)**2, U_M(2,iCF,iX1+1,iX2,iX3)**2 )

        CALL ComputeSmoothnessIndicator_Order2 &
               ( pCoeffs(1:2,1:2,iCF), &
                 OrthonormalBasis, &
                 U_M(1:2,iCF,iX1,iX2,iX3), &
                 SmoothnessIndicators(2,iCF) )

        ! --- Step 1.4: Compute nonlinear weights ---

        CALL ComputeNonLinearWeights &
               ( SmoothnessIndicators(1:2,iCF), LinearWeights(1:2), &
                 NonLinearWeights(1:2,iCF) )

      END DO ! End of loop over fields

      ! --- 3rd order ---

      IF( nNodesX(1) .GT. 2 )THEN

        DO iCF = 1, nCF

          ! --- Step 1.2.2 ---

          CALL ComputeCoefficients_Order3 &
                 ( LinearWeights(3:4), NonLinearWeights(1:2,iCF), &
                   pCoeffs(:,:,iCF) )

          DO iNodeX = 1, nDOFX ! Loop over grid points

            p(iNodeX,3,iCF) &
              = DOT_PRODUCT( pCoeffs(3,1:3,iCF), q(iNodeX,1:3,iCF) )

            p(iNodeX,4,iCF) &
              = DOT_PRODUCT( pCoeffs(4,1:3,iCF), q(iNodeX,1:3,iCF) )

          END DO

          ! --- Step 1.3: Compute smoothness indicators ---

          CALL ComputeSmoothnessIndicators_Order3 &
                 ( pCoeffs(:,:,iCF), OrthonormalBasis, U_M(:,iCF,iX1,iX2,iX3), &
                   SmoothnessIndicators(3:4,iCF) )

          ! --- Step 1.4: Compute nonlinear weights ---

          CALL ComputeNonLinearWeights &
                 ( SmoothnessIndicators(3:4,iCF), LinearWeights(3:4), &
                   NonLinearWeights(3:4,iCF) )

        END DO

      END IF

      ! --- Step 1.5: Reconstruct solution polynomial ---

      DO iCF = 1, nCF

        DO iNodeX = 1, nDOFX ! Loop over grid points

          UU(iNodeX,iX1,iX2,iX3,iCF) &
            = SUM( NonLinearWeights(2*(nNodesX(1)-1)-1:2*(nNodesX(1)-1),iCF) &
                     * p(iNodeX,2*(nNodesX(1)-1)-1:2*(nNodesX(1)-1),iCF) )

        END DO

      END DO

      ! --- Step 2.4: Project new polynomial back into physical space ---

      IF( UseCharacteristicLimiting )THEN

        DO iNodeX = 1, nDOFX

          UU(iNodeX,iX1,iX2,iX3,:) &
            = MATMUL( R_X1(:,:,0), UU(iNodeX,iX1,iX2,iX3,:) )

        END DO

      END IF

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( LimitedCell(iX1,iX2,iX3) )THEN

        DO iCF = 1, nCF

          U(:,iX1,iX2,iX3,iCF) = UU(:,iX1,iX2,iX3,iCF)

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO


  SUBROUTINE ComputeSmoothnessIndicator_Order2 &
    ( pCoeffs, OrthonormalBasis, U_M, &
      SmoothnessIndicator )

    REAL(DP), INTENT(in)  :: pCoeffs(2,2)
    REAL(DP), INTENT(in)  :: OrthonormalBasis(nDOFX,nDOFX,nDOFX)
    REAL(DP), INTENT(in)  :: U_M(2)
    REAL(DP), INTENT(out) :: SmoothnessIndicator

    INTEGER  :: iNodeX
    REAL(DP) :: dPx(nDOFX), dq(nDOFX,2)

    DO iNodeX = 1, nDOFX ! Loop over quadrature points

      dq(iNodeX,1) = U_M(1) * OrthonormalBasis(iNodeX,1,2)

      dq(iNodeX,2) = SUM( U_M * OrthonormalBasis(iNodeX,1:2,2) )

      dPx(iNodeX) = DOT_PRODUCT( pCoeffs(2,:), dq(iNodeX,:) )

    END DO

    SmoothnessIndicator = SUM( WeightsX_q * dPx**2 )

  END SUBROUTINE ComputeSmoothnessIndicator_Order2


  SUBROUTINE ComputeSmoothnessIndicators_Order3 &
    ( pCoeffs, OrthonormalBasis, U_M, &
      SmoothnessIndicators )

    REAL(DP), INTENT(in)  :: pCoeffs(:,:)
    REAL(DP), INTENT(in)  :: OrthonormalBasis(nDOFX,nDOFX,nDOFX)
    REAL(DP), INTENT(in)  :: U_M(nDOFX)
    REAL(DP), INTENT(out) :: SmoothnessIndicators(2)

    INTEGER  :: iNodeX
    REAL(DP) :: dPx(nDOFX), dq(nDOFX,nNodesX(1))
    REAL(DP) :: ddPx(nDOFX,2), ddq(nDOFX,nNodesX(1))

    DO iNodeX = 1, nDOFX ! Loop over quadrature points

      ! --- 1st derivatives ---

      dq(iNodeX,1) = U_M(1) * OrthonormalBasis(iNodeX,1,2)
      dq(iNodeX,2) = SUM( U_M(1:2) * OrthonormalBasis(iNodeX,1:2,2) )
      dq(iNodeX,3) = SUM( U_M(1:3) * OrthonormalBasis(iNodeX,1:3,2) )

      dPx(iNodeX) = DOT_PRODUCT( pCoeffs(3,:), dq(iNodeX,:) )

      ! --- 2nd derivatives ---

      ddq(iNodeX,1) = U_M(1) * OrthonormalBasis(iNodeX,1,3)
      ddq(iNodeX,2) = SUM( U_M(1:2) * OrthonormalBasis(iNodeX,1:2,3) )
      ddq(iNodeX,3) = SUM( U_M(1:3) * OrthonormalBasis(iNodeX,1:3,3) )

      ddPx(iNodeX,1) = DOT_PRODUCT( pCoeffs(3,:), ddq(iNodeX,:) )
      ddPx(iNodeX,2) = DOT_PRODUCT( pCoeffs(4,:), ddq(iNodeX,:) )

    END DO

    SmoothnessIndicators(1) = SUM( WeightsX_q * dPx**2 )

    SmoothnessIndicators(2) &
      = SUM( WeightsX_q * ddPx(:,1)**2 ) + SUM( WeightsX_q * ddPx(:,2)**2 )

  END SUBROUTINE ComputeSmoothnessIndicators_Order3


  SUBROUTINE ComputeCoefficients_Order2 &
    ( LinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(2)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(1,1) = One

    PolynomialCoefficients(1,2) = Zero

    PolynomialCoefficients(2,1) = -LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(2,2) = One / LinearWeights(2)

  END SUBROUTINE ComputeCoefficients_Order2


  SUBROUTINE ComputeCoefficients_Order3 &
    ( LinearWeights, NonLinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(2), NonLinearWeights(:)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(3,1) &
      = NonLinearWeights(1) &
          - NonLinearWeights(2) * LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(3,2) &
      = NonLinearWeights(2) / LinearWeights(2)

    PolynomialCoefficients(3,3) &
      = Zero

    PolynomialCoefficients(4,1) &
      = -LinearWeights(1) / LinearWeights(2) &
          * ( NonLinearWeights(1) &
                - NonLinearWeights(2) * LinearWeights(1) / LinearWeights(2) )

    PolynomialCoefficients(4,2) &
      = -LinearWeights(1) / LinearWeights(2) &
          *  NonLinearWeights(2) / LinearWeights(2)

    PolynomialCoefficients(4,3) &
      = One / LinearWeights(2)

  END SUBROUTINE ComputeCoefficients_Order3


  SUBROUTINE ComputeNonLinearWeights &
    ( SmoothnessIndicators, LinearWeights, NonLinearWeights )

    REAL(DP), INTENT(in)  :: SmoothnessIndicators(2), LinearWeights(2)
    REAL(DP), INTENT(out) :: NonLinearWeights(2)

    REAL(DP) :: Tau, Normalization, EPS
    EPS = 1.0d-10

    Tau = ( SmoothnessIndicators(2) - SmoothnessIndicators(1) )**2

    Normalization &
      = SUM( LinearWeights * ( One + Tau / ( EPS + SmoothnessIndicators ) ) )

    NonLinearWeights(1) &
      = LinearWeights(1) * ( One + Tau / ( EPS + SmoothnessIndicators(1) ) ) &
          / Normalization

    NonLinearWeights(2) &
      = LinearWeights(2) * ( One + Tau / ( EPS + SmoothnessIndicators(2) ) ) &
          / Normalization

  END SUBROUTINE ComputeNonLinearWeights


  SUBROUTINE MapNodalToModal_Fluid_WENO( VandermondeMatrix, uN, uM )

    REAL(DP), INTENT(in)  :: VandermondeMatrix(nDOFX,nDOFX)
    REAL(DP), INTENT(in)  :: uN(nDOFX)
    REAL(DP), INTENT(out) :: uM(nDOFX)

    uM = MATMUL( VandermondeMatrix, uN )

  END SUBROUTINE MapNodalToModal_Fluid_WENO


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


  SUBROUTINE DetectTroubledCells( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: V_K (0:2*nDimsX)
    REAL(DP) :: U_K (0:2*nDimsX,nCF)
    REAL(DP) :: U_K0(0:2*nDimsX,nCF)

    D(:,:,:,:,iDF_Sh) = Zero

    IF( .NOT. UseTroubledCellIndicator )THEN

      D(:,:,:,:,iDF_Sh) = 1.1_DP * LimiterThreshold
      RETURN

    END IF

    CALL TimersStart_Euler( Timer_Euler_TroubledCellIndicator )

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

      D(:,iX1,iX2,iX3,iDF_Sh) &
        = SUM( ABS( U_K(0,iCF_D) - U_K0(1:2*nDimsX,iCF_D) ) ) &
            / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_D) ) )

      ! --- Use Conserved Energy  to Detect Troubled Cell ---

      D(:,iX1,iX2,iX3,iDF_Sh) &
        = MAX( MAXVAL(D(:,iX1,iX2,iX3,iDF_Sh) ), &
               SUM( ABS( U_K(0,iCF_E) - U_K0(1:2*nDimsX,iCF_E) ) ) &
                 / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_E) ) ) )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_TroubledCellIndicator )

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
