MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL

  USE KindModule,                                                 ONLY: &
    DP,   &
    Zero, &
    One,  &
    Two
  USE ProgramHeaderModule,                                        ONLY: &
    nDOFX,  &
    nDimsX, &
    nNodes, &
    nNodesX
  USE ReferenceElementModuleX,                                    ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule,                                            ONLY: &
    MinModB, &
    NodeNumberX
  USE PolynomialBasisModule_Legendre,                             ONLY: &
    P_X1, &
    P_X2, &
    P_X3, &
    IndPX_Q
  USE PolynomialBasisMappingModule,                               ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE MeshModule,                                                 ONLY: &
    MeshX
  USE GeometryFieldsModule,                                       ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_SqrtGm
  USE FluidFieldsModule,                                          ONLY: &
    nCF,   &
    iCF_D, &
    iCF_E, &
    iDF_TCI
  USE Euler_BoundaryConditionsModule,                             ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL, ONLY: &
    ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL
  USE Euler_DiscontinuityDetectionModule,                         ONLY: &
    InitializeTroubledCellIndicator_Euler, &
    FinalizeTroubledCellIndicator_Euler,   &
    DetectTroubledCells_Euler,             &
    DetectShocks_Euler
  USE TimersModule_Euler,                                         ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_SlopeLimiter

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplySlopeLimiter_Euler_Relativistic_IDEAL

  REAL(DP), PUBLIC :: LimiterThreshold

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
  REAL(DP) :: LimiterThresholdParameter
  REAL(DP) :: I_6x6(1:6,1:6)

  LOGICAL :: DEBUG = .FALSE.

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
    ( UseSlopeLimiter_Option,           &
      SlopeLimiterMethod_Option,        &
      BetaTVD_Option,                   &
      BetaTVB_Option,                   &
      SlopeTolerance_Option,            &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option,  &
      LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in),     OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option, &
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

    UseSlopeLimiter = .TRUE.
    IF( PRESENT( UseSlopeLimiter_Option ) ) &
      UseSlopeLimiter = UseSlopeLimiter_Option

    SlopeLimiterMethod = 'TVD'
    IF( PRESENT( SlopeLimiterMethod_Option ) ) &
      SlopeLimiterMethod = SlopeLimiterMethod_Option

    BetaTVD = One
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    BetaTVB = Zero
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    SlopeTolerance = 1.0d-3
    IF( PRESENT( SlopeTolerance_Option ) ) &
      SlopeTolerance = SlopeTolerance_Option

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
      WRITE(*,'(A)') &
        '    INFO: InitializeSlopeLimiter_Euler_Relativistic_IDEAL'
      WRITE(*,'(A)') &
        '    -----------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,A)')         '', 'SlopeLimiterMethod: ', &
        TRIM( SlopeLimiterMethod )
      WRITE(*,*)

      IF( TRIM( SlopeLimiterMethod ) .EQ. 'TVD' )THEN

        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVD: ' , &
          BetaTVD
        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVB: ' , &
          BetaTVB
        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'SlopeTolerance: ' , &
          SlopeTolerance
        WRITE(*,*)

      END IF

      WRITE(*,'(A4,A27,L1)'       ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'LimiterThreshold: ' , &
        LimiterThreshold
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseConservativeCorrection: ' , &
        UseConservativeCorrection
    END IF

    CALL InitializeTroubledCellIndicator_Euler &
           ( UseTroubledCellIndicator_Option = UseTroubledCellIndicator, &
             LimiterThreshold_Option = LimiterThreshold )

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

        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

      CASE( 'WENO' )

        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

      CASE DEFAULT

        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    END SELECT

    CALL DetectShocks_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
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

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) ) CYCLE

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
              .GE. SlopeTolerance * ABS( U_M(iCF,0,1,iX1,iX2,iX3) ) )THEN

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

    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iGridPt
    REAL(DP) :: OrthonormalBasis1D(1:nNodesX(1),0:nNodesX(1)-1,0:nNodesX(1)-1)
    REAL(DP) :: eta

    IF( nDimsX .GT. 2 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nDimsX > 2. Stopping...'
      STOP

    END  IF

    IF( nNodesX(1) .GT. 3 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nNodesX(1) > 3. Stopping...'
      STOP

    END  IF

    IF( nNodesX(1) .EQ. 2 )THEN

      DO iNodeX = 1, nNodesX(1)

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis1D(iNodeX,0,0) &
          = One
        OrthonormalBasis1D(iNodeX,1,0) &
          = SQRT(12.0_DP) * eta

        ! --- 1st derivative ---
        OrthonormalBasis1D(iNodeX,0,1) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,1) &
          = SQRT(12.0_DP)

      END DO

    ELSE IF( nNodesX(1) .EQ. 3 )THEN

      DO iNodeX = 1, nNodesX(1)

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis1D(iNodeX,0,0) &
          = One
        OrthonormalBasis1D(iNodeX,1,0) &
          = SQRT(12.0_DP) * eta
        OrthonormalBasis1D(iNodeX,2,0) &
          = SQRT(180.0_DP) * ( eta**2 - One / 12.0_DP )

        ! --- 1st derivative ---
        OrthonormalBasis1D(iNodeX,0,1) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,1) &
          = SQRT(12.0_DP)
        OrthonormalBasis1D(iNodeX,2,1) &
          = SQRT(180.0_DP) * Two * eta

        ! --- 2nd derivative ---
        OrthonormalBasis1D(iNodeX,0,2) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,2) &
          = Zero
        OrthonormalBasis1D(iNodeX,2,2) &
          = SQRT(180.0_DP) * Two

      END DO

    END IF

    ALLOCATE( VandermondeMatrix(nDOFX,nDOFX) )

    IF( nDimsX .EQ. 1 )THEN

      ALLOCATE( OrthonormalBasis(1:nNodesX(1),0:nNodesX(1)-1,0:nNodesX(1)-1) )

      OrthonormalBasis = OrthonormalBasis1D

    ELSE IF( nDimsX .EQ. 2 )THEN

      ALLOCATE( OrthonormalBasis(1:nDOFX, &
                                 0:nDOFX-1, &
                                 0:nNodesX(1)*(nNodesX(1)+1)/2-1) )

      OrthonormalBasis = Zero

      IF     ( nNodesX(1) .EQ. 2 )THEN

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          ! --- 0th derivative ---
          OrthonormalBasis(iNodeX,0,0) = OrthonormalBasis1D(iNodeX1,0,0)
          OrthonormalBasis(iNodeX,1,0) = OrthonormalBasis1D(iNodeX1,1,0)
          OrthonormalBasis(iNodeX,2,0) = OrthonormalBasis1D(iNodeX2,1,0)

          ! --- 1st eta1 derivative ---
          OrthonormalBasis(iNodeX,0,1) = OrthonormalBasis1D(iNodeX1,0,1)
          OrthonormalBasis(iNodeX,1,1) = OrthonormalBasis1D(iNodeX1,1,1)
          OrthonormalBasis(iNodeX,2,1) = Zero

          ! --- 1st eta2 derivative ---
          OrthonormalBasis(iNodeX,0,2) = OrthonormalBasis1D(iNodeX2,0,1)
          OrthonormalBasis(iNodeX,1,2) = Zero
          OrthonormalBasis(iNodeX,2,2) = OrthonormalBasis1D(iNodeX2,1,1)

        END DO

      ELSE IF( nNodesX(1) .EQ. 3 )THEN

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          ! --- 0th derivative ---
          OrthonormalBasis(iNodeX,0,0) = OrthonormalBasis1D(iNodeX1,0,0)
          OrthonormalBasis(iNodeX,1,0) = OrthonormalBasis1D(iNodeX1,1,0)
          OrthonormalBasis(iNodeX,2,0) = OrthonormalBasis1D(iNodeX2,1,0)
          OrthonormalBasis(iNodeX,3,0) = OrthonormalBasis1D(iNodeX1,2,0)
          OrthonormalBasis(iNodeX,4,0) = OrthonormalBasis1D(iNodeX1,1,0) &
                                           * OrthonormalBasis1D(iNodeX2,1,0)
          OrthonormalBasis(iNodeX,5,0) = OrthonormalBasis1D(iNodeX2,2,0)

          ! --- 1st eta1 derivative ---
          OrthonormalBasis(iNodeX,0,1) = OrthonormalBasis1D(iNodeX1,0,1)
          OrthonormalBasis(iNodeX,1,1) = OrthonormalBasis1D(iNodeX1,1,1)
          OrthonormalBasis(iNodeX,2,1) = Zero
          OrthonormalBasis(iNodeX,3,1) = OrthonormalBasis1D(iNodeX1,2,1)
          OrthonormalBasis(iNodeX,4,1) = OrthonormalBasis1D(iNodeX1,1,1) &
                                           * OrthonormalBasis1D(iNodeX2,1,0)

          OrthonormalBasis(iNodeX,5,1) = Zero

          ! --- 1st eta2 derivative ---
          OrthonormalBasis(iNodeX,0,2) = OrthonormalBasis1D(iNodeX2,0,1)
          OrthonormalBasis(iNodeX,1,2) = Zero
          OrthonormalBasis(iNodeX,2,2) = OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,3,2) = Zero
          OrthonormalBasis(iNodeX,4,2) = OrthonormalBasis1D(iNodeX1,1,0) &
                                           * OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,5,2) = OrthonormalBasis1D(iNodeX2,2,1)

          ! --- 2nd eta1 derivative ---
          OrthonormalBasis(iNodeX,0,3) = OrthonormalBasis1D(iNodeX1,0,2)
          OrthonormalBasis(iNodeX,1,3) = OrthonormalBasis1D(iNodeX1,1,2)
          OrthonormalBasis(iNodeX,2,3) = Zero
          OrthonormalBasis(iNodeX,3,3) = OrthonormalBasis1D(iNodeX1,2,2)
          OrthonormalBasis(iNodeX,4,3) = Zero
          OrthonormalBasis(iNodeX,5,3) = Zero

          ! --- 1st mixed derivative ---
          OrthonormalBasis(iNodeX,0,4) = Zero
          OrthonormalBasis(iNodeX,1,4) = Zero
          OrthonormalBasis(iNodeX,2,4) = Zero
          OrthonormalBasis(iNodeX,3,4) = Zero
          OrthonormalBasis(iNodeX,4,4) = OrthonormalBasis1D(iNodeX1,1,1) &
                                           * OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,5,4) = Zero

          ! --- 2nd eta2 derivative ---
          OrthonormalBasis(iNodeX,0,5) = OrthonormalBasis1D(iNodeX2,0,2)
          OrthonormalBasis(iNodeX,1,5) = Zero
          OrthonormalBasis(iNodeX,2,5) = OrthonormalBasis1D(iNodeX2,1,2)
          OrthonormalBasis(iNodeX,3,5) = Zero
          OrthonormalBasis(iNodeX,4,5) = Zero
          OrthonormalBasis(iNodeX,5,5) = OrthonormalBasis1D(iNodeX2,2,2)

        END DO

      END IF ! Order of accuracy

    END IF ! Spatial dimensions

    DO iNodeX = 1, nDOFX ! Polynomial

      DO iGridPt = 1, nDOFX

        VandermondeMatrix(iNodeX,iGridPt) &
          = WeightsX_q(iGridPt) * OrthonormalBasis(iGridPt,iNodeX-1,0)

      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter_Euler_WENO


  SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

    IF( UseTroubledCellIndicator ) &
      CALL FinalizeTroubledCellIndicator_Euler

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
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    ! --- Currently assuming dX1 = dX2 = dX3 ---

    LOGICAL  :: SuppressBC
    INTEGER  :: iX1, iX2, iX3, iNodeX, iGridPt, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: U_M(0:nDOFX-1,nCF,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))
    REAL(DP) :: UU(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),nCF)
    LOGICAL  :: LimitedCell(iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    ! --- WENO Limiter ---
    REAL(DP) :: q  (nDOFX,0:nNodesX(1)-1,nCF)
    REAL(DP) :: qX1(nDOFX,0:nNodesX(1)-1,nCF)
    REAL(DP) :: qX2(nDOFX,0:nNodesX(1)-1,nCF)
    REAL(DP) :: LinearWeights(2*(nNodesX(1)-1))
    REAL(DP) :: NonLinearWeights(2*(nNodesX(1)-1),nCF)
    REAL(DP) :: SmoothnessIndicators(2*(nNodesX(1)-1),nCF)
    REAL(DP) :: pX1(nDOFX,2*(nNodesX(1)-1),nCF)
    REAL(DP) :: pX2(nDOFX,2*(nNodesX(1)-1),nCF)
    REAL(DP) :: pX1_new(nDOFX,nCF)
    REAL(DP) :: pX2_new(nDOFX,nCF)
    REAL(DP) :: pCoeffs(2*(nNodesX(1)-1),nNodesX(1),nCF)
    INTEGER  :: nPspace, k

    ! --- Characteristic limiting ---
    INTEGER  :: iGF
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)
    REAL(DP) :: R_X2(nCF,nCF), invR_X2(nCF,nCF)
    REAL(DP) :: G_K(nGF)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    k = nNodesX(1) - 1

    ! --- Get linear weights ---
    DO iNodeX = 1, 2*k

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

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) ) CYCLE

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

        IF( nDimsX .GT. 1 )THEN

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2-1,iX3,iCF), &
                   U_M(:,iCF,iX1,iX2-1,iX3) )

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2+1,iX3,iCF), &
                   U_M(:,iCF,iX1,iX2+1,iX3) )

        END IF

      END DO

      IF( UseCharacteristicLimiting )THEN

        DO iGF = iGF_Gm_dd_11, iGF_Beta_3

          G_K(iGF) = SUM( WeightsX_q * G(:,iX1,iX2,iX3,iGF) )

        END DO

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 1, G_K, U_M(0,:,iX1,iX2,iX3), R_X1, invR_X1 )

        IF( nDimsX .GT. 1 ) &
          CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
                 ( 2, G_K, U_M(0,:,iX1,iX2,iX3), R_X2, invR_X2 )

      END IF

      DO iCF = 1, nCF

        DO iNodeX = 0, nNodesX(1)-1 ! Loop over basis polynomials

          nPspace = UpperLimit_pSpace( iNodeX )

          DO iGridPt = 1, nDOFX

            q(iGridPt,iNodeX,iCF) &
              = SUM( U_M(0:nPspace,iCF,iX1,iX2,iX3) &
                       * OrthonormalBasis(iGridPt,0:nPspace,0) )

          END DO

        END DO

      END DO

      IF( UseCharacteristicLimiting )THEN

        DO iNodeX = 0, nNodesX(1)-1

          DO iGridPt = 1, nDOFX

            qX1(iGridPt,iNodeX,:) = MATMUL( invR_X1, q(iGridPt,iNodeX,:) )

          END DO

        END DO

        IF( nDimsX .GT. 1 )THEN

          DO iNodeX = 0, nNodesX(1)-1

            DO iGridPt = 1, nDOFX

              qX2(iGridPt,iNodeX,:) = MATMUL( invR_X2, q(iGridPt,iNodeX,:) )

            END DO

          END DO

        END IF

      ELSE

        qX1 = q

      END IF

      ! --- 2nd order ---

      DO iCF = 1, nCF

        CALL ComputeCoefficients_Order2 &
               ( LinearWeights, pCoeffs(:,:,iCF) )

        DO iGridPt = 1, nDOFX

          pX1(iGridPt,1:2,iCF) &
            = MATMUL( pCoeffs(1:2,1:2,iCF), qX1(iGridPt,0:1,iCF) )

        END DO

      END DO

      IF( UseCharacteristicLimiting .AND. nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          DO iGridPt = 1, nDOFX

            pX2(iGridPt,1:2,iCF) &
              = MATMUL( pCoeffs(1:2,1:2,iCF), qX2(iGridPt,0:1,iCF) )

          END DO

        END DO

      END IF

      DO iCF = 1, nCF

        ! --- Hard-code smoothness indicator beta_01 ---

        IF     ( nDimsX .EQ. 1 )THEN

          SmoothnessIndicators(1,iCF) &
            = OrthonormalBasis(1,1,1)**2 &
                * MIN( U_M(1,iCF,iX1-1,iX2,iX3)**2, &
                       U_M(1,iCF,iX1+1,iX2,iX3)**2 )

        ELSE IF( nDimsX .EQ. 2 )THEN

          ! --- Ignoring mesh size ratios because ignoring SqrtGm ---

          SmoothnessIndicators(1,iCF) &
            = OrthonormalBasis(1,1,1)**2 &
                * MIN( U_M(1,iCF,iX1-1,iX2,iX3)**2 &
                         + U_M(2,iCF,iX1-1,iX2,iX3)**2, &
                       U_M(1,iCF,iX1+1,iX2,iX3)**2 &
                         + U_M(2,iCF,iX1+1,iX2,iX3)**2, &
                       U_M(1,iCF,iX1,iX2-1,iX3)**2 &
                         + U_M(2,iCF,iX1,iX2-1,iX3)**2, &
                       U_M(1,iCF,iX1,iX2+1,iX3)**2 &
                         + U_M(2,iCF,iX1,iX2+1,iX3)**2 )

        END IF

        CALL ComputeSmoothnessIndicator_Order2 &
               ( LinearWeights(2), &
                 OrthonormalBasis(1,1,1), &
                 U_M(1:nDimsX,iCF,iX1,iX2,iX3), &
                 SmoothnessIndicators(2,iCF) )

        CALL ComputeNonLinearWeights &
               ( SmoothnessIndicators(1:2,iCF), LinearWeights(1:2), &
                 NonLinearWeights(1:2,iCF) )

      END DO ! End of loop over fields

      ! --- 3rd order ---

      IF( nNodesX(1) .GT. 2 )THEN

        DO iCF = 1, nCF

          CALL ComputeCoefficients_Order3 &
                 ( LinearWeights, NonLinearWeights(1:2,iCF), &
                   pCoeffs(:,:,iCF) )

          DO iGridPt = 1, nDOFX

            pX1(iGridPt,3:4,iCF) &
              = MATMUL( pCoeffs(3:4,1:3,iCF), qX1(iGridPt,0:2,iCF) )

          END DO

        END DO

        IF( UseCharacteristicLimiting .AND. nDimsX .GT. 1 )THEN

          DO iCF = 1, nCF

            DO iGridPt = 1, nDOFX

              pX2(iGridPt,3:4,iCF) &
                = MATMUL( pCoeffs(3:4,1:3,iCF), qX2(iGridPt,0:2,iCF) )

            END DO

          END DO

        END IF

        DO iCF = 1, nCF

          CALL ComputeSmoothnessIndicators_Order3 &
                 ( LinearWeights(2), NonLinearWeights(2,iCF), &
                   U_M(:,iCF,iX1,iX2,iX3), &
                   pCoeffs(4,:,iCF), SmoothnessIndicators(:,iCF) )

          CALL ComputeNonLinearWeights &
                 ( SmoothnessIndicators(3:4,iCF), LinearWeights(3:4), &
                   NonLinearWeights(3:4,iCF) )

        END DO

      END IF

      DO iCF = 1, nCF

        DO iGridPt = 1, nDOFX

          pX1_New(iGridPt,iCF) &
            = SUM( NonLinearWeights(2*k-1:2*k,iCF) &
                     * pX1(iGridPt,2*k-1:2*k,iCF) )

        END DO

      END DO

      IF( UseCharacteristicLimiting .AND. nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          DO iGridPt = 1, nDOFX

            pX2_New(iGridPt,iCF) &
              = SUM( NonLinearWeights(2*k-1:2*k,iCF) &
                       * pX2(iGridPt,2*k-1:2*k,iCF) )

          END DO

        END DO

      END IF

      IF( UseCharacteristicLimiting )THEN

        DO iGridPt = 1, nDOFX

          pX1_New(iGridPt,:) = MATMUL( R_X1, pX1_New(iGridPt,:) )

        END DO

        IF( nDimsX .GT. 1 )THEN

          DO iGridPt = 1, nDOFX

            pX2_New(iGridPt,:) = MATMUL( R_X2, pX2_New(iGridPt,:) )

          END DO

        ELSE

          pX2_New = Zero

        END IF

      ELSE ! Component-Wise Limiting

        IF( nDimsX .GT. 1 )THEN

          pX2_New = pX1_New

        ELSE

          pX2_New = Zero

        END IF

      END IF

      DO iCF = 1, nCF

        UU(:,iX1,iX2,iX3,iCF) &
          = ( pX1_New(:,iCF) + pX2_New(:,iCF) ) / DBLE( nDimsX )

      END DO

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


  INTEGER FUNCTION UpperLimit_pSpace( ell )

    INTEGER, INTENT(in) :: ell

    IF( nDimsX .EQ. 1 )THEN

      UpperLimit_pSpace = ell

    ELSE IF( nDimsX .EQ. 2 )THEN

      UpperLimit_pSpace = ( ell + 1 ) * ( ell + 2 ) / 2 - 1

    END IF

    RETURN
  END FUNCTION UpperLimit_pSpace

  SUBROUTINE ComputeSmoothnessIndicator_Order2 &
    ( LinearWeight_11, Coeff_v1, U_M, SmoothnessIndicator )

    REAL(DP), INTENT(in)  :: LinearWeight_11, Coeff_v1, U_M(:)
    REAL(DP), INTENT(out) :: SmoothnessIndicator

    ! --- beta_11 ---

    SmoothnessIndicator = ( U_M(1) / LinearWeight_11 * Coeff_v1 )**2

    IF( nDimsX .GT. 1 ) &
      SmoothnessIndicator &
        = SmoothnessIndicator + ( U_M(2) / LinearWeight_11 * Coeff_v1 )**2

  END SUBROUTINE ComputeSmoothnessIndicator_Order2


  SUBROUTINE ComputeSmoothnessIndicators_Order3 &
    ( LinearWeight_11, NonLinearWeight_11, U_M, &
      pCoeffs, SmoothnessIndicators )

    REAL(DP), INTENT(in)  :: LinearWeight_11, NonLinearWeight_11
    REAL(DP), INTENT(in)  :: U_M(0:nDOFX-1), pCoeffs(:)
    REAL(DP), INTENT(out) :: SmoothnessIndicators(:)

    REAL(DP) :: Coeff_v1
    REAL(DP) :: dq_X1   (nDOFX,0:2), &
                dq_X2   (nDOFX,0:2), &
                ddq_X1  (nDOFX,0:2), &
                ddq_X2  (nDOFX,0:2), &
                ddq_X1X2(nDOFX,0:2)
    REAL(DP) :: dP_X1   (nDOFX), &
                dP_X2   (nDOFX), &
                ddP_X1  (nDOFX), &
                ddP_X2  (nDOFX), &
                ddP_X1X2(nDOFX)
    INTEGER  :: iGridPt

    ! --- beta_12 ---

    Coeff_v1 = OrthonormalBasis(1,1,1)

    SmoothnessIndicators(3) &
      = ( Coeff_v1 * NonLinearWeight_11 / LinearWeight_11 * U_M(1) )**2

    IF( nDimsX .GT. 1 ) &
      SmoothnessIndicators(3) &
        = SmoothnessIndicators(3) &
            + ( Coeff_v1 * NonLinearWeight_11 / LinearWeight_11 * U_M(2) )**2

    ! --- beta_22 ---

    IF( nDimsX .EQ. 1 )THEN

      DO iGridPt = 1, nDOFX

        dq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,1) )
        dq_X1(iGridPt,1) = SUM( U_M(0:1) * OrthonormalBasis(iGridPt,0:1,1) )
        dq_X1(iGridPt,2) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,1) )
        dp_X1 (iGridPt)  = DOT_PRODUCT( pCoeffs, dq_X1(iGridPt,:) )

        ddq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,2) )
        ddq_X1(iGridPt,1) = SUM( U_M(0:1) * OrthonormalBasis(iGridPt,0:1,2) )
        ddq_X1(iGridPt,2) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,2) )
        ddp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1(iGridPt,:) )

      END DO

      SmoothnessIndicators(4) &
        = SUM( WeightsX_q * ( dp_X1**2 + ddp_X1**2 ) )

    ELSE IF( nDimsX .GT. 1 )THEN

      DO iGridPt = 1, nDOFX

        dq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,1) )
        dq_X1(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,1) )
        dq_X1(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,1) )
        dp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, dq_X1(iGridPt,:) )

        dq_X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,2) )
        dq_X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,2) )
        dq_X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,2) )
        dp_X2(iGridPt)   = DOT_PRODUCT( pCoeffs, dq_X2(iGridPt,:) )

        ddq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,3) )
        ddq_X1(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,3) )
        ddq_X1(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,3) )
        ddp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1(iGridPt,:) )

        ddq_X1X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,4) )
        ddq_X1X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,4) )
        ddq_X1X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,4) )
        ddp_X1X2(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1X2(iGridPt,:) )

        ddq_X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,5) )
        ddq_X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,5) )
        ddq_X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,5) )
        ddp_X2(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X2(iGridPt,:) )

      END DO

      SmoothnessIndicators(4) &
        = SUM( WeightsX_q &
                 * ( dp_X1**2 + ddp_X1**2 + dp_X2**2 &
                       + ddp_X2**2 + ddp_X1X2**2 ) )

    END IF

  END SUBROUTINE ComputeSmoothnessIndicators_Order3


  SUBROUTINE ComputeCoefficients_Order2 &
    ( LinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(:)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(1,1) = One

    PolynomialCoefficients(1,2) = Zero

    PolynomialCoefficients(2,1) = -LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(2,2) = One / LinearWeights(2)

  END SUBROUTINE ComputeCoefficients_Order2


  SUBROUTINE ComputeCoefficients_Order3 &
    ( LinearWeights, NonLinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(:), NonLinearWeights(:)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(1,1) &
      = One

    PolynomialCoefficients(1,2) &
      = Zero

    PolynomialCoefficients(1,3) &
      = Zero

    PolynomialCoefficients(2,1) &
      = -LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(2,2) &
      = One / LinearWeights(2)

    PolynomialCoefficients(2,3) &
      = Zero

    PolynomialCoefficients(3,1) &
      = NonLinearWeights(1) &
          - NonLinearWeights(2) * LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(3,2) &
      = NonLinearWeights(2) / LinearWeights(2)

    PolynomialCoefficients(3,3) &
      = Zero

    PolynomialCoefficients(4,1) &
      = -LinearWeights(3) / LinearWeights(4) * PolynomialCoefficients(3,1)

    PolynomialCoefficients(4,2) &
      = -LinearWeights(3) / LinearWeights(4) * PolynomialCoefficients(3,2)

    PolynomialCoefficients(4,3) &
      = One / LinearWeights(4)

  END SUBROUTINE ComputeCoefficients_Order3


  SUBROUTINE ComputeNonLinearWeights &
    ( SmoothnessIndicators, LinearWeights, NonLinearWeights )

    REAL(DP), INTENT(in)  :: SmoothnessIndicators(2), LinearWeights(2)
    REAL(DP), INTENT(out) :: NonLinearWeights(2)

    REAL(DP) :: Tau, Normalization
    REAL(DP), PARAMETER :: EPS = 1.0d-10

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
