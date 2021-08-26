MODULE Euler_SlopeLimiterModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX, &
    nNodes, &
    nNodesX, &
    bcX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE UtilitiesModule, ONLY: &
    MinModB, &
    NodeNumberX
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, &
    P_X2, &
    P_X3, &
    IndPX_Q
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_Ne, &
    iDF_TCI
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyInnerBC_Euler,  &
    ApplyOuterBC_Euler,  &
    iApplyBC_Euler_Both, &
    ApplyBoundaryConditions_Euler
  USE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE, ONLY: &
    ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    InitializeTroubledCellIndicator_Euler, &
    FinalizeTroubledCellIndicator_Euler, &
    LimiterThreshold, &
    DetectTroubledCells_Euler
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_SlopeLimiter

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: ApplySlopeLimiter_Euler_NonRelativistic_TABLE

  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseCharacteristicLimiting
  LOGICAL  :: UseConservativeCorrection
  LOGICAL  :: Verbose
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: I_6x6(1:6,1:6)


CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
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
      UseTroubledCellIndicator_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option

    INTEGER  :: i
    LOGICAL  :: UseTroubledCellIndicator
    REAL(DP) :: LimiterThresholdParameter

    IF( PRESENT( BetaTVD_Option ) )THEN
      BetaTVD = BetaTVD_Option
    ELSE
      BetaTVD = One
    END IF

    IF( PRESENT( BetaTVB_Option ) )THEN
      BetaTVB = BetaTVB_Option
    ELSE
      BetaTVB = Zero
    END IF

    IF( PRESENT( SlopeTolerance_Option ) )THEN
      SlopeTolerance = SlopeTolerance_Option
    ELSE
      SlopeTolerance = 1.0d-3
    END IF

    IF( PRESENT( UseSlopeLimiter_Option ) )THEN
      UseSlopeLimiter = UseSlopeLimiter_Option
    ELSE
      UseSlopeLimiter = .TRUE.
    END IF

    IF( PRESENT( UseCharacteristicLimiting_Option ) )THEN
      UseCharacteristicLimiting = UseCharacteristicLimiting_Option
    ELSE
      UseCharacteristicLimiting = .FALSE.
    END IF

    IF( PRESENT( UseTroubledCellIndicator_Option ) )THEN
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option
    ELSE
      UseTroubledCellIndicator = .TRUE.
    END IF

    IF( PRESENT( LimiterThresholdParameter_Option ) )THEN
      LimiterThresholdParameter = LimiterThresholdParameter_Option
    ELSE
      LimiterThresholdParameter = 0.03_DP
    END IF

    IF( PRESENT( UseConservativeCorrection_Option ) )THEN
      UseConservativeCorrection = UseConservativeCorrection_Option
    ELSE
      UseConservativeCorrection = .TRUE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    CALL InitializeTroubledCellIndicator_Euler &
           ( UseTroubledCellIndicator_Option = UseTroubledCellIndicator, &
             LimiterThresholdParameter_Option = LimiterThresholdParameter )

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: InitializeSlopeLimiter_Euler_NonRelativistic_TABLE:'
      WRITE(*,'(A)') '  ---------------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'BetaTVD: ' , &
        BetaTVD
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'BetaTVB: ' , &
        BetaTVB
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'SlopeTolerance: ' , &
        SlopeTolerance
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'LimiterThreshold: ' , &
        LimiterThreshold
    END IF

    I_6x6 = Zero
    DO i = 1, 6
      I_6x6(i,i) = One
    END DO

  END SUBROUTINE InitializeSlopeLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeTroubledCellIndicator_Euler

  END SUBROUTINE FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

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
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL  :: LimitedCell(nCF,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    LOGICAL  :: SuppressBC
    LOGICAL  :: LimitField(nCF)
    LOGICAL  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    INTEGER  :: iX1, iX2, iX3, iGF, iCF
    INTEGER  :: iApplyBC(3)
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: dY(nDimsX)
    REAL(DP) :: SlopeDifference(nCF)
    REAL(DP) :: a(nCF), b(nCF), c(nCF) ! --- Arguments for MinMod (fluid)
    REAL(DP) :: aY    , bY    , cY     ! --- Arguments for MinMod (Ye)
    REAL(DP) :: G_K(nGF)
    REAL(DP) :: dU (nCF,nDimsX)
    REAL(DP) :: U_M(nCF,0:2*nDimsX,nDOFX)
    REAL(DP) :: Y_N(    0:2*nDimsX,nDOFX) ! --- Nodal Ye
    REAL(DP) :: Y_M(    0:2*nDimsX,nDOFX) ! --- Modal Ye
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)
    REAL(DP) :: R_X2(nCF,nCF), invR_X2(nCF,nCF)
    REAL(DP) :: R_X3(nCF,nCF), invR_X3(nCF,nCF)
    REAL(DP) :: V_K(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    U_M = Zero

    iApplyBC = iApplyBC_Euler_Both
    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    LimitedCell = .FALSE.

    ExcludeInnerGhostCell = .FALSE.
    ExcludeOuterGhostCell = .FALSE.
    IF( ApplyInnerBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeInnerGhostCell(1) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeOuterGhostCell(1) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeInnerGhostCell(2) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeOuterGhostCell(2) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeInnerGhostCell(3) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeOuterGhostCell(3) = .TRUE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_TCI) < LimiterThreshold ) ) CYCLE

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      ! --- Cell Volume ---

      V_K(iX1,iX2,iX3) &
        = DOT_PRODUCT( WeightsX_q(:), G(:,iX1,iX2,iX3,iGF_SqrtGm) )

      ! --- Cell Average of Conserved Fluid ---

      DO iCF = 1, nCF

        U_K(iCF,iX1,iX2,iX3) &
          = SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                   * U(:,iX1,iX2,iX3,iCF) ) / V_K(iX1,iX2,iX3)

      END DO

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

      IF( UseCharacteristicLimiting )THEN

        ! --- Cell Average of Geometry (Spatial Metric Only) ---

        DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33

          G_K(iGF) = DOT_PRODUCT( WeightsX_q(:), G(:,iX1,iX2,iX3,iGF) )

        END DO

        ! --- Compute Eigenvectors ---

        CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
               ( 1, G_K(:), U_M(:,0,1), R_X1, invR_X1 )

        IF( nDimsX > 1 )THEN

          CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
                 ( 2, G_K(:), U_M(:,0,1), R_X2, invR_X2 )

        END IF

        IF( nDimsX > 2 )THEN

          CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
                 ( 3, G_K(:), U_M(:,0,1), R_X3, invR_X3 )

        END IF

      ELSE

        ! --- Componentwise Limiting ---

        R_X1 = I_6x6; invR_X1 = I_6x6
        R_X2 = I_6x6; invR_X2 = I_6x6
        R_X3 = I_6x6; invR_X3 = I_6x6

      END IF

      ! --- Compute Limited Slopes ---
      IF( UseCharacteristicLimiting ) &
        U_M(iCF_Ne,:,:) = U_M(iCF_Ne,:,:) * AtomicMassUnit

      a = MATMUL( invR_X1, U_M(:,0,2) )

      IF     ( iX1 .EQ. iX_B0(1) .AND. ExcludeInnerGhostCell(1) )THEN

        c = BetaTVD * MATMUL( invR_X1, ( U_M(:,2,1) - U_M(:,0,1) ) )
        b = c

      ELSE IF( iX1 .EQ. iX_E0(1) .AND. ExcludeOuterGhostCell(1) )THEN

        b = BetaTVD * MATMUL( invR_X1, ( U_M(:,0,1) - U_M(:,1,1) ) )
        c = b

      ELSE

        b = BetaTVD * MATMUL( invR_X1, ( U_M(:,0,1) - U_M(:,1,1) ) )
        c = BetaTVD * MATMUL( invR_X1, ( U_M(:,2,1) - U_M(:,0,1) ) )

      END IF

      dU(:,1) = MinModB( a, b, c, dX1, BetaTVB )

      IF( nDimsX > 1 )THEN

        a = MATMUL( invR_X2, U_M(:,0,3) )

        IF     ( iX2 .EQ. iX_B0(2) .AND. ExcludeInnerGhostCell(2) )THEN

          c = BetaTVD * MATMUL( invR_X2, ( U_M(:,4,1) - U_M(:,0,1) ) )
          b = c

        ELSE IF( iX2 .EQ. iX_E0(2) .AND. ExcludeOuterGhostCell(2) )THEN

          b = BetaTVD * MATMUL( invR_X2, ( U_M(:,0,1) - U_M(:,3,1) ) )
          c = b

        ELSE

          b = BetaTVD * MATMUL( invR_X2, ( U_M(:,0,1) - U_M(:,3,1) ) )
          c = BetaTVD * MATMUL( invR_X2, ( U_M(:,4,1) - U_M(:,0,1) ) )

        END IF

        dU(:,2) = MinModB( a, b, c, dX2, BetaTVB )

      END IF

      IF( nDimsX > 2 )THEN

        a = MATMUL( invR_X3, U_M(:,0,4) )

        IF     ( iX3 .EQ. iX_B0(3) .AND. ExcludeInnerGhostCell(3) )THEN

          c = BetaTVD * MATMUL( invR_X3, ( U_M(:,6,1) - U_M(:,0,1) ) )
          b = c

        ELSE IF( iX3 .EQ. iX_E0(3) .AND. ExcludeOuterGhostCell(3) )THEN

          b = BetaTVD * MATMUL( invR_X3, ( U_M(:,0,1) - U_M(:,5,1) ) )
          c = b

        ELSE

          b = BetaTVD * MATMUL( invR_X3, ( U_M(:,0,1) - U_M(:,5,1) ) )
          c = BetaTVD * MATMUL( invR_X3, ( U_M(:,6,1) - U_M(:,0,1) ) )

        END IF

        dU(:,3) = MinModB( a, b, c, dX3, BetaTVB )

      END IF

      IF( UseCharacteristicLimiting )THEN

        ! --- Transform Back from Characteristic Variables ---

        dU(:,1) = MATMUL( R_X1, dU(:,1) )

        dU(iCF_Ne,1) = dU(iCF_Ne,1) / AtomicMassUnit
        U_M(iCF_Ne,:,:) = U_M(iCF_Ne,:,:) / AtomicMassUnit !hack

        IF( nDimsX > 1 )THEN

          dU(:,2) = MATMUL( R_X2, dU(:,2) )

        END IF

        IF( nDimsX > 2 )THEN

          dU(:,3) = MATMUL( R_X3, dU(:,3) )

        END IF

      END IF

      ! --- Compute Ye Slopes ---

      Y_N(0,:) = U(:,iX1,iX2,iX3,iCF_Ne) * AtomicMassUnit &
                 / U(:,iX1,iX2,iX3,iCF_D)

      CALL MapNodalToModal_Fluid( Y_N(0,:), Y_M(0,:) )

      Y_N(1,:) = U(:,iX1-1,iX2,iX3,iCF_Ne) * AtomicMassUnit &
                   / U(:,iX1-1,iX2,iX3,iCF_D)
      Y_N(2,:) = U(:,iX1+1,iX2,iX3,iCF_Ne) * AtomicMassUnit &
                   / U(:,iX1+1,iX2,iX3,iCF_D)

      Y_M(1,1) = DOT_PRODUCT( WeightsX_q, Y_N(1,:) )
      Y_M(2,1) = DOT_PRODUCT( WeightsX_q, Y_N(2,:) )

      aY = Y_M(0,2)

      IF( .NOT. ExcludeOuterGhostCell(1) &
            .AND. .NOT. ExcludeInnerGhostCell(1) )THEN

        bY = BetaTVD * ( Y_M(0,1) - Y_M(1,1) )
        cY = BetaTVD * ( Y_M(2,1) - Y_M(0,1) )

      ELSE IF( ExcludeInnerGhostCell(1) )THEN

        cY = BetaTVD * ( Y_M(2,1) - Y_M(0,1) )
        bY = cY

      ELSE IF( ExcludeOuterGhostCell(1) )THEN

        cY = BetaTVD * ( Y_M(0,1) - Y_M(1,1) )
        bY = cY

      END IF

      dY(1) = MinModB( aY, bY, cY, dX1, BetaTVB )

      IF( nDimsX > 1 )THEN

        Y_N(3,:) = U(:,iX1,iX2-1,iX3,iCF_Ne) * AtomicMassUnit &
                     / U(:,iX1,iX2-1,iX3,iCF_D)
        Y_N(4,:) = U(:,iX1,iX2+1,iX3,iCF_Ne) * AtomicMassUnit &
                     / U(:,iX1,iX2+1,iX3,iCF_D)

        Y_M(3,1) = DOT_PRODUCT( WeightsX_q, Y_N(3,:) )
        Y_M(4,1) = DOT_PRODUCT( WeightsX_q, Y_N(4,:) )

        aY = Y_M(0,3)

        IF( .NOT. ExcludeOuterGhostCell(2) &
              .AND. .NOT. ExcludeInnerGhostCell(2) )THEN

          bY = BetaTVD * ( Y_M(0,1) - Y_M(3,1) )
          cY = BetaTVD * ( Y_M(4,1) - Y_M(0,1) )

        ELSE IF( ExcludeInnerGhostCell(2) )THEN

          cY = BetaTVD * ( Y_M(4,1) - Y_M(0,1) )
          bY = cY

        ELSE IF( ExcludeOuterGhostCell(1) )THEN

          bY = BetaTVD * ( Y_M(0,1) - Y_M(3,1) )
          cY = bY

        END IF

        dY(2) = MinModB( aY, bY, cY, dX2, BetaTVB )

      END IF

      IF( nDimsX > 2 )THEN

        Y_N(5,:) = U(:,iX1,iX2,iX3-1,iCF_Ne) * AtomicMassUnit &
                     / U(:,iX1,iX2,iX3-1,iCF_D)
        Y_N(6,:) = U(:,iX1,iX2,iX3+1,iCF_Ne) * AtomicMassUnit &
                     / U(:,iX1,iX2,iX3+1,iCF_D)

        Y_M(5,1) = DOT_PRODUCT( WeightsX_q, Y_N(5,:) )
        Y_M(6,1) = DOT_PRODUCT( WeightsX_q, Y_N(6,:) )

        aY = Y_M(0,4)

        IF( .NOT. ExcludeOuterGhostCell(3) &
              .AND. .NOT. ExcludeInnerGhostCell(3) )THEN

          bY = BetaTVD * ( Y_M(0,1) - Y_M(5,1) )
          cY = BetaTVD * ( Y_M(6,1) - Y_M(0,1) )

        ELSE IF( ExcludeInnerGhostCell(3) )THEN

          cY = BetaTVD * ( Y_M(6,1) - Y_M(0,1) )
          bY = cY

        ELSE IF( ExcludeOuterGhostCell(3) )THEN

          bY = BetaTVD * ( Y_M(0,1) - Y_M(5,1) )
          cY = bY

        END IF

        dY(3) = MinModB( aY, bY, cY, dX3, BetaTVB )

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

      LimitField = SlopeDifference > SlopeTolerance * ABS( U_M(:,0,1) )

      IF( LimitField(iCF_D) .NEQV. LimitField(iCF_Ne) )THEN

        LimitField(iCF_D)  = .TRUE.
        LimitField(iCF_Ne) = .TRUE.

      END IF

      IF( ABS( Y_M(0,2) - dY(1) ) > SlopeTolerance * ABS( Y_M(0,1) ) )THEN

        LimitField(iCF_D)  = .TRUE.
        LimitField(iCF_Ne) = .TRUE.

      END IF

      IF( nDimsX > 1 )THEN

        IF( ABS( Y_M(0,3) - dY(2) ) > SlopeTolerance * ABS( Y_M(0,1) ) )THEN

          LimitField(iCF_D)  = .TRUE.
          LimitField(iCF_Ne) = .TRUE.

        END IF

      END IF

      IF( nDimsX > 2 )THEN

        IF( ABS( Y_M(0,4) - dY(3) ) > SlopeTolerance * ABS( Y_M(0,1) ) )THEN

          LimitField(iCF_D)  = .TRUE.
          LimitField(iCF_Ne) = .TRUE.

        END IF

      END IF

      DO iCF = 1, nCF

        IF( LimitField(iCF) )THEN

          U_M(iCF,0,2:nDOFX) = Zero

          IF( iCF .NE. iCF_Ne )THEN

            U_M(iCF,0,2) = dU(iCF,1)

            IF( nDimsX > 1 ) &
              U_M(iCF,0,3) = dU(iCF,2)

            IF( nDimsX > 2 ) &
              U_M(iCF,0,4) = dU(iCF,3)

          ELSE

            ! --- Electron Density ---

            U_M(iCF,0,2) &
              = ( U_M(iCF_D,0,1) * dY(1) + Y_M(0,1) * dU(iCF_D,1) ) &
                  / AtomicMassUnit

            IF( nDimsX > 1 ) &
              U_M(iCF,0,3) &
                = ( U_M(iCF_D,0,1) * dY(2) + Y_M(0,1) * dU(iCF_D,2) ) &
                    / AtomicMassUnit

            IF( nDimsX > 2 ) &
              U_M(iCF,0,4) &
                = ( U_M(iCF_D,0,1) * dY(3) + Y_M(0,1) * dU(iCF_D,3) ) &
                    / AtomicMassUnit

          END IF

          CALL MapModalToNodal_Fluid &
                 ( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

          LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

        END IF

      END DO

    END DO
    END DO
    END DO

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, V_K, U, U_K, LimitedCell )

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_NonRelativistic_TABLE


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


END MODULE Euler_SlopeLimiterModule_NonRelativistic_TABLE
