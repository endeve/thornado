MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX, &
    nNodesX, &
    bcX
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply, &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodesX1, &
    NodesX2, &
    NodesX3, &
    WeightsX1, &
    WeightsX2, &
    WeightsX3
  USE UtilitiesModule, ONLY: &
    MinModB, &
    NodeNumberX
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    L_X1, &
    L_X2, &
    L_X3, &
    IndLX_Q
  USE PolynomialBasisModuleX_Legendre, ONLY: &
    P_X1, &
    P_X2, &
    P_X3, &
    IndPX_Q, &
    MassPX
  USE PolynomialBasisMappingModule, ONLY: &
    Pij_X
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
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
    iDF_TCI
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyInnerBC_Euler, &
    ApplyOuterBC_Euler, &
    iApplyBC_Euler_Both, &
    ApplyBoundaryConditions_Euler
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL, ONLY: &
    ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    InitializeTroubledCellIndicator_Euler, &
    FinalizeTroubledCellIndicator_Euler, &
    LimiterThreshold, &
    DetectTroubledCells_Euler
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_SlopeLimiter, &
    Timer_Euler_SL_CharDecomp, &
    Timer_Euler_SL_CopyIn, &
    Timer_Euler_SL_CopyOut, &
    Timer_Euler_SL_Integrate, &
    Timer_Euler_SL_Permute, &
    Timer_Euler_SL_LimitCells, &
    Timer_Euler_SL_ConsCorr, &
    Timer_Euler_SL_Mapping

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplySlopeLimiter_Euler_Relativistic_IDEAL

  LOGICAL      :: UseSlopeLimiter
  LOGICAL      :: UseCharacteristicLimiting
  LOGICAL      :: UseConservativeCorrection
  CHARACTER(4) :: SlopeLimiterMethod

  ! --- TVD Limiter ---

  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: I_6x6(1:6,1:6)

  ! --- Conservative Correction ---

  REAL(DP), ALLOCATABLE :: LegendreX(:,:)
  REAL(DP), ALLOCATABLE :: Kij_X    (:,:)

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
  !$OMP DECLARE TARGET &
  !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
  !$OMP     UseConservativeCorrection, &
  !$OMP     BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
  !$ACC DECLARE CREATE &
  !$ACC   ( UseSlopeLimiter, UseCharacteristicLimiting, &
  !$ACC     UseConservativeCorrection, &
  !$ACC     BetaTVD, BetaTVB, SlopeTolerance )
#endif


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

    INTEGER  :: i, j, iPol, iNX, iNX1, iNX2, iNX3, qX1, qX2, qX3
    LOGICAL  :: Verbose, UseTroubledCellIndicator
    REAL(DP) :: LimiterThresholdParameter

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

    SlopeTolerance = 1.0e-3_DP
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

    UseConservativeCorrection = .TRUE.
    IF( PRESENT( UseConservativeCorrection_Option ) ) &
      UseConservativeCorrection = UseConservativeCorrection_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeTroubledCellIndicator_Euler &
           ( UseTroubledCellIndicator_Option = UseTroubledCellIndicator, &
             LimiterThresholdParameter_Option = LimiterThresholdParameter )

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Slope Limiter (Euler, Relativistic, IDEAL)'
      WRITE(*,'(A)') &
        '    ------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A27,L1)'       ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A27,A)')         '', 'SlopeLimiterMethod: ', &
        TRIM( SlopeLimiterMethod )
      WRITE(*,*)

      IF( TRIM( SlopeLimiterMethod ) .EQ. 'TVD' )THEN

        WRITE(*,'(A6,A27,ES10.3E3)' ) '', 'BetaTVD: ' , &
          BetaTVD
        WRITE(*,'(A6,A27,ES10.3E3)' ) '', 'BetaTVB: ' , &
          BetaTVB
        WRITE(*,'(A6,A27,ES10.3E3)' ) '', 'SlopeTolerance: ' , &
          SlopeTolerance
        WRITE(*,*)

      END IF

      WRITE(*,'(A6,A27,L1)'       ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A6,A27,L1)'       ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,'(A6,A27,ES10.3E3)' ) '', 'LimiterThreshold: ' , &
        LimiterThreshold
      WRITE(*,*)
      WRITE(*,'(A6,A27,L1)'       ) '', 'UseConservativeCorrection: ' , &
        UseConservativeCorrection

    END IF

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET UPDATE TO &
    !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
    !$OMP     UseConservativeCorrection, &
    !$OMP     BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC UPDATE DEVICE &
    !$ACC   ( UseSlopeLimiter, UseCharacteristicLimiting, &
    !$ACC     UseConservativeCorrection, &
    !$ACC     BetaTVD, BetaTVB, SlopeTolerance )
#endif

    I_6x6 = Zero
    DO i = 1, 6
      I_6x6(i,i) = One
    END DO

    ALLOCATE( LegendreX(1:nDOFX,1:nDOFX) )

    DO iPol = 1, nDOFX ! Only need for iPol = 2,3,4 (FIXME)

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        iNX = NodeNumberX( iNX1, iNX2, iNX3 )

        LegendreX(iNX,iPol) &
          = P_X1  (IndPX_Q(1,iPol)) % P( NodesX1(iNX1) ) &
            * P_X2(IndPX_Q(2,iPol)) % P( NodesX2(iNX2) ) &
            * P_X3(IndPX_Q(3,iPol)) % P( NodesX3(iNX3) )

      END DO
      END DO
      END DO

    END DO

    ALLOCATE( Kij_X(1:nDOFX,1:nDOFX) )

    Kij_X = Zero
    DO j = 1, nDOFX
    DO i = 1, nDOFX

      DO qX3 = 1, nNodesX(3)
      DO qX2 = 1, nNodesX(2)
      DO qX1 = 1, nNodesX(1)

        Kij_X(i,j) &
          = Kij_X(i,j) &
              + MassPX(i) * WeightsX1(qX1) * WeightsX2(qX2) * WeightsX3(qX3) &
                  * P_X1(IndPX_Q(1,i)) % P( NodesX1(qX1) ) &
                  * P_X2(IndPX_Q(2,i)) % P( NodesX2(qX2) ) &
                  * P_X3(IndPX_Q(3,i)) % P( NodesX3(qX3) ) &
                  * L_X1(IndLX_Q(1,j)) % P( NodesX1(qX1) ) &
                  * L_X2(IndLX_Q(2,j)) % P( NodesX2(qX2) ) &
                  * L_X3(IndLX_Q(3,j)) % P( NodesX3(qX3) )

      END DO
      END DO
      END DO

    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: I_6x6, LegendreX, Kij_X, Pij_X )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  I_6x6, LegendreX, Kij_X, Pij_X )
#endif

  END SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL &
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

    LOGICAL  :: SuppressBC
    LOGICAL  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    INTEGER  :: iNX, iX1, iX2, iX3, iCF
    INTEGER  :: iApplyBC(3)

    INTEGER  :: nX_BE0, nX_BE1, nCF_BE0, nCF_BE1, nGF_BE0

    REAL(DP) :: U_N(1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: U_X(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: SqrtGm(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol(              iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCF        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: U_M(1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: a1(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b1(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c1(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: a2(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b2(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c2(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: a3(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b3(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c3(1:nCF,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))

    REAL(DP) :: dU_X1(1:nCF,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X2(1:nCF,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X3(1:nCF,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    REAL(DP) :: SlopeDifference(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    LOGICAL :: LimitedCell(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    REAL(DP) :: G_X(1:nDOFX,1:8  ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: G_K(1:8          ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X1(1:nCF,nCF   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X1(1:nCF,nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X2(1:nCF,nCF   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X2(1:nCF,nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X3(1:nCF,nCF   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X3(1:nCF,nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nX_BE0  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nX_BE1  = PRODUCT( iX_E1 - iX_B1 + 1 )
    nCF_BE0 = nCF * nX_BE0
    nCF_BE1 = nCF * nX_BE1
    nGF_BE0 = 8   * nX_BE0

    iApplyBC = iApplyBC_Euler_Both
    IF( PRESENT( iApplyBC_Option ) ) &
       iApplyBC = iApplyBC_Option

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

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

    CALL TimersStart_Euler( Timer_Euler_SL_CopyIn )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3 ) &
    !$OMP MAP( alloc: U_N, U_X, SqrtGm, Vol, U_K, U_M, &
    !$OMP             a1, b1, c1, a2, b2, c2, a3, b3, c3, &
    !$OMP             dU_X1, dU_X2, dU_X3, SlopeDifference, LimitedCell, &
    !$OMP             G_X, G_K, R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3 ) &
    !$ACC CREATE(     U_N, U_X, SqrtGm, Vol, U_K, U_M, &
    !$ACC             a1, b1, c1, a2, b2, c2, a3, b3, c3, &
    !$ACC             dU_X1, dU_X2, dU_X3, SlopeDifference, LimitedCell, &
    !$ACC             G_X, G_K, R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_N(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, SqrtGm, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, U_X, U_N, SqrtGm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_X(iNX,iCF,iX1,iX2,iX3) &
        = U_N(iNX,iCF,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_Permute )

    CALL TimersStart_Euler( Timer_Euler_SL_Integrate )

    ! --- Compute volumes of compute cells ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nX_BE0 , One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    ! --- Compute cell integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_BE0, One, U_X   , nDOFX, &
             WeightsX_q, 1, Zero, U_K, 1 )

    CALL TimersStop_Euler( Timer_Euler_SL_Integrate )

    ! --- Form cell averages ---

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K, Vol )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      U_K(iCF,iX1,iX2,iX3) = U_K(iCF,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_Permute )

    ! --- Map Nodal to Modal ---

    CALL TimersStart_Euler( Timer_Euler_SL_Mapping )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_BE1, nDOFX, One, Kij_X, nDOFX, &
             U_N, nDOFX, Zero, U_M, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_SL_Mapping )

    CALL ComputeMinModArguments &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_M, &
             ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
             a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    IF( UseCharacteristicLimiting )THEN

      CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, G_X, G )
#elif defined(THORNADO_OMP) && !defined(THORNADO_EULER_NOGPU)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        G_X(iNX,1,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)
        G_X(iNX,2,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)
        G_X(iNX,3,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)
        G_X(iNX,4,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm  )**2
        G_X(iNX,5,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Alpha   )
        G_X(iNX,6,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_1  )
        G_X(iNX,7,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_2  )
        G_X(iNX,8,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_3  )

      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_SL_Permute )

      CALL TimersStart_Euler( Timer_Euler_SL_Integrate )

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nGF_BE0, One, G_X, nDOFX, &
               WeightsX_q, 1, Zero, G_K, 1 )

      CALL TimersStop_Euler( Timer_Euler_SL_Integrate )

      CALL ComputeEigenvectorMatrices &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G_K, U_M, &
               R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )

      CALL MultiplyWithInverseEigenvectorMatrices &
             ( iX_B0, iX_E0, invR_X1, invR_X2, invR_X3, &
               a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    END IF

    CALL TimersStart_Euler( Timer_Euler_SL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X1, a1, b1, c1, dX1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      dU_X1(iCF,iX1,iX2,iX3) &
        = MinModB( a1(iCF,iX1,iX2,iX3), &
                   b1(iCF,iX1,iX2,iX3), &
                   c1(iCF,iX1,iX2,iX3), &
                   dX1(iX1), BetaTVB )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X2, a2, b2, c2, dX2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        dU_X2(iCF,iX1,iX2,iX3) &
          = MinModB( a2(iCF,iX1,iX2,iX3), &
                     b2(iCF,iX1,iX2,iX3), &
                     c2(iCF,iX1,iX2,iX3), &
                     dX2(iX2), BetaTVB )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X3, a3, b3, c3, dX3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        dU_X3(iCF,iX1,iX2,iX3) &
          = MinModB( a3(iCF,iX1,iX2,iX3), &
                     b3(iCF,iX1,iX2,iX3), &
                     c3(iCF,iX1,iX2,iX3), &
                     dX3(iX3), BetaTVB )

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_SL_LimitCells )

    IF( UseCharacteristicLimiting ) &
      CALL MultiplyWithEigenvectorMatrices &
             ( iX_B0, iX_E0, R_X1, R_X2, R_X3, dU_X1, dU_X2, dU_X3 )

    CALL TimersStart_Euler( Timer_Euler_SL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X1, U_M, SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      SlopeDifference(iCF,iX1,iX2,iX3) &
        = ABS( U_M(2,iCF,iX1,iX2,iX3) - dU_X1(iCF,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X2, U_M, SlopeDifference )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(3,iCF,iX1,iX2,iX3) - dU_X2(iCF,iX1,iX2,iX3) ) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X3, U_M, SlopeDifference )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(4,iCF,iX1,iX2,iX3) - dU_X3(iCF,iX1,iX2,iX3) ) )

      END DO
      END DO
      END DO
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, LimitedCell, U_M, dU_X1, dU_X2, dU_X3, D, &
    !$ACC          SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      LimitedCell(iCF,iX1,iX2,iX3) = .FALSE.

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      IF( SlopeDifference(iCF,iX1,iX2,iX3) &
            .GT. SlopeTolerance * ABS( U_M(1,iCF,iX1,iX2,iX3) ) )THEN

        DO iNX = 2, nDOFX

          U_M(iNX,iCF,iX1,iX2,iX3) = Zero

        END DO

        U_M(2,iCF,iX1,iX2,iX3) = dU_X1(iCF,iX1,iX2,iX3)

        IF( nDimsX .GT. 1 ) U_M(3,iCF,iX1,iX2,iX3) = dU_X2(iCF,iX1,iX2,iX3)

        IF( nDimsX .GT. 2 ) U_M(4,iCF,iX1,iX2,iX3) = dU_X3(iCF,iX1,iX2,iX3)

        LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_LimitCells )

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL TimersStart_Euler( Timer_Euler_SL_Mapping )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_BE1, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_SL_Mapping )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCF) = U_N(iNX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_Permute )

    CALL TimersStart_Euler( Timer_Euler_SL_CopyOut )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$OMP               U_N, U_X, SqrtGm, Vol, U_K, U_M, &
    !$OMP               a1, b1, c1, a2, b2, c2, a3, b3, c3, &
    !$OMP               dU_X1, dU_X2, dU_X3, SlopeDifference, LimitedCell, &
    !$OMP               G_X, G_K, R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$ACC               U_N, U_X, SqrtGm, Vol, U_K, U_M, &
    !$ACC               a1, b1, c1, a2, b2, c2, a3, b3, c3, &
    !$ACC               dU_X1, dU_X2, dU_X3, SlopeDifference, LimitedCell, &
    !$ACC               G_X, G_K, R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyOut )

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: I_6x6, LegendreX, Kij_X, Pij_X )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(       I_6x6, LegendreX, Kij_X, Pij_X )
#endif

    DEALLOCATE( Kij_X )
    DEALLOCATE( LegendreX )

    CALL FinalizeTroubledCellIndicator_Euler

  END SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplyConservativeCorrection &
    ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      SqrtGm(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Vol   (        iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K   (1:nCF  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    LOGICAL, INTENT(in)     :: &
      LimitedCell(1:nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: &
      U_M(1:nDOFX,1:nCF,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, iDimX
    REAL(DP) :: Correction, Term

    IF( .NOT. UseConservativeCorrection ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SL_ConsCorr )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell, U_M )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell, U_M )
#endif

    ! --- Applies a correction to the 0-th order ---
    ! --- mode to maintain the cell average.     ---

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, LimitedCell, WeightsX_q, LegendreX, SqrtGm, &
    !$ACC          U_M, Vol, U_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      IF( LimitedCell(iCF,iX1,iX2,iX3) )THEN

        Correction = Zero

        DO iDimX = 1, nDimsX

          Term = Zero

          DO iNX = 1, nDOFX

            Term &
              = Term &
                  + WeightsX_q(iNX) &
                      * LegendreX(iNX,iDimX+1) * SqrtGm(iNX,iX1,iX2,iX3)

          END DO

          Correction &
            = Correction &
                + U_M(iDimX+1,iCF,iX1,iX2,iX3) * Term / Vol(iX1,iX2,iX3)

        END DO

        U_M(1,iCF,iX1,iX2,iX3) = U_K(iCF,iX1,iX2,iX3) - Correction

      END IF

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_M ) &
    !$OMP MAP( release: iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_M ) &
    !$ACC DELETE(       iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_ConsCorr )

  END SUBROUTINE ApplyConservativeCorrection


  SUBROUTINE ComputeMinModArguments &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U_M, &
      ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
      a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    LOGICAL,  INTENT(in)  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    REAL(DP), INTENT(in)  :: U_M(nDOFX,nCF,iX_B1(1):iX_E1(1), &
                                           iX_B1(2):iX_E1(2), &
                                           iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(out) :: a1(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b1(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c1(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: a2(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b2(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c2(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: a3(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b3(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c3(nCF,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    INTEGER :: iX1, iX2, iX3, iCF

    CALL TimersStart_Euler( Timer_Euler_SL_LimitCells )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U_M, a1, b1, c1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      a1(iCF,iX1,iX2,iX3) = U_M(2,iCF,iX1,iX2,iX3)

      b1(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1  ,iX2,iX3) &
                                       - U_M(1,iCF,iX1-1,iX2,iX3) )
      c1(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1+1,iX2,iX3) &
                                       - U_M(1,iCF,iX1  ,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, U_M, a2, b2, c2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        a2(iCF,iX1,iX2,iX3) = U_M(3,iCF,iX1,iX2,iX3)

        b2(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1,iX2  ,iX3) &
                                         - U_M(1,iCF,iX1,iX2-1,iX3) )
        c2(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1,iX2+1,iX3) &
                                         - U_M(1,iCF,iX1,iX2  ,iX3) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, U_M, a3, b3, c3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        a3(iCF,iX1,iX2,iX3) = U_M(4,iCF,iX1,iX2,iX3)

        b3(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1,iX2,iX3  ) &
                                         - U_M(1,iCF,iX1,iX2,iX3-1) )
        c3(iCF,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCF,iX1,iX2,iX3+1) &
                                         - U_M(1,iCF,iX1,iX2,iX3  ) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( ExcludeInnerGhostCell(1) )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b1, c1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iCF = 1, nCF

        b1(iCF,iX_B0(1),iX2,iX3) = c1(iCF,iX_B0(1),iX2,iX3)

      END DO
      END DO
      END DO

    END IF

    IF( ExcludeOuterGhostCell(1) )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b1, c1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iCF = 1, nCF

        c1(iCF,iX_E0(1),iX2,iX3) = b1(iCF,iX_E0(1),iX2,iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 1 .AND. ExcludeInnerGhostCell(2) )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b2, c2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        b2(iCF,iX1,iX_B0(2),iX3) = c2(iCF,iX1,iX_B0(2),iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 1 .AND. ExcludeOuterGhostCell(2) )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b2, c2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        c2(iCF,iX1,iX_E0(2),iX3) = b2(iCF,iX1,iX_E0(2),iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 .AND. ExcludeInnerGhostCell(3) )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b3, c3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        b3(iCF,iX1,iX2,iX_B0(3)) = c3(iCF,iX1,iX2,iX_B0(3))

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 .AND. ExcludeOuterGhostCell(3) )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, b3, c3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        c3(iCF,iX1,iX2,iX_E0(3)) = b3(iCF,iX1,iX2,iX_E0(3))

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_SL_LimitCells )

  END SUBROUTINE ComputeMinModArguments


  SUBROUTINE ComputeEigenvectorMatrices &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G_K, U_M, &
      R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: G_K(8        ,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in)  :: U_M(nDOFX,nCF,iX_B1(1):iX_E1(1), &
                                           iX_B1(2):iX_E1(2), &
                                           iX_B1(3):iX_E1(3))

    REAL(DP), INTENT(out) :: R_X1(   nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: invR_X1(nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: R_X2(   nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: invR_X2(nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: R_X3(   nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: invR_X3(nCF,nCF,iX_B0(1):iX_E0(1), &
                                             iX_B0(2):iX_E0(2), &
                                             iX_B0(3):iX_E0(3))

    INTEGER :: iX1, iX2, iX3, iGF, iCF, jCF

    REAL(DP) :: R(nCF,nCF), invR(nCF,nCF), UK(nCF), GK(8)

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, U_M, R_X1, invR_X1 ) &
    !$ACC PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( R, invR, UK, GK )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iGF = 1, 8

        GK(iGF) = G_K(iGF,iX1,iX2,iX3)

      END DO

      DO iCF = 1, nCF

        UK(iCF) = U_M(1,iCF,iX1,iX2,iX3)

      END DO

      CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
             ( 1, GK, UK, R, invR )

      DO iCF = 1, nCF
      DO jCF = 1, nCF

        R_X1   (jCF,iCF,iX1,iX2,iX3) = R   (iCF,jCF)
        invR_X1(jCF,iCF,iX1,iX2,iX3) = invR(iCF,jCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K, U_M, R_X2, invR_X2 ) &
      !$ACC PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( R, invR, UK, GK )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        DO iGF = 1, 8

          GK(iGF) = G_K(iGF,iX1,iX2,iX3)

        END DO

        DO iCF = 1, nCF

          UK(iCF) = U_M(1,iCF,iX1,iX2,iX3)

        END DO

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 2, GK, UK, R, invR )

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          R_X2   (jCF,iCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X2(jCF,iCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K, U_M, R_X3, invR_X3 ) &
      !$ACC PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( R, invR, UK, GK )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        DO iGF = 1, 8

          GK(iGF) = G_K(iGF,iX1,iX2,iX3)

        END DO

        DO iCF = 1, nCF

          UK(iCF) = U_M(1,iCF,iX1,iX2,iX3)

        END DO

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 3, GK, UK, R, invR )

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          R_X3   (jCF,iCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X3(jCF,iCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

  END SUBROUTINE ComputeEigenvectorMatrices


  SUBROUTINE MultiplyWithInverseEigenvectorMatrices &
    ( iX_B0, iX_E0, invR_X1, invR_X2, invR_X3, &
      a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3)

    REAL(DP), INTENT(in) :: invR_X1(nCF,nCF,iX_B0(1):iX_E0(1), &
                                            iX_B0(2):iX_E0(2), &
                                            iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in) :: invR_X2(nCF,nCF,iX_B0(1):iX_E0(1), &
                                            iX_B0(2):iX_E0(2), &
                                            iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in) :: invR_X3(nCF,nCF,iX_B0(1):iX_E0(1), &
                                            iX_B0(2):iX_E0(2), &
                                            iX_B0(3):iX_E0(3))

    REAL(DP), INTENT(inout) :: a1(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: b1(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: c1(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: a2(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: b2(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: c2(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: a3(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: b3(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: c3(nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: a1t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: b1t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: c1t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: a2t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: b2t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: c2t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: a3t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: b3t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))
    REAL(DP) :: c3t(nCF,iX_B0(1):iX_E0(1), &
                        iX_B0(2):iX_E0(2), &
                        iX_B0(3):iX_E0(3))

    INTEGER :: iX1, iX2, iX3, iCF, jCF

    REAL(DP) :: Sum_a, Sum_b, Sum_c

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: a1t, b1t, c1t, a2t, b2t, c2t, a3t, b3t, c3t )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC CREATE(     a1t, b1t, c1t, a2t, b2t, c2t, a3t, b3t, c3t )
#endif

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, a1t, b1t, c1t, a2t, b2t, c2t, &
    !$ACC          a3t, b3t, c3t, a1, b1, c1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      a1t(iCF,iX1,iX2,iX3) = a1(iCF,iX1,iX2,iX3)
      b1t(iCF,iX1,iX2,iX3) = b1(iCF,iX1,iX2,iX3)
      c1t(iCF,iX1,iX2,iX3) = c1(iCF,iX1,iX2,iX3)

      a2t(iCF,iX1,iX2,iX3) = a2(iCF,iX1,iX2,iX3)
      b2t(iCF,iX1,iX2,iX3) = b2(iCF,iX1,iX2,iX3)
      c2t(iCF,iX1,iX2,iX3) = c2(iCF,iX1,iX2,iX3)

      a3t(iCF,iX1,iX2,iX3) = a2(iCF,iX1,iX2,iX3)
      b3t(iCF,iX1,iX2,iX3) = b2(iCF,iX1,iX2,iX3)
      c3t(iCF,iX1,iX2,iX3) = c2(iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Sum_a, Sum_b, Sum_c ) &
    !$ACC PRESENT( iX_B0, iX_E0, a1t, b1t, c1t, a1, b1, c1, invR_X1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      Sum_a = Zero
      Sum_b = Zero
      Sum_c = Zero

      DO jCF = 1, nCF

        Sum_a = Sum_a + invR_X1(jCF,iCF,iX1,iX2,iX3) * a1t(jCF,iX1,iX2,iX3)
        Sum_b = Sum_b + invR_X1(jCF,iCF,iX1,iX2,iX3) * b1t(jCF,iX1,iX2,iX3)
        Sum_c = Sum_c + invR_X1(jCF,iCF,iX1,iX2,iX3) * c1t(jCF,iX1,iX2,iX3)

      END DO

      a1(iCF,iX1,iX2,iX3) = Sum_a
      b1(iCF,iX1,iX2,iX3) = Sum_b
      c1(iCF,iX1,iX2,iX3) = Sum_c

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( Sum_a, Sum_b, Sum_c ) &
      !$ACC PRESENT( iX_B0, iX_E0, a2t, b2t, c2t, a2, b2, c2, invR_X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4) &
      !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        Sum_a = Zero
        Sum_b = Zero
        Sum_c = Zero

        DO jCF = 1, nCF

          Sum_a = Sum_a + invR_X2(jCF,iCF,iX1,iX2,iX3) * a2t(jCF,iX1,iX2,iX3)
          Sum_b = Sum_b + invR_X2(jCF,iCF,iX1,iX2,iX3) * b2t(jCF,iX1,iX2,iX3)
          Sum_c = Sum_c + invR_X2(jCF,iCF,iX1,iX2,iX3) * c2t(jCF,iX1,iX2,iX3)

        END DO

        a2(iCF,iX1,iX2,iX3) = Sum_a
        b2(iCF,iX1,iX2,iX3) = Sum_b
        c2(iCF,iX1,iX2,iX3) = Sum_c

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( Sum_a, Sum_b, Sum_c ) &
      !$ACC PRESENT( iX_B0, iX_E0, a3t, b3t, c3t, a3, b3, c3, invR_X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4) &
      !$OMP PRIVATE( Sum_a, Sum_b, Sum_c )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        Sum_a = Zero
        Sum_b = Zero
        Sum_c = Zero

        DO jCF = 1, nCF

          Sum_a = Sum_a + invR_X3(jCF,iCF,iX1,iX2,iX3) * a3t(jCF,iX1,iX2,iX3)
          Sum_b = Sum_b + invR_X3(jCF,iCF,iX1,iX2,iX3) * b3t(jCF,iX1,iX2,iX3)
          Sum_c = Sum_c + invR_X3(jCF,iCF,iX1,iX2,iX3) * c3t(jCF,iX1,iX2,iX3)

        END DO

        a3(iCF,iX1,iX2,iX3) = Sum_a
        b3(iCF,iX1,iX2,iX3) = Sum_b
        c3(iCF,iX1,iX2,iX3) = Sum_c

      END DO
      END DO
      END DO
      END DO

    END IF

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: a1t, b1t, c1t, a2t, b2t, c2t, a3t, b3t, c3t )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC EXIT DATA &
    !$ACC DELETE(       a1t, b1t, c1t, a2t, b2t, c2t, a3t, b3t, c3t )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

  END SUBROUTINE MultiplyWithInverseEigenvectorMatrices


  SUBROUTINE MultiplyWithEigenvectorMatrices &
    ( iX_B0, iX_E0, R_X1, R_X2, R_X3, dU_X1, dU_X2, dU_X3 )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3)

    REAL(DP), INTENT(in) :: R_X1(nCF,nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in) :: R_X2(nCF,nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in) :: R_X3(nCF,nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))

    REAL(DP), INTENT(inout) :: dU_X1(nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: dU_X2(nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: dU_X3(nCF,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))

    REAL(DP) :: dU_X1t(nCF,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X2t(nCF,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X3t(nCF,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))

    INTEGER :: iX1, iX2, iX3, iCF, jCF

    REAL(DP) :: Sum_X1, Sum_X2, Sum_X3

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: dU_X1t, dU_X2t, dU_X3t )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC CREATE(     dU_X1t, dU_X2t, dU_X3t )
#endif

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X1, dU_X2, dU_X3, dU_X1t, dU_X2t, dU_X3t )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      dU_X1t(iCF,iX1,iX2,iX3) = dU_X1(iCF,iX1,iX2,iX3)
      dU_X2t(iCF,iX1,iX2,iX3) = dU_X2(iCF,iX1,iX2,iX3)
      dU_X3t(iCF,iX1,iX2,iX3) = dU_X3(iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Sum_X1 )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Sum_X1 ) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X1, dU_X1t, R_X1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Sum_X1 )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      Sum_X1 = Zero

      DO jCF = 1, nCF

        Sum_X1 &
          = Sum_X1 + R_X1(jCF,iCF,iX1,iX2,iX3) * dU_X1t(jCF,iX1,iX2,iX3)

      END DO

      dU_X1(iCF,iX1,iX2,iX3) = Sum_X1

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( Sum_X2 )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( Sum_X2 ) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X2, dU_X2t, R_X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4) &
      !$OMP PRIVATE( Sum_X2 )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        Sum_X2 = Zero

        DO jCF = 1, nCF

          Sum_X2 &
            = Sum_X2 + R_X2(jCF,iCF,iX1,iX2,iX3) * dU_X2t(jCF,iX1,iX2,iX3)

        END DO

        dU_X2(iCF,iX1,iX2,iX3) = Sum_X2

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( Sum_X3 )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( Sum_X3 ) &
      !$ACC PRESENT( iX_B0, iX_E0, dU_X3, dU_X3t, R_X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4) &
      !$OMP PRIVATE( Sum_X3 )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCF = 1, nCF

        Sum_X3 = Zero

        DO jCF = 1, nCF

          Sum_X3 &
            = Sum_X3 + R_X3(jCF,iCF,iX1,iX2,iX3) * dU_X3t(jCF,iX1,iX2,iX3)

        END DO

        dU_X3(iCF,iX1,iX2,iX3) = Sum_X3

      END DO
      END DO
      END DO
      END DO

    END IF

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dU_X1t, dU_X2t, dU_X3t )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC EXIT DATA &
    !$ACC DELETE(       dU_X1t, dU_X2t, dU_X3t )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

  END SUBROUTINE MultiplyWithEigenvectorMatrices


END MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL
