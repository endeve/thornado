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

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
  !$OMP     UseConservativeCorrection, &
  !$OMP     BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
    !$OMP     UseConservativeCorrection, &
    !$OMP     BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: I_6x6, LegendreX, Kij_X, Pij_X )
#elif defined(THORNADO_OACC)
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

    IF( .NOT. UseCharacteristicLimiting )THEN

      CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_Componentwise &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
               SuppressBC_Option, iApplyBC_Option )

    ELSE

      CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_Characteristic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
                 SuppressBC_Option, iApplyBC_Option )

    END IF

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_Componentwise &
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
    REAL(DP) :: SqrtGamma

    INTEGER  :: nX_BE0, nX_BE1, nCF_BE0, nCF_BE1

    REAL(DP) :: dU(1:nCF,1:nDimsX    ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: SqrtGm(1:nDOFX       ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol(                  iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_X(1:nDOFX,1:nCF    ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCF            ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: U_N(1:nDOFX,1:nCF    ,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))
    REAL(DP) :: U_M(1:nDOFX,1:nCF    ,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))

    REAL(DP) :: a(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: b(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: c(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: SlopeDifference(1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    LOGICAL  :: LimitedCell(1:nCF    ,iX_B0(1):iX_E0(1), &
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$OMP             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$OMP MAP( alloc: dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$OMP             a, b, c, SlopeDifference, LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$ACC             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC CREATE(     dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$ACC             a, b, c, SlopeDifference, LimitedCell )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, SqrtGm, U_X, G, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGamma = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      SqrtGm(iNX,iX1,iX2,iX3) = SqrtGamma

      U_X(iNX,iCF_D ,iX1,iX2,iX3) = U_N(iNX,iCF_D ,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S1,iX1,iX2,iX3) = U_N(iNX,iCF_S1,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S2,iX1,iX2,iX3) = U_N(iNX,iCF_S2,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S3,iX1,iX2,iX3) = U_N(iNX,iCF_S3,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_E ,iX1,iX2,iX3) = U_N(iNX,iCF_E ,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_Ne,iX1,iX2,iX3) = U_N(iNX,iCF_Ne,iX1,iX2,iX3) * SqrtGamma

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K, Vol )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
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

    CALL TimersStart_Euler( Timer_Euler_SL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_M, dX1, dX2, dX3, LimitedCell, &
    !$ACC          ExcludeInnerGhostCell, ExcludeOuterGhostCell, dU, &
    !$ACC          a, b, c, SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      LimitedCell(iCF,iX1,iX2,iX3) = .FALSE.

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      a(1,iCF,iX1,iX2,iX3) = U_M(2,iCF,iX1,iX2,iX3)

      IF     ( iX1 .EQ. iX_B0(1) .AND. ExcludeInnerGhostCell(1) )THEN

        c(1,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1+1,iX2,iX3) &
                                         - U_M(1,iCF,iX1  ,iX2,iX3) )

        b(1,iCF,iX1,iX2,iX3) = c(1,iCF,iX1,iX2,iX3)

      ELSE IF( iX1 .EQ. iX_E0(1) .AND. ExcludeOuterGhostCell(1) )THEN

        b(1,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1  ,iX2,iX3) &
                                         - U_M(1,iCF,iX1-1,iX2,iX3) )

        c(1,iCF,iX1,iX2,iX3) = b(1,iCF,iX1,iX2,iX3)


      ELSE

        b(1,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1  ,iX2,iX3) &
                                         - U_M(1,iCF,iX1-1,iX2,iX3) )
        c(1,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1+1,iX2,iX3) &
                                         - U_M(1,iCF,iX1  ,iX2,iX3) )

      END IF

      dU(iCF,1,iX1,iX2,iX3) &
        = MinModB( a(1,iCF,iX1,iX2,iX3), &
                   b(1,iCF,iX1,iX2,iX3), &
                   c(1,iCF,iX1,iX2,iX3), &
                   dX1(iX1), BetaTVB )

      IF( nDimsX .GT. 1 )THEN

        a(2,iCF,iX1,iX2,iX3) = U_M(3,iCF,iX1,iX2,iX3)

        IF     ( iX2 .EQ. iX_B0(2) .AND. ExcludeInnerGhostCell(2) )THEN

          c(2,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2+1,iX3) &
                                           - U_M(1,iCF,iX1,iX2  ,iX3) )

          b(2,iCF,iX1,iX2,iX3) = c(2,iCF,iX1,iX2,iX3)

        ELSE IF( iX2 .EQ. iX_E0(2) .AND. ExcludeOuterGhostCell(2) )THEN

          b(2,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2  ,iX3) &
                                           - U_M(1,iCF,iX1,iX2-1,iX3) )

          c(2,iCF,iX1,iX2,iX3) = b(2,iCF,iX1,iX2,iX3)


        ELSE

          b(2,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2  ,iX3) &
                                           - U_M(1,iCF,iX1,iX2-1,iX3) )

          c(2,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2+1,iX3) &
                                           - U_M(1,iCF,iX1,iX2  ,iX3) )

        END IF

        dU(iCF,2,iX1,iX2,iX3) &
          = MinModB( a(2,iCF,iX1,iX2,iX3), &
                     b(2,iCF,iX1,iX2,iX3), &
                     c(2,iCF,iX1,iX2,iX3), &
                     dX2(iX2), BetaTVB )

      END IF

      IF( nDimsX .GT. 2 )THEN

        a(3,iCF,iX1,iX2,iX3) = U_M(4,iCF,iX1,iX2,iX3)

        IF     ( iX3 .EQ. iX_B0(3) .AND. ExcludeInnerGhostCell(3) )THEN

          c(3,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3+1) &
                                           - U_M(1,iCF,iX1,iX2,iX3  ) )

          b(3,iCF,iX1,iX2,iX3) = c(3,iCF,iX1,iX2,iX3)

        ELSE IF( iX3 .EQ. iX_E0(3) .AND. ExcludeOuterGhostCell(3) )THEN

          b(3,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3  ) &
                                           - U_M(1,iCF,iX1,iX2,iX3-1) )

          c(3,iCF,iX1,iX2,iX3) = b(3,iCF,iX1,iX2,iX3)


        ELSE

          b(3,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3  ) &
                                           - U_M(1,iCF,iX1,iX2,iX3-1) )

          c(3,iCF,iX1,iX2,iX3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3+1) &
                                           - U_M(1,iCF,iX1,iX2,iX3  ) )

        END IF

        dU(iCF,3,iX1,iX2,iX3) &
          = MinModB( a(3,iCF,iX1,iX2,iX3), &
                     b(3,iCF,iX1,iX2,iX3), &
                     c(3,iCF,iX1,iX2,iX3), &
                     dX3(iX3), BetaTVB )

      END IF

      ! --- Compare Limited Slopes to Original Slopes ---

      SlopeDifference(iCF,iX1,iX2,iX3) &
        = ABS( U_M(2,iCF,iX1,iX2,iX3) - dU(iCF,1,iX1,iX2,iX3) )

      IF( nDimsX .GT. 1 ) &
        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(3,iCF,iX1,iX2,iX3) - dU(iCF,2,iX1,iX2,iX3) ) )

      IF( nDimsX .GT. 2 ) &
        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(4,iCF,iX1,iX2,iX3) - dU(iCF,3,iX1,iX2,iX3) ) )

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      IF( SlopeDifference(iCF,iX1,iX2,iX3) &
            .GT. SlopeTolerance * ABS( U_M(1,iCF,iX1,iX2,iX3) ) )THEN

        DO iNX = 2, nDOFX

          U_M(iNX,iCF,iX1,iX2,iX3) = Zero

        END DO

        U_M(2,iCF,iX1,iX2,iX3) = dU(iCF,1,iX1,iX2,iX3)

        IF( nDimsX .GT. 1 ) U_M(3,iCF,iX1,iX2,iX3) = dU(iCF,2,iX1,iX2,iX3)

        IF( nDimsX .GT. 2 ) U_M(4,iCF,iX1,iX2,iX3) = dU(iCF,3,iX1,iX2,iX3)

        LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_LimitCells )

    ! --- Apply Conservative Correction ---

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL TimersStart_Euler( Timer_Euler_SL_Mapping )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_BE1, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_SL_Mapping )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$OMP               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$OMP               dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$OMP               a, b, c, SlopeDifference, LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$ACC               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$ACC               dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$ACC               a, b, c, SlopeDifference, LimitedCell )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyOut )

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_Componentwise


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_Characteristic &
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
    INTEGER  :: iNX, iX1, iX2, iX3, iCF, jCF, iGF
    INTEGER  :: iApplyBC(3)
    REAL(DP) :: SqrtGamma

    INTEGER  :: nX_BE0, nX_BE1, nCF_BE0, nCF_BE1, nGF_BE0

    REAL(DP) :: dU(1:nCF,1:nDimsX    ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: SqrtGm(1:nDOFX       ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol(                  iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_X(1:nDOFX,1:nCF    ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCF            ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: U_N(1:nDOFX,1:nCF    ,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))
    REAL(DP) :: U_M(1:nDOFX,1:nCF    ,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))

    REAL(DP) :: a(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: b(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: c(1:nDimsX,1:nCF     ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: SlopeDifference(1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    LOGICAL  :: LimitedCell(1:nCF    ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    ! --- For characteristic limiting ---

    REAL(DP) :: R_X1   (1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X1(1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X2   (1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X2(1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X3   (1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X3(1:nCF,1:nCF  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: G_X(1:nDOFX,1:8      ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: G_K(1:8              ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_C(1:nCF,1:nDimsX  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: R(1:nCF,1:nCF), invR(1:nCF,1:nCF)
    REAL(DP) :: GK(8), UK(nCF)

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$OMP             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$OMP MAP( alloc: dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$OMP             a, b, c, SlopeDifference, LimitedCell, &
    !$OMP             R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$OMP             G_X, G_K, dU_C )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$ACC             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC CREATE(     dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$ACC             a, b, c, SlopeDifference, LimitedCell, &
    !$ACC             R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$ACC             G_X, G_K, dU_C )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, SqrtGm, U_X, G_X, G, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGamma = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      SqrtGm(iNX,iX1,iX2,iX3) = SqrtGamma

      U_X(iNX,iCF_D ,iX1,iX2,iX3) = U_N(iNX,iCF_D ,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S1,iX1,iX2,iX3) = U_N(iNX,iCF_S1,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S2,iX1,iX2,iX3) = U_N(iNX,iCF_S2,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_S3,iX1,iX2,iX3) = U_N(iNX,iCF_S3,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_E ,iX1,iX2,iX3) = U_N(iNX,iCF_E ,iX1,iX2,iX3) * SqrtGamma
      U_X(iNX,iCF_Ne,iX1,iX2,iX3) = U_N(iNX,iCF_Ne,iX1,iX2,iX3) * SqrtGamma

      G_X(iNX,1,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)
      G_X(iNX,2,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)
      G_X(iNX,3,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)
      G_X(iNX,4,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm  )
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

    ! --- Compute volumes of compute cells ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nX_BE0 , One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    ! --- Compute cell integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_BE0, One, U_X   , nDOFX, &
             WeightsX_q, 1, Zero, U_K, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nGF_BE0, One, G_X   , nDOFX, &
             WeightsX_q, 1, Zero, G_K, 1 )

    CALL TimersStop_Euler( Timer_Euler_SL_Integrate )

    ! --- Form cell averages ---

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K, Vol )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
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

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, D, G_K, U_M, &
    !$ACC          R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 ) &
    !$ACC PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( R, invR, UK, GK )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      DO iGF = 1, 8

        GK(iGF) = G_K(iGF,iX1,iX2,iX3)

      END DO

      DO iCF = 1, nCF

        UK(iCF) = U_M(1,iCF,iX1,iX2,iX3)

      END DO

      CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
             ( 1, GK, UK, R, invR )

      DO jCF = 1, nCF
      DO iCF = 1, nCF

        R_X1   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
        invR_X1(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

      END DO
      END DO

      IF( nDimsX .GT. 1 )THEN

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 2, GK, UK, R, invR )

        DO jCF = 1, nCF
        DO iCF = 1, nCF

          R_X2   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X2(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END IF

      IF( nDimsX .GT. 2 )THEN

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 3, GK, UK, R, invR )

        DO jCF = 1, nCF
        DO iCF = 1, nCF

          R_X3   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X3(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

    CALL TimersStart_Euler( Timer_Euler_SL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_M, dX1, dX2, dX3, &
    !$ACC          ExcludeInnerGhostCell, ExcludeOuterGhostCell, dU, &
    !$ACC          a, b, c, SlopeDifference, &
    !$ACC          R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, dU_C )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      DO iCF = 1, nCF

        a(1,iCF,iX1,iX2,iX3) = Zero
        b(1,iCF,iX1,iX2,iX3) = Zero
        c(1,iCF,iX1,iX2,iX3) = Zero

        DO jCF = 1, nCF

          a(1,iCF,iX1,iX2,iX3) &
            = a(1,iCF,iX1,iX2,iX3) &
                + invR_X1(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

          b(1,iCF,iX1,iX2,iX3) &
            = b(1,iCF,iX1,iX2,iX3) &
                + BetaTVD * invR_X1(iCF,jCF,iX1,iX2,iX3) &
                    * ( U_M(1,jCF,iX1  ,iX2,iX3) &
                      - U_M(1,jCF,iX1-1,iX2,iX3) )

          c(1,iCF,iX1,iX2,iX3) &
            = c(1,iCF,iX1,iX2,iX3) &
                + BetaTVD * invR_X1(iCF,jCF,iX1,iX2,iX3) &
                    * ( U_M(1,jCF,iX1+1,iX2,iX3) &
                      - U_M(1,jCF,iX1  ,iX2,iX3) )

        END DO

      END DO

      IF( iX1 .EQ. iX_B0(1) .AND. ExcludeInnerGhostCell(1) )THEN

        DO iCF = 1, nCF

          b(1,iCF,iX1,iX2,iX3) = c(1,iCF,iX1,iX2,iX3)

        END DO

      END IF

      IF( iX1 .EQ. iX_E0(1) .AND. ExcludeOuterGhostCell(1) )THEN

        DO iCF = 1, nCF

          c(1,iCF,iX1,iX2,iX3) = b(1,iCF,iX1,iX2,iX3)

        END DO

      END IF

      DO iCF = 1, nCF

        dU_C(iCF,1,iX1,iX2,iX3) &
          = MinModB( a(1,iCF,iX1,iX2,iX3), &
                     b(1,iCF,iX1,iX2,iX3), &
                     c(1,iCF,iX1,iX2,iX3), &
                     dX1(iX1), BetaTVB )

        dU(iCF,1,iX1,iX2,iX3) = Zero

      END DO

      DO jCF = 1, nCF
      DO iCF = 1, nCF

        dU(iCF,1,iX1,iX2,iX3) &
          = dU(iCF,1,iX1,iX2,iX3) &
              + R_X1(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF,1,iX1,iX2,iX3)

      END DO
      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          a(2,iCF,iX1,iX2,iX3) = Zero
          b(2,iCF,iX1,iX2,iX3) = Zero
          c(2,iCF,iX1,iX2,iX3) = Zero

          DO jCF = 1, nCF

            a(2,iCF,iX1,iX2,iX3) &
              = a(2,iCF,iX1,iX2,iX3) &
                  + invR_X2(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

            b(2,iCF,iX1,iX2,iX3) &
              = b(2,iCF,iX1,iX2,iX3) &
                  + BetaTVD * invR_X2(iCF,jCF,iX1,iX2,iX3) &
                      * ( U_M(1,jCF,iX1,iX2  ,iX3) &
                        - U_M(1,jCF,iX1,iX2-1,iX3) )

            c(2,iCF,iX1,iX2,iX3) &
              = c(2,iCF,iX1,iX2,iX3) &
                  + BetaTVD * invR_X2(iCF,jCF,iX1,iX2,iX3) &
                      * ( U_M(1,jCF,iX1,iX2+1,iX3) &
                        - U_M(1,jCF,iX1,iX2  ,iX3) )

          END DO

        END DO

        IF( iX2 .EQ. iX_B0(2) .AND. ExcludeInnerGhostCell(2) )THEN

          DO iCF = 1, nCF

            b(2,iCF,iX1,iX2,iX3) = c(2,iCF,iX1,iX2,iX3)

          END DO

        END IF

        IF( iX2 .EQ. iX_E0(2) .AND. ExcludeOuterGhostCell(2) )THEN

          DO iCF = 1, nCF

            c(2,iCF,iX1,iX2,iX3) = b(2,iCF,iX1,iX2,iX3)

          END DO

        END IF

        DO iCF = 1, nCF

          dU_C(iCF,2,iX1,iX2,iX3) &
            = MinModB( a(2,iCF,iX1,iX2,iX3), &
                       b(2,iCF,iX1,iX2,iX3), &
                       c(2,iCF,iX1,iX2,iX3), &
                       dX2(iX2), BetaTVB )

          dU(iCF,2,iX1,iX2,iX3) = Zero

        END DO

        DO jCF = 1, nCF
        DO iCF = 1, nCF

          dU(iCF,2,iX1,iX2,iX3) &
            = dU(iCF,2,iX1,iX2,iX3) &
                + R_X2(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF,2,iX1,iX2,iX3)

        END DO
        END DO

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iCF = 1, nCF

          a(3,iCF,iX1,iX2,iX3) = Zero
          b(3,iCF,iX1,iX2,iX3) = Zero
          c(3,iCF,iX1,iX2,iX3) = Zero

          DO jCF = 1, nCF

            a(3,iCF,iX1,iX2,iX3) &
              = a(3,iCF,iX1,iX2,iX3) &
                  + invR_X3(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

            b(3,iCF,iX1,iX2,iX3) &
              = b(3,iCF,iX1,iX2,iX3) &
                  + BetaTVD * invR_X3(iCF,jCF,iX1,iX2,iX3) &
                      * ( U_M(1,jCF,iX1,iX2,iX3  ) &
                        - U_M(1,jCF,iX1,iX2,iX3-1) )

            c(3,iCF,iX1,iX2,iX3) &
              = c(3,iCF,iX1,iX2,iX3) &
                  + BetaTVD * invR_X3(iCF,jCF,iX1,iX2,iX3) &
                      * ( U_M(1,jCF,iX1,iX2,iX3+1) &
                        - U_M(1,jCF,iX1,iX2,iX3  ) )

          END DO

        END DO

        IF( iX3 .EQ. iX_B0(3) .AND. ExcludeInnerGhostCell(3) )THEN

          DO iCF = 1, nCF

            b(3,iCF,iX1,iX2,iX3) = c(3,iCF,iX1,iX2,iX3)

          END DO

        END IF

        IF( iX3 .EQ. iX_E0(3) .AND. ExcludeOuterGhostCell(3) )THEN

          DO iCF = 1, nCF

            c(3,iCF,iX1,iX2,iX3) = b(3,iCF,iX1,iX2,iX3)

          END DO

        END IF

        DO iCF = 1, nCF

          dU_C(iCF,3,iX1,iX2,iX3) &
            = MinModB( a(3,iCF,iX1,iX2,iX3), &
                       b(3,iCF,iX1,iX2,iX3), &
                       c(3,iCF,iX1,iX2,iX3), &
                       dX3(iX3), BetaTVB )

          dU(iCF,3,iX1,iX2,iX3) = Zero

        END DO

        DO jCF = 1, nCF
        DO iCF = 1, nCF

          dU(iCF,3,iX1,iX2,iX3) &
            = dU(iCF,3,iX1,iX2,iX3) &
                + R_X3(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF,3,iX1,iX2,iX3)

        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_M, LimitedCell, dU, SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      LimitedCell(iCF,iX1,iX2,iX3) = .FALSE.

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      ! --- Compare Limited Slopes to Original Slopes ---

      SlopeDifference(iCF,iX1,iX2,iX3) &
        = ABS( U_M(2,iCF,iX1,iX2,iX3) - dU(iCF,1,iX1,iX2,iX3) )

      IF( nDimsX .GT. 1 ) &
        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(3,iCF,iX1,iX2,iX3) - dU(iCF,2,iX1,iX2,iX3) ) )

      IF( nDimsX .GT. 2 ) &
        SlopeDifference(iCF,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCF,iX1,iX2,iX3), &
                 ABS( U_M(4,iCF,iX1,iX2,iX3) - dU(iCF,3,iX1,iX2,iX3) ) )

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      IF( SlopeDifference(iCF,iX1,iX2,iX3) &
            .GT. SlopeTolerance * ABS( U_M(1,iCF,iX1,iX2,iX3) ) )THEN

        DO iNX = 2, nDOFX

          U_M(iNX,iCF,iX1,iX2,iX3) = Zero

        END DO

        U_M(2,iCF,iX1,iX2,iX3) = dU(iCF,1,iX1,iX2,iX3)

        IF( nDimsX .GT. 1 ) U_M(3,iCF,iX1,iX2,iX3) = dU(iCF,2,iX1,iX2,iX3)

        IF( nDimsX .GT. 2 ) U_M(4,iCF,iX1,iX2,iX3) = dU(iCF,3,iX1,iX2,iX3)

        LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SL_LimitCells )

    ! --- Apply Conservative Correction ---

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL TimersStart_Euler( Timer_Euler_SL_Mapping )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_BE1, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_SL_Mapping )

    CALL TimersStart_Euler( Timer_Euler_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$OMP               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$OMP               dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$OMP               a, b, c, SlopeDifference, LimitedCell, &
    !$OMP               R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$OMP               G_X, G_K, dU_C )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, dX1, dX2, dX3, &
    !$ACC               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$ACC               dU, SqrtGm, Vol, U_X, U_K, U_N, U_M, &
    !$ACC               a, b, c, SlopeDifference, LimitedCell, &
    !$ACC               R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$ACC               G_X, G_K, dU_C )
#endif

    CALL TimersStop_Euler( Timer_Euler_SL_CopyOut )

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_Characteristic


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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell, U_M )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, SqrtGm, Vol, U_K, LimitedCell, U_M )
#endif

    ! --- Applies a correction to the 0-th order ---
    ! --- mode to maintain the cell average.     ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, LimitedCell, WeightsX_q, LegendreX, SqrtGm, &
    !$ACC          U_M, Vol, U_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
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


END MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL
