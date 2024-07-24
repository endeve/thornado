MODULE MHD_SlopeLimiterModule_Relativistic_IDEAL

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
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne
  USE MHD_BoundaryConditionsModule, ONLY: &
    ApplyInnerBC_MHD, &
    ApplyOuterBC_MHD, &
    iApplyBC_MHD_Both, &
    ApplyBoundaryConditions_MHD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_MHD_Relativistic_IDEAL
  PUBLIC :: FinalizeSlopeLimiter_MHD_Relativistic_IDEAL
  PUBLIC :: ApplySlopeLimiter_MHD_Relativistic_IDEAL

  LOGICAL      :: UseSlopeLimiter
  LOGICAL      :: UseConservativeCorrection
  CHARACTER(4) :: SlopeLimiterMethod

  ! --- TVD Limiter ---

  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: I_6x6(1:6,1:6)

  ! --- Conservative Correction ---

  REAL(DP), ALLOCATABLE :: LegendreX(:,:)
  REAL(DP), ALLOCATABLE :: Kij_X    (:,:)


CONTAINS


  SUBROUTINE InitializeSlopeLimiter_MHD_Relativistic_IDEAL &
    ( UseSlopeLimiter_Option,           &
      SlopeLimiterMethod_Option,        &
      BetaTVD_Option,                   &
      BetaTVB_Option,                   &
      SlopeTolerance_Option,            &
      UseConservativeCorrection_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in),     OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option, &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in),     OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: &
      SlopeLimiterMethod_Option

    INTEGER  :: i, j, iPol, iNX, iNX1, iNX2, iNX3, qX1, qX2, qX3
    LOGICAL  :: Verbose

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

    UseConservativeCorrection = .TRUE.
    IF( PRESENT( UseConservativeCorrection_Option ) ) &
      UseConservativeCorrection = UseConservativeCorrection_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Slope Limiter (MHD, Relativistic, IDEAL)'
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

      WRITE(*,*)
      WRITE(*,'(A6,A27,L1)'       ) '', 'UseConservativeCorrection: ' , &
        UseConservativeCorrection

    END IF

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

  END SUBROUTINE InitializeSlopeLimiter_MHD_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_MHD_Relativistic_IDEAL &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    REAL(DP), INTENT(in)           :: t
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
    INTEGER  :: iNX, iX1, iX2, iX3, iCM
    INTEGER  :: iApplyBC(3)

    INTEGER  :: nX_BE0, nX_BE1, nCM_BE0, nCM_BE1, nGF_BE0

    REAL(DP) :: U_N(1:nDOFX,1:nCM,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: U_X(1:nDOFX,1:nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: SqrtGm(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol(              iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCM        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: U_M(1:nDOFX,1:nCM,iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3))

    REAL(DP) :: a1(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b1(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c1(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: a2(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b2(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c2(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: a3(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: b3(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: c3(1:nCM,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))

    REAL(DP) :: dU_X1(1:nCM,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X2(1:nCM,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: dU_X3(1:nCM,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    REAL(DP) :: SlopeDifference(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    LOGICAL :: LimitedCell(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    REAL(DP) :: G_X(1:nDOFX,1:8  ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: G_K(1:8          ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X1(1:nCM,nCM   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X1(1:nCM,nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X2(1:nCM,nCM   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X2(1:nCM,nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X3(1:nCM,nCM   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X3(1:nCM,nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nX_BE0  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nX_BE1  = PRODUCT( iX_E1 - iX_B1 + 1 )
    nCM_BE0 = nCM * nX_BE0
    nCM_BE1 = nCM * nX_BE1
    nGF_BE0 = 8   * nX_BE0

    iApplyBC = iApplyBC_MHD_Both
    IF( PRESENT( iApplyBC_Option ) ) &
       iApplyBC = iApplyBC_Option

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_MHD &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U )

    ExcludeInnerGhostCell = .FALSE.
    ExcludeOuterGhostCell = .FALSE.
    IF( ApplyInnerBC_MHD( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeInnerGhostCell(1) = .TRUE.
    IF( ApplyOuterBC_MHD( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeOuterGhostCell(1) = .TRUE.
    IF( ApplyInnerBC_MHD( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeInnerGhostCell(2) = .TRUE.
    IF( ApplyOuterBC_MHD( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeOuterGhostCell(2) = .TRUE.
    IF( ApplyInnerBC_MHD( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeInnerGhostCell(3) = .TRUE.
    IF( ApplyOuterBC_MHD( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeOuterGhostCell(3) = .TRUE.

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCM = 1, nCM
    DO iNX = 1, nDOFX

      U_N(iNX,iCM,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM
    DO iNX = 1, nDOFX

      U_X(iNX,iCM,iX1,iX2,iX3) &
        = U_N(iNX,iCM,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute volumes of compute cells ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nX_BE0 , One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    ! --- Compute cell integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCM_BE0, One, U_X   , nDOFX, &
             WeightsX_q, 1, Zero, U_K, 1 )

    ! --- Form cell averages ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      U_K(iCM,iX1,iX2,iX3) = U_K(iCM,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    ! --- Map Nodal to Modal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCM_BE1, nDOFX, One, Kij_X, nDOFX, &
             U_N, nDOFX, Zero, U_M, nDOFX )

    CALL ComputeMinModArguments &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_M, &
             ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
             a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      dU_X1(iCM,iX1,iX2,iX3) &
        = MinModB( a1(iCM,iX1,iX2,iX3), &
                   b1(iCM,iX1,iX2,iX3), &
                   c1(iCM,iX1,iX2,iX3), &
                   dX1(iX1), BetaTVB )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        dU_X2(iCM,iX1,iX2,iX3) &
          = MinModB( a2(iCM,iX1,iX2,iX3), &
                     b2(iCM,iX1,iX2,iX3), &
                     c2(iCM,iX1,iX2,iX3), &
                     dX2(iX2), BetaTVB )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        dU_X3(iCM,iX1,iX2,iX3) &
          = MinModB( a3(iCM,iX1,iX2,iX3), &
                     b3(iCM,iX1,iX2,iX3), &
                     c3(iCM,iX1,iX2,iX3), &
                     dX3(iX3), BetaTVB )

      END DO
      END DO
      END DO
      END DO

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      SlopeDifference(iCM,iX1,iX2,iX3) &
        = ABS( U_M(2,iCM,iX1,iX2,iX3) - dU_X1(iCM,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        SlopeDifference(iCM,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCM,iX1,iX2,iX3), &
                 ABS( U_M(3,iCM,iX1,iX2,iX3) - dU_X2(iCM,iX1,iX2,iX3) ) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        SlopeDifference(iCM,iX1,iX2,iX3) &
          = MAX( SlopeDifference(iCM,iX1,iX2,iX3), &
                 ABS( U_M(4,iCM,iX1,iX2,iX3) - dU_X3(iCM,iX1,iX2,iX3) ) )

      END DO
      END DO
      END DO
      END DO

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      LimitedCell(iCM,iX1,iX2,iX3) = .FALSE.

      IF( SlopeDifference(iCM,iX1,iX2,iX3) &
            .GT. SlopeTolerance * ABS( U_M(1,iCM,iX1,iX2,iX3) ) )THEN

        DO iNX = 2, nDOFX

          U_M(iNX,iCM,iX1,iX2,iX3) = Zero

        END DO

        U_M(2,iCM,iX1,iX2,iX3) = dU_X1(iCM,iX1,iX2,iX3)

        IF( nDimsX .GT. 1 ) U_M(3,iCM,iX1,iX2,iX3) = dU_X2(iCM,iX1,iX2,iX3)

        IF( nDimsX .GT. 2 ) U_M(4,iCM,iX1,iX2,iX3) = dU_X3(iCM,iX1,iX2,iX3)

        LimitedCell(iCM,iX1,iX2,iX3) = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCM_BE1, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCM = 1, nCM
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCM) = U_N(iNX,iCM,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ApplySlopeLimiter_MHD_Relativistic_IDEAL


  SUBROUTINE FinalizeSlopeLimiter_MHD_Relativistic_IDEAL

    DEALLOCATE( Kij_X )
    DEALLOCATE( LegendreX )

  END SUBROUTINE FinalizeSlopeLimiter_MHD_Relativistic_IDEAL


  SUBROUTINE ApplyConservativeCorrection &
    ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      SqrtGm(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Vol   (        iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K   (1:nCM  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    LOGICAL, INTENT(in)     :: &
      LimitedCell(1:nCM,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: &
      U_M(1:nDOFX,1:nCM,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER  :: iNX, iX1, iX2, iX3, iCM, iDimX
    REAL(DP) :: Correction, Term

    IF( .NOT. UseConservativeCorrection ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      IF( LimitedCell(iCM,iX1,iX2,iX3) )THEN

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
                + U_M(iDimX+1,iCM,iX1,iX2,iX3) * Term / Vol(iX1,iX2,iX3)

        END DO

        U_M(1,iCM,iX1,iX2,iX3) = U_K(iCM,iX1,iX2,iX3) - Correction

      END IF

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyConservativeCorrection


  SUBROUTINE ComputeMinModArguments &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U_M, &
      ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
      a1, b1, c1, a2, b2, c2, a3, b3, c3 )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    LOGICAL,  INTENT(in)  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    REAL(DP), INTENT(in)  :: U_M(nDOFX,nCM,iX_B1(1):iX_E1(1), &
                                           iX_B1(2):iX_E1(2), &
                                           iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(out) :: a1(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b1(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c1(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: a2(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b2(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c2(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: a3(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: b3(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: c3(nCM,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3))

    INTEGER :: iX1, iX2, iX3, iCM

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      a1(iCM,iX1,iX2,iX3) = U_M(2,iCM,iX1,iX2,iX3)

      b1(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1  ,iX2,iX3) &
                                       - U_M(1,iCM,iX1-1,iX2,iX3) )
      c1(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1+1,iX2,iX3) &
                                       - U_M(1,iCM,iX1  ,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        a2(iCM,iX1,iX2,iX3) = U_M(3,iCM,iX1,iX2,iX3)

        b2(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1,iX2  ,iX3) &
                                         - U_M(1,iCM,iX1,iX2-1,iX3) )
        c2(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1,iX2+1,iX3) &
                                         - U_M(1,iCM,iX1,iX2  ,iX3) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        a3(iCM,iX1,iX2,iX3) = U_M(4,iCM,iX1,iX2,iX3)

        b3(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1,iX2,iX3  ) &
                                         - U_M(1,iCM,iX1,iX2,iX3-1) )
        c3(iCM,iX1,iX2,iX3) = BetaTVD * (  U_M(1,iCM,iX1,iX2,iX3+1) &
                                         - U_M(1,iCM,iX1,iX2,iX3  ) )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( ExcludeInnerGhostCell(1) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iCM = 1, nCM

        b1(iCM,iX_B0(1),iX2,iX3) = c1(iCM,iX_B0(1),iX2,iX3)

      END DO
      END DO
      END DO

    END IF

    IF( ExcludeOuterGhostCell(1) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iCM = 1, nCM

        c1(iCM,iX_E0(1),iX2,iX3) = b1(iCM,iX_E0(1),iX2,iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 1 .AND. ExcludeInnerGhostCell(2) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        b2(iCM,iX1,iX_B0(2),iX3) = c2(iCM,iX1,iX_B0(2),iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 1 .AND. ExcludeOuterGhostCell(2) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        c2(iCM,iX1,iX_E0(2),iX3) = b2(iCM,iX1,iX_E0(2),iX3)

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 .AND. ExcludeInnerGhostCell(3) )THEN

      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        b3(iCM,iX1,iX2,iX_B0(3)) = c3(iCM,iX1,iX2,iX_B0(3))

      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 .AND. ExcludeOuterGhostCell(3) )THEN

      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iCM = 1, nCM

        c3(iCM,iX1,iX2,iX_E0(3)) = b3(iCM,iX1,iX2,iX_E0(3))

      END DO
      END DO
      END DO

    END IF

  END SUBROUTINE ComputeMinModArguments


END MODULE MHD_SlopeLimiterModule_Relativistic_IDEAL
