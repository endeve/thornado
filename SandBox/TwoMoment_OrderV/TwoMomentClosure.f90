PROGRAM TwoMomentClosure

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment, &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeEddingtonTensorComponents_ud, &
    ComputeHeatFluxTensorComponents_udd, &
    Flux_E

  IMPLICIT NONE

  INTEGER,  PARAMETER :: nPoints = 2**14
  REAL(DP), PARAMETER :: MinD = 1.0d-02
  REAL(DP), PARAMETER :: MaxD = 1.0d-00
  REAL(DP), PARAMETER :: MaxV = 5.0d-01 ! Max Velocity
  REAL(DP), PARAMETER :: MaxA = 1.0d-00 ! Max Velocity Derivative

  INTEGER  :: i, INFO
  REAL(DP) :: absI, nVec(3), absV
  REAL(DP) :: D(nPoints), I1(nPoints), I2(nPoints), I3(nPoints)
  REAL(DP) :: h(nPoints), k(nPoints), q(nPoints)
  REAL(DP) :: k_11(nPoints), k_12(nPoints), k_13(nPoints)
  REAL(DP) ::                k_22(nPoints), k_23(nPoints)
  REAL(DP) ::                               k_33(nPoints)
  REAL(DP) :: q_111(nPoints), q_121(nPoints), q_131(nPoints)
  REAL(DP) ::                 q_221(nPoints), q_231(nPoints)
  REAL(DP) ::                                 q_331(nPoints)
  REAL(DP) :: q_222(nPoints), q_223(nPoints)
  REAL(DP) :: q_332(nPoints), q_333(nPoints)
  REAL(DP) :: Tr_q1(nPoints), Tr_q2(nPoints), Tr_q3(nPoints)
  REAL(DP) :: V1(nPoints), V2(nPoints), V3(nPoints)
  REAL(DP) :: dVdX(9)
  REAL(DP) :: dV1dX1, dV1dX2, dV1dX3
  REAL(DP) :: dV2dX1, dV2dX2, dV2dX3
  REAL(DP) :: dV3dX1, dV3dX2, dV3dX3
  REAL(DP) :: A(3,3), Lambda(3), WORK(8)
  REAL(DP) :: LambdaMin(nPoints), LambdaMax(nPoints)
  REAL(DP) :: dVdX_Max(nPoints), VdotN(nPoints)
  REAL(DP) :: Alpha(nPoints), AlphaK(nPoints)
  REAL(DP) :: F_E(4,nPoints), S1(nPoints), S2(nPoints), S3(nPoints)

  CALL InitializeClosure_TwoMoment( Verbose_Option = .TRUE. )

  CALL RANDOM_SEED

  DO i = 1, nPoints

    ! --- Number Density ---

    CALL RANDOM_NUMBER( D(i) )

    D(i) = 10**( D(i) * LOG10( MinD ) + ( One - D(i) ) * LOG10( MaxD ) )

    ! --- Number Flux Components ---

    CALL RANDOM_NUMBER( absI )

    absI = absI * D(i)

    CALL RANDOM_NUMBER( nVec )

    nVec = 2.0_DP * ( nVec - 0.5_DP )
    nVec = nVec / SQRT( DOT_PRODUCT( nVec, nVec ) )

    I1(i) = absI * nVec(1)
    I2(i) = absI * nVec(2)
    I3(i) = absI * nVec(3)

    ! --- Eddington and Heat Flux Tensors ---

    h(i) = FluxFactor &
             ( D(i), I1(i), I2(i), I3(i), One, One, One )

    k(i) = EddingtonFactor &
             ( D(i), h(i) )

    q(i) = HeatFluxFactor &
             ( D(i), h(i) )

    CALL ComputeEddingtonTensorComponents_ud &
           ( D(i), I1(i), I2(i), I3(i), One, One, One, &
             k_11(i), k_12(i), k_13(i), k_22(i), k_23(i), k_33(i) )

    CALL ComputeHeatFluxTensorComponents_udd &
           ( D(i), I1(i), I2(i), I3(i), One, One, One, &
             q_111(i), q_121(i), q_131(i), q_221(i), q_231(i), q_331(i), &
             q_222(i), q_223(i), q_333(i), q_332(i) )

    ! --- Trace Conditions ---

    Tr_q1(i) = ( q_111(i) + q_221(i) + q_331(i) ) * D(i)

    Tr_q2(i) = ( q_121(i) + q_222(i) + q_332(i) ) * D(i)

    Tr_q3(i) = ( q_131(i) + q_223(i) + q_333(i) ) * D(i)

    ! --- Velocity Components ---

    CALL RANDOM_NUMBER( absV )

    absV = absV * MaxV

    CALL RANDOM_NUMBER( nVec )

    nVec = 2.0_DP * ( nVec - 0.5_DP )
    nVec = nVec / SQRT( DOT_PRODUCT( nVec, nVec ) )

    V1(i) = absV * nVec(1)
    V2(i) = absV * nVec(2)
    V3(i) = absV * nVec(3)

    ! --- Velocity Derivatives ---

    CALL RANDOM_NUMBER( dVdX )

    dV1dX1 = 2.0_DP * ( dVdX(1) - 0.5_DP ) * MaxA
    dV1dX2 = 2.0_DP * ( dVdX(2) - 0.5_DP ) * MaxA
    dV1dX3 = 2.0_DP * ( dVdX(3) - 0.5_DP ) * MaxA

    dV2dX1 = 2.0_DP * ( dVdX(4) - 0.5_DP ) * MaxA
    dV2dX2 = 2.0_DP * ( dVdX(5) - 0.5_DP ) * MaxA
    dV2dX3 = 2.0_DP * ( dVdX(6) - 0.5_DP ) * MaxA

    dV3dX1 = 2.0_DP * ( dVdX(7) - 0.5_DP ) * MaxA
    dV3dX2 = 2.0_DP * ( dVdX(8) - 0.5_DP ) * MaxA
    dV3dX3 = 2.0_DP * ( dVdX(9) - 0.5_DP ) * MaxA

    ! --- Alpha ---

    Alpha(i) &
      =   ( - dV1dX1 ) * k_11(i) &
        + ( - dV2dX2 ) * k_22(i) &
        + ( - dV3dX3 ) * k_33(i) &
        + ( - dV2dX1 - dV1dX2 ) * k_12(i) &
        + ( - dV3dX1 - dV1dX3 ) * k_13(i) &
        + ( - dV3dX2 - dV2dX3 ) * k_23(i)

    ! --- Eigenvalues of Quadratic Form Matrix ---

    A(:,1) = - Half * [    Two * dV1dX1, dV1dX2 + dV2dX1, dV1dX3 + dV3dX1 ]
    A(:,2) = - Half * [ dV2dX1 + dV1dX2,    Two * dV2dX2, dV2dX3 + dV3dX2 ]
    A(:,3) = - Half * [ dV3dX1 + dV1dX3, dV3dX2 + dV2dX3,    Two * dV3dX3 ]

    CALL DSYEV( 'N', 'U', 3, A, 3, Lambda, WORK, 8, INFO )

    LambdaMin(i) = MINVAL( Lambda )
    LambdaMax(i) = MAXVAL( Lambda )

    ! --- Quadratic Form ---

    CALL RANDOM_NUMBER( nVec )

    nVec = 2.0_DP * ( nVec - 0.5_DP )
    nVec = nVec / SQRT( DOT_PRODUCT( nVec, nVec ) )

    AlphaK(i) &
      =   ( - dV1dX1 ) * nVec(1) * nVec(1) &
        + ( - dV2dX2 ) * nVec(2) * nVec(2) &
        + ( - dV3dX3 ) * nVec(3) * nVec(3) &
        + ( - dV2dX1 - dV1dX2 ) * nVec(1) * nVec(2) &
        + ( - dV3dX1 - dV1dX3 ) * nVec(1) * nVec(3) &
        + ( - dV3dX2 - dV2dX3 ) * nVec(2) * nVec(3)

    ! --- Maximum Velocity Derivative ---

    dVdX_Max(i) &
      = MAXVAL( ABS( [ dV1dX1, dV1dX2, dV1dX3, &
                       dV2dX1, dV2dX2, dV2dX3, &
                       dV3dX1, dV3dX2, dV3dX3 ] ) )

    ! --- V dot N ---

    VdotN(i) = V1(i) * nVec(1) + V2(i) * nVec(2) + V3(i) * nVec(3)

    ! --- Energy Space Flux Components ---

    F_E(:,i) &
      = Flux_E( D(i), I1(i), I2(i), I3(i), &
                V1(i), V2(i), V3(i), &
                dV1dX1, dV2dX1, dV3dX1, &
                dV1dX2, dV2dX2, dV3dX2, &
                dV1dX3, dV2dX3, dV3dX3, &
                One, One, One, &
                Zero, Zero, Zero, &
                Zero, Zero, Zero, &
                Zero, Zero, Zero )

    ! --- Number Flux Source ---

    S1(i) = - F_E(2,i) - ( I1(i) * dV1dX1 + I2(i) * dV1dX2 + I3(i) * dV1dX3 )
    S2(i) = - F_E(3,i) - ( I1(i) * dV2dX1 + I2(i) * dV2dX2 + I3(i) * dV2dX3 )
    S3(i) = - F_E(4,i) - ( I1(i) * dV3dX1 + I2(i) * dV3dX2 + I3(i) * dV3dX3 )

  END DO

  CALL WriteVector( nPoints, D , 'D.dat'  )
  CALL WriteVector( nPoints, I1, 'I1.dat' )
  CALL WriteVector( nPoints, I2, 'I2.dat' )
  CALL WriteVector( nPoints, I3, 'I3.dat' )

  CALL WriteVector( nPoints, h, 'h.dat' )
  CALL WriteVector( nPoints, k, 'k.dat' )
  CALL WriteVector( nPoints, q, 'q.dat' )

  CALL WriteVector( nPoints, k_11, 'k_11.dat' )
  CALL WriteVector( nPoints, k_12, 'k_12.dat' )
  CALL WriteVector( nPoints, k_13, 'k_13.dat' )
  CALL WriteVector( nPoints, k_22, 'k_22.dat' )
  CALL WriteVector( nPoints, k_23, 'k_23.dat' )
  CALL WriteVector( nPoints, k_33, 'k_33.dat' )

  CALL WriteVector( nPoints, q_111, 'q_111.dat' )
  CALL WriteVector( nPoints, q_222, 'q_222.dat' )
  CALL WriteVector( nPoints, q_333, 'q_333.dat' )

  CALL WriteVector( nPoints, q_121, 'q_121.dat' )
  CALL WriteVector( nPoints, q_131, 'q_131.dat' )
  CALL WriteVector( nPoints, q_221, 'q_221.dat' )
  CALL WriteVector( nPoints, q_231, 'q_231.dat' )
  CALL WriteVector( nPoints, q_331, 'q_331.dat' )

  CALL WriteVector( nPoints, q_223, 'q_223.dat' )
  CALL WriteVector( nPoints, q_332, 'q_332.dat' )

  CALL WriteVector( nPoints, Tr_q1, 'Tr_q1.dat' )
  CALL WriteVector( nPoints, Tr_q2, 'Tr_q2.dat' )
  CALL WriteVector( nPoints, Tr_q3, 'Tr_q3.dat' )

  CALL WriteVector( nPoints, V1, 'V1.dat' )
  CALL WriteVector( nPoints, V2, 'V2.dat' )
  CALL WriteVector( nPoints, V3, 'V3.dat' )

  CALL WriteVector( nPoints, Alpha , 'Alpha.dat'  )
  CALL WriteVector( nPoints, AlphaK, 'AlphaK.dat' )

  CALL WriteVector( nPoints, LambdaMin, 'LambdaMin.dat' )
  CALL WriteVector( nPoints, LambdaMax, 'LambdaMax.dat' )

  CALL WriteVector( nPoints, dVdX_Max, 'dVdX_Max.dat' )
  CALL WriteVector( nPoints, VdotN, 'VdotN.dat' )

  CALL WriteVector( nPoints, F_E(1,:), 'F_E0.dat' )
  CALL WriteVector( nPoints, F_E(2,:), 'F_E1.dat' )
  CALL WriteVector( nPoints, F_E(3,:), 'F_E2.dat' )
  CALL WriteVector( nPoints, F_E(4,:), 'F_E3.dat' )

  CALL WriteVector( nPoints, S1, 'S1.dat' )
  CALL WriteVector( nPoints, S2, 'S2.dat' )
  CALL WriteVector( nPoints, S3, 'S3.dat' )

END PROGRAM TwoMomentClosure
