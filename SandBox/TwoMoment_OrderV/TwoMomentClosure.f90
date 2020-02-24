PROGRAM TwoMomentClosure

  USE KindModule, ONLY: &
    DP, Zero, One
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment, &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeEddingtonTensorComponents_ud, &
    ComputeHeatFluxTensorComponents_udd

  IMPLICIT NONE

  INTEGER,  PARAMETER :: nPoints = 2**14
  REAL(DP), PARAMETER :: MinD = 1.0d-02
  REAL(DP), PARAMETER :: MaxD = 1.0d-00

  INTEGER  :: i
  REAL(DP) :: absI, nVec(3)
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

    ! --- 

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


END PROGRAM TwoMomentClosure
