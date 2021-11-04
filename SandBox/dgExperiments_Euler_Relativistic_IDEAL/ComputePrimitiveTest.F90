PROGRAM ComputePrimitiveTest

  USE KindModule, ONLY: &
    DP,   &
    Zero, &
    Half, &
    One,  &
    Four, &
    Pi,   &
    SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer,  &
    Gram,       &
    Second,     &
    Erg,        &
    Centimeter, &
    SpeedOfLight
  USE GeometryFieldsModule, ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState,   &
    ComputePressureFromPrimitive
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic, &
    ComputePrimitive_Euler_Relativistic

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER,  PARAMETER :: MaxIter = 10000000
  INTEGER             :: i, ITER
  REAL(DP)            :: U(nCF), P(0:1,nPF), G(nGF)
  REAL(DP)            :: Pressure(0:1), X1, X2, Time
  REAL(DP), PARAMETER :: Gamma_IDEAL = 4.0_DP / 3.0_DP
  REAL(DP), PARAMETER :: UnitsD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitsP    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsE    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsNe   = One  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsV1   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsV2   = One / Second
  REAL(DP), PARAMETER :: UnitsV3   = One / Second
  REAL(DP), PARAMETER :: UnitsX1   = Kilometer
  REAL(DP), PARAMETER :: UnitsX2   = One
  REAL(DP), PARAMETER :: UnitsGm11 = One
  REAL(DP), PARAMETER :: UnitsGm22 = Kilometer**2
  REAL(DP), PARAMETER :: UnitsGm33 = Kilometer**2


  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  X1 = 5.0e2_DP
  X2 = 23.0e0_DP * Pi / 1.80e2_DP

  G(iGF_Gm_dd_11) = One * UnitsGm11
  G(iGF_Gm_dd_22) = X1**2 * UnitsGm22
  G(iGF_Gm_dd_33) = X1**2 * SIN( X2 ) * UnitsGm33

  Pressure(0) = 1.0e33_DP * UnitsP
  P(0,iPF_D ) = 1.0e14_DP * UnitsD
  P(0,iPF_V1) = -1.0e5     * UnitsV1
  P(0,iPF_V2) = -3.0e2_DP  * UnitsV2
  P(0,iPF_V3) = 3.0e2_DP  * UnitsV3
  P(0,iPF_E ) = Pressure(0) / ( Gamma_IDEAL - One )
  P(0,iPF_Ne) = Zero      * UnitsNe

  WRITE(*,*)
  WRITE(*,'(A,ES10.3E3,A)') 'X1 = ', X1, ' km'
  WRITE(*,'(A,ES10.3E3,A)') 'X2 = ', X2 * 1.80e2_DP / Pi, ' deg'
  WRITE(*,*)
  WRITE(*,'(A,ES23.16E3)') 'GF_Gm11 = ', G(iGF_Gm_dd_11) / UnitsGm11
  WRITE(*,'(A,ES23.16E3)') 'GF_Gm22 = ', G(iGF_Gm_dd_22) / UnitsGm22
  WRITE(*,'(A,ES23.16E3)') 'GF_Gm33 = ', G(iGF_Gm_dd_33) / UnitsGm33

  CALL ComputeConserved_Euler_Relativistic &
         ( P(0,iPF_D ), &
           P(0,iPF_V1), &
           P(0,iPF_V2), &
           P(0,iPF_V3), &
           P(0,iPF_E ), &
           P(0,iPF_Ne), &
           U(  iCF_D ), &
           U(  iCF_S1), &
           U(  iCF_S2), &
           U(  iCF_S3), &
           U(  iCF_E ), &
           U(  iCF_Ne), &
           G(  iGF_Gm_dd_11), &
           G(  iGF_Gm_dd_22), &
           G(  iGF_Gm_dd_33), &
           Pressure(0) )

  WRITE(*,*)
  WRITE(*,'(A)') 'Inputs:'
  WRITE(*,'(A,ES24.16E3)') 'Pressure = ', Pressure(0) / UnitsP
  WRITE(*,'(A,ES24.16E3)') 'PF_D     = ', P(0,iPF_D ) / UnitsD
  WRITE(*,'(A,ES24.16E3)') 'PF_V1    = ', P(0,iPF_V1) / UnitsV1
  WRITE(*,'(A,ES24.16E3)') 'PF_V2    = ', P(0,iPF_V2) / UnitsV2
  WRITE(*,'(A,ES24.16E3)') 'PF_V3    = ', P(0,iPF_V3) / UnitsV3
  WRITE(*,'(A,ES24.16E3)') 'PF_E     = ', P(0,iPF_E ) / UnitsE
  WRITE(*,'(A,ES24.16E3)') 'PF_Ne    = ', P(0,iPF_Ne) / UnitsNe
  WRITE(*,*)
  WRITE(*,'(A,ES23.16E3,A)') &
    '|V| = ', SQRT( G(iGF_Gm_dd_11) * P(0,iPF_V1)**2 &
                      + G(iGF_Gm_dd_22) * P(0,iPF_V2)**2 &
                      + G(iGF_Gm_dd_33) * P(0,iPF_V3)**2 ) / UnitsV1, &
    ' km/s'

  WRITE(*,*)

  Time = MPI_WTIME()

  DO ITER = 1, MaxIter

    CALL ComputePrimitive_Euler_Relativistic &
           ( U(  iCF_D ), &
             U(  iCF_S1), &
             U(  iCF_S2), &
             U(  iCF_S3), &
             U(  iCF_E ), &
             U(  iCF_Ne), &
             P(1,iPF_D ), &
             P(1,iPF_V1), &
             P(1,iPF_V2), &
             P(1,iPF_V3), &
             P(1,iPF_E ), &
             P(1,iPF_Ne), &
             G(  iGF_Gm_dd_11), &
             G(  iGF_Gm_dd_22), &
             G(  iGF_Gm_dd_33) )

    CALL ComputePressureFromPrimitive &
           ( P(1,iPF_D), P(1,iPF_E), P(1,iPF_Ne), Pressure(1) )

  END DO

  Time = MPI_WTIME() - Time

  WRITE(*,'(A)') 'Outputs:'
  WRITE(*,'(A,ES24.16E3)') &
    'dP/P   = ', ( Pressure(1) - Pressure(0) ) / Pressure(0)
  WRITE(*,'(A,ES24.16E3)') &
    'dD/D   = ', ( P(1,iPF_D ) - P(0,iPF_D ) ) / P(0,iPF_D )
  WRITE(*,'(A,ES24.16E3)') &
    'dV1/V1 = ', ( P(1,iPF_V1) - P(0,iPF_V1) ) / P(0,iPF_V1)
  WRITE(*,'(A,ES24.16E3)') &
    'dV2/V2 = ', ( P(1,iPF_V2) - P(0,iPF_V2) ) / P(0,iPF_V2)
  WRITE(*,'(A,ES24.16E3)') &
    'dV3/V3 = ', ( P(1,iPF_V3) - P(0,iPF_V3) ) / P(0,iPF_V3)
  WRITE(*,'(A,ES24.16E3)') &
    'dE/E   = ', ( P(1,iPF_E ) - P(0,iPF_E ) ) / P(0,iPF_E )
  WRITE(*,'(A,ES24.16E3)') &
    'dNe/Ne = ', ( P(1,iPF_Ne) - P(0,iPF_Ne) ) / MAX( P(0,iPF_Ne), SqrtTiny )
  WRITE(*,*)
  WRITE(*,'(A,ES24.16E3,A)') 'Time per iteration: ', Time / MaxIter, ' s'
  WRITE(*,*)

  CALL FinalizeEquationOfState

END PROGRAM ComputePrimitiveTest
