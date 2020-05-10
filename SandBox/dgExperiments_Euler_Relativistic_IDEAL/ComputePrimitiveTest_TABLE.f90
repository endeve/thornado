PROGRAM ComputePrimitiveTest_TABLE

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
    ComputePrimitive_Euler_Relativistic_TABLE_Scalar, &
    ComputeConserved_Euler_Relativistic

  IMPLICIT NONE

  INTEGER,  PARAMETER :: nNodes = 1
  INTEGER             :: i, ITER
  REAL(DP)            :: U(nNodes,nCF), P(nNodes,nPF), G(nNodes,nGF)
  REAL(DP)            :: Pressure(nNodes), X1, X2, X3
  REAL(DP), PARAMETER :: Gamma_IDEAL = 4.0_DP / 3.0_DP
  REAL(DP), PARAMETER :: UnitsD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitsP    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsE    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsNe   = One  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsV1   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsV2   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsV3   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsX1   = Kilometer
  REAL(DP), PARAMETER :: UnitsX2   = Kilometer
  REAL(DP), PARAMETER :: UnitsX3   = Kilometer
  REAL(DP), PARAMETER :: UnitsGm11 = Kilometer
  REAL(DP), PARAMETER :: UnitsGm22 = Kilometer
  REAL(DP), PARAMETER :: UnitsGm33 = Kilometer


  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  i = 1

  X1 = 5.0e2_DP
  X2 = 5.0e2_DP
  X3 = 5.0e2_DP

  G(i,iGF_Gm_dd_11) = One * UnitsGm11
  G(i,iGF_Gm_dd_22) = One * UnitsGm22
  G(i,iGF_Gm_dd_33) = One * UnitsGm33

  Pressure(i) = 1.0_DP    * UnitsP
  P(i,iPF_D ) = 1.0e-6_DP * UnitsD
  P(i,iPF_V1) = -0.9_DP    * UnitsV1
  P(i,iPF_V2) = 0.1_DP    * UnitsV2
  P(i,iPF_V3) = 0.1_DP    * UnitsV3
  P(i,iPF_E ) = Pressure(i) / ( Gamma_IDEAL - One )
  P(i,iPF_Ne) = Zero      * UnitsNe

  WRITE(*,*)
  WRITE(*,*) '|V| = ', SQRT( G(i,iGF_Gm_dd_11) * P(i,iPF_V1)**2 &
                               + G(i,iGF_Gm_dd_22) * P(i,iPF_V2)**2 &
                               + G(i,iGF_Gm_dd_33) * P(i,iPF_V3)**2 ) / UnitsV1
  WRITE(*,*)
  WRITE(*,'(A,ES24.16E3)') 'GF_Gm11 = ', G(i,iGF_Gm_dd_11) / UnitsGm11
  WRITE(*,'(A,ES24.16E3)') 'GF_Gm22 = ', G(i,iGF_Gm_dd_22) / UnitsGm22
  WRITE(*,'(A,ES24.16E3)') 'GF_Gm33 = ', G(i,iGF_Gm_dd_33) / UnitsGm33

  CALL ComputeConserved_Euler_Relativistic &
         ( P(i,iPF_D ), &
           P(i,iPF_V1), &
           P(i,iPF_V2), &
           P(i,iPF_V3), &
           P(i,iPF_E ), &
           P(i,iPF_Ne), &
           U(i,iCF_D ), &
           U(i,iCF_S1), &
           U(i,iCF_S2), &
           U(i,iCF_S3), &
           U(i,iCF_E ), &
           U(i,iCF_Ne), &
           G(i,iGF_Gm_dd_11), &
           G(i,iGF_Gm_dd_22), &
           G(i,iGF_Gm_dd_33), &
           Pressure(i) )

  WRITE(*,*)
  WRITE(*,*) 'Inputs:'
  WRITE(*,'(A,ES24.16E3)') 'PF_D  = ', P(i,iPF_D ) / UnitsD
  WRITE(*,'(A,ES24.16E3)') 'PF_V1 = ', P(i,iPF_V1) / UnitsV1
  WRITE(*,'(A,ES24.16E3)') 'PF_V2 = ', P(i,iPF_V2) / UnitsV2
  WRITE(*,'(A,ES24.16E3)') 'PF_V3 = ', P(i,iPF_V3) / UnitsV3
  WRITE(*,'(A,ES24.16E3)') 'PF_E  = ', P(i,iPF_E ) / UnitsE
  WRITE(*,'(A,ES24.16E3)') 'PF_Ne = ', P(i,iPF_Ne) / UnitsNe
  WRITE(*,*)

  DO ITER = 1, 1

    CALL ComputePrimitive_Euler_Relativistic_TABLE_Scalar &
           ( U(i,iCF_D ), &
             U(i,iCF_S1), &
             U(i,iCF_S2), &
             U(i,iCF_S3), &
             U(i,iCF_E ), &
             U(i,iCF_Ne), &
             G(i,iGF_Gm_dd_11), &
             G(i,iGF_Gm_dd_22), &
             G(i,iGF_Gm_dd_33), &
             P(i,iPF_D ), &
             P(i,iPF_V1), &
             P(i,iPF_V2), &
             P(i,iPF_V3), &
             P(i,iPF_E ), &
             P(i,iPF_Ne) )

    WRITE(*,*) 'Outputs:'
    WRITE(*,'(A,ES24.16E3)') 'PF_D  = ', P(i,iPF_D ) / UnitsD
    WRITE(*,'(A,ES24.16E3)') 'PF_V1 = ', P(i,iPF_V1) / UnitsV1
    WRITE(*,'(A,ES24.16E3)') 'PF_V2 = ', P(i,iPF_V2) / UnitsV2
    WRITE(*,'(A,ES24.16E3)') 'PF_V3 = ', P(i,iPF_V3) / UnitsV3
    WRITE(*,'(A,ES24.16E3)') 'PF_E  = ', P(i,iPF_E ) / UnitsE
    WRITE(*,'(A,ES24.16E3)') 'PF_Ne = ', P(i,iPF_Ne) / UnitsNe

  END DO

  CALL FinalizeEquationOfState

END PROGRAM ComputePrimitiveTest_TABLE
