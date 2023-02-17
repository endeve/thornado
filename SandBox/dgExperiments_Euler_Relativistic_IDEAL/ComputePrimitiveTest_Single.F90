PROGRAM ComputePrimitiveTest_Single

  USE KindModule, ONLY: &
    DP, &
    One
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    nGF
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF, &
    iAF_P, &
    nAF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE UnitsModule, ONLY: &
    Kilometer, &
    Gram, &
    Second, &
    Erg, &
    Centimeter

  IMPLICIT NONE

  INCLUDE 'mpif.h'

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
  REAL(DP), PARAMETER :: UnitsGm11 = One
  REAL(DP), PARAMETER :: UnitsGm22 = One
  REAL(DP), PARAMETER :: UnitsGm33 = One

  INTEGER  :: iErr
  REAL(DP) :: Gamma_IDEAL
  REAL(DP) :: G(nGF), U(nCF), P(nPF)

  iErr = 0
  Gamma_IDEAL = 4.0_DP / 3.0_DP

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL, &
           Verbose_Option = .TRUE. )

  U(iCF_D       ) =  4.7052933109437282E-017_DP
  U(iCF_S1      ) = -1.5254676338844351E-017_DP
  U(iCF_S2      ) =  0.0000000000000000E+000_DP
  U(iCF_S3      ) =  0.0000000000000000E+000_DP
  U(iCF_E       ) =  2.2061209954548441E-018_DP
  U(iCF_Ne      ) =  0.0000000000000000E+000_DP
  G(iGF_Gm_dd_11) =  1.0950993657030816E+000_DP
  G(iGF_Gm_dd_22) =  8.8703048621949596E+009_DP
  G(iGF_Gm_dd_33) =  8.8703048621949596E+009_DP

  CALL ComputePrimitive_Euler_Relativistic &
    ( U(iCF_D ), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E), U(iCF_Ne), &
      P(iPF_D ), P(iPF_V1), P(iPF_V2), P(iPF_V3), P(iPF_E), P(iPF_Ne), &
      G(iGF_Gm_dd_11), G(iGF_Gm_dd_22), G(iGF_Gm_dd_33), iErr )

  CALL DescribeError_Euler &
    ( iErr, &
      Int_Option = [ 1 ], &
      Real_Option = [ U(iCF_D ), &
                      U(iCF_S1), &
                      U(iCF_S2), &
                      U(iCF_S3), &
                      U(iCF_E ), &
                      U(iCF_Ne), &
                      G(iGF_Gm_dd_11), &
                      G(iGF_Gm_dd_22), &
                      G(iGF_Gm_dd_33) ] )

  PRINT '(A,ES24.16E3)', 'PF_D : ', P(iPF_D ) / UnitsD
  PRINT '(A,ES24.16E3)', 'PF_V1: ', P(iPF_V1)
  PRINT '(A,ES24.16E3)', 'PF_E : ', P(iPF_E ) / UnitsE

  CALL FinalizeEquationOfState

END PROGRAM ComputePrimitiveTest_Single
