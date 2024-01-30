PROGRAM ComputePrimitiveTest_Single

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
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

  INTEGER      :: iErr, ITERATION, iX_B0(3), iX_E0(3), iNX, iX1, iX2, iX3
  REAL(DP)     :: Gamma_IDEAL, X1_C, X2_C, X3_C, dX1, dX2, dX3
  REAL(DP)     :: G(nGF), U(nCF), P(nPF)
  CHARACTER(2) :: iDimX

  CHARACTER(256) :: ProgramName      = 'ComputePrimitiveTest_Single'
  CHARACTER(256) :: CoordinateSystem = 'SPHERICAL'
  LOGICAL        :: ActivateUnits    = .TRUE.

  iErr         = 0
  ITERATION    = 0
  Gamma_IDEAL = 1.30_DP!4.0_DP / 3.0_DP

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = ActivateUnits, &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL, &
           Verbose_Option = .TRUE. )

  ITERATION=           00000026
  iNX=  00000002
  iDimX= 'X1'
  iX1=1
  iX2=1
  iX3=1
  iX_B0=[                00000,  00001,  00001]
  iX_E0=[                00031,  00001,  00001]
  X1_C=  +2.929688E+002_DP * Kilometer
  X2_C=  +1.570796E+000_DP
  X3_C=  +3.141593E+000_DP
  dX1=   +1.953125E+002_DP * Kilometer
  dX2=   +3.141593E+000_DP
  dX3=   +6.283185E+000_DP
  U(iCF_D       ) = +7.4207342150334740E+009_DP * ( Gram / Centimeter**3 )
  U(iCF_S1      ) = -1.1594721657221886E+019_DP * ( Gram / Centimeter**2 / Second )
  U(iCF_S2      ) = +0.0000000000000000E+000_DP
  U(iCF_S3      ) = +0.0000000000000000E+000_DP
  U(iCF_E       ) = +8.8927190722882423E+027_DP * ( Erg / Centimeter**3 )
  U(iCF_Ne      ) = +0.0000000000000000E+000_DP
  G(iGF_Gm_dd_11) = +1.0179343532647396E+000_DP
  G(iGF_Gm_dd_22) = +3.8864700701499372E+004_DP * Kilometer**2
  G(iGF_Gm_dd_33) = +3.8864700701499372E+004_DP * Kilometer**2

  CALL ComputePrimitive_Euler_Relativistic &
    ( U(iCF_D ), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E), U(iCF_Ne), &
      P(iPF_D ), P(iPF_V1), P(iPF_V2), P(iPF_V3), P(iPF_E), P(iPF_Ne), &
      G(iGF_Gm_dd_11), G(iGF_Gm_dd_22), G(iGF_Gm_dd_33), &
      ITERATION_Option = Iteration, iErr_Option = iErr )

  CALL DescribeError_Euler &
    ( iErr, &
      Int_Option = [ ITERATION, iNX, iX_B0(1), iX_B0(2), iX_B0(3), &
                                     iX_E0(1), iX_E0(2), iX_E0(3), &
                                iNX, iX1     , iX2     , iX3 ], &
      Real_Option = [ X1_C, X2_C, X3_C, dX1, dX2, dX3, &
                      U(iCF_D ), &
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

  CALL FinalizeProgram

END PROGRAM ComputePrimitiveTest_Single
