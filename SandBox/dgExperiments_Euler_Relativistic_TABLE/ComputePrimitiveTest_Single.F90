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
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    rhoMin_Euler_GR
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
  REAL(DP)     :: X1_C, X2_C, X3_C, dX1, dX2, dX3
  REAL(DP)     :: G(nGF), U(nCF), P(nPF)
  CHARACTER(2) :: iDimX

  CHARACTER(256) :: ProgramName      = 'ComputePrimitiveTest_Single'
  CHARACTER(256) :: CoordinateSystem = 'SPHERICAL'
  CHARACTER(256) :: EosTableName
  LOGICAL        :: ActivateUnits    = .TRUE.

  iErr         = 0
  ITERATION    = 0

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = ActivateUnits, &
           BasicInitialization_Option &
             = .TRUE. )

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'TABLE', &
           EquationOfStateTableName_Option &
             = TRIM( EosTableName ), &
           Verbose_Option = .TRUE. )

  rhoMin_Euler_GR = Min_D

  ! --- Dummies ---

  ITERATION = 00000044
  iNX       = 000000021
  iDimX     = 'X1'
  iX1       = 1
  iX2       = 1
  iX3       = 1
  iX_B0     = [00000,00001,00001]
  iX_E0     = [00001,00001,00001]
  X1_C      =  Zero
  X2_C      =  Zero
  X3_C      =  Zero
  dX1       =  Zero
  dX2       =  Zero
  dX3       =  Zero
  U(iCF_D       ) = +9.3419877783595740E+011_DP * ( Gram / Centimeter**3 )
  U(iCF_S1      ) = -3.0052578151128483E+023_DP * ( Gram / Centimeter**2 / Second )
  U(iCF_S2      ) = +0.0000000000000000E+000_DP
  U(iCF_S3      ) = +0.0000000000000000E+000_DP
  U(iCF_E       ) = +6.3220455732175546E+031_DP * ( Erg / Centimeter**3 )
  U(iCF_Ne      ) = +2.6496152478424499E+035_DP * ( One / Centimeter**3 )
  G(iGF_Gm_dd_11) = +1.3573939057635909E+000_DP
  G(iGF_Gm_dd_22) = +1.2946440572402196E+002_DP * Kilometer**2
  G(iGF_Gm_dd_33) = +1.2946440572402196E+002_DP * Kilometer**2

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
