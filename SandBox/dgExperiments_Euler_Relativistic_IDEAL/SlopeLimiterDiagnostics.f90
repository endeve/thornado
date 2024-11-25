PROGRAM SlopeLimiterDiagnostics

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, TwoPi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDimsX, nDOFX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_Relativistic, ONLY: &
    InitializeFields_Relativistic
  USE Euler_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_Relativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_Relativistic_IDEAL, &
    ApplySlopeLimiter_Euler_Relativistic_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_Euler_Relativistic
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF, &
    uCF, uPF, uAF, &
    uDF
  USE GeometryFieldsModule, ONLY: &
    nGF, uGF

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  INTEGER       :: nX(3), bcX(3), nNodes
  REAL(DP)      :: xL(3), xR(3), Gamma

  LOGICAL  :: WriteGF = .FALSE., WriteFF = .TRUE.
  REAL(DP) :: Timer_Evolution

  ProgramName = 'SlopeLimiterTest'

  Gamma = 5.0_DP / 3.0_DP
  bcX = [ 2, 0, 0 ]

  nX = [ 7, 1, 1 ]
  xL = [ -1.0_DP, 0.0_DP, 0.0_DP ]
  xR = [ +1.0_DP, 1.0_DP, 1.0_DP ]

  nNodes = 2

  UseSlopeLimiter           = .TRUE.
  UseCharacteristicLimiting = .FALSE.
  UseTroubledCellIndicator  = .FALSE.

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = 0.0_DP )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeFields_Relativistic

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL InitializeFields_Relativistic

  CALL InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
         ( UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           SlopeLimiterMethod_Option &
             = 'TVD', &
           BetaTVD_Option &
             = 1.75_DP, &
           BetaTVB_Option &
             = 0.00_DP, &
           SlopeTolerance_Option &
             = 0.0_DP, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = 0.0_DP, &
           UseConservativeCorrection_Option &
             = .FALSE. )

  CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL InitializeFields_Relativistic

  CALL InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
         ( UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           SlopeLimiterMethod_Option &
             = 'WENO', &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator )

  CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM SlopeLimiterDiagnostics
