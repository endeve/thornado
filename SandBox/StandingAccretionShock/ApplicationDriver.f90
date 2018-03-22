PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP) :: wTime
  REAL(DP), PARAMETER :: Gamma = 4.0_DP / 3.0_DP

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StandingAccretionShock', &
           nX_Option &
             = [ 128, 32, 1 ], &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.2d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 2.0d0, Pi,    TwoPi ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeFields &
         ( 4.0_DP * Pi, 0.5_DP, 1.0_DP, Gamma, 3.0d2 ) 

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

!!$  ! --- Main Part of Code Will Go Here

  ! --- Finalize ---

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
