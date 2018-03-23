PROGRAM DeleptonizationWave

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE_Beta, ONLY: &
    ComputeGeometryE
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR, nSpecies
  USE ClosureModule_M1, ONLY: &
    InitializeMomentClosure
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule, ONLY: &
    InitializeOpacities, &
    FinalizeOpacities
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter, &
    FinalizePositivityLimiter
  USE TimeSteppingModule_IMEX_RK, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE dgDiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_M1_DG_Implicit

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'DeleptonizationWave', &
           nX_Option &
             = [ 64, 64, 01 ], &
           swX_Option &
             = [ 01, 01, 01 ], &
           bcX_Option &
             = [ 32, 32, 01 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, - 0.5d0 ] * Kilometer, &
           xR_Option &
             = [ 1.0d2, 1.0d2, + 0.5d0 ] * Kilometer, &
           nE_Option &
             = 16, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 1.0d2 * MeV, &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Moment Closure ---

  CALL InitializeMomentClosure

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5' )

  ! --- Initialize Opacities ---

  CALL InitializeOpacities &
         ( Opacity_Option &
             = 'TABLE', &
           OpacityTableName_Option &
             = 'OpacityTable.h5' )

  ! --- Create Neutrino Opacities ---

  CALL CreateNeutrinoOpacities &
         ( nZ, [ nNodesZ(1), 1, 1, 1 ], nSpecies )

  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter &
         ( Min_1_Option = 0.0d-00, &
           Max_1_Option = 1.0d-00, &
           Min_2_Option = 0.0d-00, &
           UsePositivityLimiter_Option &
             = .TRUE. )

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK &
         ( Scheme = 'IMEX_PC2' )

  ! --- Set Initial Condition ---

  CALL InitializeFields

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- 

  CALL ComputeIncrement_M1_DG_Implicit &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, 1.0d-4, uGE, uGF, uCR, rhsCR )

  ! --- Finalize ---

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeOpacities

  CALL DestroyNeutrinoOpacities

  CALL FinalizePositivityLimiter

  CALL Finalize_IMEX_RK

  CALL FinalizeProgram

END PROGRAM DeleptonizationWave
