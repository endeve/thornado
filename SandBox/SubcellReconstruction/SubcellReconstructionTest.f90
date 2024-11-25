PROGRAM SubcellReconstructionTest

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, nDimsX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE SubcellReconstructionModule, ONLY: &
    InitializeSubcellReconstruction, &
    FinalizeSubcellReconstruction, &
    ProjectionMatrix, &
    ReconstructionMatrix
  USE UtilitiesModule, ONLY: &
    WriteMatrix

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'SubcellReconstructionTest', &
           nX_Option &
             = [ 2, 2, 1 ], &
           swX_Option &
             = [ 0, 0, 0 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0_DP, &
           eR_Option &
             = 1.0_DP, &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL InitializeSubcellReconstruction

  CALL WriteMatrix &
         ( nDOFX, nDOFX, ProjectionMatrix, 'ProjectionMatrix.dat' )

  CALL WriteMatrix &
         ( nDOFX, nDOFX, ReconstructionMatrix, 'ReconstructionMatrix.dat' )

  CALL FinalizeSubcellReconstruction

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

END PROGRAM SubcellReconstructionTest
