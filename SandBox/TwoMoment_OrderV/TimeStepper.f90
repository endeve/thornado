PROGRAM TimeStepper

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    nDOF, nDOFX, nDOFE
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers
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
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeConserved_TwoMoment
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK

  IMPLICIT NONE

  CHARACTER(16) :: Scheme
  INTEGER       :: nNodes, nE, nX(3)
  INTEGER       :: iS, iZ1, iZ2, iZ3, iZ4, iNode, iNodeX
  REAL(DP)      :: eL, eR, xL(3), xR(3)

  nNodes = 2

  nX = [ 08, 08, 08 ]
  xL = [ Zero, Zero, Zero ]
  xR = [ One,  One,  One  ]

  nE = 16
  eL = Zero
  eR = One

  Scheme = 'IMEX_ARS_111'

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'TimeStepper', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nE_Option &
             = nE, &
           swE_Option &
             = 1, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           nSpecies_Option &
             = 1, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Initialize Timers ---

  CALL InitializeTimers

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

  CALL InitializeClosure_TwoMoment

  ! --- Initialize Dummy Radiation Field ---

  DO iS  = 1, nSpecies
  DO iZ4 = iZ_B1(4), iZ_E1(4)
  DO iZ3 = iZ_B1(3), iZ_E1(3)
  DO iZ2 = iZ_B1(2), iZ_E1(2)
  DO iZ1 = iZ_B1(1), iZ_E1(1)

    DO iNode = 1, nDOF

      iNodeX = MOD( (iNode-1) / nDOFE, nDOFX ) + 1

      uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_D, iS) = 0.5_DP
      uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
      uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
      uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP

      CALL ComputeConserved_TwoMoment &
             ( uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_D, iS), &
               uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
               uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
               uPR(iNode,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
               uCR(iNode,iZ1,iZ2,iZ3,iZ4,iCR_N, iS), &
               uCR(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
               uCR(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
               uCR(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
               Zero, Zero, Zero, &
               uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
               uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
               uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

    END DO

  END DO
  END DO
  END DO
  END DO
  END DO

  CALL Initialize_IMEX_RK( Scheme )

  CALL Update_IMEX_RK( One, uGE, uGF, uCF, uCR )

  CALL Finalize_IMEX_RK

  ! --- Finalize ---

  CALL FinalizeTimers

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM TimeStepper
