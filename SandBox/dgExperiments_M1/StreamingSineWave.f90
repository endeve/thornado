PROGRAM StreamingSineWave

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
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
  USE ReferenceElementModule_Beta, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE TimeSteppingModule_IMEX_RK, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE ClosureModule_M1, ONLY: &
    InitializeMomentClosure
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE_Beta, ONLY: &
    ComputeGeometryE
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE InitializationModule, ONLY: &
    InitializeFields
  USE dgDiscretizationModule, ONLY: &
    ComputeIncrement_M1_DG_Explicit
  USE dgDiscretizationModule_Collisions, ONLY: &
    InitializeCollisions, &
    FinalizeCollisions, &
    ComputeIncrement_M1_DG_Implicit, &
    ComputeCorrection_M1_DG_Implicit

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER, PARAMETER :: nE = 1, nX(3) = [ 32, 1, 1 ], nNodes = 2
  REAL(DP) :: t, dt, t_end, wTime

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StreamingSineWave', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 1.0d0, 1.0d0, 1.0d0 ], &
           nE_Option &
             = nE, &
           eL_Option &
             = 0.0d0, &
           eR_Option &
             = 1.0d0, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           Opacity_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = 0.0_DP )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK( 'IMEX_P_A2' )

  ! --- Initialize Moment Closure ---

  CALL InitializeMomentClosure

  ! --- Initialize Implicit Solver ---

  CALL InitializeCollisions &
         ( N0_Option = 0.0d0, SigmaA_Option = 0.0d0, SigmaS_Option = 0.0d0 )

  ! --- Set Initial Condition ---

  CALL InitializeFields &
         ( Direction_Option = 'X' )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF( Time = 0.0_DP, WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t         = 0.0d-0
  dt        = 0.2_DP / (2.0_DP*DBLE(nNodes-1)+1.0_DP) / DBLE(nX(1))
  t_end     = 1.0d+0
  iCycleD   = 10
  iCycleW   = 10
  maxCycles = 10000

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCR, &
             ComputeIncrement_M1_DG_Explicit, &
             ComputeIncrement_M1_DG_Implicit, &
             ComputeCorrection_M1_DG_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

    END IF

  END DO

  CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL Finalize_IMEX_RK

  CALL FinalizeCollisions

  CALL FinalizeProgram

END PROGRAM StreamingSineWave
