PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Half, One, Two, &
    Pi, TwoPi, FourPi
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
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GravitySolutionModule_Newtonian_PointMass, ONLY: &
    ComputeGravitationalPotential
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_SlopeLimiterModule, ONLY: &
    Euler_InitializeSlopeLimiter, &
    Euler_FinalizeSlopeLimiter, &
    Euler_ApplySlopeLimiter
  USE Euler_PositivityLimiterModule, ONLY: &
    Euler_InitializePositivityLimiter, &
    Euler_FinalizePositivityLimiter, &
    Euler_ApplyPositivityLimiter
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields, ApplyPerturbations
  USE Euler_UtilitiesModule, ONLY: &
    Euler_ComputeFromConserved, &
    Euler_ComputeTimeStep
  USE Euler_dgDiscretizationModule, ONLY: &
    Euler_ComputeIncrement_DG_Explicit
  USE Euler_TallyModule, ONLY: &
    Euler_InitializeTally, &
    Euler_FinalizeTally, &
    Euler_ComputeTally

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP), PARAMETER :: mDot   = FourPi
  REAL(DP), PARAMETER :: Mass   = Half
  REAL(DP), PARAMETER :: rShock = One
  REAL(DP), PARAMETER :: Gamma  = 4.0_DP / 3.0_DP
  REAL(DP), PARAMETER :: Mach   = 1.0d2

  LOGICAL  :: wrt
  INTEGER  :: iCycle, iCycleD
  INTEGER  :: nX(3), nNodes
  REAL(DP) :: t, dt, t_end, dt_wrt, t_wrt, wTime

  nX     = [ 256, 1, 1 ]
  nNodes = 2

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StandingAccretionShock', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 11, 3, 1 ], &
           xL_Option &
             = [ 0.2d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 2.0d0, Pi,    TwoPi ], &
           zoomX_Option &
             = [ 1.0025_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           BasicInitialization_Option &
             = .TRUE. )

  t_end   = 5.0d+1
  dt_wrt  = 2.5d-1
  iCycleD = 100

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL ComputeGravitationalPotential &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass )

  CALL InitializeFluid_SSPRK( nStages = 2 )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL Euler_InitializeSlopeLimiter &
         ( BetaTVD_Option                   = 1.15_DP, &
           BetaTVB_Option                   = 0.00_DP, &
           SlopeTolerance_Option            = 1.00d-6, &
           UseSlopeLimiter_Option           = .TRUE., &
           UseCharacteristicLimiting_Option = .TRUE., &
           UseTroubledCellIndicator_Option  = .TRUE., &
           UseConservativeCorrection_Option = .TRUE., &
           LimiterThresholdParameter_Option = 0.015_DP )

  CALL Euler_InitializePositivityLimiter &
         ( Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12, &
           UsePositivityLimiter_Option = .TRUE. )

  CALL InitializeFields &
         ( mDot, Mass, rShock, Gamma, Mach )

  CALL ApplyPerturbations( 1.4_DP, 1.6_DP, 0, 2.0_DP )

  CALL Euler_ApplySlopeLimiter &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL Euler_ApplyPositivityLimiter &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL Euler_ComputeFromConserved &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  CALL Euler_InitializeTally &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL Euler_ComputeTimeStep &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             CFL = 0.3_DP / ( Two * DBLE( nNodes - 1 ) + One ), TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, Euler_ComputeIncrement_DG_Explicit )

    t = t + dt

    IF( wrt )THEN

      CALL Euler_ComputeFromConserved &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL Euler_ComputeTally &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL Euler_ComputeFromConserved &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL Euler_ComputeTally &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL Euler_FinalizeTally

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL Euler_FinalizePositivityLimiter

  CALL Euler_FinalizeSlopeLimiter

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
