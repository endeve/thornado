PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Half, One, Two, &
    Pi, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDimsX, nDOFX
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
    nCF, nPF, nAF, &
    uCF, uPF, uAF, &
    uDF
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_SlopeLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    ApplySlopeLimiter_Euler_NonRelativistic_IDEAL
  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL
  USE GeometryFieldsModule, ONLY: &
    nGF, uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields, ApplyPerturbations
  USE Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler, &
    ComputeTimeStep_Euler
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE Euler_TallyModule_NonRelativistic, ONLY: &
    InitializeTally_Euler_NonRelativistic, &
    FinalizeTally_Euler_NonRelativistic, &
    ComputeTally_Euler_NonRelativistic

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
  INTEGER  :: iErr

  nX     = [ 256, 1, 1 ]
  nNodes = 2

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StandingAccretionShock', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 0, 0 ], &
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

  CALL InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( BetaTVD_Option                   = 1.15_DP, &
           BetaTVB_Option                   = 0.00_DP, &
           SlopeTolerance_Option            = 1.00d-6, &
           UseSlopeLimiter_Option           = .TRUE., &
           UseCharacteristicLimiting_Option = .TRUE., &
           UseTroubledCellIndicator_Option  = .TRUE., &
           UseConservativeCorrection_Option = .TRUE., &
           LimiterThresholdParameter_Option = 0.015_DP )

  CALL InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( UsePositivityLimiter_Option = .TRUE., &
           Verbose_Option = .TRUE., &
           Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12 )

  CALL InitializeFields &
         ( mDot, Mass, rShock, Gamma, Mach )

  CALL ApplyPerturbations( 1.4_DP, 1.6_DP, 0, 2.0_DP )

  CALL ApplySlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeFromConserved_Euler &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeTally_Euler_NonRelativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .TRUE. )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
             CFL = 0.3_DP / ( nDimsX * DBLE( nNodes - 1 ) + One ), &
             TimeStep = dt )

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
           ( t, dt, uGF, uCF, uDF, ComputeIncrement_Euler_DG_Explicit )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
               Verbose_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_Euler &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
           Verbose_Option = .TRUE. )

  CALL FinalizeTally_Euler_NonRelativistic

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git info'
  WRITE(*,'(2x,A)') '--------'
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git branch:'
  CALL EXECUTE_COMMAND_LINE( 'git branch' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git describe --tags:'
  CALL EXECUTE_COMMAND_LINE( 'git describe --tags' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git rev-parse HEAD:'
  CALL EXECUTE_COMMAND_LINE( 'git rev-parse HEAD' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'date:'
  CALL EXECUTE_COMMAND_LINE( 'date' )
  WRITE(*,*)

END PROGRAM ApplicationDriver
