PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, TwoPi, SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Evolve
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, uDF, &
    iCF_D
  USE RadiationFieldsModule, ONLY: &
    uCR
  USE InitializationModule, ONLY: &
    InitializeFields
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE TwoMoment_PositivityLimiterModule_Old, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    SolveGravity_Newtonian_Poseidon
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeTimeStep_Euler_NonRelativistic, &
    ComputeFromConserved_Euler_NonRelativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeTimeStep_TwoMoment
  USE TimeSteppingModule_CCSN, ONLY: &
    UpdateFields

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: EosTableName
  CHARACTER(64) :: OpacityTableName_AbEm
  CHARACTER(64) :: OpacityTableName_Iso
  CHARACTER(64) :: OpacityTableName_NES
  CHARACTER(64) :: OpacityTableName_Pair
  CHARACTER(64) :: ProgenitorFileName
  LOGICAL       :: wrt
  LOGICAL       :: EvolveEuler
  LOGICAL       :: EvolveTwoMoment
  LOGICAL       :: UseSlopeLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_TwoMoment
  LOGICAL       :: RampTimeStep
  INTEGER       :: RestartFileNumber
  INTEGER       :: nNodes, nSpecies
  INTEGER       :: nX(3), bcX(3), nE, bcE
  INTEGER       :: iCycle, iCycleD, nEquidistantX
  REAL(DP)      :: xL(3), xR(3), eL, eR
  REAL(DP)      :: dEquidistantX, zoomE, CFL, RampFactor
  REAL(DP)      :: t, t_wrt, t_end, dt, dt_wrt
  REAL(DP)      :: dt_Fluid, dt_Neutrinos
  REAL(DP)      :: dt_Initial, dt_Ramp

  ProgramName = 'CoreCollapseSupernova'

  CoordinateSystem = 'SPHERICAL'

  EosTableName          = 'wl-EOS-SFHo-15-25-50.h5'
  OpacityTableName_AbEm = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
  OpacityTableName_Iso  = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
  OpacityTableName_NES  = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
  OpacityTableName_Pair = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'

  ProgenitorFileName = '../Progenitors/WH07_15M_Sun.h5'

  RestartFileNumber = - 1

  ! --- Evolution Parameters ---

  t_end   = 3.50d+2 * Millisecond
  dt_wrt  = 1.00d-0 * Millisecond
  wrt     = .FALSE.
  iCycleD = 10

  ! --- Position Space Grid Parameters ---

  nX  = [ 400, 1, 1 ]
  xL  = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ]
  xR  = [ 8.0d3 * Kilometer, Pi   , TwoPi ]
  bcX = [ 30, 0, 0 ]

  nEquidistantX = 100
  dEquidistantX = 0.5d0 * Kilometer

  ! --- Energy Space Grid Parameters ---

  nE  = 16
  eL  = 0.0d0 * MeV
  eR  = 3.0d2 * MeV
  bcE = 10

  zoomE = 1.25_DP

  ! --- Time Step Control ---

  RampTimeStep = .TRUE.
  RampFactor   = 1.1_DP
  dt_Initial   = 1.0d-6 * Millisecond

  ! --- Solvers Parameters ---

  nNodes          = 2
  nSpecies        = 2
  CFL             = 0.5_DP / ( Two * DBLE( nNodes - 1 ) + One )

  EvolveEuler                = .TRUE.
  UseSlopeLimiter_Euler      = .TRUE.
  UsePositivityLimiter_Euler = .TRUE.

  EvolveTwoMoment                = .TRUE.
  UsePositivityLimiter_TwoMoment = .TRUE.

  TimeSteppingScheme = 'IMEX_PDARS'

  ! --- Auxiliary Initialization ---

  CALL InitializeDriver

  ! --- Initialize Fields ---

  CALL InitializeFields( ProgenitorFileName )

  IF( RestartFileNumber .LT. 0 )THEN

    t = Zero

    ! --- Apply Limiters to Initial Condition ---

    CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

    ! --- Solve for Gravitational Potential ---

    CALL SolveGravity_Newtonian_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
             uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,iCF_D) )

    ! --- Write Initial Condition ---

    CALL ComputeFromConserved_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    CALL WriteFieldsHDF &
           ( Time = 0.0_DP, &
             WriteGF_Option = .TRUE., &
             WriteFF_Option = .TRUE., &
             WriteRF_Option = .TRUE. )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadGF_Option = .TRUE., &
             ReadFF_Option = .TRUE., &
             ReadRF_Option = .TRUE. )

    ! --- Solve for Gravitational Potential (To Fill Geometry Ghost Cells) ---

    CALL SolveGravity_Newtonian_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
             uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,iCF_D) )

  END IF

  ! --- Evolve ---

  t_wrt  = t + dt_wrt
  wrt    = .FALSE.
  iCycle = 0

  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    ! --- Compute Time Step ---

    IF( RampTimeStep )THEN

      IF( iCycle == 1 )THEN

        dt_Ramp = dt_Initial

      ELSE

        dt_Ramp = RampFactor * dt

      END IF

    ELSE

      dt_Ramp = HUGE( One )

    END IF

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, CFL, dt_Fluid )

    CALL ComputeTimeStep_TwoMoment &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, CFL, dt_Neutrinos )

    IF( EvolveTwoMoment )THEN

      dt = MIN( dt_Fluid, dt_Neutrinos, dt_Ramp )

    ELSE

      dt = MIN( dt_Fluid, dt_Ramp )

    END IF

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
        '', 'Cycle = ', iCycle, &
        '', 't = ',  t / Millisecond, ' ms ', &
        'dt = ', dt / Millisecond, ' ms'
      WRITE(*,*)
      WRITE(*,'(A10,A11,ES13.6E3,A4,A15,ES13.6E3,A4,A10,ES13.6E3,A3)') &
        '', 'dt_Fluid = ', dt_Fluid / Millisecond, ' ms ', &
        'dt_Neutrinos = ', dt_Neutrinos / Millisecond, ' ms ', &
        'dt_Ramp = ', dt_Ramp / Millisecond, ' ms'
      WRITE(*,*)

    END IF

    CALL TimersStart( Timer_Evolve )

    CALL UpdateFields( dt, uGE, uGF, uCF, uCR )

    CALL TimersStop( Timer_Evolve )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL WriteFieldsHDF &
           ( Time = t, &
             WriteGF_Option = .TRUE., &
             WriteFF_Option = .TRUE., &
             WriteRF_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Auxiliary Finalization ---

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE ProgramInitializationModule, ONLY: &
         InitializeProgram
    USE TimersModule, ONLY: &
      InitializeTimers, &
      FinalizeTimers
    USE MeshModule, ONLY: &
      MeshX, &
      CreateMesh_Custom
    USE ReferenceElementModuleX, ONLY: &
      InitializeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      InitializeReferenceElementX_Lagrange
    USE GeometryComputationModule, ONLY: &
      ComputeGeometryX
    USE ReferenceElementModuleE, ONLY: &
      InitializeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      InitializeReferenceElementE_Lagrange
    USE GeometryComputationModuleE, ONLY: &
      ComputeGeometryE
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange
    USE EquationOfStateModule_TABLE, ONLY: &
      InitializeEquationOfState_TABLE, &
      Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y
    USE OpacityModule_TABLE, ONLY: &
      InitializeOpacities_TABLE
    USE TwoMoment_ClosureModule, ONLY: &
      InitializeClosure_TwoMoment
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_PositivityLimiterModule_Old, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
      InitializeGravitySolver_Newtonian_Poseidon
    USE TimeSteppingModule_CCSN, ONLY: &
      InitializeTimeStepping

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 0, 0 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             nE_Option &
               = nE, &
             swE_Option &
               = 1, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             zoomE_Option &
               = zoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             ActivateUnits_Option &
               = .TRUE., &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

    CALL InitializeTimers

    CALL CreateMesh_Custom &
           ( MeshX(1), nX(1), nNodes, 1, xL(1), xR(1), &
             nEquidistantX, dEquidistantX, Verbose_Option = .TRUE. )

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

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

    ! --- Initialize Equation of State ---

    CALL InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = EosTableName, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Opacities ---

    CALL InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option &
               = TRIM( OpacityTableName_AbEm ), &
             OpacityTableName_Iso_Option  &
               = TRIM( OpacityTableName_Iso ), &
             OpacityTableName_NES_Option &
               = TRIM( OpacityTableName_NES ), &
             OpacityTableName_Pair_Option &
               = TRIM( OpacityTableName_Pair ), &
             EquationOfStateTableName_Option &
               = TRIM( EosTableName ), &
             Verbose_Option = .TRUE. )

    ! --- Initialize Moment Closure ---

    CALL InitializeClosure_TwoMoment

    ! --- Initialize Slope Limiter (Euler) ---

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
           ( BetaTVD_Option &
               = 1.75_DP, &
             BetaTVB_Option &
               = Zero, &
             SlopeTolerance_Option &
               = 1.0d-6, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_Euler, &
             UseCharacteristicLimiting_Option &
               = .FALSE., &
             UseTroubledCellIndicator_Option &
               = .FALSE., &
             LimiterThresholdParameter_Option &
               = Zero )

    ! --- Initialize Positivity Limiter (Euler) ---

    CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option &
               = UsePositivityLimiter_Euler, &
             Verbose_Option &
               = .TRUE., &
             Min_1_Option &
               = ( One + EPSILON( One ) ) * Min_D, &
             Min_2_Option &
               = ( One + EPSILON( One ) ) * Min_T, &
             Min_3_Option &
               = ( One + EPSILON( One ) ) * Min_Y, &
             Max_1_Option &
               = ( One - EPSILON( One ) ) * Max_D, &
             Max_2_Option &
               = ( One - EPSILON( One ) ) * Max_T, &
             Max_3_Option &
               = ( One - EPSILON( One ) ) * Max_Y )

    ! --- Initialize Positivity Limiter (Two-Moment) ---

    CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option &
             = SqrtTiny, &
           Max_1_Option &
             = HUGE( One ), &
           Min_2_Option &
             = SqrtTiny, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter_TwoMoment )

    ! --- Initialize Gravity Solver ---

    CALL InitializeGravitySolver_Newtonian_Poseidon

    ! --- Initialize Time Stepping ---

    CALL InitializeTimeStepping &
         ( Scheme_Option &
             = TimeSteppingScheme, &
           EvolveEuler_Option &
             = EvolveEuler, &
           EvolveTwoMoment_Option &
             = EvolveTwoMoment )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE ReferenceElementModuleX, ONLY: &
      FinalizeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      FinalizeReferenceElementX_Lagrange
    USE ReferenceElementModuleE, ONLY: &
      FinalizeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      FinalizeReferenceElementE_Lagrange
    USE ReferenceElementModuleZ, ONLY: &
      FinalizeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      FinalizeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      FinalizeReferenceElement_Lagrange
    USE EquationOfStateModule_TABLE, ONLY: &
      FinalizeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      FinalizeOpacities_TABLE
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_PositivityLimiterModule_Old, ONLY: &
      FinalizePositivityLimiter_TwoMoment
    USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
      FinalizeGravitySolver_Newtonian_Poseidon
    USE TimeSteppingModule_CCSN, ONLY: &
      FinalizeTimeStepping
    USE TimersModule, ONLY: &
      FinalizeTimers

    CALL FinalizeTimers

    CALL FinalizeTimeStepping

    CALL FinalizeGravitySolver_Newtonian_Poseidon

    CALL FinalizePositivityLimiter_TwoMoment

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeEquationOfState_TABLE

    CALL FinalizeOpacities_TABLE

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver
