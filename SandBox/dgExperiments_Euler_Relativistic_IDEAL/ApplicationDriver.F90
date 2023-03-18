PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Pi, &
    TwoPi
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
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1, &
    nDimsX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_Relativistic, ONLY: &
    InitializeFields_Relativistic
  USE Euler_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_Relativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_Relativistic_IDEAL, &
    ApplySlopeLimiter_Euler_Relativistic_IDEAL
  USE Euler_PositivityLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_Relativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_Relativistic_IDEAL, &
    ApplyPositivityLimiter_Euler_Relativistic_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_Euler_Relativistic, &
    ComputeTimeStep_Euler_Relativistic
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, &
    uAF, &
    uDF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE Euler_TallyModule_Relativistic, ONLY: &
    InitializeTally_Euler_Relativistic, &
    FinalizeTally_Euler_Relativistic, &
    ComputeTally_Euler_Relativistic
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler, &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize, &
    Timer_Euler_Finalize

  USE MPI

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: AdvectionProfile
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: SuppressTally
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseConservativeCorrection
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: Mass = Zero
  REAL(DP)      :: ZoomX(3)
  REAL(DP)      :: Timer_Evolution

  ! --- Sedov--Taylor blast wave ---
  REAL(DP) :: Eblast
  INTEGER  :: nDetCells

  LOGICAL  :: WriteGF = .TRUE., WriteFF = .TRUE.
  LOGICAL  :: ActivateUnits = .FALSE.

  SuppressTally = .FALSE.

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  ProgramName = 'Advection'
!!$  ProgramName = 'Advection2D'
!!$  ProgramName = 'RiemannProblem'
!!$  ProgramName = 'RiemannProblem2D'
!!$  ProgramName = 'RiemannProblemSpherical'
!!$  ProgramName = 'SedovTaylorBlastWave'
!!$  ProgramName = 'KelvinHelmholtzInstability'

  swX               = [ 0, 0, 0 ]
  RestartFileNumber = -1
  t                 = 0.0_DP
  ZoomX             = 1.0_DP

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'SineWave'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 64, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'Advection2D' )

      AdvectionProfile = 'SineWaveX1X2'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 1, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 32, 32, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ One / SQRT( Two ), One / SQRT( Two ), 1.0_DP ]

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'Sod' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.2_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'IsolatedShock' )

          Gamma = 4.0_DP / 3.0_DP
          t_end = 2.0e1_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'IsolatedContact' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 2.0e1_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem1' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.4_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem4' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'PerturbedShockTube' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.35_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'ShockReflection' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.75_DP
          bcX   = [ 23, 0, 0 ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid RiemannProblemName: ', RiemannProblemName
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  Sod'
          WRITE(*,'(A)')     '  IsolatedShock'
          WRITE(*,'(A)')     '  IsolatedContact'
          WRITE(*,'(A)')     '  MBProblem1'
          WRITE(*,'(A)')     '  MBProblem4'
          WRITE(*,'(A)')     '  PerturbedShockTube'
          WRITE(*,'(A)')     '  ShockReflection'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 128, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblem2D' )

      RiemannProblemName = 'DzB2002'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'DzB2002' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4_DP
          bcX   = [ 2, 2, 0 ]

        CASE( 'IsolatedShock' )

          Gamma = 4.0_DP / 3.0_DP
          t_end = 25.0_DP
          bcX   = [ 2, 2, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 32, 32, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblemSpherical' )

      RiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 0.5_DP
      bcX = [ 2, 0, 0 ]

      nX  = [ 256, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 2.0_DP, Pi, TwoPi ]

    CASE( 'SedovTaylorBlastWave' )

      nDetCells = 1
      Eblast    = 1.0e-3_DP

      CoordinateSystem = 'SPHERICAL'

      Gamma = 4.0_DP / 3.0_DP
      t_end = 1.0_DP
      bcX = [ 3, 0, 0 ]

      nX  = [ 256, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.2_DP, Pi, TwoPi ]

    CASE( 'KelvinHelmholtzInstability' )

       CoordinateSystem = 'CARTESIAN'

       Gamma = 4.0_DP / 3.0_DP
       t_end = 0.1_DP
       bcX = [ 1, 1, 1 ]

       nX = [ 16, 32, 16 ]
      swX = [ 1, 1, 1 ]
       xL = [ -0.5_DP, -1.0_DP, -0.5_DP ]
       xR = [  0.5_DP,  1.0_DP,  0.5_DP ]

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  Advection'
      WRITE(*,'(A)')     '  Advection2D'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblem2D'
      WRITE(*,'(A)')     '  RiemannProblemSpherical'
      WRITE(*,'(A)')     '  SedovTaylorBlastWave'
      WRITE(*,'(A)')     '  KelvinHelmholtzInstability'
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  ! --- DG ---

  nNodes = 1
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  ! --- Time Stepping ---

  nStagesSSPRK = 1
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  CFL = 0.5_DP ! Cockburn & Shu, (2001), JSC, 16, 173

  ! --- Slope Limiter ---

  UseSlopeLimiter           = .FALSE.
  SlopeLimiterMethod        = 'TVD'
  BetaTVD                   = 1.75_DP
  BetaTVB                   = 0.0_DP
  SlopeTolerance            = 1.0e-6_DP
  UseCharacteristicLimiting = .FALSE.
  UseTroubledCellIndicator  = .FALSE.
  LimiterThresholdParameter = 0.015_DP
  UseConservativeCorrection = .FALSE.

  ! --- Positivity Limiter ---

  UsePositivityLimiter = .FALSE.
  Min_1                = 1.0e-13_DP
  Min_2                = 1.0e-13_DP

  ! === End of User Input ===

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = ActivateUnits, &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = Mass )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
         ( UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           SlopeLimiterMethod_Option &
             = TRIM( SlopeLimiterMethod ), &
           BetaTVD_Option &
             = BetaTVD, &
           BetaTVB_Option &
             = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter, &
           UseConservativeCorrection_Option &
             = UseConservativeCorrection )

  CALL InitializePositivityLimiter_Euler_Relativistic_IDEAL &
         ( UsePositivityLimiter_Option = UsePositivityLimiter, &
           Verbose_Option = .TRUE., &
           Min_1_Option = Min_1, &
           Min_2_Option = Min_2 )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )
  WRITE(*,*)
  WRITE(*,'(A6,A,ES11.3E3)') '', 'CFL: ', CFL

  uCF = Zero ! Without this, crashes when copying data in TimeStepper
  uDF = Zero ! Without this, crashes in IO

  CALL InitializeFields_Relativistic &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ), &
           nDetCells_Option = nDetCells, &
           Eblast_Option    = Eblast )

  IF( RestartFileNumber .LT. 0 )THEN

    CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE FROM( uGF, uCF, uPF, uAF )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE HOST       ( uGF, uCF, uPF, uAF )
#endif

    CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadFF_Option = .TRUE., ReadGF_Option = .TRUE. )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE TO( uGF, uCF )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE DEVICE   ( uGF, uCF )
#endif

  END IF

  iCycleD = 10
!!$  iCycleW = 1; dt_wrt = -1.0_DP
  dt_wrt = 1.0e-2_DP * ( t_end - t ); iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'
  WRITE(*,*)

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
           SuppressTally_Option = SuppressTally )

  CALL ComputeTally_Euler_Relativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .FALSE. )

  CALL TimersStop_Euler( Timer_Euler_Initialize )

  iCycle = 0
  Timer_Evolution = MPI_WTIME()
  DO WHILE( t .LT. t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF, uCF, &
             CFL / ( nDimsX * ( Two * DBLE( nNodes ) - One ) ), &
             dt )

    IF( t + dt .LT. t_end )THEN

      t = t + dt

    ELSE

      dt = t_end - t
      t  = t_end

    END IF

    IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN

      WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
        'Cycle: ', iCycle, ' t = ', t / UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel ), &
        ' dt = ', dt /  UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel )

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, &
             ComputeIncrement_Euler_DG_Explicit )

    IF( iCycleW .GT. 0 )THEN

      IF( MOD( iCycle, iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.

    ELSE

      IF( t + dt .GT. t_wrt )THEN

        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.

      END IF

    END IF

    IF( wrt )THEN

      CALL TimersStart_Euler( Timer_Euler_InputOutput )

      CALL ComputeFromConserved_Euler_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE FROM( uGF, uCF, uPF, uAF )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE HOST       ( uGF, uCF, uPF, uAF )
#endif

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

      CALL ComputeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
             Verbose_Option = .FALSE. )

      wrt = .FALSE.

      CALL TimersStop_Euler( Timer_Euler_InputOutput )

    END IF

  END DO

  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,I8.8,A,ES10.3E3,A)') &
    'Finished ', iCycle, ' cycles in ', Timer_Evolution, ' s'
  WRITE(*,*)

  CALL TimersStart_Euler( Timer_Euler_Finalize )

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE FROM( uGF, uCF, uPF, uAF )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE HOST       ( uGF, uCF, uPF, uAF )
#endif

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL ComputeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

  CALL FinalizeTally_Euler_Relativistic

  CALL FinalizeFluid_SSPRK

  CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeEquationOfState

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

  CALL TimersStop_Euler( Timer_Euler_Finalize )

  CALL FinalizeTimers_Euler

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
