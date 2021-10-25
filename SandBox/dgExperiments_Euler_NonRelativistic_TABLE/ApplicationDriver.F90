PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Two, Pi, TwoPi
  USE UnitsModule, ONLY: &
    Second, &
    Millisecond, &
    Microsecond, &
    Kilometer
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, swX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE MeshModule, ONLY: &
    MeshX, &
    CreateMesh_Custom
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE EquationOfStateModule_TABLE, ONLY: &
    MinD, MaxD, MinT, MaxT, MinY, MaxY
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, uPF, uAF, uDF, &
    ResetFluidFields_Diagnostic
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    InitializeGravitySolver_Newtonian_Poseidon, &
    FinalizeGravitySolver_Newtonian_Poseidon, &
    SolveGravity_Newtonian_Poseidon
  USE InitializationModule, ONLY: &
    InitializeFields
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    InitializeSlopeLimiter_Euler_NonRelativistic_TABLE, &
    FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE, &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    InitializePositivityLimiter_Euler_NonRelativistic_TABLE, &
    FinalizePositivityLimiter_Euler_NonRelativistic_TABLE, &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic, &
    ComputeTimeStep_Euler_NonRelativistic
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE Euler_TallyModule_NonRelativistic, ONLY: &
    InitializeTally_Euler_NonRelativistic, &
    FinalizeTally_Euler_NonRelativistic, &
    ComputeTally_Euler_NonRelativistic
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32)  :: ProgramName
  CHARACTER(32)  :: AdvectionProfile
  CHARACTER(32)  :: RiemannProblemName
  CHARACTER(32)  :: CoordinateSystem
  CHARACTER(128) :: EosTableName
  CHARACTER(32)  :: ProgenitorFile
  LOGICAL        :: CustomRadialGrid
  LOGICAL        :: wrt
  LOGICAL        :: UseSlopeLimiter
  LOGICAL        :: UseCharacteristicLimiting
  LOGICAL        :: UseTroubledCellIndicator
  LOGICAL        :: UsePositivityLimiter
  LOGICAL        :: SelfGravity
  INTEGER        :: iCycle, iCycleD
  INTEGER        :: RestartFileNumber
  INTEGER        :: nX(3), bcX(3), nNodes, nStages
  INTEGER        :: nEquidistantElements
  REAL(DP)       :: t, dt, t_end, dt_wrt, t_wrt, wTime
  REAL(DP)       :: xL(3), xR(3), zoomX(3), dxEquidistant
  REAL(DP)       :: BetaTVD, BetaTVB
  REAL(DP)       :: LimiterThresholdParameter

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler

  ProgramName = 'RiemannProblem'

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ProgenitorFile = '../Progenitors/WH07_15M_Sun.h5'

  CustomRadialGrid = .FALSE.

  SelfGravity = .FALSE.

  RestartFileNumber = -1

  t = 0.0_DP

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'SineWave'

      CoordinateSystem = 'CARTESIAN'

      nX = [ 16, 01, 01 ]
      xL = [ -1.0d2, 0.0d0, 0.0d0 ] * Kilometer
      xR = [  1.0d2, 1.0d2, 1.0d2 ] * Kilometer
      zoomX = One

      bcX = [ 1, 1, 1 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .FALSE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.0d-0
      UsePositivityLimiter      = .FALSE.

      iCycleD = 10
      t_end   = 1.0d1 * ( 1.0d5 / SpeedOfLightMKS ) * Second
      dt_wrt  = 1.0d-0 * t_end

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Pochik'

      CoordinateSystem = 'CARTESIAN'

      nX = [ 256, 1, 1 ]
      xL = [ - 5.0_DP,   0.0_DP, 0.0_DP ] * Kilometer
      xR = [ + 5.0_DP, + 1.0_DP, 1.0_DP ] * Kilometer
      zoomX = One

      bcX = [ 2, 0, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .FALSE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 1.0d-2
      UsePositivityLimiter      = .TRUE.

      iCycleD = 1

      SELECT CASE( TRIM( RiemannProblemName ) )
      CASE( 'Sod' )

        t_end  = 2.5d-2 * Millisecond
        dt_wrt = 5.0d-2 * t_end

      CASE( 'Pochik' )

        t_end  = 2.0d-1 * Millisecond
        dt_wrt = 5.0d-2 * t_end

      END SELECT

   CASE( 'RiemannProblemSpherical' )

      RiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      nX = [ 128, 16, 1 ]
      xL = [ 1.0d-3 * Kilometer, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP * Kilometer, Pi,     4.0_DP ]
      zoomX = One

      bcX = [ 3, 3, 0 ]

      nNodes = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 5.0d-1 * Millisecond
      dt_wrt  = 2.5d-2 * Millisecond

    CASE( 'RiemannProblemCylindrical' )

      RiemannProblemName = 'CylindricalSod'

      CoordinateSystem = 'CYLINDRICAL'

      nX = [ 100, 100, 1 ]
      xL = [    0.0_DP * Kilometer,-10.0d0 * Kilometer, 0.0_DP * Kilometer ]
      xR = [   10.0_DP * Kilometer, 10.0d0 * Kilometer, TwoPi ]
      zoomX = One

      bcX = [ 3, 3, 0 ]

      nNodes = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.01_DP

      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 2.5d-2 * Millisecond
      dt_wrt  = 1.25d-4 * Millisecond

    CASE( 'Jet' )

      CoordinateSystem = 'CARTESIAN'

      nX = [ 100, 100, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      zoomX = One

      bcX = [ 2, 2, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 1.5d-0
      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 2.5d-1 * Millisecond
      dt_wrt  = 2.5d-6 * Millisecond

    CASE( 'Implosion' )

      CoordinateSystem = 'CARTESIAN'

      nX = [ 64, 64, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR = [ 0.3_DP, 0.3_DP, 1.0_DP ] * Kilometer
      zoomX = One

      bcX = [ 3, 3, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.50_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.
      UsePositivityLimiter      = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.30_DP

      iCycleD = 1
      t_end   = 2.500_DP * Millisecond
      dt_wrt  = 0.045_DP * Millisecond

    CASE( 'GravitationalCollapse' )

      CoordinateSystem = 'SPHERICAL'

      nX = [ 256, 1, 1 ]
      xL = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ]
      xR = [ 8.0d3 * Kilometer, Pi,      TwoPi ]
      zoomX = [ 1.020059256924853_DP, 1.0_DP, 1.0_DP ]

      CustomRadialGrid     = .FALSE.
      nEquidistantElements = 100
      dxEquidistant        = 0.5_DP * Kilometer

      bcX = [30, 0, 0]

      nNodes  = 2
      nStages = 2

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.
      UsePositivityLimiter      = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.01_DP

      SelfGravity = .TRUE.

      iCycleD = 10
      t_end   = 301.5d0  * Millisecond
      dt_wrt  = 5.0d-1 * Millisecond

    CASE( 'ShockEntropyWave' )

      CoordinateSystem = 'CARTESIAN'

      nX = [ 256, 1, 1 ]
      xL = [ - 5.0_DP,   0.0_DP, 0.0_DP ] * Kilometer
      xR = [ + 5.0_DP, + 1.0_DP, 1.0_DP ] * Kilometer
      zoomX = One

      bcX = [ 2, 0, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 2.0_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .FALSE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 1.5d-0
      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 7.5d-2 * Millisecond
      dt_wrt  = 1.25d-4 * Millisecond

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  Advection'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblemSpherical'
      WRITE(*,'(A)')     '  RiemannProblemCylindrical'
      WRITE(*,'(A)')     '  Jet'
      WRITE(*,'(A)')     '  Implosion'
      WRITE(*,'(A)')     '  GravitationalCollapse'
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 0, 0 ], & ! 1 1 0
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           zoomX_Option &
             = zoomX, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  IF( CustomRadialGrid )THEN

    CALL CreateMesh_Custom &
           ( MeshX(1), nX(1), nNodes, 1, xL(1), xR(1), &
             nEquidistantElements, dxEquidistant, Verbose_Option = .TRUE. )

  END IF

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = TRIM( EosTableName ) )

  CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
         ( BetaTVD_Option = BetaTVD, &
           BetaTVB_Option = BetaTVB, &
           SlopeTolerance_Option &
             = 1.0d-6, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter )

  CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
         ( UsePositivityLimiter_Option = UsePositivityLimiter, &
           Verbose_Option = .TRUE., &
           Min_1_Option = ( One + EPSILON(One) ) * MinD, &
           Min_2_Option = ( One + EPSILON(One) ) * MinT, &
           Min_3_Option = ( One + EPSILON(One) ) * MinY, &
           Max_1_Option = ( One - EPSILON(One) ) * MaxD, &
           Max_2_Option = ( One - EPSILON(One) ) * MaxT, &
           Max_3_Option = ( One - EPSILON(One) ) * MaxY )

  CALL InitializeFields &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ), &
           ProgenitorFileName_Option &
             = TRIM( ProgenitorFile ) )

  IF( RestartFileNumber .GE. 0 )THEN

    CALL ReadFieldsHDF( RestartFileNumber, t, ReadFF_Option = .TRUE. )

  END IF

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET UPDATE TO( uCF )
#elif defined( THORNADO_OACC   )
      !$ACC UPDATE DEVICE( uCF )
#endif

  CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  IF( SelfGravity )THEN

    CALL InitializeGravitySolver_Newtonian_Poseidon

    CALL SolveGravity_Newtonian_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
             uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,iCF_D) )
  END IF

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeTally_Euler_NonRelativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .TRUE. )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
             CFL = 0.5_DP / ( Two * DBLE( nNodes - 1 ) + One ), TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t / Millisecond, ' ms ', &
          'dt = ', dt / Millisecond, ' ms'

    END IF

    IF( SelfGravity )THEN

      CALL UpdateFluid_SSPRK &
            ( t, dt, uGF, uCF, uDF, &
              ComputeIncrement_Euler_DG_Explicit, &
              SolveGravity_Newtonian_Poseidon )

    ELSE

      CALL UpdateFluid_SSPRK &
            ( t, dt, uGF, uCF, uDF, ComputeIncrement_Euler_DG_Explicit )

    END IF

    t = t + dt

    IF( wrt )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET UPDATE FROM( uGF, uCF )
#elif defined( THORNADO_OACC   )
      !$ACC UPDATE HOST( uGF, uCF )
#endif

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
               Verbose_Option = .TRUE. )

      CALL ResetFluidFields_Diagnostic( nX, swX, uDF )

      wrt = .FALSE.

    END IF

  END DO

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET UPDATE FROM( uGF, uCF )
#elif defined( THORNADO_OACC   )
      !$ACC UPDATE HOST( uGF, uCF )
#endif

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
           Verbose_Option = .TRUE. )

  CALL FinalizeTally_Euler_NonRelativistic

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

  CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

  IF( SelfGravity )THEN

    CALL FinalizeGravitySolver_Newtonian_Poseidon

  END IF

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

  CALL FinalizeTimers_Euler

!  WRITE(*,*)
!  WRITE(*,'(2x,A)') 'git info'
!  WRITE(*,'(2x,A)') '--------'
!  WRITE(*,*)
!  WRITE(*,'(2x,A)') 'git branch:'
!  CALL EXECUTE_COMMAND_LINE( 'git branch' )
!  WRITE(*,*)
!  WRITE(*,'(2x,A)') 'git describe --tags:'
!  CALL EXECUTE_COMMAND_LINE( 'git describe --tags' )
!  WRITE(*,*)
!  WRITE(*,'(2x,A)') 'git rev-parse HEAD:'
!  CALL EXECUTE_COMMAND_LINE( 'git rev-parse HEAD' )
!  WRITE(*,*)
!  WRITE(*,'(2x,A)') 'date:'
!  CALL EXECUTE_COMMAND_LINE( 'date' )
!  WRITE(*,*)

END PROGRAM ApplicationDriver
