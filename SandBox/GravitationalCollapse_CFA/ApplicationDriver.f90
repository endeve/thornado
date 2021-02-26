PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP,   &
    Zero, &
    One,  &
    Two,  &
    Pi,   &
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
  USE EquationOfStateModule_TABLE, ONLY: &
    MinD, &
    MaxD, &
    MinT, &
    MaxT, &
    MinY, &
    MaxY
  USE ProgramHeaderModule, ONLY: &
    iX_B0,  &
    iX_B1,  &
    iX_E0,  &
    iX_E1,  &
    nDimsX, &
    nDOFX
  USE InitializationModule_Relativistic, ONLY: &
    InitializeFields_Relativistic
  USE Euler_SlopeLimiterModule_Relativistic_TABLE, ONLY: &
    InitializeSlopeLimiter_Euler_Relativistic_TABLE, &
    FinalizeSlopeLimiter_Euler_Relativistic_TABLE,   &
    ApplySlopeLimiter_Euler_Relativistic_TABLE
  USE Euler_PositivityLimiterModule_Relativistic_TABLE, ONLY: &
    InitializePositivityLimiter_Euler_Relativistic_TABLE, &
    FinalizePositivityLimiter_Euler_Relativistic_TABLE,   &
    ApplyPositivityLimiter_Euler_Relativistic_TABLE
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
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    InitializeGravitySolver_CFA_Poseidon, &
    FinalizeGravitySolver_CFA_Poseidon,   &
    SolveGravity_CFA_Poseidon
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    Time
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK,   &
    UpdateFluid_SSPRK, &
    WriteSourceTerms2
  USE UnitsModule, ONLY: &
    Kilometer,   &
    Millisecond, &
    UnitsDisplay
  USE Euler_TallyModule_Relativistic, ONLY: &
    InitializeTally_Euler_Relativistic, &
    FinalizeTally_Euler_Relativistic,   &
    ComputeTally_Euler_Relativistic
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler,            &
    InitializeTimers_Euler,  &
    FinalizeTimers_Euler,    &
    TimersStart_Euler,       &
    TimersStop_Euler,        &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize,  &
    Timer_Euler_Finalize
  USE Poseidon_UtilitiesModule, ONLY: &
    ComputeSourceTerms_Poseidon

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: EosTableName
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseConservativeCorrection
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: xL(3), xR(3)
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, t_wrt2, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: ZoomX(3)

  CHARACTER(32) :: ProgenitorFileName

  LOGICAL  :: WriteGF = .TRUE., WriteFF = .TRUE.
  LOGICAL  :: ActivateUnits = .TRUE.

  REAL(DP) :: Timer_Evolution

  REAL(DP), ALLOCATABLE :: SourceTerms_Poseidon(:,:,:,:,:)

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  ProgramName = 'GravitationalCollapse'

  ProgenitorFileName = '../Progenitors/WH07_15M_Sun.h5'

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  RestartFileNumber = -1
  t                 = 0.0_DP

  CoordinateSystem = 'SPHERICAL'

  t_end = 3.0e2_DP * Millisecond

  bcX = [ 30, 0, 0 ]

  nX    = [ 256                 , 1     , 1      ]
  swX   = [ 1                   , 0     , 0      ]
  xL    = [ Zero     * Kilometer, Zero  , Zero   ]
  xR    = [ 8.0e3_DP * Kilometer, Pi    , TwoPi  ]
  ZoomX = [ 1.020059256924852_DP, 1.0_DP, 1.0_DP ]

  ! --- DG ---

  nNodes = 2
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  ! --- Time Stepping ---

  nStagesSSPRK = 2
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  CFL = 0.5_DP ! Cockburn & Shu, (2001), JSC, 16, 173

  ! --- Slope Limiter ---

  UseSlopeLimiter           = .TRUE.
  SlopeLimiterMethod        = 'TVD'
  BetaTVD                   = 1.75e0_DP
  BetaTVB                   = 0.0d0
  SlopeTolerance            = 1.0e-6_DP
  UseCharacteristicLimiting = .FALSE.
  UseTroubledCellIndicator  = .FALSE.
  LimiterThresholdParameter = 1.5e-2_DP
  UseConservativeCorrection = .TRUE.

  ! --- Positivity Limiter ---

  UsePositivityLimiter = .TRUE.

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

  ALLOCATE( SourceTerms_Poseidon(1:nDOFX,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3),1:6) )

  CALL InitializeGravitySolver_CFA_Poseidon

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = TRIM( EosTableName ) )

  CALL InitializeSlopeLimiter_Euler_Relativistic_TABLE &
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

  CALL InitializePositivityLimiter_Euler_Relativistic_TABLE &
         ( UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           Verbose_Option = .TRUE., &
           Min_1_Option = ( One + EPSILON(One) ) * MinD, &
           Min_2_Option = ( One + EPSILON(One) ) * MinT, &
           Min_3_Option = ( One + EPSILON(One) ) * MinY, &
           Max_1_Option = ( One - EPSILON(One) ) * MaxD, &
           Max_2_Option = ( One - EPSILON(One) ) * MaxT, &
           Max_3_Option = ( One - EPSILON(One) ) * MaxY )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )
  WRITE(*,*)
  WRITE(*,'(A6,A,ES11.3E3)') '', 'CFL: ', CFL

  CALL InitializeFields_Relativistic &
         ( ProgenitorFileName_Option = TRIM( ProgenitorFileName ) )

  IF( RestartFileNumber .LT. 0 )THEN

    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

    CALL ComputeSourceTerms_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, SourceTerms_Poseidon )

    CALL SolveGravity_CFA_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, SourceTerms_Poseidon )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadFF_Option = .TRUE., ReadGF_Option = .TRUE. )

  END IF

  iCycleD = 10
!!$  iCycleW = 1; dt_wrt = -1.0d0
  dt_wrt = 0.1_DP * Millisecond; iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'
  WRITE(*,*)

  t_wrt  = t + dt_wrt
  t_wrt2 = t + dt_wrt
  wrt    = .FALSE.

  CALL InitializeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

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

    IF( t + dt .GT. t_wrt2 )THEN

      t_wrt2 = t_wrt2 + dt_wrt
      WriteSourceTerms2 = .TRUE.
      Time = t

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, &
             ComputeIncrement_Euler_DG_Explicit, &
             SolveGravity_CFA_Poseidon )

    WriteSourceTerms2 = .FALSE.

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

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

      CALL ComputeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

      wrt = .FALSE.

      CALL TimersStop_Euler( Timer_Euler_InputOutput )

    END IF

  END DO

  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,ES13.6E3,A)') 'Total evolution time: ', Timer_Evolution, ' s'
  WRITE(*,*)

  CALL TimersStart_Euler( Timer_Euler_Finalize )

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL ComputeSourceTerms_Poseidon &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, SourceTerms_Poseidon )

  CALL SolveGravity_CFA_Poseidon &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, SourceTerms_Poseidon )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL ComputeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

  CALL FinalizeTally_Euler_Relativistic

  CALL FinalizePositivityLimiter_Euler_Relativistic_TABLE

  CALL FinalizeSlopeLimiter_Euler_Relativistic_TABLE

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  DEALLOCATE( SourceTerms_Poseidon )

  CALL FinalizeGravitySolver_CFA_Poseidon

  CALL FinalizeEquationOfState

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
