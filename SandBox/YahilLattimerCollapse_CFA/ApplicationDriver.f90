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
    nDimsX, &
    nDOFX
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
    ReadFieldsHDF, &
    WriteSourceTermDiagnosticsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, &
    iPF_D, &
    uAF, &
    uDF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    InitializeGravitySolver_CFA_Poseidon, &
    FinalizeGravitySolver_CFA_Poseidon, &
    SolveGravity_CFA_Poseidon
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    Time
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK, &
    WriteSourceTerms2
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    Centimeter, &
    Gram, &
    Erg, &
    UnitsDisplay
  USE Euler_TallyModule_Relativistic, ONLY: &
    InitializeTally_Euler_Relativistic, &
    FinalizeTally_Euler_Relativistic, &
    ComputeTally_Euler_Relativistic
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler,  &
    FinalizeTimers_Euler, &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize, &
    Timer_Euler_Finalize
  USE Poseidon_UtilitiesModule, ONLY: &
    ComputeSourceTerms_Poseidon

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: FileName
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
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, t_wrt2, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: ZoomX(3)

  LOGICAL :: Skip10 = .FALSE.
  LOGICAL :: Skip11 = .FALSE.
  LOGICAL :: Skip12 = .FALSE.
  LOGICAL :: Skip13 = .FALSE.
  LOGICAL :: Skip14 = .FALSE.

  LOGICAL  :: WriteGF = .TRUE., WriteFF = .TRUE.
  LOGICAL  :: ActivateUnits = .TRUE.

  REAL(DP) :: Timer_Evolution

  REAL(DP), ALLOCATABLE :: SourceTerms_Poseidon(:,:,:,:,:)

  REAL(DP) :: CentralDensity
  REAL(DP) :: CentralPressure
  REAL(DP) :: CoreRadius
  REAL(DP) :: CollapseTime
  REAL(DP) :: D0
  LOGICAL  :: ReadFromFile

  REAL(DP), ALLOCATABLE :: Sources(:,:,:,:,:)

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  ProgramName = 'YahilCollapse'

  ReadFromFile = .FALSE.

  FileName = 'YahilHomologousCollapse_Gm1.30_t1.500E+000ms.dat'

  RestartFileNumber = -1
  t                 = 0.0_DP

  CoordinateSystem = 'SPHERICAL'

  CentralDensity  = 7.0e9_DP  * ( Gram / Centimeter**3 )
  CentralPressure = 6.0e27_DP * ( Erg  / Centimeter**3 )
  CoreRadius      = 1.0e5_DP  * Kilometer
  CollapseTime    = 1.50e2_DP * Millisecond

  ! --- These values come from Table 2 in the Yahil paper ---
  Gamma = 1.30_DP
  D0    = 1.75_DP

  t_end = CollapseTime - 0.5_DP * Millisecond
  bcX = [ 30, 0, 0 ]

  nX    = [ 256                 , 1     , 1      ]
  swX   = [ 1                   , 0     , 0      ]
  xL    = [ Zero                , Zero  , Zero   ]
  xR    = [ CoreRadius          , Pi    , TwoPi  ]
  ZoomX = [ 1.033358317557642_DP, 1.0_DP, 1.0_DP ]

  ! --- DG ---

  nNodes = 2
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  ! --- Time Stepping ---

  nStagesSSPRK = 3
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  CFL = 0.5_DP ! Cockburn & Shu, (2001), JSC, 16, 173

  ! --- Slope Limiter ---

  UseSlopeLimiter           = .FALSE.
  SlopeLimiterMethod        = 'TVD'
  BetaTVD                   = 1.75e0_DP
  BetaTVB                   = 0.0d0
  SlopeTolerance            = 1.0e-6_DP
  UseCharacteristicLimiting = .TRUE.
  UseTroubledCellIndicator  = .TRUE.
  LimiterThresholdParameter = 1.5e-2_DP
  UseConservativeCorrection = .TRUE.

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

  ALLOCATE( SourceTerms_Poseidon(1:nDOFX,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3),1:6) )

  ALLOCATE( Sources(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:6) )
  Sources = 0.0_DP

  CALL InitializeGravitySolver_CFA_Poseidon

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

  CALL InitializeFields_Relativistic &
         ( ReadFromFile_Option    = ReadFromFile,    &
           D0_Option              = D0,              &
           CentralDensity_Option  = CentralDensity,  &
           CentralPressure_Option = CentralPressure, &
           CoreRadius_Option      = CoreRadius,      &
           CollapseTime_Option    = CollapseTime )

  IF( RestartFileNumber .LT. 0 )THEN

    CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

    CALL ComputeSourceTerms_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, SourceTerms_Poseidon )

    CALL SolveGravity_CFA_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, SourceTerms_Poseidon )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

    CALL WriteSourceTermDiagnosticsHDF( 0.0_DP, Sources )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadFF_Option = .TRUE., ReadGF_Option = .TRUE. )

  END IF

  iCycleD = 10
!!$  iCycleW = 1; dt_wrt = -1.0d0
!!$  dt_wrt = 1.0d-2 * ( t_end - t ); iCycleW = -1
  dt_wrt = 0.1_DP * Millisecond; iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'
  WRITE(*,*)

  t_wrt = t + dt_wrt
  t_wrt2 = t + dt_wrt
  wrt   = .FALSE.

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

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e15_DP * Gram / Centimeter**3 ) )THEN

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

      CALL WriteSourceTermDiagnosticsHDF( t, Sources )

      CALL ComputeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

      EXIT

   END IF

!!$   IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e14_DP * Gram / Centimeter**3 ) &
!!$         .AND. .NOT. Skip14 )THEN
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )
!!$
!!$      CALL ComputeTally_Euler_Relativistic &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )
!!$
!!$      Skip14 = .TRUE.
!!$
!!$    END IF
!!$
!!$   IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e13_DP * Gram / Centimeter**3 ) &
!!$         .AND. .NOT. Skip13 )THEN
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )
!!$
!!$      CALL ComputeTally_Euler_Relativistic &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )
!!$
!!$      Skip13 = .TRUE.
!!$
!!$    END IF
!!$
!!$   IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e12_DP * Gram / Centimeter**3 ) &
!!$         .AND. .NOT. Skip12 )THEN
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )
!!$
!!$      CALL ComputeTally_Euler_Relativistic &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )
!!$
!!$      Skip12 = .TRUE.
!!$
!!$    END IF
!!$
!!$   IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e11_DP * Gram / Centimeter**3 ) &
!!$         .AND. .NOT. Skip11 )THEN
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )
!!$
!!$      CALL ComputeTally_Euler_Relativistic &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )
!!$
!!$      Skip11 = .TRUE.
!!$
!!$    END IF
!!$
!!$   IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e10_DP * Gram / Centimeter**3 ) &
!!$         .AND. .NOT. Skip10 )THEN
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )
!!$
!!$      CALL ComputeTally_Euler_Relativistic &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )
!!$
!!$      Skip10 = .TRUE.
!!$
!!$    END IF

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

  CALL WriteSourceTermDiagnosticsHDF( t, Sources )

  CALL ComputeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

  CALL FinalizeTally_Euler_Relativistic

  CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  DEALLOCATE( Sources )

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
