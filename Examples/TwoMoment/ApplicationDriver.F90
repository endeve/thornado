PROGRAM DeleptonizationWave

  USE KindModule, ONLY: &
    DP, SqrtTiny, Third, Zero, One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV, &
    SpeedOfLight
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_Initialize, &
    Timer_InputOutput, &
    Timer_Evolve, &
    Timer_PositivityLimiter, &
    Timer_PL_In, &
    Timer_PL_Points, &
    Timer_PL_CellAverage, &
    Timer_PL_Theta_1, &
    Timer_PL_Theta_2, &
    Timer_PL_Out
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE RadiationFieldsModule, ONLY: &
    uCR  
  USE TimeSteppingModule_Flash, ONLY: &
    Update_IMEX_PDARS
  USE InitializationModule, ONLY: &
    InitializeFields, &
    ComputeError
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  LOGICAL       :: wrt
  LOGICAL       :: Explicit
  LOGICAL       :: Implicit
  LOGICAL       :: SingleStage
  INTEGER       :: iCycle, iCycleD
  INTEGER       :: nE, nX(3), bcX(3), nNodes, nSpecies
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt
  REAL(DP)      :: eL, eR, ZoomE
  REAL(DP)      :: xL(3), xR(3)

!  ProgramName = 'SineWaveStreaming' ! --- Test Explixit Part ---
  ProgramName = 'Relaxation'        ! --- Test Implicit Part ---

  SELECT CASE( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      nNodes   = 2
      nSpecies = 2

      nX  = [ 16, 2, 2 ]
      xL  = [ - 0.0d2, - 0.0d2, - 0.0d2 ] * Kilometer
      xR  = [ + 1.0d2, + 1.0d2, + 1.0d2 ] * Kilometer
      bcX = [ 1, 1, 1 ]

      nE = 16
      eL = 0.0d0 * MeV
      eR = 3.0d2 * MeV
      ZoomE = One

      t       = 0.0_DP
      t_end   = ( xR(1) - xL(1) ) / SpeedOfLight
      dt_wrt  = 0.5_DP * t_end
      iCycleD = 1

      Explicit    = .TRUE.
      Implicit    = .FALSE.
      SingleStage = .FALSE.

    CASE( 'Relaxation' )

      nNodes   = 2
      nSpecies = 2

      nX  = [ 2, 2, 2 ]
      xL  = [         0.0d0,         0.0d0,         0.0d0 ] * Kilometer
      xR  = [ REAL( nX(1) ), REAL( nX(2) ), REAL( nX(3) ) ] * Kilometer
      bcX = [ 0, 0, 0 ]
       
      nE = 16
      eL = 0.0d0 * MeV
      eR = 3.0d2 * MeV
      ZoomE = 1.183081754893913_DP

      t       = 0.0_DP
      t_end   = 1.0_DP * Millisecond
      dt_wrt  = 0.1_DP * t_end
      iCycleD = 1

      Explicit    = .FALSE.
      Implicit    = .TRUE.
      SingleStage = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'Invalid ProgramName'
      WRITE(*,*)
      STOP

  END SELECT

  ! --- Initialize Timers ---

  CALL InitializeTimers

  CALL TimersStart( Timer_Initialize )

  ! --- Auxiliary Initialization ---

  CALL InitializeDriver

  ! --- Set Initial Condition ---

  CALL InitializeFields

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! --- Reset These Timers ---

  Timer_PositivityLimiter = Zero
  Timer_PL_In             = Zero
  Timer_PL_Points         = Zero
  Timer_PL_CellAverage    = Zero
  Timer_PL_Theta_1        = Zero
  Timer_PL_Theta_2        = Zero
  Timer_PL_Out            = Zero

  ! --- Write Initial Condition ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  CALL TimersStop( Timer_Initialize )

  ! --- Evolve ---

  t_wrt   = dt_wrt
  wrt     = .FALSE.

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, &
    ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    dt = Third * MINVAL( (xR-xL) / DBLE( nX ) ) &
           / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t / Millisecond, &
          '', 'dt = ', dt / Millisecond

    END IF

    CALL TimersStart( Timer_Evolve )

    CALL Update_IMEX_PDARS &
           ( dt, uCF, uCR, &
             Explicit_Option &
               = Explicit, &
             Implicit_Option &
               = Implicit, &
             SingleStage_Option &
               = SingleStage, &
             CallFromThornado_Option &
               = .TRUE. )

    t = t + dt

    CALL TimersStop( Timer_Evolve )

    IF( wrt )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( uCF, uCR )
#endif

      IF( TRIM( ProgramName ) == 'Relaxation' )THEN

        CALL ComputeFromConserved_Euler_NonRelativistic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      END IF

      CALL TimersStart( Timer_InputOutput )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL TimersStop( Timer_InputOutput )

      wrt = .FALSE.

    END IF

  END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( uCF, uCR )
#endif

  ! --- Write Final Solution ---

  IF( TRIM( ProgramName ) == 'Relaxation' )THEN

    CALL ComputeFromConserved_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  END IF

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  CALL ComputeError( t )

  ! --- Finalize ---

  CALL FinalizeDriver

  CALL FinalizeTimers

CONTAINS


  SUBROUTINE InitializeDriver

    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
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
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange
    USE TwoMoment_ClosureModule, ONLY: &
      InitializeClosure_TwoMoment
    USE EquationOfStateModule_TABLE, ONLY: &
      InitializeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      InitializeOpacities_TABLE
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      InitializePositivityLimiter_TwoMoment

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 01, 01, 01 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             nE_Option &
               = nE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             ZoomE_Option &
               = ZoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = 'CARTESIAN', &
             ActivateUnits_Option &
               = .TRUE., &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

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

    ! --- Initialize Equation of State ---

    CALL InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = 'EquationOfStateTable.h5', &
             Verbose_Option = .TRUE. )

    IF( Implicit )THEN

      ! --- Initialize Opacities ---

      CALL InitializeOpacities_TABLE &
             ( OpacityTableName_EmAb_Option &
                 = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
               OpacityTableName_Iso_Option  &
                 = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &
               OpacityTableName_NES_Option &
                 = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &
               OpacityTableName_Pair_Option &
                 = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
               Verbose_Option = .TRUE. )

    END IF

    ! --- Initialize Positivity Limiter ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = 0.0d0 + SqrtTiny, &
             Min_2_Option &
               = 0.0d0 + SqrtTiny, &
             UsePositivityLimiter_Option &
               = .TRUE. )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE ReferenceElementModuleX, ONLY: &
      FinalizeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      FinalizeReferenceElementX_Lagrange
    USE ReferenceElementModuleE, ONLY: &
      FinalizeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      FinalizeReferenceElementE_Lagrange
    USE ReferenceElementModule, ONLY: &
      FinalizeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      FinalizeReferenceElement_Lagrange
    USE EquationOfStateModule_TABLE, ONLY: &
      FinalizeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      FinalizeOpacities_TABLE
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      FinalizePositivityLimiter_TwoMoment
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeEquationOfState_TABLE

    IF( Implicit )THEN

      CALL FinalizeOpacities_TABLE

    END IF

    CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

  END SUBROUTINE FinalizeDriver


END PROGRAM DeleptonizationWave
