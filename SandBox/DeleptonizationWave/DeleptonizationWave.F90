PROGRAM DeleptonizationWave

  USE KindModule, ONLY: &
    DP, SqrtTiny, Third, Zero
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
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
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, iPF_D, &
    uAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities
  USE TimeSteppingModule_Flash, ONLY: &
    Update_IMEX_PDARS
  USE InitializationModule, ONLY: &
    InitializeFields_DeleptonizationWave
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment, &
    TallyPositivityLimiter_TwoMoment

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL  :: wrt
  INTEGER  :: iCycle, iCycleD
  INTEGER  :: nE, nX(3), nNodes, nSpecies
  REAL(DP) :: t, dt, t_end, dt_wrt, t_wrt
  REAL(DP) :: eL, eR, ZoomE
  REAL(DP) :: xL(3), xR(3)

  nNodes   = 2
  nSpecies = 2

  nX = [ 12, 12, 12 ]
  xL = [ - 0.0d2, - 0.0d2, - 0.0d2 ] * Kilometer
  xR = [ + 1.0d2, + 1.0d2, + 1.0d2 ] * Kilometer

  nE = 16
  eL = 0.0d0 * MeV
  eR = 3.0d2 * MeV
  ZoomE = 1.183081754893913_DP

  t       = 0.0_DP
  t_end   = 1.0d-2 * Millisecond
  dt_wrt  = 1.0d-2 * Millisecond
  iCycleD = 1

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'DeleptonizationWave', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 01, 01, 01 ], &
           bcX_Option &
             = [ 32, 32, 32 ], &
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

  ! --- Initialize Timers ---

  CALL InitializeTimers

  CALL TimersStart( Timer_Initialize )

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

  ! --- Create Neutrino Opacities ---

  CALL CreateNeutrinoOpacities( nZ, nNodesZ, nSpecies )

  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = 0.0d0 + SqrtTiny, &
           Max_1_Option = 1.0d0 - EPSILON(1.0d0), &
           Min_2_Option = 0.0d0 + SqrtTiny, &
           UsePositivityLimiter_Option &
             = .TRUE., &
           UsePositivityLimiterTally_Option &
             = .TRUE. )

  ! --- Set Initial Condition ---

  CALL InitializeFields_DeleptonizationWave

  ! --- Write Initial Condition Before Limiter ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! Reset these timers
  Timer_PositivityLimiter = Zero
  Timer_PL_In             = Zero
  Timer_PL_Points         = Zero
  Timer_PL_CellAverage    = Zero
  Timer_PL_Theta_1        = Zero
  Timer_PL_Theta_2        = Zero
  Timer_PL_Out            = Zero

  CALL TallyPositivityLimiter_TwoMoment( 0.0_DP )

  CALL ComputeNeutrinoOpacities &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
           uPF(:,:,:,:,iPF_D), &
           uAF(:,:,:,:,iAF_T), &
           uAF(:,:,:,:,iAF_Ye) )

  ! --- Write Initial Condition After Limiter ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

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
             Explicit_Option = .TRUE., &
             Implicit_Option = .TRUE., &
             SingleStage_Option = .FALSE., &
             CallFromThornado_Option = .TRUE. )

    t = t + dt

    CALL TimersStop( Timer_Evolve )

    CALL TallyPositivityLimiter_TwoMoment( t )

    IF( wrt )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( uCF, uCR )
#endif

      CALL TimersStart( Timer_InputOutput )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE., &
               WriteOP_Option = .TRUE. )

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

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizeTimers

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  CALL DestroyNeutrinoOpacities

  CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

END PROGRAM DeleptonizationWave
