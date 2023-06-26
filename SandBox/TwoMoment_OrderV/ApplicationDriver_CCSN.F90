PROGRAM ApplicationDriver_CCSN

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
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, &
    uPF, uAF, uDF
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic, &
    ComputeTimeStep_Euler_NonRelativistic
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, uAR, uGR
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment, &
    ComputeTimeStep_TwoMoment
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    InitializeNeutrinoMatterSolverParameters
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    SolveGravity_Newtonian_Poseidon
  USE TwoMoment_TimeSteppingModule, ONLY: &
    Update_IMEX_RK
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE InitializationModule_CCSN, ONLY: &
    InitializeFields
  USE TwoMoment_TallyModule, ONLY: &
    ComputeTally

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: EosTableName
  CHARACTER(64) :: OpacityTableName_AbEm
  CHARACTER(64) :: OpacityTableName_Iso
  CHARACTER(64) :: OpacityTableName_NES
  CHARACTER(64) :: OpacityTableName_Pair
  CHARACTER(64) :: OpacityTableName_Brem
  CHARACTER(64) :: ProgenitorFileName
  LOGICAL       :: wrt
  LOGICAL       :: EvolveEuler
  LOGICAL       :: EvolveTwoMoment
  LOGICAL       :: UseSlopeLimiter_Euler
  LOGICAL       :: UseSlopeLimiter_TwoMoment
  LOGICAL       :: UsePositivityLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_TwoMoment
  LOGICAL       :: UseEnergyLimiter_TwoMoment
  LOGICAL       :: RampTimeStep
  LOGICAL       :: Include_NES
  LOGICAL       :: Include_Pair
  LOGICAL       :: Include_Brem
  LOGICAL       :: Include_LinCorr
  INTEGER       :: RestartFileNumber
  INTEGER       :: nNodes, nSpecies
  INTEGER       :: nX(3), bcX(3), nE, bcE
  INTEGER       :: iCycle, iCycleD, nEquidistantX
  INTEGER       :: M_outer, MaxIter_outer
  INTEGER       :: M_inner, MaxIter_inner
  REAL(DP)      :: xL(3), xR(3), eL, eR
  REAL(DP)      :: dEquidistantX, zoomE, CFL, RampFactor
  REAL(DP)      :: t, t_wrt, t_end, dt, dt_wrt
  REAL(DP)      :: dt_Fluid, dt_Neutrinos
  REAL(DP)      :: dt_Initial, dt_Ramp
  REAL(DP)      :: Rtol_outer, Rtol_inner
  REAL(DP)      :: wMatterRHS(5)

  ProgramName = 'CCSN'

  CoordinateSystem = 'SPHERICAL'

  EosTableName          = 'wl-EOS-SFHo-15-25-50.h5'
  OpacityTableName_AbEm = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
  OpacityTableName_Iso  = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
  OpacityTableName_NES  = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
  OpacityTableName_Pair = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'
  OpacityTableName_Brem = 'wl-Op-SFHo-15-25-50-E40-HR98-Brem.h5'

  ProgenitorFileName = 'WH07_15M_Sun.h5'

  RestartFileNumber = - 1

  ! --- Evolution Parameters ---

  t_end   = 3.50d+2 * Millisecond
  dt_wrt  = 1.00d-0 * Millisecond
  wrt     = .FALSE.
  iCycleD = 10

  ! --- Position Space Grid Parameters ---

  nX  = [ 512, 1, 1 ]
  xL  = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ]
  xR  = [ 8.0d3 * Kilometer, Pi   , TwoPi ]
  bcX = [ 30, 0, 0 ]

  nEquidistantX = 100
  dEquidistantX = 0.75d0 * Kilometer

  ! --- Energy Space Grid Parameters ---

  nE    = 16
  eL    = 0.0d0 * MeV
  eR    = 3.0d2 * MeV
  bcE   = 10
  zoomE = 1.266038160710160_DP

  ! --- Time Step Control ---

  IF( RestartFileNumber .LT. 0 )THEN
    RampTimeStep = .TRUE.    
  ELSE
    RampTimeStep = .FALSE.
  END IF
  RampFactor   = 1.1_DP
  dt_Initial   = 1.0d-6 * Millisecond

  ! --- Solvers Parameters ---

  nNodes   = 2
  nSpecies = 2
  CFL      = 0.5_DP / ( Two * DBLE( nNodes - 1 ) + One )

  EvolveEuler                    = .TRUE.
  UseSlopeLimiter_Euler          = .TRUE.
  UsePositivityLimiter_Euler     = .TRUE.

  EvolveTwoMoment                = .TRUE.
  UseSlopeLimiter_TwoMoment      = .FALSE.
  UsePositivityLimiter_TwoMoment = .TRUE.
  UseEnergyLimiter_TwoMoment     = .TRUE.

  IF( EvolveEuler .AND. .NOT. EvolveTwoMoment )THEN
    TimeSteppingScheme = 'SSPRK2'
  ELSE
    TimeSteppingScheme = 'IMEX_PDARS'
  END IF

  ! --- Neutrino-Matter Solvers Parameters ---

  M_outer         = 3
  MaxIter_outer   = 100
  Rtol_outer      = 1.0d-8
  M_inner         = 2
  MaxIter_inner   = 100
  Rtol_inner      = 1.0d-8
  Include_NES     = .TRUE.
  Include_Pair    = .TRUE.
  Include_Brem    = .TRUE.
  Include_LinCorr = .FALSE.
  wMatterRHS      = [ One, One, One, One, One ]

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
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

    ! --- Solve for Gravitational Potential ---

    CALL SolveGravity_Newtonian_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF(:,:,:,:,iCF_D) )
    
    ! --- Write Initial Condition ---

    CALL ComputeFromConserved_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    CALL ComputeFromConserved_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

    CALL WriteFieldsHDF &
           ( Time = t, &
             WriteGF_Option = .TRUE., &
             WriteFF_Option = .TRUE., &
             WriteRF_Option = .TRUE. )

    CALL ComputeTally &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR, &
             SetInitialValues_Option = .TRUE. )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadGF_Option = .TRUE., &
             ReadFF_Option = .TRUE., &
             ReadRF_Option = .TRUE. )

    ! --- Solve for Gravitational Potential (To Fill Geometry Ghost Cells) ---

    CALL SolveGravity_Newtonian_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF(:,:,:,:,iCF_D) )

    CALL ComputeTally &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR, &
             SetInitialValues_Option = .TRUE. )

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

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCF, uCR, &
             ComputeIncrement_TwoMoment_Implicit, &
             SolveGravity_Newtonian_Poseidon )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL ComputeTally &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

      wrt = .FALSE.

    END IF

  END DO
  
  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ComputeTally &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

  ! --- Auxiliary Finalization ---

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule, ONLY: &
      InitializeTimers
    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
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
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      InitializeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule, ONLY: &
      InitializeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
      InitializeGravitySolver_Newtonian_Poseidon
    USE TwoMoment_TallyModule, ONLY: &
      InitializeTally
    USE TwoMoment_TimeSteppingModule, ONLY: &
      Initialize_IMEX_RK

    CALL InitializeTimers

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
               = TRIM( OpacityTableName_Iso ) , &
             OpacityTableName_NES_Option &
               = TRIM( OpacityTableName_NES ) , &
             OpacityTableName_Pair_Option &
               = TRIM( OpacityTableName_Pair ), &
             OpacityTableName_Brem_Option &
               = TRIM( OpacityTableName_Brem ), &
             EquationOfStateTableName_Option &
               = TRIM( EosTableName ), &
             Verbose_Option = .TRUE. )

    ! --- Initialize Moment Closure ---

    CALL InitializeClosure_TwoMoment

    ! --- Initialize Slope Limiter (Euler) ---

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
           ( BetaTVD_Option &
               = 1.75_DP, &
             SlopeTolerance_Option &
               = 1.0d-6, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_Euler, &
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
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_D, &
             Min_2_Option &
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_T, &
             Min_3_Option &
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_Y, &
             Max_1_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_D, &
             Max_2_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_T, &
             Max_3_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_Y )

    ! --- Initialize Troubled Cell Indicator (Two-Moment) ---

    CALL InitializeTroubledCellIndicator_TwoMoment &
           ( UseTroubledCellIndicator_Option &
               = .FALSE., &
             C_TCI_Option &
               = Zero, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Slope Limiter (Two-Moment) ---

    CALL InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Positivity Limiter (Two-Moment) ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = SqrtTiny, &
             Min_2_Option &
               = SqrtTiny, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter_TwoMoment, &
             UseEnergyLimiter_Option &
               = UseEnergyLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Set Neutrino-Matter Solver Parameters ---

    CALL InitializeNeutrinoMatterSolverParameters &
           ( M_outer_Option &
               = M_outer, &
             M_inner_Option &
               = M_inner, &
             MaxIter_outer_Option &
               = MaxIter_outer, &
             MaxIter_inner_Option &
               = MaxIter_inner, &
             Rtol_inner_Option &
               = Rtol_inner, &
             Rtol_outer_Option &
               = Rtol_outer, &
             Include_NES_Option &
               = Include_NES, &
             Include_Pair_Option &
               = Include_Pair, &
             Include_Brem_Option &
               = Include_Brem, &
             Include_LinCorr_Option &
               = Include_LinCorr, &
             wMatrRHS_Option &
               = wMatterRHS, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Gravity Solver ---

    CALL InitializeGravitySolver_Newtonian_Poseidon

    ! --- Initialize Tally ---

    CALL InitializeTally

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK &
           ( TRIM( TimeSteppingScheme ), &
             EvolveEuler_Option = EvolveEuler, &
             EvolveTwoMoment_Option = EvolveTwoMoment )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE TwoMoment_TimeSteppingModule, ONLY: &
      Finalize_IMEX_RK
    USE TwoMoment_TallyModule, ONLY: &
      FinalizeTally
    USE EquationOfStateModule_TABLE, ONLY: &
      FinalizeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      FinalizeOpacities_TABLE
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      FinalizeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule, ONLY: &
      FinalizeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      FinalizePositivityLimiter_TwoMoment
    USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
      FinalizeGravitySolver_Newtonian_Poseidon
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
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE TwoMoment_TimersModule, ONLY: &
      FinalizeTimers

    CALL Finalize_IMEX_RK

    CALL FinalizeTally

    CALL FinalizeEquationOfState_TABLE

    CALL FinalizeOpacities_TABLE

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeTroubledCellIndicator_TwoMoment

    CALL FinalizeSlopeLimiter_TwoMoment

    CALL FinalizePositivityLimiter_TwoMoment

    CALL FinalizeGravitySolver_Newtonian_Poseidon

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver_CCSN
