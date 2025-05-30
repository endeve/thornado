PROGRAM ApplicationDriver_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Pi, TwoPi, SqrtTiny
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
    uCF, uPF, uAF, uDF
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, uAR, uGR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
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
  USE TwoMoment_TimeSteppingModule, ONLY: &
    Update_IMEX_RK
  USE InitializationModule_Neutrinos, ONLY: &
    InitializeFields, &
    ComputeError
  USE TwoMoment_TallyModule, ONLY: &
    ComputeTally

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: ProfileName
  CHARACTER(64) :: EosTableName
  CHARACTER(64) :: OpacityTableName_EmAb
  CHARACTER(64) :: OpacityTableName_Iso
  CHARACTER(64) :: OpacityTableName_NES
  CHARACTER(64) :: OpacityTableName_Pair
  CHARACTER(64) :: OpacityTableName_Brem
  LOGICAL       :: EvolveEuler
  LOGICAL       :: UseSlopeLimiter_Euler
  LOGICAL       :: UseSlopeLimiter_TwoMoment
  LOGICAL       :: UsePositivityLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_TwoMoment
  LOGICAL       :: UseEnergyLimiter_TwoMoment
  LOGICAL       :: PrescribedTimeStep
  LOGICAL       :: Include_NES
  LOGICAL       :: Include_Pair
  LOGICAL       :: Include_NuPair
  LOGICAL       :: Include_Brem
  LOGICAL       :: Include_LinCorr
  LOGICAL       :: FreezeOpacities
  INTEGER       :: RestartFileNumber
  INTEGER       :: nSpecies
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: nEquidistantX
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER       :: M_outer, MaxIter_outer
  INTEGER       :: M_inner, MaxIter_inner
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: dEquidistantX
  REAL(DP)      :: t, dt, dt_CFL, dt_0, dt_MAX, dt_RATE, t_end
  REAL(DP)      :: Rtol_outer, Rtol_inner
  REAL(DP)      :: wMatterRHS(5)
  REAL(DP)      :: DnuMax
  LOGICAL       :: Relaxation_restart_from_file

  ProgramName = 'Relaxation'
  Relaxation_restart_from_file = .FALSE. 
  !ProgramName = 'DeleptonizationWave1D'

  CoordinateSystem = 'CARTESIAN'

  EosTableName          = 'wl-EOS-SFHo-15-25-50.h5'
  OpacityTableName_EmAb = 'wl-Op-SFHo-15-25-50-E40-EmAb.h5'
  OpacityTableName_Iso  = 'wl-Op-SFHo-15-25-50-E40-Iso.h5'
  OpacityTableName_NES  = 'wl-Op-SFHo-15-25-50-E40-NES.h5'
  OpacityTableName_Pair = 'wl-Op-SFHo-15-25-50-E40-Pair.h5'
  OpacityTableName_Brem = 'wl-Op-SFHo-15-25-50-E40-Brem.h5'

  PrescribedTimeStep = .FALSE.

  RestartFileNumber = - 1

  M_outer         = 2
  MaxIter_outer   = 100
  Rtol_outer      = 1.0d-8
  M_inner         = 2
  MaxIter_inner   = 100
  Rtol_inner      = 1.0d-8
  Include_NES     = .TRUE.
  Include_Pair    = .TRUE.
  Include_NuPair  = .FALSE.
  Include_Brem    = .TRUE.
  Include_LinCorr = .FALSE.
  wMatterRHS      = [ One, One, One, One, One ]
  DnuMax          = HUGE( One )
  FreezeOpacities = .FALSE.

  nEquidistantX = 1
  dEquidistantX = 1.0_DP * Kilometer

  SELECT CASE( TRIM( ProgramName ) )

    CASE( 'Relaxation' )

      nSpecies = 6
      nNodes   = 2
      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 0
      ZoomE = 1.266038160710160_DP

      IF( Relaxation_restart_from_file ) THEN
        CoordinateSystem = 'SPHERICAL'

        nX    = [ 1, 1, 1 ]
        xL    = [ 0.0_DP            , 0.0_DP, 0.0_DP ]
        xR    = [ 8.0_DP * Kilometer, Pi    , TwoPi  ]
        bcX   = [ 0, 0, 0 ]
        ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        nEquidistantX = 0
        dEquidistantX = 1.0_DP * Kilometer

        nE    = 16
        eL    = 0.0d0 * MeV
        eR    = 3.0d2 * MeV
        bcE   = 10
        ZoomE = 1.266038160710160_DP

        !TimeSteppingScheme = 'IMEX_PDARS'
        TimeSteppingScheme = 'BackwardEuler'

        t_end = 1.0d0 * Millisecond

        PrescribedTimeStep = .TRUE. ! If .FALSE., explicit CFL will be used.
        dt_0               = 8.895d-3 * Millisecond
        dt_MAX             = 8.895d-3 * Millisecond
        dt_RATE            = 1.01_DP

        iCycleD = 1
        iCycleW = 100
        maxCycles = 1000000

        EvolveEuler                    = .FALSE.
        UseSlopeLimiter_Euler          = .FALSE.
        UseSlopeLimiter_TwoMoment      = .FALSE.
        UsePositivityLimiter_Euler     = .FALSE.
        UsePositivityLimiter_TwoMoment = .FALSE.
        UseEnergyLimiter_TwoMoment     = .FALSE.

        Include_NES     = .TRUE.
        Include_Pair    = .TRUE.
        Include_NuPair  = .FALSE.
        Include_Brem    = .TRUE.
        Include_LinCorr = .FALSE.
        FreezeOpacities = .FALSE. ! --- Keep opacities fixed during iterations?
        DnuMax          = One
        M_outer         = 2
        MaxIter_outer   = 100
        Rtol_outer      = 1.0d-8
        M_inner         = 2
        MaxIter_inner   = 100
        Rtol_inner      = 1.0d-8
        Include_NES     = .TRUE.
        Include_Pair    = .TRUE.
        Include_NuPair  = .FALSE.
        Include_Brem    = .TRUE.
        Include_LinCorr = .FALSE.
        wMatterRHS      = [ One, One, One, One, One ]
        DnuMax          = One
        FreezeOpacities = .FALSE.

      ELSE

        nX  = [ 1, 1, 1 ]
        xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
        xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
        bcX = [ 0, 0, 0 ]

        TimeSteppingScheme = 'BackwardEuler'

        t_end = 1.0d2 * Millisecond

        PrescribedTimeStep = .TRUE.
        dt_0               = 1.0d-4 * Millisecond
        dt_MAX             = 1.0d+1 * Millisecond
        dt_RATE            = 1.04_DP
        iCycleD            = 1
        iCycleW            = 1
        maxCycles          = 100000

        EvolveEuler                    = .FALSE.
        UseSlopeLimiter_Euler          = .FALSE.
        UseSlopeLimiter_TwoMoment      = .FALSE.
        UsePositivityLimiter_Euler     = .FALSE.
        UsePositivityLimiter_TwoMoment = .FALSE.
        UseEnergyLimiter_TwoMoment     = .FALSE.
      
      ENDIF

    CASE( 'DeleptonizationWave1D' )

      CoordinateSystem = 'SPHERICAL'

      nSpecies = 6
      nNodes   = 2

      nX    = [ 256, 1, 1 ]
      xL    = [ 0.0_DP           , 0.0_DP, 0.0_DP ]
      xR    = [ 5.12d2 * Kilometer, Pi    , TwoPi  ]
      bcX   = [ 31, 1, 1 ]
      ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      nEquidistantX = 0
      dEquidistantX = 1.0_DP * Kilometer

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end = 1.0d2 * Millisecond

      PrescribedTimeStep = .TRUE. ! If .FALSE., explicit CFL will be used.
      dt_0               = 1.0d-3 * Millisecond
      dt_MAX             = 1.0d-3 * Millisecond
      dt_RATE            = 1.01_DP

      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .TRUE.

      ProfileName = 'input_thornado_VX_100ms.dat'

      wMatterRHS = [ One, One, Zero, Zero, Zero ] ! --- Keep Velocity Fixed

      Include_NES     = .TRUE.
      Include_Pair    = .TRUE.
      Include_NuPair  = .TRUE.
      Include_Brem    = .TRUE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .FALSE. ! --- Keep opacities fixed during iterations?

    CASE( 'EquilibriumAdvection' )

      nSpecies = 2
      nNodes   = 2

      nX  = [ 8, 1, 1 ]
      xL  = [ - 5.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ + 5.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 1, 1, 1 ]

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end = 1.0d1 * Millisecond

      iCycleD = 1
      iCycleW = 100
      maxCycles = 100000

      EvolveEuler                    = .TRUE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .TRUE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  CALL InitializeDriver

  IF (Relaxation_restart_from_file) THEN
    CALL InitializeFields( ProfileName, 'Relaxation_input.dat' )
  ELSE
    CALL InitializeFields( ProfileName )
  ENDIF


  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  IF( RestartFileNumber .LT. 0 )THEN

    t = Zero

    ! --- Apply Slope Limiter to Initial Data ---

    CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplySlopeLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

    ! --- Apply Positivity Limiter to Initial Data ---

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

    CALL WriteFieldsHDF &
           ( Time = 0.0_DP, &
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

  END IF

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( PrescribedTimeStep )THEN

      dt = MIN( ( dt_RATE )**( iCycle - 1 ) * dt_0, dt_MAX )

    ELSE

      CALL ComputeTimeStep_TwoMoment &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, &
               0.3_DP/(Two*DBLE(nNodes-1)+One), dt_CFL )

      dt = dt_CFL

    END IF

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = '    ,  t / Millisecond, &
          '', 'dt = '   , dt / Millisecond

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

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

  CALL ComputeError( t )

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule, ONLY: &
      InitializeTimers
    USE TimersModule_Euler, ONLY: &
      InitializeTimers_Euler
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
    USE TwoMoment_TallyModule, ONLY: &
      InitializeTally
    USE TwoMoment_TimeSteppingModule, ONLY: &
      Initialize_IMEX_RK

    CALL InitializeTimers

    CALL InitializeTimers_Euler

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
             zoomX_Option &
               = zoomX, &
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

    IF( nEquidistantX > 1 )THEN

      CALL CreateMesh_Custom &
             ( MeshX(1), nX(1), nNodes, 1, xL(1), xR(1), &
               nEquidistantX, dEquidistantX, Verbose_Option = .TRUE. )

    END IF

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
             UseChemicalPotentialShift_Option &
               = .TRUE., &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Opacities ---

    CALL InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option &
               = TRIM( OpacityTableName_EmAb ), &
             OpacityTableName_Iso_Option  &
               = TRIM( OpacityTableName_Iso ), &
             OpacityTableName_NES_Option &
               = TRIM( OpacityTableName_NES ), &
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
             Include_NuPair_Option &
               = Include_NuPair, &
             Include_Brem_Option &
               = Include_Brem, &
             Include_LinCorr_Option &
               = Include_LinCorr, &
             wMatrRHS_Option &
               = wMatterRHS, &
             DnuMax_Option &
               = DnuMax, &
             FreezeOpacities_Option &
               = FreezeOpacities, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Tally ---

    CALL InitializeTally

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK &
           ( TRIM( TimeSteppingScheme ), EvolveEuler_Option = EvolveEuler )

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
    USE TimersModule_Euler, ONLY: &
      FinalizeTimers_Euler
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

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

    CALL FinalizeTimers_Euler

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver_Neutrinos
