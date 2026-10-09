PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Three, Half, &
    Pi, TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFE, nDOFX, nDOFZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE, nGE
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, uAR, uGR, nCR, nDR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE MeshModule, ONLY: &
    MeshX, NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    WeightsX1
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment, &
    ComputeTimeStep_TwoMoment, &
    ComputeTimeStep_TwoMoment_Realizability
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_Legacy &
      => ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    OffGridFlux_TwoMoment
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    nQ, nC
  USE TwoMoment_OpacityModule, ONLY: &
    SetOpacities
  USE TwoMoment_TimeSteppingModule, ONLY: &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields, &
    ComputeError
  USE TwoMoment_TallyModule, ONLY: &
    ComputeTally

  IMPLICIT NONE

  CHARACTER(2)  :: Direction
  CHARACTER(32) :: Spectrum
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(64) :: ProgramNameEnvironment
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseEnergyLimiter
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: UseRealizabilityTimeStep
  LOGICAL       :: RunRelaxationDiffusion
  INTEGER       :: nNodes
  INTEGER       :: nSpecies
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER       :: EnvironmentLength, EnvironmentStatus
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t, dt, t_end, dt_CFL, dt_0, dt_grw, V_0(3)
  REAL(DP)      :: D_0, Chi, Sigma, C_TCI, C_CFL
  REAL(DP)      :: LengthScale
  REAL(DP)      :: RelaxationTime, RelaxationAdvectionSpeed
  REAL(DP)      :: RelaxationWaveSpeedSquared, RelaxationAmplitude
  REAL(DP)      :: RelaxationBackground

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'SineWaveStreaming'

  ProgramNameEnvironment = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_PROGRAM_NAME', ProgramNameEnvironment, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 ) &
    ProgramName = TRIM( ProgramNameEnvironment )

  nSpecies = 1

  C_TCI = 1.0_DP
  UseTroubledCellIndicator = .FALSE.

  UseRealizabilityTimeStep = .FALSE.

  dt_0   = HUGE( One )
  dt_grw = One

  RunRelaxationDiffusion = .FALSE.
  RelaxationTime = 1.0d-3
  RelaxationAdvectionSpeed = Zero
  RelaxationWaveSpeedSquared = One / 3.0_DP
  RelaxationAmplitude = 0.1_DP
  RelaxationBackground = One

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'RelaxationDiffusionStatic', 'RelaxationDiffusionDynamic' )

      ! d_t D + d_x I = 0
      ! d_t I + a^2 d_x D = - ( I - v D ) / epsilon

      RunRelaxationDiffusion = .TRUE.

      nX  = [ 64, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nSpecies = 1
      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end = 1.0_DP
      iCycleD = 100
      iCycleW = 1000000
      maxCycles = 1000000

      IF( TRIM( ProgramName ) == 'RelaxationDiffusionDynamic' ) &
        RelaxationAdvectionSpeed = 0.2_DP

      V_0 = Zero
      D_0 = Zero
      Chi = Zero
      Sigma = Zero

      UseSlopeLimiter = .FALSE.
      UsePositivityLimiter = .FALSE.
      UseEnergyLimiter = .FALSE.

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 64, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 8
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nSpecies = 1
      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 1.0d-0
      iCycleD = 1
      iCycleW = 100
      maxCycles = 10000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      Direction = 'X'

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      C_TCI = 1.0_DP
      UseTroubledCellIndicator = .FALSE.
      UseSlopeLimiter          = .FALSE.
      UsePositivityLimiter     = .FALSE.
      UseEnergyLimiter         = .FALSE.

    CASE( 'SineWaveDiffusion' )

      nX  = [ 32, 1, 1 ]
      xL  = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ + 3.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      !t_end   = 7.2951_DP !sig small old
      !t_end   = 1459.02504445_DP !sig large old

      !t_end   = 60.0_DP !sig small
      t_end   = 1500.0_DP !sig large

      iCycleD = 1000
      iCycleW = 1000
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      !Sigma = 2.6666666_DP 
      Sigma = 533.33333333_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .FALSE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'SphericalDiffusion' )

      CoordinateSystem = 'SPHERICAL'

      nX  = [ 50, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP,     Pi,  TwoPi ]
      bcX = [ 31, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 5.0d+0
      iCycleD = 50
      iCycleW = 50
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 1.0d+2

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'IsotropicRadiation' )

      nX  = [ 16, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 0, 0 ]

      nE  = 10
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 2

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 1.0d1
      iCycleD = 1
      iCycleW = 10
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .FALSE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'StreamingDopplerShift' )

      Spectrum = 'Fermi-Dirac'

      Direction = 'X' ! --- (X,Y, or Z)

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        nX  = [ 100, 1, 1 ]
        xL  = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR  = [ 1.0d1, 1.0d0, 1.0d0 ]
        bcX = [ 12, 1, 1 ]

        V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

        nX  = [ 1, 32, 1 ]
        xL  = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR  = [ 1.0d0, 1.0d1, 1.0d0 ]
        bcX = [ 1, 12, 1 ]

        V_0 = [ 0.0_DP, 0.1_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

        nX  = [ 1, 1, 32 ]
        xL  = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR  = [ 1.0d0, 1.0d0, 1.0d1 ]
        bcX = [ 1, 1, 12 ]

        V_0 = [ 0.0_DP, 0.0_DP, 0.1_DP ]

      ELSE

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'StreamingDopplerShift.  Direction must be X, Y, or Z'
        WRITE(*,*)
        STOP

      END IF

      nE    = 32
      eL    = 0.0d0
      eR    = 5.0d1
      bcE   = 11
      zoomE = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 2.0d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

      UseRealizabilityTimeStep = .TRUE.

    CASE( 'TransparentTurbulence' )

      Direction = 'X' ! --- (X,Y, or Z)

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        nX  = [ 64, 1, 1 ]
        xL  = [ - 1.0_DP, 0.0_DP, 0.0_DP ]
        xR  = [ + 1.0_DP, 1.0_DP, 1.0_DP ]
        bcX = [ 12, 1, 1 ]

        V_0 = [ 0.01_DP, 0.0_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

        nX  = [ 2, 80, 1 ]
        xL  = [ 0.0_DP, - 1.0_DP, 0.0_DP ]
        xR  = [ 1.0_DP, + 1.0_DP, 1.0_DP ]
        bcX = [ 1, 12, 1 ]

        V_0 = [ 0.0_DP, 0.1_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

        nX  = [ 2, 2, 80 ]
        xL  = [ 0.0_DP, 0.0_DP, - 1.0_DP ]
        xR  = [ 1.0_DP, 1.0_DP, + 1.0_DP ]
        bcX = [ 1, 1, 12 ]

        V_0 = [ 0.0_DP, 0.0_DP, 0.1_DP ]

      ELSE

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'TransparentTurbulence.  Direction must be X, Y, or Z'
        WRITE(*,*)
        STOP

      END IF

      nE  = 16
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 10

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 5.0d0
      iCycleD = 1
      iCycleW = 250
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'TransparentShock' )

      Direction = 'X' ! --- (X,Y, or Z)

      LengthScale = 1.0d-2 ! --- Shock Width

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        nX  = [ 80, 1, 1 ]
        xL  = [ 0.0d0, 0.0_DP, 0.0_DP ]
        xR  = [ 2.0d0, 1.0_DP, 1.0_DP ]
        bcX = [ 12, 1, 1 ]

        V_0 = [ - 0.1_DP, 0.0_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

        nX  = [ 2, 80, 1 ]
        xL  = [ 0.0d0, 0.0_DP, 0.0_DP ]
        xR  = [ 1.0d0, 2.0_DP, 1.0_DP ]
        bcX = [ 1, 12, 1 ]

        V_0 = [ 0.0_DP, - 0.1_DP, 0.0_DP ]

      ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

        nX  = [ 2, 2, 80 ]
        xL  = [ 0.0d0, 0.0_DP, 0.0_DP ]
        xR  = [ 1.0d0, 1.0_DP, 2.0_DP ]
        bcX = [ 1, 1, 12 ]

        V_0 = [ 0.0_DP, 0.0_DP, - 0.1_DP ]

      ELSE

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'TransparentShock.  Direction must be X, Y, or Z'
        WRITE(*,*)
        STOP

      END IF

      nE  = 32
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 11

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 3.0d0
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

      UseRealizabilityTimeStep = .TRUE.

    CASE( 'TransparentShock2D' )

      LengthScale = 3.0d-2 ! --- Shock Width

      nX  = [ 80, 10, 1 ]
      xL  = [ 0.0d0, 0.00_DP, 0.0_DP ]
      xR  = [ 2.0d0, 0.25_DP, 1.0_DP ]
      bcX = [ 12, 2, 0 ]

      nE  = 32
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 11

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 3.0d0
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

      UseRealizabilityTimeStep = .TRUE.

    CASE( 'TransparentVortex' )

      Direction = 'X' ! --- (X or Y)

      nX  = [ 48, 48, 1 ]
      xL  = [ - 5.0_DP, - 5.0_DP, - 0.5_DP ]
      xR  = [ + 5.0_DP, + 5.0_DP, + 0.5_DP ]

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        bcX = [ 12, 3, 1 ]

      ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

        bcX = [ 3, 12, 1 ]

      ELSE

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'TransparentVortex.  Direction must be X or Y'
        WRITE(*,*)
        STOP

      END IF

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      nE  = 32
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 11

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 2.0d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

      UseRealizabilityTimeStep = .TRUE.

    CASE( 'TransparentVortex_Spherical' )

      CoordinateSystem = 'SPHERICAL'

      nX  = [ 48, 48, 1 ]
      xL  = [  2.0_DP, 0.0_DP     , 0.0_DP ]
      xR  = [ 10.0_DP, 0.5_DP * Pi,  TwoPi ]
      bcX = [ 12, 3, 1 ]

      V_0 = [ SQRT(10.0_DP)*0.001_DP, 0.0_DP, 0.0_DP ]

      nE  = 32
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 11

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 1.5d+1
      iCycleD = 1
      iCycleW = 125
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

    CASE( 'RadiatingSphere' )

      CoordinateSystem = 'SPHERICAL'

      nX    = [ 200, 1, 1 ]
      xL    = [ 1.0d1, Zero,  Zero ]
      xR    = [ 1.0d4,   Pi, TwoPi ]
      bcX   = [ 12, 1, 1 ]

      IF     ( nX(1) == 200 )THEN
         ZoomX = [ 1.028717467149282_DP, One, One ] ! dX1(1)=1 with nX(1)=200
      ELSE IF( nX(1) == 400 )THEN
         ZoomX = [ 1.012074699662899_DP, One, One ] ! dX1(1)=1 with nX(1)=400
      ELSE
         ZoomX = One
      END IF

      nE    = 32
      eL    = 0.0d0
      eR    = 3.0d2
      bcE   = 11
      ZoomE = 1.119237083677839_DP

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end = 2.0d+4
      iCycleD = 1
      iCycleW = 2000
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      C_TCI = 0.1_DP
      UseTroubledCellIndicator = .FALSE.
      UseSlopeLimiter          = .FALSE.
      UsePositivityLimiter     = .TRUE.
      UseEnergyLimiter         = .TRUE.

    CASE( 'GaussianDiffusion1D' )

      nX  = [ 96, 1, 1 ]
      xL  = [ 0.0_DP, - 0.5_DP, - 0.5_DP ]
      xR  = [ 3.0_DP, + 0.5_DP, + 0.5_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 1

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end     = 5.0_DP
      iCycleD   = 10
      iCycleW   = 1000
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 3200.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'GaussianDiffusion' )

      nX  = [ 96, 64, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, - 0.5_DP ]
      xR  = [ 3.0_DP, 2.0_DP, + 0.5_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 1

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 5.0_DP
      iCycleD = 10
      iCycleW = 10
      maxCycles = 1000000

      V_0 = [ 0.0_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 3200.0_DP
      !Sigma = 1.0d+2

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'ExpandingAtmosphere' )

      CoordinateSystem = 'SPHERICAL'

      nX  = [ 200, 1, 1 ]
      xL  = [  0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 15.0_DP,     Pi,  TwoPi ]
      bcX = [ 31, 1, 1 ]

      nE  = 50
      eL  = 00.0_DP
      eR  = 12.0_DP
      bcE = 11

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 1.5d+1
      iCycleD = 100
      iCycleW = 100
      maxCycles = 1000000

      V_0 = [ 0.3_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

    CASE( 'HomogeneousSphere1D' )

      CoordinateSystem = 'SPHERICAL'

      nX  = [ 100, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 5.0_DP,     Pi,  TwoPi ]
      bcX = [ 30, 1, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 1

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 2.0d+1
      dt_0    = 1.0d-5
      dt_grw  = 1.01_DP
      iCycleD = 1
      iCycleW = 500
      maxCycles = 1000000

      V_0 = [ 0.0_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.8d0
      Chi   = 4.0d0
      Sigma = 0.0d0

      C_TCI = 0.1_DP
      UseTroubledCellIndicator = .TRUE.
      UseSlopeLimiter          = .TRUE.
      UsePositivityLimiter     = .TRUE.
      UseEnergyLimiter         = .FALSE.

    CASE( 'HomogeneousSphere2D' )

      nX  = [ 48, 48, 1 ]
      xL  = [ - 3.0_DP, - 3.0_DP, - 1.0_DP ]
      xR  = [ + 3.0_DP, + 3.0_DP, + 1.0_DP ]
      bcX = [ 2, 2, 1 ]

      nE  = 12
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 10

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 1.0d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.8_DP
      Chi   = 4.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'ShadowCasting2D_Cylindrical' )

      CoordinateSystem = 'CYLINDRICAL'

      nX  = [ 168, 140, 1 ]
      xL  = [ 00.0_DP, - 5.0_DP, 0.0_DP ]
      xR  = [ 12.0_DP, + 5.0_DP,  TwoPi ]
      bcX = [ 30, 2, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 1

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 1.5d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      V_0 = [ 0.0_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      C_TCI = 0.1_DP
      UseTroubledCellIndicator = .FALSE.
      UseSlopeLimiter          = .FALSE.
      UsePositivityLimiter     = .TRUE.
      UseEnergyLimiter         = .FALSE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  ! --- Auxiliary Initialization ---

  CALL InitializeDriver

  CALL SetOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D_0, Chi, Sigma )

  ! --- Set Initial Condition ---

  IF( RunRelaxationDiffusion )THEN
    CALL InitializeRelaxationDiffusion
  ELSE
    CALL InitializeFields( V_0, LengthScale, Direction, Spectrum )
  END IF

  t = 0.0_DP

  ! --- Apply Slope Limiter to Initial Data ---

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Apply Positivity Limiter to Initial Data ---

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Write Initial Condition ---

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ComputeTally &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR, &
           SetInitialValues_Option = .TRUE. )

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t, ' to t = ', t_end
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( UseRealizabilityTimeStep )THEN

      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uCF )

      CALL ComputeTimeStep_TwoMoment_Realizability &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, One, dt_CFL )

    ELSE

      C_CFL  = 0.3_DP / ( Two * DBLE(nNodes-1) + One )

      CALL ComputeTimeStep_TwoMoment &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, C_CFL, dt_CFL )

    END IF

    IF( dt_grw > One )THEN

      IF( iCycle == 1 )THEN

        dt = MIN( dt_CFL, dt_0 )

      ELSE

        dt = MIN( dt_CFL, dt_grw * dt )

      END IF

    ELSE

      dt = dt_CFL

    END IF

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    IF( RunRelaxationDiffusion )THEN
      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_RelaxationDiffusion_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_RelaxationDiffusion_Explicit )
    ELSE
      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_TwoMoment_Implicit )
    END IF

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL ComputeTally( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

    END IF

  END DO

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ComputeTally( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

  IF( RunRelaxationDiffusion )THEN
    CALL ComputeError_RelaxationDiffusion( t )
  ELSE
    CALL ComputeError( t )
  END IF

  ! --- Auxiliary Finalization ---

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeRelaxationDiffusion

    INTEGER :: iE, iNodeE, iNodeX, iNodeZ, iS, iX1
    REAL(DP) :: D, I, X

    IF( nDOFX /= 2 ) &
      ERROR STOP 'Relaxation-diffusion test requires two spatial nodes.'
    IF( iX_E0(2) /= iX_B0(2) .OR. iX_E0(3) /= iX_B0(3) ) &
      ERROR STOP 'Relaxation-diffusion test is one-dimensional.'

    uCF = Zero
    uCF(:,:,:,:,iCF_D ) = One
    uCF(:,:,:,:,iCF_E ) = One
    uCF(:,:,:,:,iCF_S1) = Zero
    uCF(:,:,:,:,iCF_S2) = Zero
    uCF(:,:,:,:,iCF_S3) = Zero
    uCF(:,:,:,:,iCF_Ne) = Zero
    uCR = Zero

    DO iS = 1, nSpecies
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      X = NodeCoordinate( MeshX(1), iX1, iNodeX )
      CALL RelaxationDiffusionExact( X, Zero, D, I )

      DO iE = iE_B0, iE_E0
      DO iNodeE = 1, nDOFE
        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
        uCR(iNodeZ,iE,iX1,iX_B0(2),iX_B0(3),iCR_N ,iS) = D
        uCR(iNodeZ,iE,iX1,iX_B0(2),iX_B0(3),iCR_G1,iS) = I
      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRelaxationDiffusion


  SUBROUTINE ComputeIncrement_RelaxationDiffusion_Explicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    INTEGER :: iE, iNodeE, iNodeZ1, iNodeZ2, iPrev, iNext, iS, iX1
    REAL(DP) :: a, h
    REAL(DP) :: DL, DR, IL, IR, DMean, IMean
    REAL(DP) :: FluxNLeft, FluxNRight, FluxILeft, FluxIRight
    REAL(DP) :: MeanRateN, MeanRateI, SlopeRateN, SlopeRateI
    REAL(DP) :: D1, D2, I1, I2

    IF( nDOFX /= 2 ) &
      ERROR STOP 'Relaxation-diffusion test requires two spatial nodes.'

    dU_R = Zero
    OffGridFlux_TwoMoment = Zero
    a = SQRT( RelaxationWaveSpeedSquared )

    DO iS = 1, nSpecies
    DO iE = iZ_B0_Input(1), iZ_E0_Input(1)
    DO iNodeE = 1, nDOFE
    DO iX1 = iZ_B0_Input(2), iZ_E0_Input(2)

      iPrev = iX1 - 1
      IF( iPrev < iZ_B0_Input(2) ) iPrev = iZ_E0_Input(2)
      iNext = iX1 + 1
      IF( iNext > iZ_E0_Input(2) ) iNext = iZ_B0_Input(2)

      iNodeZ1 = iNodeE
      iNodeZ2 = nDOFE + iNodeE
      h = MeshX(1) % Width(iX1)

      DL = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iPrev,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               U_R(iNodeZ2,iE,iPrev,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               +1 )
      DR = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               -1 )
      IL = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iPrev,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               U_R(iNodeZ2,iE,iPrev,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               +1 )
      IR = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               -1 )
      FluxNLeft = Half * ( IL + IR ) - Half * a * ( DR - DL )
      FluxILeft = Half * RelaxationWaveSpeedSquared * ( DL + DR ) &
        - Half * a * ( IR - IL )

      DL = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               +1 )
      DR = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iNext,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               U_R(iNodeZ2,iE,iNext,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS), &
               -1 )
      IL = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               +1 )
      IR = LinearFaceValue &
             ( U_R(iNodeZ1,iE,iNext,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               U_R(iNodeZ2,iE,iNext,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS), &
               -1 )
      FluxNRight = Half * ( IL + IR ) - Half * a * ( DR - DL )
      FluxIRight = Half * RelaxationWaveSpeedSquared * ( DL + DR ) &
        - Half * a * ( IR - IL )

      D1 = U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS)
      D2 = U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS)
      I1 = U_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS)
      I2 = U_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS)
      DMean = Half * ( D1 + D2 )
      IMean = Half * ( I1 + I2 )

      MeanRateN = - ( FluxNRight - FluxNLeft ) / h
      MeanRateI = - ( FluxIRight - FluxILeft ) / h
      SlopeRateN = 12.0_DP / h &
        * ( IMean - Half * ( FluxNLeft + FluxNRight ) )
      SlopeRateI = 12.0_DP / h &
        * ( RelaxationWaveSpeedSquared * DMean &
            - Half * ( FluxILeft + FluxIRight ) )

      dU_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS) &
        = MeanRateN - SlopeRateN / ( Two * SQRT(Three) )
      dU_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS) &
        = MeanRateN + SlopeRateN / ( Two * SQRT(Three) )
      dU_R(iNodeZ1,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS) &
        = MeanRateI - SlopeRateI / ( Two * SQRT(Three) )
      dU_R(iNodeZ2,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS) &
        = MeanRateI + SlopeRateI / ( Two * SQRT(Three) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_RelaxationDiffusion_Explicit


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Implicit, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, MatterSource_Rate_Option, &
      MatterIncrement_Rate_Option, MatterSource_Weight_Option, &
      MatterSource_Rate_Unweighted_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Implicit
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option(iZ_B1_Input(2):iZ_E1_Input(2), &
                 iZ_B1_Input(3):iZ_E1_Input(3), &
                 iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option(iZ_B1_Input(2):iZ_E1_Input(2), &
        iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterSource_Rate_Option(1:nQ,1:1+nSpecies+nC,1:nDOFX, &
        iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
        iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterIncrement_Rate_Option(1:nQ,1:nDOFX, &
        iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
        iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: MatterSource_Weight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterSource_Rate_Unweighted_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    CALL ComputeIncrement_TwoMoment_Implicit_Legacy &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             dt_Implicit, GE, GX, U_F, dU_F, U_R, dU_R, &
             uDR_Option = uDR_Option )

    IF( PRESENT(MatterCoupling_Mask_Option) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT(MatterSource_Rate_Option) ) &
      MatterSource_Rate_Option = Zero
    IF( PRESENT(MatterIncrement_Rate_Option) ) &
      MatterIncrement_Rate_Option = Zero
    IF( PRESENT(MatterSource_Rate_Unweighted_Option) ) &
      MatterSource_Rate_Unweighted_Option = Zero

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE ComputeIncrement_RelaxationDiffusion_Implicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Implicit, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, MatterSource_Rate_Option, &
      MatterIncrement_Rate_Option, MatterSource_Weight_Option, &
      MatterSource_Rate_Unweighted_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Implicit
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX,iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ,iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option(iZ_B1_Input(2):iZ_E1_Input(2), &
                 iZ_B1_Input(3):iZ_E1_Input(3), &
                 iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option(iZ_B1_Input(2):iZ_E1_Input(2), &
        iZ_B1_Input(3):iZ_E1_Input(3),iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterSource_Rate_Option(1:nQ,1:1+nSpecies+nC,1:nDOFX, &
        iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
        iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterIncrement_Rate_Option(1:nQ,1:nDOFX, &
        iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
        iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: MatterSource_Weight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterSource_Rate_Unweighted_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    INTEGER :: iE, iNodeZ, iS, iX1

    dU_F = Zero
    dU_R = Zero

    DO iS = 1, nSpecies
    DO iX1 = iZ_B0_Input(2), iZ_E0_Input(2)
    DO iE = iZ_B0_Input(1), iZ_E0_Input(1)
    DO iNodeZ = 1, nDOFZ
      dU_R(iNodeZ,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS) &
        = - ( U_R(iNodeZ,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_G1,iS) &
              - RelaxationAdvectionSpeed &
                * U_R(iNodeZ,iE,iX1,iZ_B0_Input(3),iZ_B0_Input(4),iCR_N,iS) ) &
            / ( RelaxationTime + dt_Implicit )
    END DO
    END DO
    END DO
    END DO

    IF( PRESENT(uDR_Option) ) uDR_Option = Zero
    IF( PRESENT(MatterCoupling_Mask_Option) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT(MatterSource_Rate_Option) ) &
      MatterSource_Rate_Option = Zero
    IF( PRESENT(MatterIncrement_Rate_Option) ) &
      MatterIncrement_Rate_Option = Zero
    IF( PRESENT(MatterSource_Rate_Unweighted_Option) ) &
      MatterSource_Rate_Unweighted_Option = Zero

  END SUBROUTINE ComputeIncrement_RelaxationDiffusion_Implicit


  SUBROUTINE ComputeError_RelaxationDiffusion( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: iNodeX, iNodeZ, iX1, Unit
    REAL(DP) :: DExact, IExact, DNumerical, INumerical, Weight, X
    REAL(DP) :: DErrorL1, DErrorL2, DErrorLInf, DNormL1, DNormL2
    REAL(DP) :: IErrorL1, IErrorL2, IErrorLInf, INormL1, INormL2

    DErrorL1 = Zero; DErrorL2 = Zero; DErrorLInf = Zero
    IErrorL1 = Zero; IErrorL2 = Zero; IErrorLInf = Zero
    DNormL1 = Zero; DNormL2 = Zero
    INormL1 = Zero; INormL2 = Zero

    OPEN( NEWUNIT = Unit, FILE = 'relaxation_diffusion_profile.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') '# x D_numerical D_exact I_numerical I_exact'

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX
      iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
      X = NodeCoordinate( MeshX(1), iX1, iNodeX )
      CALL RelaxationDiffusionExact( X, Time, DExact, IExact )
      DNumerical = uCR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
      INumerical = uCR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)
      Weight = MeshX(1) % Width(iX1) * WeightsX1(iNodeX)
      DErrorL1 = DErrorL1 + Weight * ABS( DNumerical - DExact )
      DErrorL2 = DErrorL2 + Weight * ( DNumerical - DExact )**2
      DErrorLInf = MAX( DErrorLInf, ABS( DNumerical - DExact ) )
      IErrorL1 = IErrorL1 + Weight * ABS( INumerical - IExact )
      IErrorL2 = IErrorL2 + Weight * ( INumerical - IExact )**2
      IErrorLInf = MAX( IErrorLInf, ABS( INumerical - IExact ) )
      DNormL1 = DNormL1 + Weight * ABS( DExact )
      DNormL2 = DNormL2 + Weight * DExact**2
      INormL1 = INormL1 + Weight * ABS( IExact )
      INormL2 = INormL2 + Weight * IExact**2
      WRITE(Unit,'(5(ES24.16E3,1X))') &
        X, DNumerical, DExact, INumerical, IExact
    END DO
    END DO
    CLOSE(Unit)

    DErrorL1 = DErrorL1 / MAX( DNormL1, TINY(One) )
    DErrorL2 = SQRT( DErrorL2 / MAX(DNormL2,TINY(One)) )
    DErrorLInf = DErrorLInf &
      / MAX( RelaxationBackground + ABS(RelaxationAmplitude), TINY(One) )
    IErrorL1 = IErrorL1 / MAX( INormL1, TINY(One) )
    IErrorL2 = SQRT( IErrorL2 / MAX(INormL2,TINY(One)) )
    IErrorLInf = IErrorLInf &
      / MAX( ABS(RelaxationAdvectionSpeed) * RelaxationBackground &
             + ABS(RelaxationAmplitude), TINY(One) )

    WRITE(*,*)
    WRITE(*,'(A)') 'Relaxation-diffusion error:'
    WRITE(*,'(A,3ES14.6E3)') '  D relative L1/L2/Linf = ', &
      DErrorL1, DErrorL2, DErrorLInf
    WRITE(*,'(A,3ES14.6E3)') '  I relative L1/L2/Linf = ', &
      IErrorL1, IErrorL2, IErrorLInf
    WRITE(*,'(A)') '  profile: relaxation_diffusion_profile.dat'

  END SUBROUTINE ComputeError_RelaxationDiffusion


  SUBROUTINE RelaxationDiffusionExact( X, Time, D, I )

    REAL(DP), INTENT(in) :: X, Time
    REAL(DP), INTENT(out) :: D, I

    REAL(DP) :: Amplitude, Phase, WaveNumber
    COMPLEX(DP) :: FluxEigenvector, Lambda

    WaveNumber = TwoPi / ( xR(1) - xL(1) )
    ! Slow Fourier eigenvalue of
    ! lambda^2 + lambda / epsilon + a^2 k^2 + i k v / epsilon = 0.
    Lambda = ( - One &
      + SQRT( CMPLX &
          ( One - 4.0_DP * RelaxationTime**2 &
                    * RelaxationWaveSpeedSquared * WaveNumber**2, &
            - 4.0_DP * RelaxationTime * WaveNumber &
                    * RelaxationAdvectionSpeed, KIND = DP ) ) ) &
      / ( Two * RelaxationTime )
    FluxEigenvector = CMPLX(Zero,One,KIND=DP) * Lambda / WaveNumber
    Amplitude = RelaxationAmplitude * EXP( REAL(Lambda,DP) * Time )
    Phase = WaveNumber * ( X - xL(1) ) + AIMAG(Lambda) * Time

    D = RelaxationBackground + Amplitude * COS( Phase )
    I = RelaxationAdvectionSpeed * RelaxationBackground &
      + Amplitude * ( REAL(FluxEigenvector,DP) * COS(Phase) &
                      - AIMAG(FluxEigenvector) * SIN(Phase) )

  END SUBROUTINE RelaxationDiffusionExact


  REAL(DP) FUNCTION LinearFaceValue( Q1, Q2, Side )

    REAL(DP), INTENT(in) :: Q1, Q2
    INTEGER, INTENT(in) :: Side

    LinearFaceValue = Half * ( Q1 + Q2 ) &
      + Half * DBLE(Side) * SQRT(Three) * ( Q2 - Q1 )

  END FUNCTION LinearFaceValue


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule, ONLY: &
      InitializeTimers
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
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange
    USE EquationOfStateModule, ONLY: &
      InitializeEquationOfState
    USE TwoMoment_ClosureModule, ONLY: &
      InitializeClosure_TwoMoment
    USE TwoMoment_OpacityModule, ONLY: &
      CreateOpacities
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

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 1, 1 ], &
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

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

    ! --- Initialize Equation of State ---

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = 'IDEAL', &
             Gamma_IDEAL_Option = 4.0_DP / 3.0_DP, &
             Verbose_Option = .TRUE. )

    ! --- Initialize Moment Closure ---

    CALL InitializeClosure_TwoMoment

    ! --- Initialize Opacities ---

    CALL CreateOpacities &
           ( nX, [ 1, 1, 1 ], nE, 1, Verbose_Option = .TRUE. )

    ! --- Initialize Troubled Cell Indicator ---

    CALL InitializeTroubledCellIndicator_TwoMoment &
           ( UseTroubledCellIndicator_Option &
               = UseTroubledCellIndicator, &
             C_TCI_Option &
               = C_TCI, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Slope Limiter ---

    CALL InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Positivity Limiter ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = SqrtTiny, &
             Min_2_Option &
               = SqrtTiny, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter, &
             UseEnergyLimiter_Option &
               = UseEnergyLimiter, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Tally ---

    CALL InitializeTally

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE TwoMoment_TimeSteppingModule, ONLY: &
      Finalize_IMEX_RK
    USE TwoMoment_TallyModule, ONLY: &
      FinalizeTally
    USE TwoMoment_OpacityModule, ONLY: &
      DestroyOpacities
    USE EquationOfStateModule, ONLY: &
      FinalizeEquationOfState
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
    USE TwoMoment_TimersModule, ONLY: &
      FinalizeTimers

    CALL Finalize_IMEX_RK

    CALL FinalizeTally

    CALL DestroyOpacities

    CALL FinalizeEquationOfState

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

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver
