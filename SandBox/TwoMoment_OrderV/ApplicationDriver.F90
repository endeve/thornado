PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Pi, TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeFromConserved_TwoMoment
  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    SetOpacities
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields, &
    ComputeError
  USE TwoMoment_TallyModule_OrderV, ONLY: &
    ComputeTally

  IMPLICIT NONE

  CHARACTER(2)  :: Direction
  CHARACTER(32) :: Spectrum
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: TimeSteppingScheme
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseEnergyLimiter
  LOGICAL       :: UseTroubledCellIndicator
  INTEGER       :: nNodes
  INTEGER       :: nSpecies
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t, dt, t_end, dt_CFL, dt_0, dt_grw, V_0(3)
  REAL(DP)      :: D_0, Chi, Sigma, C_TCI
  REAL(DP)      :: LengthScale

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'SineWaveStreaming'

  nSpecies = 1

  C_TCI = 1.0_DP
  UseTroubledCellIndicator = .FALSE.

  dt_0   = HUGE( One )
  dt_grw = One

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 8, 8, 8 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 8
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nSpecies = 6
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

    CASE( 'SineWaveDiffusion' )

      nX  = [ 16, 1, 1 ]
      xL  = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ + 3.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 1.0d-1
      iCycleD = 10
      iCycleW = 10
      maxCycles = 1000000

      V_0 = [ 0.3_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 1.0d+2

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

      nE    = 16
      eL    = 0.0d0
      eR    = 5.0d1
      bcE   = 10
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

      LengthScale = 1.0d-1 ! --- Shock Width

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
      bcE = 10

      nNodes = 1

      TimeSteppingScheme = 'SSPRK1'

      t_end   = 5.0d0
      iCycleD = 1
      iCycleW = 10
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .TRUE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'TransparentVortex' )

      Direction = 'X' ! --- (X or Y)

      nX  = [ 16, 16, 1 ]
      xL  = [ - 5.0_DP, - 5.0_DP, - 0.5_DP ]
      xR  = [ + 5.0_DP, + 5.0_DP, + 0.5_DP ]

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        bcX = [ 12, 1, 1 ]

      ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

        bcX = [ 1, 12, 1 ]

      ELSE

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'TransparentVortex.  Direction must be X or Y'
        WRITE(*,*)
        STOP

      END IF

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      nE  = 16
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 10

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 4.0d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .FALSE.

    CASE( 'RadiatingSphere' )

      CoordinateSystem = 'SPHERICAL'

      nX    = [ 200, 1, 1 ]
      xL    = [ 1.0d1, Zero,  Zero ]
      xR    = [ 1.0d4,   Pi, TwoPi ]
      bcX   = [ 12, 1, 1 ]
      ZoomX = [ 1.024333847373375_DP, One, One ]

      nE    = 16
      eL    = 0.0d0
      eR    = 3.0d2
      bcE   = 2
      ZoomE = 1.310262775587271_DP

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end = 2.0d+4
      iCycleD = 1
      iCycleW = 2000
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

    CASE( 'GaussianDiffusion' )

      nX  = [ 48, 32, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, - 0.5_DP ]
      xR  = [ 3.0_DP, 2.0_DP, + 0.5_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 1

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 5.0d0
      iCycleD = 10
      iCycleW = 10
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 1.0d+2

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
      bcE = 1

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
      UseEnergyLimiter     = .FALSE.

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

  CALL InitializeFields( V_0, LengthScale, Direction, Spectrum )

  t = 0.0_DP

  ! --- Apply Slope Limiter to Initial Data ---

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Apply Positivity Limiter to Initial Data ---

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Write Initial Condition ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( uGE, uGF, uCF, uCR )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( uGE, uGF, uCF, uCR )
#endif

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

    dt_CFL = 0.3_DP * MINVAL( (xR-xL)/DBLE(nX) ) / ( Two*DBLE(nNodes-1)+One )

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

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( uGF, uCF, uCR )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( uGF, uCF, uCR )
#endif

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL ComputeTally( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

    END IF

  END DO

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( uGF, uCF, uCR )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( uGF, uCF, uCR )
#endif

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ComputeTally( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

  CALL ComputeError( t )

  ! --- Auxiliary Finalization ---

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule_OrderV, ONLY: &
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
    USE TwoMoment_OpacityModule_OrderV, ONLY: &
      CreateOpacities
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      InitializeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
      InitializeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE TwoMoment_TallyModule_OrderV, ONLY: &
      InitializeTally
    USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
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

    USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
      Finalize_IMEX_RK
    USE TwoMoment_TallyModule_OrderV, ONLY: &
      FinalizeTally
    USE TwoMoment_OpacityModule_OrderV, ONLY: &
      DestroyOpacities
    USE EquationOfStateModule, ONLY: &
      FinalizeEquationOfState
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      FinalizeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
      FinalizeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
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
    USE TwoMoment_TimersModule_OrderV, ONLY: &
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
