PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Three, &
    Four, &
    Pi, &
    TwoPi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram_Basic, &
    FinalizeProgram_Basic
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
  USE InitializationModule, ONLY: &
    InitializeFields_Relativistic_MHD
  USE MHD_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_MHD_Relativistic_IDEAL, &
    FinalizeSlopeLimiter_MHD_Relativistic_IDEAL, &
    ApplySlopeLimiter_MHD_Relativistic_IDEAL
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_MHD_Relativistic, &
    ComputeTimeStep_MHD_Relativistic, &
    ComputeMagneticDivergence_MHD_Relativistic
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE MagnetofluidFieldsModule, ONLY: &
    uCM, &
    uPM, &
    uAM, &
    uDM
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE MHD_DiscretizationModule_Relativistic, ONLY: &
    ComputeIncrement_MHD_DG_Explicit
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeMagnetofluid_SSPRK, &
    FinalizeMagnetofluid_SSPRK, &
    UpdateMagnetofluid_SSPRK
  USE UnitsModule, ONLY: &
    Centimeter, &
    Millisecond, &
    UnitsDisplay

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: AdvectionProfile
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: SmoothProfile, ConstantDensity
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UseConservativeCorrection
  REAL(DP)      :: SlopeTolerance
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: ZoomX(3)
  REAL(DP)      :: Timer_Evolution

  LOGICAL  :: WriteGF = .TRUE., WriteMF = .TRUE.
  LOGICAL  :: ActivateUnits = .FALSE.
  LOGICAL  :: UseMHD = .TRUE.
  LOGICAL  :: EvolveOnlyMagnetic = .FALSE.
  LOGICAL  :: UseDivergenceCleaning = .FALSE.
  LOGICAL  :: UsePowellSource = .FALSE.

  REAL(DP) :: DampingParameter = 0.0_DP
  REAL(DP) :: Angle = 0.0_DP

  REAL(DP) :: MMBlastWaveB0  = 0.5_DP
  REAL(DP) :: MMBlastWavePhi = 0.0_DP

  REAL(DP) :: OTScaleFactor = 100.0_DP

  ProgramName = 'Advection1D'
  AdvectionProfile = 'HydroSineWaveX1'
  RiemannProblemName = 'IsolatedContact'

  swX               = [ 0, 0, 0 ]
  RestartFileNumber = -1
  t                 = 0.0_DP
  ZoomX             = 1.0_DP

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection1D' )

      SELECT CASE ( TRIM( AdvectionProfile ) )

        CASE( 'HydroSineWaveX1' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'HydroSineWaveX1'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 0, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 1, 1 ]
          swX = [ 1, 0, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'MagneticSineWaveX1' )

          EvolveOnlyMagnetic = .TRUE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'MagneticSineWaveX1'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 0, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 1, 1 ]
          swX = [ 1, 0, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'CPAlfvenX1' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'CPAlfvenX1'

          UsePowellSource = .FALSE.

          Gamma = 4.0_DP / 3.0_DP
          t_end = 16.449592691810107_DP
          bcX = [ 1, 0, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 1, 1 ]
          swX = [ 1, 0, 0 ]
          xL  = [ 0.0_DP,   0.0_DP, 0.0_DP ]
          xR  = [ Two * Pi, 1.0_DP, 1.0_DP ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid AdvectionProfile: ', AdvectionProfile
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  HydroSineWaveX1'
          WRITE(*,'(A)')     '  MagneticSineWaveX1'
          WRITE(*,'(A)')     '  CPAlfvenX1'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

    CASE( 'Advection2D' )

      SELECT CASE( TRIM( AdvectionProfile ) )

        CASE( 'HydroSineWaveX2' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'HydroSineWaveX2'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'HydroSineWaveX1X2' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'HydroSineWaveX1X2'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP / SQRT( Two ), 1.0_DP / SQRT( Two ), 1.0_DP ]

        CASE( 'MagneticSineWaveX2' )

          EvolveOnlyMagnetic = .TRUE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'MagneticSineWaveX2'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'MagneticSineWaveX1X2' )

          EvolveOnlyMagnetic = .TRUE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'MagneticSineWaveX1X2'


          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP / SQRT( Two ), 1.0_DP / SQRT( Two ), 1.0_DP ]

        CASE( 'CPAlfvenX2' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'CPAlfvenX2'

          Gamma = 4.0_DP / 3.0_DP
          t_end = 16.449592691810107_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ 0.0_DP,   0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, Two * Pi, 1.0_DP ]

        CASE( 'CPAlfvenOblique' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          Angle = ATAN(TwoPi)

          AdvectionProfile = 'CPAlfvenOblique'

          Gamma = 4.0_DP / 3.0_DP
          t_end = 16.449592691810107_DP / TwoPi
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [    128,                  128,                  1 ]
          swX = [    1,                  1,                  0 ]
          xL  = [ Zero,               Zero,               -One ]
          xR  = [ One / COS( Angle ), One / SIN( Angle ),  One ]

        CASE( 'LoopAdvection' )

          EvolveOnlyMagnetic = .TRUE.
          UseDivergenceCleaning = .FALSE.
          DampingParameter = 0.0_DP

          UsePowellSource = .FALSE.

          AdvectionProfile = 'LoopAdvection'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 2.0_DP
          bcX = [ 1, 1, 0 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 128, 128, 1 ]
          swX = [ 1, 1, 0 ]
          xL  = [ -0.5_DP, -0.5_DP, 0.0_DP ]
          xR  = [  0.5_DP,  0.5_DP, 1.0_DP ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid AdvectionProfile: ', AdvectionProfile
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  HydroSineWaveX2'
          WRITE(*,'(A)')     '  HydroSineWaveX1X2'
          WRITE(*,'(A)')     '  MagneticSineWaveX1X2'
          WRITE(*,'(A)')     '  CPAlfvenX2'
          WRITE(*,'(A)')     '  CPAlfvenOblique'
          WRITE(*,'(A)')     '  LoopAdvection'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

    CASE( 'Advection3D' )

      SELECT CASE ( TRIM( AdvectionProfile ) )

        CASE( 'HydroSineWaveX3' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'HydroSineWaveX3'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 1 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 2, 128 ]
          swX = [ 1, 1, 1 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'MagneticSineWaveX3' )

          EvolveOnlyMagnetic = .TRUE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'MagneticSineWaveX3'

          Gamma = 5.0_DP / 3.0_DP
          t_end = 10.0_DP
          bcX = [ 1, 1, 1 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 2, 128 ]
          swX = [ 1, 1, 1 ]
          xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
          xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        CASE( 'CPAlfvenX3' )

          EvolveOnlyMagnetic = .FALSE.
          UseDivergenceCleaning = .FALSE.

          AdvectionProfile = 'CPAlfvenX3'

          Gamma = 4.0_DP / 3.0_DP
          t_end = 16.449592691810107_DP
          bcX = [ 1, 1, 1 ]

          CoordinateSystem = 'CARTESIAN'

          nX  = [ 2, 2, 8 ]
          swX = [ 1, 1, 1 ]
          xL  = [ 0.0_DP,   0.0_DP, 0.0_DP ]
          xR  = [ Two * Pi, Two * Pi, Two * Pi ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid AdvectionProfile: ', AdvectionProfile
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  HydroSineWaveX3'
          WRITE(*,'(A)')     '  MagneticSineWaveX3'
          WRITE(*,'(A)')     '  CPAlfvenX3'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

    CASE( 'Cleaning1D' )

      EvolveOnlyMagnetic = .TRUE.
      UseDivergenceCleaning = .TRUE.
      DampingParameter = 0.0_DP

      SmoothProfile = .TRUE.

      Gamma = 1.4_DP
      t_end = 10.0_DP
      bcX = [ 1, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 128, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ -1.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [  1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'Cleaning2D' )

      EvolveOnlyMagnetic = .TRUE.

      UseDivergenceCleaning = .TRUE.
      DampingParameter = 0.0_DP

      ConstantDensity = .TRUE.

      Gamma = 5.0_DP / 3.0_DP
      t_end = 1.0_DP
      bcX = [ 1, 1, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 128, 128, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ -0.5_DP, -0.5_DP, 0.0_DP ]
      xR  = [  1.5_DP,  1.5_DP, 1.0_DP ]

    CASE( 'Riemann1D' )

      EvolveOnlyMagnetic = .FALSE.

      UseDivergenceCleaning = .FALSE.
      DampingParameter = 0.0_DP

      UsePowellSource = .FALSE.

      bcX = [ 2, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'IsolatedContact' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = One

          nX  = [ 40, 1, 1 ]

        CASE( 'RotationalWave' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = One

          nX  = [ 400, 1, 1 ]

        CASE( 'ShockTube1' )

          Gamma = Two
          t_end = 0.4_DP

          nX  = [ 400, 1, 1 ]

        CASE( 'ShockTube2' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.55_DP

          nX  = [ 800, 1, 1 ]

        CASE( 'HydroShockTube3' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4_DP

          nX = [ 400, 1, 1 ]

        CASE( 'ShockTube3' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4_DP

          nX  = [ 400, 1, 1 ]

        CASE( 'ShockTube4' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.55_DP

          nX  = [ 800, 1, 1 ]

        CASE( 'HydroMBProblem1' )

          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.4_DP

          nX  = [ 128, 1, 1 ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid RiemannProblemName: ', RiemannProblemName
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  IsolatedContact'
          WRITE(*,'(A)')     '  RotationalWave'
          WRITE(*,'(A)')     '  ShockTube1'
          WRITE(*,'(A)')     '  ShockTube2'
          WRITE(*,'(A)')     '  HydroShockTube3'
          WRITE(*,'(A)')     '  ShockTube3'
          WRITE(*,'(A)')     '  ShockTube4'
          WRITE(*,'(A)')     '  HydroMBProblem1'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

    CASE( 'MMBlastWave2D' )

      EvolveOnlyMagnetic = .FALSE.

      UseDivergenceCleaning = .FALSE.
      DampingParameter = 0.0_DP

      MMBlastWaveB0  = 0.5_DP
      MMBlastWavePhi = 0.0_DP

      Gamma = 4.0_DP / 3.0_DP
      t_end = 4.0_DP
      bcX = [ 2, 2, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 200, 200, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ -6.0_DP, -6.0_DP, 0.0_DP ]
      xR  = [  6.0_DP,  6.0_DP, 1.0_DP ]

    CASE( 'OrszagTang2D' )

      EvolveOnlyMagnetic = .FALSE.

      UseDivergenceCleaning = .FALSE.
      DampingParameter = 0.0_DP

      OTScaleFactor = 100.0_DP

      Gamma = 5.0_DP / 3.0_DP
      t_end = 50.0_DP
      bcX = [ 1, 1, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 256, 256, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'ShearingDisk_Unperturbed' )

      ActivateUnits = .TRUE.

      EvolveOnlyMagnetic = .FALSE.

      UseDivergenceCleaning = .FALSE.
      DampingParameter = 0.0_DP

      Gamma = 4.0_DP / 3.0_DP
      t_end = 10e0 * Millisecond
      bcX = [ 0, 0, 0 ]

      CoordinateSystem = 'CYLINDRICAL'

      nX  = [ 128, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 1.45d6 * Centimeter, -0.5d5 * Centimeter, Zero ]
      xR  = [ 1.65d6 * Centimeter,  0.5d5 * Centimeter, TwoPi]

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  Advection1D'
      WRITE(*,'(A)')     '  Advection2D'
      WRITE(*,'(A)')     '  Advection3D'
      WRITE(*,'(A)')     '  Cleaning1D'
      WRITE(*,'(A)')     '  Cleaning2D'
      WRITE(*,'(A)')     '  Riemann1D'
      WRITE(*,'(A)')     '  MMBlastWave2D'
      WRITE(*,'(A)')     '  OrszagTang2D'
      WRITE(*,'(A)')     '  ShearingDisk_Unperturbed'
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
  UseConservativeCorrection = .FALSE.

  ! === End of User Input ===

  CALL InitializeProgram_Basic &
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
           UseMHD_Option & 
             = UseMHD )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_MHD_Relativistic_IDEAL &
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
           UseConservativeCorrection_Option &
             = UseConservativeCorrection )

  CALL InitializeMagnetofluid_SSPRK &
         ( nStages = nStagesSSPRK, EvolveOnlyMagnetic_Option = EvolveOnlyMagnetic, &
           UseDivergenceCleaning_Option = UseDivergenceCleaning, &
           DampingParameter_Option = DampingParameter, &
           UsePowellSource_Option = UsePowellSource )

  WRITE(*,*)
  WRITE(*,'(A6,A,ES11.3E3)') '', 'CFL: ', CFL

  uCM = Zero ! Without this, crashes when copying data in TimeStepper
  uDM = Zero ! Without this, crashes in IO

  CALL InitializeFields_Relativistic_MHD &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ), &
           SmoothProfile_Option &
             = SmoothProfile, &
           ConstantDensity_Option &
             = ConstantDensity, &
           Angle_Option &
             = Angle, &
           MMBlastWaveB0_Option &
             = MMBlastWaveB0, &
           MMBlastWavePhi_Option &
             = MMBlastWavePhi, &
           OTScaleFactor_Option &
             = OTScaleFactor, &
           EvolveOnlyMagnetic_Option &
             = EvolveOnlyMagnetic )

  IF( RestartFileNumber .LT. 0 )THEN

    CALL ApplySlopeLimiter_MHD_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uDM )

    CALL ComputeFromConserved_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uPM, uAM, &
             EvolveOnlyMagnetic )

    CALL ComputeMagneticDivergence_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uDM )

    CALL WriteFieldsHDF &
           ( t, WriteGF_Option = WriteGF, WriteMF_Option = WriteMF )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadMF_Option = .TRUE., ReadGF_Option = .TRUE. )

  END IF

  iCycleD = 10
  dt_wrt = 1.0e-2 * ( t_end - t ); iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'

  t_wrt = t + dt_wrt
  wrt = .FALSE.

  iCycle = 0
  Timer_Evolution = MPI_WTIME()
  DO WHILE( t .LT. t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF, uCM, &
             CFL / ( nDimsX * ( Two * DBLE( nNodes ) - One ) ), &
             dt, UseDivergenceCleaning, EvolveOnlyMagnetic )

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

    CALL UpdateMagnetoFluid_SSPRK &
           ( t, dt, uGF, uCM, uDM, &
             ComputeIncrement_MHD_DG_Explicit )

    IF( iCycleW .GT. 0 )THEN

      IF( MOD( iCycle, iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.

    ELSE

      IF( t + dt .GT. t_wrt )THEN

        t_wrt = t_wrt + dt_wrt
        wrt = .TRUE.

      END IF

    END IF

    IF( wrt )THEN

      CALL ComputeFromConserved_MHD_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uPM, uAM, EvolveOnlyMagnetic )

      CALL ComputeMagneticDivergence_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uDM )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteMF_Option = WriteMF )

      wrt = .FALSE.

    END IF

  END DO

  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,I8.8,A,ES10.3E3,A)') &
    'Finished ', iCycle, ' cycles in ', Timer_Evolution, ' s'
  WRITE(*,*)

  CALL ComputeFromConserved_MHD_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uPM, uAM, EvolveOnlyMagnetic )

  CALL ComputeMagneticDivergence_MHD_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uDM )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteMF_Option = WriteMF )

  CALL FinalizeMagnetofluid_SSPRK

  CALL FinalizeEquationOfState

  CALL FinalizeSlopeLimiter_MHD_Relativistic_IDEAL

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram_Basic( UseMHD_Option = UseMHD )

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
