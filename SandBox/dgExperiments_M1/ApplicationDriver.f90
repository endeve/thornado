PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, SqrtTiny, Zero, One, Pi, TwoPi, Third
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
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
  USE TimeSteppingModule_IMEX_RK, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE InitializationModule, ONLY: &
    InitializeFields, &
    ComputeError
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment, &
    TallyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE TwoMoment_TallyModule, ONLY: &
    InitializeTally_TwoMoment, &
    ComputeTally_TwoMoment
  USE dgDiscretizationModule_Collisions, ONLY: &
    InitializeCollisions, &
    FinalizeCollisions, &
    ComputeIncrement_M1_DG_Implicit, &
    ComputeCorrection_M1_DG_Implicit

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(8)  :: Direction
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: iCycle, iCycleD, iCycleW, iCycleT, maxCycles
  INTEGER       :: nE, nX(3), bcX(3), nNodes
  REAL(DP)      :: t, dt, t_end, wTime, dt_para
  REAL(DP)      :: xL(3), xR(3), ZoomX(3)
  REAL(DP)      :: eL,    eR
  REAL(DP)      :: N0, SigmaA, SigmaS
  REAL(DP)      :: Radius = 1.0d16
  REAL(DP)      :: Min_1, Max_1, Min_2

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'SineWaveDiffusion'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      Direction = 'X'

      nX = [ 16, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS_3'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+1
      iCycleD   = 10
      iCycleW   = 10
      iCycleT   = 10
      maxCycles = 10000


    CASE( 'SquareWaveStreaming' )

      ! --- Minerbo Closure Only ---

      Direction = 'X'

      nX = [ 64, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+0
      iCycleD   = 100
      iCycleW   = 100
      iCycleT   = 100
      maxCycles = 1000000

    CASE( 'SineWaveDamping' )

      ! --- Minerbo Closure Only ---

      nX = [ 32, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS_3'

      N0     = 0.0_DP
      SigmaA = 1.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+1
      iCycleD   = 10
      iCycleW   = 100
      iCycleT   = 10
      maxCycles = 100000

    CASE( 'SineWaveDiffusion' )

      nX = [ 32, 1, 1 ]
      xL = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ + 3.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS_3'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 1.0d+2

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+2
      iCycleD   = 10
      iCycleW   = 2000
      iCycleT   = 10
      maxCycles = 1000000

    CASE( 'PackedBeam' )

      nX = [ 400, 1, 1 ]
      xL = [ - 1.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ + 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PC2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero ! --- Min Density
      Max_1 = One  ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 8.5d-1
      iCycleD   = 10
      iCycleW   = 10
      iCycleT   = 10
      maxCycles = 10000

    CASE( 'LineSource' )

      nX = [ 512, 512, 1 ]
      xL = [ - 1.28_DP, - 1.28_DP, 0.0_DP ]
      xR = [ + 1.28_DP, + 1.28_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero ! --- Min Density
      Max_1 = One  ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 1.0d+0
      iCycleD   = 10
      iCycleW   = 50
      iCycleT   = 10
      maxCycles = 1000000

    CASE( 'RiemannProblem' )

    ! --- Ref: Olbrant et al. (2012) ---
    ! --- JCP 231(17)  -----------------

      Direction = 'X'

      SELECT CASE ( Direction )
      CASE ( 'X' )
  
        nX = [ 240, 1, 1 ]
        xL = [ - 0.05_DP,  0.0_DP, 0.0_DP ]
        xR = [ + 0.1_DP,   1.0_DP, 1.0_DP ]
  
        bcX = [ 2, 1, 1 ]

      CASE ( 'Y' )
  
        nX = [ 2, 240, 1 ] ! need nonzero nX(1) to trigger limiter in X2
        xL = [   0.0_DP, - 0.05_DP, 0.0_DP ]
        xR = [ + 1.0_DP, + 0.1_DP, 1.0_DP ]
  
        bcX = [ 1, 2, 1 ]

      END SELECT

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero         ! --- Min Density
      Max_1 = HUGE( ONE )  ! --- Max Density !! not done
      Min_2 = Zero         ! --- Min "Gamma"

      t_end     = 1.0d-1
      iCycleD   = 1000
      iCycleW   = 1000
      iCycleT   = 100
      maxCycles = 1000000

    CASE( 'HomogeneousSphere' )

      nX = [ 64, 64, 64 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, 2.0_DP, 2.0_DP ]

      bcX = [ 32, 32, 32 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS_3'

      N0     = 1.00_DP - 1.0d-12
      SigmaA = 1000.0_DP
      SigmaS = 0.00_DP
      Radius = 1.00_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero + 1.0d-14 ! --- Min Density
      Max_1 = One  - 1.0d-14 ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 5.0d-0
      iCycleD   = 10
      iCycleW   = 500
      iCycleT   = 10
      maxCycles = 100000

    CASE( 'HomogeneousSphere_Spherical' )

      CoordinateSystem = 'SPHERICAL'

      nX = [ 128, 8, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 5.0_DP, Pi,     TwoPi  ]

      bcX = [ 32, 3, 0 ]
 
      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'IMEX_PDARS_3'

      N0     = 0.80_DP
      SigmaA = 4.00_DP
      SigmaS = 0.00_DP
      Radius = 1.00_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero + SqrtTiny ! --- Min Density
      Max_1 = One  - SqrtTiny ! --- Max Density
      Min_2 = Zero + SqrtTiny ! --- Min "Gamma"

      t_end     = 2.0d+1 
      iCycleD   = 10
      iCycleW   = 500
      iCycleT   = 10
      maxCycles = 100000

  END SELECT

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
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
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

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Initialize Implicit Solver ---

  CALL InitializeCollisions &
         ( N0_Option = N0, SigmaA0_Option = SigmaA, &
           SigmaS0_Option = SigmaS, Radius_Option = Radius )

  ! --- Set Initial Condition ---

  CALL InitializeFields &
         ( Direction_Option = TRIM( Direction ), &
           SigmaA_Option = SigmaA, &
           SigmaS_Option = SigmaS )

  ! --- Initialize Slope Limiter ---

  CALL InitializeSlopeLimiter_TwoMoment                &
         ( BetaTVD_Option = 2.0_DP,                    &
           BetaTVB_Option = 0.d0,                      &
           SlopeTolerance_Option = 1.0d-6,             &
           UseSlopeLimiter_Option = .FALSE.,           &
           UseCharacteristicLimiting_Option = .FALSE., &
           Verbose_Option = .TRUE. )
 
  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = Min_1, &
           Max_1_Option = Max_1, &
           Min_2_Option = Min_2, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           UsePositivityLimiterTally_Option &
             = .TRUE. )

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )  

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Tally ---

  CALL InitializeTally_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t  = 0.0d-0

  IF ( CoordinateSystem == 'SPHERICAL') THEN
     dt = 0.1_DP * MINVAL( MeshX(1) % Width(1:nX(1)) ) &
            / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP ) 
  ELSE IF ( CoordinateSystem == 'CARTESIAN' ) THEN
    dt  = 0.1_DP * MINVAL( (xR-xL) / DBLE( nX ) ) &
            / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )
  END IF

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t, ' to t = ', t_end
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCR, &
             ComputeIncrement_TwoMoment_Explicit, &
             ComputeIncrement_M1_DG_Implicit, &
             ComputeCorrection_M1_DG_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleT ) == 0 )THEN

      CALL TallyPositivityLimiter_TwoMoment( t )

      CALL ComputeTally_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR, t, &
             iState_Option = 1, DisplayTally_Option = .TRUE. )

    END IF

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

    END IF

  END DO

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .FALSE., &
           WriteRF_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL ComputeError( Time = t, SigmaA = SigmaA, SigmaS = SigmaS )

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL Finalize_IMEX_RK

  CALL FinalizeCollisions

  CALL FinalizeSlopeLimiter_TwoMoment

  CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
