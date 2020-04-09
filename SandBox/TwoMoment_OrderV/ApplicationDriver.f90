PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers
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
  USE ReferenceElementModuleZ, ONLY: &
    InitializeReferenceElementZ, &
    FinalizeReferenceElementZ
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
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
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeFromConserved_TwoMoment
  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
    InitializeTroubledCellIndicator_TwoMoment, &
    FinalizeTroubledCellIndicator_TwoMoment
  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    CreateOpacities, &
    SetOpacities, &
    DestroyOpacities
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(2)  :: Direction
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: TimeSteppingScheme
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  REAL(DP)      :: eL, eR, xL(3), xR(3)
  REAL(DP)      :: t, dt, t_end, V_0(3)
  REAL(DP)      :: D_0, Chi, Sigma
  REAL(DP)      :: LengthScale

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'TransparentVortex'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 2, 4, 16 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 0

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 1.0d-0
      iCycleD = 1
      iCycleW = 50
      maxCycles = 10000

      V_0 = [ 0.0_DP, 0.0_DP, 0.1_DP ]

      Direction = 'Z'

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .FALSE.

    CASE( 'SineWaveDiffusion' )

      nX  = [ 16, 1, 1 ]
      xL  = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ + 3.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 0, 0 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 0

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 2.0d1
      iCycleD = 10
      iCycleW = 10
      maxCycles = 1000000

      V_0 = [ 0.3_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 1.0d+2

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .FALSE.

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

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .FALSE.

    CASE( 'StreamingDopplerShift' )

      Direction = 'X' ! --- (X,Y, or Z)

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        nX  = [ 32, 1, 1 ]
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

      nE  = 16
      eL  = 0.0d0
      eR  = 5.0d1
      bcE = 10

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 2.5d+1
      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .TRUE.

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

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .TRUE.

    CASE( 'TransparentShock' )

      Direction = 'Y' ! --- (X,Y, or Z)

      LengthScale = 2.5d-4 ! --- Shock Width

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

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      t_end   = 5.0d0
      iCycleD = 1
      iCycleW = 250
      maxCycles = 1000000

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      UseSlopeLimiter = .TRUE.

      UsePositivityLimiter = .TRUE.

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

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .TRUE.

    CASE( 'GaussianDiffusion' )

      nX  = [ 48, 32, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, - 0.5_DP ]
      xR  = [ 3.0_DP, 2.0_DP, + 0.5_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0d0
      eR  = 1.0d0
      bcE = 0

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

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .TRUE.

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

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

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
           swE_Option &
             = 1, &
           bcE_Option &
             = bcE, &
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

  ! --- Initialize Timers ---

  CALL InitializeTimers

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

  ! --- Initialize Troubled Cell Indicator ---

  CALL InitializeTroubledCellIndicator_TwoMoment &
         ( UseTroubledCellIndicator_Option &
             = .FALSE., &
           C_TCI_Option &
             = 0.1_DP, &
           Verbose_Option &
             = .TRUE. )

  ! --- Initialize Slope Limiter ---

  CALL InitializeSlopeLimiter_TwoMoment &
         ( BetaTVD_Option = 2.0_DP, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           Verbose_Option &
             = .TRUE. )

  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = EPSILON( One ), &
           Min_2_Option = EPSILON( One ), &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           Verbose_Option = .TRUE. )

  ! --- Initialize Opacities ---

  CALL CreateOpacities &
         ( nX, [ 1, 1, 1 ], nE, 1, Verbose_Option = .TRUE. )

  CALL SetOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D_0, Chi, Sigma )

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  ! --- Set Initial Condition ---

  CALL InitializeFields( V_0, LengthScale, Direction )

  ! --- Apply Slope Limiter to Initial Data ---

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Apply Positivity Limiter to Initial Data ---

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  t = 0.0_DP
  dt = 0.3_DP * MINVAL( (xR-xL) / DBLE(nX) ) &
       / ( Two * DBLE(nNodes-1) + One )

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

    CALL Update_IMEX_RK( dt, uGE, uGF, uCF, uCR )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

    END IF

  END DO

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Finalize ---

  CALL FinalizeTroubledCellIndicator_TwoMoment

  CALL FinalizeSlopeLimiter_TwoMoment

  CALL FinalizePositivityLimiter_TwoMoment

  CALL DestroyOpacities

  CALL Finalize_IMEX_RK

  CALL FinalizeEquationOfState

  CALL FinalizeTimers

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElementZ

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
