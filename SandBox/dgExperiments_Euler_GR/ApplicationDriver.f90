PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, Four, TwoPi
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
    iX_B0, iX_B1, iX_E0, iX_E1
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_GR, ONLY: &
    InitializeFields_GR, ReadParameters
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    InitializeSlopeLimiter_Euler_GR, &
    FinalizeSlopeLimiter_Euler_GR, &
    ApplySlopeLimiter_Euler_GR
  USE PositivityLimiterModule_Euler_GR, ONLY: &
    InitializePositivityLimiter_Euler_GR, &
    FinalizePositivityLimiter_Euler_GR, &
    ApplyPositivityLimiter_Euler_GR
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeFromConserved_GR, &
    ComputeTimeStep_GR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeIncrement_Euler_GR_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    Millisecond
  USE Euler_TallyModule_GR, ONLY: &
    InitializeTally_Euler_GR, &
    FinalizeTally_Euler_GR, &
    ComputeTally_Euler_GR

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: RiemannProblemName, RiemannProblem2dName
  CHARACTER(32) :: SphericalRiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseFixed_dt, UseSourceTerm
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: iCycle, iCycleD, iCycleW = 0
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, wTime, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter

  ! --- Sedov blast wave ---
  REAL(DP) :: Eblast
  INTEGER  :: nDetCells
  REAL(DP) :: Vmax, LorentzFactor

  ! --- Standing accretion shock ---
  REAL(DP), ALLOCATABLE :: FluidFieldParameters(:)
  REAL(DP)              :: M_PNS = Zero, Ri, R_PNS, R_shock, Rf

  LOGICAL :: DEBUG = .FALSE.
  REAL(DP) :: wTime_UF, wTime_CTS
  
!  ProgramName = 'RiemannProblem'
!  ProgramName = 'RiemannProblem2d'
!  ProgramName = 'SphericalRiemannProblem'
  ProgramName = 'SphericalSedov'
!  ProgramName = 'StandingAccretionShock'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'Sod' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem1' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem4' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'PerturbedShockTube' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.35d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'CartesianSedov' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.2d0
          bcX   = [ 32, 0, 0 ]
          nDetCells = 1
          Eblast    = 1.0d-3

        CASE( 'ShockReflection' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.75d0
          bcX   = [ 23, 0, 0 ]

        CASE DEFAULT ! Sod
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 256, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblem2d' )

      RiemannProblem2dName = 'DzB2002'

      SELECT CASE ( TRIM( RiemannProblem2dName ) )

        CASE( 'DzB2002' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 2, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 64, 64, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'SphericalRiemannProblem' )

      SphericalRiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 128, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, Pi, TwoPi ]

      swX = [ 1, 0, 0 ]
      bcX = [ 2, 0, 0 ]

      t_end = 5.0d-1

    CASE( 'SphericalSedov' )

      nDetCells = 1
      Eblast    = 1.0d-1

      CoordinateSystem = 'SPHERICAL'

      Gamma = 1.4_DP!4.0_DP / 3.0_DP

      nX = [ 256, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.2_DP, Pi, TwoPi ]

      swX = [ 1, 0, 0 ]
      bcX = [ 2, 0, 0 ]

      t_end = 1.0d0

    CASE( 'StandingAccretionShock' )

      CoordinateSystem = 'SPHERICAL'

      CALL ReadParameters &
             ( '../StandingAccretionShock_Parameters.dat', &
                 FluidFieldParameters )

      M_PNS   = FluidFieldParameters(1)
      Gamma   = FluidFieldParameters(2)
      Ri      = FluidFieldParameters(3)
      R_PNS   = FluidFieldParameters(4)
      R_shock = FluidFieldParameters(5)

      nX = [ 256, 1, 1 ]
      xL = [ R_PNS, 0.0_DP, 0.0_DP ]
      xR = [ Two * R_shock, Pi, Four ]

      swX = [ 1, 0, 0 ]
      bcX = [ 11, 0, 0 ]

      t_end = 1.0d1 * Millisecond

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblem2d'
      WRITE(*,'(A)')     '  SphericalRiemannProblem'
      WRITE(*,'(A)')     '  SphericalSedov'
      WRITE(*,'(A)')     '  StandingAccretionShock'
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  nNodes = 3
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  BetaTVD = 1.75d0
  BetaTVB = 0.0d0

  UseSlopeLimiter           = .TRUE.
  SlopeTolerance            = 1.0d-6
  UseCharacteristicLimiting = .TRUE.

  UseTroubledCellIndicator  = .TRUE.
  LimiterThresholdParameter = 0.015_DP

  UsePositivityLimiter = .TRUE.
  Min_1 = Zero
  Min_2 = Zero

  UseFixed_dt   = .FALSE.
  UseSourceTerm = .TRUE.

  iCycleD = 100
!!$  iCycleW = 1; dt_wrt = -1.0d0
  dt_wrt  = 1.0d-2 * t_end;       iCycleW = -1
!!$  dt_wrt  = 1.0d-1 * Millisecond; iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 )THEN
    STOP 'dt_wrt and iCycleW cannot both be present'
  END IF

  nStagesSSPRK = 3
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  ! --- Cockburn & Shu, (2001), JSC, 16, 173 ---
  CFL = 0.5d0 / ( 2.0d0 * DBLE( nNodes - 1 ) + 1.0d0 )

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
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = M_PNS )

  CALL InitializeFields_GR &
         ( RiemannProblemName_Option = TRIM( RiemannProblemName ), &
           RiemannProblem2dName_Option = TRIM( RiemannProblem2dName ), &
           SphericalRiemannProblemName_Option &
             = TRIM( SphericalRiemannProblemName ), &
           nDetCells_Option = nDetCells, Eblast_Option = Eblast )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )

  CALL InitializeSlopeLimiter_Euler_GR &
         ( BetaTVD_Option = BetaTVD, &
           BetaTVB_Option = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter )

  CALL InitializePositivityLimiter_Euler_GR &
         ( Min_1_Option = Min_1, &
           Min_2_Option = Min_2, &
           UsePositivityLimiter_Option = UsePositivityLimiter )

  WRITE(*,*)
  WRITE(*,'(A2,A,F4.2)') '', 'CFL: ', CFL
  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Evolving Fields...'
  WRITE(*,*)

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL InitializeTally_Euler_GR'
  CALL InitializeTally_Euler_GR &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  iCycle = 0
  DO WHILE( t .LT. t_end )

    iCycle = iCycle + 1

    IF( DEBUG ) wTime_CTS = MPI_WTIME( )
    IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL ComputeTimeStep_GR'
    IF( .NOT. UseFixed_dt )THEN
      CALL ComputeTimeStep_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               CFL = CFL, TimeStep = dt, UseSourceTerm_Option = UseSourceTerm )
    ELSE
      dt = CFL * ( xR(1) - xL(1) ) / DBLE( nX(1) )
    END IF

    IF( DEBUG ) wTime_CTS = MPI_WTIME( ) - wTime_CTS

    IF( t + dt .LT. t_end )THEN
      t = t + dt
    ELSE
      dt = t_end - t
      t  = t_end
    END IF

    IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN

      IF( ProgramName .EQ. 'StandingAccretionShock' )THEN
        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t / Millisecond, ' ms ', &
          'dt = ', dt / Millisecond, ' ms'
      ELSE
        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt
      END IF

      IF( DEBUG )THEN
        WRITE(*,'(A8,A,ES10.3E3,A)') &
          '', 'Time to call ComputeTimeStep:   ', wTime_CTS, ' s'
        WRITE(*,'(A8,A,ES10.3E3,A)') &
          '', 'Time to call UpdateFluid_SSPRK: ', wTime_UF,  ' s'
      END IF

    END IF

    IF( DEBUG ) wTime_UF = MPI_WTIME()
    IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL UpdateFluid_SSPRK'
    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )
    IF( DEBUG ) wTime_UF = MPI_WTIME() - wTime_UF

    IF( DEBUG )THEN
      WRITE(*,'(A)') 'AD: CALL ComputeFromConserved_GR'
      CALL ComputeFromConserved_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      Vmax = MAXVAL( ABS( &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),2)))
      LorentzFactor = One / SQRT( One - Vmax**2 )

      IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN
        WRITE(*,'(A8,A4,F12.10)') '', 'V = ', Vmax
        WRITE(*,'(A8,A4,F10.5)')  '', 'W = ', LorentzFactor
        WRITE(*,*)
      END IF
   END IF

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

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL ComputeFromConserved_GR'
      CALL ComputeFromConserved_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL WriteFieldsHDF'
      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL ComputeTally_Euler_GR'
      CALL ComputeTally_Euler_GR &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_GR &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler_GR &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL FinalizePositivityLimiter_Euler_GR

  CALL FinalizeSlopeLimiter_Euler_GR

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX 

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
