PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Half, One, Two, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    Gram, Erg
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    InitializeGravitySolver_Newtonian_Poseidon, &
    FinalizeGravitySolver_Newtonian_Poseidon, &
    SolveGravity_Newtonian_Poseidon
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, uDF, &
    iCF_D
  USE Euler_SlopeLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    ApplySlopeLimiter_Euler_NonRelativistic_IDEAL
  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic, &
    ComputeTimeStep_Euler_NonRelativistic
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE Euler_TallyModule_NonRelativistic, ONLY: &
    InitializeTally_Euler_NonRelativistic, &
    FinalizeTally_Euler_NonRelativistic, &
    ComputeTally_Euler_NonRelativistic

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(64), PARAMETER :: FileName = 'YahilHomologousCollapse_Gm_130.dat'
  REAL(DP),      PARAMETER :: Gamma           = 1.3_DP
  REAL(DP),      PARAMETER :: CollapseTime    = 150_DP * Millisecond
  REAL(DP),      PARAMETER :: CentralDensity  = 7.0d9 * Gram / Centimeter**3
  REAL(DP),      PARAMETER :: CentralPressure = 6.0d27 * Erg / Centimeter**3
  REAL(DP),      PARAMETER :: CoreRadius      = 1.0d5 * Kilometer

  LOGICAL  :: wrt
  INTEGER  :: iCycle, iCycleD
  INTEGER  :: nX(3), nNodes
  REAL(DP) :: xL(3), xR(3)
  REAL(DP) :: t, dt, t_end, dt_wrt, t_wrt, wTime

  nX     = [ 128, 1, 1 ]
  nNodes = 3
  xL     = [ 0.0d0,      0.0d0, 0.0d0 ]
  xR     = [ CoreRadius, Pi,    TwoPi ]

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'YahilLattimerCollapse', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 30, 3, 1 ], &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           zoomX_Option &
             = [ 1.071835456828339_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  t_end   = CollapseTime - 0.5_DP * Millisecond
  dt_wrt  = 0.1_DP * Millisecond
  iCycleD = 100

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( BetaTVD_Option = 1.15_DP, &
           BetaTVB_Option = 0.00_DP, &
           SlopeTolerance_Option = 1.0d-6, &
           UseSlopeLimiter_Option = .FALSE., &
           UseTroubledCellIndicator_Option = .FALSE., &
           LimiterThresholdParameter_Option = 0.12_DP )

  CALL InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12, &
           UsePositivityLimiter_Option = .FALSE. )

  CALL InitializeFields &
         ( FileName, Gamma, CollapseTime, CentralDensity, CentralPressure )

  CALL ApplySlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL InitializeGravitySolver_Newtonian_Poseidon

  CALL SolveGravity_Newtonian_Poseidon &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,iCF_D) )

!  CALL ComputeGravitationalPotential &
!         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 3 )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeTally_Euler_NonRelativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .TRUE. )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
             CFL = 0.3_DP / ( Two * DBLE( nNodes - 1 ) + One ), &
             TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A9,ES13.6E3,A1,A10,ES13.6E3)') &
          '', 'Cycle = ', iCycle, &
          '', 't [ms] = ',  t / Millisecond, &
          '', 'dt [ms] = ', dt / Millisecond

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, &
             ComputeIncrement_Euler_DG_Explicit, &
             SolveGravity_Newtonian_Poseidon )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
               Verbose_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
           Verbose_Option = .TRUE. )

  CALL FinalizeTally_Euler_NonRelativistic

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeGravitySolver_Newtonian_Poseidon

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
