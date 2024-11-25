PROGRAM DeleptonizationWave

  USE KindModule, ONLY: &
    DP, SqrtTiny, Third, Pi, TwoPi
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
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, &
    uAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR, nSpecies
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
  USE TimeSteppingModule_IMEX_RK, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE dgDiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_M1_DG_Implicit, &
    ComputeCorrection_M1_DG_Implicit
  USE ProgenitorModule

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL  :: TEST_DEBUG
  LOGICAL  :: UsePositivityLimiter
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  INTEGER  :: iCycle, iCycleD, iCycleW
  INTEGER  :: nE, nX(3), bcX(3), nNodes
  REAL(DP) :: t, dt, t_end, wTime
  REAL(DP) :: eL, eR
  REAL(DP) :: xL(3), xR(3), ZoomX(3)

  CoordinateSystem = 'SPHERICAL'

  TEST_DEBUG = .TRUE.

  SELECT CASE( CoordinateSystem )

    CASE( 'CARTESIAN' )

      ProgramName = 'DeleptonizationWave'

      nX = [  128, 128, 1 ]  ! 96 / 128 [96, 96, 1]
      xL = [  0.0_DP, 0.0_DP, - 1.0d2 ] * Kilometer
      xR = [  4.0d2,  4.0d2,  + 1.0d2 ] * Kilometer
      
      bcX = [ 32, 32, 32 ]

      zoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]
 
    CASE( 'SPHERICAL' )

      ProgramName = 'DeleptonizationWave_Spherical'

      nX = [ 181, 1, 1 ] ! 512
      xL = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ]
      xR = [ 5.66d2 * Kilometer, Pi,     TwoPi  ]
      !xR = [ 3.0d3 * Kilometer, Pi,     TwoPi  ]

      bcX = [ 32, 0, 0 ]

      zoomX = [ 1.007_DP, 1.0_DP, 1.0_DP ]

  END SELECT
 
  nNodes = 2

  UsePositivityLimiter = .TRUE.

  IF( nNodes == 1 ) UsePositivityLimiter = .FALSE.

  nE = 10
  eL = 0.0d0 * MeV
  eR = 300.0d0 * MeV

  CALL InitializeProgram &
         ( ProgramName_Option &
             = ProgramName, &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 01, 01, 01 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           ZoomE_Option &
             = 1.414649066384746_DP, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = CoordinateSystem, &
           ActivateUnits_Option &
             = .TRUE., &
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

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option = 'EquationOfStateTable.h5' )

  ! --- Initialize Opacities ---

  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_EmAb_Option = 'OpacityTable_EmAb.h5', &
           OpacityTableName_Iso_Option  = 'OpacityTable_Iso.h5',  &
           Verbose_Option = .TRUE. )

  ! --- Create Neutrino Opacities ---

  CALL CreateNeutrinoOpacities( nZ, nNodesZ, nSpecies )

  ! --- Initialize Slope Limiter ---

  CALL InitializeSlopeLimiter_TwoMoment                &
         ( BetaTVD_Option = 2.0_DP,                    &
           BetaTVB_Option = 0.0d0,                      &
           SlopeTolerance_Option = 1.0d-6,             &
           UseSlopeLimiter_Option = .FALSE.,            &
           UseCharacteristicLimiting_Option = .FALSE., &
           Verbose_Option = .TRUE. )

  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = 0.0d0 + SqrtTiny, &
           Max_1_Option = 1.0d0 - SqrtTiny, &
           Min_2_Option = 0.0d0 + SqrtTiny, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter )

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK &
         ( Scheme = 'IMEX_PDARS_3' )

  ! --- Set Initial Condition ---

  CALL InitializeFields
  !CALL InitializeFields( Profile_Option = 'ChimeraBounce_fined.d')

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  CALL ComputeNeutrinoOpacities &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
           uPF(:,:,:,:,iPF_D), &
           uAF(:,:,:,:,iAF_T), &
           uAF(:,:,:,:,iAF_Ye) )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP
  t_end = 5.0d1 * Millisecond

  IF ( CoordinateSystem == 'SPHERICAL') THEN
    ! Spherical Coordinate time step
    !!! NEED FIX for 2D
    dt  = 0.2_DP * MINVAL( MeshX(1) % Width(1:nX(1)) ) &
            / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )
    WRITE(*,*)
    WRITE(*,*) 'dr [km]', MINVAL( MeshX(1) % Width ) / Kilometer
    WRITE(*,*) 'dr [km]', MINLOC( MeshX(1) % Width )
    WRITE(*,*) 'dr [km]', MeshX(1) % Width(1),MeshX(1) % Width(nX(1))
    WRITE(*,*) 'para   ', 0.1_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )
    WRITE(*,*) 'dt [dimless]', dt
    WRITE(*,*) 'Millisecond', Millisecond
    WRITE(*,*)
  ELSE IF ( CoordinateSystem == 'CARTESIAN' ) THEN
    ! Cartesian Coordinate time step
    ! dt  = 0.49_DP * MINVAL( (xR-xL) / DBLE( nX ) ) &
    dt  = 0.2_DP * MINVAL( (xR-xL) / DBLE( nX ) ) &
            / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )
    WRITE(*,*)
    WRITE(*,*) 'dx [km]', MINVAL( (xR-xL) / DBLE( nX ) ) / Kilometer
    WRITE(*,*)
  END IF

  iCycleD = 10
  IF( TEST_DEBUG ) THEN
    iCycleW = INT(t_end / dt / 40)
    !t_end = 6602 * dt
  ELSE
    iCycleW = 200 ! 200 -> 128, 150 -> 96
  END IF

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A3,A8,ES8.2E2,A3)') &
    '', 'Evolving from t = ', t / Millisecond, 'ms', &
    ' to t = ', t_end / Millisecond, 'ms'
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF
    
    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t / Millisecond, &
          '', 'dt = ', dt / Millisecond

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCR, &
             ComputeIncrement_TwoMoment_Explicit, &
             ComputeIncrement_M1_DG_Implicit, &
             ComputeCorrection_M1_DG_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE., &
               WriteOP_Option = .TRUE. )

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  CALL DestroyNeutrinoOpacities

  CALL FinalizeSlopeLimiter_TwoMoment

  CALL FinalizePositivityLimiter_TwoMoment

  CALL Finalize_IMEX_RK

  CALL FinalizeProgram

END PROGRAM DeleptonizationWave
