PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, SqrtTiny, One, Zero
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  Use FluidFieldsModule, ONLY: &
    uPF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE TwoMoment_DiscretizationModule_Collisions_FMC, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    uCM, uPM, uAM, uGM
  USE TwoMoment_InputOutputModule_FMC, ONLY: &
    WriteTwoMomentFieldsHDF
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    ComputeFromConserved_TwoMoment_FMC, &
    HeatFluxTensorComponents_uuu, &
    ComputeHeatFluxTensorComponents_ddd_Lagrangian, &
    ComputeHeatFluxTensorComponents_uud_Lagrangian, &
    Flux_X1, Flux_X2, Flux_X3, ComputeTimeStep_TwoMoment, &
    ComputeTimeStep_TwoMoment_Realizable
  USE TwoMoment_OpacityModule_FMC, ONLY: &
    SetOpacities
  USE TwoMoment_TimeSteppingModule_FMC, ONLY: &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(2)  :: Direction
  CHARACTER(32) :: Spectrum
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem = 'CARTESIAN'
  CHARACTER(32) :: TimeSteppingScheme
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseEnergyLimiter
  INTEGER       :: nNodes
  INTEGER       :: nSpecies = 1
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t, t_end, dt, dt_CFL, V_0(3), CFL
  REAL(DP)      :: J_0, Chi, Sigma
  REAL(DP) :: l_uuu_munurho(0:3,0:3,0:3), l_ddd_ijk(1:3,1:3,1:3), l_uud_munurho(0:3,0:3,0:3)

  ProgramName = 'StreamingDopplerShift'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 32, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 1

      TimeSteppingScheme = 'SSPRK2'
      ! TimeSteppingScheme = 'IMEX_PDARS'

      t_end   = 1.0d-0
      iCycleD = 1
      iCycleW = 100
      maxCycles = 10000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      J_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

    CASE( 'SineWaveDiffusion' )

      nX  = [ 32, 1, 1 ]
      xL  = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ + 3.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 1

      TimeSteppingScheme = 'SSPRK2'

      !t_end   = 7.2951_DP !sig small old
      !t_end   = 1459.02504445_DP !sig large old

      t_end   = 60.0_DP !sig small
      !t_end   = 1500.0_DP !sig large

      iCycleD = 1000
      iCycleW = 1000
      maxCycles = 1000000

      V_0 = [ 0.1_DP, 0.0_DP, 0.0_DP ]

      J_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 2.6666666_DP 
      !Sigma = 533.33333333_DP

    CASE( 'StreamingDopplerShift' )

      Spectrum = 'Fermi-Dirac'

      Direction = 'X' ! --- (X,Y, or Z)

      IF(     TRIM( Direction ) .EQ. 'X' )THEN

        nX  = [ 128, 1, 1 ]
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
      iCycleW = 200
      maxCycles = 15000000

      J_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

      ! UseSlopeLimiter      = .FALSE.
      UsePositivityLimiter = .TRUE.
      UseEnergyLimiter     = .TRUE.

      ! UseRealizabilityTimeStep = .TRUE.

  CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  ! --- Auxiliary Initialization ---

  CALL InitializeDriver

  CALL SetOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, J_0, Chi, Sigma )

  ! --- Set Initial Condition ---

  CALL InitializeFields( V_0 , Direction, Spectrum )

  t = 0.0_DP

  ! --- Write Initial Condition ---

  CALL ComputeFromConserved_TwoMoment_FMC &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uPF, uCM, uPM, uAM, uGM )

  CALL WriteTwoMomentFieldsHDF( t )

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t, ' to t = ', t_end
  WRITE(*,*)

  iCycle = 0
  CFL = 0.3_DP
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    ! --- Compute Timestep ---

    ! CALL ComputeTimeStep_TwoMoment &
    !        ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, CFL, dt )

    CALL ComputeTimeStep_TwoMoment_Realizable &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uPF, CFL, dt )

    IF ( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    ! --- IMEX updating ---

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uPF, uCM, ComputeIncrement_TwoMoment_Implicit )

    t = t + dt

    ! --- Write updated values ---

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_TwoMoment_FMC &
              ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uPF, uCM, uPM, uAM, uGM )

      CALL WriteTwoMomentFieldsHDF( t )

    END IF

  END DO

  CALL ComputeFromConserved_TwoMoment_FMC &
  ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uPF, uCM, uPM, uAM, uGM )

  CALL WriteTwoMomentFieldsHDF( t )

  CALL FinalizeDriver

  ! --- Very Litte (If Any) Code Will Go Here ---

CONTAINS


  SUBROUTINE InitializeDriver

    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
    USE TwoMoment_FieldsModule_FMC, ONLY: &
      CreateTwoMomentFields
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
    USE TwoMoment_OpacityModule_FMC, ONLY: &
      CreateOpacities
    USE TwoMoment_PositivityLimiterModule_FMC, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE TwoMoment_TimeSteppingModule_FMC, ONLY: &
      Initialize_IMEX_RK

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

    CALL CreateTwoMomentFields( nX, [ 1, 1, 1 ], nE, 1, nSpecies, .TRUE. )

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

    ! --- Initialize Opacities ---

    CALL CreateOpacities &
           ( nx, [1, 1, 1], nE, 1, Verbose_Option = .TRUE.)

    ! --- Initialize Positivity Limiter ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = SqrtTiny, &
             Min_2_Option &
               = SqrtTiny, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE TwoMoment_FieldsModule_FMC, ONLY: &
      DestroyTwoMomentFields
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE TwoMoment_TimeSteppingModule_FMC, ONLY: &
      Finalize_IMEX_RK
    USE TwoMoment_OpacityModule_FMC, ONLY: &
      DestroyOpacities
    USE TwoMoment_PositivityLimiterModule_FMC, ONLY: &
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
    
    CALL Finalize_IMEX_RK

    CALL DestroyOpacities

    CALL FinalizePositivityLimiter_TwoMoment

    CALL DestroyTwoMomentFields

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange


    CALL FinalizeProgram

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver
