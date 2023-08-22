PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Zero
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
    Flux_X1, Flux_X2, Flux_X3, ComputeTimeStep_TwoMoment
  USE TwoMoment_TimeSteppingModule_FMC, ONLY: &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem = 'CARTESIAN'
  CHARACTER(32) :: TimeSteppingScheme
  INTEGER       :: nNodes
  INTEGER       :: nSpecies = 1
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t, t_end, dt, dt_CFL, V_0(3), CFL
  REAL(DP)      :: J_0, Chi, Sigma
  REAL(DP) :: l_uuu_munurho(0:3,0:3,0:3), l_ddd_ijk(1:3,1:3,1:3), l_uud_munurho(0:3,0:3,0:3)

  ProgramName = 'SineWaveStreaming'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 64, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 1, 1 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 1

      nNodes = 2

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

  CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  CALL InitializeDriver

  CALL InitializeFields( V_0 )

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
  CFL = 0.1_DP
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    ! --- Compute Timestep ---
    CALL ComputeTimeStep_TwoMoment &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, CFL, dt )

    IF ( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    ! --- IMEX updating ---

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCM, ComputeIncrement_TwoMoment_Implicit )

    t = t + dt

    ! --- Write updated values ---

    CALL ComputeFromConserved_TwoMoment_FMC &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uPF, uCM, uPM, uAM, uGM )

    CALL WriteTwoMomentFieldsHDF( t )

  END DO

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
    
    CALL Finalize_IMEX_RK

    CALL DestroyTwoMomentFields

    CALL FinalizeProgram

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver
