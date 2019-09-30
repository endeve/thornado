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
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    CreateOpacities, &
    SetConstantOpacities, &
    DestroyOpacities
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: TimeSteppingScheme
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  REAL(DP)      :: eL, eR, xL(3), xR(3)
  REAL(DP)      :: t, dt, t_end, V_0(3)
  REAL(DP)      :: D_0, Chi, Sigma

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'SineWaveDiffusion'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      nX  = [ 16, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]
      bcX = [ 1, 0, 0 ]

      nE  = 1
      eL  = 0.0_DP
      eR  = 1.0_DP
      bcE = 0

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      t_end   = 1.0d0
      iCycleD = 1
      iCycleW = 1
      maxCycles = 10000

      V_0 = [ 0.3_DP, 0.0_DP, 0.0_DP ]

      D_0   = 0.0_DP
      Chi   = 0.0_DP
      Sigma = 0.0_DP

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

  CALL SetConstantOpacities( D_0, Chi, Sigma )

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  ! --- Set Initial Condition ---

  CALL InitializeFields( V_0 )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  t = 0.0_DP
  dt = 0.5_DP * MINVAL( (xR-xL) / DBLE(nX) ) &
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

  CALL DestroyOpacities

  CALL Finalize_IMEX_RK

  CALL FinalizeEquationOfState

  CALL FinalizeTimers

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
