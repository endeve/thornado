PROGRAM StreamingSineWave

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    nDOF
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
  USE ClosureModule_M1, ONLY: &
    InitializeMomentClosure
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE RadiationFieldsModule, ONLY: &
    uCR, iCR_N
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE InitializationModule, ONLY: &
    InitializeFields
  USE dgDiscretizationModule, ONLY: &
    ComputeIncrement_M1_DG_Explicit
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
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER       :: nE, nX(3), nNodes
  REAL(DP)      :: t, dt, t_end, wTime, Error
  REAL(DP)      :: xL(3), xR(3)
  REAL(DP)      :: eL,    eR
  REAL(DP)      :: N0, SigmaA, SigmaS
  
  REAL(DP), ALLOCATABLE :: uCR_N(:,:,:,:,:)

  ProgramName = 'SineWaveStreaming'
  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      Direction = 'X'

      nX = [ 16, 16, 16 ]
      nE = 20

      nNodes = 2

      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      eL = 0.0_DP
      eR = 1.0_DP

      TimeSteppingScheme = 'SSPRK2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      t_end     = 1.0d+0
      iCycleD   = 10
      iCycleW   = 10
      maxCycles = 10000

    CASE( 'SineWaveDamping' )

      ! --- Minerbo Closure Only ---

      nX = [ 16, 1, 1 ]
      nE = 1

      nNodes = 3

      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      eL = 0.0_DP
      eR = 1.0_DP

      TimeSteppingScheme = 'IMEX_P_A2'

      N0     = 0.0_DP
      SigmaA = 1.0_DP
      SigmaS = 0.0_DP

      t_end     = 1.0d+1
      iCycleD   = 10
      iCycleW   = 100
      maxCycles = 10000

    CASE( 'SineWaveDiffusion' )

      nX = [ 16, 1, 1 ]
      nE = 1

      nNodes = 3

      xL = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ + 3.0_DP, 1.0_DP, 1.0_DP ]

      eL = 0.0_DP
      eR = 1.0_DP

      TimeSteppingScheme = 'IMEX_P_A2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 1.0d+2

      t_end     = 1.0d+2
      iCycleD   = 10
      iCycleW   = 1000
      maxCycles = 100000

    CASE( 'LineSource' )

      nX = [ 128, 128, 1 ]
      nE = 1

      nNodes = 1

      xL = [ - 1.25_DP, - 1.25_DP, 0.0_DP ]
      xR = [ + 1.25_DP, + 1.25_DP, 1.0_DP ]

      eL = 0.0_DP
      eR = 1.0_DP

      TimeSteppingScheme = 'SSPRK2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      t_end     = 1.0d+0
      iCycleD   = 10
      iCycleW   = 10
      maxCycles = 10000

    CASE( 'HomogeneousSphere' )

  END SELECT

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
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
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           Opacity_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

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

  CALL InitializeMomentClosure

  ! --- Initialize Implicit Solver ---

  CALL InitializeCollisions &
         ( N0_Option = N0, SigmaA_Option = SigmaA, SigmaS_Option = SigmaS )

  ! --- Set Initial Condition ---

  CALL InitializeFields &
         ( Direction_Option = TRIM( Direction ) )

  ALLOCATE( uCR_N(nDOF,nE,nX(1),nX(2),nX(3)) )

  uCR_N = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N,1)

  ! --- Write Initial Condition ---
  ! --- Temporarily Removed Due to HDF5 ifort Incompatibility --- 
  !CALL WriteFieldsHDF( Time = 0.0_DP, WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t  = 0.0d-0
  dt = 0.2_DP / (2.0_DP*DBLE(nNodes-1)+1.0_DP) / DBLE(nX(1))

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
             ComputeIncrement_M1_DG_Explicit, &
             ComputeIncrement_M1_DG_Implicit, &
             ComputeCorrection_M1_DG_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN
      ! --- Temporarily Removed Due to HDF5 ifort Incompatibility --- 
      !CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

    END IF

  END DO

  ! --- Temporarily Removed Due to HDF5 ifort Incompatibility --- 
  !CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  Error = MAXVAL( ABS( uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N,1) - uCR_N ) )

  WRITE(*,'(A8,A,ES16.10E2)') '', 'Error = ', Error
  WRITE(*,*)

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL Finalize_IMEX_RK

  CALL FinalizeCollisions

  CALL FinalizeProgram

END PROGRAM StreamingSineWave
