PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
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
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_MHD_Relativistic
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
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeMagnetofluid_SSPRK, &
    FinalizeMagnetofluid_SSPRK, &
    UpdateMagnetofluid_SSPRK 

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: AdvectionProfile
  CHARACTER(32) :: CoordinateSystem
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, CFL
  REAL(DP)      :: ZoomX(3)

  LOGICAL  :: WriteGF = .TRUE., WriteMF = .TRUE.
  LOGICAL  :: ActivateUnits = .FALSE.
  LOGICAL  :: UseMHD = .TRUE.

  ProgramName = 'Advection'

  swX               = [ 0, 0, 0 ]
  RestartFileNumber = -1
  t                 = 0.0_DP
  ZoomX             = 1.0_DP

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'SineWave'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 500, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 5.0_DP, 1.0_DP, 1.0_DP ]

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  Advection'
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

  CALL InitializeMagnetofluid_SSPRK( nStages = nStagesSSPRK )
  WRITE(*,*)
  WRITE(*,'(A6,A,ES11.3E3)') '', 'CFL: ', CFL

  uCM = Zero ! Without this, crashes when copying data in TimeStepper
  uDM = Zero ! Without this, crashes in IO

  PRINT*, 'Initializing.'

  CALL InitializeFields_Relativistic_MHD &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ) )

  PRINT*, 'Computing from conserved.'

  CALL ComputeFromConserved_MHD_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uPM, uAM )

  PRINT*, 'Writing.'

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteMF_Option = WriteMF )

  CALL FinalizeMagnetofluid_SSPRK

  CALL FinalizeEquationOfState

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
