PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE UnitsModule, ONLY: &
    Second
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
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
    USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP)            :: wTime
  REAL(DP), PARAMETER :: Gamma = 4.0_DP / 3.0_DP
  INTEGER             :: iCycle, iCycleD, iCycleW
  REAL(DP)            :: t, dt, t_end

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StandingAccretionShock', &
           nX_Option &
             = [ 128, 32, 1 ], &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.2d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 2.0d0, Pi,    TwoPi ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           nStages_SSP_RK_Option &
             = 1, &  
           BasicInitialization_Option &
             = .TRUE. )
 
  t_end = 1.0d-3 * Second
  iCycleD = 10
  iCycleW = 10

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeFluid_SSPRK( nStages = 1 )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeFields &
         ( 4.0_DP * Pi, 0.5_DP, 1.0_DP, Gamma, 3.0d2 ) 

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

!!$  ! --- Main Part of Code Will Go Here

  iCycle = 0

  DO WHILE ( t < t_end )

    IF( t + dt < t_end )THEN
       t = t + dt
    ELSE
       dt = t_end - t
       t = t_end
    END IF

    iCycle = iCycle + 1

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

!    CALL UpdateFluid_SSPRK &
!           ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )

    ! --- Update primitive fluid variables, pressure, and sound speed
!    CALL ComputeFromConserved( iX_B0, iX_E0, uGF, uCF, uPF, uAF )

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

    END IF

  END DO


  ! --- Finalize ---

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
