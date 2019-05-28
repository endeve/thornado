PROGRAM Relaxation

  USE KindModule, ONLY: &
    DP, Third
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    MeV, &
    Kelvin
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
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    uCF, rhsCF, &
    uPF, iPF_D, &
    uAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE InitializationModule, ONLY: &
    InitializeFields_Relaxation
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_New

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL  :: wrt
  INTEGER  :: iCycle, iCycleD
  INTEGER  :: nE, nX(3), nNodes, nSpecies
  REAL(DP) :: t, dt, dt_0, t_end, dt_wrt, t_wrt, wTime
  REAL(DP) :: eL, eR
  REAL(DP) :: xL(3), xR(3)
  REAL(DP) :: D_0, T_0, Y_0

  nNodes   = 2
  nSpecies = 2

  nX = [ 1, 1, 1 ]
  xL = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
  xR = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer

  nE = 16
  eL = 0.0d0 * MeV
  eR = 3.0d2 * MeV

!!$  D_0 = 6.233d09 * Gram / Centimeter**3
!!$  T_0 = 3.021d10 * Kelvin
!!$  Y_0 = 0.3178_DP

  D_0 = 1.0520d+12 * Gram / Centimeter**3
  T_0 = 8.9670d+10 * Kelvin
  Y_0 = 0.1352_DP

  dt_0    = 1.0d-3 * Millisecond
  t       = 0.0d-0 * Millisecond
  t_end   = 1.0d+0 * Millisecond
  dt_wrt  = 1.0d-2 * Millisecond

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'Relaxation', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 0, 0, 0 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
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
           ZoomE_Option &
             = 1.2660_DP, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
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
         ( OpacityTableName_EmAb_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &
           OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &
           OpacityTableName_Pair_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
           Verbose_Option = .TRUE. )

  ! --- Set Initial Condition ---

  CALL InitializeFields_Relaxation( D_0, T_0, Y_0 )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t_wrt   = dt_wrt
  wrt     = .FALSE.
  iCycleD = 1

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A,ES8.2E2)') &
    '', 'Evolving from t [ms] = ', t / Millisecond, &
    ' to t [ms] = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    dt = dt_0

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A,ES12.6E2,A2,A,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't [ms] = ',  t / Millisecond, &
          '', 'dt [ms] = ', dt / Millisecond

    END IF

    IF( dt > 0.0_DP )THEN

      CALL ComputeIncrement_TwoMoment_Implicit_New &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, uGE, uGF, uCF, rhsCF, uCR, rhsCR )

      uCF = uCF + dt * rhsCF

      uCR = uCR + dt * rhsCR

    END IF

    t = t + dt

    IF( wrt )THEN

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

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

  CALL FinalizeProgram

END PROGRAM Relaxation
