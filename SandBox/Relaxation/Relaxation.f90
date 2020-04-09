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
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_Initialize, &
    Timer_InputOutput, &
    Timer_Evolve
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
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    rhsCF, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_T, iAF_Ye, iAF_E, iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE, &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeElectronChemicalPotential_Table, &
    ComputeProtonChemicalPotential_Table, &
    ComputeNeutronChemicalPotential_Table
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE InitializationModule, ONLY: &
    InitializeFields_Relaxation
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_New, &
    InitializeNonlinearSolverTally, &
    FinalizeNonlinearSolverTally, &
    WriteNonlinearSolverTally

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: Configuration
  INTEGER  :: iCycle, iCycleD, iCycleW
  INTEGER  :: nE, nX(3), nNodes, nSpecies
  REAL(DP) :: t, dt, dt_0, t_end, wTime
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

  Configuration = 1

  SELECT CASE( Configuration )

    CASE ( 0 )
      ! from new profile, Nested_Newton fails
      D_0 = 3.999d14 * Gram / Centimeter**3
      T_0 = 1.356d1 * MeV
      Y_0 = 0.2783_DP

    CASE( 1 )

      !D_0 = 1.232d+14 * Gram / Centimeter**3
      !T_0 = 2.508d+11 * Kelvin
      !Y_0 = 0.2747_DP

      ! from new profile, right above D = 1e14
      D_0 = 1.084d14 * Gram / Centimeter**3
      T_0 = 1.590d1 * MeV
      Y_0 = 0.2728_DP

    CASE( 2 )

      !D_0 = 1.0520d+12 * Gram / Centimeter**3
      !T_0 = 8.9670d+10 * Kelvin
      !Y_0 = 0.1352_DP

      ! from new profile, right above D = 1e12
      D_0 = 1.032d12 * Gram / Centimeter**3
      T_0 = 7.588d0 * MeV
      Y_0 = 0.1347_DP

    CASE( 3 )

      !D_0 = 6.233d09 * Gram / Centimeter**3
      !T_0 = 3.021d10 * Kelvin
      !Y_0 = 0.3178_DP

      ! from new profile, right above D = 1e10
      D_0 = 1.104d10 * Gram / Centimeter**3
      T_0 = 2.054d0 * MeV
      Y_0 = 0.2818_DP

    CASE ( 4 )
      ! from new profile, higher Antineutrino density (max = 0.1)
      D_0 = 2.564d13 * Gram / Centimeter**3
      T_0 = 1.978d1 * MeV
      Y_0 = 0.2876_DP

  END SELECT

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'Relaxation'
  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2)') '', 'D [g/cm3] = ', D_0 / ( Gram / Centimeter**3 )
  WRITE(*,'(A6,A,ES8.2E2)') '', 'T     [K] = ', T_0 / Kelvin
  WRITE(*,'(A6,A,ES8.2E2)') '', 'Y         = ', Y_0
  WRITE(*,*)


  dt_0    = 1.0d-2 * Millisecond
  t       = 0.0d-0 * Millisecond
  t_end   = 5.0d-1 * Millisecond
  iCycleD = 1
  iCycleW = 1

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
             = 1.266038160710160_DP, &
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

   ! --- Initialize Timers ---

   CALL InitializeTimers

   CALL TimersStart( Timer_Initialize )

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

  CALL ComputeFromConserved_Fluid

  CALL ComputeFromConserved_Radiation

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_Initialize )

  CALL InitializeNonlinearSolverTally

  ! --- Evolve ---

  wTime = MPI_WTIME( )
  CALL TimersStart( Timer_Evolve )

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

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A,ES12.6E2,A2,A,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't [ms] = ',  t / Millisecond, &
          '', 'dt [ms] = ', dt / Millisecond

    END IF

    IF( dt > 0.0_DP )THEN

      CALL ComputeIncrement_TwoMoment_Implicit_New &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, uGE, uGF, &
               uCF, rhsCF, uCR, rhsCR )

      uCF = uCF + dt * rhsCF

      uCR = uCR + dt * rhsCR

      t = t + dt

    END IF

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_Fluid

      CALL ComputeFromConserved_Radiation

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL WriteNonlinearSolverTally( t )

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL ComputeFromConserved_Fluid

  CALL ComputeFromConserved_Radiation

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime
  CALL TimersStop( Timer_Evolve )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizeTimers

  CALL FinalizeNonlinearSolverTally

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  CALL FinalizeProgram

CONTAINS


  SUBROUTINE ComputeFromConserved_Fluid

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputePrimitive_Euler &
             ( uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye) )

      CALL ComputeElectronChemicalPotential_Table &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_Me) )

      CALL ComputeProtonChemicalPotential_Table &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_Mp) )

      CALL ComputeNeutronChemicalPotential_Table &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_Mn) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_Fluid


  SUBROUTINE ComputeFromConserved_Radiation

  END SUBROUTINE ComputeFromConserved_Radiation


END PROGRAM Relaxation
