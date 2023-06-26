#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, FourPi
  USE UnitsModule, ONLY: &
    SpeedOfLight, &
    PlanckConstant, &
    BoltzmannConstant, &
    AtomicMassUnit, &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV, Erg
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    nNodesX, &
    nNodesZ, &
    nDOFX,   &
    nDOFE,   &
    nDOF
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Implicit, &
    Timer_Im_In, &
    Timer_Im_MapForward, &
    Timer_Im_EosIn, &
    Timer_Im_Solve, &
    Timer_Im_ComputeOpacity, &
    Timer_Im_ComputeRate, &
    Timer_Im_ComputeLS, &
    Timer_Im_UpdateFP, &
    Timer_Im_CoupledAA, &
    Timer_Im_NestedAA, &
    Timer_Im_NestedNewton, &
    Timer_Im_Newton, &
    Timer_Im_NestInner, &
    Timer_Im_EmAb_FP, &
    Timer_Im_Presolve, &
    Timer_Im_Increment, &
    Timer_Im_EosOut, &
    Timer_Im_MapBackward, &
    Timer_Im_Out
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX3D, &
    WeightsX_q
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable, &
    NodeNumberTable4D
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_T, iAF_Ye, iAF_E
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies, iNuE, iNuE_Bar
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_Point, &
    ComputeEquilibriumDistributionAndDerivatives_Point, &
    ComputeNeutrinoOpacities_EC_Point, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Point, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacitiesAndDerivatives_NES_Point, &
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit
  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_New
  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_DGFV

  PUBLIC :: InitializeNonlinearSolverTally
  PUBLIC :: FinalizeNonlinearSolverTally
  PUBLIC :: WriteNonlinearSolverTally

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV
  REAL(DP), PARAMETER :: Unit_E = Erg / Gram

  LOGICAL, PARAMETER :: ReportConvergenceData = .FALSE.
  INTEGER  :: Iterations_Min
  INTEGER  :: Iterations_Max
  INTEGER  :: Iterations_Ave

  ! --- Solver Tally ---

  LOGICAL              :: TallyNonlinearSolver = .FALSE.
  INTEGER              :: TallyFileNumber
  INTEGER, ALLOCATABLE :: MinIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: MaxIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: AveIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: AveIterationsInner_K(:,:,:)

  LOGICAL, PARAMETER :: SolveMatter = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerEmAb = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerPair = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPairLagAllButJ0 = .FALSE.

  INTEGER  :: nE_G, nX_G
  INTEGER  :: iE_B0,    iE_E0
  INTEGER  :: iE_B1,    iE_E1
  INTEGER  :: iX_B0(3), iX_E0(3)
  INTEGER  :: iX_B1(3), iX_E1(3)
  REAL(DP) :: wTime
  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)     ! --- Spatial Geometry
  REAL(DP), ALLOCATABLE :: CF_N(:,:)     ! --- Conserved Fluid
  REAL(DP), ALLOCATABLE :: dF_N(:,:)     ! --- Conserved Fluid Increment
  REAL(DP), ALLOCATABLE :: CR_K(:,:,:)   ! --- Conserved Radiation (Element)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:) ! --- Conserved Radiation (Node)
  REAL(DP), ALLOCATABLE :: dR_N(:,:,:,:) ! --- Conserved Radiation Increment

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit_New &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iX1, iX2, iX3, iGF, iCF, iCR, iS, iNodeX, iE
    INTEGER  :: nIterations
    INTEGER  :: nIterations_Inner
    INTEGER  :: nIterations_Outer
    REAL(DP) :: CF_N(1:nDOFX,1:nCF)
    REAL(DP) :: PF_N(1:nDOFX,1:nPF)
    REAL(DP) :: AF_N(1:nDOFX,1:nAF)
    REAL(DP), ALLOCATABLE :: Kappa(:)
    REAL(DP), ALLOCATABLE :: Chi_T(:)
    REAL(DP), ALLOCATABLE :: Eta_T(:)
    REAL(DP), ALLOCATABLE :: Chi(:,:,:)
    REAL(DP), ALLOCATABLE :: fEQ(:,:,:)
    REAL(DP), ALLOCATABLE :: Sig(:,:,:)
    REAL(DP), ALLOCATABLE :: Chi_NES(:,:,:)
    REAL(DP), ALLOCATABLE :: Eta_NES(:,:,:)
    REAL(DP), ALLOCATABLE :: Chi_Pair(:,:,:)
    REAL(DP), ALLOCATABLE :: Eta_Pair(:,:,:)


    CALL TimersStart( Timer_Implicit )

    CALL TimersStart( Timer_Im_In )

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    CALL InitializeCollisions_New( iE_B0, iE_E0 )

    ALLOCATE( Kappa   (nE_G) )
    ALLOCATE( Chi_T   (nE_G) )
    ALLOCATE( Eta_T   (nE_G) )
    ALLOCATE( Chi     (nE_G,nSpecies,nDOFX) )
    ALLOCATE( fEQ     (nE_G,nSpecies,nDOFX) )
    ALLOCATE( Sig     (nE_G,nSpecies,nDOFX) )
    ALLOCATE( Chi_NES (nE_G,nSpecies,nDOFX) )
    ALLOCATE( Eta_NES (nE_G,nSpecies,nDOFX) )
    ALLOCATE( Chi_Pair(nE_G,nSpecies,nDOFX) )
    ALLOCATE( Eta_Pair(nE_G,nSpecies,nDOFX) )



    Iterations_Min = + HUGE( 1 )
    Iterations_Max = - HUGE( 1 )
    Iterations_Ave = 0

    IF( TallyNonlinearSolver )THEN

      MinIterations_K      = + HUGE( 1 )
      MaxIterations_K      = - HUGE( 1 )
      AveIterations_K      = 0
      AveIterationsInner_K = 0

    END IF

    CALL TimersStop( Timer_Im_In )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ! --- Copy inputs to local arrays ---

      CALL TimersStart( Timer_Im_MapForward )

      DO iCF = 1, nCF

        CF_N(:,            iCF) = U_F (:,iX1,iX2,iX3,iCF)
        dU_F(:,iX1,iX2,iX3,iCF) = CF_N(:,            iCF)

      END DO

      ! --- Rearrange data to group energy nodes together ---

      CALL MapForward_R_New &
             ( iE_B0, iE_E0, U_R(:,iE_B0:iE_E0,iX1,iX2,iX3,:,:), CR_N )

      !!! NEW: copy old moments to dR_N (will be used when assigning dR_N)
      CALL MapForward_R_New &
             ( iE_B0, iE_E0, U_R(:,iE_B0:iE_E0,iX1,iX2,iX3,:,:), dR_N )
      !!! END: copy old moments to dR_N (will be used when assigning dR_N)

      ! --- Compute Primitive Quantities ---

      CALL ComputePrimitive_Euler &
             ( CF_N(:,iCF_D ), CF_N(:,iCF_S1), CF_N(:,iCF_S2), &
               CF_N(:,iCF_S3), CF_N(:,iCF_E ), CF_N(:,iCF_Ne), &
               PF_N(:,iPF_D ), PF_N(:,iPF_V1), PF_N(:,iPF_V2), &
               PF_N(:,iPF_V3), PF_N(:,iPF_E ), PF_N(:,iPF_Ne), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL TimersStop( Timer_Im_MapForward )

      ! --- EOS Table Lookup ---

      CALL TimersStart( Timer_Im_EosIn )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), &
               AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye) )

      CALL TimersStop( Timer_Im_EosIn )

      CALL TimersStart( Timer_Im_Solve )

      ! --- Opacity Table Lookup ---

      CALL TimersStart( Timer_Im_ComputeOpacity )

      DO iS = 1, nSpecies

        CALL ComputeNeutrinoOpacities_EC_Points &
               ( 1, nE_G, 1, nDOFX, E_N, PF_N(:,iPF_D), &
                 AF_N(:,iAF_T), AF_N(:,iAF_Ye), iS, Chi(:,iS,:) )

      END DO

      DO iS = 1, nSpecies

        CALL ComputeNeutrinoOpacities_ES_Points &
               ( 1, nE_G, 1, nDOFX, E_N, PF_N(:,iPF_D), &
                 AF_N(:,iAF_T), AF_N(:,iAF_Ye), iS, 1, Sig(:,iS,:) )

        DO iX1 = 1, nDOFX

          Sig(1:nE_G,iS,iX1) &
            = Sig(1:nE_G,iS,iX1) * FourPi * E_N(1:nE_G) * E_N(1:nE_G)

        END DO

      END DO

      Chi_NES  = Zero
      Eta_NES  = Zero
      Chi_Pair = Zero
      Eta_Pair = Zero

      CALL TimersStop( Timer_Im_ComputeOpacity )

      IF( nSpecies .EQ. 1 )THEN

        ! --- Single Species (Electron Neutrinos) ---

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_EmAb_NuE &
                 ( CR_N    (:,iCR_N,iNuE,iNodeX), &
                   dt * Chi(:,      iNuE,iNodeX), &
                   fEQ     (:,      iNuE,iNodeX), &
                   PF_N(iNodeX,iPF_D ), AF_N(iNodeX,iAF_T), &
                   AF_N(iNodeX,iAF_Ye), AF_N(iNodeX,iAF_E) )

        END DO

      ELSE

        ! --- Electron Neutrinos and Antineutrinos ---

#ifdef NEUTRINO_MATTER_SOLVER_EMAB

        CALL TimersStart( Timer_Im_EmAb_FP )

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_EmAb_FP &
                  ( dt, iNuE, iNuE_Bar, &
                    CR_N(:,iCR_N,1:2,iNodeX), &
                    Chi (:,      1:2,iNodeX), &
                    fEQ (:,      1:2,iNodeX), &
                    PF_N(iNodeX,iPF_D), &
                    AF_N(iNodeX,iAF_T), &
                    AF_N(iNodeX,iAF_Ye), &
                    AF_N(iNodeX,iAF_E), &
                    nIterations_Out = nIterations )

          Iterations_Min = MIN( Iterations_Min, nIterations )
          Iterations_Max = MAX( Iterations_Max, nIterations )
          Iterations_Ave = Iterations_Ave + nIterations

        END DO

        CALL TimersStop( Timer_Im_EmAb_FP )

#elif NEUTRINO_MATTER_SOLVER_FIXED_POINT_COUPLED

        CALL TimersStart( Timer_Im_CoupledAA )


        DO iNodeX = 1, nDOFX

          ! PRINT*, "  INPUT MATTER STATE "
          ! PRINT*, "  D = ", PF_N(iNodeX,iPF_D)
          ! PRINT*, "  T = ", AF_N(iNodeX,iAF_T)
          ! PRINT*, "  Ye = ", AF_N(iNodeX,iAF_Ye)
          ! PRINT*, "  E = ", AF_N(iNodeX,iAF_E)
          ! PRINT*
          ! PRINT*, "  INPUT RADIATION FIELD "
          ! PRINT*, "  JNein = ", CR_N(:,iCR_N,1,iNodeX)
          ! PRINT*, "  JANein = ", CR_N(:,iCR_N,2,iNodeX)
          ! PRINT*


          CALL SolveMatterEquations_FP_Coupled &
                 ( dt, iNuE, iNuE_Bar, &
                   CR_N(:,iCR_N,1:2,iNodeX), &
                   Chi (:,      1:2,iNodeX), &
                   fEQ (:,      1:2,iNodeX), &
                   Chi_NES (:,  1:2,iNodeX), &
                   Eta_NES (:,  1:2,iNodeX), &
                   Chi_Pair(:,  1:2,iNodeX), &
                   Eta_Pair(:,  1:2,iNodeX), &
                   PF_N(iNodeX,iPF_D),  &
                   AF_N(iNodeX,iAF_T),  &
                   AF_N(iNodeX,iAF_Ye), &
                   AF_N(iNodeX,iAF_E), &
                   nIterations )

           ! PRINT*, "  OUTPUT MATTER STATE "
           ! PRINT*, "  D = ", PF_N(iNodeX,iPF_D)
           ! PRINT*, "  T = ", AF_N(iNodeX,iAF_T)
           ! PRINT*, "  Ye = ", AF_N(iNodeX,iAF_Ye)
           ! PRINT*, "  E = ", AF_N(iNodeX,iAF_E)
           ! PRINT*
           ! PRINT*, "  OUTPUT RADIATION FIELD "
           ! PRINT*, "  JNeout = ", CR_N(:,iCR_N,1,iNodeX)
           ! PRINT*, "  JANeout = ", CR_N(:,iCR_N,2,iNodeX)
           ! PRINT*

          IF( TallyNonlinearSolver )THEN

            MinIterations_K(iX1,iX2,iX2) &
              = MIN( nIterations, MinIterations_K(iX1,iX2,iX2) )
            MaxIterations_K(iX1,iX2,iX3) &
              = MAX( nIterations, MaxIterations_K(iX1,iX2,iX3) )
            AveIterations_K(iX1,iX2,iX3) &
              = AveIterations_K(iX1,iX2,iX3) + nIterations

          END IF

        END DO

        CALL TimersStop( Timer_Im_CoupledAA )

#elif NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_AA

        CALL TimersStart( Timer_Im_NestedAA )

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_FP_NestedAA &
                 ( dt, iNuE, iNuE_Bar, &
                   CR_N(:,iCR_N,1:2,iNodeX), &
                   Chi (:,      1:2,iNodeX), &
                   fEQ (:,      1:2,iNodeX), &
                   Chi_NES (:,  1:2,iNodeX), &
                   Eta_NES (:,  1:2,iNodeX), &
                   Chi_Pair(:,  1:2,iNodeX), &
                   Eta_Pair(:,  1:2,iNodeX), &
                   PF_N(iNodeX,iPF_D),  &
                   AF_N(iNodeX,iAF_T),  &
                   AF_N(iNodeX,iAF_Ye), &
                   AF_N(iNodeX,iAF_E), &
                   nIterations_Inner, &
                   nIterations_Outer )

          IF( TallyNonlinearSolver )THEN

            MinIterations_K(iX1,iX2,iX2) &
              = MIN( nIterations_Outer, MinIterations_K(iX1,iX2,iX2) )
            MaxIterations_K(iX1,iX2,iX3) &
              = MAX( nIterations_Outer, MaxIterations_K(iX1,iX2,iX3) )
            AveIterations_K(iX1,iX2,iX3) &
              = AveIterations_K(iX1,iX2,iX3)      + nIterations_Outer
            AveIterationsInner_K(iX1,iX2,iX3) &
              = AveIterationsInner_K(iX1,iX2,iX3) + nIterations_Inner

          END IF

        END DO

        CALL TimersStop( Timer_Im_NestedAA )

#elif NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON

        CALL TimersStart( Timer_Im_NestedNewton )

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_FP_NestedNewton &
                 ( dt, iNuE, iNuE_Bar, &
                   CR_N(:,iCR_N,1:2,iNodeX), &
                   Chi (:,      1:2,iNodeX), &
                   fEQ (:,      1:2,iNodeX), &
                   Chi_NES (:,  1:2,iNodeX), &
                   Eta_NES (:,  1:2,iNodeX), &
                   Chi_Pair(:,  1:2,iNodeX), &
                   Eta_Pair(:,  1:2,iNodeX), &
                   PF_N(iNodeX,iPF_D),  &
                   AF_N(iNodeX,iAF_T),  &
                   AF_N(iNodeX,iAF_Ye), &
                   AF_N(iNodeX,iAF_E), &
                   nIterations_Inner, &
                   nIterations_Outer )

          IF( TallyNonlinearSolver )THEN

            MinIterations_K(iX1,iX2,iX2) &
              = MIN( nIterations_Outer, MinIterations_K(iX1,iX2,iX2) )
            MaxIterations_K(iX1,iX2,iX3) &
              = MAX( nIterations_Outer, MaxIterations_K(iX1,iX2,iX3) )
            AveIterations_K(iX1,iX2,iX3) &
              = AveIterations_K(iX1,iX2,iX3)      + nIterations_Outer
            AveIterationsInner_K(iX1,iX2,iX3) &
              = AveIterationsInner_K(iX1,iX2,iX3) + nIterations_Inner

          END IF

        END DO

        CALL TimersStop( Timer_Im_NestedNewton )

#elif NEUTRINO_MATTER_SOLVER_NEWTON

        CALL TimersStart( Timer_Im_Newton )

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_Newton &
                 ( dt, iNuE, iNuE_Bar, &
                   CR_N(:,iCR_N,1:2,iNodeX), &
                   Chi (:,      1:2,iNodeX), &
                   fEQ (:,      1:2,iNodeX), &
                   Chi_NES (:,  1:2,iNodeX), &
                   Eta_NES (:,  1:2,iNodeX), &
                   Chi_Pair(:,  1:2,iNodeX), &
                   Eta_Pair(:,  1:2,iNodeX), &
                   PF_N(iNodeX,iPF_D),  &
                   AF_N(iNodeX,iAF_T),  &
                   AF_N(iNodeX,iAF_Ye), &
                   AF_N(iNodeX,iAF_E), &
                   nIterations )

          IF( TallyNonlinearSolver )THEN

            MinIterations_K(iX1,iX2,iX2) &
              = MIN( nIterations, MinIterations_K(iX1,iX2,iX2) )
            MaxIterations_K(iX1,iX2,iX3) &
              = MAX( nIterations, MaxIterations_K(iX1,iX2,iX3) )
            AveIterations_K(iX1,iX2,iX3) &
              = AveIterations_K(iX1,iX2,iX3) + nIterations

          END IF

        END DO

        CALL TimersStop( Timer_Im_Newton )

#else

        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'ComputeIncrement_TwoMoment_Implicit_New'
        WRITE(*,'(A6,A)') &
          '', 'in TwoMoment_DiscretizationModule_Collisions_Neutrinos'
        WRITE(*,'(A6,A)') &
          '', 'Invalid NEUTRINO_MATTER_SOLVER'
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A8,A)') '', 'EMAB'
        WRITE(*,'(A8,A)') '', 'FIXED_POINT_COUPLED'
        WRITE(*,'(A8,A)') '', 'FIXED_POINT_NESTED_AA'
        WRITE(*,'(A8,A)') '', 'FIXED_POINT_NESTED_NEWTON'
        WRITE(*,'(A8,A)') '', 'NEWTON'
        WRITE(*,*)
        STOP

#endif

        IF( TallyNonlinearSolver )THEN

          AveIterations_K(iX1,iX2,iX3) &
            = FLOOR( DBLE(AveIterations_K(iX1,iX2,iX3)     )/DBLE(nDOFX) )
          AveIterationsInner_K(iX1,iX2,iX3) &
            = FLOOR( DBLE(AveIterationsInner_K(iX1,iX2,iX3))/DBLE(nDOFX) )

        END IF

      END IF

      CALL TimersStop( Timer_Im_Solve )

      CALL TimersStart( Timer_Im_Increment )

      ! --- Update Radiation Fields ---

      DO iNodeX = 1, nDOFX
      DO iS     = 1, nSpecies

        Chi_T = Chi(:,iS,iNodeX) &
                  + Chi_NES(:,iS,iNodeX) + Chi_Pair(:,iS,iNodeX)

        Kappa = Chi_T + Sig(:,iS,iNodeX)

        Eta_T = Chi(:,iS,iNodeX) * fEQ(:,iS,iNodeX) &
                  + Eta_NES(:,iS,iNodeX) + Eta_Pair(:,iS,iNodeX)

        ! --- Number Density Updated in Nonlinear Solvers ---

        ! --- Number Flux (1) ---

        CR_N(:,iCR_G1,iS,iNodeX) &
          = CR_N(:,iCR_G1,iS,iNodeX) / ( One + dt * Kappa )

        ! --- Number Flux (2) ---

        CR_N(:,iCR_G2,iS,iNodeX) &
          = CR_N(:,iCR_G2,iS,iNodeX) / ( One + dt * Kappa )

        ! --- Number Flux (3) ---

        CR_N(:,iCR_G3,iS,iNodeX) &
          = CR_N(:,iCR_G3,iS,iNodeX) / ( One + dt * Kappa )

        ! --- Increments ---

        ! dR_N(:,iCR_N,iS,iNodeX) &
        !   = Eta_T - Chi_T * CR_N(:,iCR_N,iS,iNodeX)

        !!! NEW: make the new number density consistent to FP output
        ! increment = (new_moment - old_moment) / dt
        ! thus new_moment = old_moment + dt * increment
        dR_N(:,iCR_N,iS,iNodeX) &
          = ( CR_N(:,iCR_N,iS,iNodeX) - dR_N(:,iCR_N,iS,iNodeX) ) / dt
        !!! END: make the new number density consistent to FP output

        dR_N(:,iCR_G1,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G1,iS,iNodeX)

        dR_N(:,iCR_G2,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G2,iS,iNodeX)

        dR_N(:,iCR_G3,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G3,iS,iNodeX)

      END DO
      END DO





      CALL TimersStop( Timer_Im_Increment )

      CALL TimersStart( Timer_Im_EosOut )

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), &
               PF_N(:,iPF_E), AF_N(:,iAF_E), PF_N(:,iPF_Ne) )

      CALL TimersStop( Timer_Im_EosOut )

      CALL ComputeConserved_Euler &
             ( PF_N(:,iPF_D ), PF_N(:,iPF_V1), PF_N(:,iPF_V2), &
               PF_N(:,iPF_V3), PF_N(:,iPF_E ), PF_N(:,iPF_Ne), &
               CF_N(:,iCF_D ), CF_N(:,iCF_S1), CF_N(:,iCF_S2), &
               CF_N(:,iCF_S3), CF_N(:,iCF_E ), CF_N(:,iCF_Ne), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL TimersStart( Timer_Im_MapBackward )

      ! --- Conserved Fluid Increment ---

      DO iCF = 1, nCF

        dU_F(:,iX1,iX2,iX3,iCF) &
          = ( CF_N(:,iCF) - dU_F(:,iX1,iX2,iX3,iCF) ) / dt

      END DO

      CALL MapBackward_R_New &
             ( iE_B0, iE_E0, dU_R(:,iE_B0:iE_E0,iX1,iX2,iX3,:,:), dR_N )

      CALL TimersStop( Timer_Im_MapBackward )

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Im_Out )

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', Iterations_Min
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', Iterations_Max
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', &
        DBLE( Iterations_Ave ) / DBLE( PRODUCT(iX_E0-iX_B0+1)*nDOFX )
      WRITE(*,*)

    END IF

    DEALLOCATE &
      ( Kappa, Chi_T, Eta_T, Chi, fEQ, Sig, Chi_NES, Eta_NES, &
        Chi_Pair, Eta_Pair )


    CALL FinalizeCollisions_New

    CALL TimersStop( Timer_Im_Out )

    CALL TimersStop( Timer_Implicit )

#ifdef THORNADO_DEBUG_IMPLICIT
    WRITE(*,'(a,8x,5i4,es23.15)') 'MINLOC(dU_F), MINVAL(dU_F)', MINLOC(dU_F), MINVAL(dU_F)
    WRITE(*,'(a,8x,5i4,es23.15)') 'MAXLOC(dU_F), MAXVAL(dU_F)', MAXLOC(dU_F), MAXVAL(dU_F)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(dU_R), MINVAL(dU_F)', MINLOC(dU_R), MINVAL(dU_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(dU_R), MAXVAL(dU_F)', MAXLOC(dU_R), MAXVAL(dU_R)
#endif

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_New


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER               :: iCR, iS, iCF, iGF, iX_G
    REAL(DP)              :: PF_N(1:nPF)
    REAL(DP)              :: AF_N(1:nAF)
    REAL(DP), ALLOCATABLE :: Chi(:,:) ! --- Absorption Opacity
    REAL(DP), ALLOCATABLE :: Sig(:,:) ! --- Scattering Opacity
    REAL(DP), ALLOCATABLE :: fEQ(:,:) ! --- Equilibrium Distribution

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    wTime = MPI_WTIME( )

    CALL InitializeCollisions( iZ_B0, iZ_E0 )

    wTime = MPI_WTIME( ) - wTime

    !    WRITE(*,'(A4,A32,ES10.4E2)') '', 'InitializeCollisions: ', wTime

    ! --- Energy and Integration Weights ---

    wTime = MPI_WTIME( )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    wTime = MPI_WTIME( ) - wTime

    !    WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputePointsAndWeightsE: ', wTime

    ! --- Map Spatial Geometry Data ---

    DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33

      CALL MapForward_GX &
             ( iX_B0, iX_E0, &
               GX(:,iX_B0(1):iX_E0(1), &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3),iGF), &
               GX_N(iGF,1:nX_G) )

    END DO

    ! --- Map Fluid Data for Collision Update ---

    DO iCF = 1, nCF

      CALL MapForward_F &
             ( iX_B0, iX_E0, &
               U_F(:,iX_B0(1):iX_E0(1), &
                     iX_B0(2):iX_E0(2), &
                     iX_B0(3):iX_E0(3),iCF), &
               CF_N(iCF,1:nX_G) )

    END DO

    ! --- Map Radiation Data for Collision Update ---

    wTime = MPI_WTIME( )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U_R(:,iZ_B0(1):iZ_E0(1), &
                       iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3), &
                       iZ_B0(4):iZ_E0(4),iCR,iS), &
                 CR_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    wTime = MPI_WTIME( ) - wTime

    !    WRITE(*,'(A4,A32,ES10.4E2)') '', 'MapForward_R: ', wTime

    ! --- Allocate Local Opacities ---

    ALLOCATE( Chi(nE_G,nSpecies), Sig(nE_G,nSpecies), fEQ(nE_G,nSpecies) )

    ! --- Implicit Update ---

    Iterations_Min = + HUGE( 1 )
    Iterations_Max = - HUGE( 1 )
    Iterations_Ave = 0

    DO iX_G = 1, nX_G

      dF_N(1:nCF,iX_G) = CF_N(1:nCF,iX_G)

      wTime = MPI_WTIME( )

      CALL ComputePrimitive_Euler &
             ( CF_N(iCF_D ,iX_G), CF_N(iCF_S1,iX_G), CF_N(iCF_S2,iX_G), &
               CF_N(iCF_S3,iX_G), CF_N(iCF_E ,iX_G), CF_N(iCF_Ne,iX_G), &
               PF_N(iPF_D ),      PF_N(iPF_V1),      PF_N(iPF_V2),      &
               PF_N(iPF_V3),      PF_N(iPF_E ),      PF_N(iPF_Ne),      &
               GX_N(iGF_Gm_dd_11,iX_G), &
               GX_N(iGF_Gm_dd_22,iX_G), &
               GX_N(iGF_Gm_dd_33,iX_G) )

      wTime = MPI_WTIME( ) - wTime

      !      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputePrimitive_Euler: ', wTime

      wTime = MPI_WTIME( )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( PF_N(iPF_D:iPF_D), PF_N(iPF_E:iPF_E), PF_N(iPF_Ne:iPF_Ne), &
               AF_N(iAF_T:iAF_T), AF_N(iAF_E:iAF_E), AF_N(iAF_Ye:iAF_Ye) )

      wTime = MPI_WTIME( ) - wTime

      !      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputeTS_Auxiliary: ', wTime

      wTime = MPI_WTIME( )

      ! --- Electron Capture Opacities ---

      CALL ComputeNeutrinoOpacities_EC_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), &
               iSpecies = 1, opEC_Point = Chi(1:nE_G,1) )

      ! --- Elastic Scattering Opacities ---

      CALL ComputeNeutrinoOpacities_ES_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), &
               iSpecies = 1, iMoment = 1, opES_Point = Sig(1:nE_G,1) )

      Sig(1:nE_G,1) = Sig(1:nE_G,1) * FourPi * E_N(1:nE_G) * E_N(1:nE_G)

      wTime = MPI_WTIME( ) - wTime

      !      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputeNuOp_Point: ', wTime

      wTime = MPI_WTIME( )

      CALL SolveMatterEquations_EmAb_NuE &
             ( CR_N(:,iCR_N,iNuE,iX_G), dt * Chi(:,iNuE), fEQ(:,iNuE), &
               PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), AF_N(iAF_E) )

      wTime = MPI_WTIME( ) - wTime

      !      WRITE(*,'(A4,A32,ES10.4E2)') '', 'SolveMatterEquations_EmAb: ', wTime

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( PF_N(iPF_D:iPF_D), AF_N(iAF_T:iAF_T), AF_N(iAF_Ye:iAF_Ye), &
               PF_N(iPF_E:iPF_E), AF_N(iAF_E:iAF_E), PF_N(iPF_Ne:iPF_Ne) )

      CALL ComputeConserved_Euler &
             ( PF_N(iPF_D ),      PF_N(iPF_V1),      PF_N(iPF_V2),      &
               PF_N(iPF_V3),      PF_N(iPF_E ),      PF_N(iPF_Ne),      &
               CF_N(iCF_D ,iX_G), CF_N(iCF_S1,iX_G), CF_N(iCF_S2,iX_G), &
               CF_N(iCF_S3,iX_G), CF_N(iCF_E ,iX_G), CF_N(iCF_Ne,iX_G), &
               GX_N(iGF_Gm_dd_11,iX_G), &
               GX_N(iGF_Gm_dd_22,iX_G), &
               GX_N(iGF_Gm_dd_33,iX_G) )

      ! --- Conserved Fluid Increment ---

      dF_N(1:nCF,iX_G) = ( CF_N(1:nCF,iX_G) - dF_N(1:nCF,iX_G) ) / dt

      ! --- Update Radiation Fields ---

      DO iS = 1, nSpecies

        ! --- Number Density ---

        CR_N(:,iCR_N,iS,iX_G) &
          = ( dt * Chi(:,iS) * fEQ(:,iS) &
              + CR_N(:,iCR_N,iS,iX_G) ) / ( One + dt * Chi(:,iS) )

        ! --- Number Flux (1) ---

        CR_N(:,iCR_G1,iS,iX_G) &
          = CR_N(:,iCR_G1,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (2) ---

        CR_N(:,iCR_G2,iS,iX_G) &
          = CR_N(:,iCR_G2,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (3) ---

        CR_N(:,iCR_G3,iS,iX_G) &
          = CR_N(:,iCR_G3,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Increments ---

        dR_N(:,iCR_N,iS,iX_G) &
          = Chi(:,iS) * ( fEQ(:,iS) - CR_N(:,iCR_N,iS,iX_G) )

        dR_N(:,iCR_G1,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G1,iS,iX_G)

        dR_N(:,iCR_G2,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G2,iS,iX_G)

        dR_N(:,iCR_G3,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G3,iS,iX_G)

      END DO ! iS

    END DO ! iX_G

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', Iterations_Min
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', Iterations_Max
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', &
        DBLE( Iterations_Ave ) / DBLE( nX_G )
      WRITE(*,*)

    END IF

    ! --- Deallocate Local Opacities ---

    DEALLOCATE( Chi, Sig, fEQ )

    ! --- Map Increments Back ---

    ! --- Fluid ---

    DO iCF = 1, nCF

      CALL MapBackward_F &
             ( iX_B0, iX_E0, &
               dU_F(:,iX_B0(1):iX_E0(1), &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(3):iX_E0(3),iCF), &
               dF_N(iCF,1:nX_G) )

    END DO

    ! --- Radiation ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU_R(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                        iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dR_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    CALL FinalizeCollisions

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit_DGFV &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER               :: iX1, iX2, iX3, iGF, iCF, iCR, iS, iNX, iE
    REAL(DP)              :: wTimeTotal
    REAL(DP)              :: GX_K(1:nGF)
    REAL(DP)              :: CF_K(1:nCF)
    REAL(DP)              :: PF_K(1:nPF)
    REAL(DP)              :: AF_K(1:nAF)
    REAL(DP), ALLOCATABLE :: Chi(:,:) ! --- Absorption Opacity
    REAL(DP), ALLOCATABLE :: Sig(:,:) ! --- Scattering Opacity
    REAL(DP), ALLOCATABLE :: fEQ(:,:) ! --- Equilibrium Distribution

    wTimeTotal = MPI_WTIME( )

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    CALL InitializeCollisions_DGFV( iZ_B0, iZ_E0 )

    ! --- Energy and Integration Weights ---

    CALL ComputePointsAndWeightsE_DGFV( E_N, W2_N, W3_N )

    ! --- Allocate Local Opacities ---

    ALLOCATE( Chi(nE_G,nSpecies), Sig(nE_G,nSpecies), fEQ(nE_G,nSpecies) )

    ! --- Implicit Update ---

    Iterations_Min = + HUGE( 1 )
    Iterations_Max = - HUGE( 1 )
    Iterations_Ave = 0

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33

        GX_K(iGF) = DOT_PRODUCT( WeightsX_q(:), GX(:,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CF_K(iCF) = DOT_PRODUCT( WeightsX_q(:), U_F(:,iX1,iX2,iX3,iCF) )

        dU_F(:,iX1,iX2,iX3,iCF) = CF_K(iCF)

      END DO

      CALL ComputePrimitive_Euler &
             ( CF_K(iCF_D ), CF_K(iCF_S1), CF_K(iCF_S2), &
               CF_K(iCF_S3), CF_K(iCF_E ), CF_K(iCF_Ne), &
               PF_K(iPF_D ), PF_K(iPF_V1), PF_K(iPF_V2), &
               PF_K(iPF_V3), PF_K(iPF_E ), PF_K(iPF_Ne), &
               GX_K(iGF_Gm_dd_11), GX_K(iGF_Gm_dd_22), GX_K(iGF_Gm_dd_33) )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( PF_K(iPF_D:iPF_D), PF_K(iPF_E:iPF_E), PF_K(iPF_Ne:iPF_Ne), &
               AF_K(iAF_T:iAF_T), AF_K(iAF_E:iAF_E), AF_K(iAF_Ye:iAF_Ye) )

      ! --- Electron Capture Opacities ---

      CALL ComputeNeutrinoOpacities_EC_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_K(iPF_D), AF_K(iAF_T), AF_K(iAF_Ye), &
               iSpecies = 1, opEC_Point = Chi(1:nE_G,1) )

      ! --- Elastic Scattering Opacities ---

      CALL ComputeNeutrinoOpacities_ES_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_K(iPF_D), AF_K(iAF_T), AF_K(iAF_Ye), &
               iSpecies = 1, iMoment = 1, opES_Point = Sig(1:nE_G,1) )

      Sig(1:nE_G,1) = Sig(1:nE_G,1) * FourPi * E_N(1:nE_G) * E_N(1:nE_G)

      ! --- Map Radiation Data for Collision Update ---

      DO iS  = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R_DGFV &
               ( iE_B0, iE_E0, U_R(:,iE_B0:iE_E0,iX1,iX2,iX3,iCR,iS), &
                 CR_N(1:nE_G,iCR,iS,1:nDOFX), CR_K(1:nE_G,iCR,iS) )

      END DO
      END DO

      ! --- Update Fluid ---

      CALL SolveMatterEquations_EmAb_NuE &
             ( CR_K(:,iCR_N,iNuE), dt * Chi(:,iNuE), fEQ(:,iNuE), &
               PF_K(iPF_D), AF_K(iAF_T), AF_K(iAF_Ye), AF_K(iAF_E) )

      ! --- Compute Primitive Fluid ---

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( PF_K(iPF_D:iPF_D), AF_K(iAF_T:iAF_T), AF_K(iAF_Ye:iAF_Ye), &
               PF_K(iPF_E:iPF_E), AF_K(iAF_E:iAF_E), PF_K(iPF_Ne:iPF_Ne) )

      ! --- Compute Conserved Fluid ---

      CALL ComputeConserved_Euler &
             ( PF_K(iPF_D ), PF_K(iPF_V1), PF_K(iPF_V2), &
               PF_K(iPF_V3), PF_K(iPF_E ), PF_K(iPF_Ne), &
               CF_K(iCF_D ), CF_K(iCF_S1), CF_K(iCF_S2), &
               CF_K(iCF_S3), CF_K(iCF_E ), CF_K(iCF_Ne), &
               GX_K(iGF_Gm_dd_11), GX_K(iGF_Gm_dd_22), GX_K(iGF_Gm_dd_33) )

      ! --- Conserved Fluid Increment ---

      DO iCF = 1, nCF

        dU_F(:,iX1,iX2,iX3,iCF) = ( CF_K(iCF) - dU_F(:,iX1,iX2,iX3,iCF) ) / dt

      END DO

      ! --- Update Radiation Fields ---

      DO iNX = 1, nDOFX
      DO iS  = 1, nSpecies

        ! --- Number Density ---

        CR_N(:,iCR_N,iS,iNX) &
          = ( dt * Chi(:,iS) * fEQ(:,iS) &
              + CR_N(:,iCR_N,iS,iNX) ) / ( One + dt * Chi(:,iS) )

        ! --- Number Flux (1) ---

        CR_N(:,iCR_G1,iS,iNX) &
          = CR_N(:,iCR_G1,iS,iNX) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (2) ---

        CR_N(:,iCR_G2,iS,iNX) &
          = CR_N(:,iCR_G2,iS,iNX) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (3) ---

        CR_N(:,iCR_G3,iS,iNX) &
          = CR_N(:,iCR_G3,iS,iNX) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Increments ---

        dR_N(:,iCR_N,iS,iNX) &
          = Chi(:,iS) * ( fEQ(:,iS) - CR_N(:,iCR_N,iS,iNX) )

        dR_N(:,iCR_G1,iS,iNX) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G1,iS,iNX)

        dR_N(:,iCR_G2,iS,iNX) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G2,iS,iNX)

        dR_N(:,iCR_G3,iS,iNX) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G3,iS,iNX)

      END DO
      END DO

      ! --- Map Radiation Increments Back ---

      DO iS  = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R_DGFV &
               ( iE_B0, iE_E0, dU_R(:,iE_B0:iE_E0,iX1,iX2,iX3,iCR,iS), &
                 dR_N(1:nE_G,iCR,iS,1:nDOFX) )

      END DO
      END DO

    END DO
    END DO
    END DO

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', Iterations_Min
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', Iterations_Max
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', &
        DBLE( Iterations_Ave ) / DBLE( nX_G )
      WRITE(*,*)

    END IF

    ! --- Deallocate Local Opacities ---

    DEALLOCATE( Chi, Sig, fEQ )

    CALL FinalizeCollisions_DGFV

    wTimeTotal = MPI_WTIME( ) - wTimeTotal

    !    WRITE(*,'(A4,A32,ES10.4E2)') '', 'wTimeTotal: ', wTimeTotal

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_DGFV


  SUBROUTINE SolveMatterEquations_EmAb_NuE( J, Chi, J0, D, T, Y, E )

    REAL(DP), INTENT(inout) :: J  (1:nE_G)
    REAL(DP), INTENT(in)    :: Chi(1:nE_G)
    REAL(DP), INTENT(inout) :: J0 (1:nE_G)
    REAL(DP), INTENT(inout) :: D, T, Y, E

    INTEGER,  PARAMETER :: iY = 1, iE = 2
    INTEGER,  PARAMETER :: MaxIter = 20
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    LOGICAL  :: CONVERGED
    INTEGER  :: k
    REAL(DP) :: Yold, Eold, N_B
    REAL(DP) :: W2_S(1:nE_G), W3_S(1:nE_G)
    REAL(DP) :: Theta2_N(1:nE_G), Theta3_N(1:nE_G)
    REAL(DP) :: Mnu, dMnudT, dMnudY
    REAL(DP) :: dEdT, dEdY
    REAL(DP) :: TMP(1), dTMPdT(1), dTMPdY(1)
    REAL(DP) :: dJ0dT_Y(1:nE_G), dJ0dY_T(1:nE_G)
    REAL(DP) :: dJ0dE_Y(1:nE_G), dJ0dY_E(1:nE_G)
    REAL(DP) :: U(2), dU(2), C(2), FVEC(2), FVEC0(2)
    REAL(DP) :: FJAC(2,2), IJAC(2,2), DJAC

    Yold = Y; Eold = E

    ! --- Auxiliary Variables ---

    N_B = D / AtomicMassUnit

    IF( SolveMatter )THEN

      W2_S = W2_N / PlanckConstant**3
      W3_S = W3_N / PlanckConstant**3

      Theta2_N = FourPi * W2_S * Chi / ( One + Chi )
      Theta3_N = FourPi * W3_S * Chi / ( One + Chi )

    ELSE

      Theta2_N = Zero
      Theta3_N = Zero

    END IF

    ! --- Neutrino Chemical Potential and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu, dMnudT, dMnudY, iSpecies = iNuE )

    ! --- Equilibrium Distribution ---

    J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T )

    ! --- Initial Guess ---

    U(iY) = Yold; U(iE) = Eold

    ! --- Old States (Constant) ---

    C(iY) = DOT_PRODUCT( Theta2_N(:), J(:) ) + N_B * U(iY)
    C(iE) = DOT_PRODUCT( Theta3_N(:), J(:) ) + D   * U(iE)

    ! --- Electron Fraction Equation ---

    FVEC(iY) = DOT_PRODUCT( Theta2_N(:), J0(:) ) + N_B * U(iY) - C(iY)

    ! --- Internal Energy Equation ---

    FVEC(iE) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + D   * U(iE) - C(iE)

    ! --- Scale Equations and Save Initial Evaluation ---

    FVEC(:) = FVEC(:) / C(:); FVEC0(:) = FVEC(:);

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED )

      k = k + 1

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( [ D ], [ T ], [ Y ], E = TMP, &
               dEdT_Option = dTMPdT, dEdY_Option = dTMPdY )

      dEdT = dTMPdT(1); dEdY = dTMPdY(1)

      ! --- Derivative of J0 wrt. T (Constant Y) ---

      dJ0dT_Y = dFermiDiracdT( E_N, Mnu, BoltzmannConstant * T, dMnudT, T )

      ! --- Derivative of J0 wrt. T (Constant T) ---

      dJ0dY_T = dFermiDiracdY( E_N, Mnu, BoltzmannConstant * T, dMnudY, T )

      ! --- Derivative of J0 wrt. E (Constant Y) ---

      dJ0dE_Y = dJ0dT_Y / dEdT

      ! --- Derivative of J0 wrt. Y (Constant E) ---

      dJ0dY_E = dJ0dY_T - dJ0dT_Y * dEdY / dEdT

      ! --- Jacobian ---

      FJAC(1,1) = DOT_PRODUCT( Theta2_N(:), dJ0dY_E(:) ) + N_B

      FJAC(1,2) = DOT_PRODUCT( Theta2_N(:), dJ0dE_Y(:) )

      FJAC(2,1) = DOT_PRODUCT( Theta3_N(:), dJ0dY_E(:) )

      FJAC(2,2) = DOT_PRODUCT( Theta3_N(:), dJ0dE_Y(:) ) + D

      ! --- Scale Jacobian ---

      FJAC(:,1) = FJAC(:,1) / C(:)

      FJAC(:,2) = FJAC(:,2) / C(:)

      ! --- Determinant of Jacobian ---

      DJAC = FJAC(1,1) * FJAC(2,2) - FJAC(2,1) * FJAC(1,2)

      ! --- Invert Jacobian ---

      IJAC(1,1) =   FJAC(2,2) / DJAC
      IJAC(2,1) = - FJAC(2,1) / DJAC
      IJAC(1,2) = - FJAC(1,2) / DJAC
      IJAC(2,2) =   FJAC(1,1) / DJAC

      ! --- Correction ---

      dU = - MATMUL( IJAC, FVEC )

      ! --- Apply Correction ---

      U = U + dU

      Y = U(1); E = U(2)

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

      ! --- Neutrino Chemical Potential and Derivatives ---

      CALL ComputeNeutrinoChemicalPotentials &
             ( D, T, Y, Mnu, dMnudT, dMnudY, iSpecies = 1 )

      ! --- Equilibrium Distribution ---

      J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T )

      ! --- Electron Fraction Equation ---

      FVEC(1) = DOT_PRODUCT( Theta2_N(:), J0(:) ) + N_B * U(1) - C(1)

      ! --- Internal Energy Equation ---

      FVEC(2) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + D   * U(2) - C(2)

      ! --- Scale Equations ---

      FVEC(:) = FVEC(:) / C(:)

      ! --- Check for Convergence ---

      IF( ENORM( FVEC ) < Rtol * ENORM( FVEC0 ) .OR. &
          ENORM( dU/U ) < Utol )THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

      END IF

      IF( ( k .EQ. MaxIter ) .AND. ( .NOT. CONVERGED ) )THEN

        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'SolveMatterEquations_EmAb:'
        WRITE(*,'(A6,A20,I4.4,A11)') &
          '', 'Did not converge in ', k, ' iterations'
        WRITE(*,'(A6,A)') &
          '', 'Exiting with unconverged result'
        WRITE(*,*)
        WRITE(*,'(A4,A12,ES10.4E2,A2,A10,ES10.4E2,A2,A4,ES10.4E2)') &
          '', 'D [g/ccm] = ', D / Unit_D, &
          '', 'T [MeV] = ', T / Unit_T, &
          '', 'Y = ', Y
        WRITE(*,*)
        WRITE(*,'(A4,A24,3ES12.4E2)') &
          '', '|F(Y)|, |F0(Y)|, Rtol = ', ABS( FVEC(1) ), ABS( FVEC0(1) ), Rtol
        WRITE(*,'(A4,A24,2ES12.4E2)') &
          '', '|dY/Y|, Ytol = ', ABS( dU(1) / U(1) ), Utol
        WRITE(*,'(A4,A24,3ES12.4E2)') &
          '', '|F(E)|, |F0(E)|, Rtol = ', ABS( FVEC(2) ), ABS( FVEC0(2) ), Rtol
        WRITE(*,'(A4,A24,2ES12.4E2)') &
          '', '|dE/E|, Etol = ', ABS( dU(2) / U(2) ), Utol
        WRITE(*,*)

        CONVERGED = .TRUE.

      END IF

    END DO

    J = ( J + Chi * J0 ) / ( One + Chi )

  END SUBROUTINE SolveMatterEquations_EmAb_NuE


  SUBROUTINE SolveMatterEquations_EmAb &
    ( J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

    ! --- Electron Neutrinos (1) and Electron Antineutrinos (2) ---

    REAL(DP), INTENT(in)    :: J_1  (1:nE_G), J_2  (1:nE_G)
    REAL(DP), INTENT(in)    :: Chi_1(1:nE_G), Chi_2(1:nE_G)
    REAL(DP), INTENT(inout) :: J0_1 (1:nE_G), J0_2 (1:nE_G)
    REAL(DP), INTENT(inout) :: D, T, Y, E

    REAL(DP) :: Mnu_1, dMnudT_1, dMnudY_1
    REAL(DP) :: Mnu_2, dMnudT_2, dMnudY_2

    ! --- Neutrino Chemical Potentials and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_1, dMnudT_1, dMnudY_1, iSpecies = iNuE )

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_2, dMnudT_2, dMnudY_2, iSpecies = iNuE_Bar )

    ! --- Equilibrium Distributions ---

    J0_1 = FermiDirac( E_N, Mnu_1, BoltzmannConstant * T )

    J0_2 = FermiDirac( E_N, Mnu_2, BoltzmannConstant * T )

  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE SolveMatterEquations_EmAb_FP &
    ( dt, iS_1, iS_2, J, Chi, J0, D, T, Y, E, nIterations_Out, TOL )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    INTEGER,  INTENT(out), OPTIONAL :: nIterations_Out
    REAL(DP), INTENT(in),  OPTIONAL :: TOL

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * M




    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    INTEGER  :: OS_1, OS_2
    REAL(DP) :: h3, N_B
    REAL(DP) :: S_Y, S_E, Yold, Eold
    REAL(DP) :: TMP(1)
    REAL(DP) :: C(2), Unew(2), Jnorm(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: Alpha(1:M)
    REAL(DP) :: WORK(1:LWORK)
    REAL(DP) :: GVECm(1:2*(1+nE_G))
    REAL(DP) :: FVECm(1:2*(1+nE_G))
    REAL(DP) :: GVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: FVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta(1:nE_G,1:2)
    REAL(DP) :: BVEC(1:2*(1+nE_G))
    REAL(DP) :: AMAT(1:2*(1+nE_G),1:M)

    INTEGER  :: nIterations
    REAL(DP) :: Rtol

    IF(PRESENT(TOL)) THEN
       Rtol = TOL
     ELSE
       Rtol = 1.0d-08
    END IF

    OS_1 = 2
    OS_2 = 2 + nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J


    S_Y = N_B * Yold
    S_E = D   * Eold

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E

    Jnorm(1) = WNORM(Jold(:,1), W2_N)
    Jnorm(2) = WNORM(Jold(:,2), W2_N)

    Unew(iY) = Y / Yold ! --- Initial Guess
    Unew(iE) = E / Eold ! --- Initial Guess
    Jnew     = J        ! --- Initial Guess

    ! CALL TimersStart(Timer_Im_ComputeOpacity)
    !
    ! CALL ComputeEquilibriumDistributions_Point &
    !        ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )
    !
    ! CALL ComputeEquilibriumDistributions_Point &
    !        ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )
    !
    ! CALL TimersStop(Timer_Im_ComputeOpacity)
    !
    ! Eta(:,1) = Chi(:,1) * J0(:,1)
    ! Eta(:,2) = Chi(:,2) * J0(:,2)
    !
    !
    ! ! --- Update Neutrino Densities ---
    !
    ! Jnew(:,1) = ( Jold(:,1) + dt * Eta(:,1) ) / ( One + dt * Chi(:,1) )
    !
    ! Jnew(:,2) = ( Jold(:,2) + dt * Eta(:,2) ) / ( One + dt * Chi(:,2) )

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )

      CALL TimersStart(Timer_Im_ComputeOpacity)

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

      CALL TimersStop(Timer_Im_ComputeOpacity)

      CALL TimersStart( Timer_Im_ComputeRate )

      Eta(:,1) = Chi(:,1) * J0(:,1)
      Eta(:,2) = Chi(:,2) * J0(:,2)

      CALL TimersStop( Timer_Im_ComputeRate )

      ! --- Right-Hand Side Vectors ---

      GVEC(iY,mk) = One + C(iY) &
                      - DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y
      GVEC(iE,mk) = One + C(iE) &
                      - DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      GVEC(OS_1+1:OS_1+nE_G,mk) &
        = ( Jold(:,1) + dt * Eta(:,1) ) / ( One + dt * Chi(:,1) )

      GVEC(OS_2+1:OS_2+nE_G,mk) &
        = ( Jold(:,2) + dt * Eta(:,2) ) / ( One + dt * Chi(:,2) )

      ! --- Residuals ---

      FVEC(iY,mk) = GVEC(iY,mk) - Unew(iY)

      FVEC(iE,mk) = GVEC(iE,mk) - Unew(iE)


      DO i = 1, nE_G
        FVEC(OS_1+i,mk) = GVEC(OS_1+i,mk) - Jnew(i,1)
      END DO

      DO i = 1, nE_G
        FVEC(OS_2+i,mk) = GVEC(OS_2+i,mk) - Jnew(i,2)
      END DO

      IF( ENORM( [ FVEC(iY,mk) ] ) <= Rtol .AND. &
          ENORM( [ FVEC(iE,mk) ] ) <= Rtol .AND. &
          WNORM( FVEC(OS_1+1:OS_1+nE_G,mk) * ( One + dt * Chi(:,1) ), W2_N ) &
          <= Rtol * Jnorm(1) .AND. &
          WNORM( FVEC(OS_2+1:OS_2+nE_G,mk) * ( One + dt * Chi(:,2) ), W2_N ) &
          <= Rtol * Jnorm(2) ) &
      THEN

        CONVERGED = .TRUE.
        nIterations = k

      END IF

      IF( .NOT. CONVERGED )THEN
        CALL TimersStart( Timer_Im_ComputeLS )

        IF( mk == 1 )THEN

          GVECm = GVEC(:,mk)

        ELSE

          BVEC(:) &
            = - FVEC(:,mk)
          AMAT(:,1:mk-1) &
            = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

          CALL DGELS( 'N', 2*(nE_G+1), mk-1, 1, AMAT(:,1:mk-1), 2*(nE_G+1), &
                      BVEC, 2*(nE_G+1), WORK, LWORK, INFO )

          Alpha(1:mk-1) = BVEC(1:mk-1)

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        END IF

        CALL TimersStop( Timer_Im_ComputeLS )

        CALL TimersStart( Timer_Im_UpdateFP )

        Unew(iY)  = GVECm(iY)
        Unew(iE)  = GVECm(iE)
        Jnew(:,1) = GVECm(OS_1+1:OS_1+nE_G)
        Jnew(:,2) = GVECm(OS_2+1:OS_2+nE_G)

        IF( mk == M )THEN

          GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
          FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

        END IF

        ! --- Update Matter ---

        Y = Unew(iY) * Yold
        E = Unew(iE) * Eold

        CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
               ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

        CALL TimersStop( Timer_Im_UpdateFP )
      END IF
    END DO

    ! output J for preconditioning purpose
     J = Jnew

    IF(PRESENT(nIterations_Out)) THEN
      nIterations_Out = nIterations
    END IF

  END SUBROUTINE SolveMatterEquations_EmAb_FP


  SUBROUTINE SolveMatterEquations_Presolve &
    ( dt, iS_1, iS_2, J, Chi, J0, Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, &
      Phi_0_Ot_Pair, D, T, Y, E, LagAllButJ0, TOL)

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Phi_0_In_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Phi_0_Ot_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Phi_0_In_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    LOGICAL,  INTENT(in)    :: LagAllButJ0
    REAL(DP), INTENT(in),  OPTIONAL :: TOL


    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * M


    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    INTEGER  :: OS_1, OS_2
    REAL(DP) :: h3, N_B
    REAL(DP) :: S_Y, S_E, Yold, Eold
    REAL(DP) :: TMP(1)
    REAL(DP) :: C(2), Unew(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: Alpha(1:M)
    REAL(DP) :: WORK(1:LWORK)
    REAL(DP) :: GVECm(1:2*(1+nE_G))
    REAL(DP) :: FVECm(1:2*(1+nE_G))
    REAL(DP) :: GVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: FVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta(1:nE_G,1:2)
    REAL(DP) :: BVEC(1:2*(1+nE_G))
    REAL(DP) :: AMAT(1:2*(1+nE_G),1:M)
    REAL(DP) :: Chi_NES (1:nE_G,1:2)
    REAL(DP) :: Eta_NES (1:nE_G,1:2)
    REAL(DP) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP) :: Eta_Pair(1:nE_G,1:2)

    REAL(DP) :: Rtol

    IF(PRESENT(TOL)) THEN
       Rtol = TOL
     ELSE
       Rtol = 1.0d-08
    END IF

    OS_1 = 2
    OS_2 = 2 + nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J


    Unew(iY) = Y / Yold ! --- Initial Guess
    Unew(iE) = E / Eold ! --- Initial Guess
    Jnew     = J        ! --- Initial Guess

    S_Y = N_B * Yold
    S_E = D   * Eold

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E

    CALL TimersStart( Timer_Im_ComputeRate )

    ! --- Emissivities and Opacities ---

    Eta(:,1) = Chi(:,1) * J0(:,1)
    Eta(:,2) = Chi(:,2) * J0(:,2)

    ! --- NES Emissivities and Opacities ---

    ! --- Neutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                W2_N * Jnew(:,1),       1, Zero, Eta_NES(:,1), 1 )

    Chi_NES(:,1) = Eta_NES(:,1)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                W2_N * (One-Jnew(:,1)), 1, One,  Chi_NES(:,1), 1 )

    ! --- Antineutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                W2_N * Jnew(:,2),       1, Zero, Eta_NES(:,2), 1 )

    Chi_NES(:,2) = Eta_NES(:,2)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                W2_N * (One-Jnew(:,2)), 1, One,  Chi_NES(:,2), 1 )

    ! --- Pair Emissivities and Opacities ---

    ! --- Neutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

    Chi_Pair(:,1) = Eta_Pair(:,1)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                W2_N * Jnew(:,2),       1, One,  Chi_Pair(:,1), 1 )

    ! --- Antineutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

    Chi_Pair(:,2) = Eta_Pair(:,2)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                W2_N * Jnew(:,1),       1, One,  Chi_Pair(:,2), 1 )

    CALL TimersStop( Timer_Im_ComputeRate )


    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )


      ! --- Right-Hand Side Vectors ---

      GVEC(iY,mk) = One + C(iY) &
                      - DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y
      GVEC(iE,mk) = One + C(iE) &
                      - DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      GVEC(OS_1+1:OS_1+nE_G,mk) &
        = ( Jold(:,1) + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) &
              / ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) )

      GVEC(OS_2+1:OS_2+nE_G,mk) &
        = ( Jold(:,2) + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) &
              / ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) )

      ! --- Residuals ---

      FVEC(iY,mk) = GVEC(iY,mk) - Unew(iY)

      FVEC(iE,mk) = GVEC(iE,mk) - Unew(iE)

      DO i = 1, nE_G
        FVEC(OS_1+i,mk) = GVEC(OS_1+i,mk) - Jnew(i,1)
      END DO

      DO i = 1, nE_G
        FVEC(OS_2+i,mk) = GVEC(OS_2+i,mk) - Jnew(i,2)
      END DO

      CALL TimersStart( Timer_Im_ComputeLS )

      IF( mk == 1 )THEN

        GVECm = GVEC(:,mk)

      ELSE

        BVEC(:) &
          = - FVEC(:,mk)
        AMAT(:,1:mk-1) &
          = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 2*(nE_G+1), mk-1, 1, AMAT(:,1:mk-1), 2*(nE_G+1), &
                    BVEC, 2*(nE_G+1), WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVEC(1:mk-1)

        Alpha(mk) = One - SUM( Alpha(1:mk-1) )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      CALL TimersStop( Timer_Im_ComputeLS )

      CALL TimersStart( Timer_Im_UpdateFP )

      FVECm(iY)               = GVECm(iY)               - Unew(iY)
      FVECm(iE)               = GVECm(iE)               - Unew(iE)
      FVECm(OS_1+1:OS_1+nE_G) = GVECm(OS_1+1:OS_1+nE_G) - Jnew(:,1)
      FVECm(OS_2+1:OS_2+nE_G) = GVECm(OS_2+1:OS_2+nE_G) - Jnew(:,2)

      IF( ENORM( [ FVECm(iY) ] ) <= Rtol .AND. &
          ENORM( [ FVECm(iE) ] ) <= Rtol .AND. &
          ENORM( FVECm(OS_1+1:OS_1+nE_G) ) <= Rtol * ENORM( Jold(:,1) ) .AND. &
          ENORM( FVECm(OS_2+1:OS_2+nE_G) ) <= Rtol * ENORM( Jold(:,2) ) ) &
      THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k


      END IF

      Unew(iY)  = GVECm(iY)
      Unew(iE)  = GVECm(iE)
      Jnew(:,1) = GVECm(OS_1+1:OS_1+nE_G)
      Jnew(:,2) = GVECm(OS_2+1:OS_2+nE_G)

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF

      ! --- Update Matter ---

      Y = Unew(iY) * Yold
      E = Unew(iE) * Eold

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

      CALL TimersStart( Timer_Im_UpdateFP )

      IF( .NOT. CONVERGED )THEN

        ! --- Recompute Equilibrium Distributions and Emissivities ---
        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeEquilibriumDistributions_Point &
               ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

        CALL ComputeEquilibriumDistributions_Point &
               ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

        CALL TimersStop( Timer_Im_ComputeOpacity )
        ! --- Emissivities and Opacities ---

        Eta(:,1) = Chi(:,1) * J0(:,1)
        Eta(:,2) = Chi(:,2) * J0(:,2)

        IF( .NOT. LagAllButJ0 )THEN

          CALL TimersStart( Timer_Im_ComputeRate )

          ! --- NES Emissivities and Opacities ---

          ! --- Neutrino ---

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                      W2_N * Jnew(:,1),       1, Zero, Eta_NES(:,1), 1 )

          Chi_NES(:,1) = Eta_NES(:,1)

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                      W2_N * (One-Jnew(:,1)), 1, One,  Chi_NES(:,1), 1 )

          ! --- Antineutrino ---

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                      W2_N * Jnew(:,2),       1, Zero, Eta_NES(:,2), 1 )

          Chi_NES(:,2) = Eta_NES(:,2)

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                      W2_N * (One-Jnew(:,2)), 1, One,  Chi_NES(:,2), 1 )

          ! --- Pair Emissivities and Opacities ---

          ! --- Neutrino ---

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                      W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

          Chi_Pair(:,1) = Eta_Pair(:,1)

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                      W2_N * Jnew(:,2),       1, One,  Chi_Pair(:,1), 1 )

          ! --- Antineutrino ---

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                      W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

          Chi_Pair(:,2) = Eta_Pair(:,2)

          CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                      W2_N * Jnew(:,1),       1, One,  Chi_Pair(:,2), 1 )

          CALL TimersStop( Timer_Im_ComputeRate )

        END IF
      END IF

    END DO

    J = Jnew


  END SUBROUTINE SolveMatterEquations_Presolve


  SUBROUTINE SolveMatterEquations_FP_Coupled &
    ( dt, iS_1, iS_2, J, Chi, J0, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    INTEGER,  INTENT(out)   :: nIterations

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * M
    REAL(DP), PARAMETER :: Rtol = 1.0d-8
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    INTEGER  :: OS_1, OS_2
    REAL(DP) :: h3, N_B
    REAL(DP) :: S_Y, S_E, Yold, Eold, Told
    REAL(DP) :: TMP(1)
    REAL(DP) :: C(2), Unew(2), Jnorm(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: Alpha(1:M)
    REAL(DP) :: WORK(1:LWORK)
    REAL(DP) :: GVECm(1:2*(1+nE_G))
    ! REAL(DP) :: FVECm(1:2*(1+nE_G)) !!! NEW: removed vector for storing error
    REAL(DP) :: GVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: FVEC (1:2*(1+nE_G),1:M)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta(1:nE_G,1:2)
    REAL(DP) :: BVEC(1:2*(1+nE_G))
    REAL(DP) :: AMAT(1:2*(1+nE_G),1:M)
    REAL(DP) :: Phi_0_In_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_In_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:2)

    OS_1 = 2
    OS_2 = 2 + nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J


    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J, Chi, J0, D, T, Y, E)

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF



    S_Y = N_B * Yold
    S_E = D   * Eold

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E

    Jnorm(1) = WNORM(Jold(:,1), W2_N)
    Jnorm(2) = WNORM(Jold(:,2), W2_N)

    Unew(iY) = Y / Yold ! --- Initial Guess
    Unew(iE) = E / Eold ! --- Initial Guess
    Jnew     = J        ! --- Initial Guess


    



    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )

      CALL TimersStart( Timer_Im_ComputeOpacity )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

      ! --- NES Kernels ---

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
               Phi_0_In_NES(:,:,1), Phi_0_Ot_NES(:,:,1) )

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
               Phi_0_In_NES(:,:,2), Phi_0_Ot_NES(:,:,2) )

      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacities_Pair_Point &
            ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
              Phi_0_In_Pair(:,:,1), Phi_0_Ot_Pair(:,:,1) )

      CALL ComputeNeutrinoOpacities_Pair_Point &
            ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
              Phi_0_In_Pair(:,:,2), Phi_0_Ot_Pair(:,:,2) )

      CALL TimersStop( Timer_Im_ComputeOpacity )

      ! --- Emissivities and Opacities ---

      CALL TimersStart( Timer_Im_ComputeRate )

      Eta(:,1) = Chi(:,1) * J0(:,1)
      Eta(:,2) = Chi(:,2) * J0(:,2)

      ! --- NES Emissivities and Opacities ---

      ! --- Neutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                  W2_N * Jnew(:,1),       1, Zero, Eta_NES(:,1), 1 )

      Chi_NES(:,1) = Eta_NES(:,1)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, One,  Chi_NES(:,1), 1 )

      ! --- Antineutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                  W2_N * Jnew(:,2),       1, Zero, Eta_NES(:,2), 1 )

      Chi_NES(:,2) = Eta_NES(:,2)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, One,  Chi_NES(:,2), 1 )

      ! --- Pair Emissivities and Opacities ---

      ! --- Neutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

      Chi_Pair(:,1) = Eta_Pair(:,1)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                  W2_N * Jnew(:,2),       1, One,  Chi_Pair(:,1), 1 )

      ! --- Antineutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

      Chi_Pair(:,2) = Eta_Pair(:,2)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                  W2_N * Jnew(:,1),       1, One,  Chi_Pair(:,2), 1 )

      CALL TimersStop( Timer_Im_ComputeRate )

      ! --- Right-Hand Side Vectors ---

      GVEC(iY,mk) = One + C(iY) &
                      - DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y
      GVEC(iE,mk) = One + C(iE) &
                      - DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      GVEC(OS_1+1:OS_1+nE_G,mk) &
        = ( Jold(:,1) + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) &
              / ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) )

      GVEC(OS_2+1:OS_2+nE_G,mk) &
        = ( Jold(:,2) + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) &
              / ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) )

      ! --- Residuals ---

      FVEC(iY,mk) = GVEC(iY,mk) - Unew(iY)

      FVEC(iE,mk) = GVEC(iE,mk) - Unew(iE)

      DO i = 1, nE_G
        FVEC(OS_1+i,mk) = GVEC(OS_1+i,mk) - Jnew(i,1)
      END DO

      DO i = 1, nE_G
        FVEC(OS_2+i,mk) = GVEC(OS_2+i,mk) - Jnew(i,2)
      END DO

      IF( ENORM( [ FVEC(iY,mk) ] ) <= Rtol .AND. &
          ENORM( [ FVEC(iE,mk) ] ) <= Rtol .AND. &
          !WNORM( FVEC(OS_1+1:OS_1+nE_G,mk), W2_N ) <= Rtol * Jnorm(1) .AND. &
          !WNORM( FVEC(OS_2+1:OS_2+nE_G,mk), W2_N ) <= Rtol * Jnorm(2) ) &
          WNORM( FVEC(OS_1+1:OS_1+nE_G,mk) * &
          ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) ), W2_N ) &
          <= Rtol * Jnorm(1) .AND. &
          WNORM( FVEC(OS_2+1:OS_2+nE_G,mk) * &
          ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) ), W2_N ) &
          <= Rtol * Jnorm(2) ) &
      THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

        nIterations = k

      END IF

      IF( .NOT. CONVERGED )THEN

        CALL TimersStart( Timer_Im_ComputeLS )

        IF( mk == 1 )THEN

          GVECm = GVEC(:,mk)

        ELSE

          BVEC(:) &
            = - FVEC(:,mk)
          AMAT(:,1:mk-1) &
            = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

          CALL DGELS( 'N', 2*(nE_G+1), mk-1, 1, AMAT(:,1:mk-1), 2*(nE_G+1), &
                      BVEC, 2*(nE_G+1), WORK, LWORK, INFO )

          Alpha(1:mk-1) = BVEC(1:mk-1)

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        END IF

        CALL TimersStop( Timer_Im_ComputeLS )




        CALL TimersStart( Timer_Im_UpdateFP )

        Unew(iY)  = GVECm(iY)
        Unew(iE)  = GVECm(iE)
        Jnew(:,1) = GVECm(OS_1+1:OS_1+nE_G)
        Jnew(:,2) = GVECm(OS_2+1:OS_2+nE_G)

        IF( mk == M )THEN

          GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
          FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

        END IF


        ! --- Update Matter ---

        Y = Unew(iY) * Yold
        E = Unew(iE) * Eold

        CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
               ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

        CALL TimersStop( Timer_Im_UpdateFP )

      END IF



    END DO

    J = Jnew

         ! PRINT*
         ! PRINT*, "  dT = ", (T-Told) / Told
         ! PRINT*, "  dE = ", (E-Eold) / Eold
         ! PRINT*, "  dY = ", (Y-Yold) / Yold
         ! !PRINT*, "  dJNe = ", ENORM(J(:,1)-Jold(:,1)) / ENORM(Jold(:,1))
         ! !PRINT*, "  dJANe = ", ENORM(J(:,2)-Jold(:,2)) / ENORM(Jold(:,2))
         ! PRINT*, "  dJNe = ", WNORM(J(:,1)-Jold(:,1), W2_N) / WNORM(Jold(:,1), W2_N)
         ! PRINT*, "  dJANe = ", WNORM(J(:,2)-Jold(:,2), W2_N) / WNORM(Jold(:,2), W2_N)
         ! PRINT*, "  JANe0 = ", Jold(:,2)
         ! PRINT*, "  JANe = ", J(:,2)
         ! PRINT*, "  diffJANe0 = ", J(:,2) - Jold(:,2)
         ! PRINT*, "ChiAN =", Chi(:,2)
         ! PRINT*, "ChiNESAN =", Chi_NES(:,2)
         ! PRINT*, "ChiPairAN =", Chi_Pair(:,2)
         ! PRINT*, "SUMChiAN =", ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) )
         !
         ! !PRINT*, "  JNe0 = ", WNORM(Jold(:,1), W2_N)
         ! !PRINT*, "  JANe0 = ", WNORM(Jold(:,2), W2_N)
         ! !PRINT*, "  JNe0 = ", Jold(:,1)
         ! !PRINT*, "  W2_N = ", W2_N
         ! !PRINT*, "  W2_N*JNe0 = ", Jold(:,1) * W2_N
         ! ! PRINT*, "  JNe0 = ", ENORM(Jold(:,1))
         ! ! PRINT*, "  JANe0 = ", ENORM(Jold(:,2))
         ! PRINT*
         ! PRINT*, "SolveMatterEquations_FP_Coupled Done"
         ! PRINT*

    !
    ! PRINT*, "  J_Ne = [", Jnew(:,1), "  ];"
    ! PRINT*, "  J_ANe = [", Jnew(:,2), "  ];"

  END SUBROUTINE SolveMatterEquations_FP_Coupled


  SUBROUTINE SolveMatterEquations_FP_NestedAA &
    ( dt, iS_1, iS_2, J, Chi, J0, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      D, T, Y, E, nIterations_Inner, nIterations_Outer )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    INTEGER,  INTENT(out)   :: nIterations_Inner
    INTEGER,  INTENT(out)   :: nIterations_Outer

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M_inner = 3
    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * MAX(M_inner,2)
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED, CONVERGED_INNER
    INTEGER  :: i, k, k_inner, mk, mk_inner, INFO
    INTEGER  :: OS_1, OS_2
    REAL(DP) :: h3, N_B
    REAL(DP) :: S_Y, S_E, Yold, Eold
    REAL(DP) :: TMP(1)
    REAL(DP) :: C(2), Unew(2), Jnorm(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: Alpha(1:M_inner)
    REAL(DP) :: WORK(1:LWORK)
    REAL(DP) :: GVECm(2)
    ! REAL(DP) :: FVECm(2) !NEW: removed error vector
    REAL(DP) :: GVEC (1:2,1:M)
    REAL(DP) :: FVEC (1:2,1:M)
    REAL(DP) :: GVECm_inner (1:2*nE_G)
    ! REAL(DP) :: FVECm_inner (1:2*nE_G) !NEW: removed error vector
    REAL(DP) :: GVEC_inner (1:2*nE_G,1:M_inner)
    REAL(DP) :: FVEC_inner (1:2*nE_G,1:M_inner)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta(1:nE_G,1:2)
    REAL(DP) :: BVEC(2)
    REAL(DP) :: AMAT(1:2,1:M)
    REAL(DP) :: BVEC_inner(1:2*nE_G)
    REAL(DP) :: AMAT_inner(1:2*nE_G,1:M_inner)
    REAL(DP) :: Phi_0_In_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_In_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:2)

    OS_1 = 0
    OS_2 = nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J, Chi, J0, D, T, Y, E)

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF

    S_Y = N_B * Yold
    S_E = D   * Eold

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E

    Unew(iY) = Y / Yold ! --- Initial Guess
    Unew(iE) = E / Eold ! --- Initial Guess
    Jnew     = J        ! --- Initial Guess



    k = 0
    nIterations_Inner = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )

      k_inner = 0
      CONVERGED_INNER = .FALSE.

      CALL TimersStart( Timer_Im_ComputeOpacity )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

      ! --- NES Kernels ---

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
               Phi_0_In_NES(:,:,1), Phi_0_Ot_NES(:,:,1) )

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
               Phi_0_In_NES(:,:,2), Phi_0_Ot_NES(:,:,2) )

      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacities_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
               Phi_0_In_Pair(:,:,1), Phi_0_Ot_Pair(:,:,1) )

      CALL ComputeNeutrinoOpacities_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
               Phi_0_In_Pair(:,:,2), Phi_0_Ot_Pair(:,:,2) )

      CALL TimersStop( Timer_Im_ComputeOpacity )

      ! --- Emissivities and Opacities ---

      Eta(:,1) = Chi(:,1) * J0(:,1)
      Eta(:,2) = Chi(:,2) * J0(:,2)

      ! Calculate norm of the initial J for inner loop
      Jnorm(1) = WNORM( Jnew(:,1), W2_N )
      Jnorm(2) = WNORM( Jnew(:,2), W2_N )

      ! START INNER LOOP
      CALL TimersStart(Timer_Im_NestInner)

      DO WHILE( .NOT. CONVERGED_INNER .AND. k_inner < MaxIter )

        k_inner  = k_inner + 1
        mk_inner = MIN( M_inner, k_inner )

        CALL TimersStart( Timer_Im_ComputeRate )

        ! --- NES Emissivities and Opacities ---

        ! --- Neutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                    W2_N * Jnew(:,1),       1, Zero, Eta_NES(:,1), 1 )

        Chi_NES(:,1) = Eta_NES(:,1)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                    W2_N * (One-Jnew(:,1)), 1, One,  Chi_NES(:,1), 1 )

        ! --- Antineutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                    W2_N * Jnew(:,2),       1, Zero, Eta_NES(:,2), 1 )

        Chi_NES(:,2) = Eta_NES(:,2)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                    W2_N * (One-Jnew(:,2)), 1, One,  Chi_NES(:,2), 1 )

        ! --- Pair Emissivities and Opacities ---

        ! --- Neutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                    W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

        Chi_Pair(:,1) = Eta_Pair(:,1)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                    W2_N * Jnew(:,2),       1, One,  Chi_Pair(:,1), 1 )

        ! --- Antineutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                    W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

        Chi_Pair(:,2) = Eta_Pair(:,2)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                    W2_N * Jnew(:,1),       1, One,  Chi_Pair(:,2), 1 )

        CALL TimersStop( Timer_Im_ComputeRate )


        ! --- Right-Hand Side Vectors (inner) ---

        GVEC_inner(OS_1+1:OS_1+nE_G,mk_inner) &
          = ( Jold(:,1) + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) &
                / ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) )

        GVEC_inner(OS_2+1:OS_2+nE_G,mk_inner) &
          = ( Jold(:,2) + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) &
                / ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) )

        ! --- Residuals (inner)---

        DO i = 1, nE_G
          FVEC_inner(OS_1+i,mk_inner) = GVEC_inner(OS_1+i,mk_inner) - Jnew(i,1)
        END DO

        DO i = 1, nE_G
          FVEC_inner(OS_2+i,mk_inner) = GVEC_inner(OS_2+i,mk_inner) - Jnew(i,2)
        END DO


        IF( WNORM( FVEC_inner(OS_1+1:OS_1+nE_G, mk_inner) * &
            ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) ), W2_N ) &
            <= Rtol * Jnorm(1) .AND. &
            WNORM( FVEC_inner(OS_2+1:OS_2+nE_G, mk_inner) * &
            ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) ), W2_N ) &
            <= Rtol * Jnorm(2) ) &
        THEN

          CONVERGED_INNER = .TRUE.

          nIterations_Inner = nIterations_Inner + k_inner

        END IF

        IF ( .NOT. CONVERGED_INNER )THEN

          CALL TimersStart( Timer_Im_ComputeLS )

          IF( mk_inner == 1 )THEN

            GVECm_inner = GVEC_inner(:,mk_inner)

          ELSE

            BVEC_inner(:) &
              = - FVEC_inner(:,mk_inner)
            AMAT_inner(:,1:mk_inner-1) &
              = FVEC_inner(:,1:mk_inner-1) - SPREAD( FVEC_inner(:,mk_inner), DIM = 2, NCOPIES = mk_inner-1 )

            CALL DGELS( 'N', 2*nE_G, mk_inner-1, 1, AMAT_inner(:,1:mk_inner-1), 2*nE_G, &
                        BVEC_inner, 2*nE_G, WORK, LWORK, INFO )

            Alpha(1:mk_inner-1) = BVEC_inner(1:mk_inner-1)

            Alpha(mk_inner) = One - SUM( Alpha(1:mk_inner-1) )
            ! PRINT*, "Alpha = ", Alpha

            GVECm_inner = Zero
            DO i = 1, mk_inner

              GVECm_inner = GVECm_inner + Alpha(i) * GVEC_inner(:,i)

            END DO

          END IF

          CALL TimersStop( Timer_Im_ComputeLS )

          CALL TimersStart( Timer_Im_UpdateFP )

          Jnew(:,1) = GVECm_inner(OS_1+1:OS_1+nE_G)
          Jnew(:,2) = GVECm_inner(OS_2+1:OS_2+nE_G)

          IF( mk_inner == M_inner )THEN

            GVEC_inner = CSHIFT( GVEC_inner, SHIFT = + 1, DIM = 2 )
            FVEC_inner = CSHIFT( FVEC_inner, SHIFT = + 1, DIM = 2 )

          END IF

          CALL TimersStop( Timer_Im_UpdateFP )
        END IF

      END DO ! END OF INNER LOOP

      CALL TimersStop(Timer_Im_NestInner)

      ! --- Right-Hand Side Vectors ---

      GVEC(iY,mk) = One + C(iY) &
                      - DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y
      GVEC(iE,mk) = One + C(iE) &
                      - DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      ! --- Residuals ---

      FVEC(iY,mk) = GVEC(iY,mk) - Unew(iY)

      FVEC(iE,mk) = GVEC(iE,mk) - Unew(iE)

      IF( ENORM( [ FVEC(iY,mk) ] ) <= Rtol .AND. &
          ENORM( [ FVEC(iE,mk) ] ) <= Rtol ) &
      THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

        nIterations_Outer = k
        nIterations_Inner &
          = FLOOR( DBLE( nIterations_Inner ) / DBLE( k ) )

      END IF

      IF ( .NOT. CONVERGED )THEN

      CALL TimersStart( Timer_Im_ComputeLS )

        IF( mk == 1 )THEN

          GVECm = GVEC(:,mk)

        ELSEIF (mk == 2)THEN
          BVEC(:) = - FVEC(:,mk)
          AMAT(:,1) = FVEC(:,1) - FVEC(:,mk)

          Alpha(1) = DOT_PRODUCT( AMAT(:,1), BVEC ) / DOT_PRODUCT( AMAT(:,1), AMAT(:,1) )
          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        ELSEIF (mk == 3)THEN
          BVEC(:) = - FVEC(:,mk)
          AMAT(:,1:mk-1) = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = 2 )


          Alpha(1) = AMAT(2,mk-1) * BVEC(1) - AMAT(1,mk-1) * BVEC(2)
          Alpha(mk-1) = - AMAT(2,1) * BVEC(1) + AMAT(1,1) * BVEC(2)
          Alpha(1:mk-1) = Alpha(1:mk-1) / (AMAT(1,1)*AMAT(2,mk-1) - AMAT(1,mk-1)*AMAT(2,1))

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        ELSE
          ! SOLVING AN UNDERDETERMINED SYSTEM
          BVEC(:) &
            = - FVEC(:,mk)
          AMAT(:,1:mk-1) &
            = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

          CALL DGELS( 'N', 2, mk-1, 1, AMAT(:,1:mk-1), 2, &
                      BVEC, mk, WORK, LWORK, INFO )

          Alpha(1:mk-1) = BVEC(1:mk-1)

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        END IF

        CALL TimersStop( Timer_Im_ComputeLS )

        CALL TimersStart( Timer_Im_UpdateFP )

        Unew(iY)  = GVECm(iY)
        Unew(iE)  = GVECm(iE)

        IF( mk == M )THEN

          GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
          FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

        END IF

        ! --- Update Matter ---

        Y = Unew(iY) * Yold
        E = Unew(iE) * Eold

        ! PRINT*, "U = ", Unew

        CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
               ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

        CALL TimersStop( Timer_Im_UpdateFP )
      END IF

    END DO

    J = Jnew


  END SUBROUTINE SolveMatterEquations_FP_NestedAA


  SUBROUTINE SolveMatterEquations_FP_NestedNewton &
    ( dt, iS_1, iS_2, J, Chi, J0, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      D, T, Y, E, nIterations_Inner, nIterations_Outer )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    INTEGER,  INTENT(out)   :: nIterations_Inner
    INTEGER,  INTENT(out)   :: nIterations_Outer

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * 2
    REAL(DP), PARAMETER :: Rtol = 1.0d-8
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED, CONVERGED_INNER
    INTEGER  :: i, l, k, k_inner, mk, INFO
    INTEGER  :: OS_1, OS_2
    REAL(DP) :: h3, N_B
    REAL(DP) :: S_Y, S_E, Yold, Eold
    REAL(DP) :: TMP(1)
    REAL(DP) :: C(2), Unew(2), Jnorm(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: Alpha(1:M)
    REAL(DP) :: WORK(1:LWORK)
    REAL(DP) :: IPIV(1:nE_G)
    REAL(DP) :: GVECm(2)
    ! REAL(DP) :: FVECm(2) !NEW: removed error vector
    REAL(DP) :: GVEC (1:2,1:M)
    REAL(DP) :: FVEC (1:2,1:M)
    REAL(DP) :: GVECm_inner(1:(2*nE_G))
    ! REAL(DP) :: FVECm_inner (1:2*nE_G) !NEW: removed error vector
    REAL(DP) :: GVEC_inner (1:(2*nE_G))
    REAL(DP) :: GJAC_inner (1:(2*nE_G),1:(2*nE_G))
    REAL(DP) :: S_J (1:nE_G,1:2)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta(1:nE_G,1:2)
    REAL(DP) :: BVEC(2)
    REAL(DP) :: AMAT(1:2,1:M)
    REAL(DP) :: Phi_0_In_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_In_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:2)

    OS_1 = 0
    OS_2 = nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J, Chi, J0, D, T, Y, E)

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF


    S_Y = N_B * Yold
    S_E = D   * Eold
    S_J = One / ( One + dt * Chi )

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E


    Unew(iY) = Y / Yold ! --- Initial Guess
    Unew(iE) = E / Eold ! --- Initial Guess
    Jnew     = J        ! --- Initial Guess


    k = 0
    nIterations_Inner = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )

      k_inner = 0
      CONVERGED_INNER = .FALSE.

      CALL TimersStart( Timer_Im_ComputeOpacity )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

      CALL ComputeEquilibriumDistributions_Point &
             ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

      ! --- NES Kernels ---

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
               Phi_0_In_NES(:,:,1), Phi_0_Ot_NES(:,:,1) )

      CALL ComputeNeutrinoOpacities_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
               Phi_0_In_NES(:,:,2), Phi_0_Ot_NES(:,:,2) )

      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacities_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
               Phi_0_In_Pair(:,:,1), Phi_0_Ot_Pair(:,:,1) )

      CALL ComputeNeutrinoOpacities_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
               Phi_0_In_Pair(:,:,2), Phi_0_Ot_Pair(:,:,2) )

      CALL TimersStop( Timer_Im_ComputeOpacity )

      ! --- Emissivities and Opacities ---

      Eta(:,1) = Chi(:,1) * J0(:,1)
      Eta(:,2) = Chi(:,2) * J0(:,2)

      ! Calculate norm of the initial J for inner loop
      Jnorm(1) = WNORM( Jnew(:,1), W2_N )
      Jnorm(2) = WNORM( Jnew(:,2), W2_N )

      ! START INNER LOOP
      CALL TimersStart(Timer_Im_NestInner)

      DO WHILE( .NOT. CONVERGED_INNER .AND. k_inner < MaxIter )

        k_inner = k_inner + 1

        CALL TimersStart( Timer_Im_ComputeRate )

        ! --- NES Emissivities and Opacities ---

        ! --- Neutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                   W2_N * Jnew(:,1),       1, Zero, Eta_NES(:,1), 1 )

        Chi_NES(:,1) = Eta_NES(:,1)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                   W2_N * (One-Jnew(:,1)), 1, One,  Chi_NES(:,1), 1 )

        ! --- Antineutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                   W2_N * Jnew(:,2),       1, Zero, Eta_NES(:,2), 1 )

        Chi_NES(:,2) = Eta_NES(:,2)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                   W2_N * (One-Jnew(:,2)), 1, One,  Chi_NES(:,2), 1 )

        ! --- Pair Emissivities and Opacities ---

        ! --- Neutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                   W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

        Chi_Pair(:,1) = Eta_Pair(:,1)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                   W2_N * Jnew(:,2),       1, One,  Chi_Pair(:,1), 1 )

        ! --- Antineutrino ---

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                   W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

        Chi_Pair(:,2) = Eta_Pair(:,2)

        CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                   W2_N * Jnew(:,1),       1, One,  Chi_Pair(:,2), 1 )

        CALL TimersStop( Timer_Im_ComputeRate )

        ! IF( k == 1 )THEN
        !
        !   ! TEST: --- Update Initial Guess for Neutrinos and Antineutrinos ---
        !
        !   Jnew = ( Jold + dt * ( Eta + Eta_NES + Eta_Pair ) ) &
        !          / ( One + dt * ( Chi + Chi_NES + Chi_Pair ) )
        !
        ! END IF

        ! --- Right-Hand Side Vectors (inner) ---

        GVEC_inner(OS_1+1:OS_1+nE_G) &
          = ( One + dt * ( Chi_NES(:,1) + Chi_Pair(:,1) ) * S_J(:,1) ) &
              * Jnew(:,1) - ( Jold(:,1) &
              + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) * S_J(:,1)

        GVEC_inner(OS_2+1:OS_2+nE_G) &
          = ( One + dt * ( Chi_NES(:,2) + Chi_Pair(:,2) ) * S_J(:,2) ) &
              * Jnew(:,2) - ( Jold(:,2) &
              + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) * S_J(:,2)

        IF( WNORM( GVEC_inner(OS_1+1:OS_1+nE_G) / S_J(:,1), W2_N ) &
            <= Rtol * Jnorm(1) .AND. &
            WNORM( GVEC_inner(OS_2+1:OS_2+nE_G) / S_J(:,2), W2_N ) &
            <= Rtol * Jnorm(2) ) &
        THEN

          CONVERGED_INNER = .TRUE.

          nIterations_Inner = nIterations_Inner + k_inner

        END IF
        IF ( .NOT. CONVERGED_INNER )THEN

          CALL TimersStart( Timer_Im_ComputeLS )

          ! --- Jacobian matrix ---

          ! --- NES terms ---

          DO i = 1, nE_G
          DO l = 1, nE_G
            GJAC_inner(OS_1+i,OS_1+l) &
              = - dt * ( Phi_0_In_NES(l,i,1) * ( One - Jnew(i,1) ) * W2_N(l) &
                         + Phi_0_Ot_NES(l,i,1) * Jnew(i,1) * W2_N(l) ) * S_J(i,1)
          END DO
          END DO

          DO i = 1, nE_G
          DO l = 1, nE_G
            GJAC_inner(OS_2+i,OS_2+l) &
              = - dt * ( Phi_0_In_NES(l,i,2) * ( One - Jnew(i,2) ) * W2_N(l) &
                         + Phi_0_Ot_NES(l,i,2) * Jnew(i,2) * W2_N(l) ) * S_J(i,2)
          END DO
          END DO

          ! --- Pair terms ---

          DO i = 1, nE_G
          DO l = 1, nE_G
            GJAC_inner(OS_1+i,OS_2+l) &
              = dt * ( Phi_0_In_Pair(l,i,1) * ( One - Jnew(i,1) ) * W2_N(l) &
                       + Phi_0_Ot_Pair(l,i,1) * Jnew(i,1) * W2_N(l) ) * S_J(i,1)
          END DO
          END DO

          DO i = 1, nE_G
          DO l = 1, nE_G
            GJAC_inner(OS_2+i,OS_1+l) &
              = dt * ( Phi_0_In_Pair(l,i,2) * ( One - Jnew(i,2) ) * W2_N(l) &
                       + Phi_0_Ot_Pair(l,i,2) * Jnew(i,2) * W2_N(l) ) * S_J(i,2)
          END DO
          END DO

          ! --- diagonal terms ---

          DO i = 1, nE_G
            GJAC_inner(OS_1+i,OS_1+i) &
              = GJAC_inner(OS_1+i,OS_1+i) &
                  + One + dt * ( Chi_NES(i,1) + Chi_Pair(i,1) ) * S_J(i,1)
          END DO

          DO i = 1, nE_G
            GJAC_inner(OS_2+i,OS_2+i) &
              = GJAC_inner(OS_2+i,OS_2+i) &
                  + One + dt * ( Chi_NES(i,2) + Chi_Pair(i,2) ) * S_J(i,2)
          END DO

          CALL DGESV( 2*nE_G, 1, GJAC_inner, 2*nE_G, IPIV, GVEC_inner, 2*nE_G, INFO )

          CALL TimersStop( Timer_Im_ComputeLS )

          GVECm_inner(OS_1+1:OS_1+nE_G) = Jnew(:,1) - GVEC_inner(OS_1+1:OS_1+nE_G)
          GVECm_inner(OS_2+1:OS_2+nE_G) = Jnew(:,2) - GVEC_inner(OS_2+1:OS_2+nE_G)


          Jnew(:,1) = GVECm_inner(OS_1+1:OS_1+nE_G)
          Jnew(:,2) = GVECm_inner(OS_2+1:OS_2+nE_G)
        END IF

      END DO ! END OF INNER LOOP

      CALL TimersStop(Timer_Im_NestInner)

      ! --- Right-Hand Side Vectors ---

      GVEC(iY,mk) = One + C(iY) &
                      - DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y
      GVEC(iE,mk) = One + C(iE) &
                      - DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      ! --- Residuals ---

      FVEC(iY,mk) = GVEC(iY,mk) - Unew(iY)

      FVEC(iE,mk) = GVEC(iE,mk) - Unew(iE)

      IF( ENORM( [ FVEC(iY,mk) ] ) <= Rtol .AND. &
          ENORM( [ FVEC(iE,mk) ] ) <= Rtol ) &
      THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

        nIterations_Outer = k

        nIterations_Inner &
          = FLOOR( DBLE( nIterations_Inner ) / DBLE( k ) )

      END IF

      IF ( .NOT. CONVERGED )THEN

        CALL TimersStart( Timer_Im_ComputeLS )

        IF( mk == 1 )THEN

          GVECm = GVEC(:,mk)

        ELSEIF( mk == 2 )THEN

          BVEC(:) = - FVEC(:,mk)
          AMAT(:,1) = FVEC(:,1) - FVEC(:,mk)

          Alpha(1) = DOT_PRODUCT( AMAT(:,1), BVEC ) / DOT_PRODUCT( AMAT(:,1), AMAT(:,1) )
          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        ELSEIF( mk == 3 )THEN

          BVEC(:) = - FVEC(:,mk)
          AMAT(:,1:mk-1) = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = 2 )

          Alpha(1) = AMAT(2,mk-1) * BVEC(1) - AMAT(1,mk-1) * BVEC(2)
          Alpha(mk-1) = - AMAT(2,1) * BVEC(1) + AMAT(1,1) * BVEC(2)
          Alpha(1:mk-1) = Alpha(1:mk-1) / (AMAT(1,1)*AMAT(2,mk-1) - AMAT(1,mk-1)*AMAT(2,1))

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        ELSE

          ! SOLVING AN UNDERDETERMINED SYSTEM
          BVEC(:) &
            = - FVEC(:,mk)
          AMAT(:,1:mk-1) &
            = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

          CALL DGELS( 'N', 2, mk-1, 1, AMAT(:,1:mk-1), 2, BVEC, mk, WORK, LWORK, INFO )

          Alpha(1:mk-1) = BVEC(1:mk-1)

          Alpha(mk) = One - SUM( Alpha(1:mk-1) )

          GVECm = Zero
          DO i = 1, mk

            GVECm = GVECm + Alpha(i) * GVEC(:,i)

          END DO

        END IF

        CALL TimersStop( Timer_Im_ComputeLS )

        CALL TimersStart( Timer_Im_UpdateFP )


        Unew(iY)  = GVECm(iY)
        Unew(iE)  = GVECm(iE)

        IF( mk == M )THEN

          GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
          FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

        END IF

        ! --- Update Matter ---

        Y = Unew(iY) * Yold
        E = Unew(iE) * Eold

        CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
               ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

        CALL TimersStop( Timer_Im_UpdateFP )
      END IF

     END DO

     J = Jnew

     ! PRINT*
     ! PRINT*, "SolveMatterEquations_Nested_Newton Output"
     ! PRINT*
     ! PRINT*, "  D = ", D / Unit_D
     ! PRINT*, "  T = ", T / Kelvin
     ! PRINT*, "  E = ", E / Unit_E
     ! PRINT*, "  Y = ", Y
     !
     !
     ! PRINT*, "  J_Ne = [", Jnew(:,1), "  ];"
     ! PRINT*, "  J_ANe = [", Jnew(:,2), "  ];"
     ! stop "stop after printing"

  END SUBROUTINE SolveMatterEquations_FP_NestedNewton


  SUBROUTINE SolveMatterEquations_Newton &
    ( dt, iS_1, iS_2, J, Chi, J0, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E
    INTEGER,  INTENT(out)   :: nIterations

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1, iE = 2
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, l, m
    INTEGER  :: OS_1, OS_2
    INTEGER  :: INFO
    INTEGER  :: IPIV(2+2*nE_G)
    REAL(DP) :: h3, N_B
    REAL(DP) :: Yold, Eold, S_Y, S_E, C_Y, C_E
    REAL(DP) :: Ynew, Enew, TMP(1)
    REAL(DP) :: Jnorm(2)
    REAL(DP) :: W2_S(1:nE_G)
    REAL(DP) :: W3_S(1:nE_G)
    REAL(DP) :: UVEC(2+2*nE_G)
    REAL(DP) :: FVEC(2+2*nE_G)
    REAL(DP) :: DVEC(2+2*nE_G)
    REAL(DP) :: FJAC(2+2*nE_G,2+2*nE_G)
    REAL(DP) :: S_J (1:nE_G,1:2)
    REAL(DP) :: Jold(1:nE_G,1:2)
    REAL(DP) :: Jnew(1:nE_G,1:2)
    REAL(DP) :: Eta (1:nE_G,1:2)
    REAL(DP) :: dEtadY      (1:nE_G,1:2), dEtadE      (1:nE_G,1:2)
    REAL(DP) :: dJ0dY       (1:nE_G,1:2), dJ0dE       (1:nE_G,1:2)
    REAL(DP) :: dChi_NES_dY (1:nE_G,1:2), dChi_NES_dE (1:nE_G,1:2)
    REAL(DP) :: dEta_NES_dY (1:nE_G,1:2), dEta_NES_dE (1:nE_G,1:2)
    REAL(DP) :: dChi_Pair_dY(1:nE_G,1:2), dChi_Pair_dE(1:nE_G,1:2)
    REAL(DP) :: dEta_Pair_dY(1:nE_G,1:2), dEta_Pair_dE(1:nE_G,1:2)
    REAL(DP) ::  Phi_0_In_NES    (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_In_NES_dY (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_In_NES_dE (1:nE_G,1:nE_G,1:2)
    REAL(DP) ::  Phi_0_Ot_NES    (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_Ot_NES_dY (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_Ot_NES_dE (1:nE_G,1:nE_G,1:2)
    REAL(DP) ::  Phi_0_In_Pair   (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_In_Pair_dY(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_In_Pair_dE(1:nE_G,1:nE_G,1:2)
    REAL(DP) ::  Phi_0_Ot_Pair   (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_Ot_Pair_dY(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: dPhi_0_Ot_Pair_dE(1:nE_G,1:nE_G,1:2)

    OS_1 = 2
    OS_2 = 2 + nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E
    Jold = J

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J, Chi, J0, D, T, Y, E, TOL = Rtol )

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF

    S_Y = N_B * Yold
    S_E = D   * Eold

    S_J = One / ( One + dt * Chi )

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    C_Y = - ( One + DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y )
    C_E = - ( One + DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E )

    Jnorm(1) = WNORM(Jold(:,1), W2_N)
    Jnorm(2) = WNORM(Jold(:,2), W2_N)

    Ynew = Y ! --- Initial Guess
    Enew = E ! --- Initial Guess
    Jnew = J ! --- Initial Guess


    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k = k + 1

      CALL TimersStart( Timer_Im_ComputeOpacity )

      ! --- Equilibrium Distributions ---

      CALL ComputeEquilibriumDistributionAndDerivatives_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, J0(:,1), dJ0dY(:,1), dJ0dE(:,1) )

      CALL ComputeEquilibriumDistributionAndDerivatives_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, J0(:,2), dJ0dY(:,2), dJ0dE(:,2) )

      ! PRINT*, " J0 done "
      ! PRINT*, "  T = ", (T/Kelvin)
      ! PRINT*, "  E = ", (E/MeV)
      ! PRINT*, "  Y = ", (Y)
      ! PRINT*, "  Jinit = ", Jold(:,1)
      ! PRINT*, "  JANinit = ", Jold(:,2)
      ! PRINT*, "  J = ", Jnew(:,1)
      ! PRINT*, "  JAN = ", Jnew(:,2)
      ! --- NES Kernels ---

      CALL ComputeNeutrinoOpacitiesAndDerivatives_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
                Phi_0_In_NES   (:,:,1),  Phi_0_Ot_NES   (:,:,1), &
               dPhi_0_In_NES_dY(:,:,1), dPhi_0_In_NES_dE(:,:,1), &
               dPhi_0_Ot_NES_dY(:,:,1), dPhi_0_Ot_NES_dE(:,:,1) )

      ! PRINT*, " NES1 done "
      ! PRINT*

      CALL ComputeNeutrinoOpacitiesAndDerivatives_NES_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
                Phi_0_In_NES   (:,:,2),  Phi_0_Ot_NES   (:,:,2), &
               dPhi_0_In_NES_dY(:,:,2), dPhi_0_In_NES_dE(:,:,2), &
               dPhi_0_Ot_NES_dY(:,:,2), dPhi_0_Ot_NES_dE(:,:,2) )

      ! PRINT*, " NES2 done "
      ! PRINT*
      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
                Phi_0_In_Pair   (:,:,1),  Phi_0_Ot_Pair   (:,:,1), &
               dPhi_0_In_Pair_dY(:,:,1), dPhi_0_In_Pair_dE(:,:,1), &
               dPhi_0_Ot_Pair_dY(:,:,1), dPhi_0_Ot_Pair_dE(:,:,1) )

      CALL ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point &
             ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
                Phi_0_In_Pair   (:,:,2),  Phi_0_Ot_Pair   (:,:,2), &
               dPhi_0_In_Pair_dY(:,:,2), dPhi_0_In_Pair_dE(:,:,2), &
               dPhi_0_Ot_Pair_dY(:,:,2), dPhi_0_Ot_Pair_dE(:,:,2) )

      ! PRINT*, " PAIR done "
      ! PRINT*

      CALL TimersStop( Timer_Im_ComputeOpacity )


      CALL TimersStart( Timer_Im_ComputeRate )

      ! --- Electron Capture Emissivities ---

      Eta    = Chi * J0
      dEtadY = Chi * dJ0dY
      dEtadE = Chi * dJ0dE

      ! --- NES Emissivities and Opacities ---

      ! --- Neutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                  W2_N * Jnew(:,1), 1, Zero, Eta_NES(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_NES_dY(:,:,1), nE_G, &
                  W2_N * Jnew(:,1), 1, Zero, dEta_NES_dY(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_NES_dE(:,:,1), nE_G, &
                  W2_N * Jnew(:,1), 1, Zero, dEta_NES_dE(:,1), 1 )

       Chi_NES   (:,1) = Eta_NES    (:,1)
      dChi_NES_dY(:,1) = dEta_NES_dY(:,1)
      dChi_NES_dE(:,1) = dEta_NES_dE(:,1)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, One, Chi_NES(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_NES_dY(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, One, dChi_NES_dY(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_NES_dE(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, One, dChi_NES_dE(:,1), 1 )

      ! --- Antineutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                  W2_N * Jnew(:,2), 1, Zero, Eta_NES(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_NES_dY(:,:,2), nE_G, &
                  W2_N * Jnew(:,2), 1, Zero, dEta_NES_dY(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_NES_dE(:,:,2), nE_G, &
                  W2_N * Jnew(:,2), 1, Zero, dEta_NES_dE(:,2), 1 )

       Chi_NES   (:,2) = Eta_NES    (:,2)
      dChi_NES_dY(:,2) = dEta_NES_dY(:,2)
      dChi_NES_dE(:,2) = dEta_NES_dE(:,2)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, One, Chi_NES(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_NES_dY(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, One, dChi_NES_dY(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_NES_dE(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, One, dChi_NES_dE(:,2), 1 )

      ! --- Pair Emissivities and Opacities ---

      ! --- Neutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_Pair_dY(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, Zero, dEta_Pair_dY(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_Pair_dE(:,:,1), nE_G, &
                  W2_N * (One-Jnew(:,2)), 1, Zero, dEta_Pair_dE(:,1), 1 )

       Chi_Pair   (:,1) =  Eta_Pair   (:,1)
      dChi_Pair_dY(:,1) = dEta_Pair_dY(:,1)
      dChi_Pair_dE(:,1) = dEta_Pair_dE(:,1)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                  W2_N * Jnew(:,2), 1, One, Chi_Pair(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_Pair_dY(:,:,1), nE_G, &
                  W2_N * Jnew(:,2), 1, One, dChi_Pair_dY(:,1), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_Pair_dE(:,:,1), nE_G, &
                  W2_N * Jnew(:,2), 1, One, dChi_Pair_dE(:,1), 1 )

      ! --- Antineutrino ---

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_Pair_dY(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, Zero, dEta_Pair_dY(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_In_Pair_dE(:,:,2), nE_G, &
                  W2_N * (One-Jnew(:,1)), 1, Zero, dEta_Pair_dE(:,2), 1 )

       Chi_Pair   (:,2) =  Eta_Pair   (:,2)
      dChi_Pair_dY(:,2) = dEta_Pair_dY(:,2)
      dChi_Pair_dE(:,2) = dEta_Pair_dE(:,2)

      CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                  W2_N * Jnew(:,1), 1, One, Chi_Pair(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_Pair_dY(:,:,2), nE_G, &
                  W2_N * Jnew(:,1), 1, One, dChi_Pair_dY(:,2), 1 )

      CALL DGEMV( 'T', nE_G, nE_G, One, dPhi_0_Ot_Pair_dE(:,:,2), nE_G, &
                  W2_N * Jnew(:,1), 1, One, dChi_Pair_dE(:,2), 1 )

      CALL TimersStop( Timer_Im_ComputeRate )

      ! IF( k == 1 )THEN
      !
      !   ! TEST:--- Update Initial Guess for Neutrinos and Antineutrinos ---
      !
      !   Jnew = ( Jold + dt * ( Eta + Eta_NES + Eta_Pair ) ) &
      !          / ( One + dt * ( Chi + Chi_NES + Chi_Pair ) )
      !
      ! END IF

      CALL TimersStart( Timer_Im_ComputeLS )

      ! --- Solution Vector ---

      UVEC = Zero

      UVEC(iY) = Ynew
      UVEC(iE) = Enew
      UVEC(OS_1+1:OS_1+nE_G) = Jnew(:,1)
      UVEC(OS_2+1:OS_2+nE_G) = Jnew(:,2)

      ! --- Equation Vector ---

      FVEC = Zero

      ! --- Electron Fraction Equation ---

      FVEC(iY) &
        = Ynew / Yold + C_Y &
            + DOT_PRODUCT( W2_S, Jnew(:,1) - Jnew(:,2) ) / S_Y

      ! --- Internal Energy Equation ---

      FVEC(iE) &
        = Enew / Eold + C_E &
            + DOT_PRODUCT( W3_S, Jnew(:,1) + Jnew(:,2) ) / S_E

      ! --- Neutrino Equation ---

      FVEC(OS_1+1:OS_1+nE_G) &
        = ( One + dt * ( Chi_NES(:,1) + Chi_Pair(:,1) ) * S_J(:,1) ) * Jnew(:,1) &
          - ( Jold(:,1) + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) * S_J(:,1)

      ! --- Antineutrino Equation ---

      FVEC(OS_2+1:OS_2+nE_G) &
        = ( One + dt * ( Chi_NES(:,2) + Chi_Pair(:,2) ) * S_J(:,2) ) * Jnew(:,2) &
          - ( Jold(:,2) + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) * S_J(:,2)

      CALL TimersStop( Timer_Im_ComputeLS )

      IF( ENORM( [ FVEC(iY) ] ) <= Rtol .AND. &
          ENORM( [ FVEC(iE) ] ) <= Rtol .AND. &
          WNORM( FVEC(OS_1+1:OS_1+nE_G) / S_J(:,1), W2_N ) &
          <= Rtol * Jnorm(1) .AND. &
          WNORM( FVEC(OS_2+1:OS_2+nE_G) / S_J(:,2), W2_N ) &
          <= Rtol * Jnorm(2) ) &
      THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

        nIterations = k

      END IF

      IF ( .NOT. CONVERGED )THEN

        CALL TimersStart( Timer_Im_ComputeLS )

        ! --- Jacobian Matrix (Transpose) ---

        FJAC = Zero

        ! --- Electron Fraction Equation Jacobian ---

        FJAC(iY,iY) = One / Yold
        FJAC(iE,iY) = Zero
        FJAC(OS_1+1:OS_1+nE_G,iY) =   W2_S / S_Y
        FJAC(OS_2+1:OS_2+nE_G,iY) = - W2_S / S_Y

        ! --- Internal Energy Equation Jacobian ---

        FJAC(iY,iE) = Zero
        FJAC(iE,iE) = One / Eold
        FJAC(OS_1+1:OS_1+nE_G,iE) = W3_S / S_Y
        FJAC(OS_2+1:OS_2+nE_G,iE) = W3_S / S_Y

        ! --- Neutrino Equation Jacobian ---

        DO l = 1, nE_G

          ! --- Derivative wrt. Y ---

          FJAC(iY,OS_1+l) &
            = FJAC(iY,OS_1+l) &
              + dt * ( ( dChi_NES_dY(l,1) + dChi_Pair_dY(l,1) ) * Jnew(l,1) &
                       - ( dEtadY(l,1) + dEta_NES_dY(l,1) + dEta_Pair_dY(l,1) ) ) * S_J(l,1)

          ! --- Derivative wrt. E ---

          FJAC(iE,OS_1+l) &
            = FJAC(iE,OS_1+l) &
              + dt * ( ( dChi_NES_dE(l,1) + dChi_Pair_dE(l,1) ) * Jnew(l,1) &
                       - ( dEtadE(l,1) + dEta_NES_dE(l,1) + dEta_Pair_dE(l,1) ) ) * S_J(l,1)

          ! --- Derivative wrt. J_1 ---

          FJAC(OS_1+l,OS_1+l) &
            = FJAC(OS_1+l,OS_1+l) &
              + One + dt * ( Chi_NES(l,1) + Chi_Pair(l,1) ) * S_J(l,1)

          DO m = 1, nE_G

            FJAC(OS_1+m,OS_1+l) &
              = FJAC(OS_1+m,OS_1+l) &
                + dt * W2_N(m) * ( ( Phi_0_In_NES(m,l,1) &
                                     - Phi_0_Ot_NES(m,l,1) ) * Jnew(l,1) &
                                   - Phi_0_In_NES(m,l,1) ) * S_J(l,1)

          END DO

          ! --- Derivative wrt. J_2 ---

          DO m = 1, nE_G

            FJAC(OS_2+m,OS_1+l) &
              = FJAC(OS_2+m,OS_1+l) &
                + dt * W2_N(m) * ( ( Phi_0_Ot_Pair(m,l,1) &
                                     - Phi_0_In_Pair(m,l,1) ) * Jnew(l,1) &
                                   + Phi_0_In_Pair(m,l,1) ) * S_J(l,1)

          END DO

        END DO

        ! --- Antineutrino Equation Jacobian ---

        DO l = 1, nE_G

          ! --- Derivative wrt. Y ---

          FJAC(iY,OS_2+l) &
            = FJAC(iY,OS_2+l) &
              + dt * ( ( dChi_NES_dY(l,2) + dChi_Pair_dY(l,2) ) * Jnew(l,2) &
                       - ( dEtadY(l,2) + dEta_NES_dY(l,2) + dEta_Pair_dY(l,2) ) ) * S_J(l,2)

          ! --- Derivative wrt. E ---

          FJAC(iE,OS_2+l) &
            = FJAC(iE,OS_2+l) &
              + dt * ( ( dChi_NES_dE(l,2) + dChi_Pair_dE(l,2) ) * Jnew(l,2) &
                       - ( dEtadE(l,2) + dEta_NES_dE(l,2) + dEta_Pair_dE(l,2) ) ) * S_J(l,2)

          ! --- Derivative wrt. J_2 ---

          FJAC(OS_2+l,OS_2+l) &
            = FJAC(OS_2+l,OS_2+l) &
              + One + dt * ( Chi_NES(l,2) + Chi_Pair(l,2) ) * S_J(l,2)

          DO m = 1, nE_G

            FJAC(OS_2+m,OS_2+l) &
              = FJAC(OS_2+m,OS_2+l) &
                + dt * W2_N(m) * ( ( Phi_0_In_NES(m,l,2) &
                                     - Phi_0_Ot_NES(m,l,2) ) * Jnew(l,2) &
                                   - Phi_0_In_NES(m,l,2) ) * S_J(l,2)

          END DO

          ! --- Derivative wrt. J_1 ---

          DO m = 1, nE_G

            FJAC(OS_1+m,OS_2+l) &
              = FJAC(OS_1+m,OS_2+l) &
                + dt * W2_N(m) * ( ( Phi_0_Ot_Pair(m,l,2) &
                                     - Phi_0_In_Pair(m,l,2) ) * Jnew(l,2) &
                                   + Phi_0_In_Pair(m,l,2) ) * S_J(l,2)

          END DO

        END DO

        FJAC = TRANSPOSE( FJAC )



        ! --- Solve Linear System ---

        DVEC = - FVEC

        CALL DGESV( 2+2*nE_G, 1, FJAC, 2+2*nE_G, IPIV, DVEC, 2+2*nE_G, INFO )

        CALL TimersStop( Timer_Im_ComputeLS )

        IF( INFO .NE. 0 )THEN
          PRINT*, "INFO = ", INFO
          STOP
        END IF


        CALL TimersStart( Timer_Im_UpdateFP )

        ! IF( ENORM( [ DVEC(iY) ] ) <= Rtol * ENORM( [ Yold ] ) .AND. &
        !     ENORM( [ DVEC(iE) ] ) <= Rtol * ENORM( [ Eold ] ) .AND. &
        !     ENORM( DVEC(OS_1+1:OS_1+nE_G) ) <= Rtol * ENORM( Jold(:,1) ) .AND. &
        !     ENORM( DVEC(OS_2+1:OS_2+nE_G) ) <= Rtol * ENORM( Jold(:,2) ) ) &
        ! THEN
        !
        !   CONVERGED = .TRUE.
        !
        !   Iterations_Min = MIN( Iterations_Min, k )
        !   Iterations_Max = MAX( Iterations_Max, k )
        !   Iterations_Ave = Iterations_Ave + k
        !
        !   nIterations = k
        !
        ! END IF

        UVEC = UVEC + DVEC

        Ynew = UVEC(iY)
        Enew = UVEC(iE)
        Jnew(:,1) = UVEC(OS_1+1:OS_1+nE_G)
        Jnew(:,2) = UVEC(OS_2+1:OS_2+nE_G)

        Y = Ynew
        E = Enew

        CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
               ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

        CALL TimersStop( Timer_Im_UpdateFP )

      END IF



    END DO

    J = Jnew

  END SUBROUTINE SolveMatterEquations_Newton


  SUBROUTINE ComputeNeutrinoChemicalPotentials &
    ( D, T, Y, M, dMdT, dMdY, iSpecies )

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: M, dMdT, dMdY
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: Me, dMedT, dMedY
    REAL(DP) :: Mp, dMpdT, dMpdY
    REAL(DP) :: Mn, dMndT, dMndY
    REAL(DP) :: TMP(1), dTMPdT(1), dTMPdY(1)

    ! --- Matter Chemical Potentials and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Me = TMP(1); dMedT = dTMPdT(1); dMedY = dTMPdY(1)

    CALL ComputeProtonChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Mp = TMP(1); dMpdT = dTMPdT(1); dMpdY = dTMPdY(1)

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Mn = TMP(1); dMndT = dTMPdT(1); dMndY = dTMPdY(1)

    ! --- Neutrino Chemical Potential and Derivatives ---

    IF( iSpecies .EQ. iNuE )THEN

      M = ( Me + Mp ) - Mn
      dMdT = ( dMedT + dMpdT ) - dMndT
      dMdY = ( dMedY + dMpdY ) - dMndY

    ELSEIF( iSpecies .EQ. iNuE_Bar )THEN

      M = Mn - ( Me + Mp )
      dMdT = dMndT - ( dMedT + dMpdT )
      dMdY = dMndY - ( dMedY + dMpdY )

    END IF

  END SUBROUTINE ComputeNeutrinoChemicalPotentials


  SUBROUTINE InitializeCollisions_New( iE_B, iE_E )

    INTEGER, INTENT(in) :: iE_B, iE_E

    nE_G = (iE_E-iE_B+1) * nNodesZ(1)

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( CR_N(nE_G,nCR,nSpecies,nDOFX) )
    ALLOCATE( dR_N(nE_G,nCR,nSpecies,nDOFX) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

  END SUBROUTINE InitializeCollisions_New


  SUBROUTINE InitializeCollisions( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    INTEGER :: nZ(4)

    nZ = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( GX_N(nGF,nX_G) )
    ALLOCATE( CF_N(nCF,nX_G) )
    ALLOCATE( dF_N(nCF,nX_G) )
    ALLOCATE( CR_N(nE_G,nCR,nSpecies,nX_G) )
    ALLOCATE( dR_N(nE_G,nCR,nSpecies,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE InitializeCollisions_DGFV( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    INTEGER :: nZ(4)

    nZ = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) )

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( CR_K(nE_G,nCR,nSpecies) )
    ALLOCATE( CR_N(nE_G,nCR,nSpecies,nDOFX) )
    ALLOCATE( dR_N(nE_G,nCR,nSpecies,nDOFX) )

  END SUBROUTINE InitializeCollisions_DGFV


  SUBROUTINE FinalizeCollisions_New

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( CR_N, dR_N )

  END SUBROUTINE FinalizeCollisions_New


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( GX_N )
    DEALLOCATE( CF_N, dF_N )
    DEALLOCATE( CR_N, dR_N )

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE FinalizeCollisions_DGFV

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( CR_K, CR_N, dR_N )

  END SUBROUTINE FinalizeCollisions_DGFV


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

    REAL(DP), INTENT(out) :: &
      E(:), W2(:), W3(:)

    INTEGER  :: iE_G, iE, iN

    ASSOCIATE( dE => MeshE % Width(iE_B0:iE_E0) )

    iE_G = 0
    DO iE = iE_B0, iE_E0
    DO iN = 1, nNodesE

      iE_G = iE_G + 1

      E (iE_G) = NodeCoordinate( MeshE, iE, iN )

      W2(iE_G) = WeightsE(iN) * dE(iE) * E(iE_G)**2
      W3(iE_G) = WeightsE(iN) * dE(iE) * E(iE_G)**3

    END DO
    END DO

    END ASSOCIATE ! -- dE

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE ComputePointsAndWeightsE_DGFV( E, W2, W3 )

    REAL(DP), INTENT(out) :: &
      E(:), W2(:), W3(:)

    INTEGER  :: iE_G, iE, iN
    REAL(DP) :: hc

    hc = PlanckConstant * SpeedOfLight

    ASSOCIATE( dE => MeshE % Width(iE_B0:iE_E0) )

    iE_G = 0
    DO iE = iE_B0, iE_E0
    DO iN = 1, nNodesE

      iE_G = iE_G + 1

      E (iE_G) = NodeCoordinate( MeshE, iE, iN )
      W2(iE_G) = WeightsE(iN) * ( dE(iE) / hc ) * ( E(iE_G) / hc )**2
      W3(iE_G) = W2(iE_G) * ( E(iE_G) / AtomicMassUnit )

    END DO
    END DO

    END ASSOCIATE ! -- dE

  END SUBROUTINE ComputePointsAndWeightsE_DGFV


  SUBROUTINE MapForward_GX( iX_B, iX_E, GX, GX_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iX_B(1):iX_E(1),iX_B(2):iX_E(2),iX_B(3):iX_E(3))
    REAL(DP), INTENT(out) :: &
      GX_N(1:nX_G)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      GX_N(iX) = GX(iN,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_GX


  SUBROUTINE MapForward_F( iX_B, iX_E, FF, FF_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      FF(1:nDOFX,iX_B(1):iX_E(1),iX_B(2):iX_E(2),iX_B(3):iX_E(3))
    REAL(DP), INTENT(out) :: &
      FF_N(1:nX_G)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      FF_N(iX) = FF(iN,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_F


  SUBROUTINE MapBackward_F( iX_B, iX_E, FF, FF_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      FF(1:nDOFX,iX_B(1):iX_E(1),iX_B(2):iX_E(2),iX_B(3):iX_E(3))
    REAL(DP), INTENT(in)  :: &
      FF_N(1:nX_G)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      FF(iN,iX1,iX2,iX3) = FF_N(iX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_F


  SUBROUTINE MapForward_R_New( iE_B, iE_E, RF, RF_K )

    INTEGER,  INTENT(in)  :: &
      iE_B, iE_E
    REAL(DP), INTENT(in)  :: &
      RF(1:nDOF,iE_B:iE_E,1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      RF_K(1:nE_G,1:nCR,1:nSpecies,1:nDOFX)

    INTEGER :: iE, iN_E, iN, iN_X, iCR, iS, iNodeE, iNodeX(3)

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iE  = iE_B, iE_E
    DO iN  = 1, nDOF

      iNodeE = NodeNumberTable(1,  iN)
      iNodeX = NodeNumberTable(2:4,iN)

      iN_E = (iE-1)*nNodesE+iNodeE
      iN_X = NodeNumberTableX3D(iNodeX(1),iNodeX(2),iNodeX(3))

      RF_K(iN_E,iCR,iS,iN_X) = RF(iN,iE,iCR,iS)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_R_New


  SUBROUTINE MapBackward_R_New( iE_B, iE_E, RF, RF_K )

    INTEGER,  INTENT(in)  :: &
      iE_B, iE_E
    REAL(DP), INTENT(out) :: &
      RF(1:nDOF,iE_B:iE_E,1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      RF_K(1:nE_G,1:nCR,1:nSpecies,1:nDOFX)

    INTEGER :: iE, iN_E, iN, iN_X, iCR, iS, iNodeE, iNodeX(3)

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iE  = iE_B, iE_E
    DO iN  = 1, nDOF

      iNodeE = NodeNumberTable(1,  iN)
      iNodeX = NodeNumberTable(2:4,iN)

      iN_E = (iE-1)*nNodesE+iNodeE
      iN_X = NodeNumberTableX3D(iNodeX(1),iNodeX(2),iNodeX(3))

      RF(iN,iE,iCR,iS) = RF_K(iN_E,iCR,iS,iN_X)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_R_New


  SUBROUTINE MapForward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      RF(1:nDOF,iZ_B(1):iZ_E(1),iZ_B(2):iZ_E(2),iZ_B(3):iZ_E(3),iZ_B(4):iZ_E(4))
    REAL(DP), INTENT(out) :: &
      RF_N(1:nE_G,1:nX_G)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    iX = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = iX + 1

      iE = 0
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF_N(iE,iX) = RF(iN,iZ1,iZ2,iZ3,iZ4)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_R


  SUBROUTINE MapBackward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(out) :: &
      RF(1:nDOF,iZ_B(1):iZ_E(1),iZ_B(2):iZ_E(2),iZ_B(3):iZ_E(3),iZ_B(4):iZ_E(4))
    REAL(DP), INTENT(in)  :: &
      RF_N(1:nE_G,1:nX_G)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    iX = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = iX + 1

      iE = 0
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF(iN,iZ1,iZ2,iZ3,iZ4) = RF_N(iE,iX)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_R


  SUBROUTINE MapForward_R_DGFV( iZ1_B, iZ1_E, RF, RF_N, RF_K )

    INTEGER,  INTENT(in)  :: &
      iZ1_B, iZ1_E
    REAL(DP), INTENT(in)  :: &
      RF(1:nDOF,iZ1_B:iZ1_E)
    REAL(DP), INTENT(out) :: &
      RF_N(1:nE_G,1:nDOFX), RF_K(1:nE_G)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    RF_K = Zero
    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = NodeNumberTableX3D(iN2,iN3,iN4)

      iE = 0
      DO iZ1 = iZ1_B, iZ1_E
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF_N(iE,iX) = RF(iN,iZ1)
        RF_K(iE)    = RF_K(iE) + WeightsX_q(iX) * RF_N(iE,iX)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_R_DGFV


  SUBROUTINE MapBackward_R_DGFV( iZ1_B, iZ1_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ1_B, iZ1_E
    REAL(DP), INTENT(out) :: &
      RF(1:nDOF,iZ1_B:iZ1_E)
    REAL(DP), INTENT(in)  :: &
      RF_N(1:nE_G,1:nDOFX)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = NodeNumberTableX3D(iN2,iN3,iN4)

      iE = 0
      DO iZ1 = iZ1_B, iZ1_E
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF(iN,iZ1) = RF_N(iE,iX)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_R_DGFV


  ELEMENTAL SUBROUTINE ComputePrimitive_Euler &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), INTENT(out) :: D, V_1, V_2, V_3, E, De
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    D   = N
    V_1 = S_1 / ( Gm_dd_11 * N )
    V_2 = S_2 / ( Gm_dd_22 * N )
    V_3 = S_3 / ( Gm_dd_33 * N )
    E   = G - Half * ( S_1**2 / Gm_dd_11 &
                       + S_2**2 / Gm_dd_22 &
                       + S_3**2 / Gm_dd_33 ) / N
    De  = Ne

  END SUBROUTINE ComputePrimitive_Euler


  ELEMENTAL SUBROUTINE ComputeConserved_Euler &
    ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: D, V_1, V_2, V_3, E, De
    REAL(DP), INTENT(out) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    N   = D
    S_1 = D * Gm_dd_11 * V_1
    S_2 = D * Gm_dd_22 * V_2
    S_3 = D * Gm_dd_33 * V_3
    G   = E + Half * D * (   Gm_dd_11 * V_1**2 &
                           + Gm_dd_22 * V_2**2 &
                           + Gm_dd_33 * V_3**2 )
    Ne  = De

  END SUBROUTINE ComputeConserved_Euler


  PURE REAL(DP) FUNCTION ENORM( X )

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM

  PURE REAL(DP) FUNCTION WNORM( X, W )

    REAL(DP), DIMENSION(:), INTENT(in) :: X, W

    WNORM = SQRT( DOT_PRODUCT( W * X, X ) )

    RETURN
  END FUNCTION WNORM

  SUBROUTINE InitializeNonlinearSolverTally

    USE ProgramHeaderModule, ONLY: nX

    TallyNonlinearSolver = .TRUE.
    TallyFileNumber = 0

    ALLOCATE( MinIterations_K(nX(1),nX(2),nX(3)) )
    ALLOCATE( MaxIterations_K(nX(1),nX(2),nX(3)) )
    ALLOCATE( AveIterations_K(nX(1),nX(2),nX(3)) )
    ALLOCATE( AveIterationsInner_K(nX(1),nX(2),nX(3)) )

    MinIterations_K = 0
    MaxIterations_K = 0
    AveIterations_K = 0
    AveIterationsInner_K = 0

  END SUBROUTINE InitializeNonlinearSolverTally


  SUBROUTINE FinalizeNonlinearSolverTally

    DEALLOCATE( MinIterations_K )
    DEALLOCATE( MaxIterations_K )
    DEALLOCATE( AveIterations_K )
    DEALLOCATE( AveIterationsInner_K )

  END SUBROUTINE FinalizeNonlinearSolverTally


  SUBROUTINE WriteNonlinearSolverTally( t )

    USE HDF5
    USE UnitsModule, ONLY: &
      Millisecond, &
      Kilometer
    USE ProgramHeaderModule, ONLY: &
      nX
    USE MeshModule, ONLY: &
      MeshX

    REAL(DP), INTENT(in) :: t

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    CHARACTER(256) :: GroupName
    INTEGER        :: HDFERR
    INTEGER(HID_T) :: FILE_ID

    IF( .NOT. TallyNonlinearSolver ) RETURN

    WRITE( FileNumberString, FMT='(i6.6)') TallyFileNumber

    FileName = '../Output/' // 'NonlinearSolverTally_' // &
               FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1D_REAL &
           ( [ t ] / Millisecond, DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1D_REAL &
           ( MeshX(1) % Center(1:nX(1)) / Kilometer, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1D_REAL &
           ( MeshX(2) % Center(1:nX(2)) / Kilometer, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1D_REAL &
           ( MeshX(3) % Center(1:nX(3)) / Kilometer, DatasetName, FILE_ID )

    ! --- Write Iteration Counts ---

    GroupName = 'Iteration Counts'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Max Iterations'

    CALL WriteDataset3D_INTEGER &
           ( MaxIterations_K, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Min Iterations'

    CALL WriteDataset3D_INTEGER &
           ( MinIterations_K, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Average Iterations'

    CALL WriteDataset3D_INTEGER &
           ( AveIterations_K, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Average Iterations (Inner)'

    CALL WriteDataset3D_INTEGER &
           ( AveIterationsInner_K, DatasetName, FILE_ID )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

    TallyFileNumber = TallyFileNumber + 1

  END SUBROUTINE WriteNonlinearSolverTally


  SUBROUTINE CreateGroupHDF( FileName, GroupName, FILE_ID )

    USE HDF5

    CHARACTER(len=*), INTENT(in) :: FileName
    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER        :: HDFERR
    INTEGER(HID_T) :: GROUP_ID

    CALL H5GCREATE_F( FILE_ID, TRIM( GroupName ), GROUP_ID, HDFERR )

    CALL H5GCLOSE_F( GROUP_ID, HDFERR )

  END SUBROUTINE CreateGroupHDF


  SUBROUTINE WriteDataset1D_REAL( Dataset, DatasetName, FILE_ID )

    USE HDF5

    REAL(DP),         INTENT(in) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER          :: HDFERR
    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SIZE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 1, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset1D_REAL


  SUBROUTINE WriteDataset3D_INTEGER( Dataset, DatasetName, FILE_ID )

    USE HDF5

    INTEGER,          INTENT(in) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER          :: HDFERR
    INTEGER(HSIZE_T) :: DATASIZE(3)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 3, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_INTEGER, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_INTEGER, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset3D_INTEGER


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
