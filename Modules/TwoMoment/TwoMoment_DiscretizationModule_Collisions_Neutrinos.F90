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
    Timer_Im_ComputeTS_Aux, &
    Timer_Im_ComputeOpacity, &
    Timer_Im_MapForward, &
    Timer_Im_Solve, &
    Timer_Im_Out, &
    Timer_Im_ComputeTS_Prim, &
    Timer_Im_Increment, &
    Timer_Im_MapBackward
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
    ComputeNeutrinoOpacities_EC_Point, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Point, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacities_Pair_Point, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit
  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_New
  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_DGFV

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV
  REAL(DP), PARAMETER :: Unit_E = Erg / Gram

  LOGICAL, PARAMETER :: ReportConvergenceData = .FALSE.
  INTEGER  :: Iterations_Min
  INTEGER  :: Iterations_Max
  INTEGER  :: Iterations_Ave

  LOGICAL, PARAMETER :: SolveMatter = .TRUE.

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
      dU_F(1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies)

    INTEGER  :: iX1, iX2, iX3, iGF, iCF, iCR, iS, iNodeX, iE
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

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

!!$    PRINT*, "ComputeIncrement_TwoMoment_Implicit_New"

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

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

!!$      PRINT*, "iX1,iX2,iX3 = ", iX1,iX2,iX3

      CALL TimersStart( Timer_Im_In )

      DO iCF = 1, nCF

        CF_N(:,            iCF) = U_F (:,iX1,iX2,iX3,iCF)
        dU_F(:,iX1,iX2,iX3,iCF) = CF_N(:,            iCF)

      END DO

      CALL ComputePrimitive_Euler &
             ( CF_N(:,iCF_D ), CF_N(:,iCF_S1), CF_N(:,iCF_S2), &
               CF_N(:,iCF_S3), CF_N(:,iCF_E ), CF_N(:,iCF_Ne), &
               PF_N(:,iPF_D ), PF_N(:,iPF_V1), PF_N(:,iPF_V2), &
               PF_N(:,iPF_V3), PF_N(:,iPF_E ), PF_N(:,iPF_Ne), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

!!$      WRITE(*,'(A4,A32,ES10.4E2)') '', 'Copy to CF_N: ', Timer_Implicit_In

!!$      PRINT*, "D  = ", PF_N(:,iPF_D  )
!!$      PRINT*, "V1 = ", PF_N(:,iPF_V1 )
!!$      PRINT*, "V2 = ", PF_N(:,iPF_V2 )
!!$      PRINT*, "V3 = ", PF_N(:,iPF_V3 )
!!$      PRINT*, "E  = ", PF_N(:,iPF_E  )
!!$      PRINT*, "Ne = ", PF_N(:,iPF_Ne )

      CALL TimersStart( Timer_Im_ComputeTS_Aux )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), &
               AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye) )

!!$      PRINT*, "T  = ", AF_N(:,iAF_T)
!!$      PRINT*, "E  = ", AF_N(:,iAF_E)
!!$      PRINT*, "Ye = ", AF_N(:,iAF_Ye)

      CALL TimersStop( Timer_Im_ComputeTS_Aux )

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

      END DO

      Chi_NES  = Zero
      Eta_NES  = Zero
      Chi_Pair = Zero
      Eta_Pair = Zero

      CALL TimersStop( Timer_Im_ComputeOpacity )

!!$      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputeNeutrinoOpacities: ', Timer_Im_ComputeOpacity

      CALL TimersStart( Timer_Im_MapForward )

      CALL MapForward_R_New &
             ( iE_B0, iE_E0, U_R(:,iE_B0:iE_E0,iX1,iX2,iX3,:,:), CR_N )

      CALL TimersStop( Timer_Im_MapForward )

      CALL TimersStop( Timer_Im_In )

!!$      WRITE(*,'(A4,A32,ES10.4E2)') '', 'MapForward_R: ', Timer_Im_MapR

      CALL TimersStart( Timer_Im_Solve )

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

!!$        DO iNodeX = 1, nDOFX
!!$
!!$          CALL SolveMatterEquations_EmAb &
!!$                 ( CR_N    (:,iCR_N,iNuE,    iNodeX), &
!!$                   CR_N    (:,iCR_N,iNuE_Bar,iNodeX), &
!!$                   dt * Chi(:,      iNuE,    iNodeX), &
!!$                   dt * Chi(:,      iNuE_Bar,iNodeX), &
!!$                   fEQ     (:,      iNuE,    iNodeX), &
!!$                   fEQ     (:,      iNuE_Bar,iNodeX), &
!!$                   PF_N(iNodeX,iPF_D ), AF_N(iNodeX,iAF_T), &
!!$                   AF_N(iNodeX,iAF_Ye), AF_N(iNodeX,iAF_E) )
!!$
!!$        END DO

        DO iNodeX = 1, nDOFX

          CALL SolveMatterEquations_FP_Coupled &
                 ( dt, iNuE, iNuE_Bar, &
                   CR_N(:,iCR_N,1:2,iNodeX), &
                   Chi (:,      1:2,iNodeX), &
                   fEQ (:,      1:2,iNodeX), &
                   Chi_NES (:,  1:2,iNodeX), &
                   Eta_NES (:,  1:2,iNodeX), &
                   Chi_Pair(:,  1:2,iNodeX), &
                   Eta_Pair(:,  1:2,iNodeX), &
                   PF_N(iNodeX,iPF_D ), &
                   AF_N(iNodeX,iAF_T ), &
                   AF_N(iNodeX,iAF_Ye), &
                   AF_N(iNodeX,iAF_E ) )

        END DO

      END IF

      CALL TimersStop( Timer_Im_Solve )

!!$      WRITE(*,'(A4,A32,ES10.4E2)') '', 'SolveMatterEquations_EmAb: ', Timer_Im_Solve

      CALL TimersStart( Timer_Im_Out )

      CALL TimersStart( Timer_Im_ComputeTS_Prim )

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), &
               PF_N(:,iPF_E), AF_N(:,iAF_E), PF_N(:,iPF_Ne) )

      CALL TimersStop( Timer_Im_ComputeTS_Prim )

      CALL ComputeConserved_Euler &
             ( PF_N(:,iPF_D ), PF_N(:,iPF_V1), PF_N(:,iPF_V2), &
               PF_N(:,iPF_V3), PF_N(:,iPF_E ), PF_N(:,iPF_Ne), &
               CF_N(:,iCF_D ), CF_N(:,iCF_S1), CF_N(:,iCF_S2), &
               CF_N(:,iCF_S3), CF_N(:,iCF_E ), CF_N(:,iCF_Ne), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               GX(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL TimersStart( Timer_Im_Increment )

      ! --- Conserved Fluid Increment ---

      DO iCF = 1, nCF

        dU_F(:,iX1,iX2,iX3,iCF) &
          = ( CF_N(:,iCF) - dU_F(:,iX1,iX2,iX3,iCF) ) / dt

      END DO

      ! --- Update Radiation Fields ---

      DO iNodeX = 1, nDOFX
      DO iS     = 1, nSpecies

        Chi_T = Chi(:,iS,iNodeX) &
                  + Chi_NES(:,iS,iNodeX) + Chi_Pair(:,iS,iNodeX)

        Kappa = Chi_T + Sig(:,iS,iNodeX)

        Eta_T = Chi(:,iS,iNodeX) * fEQ(:,iS,iNodeX) &
                  + Eta_NES(:,iS,iNodeX) + Eta_Pair(:,iS,iNodeX)

        ! --- Number Density ---

        CR_N(:,iCR_N,iS,iNodeX) &
          = ( CR_N(:,iCR_N,iS,iNodeX) + dt * Eta_T ) / ( One + dt * Chi_T )

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

        dR_N(:,iCR_N,iS,iNodeX) &
          = Eta_T - Chi_T * CR_N(:,iCR_N,iS,iNodeX)

        dR_N(:,iCR_G1,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G1,iS,iNodeX)

        dR_N(:,iCR_G2,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G2,iS,iNodeX)

        dR_N(:,iCR_G3,iS,iNodeX) &
          = - Kappa * CR_N(:,iCR_G3,iS,iNodeX)
        
      END DO
      END DO

      CALL TimersStop( Timer_Im_Increment )

      CALL TimersStart( Timer_Im_MapBackward )

      CALL MapBackward_R_New &
             ( iE_B0, iE_E0, dU_R(:,iE_B0:iE_E0,iX1,iX2,iX3,:,:), dR_N )

      CALL TimersStop( Timer_Im_MapBackward )

      CALL TimersStop( Timer_Im_Out )

    END DO
    END DO
    END DO

    DEALLOCATE &
      ( Kappa, Chi_T, Eta_T, Chi, fEQ, Sig, Chi_NES, Eta_NES, &
        Chi_Pair, Eta_Pair )

    CALL FinalizeCollisions_New

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
      dU_F(1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies)

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
      dU_F(1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies)

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

    REAL(DP), INTENT(in)    :: J  (1:nE_G)
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
      W3_S = W3_N / PlanckConstant**3 / AtomicMassUnit

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
    C(iE) = DOT_PRODUCT( Theta3_N(:), J(:) ) + N_B * U(iE)

    ! --- Electron Fraction Equation ---

    FVEC(iY) = DOT_PRODUCT( Theta2_N(:), J0(:) ) + N_B * U(iY) - C(iY)

    ! --- Internal Energy Equation ---

    FVEC(iE) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + N_B * U(iE) - C(iE)

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

      FJAC(2,2) = DOT_PRODUCT( Theta3_N(:), dJ0dE_Y(:) ) + N_B

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

      FVEC(2) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + N_B * U(2) - C(2)

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


  SUBROUTINE SolveMatterEquations_FP_Coupled &
    ( dt, iS_1, iS_2, J, Chi, J0, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      D, T, Y, E )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(in)    :: J       (1:nE_G,1:2)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: J0      (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_NES (1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Chi_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: Eta_Pair(1:nE_G,1:2)
    REAL(DP), INTENT(inout) :: D, T, Y, E

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: M = 5
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: LWORK = 2 * M
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

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
    REAL(DP) :: Phi_0_In_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_NES (1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_In_Pair(1:nE_G,1:nE_G,1:2)
    REAL(DP) :: Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:2)

!!$    PRINT*
!!$    PRINT*, "SolveMatterEquations_FP_Coupled Start"
!!$    PRINT*

    OS_1 = 2
    OS_2 = 2 + nE_G

    h3  = PlanckConstant**3
    N_B = D / AtomicMassUnit

    Yold = Y
    Eold = E

    Jold(:,1) = J(:,1)
    Jold(:,2) = J(:,2)

    S_Y = N_B * Yold
    S_E = D   * Eold

    W2_S = FourPi * W2_N / h3
    W3_S = FourPi * W3_N / h3

    Unew = One ! --- Initial Guess

    C(iY) = DOT_PRODUCT( W2_S, Jold(:,1) - Jold(:,2) ) / S_Y
    C(iE) = DOT_PRODUCT( W3_S, Jold(:,1) + Jold(:,2) ) / S_E

    CALL ComputeEquilibriumDistributions_Point &
           ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

    CALL ComputeEquilibriumDistributions_Point &
           ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

    Eta(:,1) = Chi(:,1) * J0(:,1)
    Eta(:,2) = Chi(:,2) * J0(:,2)

    ! --- NES Kernels ---

    CALL ComputeNeutrinoOpacities_NES_Point &
           ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_NES(:,:,1), Phi_0_Ot_NES(:,:,1) )

    CALL ComputeNeutrinoOpacities_NES_Point &
           ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_NES(:,:,2), Phi_0_Ot_NES(:,:,2) )

    ! --- NES Emissivities and Opacities ---

    ! --- Neutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,1), nE_G, &
                W2_N * Jold(:,1),       1, Zero, Eta_NES(:,1), 1 )

    Chi_NES(:,1) = Eta_NES(:,1)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,1), nE_G, &
                W2_N * (One-Jold(:,1)), 1, One,  Chi_NES(:,1), 1 )

    ! --- Antineutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_NES(:,:,2), nE_G, &
                W2_N * Jold(:,2),       1, Zero, Eta_NES(:,2), 1 )

    Chi_NES(:,2) = Eta_NES(:,2)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_NES(:,:,2), nE_G, &
                W2_N * (One-Jold(:,2)), 1, One,  Chi_NES(:,2), 1 )

    ! --- Pair Kernels ---

    CALL ComputeNeutrinoOpacities_Pair_Point &
           ( 1, nE_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_Pair(:,:,1), Phi_0_Ot_Pair(:,:,1) )

    CALL ComputeNeutrinoOpacities_Pair_Point &
           ( 1, nE_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_Pair(:,:,2), Phi_0_Ot_Pair(:,:,2) )

    ! --- Pair Emissivities and Opacities ---

    ! --- Neutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,1), nE_G, &
                W2_N * (One-Jold(:,2)), 1, Zero, Eta_Pair(:,1), 1 )

    Chi_Pair(:,1) = Eta_Pair(:,1)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,1), nE_G, &
                W2_N * Jold(:,2),       1, One,  Chi_Pair(:,1), 1 )

    ! --- Antineutrino ---

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_In_Pair(:,:,2), nE_G, &
                W2_N * (One-Jold(:,1)), 1, Zero, Eta_Pair(:,2), 1 )

    Chi_Pair(:,2) = Eta_Pair(:,2)

    CALL DGEMV( 'T', nE_G, nE_G, One, Phi_0_Ot_Pair(:,:,2), nE_G, &
                W2_N * Jold(:,1),       1, One,  Chi_Pair(:,2), 1 )

    ! --- Update Neutrino Densities ---

    Jnew(:,1) &
      = ( Jold(:,1) + dt * ( Eta(:,1) + Eta_NES(:,1) + Eta_Pair(:,1) ) ) &
            / ( One + dt * ( Chi(:,1) + Chi_NES(:,1) + Chi_Pair(:,1) ) )

    Jnew(:,2) &
      = ( Jold(:,2) + dt * ( Eta(:,2) + Eta_NES(:,2) + Eta_Pair(:,2) ) ) &
            / ( One + dt * ( Chi(:,2) + Chi_NES(:,2) + Chi_Pair(:,2) ) )

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k  = k + 1
      mk = MIN( M, k )

!!$      PRINT*, "  k, mk = ", k, mk

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

      IF( .NOT. CONVERGED )THEN

        ! --- Recompute Equilibrium Distributions and Emissivities ---

        CALL ComputeEquilibriumDistributions_Point &
               ( 1, nE_G, E_N, D, T, Y, J0(:,1), iS_1 )

        CALL ComputeEquilibriumDistributions_Point &
               ( 1, nE_G, E_N, D, T, Y, J0(:,2), iS_2 )

        Eta(:,1) = Chi(:,1) * J0(:,1)
        Eta(:,2) = Chi(:,2) * J0(:,2)

        ! --- Recompute Kernels ---

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

      END IF

    END DO

!!$    PRINT*
!!$    PRINT*, "  D = ", D / Unit_D
!!$    PRINT*, "  T = ", T / Kelvin
!!$    PRINT*, "  E = ", E / Unit_E
!!$    PRINT*, "  Y = ", Y
!!$    PRINT*
!!$    PRINT*, "SolveMatterEquations_FP_Coupled Done"
!!$    PRINT*

  END SUBROUTINE SolveMatterEquations_FP_Coupled


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


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
