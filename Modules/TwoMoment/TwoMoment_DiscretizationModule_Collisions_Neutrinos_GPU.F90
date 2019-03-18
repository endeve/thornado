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
    MeV
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
    ComputeNeutrinoOpacities_EC_Point, &
    ComputeNeutrinoOpacities_ES_Point, &
    ComputeNeutrinoOpacities_EC_Points, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_New

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

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
  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)

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
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies)

    INTEGER  :: iX1, iX2, iX3, iGF, iCF, iCR, iS, iE, iN_E
    INTEGER  :: iNode, iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: CF_N(1:nDOFX,1:nCF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: PF_N(1:nDOFX,1:nPF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: AF_N(1:nDOFX,1:nAF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: GX_N(1:nDOFX,1:nGF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: Kappa
    REAL(DP) :: Chi (1:nNodesZ(1)*(iZ_E0(1)-iZ_B0(1)+1),1:nSpecies, &
                     1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: Sig (1:nNodesZ(1)*(iZ_E0(1)-iZ_B0(1)+1),1:nSpecies, &
                     1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: fEQ (1:nNodesZ(1)*(iZ_E0(1)-iZ_B0(1)+1),1:nSpecies, &
                     1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: CR_N(1:nNodesZ(1)*(iZ_E0(1)-iZ_B0(1)+1),1:nCR,1:nSpecies, &
                     1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: dR_N(1:nNodesZ(1)*(iZ_E0(1)-iZ_B0(1)+1),1:nCR,1:nSpecies, &
                     1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    CALL TimersStart( Timer_Implicit )

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    CALL InitializeCollisions_New( iE_B0, iE_E0 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GX, U_F, U_R, iX_B0, iX_E0, iX_B1, iX_E1 ) &
    !$OMP MAP( alloc: dU_F, dU_R, CF_N, PF_N, AF_N, GX_N, CR_N, dR_N, Chi, Sig, fEQ )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U_F, U_R ) &
    !$ACC CREATE( dU_F, dU_R, CF_N, PF_N, AF_N, GX_N, CR_N, dR_N, Chi, Sig, fEQ )
#endif

    CALL TimersStart( Timer_Im_In )

    ! --- Copy inputs to locals ---

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX
              CF_N(iNodeX,iCF,iX1,iX2,iX3) = U_F(iNodeX,iX1,iX2,iX3,iCF)
              dU_F(iNodeX,iX1,iX2,iX3,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)
            END DO
          END DO
        END DO
      END DO
    END DO

    DO iGF = 1, nGF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX
              GX_N(iNodeX,iGF,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    ! --- Compute Primitive Quantities ---

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iNodeX = 1, nDOFX

            CALL ComputePrimitive_Euler &
                   ( CF_N(iNodeX,iCF_D ,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S1,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S2,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S3,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_E ,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_Ne,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_D ,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V1,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V2,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V3,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_E ,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_Ne,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3) )

          END DO
        END DO
      END DO
    END DO

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeTS_Aux )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputeThermodynamicStates_Auxiliary_TABLE &
                 ( PF_N(:,iPF_D ,iX1,iX2,iX3), &
                   PF_N(:,iPF_E ,iX1,iX2,iX3), &
                   PF_N(:,iPF_Ne,iX1,iX2,iX3), &
                   AF_N(:,iAF_T ,iX1,iX2,iX3), &
                   AF_N(:,iAF_E ,iX1,iX2,iX3), &
                   AF_N(:,iAF_Ye,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_ComputeTS_Aux )

    ! --- Opacity Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iS = 1, nSpecies

            CALL ComputeNeutrinoOpacities_EC_Points &
                   ( 1, nE_G, 1, nDOFX, E_N, &
                     PF_N(:,iPF_D ,iX1,iX2,iX3), &
                     AF_N(:,iAF_T ,iX1,iX2,iX3), &
                     AF_N(:,iAF_Ye,iX1,iX2,iX3), &
                     iS, Chi(:,iS,:,iX1,iX2,iX3) )

          END DO
        END DO
      END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iNodeX = 1, nDOFX
            DO iS = 1, nSpecies
              DO iN_E = 1, nE_G
                Sig(iN_E,iS,iNodeX,iX1,iX2,iX3) = Zero
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Rearrange data to group energy nodes together ---

    CALL TimersStart( Timer_Im_MapForward )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
            DO iX1 = iX_B0(1), iX_E0(1)
              DO iE = iE_B0, iE_E0
                DO iNodeX3 = 1, nNodesX(3)
                  DO iNodeX2 = 1, nNodesX(2)
                    DO iNodeX1 = 1, nNodesX(1)
                      DO iNodeE = 1, nNodesE

                        iNode  = + ( iNodeX3 - 1 ) * nNodesE * nNodesX(1) * nNodesX(2) &
                                 + ( iNodeX2 - 1 ) * nNodesE * nNodesX(1) &
                                 + ( iNodeX1 - 1 ) * nNodesE &
                                 +   iNodeE

                        iNodeX = + ( iNodeX3 - 1 ) * nNodesX(1) * nNodesX(2) &
                                 + ( iNodeX2 - 1 ) * nNodesX(1) &
                                 +   iNodeX1

                        iN_E   = + ( iE - iE_B0 ) * nNodesE &
                                 +   iNodeE

                        CR_N(iN_E,iCR,iS,iNodeX,iX1,iX2,iX3) = U_R(iNode,iE,iX1,iX2,iX3,iCR,iS)

                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_MapForward )

    CALL TimersStop( Timer_Im_In )

    CALL TimersStart( Timer_Im_Solve )

!    IF( nSpecies .EQ. 1 )THEN
!
!      ! --- Single Species (Electron Neutrinos) ---
!
!      DO iX3 = iX_B0(3), iX_E0(3)
!        DO iX2 = iX_B0(2), iX_E0(2)
!          DO iX1 = iX_B0(1), iX_E0(1)
!            DO iNodeX = 1, nDOFX
!
!              CALL SolveMatterEquations_EmAb_NuE &
!                     ( CR_N    (:,iCR_N,iNuE,iNodeX), &
!                       dt * Chi(:,      iNuE,iNodeX), &
!                       fEQ     (:,      iNuE,iNodeX), &
!                       PF_N(iNodeX,iPF_D ), AF_N(iNodeX,iAF_T), &
!                       AF_N(iNodeX,iAF_Ye), AF_N(iNodeX,iAF_E) )
!
!            END DO
!          END DO
!        END DO
!      END DO
!
!    ELSE

      ! --- Electron Neutrinos and Antineutrinos ---

      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX

              CALL SolveMatterEquations_EmAb &
                     ( CR_N    (:,iCR_N,iNuE,    iNodeX,iX1,iX2,iX3), &
                       CR_N    (:,iCR_N,iNuE_Bar,iNodeX,iX1,iX2,iX3), &
                       dt * Chi(:,      iNuE,    iNodeX,iX1,iX2,iX3), &
                       dt * Chi(:,      iNuE_Bar,iNodeX,iX1,iX2,iX3), &
                       fEQ     (:,      iNuE,    iNodeX,iX1,iX2,iX3), &
                       fEQ     (:,      iNuE_Bar,iNodeX,iX1,iX2,iX3), &
                       PF_N(iNodeX,iPF_D ,iX1,iX2,iX3), &
                       AF_N(iNodeX,iAF_T ,iX1,iX2,iX3), &
                       AF_N(iNodeX,iAF_Ye,iX1,iX2,iX3), &
                       AF_N(iNodeX,iAF_E ,iX1,iX2,iX3) )

            END DO
          END DO
        END DO
      END DO

!    END IF

    CALL TimersStop( Timer_Im_Solve )

    CALL TimersStart( Timer_Im_Out )

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeTS_Prim )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputeThermodynamicStates_Primitive_TABLE &
                 ( PF_N(:,iPF_D ,iX1,iX2,iX3), &
                   AF_N(:,iAF_T ,iX1,iX2,iX3), &
                   AF_N(:,iAF_Ye,iX1,iX2,iX3), &
                   PF_N(:,iPF_E ,iX1,iX2,iX3), &
                   AF_N(:,iAF_E ,iX1,iX2,iX3), &
                   PF_N(:,iPF_Ne,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_ComputeTS_Prim )

    ! --- Compute Conserved Quantities ---

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iNodeX = 1, nDOFX

            CALL ComputeConserved_Euler &
                   ( PF_N(iNodeX,iPF_D ,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V1,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V2,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_V3,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_E ,iX1,iX2,iX3), &
                     PF_N(iNodeX,iPF_Ne,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_D ,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S1,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S2,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_S3,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_E ,iX1,iX2,iX3), &
                     CF_N(iNodeX,iCF_Ne,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3), &
                     GX_N(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3) )

          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Im_Increment )

    ! --- Conserved Fluid Increment ---

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX

              dU_F(iNodeX,iX1,iX2,iX3,iCF) &
                = ( CF_N(iNodeX,iCF,iX1,iX2,iX3) - dU_F(iNodeX,iX1,iX2,iX3,iCF) ) / dt

            END DO
          END DO
        END DO
      END DO
    END DO

    ! --- Update Radiation Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iNodeX = 1, nDOFX
            DO iS = 1, nSpecies
              DO iN_E = 1, nE_G

                Kappa = Chi(iN_E,iS,iNodeX,iX1,iX2,iX3) + Sig(iN_E,iS,iNodeX,iX1,iX2,iX3)

                ! --- Number Density ---

                CR_N(iN_E,iCR_N,iS,iNodeX,iX1,iX2,iX3) &
                  = ( dt * Chi(iN_E,iS,iNodeX,iX1,iX2,iX3) * fEQ(iN_E,iS,iNodeX,iX1,iX2,iX3) &
                      + CR_N(iN_E,iCR_N,iS,iNodeX,iX1,iX2,iX3) ) / ( One + dt * Chi(iN_E,iS,iNodeX,iX1,iX2,iX3) )

                ! --- Number Flux (1) ---

                CR_N(iN_E,iCR_G1,iS,iNodeX,iX1,iX2,iX3) &
                  = CR_N(iN_E,iCR_G1,iS,iNodeX,iX1,iX2,iX3) / ( One + dt * Kappa )

                ! --- Number Flux (2) ---

                CR_N(iN_E,iCR_G2,iS,iNodeX,iX1,iX2,iX3) &
                  = CR_N(iN_E,iCR_G2,iS,iNodeX,iX1,iX2,iX3) / ( One + dt * Kappa )

                ! --- Number Flux (3) ---

                CR_N(iN_E,iCR_G3,iS,iNodeX,iX1,iX2,iX3) &
                  = CR_N(iN_E,iCR_G3,iS,iNodeX,iX1,iX2,iX3) / ( One + dt * Kappa )

                ! --- Increments ---

                dR_N(iN_E,iCR_N,iS,iNodeX,iX1,iX2,iX3) &
                  = Chi(iN_E,iS,iNodeX,iX1,iX2,iX3) * ( fEQ(iN_E,iS,iNodeX,iX1,iX2,iX3) - CR_N(iN_E,iCR_N,iS,iNodeX,iX1,iX2,iX3) )

                dR_N(iN_E,iCR_G1,iS,iNodeX,iX1,iX2,iX3) &
                  = - Kappa * CR_N(iN_E,iCR_G1,iS,iNodeX,iX1,iX2,iX3)

                dR_N(iN_E,iCR_G2,iS,iNodeX,iX1,iX2,iX3) &
                  = - Kappa * CR_N(iN_E,iCR_G2,iS,iNodeX,iX1,iX2,iX3)

                dR_N(iN_E,iCR_G3,iS,iNodeX,iX1,iX2,iX3) &
                  = - Kappa * CR_N(iN_E,iCR_G3,iS,iNodeX,iX1,iX2,iX3)
            
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_Increment )

    ! --- Rearrange data back to original layout ---

    CALL TimersStart( Timer_Im_MapBackward )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
            DO iX1 = iX_B0(1), iX_E0(1)
              DO iE = iE_B0, iE_E0
                DO iNodeX3 = 1, nNodesX(3)
                  DO iNodeX2 = 1, nNodesX(2)
                    DO iNodeX1 = 1, nNodesX(1)
                      DO iNodeE = 1, nNodesE

                        iNode  = + ( iNodeX3 - 1 ) * nNodesE * nNodesX(1) * nNodesX(2) &
                                 + ( iNodeX2 - 1 ) * nNodesE * nNodesX(1) &
                                 + ( iNodeX1 - 1 ) * nNodesE &
                                 +   iNodeE

                        iNodeX = + ( iNodeX3 - 1 ) * nNodesX(1) * nNodesX(2) &
                                 + ( iNodeX2 - 1 ) * nNodesX(1) &
                                 +   iNodeX1

                        iN_E   = + ( iE - iE_B0 ) * nNodesE &
                                 +   iNodeE

                        dU_R(iNode,iE,iX1,iX2,iX3,iCR,iS) = dR_N(iN_E,iCR,iS,iNodeX,iX1,iX2,iX3)

                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_MapBackward )

    CALL TimersStop( Timer_Im_Out )

    CALL FinalizeCollisions_New

    CALL TimersStop( Timer_Implicit )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: GX, U_F, U_R, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               CF_N, PF_N, AF_N, CR_N, dR_N, Kappa, Chi, Sig, fEQ ) &
    !$OMP MAP( from: dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( GX, U_F, U_R, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$ACC         CF_N, PF_N, AF_N, GX_N, CR_N, dR_N, Kappa, Chi, Sig, fEQ ) &
    !$ACC COPYOUT( dU_F, dU_R )
#endif

#ifdef THORNADO_DEBUG_IMPLICIT
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_F)', MAXLOC(dU_F)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_F)', MAXVAL(dU_F)
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_R)', MAXLOC(dU_R)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_R)', MAXVAL(dU_R)
#endif

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_New


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

      Theta2_N = FourPi * W2_N * Chi / ( One + Chi )
      Theta3_N = FourPi * W3_N * Chi / ( One + Chi )

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: E_N, W2_N, W3_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N )
#endif

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

  END SUBROUTINE InitializeCollisions_New


  SUBROUTINE FinalizeCollisions_New

    DEALLOCATE( E_N, W2_N, W3_N )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC DELETE( E_N, W2_N, W3_N )
#endif

  END SUBROUTINE FinalizeCollisions_New


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

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

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE MapForward_R_New( iE_B, iE_E, RF, RF_K )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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


  ELEMENTAL SUBROUTINE ComputePrimitive_Euler &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
