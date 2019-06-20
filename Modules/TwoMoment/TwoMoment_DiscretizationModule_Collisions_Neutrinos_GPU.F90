#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE, INTRINSIC :: ieee_arithmetic, ONLY: &
    IEEE_IS_NAN

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
    Timer_Im_MapBackward, &
    Timer_FP_FD, &
    Timer_FP_Op, &
    Timer_FP_Rate, &
    Timer_FP_RHS, &
    Timer_FP_RateK, &
    Timer_FP_RHSK, &
    Timer_FP_AA, &
    Timer_FP_LS, &
    Timer_FP_Update, &
    Timer_FP_FDK, &
    Timer_FP_OpK, &
    Timer_FP_Shift
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply, &
    LinearLeastSquares, &
    LinearLeastSquares_LWORK
    !VectorNorm2
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
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
    ComputeTemperatureFromSpecificInternalEnergyPoint_TABLE, &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_Point, &
    ComputeEquilibriumDistributions_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacitiesRates_NES_Points, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacitiesRates_Pair_Points

  USE wlEOSInversionModule, ONLY: &
    DescribeEOSInversionError

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

  INTEGER,  PARAMETER :: M_FP = 3
  REAL(DP), PARAMETER :: WFactor_FP = FourPi / PlanckConstant**3

  INTEGER  :: nE_G, nX_G, nZ(4), nX(3), n_FP
  INTEGER  :: iE_B0,    iE_E0
  INTEGER  :: iE_B1,    iE_E1
  INTEGER  :: iX_B0(3), iX_E0(3)
  INTEGER  :: iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:)
  REAL(DP), ALLOCATABLE :: W3_S(:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: dF_N(:,:)
  REAL(DP), ALLOCATABLE :: Chi(:,:,:)
  REAL(DP), ALLOCATABLE :: Sig(:,:,:)
  REAL(DP), ALLOCATABLE :: fEQ(:,:,:)
  REAL(DP), ALLOCATABLE :: Chi_NES(:,:,:)
  REAL(DP), ALLOCATABLE :: Eta_NES(:,:,:)
  REAL(DP), ALLOCATABLE :: Chi_Pair(:,:,:)
  REAL(DP), ALLOCATABLE :: Eta_Pair(:,:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dR_N(:,:,:,:)

  REAL(DP), ALLOCATABLE :: AMAT(:,:,:)
  REAL(DP), ALLOCATABLE :: BVEC(:,:)
  REAL(DP), ALLOCATABLE :: TAU (:,:)
  REAL(DP), ALLOCATABLE :: WORK(:,:)
  INTEGER,  ALLOCATABLE :: INFO(:)  
  INTEGER               :: LWORK

  INTERFACE ComputePrimitive_Euler
    MODULE PROCEDURE ComputePrimitive_Euler_Scalar
    MODULE PROCEDURE ComputePrimitive_Euler_Vector
  END INTERFACE

  INTERFACE ComputeConserved_Euler
    MODULE PROCEDURE ComputeConserved_Euler_Scalar
    MODULE PROCEDURE ComputeConserved_Euler_Vector
  END INTERFACE

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
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOF ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies)

    INTEGER  :: iX1, iX2, iX3, iGF, iCF, iCR, iS, iE, iN_E, iN_X
    INTEGER  :: iNode, iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: Chi_T, Eta_T, Eta, Kappa

    CALL TimersStart( Timer_Implicit )

    CALL InitializeCollisions_New( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( dU_F, dU_R )
#endif

    CALL TimersStart( Timer_Im_In )

    ! --- Copy inputs to locals ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, dF_N, U_F )
#elif defined(THORNADO_OMP)
#endif
    DO iCF = 1, nCF
      DO iN_X = 1, nX_G

        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        CF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)
        dF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)

      END DO
    END DO


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, GX_N, GX )
#elif defined(THORNADO_OMP)
#endif
    DO iGF = 1, nGF
      DO iN_X = 1, nX_G

        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

      END DO
    END DO

    ! --- Compute Primitive Quantities ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined(THORNADO_OMP)
#endif
    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeTS_Aux )

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D ), &
             PF_N(:,iPF_E ), &
             PF_N(:,iPF_Ne), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_E ), &
             AF_N(:,iAF_Ye) )

    CALL TimersStop( Timer_Im_ComputeTS_Aux )

    ! --- Opacity Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    DO iS = 1, nSpecies

      CALL ComputeNeutrinoOpacities_EC_Points &
             ( 1, nE_G, 1, nX_G, &
               E_N (:), &
               PF_N(:,iPF_D ), &
               AF_N(:,iAF_T ), &
               AF_N(:,iAF_Ye), &
               iS, Chi(:,:,iS) )

    END DO

    DO iS = 1, nSpecies

      CALL ComputeNeutrinoOpacities_ES_Points &
             ( 1, nE_G, 1, nX_G, &
               E_N (:), &
               PF_N(:,iPF_D ), &
               AF_N(:,iAF_T ), &
               AF_N(:,iAF_Ye), &
               iS, 1, Sig(:,:,iS) )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Chi_NES (iN_E,iN_X,iS) = Zero
          Eta_NES (iN_E,iN_X,iS) = Zero
          Chi_Pair(iN_E,iN_X,iS) = Zero
          Eta_Pair(iN_E,iN_X,iS) = Zero

        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Rearrange data to group energy nodes together ---

    CALL TimersStart( Timer_Im_MapForward )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNode, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNode, iNodeX, iNodeE, iE, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iN_X = 1, nX_G
          DO iN_E = 1, nE_G

            iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
            iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

            iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
            iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
            iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
            iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

            iNode  = iNodeE &
                     + ( iNodeX - 1 ) * nDOFE

            CR_N(iN_E,iN_X,iCR,iS) = U_R(iNode,iE,iX1,iX2,iX3,iCR,iS)

          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_MapForward )

    CALL TimersStop( Timer_Im_In )

    CALL TimersStart( Timer_Im_Solve )

    ! --- Electron Neutrinos and Antineutrinos ---

    !CALL SolveMatterEquations_EmAb &
    !       ( CR_N    (:,:,iCR_N,iNuE    ), &
    !         CR_N    (:,:,iCR_N,iNuE_Bar), &
    !         dt * Chi(:,:,      iNuE    ), &
    !         dt * Chi(:,:,      iNuE_Bar), &
    !         fEQ     (:,:,      iNuE    ), &
    !         fEQ     (:,:,      iNuE_Bar), &
    !         PF_N    (:,iPF_D ), &
    !         AF_N    (:,iAF_T ), &
    !         AF_N    (:,iAF_Ye), &
    !         AF_N    (:,iAF_E ) )

    CALL SolveMatterEquations_FP_Coupled &
           ( dt, iNuE, iNuE_Bar, &
             CR_N    (:,:,iCR_N,iNuE    ), &
             CR_N    (:,:,iCR_N,iNuE_Bar), &
             Chi     (:,:,      iNuE    ), &
             Chi     (:,:,      iNuE_Bar), &
             fEQ     (:,:,      iNuE    ), &
             fEQ     (:,:,      iNuE_Bar), &
             Chi_NES (:,:,      iNuE    ), &
             Chi_NES (:,:,      iNuE_Bar), &
             Eta_NES (:,:,      iNuE    ), &
             Eta_NES (:,:,      iNuE_Bar), &
             Chi_Pair(:,:,      iNuE    ), &
             Chi_Pair(:,:,      iNuE_Bar), &
             Eta_Pair(:,:,      iNuE    ), &
             Eta_Pair(:,:,      iNuE_Bar), &
             PF_N    (:,iPF_D ), &
             AF_N    (:,iAF_T ), &
             AF_N    (:,iAF_Ye), &
             AF_N    (:,iAF_E ) )

    CALL TimersStop( Timer_Im_Solve )

    CALL TimersStart( Timer_Im_Out )

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeTS_Prim )

    CALL ComputeThermodynamicStates_Primitive_TABLE &
           ( PF_N(:,iPF_D ), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_Ye), &
             PF_N(:,iPF_E ), &
             AF_N(:,iAF_E ), &
             PF_N(:,iPF_Ne) )

    CALL TimersStop( Timer_Im_ComputeTS_Prim )

    ! --- Compute Conserved Quantities ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( PF_N, CF_N, GX_N )
#elif defined(THORNADO_OMP)
#endif
    DO iN_X = 1, nX_G

      CALL ComputeConserved_Euler &
             ( PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

    CALL TimersStart( Timer_Im_Increment )

    ! --- Conserved Fluid Increment ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( nX, iX_B0, dU_F, CF_N, dF_N )
#elif defined(THORNADO_OMP)
#endif
    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX

              iN_X = iNodeX &
                     + ( iX1 - iX_B0(1) ) * nDOFX &
                     + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
                     + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

              dU_F(iNodeX,iX1,iX2,iX3,iCF) = ( CF_N(iN_X,iCF) - dF_N(iN_X,iCF) ) / dt

            END DO
          END DO
        END DO
      END DO
    END DO

    ! --- Update Radiation Fields ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Chi_T, Eta_T, Eta, Kappa )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Chi_T, Eta_T, Eta, Kappa ) &
    !$ACC PRESENT( Chi, Chi_NES, Chi_Pair, Eta_NES, Eta_Pair, &
    !$ACC          Sig, fEQ, CR_N, dR_N )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Chi_T = Chi(iN_E,iN_X,iS) + Chi_NES(iN_E,iN_X,iS) + Chi_Pair(iN_E,iN_X,iS)
          Kappa = Chi_T + Sig(iN_E,iN_X,iS)

          Eta   = Chi(iN_E,iN_X,iS) * fEQ(iN_E,iN_X,iS)
          Eta_T = Eta + Eta_NES(iN_E,iN_X,iS) + Eta_Pair(iN_E,iN_X,iS)

          ! --- Number Density ---

          CR_N(iN_E,iN_X,iCR_N,iS) &
            = ( CR_N(iN_E,iN_X,iCR_N,iS) + dt * Eta_T ) / ( One + dt * Chi_T )

          ! --- Number Flux (1) ---

          CR_N(iN_E,iN_X,iCR_G1,iS) &
            = CR_N(iN_E,iN_X,iCR_G1,iS) / ( One + dt * Kappa )

          ! --- Number Flux (2) ---

          CR_N(iN_E,iN_X,iCR_G2,iS) &
            = CR_N(iN_E,iN_X,iCR_G2,iS) / ( One + dt * Kappa )

          ! --- Number Flux (3) ---

          CR_N(iN_E,iN_X,iCR_G3,iS) &
            = CR_N(iN_E,iN_X,iCR_G3,iS) / ( One + dt * Kappa )

          ! --- Increments ---

          dR_N(iN_E,iN_X,iCR_N,iS) &
            = Eta_T - Chi_T * CR_N(iN_E,iN_X,iCR_N,iS)

          dR_N(iN_E,iN_X,iCR_G1,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iCR_G1,iS)

          dR_N(iN_E,iN_X,iCR_G2,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iCR_G2,iS)

          dR_N(iN_E,iN_X,iCR_G3,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iCR_G3,iS)
            
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_Increment )

    ! --- Rearrange data back to original layout ---

    CALL TimersStart( Timer_Im_MapBackward )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, dU_R, dR_N )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
            DO iX1 = iX_B0(1), iX_E0(1)
              DO iE = iE_B0, iE_E0
                DO iNode = 1, nDOF

                  iNodeX = MOD( (iNode-1) / nNodesE, nDOFX   ) + 1
                  iNodeE = MOD( (iNode-1)          , nNodesE ) + 1

                  iN_X = iNodeX &
                         + ( iX1 - iX_B0(1) ) * nDOFX &
                         + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
                         + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
                  iN_E = iNodeE &
                         + ( iE  - iE_B0    ) * nDOFE

                  dU_R(iNode,iE,iX1,iX2,iX3,iCR,iS) = dR_N(iN_E,iN_X,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_MapBackward )

    CALL TimersStop( Timer_Im_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_F, dU_R ) &
    !$OMP MAP( release: GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_F, dU_R ) &
    !$ACC DELETE( GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

    CALL FinalizeCollisions_New

    CALL TimersStop( Timer_Implicit )

#ifdef THORNADO_DEBUG_IMPLICIT
    WRITE(*,'(a,8x,5i4,es23.15)') 'MINLOC(dU_F), MINVAL(dU_F)', MINLOC(dU_F), MINVAL(dU_F)
    WRITE(*,'(a,8x,5i4,es23.15)') 'MAXLOC(dU_F), MAXVAL(dU_F)', MAXLOC(dU_F), MAXVAL(dU_F)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(dU_R), MINVAL(dU_R)', MINLOC(dU_R), MINVAL(dU_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(dU_R), MAXVAL(dU_R)', MAXLOC(dU_R), MAXVAL(dU_R)
#endif

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_New


  SUBROUTINE SolveMatterEquations_EmAb &
    ( J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

    ! --- Electron Neutrinos (1) and Electron Antineutrinos (2) ---

    REAL(DP), INTENT(in)    :: J_1  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J_2  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: D    (1:nX_G)
    REAL(DP), INTENT(in)    :: T    (1:nX_G)
    REAL(DP), INTENT(in)    :: Y    (1:nX_G)
    REAL(DP), INTENT(in)    :: E    (1:nX_G)

    REAL(DP), PARAMETER :: Log1d100 = LOG( 1.0d100 )
    REAL(DP) :: Mnu_1(1:nX_G), dMnudT_1(1:nX_G), dMnudY_1(1:nX_G)
    REAL(DP) :: Mnu_2(1:nX_G), dMnudT_2(1:nX_G), dMnudY_2(1:nX_G)
    REAL(DP) :: FD1_Exp, FD2_Exp

    INTEGER  :: iN_X, iN_E

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#endif

    ! --- Neutrino Chemical Potentials and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_1, dMnudT_1, dMnudY_1, iSpecies = iNuE )

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_2, dMnudT_2, dMnudY_2, iSpecies = iNuE_Bar )

    ! --- Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( FD1_Exp, FD2_Exp )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2), &
    !$ACC PRIVATE( FD1_Exp, FD2_Exp ) &
    !$ACC PRESENT( E_N, T, Mnu_1, Mnu_2, J0_1, J0_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( FD1_Exp, FD2_Exp )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        FD1_Exp = ( E_N(iN_E) - Mnu_1(iN_X) ) / ( BoltzmannConstant * T(iN_X) )
        FD1_Exp = MIN( MAX( FD1_Exp, - Log1d100 ), + Log1d100 )
        J0_1(iN_E,iN_X) = One / ( EXP( FD1_Exp ) + One )

        FD2_Exp = ( E_N(iN_E) - Mnu_2(iN_X) ) / ( BoltzmannConstant * T(iN_X) )
        FD2_Exp = MIN( MAX( FD2_Exp, - Log1d100 ), + Log1d100 )
        J0_2(iN_E,iN_X) = One / ( EXP( FD2_Exp ) + One )

      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#endif


  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE SolveMatterEquations_FP_Coupled &
    ( dt, iS_1, iS_2, Jold_1, Jold_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      D, T, Y, E )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(in)    :: Jold_1    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jold_2    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_1      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_2      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_Pair_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_Pair_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: D         (1:nX_G)
    REAL(DP), INTENT(inout) :: T         (1:nX_G)
    REAL(DP), INTENT(inout) :: Y         (1:nX_G)
    REAL(DP), INTENT(inout) :: E         (1:nX_G)

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = iE
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(              1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, FVECm_Y, GVECm_Y
    REAL(DP), DIMENSION(              1:nX_G) :: Eold, S_E, C_E, Unew_E, FVECm_E, GVECm_E
    REAL(DP), DIMENSION(1:nE_G,       1:nX_G) :: Jnew_1, FVECm_1, GVECm_1
    REAL(DP), DIMENSION(1:nE_G,       1:nX_G) :: Jnew_2, FVECm_2, GVECm_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

    INTEGER,  DIMENSION(       1:nX_G) :: PackedToUnpackedTable, UnpackedToPackedTable
    REAL(DP), DIMENSION(       1:nX_G) :: D_P, Y_P, T_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: J0_1_P, Jnew_1_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: J0_2_P, Jnew_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Eta_NES_1_P, Eta_NES_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Chi_NES_1_P, Chi_NES_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Eta_Pair_1_P, Eta_Pair_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Chi_Pair_1_P, Chi_Pair_2_P

    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    REAL(DP), DIMENSION(    1:nX_G) :: JNRM_1, JNRM_2
    LOGICAL,  DIMENSION(    1:nX_G) :: CONVERGED, ITERATE
    INTEGER,  DIMENSION(    1:nX_G) :: Error, K_FP

    REAL(DP) :: AERR_Y, AERR_E, AERR_J1, AERR_J2
    REAL(DP) :: RERR_Y, RERR_E, RERR_J1, RERR_J2
    REAL(DP) :: SUM1, SUM2
    REAL(DP) :: Eta_1, Eta_2
    REAL(DP) :: Eta_Total_1, Eta_Total_2, Chi_Total_1, Chi_Total_2
    INTEGER  :: i, k, iFP, iM, Mk, iN_X, iN_E, iX_P, nX_P
    INTEGER  :: OS_2

    OS_2 = OS_1 + nE_G
    ITERATE(:) = .TRUE.
    CONVERGED(:) = .FALSE.
    K_FP(:) = 0

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: CONVERGED, ITERATE, K_FP ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, FVECm_Y, GVECm_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, FVECm_E, GVECm_E, &
    !$OMP             Jnew_1, FVECm_1, GVECm_1, JNRM_1, &
    !$OMP             Jnew_2, FVECm_2, GVECm_2, JNRM_2, &
    !$OMP             Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP             Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP             PackedToUnpackedTable, UnpackedToPackedTable, &
    !$OMP             D_P, Y_P, T_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$OMP             Eta_NES_1_P, Eta_NES_2_P, &
    !$OMP             Chi_NES_1_P, Chi_NES_2_P, &
    !$OMP             Eta_Pair_1_P, Eta_Pair_2_P, &
    !$OMP             Chi_Pair_1_P, Chi_Pair_2_P, &
    !$OMP             GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( CONVERGED, ITERATE, K_FP ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, FVECm_Y, GVECm_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, FVECm_E, GVECm_E, &
    !$ACC         Jnew_1, FVECm_1, GVECm_1, JNRM_1, &
    !$ACC         Jnew_2, FVECm_2, GVECm_2, JNRM_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackedToUnpackedTable, UnpackedToPackedTable, &
    !$ACC         D_P, Y_P, T_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$ACC         Eta_NES_1_P, Eta_NES_2_P, &
    !$ACC         Chi_NES_1_P, Chi_NES_2_P, &
    !$ACC         Eta_Pair_1_P, Eta_Pair_2_P, &
    !$ACC         Chi_Pair_1_P, Chi_Pair_2_P, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM1, SUM2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM1, SUM2 ) &
    !$ACC PRESENT( JNRM_1, JNRM_2, Jold_1, Jold_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1, SUM2 )
#endif
    DO iN_X = 1, nX_G

      SUM1 = Zero
      SUM2 = Zero
#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP REDUCTION( +: SUM1, SUM2 )
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR &
      !$ACC REDUCTION( +: SUM1, SUM2 )
#endif
      DO iN_E = 1, nE_G
        SUM1 = SUM1 + Jold_1(iN_E,iN_X) * Jold_1(iN_E,iN_X)
        SUM2 = SUM2 + Jold_2(iN_E,iN_X) * Jold_2(iN_E,iN_X)
      END DO
      JNRM_1(iN_X) = SQRT( SUM1 )
      JNRM_2(iN_X) = SQRT( SUM2 )
      !CALL VectorNorm2( nE_G, Jold_1(1,iN_X), 1, JNRM_1(iN_X) )
      !CALL VectorNorm2( nE_G, Jold_2(1,iN_X), 1, JNRM_2(iN_X) )

    END DO

    ! --- Equilibrium Distributions ---

    CALL TimersStart( Timer_FP_FD )

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, J0_1, iS_1 )

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, J0_2, iS_2 )

    CALL TimersStop( Timer_FP_FD )

    ! --- Compute Kernels ---

    CALL TimersStart( Timer_FP_Op )

    ! --- NES Kernels ---

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1 )

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_NES_2, Phi_0_Ot_NES_2 )

    ! --- Pair Kernels ---

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1 )

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

    CALL TimersStop( Timer_FP_Op )

    ! --- NES Emissivities and Opacities ---

    CALL TimersStart( Timer_FP_Rate )

    CALL ComputeNeutrinoOpacitiesRates_NES_Points &
           ( 1, nE_G, 1, nX_G, W2_N, Jold_1, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1, Eta_NES_1, Chi_NES_1 )

    CALL ComputeNeutrinoOpacitiesRates_NES_Points &
           ( 1, nE_G, 1, nX_G, W2_N, Jold_2, &
             Phi_0_In_NES_2, Phi_0_Ot_NES_2, Eta_NES_2, Chi_NES_2 )

    ! --- Pair Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
           ( 1, nE_G, 1, nX_G, W2_N, Jold_2, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Eta_Pair_1, Chi_Pair_1 )

    CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
           ( 1, nE_G, 1, nX_G, W2_N, Jold_1, &
             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, Eta_Pair_2, Chi_Pair_2 )

    CALL TimersStop( Timer_FP_Rate )

    ! --- Right-Hand Side Vectors ---

    CALL TimersStart( Timer_FP_RHS )

    ! --- Fixed Electron Fraction and Internal Energy RHS ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W2_S, 1, Zero, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, Jold_2, nE_G, W2_S, 1,  One, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W3_S, 1, Zero, C_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_2, nE_G, W3_S, 1,  One, C_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( Y, Yold, E, Eold, D, S_Y, S_E, C_Y, C_E, Unew_Y, Unew_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      Yold(iN_X) = Y(iN_X)
      Eold(iN_X) = E(iN_X)

      S_Y(iN_X) = One / ( D(iN_X) * Yold(iN_X) / AtomicMassUnit )
      S_E(iN_X) = One / ( D(iN_X) * Eold(iN_X) )

      C_Y(iN_X) = C_Y(iN_X) * S_Y(iN_X)
      C_E(iN_X) = C_E(iN_X) * S_E(iN_X)

      Unew_Y(iN_X) = One ! --- Initial Guess
      Unew_E(iN_X) = One ! --- Initial Guess

    END DO

    ! --- Update Neutrino Densities ---

    CALL BuildRHS_J_FP &
           ( ITERATE, dt, Jold_1, Jnew_1, J0_1, Chi_1, Chi_NES_1, &
             Eta_NES_1, Chi_Pair_1, Eta_Pair_1, FVECm_1, Jnew_1 )

    CALL BuildRHS_J_FP &
           ( ITERATE, dt, Jold_2, Jnew_2, J0_2, Chi_2, Chi_NES_2, &
             Eta_NES_2, Chi_Pair_2, Eta_Pair_2, FVECm_2, Jnew_2 )

    CALL TimersStop( Timer_FP_RHS )

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1
      Mk = MIN( M_FP, k )
      iM = Mk
      !iM = 1 + MOD( k-1, M_FP )

      CALL CreatePackTable_FP &
             ( ITERATE, nX_P, PackedToUnpackedTable, UnpackedToPackedTable )

      IF ( k > 1 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackedToPackedTable, D, T, Y, D_P, T_P, Y_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#endif
        DO iX_P = 1, nX_P
          iN_X = UnpackedToPackedTable(iX_P)
          D_P(iX_P) = D(iN_X)
          T_P(iX_P) = T(iN_X)
          Y_P(iX_P) = Y(iN_X)
        END DO

        ! --- Recompute Equilibrium Distributions ---

        CALL TimersStart( Timer_FP_FDK )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackedToPackedTable, J0_1, J0_2, J0_1_P, J0_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#endif
        DO iX_P = 1, nX_P
          DO iN_E = 1, nE_G
            iN_X = UnpackedToPackedTable(iX_P)
            J0_1_P(iN_E,iX_P) = J0_1(iN_E,iN_X)
            J0_2_P(iN_E,iX_P) = J0_2(iN_E,iN_X)
          END DO
        END DO

        CALL ComputeEquilibriumDistributions_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), &
                 J0_1_P(:,1:nX_P), iS_1 )

        CALL ComputeEquilibriumDistributions_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), &
                 J0_2_P(:,1:nX_P), iS_2 )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iX_P ) &
      !$ACC PRESENT( ITERATE, PackedToUnpackedTable, J0_1, J0_2, J0_1_P, J0_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#endif
        DO iN_X = 1, nX_G
          DO iN_E = 1, nE_G
            IF ( ITERATE(iN_X) ) THEN
              iX_P = PackedToUnpackedTable(iN_X)
              J0_1(iN_E,iN_X) = J0_1_P(iN_E,iX_P)
              J0_2(iN_E,iN_X) = J0_2_P(iN_E,iX_P)
            END IF
          END DO
        END DO

        CALL TimersStop( Timer_FP_FDK )

        ! --- Recompute Kernels ---

        CALL TimersStart( Timer_FP_OpK )

        ! --- NES Kernels ---

        CALL ComputeNeutrinoOpacities_NES_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_1, 1, &
                 Phi_0_In_NES_1(:,:,1:nX_P), Phi_0_Ot_NES_1(:,:,1:nX_P) )

        CALL ComputeNeutrinoOpacities_NES_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_2, 1, &
                 Phi_0_In_NES_2(:,:,1:nX_P), Phi_0_Ot_NES_2(:,:,1:nX_P) )

        ! --- Pair Kernels ---

        CALL ComputeNeutrinoOpacities_Pair_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_1, 1, &
                 Phi_0_In_Pair_1(:,:,1:nX_P), Phi_0_Ot_Pair_1(:,:,1:nX_P) )

        CALL ComputeNeutrinoOpacities_Pair_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_2, 1, &
                 Phi_0_In_Pair_2(:,:,1:nX_P), Phi_0_Ot_Pair_2(:,:,1:nX_P) )

        CALL TimersStop( Timer_FP_OpK )

      END IF

      ! --- NES Emissivities and Opacities ---

      CALL TimersStart( Timer_FP_RateK )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackedToPackedTable, Jnew_1, Jnew_2, Jnew_1_P, Jnew_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#endif
      DO iX_P = 1, nX_P
        DO iN_E = 1, nE_G
          iN_X = UnpackedToPackedTable(iX_P)
          Jnew_1_P(iN_E,iX_P) = Jnew_1(iN_E,iN_X)
          Jnew_2_P(iN_E,iX_P) = Jnew_2(iN_E,iN_X)
        END DO
      END DO

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_1_P(:,1:nX_P), &
               Phi_0_In_NES_1(:,:,1:nX_P), Phi_0_Ot_NES_1(:,:,1:nX_P), &
               Eta_NES_1_P(:,1:nX_P), Chi_NES_1_P(:,1:nX_P) )

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_2_P(:,1:nX_P), &
               Phi_0_In_NES_2(:,:,1:nX_P), Phi_0_Ot_NES_2(:,:,1:nX_P), &
               Eta_NES_2_P(:,1:nX_P), Chi_NES_2_P(:,1:nX_P) )

      ! --- Pair Emissivities and Opacities ---

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_2_P(:,1:nX_P), &
               Phi_0_In_Pair_1(:,:,1:nX_P), Phi_0_Ot_Pair_1(:,:,1:nX_P), &
               Eta_Pair_1_P(:,1:nX_P), Chi_Pair_1_P(:,1:nX_P) )

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_1_P(:,1:nX_P), &
               Phi_0_In_Pair_2(:,:,1:nX_P), Phi_0_Ot_Pair_2(:,:,1:nX_P), &
               Eta_Pair_2_P(:,1:nX_P), Chi_Pair_2_P(:,1:nX_P) )

      CALL TimersStop( Timer_FP_RateK )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iX_P ) &
      !$ACC PRESENT( ITERATE, PackedToUnpackedTable, &
      !$ACC          Chi_NES_1, Chi_NES_2, Chi_NES_1_P, Chi_NES_2_P, &
      !$ACC          Eta_NES_1, Eta_NES_2, Eta_NES_1_P, Eta_NES_2_P, &
      !$ACC          Chi_Pair_1, Chi_Pair_2, Chi_Pair_1_P, Chi_Pair_2_P, &
      !$ACC          Eta_Pair_1, Eta_Pair_2, Eta_Pair_1_P, Eta_Pair_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#endif
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G
          IF ( ITERATE(iN_X) ) THEN
            iX_P = PackedToUnpackedTable(iN_X)
            Chi_NES_1(iN_E,iN_X) = Chi_NES_1_P(iN_E,iX_P)
            Chi_NES_2(iN_E,iN_X) = Chi_NES_2_P(iN_E,iX_P)
            Eta_NES_1(iN_E,iN_X) = Eta_NES_1_P(iN_E,iX_P)
            Eta_NES_2(iN_E,iN_X) = Eta_NES_2_P(iN_E,iX_P)
            Chi_Pair_1(iN_E,iN_X) = Chi_Pair_1_P(iN_E,iX_P)
            Chi_Pair_2(iN_E,iN_X) = Chi_Pair_2_P(iN_E,iX_P)
            Eta_Pair_1(iN_E,iN_X) = Eta_Pair_1_P(iN_E,iX_P)
            Eta_Pair_2(iN_E,iN_X) = Eta_Pair_2_P(iN_E,iX_P)
          END IF
        END DO
      END DO

      ! --- Right-Hand Side Vectors and Residuals ---

      CALL TimersStart( Timer_FP_RHSK )

      CALL BuildRHS_Y_FP &
             ( ITERATE, Jnew_1, Jnew_2, Unew_Y, C_Y, S_Y, FVECm_Y, GVECm_Y )
                                                 
      CALL BuildRHS_E_FP &                      
             ( ITERATE, Jnew_1, Jnew_2, Unew_E, C_E, S_E, FVECm_E, GVECm_E )

      CALL BuildRHS_J_FP &
             ( ITERATE, dt, Jold_1, Jnew_1, J0_1, Chi_1, &
               Chi_NES_1, Eta_NES_1, Chi_Pair_1, Eta_Pair_1, &
               FVECm_1, GVECm_1 )

      CALL BuildRHS_J_FP &
             ( ITERATE, dt, Jold_2, Jnew_2, J0_2, Chi_2, &
               Chi_NES_2, Eta_NES_2, Chi_Pair_2, Eta_Pair_2, &
               FVECm_2, GVECm_2 )

      CALL CombineRHS_FP &
             ( ITERATE, Mk, iY, iE, OS_1, OS_2, &
               FVECm_Y, FVECm_E, FVECm_1, FVECm_2, &
               GVECm_Y, GVECm_E, GVECm_1, GVECm_2, &
               FVECm, GVECm, FVEC, GVEC )

      CALL TimersStop( Timer_FP_RHSK )

      ! --- Update Solution ---

      CALL TimersStart( Timer_FP_AA )

      IF ( Mk > 1 ) THEN

        ! --- Anderson Acceleration ---

        CALL BuildLS_FP &
               ( ITERATE, Mk, FVECm, FVEC, AMAT, BVEC )

        CALL TimersStart( Timer_FP_LS )

        CALL SolveLS_FP &
              ( ITERATE, Mk, AMAT, BVEC, Alpha )

        CALL TimersStop( Timer_FP_LS )

        CALL UpdateRHS_AA_FP &
               ( ITERATE, Mk, Alpha, GVEC, GVECm )

      END IF

      CALL TimersStop( Timer_FP_AA )

      CALL TimersStart( Timer_FP_Update )

      ! --- Update Residuals and Solution Vectors---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( ITERATE, FVECm, GVECm, Y, E, Yold, Eold, &
      !$ACC          Unew_Y, Unew_E, Jnew_1, Jnew_2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          FVECm(iY,iN_X) = GVECm(iY,iN_X) - Unew_Y(iN_X)
          FVECm(iE,iN_X) = GVECm(iE,iN_X) - Unew_E(iN_X)

          DO iFP = OS_1+1, OS_1+nE_G
            FVECm(iFP,iN_X) = GVECm(iFP,iN_X) - Jnew_1(iFP-OS_1,iN_X)
          END DO

          DO iFP = OS_2+1, OS_2+nE_G
            FVECm(iFP,iN_X) = GVECm(iFP,iN_X) - Jnew_2(iFP-OS_2,iN_X)
          END DO

          Unew_Y(iN_X) = GVECm(iY,iN_X)
          Y(iN_X) = Unew_Y(iN_X) * Yold(iN_X)

          Unew_E(iN_X) = GVECm(iE,iN_X)
          E(iN_X) = Unew_E(iN_X) * Eold(iN_X)

          DO iN_E = 1, nE_G
            Jnew_1(iN_E,iN_X) = GVECm(iN_E+OS_1,iN_X)
            Jnew_2(iN_E,iN_X) = GVECm(iN_E+OS_2,iN_X)
          END DO

        END IF
      END DO

      ! --- Update Temperature ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !$OMP MAP( from: Error )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYOUT( Error ) &
      !$ACC PRESENT( ITERATE, D, E, Y, T )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          CALL ComputeTemperatureFromSpecificInternalEnergyPoint_TABLE &
                 ( D(iN_X), E(iN_X), Y(iN_X), T(iN_X), &
                   Guess_Option = T(iN_X), Error_Option = Error(iN_X) )

        ELSE

          Error(iN_X) = 0

        END IF
      END DO

      IF ( ANY( Error(:) > 0 ) ) THEN
        DO iN_X = 1, nX_G
          IF ( Error(iN_X) > 0 ) THEN
#if defined(THORNADO_OMP_OL)
            !$OMP TARGET UPDATE FROM( Eold(iN_X), E(iN_X), Yold(iN_X), Y(iN_X) )
#elif defined(THORNADO_OACC)
            !$ACC UPDATE HOST( Eold(iN_X), E(iN_X), Yold(iN_X), Y(iN_X) )
#endif
            CALL DescribeEOSInversionError( Error(iN_X) )
            WRITE(*,'(2i4,5es12.5)') k, iN_X, D(iN_X), E(iN_X), Y(iN_X), Eold(iN_X), Yold(iN_X)
          END IF
        END DO
        STOP
      END IF
!      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
!             ( D, E, Y, T )

      ! --- Check Convergence ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !$OMP PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$OMP          RERR_Y, RERR_E, RERR_J1, RERR_J2 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$ACC          RERR_Y, RERR_E, RERR_J1, RERR_J2 ) &
      !$ACC PRESENT( ITERATE, CONVERGED, K_FP, JNRM_1, JNRM_2, FVECm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$OMP          RERR_Y, RERR_E, RERR_J1, RERR_J2 )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          AERR_Y = ABS( FVECm(iY,iN_X) )
          AERR_E = ABS( FVECm(iE,iN_X) )

          AERR_J1 = Zero
          AERR_J2 = Zero
          DO iN_E = 1, nE_G
            AERR_J1 = AERR_J1 + FVECm(OS_1+iN_E,iN_X)**2
            AERR_J2 = AERR_J2 + FVECm(OS_2+iN_E,iN_X)**2
          END DO
          AERR_J1 = SQRT( AERR_J1 )
          AERR_J2 = SQRT( AERR_J2 )

          RERR_Y  = AERR_Y
          RERR_E  = AERR_E
          RERR_J1 = AERR_J1 / JNRM_1(iN_X)
          RERR_J2 = AERR_J2 / JNRM_2(iN_X)

          CONVERGED(iN_X) = RERR_Y  <= Rtol &
                      .AND. RERR_E  <= Rtol &
                      .AND. RERR_J1 <= Rtol &
                      .AND. RERR_J2 <= Rtol

          ITERATE(iN_X) = .NOT. CONVERGED(iN_X)

          IF ( CONVERGED(iN_X) ) K_FP(iN_X) = k

        END IF
      END DO

      CALL TimersStop( Timer_FP_Update )

      ! --- Shift History Arrays ---

      IF ( Mk == M_FP ) THEN

        CALL TimersStart( Timer_FP_Shift )

        CALL ShiftRHS_FP &
               ( ITERATE, Mk, FVEC, GVEC )

        CALL TimersStop( Timer_FP_Shift )

      END IF

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( ITERATE, CONVERGED )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( ITERATE, CONVERGED )
#endif

    END DO

#ifdef THORNADO_DEBUG_IMPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( K_FP )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( K_FP )
#endif
    WRITE(*,'(a,2i4)') 'MINLOC(K_FP), MINVAL(K_FP)', MINLOC(K_FP), MINVAL(K_FP)
    WRITE(*,'(a,2i4)') 'MAXLOC(K_FP), MAXVAL(K_FP)', MAXLOC(K_FP), MAXVAL(K_FP)
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: CONVERGED, ITERATE, K_FP, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, FVECm_Y, GVECm_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, FVECm_E, GVECm_E, &
    !$OMP               Jnew_1, FVECm_1, GVECm_1, JNRM_1, &
    !$OMP               Jnew_2, FVECm_2, GVECm_2, JNRM_2, &
    !$OMP               Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP               Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP               PackedToUnpackedTable, UnpackedToPackedTable, &
    !$OMP               D_P, Y_P, T_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$OMP               Eta_NES_1_P, Eta_NES_2_P, &
    !$OMP               Chi_NES_1_P, Chi_NES_2_P, &
    !$OMP               Eta_Pair_1_P, Eta_Pair_2_P, &
    !$OMP               Chi_Pair_1_P, Chi_Pair_2_P, &
    !$OMP               GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( CONVERGED, ITERATE, K_FP, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, FVECm_Y, GVECm_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, FVECm_E, GVECm_E, &
    !$ACC         Jnew_1, FVECm_1, GVECm_1, JNRM_1, &
    !$ACC         Jnew_2, FVECm_2, GVECm_2, JNRM_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackedToUnpackedTable, UnpackedToPackedTable, &
    !$ACC         D_P, Y_P, T_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$ACC         Eta_NES_1_P, Eta_NES_2_P, &
    !$ACC         Chi_NES_1_P, Chi_NES_2_P, &
    !$ACC         Eta_Pair_1_P, Eta_Pair_2_P, &
    !$ACC         Chi_Pair_1_P, Chi_Pair_2_P, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

  END SUBROUTINE SolveMatterEquations_FP_Coupled


  SUBROUTINE CreatePackTable_FP &
    ( MASK, nX_P, PackedToUnpackedTable, UnpackedToPackedTable )

    LOGICAL,  INTENT(in)  :: MASK(1:nX_G)
    INTEGER,  INTENT(out) :: nX_P
    INTEGER,  INTENT(out) :: PackedToUnpackedTable(1:nX_G)
    INTEGER,  INTENT(out) :: UnpackedToPackedTable(1:nX_G)

    INTEGER  :: iN_X, iX_P, iX_U

    ! --- Build Lookup Tables ---

    nX_P = COUNT( MASK(:) )
    iX_P = 0
    iX_U = nX_P
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN
        iX_P = iX_P + 1
        PackedToUnpackedTable(iN_X) = iX_P
        UnpackedToPackedTable(iX_P) = iN_X
      ELSE
        iX_U = iX_U + 1
        PackedToUnpackedTable(iN_X) = iX_U
        UnpackedToPackedTable(iX_U) = iN_X
      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( PackedToUnpackedTable, UnpackedToPackedTable )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( PackedToUnpackedTable, UnpackedToPackedTable )
#endif

  END SUBROUTINE CreatePackTable_FP


  SUBROUTINE BuildRHS_Y_FP &
    ( MASK, J_1, J_2, U_Y, C_Y, S_Y, F_Y, G_Y )

    LOGICAL,  INTENT(in)    :: MASK      (1:nX_G)
    REAL(DP), INTENT(in)    :: J_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: U_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: C_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: S_Y       (1:nX_G)
    REAL(DP), INTENT(inout) :: F_Y       (1:nX_G)
    REAL(DP), INTENT(inout) :: G_Y       (1:nX_G)

    INTEGER :: iN_X

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_1, nE_G, W2_S, 1, Zero, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, J_2, nE_G, W2_S, 1,  One, G_Y, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, U_Y, C_Y, S_Y, F_Y, G_Y )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_Y(iN_X) = One + C_Y(iN_X) - G_Y(iN_X) * S_Y(iN_X)
        F_Y(iN_X) = G_Y(iN_X) - U_Y(iN_X)

      END IF
    END DO

  END SUBROUTINE BuildRHS_Y_FP


  SUBROUTINE BuildRHS_E_FP &
    ( MASK, J_1, J_2, U_E, C_E, S_E, F_E, G_E )

    LOGICAL,  INTENT(in)    :: MASK      (1:nX_G)
    REAL(DP), INTENT(in)    :: J_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: U_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: C_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: S_E       (1:nX_G)
    REAL(DP), INTENT(inout) :: F_E       (1:nX_G)
    REAL(DP), INTENT(inout) :: G_E       (1:nX_G)

    INTEGER :: iN_X

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_1, nE_G, W3_S, 1, Zero, G_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_2, nE_G, W3_S, 1,  One, G_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, U_E, C_E, S_E, F_E, G_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_E(iN_X) = One + C_E(iN_X) - G_E(iN_X) * S_E(iN_X)
        F_E(iN_X) = G_E(iN_X) - U_E(iN_X)

      END IF
    END DO

  END SUBROUTINE BuildRHS_E_FP


  SUBROUTINE BuildRHS_J_FP &
    ( MASK, dt, Jold, Jnew, J0, Chi, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, F_J, G_J )

    LOGICAL,  INTENT(in)    :: MASK           (1:nX_G)
    REAL(DP), INTENT(in)    :: dt
    REAL(DP), INTENT(in)    :: Jold    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jnew    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J0      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_NES (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_NES (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_Pair(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_Pair(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: F_J     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: G_J     (1:nE_G,1:nX_G)

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T ) &
    !$ACC PRESENT( MASK, Jold, Jnew, J0, Chi, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, F_J, G_J )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Eta = Chi(iN_E,iN_X) * J0(iN_E,iN_X)

          Eta_T = Eta + Eta_NES(iN_E,iN_X) + Eta_Pair(iN_E,iN_X)
          Chi_T = Chi(iN_E,iN_X) + Chi_NES(iN_E,iN_X) + Chi_Pair(iN_E,iN_X)

          G_J(iN_E,iN_X) = ( Jold(iN_E,iN_X) + dt * Eta_T ) / ( One + dt * Chi_T )
          F_J(iN_E,iN_X) = G_J(iN_E,iN_X) - Jnew(iN_E,iN_X)

        END IF
      END DO
    END DO

  END SUBROUTINE BuildRHS_J_FP


  SUBROUTINE CombineRHS_FP &
    ( MASK, Mk, iY, iE, OS_1, OS_2, F_Y, F_E, F_1, F_2, G_Y, G_E, G_1, G_2, Fm, Gm, F, G )

    LOGICAL,  INTENT(in)    :: MASK      (1:nX_G)
    INTEGER,  INTENT(in)    :: Mk, iY, iE, OS_1, OS_2
    REAL(DP), INTENT(in)    :: F_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: F_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: F_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: F_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: G_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: G_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: G_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: G_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Fm (1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm (1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: F  (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: G  (1:n_FP,1:M_FP,1:nX_G)

    INTEGER  :: iN_E, iN_X, iFP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, F_Y, F_E, G_Y, G_E, Fm, Gm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fm(iY,iN_X) = F_Y(iN_X)
        Fm(iE,iN_X) = F_E(iN_X)

        Gm(iY,iN_X) = G_Y(iN_X)
        Gm(iE,iN_X) = G_E(iN_X)

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, F_1, F_2, G_1, G_2, Fm, Gm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Fm(OS_1+iN_E,iN_X) = F_1(iN_E,iN_X)
          Fm(OS_2+iN_E,iN_X) = F_2(iN_E,iN_X)

          Gm(OS_1+iN_E,iN_X) = G_1(iN_E,iN_X)
          Gm(OS_2+iN_E,iN_X) = G_2(iN_E,iN_X)

        END IF
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, Gm, F, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN
          F(iFP,Mk,iN_X) = Fm(iFP,iN_X)
          G(iFP,Mk,iN_X) = Gm(iFP,iN_X)
        END IF
      END DO
    END DO

  END SUBROUTINE CombineRHS_FP


  SUBROUTINE BuildLS_FP &
    ( MASK, Mk, Fm, F, AMAT, BVEC )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(in)    :: Fm  (1:n_FP,       1:nX_G)
    REAL(DP), INTENT(in)    :: F   (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: AMAT(1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: BVEC(1:n_FP,       1:nX_G)

    INTEGER  :: iN_X, iFP, i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, BVEC )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN
          BVEC(iFP,iN_X) = - Fm(iFP,iN_X)
        END IF
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( MASK, F, Fm, AMAT )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
      DO i = 1, Mk-1
        DO iFP = 1, n_FP
          IF ( MASK(iN_X) ) THEN
            AMAT(iFP,i,iN_X) = F(iFP,i,iN_X) - Fm(iFP,iN_X)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE BuildLS_FP


  SUBROUTINE SolveLS_FP &
    ( MASK, Mk, AMAT, BVEC, Alpha )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: AMAT (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: BVEC (1:n_FP,       1:nX_G)
    REAL(DP), INTENT(inout) :: Alpha(       1:M_FP,1:nX_G)

    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1
    INTEGER  :: iN_X, iFP, i

    IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( AA11, AB1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( AA11, AB1 ) &
      !$ACC PRESENT( MASK, AMAT, BVEC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( AA11, AB1 )
#endif
      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          AA11 = Zero
          AB1 = Zero

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD &
          !$OMP REDUCTION( +: AA11, AB1 )
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR &
          !$ACC REDUCTION( +: AA11, AB1 )
#endif
          DO iFP = 1, n_FP

            AA11 = AA11 + AMAT(iFP,1,iN_X) * AMAT(iFP,1,iN_X)
            AB1  = AB1  + AMAT(iFP,1,iN_X) * BVEC(iFP,iN_X)

          END DO

          BVEC(1,iN_X) = AB1 / AA11

        END IF
      END DO

    ELSE IF ( Mk == 3 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA ) &
      !$ACC PRESENT( MASK, AMAT, BVEC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#endif
      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          AA11 = Zero
          AA12 = Zero
          AA22 = Zero
          AB1  = Zero
          AB2  = Zero

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD &
          !$OMP REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR &
          !$ACC REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#endif
          DO iFP = 1, n_FP

            AA11 = AA11 + AMAT(iFP,1,iN_X) * AMAT(iFP,1,iN_X)
            AA12 = AA12 + AMAT(iFP,1,iN_X) * AMAT(iFP,2,iN_X)
            AA22 = AA22 + AMAT(iFP,2,iN_X) * AMAT(iFP,2,iN_X)

            AB1  = AB1  + AMAT(iFP,1,iN_X) * BVEC(iFP,iN_X)
            AB2  = AB2  + AMAT(iFP,2,iN_X) * BVEC(iFP,iN_X)

          END DO

          DET_AA = AA11*AA22 - AA12*AA12

          BVEC(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
          BVEC(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

        END IF
      END DO

    ELSE

      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          CALL LinearLeastSquares &
            ( 'N', n_FP, Mk-1, 1, AMAT(1,1,iN_X), n_FP, &
              BVEC(1,iN_X), n_FP, TAU(1,iN_X), WORK(1,iN_X), LWORK, INFO(iN_X) )

        END IF
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( SUM1 ) &
    !$ACC PRESENT( MASK, Alpha, BVEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        SUM1 = Zero
        DO i = 1, Mk-1
          Alpha(i,iN_X) = BVEC(i,iN_X)
          SUM1 = SUM1 + BVEC(i,iN_X)
        END DO
        Alpha(Mk,iN_X) = One - SUM1

      END IF
    END DO

  END SUBROUTINE SolveLS_FP


  SUBROUTINE UpdateRHS_AA_FP &
    ( Mask, Mk, Alpha, G, Gm )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(in)    :: Alpha(       1:M_FP,1:nX_G)
    REAL(DP), INTENT(in)    :: G    (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm   (1:n_FP,       1:nX_G)

    REAL(DP) :: SUM1
    INTEGER  :: iN_X, iFP, i

    !CALL MatrixMatrixMultiplyBatched &
    !       ( 'N', 'N', n_FP, 1, Mk, Zero, G, n_FP, n_FP*M_FP, &
    !         Alpha, M_FP, M_FP, Gm, n_FP, n_FP, nX_G )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
    !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( SUM1 ) &
    !$ACC PRESENT( MASK, G, Gm, Alpha )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( SUM1 )
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN

          SUM1 = Zero
          DO i = 1, Mk
            SUM1 = SUM1 + G(iFP,i,iN_X) * Alpha(i,iN_X)
          END DO
          Gm(iFP,iN_X) = SUM1

        END IF
      END DO
    END DO

  END SUBROUTINE UpdateRHS_AA_FP


  SUBROUTINE ShiftRHS_FP &
    ( MASK, Mk, F, G )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: F(1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: G(1:n_FP,1:M_FP,1:nX_G)

    REAL(DP) :: FTMP(1:n_FP,1:M_FP)
    REAL(DP) :: GTMP(1:n_FP,1:M_FP)
    INTEGER  :: iN_X, iFP, i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( FTMP, GTMP )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( FTMP, GTMP ) &
    !$ACC PRESENT( MASK, F, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( FTMP, GTMP )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
        !$ACC LOOP VECTOR COLLAPSE(2)
#endif
        DO i = 1, Mk-1
          DO iFP = 1, n_FP
            FTMP(iFP,i) = F(iFP,i+1,iN_X)
            GTMP(iFP,i) = G(iFP,i+1,iN_X)
          END DO
        END DO

#if defined(THORNADO_OMP_OL)
        !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
        !$ACC LOOP VECTOR COLLAPSE(2)
#endif
        DO i = 1, Mk-1
          DO iFP = 1, n_FP
            F(iFP,i,iN_X) = FTMP(iFP,i)
            G(iFP,i,iN_X) = GTMP(iFP,i)
          END DO
        END DO

      END IF
    END DO

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE ComputeNeutrinoChemicalPotentials &
    ( D, T, Y, M, dMdT, dMdY, iSpecies )

    REAL(DP), INTENT(in)  :: D   (1:nX_G)
    REAL(DP), INTENT(in)  :: T   (1:nX_G)
    REAL(DP), INTENT(in)  :: Y   (1:nX_G)
    REAL(DP), INTENT(out) :: M   (1:nX_G)
    REAL(DP), INTENT(out) :: dMdT(1:nX_G)
    REAL(DP), INTENT(out) :: dMdY(1:nX_G)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: dMdD(1:nX_G)
    REAL(DP) :: Me(1:nX_G), dMedT(1:nX_G), dMedY(1:nX_G)
    REAL(DP) :: Mp(1:nX_G), dMpdT(1:nX_G), dMpdY(1:nX_G)
    REAL(DP) :: Mn(1:nX_G), dMndT(1:nX_G), dMndY(1:nX_G)

    INTEGER :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#endif

    ! --- Matter Chemical Potentials and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, M = Me, dMdD_Option = dMdD, dMdT_Option = dMedT, dMdY_Option = dMedY )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, M = Mp, dMdD_Option = dMdD, dMdT_Option = dMpdT, dMdY_Option = dMpdY )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, M = Mn, dMdD_Option = dMdD, dMdT_Option = dMndT, dMdY_Option = dMndY )

    ! --- Neutrino Chemical Potential and Derivatives ---

    IF ( iSpecies == iNuE )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( M, dMdT, dMdY, Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO i = 1, nX_G
        M   (i) = ( Me   (i) + Mp   (i) ) - Mn   (i)
        dMdT(i) = ( dMedT(i) + dMpdT(i) ) - dMndT(i)
        dMdY(i) = ( dMedY(i) + dMpdY(i) ) - dMndY(i)
      END DO

    ELSE IF ( iSpecies == iNuE_Bar )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( M, dMdT, dMdY, Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO i = 1, nX_G
        M   (i) = Mn   (i) - ( Me   (i) + Mp   (i) )
        dMdT(i) = dMndT(i) - ( dMedT(i) + dMpdT(i) )
        dMdY(i) = dMndY(i) - ( dMedY(i) + dMpdY(i) )
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#endif

  END SUBROUTINE ComputeNeutrinoChemicalPotentials


  SUBROUTINE InitializeCollisions_New( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)

    REAL(DP) :: TMP(1)

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nX_G = nDOFX * PRODUCT( nX )
    nE_G = nNodesZ(1) * nZ(1)
    n_FP = 2 + 2*nE_G

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )
    ALLOCATE( GX_N(nX_G,nGF) )
    ALLOCATE( dF_N(nX_G,nGF) )

    ALLOCATE( Chi     (nE_G,nX_G,nSpecies) )
    ALLOCATE( Sig     (nE_G,nX_G,nSpecies) )
    ALLOCATE( fEQ     (nE_G,nX_G,nSpecies) )
    ALLOCATE( Chi_NES (nE_G,nX_G,nSpecies) )
    ALLOCATE( Eta_NES (nE_G,nX_G,nSpecies) )
    ALLOCATE( Chi_Pair(nE_G,nX_G,nSpecies) )
    ALLOCATE( Eta_Pair(nE_G,nX_G,nSpecies) )

    ALLOCATE( CR_N(nE_G,nX_G,nCR,nSpecies) )
    ALLOCATE( dR_N(nE_G,nX_G,nCR,nSpecies) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S(:) = WFactor_FP * W2_N(:)
    W3_S(:) = WFactor_FP * W3_N(:)

    ALLOCATE( AMAT(n_FP,M_FP,nX_G) )
    ALLOCATE( BVEC(n_FP,nX_G) )
    ALLOCATE( TAU (n_FP,nX_G) )
    ALLOCATE( INFO(nX_G) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX ) &
    !$OMP MAP( alloc: CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, AMAT, BVEC, TAU, INFO, &
    !$OMP             Chi, Sig, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX ) &
    !$ACC CREATE( CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, AMAT, BVEC, TAU, INFO, &
    !$ACC         Chi, Sig, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#endif

    CALL LinearLeastSquares_LWORK( 'N', n_FP, M_FP-1, 1, AMAT, n_FP, BVEC, n_FP, TMP, LWORK )
    ALLOCATE( WORK(LWORK,nX_G) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: WORK )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( WORK )
#endif

  END SUBROUTINE InitializeCollisions_New


  SUBROUTINE FinalizeCollisions_New

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$OMP               CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, AMAT, BVEC, TAU, WORK, INFO, &
    !$OMP               Chi, Sig, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$ACC         CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, AMAT, BVEC, TAU, WORK, INFO, &
    !$ACC         Chi, Sig, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( CF_N, PF_N, AF_N, GX_N, dF_N )
    DEALLOCATE( Chi, Sig, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
    DEALLOCATE( CR_N, dR_N )
    DEALLOCATE( AMAT, BVEC, TAU, WORK, INFO )

  END SUBROUTINE FinalizeCollisions_New


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


  SUBROUTINE ComputePrimitive_Euler_Scalar &
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

  END SUBROUTINE ComputePrimitive_Euler_Scalar

  SUBROUTINE ComputePrimitive_Euler_Vector &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N(:), S_1(:), S_2(:), S_3(:), G(:), Ne(:)
    REAL(DP), INTENT(out) :: D(:), V_1(:), V_2(:), V_3(:), E(:), De(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

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

  END SUBROUTINE ComputePrimitive_Euler_Vector


  SUBROUTINE ComputeConserved_Euler_Scalar &
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

  END SUBROUTINE ComputeConserved_Euler_Scalar

  SUBROUTINE ComputeConserved_Euler_Vector &
    ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D(:), V_1(:), V_2(:), V_3(:), E(:), De(:)
    REAL(DP), INTENT(out) :: N(:), S_1(:), S_2(:), S_3(:), G(:), Ne(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

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

  END SUBROUTINE ComputeConserved_Euler_Vector


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
