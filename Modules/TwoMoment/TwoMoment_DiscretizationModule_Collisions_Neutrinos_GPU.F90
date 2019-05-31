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
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Points

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

  INTEGER  :: nE_G, nX_G, nZ(4), nX(3)
  INTEGER  :: iE_B0,    iE_E0
  INTEGER  :: iE_B1,    iE_E1
  INTEGER  :: iX_B0(3), iX_E0(3)
  INTEGER  :: iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: dF_N(:,:)
  REAL(DP), ALLOCATABLE :: Chi(:,:,:)
  REAL(DP), ALLOCATABLE :: Sig(:,:,:)
  REAL(DP), ALLOCATABLE :: fEQ(:,:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dR_N(:,:,:,:)

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
    REAL(DP) :: Kappa

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Sig )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Sig(iN_E,iN_X,iS) = Zero

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

    CALL SolveMatterEquations_EmAb &
           ( CR_N    (:,:,iCR_N,iNuE    ), &
             CR_N    (:,:,iCR_N,iNuE_Bar), &
             dt * Chi(:,:,      iNuE    ), &
             dt * Chi(:,:,      iNuE_Bar), &
             fEQ     (:,:,      iNuE    ), &
             fEQ     (:,:,      iNuE_Bar), &
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
    !$OMP PRIVATE( Kappa )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Kappa ) &
    !$ACC PRESENT( Chi, Sig, fEQ, CR_N, dR_N )
#elif defined(THORNADO_OMP)
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Kappa = Chi(iN_E,iN_X,iS) + Sig(iN_E,iN_X,iS)

          ! --- Number Density ---

          CR_N(iN_E,iN_X,iCR_N,iS) &
            = ( dt * Chi(iN_E,iN_X,iS) * fEQ(iN_E,iN_X,iS) &
                + CR_N(iN_E,iN_X,iCR_N,iS) ) / ( One + dt * Chi(iN_E,iN_X,iS) )

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
            = Chi(iN_E,iN_X,iS) * ( fEQ(iN_E,iN_X,iS) - CR_N(iN_E,iN_X,iCR_N,iS) )

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
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_F)', MAXLOC(dU_F)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_F)', MAXVAL(dU_F)
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_R)', MAXLOC(dU_R)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_R)', MAXVAL(dU_R)
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

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nX_G = nDOFX * PRODUCT( nX )
    nE_G = nNodesZ(1) * nZ(1)

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )
    ALLOCATE( GX_N(nX_G,nGF) )
    ALLOCATE( dF_N(nX_G,nGF) )

    ALLOCATE( Chi(nE_G,nX_G,nSpecies) )
    ALLOCATE( Sig(nE_G,nX_G,nSpecies) )
    ALLOCATE( fEQ(nE_G,nX_G,nSpecies) )

    ALLOCATE( CR_N(nE_G,nX_G,nCR,nSpecies) )
    ALLOCATE( dR_N(nE_G,nX_G,nCR,nSpecies) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX ) &
    !$OMP MAP( alloc: CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, Chi, Sig, fEQ )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX ) &
    !$ACC CREATE( CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, Chi, Sig, fEQ )
#endif

  END SUBROUTINE InitializeCollisions_New


  SUBROUTINE FinalizeCollisions_New

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$OMP               CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, Chi, Sig, fEQ )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$ACC         CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, Chi, Sig, fEQ )
#endif

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( CF_N, PF_N, AF_N, GX_N, dF_N )
    DEALLOCATE( Chi, Sig, fEQ )
    DEALLOCATE( CR_N, dR_N )

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
