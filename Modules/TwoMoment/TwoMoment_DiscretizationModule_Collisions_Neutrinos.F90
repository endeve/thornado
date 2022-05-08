#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, FourPi
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
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
    ComputeThermodynamicStates_Primitive_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    InitializeNeutrinoMatterSolver, &
    FinalizeNeutrinoMatterSolver, &
    SolveMatterEquations_EmAb_NuE, &
    SolveMatterEquations_EmAb_FP, &
    SolveMatterEquations_FP_Coupled, &
    SolveMatterEquations_FP_NestedAA, &
    SolveMatterEquations_FP_NestedNewton, &
    E_N

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_New

  PUBLIC :: InitializeNonlinearSolverTally
  PUBLIC :: FinalizeNonlinearSolverTally
  PUBLIC :: WriteNonlinearSolverTally

  LOGICAL, PARAMETER :: ReportConvergenceData = .FALSE.
  INTEGER  :: Iterations_Min
  INTEGER  :: Iterations_Max
  REAL(DP) :: Iterations_Ave

  ! --- Solver Tally ---

  LOGICAL              :: TallyNonlinearSolver = .FALSE.
  INTEGER              :: TallyFileNumber
  INTEGER, ALLOCATABLE :: MinIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: MaxIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: AveIterations_K(:,:,:)
  INTEGER, ALLOCATABLE :: AveIterationsInner_K(:,:,:)

  INTEGER, ALLOCATABLE :: nIterations(:)
  INTEGER, ALLOCATABLE :: nIterations_Inner(:)
  INTEGER, ALLOCATABLE :: nIterations_Outer(:)

  INTEGER  :: nE_G, nX_G, nZ(4), nX(3)
  INTEGER  :: iE_B0,    iE_E0
  INTEGER  :: iE_B1,    iE_E1
  INTEGER  :: iX_B0(3), iX_E0(3)
  INTEGER  :: iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: dF_N(:,:)
  REAL(DP), ALLOCATABLE :: Chi(:,:,:)
  REAL(DP), ALLOCATABLE :: Sig_0(:,:,:)
  REAL(DP), ALLOCATABLE :: Sig_1(:,:,:)
  REAL(DP), ALLOCATABLE :: fEQ(:,:,:)
  REAL(DP), ALLOCATABLE :: Chi_NES(:,:,:)
  REAL(DP), ALLOCATABLE :: Eta_NES(:,:,:)
  REAL(DP), ALLOCATABLE :: Chi_Pair(:,:,:)
  REAL(DP), ALLOCATABLE :: Eta_Pair(:,:,:)
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
      dU_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iX1, iX2, iX3, iGF, iCF, iCR, iS, iE, iN_E, iN_X
    INTEGER  :: iNode, iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: Chi_T, Eta_T, Eta, Kappa

    CALL TimersStart( Timer_Implicit )

    CALL TimersStart( Timer_Im_In )

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

    CALL TimersStop( Timer_Im_In )

    ! --- Copy inputs to local arrays ---

    CALL TimersStart( Timer_Im_MapForward )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_F, iX_B1, iX_E1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
      DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
          DO iX1 = iX_B1(1), iX_E1(1)
            DO iNodeX = 1, nDOFX
              dU_F(iNodeX,iX1,iX2,iX3,iCF) = Zero
            END DO
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B1(4), iZ_E1(4)
          DO iZ3 = iZ_B1(3), iZ_E1(3)
            DO iZ2 = iZ_B1(2), iZ_E1(2)
              DO iZ1 = iZ_B1(1), iZ_E1(1)
                DO iNode = 1, nDOF
                  dU_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                END DO
              END DO
            END DO
          END DO
        END DO
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
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, dF_N, U_F )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
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

    ! --- Rearrange data to group energy nodes together ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNode, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNode, iNodeX, iNodeE, iE, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iNode, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
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

            CR_N(iN_E,iN_X,iS,iCR) = U_R(iNode,iE,iX1,iX2,iX3,iCR,iS)

          END DO
        END DO
      END DO
    END DO

    ! --- Compute Primitive Quantities ---

    CALL ComputePrimitive_Euler &
           ( CF_N(:,iCF_D ), &
             CF_N(:,iCF_S1), &
             CF_N(:,iCF_S2), &
             CF_N(:,iCF_S3), &
             CF_N(:,iCF_E ), &
             CF_N(:,iCF_Ne), &
             PF_N(:,iPF_D ), &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             PF_N(:,iPF_E ), &
             PF_N(:,iPF_Ne), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33) )

    CALL TimersStop( Timer_Im_MapForward )

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_EosIn )

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D ), &
             PF_N(:,iPF_E ), &
             PF_N(:,iPF_Ne), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_E ), &
             AF_N(:,iAF_Ye) )

    CALL TimersStop( Timer_Im_EosIn )

    CALL TimersStart( Timer_Im_Solve )

    ! --- Opacity Table Lookup ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    !DO iS = 1, nSpecies

    !  CALL ComputeNeutrinoOpacities_EC_Points &
    !         ( 1, nE_G, 1, nX_G, &
    !           E_N (:), &
    !           PF_N(:,iPF_D ), &
    !           AF_N(:,iAF_T ), &
    !           AF_N(:,iAF_Ye), &
    !           iS, Chi(:,:,iS) )

    !END DO

    !DO iS = 1, nSpecies

    !  ! iMoment = 1
    !  CALL ComputeNeutrinoOpacities_ES_Points &
    !         ( 1, nE_G, 1, nX_G, &
    !           E_N (:), &
    !           PF_N(:,iPF_D ), &
    !           AF_N(:,iAF_T ), &
    !           AF_N(:,iAF_Ye), &
    !           iS, 1, Sig_0(:,:,iS) )
    !  ! iMoment = 2
    !  CALL ComputeNeutrinoOpacities_ES_Points &
    !         ( 1, nE_G, 1, nX_G, &
    !           E_N (:), &
    !           PF_N(:,iPF_D ), &
    !           AF_N(:,iAF_T ), &
    !           AF_N(:,iAF_Ye), &
    !           iS, 2, Sig_1(:,:,iS) )

    !END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, Sig_0, Sig_1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Chi_NES (iN_E,iN_X,iS) = Zero
          Eta_NES (iN_E,iN_X,iS) = Zero
          Chi_Pair(iN_E,iN_X,iS) = Zero
          Eta_Pair(iN_E,iN_X,iS) = Zero

          Sig_0   (iN_E,iN_X,iS) = FourPi * E_N(iN_E)**2 * Sig_0(iN_E,iN_X,iS)
          Sig_1   (iN_E,iN_X,iS) = FourPi * E_N(iN_E)**2 * Sig_1(iN_E,iN_X,iS)

        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_ComputeOpacity )

    IF( nSpecies .EQ. 1 )THEN

      ! --- Single Species (Electron Neutrinos) ---

      CALL SolveMatterEquations_EmAb_NuE &
             ( dt, &
               CR_N(:,:,iNuE,iCR_N), &
               Chi (:,:,iNuE), &
               fEQ (:,:,iNuE), &
               PF_N(:,iPF_D ), &
               AF_N(:,iAF_T ), &
               AF_N(:,iAF_Ye), &
               AF_N(:,iAF_E ), &
               nIterations(:) )

    ELSE

      ! --- Electron Neutrinos and Antineutrinos ---

#if defined(NEUTRINO_MATTER_SOLVER_EMAB)

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iNuE, iNuE_Bar, &
               CR_N    (:,:,iNuE    ,iCR_N), &
               CR_N    (:,:,iNuE_Bar,iCR_N), &
               Chi     (:,:,iNuE    ), &
               Chi     (:,:,iNuE_Bar), &
               fEQ     (:,:,iNuE    ), &
               fEQ     (:,:,iNuE_Bar), &
               PF_N(:,iPF_D ), &
               AF_N(:,iAF_T ), &
               AF_N(:,iAF_Ye), &
               AF_N(:,iAF_E ), &
               nIterations(:), &
               TOL = 1.0d-7 )

      CALL TimersStop( Timer_Im_EmAb_FP )

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_COUPLED)

      CALL TimersStart( Timer_Im_CoupledAA )

      CALL SolveMatterEquations_FP_Coupled &
             ( dt, iNuE, iNuE_Bar, &
               CR_N    (:,:,iNuE    ,iCR_N), &
               CR_N    (:,:,iNuE_Bar,iCR_N), &
               Chi     (:,:,iNuE    ), &
               Chi     (:,:,iNuE_Bar), &
               fEQ     (:,:,iNuE    ), &
               fEQ     (:,:,iNuE_Bar), &
               Chi_NES (:,:,iNuE    ), &
               Chi_NES (:,:,iNuE_Bar), &
               Eta_NES (:,:,iNuE    ), &
               Eta_NES (:,:,iNuE_Bar), &
               Chi_Pair(:,:,iNuE    ), &
               Chi_Pair(:,:,iNuE_Bar), &
               Eta_Pair(:,:,iNuE    ), &
               Eta_Pair(:,:,iNuE_Bar), &
               PF_N    (:,iPF_D ), &
               AF_N    (:,iAF_T ), &
               AF_N    (:,iAF_Ye), &
               AF_N    (:,iAF_E ), &
               nIterations(:) )

      CALL TimersStop( Timer_Im_CoupledAA )

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_AA)

      CALL TimersStart( Timer_Im_NestedAA )

      CALL SolveMatterEquations_FP_NestedAA &
             ( dt, iNuE, iNuE_Bar, &
               CR_N    (:,:,iNuE    ,iCR_N), &
               CR_N    (:,:,iNuE_Bar,iCR_N), &
               Chi     (:,:,iNuE    ), &
               Chi     (:,:,iNuE_Bar), &
               fEQ     (:,:,iNuE    ), &
               fEQ     (:,:,iNuE_Bar), &
               Chi_NES (:,:,iNuE    ), &
               Chi_NES (:,:,iNuE_Bar), &
               Eta_NES (:,:,iNuE    ), &
               Eta_NES (:,:,iNuE_Bar), &
               Chi_Pair(:,:,iNuE    ), &
               Chi_Pair(:,:,iNuE_Bar), &
               Eta_Pair(:,:,iNuE    ), &
               Eta_Pair(:,:,iNuE_Bar), &
               PF_N    (:,iPF_D ), &
               AF_N    (:,iAF_T ), &
               AF_N    (:,iAF_Ye), &
               AF_N    (:,iAF_E ), &
               nIterations_Inner(:), &
               nIterations_Outer(:) )

      CALL TimersStop( Timer_Im_NestedAA )

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON)

      CALL TimersStart( Timer_Im_NestedNewton )

      CALL SolveMatterEquations_FP_NestedNewton &
             ( dt, iNuE, iNuE_Bar, &
               CR_N    (:,:,iNuE    ,iCR_N), &
               CR_N    (:,:,iNuE_Bar,iCR_N), &
               Chi     (:,:,iNuE    ), &
               Chi     (:,:,iNuE_Bar), &
               fEQ     (:,:,iNuE    ), &
               fEQ     (:,:,iNuE_Bar), &
               Chi_NES (:,:,iNuE    ), &
               Chi_NES (:,:,iNuE_Bar), &
               Eta_NES (:,:,iNuE    ), &
               Eta_NES (:,:,iNuE_Bar), &
               Chi_Pair(:,:,iNuE    ), &
               Chi_Pair(:,:,iNuE_Bar), &
               Eta_Pair(:,:,iNuE    ), &
               Eta_Pair(:,:,iNuE_Bar), &
               PF_N    (:,iPF_D ), &
               AF_N    (:,iAF_T ), &
               AF_N    (:,iAF_Ye), &
               AF_N    (:,iAF_E ), &
               nIterations_Inner(:), &
               nIterations_Outer(:) )

      CALL TimersStop( Timer_Im_NestedNewton )

#elif defined(NEUTRINO_MATTER_SOLVER_NEWTON)

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

    END IF

    CALL TimersStop( Timer_Im_Solve )

    CALL TimersStart( Timer_Im_Increment )

    ! --- Update Radiation Fields ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Chi_T, Eta_T, Eta, Kappa )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Chi_T, Eta_T, Eta, Kappa ) &
    !$ACC PRESENT( Chi, Chi_NES, Chi_Pair, Eta_NES, Eta_Pair, &
    !$ACC          Sig_0, Sig_1, fEQ, CR_N, dR_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Chi_T, Eta_T, Eta, Kappa )
#endif
    DO iS = 1, nSpecies
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G

          Chi_T = Chi(iN_E,iN_X,iS) + Chi_NES(iN_E,iN_X,iS) + Chi_Pair(iN_E,iN_X,iS)
          Kappa = Chi_T + Sig_0(iN_E,iN_X,iS) - Sig_1(iN_E, iN_X,iS) / 3.0_DP

          Eta   = Chi(iN_E,iN_X,iS) * fEQ(iN_E,iN_X,iS)
          Eta_T = Eta + Eta_NES(iN_E,iN_X,iS) + Eta_Pair(iN_E,iN_X,iS)

          ! --- Number Density Updated in Nonlinear Solvers ---

          ! --- Number Flux (1) ---

          CR_N(iN_E,iN_X,iS,iCR_G1) &
            = CR_N(iN_E,iN_X,iS,iCR_G1) / ( One + dt * Kappa )

          ! --- Number Flux (2) ---

          CR_N(iN_E,iN_X,iS,iCR_G2) &
            = CR_N(iN_E,iN_X,iS,iCR_G2) / ( One + dt * Kappa )

          ! --- Number Flux (3) ---

          CR_N(iN_E,iN_X,iS,iCR_G3) &
            = CR_N(iN_E,iN_X,iS,iCR_G3) / ( One + dt * Kappa )

          ! --- Increments ---

          dR_N(iN_E,iN_X,iCR_N,iS) &
            = Eta_T - Chi_T * CR_N(iN_E,iN_X,iS,iCR_N)

          dR_N(iN_E,iN_X,iCR_G1,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iS,iCR_G1)

          dR_N(iN_E,iN_X,iCR_G2,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iS,iCR_G2)

          dR_N(iN_E,iN_X,iCR_G3,iS) &
            = - Kappa * CR_N(iN_E,iN_X,iS,iCR_G3)

        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Im_Increment )

    ! --- EOS Table Lookup ---

    CALL TimersStart( Timer_Im_EosOut )

    CALL ComputeThermodynamicStates_Primitive_TABLE &
           ( PF_N(:,iPF_D ), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_Ye), &
             PF_N(:,iPF_E ), &
             AF_N(:,iAF_E ), &
             PF_N(:,iPF_Ne) )

    CALL TimersStop( Timer_Im_EosOut )

    CALL TimersStart( Timer_Im_MapBackward )

    ! --- Compute Conserved Quantities ---

    CALL ComputeConserved_Euler &
           ( PF_N(:,iPF_D ), &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             PF_N(:,iPF_E ), &
             PF_N(:,iPF_Ne), &
             CF_N(:,iCF_D ), &
             CF_N(:,iCF_S1), &
             CF_N(:,iCF_S2), &
             CF_N(:,iCF_S3), &
             CF_N(:,iCF_E ), &
             CF_N(:,iCF_Ne), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33) )

    ! --- Rearrange data back to original layout ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( nX, iX_B0, iX_E0, dU_F, CF_N, dF_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, dU_R, dR_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
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

#ifdef THORNADO_DEBUG_IMPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_F, dU_R )
#endif
    WRITE(*,'(a,8x,5i4,es23.15)') 'MINLOC(dU_F), MINVAL(dU_F)', MINLOC(dU_F), MINVAL(dU_F)
    WRITE(*,'(a,8x,5i4,es23.15)') 'MAXLOC(dU_F), MAXVAL(dU_F)', MAXLOC(dU_F), MAXVAL(dU_F)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(dU_R), MINVAL(dU_R)', MINLOC(dU_R), MINVAL(dU_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(dU_R), MAXVAL(dU_R)', MAXLOC(dU_R), MAXVAL(dU_R)
#endif

    CALL TimersStart( Timer_Im_Out )

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

    CALL TimersStop( Timer_Im_Out )

    CALL TimersStop( Timer_Implicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_New


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

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )
    ALLOCATE( GX_N(nX_G,nGF) )
    ALLOCATE( dF_N(nX_G,nGF) )

    ALLOCATE( Chi     (nE_G,nX_G,nSpecies) )
    ALLOCATE( Sig_0     (nE_G,nX_G,nSpecies) )
    ALLOCATE( Sig_1     (nE_G,nX_G,nSpecies) )
    ALLOCATE( fEQ     (nE_G,nX_G,nSpecies) )
    ALLOCATE( Chi_NES (nE_G,nX_G,nSpecies) )
    ALLOCATE( Eta_NES (nE_G,nX_G,nSpecies) )
    ALLOCATE( Chi_Pair(nE_G,nX_G,nSpecies) )
    ALLOCATE( Eta_Pair(nE_G,nX_G,nSpecies) )

    ALLOCATE( CR_N(nE_G,nX_G,nSpecies,nCR) )
    ALLOCATE( dR_N(nE_G,nX_G,nCR,nSpecies) )

    ALLOCATE( nIterations(nX_G) )
    ALLOCATE( nIterations_Inner(nX_G) )
    ALLOCATE( nIterations_Outer(nX_G) )

    nIterations(:) = 0
    nIterations_Inner(:) = 0
    nIterations_Outer(:) = 0

    CALL InitializeNeutrinoMatterSolver( iZ_B0, iZ_E0 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$OMP          nIterations, nIterations_Inner, nIterations_Outer ) &
    !$OMP MAP( alloc: CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, &
    !$OMP             Chi, Sig_0, Sig_1, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$ACC         nIterations, nIterations_Inner, nIterations_Outer ) &
    !$ACC CREATE( CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, &
    !$ACC         Chi, Sig_0, Sig_1, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#endif

  END SUBROUTINE InitializeCollisions_New


  SUBROUTINE FinalizeCollisions_New

    INTEGER :: iX1, iX2, iX3, iN_X, iNodeX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: nIterations, nIterations_Inner, nIterations_Outer ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$OMP               CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, &
    !$OMP               Chi, Sig_0, Sig_1, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( nIterations, nIterations_Inner, nIterations_Outer ) &
    !$ACC DELETE( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, &
    !$ACC         CF_N, PF_N, AF_N, GX_N, dF_N, CR_N, dR_N, &
    !$ACC         Chi, Sig_0, Sig_1, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
#endif

    IF( TallyNonlinearSolver )THEN

      MinIterations_K      = + HUGE( 1 )
      MaxIterations_K      = - HUGE( 1 )
      AveIterations_K      = 0
      AveIterationsInner_K = 0

      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iNodeX = 1, nDOFX

              iN_X = iNodeX &
                     + ( iX1 - iX_B0(1) ) * nDOFX &
                     + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
                     + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

#if defined(NEUTRINO_MATTER_SOLVER_EMAB)

              MinIterations_K(iX1,iX2,iX2) &
                = MIN( nIterations(iN_X), MinIterations_K(iX1,iX2,iX2) )
              MaxIterations_K(iX1,iX2,iX3) &
                = MAX( nIterations(iN_X), MaxIterations_K(iX1,iX2,iX3) )
              AveIterations_K(iX1,iX2,iX3) &
                = AveIterations_K(iX1,iX2,iX3) + nIterations(iN_X)

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_COUPLED)

              MinIterations_K(iX1,iX2,iX2) &
                = MIN( nIterations(iN_X), MinIterations_K(iX1,iX2,iX2) )
              MaxIterations_K(iX1,iX2,iX3) &
                = MAX( nIterations(iN_X), MaxIterations_K(iX1,iX2,iX3) )
              AveIterations_K(iX1,iX2,iX3) &
                = AveIterations_K(iX1,iX2,iX3) + nIterations(iN_X)

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_AA)

              MinIterations_K(iX1,iX2,iX2) &
                = MIN( nIterations_Outer(iN_X), MinIterations_K(iX1,iX2,iX2) )
              MaxIterations_K(iX1,iX2,iX3) &
                = MAX( nIterations_Outer(iN_X), MaxIterations_K(iX1,iX2,iX3) )
              AveIterations_K(iX1,iX2,iX3) &
                = AveIterations_K(iX1,iX2,iX3)      + nIterations_Outer(iN_X)
              AveIterationsInner_K(iX1,iX2,iX3) &
                = AveIterationsInner_K(iX1,iX2,iX3) + nIterations_Inner(iN_X)

#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON)

              MinIterations_K(iX1,iX2,iX2) &
                = MIN( nIterations_Outer(iN_X), MinIterations_K(iX1,iX2,iX2) )
              MaxIterations_K(iX1,iX2,iX3) &
                = MAX( nIterations_Outer(iN_X), MaxIterations_K(iX1,iX2,iX3) )
              AveIterations_K(iX1,iX2,iX3) &
                = AveIterations_K(iX1,iX2,iX3)      + nIterations_Outer(iN_X)
              AveIterationsInner_K(iX1,iX2,iX3) &
                = AveIterationsInner_K(iX1,iX2,iX3) + nIterations_Inner(iN_X)

#elif defined(NEUTRINO_MATTER_SOLVER_NEWTON)

              MinIterations_K(iX1,iX2,iX2) &
                = MIN( nIterations(iN_X), MinIterations_K(iX1,iX2,iX2) )
              MaxIterations_K(iX1,iX2,iX3) &
                = MAX( nIterations(iN_X), MaxIterations_K(iX1,iX2,iX3) )
              AveIterations_K(iX1,iX2,iX3) &
                = AveIterations_K(iX1,iX2,iX3) + nIterations(iN_X)

#endif

            END DO

            AveIterations_K(iX1,iX2,iX3) &
              = FLOOR( DBLE(AveIterations_K(iX1,iX2,iX3)     )/DBLE(nDOFX) )
            AveIterationsInner_K(iX1,iX2,iX3) &
              = FLOOR( DBLE(AveIterationsInner_K(iX1,iX2,iX3))/DBLE(nDOFX) )

          END DO
        END DO
      END DO

    END IF

    IF( ReportConvergenceData )THEN

#if defined(NEUTRINO_MATTER_SOLVER_EMAB)
      Iterations_Min = MINVAL( nIterations(:) )
      Iterations_Max = MAXVAL( nIterations(:) )
      Iterations_Ave = DBLE( SUM( nIterations(:) ) ) / DBLE( nX_G )
#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_COUPLED)
      Iterations_Min = MINVAL( nIterations(:) )
      Iterations_Max = MAXVAL( nIterations(:) )
      Iterations_Ave = DBLE( SUM( nIterations(:) ) ) / DBLE( nX_G )
#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_AA)
      Iterations_Min = MINVAL( nIterations_Outer(:) )
      Iterations_Max = MAXVAL( nIterations_Outer(:) )
      Iterations_Ave = DBLE( SUM( nIterations_Outer(:) ) ) / DBLE( nX_G )
#elif defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON)
      Iterations_Min = MINVAL( nIterations_Outer(:) )
      Iterations_Max = MAXVAL( nIterations_Outer(:) )
      Iterations_Ave = DBLE( SUM( nIterations_Outer(:) ) ) / DBLE( nX_G )
#elif defined(NEUTRINO_MATTER_SOLVER_NEWTON)
      Iterations_Min = MINVAL( nIterations(:) )
      Iterations_Max = MAXVAL( nIterations(:) )
      Iterations_Ave = DBLE( SUM( nIterations(:) ) ) / DBLE( nX_G )
#endif

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', Iterations_Min
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', Iterations_Max
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', Iterations_Ave
      WRITE(*,*)

    END IF

    DEALLOCATE( CF_N, PF_N, AF_N, GX_N, dF_N )
    DEALLOCATE( Chi, Sig_0, Sig_1, fEQ, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair )
    DEALLOCATE( CR_N, dR_N )
    DEALLOCATE( nIterations, nIterations_Inner, nIterations_Outer )

    CALL FinalizeNeutrinoMatterSolver

  END SUBROUTINE FinalizeCollisions_New


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
    E   = G - Half * ( + S_1**2 / Gm_dd_11 &
                       + S_2**2 / Gm_dd_22 &
                       + S_3**2 / Gm_dd_33 ) / N
    De  = Ne

  END SUBROUTINE ComputePrimitive_Euler_Scalar


  SUBROUTINE ComputePrimitive_Euler_Vector &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: N(:), S_1(:), S_2(:), S_3(:), G(:), Ne(:)
    REAL(DP), INTENT(out) :: D(:), V_1(:), V_2(:), V_3(:), E(:), De(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

    INTEGER :: iN_X, nX_P

    nX_P = SIZE(N)

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( N, S_1, S_2, S_3, G, Ne, &
    !$ACC          D, V_1, V_2, V_3, E, De, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_P

      D  (iN_X) = N(iN_X)
      V_1(iN_X) = S_1(iN_X) / ( Gm_dd_11(iN_X) * N(iN_X) )
      V_2(iN_X) = S_2(iN_X) / ( Gm_dd_22(iN_X) * N(iN_X) )
      V_3(iN_X) = S_3(iN_X) / ( Gm_dd_33(iN_X) * N(iN_X) )
      E  (iN_X) = G(iN_X) - Half * ( + S_1(iN_X)**2 / Gm_dd_11(iN_X) &
                                     + S_2(iN_X)**2 / Gm_dd_22(iN_X) &
                                     + S_3(iN_X)**2 / Gm_dd_33(iN_X) ) / N(iN_X)
      De (iN_X) = Ne(iN_X)

    END DO

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

    REAL(DP), INTENT(in)  :: D(:), V_1(:), V_2(:), V_3(:), E(:), De(:)
    REAL(DP), INTENT(out) :: N(:), S_1(:), S_2(:), S_3(:), G(:), Ne(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

    INTEGER :: iN_X, nX_P

    nX_P = SIZE(N)

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( N, S_1, S_2, S_3, G, Ne, &
    !$ACC          D, V_1, V_2, V_3, E, De, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_P

      N  (iN_X) = D(iN_X)
      S_1(iN_X) = D(iN_X) * Gm_dd_11(iN_X) * V_1(iN_X)
      S_2(iN_X) = D(iN_X) * Gm_dd_22(iN_X) * V_2(iN_X)
      S_3(iN_X) = D(iN_X) * Gm_dd_33(iN_X) * V_3(iN_X)
      G  (iN_X) = E(iN_X) + Half * D(iN_X) * ( + Gm_dd_11(iN_X) * V_1(iN_X)**2 &
                                               + Gm_dd_22(iN_X) * V_2(iN_X)**2 &
                                               + Gm_dd_33(iN_X) * V_3(iN_X)**2 )
      Ne (iN_X) = De(iN_X)

    END DO

  END SUBROUTINE ComputeConserved_Euler_Vector


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
