MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE KindModule, ONLY: &
    DP, &
    Zero
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    AtomicMassUnit
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions, &
    Timer_Collisions_PrimitiveFluid, &
    Timer_Collisions_PrimitiveTwoMoment, &
    Timer_Collisions_Solve
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_T, iAF_E , iAF_Ye, iAF_P
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nDR, iDR_iter_outer, iDR_iter_inner
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    SolveNeutrinoMatterCoupling_FP_Nested_AA, &
    InitializeNeutrinoMatterSolver, &
    FinalizeNeutrinoMatterSolver, &
    InitializeNeutrinoMatterSolverParameters
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, &
    ComputeConserved_TwoMoment
  USE OpacityModule_TABLE, ONLY: &
    QueryOpacity
#if   defined( TWOMOMENT_ORDER_1 )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeConserved_Euler_NonRelativistic

#elif defined( TWOMOMENT_ORDER_V )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeConserved_Euler_NonRelativistic

#elif defined( TWOMOMENT_RELATIVISTIC )

  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputePressureFromPrimitive_TABLE
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    ComputeConserved_Euler_Relativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER               :: nE_G, nX_G
  INTEGER               :: nZ(4), nX(3), nE
  INTEGER               :: iE_B0, iE_E0, iX_B0(3), iX_E0(3)
  INTEGER               :: iE_B1, iE_E1, iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: GE_N(:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)

  REAL(DP), ALLOCATABLE, TARGET :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: PR_N(:,:,:,:)

  REAL(DP), DIMENSION(:), CONTIGUOUS, POINTER :: N_P, G1_P, G2_P, G3_P
  REAL(DP), DIMENSION(:), CONTIGUOUS, POINTER :: J_P, H1_P, H2_P, H3_P

  INTEGER,  DIMENSION(:), ALLOCATABLE :: PositionIndexZ

  INTEGER, ALLOCATABLE :: nIterations_Inner(:)
  INTEGER, ALLOCATABLE :: nIterations_Outer(:)
  INTEGER, ALLOCATABLE :: nIterations_Prim(:)

  LOGICAL, PARAMETER :: ReportConvergenceData = .FALSE.

  REAL(DP), PARAMETER :: UnitD    = Gram / Centimeter**3

CONTAINS


  ! --- Public Subroutines ---


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R, &
      uDR_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      dt
    REAL(DP), INTENT(in) :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in) :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option &
          (iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nDR)

    INTEGER :: iN_X, iN_E, iS
    INTEGER :: iNodeX, iX1, iX2, iX3

    CALL TimersStart( Timer_Collisions )

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: dU_F, dU_R )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( dU_F, dU_R )
#endif

    CALL MapDataForCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    ! --- Compute Primitive Fluid ---

    CALL TimersStart( Timer_Collisions_PrimitiveFluid )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif



#if   defined( TWOMOMENT_ORDER_1 )

    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
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

#elif defined( TWOMOMENT_ORDER_V )

    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
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

#elif defined( TWOMOMENT_RELATIVISTIC )


    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_Relativistic &
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

#endif


    CALL TimersStop( Timer_Collisions_PrimitiveFluid )

    ! --- Compute Primitive Moments ---

    CALL TimersStart( Timer_Collisions_PrimitiveTwoMoment )

#if   defined( TWOMOMENT_ORDER_1 )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_ORDER_V )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             PositionIndexZ, &
             nIterations_Prim )

#elif defined( TWOMOMENT_RELATIVISTIC )

    CALL ComputePrimitive_TwoMoment &
           ( N_P, G1_P, G2_P, G3_P, &
             J_P, H1_P, H2_P, H3_P, &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             GX_N(:,iGF_Alpha   ), &
             GX_N(:,iGF_Beta_1  ), &
             GX_N(:,iGF_Beta_2  ), &
             GX_N(:,iGF_Beta_3  ), &
             PositionIndexZ, &
             nIterations_Prim )

#endif

    CALL TimersStop( Timer_Collisions_PrimitiveTwoMoment )

    ! --- EOS Table Lookup ---

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye) )

    ! --- Neutrino-Matter Coupling Solve ---

    CALL TimersStart( Timer_Collisions_Solve )

    CALL SolveNeutrinoMatterCoupling_FP_Nested_AA &
           ( dt, &
             PR_N(:,:,:,iCR_N ), &
             PR_N(:,:,:,iCR_G1), &
             PR_N(:,:,:,iCR_G2), &
             PR_N(:,:,:,iCR_G3), &
             PF_N(:,iPF_V1), &
             PF_N(:,iPF_V2), &
             PF_N(:,iPF_V3), &
             PF_N(:,iPF_D ), &
             AF_N(:,iAF_T ), &
             AF_N(:,iAF_Ye), &
             AF_N(:,iAF_E ), &
             GX_N(:,iGF_Gm_dd_11), &
             GX_N(:,iGF_Gm_dd_22), &
             GX_N(:,iGF_Gm_dd_33), &
             GX_N(:,iGF_Alpha   ), &
             GX_N(:,iGF_Beta_1  ), &
             GX_N(:,iGF_Beta_2  ), &
             GX_N(:,iGF_Beta_3  ), &
             CR_N(:,:,:,iCR_N ), &
             CR_N(:,:,:,iCR_G1), &
             CR_N(:,:,:,iCR_G2), &
             CR_N(:,:,:,iCR_G3), &
             nIterations_Inner, &
             nIterations_Outer )

    CALL TimersStop( Timer_Collisions_Solve )

    ! --- Store Iteration Counts (If Requested) ---

    IF( PRESENT( uDR_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( uDR_Option, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        uDR_Option(iX1,iX2,iX3,iDR_iter_outer) = 0.0_DP
        uDR_Option(iX1,iX2,iX3,iDR_iter_inner) = 0.0_DP

      END DO
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( uDR_Option, iX_B0, iX_E0, nX, &
      !$ACC          nIterations_Outer, nIterations_Inner )
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(4) &
      !$OMP PRIVATE( iN_X )
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        IF ( nIterations_Outer(iN_X) == 0 ) THEN
          nIterations_Inner(iN_X) = 0
        ELSE
          nIterations_Inner(iN_X) &
            = FLOOR(   REAL( nIterations_Inner(iN_X), DP ) &
                     / REAL( nIterations_Outer(iN_X), DP ) )

        ENDIF

        uDR_Option(iX1,iX2,iX3,iDR_iter_outer) &
          = MAX( uDR_Option(iX1,iX2,iX3,iDR_iter_outer), &
                 REAL( nIterations_Outer(iN_X), DP ) )

        uDR_Option(iX1,iX2,iX3,iDR_iter_inner) &
          = MAX( uDR_Option(iX1,iX2,iX3,iDR_iter_inner), &
                 REAL( nIterations_Inner(iN_X), DP ) )

      END DO
      END DO
      END DO
      END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( PR_N, CR_N, PF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

#if   defined( TWOMOMENT_ORDER_1 )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_ORDER_V )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

#elif defined( TWOMOMENT_RELATIVISTIC )

      CALL ComputeConserved_TwoMoment &
             ( PR_N(iN_E,iS,iN_X,iCR_N ), &
               PR_N(iN_E,iS,iN_X,iCR_G1), &
               PR_N(iN_E,iS,iN_X,iCR_G2), &
               PR_N(iN_E,iS,iN_X,iCR_G3), &
               CR_N(iN_E,iS,iN_X,iCR_N ), &
               CR_N(iN_E,iS,iN_X,iCR_G1), &
               CR_N(iN_E,iS,iN_X,iCR_G2), &
               CR_N(iN_E,iS,iN_X,iCR_G3), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33), &
               Zero, Zero, Zero, &
               GX_N(iN_X,iGF_Alpha   ), &
               GX_N(iN_X,iGF_Beta_1  ), &
               GX_N(iN_X,iGF_Beta_2  ), &
               GX_N(iN_X,iGF_Beta_3  ) )

#endif

    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( PF_N, CF_N, GX_N )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif



#if   defined( TWOMOMENT_ORDER_1 )

    DO iN_X = 1, nX_G

      PF_N(iN_X,iPF_E ) = AF_N(iN_X,iAF_E ) * PF_N(iN_X,iPF_D)
      PF_N(iN_X,iPF_Ne) = AF_N(iN_X,iAF_Ye) * PF_N(iN_X,iPF_D) / AtomicMassUnit

      CALL ComputeConserved_Euler_NonRelativistic &
             ( PF_N(iN_X,iPF_D),  &
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


#elif defined( TWOMOMENT_ORDER_V )

    DO iN_X = 1, nX_G

      PF_N(iN_X,iPF_E ) = AF_N(iN_X,iAF_E ) * PF_N(iN_X,iPF_D)
      PF_N(iN_X,iPF_Ne) = AF_N(iN_X,iAF_Ye) * PF_N(iN_X,iPF_D) / AtomicMassUnit

      CALL ComputeConserved_Euler_NonRelativistic &
             ( PF_N(iN_X,iPF_D),  &
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


#elif defined( TWOMOMENT_RELATIVISTIC )


    CALL ComputeThermodynamicStates_Primitive_TABLE &
           ( PF_N(:,iPF_D), AF_N(:,iAF_T), AF_N(:,iAF_Ye), &
             PF_N(:,iPF_E), AF_N(:,iAF_E), PF_N(:,iPF_Ne) )

    CALL ComputePressureFromPrimitive_TABLE & 
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), &
             AF_N(:,iAF_P) )


    DO iN_X = 1, nX_G

      CALL ComputeConserved_Euler_Relativistic &
             ( PF_N(iN_X,iPF_D),  &
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
               GX_N(iN_X,iGF_Gm_dd_33), &
               AF_N(iN_X,iAF_P) )

    END DO


#endif


    CALL ComputeAndMapIncrement &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, U_F, U_R, dU_F, dU_R )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_F, dU_R ) &
    !$OMP MAP( release: GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_F, dU_R ) &
    !$ACC DELETE( GE, GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

    CALL FinalizeCollisions

    CALL TimersStop( Timer_Collisions )

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  ! --- Private Subroutines ---


  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)

    INTEGER :: iN_X, iN_E, iS, iZ, nZ_G

    iE_B0 = iZ_B0(1)  ; iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1)  ; iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1
    nZ = iZ_E0 - iZ_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )
    nZ_G = nE_G * nX_G * nSpecies

    ALLOCATE( GE_N(nE_G,nGE) )
    ALLOCATE( GX_N(nX_G,nGF) )

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )

    ALLOCATE( CR_N(nE_G,nSpecies,nX_G,nCR) )
    ALLOCATE( PR_N(nE_G,nSpecies,nX_G,nCR) )

    ALLOCATE( PositionIndexZ(nZ_G) )

    iZ = 0
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      iZ = iZ + 1

      PositionIndexZ(iZ) = iN_X

    END DO
    END DO
    END DO

    ALLOCATE( nIterations_Inner(nX_G) )
    ALLOCATE( nIterations_Outer(nX_G) )
    ALLOCATE( nIterations_Prim (nZ_G) )

    nIterations_Inner(:) = 0
    nIterations_Outer(:) = 0
    nIterations_Prim (:) = 0

    N_P (1:nZ_G) => CR_N(:,:,:,iCR_N )
    G1_P(1:nZ_G) => CR_N(:,:,:,iCR_G1)
    G2_P(1:nZ_G) => CR_N(:,:,:,iCR_G2)
    G3_P(1:nZ_G) => CR_N(:,:,:,iCR_G3)

    J_P (1:nZ_G) => PR_N(:,:,:,iCR_N )
    H1_P(1:nZ_G) => PR_N(:,:,:,iCR_G1)
    H2_P(1:nZ_G) => PR_N(:,:,:,iCR_G2)
    H3_P(1:nZ_G) => PR_N(:,:,:,iCR_G3)

    ! --- Neutrino-Matter Solver Parameter Initialization ---
    ! --- ( can be moved to the program init ) --------------

    CALL InitializeNeutrinoMatterSolverParameters()

    CALL InitializeNeutrinoMatterSolver( iZ_B0, iZ_E0 )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$OMP          nIterations_Inner, nIterations_Outer, nIterations_Prim ) &
    !$OMP MAP( alloc: GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$ACC         nIterations_Inner, nIterations_Outer, nIterations_Prim ) &
    !$ACC CREATE( GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#endif

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: nIterations_Inner, nIterations_Outer, nIterations_Prim ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$OMP               GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( nIterations_Inner, nIterations_Outer, nIterations_Prim ) &
    !$ACC DELETE( iX_B0, iX_E0, iX_B1, iX_E1, nZ, nX, PositionIndexZ, &
    !$ACC         GE_N, GX_N, CF_N, PF_N, AF_N, CR_N, PR_N )
#endif

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', MINVAL( nIterations_Outer(:) )
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', MAXVAL( nIterations_Outer(:) )
      WRITE(*,'(A6,A18,ES8.2E2)') &
        '', 'Iterations (Ave): ', DBLE( SUM( nIterations_Outer(:) ) ) / DBLE( nX_G )
      WRITE(*,*)

    END IF

    DEALLOCATE( GE_N, GX_N )
    DEALLOCATE( CF_N, PF_N, AF_N )
    DEALLOCATE( CR_N, PR_N )
    DEALLOCATE( PositionIndexZ )
    DEALLOCATE( nIterations_Inner )
    DEALLOCATE( nIterations_Outer )
    DEALLOCATE( nIterations_Prim )

    NULLIFY( N_P, G1_P, G2_P, G3_P )
    NULLIFY( J_P, H1_P, H2_P, H3_P )

    CALL FinalizeNeutrinoMatterSolver

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE MapDataForCollisions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in) :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iGE, iGF, iCF, iCR, iS
    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeZ
    INTEGER :: iN_E, iN_X

    ! --- Momentum Space Geometry ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iE, iNodeE )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iE, iNodeE ) &
    !$ACC PRESENT( GE_N, GE )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE, iNodeE )
#endif
    DO iGE  = 1, nGE
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      GE_N(iN_E,iGE) = GE(iNodeE,iE,iGE)

    END DO
    END DO

    ! --- Position Space Geometry ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, GX_N, GX )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iGF  = 1, nGF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

    ! --- Fluid Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, U_F )
#elif defined(THORNADO_OMP  )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#endif
    DO iCF  = 1, nCF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Radiation Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iNodeZ, iNodeX, iNodeE, iE, iX1, iX2, iX3 )
#endif
    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

      CR_N(iN_E,iS,iN_X,iCR) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapDataForCollisions


  SUBROUTINE ComputeAndMapIncrement &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, U_F, U_R, dU_F, dU_R )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      dt
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iCF, iCR, iS
    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeZ
    INTEGER :: iN_E, iN_X

    ! --- Fluid Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iN_X ) &
    !$ACC PRESENT( nX, iX_B0, iX_E0, dU_F, CF_N, U_F )
#elif defined(THORNADO_OMP   )
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

      IF ( QueryOpacity( U_F(iNodeX,iX1,iX2,iX3,iCF_D) / UnitD ) ) THEN
        dU_F(iNodeX,iX1,iX2,iX3,iCF) &
          = ( CF_N(iN_X,iCF) - U_F(iNodeX,iX1,iX2,iX3,iCF) ) / dt
      ELSE
        dU_F(iNodeX,iX1,iX2,iX3,iCF) = Zero
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, dU_R, CR_N, U_R )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iN_X, iN_E )
#endif
    DO iS = 1, nSpecies
    DO iCR = 1, nCR
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE = iE_B0, iE_E0
    DO iNodeZ = 1, nDOFZ

      iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX   ) + 1
      iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

      iN_X = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
      iN_E = iNodeE &
             + ( iE  - iE_B0    ) * nDOFE

      IF ( QueryOpacity( U_F(iNodeX,iX1,iX2,iX3,iCF_D) / UnitD ) ) THEN
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
          = ( CR_N(iN_E,iS,iN_X,iCR) - U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) ) / dt
      ELSE
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = Zero
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeAndMapIncrement


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
