MODULE TwoMoment_DiscretizationModule_Streaming

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    SqrtTiny, Pi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Streaming, &
    Timer_Streaming_LinearAlgebra, &
    Timer_Streaming_Divergence, &
    Timer_Streaming_ObserverCorrections, &
    Timer_Streaming_Derivatives, &
    Timer_Streaming_Eigenvalues, &
    Timer_Streaming_PrimitiveTwoMoment, &
    Timer_Streaming_NumericalFlux, &
    Timer_Streaming_Sources
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    EigenvaluesSymmetric3
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q, &
    Weights_E, &
    Weights_X1, &
    Weights_X2, &
    Weights_X3
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    dLdE_q, &
    L_X1_Dn, L_X1_Up, &
    dLdX1_q, &
    L_X2_Dn, L_X2_Up, &
    dLdX2_q, &
    L_X3_Dn, L_X3_Up, &
    dLdX3_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
    iGE_Ep1, &
    iGE_Ep2, &
    iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, &
    ComputeConserved_TwoMoment, &
    Flux_E, &
    Flux_X1, &
    Flux_X2, &
    Flux_X3, &
    ComputeEddingtonTensorComponents_uu, &
    ComputeEddingtonTensorComponents_ud, &
    NumericalFlux_LLF, &
    ComputeWeakDerivatives_X1, &
    ComputeWeakDerivatives_X2, &
    ComputeWeakDerivatives_X3, &
    FaceVelocity_X1, &
    FaceVelocity_X2, &
    FaceVelocity_X3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

  ! --- Geometry Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    uGF_K, uGF_F

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K, &
    Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F

  ! --- Conserved Fluid Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    uCF_K, uCF_L, uCF_R

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    uFD_K, uS1_K, uS2_K, uS3_K, &
    uFD_L, uS1_L, uS2_L, uS3_L, &
    uFD_R, uS1_R, uS2_R, uS3_R

  ! --- Conserved Radiation Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    uCR_K, uCR_L, uCR_R

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    uN_K, uG1_K, uG2_K, uG3_K, &
    uN_L, uG1_L, uG2_L, uG3_L, &
    uN_R, uG1_R, uG2_R, uG3_R

  ! --- Primitive Radiation Fields ---

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    uD_K, uI1_K, uI2_K, uI3_K, &
    uD_L, uI1_L, uI2_L, uI3_L, &
    uD_R, uI1_R, uI2_R, uI3_R

  ! --- Primitive Fluid Velocities ---

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    uV1_F, uV2_F, uV3_F, &
    uV1_K, uV2_K, uV3_K

  ! --- Fluxes ---

  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    NumericalFlux, NumericalFlux2, Flux_q

  ! --- Increment ---

  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    dU_Z

  ! --- Eigenvalues ---

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
    Alpha

  ! --- Derivatives ---

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3

  ! ---

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    nIterations_L, nIterations_R, nIterations_K

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    PositionIndexZ_F, PositionIndexZ_K

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableZ_F, IndexTableZ_K

  INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
  INTEGER :: nZ   (4), nK_X , nK_Z , nNodesX_K, nNodesZ_K
  INTEGER :: nZ_E (4), nE_X , nE_Z , nNodesZ_E
  INTEGER :: nZ_X1(4), nX1_X, nX1_Z, nNodesX_X1, nNodesZ_X1
  INTEGER :: nZ_X2(4), nX2_X, nX2_Z, nNodesX_X2, nNodesZ_X2
  INTEGER :: nZ_X3(4), nX3_X, nX3_Z, nNodesX_X3, nNodesZ_X3

  REAL(DP), PUBLIC :: OffGridFlux_TwoMoment(2*nCR)

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(inout) :: &
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

    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_Streaming )

    ASSOCIATE ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
                dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ1, dZ2, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, GE, GX, U_F, U_R ) &
    !$OMP MAP( alloc: dU_R )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ2, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         iX_B0, iX_E0, iX_B1, iX_E1, GE, GX, U_F, U_R ) &
    !$ACC CREATE( dU_R )
#endif

    CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

    CALL InitializeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    OffGridFlux_TwoMoment = Zero

    CALL TimersStart( Timer_Streaming_Divergence )

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL ComputeIncrement_Divergence_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL ComputeIncrement_Divergence_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL TimersStop( Timer_Streaming_Divergence )

    CALL TimersStart( Timer_Streaming_ObserverCorrections )

    CALL ComputeIncrement_ObserverCorrections &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL TimersStop( Timer_Streaming_ObserverCorrections )

    ! --- Multiply Inverse Mass Matrix ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(8) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, dZ1, dZ2, dZ3, dZ4, GX, GE, dU_R, Weights_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
              / ( Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                    * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                    * dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_TwoMoment_Explicit

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_R, U_R, U_F ) &
    !$OMP MAP( release: dZ1, dZ2, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               iX_B0, iX_E0, iX_B1, iX_E1, GE, GX )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_R, U_R, U_F ) &
    !$ACC DELETE( dZ1, dZ2, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         iX_B0, iX_E0, iX_B1, iX_E1, GE, GX )
    !$ACC WAIT
#endif

    END ASSOCIATE

    CALL TimersStop( Timer_Streaming )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_X1, iNodeX_X1
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: Flux_L(nCR), uCR_X1_L(nCR)
    REAL(DP) :: Flux_R(nCR), uCR_X1_R(nCR)
    REAL(DP) :: Flux_K(nCR)
    REAL(DP) :: uV1_L, uV2_L, uV3_L
    REAL(DP) :: uV1_R, uV2_R, uV3_R

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(2) ; iZP_E0(4) = iZ_E0(2)

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X1, nDOF_X1 )

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uGF_K, GX, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    DO iGF = iGF_h_1, iGF_h_3

      ! --- Face States (Average of Left and Right States) ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               uGF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               uGF_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               uGF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               uGF_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Recompute Geometry from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( uGF_F, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ2  = iZ_B0(2), iZ_E0(2)+1
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_2)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_3)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_SqrtGm) &
          = SQRT(   uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
                  * uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
                  * uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U_F, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Fluid Fields ---

    DO iCF = iCF_D, iCF_S3

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Face Velocity Components ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R ) &
    !$ACC PRESENT( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          uFD_L, uS1_L, uS2_L, uS3_L, &
    !$ACC          uFD_R, uS1_R, uS2_R, uS3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F)
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#endif
    DO iX_F = 1, nNodesX_X1

      ! --- Left State ---

      uV1_L = uS1_L(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_L(iX_F) )
      uV2_L = uS2_L(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_L(iX_F) )
      uV3_L = uS3_L(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_L(iX_F) )

      ! --- Right State ---

      uV1_R = uS1_R(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_R(iX_F) )
      uV2_R = uS2_R(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_R(iX_F) )
      uV3_R = uS3_R(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_R(iX_F) )

      CALL FaceVelocity_X1 &
             ( uV1_L, uV2_L, uV3_L, &
               uV1_R, uV2_R, uV3_R, &
               uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F) )

    END DO

    ! --- Permute Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iCR = 1, nCR
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iZ1,iZ3,iZ4,iS,iZ2,iCR) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Radiation Fields ---

    DO iCR = 1, nCR

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X1, nX1_Z, nDOFZ, One, L_X1_Up, nDOF_X1, &
               uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)-1,iCR), nDOFZ, Zero, &
               uCR_L(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCR), nDOF_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X1, nX1_Z, nDOFZ, One, L_X1_Dn, nDOF_X1, &
               uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCR), nDOFZ, Zero, &
               uCR_R(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCR), nDOF_X1 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    ! --- Left State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_L, uG1_L, uG2_L, uG3_L, &
             uD_L, uI1_L, uI2_L, uI3_L, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_L )

    ! --- Right State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_R, uG1_R, uG2_R, uG3_R, &
             uD_R, uI1_R, uI2_R, uI3_R, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_R )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_F, iNodeZ_X1, iNodeX_X1, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, &
    !$OMP          Flux_L, uCR_X1_L, Flux_R, uCR_X1_R )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_F, iNodeZ_X1, iNodeX_X1, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, &
    !$ACC          Flux_L, uCR_X1_L, Flux_R, uCR_X1_R ) &
    !$ACC PRESENT( uD_L, uI1_L, uI2_L, uI3_L, uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F, Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          NumericalFlux, NumericalFlux2, GE, SqrtGm_F, Weights_X1, dZ1, dZ3, dZ4, &
    !$ACC          nZ_X1, iZ_B0, PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_F, iNodeZ_X1, iNodeX_X1, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, &
    !$OMP          Flux_L, uCR_X1_L, Flux_R, uCR_X1_R )
#endif
    DO iZ_F = 1, nNodesZ_X1

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeE    = IndexTableZ_F(1,iZ_F)
      iNodeX_X1 = IndexTableZ_F(2,iZ_F)
      iZ1       = IndexTableZ_F(3,iZ_F)
      iZ3       = IndexTableZ_F(4,iZ_F)
      iZ4       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ2       = IndexTableZ_F(7,iZ_F)

      iNodeZ_X1 = iNodeE + ( iNodeX_X1 - 1 ) * nDOFE

      ! --- Left State Flux ---

      Flux_L &
        = Flux_X1 &
            ( uD_L (iZ_F), uI1_L(iZ_F), &
              uI2_L(iZ_F), uI3_L(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_L (iZ_F), uI1_L(iZ_F), &
               uI2_L(iZ_F), uI3_L(iZ_F), &
               uCR_X1_L(iCR_N ), uCR_X1_L(iCR_G1), &
               uCR_X1_L(iCR_G2), uCR_X1_L(iCR_G3), &
               uV1_F(iX_F), Zero, Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Right State Flux ---

      Flux_R &
        = Flux_X1 &
            ( uD_R (iZ_F), uI1_R(iZ_F), &
              uI2_R(iZ_F), uI3_R(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_R (iZ_F), uI1_R(iZ_F), &
               uI2_R(iZ_F), uI3_R(iZ_F), &
               uCR_X1_R(iCR_N ), uCR_X1_R(iCR_G1), &
               uCR_X1_R(iCR_G2), uCR_X1_R(iCR_G3), &
               uV1_F(iX_F), Zero, Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Numerical Flux ---

      DO iCR = 1, nCR

        NumericalFlux(iNodeZ_X1,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = NumericalFlux_LLF &
              ( uCR_X1_L(iCR), uCR_X1_R(iCR), Flux_L(iCR), Flux_R(iCR), One )

        NumericalFlux(iNodeZ_X1,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_F(iX_F) &
              * NumericalFlux(iNodeZ_X1,iCR,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

      ! --- Energy and Momentum Fluxes for Conservation Tally ---

      ! --- Energy ---

      NumericalFlux2(iNodeZ_X1,iCR_N,iZ1,iZ3,iZ4,iS,iZ2) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X1,iCR_N,iZ1,iZ3,iZ4,iS,iZ2) &
                + uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2) &
                + uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2) &
                + uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2) )

      ! --- Momentum 1 ---

      NumericalFlux2(iNodeZ_X1,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X1,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2) &
                + Gm_dd_11_F(iX_F) * uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_N,iZ1,iZ3,iZ4,iS,iZ2) )

      ! --- Momentum 2 ---

      NumericalFlux2(iNodeZ_X1,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X1,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2) &
                + Gm_dd_22_F(iX_F) * uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_N,iZ1,iZ3,iZ4,iS,iZ2) )

      ! --- Momentum 3 ---

      NumericalFlux2(iNodeZ_X1,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X1,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2) &
                + Gm_dd_33_F(iX_F) * uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X1,iCR_N,iZ1,iZ3,iZ4,iS,iZ2) )

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_Z, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    CALL ComputeOffGridFlux( iZP_B0, iZP_E0, nDOF_X1, NumericalFlux, NumericalFlux2 )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( uV1_K, uV2_K, uV3_K, uFD_K, uS1_K, uS2_K, uS3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO
#endif
    DO iX_K = 1, nNodesX_K

      uV1_K(iX_K) = uS1_K(iX_K) / ( Gm_dd_11_K(iX_K) * uFD_K(iX_K) )
      uV2_K(iX_K) = uS2_K(iX_K) / ( Gm_dd_22_K(iX_K) * uFD_K(iX_K) )
      uV3_K(iX_K) = uS3_K(iX_K) / ( Gm_dd_33_K(iX_K) * uFD_K(iX_K) )

    END DO

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    CALL ComputePrimitive_TwoMoment &
           ( uN_K, uG1_K, uG2_K, uG3_K, &
             uD_K, uI1_K, uI2_K, uI3_K, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_K, nIterations_K )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, Flux_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, Flux_K ) &
    !$ACC PRESENT( uD_K, uI1_K, uI2_K, uI3_K, uV1_K, uV2_K, uV3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          Flux_q, GE, SqrtGm_K, Weights_q, dZ1, dZ3, dZ4, &
    !$ACC          nZ, iZ_B0, PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ2, iZ3, iZ4, iS, Flux_K )
#endif
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ1    = IndexTableZ_K(3,iZ_K)
      iZ3    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ2    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_X1 &
            ( uD_K (iZ_K), uI1_K(iZ_K), &
              uI2_K(iZ_K), uI3_K(iZ_K), &
              uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
              Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K) )

      DO iCR = 1, nCR

        Flux_q(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCR)

      END DO

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, dU_Z, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_Z(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence_X

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_X2, iNodeX_X2
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: Flux_L(nCR), uCR_X2_L(nCR)
    REAL(DP) :: Flux_R(nCR), uCR_X2_R(nCR)
    REAL(DP) :: Flux_K(nCR)
    REAL(DP) :: uV1_L, uV2_L, uV3_L
    REAL(DP) :: uV1_R, uV2_R, uV3_R
    REAL(DP) :: SUM_N, SUM_G1, SUM_G2, SUM_G3

    IF( iZ_E0(3) .EQ. iZ_B0(3) ) RETURN

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(2) ; iZP_E0(2) = iZ_E0(2)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(3) ; iZP_E0(4) = iZ_E0(3)

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ1, dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X2, nDOF_X2 )

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uGF_K, GX, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iZ3 = iZ_B0(3)-1, iZ_E0(3)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ2,iZ4,iZ3,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    DO iGF = iGF_h_1, iGF_h_3

      ! --- Face States (Average of Left and Right States) ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               uGF_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1,iGF), nDOFX, Zero, &
               uGF_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX_X2 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               uGF_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX, Half, &
               uGF_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX_X2 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Recompute Geometry from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( uGF_F, iZ_B0, iZ_E0, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ3  = iZ_B0(3), iZ_E0(3)+1
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2

        uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_11) &
          = MAX( uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_h_1)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_22) &
          = MAX( uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_h_2)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_33) &
          = MAX( uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_h_3)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_SqrtGm) &
          = SQRT(   uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_11) &
                  * uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_22) &
                  * uGF_F(iNodeX,iZ2,iZ4,iZ3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U_F, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ3 = iZ_B0(3)-1, iZ_E0(3)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ2,iZ4,iZ3,iCF) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Fluid Fields ---

    DO iCF = iCF_D, iCF_S3

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               uCF_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iCF), nDOFX_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               uCF_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iCF), nDOFX_X2 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Face Velocity Components ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R ) &
    !$ACC PRESENT( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          uFD_L, uS1_L, uS2_L, uS3_L, &
    !$ACC          uFD_R, uS1_R, uS2_R, uS3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F)
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#endif
    DO iX_F = 1, nNodesX_X2

      ! --- Left State ---

      uV1_L = uS1_L(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_L(iX_F) )
      uV2_L = uS2_L(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_L(iX_F) )
      uV3_L = uS3_L(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_L(iX_F) )

      ! --- Right State ---

      uV1_R = uS1_R(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_R(iX_F) )
      uV2_R = uS2_R(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_R(iX_F) )
      uV3_R = uS3_R(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_R(iX_F) )

      CALL FaceVelocity_X2 &
             ( uV1_L, uV2_L, uV3_L, &
               uV1_R, uV2_R, uV3_R, &
               uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F) )

    END DO

    ! --- Permute Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iCR = 1, nCR
    DO iZ3 = iZ_B0(3)-1, iZ_E0(3)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iZ1,iZ2,iZ4,iS,iZ3,iCR) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Radiation Fields ---

    DO iCR = 1, nCR

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X2, nX2_Z, nDOFZ, One, L_X2_Up, nDOF_X2, &
               uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)-1,iCR), nDOFZ, Zero, &
               uCR_L(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ,iCR), nDOF_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X2, nX2_Z, nDOFZ, One, L_X2_Dn, nDOF_X2, &
               uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ,iCR), nDOFZ, Zero, &
               uCR_R(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ,iCR), nDOF_X2 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    ! --- Left State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_L, uG1_L, uG2_L, uG3_L, &
             uD_L, uI1_L, uI2_L, uI3_L, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_L )

    ! --- Right State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_R, uG1_R, uG2_R, uG3_R, &
             uD_R, uI1_R, uI2_R, uI3_R, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_R )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_F, iNodeZ_X2, iNodeX_X2, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, &
    !$OMP          Flux_L, uCR_X2_L, Flux_R, uCR_X2_R )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_F, iNodeZ_X2, iNodeX_X2, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, &
    !$ACC          Flux_L, uCR_X2_L, Flux_R, uCR_X2_R ) &
    !$ACC PRESENT( uD_L, uI1_L, uI2_L, uI3_L, uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F, Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          NumericalFlux, NumericalFlux2, GE, SqrtGm_F, Weights_X2, dZ1, dZ2, dZ4, &
    !$ACC          nZ_X2, iZ_B0, PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_F, iNodeZ_X2, iNodeX_X2, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, &
    !$OMP          Flux_L, uCR_X2_L, Flux_R, uCR_X2_R )
#endif
    DO iZ_F = 1, nNodesZ_X2

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeE    = IndexTableZ_F(1,iZ_F)
      iNodeX_X2 = IndexTableZ_F(2,iZ_F)
      iZ1       = IndexTableZ_F(3,iZ_F)
      iZ2       = IndexTableZ_F(4,iZ_F)
      iZ4       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ3       = IndexTableZ_F(7,iZ_F)

      iNodeZ_X2 = iNodeE + ( iNodeX_X2 - 1 ) * nDOFE

      ! --- Left State Flux ---

      Flux_L &
        = Flux_X2 &
            ( uD_L (iZ_F), uI1_L(iZ_F), &
              uI2_L(iZ_F), uI3_L(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_L (iZ_F), uI1_L(iZ_F), &
               uI2_L(iZ_F), uI3_L(iZ_F), &
               uCR_X2_L(iCR_N ), uCR_X2_L(iCR_G1), &
               uCR_X2_L(iCR_G2), uCR_X2_L(iCR_G3), &
               Zero, uV2_F(iX_F), Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Right State Flux ---

      Flux_R &
        = Flux_X2 &
            ( uD_R (iZ_F), uI1_R(iZ_F), &
              uI2_R(iZ_F), uI3_R(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_R (iZ_F), uI1_R(iZ_F), &
               uI2_R(iZ_F), uI3_R(iZ_F), &
               uCR_X2_R(iCR_N ), uCR_X2_R(iCR_G1), &
               uCR_X2_R(iCR_G2), uCR_X2_R(iCR_G3), &
               Zero, uV2_F(iX_F), Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Numerical Flux ---

      DO iCR = 1, nCR

        NumericalFlux(iNodeZ_X2,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
          = NumericalFlux_LLF &
              ( uCR_X2_L(iCR), uCR_X2_R(iCR), Flux_L(iCR), Flux_R(iCR), One )

        NumericalFlux(iNodeZ_X2,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
          = dZ1(iZ1) * dZ2(iZ2) * dZ4(iZ4) &
              * Weights_X2(iNodeZ_X2) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_F(iX_F) &
              * NumericalFlux(iNodeZ_X2,iCR,iZ1,iZ2,iZ4,iS,iZ3)

      END DO

      ! --- Energy and Momentum Fluxes for Conservation Tally ---

      ! --- Energy ---

      NumericalFlux2(iNodeZ_X2,iCR_N,iZ1,iZ2,iZ4,iS,iZ3) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X2,iCR_N,iZ1,iZ2,iZ4,iS,iZ3) &
                + uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3) &
                + uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3) &
                + uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3) )

      ! --- Momentum 1 ---

      NumericalFlux2(iNodeZ_X2,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X2,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3) &
                + Gm_dd_11_F(iX_F) * uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_N,iZ1,iZ2,iZ4,iS,iZ3) )

      ! --- Momentum 2 ---

      NumericalFlux2(iNodeZ_X2,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X2,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3) &
                + Gm_dd_22_F(iX_F) * uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_N,iZ1,iZ2,iZ4,iS,iZ3) )

      ! --- Momentum 3 ---

      NumericalFlux2(iNodeZ_X2,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X2,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3) &
                + Gm_dd_33_F(iX_F) * uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X2,iCR_N,iZ1,iZ2,iZ4,iS,iZ3) )

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X2, + One, L_X2_Dn, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), &
             nDOF_X2, Zero, dU_Z, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X2, - One, L_X2_Up, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)+1), &
             nDOF_X2, One,  dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    CALL ComputeOffGridFlux( iZP_B0, iZP_E0, nDOF_X2, NumericalFlux, NumericalFlux2 )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( uV1_K, uV2_K, uV3_K, uFD_K, uS1_K, uS2_K, uS3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO
#endif
    DO iX_K = 1, nNodesX_K

      uV1_K(iX_K) = uS1_K(iX_K) / ( Gm_dd_11_K(iX_K) * uFD_K(iX_K) )
      uV2_K(iX_K) = uS2_K(iX_K) / ( Gm_dd_22_K(iX_K) * uFD_K(iX_K) )
      uV3_K(iX_K) = uS3_K(iX_K) / ( Gm_dd_33_K(iX_K) * uFD_K(iX_K) )

    END DO

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    CALL ComputePrimitive_TwoMoment &
           ( uN_K, uG1_K, uG2_K, uG3_K, &
             uD_K, uI1_K, uI2_K, uI3_K, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_K, nIterations_K )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, Flux_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, Flux_K ) &
    !$ACC PRESENT( uD_K, uI1_K, uI2_K, uI3_K, uV1_K, uV2_K, uV3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          Flux_q, GE, SqrtGm_K, Weights_q, dZ1, dZ2, dZ4, &
    !$ACC          nZ, iZ_B0, PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ3, iZ2, iZ4, iS, Flux_K )
#endif
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ1    = IndexTableZ_K(3,iZ_K)
      iZ2    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ3    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_X2 &
            ( uD_K (iZ_K), uI1_K(iZ_K), &
              uI2_K(iZ_K), uI3_K(iZ_K), &
              uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
              Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K) )

      DO iCR = 1, nCR

        Flux_q(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
          = dZ1(iZ1) * dZ2(iZ2) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCR)

      END DO

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOFZ, One, dLdX2_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, dU_Z, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_Z(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence_X

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ1, dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ1, dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X3,X4} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_X3, iNodeX_X3
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: Flux_L(nCR), uCR_X3_L(nCR)
    REAL(DP) :: Flux_R(nCR), uCR_X3_R(nCR)
    REAL(DP) :: Flux_K(nCR)
    REAL(DP) :: uV1_L, uV2_L, uV3_L
    REAL(DP) :: uV1_R, uV2_R, uV3_R
    REAL(DP) :: SUM_N, SUM_G1, SUM_G2, SUM_G3

    IF( iZ_E0(4) .EQ. iZ_B0(4) ) RETURN

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(2) ; iZP_E0(2) = iZ_E0(2)
    iZP_B0(3) = iZ_B0(3) ; iZP_E0(3) = iZ_E0(3)
    iZP_B0(4) = iZ_B0(4) ; iZP_E0(4) = iZ_E0(4)

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ1, dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X3, nDOF_X3 )

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uGF_K, GX, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iZ4 = iZ_B0(4)-1, iZ_E0(4)+1
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ2,iZ3,iZ4,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    DO iGF = iGF_h_1, iGF_h_3

      ! --- Face States (Average of Left and Right States) ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               uGF_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1,iGF), nDOFX, Zero, &
               uGF_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX_X3 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               uGF_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX, Half, &
               uGF_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX_X3 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Recompute Geometry from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( uGF_F, iZ_B0, iZ_E0, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ4  = iZ_B0(4), iZ_E0(4)+1
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) &
          = MAX( uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_h_1)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) &
          = MAX( uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_h_2)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) &
          = MAX( uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_h_3)**2, SqrtTiny )
        uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
          = SQRT(   uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) &
                  * uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) &
                  * uGF_F(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U_F, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ4 = iZ_B0(4)-1, iZ_E0(4)+1
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ2,iZ3,iZ4,iCF) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Fluid Fields ---

    DO iCF = iCF_D, iCF_S3

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               uCF_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iCF), nDOFX_X3 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               uCF_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iCF), nDOFX_X3 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Face Velocity Components ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R ) &
    !$ACC PRESENT( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          uFD_L, uS1_L, uS2_L, uS3_L, &
    !$ACC          uFD_R, uS1_R, uS2_R, uS3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F)
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( uV1_L, uV2_L, uV3_L, uV1_R, uV2_R, uV3_R )
#endif
    DO iX_F = 1, nNodesX_X3

      ! --- Left State ---

      uV1_L = uS1_L(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_L(iX_F) )
      uV2_L = uS2_L(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_L(iX_F) )
      uV3_L = uS3_L(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_L(iX_F) )

      ! --- Right State ---

      uV1_R = uS1_R(iX_F) / ( Gm_dd_11_F(iX_F) * uFD_R(iX_F) )
      uV2_R = uS2_R(iX_F) / ( Gm_dd_22_F(iX_F) * uFD_R(iX_F) )
      uV3_R = uS3_R(iX_F) / ( Gm_dd_33_F(iX_F) * uFD_R(iX_F) )

      CALL FaceVelocity_X3 &
             ( uV1_L, uV2_L, uV3_L, &
               uV1_R, uV2_R, uV3_R, &
               uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F) )

    END DO

    ! --- Permute Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4)-1, iZ_E0(4)+1
    DO iS  = 1, nSpecies
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iZ1,iZ2,iZ3,iS,iZ4,iCR) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Radiation Fields ---

    DO iCR = 1, nCR

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X3, nX3_Z, nDOFZ, One, L_X3_Up, nDOF_X3, &
               uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)-1,iCR), nDOFZ, Zero, &
               uCR_L(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ,iCR), nDOF_X3 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X3, nX3_Z, nDOFZ, One, L_X3_Dn, nDOF_X3, &
               uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ,iCR), nDOFZ, Zero, &
               uCR_R(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ,iCR), nDOF_X3 )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    ! --- Left State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_L, uG1_L, uG2_L, uG3_L, &
             uD_L, uI1_L, uI2_L, uI3_L, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_L )

    ! --- Right State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_R, uG1_R, uG2_R, uG3_R, &
             uD_R, uI1_R, uI2_R, uI3_R, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             PositionIndexZ_F, nIterations_R )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_F, iNodeZ_X3, iNodeX_X3, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, &
    !$OMP          Flux_L, uCR_X3_L, Flux_R, uCR_X3_R )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_F, iNodeZ_X3, iNodeX_X3, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, &
    !$ACC          Flux_L, uCR_X3_L, Flux_R, uCR_X3_R ) &
    !$ACC PRESENT( uD_L, uI1_L, uI2_L, uI3_L, uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC          uV1_F, uV2_F, uV3_F, Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
    !$ACC          NumericalFlux, NumericalFlux2, GE, SqrtGm_F, Weights_X3, dZ1, dZ2, dZ3, &
    !$ACC          nZ_X3, iZ_B0, PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_F, iNodeZ_X3, iNodeX_X3, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, &
    !$OMP          Flux_L, uCR_X3_L, Flux_R, uCR_X3_R )
#endif
    DO iZ_F = 1, nNodesZ_X3

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeE    = IndexTableZ_F(1,iZ_F)
      iNodeX_X3 = IndexTableZ_F(2,iZ_F)
      iZ1       = IndexTableZ_F(3,iZ_F)
      iZ2       = IndexTableZ_F(4,iZ_F)
      iZ3       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ4       = IndexTableZ_F(7,iZ_F)

      iNodeZ_X3 = iNodeE + ( iNodeX_X3 - 1 ) * nDOFE

        ! --- Left State Flux ---

      Flux_L &
        = Flux_X3 &
            ( uD_L (iZ_F), uI1_L(iZ_F), &
              uI2_L(iZ_F), uI3_L(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_L (iZ_F), uI1_L(iZ_F), &
               uI2_L(iZ_F), uI3_L(iZ_F), &
               uCR_X3_L(iCR_N ), uCR_X3_L(iCR_G1), &
               uCR_X3_L(iCR_G2), uCR_X3_L(iCR_G3), &
               Zero, Zero, uV3_F(iX_F), &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Right State Flux ---

      Flux_R &
        = Flux_X3 &
            ( uD_R (iZ_F), uI1_R(iZ_F), &
              uI2_R(iZ_F), uI3_R(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_R (iZ_F), uI1_R(iZ_F), &
               uI2_R(iZ_F), uI3_R(iZ_F), &
               uCR_X3_R(iCR_N ), uCR_X3_R(iCR_G1), &
               uCR_X3_R(iCR_G2), uCR_X3_R(iCR_G3), &
               Zero, Zero, uV3_F(iX_F), &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Numerical Flux ---

      DO iCR = 1, nCR

        NumericalFlux(iNodeZ_X3,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
          = NumericalFlux_LLF &
              ( uCR_X3_L(iCR), uCR_X3_R(iCR), Flux_L(iCR), Flux_R(iCR), One )

        NumericalFlux(iNodeZ_X3,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
          = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) &
              * Weights_X3(iNodeZ_X3) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_F(iX_F) &
              * NumericalFlux(iNodeZ_X3,iCR,iZ1,iZ2,iZ3,iS,iZ4)

      END DO

      ! --- Energy and Momentum Fluxes for Conservation Tally ---

      ! --- Energy ---

      NumericalFlux2(iNodeZ_X3,iCR_N,iZ1,iZ2,iZ3,iS,iZ4) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X3,iCR_N,iZ1,iZ2,iZ3,iS,iZ4) &
                + uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4) &
                + uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4) &
                + uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4) )

      ! --- Momentum 1 ---

      NumericalFlux2(iNodeZ_X3,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X3,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4) &
                + Gm_dd_11_F(iX_F) * uV1_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_N,iZ1,iZ2,iZ3,iS,iZ4) )

      ! --- Momentum 2 ---

      NumericalFlux2(iNodeZ_X3,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X3,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4) &
                + Gm_dd_22_F(iX_F) * uV2_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_N,iZ1,iZ2,iZ3,iS,iZ4) )

      ! --- Momentum 3 ---

      NumericalFlux2(iNodeZ_X3,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4) &
        = GE(iNodeE,iZ1,iGE_Ep1) &
            * ( NumericalFlux(iNodeZ_X3,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4) &
                + Gm_dd_33_F(iX_F) * uV3_F(iX_F) &
                    * NumericalFlux(iNodeZ_X3,iCR_N,iZ1,iZ2,iZ3,iS,iZ4) )

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X3, + One, L_X3_Dn, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), &
             nDOF_X3, Zero, dU_Z, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X3, - One, L_X3_Up, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)+1), &
             nDOF_X3, One,  dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    CALL ComputeOffGridFlux( iZP_B0, iZP_E0, nDOF_X3, NumericalFlux, NumericalFlux2 )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( uV1_K, uV2_K, uV3_K, uFD_K, uS1_K, uS2_K, uS3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO
#endif
    DO iX_K = 1, nNodesX_K

      uV1_K(iX_K) = uS1_K(iX_K) / ( Gm_dd_11_K(iX_K) * uFD_K(iX_K) )
      uV2_K(iX_K) = uS2_K(iX_K) / ( Gm_dd_22_K(iX_K) * uFD_K(iX_K) )
      uV3_K(iX_K) = uS3_K(iX_K) / ( Gm_dd_33_K(iX_K) * uFD_K(iX_K) )

    END DO

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    CALL ComputePrimitive_TwoMoment &
           ( uN_K, uG1_K, uG2_K, uG3_K, &
             uD_K, uI1_K, uI2_K, uI3_K, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_K, nIterations_K )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, Flux_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, Flux_K ) &
    !$ACC PRESENT( uD_K, uI1_K, uI2_K, uI3_K, uV1_K, uV2_K, uV3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          Flux_q, GE, SqrtGm_K, Weights_q, dZ1, dZ2, dZ3, &
    !$ACC          nZ, iZ_B0, PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_K, iNodeZ, iNodeX, iNodeE, iZ1, iZ4, iZ2, iZ3, iS, Flux_K )
#endif
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ1    = IndexTableZ_K(3,iZ_K)
      iZ2    = IndexTableZ_K(4,iZ_K)
      iZ3    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ4    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_X3 &
            ( uD_K (iZ_K), uI1_K(iZ_K), &
              uI2_K(iZ_K), uI3_K(iZ_K), &
              uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
              Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K) )

      DO iCR = 1, nCR

        Flux_q(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
          = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCR)

      END DO

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOFZ, One, dLdX3_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, dU_Z, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_Z(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence_X

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ1, dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ1, dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_ObserverCorrections &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_E
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: EdgeEnergyCubed, Beta
    REAL(DP) :: A(3,3), Lambda(3)
    REAL(DP) :: k_uu(3,3), S_uu_11, S_uu_22, S_uu_33
    REAL(DP) :: Flux_K(nCR), dFlux_K(nPR)
    REAL(DP) :: Flux_L(nCR), uPR_L(nPR)
    REAL(DP) :: Flux_R(nCR), uPR_R(nPR)
    REAL(DP) :: SUM_N, SUM_G1, SUM_G2, SUM_G3

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(2) ; iZP_E0(1) = iZ_E0(2)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(1) ; iZP_E0(4) = iZ_E0(1)

    ASSOCIATE &
      ( xZ1 => MeshE    % Center, &
        dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: xZ1, dZ1, dZ2, dZ3, dZ4, &
    !$OMP          iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( xZ1, dZ1, dZ2, dZ3, dZ4, &
    !$ACC         iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    CALL InitializeIncrement_ObserverCorrections &
           ( iZ_B0, iZ_E0 )

    ! --- Calculate Weak Derivatives ---

    CALL TimersStart( Timer_Streaming_Derivatives )

    CALL ComputeWeakDerivatives_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX1, dV_d_dX1, dGm_dd_dX1 )

    CALL ComputeWeakDerivatives_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX2, dV_d_dX2, dGm_dd_dX2 )

    CALL ComputeWeakDerivatives_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )

    CALL TimersStop( Timer_Streaming_Derivatives )

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uGF_K, GX, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ2,iZ3,iZ4,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U_F, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ2,iZ3,iZ4,iCF) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U_R, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iCR = 1, nCR
    DO iZ1 = iZ_B0(1)-1, iZ_E0(1)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iZ2,iZ3,iZ4,iS,iZ1,iCR) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute Primitive Fluid ---

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( uV1_K, uV2_K, uV3_K, uFD_K, uS1_K, uS2_K, uS3_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO
#endif
    DO iX_K = 1, nNodesX_K

      uV1_K(iX_K) = uS1_K(iX_K) / ( Gm_dd_11_K(iX_K) * uFD_K(iX_K) )
      uV2_K(iX_K) = uS2_K(iX_K) / ( Gm_dd_22_K(iX_K) * uFD_K(iX_K) )
      uV3_K(iX_K) = uS3_K(iX_K) / ( Gm_dd_33_K(iX_K) * uFD_K(iX_K) )

    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Radiation Fields ---

    DO iCR = 1, nCR

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_E, nE_Z, nDOFZ, One, L_E_Up, nDOF_E, &
               uCR_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)-1,iCR), nDOFZ, Zero, &
               uCR_L(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCR), nDOF_E )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_E, nE_Z, nDOFZ, One, L_E_Dn, nDOF_E, &
               uCR_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCR), nDOFZ, Zero, &
               uCR_R(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCR), nDOF_E )

    END DO

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Eigenvalues ---

    CALL TimersStart( Timer_Streaming_Eigenvalues )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( A, Lambda )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( A, Lambda ) &
    !$ACC PRESENT( dV_u_dX1, dV_u_dX2, dV_u_dX3, Alpha, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( A, Lambda )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        ! --- Quadratic Form Matrix ---

        A(:,1) = Half * [ Two * dV_u_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                                dV_u_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
                              + dV_u_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                                dV_u_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
                              + dV_u_dX1(iNodeX,3,iZ2,iZ3,iZ4) ]
        A(:,2) = Half * [       dV_u_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
                              + dV_u_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                          Two * dV_u_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                                dV_u_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
                              + dV_u_dX2(iNodeX,3,iZ2,iZ3,iZ4) ]
        A(:,3) = Half * [       dV_u_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
                              + dV_u_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                                dV_u_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
                              + dV_u_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                          Two * dV_u_dX3(iNodeX,3,iZ2,iZ3,iZ4) ]

        CALL EigenvaluesSymmetric3( A, Lambda )

        Alpha(iNodeX,iZ2,iZ3,iZ4) = MAXVAL( ABS( Lambda ) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Streaming_Eigenvalues )

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    ! --- Left State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_L, uG1_L, uG2_L, uG3_L, &
             uD_L, uI1_L, uI2_L, uI3_L, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_F, nIterations_L )

    ! --- Right State Primitive ---

    CALL ComputePrimitive_TwoMoment &
           ( uN_R, uG1_R, uG2_R, uG3_R, &
             uD_R, uI1_R, uI2_R, uI3_R, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_F, nIterations_R )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_F, iNodeZ_E, iZ1, iZ2, iZ3, iZ4, iS, uPR_L, uPR_R, &
    !$OMP          Flux_L, Flux_R, A, Lambda, EdgeEnergyCubed )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_F, iNodeZ_E, iZ1, iZ2, iZ3, iZ4, iS, uPR_L, uPR_R, &
    !$ACC          Flux_L, Flux_R, EdgeEnergyCubed ) &
    !$ACC PRESENT( dV_u_dX1, dV_u_dX2, dV_u_dX3, uV1_K, uV2_K, uV3_K, &
    !$ACC          dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          uD_L, uI1_L, uI2_L, uI3_L, uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC          NumericalFlux, NumericalFlux2, SqrtGm_K, Weights_E, xZ1, dZ2, dZ3, dZ4, &
    !$ACC          Alpha, PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_F, iNodeZ_E, iZ1, iZ2, iZ3, iZ4, iS, uPR_L, uPR_R, &
    !$OMP          Flux_L, Flux_R, A, Lambda, EdgeEnergyCubed )
#endif
    DO iZ_F = 1, nNodesZ_E

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeZ_E  = IndexTableZ_F(2,iZ_F) ! = iNodeX
      iZ2       = IndexTableZ_F(3,iZ_F)
      iZ3       = IndexTableZ_F(4,iZ_F)
      iZ4       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ1       = IndexTableZ_F(7,iZ_F)

      ! --- Left State Flux ---

      Flux_L &
        = Flux_E( uD_L (iZ_F), uI1_L(iZ_F), uI2_L(iZ_F), uI3_L(iZ_F), &
                  uV1_K(iX_F), uV2_K(iX_F), uV3_K(iX_F), &
                  dV_u_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_F), Gm_dd_22_K(iX_F), Gm_dd_33_K(iX_F), &
                  dGm_dd_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4) )

      uPR_L = [ uD_L(iZ_F), uI1_L(iZ_F), uI2_L(iZ_F), uI3_L(iZ_F) ]

      ! --- Right State Flux ---

      Flux_R &
        = Flux_E( uD_R (iZ_F), uI1_R(iZ_F), uI2_R(iZ_F), uI3_R(iZ_F), &
                  uV1_K(iX_F), uV2_K(iX_F), uV3_K(iX_F), &
                  dV_u_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_F), Gm_dd_22_K(iX_F), Gm_dd_33_K(iX_F), &
                  dGm_dd_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4) )

      uPR_R = [ uD_R(iZ_F), uI1_R(iZ_F), uI2_R(iZ_F), uI3_R(iZ_F) ]

      ! --- Numerical Flux (Local Lax-Friedrichs) ---

      EdgeEnergyCubed = ( xZ1(iZ1) - Half * dZ1(iZ1) )**3

      DO iCR = 1, nCR

        NumericalFlux(iNodeZ_E,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = NumericalFlux_LLF &
              ( uPR_L(iCR), uPR_R(iCR), Flux_L(iCR), Flux_R(iCR), &
                Alpha(iNodeZ_E,iZ2,iZ3,iZ4) )

        NumericalFlux(iNodeZ_E,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_E(iNodeZ_E) * EdgeEnergyCubed &
              * SqrtGm_K(iX_F) &
              * NumericalFlux(iNodeZ_E,iCR,iZ2,iZ3,iZ4,iS,iZ1)

      END DO

      ! --- Energy and Momentum Fluxes for Conservation Tally ---

      EdgeEnergyCubed = ( xZ1(iZ1) - Half * dZ1(iZ1) )

      ! --- Energy ---

      NumericalFlux2(iNodeZ_E,iCR_N,iZ2,iZ3,iZ4,iS,iZ1) &
        = EdgeEnergyCubed &
            * ( NumericalFlux(iNodeZ_E,iCR_N,iZ2,iZ3,iZ4,iS,iZ1) &
                + uV1_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
                + uV2_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
                + uV3_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) )

      ! --- Momentum 1 ---

      NumericalFlux2(iNodeZ_E,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
        = EdgeEnergyCubed &
            * ( NumericalFlux(iNodeZ_E,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
                + Gm_dd_11_K(iX_F) * uV1_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_N,iZ2,iZ3,iZ4,iS,iZ1) )

      ! --- Momentum 2 ---

      NumericalFlux2(iNodeZ_E,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
        = EdgeEnergyCubed &
            * ( NumericalFlux(iNodeZ_E,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
                + Gm_dd_22_K(iX_F) * uV2_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_N,iZ2,iZ3,iZ4,iS,iZ1) )

      ! --- Momentum 3 ---

      NumericalFlux2(iNodeZ_E,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
        = EdgeEnergyCubed &
            * ( NumericalFlux(iNodeZ_E,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
                + Gm_dd_33_K(iX_F) * uV3_K(iX_F) &
                    * NumericalFlux(iNodeZ_E,iCR_N,iZ2,iZ3,iZ4,iS,iZ1) )

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contributions from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_E, + One, L_E_Dn, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), &
             nDOF_E, Zero, dU_Z, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_E, - One, L_E_Up, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)+1), &
             nDOF_E, One,  dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    CALL ComputeOffGridFlux( iZP_B0, iZP_E0, nDOF_E, NumericalFlux, NumericalFlux2 )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart( Timer_Streaming_PrimitiveTwoMoment )

    CALL ComputePrimitive_TwoMoment &
           ( uN_K, uG1_K, uG2_K, uG3_K, &
             uD_K, uI1_K, uI2_K, uI3_K, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             PositionIndexZ_K, nIterations_K )

    CALL TimersStop( Timer_Streaming_PrimitiveTwoMoment )

    CALL TimersStart( Timer_Streaming_NumericalFlux )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, &
    !$OMP          Flux_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, &
    !$ACC          Flux_K ) &
    !$ACC PRESENT( dV_u_dX1, dV_u_dX2, dV_u_dX3, uV1_K, uV2_K, uV3_K, &
    !$ACC          dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          uD_K, uI1_K, uI2_K, uI3_K, &
    !$ACC          Flux_q, GE, SqrtGm_K, Weights_q, dZ2, dZ3, dZ4, &
    !$ACC          PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, &
    !$OMP          Flux_K )
#endif
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ2    = IndexTableZ_K(3,iZ_K)
      iZ3    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ1    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_E( uD_K(iZ_K), uI1_K(iZ_K), uI2_K(iZ_K), uI3_K(iZ_K), &
                  uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
                  dV_u_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K), &
                  dGm_dd_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeX,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,3,iZ2,iZ3,iZ4) )

      DO iCR = 1, nCR

        Flux_q(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)  &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep3) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCR)

      END DO

    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux )

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOFZ, One, dLdE_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    !--------------------------------------------------
    !--- Volume Sources (Number Flux Equation Only) ---
    !--------------------------------------------------

    CALL TimersStart( Timer_Streaming_Sources )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, Flux_K, &
    !$OMP          Beta, dFlux_K, k_uu, S_uu_11, S_uu_22, S_uu_33 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, Flux_K, &
    !$ACC          Beta, dFlux_K, k_uu, S_uu_11, S_uu_22, S_uu_33 ) &
    !$ACC PRESENT( dV_u_dX1, dV_u_dX2, dV_u_dX3, uV1_K, uV2_K, uV3_K, &
    !$ACC          dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
    !$ACC          uD_K, uI1_K, uI2_K, uI3_K, &
    !$ACC          dU_Z, GE, SqrtGm_K, Weights_q, dZ1, dZ2, dZ3, dZ4, &
    !$ACC          PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX_K, iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iS, Flux_K, &
    !$OMP          Beta, dFlux_K, k_uu, S_uu_11, S_uu_22, S_uu_33 )
#endif
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ2    = IndexTableZ_K(3,iZ_K)
      iZ3    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ1    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_E( uD_K(iZ_K), uI1_K(iZ_K), uI2_K(iZ_K), uI3_K(iZ_K), &
                  uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
                  dV_u_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX1(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX2(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_u_dX3(iNodeX,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K), &
                  dGm_dd_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX1(iNodeX,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX2(iNodeX,3,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                  dGm_dd_dX3(iNodeX,3,iZ2,iZ3,iZ4) )

      CALL ComputeEddingtonTensorComponents_uu &
             ( uD_K(iZ_K), uI1_K(iZ_K), uI2_K(iZ_K), uI3_K(iZ_K), &
               Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K), &
               k_uu )

      S_uu_11 = Half * k_uu(1,1) * uD_K(iZ_K) + uI1_K(iZ_K) * uV1_K(iX_K)
      S_uu_22 = Half * k_uu(2,2) * uD_K(iZ_K) + uI2_K(iZ_K) * uV2_K(iX_K)
      S_uu_33 = Half * k_uu(3,3) * uD_K(iZ_K) + uI3_K(iZ_K) * uV3_K(iX_K)

      dFlux_K(iCR_G1) &
        = + uI1_K(iZ_K) * dV_d_dX1(iNodeX,1,iZ2,iZ3,iZ4) &
          + uI2_K(iZ_K) * dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4) &
          + uI3_K(iZ_K) * dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4) &
          - S_uu_11 * dGm_dd_dX1(iNodeX,1,iZ2,iZ3,iZ4) &
          - S_uu_22 * dGm_dd_dX1(iNodeX,2,iZ2,iZ3,iZ4) &
          - S_uu_33 * dGm_dd_dX1(iNodeX,3,iZ2,iZ3,iZ4)

      dFlux_K(iCR_G2) &
        = + uI1_K(iZ_K) * dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4) &
          + uI2_K(iZ_K) * dV_d_dX2(iNodeX,2,iZ2,iZ3,iZ4) &
          + uI3_K(iZ_K) * dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4) &
          - S_uu_11 * dGm_dd_dX2(iNodeX,1,iZ2,iZ3,iZ4) &
          - S_uu_22 * dGm_dd_dX2(iNodeX,2,iZ2,iZ3,iZ4) &
          - S_uu_33 * dGm_dd_dX2(iNodeX,3,iZ2,iZ3,iZ4)

      dFlux_K(iCR_G3) &
        = + uI1_K(iZ_K) * dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4) &
          + uI2_K(iZ_K) * dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4) &
          + uI3_K(iZ_K) * dV_d_dX3(iNodeX,3,iZ2,iZ3,iZ4) &
          - S_uu_11 * dGm_dd_dX3(iNodeX,1,iZ2,iZ3,iZ4) &
          - S_uu_22 * dGm_dd_dX3(iNodeX,2,iZ2,iZ3,iZ4) &
          - S_uu_33 * dGm_dd_dX3(iNodeX,3,iZ2,iZ3,iZ4)

      Beta &
        = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)  &
          * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
          * SqrtGm_K(iX_K)

      DO iCR = iCR_G1, iCR_G3

        dU_Z(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_Z(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
            - Beta * ( Flux_K(iCR) + dFlux_K(iCR) )

      END DO

    END DO

    CALL TimersStop( Timer_Streaming_Sources )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, dU_Z, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_Z(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_ObserverCorrections

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: xZ1, dZ1, dZ2, dZ3, dZ4, &
    !$OMP               iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( xZ1, dZ1, dZ2, dZ3, dZ4, &
    !$ACC         iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_ObserverCorrections


  SUBROUTINE InitializeIncrement_TwoMoment_Explicit( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    nZ         = iZ_E0 - iZ_B0 + 1           ! Number of Elements per Phase Space Dimension
    nK_X       = PRODUCT( nZ(2:4) )          ! Number of Elements in Position Space
    nK_Z       = nSpecies * nZ(1) * nK_X     ! Number of Elements in Phase Space
    nNodesX_K  = nDOFX * nK_X                ! Number of Nodes in Elements in Position Space
    nNodesZ_K  = nDOFZ * nK_Z                ! Number of Nodes in Elements in Phase Space

    nZ_E  = nZ + [1,0,0,0]                   ! Number of E  Faces per Phase Space Dimension
    nZ_X1 = nZ + [0,1,0,0]                   ! Number of X1 Faces per Phase Space Dimension
    nZ_X2 = nZ + [0,0,1,0]                   ! Number of X2 Faces per Phase Space Dimension
    nZ_X3 = nZ + [0,0,0,1]                   ! Number of X3 Faces per Phase Space Dimension

    nE_X  = PRODUCT( nZ_E (2:4) )            ! Number of E  Faces in Position Space
    nX1_X = PRODUCT( nZ_X1(2:4) )            ! Number of X1 Faces in Position Space
    nX2_X = PRODUCT( nZ_X2(2:4) )            ! Number of X2 Faces in Position Space
    nX3_X = PRODUCT( nZ_X3(2:4) )            ! Number of X3 Faces in Position Space

    nE_Z  = nSpecies * nZ_E (1) * nE_X       ! Number of E  Faces in Phase Space
    nX1_Z = nSpecies * nZ_X1(1) * nX1_X      ! Number of X1 Faces in Phase Space
    nX2_Z = nSpecies * nZ_X2(1) * nX2_X      ! Number of X2 Faces in Phase Space
    nX3_Z = nSpecies * nZ_X3(1) * nX3_X      ! Number of X3 Faces in Phase Space

    !nNodesX_E  = nDOFX_E  * nE_X             ! Number of Nodes on X1 Faces in Position space
    nNodesX_X1 = nDOFX_X1 * nX1_X            ! Number of Nodes on X1 Faces in Position space
    nNodesX_X2 = nDOFX_X2 * nX2_X            ! Number of Nodes on X2 Faces in Position space
    nNodesX_X3 = nDOFX_X3 * nX3_X            ! Number of Nodes on X3 Faces in Position space

    nNodesZ_E  = nDOF_E  * nE_Z              ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X1 = nDOF_X1 * nX1_Z             ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X2 = nDOF_X2 * nX2_Z             ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X3 = nDOF_X3 * nX3_Z             ! Number of Nodes on X1 Faces in Phase space

    ALLOCATE( uV1_K(nNodesX_K) )
    ALLOCATE( uV2_K(nNodesX_K) )
    ALLOCATE( uV3_K(nNodesX_K) )

    ALLOCATE( uD_K (nNodesZ_K) )
    ALLOCATE( uI1_K(nNodesZ_K) )
    ALLOCATE( uI2_K(nNodesZ_K) )
    ALLOCATE( uI3_K(nNodesZ_K) )

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: MeshE % Width, MeshX(1) % Width, MeshX(2) % Width, MeshX(3) % Width, &
    !$OMP          nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3 ) &
    !$OMP MAP( alloc: uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( MeshE % Width, MeshX(1) % Width, MeshX(2) % Width, MeshX(3) % Width, &
    !$ACC         nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3 ) &
    !$ACC CREATE( uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE InitializeIncrement_TwoMoment_Explicit


  SUBROUTINE FinalizeIncrement_TwoMoment_Explicit

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: MeshE % Width, MeshX(1) % Width, MeshX(2) % Width, MeshX(3) % Width, &
    !$OMP               nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3, &
    !$OMP               uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( MeshE % Width, MeshX(1) % Width, MeshX(2) % Width, MeshX(3) % Width, &
    !$ACC         nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3,  &
    !$ACC         uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#endif

    END ASSOCIATE ! dZ1, etc.

    DEALLOCATE( uV1_K, uV2_K, uV3_K )
    DEALLOCATE( uD_K, uI1_K, uI2_K, uI3_K )

  END SUBROUTINE FinalizeIncrement_TwoMoment_Explicit


  SUBROUTINE InitializeIncrement_Divergence_X &
    ( iZP_B0, iZP_E0, nDOFX_X, nDOFZ_X )

    INTEGER, INTENT(in) :: iZP_B0(4), iZP_E0(4) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X ! nDOFX_X1, ...
    INTEGER, INTENT(in) :: nDOFZ_X ! nDOFZ_X1, ...

    INTEGER :: iZP1, iZP2, iZP3, iZP4, iS, iNodeE
    INTEGER :: iX_K, iZ_K, iNodeX, iNodeZ
    INTEGER :: iX_F, iZ_F, iNodeX_X, iNodeZ_X
    INTEGER :: nZP(4), nZP_X(4), nX_X, nX_Z, nNodesX_X, nNodesZ_X

    nZP = iZP_E0 - iZP_B0 + 1
    nZP_X = nZP + [0,0,0,1]
    nX_X = PRODUCT( nZP_X(2:4) )
    nX_Z = nSpecies * nZP_X(1) * nX_X
    nNodesX_X = nDOFX_X * nX_X
    nNodesZ_X = nDOFZ_X * nX_Z

    ! --- Geometry Fields ---

    ALLOCATE( uGF_K (1:nDOFX, &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    iZP_B0(4)-1:iZP_E0(4)+1, &
                    1:nGF) )
    ALLOCATE( uGF_F (1:nDOFX_X, &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    iZP_B0(4)  :iZP_E0(4)+1, &
                    1:nGF) )

    Gm_dd_11_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_SqrtGm  )

    Gm_dd_11_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_SqrtGm  )

    ! --- Conserved Fluid Fields ---

    ALLOCATE( uCF_K(1:nDOFX, &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    iZP_B0(4)-1:iZP_E0(4)+1, &
                    1:nCF) )
    ALLOCATE( uCF_L(1:nDOFX_X, &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    iZP_B0(4)  :iZP_E0(4)+1, &
                    1:nCF) )
    ALLOCATE( uCF_R(1:nDOFX_X, &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    iZP_B0(4)  :iZP_E0(4)+1, &
                    1:nCF) )

    uFD_K(1:nNodesX_K) => uCF_K(:,:,:,iZP_B0(4):iZP_E0(4),iCF_D )
    uS1_K(1:nNodesX_K) => uCF_K(:,:,:,iZP_B0(4):iZP_E0(4),iCF_S1)
    uS2_K(1:nNodesX_K) => uCF_K(:,:,:,iZP_B0(4):iZP_E0(4),iCF_S2)
    uS3_K(1:nNodesX_K) => uCF_K(:,:,:,iZP_B0(4):iZP_E0(4),iCF_S3)

    uFD_L(1:nNodesX_X) => uCF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_D )
    uS1_L(1:nNodesX_X) => uCF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S1)
    uS2_L(1:nNodesX_X) => uCF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S2)
    uS3_L(1:nNodesX_X) => uCF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S3)

    uFD_R(1:nNodesX_X) => uCF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_D )
    uS1_R(1:nNodesX_X) => uCF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S1)
    uS2_R(1:nNodesX_X) => uCF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S2)
    uS3_R(1:nNodesX_X) => uCF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iCF_S3)

    ! --- Conserved Radiation Fields ---

    ALLOCATE( uCR_K(1:nDOFZ, &
                    iZP_B0(1)  :iZP_E0(1)  , &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    1:nSpecies, &
                    iZP_B0(4)-1:iZP_E0(4)+1, &
                    1:nCR) )
    ALLOCATE( uCR_L(1:nDOFZ_X, &
                    iZP_B0(1)  :iZP_E0(1)  , &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    1:nSpecies, &
                    iZP_B0(4)  :iZP_E0(4)+1, &
                    1:nCR) )
    ALLOCATE( uCR_R(1:nDOFZ_X, &
                    iZP_B0(1)  :iZP_E0(1)  , &
                    iZP_B0(2)  :iZP_E0(2)  , &
                    iZP_B0(3)  :iZP_E0(3)  , &
                    1:nSpecies, &
                    iZP_B0(4)  :iZP_E0(4)+1, &
                    1:nCR) )

    uN_K (1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCR_N )
    uG1_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCR_G1)
    uG2_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCR_G2)
    uG3_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCR_G3)

    uN_L (1:nNodesZ_X) => uCR_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_N )
    uG1_L(1:nNodesZ_X) => uCR_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G1)
    uG2_L(1:nNodesZ_X) => uCR_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G2)
    uG3_L(1:nNodesZ_X) => uCR_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G3)

    uN_R (1:nNodesZ_X) => uCR_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_N )
    uG1_R(1:nNodesZ_X) => uCR_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G1)
    uG2_R(1:nNodesZ_X) => uCR_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G2)
    uG3_R(1:nNodesZ_X) => uCR_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCR_G3)

    ! --- Primitive Radiation Fields ---

    ALLOCATE( uD_L (nNodesZ_X) )
    ALLOCATE( uI1_L(nNodesZ_X) )
    ALLOCATE( uI2_L(nNodesZ_X) )
    ALLOCATE( uI3_L(nNodesZ_X) )

    ALLOCATE( uD_R (nNodesZ_X) )
    ALLOCATE( uI1_R(nNodesZ_X) )
    ALLOCATE( uI2_R(nNodesZ_X) )
    ALLOCATE( uI3_R(nNodesZ_X) )

    ! --- Primitive Fluid Velocities ---

    ALLOCATE( uV1_F(nNodesX_X) )
    ALLOCATE( uV2_F(nNodesX_X) )
    ALLOCATE( uV3_F(nNodesX_X) )

    ! --- Fluxes ---

    ALLOCATE( NumericalFlux (nDOFZ_X,nCR, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)+1) )
    ALLOCATE( NumericalFlux2(nDOFZ_X,nCR, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)+1) )
    ALLOCATE( Flux_q        (nDOFZ,nCR, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)  ) )

    ! --- Increment ---

    ALLOCATE( dU_Z          (nDOFZ,nCR, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)  ) )

    ! ---

    ALLOCATE( nIterations_L(nNodesZ_X) )
    ALLOCATE( nIterations_R(nNodesZ_X) )
    ALLOCATE( nIterations_K(nNodesZ_K) )

    ALLOCATE( PositionIndexZ_F(nNodesZ_X) )
    ALLOCATE( PositionIndexZ_K(nNodesZ_K) )

    ALLOCATE( IndexTableZ_F(7,nNodesZ_X) )
    ALLOCATE( IndexTableZ_K(7,nNodesZ_K) )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP             IndexTableZ_F, IndexTableZ_K, &
    !$OMP             uGF_K, uGF_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R, &
    !$OMP             uV1_F, uV2_F, uV3_F, &
    !$OMP             uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP             uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP             NumericalFlux, NumericalFlux2, Flux_q, dU_Z, &
    !$OMP             nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K, &
    !$ACC         uGF_K, uGF_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R, &
    !$ACC         uV1_F, uV2_F, uV3_F, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         NumericalFlux, NumericalFlux2, Flux_q, dU_Z, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP MAP( to: nZP_X ) &
    !$OMP PRIVATE( iNodeZ_X, iZ_F, iX_F )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC COPYIN( nZP_X ) &
    !$ACC PRIVATE( iNodeZ_X, iZ_F, iX_F ) &
    !$ACC PRESENT( PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ_X, iZ_F, iX_F )
#endif
    DO iZP4 = iZP_B0(4), iZP_E0(4)+1
    DO iS  = 1, nSpecies
    DO iZP3 = iZP_B0(3), iZP_E0(3)
    DO iZP2 = iZP_B0(2), iZP_E0(2)
    DO iZP1 = iZP_B0(1), iZP_E0(1)
    DO iNodeX_X = 1, nDOFX_X
    DO iNodeE = 1, nDOFE

      iNodeZ_X = iNodeE &
                 + ( iNodeX_X - 1 ) * nDOFE

      iZ_F = iNodeZ_X &
             + ( iZP1 - iZP_B0(1) ) * nDOFZ_X &
             + ( iZP2 - iZP_B0(2) ) * nDOFZ_X * nZP_X(1) &
             + ( iZP3 - iZP_B0(3) ) * nDOFZ_X * nZP_X(1) * nZP_X(2) &
             + ( iS  - 1          ) * nDOFZ_X * nZP_X(1) * nZP_X(2) * nZP_X(3) &
             + ( iZP4 - iZP_B0(4) ) * nDOFZ_X * nZP_X(1) * nZP_X(2) * nZP_X(3) * nSpecies

      iX_F = iNodeX_X &
             + ( iZP2 - iZP_B0(2) ) * nDOFX_X &
             + ( iZP3 - iZP_B0(3) ) * nDOFX_X * nZP_X(2) &
             + ( iZP4 - iZP_B0(4) ) * nDOFX_X * nZP_X(2) * nZP_X(3)

      PositionIndexZ_F(iZ_F) = iX_F

      IndexTableZ_F(:,iZ_F) = [ iNodeE, iNodeX_X, iZP1, iZP2, iZP3, iS, iZP4 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP MAP( to: nZP ) &
    !$OMP PRIVATE( iNodeZ, iZ_K, iX_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC COPYIN( nZP ) &
    !$ACC PRIVATE( iNodeZ, iZ_K, iX_K ) &
    !$ACC PRESENT( PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iZ_K, iX_K )
#endif
    DO iZP4 = iZP_B0(4), iZP_E0(4)
    DO iS  = 1, nSpecies
    DO iZP3 = iZP_B0(3), iZP_E0(3)
    DO iZP2 = iZP_B0(2), iZP_E0(2)
    DO iZP1 = iZP_B0(1), iZP_E0(1)
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = iNodeE &
               + ( iNodeX - 1 ) * nDOFE

      iZ_K = iNodeZ &
             + ( iZP1 - iZP_B0(1) ) * nDOFZ &
             + ( iZP2 - iZP_B0(2) ) * nDOFZ * nZP(1) &
             + ( iZP3 - iZP_B0(3) ) * nDOFZ * nZP(1) * nZP(2) &
             + ( iS  - 1          ) * nDOFZ * nZP(1) * nZP(2) * nZP(3) &
             + ( iZP4 - iZP_B0(4) ) * nDOFZ * nZP(1) * nZP(2) * nZP(3) * nSpecies

      iX_K = iNodeX &
             + ( iZP2 - iZP_B0(2) ) * nDOFX &
             + ( iZP3 - iZP_B0(3) ) * nDOFX * nZP(2) &
             + ( iZP4 - iZP_B0(4) ) * nDOFX * nZP(2) * nZP(3)

      PositionIndexZ_K(iZ_K) = iX_K

      IndexTableZ_K(:,iZ_K) = [ iNodeE, iNodeX, iZP1, iZP2, iZP3, iS, iZP4 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! AH: This update *may* not be necessary.
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP FROM( PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP       IndexTableZ_F, IndexTableZ_K )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC HOST( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC       IndexTableZ_F, IndexTableZ_K )
#endif

    RETURN
  END SUBROUTINE InitializeIncrement_Divergence_X


  SUBROUTINE FinalizeIncrement_Divergence_X

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP               IndexTableZ_F, IndexTableZ_K, &
    !$OMP               uGF_K, uGF_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R, &
    !$OMP               uV1_F, uV2_F, uV3_F, &
    !$OMP               uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP               uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP               NumericalFlux, NumericalFlux2, Flux_q, dU_Z, &
    !$OMP               nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K, &
    !$ACC         uGF_K, uGF_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R, &
    !$ACC         uV1_F, uV2_F, uV3_F, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         NumericalFlux, NumericalFlux2, Flux_q, dU_Z, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif

    DEALLOCATE( PositionIndexZ_F, PositionIndexZ_K )
    DEALLOCATE( IndexTableZ_F, IndexTableZ_K )
    DEALLOCATE( uGF_K, uGF_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R )
    DEALLOCATE( uD_L, uI1_L, uI2_L, uI3_L )
    DEALLOCATE( uD_R, uI1_R, uI2_R, uI3_R )
    DEALLOCATE( uV1_F, uV2_F, uV3_F )
    DEALLOCATE( NumericalFlux, NumericalFlux2, Flux_q, dU_Z )
    DEALLOCATE( nIterations_L, nIterations_R, nIterations_K )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY( uFD_K, uS1_K, uS2_K, uS3_K )
    NULLIFY( uFD_L, uS1_L, uS2_L, uS3_L )
    NULLIFY( uFD_R, uS1_R, uS2_R, uS3_R )
    NULLIFY( uN_K, uG1_K, uG2_K, uG3_K )
    NULLIFY( uN_L, uG1_L, uG2_L, uG3_L )
    NULLIFY( uN_R, uG1_R, uG2_R, uG3_R )

  END SUBROUTINE FinalizeIncrement_Divergence_X


  SUBROUTINE InitializeIncrement_ObserverCorrections &
    ( iZ_B0, iZ_E0 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iNodeE
    INTEGER :: iX_K, iZ_K, iNodeX, iNodeZ
    INTEGER :: iX_F, iZ_F, iNode_E

    ! --- Geometry Fields ---

    ALLOCATE( uGF_K(1:nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nGF) )

    Gm_dd_11_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_SqrtGm  )

    ! --- Conserved Fluid Fields ---

    ALLOCATE( uCF_K(1:nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nCF) )

    uFD_K(1:nNodesX_K) => uCF_K(:,:,:,iZ_B0(4):iZ_E0(4),iCF_D )
    uS1_K(1:nNodesX_K) => uCF_K(:,:,:,iZ_B0(4):iZ_E0(4),iCF_S1)
    uS2_K(1:nNodesX_K) => uCF_K(:,:,:,iZ_B0(4):iZ_E0(4),iCF_S2)
    uS3_K(1:nNodesX_K) => uCF_K(:,:,:,iZ_B0(4):iZ_E0(4),iCF_S3)

    ! --- Conserved Radiation Fields ---

    ALLOCATE( uCR_K(1:nDOFZ, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)-1:iZ_E0(1)+1, &
                    1:nCR) )
    ALLOCATE( uCR_L(1:nDOF_E, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)  :iZ_E0(1)+1, &
                    1:nCR) )
    ALLOCATE( uCR_R(1:nDOF_E, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)  :iZ_E0(1)+1, &
                    1:nCR) )

    uN_K (1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCR_N )
    uG1_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCR_G1)
    uG2_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCR_G2)
    uG3_K(1:nNodesZ_K) => uCR_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCR_G3)

    uN_L (1:nNodesZ_E) => uCR_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_N )
    uG1_L(1:nNodesZ_E) => uCR_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G1)
    uG2_L(1:nNodesZ_E) => uCR_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G2)
    uG3_L(1:nNodesZ_E) => uCR_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G3)

    uN_R (1:nNodesZ_E) => uCR_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_N )
    uG1_R(1:nNodesZ_E) => uCR_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G1)
    uG2_R(1:nNodesZ_E) => uCR_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G2)
    uG3_R(1:nNodesZ_E) => uCR_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCR_G3)

    ! --- Primitive Radiation Fields ---

    ALLOCATE( uD_L (nNodesZ_E) )
    ALLOCATE( uI1_L(nNodesZ_E) )
    ALLOCATE( uI2_L(nNodesZ_E) )
    ALLOCATE( uI3_L(nNodesZ_E) )

    ALLOCATE( uD_R (nNodesZ_E) )
    ALLOCATE( uI1_R(nNodesZ_E) )
    ALLOCATE( uI2_R(nNodesZ_E) )
    ALLOCATE( uI3_R(nNodesZ_E) )

    ! --- Fluxes ---

    ALLOCATE( NumericalFlux (nDOF_E,nCR, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)+1) )
    ALLOCATE( NumericalFlux2(nDOF_E,nCR, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)+1) )
    ALLOCATE( Flux_q        (nDOFZ,nCR, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)  ) )

    ! --- Increment ---

    ALLOCATE( dU_Z          (nDOFZ,nCR, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)  ) )

    ! --- Eigenvalues ---

    ALLOCATE( Alpha(nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  ) )

    ! --- Derivatives ---

    ALLOCATE( dV_u_dX1   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_u_dX2   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_u_dX3   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )

    ALLOCATE( dV_d_dX1   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_d_dX2   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_d_dX3   (nDOFX,3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )

    ALLOCATE( dGm_dd_dX1(nDOFX,3, &
                         iZ_B0(2)  :iZ_E0(2)  , &
                         iZ_B0(3)  :iZ_E0(3)  , &
                         iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dGm_dd_dX2(nDOFX,3, &
                         iZ_B0(2)  :iZ_E0(2)  , &
                         iZ_B0(3)  :iZ_E0(3)  , &
                         iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dGm_dd_dX3(nDOFX,3, &
                         iZ_B0(2)  :iZ_E0(2)  , &
                         iZ_B0(3)  :iZ_E0(3)  , &
                         iZ_B0(4)  :iZ_E0(4)  ) )

    ! ---

    ALLOCATE( nIterations_L(nNodesZ_E) )
    ALLOCATE( nIterations_R(nNodesZ_E) )
    ALLOCATE( nIterations_K(nNodesZ_K) )

    ALLOCATE( PositionIndexZ_F(nNodesZ_E) )
    ALLOCATE( PositionIndexZ_K(nNodesZ_K) )

    ALLOCATE( IndexTableZ_F(7,nNodesZ_E) )
    ALLOCATE( IndexTableZ_K(7,nNodesZ_K) )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP             IndexTableZ_F, IndexTableZ_K, &
    !$OMP             uGF_K, uCF_K, uCR_K, uCR_L, uCR_R, &
    !$OMP             uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP             uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP             NumericalFlux, NumericalFlux2, Flux_q, dU_Z, Alpha, &
    !$OMP             dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    !$OMP             dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    !$OMP             dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$OMP             nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K, &
    !$ACC         uGF_K, uCF_K, uCR_K, uCR_L, uCR_R, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         NumericalFlux, NumericalFlux2, Flux_q, dU_Z, Alpha, &
    !$ACC         dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    !$ACC         dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    !$ACC         dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif


#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iZ_F, iX_F )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iZ_F, iX_F ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, nZ_E, PositionIndexZ_F, IndexTableZ_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iZ_F, iX_F )
#endif
    DO iZ1 = iZ_B0(1), iZ_E0(1)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iNode_E = 1, nDOF_E

      iZ_F = iNode_E &
             + ( iZ2 - iZ_B0(2) ) * nDOF_E &
             + ( iZ3 - iZ_B0(3) ) * nDOF_E * nZ_E(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOF_E * nZ_E(2) * nZ_E(3) &
             + ( iS  - 1        ) * nDOF_E * nZ_E(2) * nZ_E(3) * nZ_E(4) &
             + ( iZ1 - iZ_B0(1) ) * nDOF_E * nZ_E(2) * nZ_E(3) * nZ_E(4) * nSpecies

      iX_F = iNode_E &
             + ( iZ2 - iZ_B0(2) ) * nDOF_E &
             + ( iZ3 - iZ_B0(3) ) * nDOF_E * nZ_E(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOF_E * nZ_E(2) * nZ_E(3)

      PositionIndexZ_F(iZ_F) = iX_F

      IndexTableZ_F(:,iZ_F) = [ 1, iNode_E, iZ2, iZ3, iZ4, iS, iZ1 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iZ_K, iX_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeZ, iZ_K, iX_K ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, nZ, PositionIndexZ_K, IndexTableZ_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iZ_K, iX_K )
#endif
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = iNodeE &
               + ( iNodeX - 1 ) * nDOFE

      iZ_K = iNodeZ &
             + ( iZ2 - iZ_B0(2) ) * nDOFZ &
             + ( iZ3 - iZ_B0(3) ) * nDOFZ * nZ(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOFZ * nZ(2) * nZ(3) &
             + ( iS  - 1        ) * nDOFZ * nZ(2) * nZ(3) * nZ(4) &
             + ( iZ1 - iZ_B0(1) ) * nDOFZ * nZ(2) * nZ(3) * nZ(4) * nSpecies

      iX_K = iNodeX &
             + ( iZ2 - iZ_B0(2) ) * nDOFX &
             + ( iZ3 - iZ_B0(3) ) * nDOFX * nZ(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOFX * nZ(2) * nZ(3)

      PositionIndexZ_K(iZ_K) = iX_K

      IndexTableZ_K(:,iZ_K) = [ iNodeE, iNodeX, iZ2, iZ3, iZ4, iS, iZ1 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! AH: This update *may* not be necessary.
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP FROM( PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP       IndexTableZ_F, IndexTableZ_K )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC HOST( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC       IndexTableZ_F, IndexTableZ_K )
#endif

    RETURN
  END SUBROUTINE InitializeIncrement_ObserverCorrections


  SUBROUTINE FinalizeIncrement_ObserverCorrections

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP               IndexTableZ_F, IndexTableZ_K, &
    !$OMP               uGF_K, uCF_K, uCR_K, uCR_L, uCR_R, &
    !$OMP               uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP               uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP               NumericalFlux, NumericalFlux2, Flux_q, dU_Z, Alpha, &
    !$OMP               dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    !$OMP               dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    !$OMP               dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$OMP               nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K, &
    !$ACC         uGF_K, uCF_K, uCR_K, uCR_L, uCR_R, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         NumericalFlux, NumericalFlux2, Flux_q, dU_Z, Alpha, &
    !$ACC         dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    !$ACC         dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    !$ACC         dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif

    DEALLOCATE( PositionIndexZ_F, PositionIndexZ_K )
    DEALLOCATE( IndexTableZ_F, IndexTableZ_K )
    DEALLOCATE( uGF_K, uCF_K, uCR_K, uCR_L, uCR_R )
    DEALLOCATE( uD_L, uI1_L, uI2_L, uI3_L )
    DEALLOCATE( uD_R, uI1_R, uI2_R, uI3_R )
    DEALLOCATE( NumericalFlux, NumericalFlux2, Flux_q, dU_Z, Alpha )
    DEALLOCATE( dV_u_dX1, dV_u_dX2, dV_u_dX3 )
    DEALLOCATE( dV_d_dX1, dV_d_dX2, dV_d_dX3 )
    DEALLOCATE( dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3 )
    DEALLOCATE( nIterations_L, nIterations_R, nIterations_K )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( uFD_K, uS1_K, uS2_K, uS3_K )
    NULLIFY( uN_K, uG1_K, uG2_K, uG3_K )
    NULLIFY( uN_L, uG1_L, uG2_L, uG3_L )
    NULLIFY( uN_R, uG1_R, uG2_R, uG3_R )

  END SUBROUTINE FinalizeIncrement_ObserverCorrections


  SUBROUTINE ComputeOffGridFlux( iZP_B0, iZP_E0, nDOFZ_X, NumericalFlux, NumericalFlux2 )

    INTEGER, INTENT(in) :: iZP_B0(4), iZP_E0(4) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFZ_X ! nDOFZ_X1, ...

    REAL(DP), INTENT(in) :: &
      NumericalFlux (nDOFZ_X,nCR, &
                     iZP_B0(1):iZP_E0(1)  , &
                     iZP_B0(2):iZP_E0(2)  , &
                     iZP_B0(3):iZP_E0(3)  , &
                     nSpecies, &
                     iZP_B0(4):iZP_E0(4)+1)
    REAL(DP), INTENT(in) :: &
      NumericalFlux2(nDOFZ_X,nCR, &
                     iZP_B0(1):iZP_E0(1)  , &
                     iZP_B0(2):iZP_E0(2)  , &
                     iZP_B0(3):iZP_E0(3)  , &
                     nSpecies, &
                     iZP_B0(4):iZP_E0(4)+1)

    INTEGER :: iZP1, iZP2, iZP3, iZP4, iS, iCR
    INTEGER :: iNodeZ_X

    REAL(DP) :: FluxIn1, FluxIn2, FluxOt1, FluxOt2 

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( FluxIn1, FluxIn2, FluxOt1, FluxOt2  ) &
    !$OMP MAP( tofrom: OffGridFlux_TwoMoment )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( FluxIn1, FluxIn2, FluxOt1, FluxOt2 ) &
    !$ACC COPY( OffGridFlux_TwoMoment ) &
    !$ACC PRESENT( iZP_B0, iZP_E0, LeptonNumber, NumericalFlux, NumericalFlux2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( FluxIn1, FluxIn2, FluxOt1, FluxOt2  )
#endif
    DO iCR = 1, nCR

      FluxIn1 = Zero
      FluxIn2 = Zero

#if   defined( THORNADO_OMP_OL )
      !$OMP PARALLEL DO SIMD COLLAPSE(5) &
      !$OMP REDUCTION( + : FluxIn1, FluxIn2 )
#elif defined( THORNADO_OACC   )
      !$ACC LOOP VECTOR COLLAPSE(5) &
      !$ACC REDUCTION( + : FluxIn1, FluxIn2 )
#endif
      DO iS  = 1       , nSpecies
      DO iZP3 = iZP_B0(3), iZP_E0(3)
      DO iZP2 = iZP_B0(2), iZP_E0(2)
      DO iZP1 = iZP_B0(1), iZP_E0(1)

        DO iNodeZ_X = 1, nDOFZ_X

          FluxIn1 = FluxIn1 + LeptonNumber(iS) * NumericalFlux (iNodeZ_X,iCR,iZP1,iZP2,iZP3,iS,iZP_B0(4))
          FluxIn2 = FluxIn2 +                    NumericalFlux2(iNodeZ_X,iCR,iZP1,iZP2,iZP3,iS,iZP_B0(4))

        END DO

      END DO
      END DO
      END DO
      END DO

      FluxOt1 = Zero
      FluxOt2 = Zero

#if   defined( THORNADO_OMP_OL )
      !$OMP PARALLEL DO SIMD COLLAPSE(5) &
      !$OMP REDUCTION( + : FluxOt1, FluxOt2 )
#elif defined( THORNADO_OACC   )
      !$ACC LOOP VECTOR COLLAPSE(5) &
      !$ACC REDUCTION( + : FluxOt1, FluxOt2 )
#endif
      DO iS  = 1       , nSpecies
      DO iZP3 = iZP_B0(3), iZP_E0(3)
      DO iZP2 = iZP_B0(2), iZP_E0(2)
      DO iZP1 = iZP_B0(1), iZP_E0(1)

        DO iNodeZ_X = 1, nDOFZ_X

          FluxOt1 = FluxOt1 + LeptonNumber(iS) * NumericalFlux (iNodeZ_X,iCR,iZP1,iZP2,iZP3,iS,iZP_E0(4)+1)
          FluxOt2 = FluxOt2 +                    NumericalFlux2(iNodeZ_X,iCR,iZP1,iZP2,iZP3,iS,iZP_E0(4)+1)

        END DO

      END DO
      END DO
      END DO
      END DO

      OffGridFlux_TwoMoment(iCR)     = OffGridFlux_TwoMoment(iCR)     + FluxIn1 - FluxOt1
      OffGridFlux_TwoMoment(iCR+nCR) = OffGridFlux_TwoMoment(iCR+nCR) + FluxIn2 - FluxOt2

    END DO

  END SUBROUTINE ComputeOffGridFlux


END MODULE TwoMoment_DiscretizationModule_Streaming
