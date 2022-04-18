#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EULER
#endif

MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Third, &
    Half, &
    One, &
    Two, &
    Three, &
    Four
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    dLXdX3_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    iGF_Phi_N, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nDF, &
    iDF_Sh_X1, &
    iDF_Sh_X2, &
    iDF_Sh_X3
  USE GeometryComputationModule_XCFC, ONLY: &
    ComputeChristoffelSymbols_3D_XCFC
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler, &
    Eigenvalues_Euler, &
    AlphaMiddle_Euler, &
    Flux_X1_Euler, &
    Flux_X2_Euler, &
    Flux_X3_Euler, &
    StressTensor_Diagonal_Euler, &
    NumericalFlux_Euler_X1, &
    NumericalFlux_Euler_X2, &
    NumericalFlux_Euler_X3
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE TimersModule_Euler,  ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_DG, &
    Timer_Euler_Divergence, &
    Timer_Euler_Geometry, &
    Timer_Euler_Gravity, &
    Timer_Euler_SurfaceTerm, &
    Timer_Euler_VolumeTerm, &
    Timer_Euler_Increment, &
    Timer_Euler_DG_CopyIn, &
    Timer_Euler_DG_Permute, &
    Timer_Euler_DG_Interpolate, &
    Timer_Euler_DG_ComputePrimitive, &
    Timer_Euler_DG_CopyOut, &
    Timer_Euler_DG_ErrorCheck
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_DG_Explicit

  REAL(DP), POINTER, CONTIGUOUS :: &
    Gm_dd_11_K(:), Gm_dd_22_K(:), Gm_dd_33_K(:), SqrtGm_K(:), &
    Gm_dd_11_F(:), Gm_dd_22_F(:), Gm_dd_33_F(:), SqrtGm_F(:), &
    Beta_1_K  (:), Beta_2_K  (:), Beta_3_K  (:), Alpha_K (:), &
    Beta_1_F  (:), Beta_2_F  (:), Beta_3_F  (:), Alpha_F (:)

  REAL(DP), POINTER, CONTIGUOUS :: &
    uD_K(:), uS1_K(:), uS2_K(:), uS3_K(:), uE_K(:), uNe_K(:), &
    uD_L(:), uS1_L(:), uS2_L(:), uS3_L(:), uE_L(:), uNe_L(:), &
    uD_R(:), uS1_R(:), uS2_R(:), uS3_R(:), uE_R(:), uNe_R(:)

  REAL(DP), ALLOCATABLE :: &
    pD_K(:), pV1_K(:), pV2_K(:), pV3_K(:), pE_K(:), pNe_K(:), &
    pD_L(:), pV1_L(:), pV2_L(:), pV3_L(:), pE_L(:), pNe_L(:), &
    pD_R(:), pV1_R(:), pV2_R(:), pV3_R(:), pE_R(:), pNe_R(:)

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableX_F(:,:), IndexTableX_V(:,:)

  INTEGER :: nX   (3), nX_K,  nNodesX_K
  INTEGER :: nX_X1(3), nX1_X, nNodesX_X1
  INTEGER :: nX_X2(3), nX2_X, nNodesX_X2
  INTEGER :: nX_X3(3), nX3_X, nNodesX_X3

  REAL(DP), PUBLIC :: OffGridFlux_Euler_X1_Inner(nCF), &
                      OffGridFlux_Euler_X1_Outer(nCF), &
                      OffGridFlux_Euler_X2_Inner(nCF), &
                      OffGridFlux_Euler_X2_Outer(nCF), &
                      OffGridFlux_Euler_X3_Inner(nCF), &
                      OffGridFlux_Euler_X3_Outer(nCF)


CONTAINS


  SUBROUTINE ComputeIncrement_Euler_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
      SuppressBC_Option, UseXCFC_Option, &
      SurfaceFlux_X1_Option, &
      SurfaceFlux_X2_Option, &
      SurfaceFlux_X3_Option )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)         :: &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      D (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)           :: &
      dU(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      SuppressBC_Option
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      UseXCFC_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      SurfaceFlux_X1_Option(:,:,:,:,:), &
      SurfaceFlux_X2_Option(:,:,:,:,:), &
      SurfaceFlux_X3_Option(:,:,:,:,:)

    ! --- Surface flux for coarse/fine corrections ---

    REAL(DP) :: &
      SurfaceFlux_X1(nDOFX_X1, &
          iX_B0(1):iX_E0(1)+1, &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3), &
          nCF)

    REAL(DP) :: &
      SurfaceFlux_X2(nDOFX_X2, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2)+1, &
          iX_B0(3):iX_E0(3), &
          nCF)

    REAL(DP) :: &
      SurfaceFlux_X3(nDOFX_X3, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3)+1, &
          nCF)

    INTEGER  :: iNX, iX1, iX2, iX3, iCF
    LOGICAL  :: SuppressBC, UseXCFC
    REAL(DP) :: tau(nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3))

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG )

    UseXCFC = .FALSE.
    IF( PRESENT( UseXCFC_Option ) ) &
      UseXCFC = UseXCFC_Option

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3 ) &
    !$OMP MAP( alloc: dU, tau )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3 ) &
    !$ACC CREATE(     dU, tau )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    CALL TimersStop_Euler( Timer_Euler_DG )

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

      CALL DetectShocks_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    END IF

    CALL TimersStart_Euler( Timer_Euler_DG )

    CALL InitializeIncrement_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1 )

    CALL TimersStart_Euler( Timer_Euler_Increment )

    IF( UseXCFC )THEN

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B1, iX_E1, tau, G )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        tau(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

      END DO
      END DO
      END DO
      END DO

    ELSE

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B1, iX_E1, tau )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        tau(iNX,iX1,iX2,iX3) = One

      END DO
      END DO
      END DO
      END DO

    END IF

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, dU )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCF) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    OffGridFlux_Euler_X1_Inner = Zero
    OffGridFlux_Euler_X1_Outer = Zero
    OffGridFlux_Euler_X2_Inner = Zero
    OffGridFlux_Euler_X2_Outer = Zero
    OffGridFlux_Euler_X3_Inner = Zero
    OffGridFlux_Euler_X3_Outer = Zero

    CALL TimersStop_Euler( Timer_Euler_Increment )

    CALL TimersStart_Euler( Timer_Euler_Divergence )

    CALL ComputeIncrement_Euler_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X1 )

    CALL ComputeIncrement_Euler_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X2 )

    CALL ComputeIncrement_Euler_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X3 )

    IF( PRESENT( SurfaceFlux_X1_Option ) ) &
      SurfaceFlux_X1_Option = SurfaceFlux_X1

    IF( PRESENT( SurfaceFlux_X2_Option ) ) &
      SurfaceFlux_X2_Option = SurfaceFlux_X2

    IF( PRESENT( SurfaceFlux_X3_Option ) ) &
      SurfaceFlux_X3_Option = SurfaceFlux_X3

    CALL TimersStop_Euler( Timer_Euler_Divergence )

    ! --- Multiply Inverse Mass Matrix ---

    CALL TimersStart_Euler( Timer_Euler_Increment )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, dX1, dX2, dX3, dU, G, tau, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            * tau(iNX,iX1,iX2,iX3) &
            / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * dX1(iX1) * dX2(iX2) * dX3(iX3) )

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Increment )

    CALL TimersStart_Euler( Timer_Euler_Geometry )

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, UseXCFC, dU )

    CALL TimersStop_Euler( Timer_Euler_Geometry )

    CALL TimersStart_Euler( Timer_Euler_Gravity )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Gravity )

    CALL FinalizeIncrement_Euler

#ifdef THORNADO_DEBUG_EULER
#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET UPDATE FROM( dU )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC UPDATE HOST       ( dU )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU)', MAXLOC(dU)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU)', MAXVAL(dU)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    dU, U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G, tau )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      dU, U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G, tau )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_DG )

  END SUBROUTINE ComputeIncrement_Euler_DG_Explicit


  SUBROUTINE ComputeIncrement_Euler_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X1 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X1(1:nDOFX_X1,iX_B0(1):iX_E0(1)+1, &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3),1:nCF)

    INTEGER  :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCF, iGF
    INTEGER  :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: Flux_L   (nCF), Flux_R   (nCF)
    REAL(DP) :: Flux_F   (nCF), Flux_K   (nCF)
    REAL(DP) :: uCF_L_nCF(nCF), uCF_R_nCF(nCF)

    INTEGER  :: iErr(nNodesX_X1)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)-1:iX_E0(1)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)  :iX_E0(1)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCF_K(nDOFX, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)-1:iX_E0(1)+1, &
            nCF)
    REAL(DP) :: &
      uCF_L(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCF)
    REAL(DP) :: &
      uCF_R(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCF)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDF_L(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)
    REAL(DP) :: &
      uDF_R(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X1,nCF, &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3), &
                    iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCF, &
             iX_B0(2):iX_E0(2), &
             iX_B0(3):iX_E0(3), &
             iX_B0(1):iX_E0(1))

    ! --- X1 Increment ---

    REAL(DP) :: &
      dU_X1(nDOFX,nCF, &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1))

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(2) ; iXP_E0(1) = iX_E0(2)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(1) ; iXP_E0(3) = iX_E0(1)

    ASSOCIATE( dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX2, dX3 ) &
    !$OMP MAP( alloc: EigVals_L, EigVals_R, &
    !$OMP             Flux_L, Flux_R, &
    !$OMP             Flux_F, Flux_K, &
    !$OMP             uCF_L_nCF, uCF_R_nCF, &
    !$OMP             iErr, &
    !$OMP             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP             NumericalFlux, Flux_q, dU_X1, SurfaceFlux_X1 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX2, dX3 ) &
    !$ACC CREATE(     EigVals_L, EigVals_R, &
    !$ACC             Flux_L, Flux_R, &
    !$ACC             Flux_F, Flux_K, &
    !$ACC             uCF_L_nCF, uCF_R_nCF, &
    !$ACC             iErr, &
    !$ACC             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC             NumericalFlux, Flux_q, dU_X1, SurfaceFlux_X1 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X1, &
             G_K, G_F, uCF_K, uCF_L, uCF_R )

    ! --- Permute data ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1           , nGF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX2,iX3,iX1,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, uCF_K, U )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1           , nCF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      uCF_K(iNX,iX2,iX3,iX1,iCF) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      uDF_L(iX2,iX3,iX1,1) = D(1,iX1-1,iX2,iX3,iDF_Sh_X2)
      uDF_L(iX2,iX3,iX1,2) = D(1,iX1-1,iX2,iX3,iDF_Sh_X3)
      uDF_R(iX2,iX3,iX1,1) = D(1,iX1  ,iX2,iX3,iDF_Sh_X2)
      uDF_R(iX2,iX3,iX1,2) = D(1,iX1  ,iX2,iX3,iDF_Sh_X3)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )

    ! --- Lapse function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Alpha), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX_X1 )

    ! --- Shift vector (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX_X1 )

    ! --- Compute metric and metric determinant on faces ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX1   = iX_B0(1), iX_E0(1) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iNX_X = 1       , nDOFX_X1

      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_2   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_3   ), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Alpha) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_Alpha), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCF_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCF_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCF), nDOFX_X1 )

    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL ComputePrimitive_Euler &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$ACC          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, uD_L, uS1_L, uE_L, &
    !$ACC          pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, uD_R, uS1_R, uE_R, &
    !$ACC          Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F, &
    !$ACC          Alpha_F, Beta_1_F, &
    !$ACC          iErr, &
    !$ACC          uCF_L, uCF_R, uDF_L, uDF_R, NumericalFlux, &
    !$ACC          dX2, dX3, WeightsX_X1, IndexTableX_F, SurfaceFlux_X1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_X = 1, nNodesX_X1

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_Euler &
            ( pV1_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_11_F(iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X) )

      Flux_L &
        = Flux_X1_Euler &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X) )

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_Euler &
            ( pV1_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_11_F(iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X) )

      Flux_R &
        = Flux_X1_Euler &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X) )

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iErr(iNX_X) = 0

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uD_L          (iNX_X), &
              uS1_L         (iNX_X), &
              uE_L          (iNX_X), &
              Flux_L(iCF_D)        , &
              Flux_L(iCF_S1)       , &
              Flux_L(iCF_E)        , &
              uD_R          (iNX_X), &
              uS1_R         (iNX_X), &
              uE_R          (iNX_X), &
              Flux_R(iCF_D)        , &
              Flux_R(iCF_S1)       , &
              Flux_R(iCF_E)        , &
              Gm_dd_11_F    (iNX_X), &
              AlphaPls             , &
              AlphaMns             , &
              Alpha_F       (iNX_X), &
              Beta_1_F      (iNX_X), &
              iErr          (iNX_X) )

      iNX = IndexTableX_F(1,iNX_X)
      iX2 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX1 = IndexTableX_F(4,iNX_X)

      DO iCF = 1, nCF

        uCF_L_nCF(iCF) = uCF_L(iNX,iX2,iX3,iX1,iCF)
        uCF_R_nCF(iCF) = uCF_R(iNX,iX2,iX3,iX1,iCF)

      END DO

      Flux_F &
        = NumericalFlux_Euler_X1 &
            ( uCF_L_nCF           , &
              uCF_R_nCF           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            , &
              AlphaMdl            , &
              Gm_dd_11_F(iNX_X)   , &
              pV1_L     (iNX_X)   , &
              pV1_R     (iNX_X)   , &
              P_L                 , &
              P_R                 , &
              Alpha_F   (iNX_X)   , &
              Beta_1_F  (iNX_X)   , &
              uDF_L(iX2,iX3,iX1,1), &
              uDF_R(iX2,iX3,iX1,2), &
              uDF_L(iX2,iX3,iX1,1), &
              uDF_R(iX2,iX3,iX1,2) )

      DO iCF = 1, nCF

        SurfaceFlux_X1(iNX,iX1,iX2,iX3,iCF) &
          = Flux_F(iCF) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

        NumericalFlux(iNX,iCF,iX2,iX3,iX1) &
          = Flux_F(iCF) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX2(iX2) * dX3(iX3) &
              * WeightsX_X1(iNX)

      END DO ! iCF

    END DO ! iNX_X

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X1, + One, LX_X1_Dn, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, &
             Zero, dU_X1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X1, - One, LX_X1_Up, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, &
             One,  dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, Alpha_K, Beta_1_K, &
    !$ACC          Flux_q, SqrtGm_K, dX2, dX3, WeightsX_q, IndexTableX_V )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X1_Euler &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_1_K  (iNX_K) )

      iNX = IndexTableX_V(1,iNX_K)
      iX2 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX1 = IndexTableX_V(4,iNX_K)

      DO iCF = 1, nCF

        Flux_q(iNX,iCF,iX2,iX3,iX1) &
          = Flux_K(iCF) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX2(iX2) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX, One, dLXdX1_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, dU_X1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1       , nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X1(iNX,iCF,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

#ifdef THORNADO_DEBUG_EULER
#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET UPDATE FROM( dU_X1 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC UPDATE HOST       ( dU_X1 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X1)', MAXLOC(dU_X1)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X1)', MAXVAL(dU_X1)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    NumericalFlux, iErr, SurfaceFlux_X1 ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX2, dX3, &
    !$OMP               EigVals_L, EigVals_R, &
    !$OMP               Flux_L, Flux_R, &
    !$OMP               Flux_F, Flux_K, &
    !$OMP               uCF_L_nCF, uCF_R_nCF, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP               Flux_q, dU_X1 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      NumericalFlux, iErr, SurfaceFlux_X1 ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX2, dX3, &
    !$ACC               EigVals_L, EigVals_R, &
    !$ACC               Flux_L, Flux_R, &
    !$ACC               Flux_F, Flux_K, &
    !$ACC               uCF_L_nCF, uCF_R_nCF, &
    !$ACC               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC               Flux_q, dU_X1 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iCF   = 1       , nCF
    DO iNX_X = 1       , nDOFX_X1

      OffGridFlux_Euler_X1_Inner(iCF) &
        = OffGridFlux_Euler_X1_Inner(iCF) &
            + NumericalFlux(iNX_X,iCF,iX2,iX3,iX_B0(1))

      OffGridFlux_Euler_X1_Outer(iCF) &
        = OffGridFlux_Euler_X1_Outer(iCF) &
            + NumericalFlux(iNX_X,iCF,iX2,iX3,iX_E0(1)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

#ifdef HYDRO_RELATIVISTIC

    DO iNX_X = 1, nNodesX_X1

      CALL DescribeError_Euler( iErr(iNX_X) )

    END DO

#endif

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X1


  SUBROUTINE ComputeIncrement_Euler_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X2 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X2(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2)+1, &
                                iX_B0(3):iX_E0(3),1:nCF)

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCF, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: Flux_L   (nCF), Flux_R   (nCF)
    REAL(DP) :: Flux_F   (nCF), Flux_K   (nCF)
    REAL(DP) :: uCF_L_nCF(nCF), uCF_R_nCF(nCF)

    INTEGER  :: iErr(nNodesX_X2)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)-1:iX_E0(2)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X2, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)  :iX_E0(2)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCF_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)-1:iX_E0(2)+1, &
            nCF)
    REAL(DP) :: &
      uCF_L(nDOFX_X2, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCF)
    REAL(DP) :: &
      uCF_R(nDOFX_X2, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCF)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDF_L(iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2)+1,2)
    REAL(DP) :: &
      uDF_R(iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X2,nCF, &
                    iX_B0(1):iX_E0(1), &
                    iX_B0(3):iX_E0(3), &
                    iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCF, &
             iX_B0(1):iX_E0(1), &
             iX_B0(3):iX_E0(3), &
             iX_B0(2):iX_E0(2))

    ! --- X2 Increment ---

    REAL(DP) :: &
      dU_X2(nDOFX,nCF, &
            iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2))

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(2) ; iXP_E0(3) = iX_E0(2)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX3 ) &
    !$OMP MAP( alloc: EigVals_L, EigVals_R, &
    !$OMP             Flux_L, Flux_R, &
    !$OMP             Flux_F, Flux_K, &
    !$OMP             uCF_L_nCF, uCF_R_nCF, &
    !$OMP             iErr, &
    !$OMP             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP             NumericalFlux, Flux_q, dU_X2, SurfaceFlux_X2 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX3 ) &
    !$ACC CREATE(     EigVals_L, EigVals_R, &
    !$ACC             Flux_L, Flux_R, &
    !$ACC             Flux_F, Flux_K, &
    !$ACC             uCF_L_nCF, uCF_R_nCF, &
    !$ACC             iErr, &
    !$ACC             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC             NumericalFlux, Flux_q, dU_X2, SurfaceFlux_X2 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X2, &
             G_K, G_F, uCF_K, uCF_L, uCF_R )

    ! --- Permute data ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1           , nGF
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX1,iX3,iX2,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, uCF_K, U )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1           , nCF
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCF_K(iNX,iX1,iX3,iX2,iCF) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDF_L(iX1,iX3,iX2,1) = D(1,iX1,iX2-1,iX3,iDF_Sh_X1)
      uDF_L(iX1,iX3,iX2,2) = D(1,iX1,iX2-1,iX3,iDF_Sh_X3)
      uDF_R(iX1,iX3,iX2,1) = D(1,iX1,iX2  ,iX3,iDF_Sh_X1)
      uDF_R(iX1,iX3,iX2,2) = D(1,iX1,iX2  ,iX3,iDF_Sh_X3)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )

    ! --- Lapse function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Alpha), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Alpha), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Alpha), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Alpha), nDOFX_X2 )

    ! --- Shift vector (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Beta_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX_X2 )

    ! --- Compute metric and metric determinant on faces ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX2   = iX_B0(2), iX_E0(2) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iNX_X = 1       , nDOFX_X2

      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_2   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_3   ), SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_Alpha) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_Alpha), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               uCF_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCF), nDOFX_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               uCF_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCF), nDOFX_X2 )

    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL ComputePrimitive_Euler &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$ACC          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, uD_L, uS2_L, uE_L, &
    !$ACC          pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, uD_R, uS2_R, uE_R, &
    !$ACC          Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F, &
    !$ACC          Alpha_F, Beta_2_F, &
    !$ACC          iErr, &
    !$ACC          uCF_L, uCF_R, uDF_L, uDF_R, NumericalFlux, &
    !$ACC          dX1, dX3, WeightsX_X2, IndexTableX_F, SurfaceFlux_X2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_X = 1, nNodesX_X2

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_Euler &
            ( pV2_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_22_F(iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_2_F  (iNX_X) )

      Flux_L &
        = Flux_X2_Euler &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_2_F  (iNX_X) )

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_Euler &
            ( pV2_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_22_F(iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_2_F  (iNX_X) )

      Flux_R &
        = Flux_X2_Euler &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_2_F  (iNX_X) )

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iErr(iNX_X) = 0

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uD_L          (iNX_X), &
              uS2_L         (iNX_X), &
              uE_L          (iNX_X), &
              Flux_L(iCF_D)        , &
              Flux_L(iCF_S2)       , &
              Flux_L(iCF_E)        , &
              uD_R          (iNX_X), &
              uS2_R         (iNX_X), &
              uE_R          (iNX_X), &
              Flux_R(iCF_D)        , &
              Flux_R(iCF_S2)       , &
              Flux_R(iCF_E)        , &
              Gm_dd_22_F    (iNX_X), &
              AlphaPls             , &
              AlphaMns             , &
              Alpha_F       (iNX_X), &
              Beta_2_F      (iNX_X), &
              iErr          (iNX_X) )

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX2 = IndexTableX_F(4,iNX_X)

      DO iCF = 1, nCF

        uCF_L_nCF(iCF) = uCF_L(iNX,iX1,iX3,iX2,iCF)
        uCF_R_nCF(iCF) = uCF_R(iNX,iX1,iX3,iX2,iCF)

      END DO

      Flux_F &
        = NumericalFlux_Euler_X2 &
            ( uCF_L_nCF           , &
              uCF_R_nCF           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            , &
              AlphaMdl            , &
              Gm_dd_22_F(iNX_X)   , &
              pV2_L     (iNX_X)   , &
              pV2_R     (iNX_X)   , &
              P_L                 , &
              P_R                 , &
              Alpha_F   (iNX_X)   , &
              Beta_2_F  (iNX_X)   , &
              uDF_L(iX1,iX3,iX2,1), &
              uDF_R(iX1,iX3,iX2,2), &
              uDF_L(iX1,iX3,iX2,1), &
              uDF_R(iX1,iX3,iX2,2) )

      DO iCF = 1, nCF

        SurfaceFlux_X2(iNX,iX1,iX2,iX3,iCF) &
          = Flux_F(iCF) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

        NumericalFlux(iNX,iCF,iX1,iX3,iX2) &
          = Flux_F(iCF) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX1(iX1) * dX3(iX3) &
              * WeightsX_X2(iNX)

      END DO ! iCF

    END DO ! iNX_X

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X2, + One, LX_X2_Dn, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, &
             Zero, dU_X2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X2, - One, LX_X2_Up, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, &
             One,  dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, Alpha_K, Beta_2_K, &
    !$ACC          Flux_q, SqrtGm_K, dX1, dX3, WeightsX_q, IndexTableX_V )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X2_Euler &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_2_K  (iNX_K) )

      iNX = IndexTableX_V(1,iNX_K)
      iX1 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX2 = IndexTableX_V(4,iNX_K)

      DO iCF = 1, nCF

        Flux_q(iNX,iCF,iX1,iX3,iX2) &
          = Flux_K(iCF) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX1(iX1) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX, One, dLXdX2_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, dU_X2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1       , nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X2(iNX,iCF,iX1,iX3,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

#ifdef THORNADO_DEBUG_EULER
#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET UPDATE FROM( dU_X2 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC UPDATE HOST       ( dU_X2 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X2)', MAXLOC(dU_X2)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X2)', MAXVAL(dU_X2)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    NumericalFlux, iErr, SurfaceFlux_X2 ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX3, &
    !$OMP               EigVals_L, EigVals_R, &
    !$OMP               Flux_L, Flux_R, &
    !$OMP               Flux_F, Flux_K, &
    !$OMP               uCF_L_nCF, uCF_R_nCF, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP               Flux_q, dU_X2 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      NumericalFlux, iErr, SurfaceFlux_X2 ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX3, &
    !$ACC               EigVals_L, EigVals_R, &
    !$ACC               Flux_L, Flux_R, &
    !$ACC               Flux_F, Flux_K, &
    !$ACC               uCF_L_nCF, uCF_R_nCF, &
    !$ACC               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC               Flux_q, dU_X2 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iCF   = 1       , nCF
    DO iNX_X = 1       , nDOFX_X2

      OffGridFlux_Euler_X2_Inner(iCF) &
        = OffGridFlux_Euler_X2_Inner(iCF) &
            + NumericalFlux(iNX_X,iCF,iX1,iX3,iX_B0(2))

      OffGridFlux_Euler_X2_Outer(iCF) &
        = OffGridFlux_Euler_X2_Outer(iCF) &
            + NumericalFlux(iNX_X,iCF,iX1,iX3,iX_E0(2)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

#ifdef HYDRO_RELATIVISTIC

    DO iNX_X = 1, nNodesX_X2

      CALL DescribeError_Euler( iErr(iNX_X) )

    END DO

#endif

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X2


  SUBROUTINE ComputeIncrement_Euler_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X3 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X3(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1,1:nCF)

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCF, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: Flux_L   (nCF), Flux_R   (nCF)
    REAL(DP) :: Flux_F   (nCF), Flux_K   (nCF)
    REAL(DP) :: uCF_L_nCF(nCF), uCF_R_nCF(nCF)

    INTEGER  :: iErr(nNodesX_X3)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)-1:iX_E0(3)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X3, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCF_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)-1:iX_E0(3)+1, &
            nCF)
    REAL(DP) :: &
      uCF_L(nDOFX_X3, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3)+1, &
            nCF)
    REAL(DP) :: &
      uCF_R(nDOFX_X3, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3)+1, &
            nCF)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDF_L(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3)+1,2)
    REAL(DP) :: &
      uDF_R(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X3,nCF, &
                    iX_B0(1):iX_E0(1), &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCF, &
             iX_B0(1):iX_E0(1), &
             iX_B0(2):iX_E0(2), &
             iX_B0(3):iX_E0(3))

    ! --- X3 Increment ---

    REAL(DP) :: &
      dU_X3(nDOFX,nCF, &
            iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3))

    IF( iX_E0(3) .EQ. iX_B0(3) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(2) ; iXP_E0(2) = iX_E0(2)
    iXP_B0(3) = iX_B0(3) ; iXP_E0(3) = iX_E0(3)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX2 ) &
    !$OMP MAP( alloc: EigVals_L, EigVals_R, &
    !$OMP             Flux_L, Flux_R, &
    !$OMP             Flux_F, Flux_K, &
    !$OMP             uCF_L_nCF, uCF_R_nCF, &
    !$OMP             iErr, &
    !$OMP             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP             NumericalFlux, Flux_q, dU_X3, SurfaceFlux_X3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX2 ) &
    !$ACC CREATE(     EigVals_L, EigVals_R, &
    !$ACC             Flux_L, Flux_R, &
    !$ACC             Flux_F, Flux_K, &
    !$ACC             uCF_L_nCF, uCF_R_nCF, &
    !$ACC             iErr, &
    !$ACC             G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC             NumericalFlux, Flux_q, dU_X3, SurfaceFlux_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X3, &
             G_K, G_F, uCF_K, uCF_L, uCF_R )

    ! --- Permute data ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1           , nGF
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX1,iX2,iX3,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, uCF_K, U )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1           , nCF
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCF_K(iNX,iX1,iX2,iX3,iCF) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDF_L(iX1,iX2,iX3,1) = D(1,iX1,iX2,iX3-1,iDF_Sh_X1)
      uDF_L(iX1,iX2,iX3,2) = D(1,iX1,iX2,iX3-1,iDF_Sh_X2)
      uDF_R(iX1,iX2,iX3,1) = D(1,iX1,iX2,iX3  ,iDF_Sh_X1)
      uDF_R(iX1,iX2,iX3,2) = D(1,iX1,iX2,iX3  ,iDF_Sh_X2)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_1), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_1), nDOFX_X3 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_2), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_2), nDOFX_X3 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_3), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_h_3), nDOFX_X3 )

    ! --- Lapse function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_Alpha), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Alpha), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Alpha), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Alpha), nDOFX_X3 )

    ! --- Shift vector (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_Beta_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_3), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_3), nDOFX_X3 )

    ! --- Compute metric and metric determinant on faces ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3   = iX_B0(3), iX_E0(3) + 1
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iNX_X = 1       , nDOFX_X3

      G_F             (iNX_X,iX1,iX2,iX3,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX1,iX2,iX3,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX1,iX2,iX3,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX1,iX2,iX3,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX2,iX3,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX2,iX3,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX1,iX2,iX3,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_h_1   ) &
                 * G_F(iNX_X,iX1,iX2,iX3,iGF_h_2   ) &
                 * G_F(iNX_X,iX1,iX2,iX3,iGF_h_3   ), SqrtTiny )
      G_F             (iNX_X,iX1,iX2,iX3,iGF_Alpha) &
        = MAX( G_F    (iNX_X,iX1,iX2,iX3,iGF_Alpha), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               uCF_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCF), nDOFX_X3 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               uCF_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCF), nDOFX_X3 )

    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL ComputePrimitive_Euler &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$ACC          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, uD_L, uS3_L, uE_L, &
    !$ACC          pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, uD_R, uS3_R, uE_R, &
    !$ACC          Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F, &
    !$ACC          Alpha_F, Beta_3_F, &
    !$ACC          iErr, &
    !$ACC          uCF_L, uCF_R, uDF_L, uDF_R, NumericalFlux, &
    !$ACC          dX1, dX2, WeightsX_X3, IndexTableX_F, SurfaceFlux_X3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( AlphaMns, AlphaPls, AlphaMdl, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, Flux_L, Flux_R, Flux_F, &
    !$OMP          uCF_L_nCF, uCF_R_nCF, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_X = 1, nNodesX_X3

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_Euler &
            ( pV3_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_33_F(iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_3_F  (iNX_X) )

      Flux_L &
        = Flux_X3_Euler &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_3_F  (iNX_X) )

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_Euler &
            ( pV3_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_22_F(iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_3_F  (iNX_X) )

      Flux_R &
        = Flux_X3_Euler &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_3_F  (iNX_X) )

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iErr(iNX_X) = 0

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uD_L          (iNX_X), &
              uS3_L         (iNX_X), &
              uE_L          (iNX_X), &
              Flux_L(iCF_D)        , &
              Flux_L(iCF_S3)       , &
              Flux_L(iCF_E)        , &
              uD_R          (iNX_X), &
              uS3_R         (iNX_X), &
              uE_R          (iNX_X), &
              Flux_R(iCF_D)        , &
              Flux_R(iCF_S3)       , &
              Flux_R(iCF_E)        , &
              Gm_dd_33_F    (iNX_X), &
              AlphaPls             , &
              AlphaMns             , &
              Alpha_F       (iNX_X), &
              Beta_3_F      (iNX_X), &
              iErr          (iNX_X) )

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX2 = IndexTableX_F(3,iNX_X)
      iX3 = IndexTableX_F(4,iNX_X)

      DO iCF = 1, nCF

        uCF_L_nCF(iCF) = uCF_L(iNX,iX1,iX2,iX3,iCF)
        uCF_R_nCF(iCF) = uCF_R(iNX,iX1,iX2,iX3,iCF)

      END DO

      Flux_F &
        = NumericalFlux_Euler_X3 &
            ( uCF_L_nCF           , &
              uCF_R_nCF           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            , &
              AlphaMdl            , &
              Gm_dd_33_F(iNX_X)   , &
              pV3_L     (iNX_X)   , &
              pV3_R     (iNX_X)   , &
              P_L                 , &
              P_R                 , &
              Alpha_F   (iNX_X)   , &
              Beta_3_F  (iNX_X)   , &
              uDF_L(iX1,iX2,iX3,1), &
              uDF_R(iX1,iX2,iX3,2), &
              uDF_L(iX1,iX2,iX3,1), &
              uDF_R(iX1,iX2,iX3,2) )

      DO iCF = 1, nCF

        SurfaceFlux_X3(iNX,iX1,iX2,iX3,iCF) &
          = Flux_F(iCF) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

        NumericalFlux(iNX,iCF,iX1,iX2,iX3) &
          = Flux_F(iCF) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX1(iX1) * dX2(iX2) &
              * WeightsX_X3(iNX)

      END DO ! iCF

    END DO ! iNX_X

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X3, + One, LX_X3_Dn, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, &
             Zero, dU_X3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X3, - One, LX_X3_Up, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, &
             One,  dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL ComputePrimitive_Euler &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K )

    CALL TimersStop_Euler( Timer_Euler_DG_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
    !$ACC          Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, Alpha_K, Beta_3_K, &
    !$ACC          Flux_q, SqrtGm_K, dX1, dX2, WeightsX_q, IndexTableX_V )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( P_K, Flux_K, iNX, iX1, iX2, iX3 )
#endif
    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X3_Euler &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_3_K  (iNX_K) )

      iNX = IndexTableX_V(1,iNX_K)
      iX1 = IndexTableX_V(2,iNX_K)
      iX2 = IndexTableX_V(3,iNX_K)
      iX3 = IndexTableX_V(4,iNX_K)

      DO iCF = 1, nCF

        Flux_q(iNX,iCF,iX1,iX2,iX3) &
          = Flux_K(iCF) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX1(iX1) * dX2(iX2) * WeightsX_q(iNX)

      END DO

    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX, One, dLXdX3_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, dU_X3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1       , nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X3(iNX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

#ifdef THORNADO_DEBUG_EULER
#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET UPDATE FROM( dU_X3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC UPDATE HOST       ( dU_X3 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X3)', MAXLOC(dU_X3)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X3)', MAXVAL(dU_X3)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    NumericalFlux, iErr, SurfaceFlux_X3 ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX2, &
    !$OMP               EigVals_L, EigVals_R, &
    !$OMP               Flux_L, Flux_R, &
    !$OMP               Flux_F, Flux_K, &
    !$OMP               uCF_L_nCF, uCF_R_nCF, &
    !$OMP               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$OMP               Flux_q, dU_X3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      NumericalFlux, iErr, SurfaceFlux_X3 ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, iXP_B0, iXP_E0, dX1, dX2, &
    !$ACC               EigVals_L, EigVals_R, &
    !$ACC               Flux_L, Flux_R, &
    !$ACC               Flux_F, Flux_K, &
    !$ACC               uCF_L_nCF, uCF_R_nCF, &
    !$ACC               G_K, G_F, uCF_K, uCF_L, uCF_R, uDF_L, uDF_R, &
    !$ACC               Flux_q, dU_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX2   = iX_B0(2), iX_E0(2)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iCF   = 1       , nCF
    DO iNX_X = 1       , nDOFX_X3

      OffGridFlux_Euler_X3_Inner(iCF) &
        = OffGridFlux_Euler_X3_Inner(iCF) &
            + NumericalFlux(iNX_X,iCF,iX1,iX2,iX_B0(3))

      OffGridFlux_Euler_X3_Outer(iCF) &
        = OffGridFlux_Euler_X3_Outer(iCF) &
            + NumericalFlux(iNX_X,iCF,iX1,iX2,iX_E0(3)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

#ifdef HYDRO_RELATIVISTIC

    DO iNX_X = 1, nNodesX_X3

      CALL DescribeError_Euler( iErr(iNX_X) )

    END DO

#endif

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, UseXCFC, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G  (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U  (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      tau(:,iX_B1(1):,iX_B1(2):,iX_B1(3):)
    LOGICAL,  INTENT(in)    :: &
      UseXCFC
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

#ifdef HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, UseXCFC, dU )

#else

    CALL ComputeIncrement_Geometry_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Geometry_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNX
    REAL(DP) :: dX1, dX2
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh2dX1(nDOFX), dh3dX1(nDOFX), dh3dX2(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF)
    REAL(DP) :: G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)
    REAL(DP) :: G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF)

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

!      print*,"iX1, iX2, iX3 = ", iX1, iX2, iX3

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)

        END DO

      END IF

      CALL ComputePrimitive_Euler &
             ( uCF_K(:,iCF_D ), &
               uCF_K(:,iCF_S1), &
               uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), &
               uCF_K(:,iCF_E ), &
               uCF_K(:,iCF_Ne), &
               uPF_K(:,iPF_D ), &
               uPF_K(:,iPF_V1), &
               uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), &
               uPF_K(:,iPF_E ), &
               uPF_K(:,iPF_Ne), &
               G_K(:,iGF_Gm_dd_11), &
               G_K(:,iGF_Gm_dd_22), &
               G_K(:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      DO iNX = 1, nDOFX

        Stress(iNX,1:3) &
          = StressTensor_Diagonal_Euler &
              ( uCF_K(iNX,iCF_S1), &
                uCF_K(iNX,iCF_S2), &
                uCF_K(iNX,iCF_S3), &
                uPF_K(iNX,iPF_V1), &
                uPF_K(iNX,iPF_V2), &
                uPF_K(iNX,iPF_V3), &
                P_K  (iNX) )

      END DO

      ! --- Scale Factor Derivatives wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      DO iGF = iGF_h_2, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        G_X1_Dn(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Dn(1:nDOFX_X1,iGF), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        G_X1_Up(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Up(1:nDOFX_X1,iGF), SqrtTiny )

      END DO

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

      dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

      dh3dx1 = dh3dx1 / ( WeightsX_q(:) * dX1 )

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)  &
            + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3)

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale Factor Derivatives wrt X2 ---

        ! --- Face States (Average of Left and Right States) ---

        DO iGF = iGF_h_3, iGF_h_3

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          G_X2_Dn(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Dn(1:nDOFX_X2,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          G_X2_Up(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Up(1:nDOFX_X2,iGF), SqrtTiny )

        END DO

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Dn(:,iGF_h_3), 1,  One, dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX2, 1 )

        dh3dx2 = dh3dx2 / ( WeightsX_q(:) * dX2 )

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2) &
              + ( Stress(:,3) * dh3dX2(:) ) / G_K(:,iGF_h_3)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Geometry_NonRelativistic


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, UseXCFC, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(in)    :: &
      tau(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    LOGICAL,  INTENT(in)    :: &
      UseXCFC
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER :: iX1, iX2, iX3, iNX, iCF, iGF
    INTEGER :: nK(3), nGF_K

    REAL(DP) :: DivGridVolume
    REAL(DP) :: P(nPF)
    REAL(DP) :: Pressure
    REAL(DP) :: PressureTensor(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3))

    REAL(DP) :: G_K_X1 (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_Dn_X1(nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))
    REAL(DP) :: G_Up_X1(nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))
    REAL(DP) :: dGdX1  (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))

    REAL(DP) :: G_K_X2 (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_Dn_X2(nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))
    REAL(DP) :: G_Up_X2(nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))
    REAL(DP) :: dGdX2  (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))

    REAL(DP) :: G_K_X3 (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_Dn_X3(nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))
    REAL(DP) :: G_Up_X3(nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))
    REAL(DP) :: dGdX3  (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))

    REAL(DP) :: Christoffel3D_X1(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3))
    REAL(DP) :: Christoffel3D_X2(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3))
    REAL(DP) :: Christoffel3D_X3(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3))

    INTEGER :: iErr(nDOFX,iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3))

    REAL(DP) :: GradPsiF(nDOFX)

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nK = iX_E0 - iX_B0 + 1
    nGF_K = nGF * PRODUCT( nK )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, dX1, dX2, dX3, tau ) &
    !$OMP MAP( alloc: PressureTensor, &
    !$OMP             G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$OMP             G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$OMP             G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$OMP             Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, &
    !$OMP             iErr )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, dX1, dX2, dX3, tau ) &
    !$ACC CREATE(     PressureTensor, &
    !$ACC             G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$ACC             G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$ACC             G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$ACC             Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, &
    !$ACC             iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    ! --- Permute data ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K_X1, G )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iGF = 1           , nGF
    DO iNX = 1           , nDOFX

      G_K_X1(iNX,iGF,iX2,iX3,iX1) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X2, G )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
      DO iX3 = iX_B0(3)    , iX_E0(3)
      DO iX1 = iX_B0(1)    , iX_E0(1)
      DO iGF = 1           , nGF
      DO iNX = 1           , nDOFX

        G_K_X2(iNX,iGF,iX1,iX3,iX2) = G(iNX,iX1,iX2,iX3,iGF)

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X3, G )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
      DO iX2 = iX_B0(2)    , iX_E0(2)
      DO iX1 = iX_B0(1)    , iX_E0(1)
      DO iGF = 1           , nGF
      DO iNX = 1           , nDOFX

        G_K_X3(iNX,iGF,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF)

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate to faces ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- X1 ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, One , LX_X1_Up, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX, Half, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, One , LX_X1_Up, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Half, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

      ! --- Compute metric on faces from scale factors ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_Up_X1, G_Dn_X1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1       , nDOFX_X1

      G_Up_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Up_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Up_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_Up_X1        (iNX,iGF_SqrtGm,iX2,iX3,iX1) &
        = MAX( G_Up_X1    (iNX,iGF_h_1   ,iX2,iX3,iX1) &
                 * G_Up_X1(iNX,iGF_h_2   ,iX2,iX3,iX1) &
                 * G_Up_X1(iNX,iGF_h_3   ,iX2,iX3,iX1), SqrtTiny )

      G_Dn_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Dn_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Dn_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_Dn_X1        (iNX,iGF_SqrtGm,iX2,iX3,iX1) &
        = MAX( G_Dn_X1    (iNX,iGF_h_1   ,iX2,iX3,iX1) &
                 * G_Dn_X1(iNX,iGF_h_2   ,iX2,iX3,iX1) &
                 * G_Dn_X1(iNX,iGF_h_3   ,iX2,iX3,iX1), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

      ! --- X2 ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF_K, nDOFX, One , LX_X2_Up, nDOFX_X2, &
               G_K_X2 (1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Zero, &
               G_Up_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF_K, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_K_X2 (1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX, Half, &
               G_Up_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF_K, nDOFX, One , LX_X2_Up, nDOFX_X2, &
               G_K_X2 (1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
               G_Dn_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF_K, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_K_X2 (1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Half, &
               G_Dn_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Up_X2, G_Dn_X2 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX_X2

        G_Up_X2         (iNX,iGF_Gm_dd_11,iX1,iX3,iX2) &
          = MAX( G_Up_X2(iNX,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )
        G_Up_X2         (iNX,iGF_Gm_dd_22,iX1,iX3,iX2) &
          = MAX( G_Up_X2(iNX,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )
        G_Up_X2         (iNX,iGF_Gm_dd_33,iX1,iX3,iX2) &
          = MAX( G_Up_X2(iNX,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )

        G_Up_X2        (iNX,iGF_SqrtGm,iX1,iX3,iX2) &
          = MAX( G_Up_X2    (iNX,iGF_h_1   ,iX1,iX3,iX2) &
                   * G_Up_X2(iNX,iGF_h_2   ,iX1,iX3,iX2) &
                   * G_Up_X2(iNX,iGF_h_3   ,iX1,iX3,iX2), SqrtTiny )

        G_Dn_X2         (iNX,iGF_Gm_dd_11,iX1,iX3,iX2) &
          = MAX( G_Dn_X2(iNX,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )
        G_Dn_X2         (iNX,iGF_Gm_dd_22,iX1,iX3,iX2) &
          = MAX( G_Dn_X2(iNX,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )
        G_Dn_X2         (iNX,iGF_Gm_dd_33,iX1,iX3,iX2) &
          = MAX( G_Dn_X2(iNX,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )

        G_Dn_X2        (iNX,iGF_SqrtGm,iX1,iX3,iX2) &
          = MAX( G_Dn_X2    (iNX,iGF_h_1   ,iX1,iX3,iX2) &
                   * G_Dn_X2(iNX,iGF_h_2   ,iX1,iX3,iX2) &
                   * G_Dn_X2(iNX,iGF_h_3   ,iX1,iX3,iX2), SqrtTiny )

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

      ! --- X3 ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF_K, nDOFX, One , LX_X3_Up, nDOFX_X3, &
               G_K_X3 (1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Zero, &
               G_Up_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF_K, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               G_K_X3 (1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX, Half, &
               G_Up_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF_K, nDOFX, One , LX_X3_Up, nDOFX_X3, &
               G_K_X3 (1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
               G_Dn_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF_K, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               G_K_X3 (1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Half, &
               G_Dn_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Up_X3, G_Dn_X3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX_X3

        G_Up_X3         (iNX,iGF_Gm_dd_11,iX1,iX2,iX3) &
          = MAX( G_Up_X3(iNX,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )
        G_Up_X3         (iNX,iGF_Gm_dd_22,iX1,iX2,iX3) &
          = MAX( G_Up_X3(iNX,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )
        G_Up_X3         (iNX,iGF_Gm_dd_33,iX1,iX2,iX3) &
          = MAX( G_Up_X3(iNX,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )

        G_Up_X3        (iNX,iGF_SqrtGm,iX1,iX2,iX3) &
          = MAX( G_Up_X3    (iNX,iGF_h_1   ,iX1,iX2,iX3) &
                   * G_Up_X3(iNX,iGF_h_2   ,iX1,iX2,iX3) &
                   * G_Up_X3(iNX,iGF_h_3   ,iX1,iX2,iX3), SqrtTiny )

        G_Dn_X3         (iNX,iGF_Gm_dd_11,iX1,iX2,iX3) &
          = MAX( G_Dn_X3(iNX,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )
        G_Dn_X3         (iNX,iGF_Gm_dd_22,iX1,iX2,iX3) &
          = MAX( G_Dn_X3(iNX,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )
        G_Dn_X3         (iNX,iGF_Gm_dd_33,iX1,iX2,iX3) &
          = MAX( G_Dn_X3(iNX,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )

        G_Dn_X3        (iNX,iGF_SqrtGm,iX1,iX2,iX3) &
          = MAX( G_Dn_X3    (iNX,iGF_h_1   ,iX1,iX2,iX3) &
                   * G_Dn_X3(iNX,iGF_h_2   ,iX1,iX2,iX3) &
                   * G_Dn_X3(iNX,iGF_h_3   ,iX1,iX2,iX3), SqrtTiny )

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Compute derivatives (X1) ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X1, G_Up_X1, WeightsX_X1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF
    DO iNX = 1       , nDOFX_X1

      G_Dn_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Dn_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

      G_Up_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Up_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K_X1, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF
    DO iNX = 1       , nDOFX

      G_K_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_K_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_q(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, +One, LX_X1_Up, nDOFX_X1, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, Zero, &
             dGdX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, -One, LX_X1_Dn, nDOFX_X1, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, One,  &
             dGdX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX   , -One, dLXdX1_q, nDOFX   , &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX   , One , &
             dGdX1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dGdX1, WeightsX_q, dX1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF
    DO iNX = 1       , nDOFX

      dGdX1(iNX,iGF,iX2,iX3,iX1) &
        = dGdX1(iNX,iGF,iX2,iX3,iX1) / ( WeightsX_q(iNX) * dX1(iX1) )

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    IF( nDimsX .GT. 1 )THEN

      ! --- Compute derivatives (X2) ---

      CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X2, G_Up_X2, WeightsX_X2 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX_X2

        G_Dn_X2(iNX,iGF,iX1,iX3,iX2) &
          = G_Dn_X2(iNX,iGF,iX1,iX3,iX2) * WeightsX_X2(iNX)

        G_Up_X2(iNX,iGF,iX1,iX3,iX2) &
          = G_Up_X2(iNX,iGF,iX1,iX3,iX2) * WeightsX_X2(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X2, WeightsX_q )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
      DO iX3 = iX_B0(3)    , iX_E0(3)
      DO iX1 = iX_B0(1)    , iX_E0(1)
      DO iGF = 1           , nGF
      DO iNX = 1           , nDOFX

        G_K_X2(iNX,iGF,iX1,iX3,iX2) &
          = G_K_X2(iNX,iGF,iX1,iX3,iX2) * WeightsX_q(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

      CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, +One, LX_X2_Up, nDOFX_X2, &
               G_Up_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2, Zero, &
               dGdX2, nDOFX )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX_X2, -One, LX_X2_Dn, nDOFX_X2, &
               G_Dn_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2, One,  &
               dGdX2, nDOFX )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX   , -One, dLXdX2_q, nDOFX   , &
               G_K_X2 (1,1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX   , One , &
               dGdX2, nDOFX )

      CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

      CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX2, WeightsX_q, dX2 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX

        dGdX2(iNX,iGF,iX1,iX3,iX2) &
          = dGdX2(iNX,iGF,iX1,iX3,iX2) / ( WeightsX_q(iNX) * dX2(iX2) )

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ELSE

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX2 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX

        dGdX2(iNX,iGF,iX1,iX3,iX2) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

      ! --- Compute derivatives (X3) ---

      CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X3, G_Up_X3, WeightsX_X3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX_X3

        G_Dn_X3(iNX,iGF,iX1,iX2,iX3) &
          = G_Dn_X3(iNX,iGF,iX1,iX2,iX3) * WeightsX_X3(iNX)

        G_Up_X3(iNX,iGF,iX1,iX2,iX3) &
          = G_Up_X3(iNX,iGF,iX1,iX2,iX3) * WeightsX_X3(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X3, WeightsX_q )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
      DO iX2 = iX_B0(2)    , iX_E0(2)
      DO iX1 = iX_B0(1)    , iX_E0(1)
      DO iGF = 1           , nGF
      DO iNX = 1           , nDOFX

        G_K_X3(iNX,iGF,iX1,iX2,iX3) &
          = G_K_X3(iNX,iGF,iX1,iX2,iX3) * WeightsX_q(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

      CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, +One, LX_X3_Up, nDOFX_X3, &
               G_Up_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3, Zero, &
               dGdX3, nDOFX )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX_X3, -One, LX_X3_Dn, nDOFX_X3, &
               G_Dn_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3, One,  &
               dGdX3, nDOFX )

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nGF_K, nDOFX   , -One, dLXdX3_q, nDOFX   , &
               G_K_X3 (1,1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX   , One , &
               dGdX3, nDOFX )

      CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

      CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX3, WeightsX_q, dX3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX

        dGdX3(iNX,iGF,iX1,iX2,iX3) &
          = dGdX3(iNX,iGF,iX1,iX2,iX3) / ( WeightsX_q(iNX) * dX3(iX3) )

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ELSE

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1       , nGF
      DO iNX = 1       , nDOFX

        dGdX3(iNX,iGF,iX1,iX2,iX3) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    ! --- Contributions from time-independent metric ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( P, Pressure )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( P, Pressure ) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, U, G, &
    !$ACC          PressureTensor, tau, &
    !$ACC          dGdX1, dGdX2, dGdX3, iErr )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( P, Pressure )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      iErr(iNX,iX1,iX2,iX3) = 0

      CALL ComputePrimitive_Euler &
             ( U(   iNX,iX1,iX2,iX3,iCF_D ), &
               U(   iNX,iX1,iX2,iX3,iCF_S1), &
               U(   iNX,iX1,iX2,iX3,iCF_S2), &
               U(   iNX,iX1,iX2,iX3,iCF_S3), &
               U(   iNX,iX1,iX2,iX3,iCF_E ), &
               U(   iNX,iX1,iX2,iX3,iCF_Ne), &
               P(iPF_D ), &
               P(iPF_V1), &
               P(iPF_V2), &
               P(iPF_V3), &
               P(iPF_E ), &
               P(iPF_Ne), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               iErr(iNX,iX1,iX2,iX3) )

      CALL ComputePressureFromPrimitive &
             ( P(iPF_D), P(iPF_E), P(iPF_Ne), Pressure )

      PressureTensor(1,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCF_S1) * P(iPF_V1) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)

      PressureTensor(2,1,iNX,iX1,iX2,iX3) &
        =   U(iNX,iX1,iX2,iX3,iCF_S2) * P(iPF_V1) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)

      PressureTensor(3,1,iNX,iX1,iX2,iX3) &
        =   U(iNX,iX1,iX2,iX3,iCF_S3) * P(iPF_V1) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)

      PressureTensor(1,2,iNX,iX1,iX2,iX3) &
        = PressureTensor(2,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,2,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCF_S2) * P(iPF_V2) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)

      PressureTensor(3,2,iNX,iX1,iX2,iX3) &
        =   U(iNX,iX1,iX2,iX3,iCF_S3) * P(iPF_V2) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)

      PressureTensor(1,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,2,iNX,iX1,iX2,iX3)

      PressureTensor(3,3,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCF_S3) * P(iPF_V3) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)

      ! --- X1 increments ---

      ! --- Momentum increment ---

      dU(iNX,iX1,iX2,iX3,iCF_S1) &
        = dU(iNX,iX1,iX2,iX3,iCF_S1) &
            + tau(iNX,iX1,iX2,iX3) &
                * ( G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                      * (   PressureTensor(1,1,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_1) &
                              * dGdX1(iNX,iGF_h_1,iX2,iX3,iX1) &
                          + PressureTensor(2,2,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_2) &
                              * dGdX1(iNX,iGF_h_2,iX2,iX3,iX1) &
                          + PressureTensor(3,3,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_3) &
                              * dGdX1(iNX,iGF_h_3,iX2,iX3,iX1) ) &
                  + U(iNX,iX1,iX2,iX3,iCF_S1) &
                      * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCF_S2) &
                      * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCF_S3) &
                      * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) &
                  - ( U(iNX,iX1,iX2,iX3,iCF_D) + U(iNX,iX1,iX2,iX3,iCF_E) ) &
                      * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1) )

      ! --- Energy increment ---

      dU(iNX,iX1,iX2,iX3,iCF_E) &
        = dU(iNX,iX1,iX2,iX3,iCF_E) &
            -tau(iNX,iX1,iX2,iX3) &
              * U(iNX,iX1,iX2,iX3,iCF_S1) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
              * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      ! --- X2 increments ---

      IF( nDimsX .GT. 1 )THEN

        dU(iNX,iX1,iX2,iX3,iCF_S2) &
          = dU(iNX,iX1,iX2,iX3,iCF_S2) &
              + tau(iNX,iX1,iX2,iX3) &
                  * ( G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                        * (   PressureTensor(1,1,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_1) &
                                  * dGdX2(iNX,iGF_h_1,iX1,iX3,iX2) &
                            + PressureTensor(2,2,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_2) &
                                  * dGdX2(iNX,iGF_h_2,iX1,iX3,iX2) &
                            + PressureTensor(3,3,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_3) &
                                  * dGdX2(iNX,iGF_h_3,iX1,iX3,iX2) ) &
                    + U(iNX,iX1,iX2,iX3,iCF_S1) &
                        * dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) &
                    + U(iNX,iX1,iX2,iX3,iCF_S2) &
                        * dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
                    + U(iNX,iX1,iX2,iX3,iCF_S3) &
                        * dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) &
                    - ( U(iNX,iX1,iX2,iX3,iCF_D) + U(iNX,iX1,iX2,iX3,iCF_E) ) &
                        * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2) )

        dU(iNX,iX1,iX2,iX3,iCF_E) &
          = dU(iNX,iX1,iX2,iX3,iCF_E) &
              -tau(iNX,iX1,iX2,iX3) &
                 * U(iNX,iX1,iX2,iX3,iCF_S2) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                 * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2)

      END IF

      ! --- X3 increments ---

      IF( nDimsX .GT. 2 )THEN

        dU(iNX,iX1,iX2,iX3,iCF_S3) &
          = dU(iNX,iX1,iX2,iX3,iCF_S3) &
              + tau(iNX,iX1,iX2,iX3) &
                  * ( G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                        * (   PressureTensor(1,1,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_1) &
                                  * dGdX3(iNX,iGF_h_1,iX1,iX2,iX3) &
                            + PressureTensor(2,2,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_2) &
                                  * dGdX3(iNX,iGF_h_2,iX1,iX2,iX3) &
                            + PressureTensor(3,3,iNX,iX1,iX2,iX3) &
                                * G(iNX,iX1,iX2,iX3,iGF_h_3) &
                                  * dGdX3(iNX,iGF_h_3,iX1,iX2,iX3) ) &
                    + U(iNX,iX1,iX2,iX3,iCF_S1) &
                        * dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) &
                    + U(iNX,iX1,iX2,iX3,iCF_S2) &
                        * dGdX3(iNX,iGF_Beta_2,iX1,iX2,iX3) &
                    + U(iNX,iX1,iX2,iX3,iCF_S3) &
                        * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
                    - ( U(iNX,iX1,iX2,iX3,iCF_D) + U(iNX,iX1,iX2,iX3,iCF_E) ) &
                        * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3) )

        dU(iNX,iX1,iX2,iX3,iCF_E) &
          = dU(iNX,iX1,iX2,iX3,iCF_E) &
              -tau(iNX,iX1,iX2,iX3) &
                 * U(iNX,iX1,iX2,iX3,iCF_S3) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                 * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3)

      END IF

    END DO
    END DO
    END DO
    END DO

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Contributions from time-dependent metric ---

    CALL ComputeChristoffelSymbols_3D_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX1, dGdX2, dGdX3, &
             Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3 )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( DivGridVolume )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( DivGridVolume ) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, U, G, dGdX1, dGdX2, dGdX3, &
    !$ACC          tau, PressureTensor, &
    !$ACC          Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( DivGridVolume )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      DivGridVolume &
        =   dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
          + dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
          + dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
          + (   Christoffel3D_X1(1,1,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,1,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
          + (   Christoffel3D_X1(1,2,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,2,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
          + (   Christoffel3D_X1(1,3,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,3,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,3,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_3)

      dU(iNX,iX1,iX2,iX3,iCF_E) &
        = dU(iNX,iX1,iX2,iX3,iCF_E) &
            + tau(iNX,iX1,iX2,iX3) * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                 * (    PressureTensor(1,1,iNX,iX1,iX2,iX3) &
                          * G(iNX,iX1,iX2,iX3,iGF_K_dd_11) &
                     +  PressureTensor(2,2,iNX,iX1,iX2,iX3) &
                         * G(iNX,iX1,iX2,iX3,iGF_K_dd_22) &
                     +  PressureTensor(3,3,iNX,iX1,iX2,iX3) &
                         * G(iNX,iX1,iX2,iX3,iGF_K_dd_33) &
                     + Two * PressureTensor(1,2,iNX,iX1,iX2,iX3) &
                         * G(iNX,iX1,iX2,iX3,iGF_K_dd_12) &
                     + Two * PressureTensor(1,3,iNX,iX1,iX2,iX3) &
                         * G(iNX,iX1,iX2,iX3,iGF_K_dd_13) &
                     + Two * PressureTensor(2,3,iNX,iX1,iX2,iX3) &
                         * G(iNX,iX1,iX2,iX3,iGF_K_dd_23) )

      IF( .NOT. UseXCFC )THEN

        DO iCF = 1, nCF

          dU(iNX,iX1,iX2,iX3,iCF) &
            = dU(iNX,iX1,iX2,iX3,iCF) &
                - U(iNX,iX1,iX2,iX3,iCF) * DivGridVolume

        END DO

      END IF

    END DO
    END DO
    END DO
    END DO

    ! --- GradPsiF ---

    GradPsiF = Zero

    CALL DGEMV &
           ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             dGdX1(:,iGF_Psi,1,1,iX_E0(1)), 1, Zero, GradPsiF(1:nDOFX_X1), 1 )

#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr ) &
    !$OMP MAP( release: iX_B0, iX_E0, dX1, dX2, dX3, tau, &
    !$OMP               PressureTensor, &
    !$OMP               G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$OMP               G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$OMP               G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$OMP               Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr ) &
    !$ACC DELETE(       iX_B0, iX_E0, dX1, dX2, dX3, tau, &
    !$ACC               PressureTensor, &
    !$ACC               G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$ACC               G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$ACC               G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$ACC               Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    END ASSOCIATE ! dX1, dX2, dX3

    CALL TimersStart_Euler( Timer_Euler_DG_ErrorCheck )

    IF( ANY( iErr .NE. 0 ) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'ERROR: ComputeIncrement_Geometry_Relativistic'

          WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr(iNX,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_ErrorCheck )

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic


  SUBROUTINE ComputeIncrement_Gravity &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

   INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

#ifdef HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Gravity_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#else

    CALL ComputeIncrement_Gravity_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

  END SUBROUTINE ComputeIncrement_Gravity


  SUBROUTINE ComputeIncrement_Gravity_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: Phi_P_X1(nDOFX)
    REAL(DP) :: Phi_K   (nDOFX)
    REAL(DP) :: Phi_N_X1(nDOFX)
    REAL(DP) :: dPhidX1(nDOFX)
    REAL(DP) :: Phi_X1_Dn(nDOFX_X1)
    REAL(DP) :: Phi_X1_Up(nDOFX_X1)
    REAL(DP) :: uCF_K(nDOFX,nCF)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      Phi_P_X1(:) = G(:,iX1-1,iX2,iX3,iGF_Phi_N)
      Phi_K   (:) = G(:,iX1,  iX2,iX3,iGF_Phi_N)
      Phi_N_X1(:) = G(:,iX1+1,iX2,iX3,iGF_Phi_N)

      ! --- Derivative of Potential wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      ! --- Face at X1_L ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               Phi_P_X1(:), 1, Zero, Phi_X1_Dn(:), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               Phi_K   (:), 1, Half, Phi_X1_Dn(:), 1 )

      ! --- Face at X1_H ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               Phi_K   (:), 1, Zero, Phi_X1_Up(:), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               Phi_N_X1(:), 1, Half, Phi_X1_Up(:), 1 )

      ! --- dPhidX1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * Phi_X1_Up(:), 1, Zero, dPhidX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * Phi_X1_Dn(:), 1,  One, dPhidX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * Phi_K    (:), 1,  One, dPhidX1, 1 )

      dPhidX1 = dPhidX1 / ( WeightsX_q(:) * dX1 )

      ! --- Increments ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            - uCF_K(:,iCF_D) * dPhidX1(:)

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            - uCF_K(:,iCF_S1) * dPhidX1(:)

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Gravity_NonRelativistic


  SUBROUTINE ComputeIncrement_Gravity_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
  END SUBROUTINE ComputeIncrement_Gravity_Relativistic


  SUBROUTINE InitializeIncrement_Euler( iX_B0, iX_E0, iX_B1, iX_E1 )

    INTEGER, INTENT(in) :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    nX        = iX_E0 - iX_B0 + 1 ! Number of Elements per Dimension
    nX_K      = PRODUCT( nX )     ! Number of Elements in Position Space
    nNodesX_K = nDOFX * nX_K      ! Number of Nodes in Elements

    nX_X1 = nX + [1,0,0] ! Number of X1 Faces per Dimension
    nX_X2 = nX + [0,1,0] ! Number of X2 Faces per Dimension
    nX_X3 = nX + [0,0,1] ! Number of X3 Faces per Dimension

    nX1_X = PRODUCT( nX_X1 ) ! Number of X1 Faces
    nX2_X = PRODUCT( nX_X2 ) ! Number of X2 Faces
    nX3_X = PRODUCT( nX_X3 ) ! Number of X3 Faces

    nNodesX_X1 = nDOFX_X1 * nX1_X ! Number of Nodes on X1 Faces
    nNodesX_X2 = nDOFX_X2 * nX2_X ! Number of Nodes on X2 Faces
    nNodesX_X3 = nDOFX_X3 * nX3_X ! Number of Nodes on X3 Faces

    ALLOCATE( pD_K (nNodesX_K) )
    ALLOCATE( pV1_K(nNodesX_K) )
    ALLOCATE( pV2_K(nNodesX_K) )
    ALLOCATE( pV3_K(nNodesX_K) )
    ALLOCATE( pE_K (nNodesX_K) )
    ALLOCATE( pNe_K(nNodesX_K) )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    dX1, dX2, dX3, &
    !$OMP             nX, nX_X1, nX_X2, nX_X3 ) &
    !$OMP MAP( alloc: pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     dX1, dX2, dX3, &
    !$ACC             nX, nX_X1, nX_X2, nX_X3 ) &
    !$ACC CREATE(     pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K )
#endif

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE InitializeIncrement_Euler


  SUBROUTINE FinalizeIncrement_Euler

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX1, dX2, dX3, &
    !$OMP               nX, nX_X1, nX_X2, nX_X3, &
    !$OMP               pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC DELETE(       dX1, dX2, dX3, &
    !$ACC               nX, nX_X1, nX_X2, nX_X3,  &
    !$ACC               pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K )
#endif

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K )

  END SUBROUTINE FinalizeIncrement_Euler


  SUBROUTINE InitializeIncrement_Divergence &
    ( iXP_B0, iXP_E0, nDOFX_X, G_K, G_F, uCF_K, uCF_L, uCF_R )

    INTEGER, INTENT(in) :: iXP_B0(3), iXP_E0(3) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X              ! nDOFX_X1, ...

    ! --- Geometry Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      G_K (1:nDOFX, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)-1:iXP_E0(3)+1, &
           1:nGF), &
      G_F (1:nDOFX_X, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)  :iXP_E0(3)+1, &
           1:nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCF_K(1:nDOFX, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)-1:iXP_E0(3)+1, &
            1:nCF), &
      uCF_L(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCF), &
      uCF_R(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCF)

    INTEGER :: iNX, iX1, iX2, iX3, iCF
    INTEGER :: nXP(3), nXP_X(3), nX_X, nNodesX_X, iX_F, iX_V

    nXP       = iXP_E0 - iXP_B0 + 1
    nXP_X     = nXP + [0,0,1]
    nX_X      = PRODUCT( nXP_X )
    nNodesX_X = nDOFX_X * nX_X

    Gm_dd_11_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_SqrtGm  )
    Beta_1_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_1  )
    Beta_2_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_2  )
    Beta_3_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_3  )
    Alpha_K   (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Alpha   )

    Gm_dd_11_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_SqrtGm  )
    Beta_1_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_1  )
    Beta_2_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_2  )
    Beta_3_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_3  )
    Alpha_F   (1:nNodesX_X) => G_F(:,:,:,:,iGF_Alpha   )

    uD_K (1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_D )
    uS1_K(1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_S1)
    uS2_K(1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_S2)
    uS3_K(1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_S3)
    uE_K (1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_E )
    uNe_K(1:nNodesX_K) => uCF_K(:,:,:,iXP_B0(3):iXP_E0(3),iCF_Ne)

    uD_L (1:nNodesX_X) => uCF_L(:,:,:,:,iCF_D )
    uS1_L(1:nNodesX_X) => uCF_L(:,:,:,:,iCF_S1)
    uS2_L(1:nNodesX_X) => uCF_L(:,:,:,:,iCF_S2)
    uS3_L(1:nNodesX_X) => uCF_L(:,:,:,:,iCF_S3)
    uE_L (1:nNodesX_X) => uCF_L(:,:,:,:,iCF_E )
    uNe_L(1:nNodesX_X) => uCF_L(:,:,:,:,iCF_Ne)

    uD_R (1:nNodesX_X) => uCF_R(:,:,:,:,iCF_D )
    uS1_R(1:nNodesX_X) => uCF_R(:,:,:,:,iCF_S1)
    uS2_R(1:nNodesX_X) => uCF_R(:,:,:,:,iCF_S2)
    uS3_R(1:nNodesX_X) => uCF_R(:,:,:,:,iCF_S3)
    uE_R (1:nNodesX_X) => uCF_R(:,:,:,:,iCF_E )
    uNe_R(1:nNodesX_X) => uCF_R(:,:,:,:,iCF_Ne)

    ALLOCATE( pD_L (nNodesX_X) )
    ALLOCATE( pV1_L(nNodesX_X) )
    ALLOCATE( pV2_L(nNodesX_X) )
    ALLOCATE( pV3_L(nNodesX_X) )
    ALLOCATE( pE_L (nNodesX_X) )
    ALLOCATE( pNe_L(nNodesX_X) )

    ALLOCATE( pD_R (nNodesX_X) )
    ALLOCATE( pV1_R(nNodesX_X) )
    ALLOCATE( pV2_R(nNodesX_X) )
    ALLOCATE( pV3_R(nNodesX_X) )
    ALLOCATE( pE_R (nNodesX_X) )
    ALLOCATE( pNe_R(nNodesX_X) )

    ALLOCATE( IndexTableX_F(4,nNodesX_X) )
    ALLOCATE( IndexTableX_V(4,nNodesX_K) )

    DO iX3 = iXP_B0(3), iXP_E0(3) + 1
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX_X

      iX_F = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX_X &
               + ( iX2 - iXP_B0(2) ) * nDOFX_X * nXP_X(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX_X * nXP_X(1) * nXP_X(2)

      IndexTableX_F(:,iX_F) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iXP_B0(3), iXP_E0(3)
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX

      iX_V = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX &
               + ( iX2 - iXP_B0(2) ) * nDOFX * nXP(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX * nXP(1) * nXP(2)

      IndexTableX_V(:,iX_V) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    IndexTableX_F, IndexTableX_V ) &
    !$OMP MAP( alloc: pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
    !$OMP             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     IndexTableX_F, IndexTableX_V ) &
    !$ACC CREATE(     pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
    !$ACC             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R )
#endif

  END SUBROUTINE InitializeIncrement_Divergence


  SUBROUTINE FinalizeIncrement_Divergence

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
    !$OMP               pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
    !$OMP               IndexTableX_F, IndexTableX_V )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC DELETE(       pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
    !$ACC               pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
    !$ACC               IndexTableX_F, IndexTableX_V )
#endif

    DEALLOCATE( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L )
    DEALLOCATE( pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R )
    DEALLOCATE( IndexTableX_F, IndexTableX_V )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY(   Beta_1_K,   Beta_2_K,   Beta_3_K,  Alpha_K )
    NULLIFY(   Beta_1_F,   Beta_2_F,   Beta_3_F,  Alpha_F )
    NULLIFY( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K )
    NULLIFY( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L )
    NULLIFY( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R )

  END SUBROUTINE FinalizeIncrement_Divergence


END MODULE Euler_dgDiscretizationModule
