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
    WeightsX_q, &
    NodeNumberTableX
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
    nAF, &
    iAF_P, &
    nDF, &
    iDF_Sh_X1, &
    iDF_Sh_X2, &
    iDF_Sh_X3
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
    Timer_Euler_ComputePrimitive, &
    Timer_Euler_Geometry, &
    Timer_Euler_Gravity, &
    Timer_Euler_SurfaceTerm, &
    Timer_Euler_VolumeTerm, &
    Timer_Euler_Increment, &
    Timer_Euler_DG_CopyIn, &
    Timer_Euler_DG_Permute, &
    Timer_Euler_DG_Interpolate, &
    Timer_Euler_DG_CopyOut, &
    Timer_Euler_DG_ErrorCheck
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

#ifndef USE_AMREX_TRUE

  USE InputOutputModuleHDF, ONLY: &
    WriteSourceTermDiagnosticsHDF

#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_DG_Explicit

  LOGICAL,  PUBLIC :: WriteSourceTerms
  REAL(DP), PUBLIC :: Time


CONTAINS


  SUBROUTINE ComputeIncrement_Euler_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SuppressBC_Option )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)         :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)         :: &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)           :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      SuppressBC_Option

    INTEGER :: iNX, iX1, iX2, iX3, iCF
    LOGICAL :: SuppressBC

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G, U, D ) &
    !$OMP MAP( alloc: dU )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G, U, D ) &
    !$ACC CREATE(     dU )
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

    CALL TimersStart_Euler( Timer_Euler_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, dU )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
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

    CALL TimersStop_Euler( Timer_Euler_Increment )

    CALL TimersStart_Euler( Timer_Euler_Divergence )

    CALL ComputeIncrement_Euler_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Euler_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Euler_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL TimersStop_Euler( Timer_Euler_Divergence )

    ! --- Multiply Inverse Mass Matrix ---

    CALL TimersStart_Euler( Timer_Euler_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, dX1, dX2, dX3, dU, G, WeightsX_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
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
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Geometry )

    CALL TimersStart_Euler( Timer_Euler_Gravity )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Gravity )

    CALL TimersStop_Euler( Timer_Euler_DG )

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU)', MAXLOC(dU)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU)', MAXVAL(dU)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    dU, U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      dU, U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, dX3, G )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Euler_DG_Explicit


  SUBROUTINE ComputeIncrement_Euler_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER  :: nK(3), nF_X1(3), nCF_K, nCF_F, nGF_F, nDF_F
    INTEGER  :: iNX, iNX_X1, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uCF_F_L(nCF), uCF_F_R(nCF)
    REAL(DP) :: Flux_X1_L(nCF), Flux_X1_R(nCF)
    REAL(DP) :: Flux_X1_F(nCF), Flux_X1_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L, P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_F          (nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    REAL(DP) :: uCF_L        (nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uCF_R        (nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)-1:iX_E0(1)+1)

    REAL(DP) :: uPF_L        (nDOFX_X1,nPF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uPF_R        (nDOFX_X1,nPF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uPF_K        (nDOFX,   nPF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)-1:iX_E0(1)+1)

    REAL(DP) :: uDF_L        (2           ,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: uDF_R        (2           ,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    REAL(DP) :: dU_X1        (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)  )

    REAL(DP) :: Flux_X1_q    (nDOFX,   nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)  )

    REAL(DP) :: NumericalFlux(nDOFX_X1,nCF,iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(1)  :iX_E0(1)+1)

    INTEGER :: iErr_L(1:nDOFX_X1,iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(1):iX_E0(1)+1 )
    INTEGER :: iErr_R(1:nDOFX_X1,iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(1):iX_E0(1)+1 )
    INTEGER :: iErr_M(1:nDOFX_X1,iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(1):iX_E0(1)+1 )
    INTEGER :: iErr_V(1:nDOFX   ,iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(1):iX_E0(1) )

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nF_X1 = nK + [1,0,0]           ! Number of X1 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nF_X1 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nF_X1 ) ! Number of Geometry Fields on Interfaces
    nDF_F = 2   * PRODUCT( nF_X1 ) ! Number of Diagnostic Fields on Interfaces

    ASSOCIATE( dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, dX2, dX3 ) &
    !$OMP MAP( alloc: uCF_L, uCF_R, uCF_K, &
    !$OMP             uPF_L, uPF_R, uPF_K, &
    !$OMP             uDF_L, uDF_R,        &
    !$OMP             G_K, G_F, dU_X1, Flux_X1_q, NumericalFlux, &
    !$OMP             iErr_L, iErr_R, iErr_M, iErr_V )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, dX2, dX3 ) &
    !$ACC CREATE(     uCF_L, uCF_R, uCF_K, &
    !$ACC             uPF_L, uPF_R, uPF_K, &
    !$ACC             uDF_L, uDF_R,        &
    !$ACC             G_K, G_F, dU_X1, Flux_X1_q, NumericalFlux, &
    !$ACC             iErr_L, iErr_R, iErr_M, iErr_V )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF
    DO iNX = 1, nDOFX

      G_K(iNX,iGF,iX2,iX3,iX1) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_F, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_F, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX_X1 = 1, nDOFX_X1

      G_F         (iNX_X1,iGF_h_1,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_1,iX2,iX3,iX1), SqrtTiny )

      G_F         (iNX_X1,iGF_h_2,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_2,iX2,iX3,iX1), SqrtTiny )

      G_F         (iNX_X1,iGF_h_3,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_3,iX2,iX3,iX1), SqrtTiny )

      G_F         (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F         (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F         (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_F        (iNX_X1,iGF_SqrtGm,iX2,iX3,iX1) &
        = G_F    (iNX_X1,iGF_h_1   ,iX2,iX3,iX1) &
            * G_F(iNX_X1,iGF_h_2   ,iX2,iX3,iX1) &
            * G_F(iNX_X1,iGF_h_3   ,iX2,iX3,iX1)

      G_F         (iNX_X1,iGF_Alpha,iX2,iX3,iX1) &
        = MAX( G_F(iNX_X1,iGF_Alpha,iX2,iX3,iX1), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, uCF_K, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      uCF_K(iNX,iCF,iX2,iX3,iX1) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      uDF_L(1,iX2,iX3,iX1) = D(1,iX1-1,iX2,iX3,iDF_Sh_X2)
      uDF_L(2,iX2,iX3,iX1) = D(1,iX1-1,iX2,iX3,iDF_Sh_X3)
      uDF_R(1,iX2,iX3,iX1) = D(1,iX1  ,iX2,iX3,iDF_Sh_X2)
      uDF_R(2,iX2,iX3,iX1) = D(1,iX1  ,iX2,iX3,iDF_Sh_X3)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nCF_F, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             uCF_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nCF_F, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             uCF_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F, uCF_L, uCF_R, uPF_L, uPF_R, &
    !$ACC          iErr_L, iErr_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX_X1 = 1, nDOFX_X1

      iErr_L(iNX_X1,iX2,iX3,iX1) = 0
      iErr_R(iNX_X1,iX2,iX3,iX1) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_L (iNX_X1,iCF_D       ,iX2,iX3,iX1), &
               uCF_L (iNX_X1,iCF_S1      ,iX2,iX3,iX1), &
               uCF_L (iNX_X1,iCF_S2      ,iX2,iX3,iX1), &
               uCF_L (iNX_X1,iCF_S3      ,iX2,iX3,iX1), &
               uCF_L (iNX_X1,iCF_E       ,iX2,iX3,iX1), &
               uCF_L (iNX_X1,iCF_Ne      ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_D       ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_E       ,iX2,iX3,iX1), &
               uPF_L (iNX_X1,iPF_Ne      ,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
               iErr_L(iNX_X1             ,iX2,iX3,iX1) )

      CALL ComputePrimitive_Euler &
             ( uCF_R (iNX_X1,iCF_D       ,iX2,iX3,iX1), &
               uCF_R (iNX_X1,iCF_S1      ,iX2,iX3,iX1), &
               uCF_R (iNX_X1,iCF_S2      ,iX2,iX3,iX1), &
               uCF_R (iNX_X1,iCF_S3      ,iX2,iX3,iX1), &
               uCF_R (iNX_X1,iCF_E       ,iX2,iX3,iX1), &
               uCF_R (iNX_X1,iCF_Ne      ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_D       ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_E       ,iX2,iX3,iX1), &
               uPF_R (iNX_X1,iPF_Ne      ,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
               G_F   (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
               iErr_R(iNX_X1             ,iX2,iX3,iX1) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X1_F, AlphaMns, AlphaPls, AlphaMdl )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$ACC          Flux_X1_F, AlphaMns, AlphaPls, AlphaMdl ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX2, dX3, &
    !$ACC          G_F, uCF_L, uCF_R, uPF_L, uPF_R, uDF_L, uDF_R, &
    !$ACC          NumericalFlux, WeightsX_X1, iErr_M )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X1_F, AlphaMns, AlphaPls, AlphaMdl )
#endif
    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX_X1 = 1, nDOFX_X1

      iErr_M(iNX_X1,iX2,iX3,iX1) = 0

      DO iCF = 1, nCF

        uCF_F_L(iCF) = uCF_L(iNX_X1,iCF,iX2,iX3,iX1)
        uCF_F_R(iCF) = uCF_R(iNX_X1,iCF,iX2,iX3,iX1)

      END DO

      ! --- Left State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_L(iNX_X1,iPF_D ,iX2,iX3,iX1), &
               uPF_L(iNX_X1,iPF_E ,iX2,iX3,iX1), &
               uPF_L(iNX_X1,iPF_Ne,iX2,iX3,iX1), &
               P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(iNX_X1,iPF_D ,iX2,iX3,iX1), &
               uPF_L(iNX_X1,iPF_E ,iX2,iX3,iX1), &
               uPF_L(iNX_X1,iPF_Ne,iX2,iX3,iX1), &
               Cs_L )

      EigVals_L &
        = Eigenvalues_Euler &
            ( uPF_L(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              Cs_L,                                   &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

      Flux_X1_L &
        = Flux_X1_Euler &
            ( uPF_L(iNX_X1,iPF_D       ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_E       ,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_Ne      ,iX2,iX3,iX1), &
              P_L,                                    &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

      ! --- Right State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_R(iNX_X1,iPF_D ,iX2,iX3,iX1), &
               uPF_R(iNX_X1,iPF_E ,iX2,iX3,iX1), &
               uPF_R(iNX_X1,iPF_Ne,iX2,iX3,iX1), &
               P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(iNX_X1,iPF_D ,iX2,iX3,iX1), &
               uPF_R(iNX_X1,iPF_E ,iX2,iX3,iX1), &
               uPF_R(iNX_X1,iPF_Ne,iX2,iX3,iX1), &
               Cs_R )

      EigVals_R &
        = Eigenvalues_Euler &
            ( uPF_R(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              Cs_R,                                   &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

      Flux_X1_R &
        = Flux_X1_Euler &
            ( uPF_R(iNX_X1,iPF_D       ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V2      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V3      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_E       ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_Ne      ,iX2,iX3,iX1), &
              P_R,                                    &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_22,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Gm_dd_33,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1) )

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uCF_F_L  (iCF_D), uCF_F_L  (iCF_S1), uCF_F_L  (iCF_E), &
              Flux_X1_L(iCF_D), Flux_X1_L(iCF_S1), Flux_X1_L(iCF_E), &
              uCF_F_R  (iCF_D), uCF_F_R  (iCF_S1), uCF_F_R  (iCF_E), &
              Flux_X1_R(iCF_D), Flux_X1_R(iCF_S1), Flux_X1_R(iCF_E), &
              G_F(iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), AlphaPls, AlphaMns, &
              G_F(iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F(iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1), &
              iErr_M(iNX_X1,iX2,iX3,iX1) )

      Flux_X1_F &
        = NumericalFlux_Euler_X1 &
            ( uCF_F_L, uCF_F_R, Flux_X1_L, Flux_X1_R, &
              AlphaPls, AlphaMns, AlphaMdl,           &
              G_F  (iNX_X1,iGF_Gm_dd_11,iX2,iX3,iX1), &
              uPF_L(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              uPF_R(iNX_X1,iPF_V1      ,iX2,iX3,iX1), &
              P_L, P_R,                               &
              G_F  (iNX_X1,iGF_Alpha   ,iX2,iX3,iX1), &
              G_F  (iNX_X1,iGF_Beta_1  ,iX2,iX3,iX1), &
              uDF_L(1                  ,iX2,iX3,iX1), &
              uDF_R(1                  ,iX2,iX3,iX1), &
              uDF_L(2                  ,iX2,iX3,iX1), &
              uDF_R(2                  ,iX2,iX3,iX1) )

      DO iCF = 1, nCF

        NumericalFlux    (iNX_X1,iCF       ,iX2,iX3,iX1) &
          = Flux_X1_F(iCF)                               &
              * G_F      (iNX_X1,iGF_Alpha ,iX2,iX3,iX1) &
              * G_F      (iNX_X1,iGF_SqrtGm,iX2,iX3,iX1) &
              * dX2(iX2) * dX3(iX3) * WeightsX_X1(iNX_X1)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X1, + One, LX_X1_Dn, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, Zero, &
             dU_X1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X1, - One, LX_X1_Up, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One, &
             dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, uCF_K, uPF_K, iErr_V )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      iErr_V(iNX,iX2,iX3,iX1) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_K (iNX,iCF_D       ,iX2,iX3,iX1), &
               uCF_K (iNX,iCF_S1      ,iX2,iX3,iX1), &
               uCF_K (iNX,iCF_S2      ,iX2,iX3,iX1), &
               uCF_K (iNX,iCF_S3      ,iX2,iX3,iX1), &
               uCF_K (iNX,iCF_E       ,iX2,iX3,iX1), &
               uCF_K (iNX,iCF_Ne      ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_D       ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_V1      ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_V2      ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_V3      ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_E       ,iX2,iX3,iX1), &
               uPF_K (iNX,iPF_Ne      ,iX2,iX3,iX1), &
               G_K   (iNX,iGF_Gm_dd_11,iX2,iX3,iX1), &
               G_K   (iNX,iGF_Gm_dd_22,iX2,iX3,iX1), &
               G_K   (iNX,iGF_Gm_dd_33,iX2,iX3,iX1), &
               iErr_V(iNX             ,iX2,iX3,iX1) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_K, P_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X1_K, P_K ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX2, dX3, &
    !$ACC          G_K, uCF_K, uPF_K, Flux_X1_q, WeightsX_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X1_K, P_K )
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      CALL ComputePressureFromPrimitive &
             ( uPF_K(iNX,iPF_D ,iX2,iX3,iX1), &
               uPF_K(iNX,iPF_E ,iX2,iX3,iX1), &
               uPF_K(iNX,iPF_Ne,iX2,iX3,iX1), &
               P_K )

      Flux_X1_K(1:nCF) &
        = Flux_X1_Euler &
          ( uPF_K(iNX,iPF_D       ,iX2,iX3,iX1), &
            uPF_K(iNX,iPF_V1      ,iX2,iX3,iX1), &
            uPF_K(iNX,iPF_V2      ,iX2,iX3,iX1), &
            uPF_K(iNX,iPF_V3      ,iX2,iX3,iX1), &
            uPF_K(iNX,iPF_E       ,iX2,iX3,iX1), &
            uPF_K(iNX,iPF_Ne      ,iX2,iX3,iX1), &
            P_K,                                 &
            G_K  (iNX,iGF_Gm_dd_11,iX2,iX3,iX1), &
            G_K  (iNX,iGF_Gm_dd_22,iX2,iX3,iX1), &
            G_K  (iNX,iGF_Gm_dd_33,iX2,iX3,iX1), &
            G_K  (iNX,iGF_Alpha   ,iX2,iX3,iX1), &
            G_K  (iNX,iGF_Beta_1  ,iX2,iX3,iX1) )

      DO iCF = 1, nCF

        Flux_X1_q(iNX,iCF,iX2,iX3,iX1) = Flux_X1_K(iCF)

        Flux_X1_q    (iNX,iCF       ,iX2,iX3,iX1) &
          = Flux_X1_q(iNX,iCF       ,iX2,iX3,iX1) &
              * G_K  (iNX,iGF_Alpha ,iX2,iX3,iX1) &
              * G_K  (iNX,iGF_SqrtGm,iX2,iX3,iX1) &
              * dX2(iX2) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX1_q, nDOFX, &
             Flux_X1_q, nDOFX, One, dU_X1, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, dU_X1 )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X1(iNX,iCF,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST       ( dU_X1 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X1)', MAXLOC(dU_X1)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X1)', MAXVAL(dU_X1)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$OMP MAP( release: dX2, dX3, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP               uCF_L, uCF_R, uCF_K, &
    !$OMP               uPF_L, uPF_R, uPF_K, &
    !$OMP               uDF_L, uDF_R,        &
    !$OMP               G_K, G_F, dU_X1, Flux_X1_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, dX2, dX3, &
    !$ACC               uCF_L, uCF_R, uCF_K, &
    !$ACC               uPF_L, uPF_R, uPF_K, &
    !$ACC               uDF_L, uDF_R,        &
    !$ACC               G_K, G_F, dU_X1, Flux_X1_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_DG_ErrorCheck )

    IF( ANY( iErr_L .NE. 0 ) .OR. &
        ANY( iErr_R .NE. 0 ) .OR. &
        ANY( iErr_M .NE. 0 ) .OR. &
        ANY( iErr_V .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR: Surface term (X1) (Left)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1) + 1
      DO iNX_X1 = 1, nDOFX_X1

        IF( iErr_L(iNX_X1,iX2,iX3,iX1) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X1, iX1, iX2, iX3 = ', iNX_X1, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_L(iNX_X1,iX2,iX3,iX1) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Surface term (X1) (Right)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1) + 1
      DO iNX_X1 = 1, nDOFX_X1

        IF( iErr_R(iNX_X1,iX2,iX3,iX1) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X1, iX1, iX2, iX3 = ', iNX_X1, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_R(iNX_X1,iX2,iX3,iX1) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: AlphaMiddle (X1)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1) + 1
      DO iNX_X1 = 1, nDOFX_X1

        IF( iErr_M(iNX_X1,iX2,iX3,iX1) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X1, iX1, iX2, iX3 = ', iNX_X1, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_M(iNX_X1,iX2,iX3,iX1) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Volume term (X1)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr_V(iNX,iX2,iX3,iX1) .NE. 0 )THEN

          WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_V(iNX,iX2,iX3,iX1) )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_ErrorCheck )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X1


  SUBROUTINE ComputeIncrement_Euler_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER  :: nK(3), nF_X2(3), nCF_K, nCF_F, nGF_F, nDF_F
    INTEGER  :: iNX, iNX_X2, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uCF_F_L(nCF), uCF_F_R(nCF)
    REAL(DP) :: Flux_X2_L(nCF), Flux_X2_R(nCF)
    REAL(DP) :: Flux_X2_F(nCF), Flux_X2_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L, P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_F          (nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    REAL(DP) :: uCF_L        (nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uCF_R        (nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)-1:iX_E0(2)+1)

    REAL(DP) :: uPF_L        (nDOFX_X2,nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uPF_R        (nDOFX_X2,nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uPF_K        (nDOFX,   nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)-1:iX_E0(2)+1)

    REAL(DP) :: uDF_L        (2           ,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: uDF_R        (2           ,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    REAL(DP) :: dU_X2        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)  )

    REAL(DP) :: Flux_X2_q    (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)  )

    REAL(DP) :: NumericalFlux(nDOFX_X2,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(3)  :iX_E0(3), &
                                           iX_B0(2)  :iX_E0(2)+1)

    INTEGER :: iErr_L(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(2):iX_E0(2)+1 )
    INTEGER :: iErr_R(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(2):iX_E0(2)+1 )
    INTEGER :: iErr_M(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(2):iX_E0(2)+1 )
    INTEGER :: iErr_V(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                 iX_B0(3):iX_E0(3), &
                                 iX_B0(2):iX_E0(2) )

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nF_X2 = nK + [0,1,0]           ! Number of X2 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nF_X2 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nF_X2 ) ! Number of Geometry Fields on Interfaces
    nDF_F = 2   * PRODUCT( nF_X2 ) ! Number of Diagnostic Fields on Interfaces

    ASSOCIATE( dX1 => MeshX(1) % Width, dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX3 ) &
    !$OMP MAP( alloc: uCF_L, uCF_R, uCF_K, &
    !$OMP             uPF_L, uPF_R, uPF_K, &
    !$OMP             uDF_L, uDF_R,        &
    !$OMP             G_K, G_F, dU_X2, Flux_X2_q, NumericalFlux, &
    !$OMP             iErr_L, iErr_R, iErr_M, iErr_V )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX3 ) &
    !$ACC CREATE(     uCF_L, uCF_R, uCF_K, &
    !$ACC             uPF_L, uPF_R, uPF_K, &
    !$ACC             uDF_L, uDF_R,        &
    !$ACC             G_K, G_F, dU_X2, Flux_X2_q, NumericalFlux, &
    !$ACC             iErr_L, iErr_R, iErr_M, iErr_V )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNX = 1, nDOFX

      G_K(iNX,iGF,iX1,iX3,iX2) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_F, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nGF_F, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X2 = 1, nDOFX_X2

      G_F         (iNX_X2,iGF_h_1,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_1,iX1,iX3,iX2), SqrtTiny )

      G_F         (iNX_X2,iGF_h_2,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_2,iX1,iX3,iX2), SqrtTiny )

      G_F         (iNX_X2,iGF_h_3,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_3,iX1,iX3,iX2), SqrtTiny )

      G_F         (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F         (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F         (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )

      G_F        (iNX_X2,iGF_SqrtGm,iX1,iX3,iX2) &
        = G_F    (iNX_X2,iGF_h_1   ,iX1,iX3,iX2) &
            * G_F(iNX_X2,iGF_h_2   ,iX1,iX3,iX2) &
            * G_F(iNX_X2,iGF_h_3   ,iX1,iX3,iX2)

      G_F         (iNX_X2,iGF_Alpha,iX1,iX3,iX2) &
        = MAX( G_F(iNX_X2,iGF_Alpha,iX1,iX3,iX2), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( uCF_K, U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      uCF_K(iNX,iCF,iX1,iX3,iX2) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDF_L(1,iX1,iX3,iX2) = D(1,iX1,iX2-1,iX3,iDF_Sh_X1)
      uDF_L(2,iX1,iX3,iX2) = D(1,iX1,iX2-1,iX3,iDF_Sh_X3)
      uDF_R(1,iX1,iX3,iX2) = D(1,iX1,iX2  ,iX3,iDF_Sh_X1)
      uDF_R(2,iX1,iX3,iX2) = D(1,iX1,iX2  ,iX3,iDF_Sh_X3)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nCF_F, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             uCF_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nCF_F, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             uCF_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F, uCF_L, uCF_R, uPF_L, uPF_R, &
    !$ACC          iErr_L, iErr_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X2 = 1, nDOFX_X2

      iErr_L(iNX_X2,iX1,iX3,iX2) = 0
      iErr_R(iNX_X2,iX1,iX3,iX2) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_L (iNX_X2,iCF_D       ,iX1,iX3,iX2), &
               uCF_L (iNX_X2,iCF_S1      ,iX1,iX3,iX2), &
               uCF_L (iNX_X2,iCF_S2      ,iX1,iX3,iX2), &
               uCF_L (iNX_X2,iCF_S3      ,iX1,iX3,iX2), &
               uCF_L (iNX_X2,iCF_E       ,iX1,iX3,iX2), &
               uCF_L (iNX_X2,iCF_Ne      ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_D       ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_E       ,iX1,iX3,iX2), &
               uPF_L (iNX_X2,iPF_Ne      ,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
               iErr_L(iNX_X2             ,iX1,iX3,iX2) )

      CALL ComputePrimitive_Euler &
             ( uCF_R (iNX_X2,iCF_D       ,iX1,iX3,iX2), &
               uCF_R (iNX_X2,iCF_S1      ,iX1,iX3,iX2), &
               uCF_R (iNX_X2,iCF_S2      ,iX1,iX3,iX2), &
               uCF_R (iNX_X2,iCF_S3      ,iX1,iX3,iX2), &
               uCF_R (iNX_X2,iCF_E       ,iX1,iX3,iX2), &
               uCF_R (iNX_X2,iCF_Ne      ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_D       ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_E       ,iX1,iX3,iX2), &
               uPF_R (iNX_X2,iPF_Ne      ,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
               G_F   (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
               iErr_R(iNX_X2             ,iX1,iX3,iX2) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X2_F, AlphaMns, AlphaPls, AlphaMdl )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$ACC          Flux_X2_F, AlphaMns, AlphaPls, AlphaMdl ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX1, dX3, &
    !$ACC          G_F, uCF_L, uCF_R, uPF_L, uPF_R, uDF_L, uDF_R, &
    !$ACC          NumericalFlux, WeightsX_X2, iErr_M )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X2_F, AlphaMns, AlphaPls, AlphaMdl )
#endif
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X2 = 1, nDOFX_X2

      iErr_M(iNX_X2,iX1,iX3,iX2) = 0

      DO iCF = 1, nCF

        uCF_F_L(iCF) = uCF_L(iNX_X2,iCF,iX1,iX3,iX2)
        uCF_F_R(iCF) = uCF_R(iNX_X2,iCF,iX1,iX3,iX2)

      END DO

      ! --- Left State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_L(iNX_X2,iPF_D ,iX1,iX3,iX2), &
               uPF_L(iNX_X2,iPF_E ,iX1,iX3,iX2), &
               uPF_L(iNX_X2,iPF_Ne,iX1,iX3,iX2), &
               P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(iNX_X2,iPF_D ,iX1,iX3,iX2), &
               uPF_L(iNX_X2,iPF_E ,iX1,iX3,iX2), &
               uPF_L(iNX_X2,iPF_Ne,iX1,iX3,iX2), &
               Cs_L )

      EigVals_L = &
        Eigenvalues_Euler &
          ( uPF_L(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
            Cs_L,                                   &
            G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
            uPF_L(iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
            uPF_L(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
            uPF_L(iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

      Flux_X2_L &
        = Flux_X2_Euler &
            ( uPF_L(iNX_X2,iPF_D       ,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_E       ,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_Ne      ,iX1,iX3,iX2), &
              P_L,                                    &
              G_F  (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

      ! --- Right State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_R(iNX_X2,iPF_D ,iX1,iX3,iX2), &
               uPF_R(iNX_X2,iPF_E ,iX1,iX3,iX2), &
               uPF_R(iNX_X2,iPF_Ne,iX1,iX3,iX2), &
               P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(iNX_X2,iPF_D ,iX1,iX3,iX2), &
               uPF_R(iNX_X2,iPF_E ,iX1,iX3,iX2), &
               uPF_R(iNX_X2,iPF_Ne,iX1,iX3,iX2), &
               Cs_R )

      EigVals_R = &
        Eigenvalues_Euler &
          ( uPF_R(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
            Cs_R,                                   &
            G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
            uPF_R(iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
            uPF_R(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
            uPF_R(iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
            G_F  (iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

      Flux_X2_R &
        = Flux_X2_Euler &
            ( uPF_R(iNX_X2,iPF_D       ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_V1      ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_V3      ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_E       ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_Ne      ,iX1,iX3,iX2), &
              P_R,                                    &
              G_F  (iNX_X2,iGF_Gm_dd_11,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Gm_dd_33,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2) )

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uCF_F_L  (iCF_D), uCF_F_L  (iCF_S2), uCF_F_L  (iCF_E), &
              Flux_X2_L(iCF_D), Flux_X2_L(iCF_S2), Flux_X2_L(iCF_E), &
              uCF_F_R  (iCF_D), uCF_F_R  (iCF_S2), uCF_F_R  (iCF_E), &
              Flux_X2_R(iCF_D), Flux_X2_R(iCF_S2), Flux_X2_R(iCF_E), &
              G_F(iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), AlphaPls, AlphaMns, &
              G_F(iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
              G_F(iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2), &
              iErr_M(iNX_X2,iX1,iX3,iX2) )

      Flux_X2_F &
        = NumericalFlux_Euler_X2 &
            ( uCF_F_L, uCF_F_R, Flux_X2_L, Flux_X2_R, &
              AlphaPls, AlphaMns, AlphaMdl,           &
              G_F  (iNX_X2,iGF_Gm_dd_22,iX1,iX3,iX2), &
              uPF_L(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
              uPF_R(iNX_X2,iPF_V2      ,iX1,iX3,iX2), &
              P_L, P_R,                               &
              G_F  (iNX_X2,iGF_Alpha   ,iX1,iX3,iX2), &
              G_F  (iNX_X2,iGF_Beta_2  ,iX1,iX3,iX2), &
              uDF_L(1                  ,iX1,iX3,iX2), &
              uDF_R(1                  ,iX1,iX3,iX2), &
              uDF_L(2                  ,iX1,iX3,iX2), &
              uDF_R(2                  ,iX1,iX3,iX2) )

      DO iCF = 1, nCF

        NumericalFlux    (iNX_X2,iCF       ,iX1,iX3,iX2) &
          = Flux_X2_F(iCF)                               &
              * G_F      (iNX_X2,iGF_Alpha ,iX1,iX3,iX2) &
              * G_F      (iNX_X2,iGF_SqrtGm,iX1,iX3,iX2) &
              * dX1(iX1) * dX3(iX3) * WeightsX_X2(iNX_X2)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X2, + One, LX_X2_Dn, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)), nDOFX_X2, Zero, &
             dU_X2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X2, - One, LX_X2_Up, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One, &
             dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, uCF_K, uPF_K, iErr_V )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iErr_V(iNX,iX1,iX3,iX2) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_K (iNX,iCF_D       ,iX1,iX3,iX2), &
               uCF_K (iNX,iCF_S1      ,iX1,iX3,iX2), &
               uCF_K (iNX,iCF_S2      ,iX1,iX3,iX2), &
               uCF_K (iNX,iCF_S3      ,iX1,iX3,iX2), &
               uCF_K (iNX,iCF_E       ,iX1,iX3,iX2), &
               uCF_K (iNX,iCF_Ne      ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_D       ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_V1      ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_V2      ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_V3      ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_E       ,iX1,iX3,iX2), &
               uPF_K (iNX,iPF_Ne      ,iX1,iX3,iX2), &
               G_K   (iNX,iGF_Gm_dd_11,iX1,iX3,iX2), &
               G_K   (iNX,iGF_Gm_dd_22,iX1,iX3,iX2), &
               G_K   (iNX,iGF_Gm_dd_33,iX1,iX3,iX2), &
               iErr_V(iNX             ,iX1,iX3,iX2) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_K, P_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRIVATE( Flux_X2_K, P_K ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX1, dX3, &
    !$ACC          G_K, uCF_K, uPF_K, Flux_X2_q, WeightsX_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X2_K, P_K )
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      CALL ComputePressureFromPrimitive &
             ( uPF_K(iNX,iPF_D ,iX1,iX3,iX2), &
               uPF_K(iNX,iPF_E ,iX1,iX3,iX2), &
               uPF_K(iNX,iPF_Ne,iX1,iX3,iX2), &
               P_K )

      Flux_X2_K &
        = Flux_X2_Euler &
          ( uPF_K(iNX,iPF_D       ,iX1,iX3,iX2), &
            uPF_K(iNX,iPF_V1      ,iX1,iX3,iX2), &
            uPF_K(iNX,iPF_V2      ,iX1,iX3,iX2), &
            uPF_K(iNX,iPF_V3      ,iX1,iX3,iX2), &
            uPF_K(iNX,iPF_E       ,iX1,iX3,iX2), &
            uPF_K(iNX,iPF_Ne      ,iX1,iX3,iX2), &
            P_K,                                 &
            G_K  (iNX,iGF_Gm_dd_11,iX1,iX3,iX2), &
            G_K  (iNX,iGF_Gm_dd_22,iX1,iX3,iX2), &
            G_K  (iNX,iGF_Gm_dd_33,iX1,iX3,iX2), &
            G_K  (iNX,iGF_Alpha,   iX1,iX3,iX2), &
            G_K  (iNX,iGF_Beta_2,  iX1,iX3,iX2) )

      DO iCF = 1, nCF

        Flux_X2_q(iNX,iCF,iX1,iX3,iX2) = Flux_X2_K(iCF)

        Flux_X2_q    (iNX,iCF       ,iX1,iX3,iX2) &
          = Flux_X2_q(iNX,iCF       ,iX1,iX3,iX2) &
              * G_K  (iNX,iGF_Alpha ,iX1,iX3,iX2) &
              * G_K  (iNX,iGF_SqrtGm,iX1,iX3,iX2) &
              * dX1(iX1) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX2_q, nDOFX, &
             Flux_X2_q, nDOFX, One, dU_X2, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X2, dU )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X2(iNX,iCF,iX1,iX3,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X2 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X2 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X2)', MAXLOC(dU_X2)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X2)', MAXVAL(dU_X2)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX3, &
    !$OMP               uCF_L, uCF_R, uCF_K, &
    !$OMP               uPF_L, uPF_R, uPF_K, &
    !$OMP               uDF_L, uDF_R,        &
    !$OMP               G_K, G_F, dU_X2, Flux_X2_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX3, &
    !$ACC               uCF_L, uCF_R, uCF_K, &
    !$ACC               uPF_L, uPF_R, uPF_K, &
    !$ACC               uDF_L, uDF_R,        &
    !$ACC               G_K, G_F, dU_X2, Flux_X2_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_DG_ErrorCheck )

    IF( ANY( iErr_L .NE. 0 ) .OR. &
        ANY( iErr_R .NE. 0 ) .OR. &
        ANY( iErr_M .NE. 0 ) .OR. &
        ANY( iErr_V .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR: Surface term (X2) (Left)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2) + 1
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X2 = 1, nDOFX_X2

        IF( iErr_L(iNX_X2,iX1,iX3,iX2) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X2, iX1, iX2, iX3 = ', iNX_X2, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_L(iNX_X2,iX1,iX3,iX2) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Surface term (X2) (Right)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2) + 1
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X2 = 1, nDOFX_X2

        IF( iErr_R(iNX_X2,iX1,iX3,iX2) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X2, iX1, iX2, iX3 = ', iNX_X2, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_R(iNX_X2,iX1,iX3,iX2) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: AlphaMiddle (X2)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2) + 1
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X2 = 1, nDOFX_X2

        IF( iErr_M(iNX_X2,iX1,iX3,iX2) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X2, iX1, iX2, iX3 = ', iNX_X2, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_M(iNX_X2,iX1,iX3,iX2) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Volume term (X2)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr_V(iNX,iX1,iX3,iX2) .NE. 0 )THEN

          WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_V(iNX,iX1,iX3,iX2) )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_ErrorCheck )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X2


  SUBROUTINE ComputeIncrement_Euler_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER  :: nK(3), nF_X3(3), nCF_K, nCF_F, nGF_F, nDF_F
    INTEGER  :: iNX, iNX_X3, iX1, iX2, iX3, iCF, iGF
    REAL(DP) :: AlphaMns, AlphaPls, AlphaMdl
    REAL(DP) :: uCF_F_L(nCF), uCF_F_R(nCF)
    REAL(DP) :: Flux_X3_L(nCF), Flux_X3_R(nCF)
    REAL(DP) :: Flux_X3_F(nCF), Flux_X3_K(nCF)
    REAL(DP) :: EigVals_L(nCF), EigVals_R(nCF)
    REAL(DP) :: P_L, P_R, P_K
    REAL(DP) :: Cs_L, Cs_R

    REAL(DP) :: G_K          (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_F          (nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    REAL(DP) :: uCF_L        (nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uCF_R        (nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uCF_K        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)-1:iX_E0(3)+1)

    REAL(DP) :: uPF_L        (nDOFX_X3,nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uPF_R        (nDOFX_X3,nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uPF_K        (nDOFX,   nPF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)-1:iX_E0(3)+1)

    REAL(DP) :: uDF_L        (2           ,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: uDF_R        (2           ,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    REAL(DP) :: dU_X3        (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)  )

    REAL(DP) :: Flux_X3_q    (nDOFX,   nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)  )

    REAL(DP) :: NumericalFlux(nDOFX_X3,nCF,iX_B0(1)  :iX_E0(1), &
                                           iX_B0(2)  :iX_E0(2), &
                                           iX_B0(3)  :iX_E0(3)+1)

    INTEGER :: iErr_L(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3)+1 )
    INTEGER :: iErr_R(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3)+1 )
    INTEGER :: iErr_M(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3)+1 )
    INTEGER :: iErr_V(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3) )

    IF( iX_E0(3) .EQ. iX_B0(3) ) RETURN

    nK    = iX_E0 - iX_B0 + 1      ! Number of Elements per Spatial Dimension
    nF_X3 = nK + [0,0,1]           ! Number of X3 Faces per Spatial Dimension
    nCF_K = nCF * PRODUCT( nK )    ! Number of Fluid Fields in Domain
    nCF_F = nCF * PRODUCT( nF_X3 ) ! Number of Fluid Fields on Interfaces
    nGF_F = nGF * PRODUCT( nF_X3 ) ! Number of Geometry Fields on Interfaces
    nDF_F = 2   * PRODUCT( nF_X3 ) ! Number of Diagnostic Fields on Interfaces

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2 ) &
    !$OMP MAP( alloc: uCF_L, uCF_R, uCF_K, &
    !$OMP             uPF_L, uPF_R, uPF_K, &
    !$OMP             uDF_L, uDF_R,        &
    !$OMP             G_K, G_F, dU_X3, Flux_X3_q, NumericalFlux, &
    !$OMP             iErr_L, iErr_R, iErr_M, iErr_V )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2 ) &
    !$ACC CREATE(     uCF_L, uCF_R, uCF_K, &
    !$ACC             uPF_L, uPF_R, uPF_K, &
    !$ACC             uDF_L, uDF_R,        &
    !$ACC             G_K, G_F, dU_X3, Flux_X3_q, NumericalFlux, &
    !$ACC             iErr_L, iErr_R, iErr_M, iErr_V )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF
    DO iNX = 1, nDOFX

      G_K(iNX,iGF,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_F, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             G_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nGF_F, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Half, &
             G_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X3 = 1, nDOFX_X3

      G_F         (iNX_X3,iGF_h_1,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_1,iX1,iX2,iX3), SqrtTiny )

      G_F         (iNX_X3,iGF_h_2,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_2,iX1,iX2,iX3), SqrtTiny )

      G_F         (iNX_X3,iGF_h_3,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_3,iX1,iX2,iX3), SqrtTiny )

      G_F         (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F         (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F         (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )

      G_F        (iNX_X3,iGF_SqrtGm,iX1,iX2,iX3) &
        = G_F    (iNX_X3,iGF_h_1   ,iX1,iX2,iX3) &
            * G_F(iNX_X3,iGF_h_2   ,iX1,iX2,iX3) &
            * G_F(iNX_X3,iGF_h_3   ,iX1,iX2,iX3)

      G_F         (iNX_X3,iGF_Alpha,iX1,iX2,iX3) &
        = MAX( G_F(iNX_X3,iGF_Alpha,iX1,iX2,iX3), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, uCF_K, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      uCF_K(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, uDF_L, uDF_R, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDF_L(1,iX1,iX2,iX3) = D(1,iX1,iX2,iX3-1,iDF_Sh_X1)
      uDF_L(2,iX1,iX2,iX3) = D(1,iX1,iX2,iX3-1,iDF_Sh_X2)
      uDF_R(1,iX1,iX2,iX3) = D(1,iX1,iX2,iX3  ,iDF_Sh_X1)
      uDF_R(2,iX1,iX2,iX3) = D(1,iX1,iX2,iX3  ,iDF_Sh_X2)

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nCF_F, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             uCF_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             uCF_L(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nCF_F, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             uCF_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Zero, &
             uCF_R(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_F, uCF_L, uCF_R, uPF_L, uPF_R, &
    !$ACC          iErr_L, iErr_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X3 = 1, nDOFX_X3

      iErr_L(iNX_X3,iX1,iX2,iX3) = 0
      iErr_R(iNX_X3,iX1,iX2,iX3) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_L (iNX_X3,iCF_D       ,iX1,iX2,iX3), &
               uCF_L (iNX_X3,iCF_S1      ,iX1,iX2,iX3), &
               uCF_L (iNX_X3,iCF_S2      ,iX1,iX2,iX3), &
               uCF_L (iNX_X3,iCF_S3      ,iX1,iX2,iX3), &
               uCF_L (iNX_X3,iCF_E       ,iX1,iX2,iX3), &
               uCF_L (iNX_X3,iCF_Ne      ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_D       ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_E       ,iX1,iX2,iX3), &
               uPF_L (iNX_X3,iPF_Ne      ,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
               iErr_L(iNX_X3             ,iX1,iX2,iX3) )

      CALL ComputePrimitive_Euler &
             ( uCF_R (iNX_X3,iCF_D       ,iX1,iX2,iX3), &
               uCF_R (iNX_X3,iCF_S1      ,iX1,iX2,iX3), &
               uCF_R (iNX_X3,iCF_S2      ,iX1,iX2,iX3), &
               uCF_R (iNX_X3,iCF_S3      ,iX1,iX2,iX3), &
               uCF_R (iNX_X3,iCF_E       ,iX1,iX2,iX3), &
               uCF_R (iNX_X3,iCF_Ne      ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_D       ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_E       ,iX1,iX2,iX3), &
               uPF_R (iNX_X3,iPF_Ne      ,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
               G_F   (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
               iErr_R(iNX_X3             ,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X3_F, AlphaMns, AlphaPls, AlphaMdl )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$ACC          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$ACC          Flux_X3_F, AlphaMns, AlphaPls, AlphaMdl ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX1, dX2, &
    !$ACC          G_F, uCF_L, uCF_R, uPF_L, uPF_R, uDF_L, uDF_R, &
    !$ACC          NumericalFlux, WeightsX_X3, iErr_M )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, P_L, P_R, Cs_L, Cs_R, &
    !$OMP          EigVals_L, EigVals_R, uCF_F_L, uCF_F_R, &
    !$OMP          Flux_X3_F, AlphaMns, AlphaPls, AlphaMdl )
#endif
    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX_X3 = 1, nDOFX_X3

      iErr_M(iNX_X3,iX1,iX2,iX3) = 0

      DO iCF = 1, nCF

        uCF_F_L(iCF) = uCF_L(iNX_X3,iCF,iX1,iX2,iX3)
        uCF_F_R(iCF) = uCF_R(iNX_X3,iCF,iX1,iX2,iX3)

      END DO

      ! --- Left State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_L(iNX_X3,iPF_D ,iX1,iX2,iX3), &
               uPF_L(iNX_X3,iPF_E ,iX1,iX2,iX3), &
               uPF_L(iNX_X3,iPF_Ne,iX1,iX2,iX3), &
               P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(iNX_X3,iPF_D ,iX1,iX2,iX3), &
               uPF_L(iNX_X3,iPF_E ,iX1,iX2,iX3), &
               uPF_L(iNX_X3,iPF_Ne,iX1,iX2,iX3), &
               Cs_L )

      EigVals_L &
        = Eigenvalues_Euler &
            ( uPF_L(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              Cs_L,                                   &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

      Flux_X3_L &
        = Flux_X3_Euler &
            ( uPF_L(iNX_X3,iPF_D       ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_E       ,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_Ne      ,iX1,iX2,iX3), &
              P_L,                                    &
              G_F  (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

      ! --- Right State Primitive, etc. ---

      CALL ComputePressureFromPrimitive &
             ( uPF_R(iNX_X3,iPF_D ,iX1,iX2,iX3), &
               uPF_R(iNX_X3,iPF_E ,iX1,iX2,iX3), &
               uPF_R(iNX_X3,iPF_Ne,iX1,iX2,iX3), &
               P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(iNX_X3,iPF_D ,iX1,iX2,iX3), &
               uPF_R(iNX_X3,iPF_E ,iX1,iX2,iX3), &
               uPF_R(iNX_X3,iPF_Ne,iX1,iX2,iX3), &
               Cs_R )

      EigVals_R &
        = Eigenvalues_Euler &
            ( uPF_R(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              Cs_R,                                   &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

      Flux_X3_R &
        = Flux_X3_Euler &
            ( uPF_R(iNX_X3,iPF_D       ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V1      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V2      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_E       ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_Ne      ,iX1,iX2,iX3), &
              P_R,                                    &
              G_F  (iNX_X3,iGF_Gm_dd_11,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_22,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3) )

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      AlphaMdl &
        = AlphaMiddle_Euler &
            ( uCF_F_L  (iCF_D), uCF_F_L  (iCF_S3), uCF_F_L  (iCF_E), &
              Flux_X3_L(iCF_D), Flux_X3_L(iCF_S3), Flux_X3_L(iCF_E), &
              uCF_F_R  (iCF_D), uCF_F_R  (iCF_S3), uCF_F_R  (iCF_E), &
              Flux_X3_R(iCF_D), Flux_X3_R(iCF_S3), Flux_X3_R(iCF_E), &
              G_F(iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), AlphaPls, AlphaMns, &
              G_F(iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F(iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3), &
              iErr_M(iNX_X3,iX1,iX2,iX3) )

      Flux_X3_F &
        = NumericalFlux_Euler_X3 &
            ( uCF_F_L, uCF_F_R, Flux_X3_L, Flux_X3_R,   &
              AlphaPls, AlphaMns, AlphaMdl,             &
              G_F  (iNX_X3,iGF_Gm_dd_33,iX1,iX2,iX3), &
              uPF_L(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              uPF_R(iNX_X3,iPF_V3      ,iX1,iX2,iX3), &
              P_L, P_R,   &
              G_F  (iNX_X3,iGF_Alpha   ,iX1,iX2,iX3), &
              G_F  (iNX_X3,iGF_Beta_3  ,iX1,iX2,iX3), &
              uDF_L(1                  ,iX1,iX2,iX3), &
              uDF_R(1                  ,iX1,iX2,iX3), &
              uDF_L(2                  ,iX1,iX2,iX3), &
              uDF_R(2                  ,iX1,iX2,iX3) )

      DO iCF = 1, nCF

        NumericalFlux    (iNX_X3,iCF       ,iX1,iX2,iX3) &
          = Flux_X3_F(iCF)                               &
              * G_F      (iNX_X3,iGF_Alpha ,iX1,iX2,iX3) &
              * G_F      (iNX_X3,iGF_SqrtGm,iX1,iX2,iX3) &
              * dX1(iX1) * dX2(iX2) * WeightsX_X3(iNX_X3)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    ! --- Surface Contribution ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X3, + One, LX_X3_Dn, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)), nDOFX_X3, Zero, &
             dU_X3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX_X3, - One, LX_X3_Up, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One, &
             dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K, uCF_K, uPF_K, iErr_V )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iErr_V(iNX,iX1,iX2,iX3) = 0

      CALL ComputePrimitive_Euler &
             ( uCF_K (iNX,iCF_D       ,iX1,iX2,iX3), &
               uCF_K (iNX,iCF_S1      ,iX1,iX2,iX3), &
               uCF_K (iNX,iCF_S2      ,iX1,iX2,iX3), &
               uCF_K (iNX,iCF_S3      ,iX1,iX2,iX3), &
               uCF_K (iNX,iCF_E       ,iX1,iX2,iX3), &
               uCF_K (iNX,iCF_Ne      ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_D       ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_V1      ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_V2      ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_V3      ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_E       ,iX1,iX2,iX3), &
               uPF_K (iNX,iPF_Ne      ,iX1,iX2,iX3), &
               G_K   (iNX,iGF_Gm_dd_11,iX1,iX2,iX3), &
               G_K   (iNX,iGF_Gm_dd_22,iX1,iX2,iX3), &
               G_K   (iNX,iGF_Gm_dd_33,iX1,iX2,iX3), &
               iErr_V(iNX             ,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_K, P_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4), &
    !$ACC PRIVATE( Flux_X3_K, P_K ) &
    !$ACC PRESENT( iX_B0, iX_E0, dX1, dX2, &
    !$ACC          G_K, uCF_K, uPF_K, Flux_X3_q, WeightsX_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Flux_X3_K, P_K )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      CALL ComputePressureFromPrimitive &
             ( uPF_K(iNX,iPF_D ,iX1,iX2,iX3), &
               uPF_K(iNX,iPF_E ,iX1,iX2,iX3), &
               uPF_K(iNX,iPF_Ne,iX1,iX2,iX3), &
               P_K )

      Flux_X3_K &
        = Flux_X3_Euler &
          ( uPF_K(iNX,iPF_D       ,iX1,iX2,iX3), &
            uPF_K(iNX,iPF_V1      ,iX1,iX2,iX3), &
            uPF_K(iNX,iPF_V2      ,iX1,iX2,iX3), &
            uPF_K(iNX,iPF_V3      ,iX1,iX2,iX3), &
            uPF_K(iNX,iPF_E       ,iX1,iX2,iX3), &
            uPF_K(iNX,iPF_Ne      ,iX1,iX2,iX3), &
            P_K,                                 &
            G_K  (iNX,iGF_Gm_dd_11,iX1,iX2,iX3), &
            G_K  (iNX,iGF_Gm_dd_22,iX1,iX2,iX3), &
            G_K  (iNX,iGF_Gm_dd_33,iX1,iX2,iX3), &
            G_K  (iNX,iGF_Alpha,   iX1,iX2,iX3), &
            G_K  (iNX,iGF_Beta_3,  iX1,iX2,iX3) )

      DO iCF = 1, nCF

        Flux_X3_q(iNX,iCF,iX1,iX2,iX3) = Flux_X3_K(iCF)

        Flux_X3_q    (iNX,iCF       ,iX1,iX2,iX3) &
          = Flux_X3_q(iNX,iCF       ,iX1,iX2,iX3) &
              * G_K  (iNX,iGF_Alpha ,iX1,iX2,iX3) &
              * G_K  (iNX,iGF_SqrtGm,iX1,iX2,iX3) &
              * dX1(iX1) * dX2(iX2) * WeightsX_q(iNX)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

    ! --- Contribution from Volume ---

    CALL TimersStart_Euler( Timer_Euler_DG_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nCF_K, nDOFX, One, dLXdX3_q, nDOFX, &
             Flux_X3_q, nDOFX, One, dU_X3, nDOFX )

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dU_X3, dU )
#elif defined(THORNADO_OMP)
    !$ACC PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      dU    (iNX,iX1,iX2,iX3,iCF) &
        = dU(iNX,iX1,iX2,iX3,iCF) &
            + dU_X3(iNX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X3 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X3 )
#endif
    WRITE(*,'(A20,7I4)')     'MAXLOC(dU_X3)', MAXLOC(dU_X3)
    WRITE(*,'(A20,ES23.15)') 'MAXVAL(dU_X3)', MAXVAL(dU_X3)
#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, &
    !$OMP               uCF_L, uCF_R, uCF_K, &
    !$OMP               uPF_L, uPF_R, uPF_K, &
    !$OMP               uDF_L, uDF_R,        &
    !$OMP               G_K, G_F, dU_X3, Flux_X3_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr_L, iErr_R, iErr_M, iErr_V ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, dX1, dX2, &
    !$ACC               uCF_L, uCF_R, uCF_K, &
    !$ACC               uPF_L, uPF_R, uPF_K, &
    !$ACC               uDF_L, uDF_R,        &
    !$ACC               G_K, G_F, dU_X3, Flux_X3_q, NumericalFlux )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_DG_ErrorCheck )

    IF( ANY( iErr_L .NE. 0 ) .OR. &
        ANY( iErr_R .NE. 0 ) .OR. &
        ANY( iErr_M .NE. 0 ) .OR. &
        ANY( iErr_V .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR: Surface term (X3) (Left)'

      DO iX3 = iX_B0(3), iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X3 = 1, nDOFX_X3

        IF( iErr_L(iNX_X3,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X3, iX1, iX2, iX3 = ', iNX_X3, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_L(iNX_X3,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Surface term (X3) (Right)'

      DO iX3 = iX_B0(3), iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X3 = 1, nDOFX_X3

        IF( iErr_R(iNX_X3,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X3, iX1, iX2, iX3 = ', iNX_X3, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_R(iNX_X3,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: AlphaMiddle (X3)'

      DO iX3 = iX_B0(3), iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX_X3 = 1, nDOFX_X3

        IF( iErr_M(iNX_X3,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iNX_X3, iX1, iX2, iX3 = ', iNX_X3, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_M(iNX_X3,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

      WRITE(*,*) 'ERROR: Volume term (X3)'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr_V(iNX,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr_V(iNX,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_ErrorCheck )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Euler_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

#ifdef HYDRO_RELATIVISTIC

!!$    CALL ComputeIncrement_Geometry_Relativistic_CPU &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Geometry_Relativistic_GPU &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

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


  SUBROUTINE ComputeIncrement_Geometry_Relativistic_CPU &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNX, iDim, jDim, iErr(nDOFX)
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh1dX1(nDOFX), dh2dX1(nDOFX), dh3dX1(nDOFX), &
                dh1dX2(nDOFX), dh2dX2(nDOFX), dh3dX2(nDOFX), &
                dh1dX3(nDOFX), dh2dX3(nDOFX), dh3dX3(nDOFX)
    REAL(DP) :: db1dX1(nDOFX), db2dX1(nDOFX), db3dX1(nDOFX), &
                db1dX2(nDOFX), db2dX2(nDOFX), db3dX2(nDOFX), &
                db1dX3(nDOFX), db2dX3(nDOFX), db3dX3(nDOFX)
    REAL(DP) :: dadx1(nDOFX), dadx2(nDOFX), dadx3(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF), uPF_K(nDOFX,nPF), G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF), &
                G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF), &
                G_P_X3(nDOFX,nGF), G_N_X3(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF), &
                G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF), &
                G_X3_Dn(nDOFX_X3,nGF), G_X3_Up(nDOFX_X3,nGF)

    REAL(DP) :: EnergyDensitySourceTerms(7,nDOFX,iX_B0(1):iX_E0(1), &
                                                 iX_B0(2):iX_E0(2), &
                                                 iX_B0(3):iX_E0(3))

    REAL(DP) :: DivGridVolume      (nDOFX)
    REAL(DP) :: PressureTensorTrace(nDOFX)
    REAL(DP) :: PressureTensor     (nDOFX,3,3)
    REAL(DP) :: Xij                (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X1   (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X2   (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X3   (nDOFX,3,3)
    REAL(DP) :: Christoffel_X1     (nDOFX,3,3)
    REAL(DP) :: Christoffel_X2     (nDOFX,3,3)
    REAL(DP) :: Christoffel_X3     (nDOFX,3,3)

    REAL(DP) :: GradPsi (nDOFX)
    REAL(DP) :: GradPsiF(nDOFX)
    REAL(DP) :: X1      (nDOFX)

    dadx1  = Zero
    dadx2  = Zero
    dadx3  = Zero
    dh1dX1 = Zero
    dh2dX1 = Zero
    dh3dX1 = Zero
    dh1dX2 = Zero
    dh2dX2 = Zero
    dh3dX2 = Zero
    dh1dX3 = Zero
    dh2dX3 = Zero
    dh3dX3 = Zero
    db1dX1 = Zero
    db2dX1 = Zero
    db3dX1 = Zero
    db1dX2 = Zero
    db2dX2 = Zero
    db3dX2 = Zero
    db1dX3 = Zero
    db2dX3 = Zero
    db3dX3 = Zero
    PressureTensor = Zero

    GradPsi = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX3 = MeshX(3) % Width(iX3)
      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      IF     ( nDimsX .EQ. 2 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)

        END DO

      ELSE IF( nDimsX .EQ. 3 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)
          G_P_X3(:,iGF) = G(:,iX1,iX2,iX3-1,iGF)
          G_N_X3(:,iGF) = G(:,iX1,iX2,iX3+1,iGF)

        END DO

      END IF

      CALL ComputePrimitive_Euler &
           ( uCF_K(:,iCF_D ),     &
             uCF_K(:,iCF_S1),     &
             uCF_K(:,iCF_S2),     &
             uCF_K(:,iCF_S3),     &
             uCF_K(:,iCF_E ),     &
             uCF_K(:,iCF_Ne),     &
             uPF_K(:,iPF_D ),     &
             uPF_K(:,iPF_V1),     &
             uPF_K(:,iPF_V2),     &
             uPF_K(:,iPF_V3),     &
             uPF_K(:,iPF_E ),     &
             uPF_K(:,iPF_Ne),     &
             G_K(:,iGF_Gm_dd_11), &
             G_K(:,iGF_Gm_dd_22), &
             G_K(:,iGF_Gm_dd_33), &
             iErr )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      ! --- Compute P^{ij} ---

      DO iNX = 1, nDOFX

        PressureTensor(:,1,1) &
          = ( uPF_K(:,iPF_V1) * uCF_K(:,iCF_S1) + P_K ) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,1,2) &
          = uPF_K(:,iPF_V1) * uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,1,3) &
          = uPF_K(:,iPF_V1) * uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33)

        PressureTensor(:,2,1) &
          = uPF_K(:,iPF_V2) * uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,2,2) &
          = ( uPF_K(:,iPF_V2) * uCF_K(:,iCF_S2) + P_K ) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,2,3) &
          = uPF_K(:,iPF_V2) * uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33)

        PressureTensor(:,3,1) &
          = uPF_K(:,iPF_V3) * uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,3,2) &
          = uPF_K(:,iPF_V3) * uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,3,3) &
          = ( uPF_K(:,iPF_V3) * uCF_K(:,iCF_S3) + P_K ) / G_K(:,iGF_Gm_dd_33)

      END DO

      ! --- Redundant calculation. Will modify code so this can be
      !     removed later ---

      DO iNX = 1, nDOFX

        Stress(iNX,:) &
          = StressTensor_Diagonal_Euler &
              ( uCF_K(iNX,iCF_S1), &
                uCF_K(iNX,iCF_S2), &
                uCF_K(iNX,iCF_S3), &
                uPF_K(iNX,iPF_V1), &
                uPF_K(iNX,iPF_V2), &
                uPF_K(iNX,iPF_V3), &
                P_K  (iNX) )

      END DO

      ! --- Scale factor derivatives wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_h_1), &
               G_K    (:,iGF_h_1), &
               G_N_X1 (:,iGF_h_1), &
               G_X1_Dn(:,iGF_h_1), &
               G_X1_Up(:,iGF_h_1) )

      G_X1_Dn(:,iGF_h_1) = MAX( G_X1_Dn(:,iGF_h_1), SqrtTiny )
      G_X1_Up(:,iGF_h_1) = MAX( G_X1_Up(:,iGF_h_1), SqrtTiny )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_h_2), &
               G_K    (:,iGF_h_2), &
               G_N_X1 (:,iGF_h_2), &
               G_X1_Dn(:,iGF_h_2), &
               G_X1_Up(:,iGF_h_2) )

      G_X1_Dn(:,iGF_h_2) = MAX( G_X1_Dn(:,iGF_h_2), SqrtTiny )
      G_X1_Up(:,iGF_h_2) = MAX( G_X1_Up(:,iGF_h_2), SqrtTiny )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_h_3), &
               G_K    (:,iGF_h_3), &
               G_N_X1 (:,iGF_h_3), &
               G_X1_Dn(:,iGF_h_3), &
               G_X1_Up(:,iGF_h_3) )

      G_X1_Dn(:,iGF_h_3) = MAX( G_X1_Dn(:,iGF_h_3), SqrtTiny )
      G_X1_Up(:,iGF_h_3) = MAX( G_X1_Up(:,iGF_h_3), SqrtTiny )

      ! --- Differentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_1), &
               G_X1_Dn(:,iGF_h_1), &
               G_K    (:,iGF_h_1), &
               dh1dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_2), &
               G_X1_Dn(:,iGF_h_2), &
               G_K    (:,iGF_h_2), &
               dh2dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_3), &
               G_X1_Dn(:,iGF_h_3), &
               G_K    (:,iGF_h_3), &
               dh3dX1 )

      ! --- Shift vector derivative wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_Beta_1), &
               G_K    (:,iGF_Beta_1), &
               G_N_X1 (:,iGF_Beta_1), &
               G_X1_Dn(:,iGF_Beta_1), &
               G_X1_Up(:,iGF_Beta_1) )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_Beta_2), &
               G_K    (:,iGF_Beta_2), &
               G_N_X1 (:,iGF_Beta_2), &
               G_X1_Dn(:,iGF_Beta_2), &
               G_X1_Up(:,iGF_Beta_2) )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_Beta_3), &
               G_K    (:,iGF_Beta_3), &
               G_N_X1 (:,iGF_Beta_3), &
               G_X1_Dn(:,iGF_Beta_3), &
               G_X1_Up(:,iGF_Beta_3) )

      ! --- Differentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_1), &
               G_X1_Dn(:,iGF_Beta_1), &
               G_K    (:,iGF_Beta_1), &
               db1dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_2), &
               G_X1_Dn(:,iGF_Beta_2), &
               G_K    (:,iGF_Beta_2), &
               db2dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_3), &
               G_X1_Dn(:,iGF_Beta_3), &
               G_K    (:,iGF_Beta_3), &
               db3dX1 )

      ! --- Lapse function derivative wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_Alpha), &
               G_K    (:,iGF_Alpha), &
               G_N_X1 (:,iGF_Alpha), &
               G_X1_Dn(:,iGF_Alpha), &
               G_X1_Up(:,iGF_Alpha) )

      G_X1_Dn(:,iGF_Alpha) = MAX( G_X1_Dn(:,iGF_Alpha), SqrtTiny )
      G_X1_Up(:,iGF_Alpha) = MAX( G_X1_Up(:,iGF_Alpha), SqrtTiny )

      ! --- Diffentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1,dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Alpha), &
               G_X1_Dn(:,iGF_Alpha), &
               G_K    (:,iGF_Alpha), &
               dadX1 )

      ! --- Momentum density source term (S1) ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1)                                &
            + G_K(:,iGF_Alpha)                                    &
                * ( ( Stress(:,1) * dh1dX1 ) / G_K(:,iGF_h_1)     &
                    + ( Stress(:,2) * dh2dX1 ) / G_K(:,iGF_h_2)   &
                    + ( Stress(:,3) * dh3dX1 ) / G_K(:,iGF_h_3) ) &
            + uCF_K(:,iCF_S1) * db1dX1 &
                + uCF_K(:,iCF_S2) * db2dX1 &
                + uCF_K(:,iCF_S3) * db3dX1 &
            - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx1

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale factor derivatives wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_h_1), &
                 G_K    (:,iGF_h_1), &
                 G_N_X2 (:,iGF_h_1), &
                 G_X2_Dn(:,iGF_h_1), &
                 G_X2_Up(:,iGF_h_1) )

        G_X2_Dn(:,iGF_h_1) = MAX( G_X2_Dn(:,iGF_h_1), SqrtTiny )
        G_X2_Up(:,iGF_h_1) = MAX( G_X2_Up(:,iGF_h_1), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_h_2), &
                 G_K    (:,iGF_h_2), &
                 G_N_X2 (:,iGF_h_2), &
                 G_X2_Dn(:,iGF_h_2), &
                 G_X2_Up(:,iGF_h_2) )

        G_X2_Dn(:,iGF_h_2) = MAX( G_X2_Dn(:,iGF_h_2), SqrtTiny )
        G_X2_Up(:,iGF_h_2) = MAX( G_X2_Up(:,iGF_h_2), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_h_3), &
                 G_K    (:,iGF_h_3), &
                 G_N_X2 (:,iGF_h_3), &
                 G_X2_Dn(:,iGF_h_3), &
                 G_X2_Up(:,iGF_h_3) )

        G_X2_Dn(:,iGF_h_3) = MAX( G_X2_Dn(:,iGF_h_3), SqrtTiny )
        G_X2_Up(:,iGF_h_3) = MAX( G_X2_Up(:,iGF_h_3), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_1), &
                 G_X2_Dn(:,iGF_h_1), &
                 G_K    (:,iGF_h_1), &
                 dh1dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_2), &
                 G_X2_Dn(:,iGF_h_2), &
                 G_K    (:,iGF_h_2), &
                 dh2dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_3), &
                 G_X2_Dn(:,iGF_h_3), &
                 G_K    (:,iGF_h_3), &
                 dh3dX2 )

        ! --- Shift vector derivatives wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_Beta_1), &
                 G_K    (:,iGF_Beta_1), &
                 G_N_X2 (:,iGF_Beta_1), &
                 G_X2_Dn(:,iGF_Beta_1), &
                 G_X2_Up(:,iGF_Beta_1) )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_Beta_2), &
                 G_K    (:,iGF_Beta_2), &
                 G_N_X2 (:,iGF_Beta_2), &
                 G_X2_Dn(:,iGF_Beta_2), &
                 G_X2_Up(:,iGF_Beta_2) )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_Beta_3), &
                 G_K    (:,iGF_Beta_3), &
                 G_N_X2 (:,iGF_Beta_3), &
                 G_X2_Dn(:,iGF_Beta_3), &
                 G_X2_Up(:,iGF_Beta_3) )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_1), &
                 G_X2_Dn(:,iGF_Beta_1), &
                 G_K    (:,iGF_Beta_1), &
                 db1dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_2), &
                 G_X2_Dn(:,iGF_Beta_2), &
                 G_K    (:,iGF_Beta_2), &
                 db2dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_3), &
                 G_X2_Dn(:,iGF_Beta_3), &
                 G_K    (:,iGF_Beta_3), &
                 db3dX2 )

        ! --- Lapse function derivative wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_Alpha), &
                 G_K    (:,iGF_Alpha), &
                 G_N_X2 (:,iGF_Alpha), &
                 G_X2_Dn(:,iGF_Alpha), &
                 G_X2_Up(:,iGF_Alpha) )

        G_X2_Dn(:,iGF_Alpha) = MAX( G_X2_Dn(:,iGF_Alpha), SqrtTiny )
        G_X2_Up(:,iGF_Alpha) = MAX( G_X2_Up(:,iGF_Alpha), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Alpha), &
                 G_X2_Dn(:,iGF_Alpha), &
                 G_K    (:,iGF_Alpha), &
                 dadX2 )

        ! --- Momentum density source term (S2) ---

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2) &
              + G_K(:,iGF_Alpha) &
                  * (   Stress(:,1) * dh1dX2 / G_K(:,iGF_h_1) &
                      + Stress(:,2) * dh2dX2 / G_K(:,iGF_h_2) &
                      + Stress(:,3) * dh3dX2 / G_K(:,iGF_h_3) ) &
              + uCF_K(:,iCF_S1) * db1dX2 &
              + uCF_K(:,iCF_S2) * db2dX2 &
              + uCF_K(:,iCF_S3) * db3dX2 &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx2

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Scale factor derivatives wrt X3 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_h_1), &
                 G_K    (:,iGF_h_1), &
                 G_N_X3 (:,iGF_h_1), &
                 G_X3_Dn(:,iGF_h_1), &
                 G_X3_Up(:,iGF_h_1) )

        G_X3_Dn(:,iGF_h_1) = MAX( G_X3_Dn(:,iGF_h_1), SqrtTiny )
        G_X3_Up(:,iGF_h_1) = MAX( G_X3_Up(:,iGF_h_1), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_h_2), &
                 G_K    (:,iGF_h_2), &
                 G_N_X3 (:,iGF_h_2), &
                 G_X3_Dn(:,iGF_h_2), &
                 G_X3_Up(:,iGF_h_2) )

        G_X3_Dn(:,iGF_h_2) = MAX( G_X3_Dn(:,iGF_h_2), SqrtTiny )
        G_X3_Up(:,iGF_h_2) = MAX( G_X3_Up(:,iGF_h_2), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_h_3), &
                 G_K    (:,iGF_h_3), &
                 G_N_X3 (:,iGF_h_3), &
                 G_X3_Dn(:,iGF_h_3), &
                 G_X3_Up(:,iGF_h_3) )

        G_X3_Dn(:,iGF_h_3) = MAX( G_X3_Dn(:,iGF_h_3), SqrtTiny )
        G_X3_Up(:,iGF_h_3) = MAX( G_X3_Up(:,iGF_h_3), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_1), &
                 G_X3_Dn(:,iGF_h_1), &
                 G_K    (:,iGF_h_1), &
                 dh1dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_2), &
                 G_X3_Dn(:,iGF_h_2), &
                 G_K    (:,iGF_h_2), &
                 dh2dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_3), &
                 G_X3_Dn(:,iGF_h_3), &
                 G_K    (:,iGF_h_3), &
                 dh3dX3 )

        ! --- Shift vector derivatives wrt X3 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_Beta_1), &
                 G_K    (:,iGF_Beta_1), &
                 G_N_X3 (:,iGF_Beta_1), &
                 G_X3_Dn(:,iGF_Beta_1), &
                 G_X3_Up(:,iGF_Beta_1) )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_Beta_2), &
                 G_K    (:,iGF_Beta_2), &
                 G_N_X3 (:,iGF_Beta_2), &
                 G_X3_Dn(:,iGF_Beta_2), &
                 G_X3_Up(:,iGF_Beta_2) )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_Beta_3), &
                 G_K    (:,iGF_Beta_3), &
                 G_N_X3 (:,iGF_Beta_3), &
                 G_X3_Dn(:,iGF_Beta_3), &
                 G_X3_Up(:,iGF_Beta_3) )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_1), &
                 G_X3_Dn(:,iGF_Beta_1), &
                 G_K    (:,iGF_Beta_1), &
                 db1dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_2), &
                 G_X3_Dn(:,iGF_Beta_2), &
                 G_K    (:,iGF_Beta_2), &
                 db2dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_3), &
                 G_X3_Dn(:,iGF_Beta_3), &
                 G_K    (:,iGF_Beta_3), &
                 db3dX3 )

        ! --- Lapse function derivative wrt X3 ---

        ! -- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_Alpha), &
                 G_K    (:,iGF_Alpha), &
                 G_N_X3 (:,iGF_Alpha), &
                 G_X3_Dn(:,iGF_Alpha), &
                 G_X3_Up(:,iGF_Alpha) )

        G_X3_Dn(:,iGF_Alpha) = MAX( G_X3_Dn(:,iGF_Alpha), SqrtTiny )
        G_X3_Up(:,iGF_Alpha) = MAX( G_X3_Up(:,iGF_Alpha), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Alpha), &
                 G_X3_Dn(:,iGF_Alpha), &
                 G_K    (:,iGF_Alpha), &
                 dadX3 )

        ! --- Momentum density source term (S3) ---

        dU(:,iX1,iX2,iX3,iCF_S3) &
          = dU(:,iX1,iX2,iX3,iCF_S3)                                &
              + G_K(:,iGF_Alpha)                                    &
                  * ( ( Stress(:,1) * dh1dX3 ) / G_K(:,iGF_h_1)     &
                      + ( Stress(:,2) * dh2dX3 ) / G_K(:,iGF_h_2)   &
                      + ( Stress(:,3) * dh3dX3 ) / G_K(:,iGF_h_3) ) &
              + uCF_K(:,iCF_S1) * db1dX3 &
                  + uCF_K(:,iCF_S2) * db2dX3 &
                  + uCF_K(:,iCF_S3) * db3dX3 &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx3

      END IF

      EnergyDensitySourceTerms(1,:,iX1,iX2,iX3) &
        = -(   uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11) * dadx1 &
             + uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22) * dadx2 &
             + uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33) * dadx3 )

      ! --- Compute divergence of grid volume ---

      ! --- Interpolate SqrtGm to faces (do this with scale factors instead?)---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_SqrtGm), &
               G_K    (:,iGF_SqrtGm), &
               G_N_X1 (:,iGF_SqrtGm), &
               G_X1_Dn(:,iGF_SqrtGm), &
               G_X1_Up(:,iGF_SqrtGm) )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Alpha) * G_X1_Up(:,iGF_SqrtGm) &
                 * G_X1_Up(:,iGF_Beta_1), &
               G_X1_Dn(:,iGF_Alpha) * G_X1_Dn(:,iGF_SqrtGm) &
                 * G_X1_Dn(:,iGF_Beta_1), &
               G_K    (:,iGF_Alpha) * G_K    (:,iGF_SqrtGm) &
                 * G_K    (:,iGF_Beta_1), &
               DivGridVolume, &
               Alpha_Option = One, Beta_Option = Zero )

      IF( nDimsX .GT. 1 )THEN

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2 (:,iGF_SqrtGm), &
                 G_K    (:,iGF_SqrtGm), &
                 G_N_X2 (:,iGF_SqrtGm), &
                 G_X2_Dn(:,iGF_SqrtGm), &
                 G_X2_Up(:,iGF_SqrtGm) )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Alpha) * G_X2_Up(:,iGF_SqrtGm) &
                   * G_X2_Up(:,iGF_Beta_2), &
                 G_X2_Dn(:,iGF_Alpha) * G_X2_Dn(:,iGF_SqrtGm) &
                   * G_X2_Dn(:,iGF_Beta_2), &
                 G_K    (:,iGF_Alpha) * G_K    (:,iGF_SqrtGm) &
                   * G_K    (:,iGF_Beta_2), &
                 DivGridVolume, &
                 Alpha_Option = One, Beta_Option = One )

      END IF

      IF( nDimsX .GT. 2 )THEN

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3 (:,iGF_SqrtGm), &
                 G_K    (:,iGF_SqrtGm), &
                 G_N_X3 (:,iGF_SqrtGm), &
                 G_X3_Dn(:,iGF_SqrtGm), &
                 G_X3_Up(:,iGF_SqrtGm) )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Alpha) * G_X3_Up(:,iGF_SqrtGm) &
                   * G_X3_Up(:,iGF_Beta_3), &
                 G_X3_Dn(:,iGF_Alpha) * G_X3_Dn(:,iGF_SqrtGm) &
                   * G_X3_Dn(:,iGF_Beta_3), &
                 G_K    (:,iGF_Alpha) * G_K    (:,iGF_SqrtGm) &
                   * G_K    (:,iGF_Beta_3), &
                 DivGridVolume, &
                 Alpha_Option = One, Beta_Option = One )

      END IF

      DivGridVolume &
        = One / ( G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) ) * DivGridVolume

      ! --- Compute energy increment ---

      ! --- Extrinsic curvature term ---

      EnergyDensitySourceTerms(2,:,iX1,iX2,iX3) &
        =   PressureTensor(:,1,1) * db1dX1 &
          + PressureTensor(:,1,2) * db2dX1 &
          + PressureTensor(:,1,3) * db3dX1 &
          + PressureTensor(:,2,1) * db1dX2 &
          + PressureTensor(:,2,2) * db2dX2 &
          + PressureTensor(:,2,3) * db3dX2 &
          + PressureTensor(:,3,1) * db1dX1 &
          + PressureTensor(:,3,2) * db2dX3 &
          + PressureTensor(:,3,3) * db3dX3

      ! --- Need to add Christoffel symbol term ---

      DO iNX = 1, nDOFX

        Christoffel3D_X1(iNX,1,1) &
          = One / G_K(iNX,iGF_h_1) * dh1dX1(iNX)
        Christoffel3D_X1(iNX,1,2) &
          = One / G_K(iNX,iGF_h_1) * dh1dX2(iNX)
        Christoffel3D_X1(iNX,1,3) &
          = One / G_K(iNX,iGF_h_1) * dh1dX3(iNX)
        Christoffel3D_X1(iNX,2,1) &
          = Christoffel3D_X1(iNX,1,2)
        Christoffel3D_X1(iNX,2,2) &
          = -G_K(iNX,iGF_h_2) / G_K(iNX,iGF_h_1)**2 * dh2dX1(iNX)
        Christoffel3D_X1(iNX,2,3) &
          = Zero
        Christoffel3D_X1(iNX,3,1) &
          = Christoffel3D_X1(iNX,1,3)
        Christoffel3D_X1(iNX,3,2) &
          = Christoffel3D_X1(iNX,2,3)
        Christoffel3D_X1(iNX,3,3) &
          = -G_K(iNX,iGF_h_3) / G_K(iNX,iGF_h_1)**2 * dh3dX1(iNX)

        Christoffel3D_X2(iNX,1,1) &
          = -G_K(iNX,iGF_h_1) / G_K(iNX,iGF_h_2)**2 * dh1dX2(iNX)
        Christoffel3D_X2(iNX,1,2) &
          = One / G_K(iNX,iGF_h_2) * dh2dX1(iNX)
        Christoffel3D_X2(iNX,1,3) &
          = Zero
        Christoffel3D_X2(iNX,2,1) &
          = Christoffel3D_X2(iNX,1,2)
        Christoffel3D_X2(iNX,2,2) &
          = One / G_K(iNX,iGF_h_2) * dh2dX2(iNX)
        Christoffel3D_X2(iNX,2,3) &
          = One / G_K(iNX,iGF_h_2) * dh2dX3(iNX)
        Christoffel3D_X2(iNX,3,1) &
          = Christoffel3D_X2(iNX,1,3)
        Christoffel3D_X2(iNX,3,2) &
          = Christoffel3D_X2(iNX,2,3)
        Christoffel3D_X2(iNX,3,3) &
          = -G_K(iNX,iGF_h_3) / G_K(iNX,iGF_h_2)**2 * dh3dX2(iNX)

        Christoffel3D_X3(iNX,1,1) &
          = -G_K(iNX,iGF_h_1) / G_K(iNX,iGF_h_3)**2 * dh1dX3(iNX)
        Christoffel3D_X3(iNX,1,2) &
          = Zero
        Christoffel3D_X3(iNX,1,3) &
          = One / G_K(iNX,iGF_h_3) * dh3dX1(iNX)
        Christoffel3D_X3(iNX,2,1) &
          = Christoffel3D_X3(iNX,1,2)
        Christoffel3D_X3(iNX,2,2) &
          = -G_K(iNX,iGF_h_2) / G_K(iNX,iGF_h_3)**2 * dh2dX3(iNX)
        Christoffel3D_X3(iNX,2,3) &
          = One / G_K(iNX,iGF_h_3) * dh3dX2(iNX)
        Christoffel3D_X3(iNX,3,1) &
          = Christoffel3D_X3(iNX,1,3)
        Christoffel3D_X3(iNX,3,2) &
          = Christoffel3D_X3(iNX,2,3)
        Christoffel3D_X3(iNX,3,3) &
          = One / G_K(iNX,iGF_h_3) * dh3dX3(iNX)

        Xij(iNX,1,1) &
          = G_K(iNX,iGF_Alpha)**( -2 ) &
              * ( G_K(iNX,iGF_Gm_dd_11) * db1dX1(iNX) &
                    + G_K(iNX,iGF_h_1) &
                        * (   G_K(iNX,iGF_Beta_1) * dh1dX1(iNX)   &
                            + G_K(iNX,iGF_Beta_2) * dh1dX2(iNX)   &
                            + G_K(iNX,iGF_Beta_3) * dh1dX3(iNX) ) &
                    - Third * G_K(iNX,iGF_Gm_dd_11) * DivGridVolume(iNX) )

        Xij(iNX,1,2) &
          = Half * G_K(iNX,iGF_Alpha)**( -2 ) &
              * (   G_K(iNX,iGF_Gm_dd_11) * db1dX2(iNX) &
                  + G_K(iNX,iGF_Gm_dd_22) * db2dX1(iNX) )

        Xij(iNX,1,3) &
          = Half * G_K(iNX,iGF_Alpha)**( -2 ) &
              * (   G_K(iNX,iGF_Gm_dd_11) * db1dX3(iNX) &
                  + G_K(iNX,iGF_Gm_dd_33) * db3dX1(iNX) )

        Xij(iNX,2,1) = Xij(iNX,1,2)

        Xij(iNX,2,2) &
          = G_K(iNX,iGF_Alpha)**( -2 ) &
              * ( G_K(iNX,iGF_Gm_dd_22) * db2dX2(iNX) &
                    + G_K(iNX,iGF_h_2) &
                        * (   G_K(iNX,iGF_Beta_1) * dh2dX1(iNX)   &
                            + G_K(iNX,iGF_Beta_2) * dh2dX2(iNX)   &
                            + G_K(iNX,iGF_Beta_3) * dh2dX3(iNX) ) &
                    - Third * G_K(iNX,iGF_Gm_dd_22) * DivGridVolume(iNX) )

        Xij(iNX,2,3) &
          = Half * G_K(iNX,iGF_Alpha)**( -2 ) &
              * (   G_K(iNX,iGF_Gm_dd_22) * db2dX3(iNX) &
                  + G_K(iNX,iGF_Gm_dd_33) * db3dX2(iNX) )

        Xij(iNX,3,1) = Xij(iNX,1,3)

        Xij(iNX,3,2) = Xij(iNX,2,3)

        Xij(iNX,3,3) &
          = G_K(iNX,iGF_Alpha)**( -2 ) &
              * ( G_K(iNX,iGF_Gm_dd_33) * db3dX3(iNX) &
                    + G_K(iNX,iGF_h_3) &
                        * (   G_K(iNX,iGF_Beta_1) * dh3dX1(iNX)   &
                            + G_K(iNX,iGF_Beta_2) * dh3dX2(iNX)   &
                            + G_K(iNX,iGF_Beta_3) * dh3dX3(iNX) ) &
                    - Third * G_K(iNX,iGF_Gm_dd_33) * DivGridVolume(iNX) )

      END DO

      DO iDim = 1, 3
      DO jDim = 1, 3

        DO iNX = 1, nDOFX

          Christoffel_X1(iNX,iDim,jDim) &
            = G_K(iNX,iGF_Beta_1) * Xij(iNX,iDim,jDim) &
                + Christoffel3D_X1(iNX,iDim,jDim)

          Christoffel_X2(iNX,iDim,jDim) &
            = G_K(iNX,iGF_Beta_2) * Xij(iNX,iDim,jDim) &
                + Christoffel3D_X2(iNX,iDim,jDim)

          Christoffel_X3(iNX,iDim,jDim) &
            = G_K(iNX,iGF_Beta_3) * Xij(iNX,iDim,jDim) &
                + Christoffel3D_X3(iNX,iDim,jDim)

        END DO

      END DO
      END DO

      EnergyDensitySourceTerms(3,:,iX1,iX2,iX3) = Zero

      DO iNX = 1, nDOFX

        PressureTensor(iNX,1,:) &
          = G_K(iNX,iGF_Gm_dd_11) * PressureTensor(iNX,1,:)
        PressureTensor(iNX,2,:) &
          = G_K(iNX,iGF_Gm_dd_22) * PressureTensor(iNX,2,:)
        PressureTensor(iNX,3,:) &
          = G_K(iNX,iGF_Gm_dd_33) * PressureTensor(iNX,3,:)

      END DO

      ! Get gradient of conformal factor on upper face

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1 (:,iGF_Psi), &
               G_K    (:,iGF_Psi), &
               G_N_X1 (:,iGF_Psi), &
               G_X1_Dn(:,iGF_Psi), &
               G_X1_Up(:,iGF_Psi) )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Psi), &
               G_X1_Dn(:,iGF_Psi), &
               G_K    (:,iGF_Psi), &
               GradPsi,            &
               Alpha_Option = One, Beta_Option = Zero )

      GradPsiF = Zero

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GradPsi, 1, Zero, GradPsiF(1:nDOFX_X1), 1 )

      DO iDim = 1, 3

        EnergyDensitySourceTerms(3,:,iX1,iX2,iX3) &
        = EnergyDensitySourceTerms(3,:,iX1,iX2,iX3) &
            + PressureTensor(:,iDim,1) &
                * (   Christoffel_X1(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X1(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X1(:,iDim,3) * G_K(:,iGF_Beta_3) ) &
            + PressureTensor(:,iDim,2) &
                * (   Christoffel_X2(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X2(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X2(:,iDim,3) * G_K(:,iGF_Beta_3) ) &
            + PressureTensor(:,iDim,3) &
                * (   Christoffel_X3(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X3(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X3(:,iDim,3) * G_K(:,iGF_Beta_3) )

      END DO

      PressureTensorTrace &
        =   uCF_K(:,iCF_S1) * uPF_K(:,iPF_V1) &
          + uCF_K(:,iCF_S2) * uPF_K(:,iPF_V2) &
          + uCF_K(:,iCF_S3) * uPF_K(:,iPF_V3) + Three * P_K

      EnergyDensitySourceTerms(4,:,iX1,iX2,iX3) &
        = -Third * PressureTensorTrace * DivGridVolume

      EnergyDensitySourceTerms(5,:,iX1,iX2,iX3) &
        = DivGridVolume

      EnergyDensitySourceTerms(6,:,iX1,iX2,iX3) &
        = Half / G_K(:,iGF_Alpha)**2 &
            * ( Four / Three &
                  * ( db1dX1 + Christoffel_X1(:,1,1) * G_K(:,iGF_Beta_1) )**2 &
                  + G_K(:,iGF_Gm_dd_11) * G_K(:,iGF_Beta_1)**2 &
                      * ( One / G_K(:,iGF_Gm_dd_22) * Christoffel_X1(:,2,1)**2 &
                        + One / G_K(:,iGF_Gm_dd_33) * Christoffel_X1(:,3,1)**2 &
                        ) )

      EnergyDensitySourceTerms(7,:,iX1,iX2,iX3) &
        = GradPsiF

      ! --- Add to increments ---

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            + EnergyDensitySourceTerms(1,:,iX1,iX2,iX3) &
            + EnergyDensitySourceTerms(2,:,iX1,iX2,iX3) &
            + EnergyDensitySourceTerms(3,:,iX1,iX2,iX3) &
            + EnergyDensitySourceTerms(4,:,iX1,iX2,iX3)

      DO iCF = 1, nCF

        dU(:,iX1,iX2,iX3,iCF) &
          = dU(:,iX1,iX2,iX3,iCF) &
              - U(:,iX1,iX2,iX3,iCF) * DivGridVolume

      END DO

    END DO
    END DO
    END DO

#ifndef USE_AMREX_TRUE

    IF( WriteSourceTerms )THEN

      CALL WriteSourceTermDiagnosticsHDF( Time, EnergyDensitySourceTerms )

    END IF

#endif

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic_CPU


  SUBROUTINE ComputeIncrement_Geometry_Relativistic_GPU &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER :: iX1, iX2, iX3, iNX, iCF, iGF, i, k, iDim
    INTEGER :: nK(3), nGF_K

    REAL(DP) :: P(nPF)
    REAL(DP) :: Pressure
    REAL(DP) :: PressureTensor(3,3, nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))
    REAL(DP) :: PressureTensorTrace(nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))
    INTEGER :: iErr                (nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))

    REAL(DP) :: G_K_X1 (nDOFX,   nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_Dn_X1(nDOFX_X1,nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))
    REAL(DP) :: G_Up_X1(nDOFX_X1,nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))
    REAL(DP) :: dGdX1  (nDOFX,   nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))

    REAL(DP) :: G_K_X2 (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_Dn_X2(nDOFX_X2,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))
    REAL(DP) :: G_Up_X2(nDOFX_X2,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))
    REAL(DP) :: dGdX2  (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))

    REAL(DP) :: G_K_X3 (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_Dn_X3(nDOFX_X3,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))
    REAL(DP) :: G_Up_X3(nDOFX_X3,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))
    REAL(DP) :: dGdX3  (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))

    REAL(DP) :: EnergyDensitySourceTerms(7,nDOFX,iX_B0(1):iX_E0(1), &
                                                 iX_B0(2):iX_E0(2), &
                                                 iX_B0(3):iX_E0(3))

    REAL(DP) :: Christoffel3D_X1 (3,3)
    REAL(DP) :: Christoffel3D_X2 (3,3)
    REAL(DP) :: Christoffel3D_X3 (3,3)
    REAL(DP) :: Christoffel_X1   (3,3)
    REAL(DP) :: Christoffel_X2   (3,3)
    REAL(DP) :: Christoffel_X3   (3,3)
    REAL(DP) :: PressureTensor_ud(3,3)
    REAL(DP) :: Xij              (3,3)
    REAL(DP) :: DivGridVolume

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nK = iX_E0 - iX_B0 + 1
    nGF_K = ( nGF + 1 ) * PRODUCT( nK )

    CALL TimersStart_Euler( Timer_Euler_DG_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, dX1, dX2, dX3 ) &
    !$OMP MAP( alloc: PressureTensor, PressureTensorTrace, &
    !$OMP             G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$OMP             G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$OMP             G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$OMP             EnergyDensitySourceTerms, iErr )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, dX1, dX2, dX3 ) &
    !$ACC CREATE(     PressureTensor, PressureTensorTrace, &
    !$ACC             G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$ACC             G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$ACC             G_K_X3, G_Dn_X3, G_Up_X3, dGdX3, &
    !$ACC             EnergyDensitySourceTerms, iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyIn )

    ! --- Permute data ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K_X1, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF
    DO iNX = 1, nDOFX

      G_K_X1(iNX,iGF,iX2,iX3,iX1) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K_X1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1, nDOFX

      G_K_X1(iNX,nGF+1,iX2,iX3,iX1) &
        =   G_K_X1(iNX,iGF_Alpha ,iX2,iX3,iX1) &
          * G_K_X1(iNX,iGF_SqrtGm,iX2,iX3,iX1) &
          * G_K_X1(iNX,iGF_Beta_1,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO

    IF( nDimsX .GT. 1 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X2, G )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF
      DO iNX = 1, nDOFX

        G_K_X2(iNX,iGF,iX1,iX3,iX2) = G(iNX,iX1,iX2,iX3,iGF)

      END DO
      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        G_K_X2(iNX,nGF+1,iX1,iX3,iX2) &
          =   G_K_X2(iNX,iGF_Alpha ,iX1,iX3,iX2) &
            * G_K_X2(iNX,iGF_SqrtGm,iX1,iX3,iX2) &
            * G_K_X2(iNX,iGF_Beta_2,iX1,iX3,iX2)

      END DO
      END DO
      END DO
      END DO

    END IF

    IF( nDimsX .GT. 2 )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X3, G )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF
      DO iNX = 1, nDOFX

        G_K_X3(iNX,iGF,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF)

      END DO
      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
      DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        G_K_X3(iNX,nGF+1,iX1,iX2,iX3) &
          =   G_K_X3(iNX,iGF_Alpha ,iX1,iX2,iX3) &
            * G_K_X3(iNX,iGF_SqrtGm,iX1,iX2,iX3) &
            * G_K_X3(iNX,iGF_Beta_3,iX1,iX2,iX3)

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

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_Interpolate )

    ! --- Compute derivatives (X1) ---

    CALL TimersStart_Euler( Timer_Euler_DG_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X1, G_Up_X1, WeightsX_X1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF + 1
    DO iNX = 1, nDOFX_X1

      G_Dn_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Dn_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

      G_Up_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Up_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, G_K_X1, WeightsX_q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF + 1
    DO iNX = 1, nDOFX

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dGdX1, WeightsX_q, dX1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1, nGF + 1
    DO iNX = 1, nDOFX

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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X2, G_Up_X2, WeightsX_X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX_X2

        G_Dn_X2(iNX,iGF,iX1,iX3,iX2) &
          = G_Dn_X2(iNX,iGF,iX1,iX3,iX2) * WeightsX_X2(iNX)

        G_Up_X2(iNX,iGF,iX1,iX3,iX2) &
          = G_Up_X2(iNX,iGF,iX1,iX3,iX2) * WeightsX_X2(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X2, WeightsX_q )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX2, WeightsX_q, dX2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

        dGdX2(iNX,iGF,iX1,iX3,iX2) &
          = dGdX2(iNX,iGF,iX1,iX3,iX2) / ( WeightsX_q(iNX) * dX2(iX2) )

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_Dn_X3, G_Up_X3, WeightsX_X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX_X3

        G_Dn_X3(iNX,iGF,iX1,iX2,iX3) &
          = G_Dn_X3(iNX,iGF,iX1,iX2,iX3) * WeightsX_X3(iNX)

        G_Up_X3(iNX,iGF,iX1,iX2,iX3) &
          = G_Up_X3(iNX,iGF,iX1,iX2,iX3) * WeightsX_X3(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, G_K_X3, WeightsX_q )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX3, WeightsX_q, dX3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

        dGdX3(iNX,iGF,iX1,iX2,iX3) &
          = dGdX3(iNX,iGF,iX1,iX2,iX3) / ( WeightsX_q(iNX) * dX3(iX3) )

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_Euler_DG_Permute )

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B0, iX_E0, dGdX3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iGF = 1, nGF + 1
      DO iNX = 1, nDOFX

        dGdX3(iNX,iGF,iX1,iX2,iX3) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

    ! --- Contributions from time-independent metric ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( P, Pressure )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( P, Pressure ) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, U, G, &
    !$ACC          PressureTensor, PressureTensorTrace, &
    !$ACC          dGdX1, dGdX2, dGdX3, EnergyDensitySourceTerms, iErr )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( P, Pressure )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

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

      PressureTensorTrace(iNX,iX1,iX2,iX3) &
        =   U(iNX,iX1,iX2,iX3,iCF_S1) * P(iPF_V1) &
          + U(iNX,iX1,iX2,iX3,iCF_S2) * P(iPF_V2) &
          + U(iNX,iX1,iX2,iX3,iCF_S3) * P(iPF_V3) &
          + Three * Pressure

      ! --- X1 increments ---

      ! --- Momentum increment ---

      dU(iNX,iX1,iX2,iX3,iCF_S1) &
        = dU(iNX,iX1,iX2,iX3,iCF_S1) &
            + G(iNX,iX1,iX2,iX3,iGF_Alpha) &
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
                * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      ! --- Energy increment ---

      EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
        = -U(iNX,iX1,iX2,iX3,iCF_S1) &
             / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
             * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      ! --- X2 increments ---

      IF( nDimsX .GT. 1 )THEN

        dU(iNX,iX1,iX2,iX3,iCF_S2) &
          = dU(iNX,iX1,iX2,iX3,iCF_S2) &
              + G(iNX,iX1,iX2,iX3,iGF_Alpha) &
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
                  * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2)

        EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
          = EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
              -U(iNX,iX1,iX2,iX3,iCF_S2) &
                 / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                 * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2)

      END IF

      ! --- X3 increments ---

      IF( nDimsX .GT. 2 )THEN

        dU(iNX,iX1,iX2,iX3,iCF_S3) &
          = dU(iNX,iX1,iX2,iX3,iCF_S3) &
              + G(iNX,iX1,iX2,iX3,iGF_Alpha) &
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
                  * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3)

          EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
            = EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
                -U(iNX,iX1,iX2,iX3,iCF_S3) &
                   / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                   * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3)

      END IF

      dU(iNX,iX1,iX2,iX3,iCF_E) &
        = dU(iNX,iX1,iX2,iX3,iCF_E) &
            + EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Contributions from time-dependent metric ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, &
    !$OMP          Christoffel_X1  , Christoffel_X2  , Christoffel_X3  , &
    !$OMP          Xij, PressureTensor_ud, DivGridVolume )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, &
    !$ACC          Christoffel_X1  , Christoffel_X2  , Christoffel_X3  , &
    !$ACC          Xij, PressureTensor_ud, DivGridVolume ) &
    !$ACC PRESENT( iX_B0, iX_E0, dU, U, G, dGdX1, dGdX2, dGdX3, &
    !$ACC          EnergyDensitySourceTerms, &
    !$ACC          PressureTensor, PressureTensorTrace )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, &
    !$OMP          Christoffel_X1  , Christoffel_X2  , Christoffel_X3  , &
    !$OMP          Xij, PressureTensor_ud, DivGridVolume )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      DivGridVolume &
        =  (   dGdX1(iNX,nGF+1,iX2,iX3,iX1) &
             + dGdX2(iNX,nGF+1,iX1,iX3,iX2) &
             + dGdX3(iNX,nGF+1,iX1,iX2,iX3) ) &
           / ( G(iNX,iX1,iX2,iX3,iGF_Alpha ) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) )

      ! --- Extrinsic curvature term ---

      EnergyDensitySourceTerms(2,iNX,iX1,iX2,iX3) &
        =   PressureTensor(1,1,iNX  ,iX1,iX2,iX3) &
              * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
          + PressureTensor(1,2,iNX  ,iX1,iX2,iX3) &
              * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) &
          + PressureTensor(1,3,iNX  ,iX1,iX2,iX3) &
              * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) &
          + PressureTensor(2,1,iNX  ,iX1,iX2,iX3) &
              * dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) &
          + PressureTensor(2,2,iNX  ,iX1,iX2,iX3) &
              * dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
          + PressureTensor(2,3,iNX  ,iX1,iX2,iX3) &
              * dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) &
          + PressureTensor(3,1,iNX  ,iX1,iX2,iX3) &
              * dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) &
          + PressureTensor(3,2,iNX  ,iX1,iX2,iX3) &
              * dGdX3(iNX,iGF_Beta_2,iX1,iX2,iX3) &
          + PressureTensor(3,3,iNX  ,iX1,iX2,iX3) &
              * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3)

        ! --- X1 ---

        Christoffel3D_X1(1,1) &
          = dGdX1(iNX,iGF_h_1,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(2,1) &
          = dGdX2(iNX,iGF_h_1,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(3,1) &
          = dGdX3(iNX,iGF_h_1,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(1,2) &
          = Christoffel3D_X1(2,1)

        Christoffel3D_X1(2,2) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_2) / G(iNX,iX1,iX2,iX3,iGF_h_1)**2 &
              * dGdX1(iNX,iGF_h_2,iX2,iX3,iX1)

        Christoffel3D_X1(3,2) &
          = Zero

        Christoffel3D_X1(1,3) &
          = Christoffel3D_X1(3,1)

        Christoffel3D_X1(2,3) &
          = Christoffel3D_X1(3,2)

        Christoffel3D_X1(3,3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_3) / G(iNX,iX1,iX2,iX3,iGF_h_1)**2 &
              * dGdX1(iNX,iGF_h_3,iX2,iX3,iX1)

        ! --- X2 ---

        Christoffel3D_X2(1,1) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_1) / G(iNX,iX1,iX2,iX3,iGF_h_2)**2 &
              * dGdX2(iNX,iGF_h_1,iX1,iX3,iX2)

        Christoffel3D_X2(2,1) &
          = dGdX1(iNX,iGF_h_2,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(3,1) &
          = Zero

        Christoffel3D_X2(1,2) &
          = Christoffel3D_X2(2,1)

        Christoffel3D_X2(2,2) &
          = dGdX2(iNX,iGF_h_2,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(3,2) &
          = dGdX3(iNX,iGF_h_2,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(1,3) &
          = Christoffel3D_X2(3,1)

        Christoffel3D_X2(2,3) &
          = Christoffel3D_X2(3,2)

        Christoffel3D_X2(3,3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_3) / G(iNX,iX1,iX2,iX3,iGF_h_2)**2 &
              * dGdX2(iNX,iGF_h_3,iX1,iX3,iX2)

        ! --- X3 ---

        Christoffel3D_X3(1,1) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_1) / G(iNX,iX1,iX2,iX3,iGF_h_3)**2 &
              * dGdX3(iNX,iGF_h_1,iX1,iX2,iX3)

        Christoffel3D_X3(2,1) &
          = Zero

        Christoffel3D_X3(3,1) &
          = dGdX1(iNX,iGF_h_3,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_3)

        Christoffel3D_X3(1,2) &
          = Christoffel3D_X3(2,1)

        Christoffel3D_X3(2,2) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_2) / G(iNX,iX1,iX2,iX3,iGF_h_3)**2 &
              * dGdX3(iNX,iGF_h_2,iX1,iX2,iX3)

        Christoffel3D_X3(3,2) &
          = dGdX2(iNX,iGF_h_3,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_3)

        Christoffel3D_X3(1,3) &
          = Christoffel3D_X3(3,1)

        Christoffel3D_X3(2,3) &
          = Christoffel3D_X3(3,2)

        Christoffel3D_X3(3,3) &
          = dGdX3(iNX,iGF_h_3,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_3)

        Xij(1,1) &
          = G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * dGdX1 (iNX,iGF_Beta_1,iX2,iX3,iX1) &
                    + G(iNX,iX1,iX2,iX3,iGF_h_1) &
                        * (   G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                                * dGdX1(iNX,iGF_h_1,iX2,iX3,iX1) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                                * dGdX2(iNX,iGF_h_1,iX1,iX3,iX2) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                                * dGdX3(iNX,iGF_h_1,iX1,iX2,iX3) ) &
                    - Third * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                        * DivGridVolume )

        Xij(2,1) &
          = Half * G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                      * dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) &
                  + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                      * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) )

        Xij(3,1) &
          = Half * G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                      * dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) &
                  + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                      * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) )

        Xij(1,2) = Xij(2,1)

        Xij(2,2) &
          = G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * dGdX2(iNX,iGF_Beta_2  ,iX1,iX3,iX2) &
                    + G(iNX,iX1,iX2,iX3,iGF_h_2) &
                        * (   G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                                * dGdX1(iNX,iGF_h_2,iX2,iX3,iX1) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                                * dGdX2(iNX,iGF_h_2,iX1,iX3,iX2) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                                * dGdX3(iNX,iGF_h_2,iX1,iX2,iX3) ) &
                    - Third * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                        * DivGridVolume )

        Xij(3,2) &
          = Half * G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                      * dGdX3(iNX,iGF_Beta_2,iX1,iX2,iX3) &
                  + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                      * dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) )

        Xij(1,3) = Xij(3,1)

        Xij(2,3) = Xij(3,2)

        Xij(3,3) &
          = G(iNX,iX1,iX2,iX3,iGF_Alpha)**( -2 ) &
              * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
                    + G(iNX,iX1,iX2,iX3,iGF_h_3) &
                        * (   G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                                * dGdX1(iNX,iGF_h_3,iX2,iX3,iX1) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                                * dGdX2(iNX,iGF_h_3,iX1,iX3,iX2) &
                            + G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                                * dGdX3(iNX,iGF_h_3,iX1,iX2,iX3) ) &
                    - Third * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                        * DivGridVolume )

      DO k = 1, 3
      DO i = 1, 3

        Christoffel_X1(i,k) &
          = G(iNX,iX1,iX2,iX3,iGF_Beta_1) * Xij(i,k) &
              + Christoffel3D_X1(i,k)

        Christoffel_X2(i,k) &
          = G(iNX,iX1,iX2,iX3,iGF_Beta_2) * Xij(i,k) &
              + Christoffel3D_X2(i,k)

        Christoffel_X3(i,k) &
          = G(iNX,iX1,iX2,iX3,iGF_Beta_3) * Xij(i,k) &
              + Christoffel3D_X3(i,k)

      END DO
      END DO

      EnergyDensitySourceTerms(3,iNX,iX1,iX2,iX3) = Zero

      DO iDim = 1, 3

        PressureTensor_ud(1,iDim) &
          = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
              * PressureTensor(1,iDim,iNX,iX1,iX2,iX3)

        PressureTensor_ud(2,iDim) &
          = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
              * PressureTensor(2,iDim,iNX,iX1,iX2,iX3)

        PressureTensor_ud(3,iDim) &
          = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
              * PressureTensor(3,iDim,iNX,iX1,iX2,iX3)

      END DO

      DO iDim = 1, 3

        EnergyDensitySourceTerms(3,iNX,iX1,iX2,iX3) &
          = EnergyDensitySourceTerms(3,iNX,iX1,iX2,iX3) &
              + PressureTensor_ud(iDim,1) &
                  * (   Christoffel_X1(iDim,1) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                      + Christoffel_X1(iDim,2) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                      + Christoffel_X1(iDim,3) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
              + PressureTensor_ud(iDim,2) &
                  * (   Christoffel_X2(iDim,1) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                      + Christoffel_X2(iDim,2) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                      + Christoffel_X2(iDim,3) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
              + PressureTensor_ud(iDim,3) &
                  * (   Christoffel_X3(iDim,1) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                      + Christoffel_X3(iDim,2) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                      + Christoffel_X3(iDim,3) &
                          * G(iNX,iX1,iX2,iX3,iGF_Beta_3) )

      END DO

      EnergyDensitySourceTerms(4,iNX,iX1,iX2,iX3) &
        = -Third * PressureTensorTrace(iNX,iX1,iX2,iX3) * DivGridVolume

      ! --- 5, 6, and 7 are diagnostics ---

      EnergyDensitySourceTerms(5,iNX,iX1,iX2,iX3) &
        = DivGridVolume

      EnergyDensitySourceTerms(6,iNX,iX1,iX2,iX3) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha)**2 &
            * ( Four / Three * ( dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                  + Christoffel_X1(1,1) &
                      * G(iNX,iX1,iX2,iX3,iGF_Beta_1) )**2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * G(iNX,iX1,iX2,iX3,iGF_Beta_1)**2 &
                    * ( One / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                          * Christoffel_X1(2,1)**2 &
                          + One / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                          * Christoffel_X1(3,1)**2 ) )

      ! --- 7 (GradPsi) is deferred ---

      ! --- Energy Increment ---

      dU(iNX,iX1,iX2,iX3,iCF_E) &
        = dU(iNX,iX1,iX2,iX3,iCF_E) &
            + EnergyDensitySourceTerms(2,iNX,iX1,iX2,iX3) &
            + EnergyDensitySourceTerms(3,iNX,iX1,iX2,iX3) &
            + EnergyDensitySourceTerms(4,iNX,iX1,iX2,iX3)

      DO iCF = 1, nCF

        dU(iNX,iX1,iX2,iX3,iCF) &
          = dU(iNX,iX1,iX2,iX3,iCF) &
              - U(iNX,iX1,iX2,iX3,iCF) * DivGridVolume

      END DO

    END DO
    END DO
    END DO
    END DO

#endif

    CALL TimersStart_Euler( Timer_Euler_DG_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    EnergyDensitySourceTerms, iErr ) &
    !$OMP MAP( release: iX_B0, iX_E0, dX1, dX2, dX3, &
    !$OMP               PressureTensor, PressureTensorTrace, &
    !$OMP               G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$OMP               G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$OMP               G_K_X3, G_Dn_X3, G_Up_X3, dGdX3 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      EnergyDensitySourceTerms, iErr ) &
    !$ACC DELETE(       iX_B0, iX_E0, dX1, dX2, dX3, &
    !$ACC               PressureTensor, PressureTensorTrace, &
    !$ACC               G_K_X1, G_Dn_X1, G_Up_X1, dGdX1, &
    !$ACC               G_K_X2, G_Dn_X2, G_Up_X2, dGdX2, &
    !$ACC               G_K_X3, G_Dn_X3, G_Up_X3, dGdX3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_DG_CopyOut )

    END ASSOCIATE ! dX1, dX2, dX3

    CALL TimersStart_Euler( Timer_Euler_DG_ErrorCheck )

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR: ComputeIncrement_Geometry_Relativistic'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

          CALL DescribeError_Euler( iErr(iNX,iX1,iX2,iX3) )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_DG_ErrorCheck )

#ifndef USE_AMREX_TRUE

    IF( WriteSourceTerms )THEN

      CALL WriteSourceTermDiagnosticsHDF( Time, EnergyDensitySourceTerms )

    END IF

#endif

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic_GPU


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


  SUBROUTINE InterpolateToFace &
    ( nDOFX_X, LX_Up, LX_Dn, &
      InterpolantP, InterpolantK, InterpolantN, &
      Answer_Dn, Answer_Up )

    INTEGER,  INTENT(in)  :: nDOFX_X
    REAL(DP), INTENT(in)  :: LX_Up(:,:)
    REAL(DP), INTENT(in)  :: LX_Dn(:,:)
    REAL(DP), INTENT(in)  :: InterpolantP(:), InterpolantK(:), InterpolantN(:)
    REAL(DP), INTENT(out) :: Answer_Dn(:), Answer_Up(:)

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, One,  LX_Up, nDOFX_X, &
             InterpolantP, 1, Zero, Answer_Dn, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, Half, LX_Dn, nDOFX_X, &
             InterpolantK, 1, Half, Answer_Dn, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, One,  LX_Up, nDOFX_X, &
             InterpolantK, 1, Zero, Answer_Up, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, Half, LX_Dn, nDOFX_X, &
             InterpolantN, 1, Half, Answer_Up, 1 )

  END SUBROUTINE InterpolateToFace


  SUBROUTINE ComputeDerivative &
    ( nDOFX_X, dX, LX_Up, LX_Dn, dLX, WeightsX_X, &
      UpperFaceValues, LowerFaceValues, VolumeValues, &
      Answer, Alpha_Option, Beta_Option )

    INTEGER,  INTENT(in)  :: nDOFX_X
    REAL(DP), INTENT(in)  :: dX
    REAL(DP), INTENT(in)  :: LX_Up(:,:)
    REAL(DP), INTENT(in)  :: LX_Dn(:,:)
    REAL(DP), INTENT(in)  :: dLX  (:,:)
    REAL(DP), INTENT(in)  :: WeightsX_X(:)
    REAL(DP), INTENT(in)  :: UpperFaceValues(:)
    REAL(DP), INTENT(in)  :: LowerFaceValues(:)
    REAL(DP), INTENT(in)  :: VolumeValues   (:)
    REAL(DP), INTENT(out) :: Answer(:)
    REAL(DP), INTENT(in), OPTIONAL :: Alpha_Option, Beta_Option

    REAL(DP) :: Alpha, Beta

    Alpha = One
    IF( PRESENT( Alpha_Option ) ) Alpha = Alpha_Option

    Beta = Zero
    IF( PRESENT( Beta_Option ) ) Beta = Beta_Option

    Answer = Zero

    CALL DGEMV( 'T', nDOFX_X, nDOFX, + Alpha, LX_Up, nDOFX_X, &
                WeightsX_X * UpperFaceValues, 1, Beta, Answer, 1 )

    CALL DGEMV( 'T', nDOFX_X, nDOFX, - Alpha, LX_Dn, nDOFX_X, &
                WeightsX_X * LowerFaceValues, 1, One , Answer, 1 )

    CALL DGEMV( 'T', nDOFX,    nDOFX, - Alpha, dLX, nDOFX,    &
                WeightsX_q  * VolumeValues  , 1, One , Answer, 1 )

    Answer = Answer / ( WeightsX_q * dX )

  END SUBROUTINE ComputeDerivative


END MODULE Euler_dgDiscretizationModule
