MODULE MHD_DiscretizationModule_Relativistic

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
    nNodesX, &
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
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    iGF_Phi_N, &
    CoordinateSystem
  USE MagnetofluidFieldsModule, ONLY: &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nCM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nPM, &
    iDM_Sh_X1, &
    iDM_Sh_X2, &
    iDM_Sh_X3, &
    iDM_Div, &
    nDM
  USE MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD
  USE MHD_UtilitiesModule, ONLY: &
    ComputePrimitive_MHD, &
    Eigenvalues_MHD, &
    Flux_X1_MHD, &
    Flux_X2_MHD, &
    Flux_X3_MHD, &
    NumericalFlux_MHD_X1, &
    NumericalFlux_MHD_X2, &
    NumericalFlux_MHD_X3
  USE MHD_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_MHD
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
   ComputeMagneticDivergence_MHD_Relativistic, &
   ComputeWeakMagneticDivergence_MHD_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_MHD_DG_Explicit

  LOGICAL  :: EvolveOnlyMagnetic
  LOGICAL  :: UseDivergenceCleaning
  LOGICAL  :: UsePowellSource
  REAL(DP) :: DampingParameter

  REAL(DP), POINTER, CONTIGUOUS :: &
    Gm_dd_11_K(:), Gm_dd_22_K(:), Gm_dd_33_K(:), SqrtGm_K(:), &
    Gm_dd_11_F(:), Gm_dd_22_F(:), Gm_dd_33_F(:), SqrtGm_F(:), &
    Beta_1_K  (:), Beta_2_K  (:), Beta_3_K  (:), Alpha_K (:), &
    Beta_1_F  (:), Beta_2_F  (:), Beta_3_F  (:), Alpha_F (:)

  REAL(DP), POINTER, CONTIGUOUS :: &
    uD_K(:), uS1_K(:), uS2_K(:), uS3_K(:), uE_K(:), uNe_K(:), &
    uB1_K(:), uB2_K(:), uB3_K(:), uChi_K(:), &
    uD_L(:), uS1_L(:), uS2_L(:), uS3_L(:), uE_L(:), uNe_L(:), &
    uB1_L(:), uB2_L(:), uB3_L(:), uChi_L(:), &
    uD_R(:), uS1_R(:), uS2_R(:), uS3_R(:), uE_R(:), uNe_R(:), &
    uB1_R(:), uB2_R(:), uB3_R(:), uChi_R(:)

  REAL(DP), ALLOCATABLE :: &
    pD_K(:), pV1_K(:), pV2_K(:), pV3_K(:), pE_K(:), pNe_K(:), &
    pB1_K(:), pB2_K(:), pB3_K(:), pChi_K(:), &
    pD_L(:), pV1_L(:), pV2_L(:), pV3_L(:), pE_L(:), pNe_L(:), &
    pB1_L(:), pB2_L(:), pB3_L(:), pChi_L(:), &
    pD_R(:), pV1_R(:), pV2_R(:), pV3_R(:), pE_R(:), pNe_R(:), &
    pB1_R(:), pB2_R(:), pB3_R(:), pChi_R(:)

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableX_F(:,:), IndexTableX_V(:,:)

  INTEGER :: nX   (3), nX_K,  nNodesX_K
  INTEGER :: nX_X1(3), nX1_X, nNodesX_X1
  INTEGER :: nX_X2(3), nX2_X, nNodesX_X2
  INTEGER :: nX_X3(3), nX3_X, nNodesX_X3

  REAL(DP), PUBLIC :: OffGridFlux_MHD_X1_Inner(nCM), &
                      OffGridFlux_MHD_X1_Outer(nCM), &
                      OffGridFlux_MHD_X2_Inner(nCM), &
                      OffGridFlux_MHD_X2_Outer(nCM), &
                      OffGridFlux_MHD_X3_Inner(nCM), &
                      OffGridFlux_MHD_X3_Outer(nCM)

CONTAINS


  SUBROUTINE ComputeIncrement_MHD_DG_Explicit &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
      SuppressBC_Option, &
      EvolveOnlyMagnetic_Option, &
      UseDivergenceCleaning_Option, &
      DampingParameter_Option, &
      UsePowellSource_Option, &
      SurfaceFlux_X1_Option, &
      SurfaceFlux_X2_Option, &
      SurfaceFlux_X3_Option )

    REAL(DP), INTENT(in)            :: t
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
      EvolveOnlyMagnetic_Option, &
      UseDivergenceCleaning_Option, &
      UsePowellSource_Option
    REAL(DP), INTENT(in),  OPTIONAL :: &
      DampingParameter_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      SurfaceFlux_X1_Option(:,:,:,:,:), &
      SurfaceFlux_X2_Option(:,:,:,:,:), &
      SurfaceFlux_X3_Option(:,:,:,:,:)

    ! -- Surface flux for coarse/fine corrections ---

    REAL(DP) :: &
      SurfaceFlux_X1(nDOFX_X1, &
          iX_B0(1):iX_E0(1)+1, &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3), &
          nCM)

    REAL(DP) :: &
      SurfaceFlux_X2(nDOFX_X2, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2)+1, &
          iX_B0(3):iX_E0(3), &
          nCM)

    REAL(DP) :: &
      SurfaceFlux_X3(nDOFX_X3, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3)+1, &
          nCM)

    INTEGER  :: iNX, iX1, iX2, iX3, iCM
    LOGICAL  :: SuppressBC
    REAL(DP) :: tau(nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3))

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    EvolveOnlyMagnetic = .FALSE.
    IF( PRESENT( EvolveOnlyMagnetic_Option ) ) &
      EvolveOnlyMagnetic = EvolveOnlyMagnetic_Option

    UseDivergenceCleaning = .FALSE.
    IF( PRESENT( UseDivergenceCleaning_Option ) ) &
      UseDivergenceCleaning = UseDivergenceCleaning_Option

    DampingParameter = 0.0_DP
    IF( UseDivergenceCleaning .AND. PRESENT( DampingParameter_Option ) )THEN
      DampingParameter = DampingParameter_Option
    END IF

    UsePowellSource = .FALSE.
    IF( PRESENT( UsePowellSource_Option ) ) &
      UsePowellSource = UsePowellSource_Option

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_MHD &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U )

      CALL DetectShocks_MHD &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, EvolveOnlyMagnetic, D )

    END IF

    CALL InitializeIncrement_MHD &
           ( iX_B0, iX_E0, iX_B1, iX_E1 )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      tau(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

    END DO
    END DO
    END DO
    END DO

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCM) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    OffGridFlux_MHD_X1_Inner = Zero
    OffGridFlux_MHD_X1_Outer = Zero
    OffGridFlux_MHD_X2_Inner = Zero
    OffGridFlux_MHD_X2_Outer = Zero
    OffGridFlux_MHD_X3_Inner = Zero
    OffGridFlux_MHD_X3_Outer = Zero

    CALL ComputeIncrement_MHD_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X1 )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

   !PRINT*, 'In node: ', iNX
   !PRINT*, 'In cell: ', iX1, iX2, iX3

   !PRINT*, 'Individual inc. after X1 divergence.'
   !PRINT*, 'CD: ', dU(iNX,iX1,iX2,iX3,iCM_D )
   !PRINT*, 'S1: ', dU(iNX,iX1,iX2,iX3,iCM_S1)
   !PRINT*, 'S2: ', dU(iNX,iX1,iX2,iX3,iCM_S2)
   !PRINT*, 'S3: ', dU(iNX,iX1,iX2,iX3,iCM_S3)
   !PRINT*, 'CE: ', dU(iNX,iX1,iX2,iX3,iCM_E )

    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_MHD_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X2 )

    CALL ComputeIncrement_MHD_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X3 )

    IF( PRESENT( SurfaceFlux_X1_Option ) )THEN
     SurfaceFlux_X1_Option = SurfaceFlux_X1
    END IF

    IF( PRESENT( SurfaceFlux_X2_Option ) )THEN
     SurfaceFlux_X2_Option = SurfaceFlux_X2
    END IF

    IF( PRESENT( SurfaceFlux_X3_Option ) )THEN
     SurfaceFlux_X3_Option = SurfaceFlux_X3
    END IF

    ! --- Multiply Inverse Mass Matrix ---

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

    !PRINT*, 'In cell: ', iX1, iX2, iX3
    !PRINT*, 'In node: ', iNX

    !PRINT*, 'tau: ', tau(iNX,iX1,iX2,iX3)
    !PRINT*, 'SqrtGm: ', G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

    !PRINT*, 'dU for field: ', iCM, ': ', dU(iNX,iX1,iX2,iX3,iCM)
    !PRINT*

      dU(iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            * tau(iNX,iX1,iX2,iX3) &
            / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * dX1(iX1) * dX2(iX2) * dX3(iX3) )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

   !PRINT*, 'In node: ', iNX
   !PRINT*, 'In cell: ', iX1, iX2, iX3

   !PRINT*, 'Individual inc. after matrix.'
   !PRINT*, 'CD: ', dU(iNX,iX1,iX2,iX3,iCM_D )
   !PRINT*, 'S1: ', dU(iNX,iX1,iX2,iX3,iCM_S1)
   !PRINT*, 'S2: ', dU(iNX,iX1,iX2,iX3,iCM_S2)
   !PRINT*, 'S3: ', dU(iNX,iX1,iX2,iX3,iCM_S3)
   !PRINT*, 'CE: ', dU(iNX,iX1,iX2,iX3,iCM_E )

    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

   !PRINT*, 'In node: ', iNX
   !PRINT*, 'In cell: ', iX1, iX2, iX3

   !PRINT*, 'Individual inc. after geometry.'
   !PRINT*, 'CD: ', dU(iNX,iX1,iX2,iX3,iCM_D )
   !PRINT*, 'S1: ', dU(iNX,iX1,iX2,iX3,iCM_S1)
   !PRINT*, 'S2: ', dU(iNX,iX1,iX2,iX3,iCM_S2)
   !PRINT*, 'S3: ', dU(iNX,iX1,iX2,iX3,iCM_S3)
   !PRINT*, 'CE: ', dU(iNX,iX1,iX2,iX3,iCM_E )

    END DO
    END DO
    END DO
    END DO

    IF( UsePowellSource ) THEN

      CALL ComputeIncrement_Powell &
            ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    END IF

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

   !PRINT*, 'In node: ', iNX
   !PRINT*, 'In cell: ', iX1, iX2, iX3

   !PRINT*, 'Individual inc. after gravity.'
   !PRINT*, 'CD: ', dU(iNX,iX1,iX2,iX3,iCM_D )
   !PRINT*, 'S1: ', dU(iNX,iX1,iX2,iX3,iCM_S1)
   !PRINT*, 'S2: ', dU(iNX,iX1,iX2,iX3,iCM_S2)
   !PRINT*, 'S3: ', dU(iNX,iX1,iX2,iX3,iCM_S3)
   !PRINT*, 'CE: ', dU(iNX,iX1,iX2,iX3,iCM_E )

    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_MHD

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_DG_Explicit


  SUBROUTINE ComputeIncrement_MHD_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X1 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X1(1:nDOFX_X1,iX_B0(1):iX_E0(1)+1, &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3),1:nCM)

    INTEGER  :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER  :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(2), EigVals_R(2)
    REAL(DP) :: Flux_L   (nCM), Flux_R   (nCM)
    REAL(DP) :: Flux_F   (nCM), Flux_K   (nCM)
    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

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
      uCM_K(nDOFX, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)-1:iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_L(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)
    REAL(DP) :: &
      uDM_R(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X1,nCM, &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3), &
                    iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCM, &
             iX_B0(2):iX_E0(2), &
             iX_B0(3):iX_E0(3), &
             iX_B0(1):iX_E0(1))

    ! --- X1 Increment ---

    REAL(DP) :: &
      dU_X1(nDOFX,nCM, &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1))

   !PRINT*, 'dU_X1: ', dU_X1
   !PRINT*, 'G_F:   ', G_F
   !PRINT*, 'uCM_L: ', uCM_L
   !PRINT*, 'uCM_R: ', uCM_R

   !PRINT*, 'dU_X1 after zeroing: ', dU_X1

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(2) ; iXP_E0(1) = iX_E0(2)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(1) ; iXP_E0(3) = iX_E0(1)

    ASSOCIATE( dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X1, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    G_F = Zero
    G_K = Zero
    uCM_K = Zero
    uCM_L = Zero
    uCM_R = Zero
    dU_X1 = Zero

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

    DO iCM = 1           , nCM
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX2,iX3,iX1,iCM) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      uDM_L(iX2,iX3,iX1,1) = D(1,iX1-1,iX2,iX3,iDM_Sh_X2)
      uDM_L(iX2,iX3,iX1,2) = D(1,iX1-1,iX2,iX3,iDM_Sh_X3)
      uDM_R(iX2,iX3,iX1,1) = D(1,iX1  ,iX2,iX3,iDM_Sh_X2)
      uDM_R(iX2,iX3,iX1,2) = D(1,iX1  ,iX2,iX3,iDM_Sh_X3)

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

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

    ! --- Shift vector (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX_X1 )

    ! --- Shift vector (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX_X1 )

    ! --- Compute metric and metric determinant on faces ---

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

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

    END DO

    ! --- Numerical Flux ---

    !PRINT*, 'Computing the primitive variables for the left state.'

    CALL ComputePrimitive_MHD &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             pB1_L, pB2_L, pB3_L, pChi_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    !PRINT*, 'Computing the primitive variables for the right state.'

    CALL ComputePrimitive_MHD &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             pB1_R, pB2_R, pB3_R, pChi_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    NumericalFlux  = Zero
    Flux_q         = Zero
    SurfaceFlux_X1 = Zero

    DO iNX_X = 1, nNodesX_X1

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_MHD &
            ( pV1_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_11_F(iNX_X), &
              Beta_1_F  (iNX_X), &
              pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_L &
        = Flux_X1_MHD &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .TRUE. )THEN

     !PRINT*
     !PRINT*, 'Left states and flux for iNX_X = ', iNX_X
     !PRINT*, 'D_L:      ', pD_L(iNX_X)
     !PRINT*, 'V1_L:     ', pV1_L(iNX_X)
     !PRINT*, 'V2_L:     ', pV2_L(iNX_X)
     !PRINT*, 'V3_L:     ', pV3_L(iNX_X)
     !PRINT*, 'E_L:      ', pE_L(iNX_X)
     !PRINT*, 'Ne_L:     ', pNe_L(iNX_X)
     !PRINT*, 'P_L:      ', P_L
     !PRINT*, 'Gm_dd_11_F: ', Gm_dd_11_F(iNX_X)
     !PRINT*, 'Gm_dd_22_F: ', Gm_dd_22_F(iNX_X)
     !PRINT*, 'Gm_dd_33_F: ', Gm_dd_33_F(iNX_X)
     !PRINT*, 'Alpha_F: ',    Alpha_F(iNX_X)
     !PRINT*, 'Flux_L:   ', Flux_L

      END IF

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_MHD &
            ( pV1_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_11_F(iNX_X), &
              Beta_1_F  (iNX_X), &
              pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_R &
        = Flux_X1_MHD &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .TRUE. )THEN

    !PRINT*
    !PRINT*, 'Right states and flux for iNX_X = ', iNX_X
    !PRINT*, 'D_R:      ', pD_R(iNX_X)
    !PRINT*, 'V1_R:     ', pV1_R(iNX_X)
    !PRINT*, 'V2_R:     ', pV2_R(iNX_X)
    !PRINT*, 'V3_R:     ', pV3_R(iNX_X)
    !PRINT*, 'E_R:      ', pE_R(iNX_X)
    !PRINT*, 'Ne_R:     ', pNe_R(iNX_X)
    !PRINT*, 'P_R:      ', P_R
    !PRINT*, 'Gm_dd_11_F: ', Gm_dd_11_F(iNX_X)
    !PRINT*, 'Gm_dd_22_F: ', Gm_dd_22_F(iNX_X)
    !PRINT*, 'Gm_dd_33_F: ', Gm_dd_33_F(iNX_X)
    !PRINT*, 'Alpha_F: ', Alpha_F(iNX_X)
    !PRINT*, 'Flux_R:   ', Flux_R

      END IF

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

   !PRINT*, 'AlphaMns: ', AlphaMns
   !PRINT*, 'AlphaPls: ', AlphaPls

      iNX = IndexTableX_F(1,iNX_X)
      iX2 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX1 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX2,iX3,iX1,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX2,iX3,iX1,iCM)

      END DO

      Flux_F &
        = NumericalFlux_MHD_X1 &
            ( uCM_L_nCM           , &
              uCM_R_nCM           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            , &
              [ Gm_dd_11_F(iNX_X), &
                Gm_dd_22_F(iNX_X), &
                Gm_dd_33_F(iNX_X), &
                Alpha_F(iNX_X),    &
                Beta_1_F(iNX_X),   &
                Beta_2_F(iNX_X),   &
                Beta_3_F(iNX_X)  ], &
              EvolveOnlyMagnetic  , &
              UseDivergenceCleaning )

     !PRINT*
     !PRINT*, 'uCM_L_nCM: ', uCM_L_nCM
     !PRINT*, 'uCM_R_nCM: ', uCM_R_nCM
     !PRINT*, 'Flux_L:    ', Flux_L
     !PRINT*, 'Flux_R:    ', Flux_R
     !PRINT*, 'Flux_F:    ', Flux_F
     !PRINT*

      DO iCM = 1, nCM

        SurfaceFlux_X1(iNX,iX1,iX2,iX3,iCM) &
          = Flux_F(iCM) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

       !PRINT*, 'SurfaceFlux_X1: ', SurfaceFlux_X1(iNX,iX1,iX2,iX3,iCM)

        NumericalFlux(iNX,iCM,iX2,iX3,iX1) &
          = Flux_F(iCM) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX2(iX2) * dX3(iX3) &
              * WeightsX_X1(iNX)

       !PRINT*, 'iNX, iX1, iX2, iX3, iCM: ', iNX, iX1, iX2, iX3, iCM
       !PRINT*, 'iNX_X: ', iNX_X
       !PRINT*, 'Alpha_F: ', Alpha_F(iNX_X)
       !PRINT*, 'SqrtGm_F: ', SqrtGm_F(iNX_X)
       !PRINT*, 'dX2(iX2): ', dX2(iX2)
       !PRINT*, 'dX3(iX3): ', dX3(iX3)
       !PRINT*, 'WeightsX_X1(iNX): ', WeightsX_X1(iNX)

      END DO ! iCM

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

   !PRINT*
   !PRINT*, 'iNX, iX2, iX3, iX1: ', iX2, iX3, iX1  
   !PRINT*, 'dU_X1(iCM_D ) before left: ', dU_X1(iNX,1,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S1) before left: ', dU_X1(iNX,2,iX2,iX3,iX1)
   !PRINT*
   !PRINT*, 'NumericalFlux Left:  ', NumericalFlux(iNX,:,iX2,iX3,iX1)
   !IF(iX1 + 1 .LE. iX_E1(1))!PRINT*, 'NumericalFlux Right: ', NumericalFlux(iNX,:,iX2,iX3,iX1+1)
   !PRINT*

    END DO
    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X1, + One, LX_X1_Dn, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, &
             Zero, dU_X1, nDOFX )

    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX
  
   !PRINT*
   !PRINT*, 'iNX, iX2, iX3, iX1: ', iX2, iX3, iX1  
   !PRINT*, 'dU_X1(iCM_D ) after left:  ', dU_X1(iNX,1,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S1) after left:  ', dU_X1(iNX,2,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S2) after left:  ', dU_X1(iNX,3,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S3) after left:  ', dU_X1(iNX,4,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_E ) after left:  ', dU_X1(iNX,5,iX2,iX3,iX1)
   !PRINT*

    END DO
    END DO
    END DO
    END DO

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X1, - One, LX_X1_Up, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, &
             One,  dU_X1, nDOFX )

    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

   !PRINT*
   !PRINT*, 'iNX, iX2, iX3, iX1: ', iX2, iX3, iX1
   !PRINT*, 'dU_X1(iCM_D ) after right: ', dU_X1(iNX,1,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S1) after right: ', dU_X1(iNX,2,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S2) after right: ', dU_X1(iNX,3,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_S3) after right: ', dU_X1(iNX,4,iX2,iX3,iX1)
   !PRINT*, 'dU_X1(iCM_E ) after right: ', dU_X1(iNX,5,iX2,iX3,iX1)
   !PRINT*

    END DO
    END DO
    END DO
    END DO

    !--------------------
    ! --- Volume Term ---
    !--------------------

   !PRINT*, 'Computing the primitive variables for the volume term.'

    CALL ComputePrimitive_MHD &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             pB1_K, pB2_K, pB3_K, pChi_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             Alpha_K, Beta_1_K, Beta_2_K, Beta_3_K, &
             EvolveOnlyMagnetic )

    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X1_MHD &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              pB1_K     (iNX_K), &
              pB2_K     (iNX_K), &
              pB3_K     (iNX_K), &
              pChi_K    (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_1_K  (iNX_K), &
              Beta_2_K  (iNX_K), &
              Beta_3_K  (iNX_K), &
              UseDivergenceCleaning )

      iNX = IndexTableX_V(1,iNX_K)
      iX2 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX1 = IndexTableX_V(4,iNX_K)

      DO iCM = 1, nCM

        Flux_q(iNX,iCM,iX2,iX3,iX1) &
          = Flux_K(iCM) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX2(iX2) * dX3(iX3) * WeightsX_q(iNX)

        !IF( iCM .EQ. 1 )PRINT*, 'Flux_q(iCM_D ): ', Flux_q(iNX,iCM,iX2,iX3,iX1)
        !IF( iCM .EQ. 2 )PRINT*, 'Flux_q(iCM_S1): ', Flux_q(iNX,iCM,iX2,iX3,iX1)
        !IF( iCM .EQ. 3 )PRINT*, 'Flux_q(iCM_S2): ', Flux_q(iNX,iCM,iX2,iX3,iX1)
        !IF( iCM .EQ. 4 )PRINT*, 'Flux_q(iCM_S3): ', Flux_q(iNX,iCM,iX2,iX3,iX1)
        !IF( iCM .EQ. 5 )PRINT*, 'Flux_q(iCM_E ): ', Flux_q(iNX,iCM,iX2,iX3,iX1)

       END DO

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX, One, dLXdX1_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X1, nDOFX )

    DO iCM = 1       , nCM
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            + dU_X1(iNX,iCM,iX2,iX3,iX1)

      !IF( iCM .EQ. 1 )PRINT*, 'dU(iCM_D ): ', dU(iNX,iX1,iX2,iX3,iCM)
      !IF( iCM .EQ. 2 )PRINT*, 'dU(iCM_S1): ', dU(iNX,iX1,iX2,iX3,iCM)
      !IF( iCM .EQ. 3 )PRINT*, 'dU(iCM_S2): ', dU(iNX,iX1,iX2,iX3,iCM)
      !IF( iCM .EQ. 4 )PRINT*, 'dU(iCM_S3): ', dU(iNX,iX1,iX2,iX3,iCM)
      !IF( iCM .EQ. 5 )PRINT*, 'dU(iCM_E ): ', dU(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iCM   = 1       , nCM
    DO iNX_X = 1       , nDOFX_X1

      OffGridFlux_MHD_X1_Inner(iCM) &
        = OffGridFlux_MHD_X1_Inner(iCM) &
            + NumericalFlux(iNX_X,iCM,iX2,iX3,iX_B0(1))

      OffGridFlux_MHD_X1_Outer(iCM) &
        = OffGridFlux_MHD_X1_Outer(iCM) &
            + NumericalFlux(iNX_X,iCM,iX2,iX3,iX_E0(1)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_Divergence_X1


  SUBROUTINE ComputeIncrement_MHD_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X2 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X2(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2)+1, &
                                iX_B0(3):iX_E0(3),1:nCM)

    INTEGER  :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER  :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(2), EigVals_R(2)
    REAL(DP) :: Flux_L   (nCM), Flux_R   (nCM)
    REAL(DP) :: Flux_F   (nCM), Flux_K   (nCM)
    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

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
      uCM_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)-1:iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X2, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X2, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_L(iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2)+1,2)
    REAL(DP) :: &
      uDM_R(iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X2,nCM, &
                    iX_B0(1):iX_E0(1), &
                    iX_B0(3):iX_E0(3), &
                    iX_B0(2):iX_E0(2)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCM, &
             iX_B0(1):iX_E0(1), &
             iX_B0(3):iX_E0(3), &
             iX_B0(2):iX_E0(2))

    ! --- X2 Increment ---

    REAL(DP) :: &
      dU_X2(nDOFX,nCM, &
            iX_B0(1):iX_E0(1), &
            iX_B0(3):iX_E0(3), &
            iX_B0(2):iX_E0(2))

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(2) ; iXP_E0(3) = iX_E0(2)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X2, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

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

    DO iCM = 1           , nCM
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX1,iX3,iX2,iCM) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDM_L(iX1,iX3,iX2,1) = D(1,iX1,iX2-1,iX3,iDM_Sh_X1)
      uDM_L(iX1,iX3,iX2,2) = D(1,iX1,iX2-1,iX3,iDM_Sh_X3)
      uDM_R(iX1,iX3,iX2,1) = D(1,iX1,iX2,  iX3,iDM_Sh_X1)
      uDM_R(iX1,iX3,iX2,2) = D(1,iX1,iX2,  iX3,iDM_Sh_X3)

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

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

    ! --- Shift vector (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Beta_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_1), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_1), nDOFX_X2 )

    ! --- Shift vector (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Beta_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_2), nDOFX_X2 )

    ! --- Shift vector (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_Beta_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_3), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_Beta_3), nDOFX_X2 )

    ! --- Compute metric and metric determinant on faces ---

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

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

    END DO

    ! --- Numerical Flux ---

    !PRINT*, 'Computing the primitive variables for the left state.'

    CALL ComputePrimitive_MHD &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             pB1_L, pB2_L, pB3_L, pChi_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    !PRINT*, 'Computing the primitive variables for the right state.'

    CALL ComputePrimitive_MHD &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             pB1_R, pB2_R, pB3_R, pChi_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    ! Initializations needed in debug mode
    NumericalFlux  = Zero
    Flux_q         = Zero
    SurfaceFlux_X2 = Zero

    DO iNX_X = 1, nNodesX_X2

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_MHD &
            ( pV2_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_22_F(iNX_X), &
              Beta_2_F  (iNX_X), &
              pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_L &
        = Flux_X2_MHD &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN

    !PRINT*
    !PRINT*, 'Left states and flux for iNX_X = ', iNX_X
    !PRINT*, 'D_L:      ', pD_L(iNX_X)
    !PRINT*, 'V1_L:     ', pV1_L(iNX_X)
    !PRINT*, 'V2_L:     ', pV2_L(iNX_X)
    !PRINT*, 'V3_L:     ', pV3_L(iNX_X)
    !PRINT*, 'E_L:      ', pE_L(iNX_X)
    !PRINT*, 'Ne_L:     ', pNe_L(iNX_X)
    !PRINT*, 'B1_L:     ', pB1_L(iNX_X)
    !PRINT*, 'B2_L:     ', pB2_L(iNX_X)
    !PRINT*, 'B3_L:     ', pB3_L(iNX_X)
    !PRINT*, 'Chi_L:    ', pChi_L(iNX_X)
    !PRINT*, 'P_L:      ', P_L
    !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
    !PRINT*, 'Flux_L:   ', Flux_L

      END IF

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_MHD &
            ( pV2_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_22_F(iNX_X), &
              Beta_2_F  (iNX_X), &
              pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_R &
        = Flux_X2_MHD &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN

     !PRINT*
     !PRINT*, 'Right states and flux for iNX_X = ', iNX_X
     !PRINT*, 'D_R:      ', pD_R(iNX_X)
     !PRINT*, 'V1_R:     ', pV1_R(iNX_X)
     !PRINT*, 'V2_R:     ', pV2_R(iNX_X)
     !PRINT*, 'V3_R:     ', pV3_R(iNX_X)
     !PRINT*, 'E_R:      ', pE_R(iNX_X)
     !PRINT*, 'Ne_R:     ', pNe_R(iNX_X)
     !PRINT*, 'B1_R:     ', pB1_R(iNX_X)
     !PRINT*, 'B2_R:     ', pB2_R(iNX_X)
     !PRINT*, 'B3_R:     ', pB3_R(iNX_X)
     !PRINT*, 'Chi_R:    ', pChi_R(iNX_X)
     !PRINT*, 'P_R:      ', P_R
     !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
     !PRINT*, 'Flux_R:   ', Flux_R

      END IF

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX2 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX1,iX3,iX2,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX1,iX3,iX2,iCM)

      END DO

      Flux_F &
        = NumericalFlux_MHD_X2 &
            ( uCM_L_nCM           , &
              uCM_R_nCM           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            )

      DO iCM = 1, nCM

        SurfaceFlux_X2(iNX,iX1,iX2,iX3,iCM) &
          = Flux_F(iCM) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

        NumericalFlux(iNX,iCM,iX1,iX3,iX2) &
          = Flux_F(iCM) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX1(iX1) * dX3(iX3) &
              * WeightsX_X2(iNX)

      END DO ! iCM

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X2, + One, LX_X2_Dn, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, &
             Zero, dU_X2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X2, - One, LX_X2_Up, nDOFX_X2, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, &
             One,  dU_X2, nDOFX )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    !PRINT*, 'Computing the primitive variables for the volume term.'

    CALL ComputePrimitive_MHD &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             pB1_K, pB2_K, pB3_K, pChi_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             Alpha_K, Beta_1_K, Beta_2_K, Beta_3_K, &
             EvolveOnlyMagnetic )

    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X2_MHD &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              pB1_K     (iNX_K), &
              pB2_K     (iNX_K), &
              pB3_K     (iNX_K), &
              pChi_K    (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_1_K  (iNX_K), &
              Beta_2_K  (iNX_K), &
              Beta_3_K  (iNX_K), &
              UseDivergenceCleaning )

      iNX = IndexTableX_V(1,iNX_K)
      iX1 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX2 = IndexTableX_V(4,iNX_K)

      DO iCM = 1, nCM

        Flux_q(iNX,iCM,iX1,iX3,iX2) &
          = Flux_K(iCM) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX1(iX1) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX, One, dLXdX2_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X2, nDOFX )

    DO iCM = 1       , nCM
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            + dU_X2(iNX,iCM,iX1,iX3,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iCM   = 1       , nCM
    DO iNX_X = 1       , nDOFX_X2

      OffGridFlux_MHD_X2_Inner(iCM) &
        = OffGridFlux_MHD_X2_Inner(iCM) &
            + NumericalFlux(iNX_X,iCM,iX1,iX3,iX_B0(2))

      OffGridFlux_MHD_X2_Outer(iCM) &
        = OffGridFlux_MHD_X2_Outer(iCM) &
            + NumericalFlux(iNX_X,iCM,iX1,iX3,iX_E0(2)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_Divergence_X2


  SUBROUTINE ComputeIncrement_MHD_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SurfaceFlux_X3 )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(out)   :: &
      SurfaceFlux_X3(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3)+1,1:nCM)

    INTEGER  :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER  :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(2), EigVals_R(2)
    REAL(DP) :: Flux_L   (nCM), Flux_R   (nCM)
    REAL(DP) :: Flux_F   (nCM), Flux_K   (nCM)
    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

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
      uCM_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)-1:iX_E0(3)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X3, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X3, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_L(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3)+1,2)
    REAL(DP) :: &
      uDM_R(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X3,nCM, &
                    iX_B0(1):iX_E0(1), &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCM, &
             iX_B0(1):iX_E0(1), &
             iX_B0(2):iX_E0(2), &
             iX_B0(3):iX_E0(3))

    ! --- X3 Increment ---

    REAL(DP) :: &
      dU_X3(nDOFX,nCM, &
            iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3))

    IF( iX_E0(3) .EQ. iX_B0(3) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(2) ; iXP_E0(2) = iX_E0(2)
    iXP_B0(3) = iX_B0(3) ; iXP_E0(3) = iX_E0(3)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X3, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

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

    DO iCM = 1           , nCM
    DO iX3 = iX_B0(3) - 1, iX_E0(3) + 1
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX1,iX2,iX3,iCM) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      uDM_L(iX1,iX2,iX3,1) = D(1,iX1,iX2,iX3-1,iDM_Sh_X1)
      uDM_L(iX1,iX2,iX3,2) = D(1,iX1,iX2,iX3-1,iDM_Sh_X2)
      uDM_R(iX1,iX2,iX3,1) = D(1,iX1,iX2,iX3,  iDM_Sh_X1)
      uDM_R(iX1,iX2,iX3,2) = D(1,iX1,iX2,iX3,  iDM_Sh_X2)

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

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

    ! --- Shift vector (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_Beta_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_1), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_1), nDOFX_X3 )

    ! --- Shift vector (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iGF_Beta_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_2), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             G_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iGF_Beta_2), nDOFX_X3 )

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

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               uCM_K(1,iX_B0(1),iX_B0(2),iX_B0(3)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCM), nDOFX_X3 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nX3_X, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               uCM_K(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(1),iX_B0(2),iX_B0(3)  ,iCM), nDOFX_X3 )

    END DO

    ! --- Numerical Flux ---

    !PRINT*, 'Computing the primitive variables for the left state.'

    CALL ComputePrimitive_MHD &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             pB1_L, pB2_L, pB3_L, pChi_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    !PRINT*, 'Computing the primitive variables for the right state.'

    CALL ComputePrimitive_MHD &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             pB1_R, pB2_R, pB3_R, pChi_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F, &
             EvolveOnlyMagnetic )

    ! Initializations needed in debug mode
    NumericalFlux  = Zero
    Flux_q         = Zero
    SurfaceFlux_X3 = Zero

    DO iNX_X = 1, nNodesX_X3

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_MHD &
            ( pV3_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_33_F(iNX_X), &
              Beta_3_F  (iNX_X), &
              pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_L &
        = Flux_X3_MHD &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN

    !PRINT*
    !PRINT*, 'Left states and flux for iNX_X = ', iNX_X
    !PRINT*, 'D_L:      ', pD_L(iNX_X)
    !PRINT*, 'V1_L:     ', pV1_L(iNX_X)
    !PRINT*, 'V2_L:     ', pV2_L(iNX_X)
    !PRINT*, 'V3_L:     ', pV3_L(iNX_X)
    !PRINT*, 'E_L:      ', pE_L(iNX_X)
    !PRINT*, 'Ne_L:     ', pNe_L(iNX_X)
    !PRINT*, 'B1_L:     ', pB1_L(iNX_X)
    !PRINT*, 'B2_L:     ', pB2_L(iNX_X)
    !PRINT*, 'B3_L:     ', pB3_L(iNX_X)
    !PRINT*, 'Chi_L:    ', pChi_L(iNX_X)
    !PRINT*, 'P_L:      ', P_L
    !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
    !PRINT*, 'Flux_L:   ', Flux_L

      END IF

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_MHD &
            ( pV3_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_33_F(iNX_X), &
              Beta_3_F  (iNX_X), &
              pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_R &
        = Flux_X3_MHD &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN

     !PRINT*
     !PRINT*, 'Right states and flux for iNX_X = ', iNX_X
     !PRINT*, 'D_R:      ', pD_R(iNX_X)
     !PRINT*, 'V1_R:     ', pV1_R(iNX_X)
     !PRINT*, 'V2_R:     ', pV2_R(iNX_X)
     !PRINT*, 'V3_R:     ', pV3_R(iNX_X)
     !PRINT*, 'E_R:      ', pE_R(iNX_X)
     !PRINT*, 'Ne_R:     ', pNe_R(iNX_X)
     !PRINT*, 'B1_R:     ', pB1_R(iNX_X)
     !PRINT*, 'B2_R:     ', pB2_R(iNX_X)
     !PRINT*, 'B3_R:     ', pB3_R(iNX_X)
     !PRINT*, 'Chi_R:    ', pChi_R(iNX_X)
     !PRINT*, 'P_R:      ', P_R
     !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
     !PRINT*, 'Flux_R:   ', Flux_R

      END IF

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX2 = IndexTableX_F(3,iNX_X)
      iX3 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX1,iX2,iX3,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX1,iX2,iX3,iCM)

      END DO

      Flux_F &
        = NumericalFlux_MHD_X3 &
            ( uCM_L_nCM           , &
              uCM_R_nCM           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            )

      DO iCM = 1, nCM

        SurfaceFlux_X3(iNX,iX1,iX2,iX3,iCM) &
          = Flux_F(iCM) * Alpha_F(iNX_X) * SqrtGm_F(iNX_X)

        NumericalFlux(iNX,iCM,iX1,iX2,iX3) &
          = Flux_F(iCM) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX1(iX1) * dX2(iX2) &
              * WeightsX_X3(iNX)

      END DO ! iCM

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X3, + One, LX_X3_Dn, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, &
             Zero, dU_X3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X3, - One, LX_X3_Up, nDOFX_X3, &
             NumericalFlux(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, &
             One,  dU_X3, nDOFX )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    !PRINT*, 'Computing the primitive variables for the volume term.'

    CALL ComputePrimitive_MHD &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             pB1_K, pB2_K, pB3_K, pChi_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             Alpha_K, Beta_1_K, Beta_2_K, Beta_3_K, &
             EvolveOnlyMagnetic )

    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X3_MHD &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              pB1_K     (iNX_K), &
              pB2_K     (iNX_K), &
              pB3_K     (iNX_K), &
              pChi_K    (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_1_K  (iNX_K), &
              Beta_2_K  (iNX_K), &
              Beta_3_K  (iNX_K), &
              UseDivergenceCleaning )

      iNX = IndexTableX_V(1,iNX_K)
      iX1 = IndexTableX_V(2,iNX_K)
      iX2 = IndexTableX_V(3,iNX_K)
      iX3 = IndexTableX_V(4,iNX_K)

      DO iCM = 1, nCM

        Flux_q(iNX,iCM,iX1,iX2,iX3) &
          = Flux_K(iCM) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX1(iX1) * dX2(iX2) * WeightsX_q(iNX)

      END DO

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX, One, dLXdX3_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X3, nDOFX )

    DO iCM = 1       , nCM
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            + dU_X3(iNX,iCM,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

    ! --- Off-Grid Fluxes for Conservation Tally ---

    DO iX2   = iX_B0(2), iX_E0(2)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iCM   = 1       , nCM
    DO iNX_X = 1       , nDOFX_X3

      OffGridFlux_MHD_X3_Inner(iCM) &
        = OffGridFlux_MHD_X3_Inner(iCM) &
            + NumericalFlux(iNX_X,iCM,iX1,iX2,iX_B0(3))

      OffGridFlux_MHD_X3_Outer(iCM) &
        = OffGridFlux_MHD_X3_Outer(iCM) &
            + NumericalFlux(iNX_X,iCM,iX1,iX2,iX_E0(3)+1)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G  (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      tau(:,iX_B1(1):,iX_B1(2):,iX_B1(3):)
    REAL(DP), INTENT(inout) :: &
      U  (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      dU (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(inout)        :: &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(in)           :: &
      tau(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(inout)        :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)

    INTEGER :: iNX, iX1, iX2, iX3

    REAL(DP) :: P(nPM)
    REAL(DP) :: Pressure
    REAL(DP) :: W, B0u, VdotB, bSq
    REAL(DP) :: PressureTensor(3,3, nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))

    REAL(DP) :: dGdX1  (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))
    REAL(DP) :: dGdX2  (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))
    REAL(DP) :: dGdX3  (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))

    CALL ComputeDerivatives_Geometry_Relativistic_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX1 )

    CALL ComputeDerivatives_Geometry_Relativistic_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX2 )

    CALL ComputeDerivatives_Geometry_Relativistic_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX3 )

   !PRINT*, 'iX_B0: ', iX_B0
   !PRINT*, 'iX_E0: ', iX_E0
   !PRINT*, 'iX_B1: ', iX_B1
   !PRINT*, 'iX_E1: ', iX_E1

    ! --- Contributions from time-independent metric ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      !PRINT*, 'In cell: ', iX1, iX2, iX3
      !PRINT*, 'In node: ', iNX

      !PRINT*, 'Computing the primitive variables for the time-independent metric source term.'

      CALL ComputePrimitive_MHD  &
             ( U(   iNX,iX1,iX2,iX3,iCM_D ),  &
               U(   iNX,iX1,iX2,iX3,iCM_S1),  &
               U(   iNX,iX1,iX2,iX3,iCM_S2),  &
               U(   iNX,iX1,iX2,iX3,iCM_S3),  &
               U(   iNX,iX1,iX2,iX3,iCM_E ),  &
               U(   iNX,iX1,iX2,iX3,iCM_Ne),  &
               U(   iNX,iX1,iX2,iX3,iCM_B1),  &
               U(   iNX,iX1,iX2,iX3,iCM_B2),  &
               U(   iNX,iX1,iX2,iX3,iCM_B3),  &
               U(   iNX,iX1,iX2,iX3,iCM_Chi), &
               P(iPM_D ),  &
               P(iPM_V1),  &
               P(iPM_V2),  &
               P(iPM_V3),  &
               P(iPM_E ),  &
               P(iPM_Ne),  &
               P(iPM_B1),  &
               P(iPM_B2),  &
               P(iPM_B3),  &
               P(iPM_Chi), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(   iNX,iX1,iX2,iX3,iGF_Alpha   ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_1  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_2  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_3  ), &
               EvolveOnlyMagnetic )

      CALL ComputePressureFromPrimitive &
             ( P(iPM_D), P(iPM_E), P(iPM_Ne), Pressure )

      W   = One &
            / SQRT( One &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * P(iPM_V1)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * P(iPM_V2)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * P(iPM_V3)**2 )

      B0u = ( W / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
            * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * P(iPM_V1) * U(iNX,iX1,iX2,iX3,iCM_B1) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  * P(iPM_V2) * U(iNX,iX1,iX2,iX3,iCM_B2) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  * P(iPM_V3) * U(iNX,iX1,iX2,iX3,iCM_B3) )

      VdotB = B0u * ( G(iNX,iX1,iX2,iX3,iGF_Alpha) / W )

      bSq = ( 1.0_DP / W**2 ) &
            * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * U(iNX,iX1,iX2,iX3,iCM_B1)**2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * U(iNX,iX1,iX2,iX3,iCM_B2)**2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * U(iNX,iX1,iX2,iX3,iCM_B3)**2 ) &
            + ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * P(iPM_V1) &
                * U(iNX,iX1,iX2,iX3,iCM_B1) &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                * P(iPM_V2) &
                * U(iNX,iX1,iX2,iX3,iCM_B2) &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                * P(iPM_V3) &
                * U(iNX,iX1,iX2,iX3,iCM_B3) )**2

      Pressure = Pressure + bSq / 2.0_DP

      PressureTensor(1,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S1) * P(iPM_V1) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
          - ( P(iPM_B1)**2 &
              + Two * B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_1)**2 )

      PressureTensor(2,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S2) * P(iPM_V1) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
          - ( P(iPM_B2) * P(iPM_B1) &
              + B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2) )

      PressureTensor(3,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V1) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3) * P(iPM_B1) &
              + B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_1) )

      PressureTensor(1,2,iNX,iX1,iX2,iX3) &
        = PressureTensor(2,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,2,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S2) * P(iPM_V2) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
          - ( P(iPM_B2)**2 &
              + Two * B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2)**2 )

      PressureTensor(3,2,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V2) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3) * P(iPM_B2) &
              + B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2) )

      PressureTensor(1,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,2,iNX,iX1,iX2,iX3)

      PressureTensor(3,3,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V3) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3)**2 &
              + Two * B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_3)**2 )

     !PRINT*, 'PressureTensor: ', PressureTensor(:,:,iNX,iX1,iX2,iX3)

      ! --- X1 increments ---

      ! --- Momentum increment ---

      dU(iNX,iX1,iX2,iX3,iCM_S1) &
        = dU(iNX,iX1,iX2,iX3,iCM_S1) &
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
                  + U(iNX,iX1,iX2,iX3,iCM_S1) &
                      * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCM_S2) &
                      * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCM_S3) &
                      * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) &
                  - ( U(iNX,iX1,iX2,iX3,iCM_D) + U(iNX,iX1,iX2,iX3,iCM_E) ) &
                      * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1) )

     !PRINT*, 'Momentum increment: ', dU(iNX,iX1,iX2,iX3,iCM_S1)

      ! --- Energy increment ---

      dU(iNX,iX1,iX2,iX3,iCM_E) &
        = dU(iNX,iX1,iX2,iX3,iCM_E) &
            - tau(iNX,iX1,iX2,iX3) &
                * U(iNX,iX1,iX2,iX3,iCM_S1) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      IF( UseDivergenceCleaning )THEN

       !PRINT*, 'Using divergence cleaning for source terms.'

        ! --- Eulerian Magnetic Field increment ---

        dU(iNX,iX1,iX2,iX3,iCM_B1) &
          = dU(iNX,iX1,iX2,iX3,iCM_B1) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) ) &
              + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                  * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * ( dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1) &
                        / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11 ) ) &
              + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                  * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                  * ( - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                            / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)**2 ) &
                        * dGdX1(iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
                      + ( dGdX1(iNX,iGF_SqrtGm,iX2,iX3,iX1) &
                            / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) ) )

        dU(iNX,iX1,iX2,iX3,iCM_B2) &
          = dU(iNX,iX1,iX2,iX3,iCM_B2) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) )

        dU(iNX,iX1,iX2,iX3,iCM_B3) &
          = dU(iNX,iX1,iX2,iX3,iCM_B3) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) )

        ! --- Divergence violation field increment ---

        dU(iNX,iX1,iX2,iX3,iCM_Chi) &
          = dU(iNX,iX1,iX2,iX3,iCM_Chi) &
              - DampingParameter * U(iNX,iX1,iX2,iX3,iCM_Chi) &
              + ( U(iNX,iX1,iX2,iX3,iCM_B1) / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
                * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      END IF

      ! -- X2 increments ---

      IF( nDimsX .GT. 1 )THEN

        dU(iNX,iX1,iX2,iX3,iCM_S2) &
          = dU(iNX,iX1,iX2,iX3,iCM_S2) &
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
                    + U(iNX,iX1,iX2,iX3,iCM_S1) &
                        * dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) &
                    + U(iNX,iX1,iX2,iX3,iCM_S2) &
                        * dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
                    + U(iNX,iX1,iX2,iX3,iCM_S3) &
                        * dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) &
                    - ( U(iNX,iX1,iX2,iX3,iCM_D) + U(iNX,iX1,iX2,iX3,iCM_E) ) &
                        * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2) )

        dU(iNX,iX1,iX2,iX3,iCM_E) &
          = dU(iNX,iX1,iX2,iX3,iCM_E) &
              -tau(iNX,iX1,iX2,iX3) &
                 * U(iNX,iX1,iX2,iX3,iCM_S2) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                 * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2)

        IF( UseDivergenceCleaning )THEN

          !PRINT*, 'Using divergence cleaning for source terms.'

          ! --- Eulerian Magnetic Field increment ---

          dU(iNX,iX1,iX2,iX3,iCM_B2) &
            = dU(iNX,iX1,iX2,iX3,iCM_B2) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B2) &
                      * dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) ) &
                + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * ( dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2) &
                          / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22 ) ) &
                + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                    * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                    * ( - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)**2 ) &
                          * dGdX2(iNX,iGF_Gm_dd_22,iX1,iX3,iX2) &
                        + ( dGdX2(iNX,iGF_SqrtGm,iX1,iX3,iX2) &
                              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) ) )

          dU(iNX,iX1,iX2,iX3,iCM_B1) &
            = dU(iNX,iX1,iX2,iX3,iCM_B1) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B2) &
                      * dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) )

          dU(iNX,iX1,iX2,iX3,iCM_B3) &
            = dU(iNX,iX1,iX2,iX3,iCM_B3) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B2) &
                      * dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) )

          ! --- Divergence violation field increment ---

          dU(iNX,iX1,iX2,iX3,iCM_Chi) &
            = dU(iNX,iX1,iX2,iX3,iCM_Chi) &
                + ( U(iNX,iX1,iX2,iX3,iCM_B2) / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
                  * dGdX2(iNX,iGF_Alpha,iX1,iX3,iX2)

        END IF

      END IF

      ! -- X3 increments ---

      IF( nDimsX .GT. 2 )THEN

        dU(iNX,iX1,iX2,iX3,iCM_S3) &
          = dU(iNX,iX1,iX2,iX3,iCM_S3) &
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
                    + U(iNX,iX1,iX2,iX3,iCM_S1) &
                        * dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) &
                    + U(iNX,iX1,iX2,iX3,iCM_S2) &
                        * dGdX3(iNX,iGF_Beta_2,iX1,iX2,iX3) &
                    + U(iNX,iX1,iX2,iX3,iCM_S3) &
                        * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
                    - ( U(iNX,iX1,iX2,iX3,iCM_D) + U(iNX,iX1,iX2,iX3,iCM_E) ) &
                        * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3) )

        dU(iNX,iX1,iX2,iX3,iCM_E) &
          = dU(iNX,iX1,iX2,iX3,iCM_E) &
              -tau(iNX,iX1,iX2,iX3) &
                 * U(iNX,iX1,iX2,iX3,iCM_S3) / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                 * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3)

        IF( UseDivergenceCleaning )THEN

          !PRINT*, 'Using divergence cleaning for source terms.'

          ! --- Eulerian Magnetic Field increment ---

          dU(iNX,iX1,iX2,iX3,iCM_B3) &
            = dU(iNX,iX1,iX2,iX3,iCM_B3) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B3) &
                      * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) ) &
                + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * ( dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3) &
                          / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33 ) ) &
                + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                    * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                    * ( - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)**2 ) &
                          * dGdX3(iNX,iGF_Gm_dd_33,iX1,iX2,iX3) &
                        + ( dGdX3(iNX,iGF_SqrtGm,iX1,iX2,iX3) &
                              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ) )

          dU(iNX,iX1,iX2,iX3,iCM_B1) &
            = dU(iNX,iX1,iX2,iX3,iCM_B1) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B3) &
                      * dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) )

          dU(iNX,iX1,iX2,iX3,iCM_B2) &
            = dU(iNX,iX1,iX2,iX3,iCM_B2) &
                - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(iNX,iX1,iX2,iX3,iCM_B3) &
                      * dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) )

          ! --- Divergence violation field increment ---

          dU(iNX,iX1,iX2,iX3,iCM_Chi) &
            = dU(iNX,iX1,iX2,iX3,iCM_Chi) &
                + ( U(iNX,iX1,iX2,iX3,iCM_B3) / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
                  * dGdX3(iNX,iGF_Alpha,iX1,iX2,iX3)

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO

    ! --- Contributions from time-dependent metric ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU(iNX,iX1,iX2,iX3,iCM_E) &
        = dU(iNX,iX1,iX2,iX3,iCM_E) &
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

      IF( UseDivergenceCleaning)THEN

        dU(iNX,iX1,iX2,iX3,iCM_Chi) &
          = dU(iNX,iX1,iX2,iX3,iCM_Chi) &
              - tau(iNX,iX1,iX2,iX3) * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                 * (   G(iNX,iX1,iX2,iX3,iGF_K_dd_11) &
                       / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                     + G(iNX,iX1,iX2,iX3,iGF_K_dd_22) &
                       / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                     + G(iNX,iX1,iX2,iX3,iGF_K_dd_33) &
                       / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ) &
                 * U(iNX,iX1,iX2,iX3,iCM_Chi)
      END IF

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic


  SUBROUTINE ComputeIncrement_Powell &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

   REAL(DP), INTENT(in)    :: t
   INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    CALL ComputeIncrement_Powell_Relativistic &
           ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

  END SUBROUTINE ComputeIncrement_Powell


  SUBROUTINE ComputeIncrement_Powell_Relativistic &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    REAL(DP) :: WeakDiv(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER :: iX1, iX2, iX3, iNX, iCM, iGF

    REAL(DP) :: P(nPM)
    REAL(DP) :: W, B0u, VdotB, bSq

    ! --- Compute Magnetic Divergence for Powell Sources ---

    CALL ComputeMagneticDivergence_MHD_Relativistic    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )
    !CALL ComputeWeakMagneticDivergence_MHD_Relativistic( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

    DO iNX = 1       , nDOFX

      CALL ComputePrimitive_MHD  &
             ( U(   iNX,iX1,iX2,iX3,iCM_D ),  &
               U(   iNX,iX1,iX2,iX3,iCM_S1),  &
               U(   iNX,iX1,iX2,iX3,iCM_S2),  &
               U(   iNX,iX1,iX2,iX3,iCM_S3),  &
               U(   iNX,iX1,iX2,iX3,iCM_E ),  &
               U(   iNX,iX1,iX2,iX3,iCM_Ne),  &
               U(   iNX,iX1,iX2,iX3,iCM_B1),  &
               U(   iNX,iX1,iX2,iX3,iCM_B2),  &
               U(   iNX,iX1,iX2,iX3,iCM_B3),  &
               U(   iNX,iX1,iX2,iX3,iCM_Chi), &
               P(iPM_D ),  &
               P(iPM_V1),  &
               P(iPM_V2),  &
               P(iPM_V3),  &
               P(iPM_E ),  &
               P(iPM_Ne),  &
               P(iPM_B1),  &
               P(iPM_B2),  &
               P(iPM_B3),  &
               P(iPM_Chi), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(   iNX,iX1,iX2,iX3,iGF_Alpha   ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_1  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_2  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_3  ), &
               EvolveOnlyMagnetic )

      W   = One &
            / SQRT( One &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * P(iPM_V1)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * P(iPM_V2)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * P(iPM_V3)**2 )

      VdotB = ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * P(iPM_V1) * U(iNX,iX1,iX2,iX3,iCM_B1) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  * P(iPM_V2) * U(iNX,iX1,iX2,iX3,iCM_B2) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  * P(iPM_V3) * U(iNX,iX1,iX2,iX3,iCM_B3) )

      WeakDiv(iNX,iX1,iX2,iX3) = D(iNX,iX1,iX2,iX3,iDM_Div)

      dU(iNX,iX1,iX2,iX3,iCM_S1) = dU(iNX,iX1,iX2,iX3,iCM_S1) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                                       * ( U(iNX,iX1,iX2,iX3,iCM_B1) / W**2 &
                                           + P(iPM_V1) * VdotB )

      dU(iNX,iX1,iX2,iX3,iCM_S2) = dU(iNX,iX1,iX2,iX3,iCM_S2) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                                       * ( U(iNX,iX1,iX2,iX3,iCM_B2) / W**2 &
                                           + P(iPM_V2) * VdotB )

      dU(iNX,iX1,iX2,iX3,iCM_S3) = dU(iNX,iX1,iX2,iX3,iCM_S3) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                                       * ( U(iNX,iX1,iX2,iX3,iCM_B3) / W**2 &
                                           + P(iPM_V3) * VdotB )

      dU(iNX,iX1,iX2,iX3,iCM_E)  = dU(iNX,iX1,iX2,iX3,iCM_E ) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) * VdotB

      dU(iNX,iX1,iX2,iX3,iCM_B1) = dU(iNX,iX1,iX2,iX3,iCM_B1) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * ( P(iPM_V1) &
                                         - G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                                           / G(iNX,iX1,iX2,iX3,iGF_Alpha ) )

      dU(iNX,iX1,iX2,iX3,iCM_B2) = dU(iNX,iX1,iX2,iX3,iCM_B2) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * ( P(iPM_V2) &
                                         - G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                                           / G(iNX,iX1,iX2,iX3,iGF_Alpha ) )

      dU(iNX,iX1,iX2,iX3,iCM_B3) = dU(iNX,iX1,iX2,iX3,iCM_B3) &
                                   - WeakDiv(iNX,iX1,iX2,iX3) &
                                     * ( P(iPM_V3) &
                                         - G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                                           / G(iNX,iX1,iX2,iX3,iGF_Alpha ) )

    END DO

    !PRINT*
    !PRINT*, 'Avg. Div.     : ', D(1,iX1,iX2,iX3,iDM_Div)
    !PRINT*, 'Weak Div. Sum : ', WeakDivSum
    !PRINT*, 'Diff.         : ', D(1,iX1,iX2,iX3,iDM_Div) -  WeakDivSum
    !PRINT*

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Powell_Relativistic


  SUBROUTINE ComputeIncrement_Gravity &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

   INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    CALL ComputeIncrement_Gravity_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

  END SUBROUTINE ComputeIncrement_Gravity


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


  SUBROUTINE InitializeIncrement_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1 )

    INTEGER,  INTENT(in)           :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

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

    ALLOCATE( pD_K  (nNodesX_K) )
    ALLOCATE( pV1_K (nNodesX_K) )
    ALLOCATE( pV2_K (nNodesX_K) )
    ALLOCATE( pV3_K (nNodesX_K) )
    ALLOCATE( pE_K  (nNodesX_K) )
    ALLOCATE( pNe_K (nNodesX_K) )
    ALLOCATE( pB1_K (nNodesX_K) )
    ALLOCATE( pB2_K (nNodesX_K) )
    ALLOCATE( pB3_K (nNodesX_K) )
    ALLOCATE( pChi_K(nNodesX_K) )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE InitializeIncrement_MHD


  SUBROUTINE FinalizeIncrement_MHD

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
                pB1_K, pB2_K, pB3_K, pChi_K )

  END SUBROUTINE FinalizeIncrement_MHD


  SUBROUTINE InitializeIncrement_Divergence &
    ( iXP_B0, iXP_E0, nDOFX_X, G_K, G_F, uCM_K, uCM_L, uCM_R )

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
      uCM_K(1:nDOFX, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)-1:iXP_E0(3)+1, &
            1:nCM), &
      uCM_L(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM), &
      uCM_R(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM)

    INTEGER :: iNX, iX1, iX2, iX3, iCM
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

    uD_K  (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_D )
    uS1_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S1)
    uS2_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S2)
    uS3_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S3)
    uE_K  (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_E )
    uNe_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_Ne)
    uB1_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B1)
    uB2_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B2)
    uB3_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B3)
    uChi_K(1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_Chi)

    uD_L  (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_D )
    uS1_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S1)
    uS2_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S2)
    uS3_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S3)
    uE_L  (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_E )
    uNe_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_Ne)
    uB1_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B1)
    uB2_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B2)
    uB3_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B3)
    uChi_L(1:nNodesX_X) => uCM_L(:,:,:,:,iCM_Chi)

    uD_R  (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_D )
    uS1_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S1)
    uS2_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S2)
    uS3_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S3)
    uE_R  (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_E )
    uNe_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_Ne)
    uB1_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B1)
    uB2_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B2)
    uB3_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B3)
    uChi_R(1:nNodesX_X) => uCM_R(:,:,:,:,iCM_Chi)

    ALLOCATE( pD_L  (nNodesX_X) )
    ALLOCATE( pV1_L (nNodesX_X) )
    ALLOCATE( pV2_L (nNodesX_X) )
    ALLOCATE( pV3_L (nNodesX_X) )
    ALLOCATE( pE_L  (nNodesX_X) )
    ALLOCATE( pNe_L (nNodesX_X) )
    ALLOCATE( pB1_L (nNodesX_X) )
    ALLOCATE( pB2_L (nNodesX_X) )
    ALLOCATE( pB3_L (nNodesX_X) )
    ALLOCATE( pChi_L(nNodesX_X) )

    ALLOCATE( pD_R  (nNodesX_X) )
    ALLOCATE( pV1_R (nNodesX_X) )
    ALLOCATE( pV2_R (nNodesX_X) )
    ALLOCATE( pV3_R (nNodesX_X) )
    ALLOCATE( pE_R  (nNodesX_X) )
    ALLOCATE( pNe_R (nNodesX_X) )
    ALLOCATE( pB1_R (nNodesX_X) )
    ALLOCATE( pB2_R (nNodesX_X) )
    ALLOCATE( pB3_R (nNodesX_X) )
    ALLOCATE( pChi_R(nNodesX_X) )

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

  END SUBROUTINE InitializeIncrement_Divergence


  SUBROUTINE FinalizeIncrement_Divergence

    DEALLOCATE( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
                pB1_L, pB2_L, pB3_L, pChi_L )
    DEALLOCATE( pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
                pB1_R, pB2_R, pB3_R, pChi_R )
    DEALLOCATE( IndexTableX_F, IndexTableX_V )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY(   Beta_1_K,   Beta_2_K,   Beta_3_K,  Alpha_K )
    NULLIFY(   Beta_1_F,   Beta_2_F,   Beta_3_F,  Alpha_F )
    NULLIFY( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K )
    NULLIFY( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L )
    NULLIFY( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R )

  END SUBROUTINE FinalizeIncrement_Divergence


  SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX1 )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(inout) :: &
      dGdX1(1:nDOFX,1:nGF,iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3), &
                          iX_B0(1):iX_E0(1))

    INTEGER :: jNX, iNX, iX1, iX2, iX3, iGF
    INTEGER :: nK(3), nGF_K

    REAL(DP) :: G_K_X1 (nDOFX,   nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_Dn_X1(nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))
    REAL(DP) :: G_Up_X1(nDOFX_X1,nGF,iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(1)  :iX_E0(1))

    ASSOCIATE( dX1 => MeshX(1) % Width )

!!$    dGdX1 = Zero
!!$
!!$    DO iX1 = iX_B0(1), iX_E0(1)
!!$    DO iX3 = iX_B0(3), iX_E0(3)
!!$    DO iX2 = iX_B0(2), iX_E0(2)
!!$    DO iGF = 1       , nGF
!!$    DO iNX = 1       , nDOFX
!!$
!!$      DO jNX = 1, nDOFX
!!$
!!$        dGdX1(iNX,iGF,iX2,iX3,iX1) &
!!$          = dGdX1(iNX,iGF,iX2,iX3,iX1) &
!!$              + G(jNX,iX1,iX2,iX3,iGF) * dLXdX1_q(iNX,jNX)
!!$
!!$      END DO
!!$
!!$      dGdX1(iNX,iGF,iX2,iX3,iX1) &
!!$        = dGdX1(iNX,iGF,iX2,iX3,iX1) / MeshX(1) % Width(iX1)
!!$
!!$    END DO
!!$    END DO
!!$    END DO
!!$    END DO
!!$    END DO
!!$
!!$    END ASSOCIATE
!!$
!!$    RETURN

    nK    = iX_E0 - iX_B0 + 1
    nGF_K = nGF * PRODUCT( nK )

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

    ! --- Define field at interface as average of solutions interpolated
    !     to upper and lower faces ---

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

!!$    ! --- Compute metric on faces from scale factors ---
!!$
!!$#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
!!$    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
!!$#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
!!$    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
!!$    !$ACC PRESENT( iX_B0, iX_E0, G_Up_X1, G_Dn_X1 )
!!$#elif defined( THORNADO_OMP    )
!!$    !$OMP PARALLEL DO COLLAPSE(4)
!!$#endif
!!$    DO iX1 = iX_B0(1), iX_E0(1)
!!$    DO iX3 = iX_B0(3), iX_E0(3)
!!$    DO iX2 = iX_B0(2), iX_E0(2)
!!$    DO iNX = 1       , nDOFX_X1
!!$
!!$      G_Up_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
!!$        = MAX( G_Up_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$      G_Up_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
!!$        = MAX( G_Up_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$      G_Up_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
!!$        = MAX( G_Up_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$
!!$      G_Up_X1             (iNX,iGF_SqrtGm,iX2,iX3,iX1) &
!!$        = MAX( G_Up_X1    (iNX,iGF_h_1   ,iX2,iX3,iX1) &
!!$                 * G_Up_X1(iNX,iGF_h_2   ,iX2,iX3,iX1) &
!!$                 * G_Up_X1(iNX,iGF_h_3   ,iX2,iX3,iX1), SqrtTiny )
!!$
!!$      G_Dn_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
!!$        = MAX( G_Dn_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$      G_Dn_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
!!$        = MAX( G_Dn_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$      G_Dn_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
!!$        = MAX( G_Dn_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )
!!$
!!$      G_Dn_X1             (iNX,iGF_SqrtGm,iX2,iX3,iX1) &
!!$        = MAX( G_Dn_X1    (iNX,iGF_h_1   ,iX2,iX3,iX1) &
!!$                 * G_Dn_X1(iNX,iGF_h_2   ,iX2,iX3,iX1) &
!!$                 * G_Dn_X1(iNX,iGF_h_3   ,iX2,iX3,iX1), SqrtTiny )
!!$
!!$    END DO
!!$    END DO
!!$    END DO
!!$    END DO

    ! --- Compute derivatives ---

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

    END ASSOCIATE ! dX1

  END SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X1


  SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX2 )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(inout) :: &
      dGdX2(1:nDOFX,1:nGF,iX_B0(1):iX_E0(1), &
                          iX_B0(3):iX_E0(3), &
                          iX_B0(2):iX_E0(2))

    INTEGER :: jNX, iNX, iX1, iX2, iX3, iGF
    INTEGER :: nK(3), nGF_K

    REAL(DP) :: G_K_X2 (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_Dn_X2(nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))
    REAL(DP) :: G_Up_X2(nDOFX_X2,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(3)  :iX_E0(3), &
                                     iX_B0(2)  :iX_E0(2))

    ASSOCIATE( dX2 => MeshX(2) % Width )

    IF( nDimsX .LT. 2 )THEN

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

    ELSE

      nK    = iX_E0 - iX_B0 + 1
      nGF_K = nGF * PRODUCT( nK )

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

      ! --- Define field at interface as average of solutions interpolated
      !     to upper and lower faces ---

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

!!$      ! --- Compute metric on faces from scale factors ---
!!$
!!$#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
!!$      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
!!$#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
!!$      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
!!$      !$ACC PRESENT( iX_B0, iX_E0, G_Up_X2, G_Dn_X2 )
!!$#elif defined( THORNADO_OMP    )
!!$      !$OMP PARALLEL DO COLLAPSE(4)
!!$#endif
!!$      DO iX2 = iX_B0(2), iX_E0(2)
!!$      DO iX3 = iX_B0(3), iX_E0(3)
!!$      DO iX1 = iX_B0(1), iX_E0(1)
!!$      DO iNX = 1       , nDOFX_X2
!!$
!!$        G_Up_X2         (iNX,iGF_Gm_dd_11,iX1,iX3,iX2) &
!!$          = MAX( G_Up_X2(iNX,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$        G_Up_X2         (iNX,iGF_Gm_dd_22,iX1,iX3,iX2) &
!!$          = MAX( G_Up_X2(iNX,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$        G_Up_X2         (iNX,iGF_Gm_dd_33,iX1,iX3,iX2) &
!!$          = MAX( G_Up_X2(iNX,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$
!!$        G_Up_X2             (iNX,iGF_SqrtGm,iX1,iX3,iX2) &
!!$          = MAX( G_Up_X2    (iNX,iGF_h_1   ,iX1,iX3,iX2) &
!!$                   * G_Up_X2(iNX,iGF_h_2   ,iX1,iX3,iX2) &
!!$                   * G_Up_X2(iNX,iGF_h_3   ,iX1,iX3,iX2), SqrtTiny )
!!$
!!$        G_Dn_X2         (iNX,iGF_Gm_dd_11,iX1,iX3,iX2) &
!!$          = MAX( G_Dn_X2(iNX,iGF_h_1     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$        G_Dn_X2         (iNX,iGF_Gm_dd_22,iX1,iX3,iX2) &
!!$          = MAX( G_Dn_X2(iNX,iGF_h_2     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$        G_Dn_X2         (iNX,iGF_Gm_dd_33,iX1,iX3,iX2) &
!!$          = MAX( G_Dn_X2(iNX,iGF_h_3     ,iX1,iX3,iX2)**2, SqrtTiny )
!!$
!!$        G_Dn_X2             (iNX,iGF_SqrtGm,iX1,iX3,iX2) &
!!$          = MAX( G_Dn_X2    (iNX,iGF_h_1   ,iX1,iX3,iX2) &
!!$                   * G_Dn_X2(iNX,iGF_h_2   ,iX1,iX3,iX2) &
!!$                   * G_Dn_X2(iNX,iGF_h_3   ,iX1,iX3,iX2), SqrtTiny )
!!$
!!$      END DO
!!$      END DO
!!$      END DO
!!$      END DO

      ! --- Compute derivatives ---

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

    END IF ! nDimsX .LT. 2

    END ASSOCIATE ! dX2

  END SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X2


  SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX3 )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(inout) :: &
      dGdX3(1:nDOFX,1:nGF,iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3))

    INTEGER :: jNX, iNX, iX1, iX2, iX3, iGF
    INTEGER :: nK(3), nGF_K

    REAL(DP) :: G_K_X3 (nDOFX,   nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_Dn_X3(nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))
    REAL(DP) :: G_Up_X3(nDOFX_X3,nGF,iX_B0(1)  :iX_E0(1), &
                                     iX_B0(2)  :iX_E0(2), &
                                     iX_B0(3)  :iX_E0(3))

    ASSOCIATE( dX3 => MeshX(3) % Width )

    IF( nDimsX .LT. 3 )THEN

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

    ELSE

      nK    = iX_E0 - iX_B0 + 1
      nGF_K = nGF * PRODUCT( nK )

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

      ! --- Define field at interface as average of solutions interpolated
      !     to upper and lower faces ---

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

!!$      ! --- Compute metric on faces from scale factors ---
!!$
!!$#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
!!$      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
!!$#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
!!$      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
!!$      !$ACC PRESENT( iX_B0, iX_E0, G_Up_X3, G_Dn_X3 )
!!$#elif defined( THORNADO_OMP    )
!!$      !$OMP PARALLEL DO COLLAPSE(4)
!!$#endif
!!$      DO iX3 = iX_B0(3), iX_E0(3)
!!$      DO iX2 = iX_B0(2), iX_E0(2)
!!$      DO iX1 = iX_B0(1), iX_E0(1)
!!$      DO iNX = 1       , nDOFX_X3
!!$
!!$        G_Up_X3         (iNX,iGF_Gm_dd_11,iX1,iX2,iX3) &
!!$          = MAX( G_Up_X3(iNX,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$        G_Up_X3         (iNX,iGF_Gm_dd_22,iX1,iX2,iX3) &
!!$          = MAX( G_Up_X3(iNX,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$        G_Up_X3         (iNX,iGF_Gm_dd_33,iX1,iX2,iX3) &
!!$          = MAX( G_Up_X3(iNX,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$
!!$        G_Up_X3             (iNX,iGF_SqrtGm,iX1,iX2,iX3) &
!!$          = MAX( G_Up_X3    (iNX,iGF_h_1   ,iX1,iX2,iX3) &
!!$                   * G_Up_X3(iNX,iGF_h_2   ,iX1,iX2,iX3) &
!!$                   * G_Up_X3(iNX,iGF_h_3   ,iX1,iX2,iX3), SqrtTiny )
!!$
!!$        G_Dn_X3         (iNX,iGF_Gm_dd_11,iX1,iX2,iX3) &
!!$          = MAX( G_Dn_X3(iNX,iGF_h_1     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$        G_Dn_X3         (iNX,iGF_Gm_dd_22,iX1,iX2,iX3) &
!!$          = MAX( G_Dn_X3(iNX,iGF_h_2     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$        G_Dn_X3         (iNX,iGF_Gm_dd_33,iX1,iX2,iX3) &
!!$          = MAX( G_Dn_X3(iNX,iGF_h_3     ,iX1,iX2,iX3)**2, SqrtTiny )
!!$
!!$        G_Dn_X3             (iNX,iGF_SqrtGm,iX1,iX2,iX3) &
!!$          = MAX( G_Dn_X3    (iNX,iGF_h_1   ,iX1,iX2,iX3) &
!!$                   * G_Dn_X3(iNX,iGF_h_2   ,iX1,iX2,iX3) &
!!$                   * G_Dn_X3(iNX,iGF_h_3   ,iX1,iX2,iX3), SqrtTiny )
!!$
!!$      END DO
!!$      END DO
!!$      END DO
!!$      END DO

      ! --- Compute derivatives ---

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

    END IF ! nDimsX .LT. 3

    END ASSOCIATE ! dX3

  END SUBROUTINE ComputeDerivatives_Geometry_Relativistic_X3


END MODULE MHD_DiscretizationModule_Relativistic
