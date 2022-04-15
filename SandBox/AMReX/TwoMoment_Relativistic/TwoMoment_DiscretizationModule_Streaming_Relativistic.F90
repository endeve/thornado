MODULE TwoMoment_DiscretizationModule_Streaming_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nX,    &
    nDOFZ
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
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
    NodeCoordinate, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
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
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_BoundaryConditionsModule_Relativistic, ONLY: &
    ApplyBoundaryConditions_Euler_Relativistic
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment, &
    ComputePrimitive_TwoMoment_Vector_Richardson, &
    ComputeConserved_TwoMoment, &
    Flux_X1, &
    Flux_E, &
    Source_E, &
    ComputeEddingtonTensorComponents_ud, &
    ComputeHeatFluxTensorComponents_uud_Lagrangian, &
    NumericalFlux_LLF
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MyAmrModule, ONLY: &
    Mass,     &
    R0
  USE UnitsModule,             ONLY: &
    GravitationalConstant

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: UpwindFlux = .FALSE.


  PUBLIC :: ComputeIncrement_TwoMoment_Explicit
  PUBLIC :: ComputeGeometryDerivatives_X0
  PUBLIC :: ComputeGeometryDerivatives_X1
  PUBLIC :: ComputeGeometryDerivatives_X2
  PUBLIC :: ComputeGeometryDerivatives_X3
  PUBLIC :: ComputeWeakDerivatives_X0
  PUBLIC :: ComputeWeakDerivatives_X1
  PUBLIC :: ComputeWeakDerivatives_X2
  PUBLIC :: ComputeWeakDerivatives_X3
  PUBLIC :: ComputeChristoffel

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K, &
    G_Alpha_K,  G_Beta_1_K, G_Beta_2_K, G_Beta_3_K, &
    Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F, &
    G_Alpha_F,  G_Beta_1_F, G_Beta_2_F, G_Beta_3_F


  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    uN_K, uG1_K, uG2_K, uG3_K, &
    uN_L, uG1_L, uG2_L, uG3_L, &
    uN_R, uG1_R, uG2_R, uG3_R

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    uD_K, uI1_K, uI2_K, uI3_K, &
    uD_L, uI1_L, uI2_L, uI3_L, &
    uD_R, uI1_R, uI2_R, uI3_R

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    uV1_F, uV2_F, uV3_F, &
    uV1_K, uV2_K, uV3_K

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    nIterations_L, nIterations_R, nIterations_K

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    PositionIndexZ_F, PositionIndexZ_K

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableZ_F, IndexTableZ_K

  INTEGER :: nZ   (4), nK_X , nK_Z , nNodesX_K, nNodesZ_K
  INTEGER :: nZ_E (4), nE_X , nE_Z , nNodesZ_E
  INTEGER :: nZ_X1(4), nX1_X, nX1_Z, nNodesX_X1, nNodesZ_X1
  INTEGER :: nZ_X2(4), nX2_X, nX2_Z, nNodesX_X2, nNodesZ_X2
  INTEGER :: nZ_X3(4), nX3_X, nX3_Z, nNodesX_X3, nNodesZ_X3

  REAL(DP), PUBLIC :: OffGridFlux_TwoMoment(2*nCR)



CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, &
      Verbose_Option, SuppressBC_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      Verbose_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS
    LOGICAL :: Verbose
    LOGICAL  :: SuppressBC
    REAL(DP) :: dU_R_NEW(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)


    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF
   
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)
    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )


    CALL InitializeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )


    CALL ApplyBoundaryConditions_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )



    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
        dU_R_NEW(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, &
             Verbose_Option = Verbose )
    CALL ComputeIncrement_Divergence_X1_NEW &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R_NEW)
    CALL ComputeIncrement_ObserverCorrections &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, &
             Verbose_Option = Verbose )
    ! --- Multiply Inverse Mass Matrix ---

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE


      IF (dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) .NE. 0.0_DP) THEN
      print*,   (dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)- dU_R_NEW(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS))/ dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)
      ELSE
      print*,   (dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)- dU_R_NEW(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS))
      END IF
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
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


    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit

  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, Verbose_Option )

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
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nZ(4), nZ_X1(4), nV_X1, nV, nX_X1
    INTEGER  :: nIterations
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uCR_X1_L(nCR)
    REAL(DP) :: uCR_X1_R(nCR)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      GX_F(nDOFX_X1,nGF, &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_L(nDOFX_X1,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_R(nDOFX_X1,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      V_u(3,nDOFX_X1, &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_L(nDOF_X1,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_R(nDOF_X1,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      NumericalFlux &
        (nDOF_X1,nCR, &
         iZ_B0(1):iZ_E0(1), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies, &
         iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      Flux_q(nDOFZ,nCR, &
             iZ_B0(1):iZ_E0(1), &
             iZ_B0(3):iZ_E0(3), &
             iZ_B0(4):iZ_E0(4), &
             nSpecies, &
             iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      dU_X1(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E0(2))
    LOGICAL :: Verbose

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
      PRINT*, "      ComputeIncrement_Divergence_X1"
    END IF

    nZ    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nZ_X1 = nZ + [0,1,0,0]    ! Number of X3 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nZ )
    nV_X1 = nCR * nSpecies * PRODUCT( nZ_X1 )
    nX_X1 = PRODUCT( nZ_X1(2:4) ) ! Number of X1 Faces in Position Space

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )
    ! --- Permute Geometry Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)
      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Face States (Average of Left and Right States) ---


    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )




    ! --- Permute Fluid Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iCF = 1, nCF
      DO iNodeX = 1, nDOFX
        
        uCF_K(iNodeX,iCF,iZ3,iZ4,iZ2) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)
      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nCF, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             uCF_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nCF, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             uCF_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    ! --- Compute Face Velocity Components --





    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        ! --- Left State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S1,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S2,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S3,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_E ,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_Ne,iZ3,iZ4,iZ2), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )
        ! --- Right State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S1,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S2,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S3,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_E ,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_Ne,iZ3,iZ4,iZ2), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        V_u(1:3,iNodeX,iZ3,iZ4,iZ2) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )
      END DO

    END DO
    END DO
    END DO
    ! --- Permute Radiation Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nV_X1, nDOFZ, One, L_X1_Up, nDOF_X1, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOF_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nV_X1, nDOFZ, One, L_X1_Dn, nDOF_X1, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOF_X1 )


    ! --- Numerical Flux ---

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iNodeX = 1, nDOFX_X1
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        ! --- Left State Primitive ---
        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )
        ! --- Left State Flux ---

        Flux_L &
          = Flux_X1 &
              ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                uPR_L(iPR_I2), uPR_L(iPR_I3), &
                V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), & 
                GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        CALL ComputeConserved_TwoMoment &
                (uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 uCR_X1_L(iCR_N ), uCR_X1_L(iCR_G1), &
                 uCR_X1_L(iCR_G2), uCR_X1_L(iCR_G3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP,0.0_DP,0.0_DP,                  &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP,0.0_DP,0.0_DP,                  &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )
        ! --- Right State Flux ---

        Flux_R &
          = Flux_X1 &
              ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                uPR_R(iPR_I2), uPR_R(iPR_I3), &
                V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        CALL ComputeConserved_TwoMoment &
                (uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 uCR_X1_R(iCR_N ), uCR_X1_R(iCR_G1), &
                 uCR_X1_R(iCR_G2), uCR_X1_R(iCR_G3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP,0.0_DP,0.0_DP,                  &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        ! --- Numerical Flux ---
        DO iCR = 1, nCR

          NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = NumericalFlux_LLF &
                ( uCR_X1_L(iCR), &
                  uCR_X1_R(iCR), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_X1(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_F(iNodeX,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)
        END DO


      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_X1, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_X1, nDOFZ )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX,iCF_D       ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S1      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S2      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S3      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_E       ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_Ne      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_D       ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_E       ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )
      END DO

    END DO
    END DO
    END DO

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 GX_K(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )

        Flux_K &
          = Flux_X1 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_K(iNodeX,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * Flux_K(iCR)
        END DO
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X1, nDOFZ )


    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X1(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1




  SUBROUTINE ComputeIncrement_Divergence_X1_NEW &
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

    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uPF_K(nPF)
    REAL(DP) :: Flux_L(nCR), uCR_X1_L(nCR)
    REAL(DP) :: Flux_R(nCR), uCR_X1_R(nCR)
    REAL(DP) :: Flux_K(nCR)
    REAL(DP) :: uV1_L, uV2_L, uV3_L
    REAL(DP) :: uV1_R, uV2_R, uV3_R

    ! --- Geometry Fields ---

    REAL(DP) :: &
      GX_K (nDOFX, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nGF)
    REAL(DP) :: &
      GX_F (nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCF_K(nDOFX, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nCF)
    REAL(DP) :: &
      uCF_L(nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCF)
    REAL(DP) :: &
      uCF_R(nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCF)

    ! --- Conserved Radiation Fields ---

    REAL(DP) :: &
      uCR_K(nDOFZ, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nCR)
    REAL(DP) :: &
      uCR_L(nDOF_X1, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCR)
    REAL(DP) :: &
      uCR_R(nDOF_X1, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCR)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux (nDOF_X1,nCR, &
                     iZ_B0(1)  :iZ_E0(1)  , &
                     iZ_B0(3)  :iZ_E0(3)  , &
                     iZ_B0(4)  :iZ_E0(4)  , &
                     nSpecies, &
                     iZ_B0(2)  :iZ_E0(2)+1)
    REAL(DP) :: &
      Flux_q        (nDOFZ,nCR, &
                     iZ_B0(1)  :iZ_E0(1)  , &
                     iZ_B0(3)  :iZ_E0(3)  , &
                     iZ_B0(4)  :iZ_E0(4)  , &
                     nSpecies, &
                     iZ_B0(2)  :iZ_E0(2)  )

    ! --- X1 Increment ---

    REAL(DP) :: &
      dU_X1(nDOFZ,nCR, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)  :iZ_E0(2)  )

    INTEGER :: nZP(4), nZP_X(4)


    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN
    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(2) ; iZP_E0(4) = iZ_E0(2)

    nZP = iZP_E0 - iZP_B0 + 1
    nZP_X = nZP + [0,0,0,1]

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 ) &
    !$OMP MAP( alloc: GX_K, GX_F, uCF_K, uCF_L, uCF_R, &
    !$OMP             uCR_K, uCR_L, uCR_R, NumericalFlux, Flux_q, dU_X1 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0 ) &
    !$ACC CREATE( GX_K, GX_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         uCR_K, uCR_L, uCR_R, NumericalFlux, Flux_q, dU_X1 )
#endif

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X1, nDOF_X1, &
             GX_K, GX_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R )


    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    DO iGF = iGF_h_1, iGF_h_3

      ! --- Face States (Average of Left and Right States) ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO



    ! --- Recompute Geometry from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_F, iZ_B0, iZ_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ2  = iZ_B0(2), iZ_E0(2)+1
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1)**2, SqrtTiny )
        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_2)**2, SqrtTiny )
        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_3)**2, SqrtTiny )
        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_SqrtGm) &
          = SQRT(   GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
                  * GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
                  * GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

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



    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

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


    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1


        ! --- Left State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

        ! --- Right State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )
!check to make sure this is right

      iX_F = iNodeX &
             + ( iZ3 - iZP_B0(2) ) * nDOFX_X1 &
             + ( iZ4 - iZP_B0(3) ) * nDOFX_X1 * nZP_X(2) &
             + ( iZ2 - iZP_B0(4) ) * nDOFX_X1 * nZP_X(2) * nZP_X(3)

      CALL FaceVelocity_X1_NEW &
             ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
               uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3), &
               uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F) )

      END DO

    END DO
    END DO
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



    ! --- Numerical Flux ---

    ! --- Left State Primitive ---



     CALL ComputePrimitive_TwoMoment_Vector_Richardson &
           ( uN_L, uG1_L, uG2_L, uG3_L, &
             uD_L, uI1_L, uI2_L, uI3_L, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             G_Alpha_F, G_Beta_1_F, G_Beta_2_F, G_Beta_3_F, &
             PositionIndexZ_F, &
             nIterations_L )     

    ! --- Right State Primitive ---

     CALL ComputePrimitive_TwoMoment_Vector_Richardson &
           ( uN_R, uG1_R, uG2_R, uG3_R, &
             uD_R, uI1_R, uI2_R, uI3_R, &
             uV1_F, uV2_F, uV3_F, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             G_Alpha_F, G_Beta_1_F, G_Beta_2_F, G_Beta_3_F, &
             PositionIndexZ_F, &
             nIterations_R )     

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
    !$ACC          NumericalFlux, GE, SqrtGm_F, Weights_X1, dZ1, dZ3, dZ4, &
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
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F), &
              G_Alpha_F(iX_F), G_Beta_1_F(iX_F),      &
              G_Beta_2_F(iX_F), G_Beta_3_F(iX_F) )

      CALL ComputeConserved_TwoMoment &
             ( uD_L (iZ_F), uI1_L(iZ_F), &
               uI2_L(iZ_F), uI3_L(iZ_F), &
               uCR_X1_L(iCR_N ), uCR_X1_L(iCR_G1), &
               uCR_X1_L(iCR_G2), uCR_X1_L(iCR_G3), &
               uV1_F(iX_F), Zero, Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F), &
               Zero, Zero, Zero, G_Alpha_F(iX_F), G_Beta_1_F(iX_F),      &
               G_Beta_2_F(iX_F), G_Beta_3_F(iX_F) )

      ! --- Right State Flux ---

      Flux_R &
        = Flux_X1 &
            ( uD_R (iZ_F), uI1_R(iZ_F), &
              uI2_R(iZ_F), uI3_R(iZ_F), &
              uV1_F(iX_F), uV2_F(iX_F), uV3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F), &
              G_Alpha_F(iX_F), G_Beta_1_F(iX_F),      &
              G_Beta_2_F(iX_F), G_Beta_3_F(iX_F) )


      CALL ComputeConserved_TwoMoment &
             ( uD_R (iZ_F), uI1_R(iZ_F), &
               uI2_R(iZ_F), uI3_R(iZ_F), &
               uCR_X1_R(iCR_N ), uCR_X1_R(iCR_G1), &
               uCR_X1_R(iCR_G2), uCR_X1_R(iCR_G3), &
               uV1_F(iX_F), Zero, Zero, &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F), &
               Zero, Zero, Zero, G_Alpha_F(iX_F), G_Beta_1_F(iX_F),      &
               G_Beta_2_F(iX_F), G_Beta_3_F(iX_F) )
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


    END DO



    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_X1, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_X1, nDOFZ )


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


    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        iX_K = iNodeX &
             + ( iZ3 - iZP_B0(2) ) * nDOFX &
             + ( iZ4 - iZP_B0(3) ) * nDOFX * nZP(2) &
             + ( iZ2 - iZP_B0(4) ) * nDOFX * nZP(2) * nZP(3)

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_K(iNodeX      ,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_K(iPF_D ), uV1_K(iX_K), uV2_K(iX_K),  &
                 uV3_K(iX_K), uPF_K(iPF_E), uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )
      END DO

    END DO
    END DO
    END DO


     CALL ComputePrimitive_TwoMoment_Vector_Richardson &
           ( uN_K, uG1_K, uG2_K, uG3_K, &
             uD_K, uI1_K, uI2_K, uI3_K, &
             uV1_K, uV2_K, uV3_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             G_Alpha_K, G_Beta_1_K, G_Beta_2_K, G_Beta_3_K, &
             PositionIndexZ_K, &
             nIterations_K )     
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
              Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K), &
              G_Alpha_K(iX_K), G_Beta_1_K(iX_K),      &
              G_Beta_2_K(iX_K), G_Beta_3_K(iX_K) )

      DO iCR = 1, nCR

        Flux_q(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCR)

      END DO

    END DO


    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCR, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X1, nDOFZ )


#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, dU_X1, iZ_B0, iZ_E0 )
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
            + dU_X1(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)

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
    !$OMP MAP( release: dZ1, dZ3, dZ4, &
    !$OMP               iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0, &
    !$OMP               GX_K, GX_F, uCF_K, uCF_L, uCF_R, &
    !$OMP               uCR_K, uCR_L, uCR_R, NumericalFlux, Flux_q, dU_X1 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ1, dZ3, dZ4, &
    !$ACC         iZ_B0, iZ_E0, iZ_B1, iZ_E1, iZP_B0, iZP_E0, &
    !$ACC         GX_K, GX_F, uCF_K, uCF_L, uCF_R, &
    !$ACC         uCR_K, uCR_L, uCR_R, NumericalFlux, Flux_q, dU_X1 )
#endif

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1_NEW





  SUBROUTINE ComputeIncrement_ObserverCorrections &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, Verbose_Option )

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
    LOGICAL, INTENT(in), OPTIONAL :: &
      Verbose_Option

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX, iNodeX1
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nK(4), nK_Z1(4), nV, nV_Z1
    REAL(DP) :: EdgeEnergyCubed
    REAL(DP) :: Alpha
    REAL(DP) :: AlphaM, AlphaP
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR), E, X1
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: S_E(3)
    REAL(DP) :: &
      uGF_K(nDOFX,nGF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      U_u(1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: &
      U_d(1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX0 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX1 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX2 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX3 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX0_COV &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX1_COV &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX2_COV &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX3_COV &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B1(1):iZ_E1(1))
    REAL(DP) :: uCR_L &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: uCR_R &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: NumericalFlux &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: dU_E &
                  (nDOFZ,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: Flux_q &
                  (nDOFZ,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: &
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: &
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: &
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: &
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP):: W, C

    LOGICAL :: Verbose


    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

    IF (Verbose) THEN
       PRINT*, "      ComputeIncrement_ObserverCorrections"
    END IF

    CALL ComputeGeometryDerivatives_X0 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX0, &
             Verbose_Option = Verbose )

    CALL ComputeGeometryDerivatives_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX1, &
             Verbose_Option = Verbose )

    CALL ComputeGeometryDerivatives_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX2, &
             Verbose_Option = Verbose )

    CALL ComputeGeometryDerivatives_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX3, &
             Verbose_Option = Verbose )
    CALL ComputeChristoffel( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                             dG_dd_dX0, dG_dd_dX1, dG_dd_dX2, dG_dd_dX3, Gamma_udd )

    CALL ComputeWeakDerivatives_X0 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, Gamma_udd, dU_d_dX0, dU_d_dX0_COV,  &
             Verbose_Option = Verbose )
    CALL ComputeWeakDerivatives_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, Gamma_udd, dU_d_dX1, dU_d_dX1_COV, &
             Verbose_Option = Verbose )
    
    CALL ComputeWeakDerivatives_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, Gamma_udd,  dU_d_dX2, dU_d_dX2_COV, &
             Verbose_Option = Verbose )

    CALL ComputeWeakDerivatives_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, Gamma_udd, dU_d_dX3, dU_d_dX3_COV, &
             Verbose_Option = Verbose )

    CALL ComputeFourVelocity &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, U_u, U_d, uPF_K, &
             Verbose_Option = Verbose  )




    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_Z1 = nK + [1,0,0,0]    ! Number of Z1 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nK )
    nV_Z1 = nCR * nSpecies * PRODUCT( nK_Z1 )

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    ! --- Permute Geometry Fields ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iGF = 1, nGF

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iGF,iZ2,iZ3,iZ4) &
          = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iCF = 1, nCF

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ3,iZ4) &
          = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Radiation Fields ---

    DO iZ1 = iZ_B1(1), iZ_E1(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ
        uCR_K(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute Primitive Fluid ---


    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_E, nV_Z1, nDOFZ, One, L_E_Up, nDOF_E, &
             uCR_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOF_E )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_E, nV_Z1, nDOFZ, One, L_E_Dn, nDOF_E, &
             uCR_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOF_E )

    ! --- Numerical Flux ---

    DO iZ1 = iZ_B0(1), iZ_E1(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNode = 1, nDOF_E ! = nDOFX
        ! --- Left State Primitive --
        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNode,iCR_N       ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G1      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G2      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G3      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &   
                 uGF_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_1,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_2,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_3,iZ2,iZ3,iZ4) )
        ! --- Left State Flux ---

        Flux_L &
          = Flux_E( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                    uPR_L(iPR_I2), uPR_L(iPR_I3), &
                    uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_1,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_2,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_3,iZ2,iZ3,iZ4), &
                    U_u(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX0_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3_COV(iNode,:,iZ2,iZ3,iZ4) )
        ! --- Right State Primitive ---
        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNode,iCR_N       ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G1      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G2      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G3      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &   
                 uGF_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_1,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_2,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Beta_3,iZ2,iZ3,iZ4) )

        ! --- Right State Flux ---

        Flux_R &
          = Flux_E( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                    uPR_R(iPR_I2), uPR_R(iPR_I3), &
                    uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_1,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_2,iZ2,iZ3,iZ4), &
                    uGF_K(iNode,iGF_Beta_3,iZ2,iZ3,iZ4), &
                    U_u(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX0_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2_COV(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3_COV(iNode,:,iZ2,iZ3,iZ4) )

        ! --- Numerical Flux ---

        EdgeEnergyCubed &
          = ( MeshE % Center(iZ1) - Half * MeshE % Width(iZ1) )**3

        IF( UpwindFlux )THEN

          AlphaP = 0.0_DP
          AlphaM = 1.0_DP 


          DO iCR = 1, nCR

            NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
              = AlphaP * Flux_L(iCR) + AlphaM * Flux_R(iCR)

            NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
              = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                  * EdgeEnergyCubed * Weights_E(iNode) &
                  * uGF_K(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                  * NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1)

          END DO

        ELSE

          W = 1.0_DP - ( uPF_K(iNode,iPF_V1,iZ2,iZ3,iZ4)**2 * uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
                     + uPF_K(iNode,iPF_V2,iZ2,iZ3,iZ4)**2 * uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
                     + uPF_K(iNode,iPF_V3,iZ2,iZ3,iZ4)**2 * uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )
       
          W = 1.0_DP / SQRT(W)
           
          C = 1.0_DP
          CALL ComputeAlpha( uPF_K(iNode,iPF_V1,iZ2,iZ3,iZ4), &
                             uPF_K(iNode,iPF_V2,iZ2,iZ3,iZ4), &
                             uPF_K(iNode,iPF_V3,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                             U_u(iNode,:,iZ2,iZ3,iZ4), &
                             dU_d_dX0(iNode,:,iZ2,iZ3,iZ4), &
                             dU_d_dX1(iNode,:,iZ2,iZ3,iZ4), &
                             dU_d_dX2(iNode,:,iZ2,iZ3,iZ4), &
                             dU_d_dX3(iNode,:,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Beta_1,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Beta_2,iZ2,iZ3,iZ4), &
                             uGF_K(iNode,iGF_Beta_3,iZ2,iZ3,iZ4), &
                             C, Alpha )
          DO iCR = 1, nCR


            NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
              = NumericalFlux_LLF &
                  ( W * uPR_L(iCR),  &
                    W * uPR_R(iCR), &
                    Flux_L(iCR), Flux_R(iCR), Alpha )

            NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
              = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                  * EdgeEnergyCubed * Weights_E(iNode) &
                  * uGF_K(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                  * NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1)
          END DO


        END IF

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Surface Contributions ---

    ! --- Contributions from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_E, + One, L_E_Dn, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), &
             nDOF_E, Zero, dU_E, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_E, - One, L_E_Up, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)+1), &
             nDOF_E, One,  dU_E, nDOFZ )


    !--------------------
    ! --- Volume Term ---
    !--------------------

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_K(iPR_D ), &
                 uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), &
                 uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &   
                 uGF_K(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4) )

        Flux_K &
          = Flux_E( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                    uPR_K(iPR_I2), uPR_K(iPR_I3), &
                    uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4), &
                    U_u(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX0_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3_COV(iNodeX,:,iZ2,iZ3,iZ4) )
        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
            = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) * Weights_q(iNodeZ) &
                * GE(iNodeE,iZ1,iGE_Ep3) &
                * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                * Flux_K(iCR)
        END DO
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdE_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_E, nDOFZ )


    !--------------------------------------------------
    !--- Volume Sources (Number Flux Equation Only) ---
    !--------------------------------------------------

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE




        E = NodeCoordinate( MeshE, iZ1, iNodeE ) 

        W = 1.0_DP - ( uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4)**2 * uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
                     + uPF_K(iNodeX,iPF_V2,iZ2,iZ3,iZ4)**2 * uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
                     + uPF_K(iNodeX,iPF_V3,iZ2,iZ3,iZ4)**2 * uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )
       
        W = 1.0_DP / SQRT(W)
        
        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_K(iPR_D ), &
                 uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), &
                 uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &   
                 uGF_K(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4) )

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeX1 )

         S_E = Source_E( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                    uPR_K(iPR_I2), uPR_K(iPR_I3), &
                    uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                    uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4), &
                    uGF_K(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4), &
                    U_u(iNodeX,:,iZ2,iZ3,iZ4), &
                    U_d(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX0(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX0_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3_COV(iNodeX,:,iZ2,iZ3,iZ4), &
                    dG_dd_dX1(iNodeX,:,:,iZ2,iZ3,iZ4), &
                    dG_dd_dX2(iNodeX,:,:,iZ2,iZ3,iZ4), &
                    dG_dd_dX3(iNodeX,:,:,iZ2,iZ3,iZ4), E, X1 )
        ! --- iCR_G1 ---

        dU_E(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * S_E(1)
        
        ! --- iCR_G2 ---
        
        dU_E(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * S_E(2)

        ! --- iCR_G3 ---

        dU_E(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * S_E(3)
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_E(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1)
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ1, etc    END ASSOCIATE ! dZ1, etc..

  END SUBROUTINE ComputeIncrement_ObserverCorrections 



  FUNCTION FaceVelocity_X1 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X1(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X1(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X1(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X1(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X1




  FUNCTION FaceVelocity_X2 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X2(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X2(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X2(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X2(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X2

  FUNCTION FaceVelocity_X3 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X3(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X3(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X3(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X3(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X3

  SUBROUTINE FaceVelocity_X1_NEW &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R, V1_F, V2_F, V3_F )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in)  :: V1_R, V2_R, V3_R
    REAL(DP), INTENT(out) :: V1_F, V2_F, V3_F

    ! --- Average Left and Right States ---

    V1_F = Half * ( V1_L + V1_R )
    V2_F = Half * ( V2_L + V2_R )
    V3_F = Half * ( V3_L + V3_R )

    RETURN
  END SUBROUTINE FaceVelocity_X1_NEW

  
  SUBROUTINE ComputeFourVelocity &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, U_u, U_d, uPF, Verbose_Option  )


    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      U_u(1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      U_d(1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      uPF(nDOFX,nPF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4
    REAL(DP) :: W, V_0, V1, V2, V3

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)

      DO iNodeX = 1, nDOFX
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF(iNodeX,iZ2,iZ3,iZ4,iCF_D ), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S1), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S2), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S3), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_E ), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_Ne), &
                 uPF(iNodeX,iPF_D,iZ2,iZ3,iZ4), &
                 uPF(iNodeX,iPF_V1,iZ2,iZ3,iZ4), &
                 uPF(iNodeX,iPF_V2,iZ2,iZ3,iZ4), &
                 uPF(iNodeX,iPF_V3,iZ2,iZ3,iZ4), &
                 uPF(iNodeX,iPF_E,iZ2,iZ3,iZ4), &
                 uPF(iNodeX,iPF_Ne,iZ2,iZ3,iZ4), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )
        
        V1 = uPF(iNodeX,iPF_V1,iZ2,iZ3,iZ4)
        V2 = uPF(iNodeX,iPF_V2,iZ2,iZ3,iZ4)
        V3 = uPF(iNodeX,iPF_V3,iZ2,iZ3,iZ4)
        W =1.0_DP /SQRT( 1.0_DP - ( V1**2 * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) &
           + V2**2 * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) + V3**2 * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) ) )

        V_0 = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1) * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * V1 &
            + GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2) * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * V2 &
            + GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * V3              
 
        U_u(iNodeX,0,iZ2,iZ3,iZ4) = W * ( 1.0_DP / GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) ) 
        U_u(iNodeX,1,iZ2,iZ3,iZ4) = W * ( V1 - GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1) / GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) )
        U_u(iNodeX,2,iZ2,iZ3,iZ4) = W * ( V2 - GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2) / GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) )
        U_u(iNodeX,3,iZ2,iZ3,iZ4) = W * ( V3 - GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) / GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) )


        U_d(iNodeX,0,iZ2,iZ3,iZ4) = W * ( - GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) + V_0 )
        U_d(iNodeX,1,iZ2,iZ3,iZ4) = W * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * V1
        U_d(iNodeX,2,iZ2,iZ3,iZ4) = W * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * V2
        U_d(iNodeX,3,iZ2,iZ3,iZ4) = W * GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * V3


      END DO
    END DO
    END DO
    END DO


  END SUBROUTINE



  SUBROUTINE ComputeWeakDerivatives_X0 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, Gamma_udd, dU_d_dX0, dU_d_dX0_COV, Verbose_Option  )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX0 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dU_d_dX0_COV &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    REAL(DP) :: &
      uPF_K(nPF), &
      U_d(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: V_u_K(3), V_d_K(3), A_K, B_u_K(3), B_d_K(3), W_K, Vsq_K
    LOGICAL :: Verbose
    INTEGER :: iNodeX, mu, iZ2, iZ3, iZ4


    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X0"
    END IF

    dU_d_dX0 = Zero
    dU_d_dX0_COV = Zero

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF(iNodeX,iZ2,iZ3,iZ4,iCF_D ), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S1), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S2), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_S3), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_E ), &
                 uCF(iNodeX,iZ2,iZ3,iZ4,iCF_Ne), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

        V_u_K(1) = uPF_K(iPF_V1)
 
        V_u_K(2) = uPF_K(iPF_V2) 

        V_u_K(3) = uPF_K(iPF_V3) 

        V_d_K(1) &
          = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) &
              * V_u_K(1)

        V_d_K(2) &
          = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) &
              * V_u_K(2)

        V_d_K(3) &
          = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) &
              * V_u_K(3)

        A_K = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) = GX (iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)
        
        Vsq_K = V_u_K(1) * V_d_K(1) &
                  + V_u_K(2) * V_d_K(2) &
                  + V_u_K(3) * V_d_K(3)
     
        W_K = 1.0_DP / SQRT( 1.0_DP - Vsq_K )

       

        U_d(iNodeX,1,iZ2,iZ3,iZ4) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d(iNodeX,2,iZ2,iZ3,iZ4) = W_K * V_d_K(1) 

        U_d(iNodeX,3,iZ2,iZ3,iZ4) = W_K * V_d_K(2) 

        U_d(iNodeX,4,iZ2,iZ3,iZ4) = W_K * V_d_K(3) 


      END DO

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO mu = 0,3

        dU_d_dX0_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
        = dU_d_dX0_COV(iNodeX,0,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,0,0,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)

        dU_d_dX0_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
        = dU_d_dX0_COV(iNodeX,1,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,0,1,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)
    
        dU_d_dX0_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
        = dU_d_dX0_COV(iNodeX,2,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,0,2,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)
    
        dU_d_dX0_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
        = dU_d_dX0_COV(iNodeX,3,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,0,3,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)

      END DO
      END DO

    END DO
    END DO
    END DO
  END SUBROUTINE ComputeWeakDerivatives_X0 


  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, Gamma_udd, dU_d_dX1, dU_d_dX1_COV, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX1 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dU_d_dX1_COV &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: nK(4), nK_X1(4), nX, nX_X1
    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4, iCF, iGF, mu
    REAL(DP) :: &
      uPF_K(nPF), uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: &
      GX_K(nDOFX   ,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)), &
      GX_F(nDOFX_X1,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1)
    REAL(DP) :: &
      uCF_K(nDOFX   ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B1(2):iZ_E1(2)  ,nCF), &
      uCF_L(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF), &
      uCF_R(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF)
    REAL(DP) :: Vsq_X1, Vsq_K
    REAL(DP) :: V_u_K(3), V_d_K(3), V_u_X1(3), V_d_X1(3), W_X1, W_K
    REAL(DP) :: B_u_X1(3), B_d_X1(3), A_X1, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      U_d_X1(nDOFX_X1,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      U_d_K(nDOFX,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      U_d(nDOFX,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      dU_d_dX1_Temp &
         (1:nDOFX,4, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    LOGICAL :: Verbose


    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      dU_d_dX1 = Zero
      dU_d_dX1_COV = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X1"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X1 = nK + [0,1,0,0]    ! Number of X1 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X1 = PRODUCT( nK_X1(2:4) ) ! Number of X1 Faces in Position Space
    
    ! --- Permute Geometry Fields ---

    DO iGF = 1, nGF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX


        GX_K(iNodeX,iGF,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)




      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Geometry Fields ---


      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nGF*nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
               GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nGF*nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX, Half, &
               GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )






    ! --- Compute Metric Components ---

    ! --- Permute Fluid Fields ---

    DO iCF = 1, nCF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

    END DO
  

    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        ! --- Left States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )


        V_u_X1(1:3) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )


        V_d_X1(1) &
          = GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
              * V_u_X1(1)

        V_d_X1(2) &
          = GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
              * V_u_X1(2)

        V_d_X1(3) &
          = GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) &
              * V_u_X1(3)

        Vsq_X1 = V_u_X1(1) * V_d_X1(1) &  
                    + V_u_X1(2) * V_d_X1(2) &  
                    + V_u_X1(3) * V_d_X1(3) 
    
        W_X1 = WeightsX_X1(iNodeX)  / SQRT(1.0_DP - Vsq_X1) 

        A_X1 = GX_F(iNodeX,iGF_Alpha,iZ3,iZ4,iZ2)

        B_u_X1(1) = GX_F(iNodeX,iGF_Beta_1,iZ3,iZ4,iZ2)
        B_u_X1(2) = GX_F(iNodeX,iGF_Beta_2,iZ3,iZ4,iZ2)
        B_u_X1(3) = GX_F(iNodeX,iGF_Beta_3,iZ3,iZ4,iZ2)

        B_d_X1(1) = GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) * B_u_X1(1)
        B_d_X1(2) = GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) * B_u_X1(2)
        B_d_X1(3) = GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) * B_u_X1(3)



        U_d_X1(iNodeX,1,iZ3,iZ4,iZ2) = W_X1 * ( - A_X1 + B_d_X1(1) * V_u_X1(1) &
                                     + B_d_X1(2) * V_u_X1(2) + B_d_X1(3) * V_u_X1(3) )
        
        U_d_X1(iNodeX,2,iZ3,iZ4,iZ2) = W_X1 * V_d_X1(1) 

        U_d_X1(iNodeX,3,iZ3,iZ4,iZ2) = W_X1 * V_d_X1(2) 

        U_d_X1(iNodeX,4,iZ3,iZ4,iZ2) = W_X1 * V_d_X1(3)





     END DO

    END DO
    END DO
    END DO










    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             U_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dU_d_dX1_Temp, nDOFX )


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             U_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dU_d_dX1_Temp, nDOFX )





    ! --- Volume Term ---
    ! -------------------

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        V_u_K(1) = uPF_K(iPF_V1)
 
        V_u_K(2) = uPF_K(iPF_V2) 

        V_u_K(3) = uPF_K(iPF_V3) 

        V_d_K(1) &
          = GX_K(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
              * V_u_K(1)

        V_d_K(2) &
          = GX_K(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
              * V_u_K(2)

        V_d_K(3) &
          = GX_K(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) &
              * V_u_K(3)

        A_K = GX_K(iNodeX,iGF_Alpha,iZ3,iZ4,iZ2)

        B_u_K(1) =  GX_K(iNodeX,iGF_Beta_1,iZ3,iZ4,iZ2)
        B_u_K(2) =  GX_K(iNodeX,iGF_Beta_2,iZ3,iZ4,iZ2)
        B_u_K(3) =  GX_K(iNodeX,iGF_Beta_3,iZ3,iZ4,iZ2)

        B_d_K(1) =  GX_K(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) * B_u_K(1)
        B_d_K(2) =  GX_K(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) * B_u_K(2)
        B_d_K(3) =  GX_K(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) * B_u_K(3)
        
        Vsq_K = V_u_K(1) * V_d_K(1) &
                  + V_u_K(2) * V_d_K(2) &
                  + V_u_K(3) * V_d_K(3)
     
        W_K = 1.0_DP / SQRT( 1.0_DP - Vsq_K )

       
        U_d(iNodeX,1,iZ3,iZ4,iZ2) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d(iNodeX,2,iZ3,iZ4,iZ2) = W_K * V_d_K(1) 

        U_d(iNodeX,3,iZ3,iZ4,iZ2) = W_K * V_d_K(2) 

        U_d(iNodeX,4,iZ3,iZ4,iZ2) = W_K * V_d_K(3) 



        U_d_K(iNodeX,1,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * U_d(iNodeX,1,iZ3,iZ4,iZ2)
        U_d_K(iNodeX,2,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * U_d(iNodeX,2,iZ3,iZ4,iZ2)
        U_d_K(iNodeX,3,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * U_d(iNodeX,3,iZ3,iZ4,iZ2)
        U_d_K(iNodeX,4,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * U_d(iNodeX,3,iZ3,iZ4,iZ2)









      END DO

    END DO
    END DO
    END DO






    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX1_Temp, nDOFX )






    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX





        dU_d_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
         = dU_d_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


        dU_d_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
         = dU_d_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dU_d_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
         = dU_d_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dU_d_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
         = dU_d_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )






      END DO

    END DO
    END DO
    END DO
    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        dU_d_dX1(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2)

        dU_d_dX1(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2)

        dU_d_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2)

        dU_d_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2)




        dU_d_dX1_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2)

        dU_d_dX1_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2)

        dU_d_dX1_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2)

        dU_d_dX1_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2)


      END DO

    END DO
    END DO
    END DO


    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO mu = 0,3

    
        dU_d_dX1_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
        = dU_d_dX1_COV(iNodeX,0,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,1,0,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ3,iZ4,iZ2)

        dU_d_dX1_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
        = dU_d_dX1_COV(iNodeX,1,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,1,1,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ3,iZ4,iZ2)
    
        dU_d_dX1_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
        = dU_d_dX1_COV(iNodeX,2,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,1,2,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ3,iZ4,iZ2)
    
        dU_d_dX1_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
        = dU_d_dX1_COV(iNodeX,3,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,1,3,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ3,iZ4,iZ2)

      END DO
      END DO

    END DO
    END DO
    END DO


  END SUBROUTINE ComputeWeakDerivatives_X1

  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, Gamma_udd, dU_d_dX2, dU_d_dX2_COV, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(in) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX2 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX2_COV &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4, iGF, iCF, mu
    INTEGER  :: nK(4), nK_X2(4), nX, nX_X2
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      GX_F(nDOFX_X2,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_L(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_R(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: Vsq_X2, Vsq_K
    REAL(DP) :: V_u_K(3), V_d_K(3), V_u_X2(3), V_d_X2(3), W_X2, W_K
    REAL(DP) :: B_u_X2(3), B_d_X2(3), A_X2, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      U_d_X2(nDOFX_X2,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
             iZ_B0(3):iZ_E0(3)+1), &
      U_d_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3)), &
      U_d(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      dU_d_dX2_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    LOGICAL :: Verbose

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      dU_d_dX2 = Zero
      dU_d_dX2_COV = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X2"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X2 = nK + [0,0,1,0]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X2 = PRODUCT( nK_X2(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)



      END DO
      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---





      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF*nX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
               GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nGF*nX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX, Half, &
               GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )









    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---


    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ4,iZ3) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO


    ! --- Interpolate Fluid Fields ---


    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )


    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2

        ! --- Left States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )


        V_u_X2(1:3) &
          = FaceVelocity_X2 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

        V_d_X2(1) &
          = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
              * V_u_X2(1)

        V_d_X2(2) &
          = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
              * V_u_X2(2)

        V_d_X2(3) &
          = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
              * V_u_X2(3)

        Vsq_X2 = V_u_X2(1) * V_d_X2(1) &  
               + V_u_X2(2) * V_d_X2(2) &  
               + V_u_X2(3) * V_d_X2(3) 

        A_X2 = GX_F(iNodeX,iGF_Alpha,iZ2,iZ4,iZ3)

        B_u_X2(1) = GX_F(iNodeX,iGF_Beta_1,iZ2,iZ4,iZ3)
        B_u_X2(2) = GX_F(iNodeX,iGF_Beta_2,iZ2,iZ4,iZ3)
        B_u_X2(3) = GX_F(iNodeX,iGF_Beta_3,iZ2,iZ4,iZ3)

        B_d_X2(1) = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) * B_u_X2(1)
        B_d_X2(2) = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) * B_u_X2(2)
        B_d_X2(3) = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) * B_u_X2(3)




        U_d_X2(iNodeX,1,iZ2,iZ4,iZ3) = W_X2 * ( - A_X2 + B_d_X2(1) * V_u_X2(1) &
                                     + B_d_X2(2) * V_u_X2(2) + B_d_X2(3) * V_u_X2(3) )
        
        U_d_X2(iNodeX,2,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(1) 

        U_d_X2(iNodeX,3,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(2) 

        U_d_X2(iNodeX,4,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(3)



      END DO

    END DO
    END DO
    END DO





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             U_d_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dU_d_dX2_Temp, nDOFX )


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             U_d_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dU_d_dX2_Temp, nDOFX )





    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )




        V_u_K(1) = uPF_K(iPF_V1)
 
        V_u_K(2) = uPF_K(iPF_V2) 

        V_u_K(3) = uPF_K(iPF_V3) 

        V_d_K(1) &
          = GX_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
              * V_u_K(1)

        V_d_K(2) &
          = GX_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
              * V_u_K(2)

        V_d_K(3) &
          = GX_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
              * V_u_K(3)

        Vsq_K = V_u_K(1) * V_d_K(1) &
                  + V_u_K(2) * V_d_K(2) &
                  + V_u_K(3) * V_d_K(3)
    
        A_K = GX_K(iNodeX, iGF_Alpha, iZ2,iZ4,iZ3)

        B_u_K(1) =  GX_K(iNodeX, iGF_Beta_1,iZ2,iZ4,iZ3)
        B_u_K(2) =  GX_K(iNodeX, iGF_Beta_2,iZ2,iZ4,iZ3)
        B_u_K(3) =  GX_K(iNodeX, iGF_Beta_3,iZ2,iZ4,iZ3)

        B_d_K(1) =  GX_K(iNodeX, iGF_Gm_dd_11,iZ2,iZ4,iZ3) * B_u_K(1)
        B_d_K(2) =  GX_K(iNodeX, iGF_Gm_dd_22,iZ2,iZ4,iZ3) * B_u_K(2)
        B_d_K(3) =  GX_K(iNodeX, iGF_Gm_dd_33,iZ2,iZ4,iZ3) * B_u_K(3)


 
        W_K = 1.0_DP / SQRT( 1.0_DP - Vsq_K )

       

        U_d(iNodeX,1,iZ2,iZ4,iZ3) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d(iNodeX,2,iZ2,iZ4,iZ3) = W_K * V_d_K(1) 

        U_d(iNodeX,3,iZ2,iZ4,iZ3) = W_K * V_d_K(2) 

        U_d(iNodeX,4,iZ2,iZ4,iZ3) = W_K * V_d_K(3)


        U_d_K(iNodeX,1,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * U_d(iNodeX,1,iZ2,iZ4,iZ3)
        U_d_K(iNodeX,2,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * U_d(iNodeX,2,iZ2,iZ4,iZ3)
        U_d_K(iNodeX,3,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * U_d(iNodeX,3,iZ2,iZ4,iZ3)
        U_d_K(iNodeX,4,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * U_d(iNodeX,4,iZ2,iZ4,iZ3)

      END DO

    END DO
    END DO
    END DO






    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX2_Temp, nDOFX )









    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX




        dU_d_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
         = dU_d_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )


        dU_d_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
         = dU_d_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dU_d_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
         = dU_d_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dU_d_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
         = dU_d_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )



      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        dU_d_dX2(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3)




        dU_d_dX2_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3)

        dU_d_dX2_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3)

        dU_d_dX2_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3)

        dU_d_dX2_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3)
      END DO

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO mu = 0,3

    
        dU_d_dX2_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
        = dU_d_dX2_COV(iNodeX,0,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,2,0,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ4,iZ3)

        dU_d_dX2_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
        = dU_d_dX2_COV(iNodeX,1,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,2,1,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ4,iZ3)
    
        dU_d_dX2_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
        = dU_d_dX2_COV(iNodeX,2,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,2,2,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ4,iZ3)
    
        dU_d_dX2_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
        = dU_d_dX2_COV(iNodeX,3,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,2,3,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ4,iZ3)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeWeakDerivatives_X2


  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, Gamma_udd, dU_d_dX3, dU_d_dX3_COV, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(in) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX3 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_d_dX3_COV &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4, iGF, iCF, mu
    INTEGER  :: nK(4), nK_X3(4), nX, nX_X3
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      GX_F(nDOFX_X3,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_L(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_R(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: Vsq_X3, Vsq_K
    REAL(DP) :: V_u_K(3), V_d_K(3), V_u_X3(3), V_d_X3(3), W_X3, W_K
    REAL(DP) :: B_u_X3(3), B_d_X3(3), A_X3, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      U_d_X3(nDOFX_X3,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
             iZ_B0(4):iZ_E0(4)+1), &
      U_d_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4)), &
      U_d(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX3_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 

    LOGICAL :: Verbose

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      dU_d_dX3 = Zero
      dU_d_dX3_COV = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X3"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X3 = nK + [0,0,0,1]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X3 = PRODUCT( nK_X3(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)



      END DO
      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---



      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF*nX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
               GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nGF*nX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX, Half, &
               GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )









    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---


    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ3,iZ4) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO



    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )


    ! --- Compute Face Velocity Components ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        ! --- Left States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )



        V_u_X3(1:3) &
          = FaceVelocity_X3 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )


        V_d_X3(1) &
          = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
              * V_u_X3(1)

        V_d_X3(2) &
          = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
              * V_u_X3(2)

        V_d_X3(3) &
          = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
              * V_u_X3(3)

        Vsq_X3 = V_u_X3(1) * V_d_X3(1) &  
                    + V_u_X3(2) * V_d_X3(2) &  
                    + V_u_X3(3) * V_d_X3(3) 
    
        W_X3 = WeightsX_X3(iNodeX)  / SQRT(1.0_DP - Vsq_X3) 

        A_X3 = GX_F(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4)

        B_u_X3(1) = GX_F(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4)
        B_u_X3(2) = GX_F(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4)
        B_u_X3(3) = GX_F(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4)

        B_d_X3(1) = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) * B_u_X3(1)
        B_d_X3(2) = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) * B_u_X3(2)
        B_d_X3(3) = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) * B_u_X3(3)





        U_d_X3(iNodeX,1,iZ2,iZ3,iZ4) = W_X3 * ( - A_X3 + B_d_X3(1) * V_u_X3(1) &
                                     + B_d_X3(2) * V_u_X3(2) + B_d_X3(3) * V_u_X3(3) )
        
        U_d_X3(iNodeX,2,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(1) 

        U_d_X3(iNodeX,3,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(2) 

        U_d_X3(iNodeX,4,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(3)





      END DO

    END DO
    END DO
    END DO





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             U_d_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dU_d_dX3_Temp, nDOFX )


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             U_d_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dU_d_dX3_Temp, nDOFX )






    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )






        V_u_K(1) = uPF_K(iPF_V1)
 
        V_u_K(2) = uPF_K(iPF_V2) 

        V_u_K(3) = uPF_K(iPF_V3) 

        V_d_K(1) &
          = GX_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
              * V_u_K(1)

        V_d_K(2) &
          = GX_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
              * V_u_K(2)

        V_d_K(3) &
          = GX_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
              * V_u_K(3)

        A_K = GX_K(iNodeX,iGF_Alpha,iZ2,iZ3,iZ4)

        B_u_K(1) =  GX_K(iNodeX,iGF_Beta_1,iZ2,iZ3,iZ4)
        B_u_K(2) =  GX_K(iNodeX,iGF_Beta_2,iZ2,iZ3,iZ4)
        B_u_K(3) =  GX_K(iNodeX,iGF_Beta_3,iZ2,iZ3,iZ4)

        B_d_K(1) =  GX_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) * B_u_K(1)
        B_d_K(2) =  GX_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) * B_u_K(2)
        B_d_K(3) =  GX_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) * B_u_K(3)
        
        Vsq_K = V_u_K(1) * V_d_K(1) &
                  + V_u_K(2) * V_d_K(2) &
                  + V_u_K(3) * V_d_K(3)
     
        W_K = 1.0_DP / SQRT( 1.0_DP - Vsq_K )

       

        U_d(iNodeX,1,iZ2,iZ3,iZ4) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d(iNodeX,2,iZ2,iZ3,iZ4) = W_K * V_d_K(1) 

        U_d(iNodeX,3,iZ2,iZ3,iZ4) = W_K * V_d_K(2) 

        U_d(iNodeX,4,iZ2,iZ3,iZ4) = W_K * V_d_K(3) 

        U_d_K(iNodeX,1,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * U_d(iNodeX,1,iZ2,iZ3,iZ4)
        U_d_K(iNodeX,2,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * U_d(iNodeX,2,iZ2,iZ3,iZ4)
        U_d_K(iNodeX,3,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * U_d(iNodeX,3,iZ2,iZ3,iZ4)
        U_d_K(iNodeX,4,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * U_d(iNodeX,4,iZ2,iZ3,iZ4)



      END DO

    END DO
    END DO
    END DO





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX3_Temp, nDOFX )








    ASSOCIATE( dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX





        dU_d_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
         = dU_d_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


        dU_d_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
         = dU_d_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dU_d_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
         = dU_d_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dU_d_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
         = dU_d_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        dU_d_dX3(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4)




        dU_d_dX3_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4)

        dU_d_dX3_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4)

        dU_d_dX3_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4)

        dU_d_dX3_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4)
      END DO

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO mu = 0,3

    
        dU_d_dX3_COV(iNodeX,0,iZ2,iZ3,iZ4)  &
        = dU_d_dX3_COV(iNodeX,0,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,3,0,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)

        dU_d_dX3_COV(iNodeX,1,iZ2,iZ3,iZ4)  &
        = dU_d_dX3_COV(iNodeX,1,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,3,1,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)
    
        dU_d_dX3_COV(iNodeX,2,iZ2,iZ3,iZ4)  &
        = dU_d_dX3_COV(iNodeX,2,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,3,2,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)
    
        dU_d_dX3_COV(iNodeX,3,iZ2,iZ3,iZ4)  &
        = dU_d_dX3_COV(iNodeX,3,iZ2,iZ3,iZ4) &
        - Gamma_udd(iNodeX,mu,3,3,iZ2,iZ3,iZ4) * U_d(iNodeX,mu+1,iZ2,iZ3,iZ4)

      END DO
      END DO

    END DO
    END DO
    END DO



  END SUBROUTINE ComputeWeakDerivatives_X3



  SUBROUTINE ComputeGeometryDerivatives_X0 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX0, Verbose_Option  )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    dG_dd_dX0 = Zero
  END SUBROUTINE ComputeGeometryDerivatives_X0 


  SUBROUTINE ComputeGeometryDerivatives_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX1, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: nK(4), nK_X1(4), nX, nX_X1
    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4
    REAL(DP) :: B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      G_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)), &
      G_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)), &
      H_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1)
    REAL(DP) :: & 
      dG_dd_dX1_Temp &
         (1:nDOFX,7, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    LOGICAL :: Verbose

    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      dG_dd_dX1 = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeGeometryDerivatives_X1"
    END IF
 

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X1 = nK + [0,1,0,0]    ! Number of X1 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X1 = PRODUCT( nK_X1(2:4) ) ! Number of X1 Faces in Position Space
    
    ! --- Permute Geometry Fields ---


    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX





        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)
        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ3,iZ4,iZ2) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)
       
        G_munu_K(iNodeX,2,iZ3,iZ4,iZ2) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ3,iZ4,iZ2) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ3,iZ4,iZ2) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)
      END DO

    END DO
    END DO
    END DO
    ! --- Interpolate Geometry Fields ---






      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )
    ! --- Compute Metric Components ---

    ! --- Permute Fluid Fields ---


    ! --- Interpolate Fluid Fields ---

  
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1



        H_munu_F(iNodeX,1,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,1,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,2,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,2,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,3,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,3,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,4,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,4,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,5,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,5,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,6,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,6,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,7,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,7,iZ3,iZ4,iZ2)
     END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             H_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dG_dd_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             H_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dG_dd_dX1_Temp, nDOFX )

    ! --- Volume Term ---
    ! -------------------

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX



        H_munu_K(iNodeX,1,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,2,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,3,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,4,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,5,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,6,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,7,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ3,iZ4,iZ2)


      END DO

    END DO
    END DO
    END DO


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX1_Temp, nDOFX )

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX


        dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


        dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


      END DO

    END DO
    END DO
    END DO
    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        dG_dd_dX1(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX1(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP




      END DO

    END DO
    END DO
    END DO
  END SUBROUTINE ComputeGeometryDerivatives_X1

  SUBROUTINE ComputeGeometryDerivatives_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX2, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4
    INTEGER  :: nK(4), nK_X2(4), nX, nX_X2
    REAL(DP) :: &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3)), &
      G_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)), &
      H_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1)
    REAL(DP) ::  B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: & 
      dG_dd_dX2_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    LOGICAL :: Verbose

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      dG_dd_dX2 = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeGeometryDerivatives_X2"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X2 = nK + [0,0,1,0]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X2 = PRODUCT( nK_X2(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ2,iZ4,iZ3) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)

        G_munu_K(iNodeX,2,iZ2,iZ4,iZ3) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ2,iZ4,iZ3) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ2,iZ4,iZ3) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---











      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, 7*nX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, 7*nX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )



    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---




    ! --- Interpolate Fluid Fields ---


    ! --- Interpolate Left State ---



    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2



        H_munu_F(iNodeX,1,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,1,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,2,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,2,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,3,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,3,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,4,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,4,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,5,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,5,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,6,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,6,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,7,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,7,iZ2,iZ4,iZ3)


      END DO

    END DO
    END DO
    END DO








    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dG_dd_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dG_dd_dX2_Temp, nDOFX )










    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        H_munu_K(iNodeX,1,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,2,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,3,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,4,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,5,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,6,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,7,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ2,iZ4,iZ3)


      END DO

    END DO
    END DO
    END DO










    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX2_Temp, nDOFX )




    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX








        dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )


        dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )




      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX




        dG_dd_dX2(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX2(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGeometryDerivatives_X2


  SUBROUTINE ComputeGeometryDerivatives_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX3, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: iZ2, iZ3, iZ4
    INTEGER  :: nK(4), nK_X3(4), nX, nX_X3
    REAL(DP) :: B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4)), &
      G_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)), &
      H_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1)
    REAL(DP) :: & 
      dG_dd_dX3_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    LOGICAL :: Verbose

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      dG_dd_dX3 = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeGeometryDerivatives_X3"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X3 = nK + [0,0,0,1]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X3 = PRODUCT( nK_X3(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ2,iZ3,iZ4) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)

        G_munu_K(iNodeX,2,iZ2,iZ3,iZ4) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ2,iZ3,iZ4) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ2,iZ3,iZ4) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)


      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---









      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, 7*nX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, 7*nX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )



    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---



    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---


    ! --- Interpolate Right State ---


    ! --- Compute Face Velocity Components ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        ! --- Left States ---




        H_munu_F(iNodeX,1,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,1,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,2,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,2,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,3,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,3,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,4,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,4,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,5,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,5,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,6,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,6,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,7,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,7,iZ2,iZ3,iZ4)


      END DO

    END DO
    END DO
    END DO




    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dG_dd_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dG_dd_dX3_Temp, nDOFX )


    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX




        H_munu_K(iNodeX,1,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,2,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,3,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,4,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,5,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,6,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,7,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ2,iZ3,iZ4)




      END DO

    END DO
    END DO
    END DO






    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX3_Temp, nDOFX )



    ASSOCIATE( dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX





        dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


        dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )




      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        dG_dd_dX3(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX3(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

      END DO

    END DO
    END DO
    END DO


  END SUBROUTINE ComputeGeometryDerivatives_X3




  SUBROUTINE ComputeAlpha( V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                           U_u, dU_d_dX0, dU_d_dX1, dU_d_dX2, dU_d_dX3, &
                           Alp, B_u_1, B_u_2, B_u_3, C, Alpha )

    REAL(DP), INTENT(in)  :: V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: U_u(4), dU_d_dX0(4), dU_d_dX1(4), dU_d_dX2(4), dU_d_dX3(4)
    REAL(DP), INTENT(in)  :: Alp, B_u_1, B_u_2, B_u_3, C
    REAL(DP), INTENT(out) :: Alpha


    REAL(DP) :: dU_d_dX(4,4), G_uu_munu(4,4),  W, Vsq, Alpha_Eig, Alpha_A, A(4,4), Lambda(4), WORK(11)
    INTEGER  :: mu, nu, rho, INFO  
 
    Vsq =  V_u_1**2 * Gm_dd_11 + V_u_2**2 * Gm_dd_22 + V_u_3**2 * Gm_dd_33 
 
    W = 1.0_DP / SQRT( 1.0_DP - Vsq )

    dU_d_dX(1,1) = dU_d_dX0(1)
    dU_d_dX(1,2) = dU_d_dX0(2)
    dU_d_dX(1,3) = dU_d_dX0(3)
    dU_d_dX(1,4) = dU_d_dX0(4)

    dU_d_dX(2,1) = dU_d_dX1(1)
    dU_d_dX(2,2) = dU_d_dX1(2)
    dU_d_dX(2,3) = dU_d_dX1(3)
    dU_d_dX(2,4) = dU_d_dX1(4)

    dU_d_dX(3,1) = dU_d_dX2(1)
    dU_d_dX(3,2) = dU_d_dX2(2)
    dU_d_dX(3,3) = dU_d_dX2(3)
    dU_d_dX(3,4) = dU_d_dX2(4)

    dU_d_dX(4,1) = dU_d_dX3(1)
    dU_d_dX(4,2) = dU_d_dX3(2)
    dU_d_dX(4,3) = dU_d_dX3(3)
    dU_d_dX(4,4) = dU_d_dX3(4)

    G_uu_munu = 0.0_DP 
    
    G_uu_munu(1,1) = -1.0_DP / Alp**2
    G_uu_munu(2,1) = B_u_1 / Alp**2
    G_uu_munu(3,1) = B_u_2 / Alp**2
    G_uu_munu(4,1) = B_u_3 / Alp**2
    G_uu_munu(1,2) = B_u_1 / Alp**2
    G_uu_munu(1,3) = B_u_2 / Alp**2
    G_uu_munu(1,4) = B_u_3 / Alp**2
    G_uu_munu(2,2) = 1.0_DP / Gm_dd_11 - B_u_1 * B_u_1 / Alp**2
    G_uu_munu(3,3) = 1.0_DP / Gm_dd_22 - B_u_2 * B_u_2 / Alp**2
    G_uu_munu(4,4) = 1.0_DP / Gm_dd_33 - B_u_3 * B_u_3 / Alp**2

    A(1,1) =  dU_d_dX(1,1) + dU_d_dX(1,1)
    A(1,2) =  dU_d_dX(1,2) + dU_d_dX(2,1)
    A(1,3) =  dU_d_dX(1,3) + dU_d_dX(3,1)
    A(1,4) =  dU_d_dX(1,4) + dU_d_dX(4,1)

    A(2,1) =  dU_d_dX(2,1) + dU_d_dX(1,2)
    A(2,2) =  dU_d_dX(2,2) + dU_d_dX(2,2)
    A(2,3) =  dU_d_dX(2,3) + dU_d_dX(3,2)
    A(2,4) =  dU_d_dX(2,4) + dU_d_dX(4,2)

    A(3,1) =  dU_d_dX(3,1) + dU_d_dX(1,3)
    A(3,2) =  dU_d_dX(3,2) + dU_d_dX(2,3)
    A(3,3) =  dU_d_dX(3,3) + dU_d_dX(3,3)
    A(3,4) =  dU_d_dX(3,4) + dU_d_dX(4,3)

    A(4,1) =  dU_d_dX(4,1) + dU_d_dX(1,4)
    A(4,2) =  dU_d_dX(4,2) + dU_d_dX(2,4)
    A(4,3) =  dU_d_dX(4,3) + dU_d_dX(3,4)
    A(4,4) =  dU_d_dX(4,4) + dU_d_dX(4,4)

    A = 0.5_DP * A
  
    CALL DSYEV( 'N', 'U', 4, A, 4, Lambda, WORK, 11, INFO )

    Alpha_Eig = MAXVAL( ABS( Lambda ) ) / W

    Alpha_A = 0.0_DP

    DO rho = 1, 4
    DO mu = 1, 4 
    DO nu = 1, 4

      Alpha_A = Alpha_A + U_u(nu) * dU_d_dX(nu,mu) * U_u(nu) * dU_d_dX(nu,rho) * G_uu_munu(rho,mu)

    END DO
    END DO
    END DO

    Alpha_A = SQRT(ABS(Alpha_A)) / W

    Alpha = C * ( Alpha_Eig + Alpha_A )




  END SUBROUTINE ComputeAlpha

  SUBROUTINE ComputeChristoffel( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                                 dG_dd_dX0, dG_dd_dX1, dG_dd_dX2, dG_dd_dX3, Gamma_udd )
   

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in) :: & 
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 

    REAL(DP), INTENT(out) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))



    REAL(DP) :: &
      G_uu_munu (1:nDOFX,0:3,0:3, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4) )
     REAL(DP) :: dG_munu_dXrho &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
     INTEGER :: iZ2, iZ3, iZ4, iNodeX, mu, nu, rho, sig, iNodeX1
     REAL(DP) :: G_dd_11, G_dd_22, G_dd_33, B(3), A, X1, e, f

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeX1 )


        dG_munu_dXrho(iNodeX,0,:,:,iZ2,iZ3,iZ4) = dG_dd_dX0(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,1,:,:,iZ2,iZ3,iZ4) = dG_dd_dX1(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,2,:,:,iZ2,iZ3,iZ4) = dG_dd_dX2(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,3,:,:,iZ2,iZ3,iZ4) = dG_dd_dX3(iNodeX,:,:,iZ2,iZ3,iZ4)



        G_dd_11 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
        G_dd_22 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
        G_dd_33 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)
 

        A = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)
        B(1) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B(2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B(3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        G_uu_munu(iNodeX,:,:,iZ2,iZ3,iZ4) = 0.0_DP
          
        G_uu_munu(iNodeX,0,0,iZ2,iZ3,iZ4) = -1.0_DP / A**2
        G_uu_munu(iNodeX,1,0,iZ2,iZ3,iZ4) = B(1) / A**2
        G_uu_munu(iNodeX,2,0,iZ2,iZ3,iZ4) = B(2) / A**2
        G_uu_munu(iNodeX,3,0,iZ2,iZ3,iZ4) = B(3) / A**2
        G_uu_munu(iNodeX,0,1,iZ2,iZ3,iZ4) = B(1) / A**2
        G_uu_munu(iNodeX,0,2,iZ2,iZ3,iZ4) = B(2) / A**2
        G_uu_munu(iNodeX,0,3,iZ2,iZ3,iZ4) = B(3) / A**2
        G_uu_munu(iNodeX,1,1,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_11 - B(1) * B(1) / A**2
        G_uu_munu(iNodeX,2,2,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_22 - B(2) * B(2) / A**2
        G_uu_munu(iNodeX,3,3,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_33 - B(3) * B(3) / A**2
      END DO

    END DO
    END DO
    END DO
    Gamma_udd = 0.0_DP

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        DO rho = 0,3
        DO mu = 0,3
        DO nu = 0,3
        DO sig = 0,3

        Gamma_udd(iNodeX,rho,mu,nu,iZ2,iZ3,iZ4) = Gamma_udd(iNodeX,rho,mu,nu,iZ2,iZ3,iZ4) &
                                                + 0.5_DP * G_uu_munu(iNodeX,sig,rho,iZ2,iZ3,iZ4) &
                                                * ( dG_munu_dXrho(iNodeX,nu,mu,sig,iZ2,iZ3,iZ4) &
                                                + dG_munu_dXrho(iNodeX,mu,nu,sig,iZ2,iZ3,iZ4) &
                                                - dG_munu_dXrho(iNodeX,sig,mu,nu,iZ2,iZ3,iZ4) )


        END DO
        END DO
        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeChristoffel


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
    !$OMP MAP( to: dZ1, dZ2, dZ3, dZ4, &
    !$OMP          nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3 ) &
    !$OMP MAP( alloc: uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ1, dZ2, dZ3, dZ4, &
    !$ACC          nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3 ) &
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
    !$OMP MAP( release: dZ1, dZ2, dZ3, dZ4, &
    !$OMP               nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3, &
    !$OMP               uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ1, dZ2, dZ3, dZ4, &
    !$ACC         nZ, nZ_E, nZ_X1, nZ_X2, nZ_X3,  &
    !$ACC         uV1_K, uV2_K, uV3_K, uD_K, uI1_K, uI2_K, uI3_K )
#endif

    END ASSOCIATE ! dZ1, etc.

    DEALLOCATE( uV1_K, uV2_K, uV3_K )
    DEALLOCATE( uD_K, uI1_K, uI2_K, uI3_K )

  END SUBROUTINE FinalizeIncrement_TwoMoment_Explicit


  SUBROUTINE InitializeIncrement_Divergence_X &
    ( iZP_B0, iZP_E0, nDOFX_X, nDOFZ_X, &
      GX_K, GX_F, uCF_K, uCF_L, uCF_R, uCR_K, uCR_L, uCR_R )

    INTEGER, INTENT(in) :: iZP_B0(4), iZP_E0(4) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X ! nDOFX_X1, ...
    INTEGER, INTENT(in) :: nDOFZ_X ! nDOFZ_X1, ...

    ! --- Geometry Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      GX_K (1:nDOFX, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nGF), &
      GX_F (1:nDOFX_X, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCF_K(1:nDOFX, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nCF), &
      uCF_L(1:nDOFX_X, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCF), &
      uCF_R(1:nDOFX_X, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCF)

    ! --- Conserved Radiation Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCR_K(1:nDOFZ, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nCR), &
      uCR_L(1:nDOFZ_X, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCR), &
      uCR_R(1:nDOFZ_X, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCR)

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

    Gm_dd_11_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_SqrtGm  )
    G_Alpha_K (1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Alpha   )
    G_Beta_1_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_1  )
    G_Beta_2_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_2  )
    G_Beta_3_K(1:nNodesX_K) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_3  )


    Gm_dd_11_F(1:nNodesX_X) => GX_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => GX_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => GX_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => GX_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_SqrtGm  )
    G_Alpha_F (1:nNodesX_X) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Alpha   )
    G_Beta_1_F(1:nNodesX_X) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_1  )
    G_Beta_2_F(1:nNodesX_X) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_2  )
    G_Beta_3_F(1:nNodesX_X) => GX_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Beta_3  )


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

    ALLOCATE( uD_L (nNodesZ_X) )
    ALLOCATE( uI1_L(nNodesZ_X) )
    ALLOCATE( uI2_L(nNodesZ_X) )
    ALLOCATE( uI3_L(nNodesZ_X) )

    ALLOCATE( uD_R (nNodesZ_X) )
    ALLOCATE( uI1_R(nNodesZ_X) )
    ALLOCATE( uI2_R(nNodesZ_X) )
    ALLOCATE( uI3_R(nNodesZ_X) )

    ALLOCATE( uV1_F(nNodesX_X) )
    ALLOCATE( uV2_F(nNodesX_X) )
    ALLOCATE( uV3_F(nNodesX_X) )

    ALLOCATE( nIterations_L(nNodesZ_X) )
    ALLOCATE( nIterations_R(nNodesZ_X) )
    ALLOCATE( nIterations_K(nNodesZ_K) )

    ALLOCATE( PositionIndexZ_F(nNodesZ_X) )
    ALLOCATE( PositionIndexZ_K(nNodesZ_K) )

    ALLOCATE( IndexTableZ_F(7,nNodesZ_X) )
    ALLOCATE( IndexTableZ_K(7,nNodesZ_K) )

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

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP          IndexTableZ_F, IndexTableZ_K ) &
    !$OMP MAP( alloc: uV1_F, uV2_F, uV3_F, &
    !$OMP             uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP             uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP             nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K ) &
    !$ACC CREATE( uV1_F, uV2_F, uV3_F, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif

    RETURN
  END SUBROUTINE InitializeIncrement_Divergence_X


  SUBROUTINE FinalizeIncrement_Divergence_X

#if defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: PositionIndexZ_F, PositionIndexZ_K, &
    !$OMP               IndexTableZ_F, IndexTableZ_K, &
    !$OMP               uV1_F, uV2_F, uV3_F, &
    !$OMP               uD_L, uI1_L, uI2_L, uI3_L, &
    !$OMP               uD_R, uI1_R, uI2_R, uI3_R, &
    !$OMP               nIterations_L, nIterations_R, nIterations_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( PositionIndexZ_F, PositionIndexZ_K, &
    !$ACC         IndexTableZ_F, IndexTableZ_K, &
    !$ACC         uV1_F, uV2_F, uV3_F, &
    !$ACC         uD_L, uI1_L, uI2_L, uI3_L, &
    !$ACC         uD_R, uI1_R, uI2_R, uI3_R, &
    !$ACC         nIterations_L, nIterations_R, nIterations_K )
#endif

    DEALLOCATE( PositionIndexZ_F, PositionIndexZ_K )
    DEALLOCATE( IndexTableZ_F, IndexTableZ_K )
    DEALLOCATE( uD_L, uI1_L, uI2_L, uI3_L )
    DEALLOCATE( uD_R, uI1_R, uI2_R, uI3_R )
    DEALLOCATE( uV1_F, uV2_F, uV3_F )
    DEALLOCATE( nIterations_L, nIterations_R, nIterations_K )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( G_Alpha_K, G_Beta_1_K, G_Beta_2_K, G_Beta_3_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY( G_Alpha_F, G_Beta_1_F, G_Beta_2_F, G_Beta_3_F )

    NULLIFY( uN_K, uG1_K, uG2_K, uG3_K )
    NULLIFY( uN_L, uG1_L, uG2_L, uG3_L )
    NULLIFY( uN_R, uG1_R, uG2_R, uG3_R )

  END SUBROUTINE FinalizeIncrement_Divergence_X




END MODULE TwoMoment_DiscretizationModule_Streaming_Relativistic
