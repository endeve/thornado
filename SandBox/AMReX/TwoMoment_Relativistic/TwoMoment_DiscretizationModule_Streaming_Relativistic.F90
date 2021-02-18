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
    ComputeConserved_TwoMoment, &
    Flux_X1, &
    Flux_E, &
    Source_E, &
    ComputeEddingtonTensorComponents_ud, &
    ComputeHeatFluxTensorComponents_uud_Lagrangian, &
    NumericalFlux_LLF
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: UpwindFlux = .FALSE.


  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

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
    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS, i
    LOGICAL :: Verbose
    LOGICAL  :: SuppressBC
    REAL(DP):: PD, V1, V2, V3, PE, PNe, W

    REAL(DP) :: &
      dWV_u_dX1 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dWV_d_dX1 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dW_dX1 &
        (nDOFX, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))

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


    CALL ApplyBoundaryConditions_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )


!    CALL CheckRealizability(U_R, U_F, GX, iZ_B1, iZ_E1, iZ_B0, iZ_E0)

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


    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, &
             Verbose_Option = Verbose )

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
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iPF
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

!this is where the issue is

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

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX, INFO
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nK(4), nK_Z1(4), nV, nV_Z1, i, j
    REAL(DP) :: EdgeEnergyCubed
    REAL(DP) :: Alpha, AlphaL, AlphaR
    REAL(DP) ::        AlphaM, AlphaP
    REAL(DP) :: A(4,4), Lambda(4), WORK(11)
    REAL(DP) :: k_ud_11_L, k_ud_12_L, k_ud_13_L
    REAL(DP) ::            k_ud_22_L, k_ud_23_L
    REAL(DP) ::                       k_ud_33_L
    REAL(DP) :: k_ud_11_R, k_ud_12_R, k_ud_13_R
    REAL(DP) ::            k_ud_22_R, k_ud_23_R
    REAL(DP) ::                       k_ud_33_R
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
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
      dU_u_dX0 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_u_dX1 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_u_dX2 &
        (nDOFX,0:3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_u_dX3 &
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
    REAL(DP):: W, E2, C, l_uud_munurho(0:3,0:3,0:3), k_ud_munu(0:3,0:3)

    LOGICAL :: Verbose
    CHARACTER(len=40) :: name1, name2
    CHARACTER(len=1):: nds
    CHARACTER(len=2)::nxn1
    CHARACTER(len=3)::nxn2

!    IF ( nDOFX == 1) THEN
!      nds="1"
!    ELSE IF( nDOFX == 2) THEN
!      nds="2"
!    ELSE
!      nds="3"
!    END IF 
!
!    IF ( nX(1) == 32) THEN
!      nxn1="32"
!    ELSE IF( nX(1) == 64) THEN
!      nxn1="64"
!    ELSE IF (nX(1) == 128) THEN
!      nxn2="128"
!    ELSE
!      nxn2="256"
!    END IF 
!
!    
!
!    print*, name1
!    IF (nX(1)==32 .OR. nX(1)==64) THEN
!      name1='dU0'//nds//nxn1//'.txt'
!      name2='dU1'//nds//nxn1//'.txt'
!    ELSE
!      name1='dU0'//nds//nxn2//'.txt'
!      name1='dU1'//nds//nxn2//'.txt'
!    END IF
!    name1=trim(name1)
!    name2=trim(name2)
!
!
    IF( iZ_E0(1) .EQ. iZ_B0(1) ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

    IF (Verbose) THEN
       PRINT*, "      ComputeIncrement_ObserverCorrections"
    END IF

    CALL ComputeWeakDerivatives_X0 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dU_d_dX0, dU_u_dX0, dG_dd_dX0, &
             Verbose_Option = Verbose )

    CALL ComputeWeakDerivatives_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dU_d_dX1, dU_u_dX1, dG_dd_dX1, &
             Verbose_Option = Verbose )

    CALL ComputeWeakDerivatives_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dU_d_dX2, dU_u_dX2, dG_dd_dX2, &
             Verbose_Option = Verbose )

    CALL ComputeWeakDerivatives_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dU_d_dX3, dU_u_dX3, dG_dd_dX3, &
             Verbose_Option = Verbose )

    CALL ComputeFourVelocity &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, U_u, U_d, uPF_K, &
             Verbose_Option = Verbose  )


    !open(2, file = name1, status = 'new') 
    !open(3, file = name2, status = 'new') 

    !DO iZ4 = iZ_B0(4), iZ_E0(4)
    !DO iZ3 = iZ_B0(3), iZ_E0(3)
    !DO iZ2 = iZ_B0(2), iZ_E0(2)

      !DO iNodeX = 1, nDOFX
         !write(2,*) dU_dX1(iNodeX,0,iZ2,iZ3,iZ4)
         !write(3,*) dU_dX1(iNodeX,1,iZ2,iZ3,iZ4)
      !END DO

    !END DO
    !END DO
    !END DO

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
                    dU_d_dX0(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3(iNode,:,iZ2,iZ3,iZ4) )
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
                    dU_d_dX0(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2(iNode,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3(iNode,:,iZ2,iZ3,iZ4) )

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
                             dU_u_dX0(iNode,:,iZ2,iZ3,iZ4), &
                             dU_u_dX1(iNode,:,iZ2,iZ3,iZ4), &
                             dU_u_dX2(iNode,:,iZ2,iZ3,iZ4), &
                             dU_u_dX3(iNode,:,iZ2,iZ3,iZ4), &
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
                    dU_d_dX0(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX1(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX2(iNodeX,:,iZ2,iZ3,iZ4), &
                    dU_d_dX3(iNodeX,:,iZ2,iZ3,iZ4) )

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
                    dU_d_dX3(iNodeX,:,iZ2,iZ3,iZ4) )

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

print*, "Fuck"
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
    INTEGER  :: i, iZ2, iZ3, iZ4, iCF, iGF
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
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dU_d_dX0, dU_u_dX0, dG_dd_dX0, Verbose_Option  )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_d_dX0 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dU_u_dX0 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    dU_d_dX0 = Zero
    dU_u_dX0 = Zero
    dG_dd_dX0 = Zero
  END SUBROUTINE ComputeWeakDerivatives_X0 


  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dU_d_dX1, dU_u_dX1, dG_dd_dX1, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_d_dX1 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dU_u_dX1 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP), INTENT(out) :: &
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: nK(4), nK_X1(4), nX, nX_X1
    INTEGER  :: iNodeX, iNodeX1
    INTEGER  :: i, iZ2, iZ3, iZ4, iCF, iGF
    REAL(DP) :: &
      uPF_K(nPF), uPF_L(nPF), uPF_R(nPF), X1, V0, pi3, pi23, sine, cose
    REAL(DP) :: &
      V_u_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      V_d_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      V_u_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      V_d_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      WV_u_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      WV_d_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      WV_u_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      WV_d_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      W_X1(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      W_K(nDOFX,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      dW_dX1(nDOFX,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
               iZ_B0(2):iZ_E0(2)), &
      dWV_u_dX1(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
               iZ_B0(2):iZ_E0(2)), &
      dWV_d_dX1(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
               iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      GX_K(nDOFX   ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)  ,nGF), &
      GX_F(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1,nGF), &
      GX_K_New(nDOFX   ,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)), &
      GX_F_New(nDOFX_X1,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1), &
      G_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)), &
      G_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)), &
      H_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1)
    REAL(DP) :: &
      uCF_K(nDOFX   ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B1(2):iZ_E1(2)  ,nCF), &
      uCF_L(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF), &
      uCF_R(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF)
    REAL(DP) :: Vsq_X1, Vsq_K, Vsq_X1_New, Vsq_K_New
    REAL(DP) :: V_u_K_New(3), V_d_K_New(3), V_u_X1_New(3), V_d_X1_New(3), W_X1_New, W_K_New
    REAL(DP) :: B_u_X1(3), B_d_X1(3), A_X1, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      U_u_X1(nDOFX_X1,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      U_d_X1(nDOFX_X1,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      U_u_K(nDOFX,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      U_d_K(nDOFX,4,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      dU_d_dX1_Temp &
         (1:nDOFX,4, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    REAL(DP) :: & 
      dU_u_dX1_Temp &
         (1:nDOFX,4, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    REAL(DP) :: & 
      dG_dd_dX1_Temp &
         (1:nDOFX,7, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    CHARACTER(len=40) :: name1, name2
    CHARACTER(len=1):: nds
    CHARACTER(len=2)::nxn1
    CHARACTER(len=3)::nxn2
    LOGICAL :: Verbose

    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      dU_d_dX1 = Zero
      dU_u_dX1 = Zero
      dG_dd_dX1 = Zero
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
    
!    IF ( nDOFX == 1) THEN
!      nds="1"
!    ELSE IF( nDOFX == 2) THEN
!      nds="2"
!    ELSE
!      nds="3"
!    END IF 
!
!    IF ( nX == 32) THEN
!      nxn1="32"
!    ELSE IF( nX == 64) THEN
!      nxn1="64"
!    ELSE IF (nX == 128) THEN
!      nxn2="128"
!    ELSE
!      nxn2="256"
!    END IF 
!    
!    IF (nX==32 .OR. nX==64) THEN
!      name1="dWdX"//nds//nxn1//".txt"
!      name2="dWvdX"//nds//nxn1//".txt"
!    ELSE 
!      name1='dWdX'//nds//nxn2//'.txt'
!      name2="dWvdX"//nds//nxn2//".txt"
!    END IF
!
!    print*, name1
!    print*, name2
!    name1=trim(name1)
!    name2=trim(name2)
    ! --- Permute Geometry Fields ---

    DO iGF = 1, nGF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

        GX_K_New(iNodeX,iGF,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)




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
    END DO

    ! --- Interpolate Geometry Fields ---

    DO iGF = iGF_h_1, iGF_h_3

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nGF*nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K_New(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
               GX_F_New(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nGF*nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K_New(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX, Half, &
               GX_F_New(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )






      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )
    ! --- Compute Metric Components ---

    DO i   = 1, 3
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11+i-1) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1+i-1)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO
    END DO

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
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

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
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )


        V_u_X1(iNodeX,1:3,iZ3,iZ4,iZ2) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

        V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) 

        V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) 

        V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) 

        V_d_X1(iNodeX,1,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_X1(iNodeX,1,iZ3,iZ4,iZ2)

        V_d_X1(iNodeX,2,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_X1(iNodeX,2,iZ3,iZ4,iZ2)

        V_d_X1(iNodeX,3,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_X1(iNodeX,3,iZ3,iZ4,iZ2)

        Vsq_X1 = V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) * V_d_X1(iNodeX,1,iZ3,iZ4,iZ2) &
               + V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) * V_d_X1(iNodeX,2,iZ3,iZ4,iZ2) &
               + V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) * V_d_X1(iNodeX,3,iZ3,iZ4,iZ2) 
 
        W_X1(iNodeX,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX)  / SQRT(1.0_DP - Vsq_X1) 

        WV_u_X1(iNodeX,1,iZ3,iZ4,iZ2) = V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2) 

        WV_u_X1(iNodeX,2,iZ3,iZ4,iZ2) = V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2) 

        WV_u_X1(iNodeX,3,iZ3,iZ4,iZ2) = V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2)
 
        WV_d_X1(iNodeX,1,iZ3,iZ4,iZ2) = V_d_X1(iNodeX,1,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2) 

        WV_d_X1(iNodeX,2,iZ3,iZ4,iZ2) = V_d_X1(iNodeX,2,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2) 

        WV_d_X1(iNodeX,3,iZ3,iZ4,iZ2) = V_d_X1(iNodeX,3,iZ3,iZ4,iZ2) * W_X1(iNodeX,iZ3,iZ4,iZ2) 









        V_u_X1_New(1:3) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )


        V_d_X1_New(1) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_X1_New(1)

        V_d_X1_New(2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_X1_New(2)

        V_d_X1_New(3) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_X1_New(3)

        Vsq_X1_New = V_u_X1_New(1) * V_d_X1_New(1) &  
                    + V_u_X1_New(2) * V_d_X1_New(2) &  
                    + V_u_X1_New(3) * V_d_X1_New(3) 
    
        W_X1_New = WeightsX_X1(iNodeX)  / SQRT(1.0_DP - Vsq_X1_New) 

        A_X1 = GX_F_New(iNodeX,iGF_Alpha,iZ3,iZ4,iZ2)

        B_u_X1(1) = GX_F_New(iNodeX,iGF_Beta_1,iZ3,iZ4,iZ2)
        B_u_X1(2) = GX_F_New(iNodeX,iGF_Beta_2,iZ3,iZ4,iZ2)
        B_u_X1(3) = GX_F_New(iNodeX,iGF_Beta_3,iZ3,iZ4,iZ2)

        B_d_X1(1) = GX_F_New(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) * B_u_X1(1)
        B_d_X1(2) = GX_F_New(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) * B_u_X1(2)
        B_d_X1(3) = GX_F_New(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) * B_u_X1(3)


        U_u_X1(iNodeX,1,iZ3,iZ4,iZ2) = W_X1_New / A_X1
        
        U_u_X1(iNodeX,2,iZ3,iZ4,iZ2) = W_X1_New * ( V_u_X1_New(1) - B_u_X1(1) / A_X1)

        U_u_X1(iNodeX,3,iZ3,iZ4,iZ2) = W_X1_New * ( V_u_X1_New(2) - B_u_X1(2) / A_X1)

        U_u_X1(iNodeX,4,iZ3,iZ4,iZ2) = W_X1_New * ( V_u_X1_New(3) - B_u_X1(3) / A_X1)



        U_d_X1(iNodeX,1,iZ3,iZ4,iZ2) = W_X1_New * ( - A_X1 + B_d_X1(1) * V_u_X1_New(1) &
                                     + B_d_X1(2) * V_u_X1_New(2) + B_d_X1(3) * V_u_X1_New(3) )
        
        U_d_X1(iNodeX,2,iZ3,iZ4,iZ2) = W_X1_New * V_d_X1_New(1) 

        U_d_X1(iNodeX,3,iZ3,iZ4,iZ2) = W_X1_New * V_d_X1_New(2) 

        U_d_X1(iNodeX,4,iZ3,iZ4,iZ2) = W_X1_New * V_d_X1_New(3)




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
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             WV_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dWV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             WV_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dWV_d_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             WV_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dWV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             WV_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dWV_d_dX1, nDOFX )

    ! -------------------
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             W_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dW_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             W_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dW_dX1, nDOFX )









    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             U_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dU_u_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             U_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dU_d_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             U_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dU_u_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             U_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dU_d_dX1_Temp, nDOFX )




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
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

        V_u_K(iNodeX,1,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V1) 

        V_u_K(iNodeX,2,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V2)  

        V_u_K(iNodeX,3,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V3) 

        V_d_K(iNodeX,1,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_K(iNodeX,1,iZ3,iZ4,iZ2)

        V_d_K(iNodeX,2,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_K(iNodeX,2,iZ3,iZ4,iZ2)

        V_d_K(iNodeX,3,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_K(iNodeX,3,iZ3,iZ4,iZ2)

        Vsq_K =  V_u_K(iNodeX,1,iZ3,iZ4,iZ2) * V_d_K(iNodeX,1,iZ3,iZ4,iZ2) & 
              +  V_u_K(iNodeX,2,iZ3,iZ4,iZ2) * V_d_K(iNodeX,2,iZ3,iZ4,iZ2) &
              +  V_u_K(iNodeX,3,iZ3,iZ4,iZ2) * V_d_K(iNodeX,3,iZ3,iZ4,iZ2) 

        W_K(iNodeX,iZ3,iZ4,iZ2) = ( WeightsX_q(iNodeX)  / SQRT( 1.0_DP - Vsq_K ) )  

        WV_u_K(iNodeX,1,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_u_K(iNodeX,1,iZ3,iZ4,iZ2)

        WV_u_K(iNodeX,2,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_u_K(iNodeX,2,iZ3,iZ4,iZ2)

        WV_u_K(iNodeX,3,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_u_K(iNodeX,3,iZ3,iZ4,iZ2)

        WV_d_K(iNodeX,1,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_d_K(iNodeX,1,iZ3,iZ4,iZ2)

        WV_d_K(iNodeX,2,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_d_K(iNodeX,2,iZ3,iZ4,iZ2)

        WV_d_K(iNodeX,3,iZ3,iZ4,iZ2) = W_K(iNodeX,iZ3,iZ4,iZ2) * V_d_K(iNodeX,3,iZ3,iZ4,iZ2)





        V_u_K_New(1) = uPF_K(iPF_V1)
 
        V_u_K_New(2) = uPF_K(iPF_V2) 

        V_u_K_New(3) = uPF_K(iPF_V3) 

        V_d_K_New(1) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_K_New(1)

        V_d_K_New(2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_K_New(2)

        V_d_K_New(3) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_K_New(3)

        A_K = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Alpha)

        B_u_K(1) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Beta_1)
        B_u_K(2) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Beta_2)
        B_u_K(3) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Beta_3)

        B_d_K(1) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) * B_u_K(3)
        
        Vsq_K_New = V_u_K_New(1) * V_d_K_New(1) &
                  + V_u_K_New(2) * V_d_K_New(2) &
                  + V_u_K_New(3) * V_d_K_New(3)
     
        W_K_New = WeightsX_q(iNodeX) / SQRT( 1.0_DP - Vsq_K_New )

        U_u_K(iNodeX,1,iZ3,iZ4,iZ2) = W_K_New / A_K

        U_u_K(iNodeX,2,iZ3,iZ4,iZ2) = W_K_New * ( - B_u_K(1) / A_K + V_u_K_New(1) )

        U_u_K(iNodeX,3,iZ3,iZ4,iZ2) = W_K_New * ( - B_u_K(2) / A_K + V_u_K_New(2) )

        U_u_K(iNodeX,4,iZ3,iZ4,iZ2) = W_K_New * ( - B_u_K(3) / A_K + V_u_K_New(3) )
       

        U_d_K(iNodeX,1,iZ3,iZ4,iZ2) = W_K_New * ( -A_K + B_d_K(1) * V_u_K_New(1) &
                                              +    B_d_K(2) * V_u_K_New(2) &
                                              +    B_d_K(3) * V_u_K_New(3)  )

        U_d_K(iNodeX,2,iZ3,iZ4,iZ2) = W_K_New * V_d_K_New(1) 

        U_d_K(iNodeX,3,iZ3,iZ4,iZ2) = W_K_New * V_d_K_New(2) 

        U_d_K(iNodeX,4,iZ3,iZ4,iZ2) = W_K_New * V_d_K_New(3) 









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
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             WV_u_K, nDOFX, One, dWV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             WV_d_K, nDOFX, One, dWV_d_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             W_K, nDOFX, One, dW_dX1, nDOFX )




    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             U_u_K, nDOFX, One, dU_u_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX1_Temp, nDOFX )





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX1_Temp, nDOFX )

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        dWV_u_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
         = dWV_u_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dWV_d_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
         = dWV_d_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dWV_u_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
         = dWV_u_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dWV_d_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
         = dWV_d_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dWV_u_dX1(iNodeX,3,iZ3,iZ4,iZ2) &
         = dWV_u_dX1(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dWV_d_dX1(iNodeX,3,iZ3,iZ4,iZ2) &
         = dWV_d_dX1(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dW_dX1(iNodeX,iZ3,iZ4,iZ2) &
         = dW_dX1(iNodeX,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )






        dU_u_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
         = dU_u_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


        dU_u_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
         = dU_u_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dU_u_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
         = dU_u_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dU_u_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
         = dU_u_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


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

        dU_d_dX1(iNodeX,0,iZ2,iZ3,iZ4)  &
          = -GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) * dW_dX1(iNodeX,iZ3,iZ4,iZ2) &
            + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1) * dWV_u_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
            + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2) * dWV_u_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
            + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) * dWV_u_dX1(iNodeX,3,iZ3,iZ4,iZ2) 


        dU_d_dX1(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,1,iZ3,iZ4,iZ2)

        dU_d_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,2,iZ3,iZ4,iZ2)

        dU_d_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,3,iZ3,iZ4,iZ2)

        dU_u_dX1(iNodeX,0,iZ2,iZ3,iZ4)  &
          = 1.0_DP / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) * dW_dX1(iNodeX,iZ3,iZ4,iZ2)

        dU_u_dX1(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,1,iZ3,iZ4,iZ2) &
          - GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1) / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) * dW_dX1(iNodeX,iZ3,iZ4,iZ2)
 
        dU_u_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,2,iZ3,iZ4,iZ2) &
          - GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2) / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) * dW_dX1(iNodeX,iZ3,iZ4,iZ2)
 
        dU_u_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dWV_d_dX1(iNodeX,3,iZ3,iZ4,iZ2) &
          - GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha) * dW_dX1(iNodeX,iZ3,iZ4,iZ2)


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










!         print*,iZ2, dU_d_dX1(iNodeX,0,iZ2,iZ3,iZ4) - dU_d_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2), &
!                 dU_d_dX1(iNodeX,1,iZ2,iZ3,iZ4) - dU_d_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2), &
!                 dU_d_dX1(iNodeX,2,iZ2,iZ3,iZ4) - dU_d_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2), &
!                 dU_d_dX1(iNodeX,3,iZ2,iZ3,iZ4) - dU_d_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2)
!
!
 !        print*,iZ2, dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2), dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) 
      END DO

    END DO
    END DO
    END DO
!STOP
  END SUBROUTINE ComputeWeakDerivatives_X1

  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dU_d_dX2, dU_u_dX2, dG_dd_dX2, Verbose_Option  )

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
    REAL(DP), INTENT(out) :: &
      dU_d_dX2 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_u_dX2 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
    INTEGER  :: nK(4), nK_X2(4), nX, nX_X2
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3)), &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3)), &
      G_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)), &
      H_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1)
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
      U_u_X2(nDOFX_X2,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
             iZ_B0(3):iZ_E0(3)+1), &
      U_d_X2(nDOFX_X2,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
             iZ_B0(3):iZ_E0(3)+1), &
      U_u_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3)), &
      U_d_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      dU_d_dX2_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    REAL(DP) :: & 
      dU_u_dX2_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    REAL(DP) :: & 
      dG_dd_dX2_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    LOGICAL :: Verbose

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      dU_d_dX2 = Zero
      dU_u_dX2 = Zero
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


        U_u_X2(iNodeX,1,iZ2,iZ4,iZ3) = W_X2 / A_X2
        
        U_u_X2(iNodeX,2,iZ2,iZ4,iZ3) = W_X2 * ( V_u_X2(1) - B_u_X2(1) / A_X2)

        U_u_X2(iNodeX,3,iZ2,iZ4,iZ3) = W_X2 * ( V_u_X2(2) - B_u_X2(2) / A_X2)

        U_u_X2(iNodeX,4,iZ2,iZ4,iZ3) = W_X2 * ( V_u_X2(3) - B_u_X2(3) / A_X2)



        U_d_X2(iNodeX,1,iZ2,iZ4,iZ3) = W_X2 * ( - A_X2 + B_d_X2(1) * V_u_X2(1) &
                                     + B_d_X2(2) * V_u_X2(2) + B_d_X2(3) * V_u_X2(3) )
        
        U_d_X2(iNodeX,2,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(1) 

        U_d_X2(iNodeX,3,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(2) 

        U_d_X2(iNodeX,4,iZ2,iZ4,iZ3) = W_X2 * V_d_X2(3)


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
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             U_u_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dU_u_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             U_d_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dU_d_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             U_u_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dU_u_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             U_d_X2(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dU_d_dX2_Temp, nDOFX )




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


 
        W_K = WeightsX_q(iNodeX) / SQRT( 1.0_DP - Vsq_K )

        U_u_K(iNodeX,1,iZ2,iZ4,iZ3) = W_K / A_K

        U_u_K(iNodeX,2,iZ2,iZ4,iZ3) = W_K * ( - B_u_K(1) / A_K + V_u_K(1) )

        U_u_K(iNodeX,3,iZ2,iZ4,iZ3) = W_K * ( - B_u_K(2) / A_K + V_u_K(2) )

        U_u_K(iNodeX,4,iZ2,iZ4,iZ3) = W_K * ( - B_u_K(3) / A_K + V_u_K(3) )
       

        U_d_K(iNodeX,1,iZ2,iZ4,iZ3) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d_K(iNodeX,2,iZ2,iZ4,iZ3) = W_K * V_d_K(1) 

        U_d_K(iNodeX,3,iZ2,iZ4,iZ3) = W_K * V_d_K(2) 

        U_d_K(iNodeX,4,iZ2,iZ4,iZ3) = W_K * V_d_K(3)



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
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             U_u_K, nDOFX, One, dU_u_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX2_Temp, nDOFX )





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX2_Temp, nDOFX )




    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        dU_u_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
         = dU_u_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )


        dU_u_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
         = dU_u_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dU_u_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
         = dU_u_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dU_u_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
         = dU_u_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )


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


        dU_d_dX2(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3)

        dU_d_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3)


        dU_u_dX2(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_u_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3)

        dU_u_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_u_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3)

        dU_u_dX2(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_u_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3)

        dU_u_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_u_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3)


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

  END SUBROUTINE ComputeWeakDerivatives_X2


  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dU_d_dX3, dU_u_dX3, dG_dd_dX3, Verbose_Option  )

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
    REAL(DP), INTENT(out) :: &
      dU_d_dX3 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dU_u_dX3 &
         (1:nDOFX,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(out) :: &
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
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
           iZ_B0(4):iZ_E1(4)), &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4)), &
      G_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)), &
      H_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1)
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
      U_u_X3(nDOFX_X3,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
             iZ_B0(4):iZ_E0(4)+1), &
      U_d_X3(nDOFX_X3,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
             iZ_B0(4):iZ_E0(4)+1), &
      U_u_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4)), &
      U_d_K(nDOFX,4,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_d_dX3_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: & 
      dU_u_dX3_Temp &
         (1:nDOFX,4, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    REAL(DP) :: & 
      dG_dd_dX3_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    LOGICAL :: Verbose

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      dU_d_dX3 = Zero
      dU_u_dX3 = Zero
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


        U_u_X3(iNodeX,1,iZ2,iZ3,iZ4) = W_X3 / A_X3
        
        U_u_X3(iNodeX,2,iZ2,iZ3,iZ4) = W_X3 * ( V_u_X3(1) - B_u_X3(1) / A_X3)

        U_u_X3(iNodeX,3,iZ2,iZ3,iZ4) = W_X3 * ( V_u_X3(2) - B_u_X3(2) / A_X3)

        U_u_X3(iNodeX,4,iZ2,iZ3,iZ4) = W_X3 * ( V_u_X3(3) - B_u_X3(3) / A_X3)



        U_d_X3(iNodeX,1,iZ2,iZ3,iZ4) = W_X3 * ( - A_X3 + B_d_X3(1) * V_u_X3(1) &
                                     + B_d_X3(2) * V_u_X3(2) + B_d_X3(3) * V_u_X3(3) )
        
        U_d_X3(iNodeX,2,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(1) 

        U_d_X3(iNodeX,3,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(2) 

        U_d_X3(iNodeX,4,iZ2,iZ3,iZ4) = W_X3 * V_d_X3(3)




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
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             U_u_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dU_u_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             U_d_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dU_d_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             U_u_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dU_u_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             U_d_X3(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dU_d_dX3_Temp, nDOFX )




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
     
        W_K = WeightsX_q(iNodeX) / SQRT( 1.0_DP - Vsq_K )

        U_u_K(iNodeX,1,iZ2,iZ3,iZ4) = W_K / A_K

        U_u_K(iNodeX,2,iZ2,iZ3,iZ4) = W_K * ( - B_u_K(1) / A_K + V_u_K(1) )

        U_u_K(iNodeX,3,iZ2,iZ3,iZ4) = W_K * ( - B_u_K(2) / A_K + V_u_K(2) )

        U_u_K(iNodeX,4,iZ2,iZ3,iZ4) = W_K * ( - B_u_K(3) / A_K + V_u_K(3) )
       

        U_d_K(iNodeX,1,iZ2,iZ3,iZ4) = W_K * ( -A_K + B_d_K(1) * V_u_K(1) &
                                              +    B_d_K(2) * V_u_K(2) &
                                              +    B_d_K(3) * V_u_K(3)  )

        U_d_K(iNodeX,2,iZ2,iZ3,iZ4) = W_K * V_d_K(1) 

        U_d_K(iNodeX,3,iZ2,iZ3,iZ4) = W_K * V_d_K(2) 

        U_d_K(iNodeX,4,iZ2,iZ3,iZ4) = W_K * V_d_K(3) 









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
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             U_u_K, nDOFX, One, dU_u_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             U_d_K, nDOFX, One, dU_d_dX3_Temp, nDOFX )





    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX3_Temp, nDOFX )



    ASSOCIATE( dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        dU_u_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
         = dU_u_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


        dU_u_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
         = dU_u_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dU_u_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
         = dU_u_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dU_u_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
         = dU_u_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


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

        dU_d_dX3(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4)

        dU_d_dX3(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_d_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4)


        dU_u_dX3(iNodeX,0,iZ2,iZ3,iZ4)  &
          = dU_u_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4)

        dU_u_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
          = dU_u_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4)

        dU_u_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
          = dU_u_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4)

        dU_u_dX3(iNodeX,3,iZ2,iZ3,iZ4)  &
          = dU_u_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4)

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


  END SUBROUTINE ComputeWeakDerivatives_X3

  SUBROUTINE ComputeAlpha( V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                           U_u, dU_d_dX0, dU_d_dX1, dU_d_dX2, dU_d_dX3, &
                           dU_u_dX0, dU_u_dX1, dU_u_dX2, dU_u_dX3, C, Alpha )

    REAL(DP), INTENT(in)  :: V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: U_u(4), dU_d_dX0(4), dU_d_dX1(4), dU_d_dX2(4), dU_d_dX3(4)
    REAL(DP), INTENT(in)  :: dU_u_dX0(4), dU_u_dX1(4), dU_u_dX2(4), dU_u_dX3(4), C     
    REAL(DP), INTENT(out) :: Alpha


    REAL(DP) :: dU_d_dX(4,4), dU_u_dX(4,4), W, Vsq, Alpha_Eig, Alpha_A, A(4,4), Lambda(4), WORK(11), l_mu(4)
    INTEGER  :: mu, nu, INFO, B  
 
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

    dU_u_dX(1,1) = dU_u_dX0(1)
    dU_u_dX(1,2) = dU_u_dX0(2)
    dU_u_dX(1,3) = dU_u_dX0(3)
    dU_u_dX(1,4) = dU_u_dX0(4)

    dU_u_dX(2,1) = dU_u_dX1(1)
    dU_u_dX(2,2) = dU_u_dX1(2)
    dU_u_dX(2,3) = dU_u_dX1(3)
    dU_u_dX(2,4) = dU_u_dX1(4)

    dU_u_dX(3,1) = dU_u_dX2(1)
    dU_u_dX(3,2) = dU_u_dX2(2)
    dU_u_dX(3,3) = dU_u_dX2(3)
    dU_u_dX(3,4) = dU_u_dX2(4)

    dU_u_dX(4,1) = dU_u_dX3(1)
    dU_u_dX(4,2) = dU_u_dX3(2)
    dU_u_dX(4,3) = dU_u_dX3(3)
    dU_u_dX(4,4) = dU_u_dX3(4)

    B = 0

IF (B==0) THEN

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

    DO mu = 1, 4 
    DO nu = 1, 4

      Alpha_A = Alpha_A + U_u(nu) * dU_d_dX(nu,mu) * U_u(nu) * dU_u_dX(nu,mu)

    END DO
    END DO

    Alpha_A = SQRT(ABS(Alpha_A)) / W

    Alpha = C * ( Alpha_Eig + Alpha_A )

END IF


IF (B==1) THEN

    l_mu(1) = W * V_u_1 + W * V_u_2 + W * V_u_3  
    IF (Vsq .NE. 0.0_DP) THEN

      l_mu(2) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_1**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_2 / Vsq &
              + ( W - 1.0_DP ) * V_u_1 * V_u_3 / Vsq 

      l_mu(3) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_2**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_2 / Vsq &
              + ( W - 1.0_DP ) * V_u_2 * V_u_3 / Vsq 

      l_mu(4) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_3**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_3 / Vsq &
              + ( W - 1.0_DP ) * V_u_2 * V_u_3 / Vsq 
    ELSE

      l_mu(2) = 0.0_DP
      l_mu(3) = 0.0_DP
      l_mu(4) = 0.0_DP

    END IF

    Alpha = 0.0_DP
 


    DO mu = 1, 4 
    DO nu = 1, 4

      ! Alpha = Alpha + ( l_mu(nu) * U_u(mu) + l_mu(mu) * l_mu(nu) ) * dU_dX(mu,nu)
       Alpha = Alpha + 0.5_DP * ( U_u(mu) + l_mu(mu) ) * ( dU_d_dX(mu,nu) + dU_d_dX(nu,mu) ) &
             * ( U_u(nu) + l_mu(nu))
    END DO
    END DO
 
    Alpha = C * ABS( Alpha ) / W

END IF




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
     INTEGER :: iZ2, iZ3, iZ4, iNodeX, mu, nu, rho, sig
     REAL(DP) :: G_dd_11, G_dd_22, G_dd_33, B(3), A


    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

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
                                                + dG_munu_dXrho(iNodeX,sig,mu,nu,iZ2,iZ3,iZ4) )
        END DO
        END DO
        END DO
        END DO

      END DO

    END DO
    END DO
    END DO


  END SUBROUTINE ComputeChristoffel

  SUBROUTINE CheckRealizability(uCR_K, uCF_K, GX_K, iZ_B1, iZ_E1, iZ_B0, iZ_E0)

    INTEGER,  INTENT(in)    :: &
       iZ_B1(4), iZ_E1(4), iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      uCF_K(nDOFX, &
            iZ_B1(2):iZ_E1(2), &
            iZ_B1(3):iZ_E1(3), &
            iZ_B1(4):iZ_E1(4), &
            1: nCF)

    REAL(DP), INTENT(in)  :: &
      GX_K(nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(inout) :: &
      uCR_K (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP):: D, I1, I2, I3
    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iNodeX, iNodeZ
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B1(2):iZ_E1(2), &
            iZ_B1(3):iZ_E1(3), &
            iZ_B1(4):iZ_E1(4))

    REAL(DP) :: uPR_K(nPR), W
    INTEGER:: n
    n=0
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX       ,iZ2,iZ3,iZ4, iCF_D), &
                 uCF_K(iNodeX      ,iZ2,iZ3,iZ4,iCF_S1), &
                 uCF_K(iNodeX      ,iZ2,iZ3,iZ4,iCF_S2), &
                 uCF_K(iNodeX      ,iZ2,iZ3,iZ4,iCF_S3), &
                 uCF_K(iNodeX       ,iZ2,iZ3,iZ4,iCF_E), &
                 uCF_K(iNodeX,iZ2,iZ3,iZ4,iCF_Ne), &
                 uPF_K(iNodeX,iPF_D       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_E       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iNodeZ = 1, nDOFZ
        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        W = 1.0_DP / SQRT( 1.0_DP - uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4)**2)
!        IF( iZ1 .EQ. 0 .AND. iZ2 .EQ. 0 ) THEN
 !         CYCLE
  !      END IF

   !     IF (iZ1 .EQ. iZ_E1(1) .AND. iZ2 .EQ. iZ_E1(2)) THEN
    !      CYCLE
     !   END IF

      !  IF (iZ1 .EQ. iZ_E1(1) .AND. iZ2 .EQ. 0) THEN
       !   CYCLE
       ! END IF

       ! IF (iZ1 .EQ. 0 .AND. iZ2 .EQ. iZ_E1(2)) THEN
       !   CYCLE
       ! END IF
      IF (uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS) <= 0.0_DP) THEN
        ! uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS) = 1.d-8
      END IF 
      IF((W * uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS) -ABS(uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS) )) .LT. 0.0_DP) THEN
         print*, iZ1, iZ2, W
         print*, "N ", uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS) 
         print*, "G1 ", uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS)  
         print*, "G2 ", uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G2, iS) 
         print*, "G3 ", uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G3, iS) 
         print*, "WN-ABS(G1) ", W *uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS) -ABS(uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS)) 
         n = n+1
      !   uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS) &
       !    = SIGN( W * uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
        !               uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) )
      END IF
        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS), &
                 uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS), &
                 uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G2, iS), &
                 uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G3, iS), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX_K (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 GX_K(iNodeX ,iZ2,iZ3,iZ4,iGF_Alpha  ), &
                 GX_K(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_1  ), &
                 GX_K(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_2  ), &
                 GX_K(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_3  ) )
      IF((W * uPR_K(iPR_D)-ABS(uPR_K(iPR_I1))) .LT. 0.0_DP) THEN
         print*, iZ1, iZ2
         print*, "D ", uPR_K(iPR_D)
         print*, "I1 ", uPR_K(iPR_I1) 
         print*, "I2 ", uPR_K(iPR_I2) 
         print*, "I3 ", uPR_K(iPR_I3) 
         print*, "WD-ABS(I) ", W*uPR_K(iPR_D) - ABS(uPR_K(iPR_I1))
         n = n+1
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO




  END SUBROUTINE CheckRealizability




END MODULE TwoMoment_DiscretizationModule_Streaming_Relativistic
