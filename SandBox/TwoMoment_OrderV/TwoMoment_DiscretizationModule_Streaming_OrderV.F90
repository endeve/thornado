MODULE TwoMoment_DiscretizationModule_Streaming_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Explicit, &
    Timer_Ex_Permute, &
    Timer_Ex_Interpolate, &
    Timer_Ex_Flux, &
    Timer_Ex_Increment
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, &
    Weights_q, &
    Weights_X1
  USE ReferenceElementModule_Lagrange, ONLY: &
    dLdX1_q, &
    L_X1_Dn, L_X1_Up
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
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
    Euler_ApplyBoundaryConditions
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    Euler_ComputePrimitive_NonRelativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_X1, &
    NumericalFlux_LLF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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
      U_R (1:nDOFZ ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_Explicit )

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (Start)"

    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    ASSOCIATE &
      ( dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL Euler_ApplyBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

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

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- Multiply Inverse Mass Matrix ---

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
              / ( Weights_q(iNodeZ) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                    * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ2, etc.

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (End)"

    CALL TimersStop( Timer_Explicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX, iDim
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iPF
    INTEGER  :: nZ(4), nZ_X1(4), nK, nF, nF_X1
    INTEGER  :: nIterations
    REAL(DP) :: uPF_L(nDOFX_X1,nPF), P_L(nDOFX_X1), Cs_L(nDOFX_X1)
    REAL(DP) :: uPF_R(nDOFX_X1,nPF), P_R(nDOFX_X1), Cs_R(nDOFX_X1)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: GX_K(nDOFX,   iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                              iZ_B1(2):iZ_E1(2),nGF)
    REAL(DP) :: GX_F(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                           iZ_B0(2):iZ_E1(2),nGF)
    REAL(DP) :: G_K(nDOFZ,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    iZ_B0(2):iZ_E0(2)  )
    REAL(DP) :: G_F(nDOF_X1,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                                iZ_B0(2):iZ_E1(2))
    REAL(DP) :: uCF_K(nDOFX,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                            iZ_B1(2):iZ_E1(2),nCF)
    REAL(DP) :: uCF_L(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                            iZ_B0(2):iZ_E1(2),nCF)
    REAL(DP) :: uCF_R(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                            iZ_B0(2):iZ_E1(2),nCF)
    REAL(DP) :: uPF_KX(nDOFX,nPF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       iZ_B0(2):iZ_E0(2))
    REAL(DP) :: uPF_K(nDOFZ,nPF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                      iZ_B0(2):iZ_E0(2))
    REAL(DP) :: V_u_X1(3,nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       iZ_B0(2):iZ_E1(2))
    REAL(DP) :: V_u(3,nDOF_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    iZ_B0(2):iZ_E1(2))
    REAL(DP) :: uCR_K(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B1(2):iZ_E1(2))
    REAL(DP) :: uCR_L(nDOF_X1,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: uCR_R(nDOF_X1,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: NumericalFlux(nDOF_X1,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: dU_X1(nDOFZ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E0(2)  )
    REAL(DP) :: Flux_q(nDOFZ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                       iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E0(2)  )

    PRINT*, "      ComputeIncrement_Divergence_X1"

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    nZ    = iZ_E0 - iZ_B0 + 1
    nZ_X1 = nZ + [0,1,0,0]
    nK    = nSpecies * nCR * PRODUCT( nZ )
    nF    = nSpecies * nCR * PRODUCT( nZ_X1 )
    nF_X1 = PRODUCT( nZ_X1(2:4) )

    ASSOCIATE &
      ( dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Geometry Fields ---

    DO iGF = 1, nGF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale Factors ---

    DO iGF = iGF_h_1, iGF_h_3

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    DO iDim = 1, 3
    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11+iDim-1) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1+iDim-1)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO
    END DO

    DO iGF = iGF_h_1, iGF_h_3
    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1)*nDOFE + iNodeE

        G_F(iNodeZ,iGF,iZ3,iZ4,iZ2) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF), SqrtTiny )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNode = 1, nDOF_X1

        G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
          =   G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2) &
            * G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2) &
            * G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Fluid Fields ---

    DO iCF = 1, nCF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

    END DO

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Compute Face Velocity Components ---

    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      ! --- Left States ---

      CALL Euler_ComputePrimitive_NonRelativistic &
             ( uCF_L(:,iZ3,iZ4,iZ2,iCF_D ), &
               uCF_L(:,iZ3,iZ4,iZ2,iCF_S1), &
               uCF_L(:,iZ3,iZ4,iZ2,iCF_S2), &
               uCF_L(:,iZ3,iZ4,iZ2,iCF_S3), &
               uCF_L(:,iZ3,iZ4,iZ2,iCF_E ), &
               uCF_L(:,iZ3,iZ4,iZ2,iCF_Ne), &
               uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), &
               uPF_L(:,iPF_V2), uPF_L(:,iPF_V3), &
               uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      ! --- Right States ---

      CALL Euler_ComputePrimitive_NonRelativistic &
             ( uCF_R(:,iZ3,iZ4,iZ2,iCF_D ), &
               uCF_R(:,iZ3,iZ4,iZ2,iCF_S1), &
               uCF_R(:,iZ3,iZ4,iZ2,iCF_S2), &
               uCF_R(:,iZ3,iZ4,iZ2,iCF_S3), &
               uCF_R(:,iZ3,iZ4,iZ2,iCF_E ), &
               uCF_R(:,iZ3,iZ4,iZ2,iCF_Ne), &
               uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), &
               uPF_R(:,iPF_V2), uPF_R(:,iPF_V3), &
               uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
               GX_F(:,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNode = 1, nDOFX_X1

        V_u_X1(1:3,iNode,iZ3,iZ4,iZ2) &
          = FaceVelocity_X1 &
              ( uPF_L(iNode,iPF_D ), uPF_L(iNode,iPF_V1), &
                uPF_L(iNode,iPF_V2), uPF_L(iNode,iPF_V3), &
                P_L(iNode), Cs_L(iNode), &
                uPF_R(iNode,iPF_D ), uPF_R(iNode,iPF_V1), &
                uPF_R(iNode,iPF_V2), uPF_R(iNode,iPF_V3), &
                P_R(iNode), Cs_R(iNode), &
                GX_F(iNode,iZ3,iZ4,iZ2,iGF_Gm_dd_11) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Expand Face Velocity ---

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1)*nDOFE + iNodeE

        V_u(1:3,iNodeZ,iZ3,iZ4,iZ2) = V_u_X1(1:3,iNodeX,iZ3,iZ4,iZ2)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Radiation Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iZ1,iZ3,iZ4,iCR,iS,iZ2) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nF, nDOFZ, One, L_X1_Up, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)-1), nDOFZ, Zero, &
             uCR_L(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nF, nDOFZ, One, L_X1_Dn, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOFZ, Zero, &
             uCR_R(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

    ! --- Numerical Flux ---

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOF_X1

        ! --- Left State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 V_u(1,iNode,iZ3,iZ4,iZ2), &
                 V_u(2,iNode,iZ3,iZ4,iZ2), &
                 V_u(3,iNode,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        ! --- Left State Flux ---

        Flux_L = Flux_X1 &
                   ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                     uPR_L(iPR_I2), uPR_L(iPR_I3), &
                     V_u(1,iNode,iZ3,iZ4,iZ2), &
                     V_u(2,iNode,iZ3,iZ4,iZ2), &
                     V_u(3,iNode,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 V_u(1,iNode,iZ3,iZ4,iZ2), &
                 V_u(2,iNode,iZ3,iZ4,iZ2), &
                 V_u(3,iNode,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        ! --- Right State Flux ---

        Flux_R = Flux_X1 &
                   ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                     uPR_R(iPR_I2), uPR_R(iPR_I3), &
                     V_u(1,iNode,iZ3,iZ4,iZ2), &
                     V_u(2,iNode,iZ3,iZ4,iZ2), &
                     V_u(3,iNode,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        DO iCR = 1, nCR

          NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = NumericalFlux_LLF &
                ( uCR_L(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                  uCR_R(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ3(iZ3) * dZ4(iZ4) * Weights_X1(iNode) &
                * G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

        END DO

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_X1, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_X1, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart( Timer_Ex_Permute )

    DO iGF = 1, nGF
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        G_K(iNodeZ,iGF,iZ3,iZ4,iZ2) = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      CALL Euler_ComputePrimitive_NonRelativistic &
             ( uCF_K (:,iZ3,iZ4,iZ2,iCF_D ), &
               uCF_K (:,iZ3,iZ4,iZ2,iCF_S1), &
               uCF_K (:,iZ3,iZ4,iZ2,iCF_S2), &
               uCF_K (:,iZ3,iZ4,iZ2,iCF_S3), &
               uCF_K (:,iZ3,iZ4,iZ2,iCF_E ), &
               uCF_K (:,iZ3,iZ4,iZ2,iCF_Ne), &
               uPF_KX(:,iPF_D ,iZ3,iZ4,iZ2), &
               uPF_KX(:,iPF_V1,iZ3,iZ4,iZ2), &
               uPF_KX(:,iPF_V2,iZ3,iZ4,iZ2), &
               uPF_KX(:,iPF_V3,iZ3,iZ4,iZ2), &
               uPF_KX(:,iPF_E ,iZ3,iZ4,iZ2), &
               uPF_KX(:,iPF_Ne,iZ3,iZ4,iZ2), &
               GX_K  (:,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
               GX_K  (:,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
               GX_K  (:,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Expand Primitive Fluid ---

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iPF = 1, nPF

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        uPF_K(iNodeZ,iPF,iZ3,iZ4,iZ2) = uPF_KX(iNodeX,iPF,iZ3,iZ4,iZ2)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Flux )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOFZ

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNode,iPF_V1,iZ3,iZ4,iZ2), &
                 uPF_K(iNode,iPF_V2,iZ3,iZ4,iZ2), &
                 uPF_K(iNode,iPF_V3,iZ3,iZ4,iZ2), &
                 G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        Flux_K &
          = Flux_X1 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                uPF_K(iNode,iPF_V1,iZ3,iZ4,iZ2), &
                uPF_K(iNode,iPF_V2,iZ3,iZ4,iZ2), &
                uPF_K(iNode,iPF_V3,iZ3,iZ4,iZ2), &
                G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        DO iCR = 1, nCR

          Flux_q(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ3(iZ3) * dZ4(iZ4) * Weights_q(iNode) &
                * G_K(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * Flux_K(iCR)

        END DO

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X1, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOFZ

        dU_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X1(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ3, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1


  FUNCTION FaceVelocity_X1 &
    ( D_L, V1_L, V2_L, V3_L, P_L, Cs_L, &
      D_R, V1_R, V2_R, V3_R, P_R, Cs_R, &
      Gm_dd_11 )

    REAL(DP), INTENT(in) :: D_L, V1_L, V2_L, V3_L, P_L, Cs_L
    REAL(DP), INTENT(in) :: D_R, V1_R, V2_R, V3_R, P_R, Cs_R
    REAL(DP), INTENT(in) :: Gm_dd_11
    REAL(DP)             :: FaceVelocity_X1(1:3)

    REAL(DP) :: aM, aP, D_M

    aM = MAX( Zero, Cs_L / Gm_dd_11 - V1_L, Cs_R / Gm_dd_11 - V1_R )
    aP = MAX( Zero, Cs_L / Gm_dd_11 + V1_L, Cs_R / Gm_dd_11 + V1_R )

    D_M = ( aM + V1_L ) * D_L + ( aP - V1_R ) * D_R

    FaceVelocity_X1(1) &
      = ( aM + V1_L ) * D_L * V1_L + ( aP - V1_R ) * D_R * V1_R &
        - ( P_R - P_L ) / Gm_dd_11
    FaceVelocity_X1(2) &
      = ( aM + V1_L ) * D_L * V2_L + ( aP - V1_R ) * D_R * V2_R
    FaceVelocity_X1(3) &
      = ( aM + V1_L ) * D_L * V3_L + ( aP - V1_R ) * D_R * V3_R

    FaceVelocity_X1 = FaceVelocity_X1 / D_M

    RETURN
  END FUNCTION FaceVelocity_X1


END MODULE TwoMoment_DiscretizationModule_Streaming_OrderV
