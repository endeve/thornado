#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EXPLICIT
#endif

MODULE TwoMoment_DiscretizationModule_Streaming

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, &
    nNodesE, nDOFE, &
    nDOF
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Explicit, &
    Timer_Ex_In, &
    Timer_Ex_Div, &
    Timer_Ex_Geometry, &
    Timer_Ex_Permute, &
    Timer_Ex_Interpolate, &
    Timer_Ex_Flux, &
    Timer_Ex_Increment, &
    Timer_Ex_Out
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_q, &
    Weights_X1, &
    Weights_X2, &
    Weights_X3, &
    NodeNumberTable, &
    OuterProduct1D3D, &
    NodeNumbersX
  USE ReferenceElementModule_Lagrange, ONLY: &
    dLdX1_q, &
    dLdX2_q, &
    dLdX3_q, &
    L_X1_Dn, &
    L_X1_Up, &
    L_X2_Dn, &
    L_X2_Up, &
    L_X3_Dn, &
    L_X3_Up
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
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
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
    iGE_Ep0, &
    iGE_Ep2
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_X1, &
    Flux_X2, &
    Flux_X3, &
    NumericalFlux_LLF, &
    StressTensor_Diagonal

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: Ones(16) = One

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iNodeX, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS
    REAL(DP) :: Tau

    CALL TimersStart( Timer_Explicit )

    ASSOCIATE ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
                dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GX, U, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP          dZ1, dZ2, dZ3, dZ4 ) &
    !$OMP MAP( alloc: dU )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         dZ1, dZ2, dZ3, dZ4 ) &
    !$ACC CREATE( dU )
#endif

    CALL TimersStop( Timer_Ex_In )

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU, iZ_B1, iZ_E1 )
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
                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Div )

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL TimersStop( Timer_Ex_Div )

    ! --- Multiply Inverse Mass Matrix ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, Tau )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, Tau ) &
    !$ACC PRESENT( GX, dU, dZ2, dZ3, dZ4, Weights_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, Tau )
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  iNodeX = MOD( (iNode-1) / nDOFE, nDOFX ) + 1
                  Tau = One * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                    = dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                        / ( Weights_q(iNode) * Tau * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Geometry )

    CALL ComputeIncrement_Geometry &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL TimersStop( Timer_Ex_Geometry )

#ifdef THORNADO_DEBUG_EXPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU )
#endif
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU)', MAXLOC(dU)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU)', MAXVAL(dU)
#endif

    CALL TimersStart( Timer_Ex_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU, U ) &
    !$OMP MAP( release: GX, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               dZ1, dZ2, dZ3, dZ4 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU, U ) &
    !$ACC DELETE( GX, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         dZ1, dZ2, dZ3, dZ4 )
#endif

    CALL TimersStop( Timer_Ex_Out )

    END ASSOCIATE

    CALL TimersStop( Timer_Explicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: nZ(4), nZ_X1(4), nK, nF, nF_GF
    INTEGER  :: iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF
    INTEGER  :: iNodeX, iNodeE, iNodeZ
    REAL(DP) :: FF, EF, FF_L, EF_L, FF_R, EF_R
    REAL(DP) :: alpha
    REAL(DP) :: Tau
    REAL(DP) :: Tau_X1
    REAL(DP) :: absLambda_L
    REAL(DP) :: absLambda_R
    REAL(DP), DIMENSION(nPR) :: uPR_L, uPR_R
    REAL(DP), DIMENSION(nPR) :: uPR_K
    REAL(DP), DIMENSION(nCR) :: Flux_X1_L
    REAL(DP), DIMENSION(nCR) :: Flux_X1_R
    REAL(DP), DIMENSION(nCR) :: Flux_X1_K

    REAL(DP) :: GX_K         (nDOFX       ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)-1:iZ_E0(2)+1,nGF)
    REAL(DP) :: GX_F         (nDOFX_X1    ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)+1,nGF)
    REAL(DP) :: G_K          (nDOF    ,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)  )
    REAL(DP) :: G_F          (nDOF_X1 ,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)+1)

    REAL(DP) :: uCR_K        (nDOF    ,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2)-1:iZ_E0(2)+1)
    REAL(DP) :: uCR_L        (nDOF_X1 ,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2)  :iZ_E0(2)+1)
    REAL(DP) :: uCR_R        (nDOF_X1 ,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2)  :iZ_E0(2)+1)

    REAL(DP) :: dU_X1        (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2)  :iZ_E0(2)  )
    REAL(DP) :: Flux_X1_q    (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2)  :iZ_E0(2)  )
    REAL(DP) :: NumericalFlux(nDOF_X1 ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2)  :iZ_E0(2)+1)

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    nZ = iZ_E0 - iZ_B0 + 1
    nZ_X1 = nZ + [0,1,0,0]
    nK = nSpecies * nCR * PRODUCT( nZ )
    nF = nSpecies * nCR * PRODUCT( nZ_X1 )
    nF_GF = PRODUCT( nZ_X1(2:4) )

    ASSOCIATE ( dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP             dU_X1, Flux_X1_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X1, Flux_X1_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_In )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
      DO iZ2 = iZ_B0(2) - 1, iZ_E0(2) + 1
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iNode = 1, nDOFX
              GX_K(iNode,iZ3,iZ4,iZ2,iGF) = GX(iNode,iZ2,iZ3,iZ4,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale Factors ---

    DO iGF = iGF_h_1, iGF_h_3
      
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_GF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_GF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    ! --- Lapse Function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_GF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF_Alpha), nDOFX, Zero, &
             GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_Alpha), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_GF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_Alpha), nDOFX, Half, &
             GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_Alpha), nDOFX_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iNodeX = 1, nDOFX_X1
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_F(iNodeZ,iGF,iZ3,iZ4,iZ2) = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF), SqrtTiny )
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#endif
  DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iNodeX = 1, nDOFX_X1
          DO iNodeE = 1, nDOFE
            iNodeZ = (iNodeX-1)*nDOFE + iNodeE
            G_F(iNodeZ,iGF_Alpha,iZ3,iZ4,iZ2) = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Alpha), SqrtTiny )
          END DO
        END DO
      END DO
    END DO
  END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iNode = 1, nDOF_X1
            G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2) = MAX( G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2) = MAX( G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) = MAX( G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)**2, SqrtTiny )
            G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
              =   G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2) &
                * G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2) &
                * G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iZ2 = iZ_B0(2) - 1, iZ_E0(2) + 1
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF
                  uCR_K(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Radiation Fields ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Up, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)-1), nDOF, Zero, &
             uCR_L(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Dn, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF, Zero, &
             uCR_R(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X1, alpha, absLambda_L, absLambda_R )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( Flux_X1_L, Flux_X1_R, FF_L, EF_L, FF_R, EF_R, &
    !$ACC          uPR_L, uPR_R, Tau_X1, alpha, absLambda_L, absLambda_R ) &
    !$ACC PRESENT( uCR_L, uCR_R, NumericalFlux, G_F, dZ3, dZ4, Weights_X1, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X1_L, Flux_X1_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X1, alpha, absLambda_L, absLambda_R )
#endif
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF_X1

                ! --- Left State Primitive, etc. ---

                uPR_L(iPR_D ) = uCR_L(iNode,iZ1,iZ3,iZ4,iCR_N ,iS,iZ2)
                uPR_L(iPR_I1) = uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2) / G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2)
                uPR_L(iPR_I2) = uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2) / G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2)
                uPR_L(iPR_I3) = uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2) / G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2)

                FF_L = FluxFactor &
                       ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                         uPR_L(iPR_I2), uPR_L(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                         G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                         G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )
                EF_L = EddingtonFactor( uPR_L(iPR_D), FF_L )
                Flux_X1_L(1:nCR) &
                  = Flux_X1 &
                      ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                        uPR_L(iPR_I2), uPR_L(iPR_I3), &
                        FF_L, EF_L, &
                        G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                        G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                        G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

                ! --- Right State Primitive, etc. ---

                uPR_R(iPR_D ) = uCR_R(iNode,iZ1,iZ3,iZ4,iCR_N ,iS,iZ2)
                uPR_R(iPR_I1) = uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2) / G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2)
                uPR_R(iPR_I2) = uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2) / G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2)
                uPR_R(iPR_I3) = uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2) / G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2)

                FF_R = FluxFactor &
                       ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                         uPR_R(iPR_I2), uPR_R(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                         G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                         G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )
                EF_R = EddingtonFactor( uPR_R(iPR_D), FF_R )
                Flux_X1_R(1:nCR) &
                  = Flux_X1 &
                      ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                        uPR_R(iPR_I2), uPR_R(iPR_I3), &
                        FF_R, EF_R, &
                        G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                        G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                        G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

                DO iCR = 1, nCR

                  absLambda_L = 1.0_DP
                  absLambda_R = 1.0_DP

                  Tau_X1 = One * G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2)
                  alpha = MAX( absLambda_L, absLambda_R )

                  NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
                    = NumericalFlux_LLF &
                        ( uCR_L(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                          uCR_R(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                          Flux_X1_L(iCR), &
                          Flux_X1_R(iCR), &
                          alpha )

                  NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
                    = dZ3(iZ3) * dZ4(iZ4) * Weights_X1(iNode) * Tau_X1 &
                        * G_F(iNode,iGF_Alpha,iZ3,iZ4,iZ2) &
                        * NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Surface Contribution ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOF_X1, Zero, &
             dU_X1, nDOF )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), nDOF_X1, One, &
             dU_X1, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Volume Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_K, GX_K, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iGF = 1, nGF
      DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iNodeX = 1, nDOFX
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_K(iNodeZ,iGF,iZ3,iZ4,iZ2) = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X1_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( EF, FF, uPR_K, Tau, Flux_X1_K ) &
    !$ACC PRESENT( Flux_X1_q, uCR_K, G_K, dZ3, dZ4, Weights_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X1_K )
#endif
    DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF

                uPR_K(iPR_D ) = uCR_K(iNode,iZ1,iZ3,iZ4,iCR_N ,iS,iZ2)
                uPR_K(iPR_I1) = uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2) / G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2)
                uPR_K(iPR_I2) = uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2) / G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2)
                uPR_K(iPR_I3) = uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2) / G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2)

                FF = FluxFactor &
                       ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                         uPR_K(iPR_I2), uPR_K(iPR_I3), &
                         G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                         G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                         G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

                EF = EddingtonFactor( uPR_K(iPR_D), FF )

                Flux_X1_K(1:nCR) &
                  = Flux_X1 &
                      ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                        uPR_K(iPR_I2), uPR_K(iPR_I3), &
                        FF, EF, &
                        G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                        G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                        G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

                DO iCR = 1, nCR

                  Flux_X1_q(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) = Flux_X1_K(iCR)

                  Tau = One * G_K(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2)

                  Flux_X1_q(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
                    = dZ3(iZ3) * dZ4(iZ4) * Weights_q(iNode) * Tau &
                        * G_K(iNode,iGF_Alpha,iZ3,iZ4,iZ2) &
                        * Flux_X1_q(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Contribution from Volume ---

    CALL TimersStart( Timer_Ex_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF, One, dLdX1_q, nDOF, &
             Flux_X1_q, nDOF, One, dU_X1, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_X1, dU, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                    =   dU  (iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                      + dU_X1(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

#ifdef THORNADO_DEBUG_EXPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X1 )
#endif
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_X1)', MAXLOC(dU_X1)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_X1)', MAXVAL(dU_X1)
#endif

    CALL TimersStart( Timer_Ex_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP               dU_X1, Flux_X1_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ3, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X1, Flux_X1_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_Out )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: nZ(4), nZ_X2(4), nK, nF, nF_GF
    INTEGER  :: iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF
    INTEGER  :: iNodeX, iNodeE, iNodeZ
    REAL(DP) :: FF, EF, FF_L, EF_L, FF_R, EF_R
    REAL(DP) :: alpha
    REAL(DP) :: Tau
    REAL(DP) :: Tau_X2
    REAL(DP) :: absLambda_L
    REAL(DP) :: absLambda_R
    REAL(DP), DIMENSION(nPR) :: uPR_L, uPR_R
    REAL(DP), DIMENSION(nPR) :: uPR_K
    REAL(DP), DIMENSION(nCR) :: Flux_X2_L
    REAL(DP), DIMENSION(nCR) :: Flux_X2_R
    REAL(DP), DIMENSION(nCR) :: Flux_X2_K

    REAL(DP) :: GX_K         (nDOFX       ,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)-1:iZ_E0(3)+1,nGF)
    REAL(DP) :: GX_F         (nDOFX_X2    ,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)  :iZ_E0(3)+1,nGF)
    REAL(DP) :: G_K          (nDOF    ,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)  :iZ_E0(3)  )
    REAL(DP) :: G_F          (nDOF_X2 ,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)  :iZ_E0(3)+1)

    REAL(DP) :: uCR_K        (nDOF    ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(3)-1:iZ_E0(3)+1)
    REAL(DP) :: uCR_L        (nDOF_X2 ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(3)  :iZ_E0(3)+1)
    REAL(DP) :: uCR_R        (nDOF_X2 ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(3)  :iZ_E0(3)+1)

    REAL(DP) :: dU_X2        (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(3)  :iZ_E0(3)  )
    REAL(DP) :: Flux_X2_q    (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(3)  :iZ_E0(3)  )
    REAL(DP) :: NumericalFlux(nDOF_X2 ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(3)  :iZ_E0(3)+1)

    IF( iZ_E0(3) .EQ. iZ_B0(3) ) RETURN

    nZ = iZ_E0 - iZ_B0 + 1
    nZ_X2 = nZ + [0,0,1,0]
    nK = nSpecies * nCR * PRODUCT( nZ )
    nF = nSpecies * nCR * PRODUCT( nZ_X2 )
    nF_GF = PRODUCT( nZ_X2(2:4) )

    ASSOCIATE ( dZ2 => MeshX(1) % Width, dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP             dU_X2, Flux_X2_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X2, Flux_X2_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_In )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
      DO iZ3 = iZ_B0(3) - 1, iZ_E0(3) + 1
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNode = 1, nDOFX
              GX_K(iNode,iZ2,iZ4,iZ3,iGF) = GX(iNode,iZ2,iZ3,iZ4,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale Factors ---

    DO iGF = iGF_h_1, iGF_h_3
      
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nF_GF, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               GX_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX_X2 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nF_GF, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               GX_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF), nDOFX_X2 )

    END DO

    ! --- Lapse Function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nF_GF, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             GX_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1,iGF_Alpha), nDOFX, Zero, &
             GX_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF_Alpha), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nF_GF, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             GX_K(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF_Alpha), nDOFX, Half, &
             GX_F(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF_Alpha), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNodeX = 1, nDOFX_X2
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_F(iNodeZ,iGF,iZ2,iZ4,iZ3) = MAX( GX_F(iNodeX,iZ2,iZ4,iZ3,iGF), SqrtTiny )
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#endif
  DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iNodeX = 1, nDOFX_X2
          DO iNodeE = 1, nDOFE
            iNodeZ = (iNodeX-1)*nDOFE + iNodeE
            G_F(iNodeZ,iGF_Alpha,iZ2,iZ4,iZ3) = MAX( GX_F(iNodeX,iZ2,iZ4,iZ3,iGF_Alpha), SqrtTiny )
          END DO
        END DO
      END DO
    END DO
  END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iNode = 1, nDOF_X2
            G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3) = MAX( G_F(iNode,iGF_h_1,iZ2,iZ4,iZ3)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3) = MAX( G_F(iNode,iGF_h_2,iZ2,iZ4,iZ3)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) = MAX( G_F(iNode,iGF_h_3,iZ2,iZ4,iZ3)**2, SqrtTiny )
            G_F(iNode,iGF_SqrtGm,iZ2,iZ4,iZ3) &
              =   G_F(iNode,iGF_h_1,iZ2,iZ4,iZ3) &
                * G_F(iNode,iGF_h_2,iZ2,iZ4,iZ3) &
                * G_F(iNode,iGF_h_3,iZ2,iZ4,iZ3)
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iZ3 = iZ_B0(3) - 1, iZ_E0(3) + 1
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF
                  uCR_K(iNode,iZ1,iZ2,iZ4,iCR,iS,iZ3) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Radiation Fields ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X2, nF, nDOF, One, L_X2_Up, nDOF_X2, &
             uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)-1), nDOF, Zero, &
             uCR_L(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)  ), nDOF_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X2, nF, nDOF, One, L_X2_Dn, nDOF_X2, &
             uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)  ), nDOF, Zero, &
             uCR_R(1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)  ), nDOF_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X2, alpha, absLambda_L, absLambda_R )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( Flux_X2_L, Flux_X2_R, FF_L, EF_L, FF_R, EF_R, &
    !$ACC          uPR_L, uPR_R, Tau_X2, alpha, absLambda_L, absLambda_R ) &
    !$ACC PRESENT( uCR_L, uCR_R, G_F, NumericalFlux, dZ2, dZ4, Weights_X2, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X2_L, Flux_X2_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X2, alpha, absLambda_L, absLambda_R )
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF_X2

                ! --- Left State Primitive, etc. ---

                uPR_L(iPR_D ) = uCR_L(iNode,iZ1,iZ2,iZ4,iCR_N ,iS,iZ3)
                uPR_L(iPR_I1) = uCR_L(iNode,iZ1,iZ2,iZ4,iCR_G1,iS,iZ3) / G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3)
                uPR_L(iPR_I2) = uCR_L(iNode,iZ1,iZ2,iZ4,iCR_G2,iS,iZ3) / G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3)
                uPR_L(iPR_I3) = uCR_L(iNode,iZ1,iZ2,iZ4,iCR_G3,iS,iZ3) / G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3)

                FF_L = FluxFactor &
                       ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                         uPR_L(iPR_I2), uPR_L(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                         G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                         G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )
                EF_L = EddingtonFactor( uPR_L(iPR_D), FF_L )
                Flux_X2_L(1:nCR) &
                  = Flux_X2 &
                      ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                        uPR_L(iPR_I2), uPR_L(iPR_I3), &
                        FF_L, EF_L, &
                        G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                        G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                        G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

                ! --- Right State Primitive, etc. ---

                uPR_R(iPR_D ) = uCR_R(iNode,iZ1,iZ2,iZ4,iCR_N ,iS,iZ3)
                uPR_R(iPR_I1) = uCR_R(iNode,iZ1,iZ2,iZ4,iCR_G1,iS,iZ3) / G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3)
                uPR_R(iPR_I2) = uCR_R(iNode,iZ1,iZ2,iZ4,iCR_G2,iS,iZ3) / G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3)
                uPR_R(iPR_I3) = uCR_R(iNode,iZ1,iZ2,iZ4,iCR_G3,iS,iZ3) / G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3)

                FF_R = FluxFactor &
                       ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                         uPR_R(iPR_I2), uPR_R(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                         G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                         G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )
                EF_R = EddingtonFactor( uPR_R(iPR_D), FF_R )
                Flux_X2_R(1:nCR) &
                  = Flux_X2 &
                      ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                        uPR_R(iPR_I2), uPR_R(iPR_I3), &
                        FF_R, EF_R, &
                        G_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                        G_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                        G_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

                DO iCR = 1, nCR

                  absLambda_L = 1.0_DP
                  absLambda_R = 1.0_DP

                  Tau_X2 = One * G_F(iNode,iGF_SqrtGm,iZ2,iZ4,iZ3)
                  alpha = MAX( absLambda_L, absLambda_R )

                  NumericalFlux(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
                    = NumericalFlux_LLF &
                        ( uCR_L(iNode,iZ1,iZ2,iZ4,iCR,iS,iZ3), &
                          uCR_R(iNode,iZ1,iZ2,iZ4,iCR,iS,iZ3), &
                          Flux_X2_L(iCR), &
                          Flux_X2_R(iCR), &
                          alpha )

                  NumericalFlux(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
                    = dZ2(iZ2) * dZ4(iZ4) * Weights_X2(iNode) * Tau_X2 &
                        * G_F(iNode,iGF_Alpha,iZ2,iZ4,iZ3) &
                        * NumericalFlux(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3)
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Surface Contribution ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X2, + One, L_X2_Dn, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), nDOF_X2, Zero, &
             dU_X2, nDOF )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X2, - One, L_X2_Up, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)+1), nDOF_X2, One, &
             dU_X2, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Volume Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_K, GX_K, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iGF = 1, nGF
            DO iNodeX = 1, nDOFX
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_K(iNodeZ,iGF,iZ2,iZ4,iZ3) = GX_K(iNodeX,iZ2,iZ4,iZ3,iGF)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X2_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( EF, FF, uPR_K, Tau, Flux_X2_K ) &
    !$ACC PRESENT( uCR_K, G_K, Flux_X2_q, dZ2, dZ4, Weights_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X2_K )
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF

                uPR_K(iPR_D ) = uCR_K(iNode,iZ1,iZ2,iZ4,iCR_N ,iS,iZ3)
                uPR_K(iPR_I1) = uCR_K(iNode,iZ1,iZ2,iZ4,iCR_G1,iS,iZ3) / G_K(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3)
                uPR_K(iPR_I2) = uCR_K(iNode,iZ1,iZ2,iZ4,iCR_G2,iS,iZ3) / G_K(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3)
                uPR_K(iPR_I3) = uCR_K(iNode,iZ1,iZ2,iZ4,iCR_G3,iS,iZ3) / G_K(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3)

                FF = FluxFactor &
                       ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                         uPR_K(iPR_I2), uPR_K(iPR_I3), &
                         G_K(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                         G_K(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                         G_K(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

                EF = EddingtonFactor( uPR_K(iPR_D), FF )

                Flux_X2_K(1:nCR) &
                  = Flux_X2 &
                      ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                        uPR_K(iPR_I2), uPR_K(iPR_I3), &
                        FF, EF, &
                        G_K(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                        G_K(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                        G_K(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

                DO iCR = 1, nCR

                  Flux_X2_q(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3) = Flux_X2_K(iCR)

                  Tau = One * G_K(iNode,iGF_SqrtGm,iZ2,iZ4,iZ3)

                  Flux_X2_q(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
                    = dZ2(iZ2) * dZ4(iZ4) * Weights_q(iNode) * Tau &
                        * G_K(iNode,iGF_Alpha,iZ2,iZ4,iZ3) &
                        * Flux_X2_q(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3)

                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Contribution from Volume ---

    CALL TimersStart( Timer_Ex_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF, One, dLdX2_q, nDOF, &
             Flux_X2_q, nDOF, One, dU_X2, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU, dU_X2, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                    =   dU  (iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                      + dU_X2(iNode,iCR,iZ1,iZ2,iZ4,iS,iZ3)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

#ifdef THORNADO_DEBUG_EXPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X2 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X2 )
#endif
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_X2)', MAXLOC(dU_X2)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_X2)', MAXVAL(dU_X2)
#endif

    CALL TimersStart( Timer_Ex_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP               dU_X2, Flux_X2_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ2, dZ4, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X2, Flux_X2_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_Out )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: nZ(4), nZ_X3(4), nK, nF, nF_GF
    INTEGER  :: iNode, iZ1, iZ2, iZ4, iZ3, iCR, iS, iGF
    INTEGER  :: iNodeX, iNodeE, iNodeZ
    REAL(DP) :: FF, EF, FF_L, EF_L, FF_R, EF_R
    REAL(DP) :: alpha
    REAL(DP) :: Tau
    REAL(DP) :: Tau_X3
    REAL(DP) :: absLambda_L
    REAL(DP) :: absLambda_R
    REAL(DP), DIMENSION(nPR) :: uPR_L, uPR_R
    REAL(DP), DIMENSION(nPR) :: uPR_K
    REAL(DP), DIMENSION(nCR) :: Flux_X3_L
    REAL(DP), DIMENSION(nCR) :: Flux_X3_R
    REAL(DP), DIMENSION(nCR) :: Flux_X3_K

    REAL(DP) :: GX_K         (nDOFX       ,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4)-1:iZ_E0(4)+1,nGF)
    REAL(DP) :: GX_F         (nDOFX_X3    ,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4)  :iZ_E0(4)+1,nGF)
    REAL(DP) :: G_K          (nDOF    ,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4)  :iZ_E0(4)  )
    REAL(DP) :: G_F          (nDOF_X3 ,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4)  :iZ_E0(4)+1)

    REAL(DP) :: uCR_K        (nDOF    ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nCR,nSpecies,iZ_B0(4)-1:iZ_E0(4)+1)
    REAL(DP) :: uCR_L        (nDOF_X3 ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nCR,nSpecies,iZ_B0(4)  :iZ_E0(4)+1)
    REAL(DP) :: uCR_R        (nDOF_X3 ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nCR,nSpecies,iZ_B0(4)  :iZ_E0(4)+1)

    REAL(DP) :: dU_X3        (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nSpecies,iZ_B0(4)  :iZ_E0(4)  )
    REAL(DP) :: Flux_X3_q    (nDOF    ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nSpecies,iZ_B0(4)  :iZ_E0(4)  )
    REAL(DP) :: NumericalFlux(nDOF_X3 ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),nSpecies,iZ_B0(4)  :iZ_E0(4)+1)

    IF( iZ_E0(4) .EQ. iZ_B0(4) ) RETURN

    nZ = iZ_E0 - iZ_B0 + 1
    nZ_X3 = nZ + [0,0,0,1]
    nK = nSpecies * nCR * PRODUCT( nZ )
    nF = nSpecies * nCR * PRODUCT( nZ_X3 )
    nF_GF = PRODUCT( nZ_X3(2:4) )

    ASSOCIATE ( dZ2 => MeshX(1) % Width, dZ3 => MeshX(2) % Width )

    CALL TimersStart( Timer_Ex_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP             dU_X3, Flux_X3_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X3, Flux_X3_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_In )

    ! --- Geometry Fields in Element Nodes ---

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = 1, nGF
      DO iZ4 = iZ_B0(4) - 1, iZ_E0(4) + 1
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNode = 1, nDOFX
              GX_K(iNode,iZ2,iZ3,iZ4,iGF) = GX(iNode,iZ2,iZ3,iZ4,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Geometry Fields on Shared Face ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale Factors ---

    DO iGF = iGF_h_1, iGF_h_3
 
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nF_GF, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               GX_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX_X3 )
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, nF_GF, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               GX_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF), nDOFX_X3 )

    END DO

    ! --- Lapse Function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nF_GF, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             GX_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1,iGF_Alpha), nDOFX, Zero, &
             GX_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF_Alpha), nDOFX_X3 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nF_GF, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             GX_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF_Alpha), nDOFX, Half, &
             GX_F(1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ,iGF_Alpha), nDOFX_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ4 = iZ_B0(4), iZ_E0(4) + 1
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNodeX = 1, nDOFX_X3
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_F(iNodeZ,iGF,iZ2,iZ3,iZ4) = MAX( GX_F(iNodeX,iZ2,iZ3,iZ4,iGF), SqrtTiny )
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_F, GX_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ )
#endif
  DO iZ4 = iZ_B0(4), iZ_E0(4) + 1
    DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iNodeX = 1, nDOFX_X3
          DO iNodeE = 1, nDOFE
            iNodeZ = (iNodeX-1)*nDOFE + iNodeE
            G_F(iNodeZ,iGF_Alpha,iZ2,iZ3,iZ4) = MAX( GX_F(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha), SqrtTiny )
          END DO
        END DO
      END DO
    END DO
  END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( G_F, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4) + 1
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iNode = 1, nDOF_X3
            G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4) = MAX( G_F(iNode,iGF_h_1,iZ2,iZ3,iZ4)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4) = MAX( G_F(iNode,iGF_h_2,iZ2,iZ3,iZ4)**2, SqrtTiny )
            G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) = MAX( G_F(iNode,iGF_h_3,iZ2,iZ3,iZ4)**2, SqrtTiny )
            G_F(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              =   G_F(iNode,iGF_h_1,iZ2,iZ3,iZ4) &
                * G_F(iNode,iGF_h_2,iZ2,iZ3,iZ4) &
                * G_F(iNode,iGF_h_3,iZ2,iZ3,iZ4)
          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( uCR_K, U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iZ4 = iZ_B0(4) - 1, iZ_E0(4) + 1
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF
                  uCR_K(iNode,iZ1,iZ2,iZ3,iCR,iS,iZ4) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Radiation Fields ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X3, nF, nDOF, One, L_X3_Up, nDOF_X3, &
             uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)-1), nDOF, Zero, &
             uCR_L(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)  ), nDOF_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X3, nF, nDOF, One, L_X3_Dn, nDOF_X3, &
             uCR_K(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)  ), nDOF, Zero, &
             uCR_R(1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)  ), nDOF_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X3, alpha, absLambda_L, absLambda_R )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( Flux_X3_L, Flux_X3_R, FF_L, EF_L, FF_R, EF_R, &
    !$ACC          uPR_L, uPR_R, Tau_X3, alpha, absLambda_L, absLambda_R ) &
    !$ACC PRESENT( uCR_L, uCR_R, G_F, NumericalFlux, dZ2, dZ3, Weights_X3, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Flux_X3_L, Flux_X3_R, FF_L, EF_L, FF_R, EF_R, &
    !$OMP          uPR_L, uPR_R, Tau_X3, alpha, absLambda_L, absLambda_R )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4) + 1
      DO iS = 1, nSpecies
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF_X3

                ! --- Left State Primitive, etc. ---

                uPR_L(iPR_D ) = uCR_L(iNode,iZ1,iZ2,iZ3,iCR_N ,iS,iZ4)
                uPR_L(iPR_I1) = uCR_L(iNode,iZ1,iZ2,iZ3,iCR_G1,iS,iZ4) / G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4)
                uPR_L(iPR_I2) = uCR_L(iNode,iZ1,iZ2,iZ3,iCR_G2,iS,iZ4) / G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4)
                uPR_L(iPR_I3) = uCR_L(iNode,iZ1,iZ2,iZ3,iCR_G3,iS,iZ4) / G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4)

                FF_L = FluxFactor &
                       ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                         uPR_L(iPR_I2), uPR_L(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                         G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                         G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )
                EF_L = EddingtonFactor( uPR_L(iPR_D), FF_L )
                Flux_X3_L(1:nCR) &
                  = Flux_X3 &
                      ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                        uPR_L(iPR_I2), uPR_L(iPR_I3), &
                        FF_L, EF_L, &
                        G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                        G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                        G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

                ! --- Right State Primitive, etc. ---

                uPR_R(iPR_D ) = uCR_R(iNode,iZ1,iZ2,iZ3,iCR_N ,iS,iZ4)
                uPR_R(iPR_I1) = uCR_R(iNode,iZ1,iZ2,iZ3,iCR_G1,iS,iZ4) / G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4)
                uPR_R(iPR_I2) = uCR_R(iNode,iZ1,iZ2,iZ3,iCR_G2,iS,iZ4) / G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4)
                uPR_R(iPR_I3) = uCR_R(iNode,iZ1,iZ2,iZ3,iCR_G3,iS,iZ4) / G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4)

                FF_R = FluxFactor &
                       ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                         uPR_R(iPR_I2), uPR_R(iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                         G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                         G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )
                EF_R = EddingtonFactor( uPR_R(iPR_D), FF_R )
                Flux_X3_R(1:nCR) &
                  = Flux_X3 &
                      ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                        uPR_R(iPR_I2), uPR_R(iPR_I3), &
                        FF_R, EF_R, &
                        G_F(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                        G_F(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                        G_F(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

                DO iCR = 1, nCR

                  absLambda_L = 1.0_DP
                  absLambda_R = 1.0_DP

                  Tau_X3 = One * G_F(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4)
                  alpha = MAX( absLambda_L, absLambda_R )

                  NumericalFlux(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
                    = NumericalFlux_LLF &
                        ( uCR_L(iNode,iZ1,iZ2,iZ3,iCR,iS,iZ4), &
                          uCR_R(iNode,iZ1,iZ2,iZ3,iCR,iS,iZ4), &
                          Flux_X3_L(iCR), &
                          Flux_X3_R(iCR), &
                          alpha )

                  NumericalFlux(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
                    = dZ2(iZ2) * dZ3(iZ3) * Weights_X3(iNode) * Tau_X3 &
                        * G_F(iNode,iGF_Alpha,iZ2,iZ3,iZ4) &
                        * NumericalFlux(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4)
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Surface Contribution ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X3, + One, L_X3_Dn, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), nDOF_X3, Zero, &
             dU_X3, nDOF )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X3, - One, L_X3_Up, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)+1), nDOF_X3, One, &
             dU_X3, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    !---------------------
    ! --- Volume Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( G_K, GX_K, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iGF = 1, nGF
            DO iNodeX = 1, nDOFX
              DO iNodeE = 1, nDOFE
                iNodeZ = (iNodeX-1)*nDOFE + iNodeE
                G_K(iNodeZ,iGF,iZ2,iZ3,iZ4) = GX_K(iNodeX,iZ2,iZ3,iZ4,iGF)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Numerical Flux ---

    CALL TimersStart( Timer_Ex_Flux )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X3_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( EF, FF, uPR_K, Tau, Flux_X3_K ) &
    !$ACC PRESENT( uCR_K, G_K, Flux_X3_q, dZ2, dZ3, Weights_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( EF, FF, uPR_K, Tau, Flux_X3_K )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iS = 1, nSpecies
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNode = 1, nDOF

                uPR_K(iPR_D ) = uCR_K(iNode,iZ1,iZ2,iZ3,iCR_N ,iS,iZ4)
                uPR_K(iPR_I1) = uCR_K(iNode,iZ1,iZ2,iZ3,iCR_G1,iS,iZ4) / G_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4)
                uPR_K(iPR_I2) = uCR_K(iNode,iZ1,iZ2,iZ3,iCR_G2,iS,iZ4) / G_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4)
                uPR_K(iPR_I3) = uCR_K(iNode,iZ1,iZ2,iZ3,iCR_G3,iS,iZ4) / G_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4)

                FF = FluxFactor &
                       ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                         uPR_K(iPR_I2), uPR_K(iPR_I3), &
                         G_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                         G_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                         G_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

                EF = EddingtonFactor( uPR_K(iPR_D), FF )

                Flux_X3_K(1:nCR) &
                  = Flux_X3 &
                      ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                        uPR_K(iPR_I2), uPR_K(iPR_I3), &
                        FF, EF, &
                        G_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                        G_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                        G_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

                DO iCR = 1, nCR

                  Flux_X3_q(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4) = Flux_X3_K(iCR)

                  Tau = One * G_K(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4)

                  Flux_X3_q(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
                    = dZ2(iZ2) * dZ3(iZ3) * Weights_q(iNode) * Tau &
                        * G_K(iNode,iGF_Alpha,iZ2,iZ3,iZ4) &
                        * Flux_X3_q(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4)

                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    ! --- Contribution from Volume ---

    CALL TimersStart( Timer_Ex_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF, One, dLdX3_q, nDOF, &
             Flux_X3_q, nDOF, One, dU_X3, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU, dU_X3, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                    =   dU  (iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                      + dU_X3(iNode,iCR,iZ1,iZ2,iZ3,iS,iZ4)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

#ifdef THORNADO_DEBUG_EXPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dU_X3 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dU_X3 )
#endif
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU_X3)', MAXLOC(dU_X3)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU_X3)', MAXVAL(dU_X3)
#endif

    CALL TimersStart( Timer_Ex_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$OMP               dU_X3, Flux_X3_q, NumericalFlux )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         GX_K, GX_F, G_K, G_F, uCR_K, uCR_L, uCR_R, &
    !$ACC         dU_X3, Flux_X3_q, NumericalFlux )
#endif

    CALL TimersStop( Timer_Ex_Out )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: nZ(4), nZ_X1(4), nZ_X2(4), nF_G, nF_G_X1, nF_G_X2
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iGF
    INTEGER  :: iNodeZ, iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: PR_D, PR_I1, PR_I2, PR_I3, FF, EF, Stress(3)
    REAL(DP) :: h2, h3, G11, G22, G33, dU_G1, dU_G2
    REAL(DP) :: &
      G     (nDOF   ,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4)  :iZ_E0(4)      ), &
      GX_X1 (nDOFX  ,    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)-1:iZ_E0(2)+1,nGF), &
      GX_X2 (nDOFX  ,    iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)-1:iZ_E0(3)+1,nGF)
    REAL(DP) :: &
      h2_X1(nDOFX_X1,    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)+1    ), &
      h3_X1(nDOFX_X1,    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)+1    ), &
      h3_X2(nDOFX_X2,    iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)  :iZ_E0(3)+1    )
    REAL(DP) :: &
      dh2dX1(nDOFX  ,    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)      ), &
      dh3dX1(nDOFX  ,    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2)  :iZ_E0(2)      ), &
      dh3dX2(nDOFX  ,    iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3)  :iZ_E0(3)      )

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    nZ = iZ_E0 - iZ_B0 + 1
    nZ_X1 = nZ + [0,1,0,0]
    nZ_X2 = nZ + [0,0,1,0]
    nF_G    = PRODUCT( nZ   (2:4) )
    nF_G_X1 = PRODUCT( nZ_X1(2:4) )
    nF_G_X2 = PRODUCT( nZ_X2(2:4) )

    ASSOCIATE ( dZ2 => MeshX(1) % Width, dZ3 => MeshX(2) % Width )

    CALL TimersStart( Timer_Ex_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: h2_X1, h3_X1, h3_X2, dh2dX1, dh3dX1, dh3dX2, &
    !$OMP             GX_X1, GX_X2, G )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( h2_X1, h3_X1, h3_X2, dh2dX1, dh3dX1, dh3dX2, &
    !$ACC         GX_X1, GX_X2, G )
#endif

    CALL TimersStop( Timer_Ex_In )

    ! --- X1 Face ---

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_X1, GX, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ2 = iZ_B0(2) - 1, iZ_E0(2) + 1
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iNodeX = 1, nDOFX
              GX_X1(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- h_2 and h_3 on X1 Faces ---

    CALL TimersStart( Timer_Ex_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_G_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF_h_2), nDOFX, Zero, &
             h2_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_G_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_h_2), nDOFX, Half, &
             h2_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_G_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF_h_3), nDOFX, Zero, &
             h3_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nF_G_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_h_3), nDOFX, Half, &
             h3_X1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( h2_X1, h3_X1, WeightsX_X1, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iNodeX1 = 1, nDOFX_X1

            h2_X1(iNodeX1,iZ3,iZ4,iZ2) &
              = WeightsX_X1(iNodeX1) * MAX( h2_X1(iNodeX1,iZ3,iZ4,iZ2), SqrtTiny )

            h3_X1(iNodeX1,iZ3,iZ4,iZ2) &
              = WeightsX_X1(iNodeX1) * MAX( h3_X1(iNodeX1,iZ3,iZ4,iZ2), SqrtTiny )

          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_X1, WeightsX_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iNodeX = 1, nDOFX

              GX_X1(iNodeX,iZ3,iZ4,iZ2,iGF) &
                = WeightsX_q(iNodeX) * GX_X1(iNodeX,iZ3,iZ4,iZ2,iGF)

            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Derivative of Scale Factor h_2 wrt X1 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h2_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, Zero, &
             dh2dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h2_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, One, &
             dh2dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX,    - One, dLXdX1_q, nDOFX, &
             GX_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_h_2), nDOFX, One, &
             dh2dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )

    ! --- Derivative of Scale Factor h_3 wrt X1 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h3_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, Zero, &
             dh3dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h3_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, One, &
             dh3dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX,    - One, dLXdX1_q, nDOFX, &
             GX_X1 (1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF_h_3), nDOFX, One, &
             dh3dX1(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX )

    CALL TimersStop( Timer_Ex_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( dh2dX1, dh3dX1, dZ2, WeightsX_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iNodeX = 1, nDOFX

            dh2dX1(iNodeX,iZ3,iZ4,iZ2) &
              = dh2dX1(iNodeX,iZ3,iZ4,iZ2) / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

            dh3dX1(iNodeX,iZ3,iZ4,iZ2) &
              = dh3dX1(iNodeX,iZ3,iZ4,iZ2) / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

          END DO
        END DO
      END DO
    END DO

    dh3dX2 = 0.0_DP

    IF( iZ_B1(3) < iZ_E1(3) )THEN

    ! --- X2 Face ---

    CALL TimersStart( Timer_Ex_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_X2, GX, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ3 = iZ_B1(3), iZ_E1(3)
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNodeX = 1, nDOFX
              GX_X2(iNodeX,iZ2,iZ4,iZ3,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- h_3 on X2 Faces ---

    CALL TimersStart( Timer_Ex_Interpolate )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nF_G_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             GX_X2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1,iGF_h_3), nDOFX, Zero, &
             h3_X2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nF_G_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             GX_X2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF_h_3), nDOFX, Half, &
             h3_X2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( h3_X2, WeightsX_X2, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iNodeX2 = 1, nDOFX_X2

            h3_X2(iNodeX2,iZ2,iZ4,iZ3) &
              = WeightsX_X2(iNodeX2) * MAX( h3_X2(iNodeX2,iZ2,iZ4,iZ3), SqrtTiny )

          END DO
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_X2, WeightsX_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iGF = iGF_h_1, iGF_h_3
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iNodeX = 1, nDOFX

              GX_X2(iNodeX,iZ2,iZ4,iZ3,iGF) &
                = WeightsX_q(iNodeX) * GX_X2(iNodeX,iZ2,iZ4,iZ3,iGF)

            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Derivative of Scale Factor h_3 wrt X2 ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             h3_X2 (1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, Zero, &
             dh3dX2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             h3_X2 (1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, One, &
             dh3dX2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX )
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nF_G, nDOFX,    - One, dLXdX2_q, nDOFX, &
             GX_X2 (1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ,iGF_h_3), nDOFX, One, &
             dh3dX2(1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX )

    CALL TimersStop( Timer_Ex_Interpolate )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( dh3dX2, dZ3, WeightsX_q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iNodeX = 1, nDOFX

            dh3dX2(iNodeX,iZ2,iZ4,iZ3) &
              = dh3dX2(iNodeX,iZ2,iZ4,iZ3) / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Permute )

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iNodeX )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iNodeX ) &
    !$ACC PRESENT( G, GX, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeX )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iGF = 1, nGF
            DO iNodeZ = 1, nDOF

              iNodeX = MOD( (iNodeZ-1) / nNodesE, nDOFX ) + 1

              G(iNodeZ,iGF,iZ2,iZ3,iZ4) &
                = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Add to Increments ---

    CALL TimersStart( Timer_Ex_Increment )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeX, h2, h3, G11, G22, G33, dU_G1, dU_G2, &
    !$OMP          PR_D, PR_I1, PR_I2, PR_I3, FF, EF, Stress )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNodeX, h2, h3, G11, G22, G33, dU_G1, dU_G2, &
    !$ACC          PR_D, PR_I1, PR_I2, PR_I3, FF, EF, Stress ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, dh2dX1, dh3dX1, dh3dX2, G, U, dU )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeX, h2, h3, G11, G22, G33, dU_G1, dU_G2, &
    !$OMP          PR_D, PR_I1, PR_I2, PR_I3, FF, EF, Stress )
#endif
    DO iS  = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              DO iNodeZ = 1, nDOF

                iNodeX = MOD( (iNodeZ-1) / nNodesE, nDOFX   ) + 1

                h2    = G(iNodeZ,iGF_h_2,iZ2,iZ3,iZ4)
                h3    = G(iNodeZ,iGF_h_3,iZ2,iZ3,iZ4)

                G11   = G(iNodeZ,iGF_Gm_dd_11,iZ2,iZ3,iZ4)
                G22   = G(iNodeZ,iGF_Gm_dd_22,iZ2,iZ3,iZ4)
                G33   = G(iNodeZ,iGF_Gm_dd_33,iZ2,iZ3,iZ4)

                PR_D  = U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
                PR_I1 = U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) / G11
                PR_I2 = U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) / G22
                PR_I3 = U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) / G33

                FF &
                  = FluxFactor &
                      ( PR_D, PR_I1, PR_I2, PR_I3, G11, G22, G33 )

                EF &
                  = EddingtonFactor &
                      ( PR_D, FF )

                Stress(1:3) &
                  = StressTensor_Diagonal &
                      ( PR_D, PR_I1, PR_I2, PR_I3, FF, EF, G11, G22, G33 )

                dU_G1 =   dh2dX1(iNodeX,iZ3,iZ4,iZ2) * Stress(2) / h2 &
                        + dh3dX1(iNodeX,iZ3,iZ4,iZ2) * Stress(3) / h3

                dU_G2 =   dh3dX2(iNodeX,iZ2,iZ4,iZ3) * Stress(3) / h3

                dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
                  = dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) + dU_G1

                dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
                  = dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) + dU_G2

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

#ifdef THORNADO_DEBUG_EXPLICIT
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( dh2dX1, dh3dX1, dh3dX2 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( dh2dX1, dh3dX1, dh3dX2 )
#endif
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh2dX1), MINVAL(dh2dX1)', MINLOC(dh2dX1), MINVAL(dh2dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh2dX1), MAXVAL(dh2dX1)', MAXLOC(dh2dX1), MAXVAL(dh2dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh3dX1), MINVAL(dh3dX1)', MINLOC(dh3dX1), MINVAL(dh3dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh3dX1), MAXVAL(dh3dX1)', MAXLOC(dh3dX1), MAXVAL(dh3dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh3dX2), MINVAL(dh3dX2)', MINLOC(dh3dX2), MINVAL(dh3dX2)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh3dX2), MAXVAL(dh3dX2)', MAXLOC(dh3dX2), MAXVAL(dh3dX2)
#endif

    CALL TimersStart( Timer_Ex_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               h2_X1, h3_X1, h3_X2, dh2dX1, dh3dX1, dh3dX2, &
    !$OMP               GX_X1, GX_X2, G )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( dZ2, dZ3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         h2_X1, h3_X1, h3_X2, dh2dX1, dh3dX1, dh3dX2, &
    !$ACC         GX_X1, GX_X2, G )
#endif

    CALL TimersStop( Timer_Ex_Out )

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_Geometry


END MODULE TwoMoment_DiscretizationModule_Streaming
