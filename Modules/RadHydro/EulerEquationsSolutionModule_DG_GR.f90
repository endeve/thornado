MODULE EulerEquationsSolutionModule_DG_GR

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF, iGF_Alpha, iGF_Beta_1, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    rhsCF, &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, nAF, iAF_P, iAF_Cs, iAF_Gm
  USE RiemannSolverModule, ONLY: &
     NumericalFlux_Fluid_GR
  USE EulerEquationsUtilitiesModule_GR, ONLY: &
    ComputePrimitive, Eigenvalues, &
    ComputeSoundSpeed, Flux_X1, AlphaC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeRHS_Euler_DG_GR

CONTAINS


  SUBROUTINE ComputeRHS_Euler_DG_GR( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uGF(:,iX1,iX2,iX3,1:nGF), &
                   uPF(:,iX1,iX2,iX3,1:nPF), uAF(:,iX1,iX2,iX3,1:nAF) )
          
        END DO
      END DO
    END DO

    CALL ComputeRHS_Euler_DG_GR_X1( iX_Begin, iX_End )

  END SUBROUTINE ComputeRHS_Euler_DG_GR


  SUBROUTINE ComputeRHS_Euler_DG_GR_X1( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3, iGF, iCF, iAF
    INTEGER :: iNodeX1, jNodeX1, iNodeX2, iNodeX3, iNodeX, jNodeX
    REAL(DP) :: Alpha, AlphaPls, AlphaMns, AlphaMdl, AlphaL, AlphaR, Cs
    REAL(DP), DIMENSION(1:nCF) :: Flux_L, Flux_R, Flux
    REAL(DP), DIMENSION(1:5)   :: EigVals_L, EigVals_R
    REAL(DP), DIMENSION(1:1,1:nGF) :: uGF_L, uGF_R
    REAL(DP), DIMENSION(1:1,1:nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(1:1,1:nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(1:1,1:nAF) :: uAF_L, uAF_R

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)) )

    ! -- Compute Right-Hand Side for Relativistic Euler Equations --

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          ! -- Geometry Fields ---

          ASSOCIATE &
            ( uGF_P => uGF(:,iX1-1,iX2,iX3,:), & ! Previous Element
              uGF_K => uGF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uGF_N => uGF(:,iX1+1,iX2,iX3,:) )  ! Next     Element

          ! -- Fluid Fields ---

          ! -- Conserved --

          ASSOCIATE &
            ( uCF_P => uCF(:,iX1-1,iX2,iX3,:), & ! Previous Element
              uCF_K => uCF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uCF_N => uCF(:,iX1+1,iX2,iX3,:) )  ! Next     Element

          ! -- Primitive --

          ASSOCIATE &
            ( uPF_P => uPF(:,iX1-1,iX2,iX3,:), & ! Previous Element
              uPF_K => uPF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uPF_N => uPF(:,iX1+1,iX2,iX3,:) )  ! Next     Element

          ! -- Auxiliary --

          ASSOCIATE &
            ( uAF_P => uAF(:,iX1-1,iX2,iX3,:), & ! Previous Element
              uAF_K => uAF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uAF_N => uAF(:,iX1+1,iX2,iX3,:) )  ! Next     Element

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) = 0.0_DP

                ! -- Volume Term --

! *** Deferred ***

                ! -- Left Face --

                ! -- Average Geometry Fields --

                DO iGF = 1, nGF

                  jNodeX = iNodeX
                  uGF_L(iNodeX,iGF) &
                    = 0.5_DP * ( uGF_P(jNodeX,iGF) + uGF_K(jNodeX,iGF) )

                END DO
                uGF_R = uGF_L

                ! -- Left Fluid State --

                DO iCF = 1, nCF

                  jNodeX = iNodeX
                  uCF_L(iNodeX,iCF) = uCF_P(jNodeX,iCF)

                END DO

                uAF_L(iNodeX1,iAF_P)  & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_P(:,iAF_P) )

                uAF_L(iNodeX1,iAF_Gm) & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_P(:,iAF_Gm) )

                CALL ComputePrimitive( uCF_L, uGF_L, uPF_L, uAF_L )

                Cs = ComputeSoundSpeed( uAF_L( iNodeX, iAF_P  ), &
                                        uPF_L( iNodeX, iPF_E  ), &
                                        uPF_L( iNodeX, iPF_D  ), &
                                        uAF_L( iNodeX, iAF_Gm ) )

                EigVals_L = Eigenvalues( uPF_L( iNodeX, iPF_V1 ),       &
                                         uPF_L( iNodeX, iPF_V2 ),       &
                                         uPF_L( iNodeX, iPF_V3 ),       &
                                         uGF_L( iNodeX, iGF_Gm_dd_11 ), &
                                         uGF_L( iNodeX, iGF_Gm_dd_22 ), &
                                         uGF_L( iNodeX, iGF_Gm_dd_33 ), &
                                         uGF_L( iNodeX, iGF_Alpha ),    &
                                         uGF_L( iNodeX, iGF_Beta_1 ),   &
                                         Cs )

                Flux_L &
                  = Flux_X1 &
                      ( uPF_L( iNodeX, iPF_D ),        &
                        uPF_L( iNodeX, iPF_V1 ),       &
                        uPF_L( iNodeX, iPF_V2 ),       &
                        uPF_L( iNodeX, iPF_V3 ),       &
                        uPF_L( iNodeX, iPF_E ),        &
                        uAF_L( iNodeX, iAF_P ),        &
                        uPF_L( iNodeX, iPF_Ne ),       &
                        uGF_L( iNodeX, iGF_Alpha ),    &
                        uGF_L( iNodeX, iGF_Beta_1 ),   &
                        uGF_L( iNodeX, iGF_Gm_dd_11 ), &
                        uGF_L( iNodeX, iGF_Gm_dd_22 ), &
                        uGF_L( iNodeX, iGF_Gm_dd_33 ) )

                ! -- Right Fluid State --

                DO iCF = 1, nCF

                  jNodeX = iNodeX
                  uCF_R(iNodeX,iCF) = uCF_K(jNodeX,iCF)

                END DO

                uAF_R(iNodeX1,iAF_P)  & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_K(:,iAF_Gm) )

                CALL ComputePrimitive( uCF_R, uGF_R, uPF_R, uAF_R )
                
                Cs = ComputeSoundSpeed( uAF_R( iNodeX, iAF_P  ), &
                                        uPF_R( iNodeX, iPF_E  ), &
                                        uPF_R( iNodeX, iPF_D  ), &
                                        uAF_R( iNodeX, iAF_Gm ) )

                EigVals_R = Eigenvalues( uPF_R( iNodeX, iPF_V1 ),       &
                                         uPF_R( iNodeX, iPF_V2 ),       &
                                         uPF_R( iNodeX, iPF_V3 ),       &
                                         uGF_R( iNodeX, iGF_Gm_dd_11 ), &
                                         uGF_R( iNodeX, iGF_Gm_dd_22 ), &
                                         uGF_R( iNodeX, iGF_Gm_dd_33 ), &
                                         uGF_R( iNodeX, iGF_Alpha ),    &
                                         uGF_R( iNodeX, iGF_Beta_1 ),   &
                                         Cs )

                Flux_R &
                  = Flux_X1 &
                      ( uPF_R( iNodeX, iPF_D ),        &
                        uPF_R( iNodeX, iPF_V1 ),       &
                        uPF_R( iNodeX, iPF_V2 ),       &
                        uPF_R( iNodeX, iPF_V3 ),       &
                        uPF_R( iNodeX, iPF_E ),        &
                        uAF_R( iNodeX, iAF_P ),        &
                        uPF_R( iNodeX, iPF_Ne ),       &
                        uGF_R( iNodeX, iGF_Alpha ),    &
                        uGF_R( iNodeX, iGF_Beta_1 ),   &
                        uGF_R( iNodeX, iGF_Gm_dd_11 ), &
                        uGF_R( iNodeX, iGF_Gm_dd_22 ), &
                        uGF_R( iNodeX, iGF_Gm_dd_33 ) )

                ! First Alpha will be the largest eigenvalue (absolute value) 
                ! of the left and right state
                AlphaL = MAXVAL( ABS( EigVals_L ) )
                AlphaR = MAXVAL( ABS( EigVals_R ) )

                Alpha    = MAX( AlphaL, AlphaR )
                AlphaMns = MAX( 0.0_DP, &
                     MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )
                AlphaPls = MAX( 0.0_DP, &
                     MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

                AlphaMdl = AlphaC( uCF_L( iNodeX, 1:nCF ),        &
                                   Flux_L( 1:nCF),                &
                                   uCF_R( iNodeX, 1:nCF ),        &
                                   Flux_R( 1:nCF ),               &
                                   AlphaPls,                      &
                                   AlphaMns,                      &
                                   uGF_L( iNodeX, iGF_Gm_dd_11 ), &
                                   uGF_L( iNodeX, iGF_Beta_1 ) )

                Flux &
                  = NumericalFlux_Fluid_GR                        &
                      ( uCF_L( iNodeX, 1:nCF ),                   &
                        uCF_R( iNodeX,1:nCF ),                    &
                        Flux_L( 1:nCF ),                          &
                        Flux_R( 1:nCF ),                          &
                        Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF, &
                        uPF_L( iNodeX, iPF_V1 ),                  &
                        uPF_R( iNodeX, iPF_V1 ),                  &
                        uAF_L( iNodeX, iAF_P ),                   &
                        uAF_R( iNodeX, iAF_P ),                   &
                        uGF_L( iNodeX, iGF_Beta_1 ),              &
                        uGF_L( iNodeX, iGF_Gm_dd_11 ) )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + Flux(1:nCF) / dX1(iX1)

                ! -- Right Face -- 

                ! -- Average Geometry Fields --

                DO iGF = 1, nGF

                  jNodeX = iNodeX
                  uGF_L(iNodeX,iGF) &
                    = 0.5_DP * ( uGF_K(jNodeX,iGF) + uGF_N(jNodeX,iGF) )

                END DO
                uGF_R = uGF_L

                ! -- Left Fluid State --

                DO iCF = 1, nCF

                  jNodeX = iNodeX
                  uCF_L(iNodeX,iCF) = uCF_K(jNodeX,iCF)

                END DO

                uAF_L(iNodeX1,iAF_P)  & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_K(:,iAF_P) )

                uAF_L(iNodeX1,iAF_Gm) & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_K(:,iAF_Gm) )

                CALL ComputePrimitive( uCF_L, uGF_L, uPF_L, uAF_L )
                
                Cs = ComputeSoundSpeed( uAF_L( iNodeX, iAF_P ), &
                                        uPF_L( iNodeX, iPF_E ), &
                                        uPF_L( iNodeX, iPF_D ), &
                                        uAF_L( iNodeX, iAF_Gm ) )

                EigVals_L = Eigenvalues( uPF_L( iNodeX, iPF_V1 ),       &
                                         uPF_L( iNodeX, iPF_V2 ),       &
                                         uPF_L( iNodeX, iPF_V3 ),       &
                                         uGF_L( iNodeX, iGF_Gm_dd_11 ), &
                                         uGF_L( iNodeX, iGF_Gm_dd_22 ), &
                                         uGF_L( iNodeX, iGF_Gm_dd_33 ), &
                                         uGF_L( iNodeX, iGF_Alpha ),    &
                                         uGF_L( iNodeX, iGF_Beta_1 ),   &
                                         Cs )

                Flux_L &
                  = Flux_X1 &
                      ( uPF_L( iNodeX, iPF_D ),        &
                        uPF_L( iNodeX, iPF_V1 ),       &
                        uPF_L( iNodeX, iPF_V2 ),       &
                        uPF_L( iNodeX, iPF_V3 ),       &
                        uPF_L( iNodeX, iPF_E ),        &
                        uAF_L( iNodeX, iAF_P ),        &
                        uPF_L( iNodeX, iPF_Ne ),       &
                        uGF_L( iNodeX, iGF_Alpha ),    &
                        uGF_L( iNodeX, iGF_Beta_1 ),   &
                        uGF_L( iNodeX, iGF_Gm_dd_11 ), &
                        uGF_L( iNodeX, iGF_Gm_dd_22 ), &
                        uGF_L( iNodeX, iGF_Gm_dd_33 ) )

                ! -- Right Fluid State --

                DO iCF = 1, nCF

                  jNodeX = iNodeX
                  uCF_R(iNodeX,iCF) = uCF_N(jNodeX,iCF)

                END DO

                uAF_R(iNodeX1,iAF_P)  & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_N(:,iAF_P) )

                uAF_R(iNodeX1,iAF_Gm) & ! -- Set to Cell-Average
                  = DOT_PRODUCT( WeightsX_q, uAF_N(:,iAF_Gm) )

                CALL ComputePrimitive( uCF_R, uGF_R, uPF_R, uAF_R )
                
                Cs = ComputeSoundSpeed( uAF_R( iNodeX, iAF_P ), &
                                        uPF_R( iNodeX, iPF_E ), &
                                        uPF_R( iNodeX, iPF_D ), &
                                        uAF_R( iNodeX, iAF_Gm ) )

                EigVals_R = Eigenvalues( uPF_R( iNodeX, iPF_V1 ),       &
                                         uPF_R( iNodeX, iPF_V2 ),       &
                                         uPF_R( iNodeX, iPF_V3 ),       &
                                         uGF_R( iNodeX, iGF_Gm_dd_11 ), &
                                         uGF_R( iNodeX, iGF_Gm_dd_22 ), &
                                         uGF_R( iNodeX, iGF_Gm_dd_33 ), &
                                         uGF_R( iNodeX, iGF_Alpha ),    &
                                         uGF_R( iNodeX, iGF_Beta_1 ),   &
                                         Cs )

                Flux_R &
                  = Flux_X1 &
                      ( uPF_R( iNodeX, iPF_D ),        &
                        uPF_R( iNodeX, iPF_V1 ),       &
                        uPF_R( iNodeX, iPF_V2 ),       &
                        uPF_R( iNodeX, iPF_V3 ),       &
                        uPF_R( iNodeX, iPF_E ),        &
                        uAF_R( iNodeX, iAF_P ),        &
                        uPF_R( iNodeX, iPF_Ne ),       &
                        uGF_R( iNodeX, iGF_Alpha ),    &
                        uGF_R( iNodeX, iGF_Beta_1 ),   &
                        uGF_R( iNodeX, iGF_Gm_dd_11 ), &
                        uGF_R( iNodeX, iGF_Gm_dd_22 ), &
                        uGF_R( iNodeX, iGF_Gm_dd_33 ) )

                ! First Alpha will be the largest eigenvalue (absolute value)
                ! of the left and right state
                AlphaL = MAXVAL( ABS( EigVals_L ) )
                AlphaR = MAXVAL( ABS( EigVals_R ) )

                Alpha    = MAX( AlphaL, AlphaR )
                AlphaMns = MAX( 0.0_DP, &
                     MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )
                AlphaPls = MAX( 0.0_DP, &
                     MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

                AlphaMdl = AlphaC( uCF_L( iNodeX, 1:nCF ),        &
                                   Flux_L( 1:nCF),                &
                                   uCF_R( iNodeX, 1:nCF ),        &
                                   Flux_R( 1:nCF ),               &
                                   AlphaPls,                      &
                                   AlphaMns,                      &
                                   uGF_R( iNodeX, iGF_Gm_dd_11 ), &
                                   uGF_L( iNodeX, iGF_Beta_1 ) )

                Flux &
                  = NumericalFlux_Fluid_GR                        &
                      ( uCF_L( iNodeX, 1:nCF ),                   &
                        uCF_R( iNodeX,1:nCF ),                    &
                        Flux_L( 1:nCF ),                          &
                        Flux_R( 1:nCF ),                          &
                        Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF, &
                        uPF_L( iNodeX, iPF_V1 ),                  &
                        uPF_R( iNodeX, iPF_V1 ),                  &
                        uAF_L( iNodeX, iAF_P ),                   &
                        uAF_R( iNodeX, iAF_P ),                   &
                        uGF_R( iNodeX, iGF_Beta_1 ),              &
                        uGF_R( iNodeX, iGF_Gm_dd_11 ) )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      - Flux(1:nCF) / dX1(iX1)

              END DO
            END DO
          END DO

          END ASSOCIATE ! uAF_P, etc.
          END ASSOCIATE ! uPF_P, etc.
          END ASSOCIATE ! uCF_P, etc.
          END ASSOCIATE ! uGF_P, etc.

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeRHS_Euler_DG_GR_X1


END MODULE EulerEquationsSolutionModule_DG_GR

