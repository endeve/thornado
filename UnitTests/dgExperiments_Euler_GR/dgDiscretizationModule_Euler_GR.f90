MODULE dgDiscretizationModule_Euler_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    WeightsX_q, &
    WeightsX_X1, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    LX_X1_Dn, &
    LX_X1_Up
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, iGF_Gm_uu_11, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    rhsCF, uAF, nAF, iAF_P, iAF_Gm
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY:  &
    ComputePrimitive_GR, Eigenvalues_GR, AlphaC_GR, &
    Flux_X1_GR, NumericalFlux_X1_HLLC_GR
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeRHS_Euler_GR

  LOGICAL, PARAMETER :: DisplayTimers = .TRUE.
  REAL(DP) :: Timer_RHS_GR
  REAL(DP) :: Timer_RHS_1_GR, dT_RHS_1_GR
  REAL(DP) :: Timer_RHS_2_GR, dT_RHS_2_GR
  REAL(DP) :: Timer_RHS_3_GR, dT_RHS_3_GR
  REAL(DP) :: Timer_INT_F_GR, dT_INT_F_GR
  REAL(DP) :: Timer_INT_G_GR, dT_INT_G_GR
  REAL(DP) :: Timer_FLX_N_GR, dT_FLX_N_GR

CONTAINS


  SUBROUTINE ComputeRHS_Euler_GR

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iAF, iNodeX, iNodeX_X1, iNodeX1
    REAL(DP) :: Alpha, Beta_1, Beta_2, Beta_3
    REAL(DP) :: ErrorL1, ErrorIn, Error, X1
    REAL(DP), DIMENSION(nDOFX_X1)     :: Cs
    REAL(DP), DIMENSION(nDOFX_X1,nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(nDOFX_X1,nGF) :: uGF_L, uGF_R
    REAL(DP), DIMENSION(nDOFX_X1,nAF) :: uAF_L, uAF_R    
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_L_GR, Flux_X1_R_GR
    REAL(DP), DIMENSION(nDOFX_X1) :: AlphaMns, AlphaMdl, AlphaPls, &
                                         AlphaL, AlphaR, AlphaMAX
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: EigVals_L, EigVals_R, EigVals_GR
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: NumericalFlux_GR
    REAL(DP), DIMENSION(nDOFX,   nPF) :: uPF_P, uPF_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: uCF_P, uCF_K    
    REAL(DP), DIMENSION(nDOFX,   nGF) :: uGF_P, uGF_K
    REAL(DP), DIMENSION(nDOFX,   nAF) :: uAF_P, uAF_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: Flux_X1_q_GR

    Timer_RHS_1_GR = 0.0_DP
    Timer_RHS_2_GR = 0.0_DP
    Timer_RHS_3_GR = 0.0_DP
    Timer_INT_F_GR = 0.0_DP
    Timer_INT_G_GR = 0.0_DP
    Timer_FLX_N_GR = 0.0_DP

    CALL Timer_Start( Timer_RHS_GR )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1) + 1

          DO iCF = 1, nCF

            uCF_P(:,iCF) = uCF(:,iX1-1,iX2,iX3,iCF)
            uCF_K(:,iCF) = uCF(:,iX1,  iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            uGF_P(:,iGF) = uGF(:,iX1-1,iX2,iX3,iGF)
            uGF_K(:,iGF) = uGF(:,iX1,  iX2,iX3,iGF)

          END DO

          DO iAF = 1, nAF

            uAF_P(:,iAF) = uAF(:,iX1-1,iX2,iX3,iAF)
            uAF_K(:,iAF) = uAF(:,iX1,  iX2,iX3,iAF)

          END DO

          !--------------------
          ! --- Volume Term ---
          !--------------------

          CALL ComputePrimitive_GR &
                 ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )


          DO iNodeX = 1, nDOFX

            Flux_X1_q_GR(iNodeX,1:nCF) &
              = Flux_X1_GR &
                  ( uPF_K(iNodeX,iPF_D ), &
                    uPF_K(iNodeX,iPF_V1), &
                    uPF_K(iNodeX,iPF_V2), &
                    uPF_K(iNodeX,iPF_V3), &
                    uPF_K(iNodeX,iPF_E), &                    
                    uAF_K(iNodeX,iAF_P ), &                    
                    uPF_K(iNodeX,iPF_Ne), &
                    uGF_K(iNodeX,iGF_Alpha), &
                    uGF_K(iNodeX,iGF_Beta_1), &                    
                    uGF_K(iNodeX,iGF_Gm_dd_11), &
                    uGF_K(iNodeX,iGF_Gm_dd_22), &
                    uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO

          CALL Timer_Start( dT_RHS_1_GR )

          DO iCF = 1, nCF

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = MATMUL &
                  ( TRANSPOSE( dLXdX1_q(1:nDOFX,1:nDOFX) ), &
                    WeightsX_q(1:nDOFX) * Flux_X1_q_GR(1:nDOFX,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_1_GR )

          CALL Timer_Add( Timer_RHS_1_GR, dT_RHS_1_GR )

          !------------------------
          ! --- Divergence Term ---
          !------------------------

          ! --- Interpolate Fluid Fields ---

          CALL Timer_Start( dT_INT_F_GR )

          DO iCF = 1, nCF

            ! --- Left States ---

            uCF_L(1:nDOFX_X1,iCF) &
              = MATMUL( LX_X1_Up(:,:), uCF_P(:,iCF) )

            ! --- Right States ---

            uCF_R(1:nDOFX_X1,iCF) &
              = MATMUL( LX_X1_Dn(:,:), uCF_K(:,iCF) )

          END DO

          DO iAF = 1, nAF

            ! --- Left States ---

            uAF_L(1:nDOFX_X1,iAF) &
              = MATMUL( LX_X1_Up(:,:), uAF_P(:,iAF) )

            ! --- Right States ---

            uAF_R(1:nDOFX_X1,iAF) &
              = MATMUL( LX_X1_Dn(:,:), uAF_K(:,iAF) )

          END DO

          CALL Timer_Stop( dT_INT_F_GR )

          CALL Timer_Add( Timer_INT_F_GR, dT_INT_F_GR )

          ! --- Interpolate Geometry Fields ---

          CALL Timer_Start( dT_INT_G_GR )

          ! --- Left States ---

          uGF_L(1:nDOFX_X1,iGF_Gm_dd_11) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_11) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_11) ) )
          uGF_L(1:nDOFX_X1,iGF_Gm_dd_22) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_22) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_22) ) )
          uGF_L(1:nDOFX_X1,iGF_Gm_dd_33) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_33) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_33) ) )
          uGF_L(1:nDOFX_X1,iGF_Alpha) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Alpha) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Alpha) ) )
          uGF_L(1:nDOFX_X1,iGF_Beta_1) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Beta_1) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Beta_1) ) )
          uGF_L(1:nDOFX_X1,iGF_Beta_2) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Beta_2) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Beta_2) ) )
          uGF_L(1:nDOFX_X1,iGF_Beta_3) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Beta_3) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Beta_3) ) )

          ! --- Right States (Equal to Left States) ---

          uGF_R(1:nDOFX_X1,iGF_Gm_dd_11) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_11)
          uGF_R(1:nDOFX_X1,iGF_Gm_dd_22) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_22)
          uGF_R(1:nDOFX_X1,iGF_Gm_dd_33) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_33)
          uGF_R(1:nDOFX_X1,iGF_Alpha) &
            = uGF_L(1:nDOFX_X1,iGF_Alpha)
          uGF_R(1:nDOFX_X1,iGF_Beta_1) &
            = uGF_L(1:nDOFX_X1,iGF_Beta_1)
          uGF_R(1:nDOFX_X1,iGF_Beta_2) &
            = uGF_L(1:nDOFX_X1,iGF_Beta_2)
          uGF_R(1:nDOFX_X1,iGF_Beta_3) &
            = uGF_L(1:nDOFX_X1,iGF_Beta_3)

          CALL Timer_Stop( dT_INT_G_GR )

          CALL Timer_Add( Timer_INT_G_GR, dT_INT_G_GR )

          ! --- Left State Primitive, etc. ---

          CALL ComputePrimitive_GR &
                 ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
                   uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
                   uGF_L(:,iGF_Gm_dd_11), &
                   uGF_L(:,iGF_Gm_dd_22), &
                   uGF_L(:,iGF_Gm_dd_33), &
                   uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
                   uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
                   uAF_L(:,iAF_P) )

          CALL ComputeSoundSpeed_GR &
                 ( uAF_L(:,iAF_P), uPF_L(:,iPF_E), uPF_L(:,iPF_D), &
                   uAF_L(:,iAF_Gm) )

          DO iNodeX_X1 = 1, nDOFX_X1

            CALL Eigenvalues_GR( uPF_L( iNodeX_X1, iPF_V1 ),       &
                                   uPF_L( iNodeX_X1, iPF_V2 ),       &
                                   uPF_L( iNodeX_X1, iPF_V3 ),       &
                                   uGF_L( iNodeX_X1, iGF_Gm_dd_11 ), &
                                   uGF_L( iNodeX_X1, iGF_Gm_dd_22 ), &
                                   uGF_L( iNodeX_X1, iGF_Gm_dd_33 ), &
                                   uGF_L( iNodeX_X1, iGF_Gm_uu_11 ), &
                                   uGF_L( iNodeX_X1, iGF_Alpha ),    &
                                   uGF_L( iNodeX_X1, iGF_Beta_1 ),   &
                                   Cs(    iNodeX_X1 ) )



            Flux_X1_L_GR(iNodeX_X1,1:nCF) &
              = Flux_X1_GR &
                  ( uPF_L(iNodeX_X1,iPF_D ), &
                    uPF_L(iNodeX_X1,iPF_V1), &
                    uPF_L(iNodeX_X1,iPF_V2), &
                    uPF_L(iNodeX_X1,iPF_V3), &
                    uPF_L(iNodeX_X1,iPF_E ), &
                    uAF_L(iNodeX_X1,iAF_P ), &
                    uPF_L(iNodeX_X1,iPF_Ne), &
                    uGF_L(iNodeX_X1,iGF_Alpha), &
                    uGF_L(iNodeX_X1,iGF_Beta_1), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_11), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_22), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_33) )

          END DO

          EigVals_L = Eigvals_GR

          ! --- Right State Primitive, etc. ---

          CALL ComputePrimitive_GR &
                 ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                   uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                   uGF_R(:,iGF_Gm_dd_11), &
                   uGF_R(:,iGF_Gm_dd_22), &
                   uGF_R(:,iGF_Gm_dd_33), &
                   uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                   uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                   uAF_R(:,iAF_P) )

          CALL ComputeSoundSpeed_GR &
                 ( uAF_R(:,iAF_P), uPF_R(:,iPF_E), uPF_R(:,iPF_D), &
                   uAF_R(:,iAF_Gm) )

          DO iNodeX_X1 = 1, nDOFX_X1

            CALL Eigenvalues_GR( uPF_R( iNodeX_X1, iPF_V1 ),       &
                                   uPF_R( iNodeX_X1, iPF_V2 ),       &
                                   uPF_R( iNodeX_X1, iPF_V3 ),       &
                                   uGF_R( iNodeX_X1, iGF_Gm_dd_11 ), &
                                   uGF_R( iNodeX_X1, iGF_Gm_dd_22 ), &
                                   uGF_R( iNodeX_X1, iGF_Gm_dd_33 ), &
                                   uGF_R( iNodeX_X1, iGF_Gm_uu_11 ), &
                                   uGF_R( iNodeX_X1, iGF_Alpha ),    &
                                   uGF_R( iNodeX_X1, iGF_Beta_1 ),   &
                                   Cs(    iNodeX_X1 ) )

            Flux_X1_R_GR(iNodeX_X1,1:nCF) &
              = Flux_X1_GR &
                  ( uPF_R(iNodeX_X1,iPF_D ), &
                    uPF_R(iNodeX_X1,iPF_V1), &
                    uPF_R(iNodeX_X1,iPF_V2), &
                    uPF_R(iNodeX_X1,iPF_V3), &
                    uPF_R(iNodeX_X1,iPF_E ), &
                    uAF_L(iNodeX_X1,iAF_P ), &
                    uPF_R(iNodeX_X1,iPF_Ne), &
                    uGF_R(iNodeX_X1,iGF_Alpha), &
                    uGF_R(iNodeX_X1,iGF_Beta_1), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_11), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_22), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_33) )


         END DO
          EigVals_R = EigVals_GR

          ! --- Numerical Flux ---

          CALL Timer_Start( dT_FLX_N_GR )

          DO iNodeX_X1 = 1, nDOFX_X1

                ! First Alpha will be the largest eigenvalue (absolute value)
                ! of the left and right state
            AlphaL = MAXVAL( ABS( EigVals_L( iNodeX_X1, : ) ) )
            AlphaR = MAXVAL( ABS( EigVals_R( iNodeX_X1, : ) ) )

            AlphaMAX  = MAX( AlphaL, AlphaR )
            AlphaMns  = MAX( 0.0_DP, &
                 MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )
            AlphaPls = MAX( 0.0_DP, &
                 MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

            AlphaMdl = AlphaC_GR( uCF_L( iNodeX_X1, 1:nCF ),      &
                               uCF_R( iNodeX_X1, 1:nCF ),      &
                               Flux_X1_L_GR( iNodeX_X1,1:nCF),           &
                               Flux_X1_R_GR( iNodeX_X1,1:nCF ),          &
                               AlphaPls(iNodeX_X1),                       &
                               AlphaMns(iNodeX_X1),                       &
                               uPF_L( iNodeX_X1, iPF_V1 ),     &
                               uPF_R( iNodeX_X1, iPF_V1 ),     &
                               uGF_L( iNodeX_X1, iGF_Beta_1 ), &
                               uGF_L( iNodeX_X1, iGF_Gm_dd_11 ) )

            NumericalFlux_GR( iNodeX_X1,:) &
              = NumericalFlux_X1_HLLC_GR                      &
                  ( uCF_L( iNodeX_X1, 1:nCF ),                &
                    uCF_R( iNodeX_X1,1:nCF ),                 &
                    Flux_X1_L_GR( iNodeX_X1,1:nCF ),                    &
                    Flux_X1_R_GR( iNodeX_X1,1:nCF ),                    &
                    AlphaMAX(iNodeX_X1), AlphaPls(iNodeX_X1), AlphaMns(iNodeX_X1), AlphaMdl(iNodeX_X1), nCF, &
                    uPF_L( iNodeX_X1, iPF_V1 ),               &
                    uPF_R( iNodeX_X1, iPF_V1 ),               &
                    uAF_L( iNodeX_X1, iAF_P ),                &
                    uAF_R( iNodeX_X1, iAF_P ),                &
                    uGF_R( iNodeX_X1, iGF_Beta_1 ),           &
                    uGF_R( iNodeX_X1, iGF_Gm_uu_11 ),         &
                    uGF_R( iNodeX_X1, iGF_Gm_dd_11 ) )


          END DO

          CALL Timer_Stop( dT_FLX_N_GR )

          CALL Timer_Add( Timer_FLX_N_GR, dT_FLX_N_GR )

          ! --- Contribution to This Element ---

          CALL Timer_Start( dT_RHS_2_GR )

          DO iCF = 1, nCF

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = rhsCF(:,iX1,iX2,iX3,iCF) &
                  + MATMUL &
                      ( TRANSPOSE( LX_X1_Dn(:,:) ), &
                        WeightsX_X1(:) * NumericalFlux_GR(:,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_2_GR )

          CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

          ! --- Contribution to Previous Element ---

          CALL Timer_Start( dT_RHS_3_GR )

          DO iCF = 1, nCF

            rhsCF(:,iX1-1,iX2,iX3,iCF) &
              = rhsCF(:,iX1-1,iX2,iX3,iCF) &
                  - MATMUL &
                      ( TRANSPOSE( LX_X1_Up(:,:) ), &
                        WeightsX_X1(:) * NumericalFlux_GR(:,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_3_GR )

          CALL Timer_Add( Timer_RHS_3_GR, dT_RHS_3_GR )

        END DO
      END DO
    END DO

    ! --- Multiply Inverse Mass Matrix ---

    DO iCF = 1, nCF
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = rhsCF(:,iX1,iX2,iX3,iCF) &
                  / ( WeightsX_q(:) * MeshX(1) % Width(iX1) )

          END DO
        END DO
      END DO
    END DO

    CALL Timer_Stop( Timer_RHS_GR )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A24,ES10.4E2)') &
        '', 'ComputeRHS_Euler: ', Timer_RHS_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 1: ', Timer_RHS_1_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 2: ', Timer_RHS_2_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 3: ', Timer_RHS_3_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT F: ', Timer_INT_F_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT G: ', Timer_INT_G_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'FLX N: ', Timer_FLX_N_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
           '', 'Sum: ', Timer_RHS_1_GR+Timer_RHS_2_GR+Timer_RHS_3_GR+&
           Timer_INT_F_GR &
        + Timer_INT_G_GR + Timer_FLX_N_GR

    END IF

    ! --- Compute Error ---

    ErrorL1 = 0.0_DP
    ErrorIn = 0.0_DP

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            Error &
              = ABS( - Pi * COS( TwoPi * X1 ) &
                     - rhsCF(iNodeX,iX1,iX2,iX3,1) )

            ErrorL1 = ErrorL1 + Error
            ErrorIn = MAX( ErrorIn, Error )

          END DO

        END DO
      END DO
    END DO

    ErrorL1 = ErrorL1 / REAL( nDOFX*nX(1)*nX(2)*nX(3) )

    WRITE(*,*)
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorL1: ', ErrorL1
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorIn: ', ErrorIn
    WRITE(*,*)

  END SUBROUTINE ComputeRHS_Euler_GR


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE dgDiscretizationModule_Euler_GR
