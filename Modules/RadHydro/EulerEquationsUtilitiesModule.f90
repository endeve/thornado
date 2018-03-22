MODULE EulerEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP, Zero, Half
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumberX_X1
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, WeightsX_X1
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE GravitySolutionModule, ONLY: &
    L_X1_L, L_X1_H, &
    V_X1_L, V_X1_H, &
    dVdX1, dVdX2, dVdX3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Conserved
  PUBLIC :: ComputeConserved
  PUBLIC :: Primitive
  PUBLIC :: ComputePrimitive
  PUBLIC :: Eigenvalues
  PUBLIC :: ComputeEigenvectors_L
  PUBLIC :: ComputeEigenvectors_R
  PUBLIC :: AlphaMax
  PUBLIC :: AlphaP
  PUBLIC :: AlphaM
  PUBLIC :: AlphaC
  PUBLIC :: Flux_X1
  PUBLIC :: GeometrySources
  PUBLIC :: ComputeGeometrySources_Gravity

CONTAINS


  PURE FUNCTION Conserved( Primitive )

    REAL(DP), DIMENSION(1:nPF), INTENT(in) :: Primitive
    REAL(DP), DIMENSION(1:nCF)             :: Conserved

    Conserved(iCF_D)  &
      = Primitive(iPF_D)
    Conserved(iCF_S1) &
      = Primitive(iPF_D) * Primitive(iPF_V1)
    Conserved(iCF_S2) &
      = Primitive(iPF_D) * Primitive(iPF_V2)
    Conserved(iCF_S3) &
      = Primitive(iPF_D) * Primitive(iPF_V3)
    Conserved(iCF_E)  &
      = Primitive(iPF_E) &
          + Half * Primitive(iPF_D) &
              * ( Primitive(iPF_V1)**2 + Primitive(iPF_V2)**2 &
                    + Primitive(iPF_V3)**2 )
    Conserved(iCF_Ne) &
      = Primitive(iPF_Ne)

    RETURN
  END FUNCTION Conserved


  PURE FUNCTION Primitive( Conserved )

    REAL(DP), DIMENSION(1:nCF), INTENT(in) :: Conserved
    REAL(DP), DIMENSION(1:nPF)             :: Primitive

    REAL(DP) :: CF_D

    CF_D = MAX( Conserved(iCF_D), TINY( 1.0_DP ) )

    Primitive(iPF_D)  &
      = CF_D
    Primitive(iPF_V1) &
      = Conserved(iCF_S1) / CF_D
    Primitive(iPF_V2) &
      = Conserved(iCF_S2) / CF_D
    Primitive(iPF_V3) &
      = Conserved(iCF_S3) / CF_D
    Primitive(iPF_E)  &
      = Conserved(iCF_E) &
          - Half * ( Conserved(iCF_S1)**2 + Conserved(iCF_S2)**2 &
                     + Conserved(iCF_S3)**2 ) / CF_D
    Primitive(iPF_Ne)  &
      = Conserved(iCF_Ne)

    RETURN
  END FUNCTION Primitive


  SUBROUTINE ComputeConserved( uPF, uCF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uPF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uCF

    uCF(:,iCF_D) &
      = uPF(:,iPF_D)

    uCF(:,iCF_S1) &
      = uCF(:,iCF_D) * uPF(:,iPF_V1)

    uCF(:,iCF_S2) &
      = uCF(:,iCF_D) * uPF(:,iPF_V2)

    uCF(:,iCF_S3) &
      = uCF(:,iCF_D) * uPF(:,iPF_V3)

    uCF(:,iCF_E) &
      = uPF(:,iPF_E) &
        + Half * uPF(:,iPF_D) &
            * ( uPF(:,iPF_V1)**2 + uPF(:,iPF_V2)**2 + uPF(:,iPF_V3)**2 )

    uCF(:,iCF_Ne) &
      = uPF(:,iPF_Ne)

  END SUBROUTINE ComputeConserved


  SUBROUTINE ComputePrimitive( uCF, uPF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uCF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uPF

    uPF(:,iPF_D) &
      = uCF(:,iCF_D)

    uPF(:,iPF_V1) &
      = uCF(:,iCF_S1) / uCF(:,iCF_D)

    uPF(:,iPF_V2) &
      = uCF(:,iCF_S2) / uCF(:,iCF_D)

    uPF(:,iPF_V3) &
      = uCF(:,iCF_S3) / uCF(:,iCF_D)

    uPF(:,iPF_E) &
      = uCF(:,iCF_E) &
        - Half * ( uCF(:,iCF_S1)**2 + uCF(:,iCF_S2)**2 &
                   + uCF(:,iCF_S3)**2 ) / uCF(:,iCF_D)

    uPF(:,iPF_Ne) &
      = uCF(:,iCF_Ne)

  END SUBROUTINE ComputePrimitive


  PURE FUNCTION Eigenvalues( V, Cs )

    REAL(DP), INTENT(in)     :: V, Cs
    REAL(DP), DIMENSION(nCF) :: Eigenvalues

    Eigenvalues(1:nCF) = [ V - Cs, V, V, V, V, V + Cs ]

    RETURN
  END FUNCTION Eigenvalues


  SUBROUTINE ComputeEigenvectors_L( V1, V2, V3, E, P, Cs, L1, Componentwise )

    REAL(DP),                     INTENT(in)  :: V1, V2, V3, E, P, Cs
    REAL(DP), DIMENSION(nCF,nCF), INTENT(out) :: L1
    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP) :: g, k, h, M1, M2, M3

    IF( Componentwise )THEN

      L1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      L1(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE

      g = P / ( E * Cs )
      k = 0.5_DP * ( V1**2 + V2**2 + V3**2 )
      h = g * k / Cs

      M1 = V1 / Cs; M2 = V2 / Cs; M3 = V3 / Cs

      L1(:,1) &
        = [ 0.5_DP * ( h + M1 ), &
            1.0_DP - h, &
            0.5_DP * ( h - M1 ), &
            V2, - V3, 0.0_DP ]

      L1(:,2) &
        = [ - 0.5_DP * ( g * M1 + 1.0_DP / Cs ), &
                         g * M1, &
            - 0.5_DP * ( g * M1 - 1.0_DP / Cs ), &
              0.0_DP, 0.0_DP, 0.0_DP ]

      L1(:,3) &
        = [ - 0.5_DP * g * M2, &
                       g * M2, &
            - 0.5_DP * g * M2, &
            - 1.0_DP, 0.0_DP, 0.0_DP ]

      L1(:,4) &
        = [ - 0.5_DP * g * M3, &
                       g * M3, &
            - 0.5_DP * g * M3, &
              0.0_DP, 1.0_DP, 0.0_DP ]

      L1(:,5) &
        = [ 0.5_DP * g / Cs, &
                   - g / Cs, &
            0.5_DP * g / Cs, &
            0.0_DP, 0.0_DP, 0.0_DP ]

      L1(:,6) &
        = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    END IF

  END SUBROUTINE ComputeEigenvectors_L


  SUBROUTINE ComputeEigenvectors_R( V1, V2, V3, E, P, Cs, R1, Componentwise )

    REAL(DP),                     INTENT(in)  :: V1, V2, V3, E, P, Cs
    REAL(DP), DIMENSION(nCF,nCF), INTENT(out) :: R1
    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP) :: k, h

    IF( Componentwise )THEN

      R1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      R1(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE

      k = 0.5_DP * ( V1**2 + V2**2 + V3**2 )
      h = Cs**2 * ( E / P ) + k

      R1(:,1) = [ 1.0_DP, V1 - Cs,       V2,     V3, h - Cs * V1, 0.0_DP ]
      R1(:,2) = [ 1.0_DP,      V1,       V2,     V3,           k, 0.0_DP ]
      R1(:,3) = [ 1.0_DP, V1 + Cs,       V2,     V3, h + Cs * V1, 0.0_DP ]
      R1(:,4) = [ 0.0_DP,  0.0_DP, - 1.0_DP, 0.0_DP,        - V2, 0.0_DP ]
      R1(:,5) = [ 0.0_DP,  0.0_DP,   0.0_DP, 1.0_DP,          V3, 0.0_DP ]
      R1(:,6) = [ 0.0_DP,  0.0_DP,   0.0_DP, 0.0_DP,      0.0_DP, 1.0_DP ]

    END IF

  END SUBROUTINE ComputeEigenvectors_R


  PURE REAL(DP) FUNCTION AlphaMax( V, Cs )

    REAL(DP), INTENT(in) :: V, Cs 

    AlphaMax = MAXVAL( ABS( Eigenvalues( V, Cs ) ) )

    RETURN
  END FUNCTION AlphaMax


  PURE REAL(DP) FUNCTION AlphaP( Lambda_L, Lambda_R )

    REAL(DP), DIMENSION(nCF), INTENT(in) :: Lambda_L, Lambda_R

    AlphaP = MAX( Zero, MAXVAL( + Lambda_L ), MAXVAL( + Lambda_R ) )

    RETURN
  END FUNCTION AlphaP


  PURE REAL(DP) FUNCTION AlphaM( Lambda_L, Lambda_R )

    REAL(DP), DIMENSION(nCF), INTENT(in) :: Lambda_L, Lambda_R

    AlphaM = MAX( Zero, MAXVAL( - Lambda_L ), MAXVAL( - Lambda_R ) )

    RETURN
  END FUNCTION AlphaM


  PURE REAL(DP) FUNCTION AlphaC( U_L, U_R, F_L, F_R, aP, aM )

    ! --- Middle Wavespeed as Suggested by Batten et al. (1997) ---
    ! --- (SIAM J. Sci. Comput., Vol. 18, No. 6, pp. 1553-1570) ---

    REAL(DP), DIMENSION(2), INTENT(in) :: U_L, U_R, F_L, F_R
    REAL(DP),               INTENT(in) :: aP, aM

    AlphaC &
      = ( aP * U_R(2) + aM * U_L(2) + ( F_L(2) - F_R(2) ) ) &
          / ( aP * U_R(1) + aM * U_L(1) + ( F_L(1) - F_R(1) ) )

    RETURN
  END FUNCTION AlphaC


  PURE FUNCTION Flux_X1( D, V1, V2, V3, E, P, Ne )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, P, Ne
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1

    Flux_X1(iCF_D)  = D * V1

    Flux_X1(iCF_S1) = D * V1 * V1 + P

    Flux_X1(iCF_S2) = D * V2 * V1

    Flux_X1(iCF_S3) = D * V3 * V1

    Flux_X1(iCF_E)  = ( E + Half * D * ( V1**2 + V2**2 + V3**2 ) + P ) * V1

    Flux_X1(iCF_Ne) = Ne * V1

    RETURN
  END FUNCTION Flux_X1


  FUNCTION GeometrySources( D, S_1, S_2, S_3, P, X )
   
    REAL(DP)             :: GeometrySources(1:nCF)
    REAL(DP), INTENT(in) :: D, S_1, S_2, S_3, P, X(1:3)

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'WARNING: GeometrySources Obsolete'
    WRITE(*,*)
    STOP

    GeometrySources(iCF_D) &
      = 0.0_DP

    GeometrySources(iCF_S1) &
      = 0.0_DP
!!$      = ( S_2**2 / D + P ) * dlnadX1( X ) &
!!$        + ( S_3**2 / D + P ) * dlnbdX1( X )

    GeometrySources(iCF_S2) &
      = 0.0_DP
!!$    (b2**2 * c3**2 * 1/D * S_3**2 + P)/a1 * dlc &
!!$    -(a1 * 1/D * S_1 * S_2) * dla 

    GeometrySources(iCF_S3) &
      = 0.0_DP
!!$      = -(b2 * c3 * 1/D * S_1 * S_3) * dlb &
!!$        -(b2 * c3 * 1/D * S_2 * S_3) * dlc
    
    GeometrySources(iCF_E) &
      = 0.0_DP

    GeometrySources(iCF_Ne) &
      = 0.0_DP

    RETURN
  END FUNCTION GeometrySources


  SUBROUTINE ComputeGeometrySources_Gravity &
               ( dX, uCF, uGF, uGF_P_X1, uGF_N_X1, GS )

    REAL(DP), DIMENSION(3),         INTENT(in)  :: dX
    REAL(DP), DIMENSION(nDOFX,nCF), INTENT(in)  :: uCF
    REAL(DP), DIMENSION(nDOFX,nGF), INTENT(in)  :: uGF, uGF_P_X1, uGF_N_X1
    REAL(DP), DIMENSION(nDOFX,nCF), INTENT(out) :: GS

    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP), DIMENSION(nDOFX) :: g_X1
    REAL(DP), DIMENSION(nNodesX(2)*nNodesX(3)) :: Phi_X1_L, Phi_X1_H

    Phi_X1_L = 0.0_DP
    Phi_X1_H = 0.0_DP
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)

        iNodeX = NodeNumberX_X1( iNodeX2, iNodeX3 )

        Phi_X1_L(iNodeX) &
          = 0.5_DP * ( DOT_PRODUCT &
                           ( L_X1_H(:,iNodeX), uGF_P_X1(:,iGF_Phi_N) ) &
                       + DOT_PRODUCT &
                           ( L_X1_L(:,iNodeX), uGF     (:,iGF_Phi_N) ) )

        Phi_X1_H(iNodeX) &
          = 0.5_DP * ( DOT_PRODUCT &
                           ( L_X1_H(:,iNodeX), uGF     (:,iGF_Phi_N) ) &
                       + DOT_PRODUCT &
                           ( L_X1_L(:,iNodeX), uGF_N_X1(:,iGF_Phi_N) ) )

      END DO
    END DO

    g_X1 = 0.0_DP
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

          ! --- Volume Term ---

          g_X1(iNodeX) &
            = - SUM( WeightsX_q(:) * dVdX1(:,iNodeX) * uGF(:,iGF_Phi_N) ) &
                / ( WeightsX_q(iNodeX) * dX(1) )

          ! --- Surface Terms ---

          g_X1(iNodeX) &
            = g_X1(iNodeX) &
              + ( SUM( WeightsX_X1(:) * Phi_X1_H(:) * V_X1_H(:,iNodeX) ) &
                  - SUM( WeightsX_X1(:) * Phi_X1_L(:) * V_X1_L(:,iNodeX) ) ) &
                / ( WeightsX_q(iNodeX) * dX(1) )

        END DO
      END DO
    END DO

    GS = 0.0_DP
    DO iNodeX = 1, nDOFX

      GS(iNodeX,iCF_S1) &
        = - uCF(iNodeX,iCF_D) * g_X1(iNodeX)

      GS(iNodeX,iCF_S2) &
        = 0.0_DP

      GS(iNodeX,iCF_S3) &
        = 0.0_DP

      GS(iNodeX,iCF_E)  &
        = - uCF(iNodeX,iCF_S1) * g_X1(iNodeX)

    END DO

  END SUBROUTINE ComputeGeometrySources_Gravity


END MODULE EulerEquationsUtilitiesModule
