MODULE EulerEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    a, b, c, dlnadX1, dlnbdX1, dlncdX2
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne

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
          + 0.5_DP * Primitive(iPF_D) &
              * ( Primitive(iPF_V1)**2 + Primitive(iPF_V2)**2 &
                    + Primitive(iPF_V3)**2 )
    Conserved(iCF_Ne) &
      = Primitive(iPF_Ne)

    RETURN
  END FUNCTION Conserved


  PURE FUNCTION Primitive( Conserved )

    REAL(DP), DIMENSION(1:nCF), INTENT(in) :: Conserved
    REAL(DP), DIMENSION(1:nPF)             :: Primitive

    Primitive(iPF_D)  &
      = Conserved(iCF_D)
    Primitive(iPF_V1) &
      = Conserved(iCF_S1) / Conserved(iCF_D)
    Primitive(iPF_V2) &
      = Conserved(iCF_S2) / Conserved(iCF_D)
    Primitive(iPF_V3) &
      = Conserved(iCF_S3) / Conserved(iCF_D)
    Primitive(iPF_E)  &
      = Conserved(iCF_E) &
          - 0.5_DP * ( Conserved(iCF_S1)**2 + Conserved(iCF_S2)**2 &
                       + Conserved(iCF_S3)**2 ) / Conserved(iCF_D)
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
        + 0.5_DP * uPF(:,iPF_D) &
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
        - 0.5_DP * ( uCF(:,iCF_S1)**2 + uCF(:,iCF_S2)**2 &
                     + uCF(:,iCF_S3)**2 ) / uCF(:,iCF_D)

    uPF(:,iPF_Ne) &
      = uCF(:,iCF_Ne)

  END SUBROUTINE ComputePrimitive


  PURE FUNCTION Eigenvalues( V, Cs )

    REAL(DP), INTENT(in)     :: V, Cs
    REAL(DP), DIMENSION(1:5) :: Eigenvalues

    Eigenvalues(1:5) = [ V - Cs, V, V + Cs, V, V ]

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


  PURE REAL(DP) FUNCTION AlphaP( V_L, Cs_L, V_R, Cs_R )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R

    AlphaP &
      = MAX( 0.0_DP, &
             MAXVAL( + Eigenvalues( V_L, Cs_L ) ), &
             MAXVAL( + Eigenvalues( V_R, Cs_R ) ) )

    RETURN
  END FUNCTION AlphaP


  PURE REAL(DP) FUNCTION AlphaM( V_L, Cs_L, V_R, Cs_R )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R

    AlphaM &
      = MAX( 0.0_DP, &
             MAXVAL( - Eigenvalues( V_L, Cs_L ) ), &
             MAXVAL( - Eigenvalues( V_R, Cs_R ) ) )

    RETURN
  END FUNCTION AlphaM


  PURE REAL(DP) FUNCTION AlphaC( U_L, U_R, F_L, F_R, aP, aM )

    ! --- Middle Wavespeed as Suggested by Batten et al. (1997) ---
    ! --- (SIAM J. Sci. Comput., Vol. 18, No. 6, pp. 1553-1570) ---

    REAL(DP), DIMENSION(2), INTENT(in) :: U_L, U_R, F_L, F_R
    REAL(DP),               INTENT(in) :: aP, aM

    AlphaC &
      = ( aP * U_R(2) + aM * U_L(2) - ( F_R(2) - F_L(2) ) ) &
          / ( aP * U_R(1) + aM * U_L(1) - ( F_R(1) - F_L(1) ) )

    RETURN
  END FUNCTION AlphaC


  PURE FUNCTION Flux_X1( D, V1, V2, V3, E, P, Ne )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, P, Ne
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1

    Flux_X1(iCF_D)  = D * V1

    Flux_X1(iCF_S1) = D * V1 * V1 + P

    Flux_X1(iCF_S2) = D * V2 * V1

    Flux_X1(iCF_S3) = D * V3 * V1

    Flux_X1(iCF_E)  = ( E + 0.5_DP * D * ( V1**2 + V2**2 + V3**2 ) + P ) * V1

    Flux_X1(iCF_Ne) = Ne * V1

    RETURN
  END FUNCTION Flux_X1


  PURE FUNCTION GeometrySources( D, S_1, S_2, S_3, P, X )
   
    REAL(DP)             :: GeometrySources(1:nCF)
    REAL(DP), INTENT(in) :: D, S_1, S_2, S_3, P, X(1:3)

    GeometrySources(iCF_D) &
      = 0.0_DP

    GeometrySources(iCF_S1) &
      = ( S_2**2 / D + P ) * dlnadX1( X ) &
        + ( S_3**2 / D + P ) * dlnbdX1( X )

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


END MODULE EulerEquationsUtilitiesModule
