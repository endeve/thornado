MODULE MomentEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nE, nDOF
  USE GeometryFieldsModule, ONLY: &
    a, dlnadX1, dlnbdX1, dlncdX2
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Closure_Minerbo = .TRUE.

  PUBLIC :: Conserved
  PUBLIC :: Primitive
  PUBLIC :: ComputeConserved
  PUBLIC :: ComputePrimitive
  PUBLIC :: ComputeEigenvectors_L
  PUBLIC :: ComputeEigenvectors_R
  PUBLIC :: AlphaMax
  PUBLIC :: AlphaP
  PUBLIC :: AlphaM
  PUBLIC :: Flux_X1
  PUBLIC :: GeometrySources

CONTAINS


  PURE FUNCTION Conserved( Primitive )

    REAL(DP), DIMENSION(1:nPR), INTENT(in) :: Primitive
    REAL(DP), DIMENSION(1:nCR)             :: Conserved

    Conserved(iCR_N)  = Primitive(iPR_D)
    Conserved(iCR_G1) = Primitive(iPR_I1)
    Conserved(iCR_G2) = Primitive(iPR_I2)
    Conserved(iCR_G3) = Primitive(iPR_I3)

    RETURN
  END FUNCTION Conserved


  PURE FUNCTION Primitive( Conserved )

    REAL(DP), DIMENSION(1:nPR), INTENT(in) :: Conserved
    REAL(DP), DIMENSION(1:nCR)             :: Primitive

    Primitive(iPR_D)  = Conserved(iCR_N)
    Primitive(iPR_I1) = Conserved(iCR_G1)
    Primitive(iPR_I2) = Conserved(iCR_G2)
    Primitive(iPR_I3) = Conserved(iCR_G3)

    RETURN
  END FUNCTION Primitive


  SUBROUTINE ComputeConserved( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iNode, iE, iX1, iX2, iX3, iS

    DO iS = 1, nSpecies

      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)
            DO iE = 1, nE
              DO iNode = 1, nDOF

                uCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                  = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,iS) )

              END DO
            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputeConserved


  SUBROUTINE ComputePrimitive( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iNode, iE, iX1, iX2, iX3, iS

    DO iS = 1, nSpecies

      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)
            DO iE = 1, nE
              DO iNode = 1, nDOF

                uPR(iNode,iE,iX1,iX2,iX3,1:nPR,iS) &
                  = Primitive( uCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) )

              END DO
            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputePrimitive


  PURE FUNCTION Eigenvalues( N, G_1, G_2, G_3 )

    REAL(DP), INTENT(in)     :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(1:4) :: Eigenvalues

    REAL(DP) :: h, h1, Xi, dXi, D, sqrtD

    h  = Length( ReducedFlux( N, [ G_1, G_2, G_3 ] ) )
    h1 = G_1 / Length( [ G_1, G_2, G_3 ] )

    Xi  = EddingtonFactor( h )
    dXi = EddingtonFactorDerivative( h )

    D = MAX( ( dXi - 2.0_DP * h )**2 + 4.0_DP * ( Xi - h**2 ), 0.0_DP )
    sqrtD = SQRT( D )

    Eigenvalues(1:4) &
      = [ 0.5_DP * ( h1 * dXi + sqrtD ), &
          0.5_DP * ( h1 * dXi - sqrtD ), &
          0.5_DP * ( 3.0 * Xi - 1.0_DP ) / h, &
          0.5_DP * ( 3.0 * Xi - 1.0_DP ) / h ]

    RETURN
  END FUNCTION Eigenvalues


  SUBROUTINE ComputeEigenvectors_L( N, G_1, G_2, G_3, L1, Componentwise )

    REAL(DP),                     INTENT(in)  :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(nCR,nCR), INTENT(out) :: L1
    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP)                 :: dLambda
    REAL(DP), DIMENSION(1:4) :: Lambda

    Lambda = Eigenvalues( N, G_1, G_2, G_3 )

    IF( ABS( Lambda(1) - Lambda(2) ) < 1.0d-8 .OR. Componentwise )THEN

      L1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      L1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      L1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE

      dLambda = Lambda(2) - Lambda(1)

      L1(:,1) = [ + Lambda(2) / dLambda, - 1.0_DP / dLambda, 0.0_DP, 0.0_DP ]
      L1(:,2) = [ - Lambda(1) / dLambda, + 1.0_DP / dLambda, 0.0_DP, 0.0_DP ]
      L1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      L1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    END IF

  END SUBROUTINE ComputeEigenvectors_L


  SUBROUTINE ComputeEigenvectors_R( N, G_1, G_2, G_3, R1, Componentwise )

    REAL(DP),                     INTENT(in)  :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(nCR,nCR), INTENT(out) :: R1
    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP), DIMENSION(1:4) :: Lambda

    Lambda = Eigenvalues( N, G_1, G_2, G_3 )

    IF( ABS( Lambda(1) - Lambda(2) ) < 1.0d-8 .OR. Componentwise )THEN

      R1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      R1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE

      R1(:,1) = [ 1.0_DP,    1.0_DP,    0.0_DP, 0.0_DP ]
      R1(:,2) = [ Lambda(1), Lambda(2), 0.0_DP, 0.0_DP ]
      R1(:,3) = [ 0.0_DP,    0.0_DP,    1.0_DP, 0.0_DP ]
      R1(:,4) = [ 0.0_DP,    0.0_DP,    0.0_DP, 1.0_DP ]

    END IF

  END SUBROUTINE ComputeEigenvectors_R


  PURE REAL(DP) FUNCTION AlphaMax( N, G_1, G_2, G_3 )

    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3

    AlphaMax = MAXVAL( ABS( Eigenvalues( N, G_1, G_2, G_3 ) ) )

    RETURN
  END FUNCTION AlphaMax


  PURE REAL(DP) FUNCTION AlphaP &
    ( N_L, G_1_L, G_2_L, G_3_L, N_R, G_1_R, G_2_R, G_3_R )

    REAL(DP), INTENT(in) :: N_L, G_1_L, G_2_L, G_3_L
    REAL(DP), INTENT(in) :: N_R, G_1_R, G_2_R, G_3_R

    AlphaP &
      = MAX( 0.0_DP, &
             MAXVAL( + Eigenvalues( N_L, G_1_L, G_2_L, G_3_L ) ), &
             MAXVAL( + Eigenvalues( N_R, G_1_R, G_2_R, G_3_R ) ) )

    RETURN
  END FUNCTION AlphaP


  PURE REAL(DP) FUNCTION AlphaM &
    ( N_L, G_1_L, G_2_L, G_3_L, N_R, G_1_R, G_2_R, G_3_R )

    REAL(DP), INTENT(in) :: N_L, G_1_L, G_2_L, G_3_L
    REAL(DP), INTENT(in) :: N_R, G_1_R, G_2_R, G_3_R

    AlphaM &
      = MAX( 0.0_DP, &
             MAXVAL( - Eigenvalues( N_L, G_1_L, G_2_L, G_3_L ) ), &
             MAXVAL( - Eigenvalues( N_R, G_1_R, G_2_R, G_3_R ) ) )

    RETURN
  END FUNCTION AlphaM


  PURE FUNCTION Flux_X1( N, G_1, G_2, G_3 )

    REAL(DP)             :: Flux_X1(1:4)
    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3

    REAL(DP) :: Xi, G2

    Xi = EddingtonFactor &
           ( Length( ReducedFlux( N, [ G_1, G_2, G_3 ] ) ) )

    G2 = Length( [ G_1, G_2, G_3 ] )**2

    Flux_X1(1) &
      = G_1

    Flux_X1(2) &
      = N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_1*G_1/G2 + (1.0_DP - Xi) )

    Flux_X1(3) &
      = N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_2*G_1/G2 )

    Flux_X1(4) &
      = N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_3*G_1/G2 )

    RETURN
  END FUNCTION Flux_X1


  PURE FUNCTION GeometrySources( N, G_1, G_2, G_3, X )

    REAL(DP)             :: GeometrySources(1:4)
    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3, X(1:3)

    REAL(DP) :: Xi, G2

    Xi = EddingtonFactor &
           ( Length( ReducedFlux( N, [ G_1, G_2, G_3 ] ) ) )

    G2 = Length( [ G_1, G_2, G_3 ] )**2

    GeometrySources(1) &
      = 0.0_DP

    GeometrySources(2) &
      = N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_2*G_2/G2 + (1.0_DP - Xi) ) &
          * dlnadX1( X ) &
        + N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_3*G_3/G2 + (1.0_DP - Xi) ) &
            * dlnbdX1( X )

    GeometrySources(3) &
      = 0.0_DP
!!$        N * 0.5_DP * ( (3.0_DP*Xi - 1.0_DP)*G_3*G_3/G2 + (1.0_DP - Xi) ) &
!!$            * dlncdX2( X ) / a( X ) &
!!$        - N * 0.5_DP * (3.0_DP*Xi - 1.0_DP)*G_1*G_2/G2 &
!!$            * dlnadX1( X )

    GeometrySources(4) &
      = 0.0_DP
!!$        - N * 0.5_DP * (3.0_DP*Xi - 1.0_DP)*G_1*G_3/G2 &
!!$            * dlnbdX1( X ) &
!!$        - N * 0.5_DP * (3.0_DP*Xi - 1.0_DP)*G_2*G_3/G2 &
!!$            * dlncdX2( X ) / a( X )

    RETURN
  END FUNCTION GeometrySources


  PURE REAL(DP) FUNCTION EddingtonFactor( h )

    REAL(DP), INTENT(in) :: h

    IF( Closure_Minerbo )THEN

      ! Minerbo:
      EddingtonFactor &
        = 1.0_DP / 3.0_DP &
            + 2.0_DP * h**2 * ( 3.0_DP - h + 3.0_DP * h**2 ) / 15.0_DP

    ELSE

      ! Levermore:
      eddingtonFactor &
        = ( 5.0_DP - 2.0_DP * SQRT( 4.0_DP - 3.0_DP * h**2 ) ) / 3.0_DP

    END IF

    RETURN
  END FUNCTION EddingtonFactor


  PURE REAL(DP) FUNCTION EddingtonFactorDerivative( h )

    REAL(DP), INTENT(in) :: h

    IF( Closure_Minerbo )THEN

      ! Minerbo:
      EddingtonFactorDerivative &
        = 2.0_DP * h * ( 2.0_DP - h + 4.0_DP * h**2 ) / 5.0_DP

    ELSE

      ! Levermore:
      EddingtonFactorDerivative &
        = 2.0_DP * h / SQRT( 4.0_DP - 3.0_DP * h**2 )

    END IF

    RETURN
  END FUNCTION EddingtonFactorDerivative


  PURE FUNCTION ReducedFlux( N, G )

    REAL(DP)             :: ReducedFlux(1:3)
    REAL(DP), INTENT(in) :: N, G(1:3)

    ReducedFlux(1:3) = G(1:3) / N

    RETURN
  END FUNCTION ReducedFlux


  PURE REAL(DP) FUNCTION Length( V )

    REAL(DP), INTENT(in) :: V(1:3)

    Length = SQRT( MAX( DOT_PRODUCT( V, V ), TINY( 1.0_DP ) ) )

    RETURN
  END FUNCTION Length


END MODULE MomentEquationsUtilitiesModule
