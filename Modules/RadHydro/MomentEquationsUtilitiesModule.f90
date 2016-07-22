MODULE MomentEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nE, nDOF
  USE GeometryModule, ONLY: &
    a, dlnadX1, dlnbdX1, dlncdX2
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Conserved
  PUBLIC :: Primitive
  PUBLIC :: ComputeConserved
  PUBLIC :: ComputePrimitive
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

    EddingtonFactor &
      = 1.0_DP / 3.0_DP &
          + ( 6.0_DP * h**2 - 2.0_DP * h**3 + 6.0_DP * h**4 ) / 15.0_DP

    RETURN
  END FUNCTION EddingtonFactor


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
