MODULE MomentEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP, Zero, Fifth, Third, Half, &
    One, Two, Three, Four, Five
  USE ProgramHeaderModule, ONLY: &
    nE, nDOF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Closure_Minerbo = .TRUE.

  PUBLIC :: Conserved
  PUBLIC :: Primitive
  PUBLIC :: ComputeConservedMoments
  PUBLIC :: ComputePrimitiveMoments
  PUBLIC :: Eigenvalues
  PUBLIC :: ComputeEigenvectors_L
  PUBLIC :: ComputeEigenvectors_R
  PUBLIC :: AlphaMax
  PUBLIC :: AlphaP
  PUBLIC :: AlphaM
  PUBLIC :: Flux_X1
  PUBLIC :: GeometrySources
  PUBLIC :: FluxFactor
  PUBLIC :: EddingtonFactor

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


  SUBROUTINE ComputeConservedMoments( iX_Begin, iX_End )

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

  END SUBROUTINE ComputeConservedMoments


  SUBROUTINE ComputePrimitiveMoments( iX_Begin, iX_End )

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

  END SUBROUTINE ComputePrimitiveMoments


  PURE FUNCTION Eigenvalues( N, G_1, G_2, G_3, FF, EF )

    REAL(DP), INTENT(in)     :: N, G_1, G_2, G_3, FF, EF
    REAL(DP), DIMENSION(1:4) :: Eigenvalues

    REAL(DP) :: G1, dEF, Xi, dXi, D, sqrtD

    dEF = EddingtonFactorDerivative( FF )

    sqrtD &
      = SQRT( MAX( ( dEF - Two * FF )**2 + Four * ( EF - FF**2 ), Zero ) )

    G1 = G_1 / ( MAX( N * FF, TINY( One ) ) )

    Eigenvalues(1) = Half * ( G1 * dEF + sqrtD )
    Eigenvalues(2) = Half * ( G1 * dEF - sqrtD )
    Eigenvalues(3) = Half * ( Three * EF - One ) / EF
    Eigenvalues(4) = Half * ( Three * EF - One ) / EF

    RETURN
  END FUNCTION Eigenvalues


  SUBROUTINE ComputeEigenvectors_L( N, G_1, G_2, G_3, L1, Componentwise )

    REAL(DP),                     INTENT(in)  :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(nCR,nCR), INTENT(out) :: L1
    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP)                 :: FF, EF, dLambda
    REAL(DP), DIMENSION(1:4) :: Lambda

    FF = FluxFactor( N, G_1, G_2, G_3 )

    EF = EddingtonFactor( FF )

    Lambda = Eigenvalues( N, G_1, G_2, G_3, FF, EF )

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

    REAL(DP)                 :: FF, EF
    REAL(DP), DIMENSION(1:4) :: Lambda

    FF = FluxFactor( N, G_1, G_2, G_3 )

    EF = EddingtonFactor( FF )

    Lambda = Eigenvalues( N, G_1, G_2, G_3, FF, EF )

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


  PURE REAL(DP) FUNCTION AlphaMax( Lambda )

    REAL(DP), DIMENSION(4), INTENT(in) :: Lambda

    AlphaMax = MAXVAL( ABS( Lambda ) )

    RETURN
  END FUNCTION AlphaMax


  PURE REAL(DP) FUNCTION AlphaP( Lambda_L, Lambda_R )

    REAL(DP), DIMENSION(4), INTENT(in) :: Lambda_L, Lambda_R

    AlphaP = MAX( Zero, MAXVAL( + Lambda_L ), MAXVAL( + Lambda_R ) )

    RETURN
  END FUNCTION AlphaP


  PURE REAL(DP) FUNCTION AlphaM( Lambda_L, Lambda_R )

    REAL(DP), DIMENSION(4), INTENT(in) :: Lambda_L, Lambda_R

    AlphaM = MAX( Zero, MAXVAL( - Lambda_L ), MAXVAL( - Lambda_R ) )

    RETURN
  END FUNCTION AlphaM


  PURE FUNCTION Flux_X1( N, G_1, G_2, G_3, FF, EF )

    REAL(DP)             :: Flux_X1(1:4)
    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3, FF, EF

    REAL(DP) :: G2

    G2 = MAX( ( FF * N )**2, TINY( One ) )

    Flux_X1(1) &
      = G_1

    Flux_X1(2) &
      = N * Half * ( (Three*EF - One)*G_1*G_1/G2 + (One - EF) )

    Flux_X1(3) &
      = N * Half * ( (Three*EF - One)*G_2*G_1/G2 )

    Flux_X1(4) &
      = N * Half * ( (Three*EF - One)*G_3*G_1/G2 )

    RETURN
  END FUNCTION Flux_X1


  FUNCTION GeometrySources( N, G_1, G_2, G_3, X )

    REAL(DP)             :: GeometrySources(1:4)
    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3, X(1:3)

    REAL(DP) :: FF, EF, G2

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'WARNING: GeometrySources Obsolete'
    WRITE(*,*)
    STOP

    FF = FluxFactor( N, G_1, G_2, G_3 )

    EF = EddingtonFactor( FF )

    G2 = MAX( ( FF * N )**2, TINY( One ) )

    GeometrySources(1) &
      = 0.0_DP

    GeometrySources(2) &
      = 0.0_DP
!!$      = N * Half * ( (Three*EF - One)*G_2*G_2/G2 + (One - EF) ) &
!!$          * dlnadX1( X ) &
!!$        + N * Half * ( (Three*EF - One)*G_3*G_3/G2 + (One - EF) ) &
!!$            * dlnbdX1( X )

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


  PURE REAL(DP) FUNCTION FluxFactor( N, G_1, G_2, G_3 )

    REAL(DP), INTENT(in) :: N, G_1, G_2, G_3

    FluxFactor &
      = MAX( Zero, MIN( SQRT( G_1**2 + G_2**2 + G_3**2 ) / N, One ) )
      
    RETURN
  END FUNCTION FluxFactor


  PURE REAL(DP) FUNCTION EddingtonFactor( FF )

    REAL(DP), INTENT(in) :: FF

    IF( Closure_Minerbo )THEN

      ! Minerbo:
      EddingtonFactor &
        = Third + Two * FF**2 * ( One - Third * FF + FF**2 ) * Fifth

    ELSE

      ! Levermore:
      eddingtonFactor &
        = ( Five - Two * SQRT( Four - Three * FF**2 ) ) * Third

    END IF

    RETURN
  END FUNCTION EddingtonFactor


  PURE REAL(DP) FUNCTION EddingtonFactorDerivative( FF )

    REAL(DP), INTENT(in) :: FF

    IF( Closure_Minerbo )THEN

      ! Minerbo:
      EddingtonFactorDerivative &
        = Four * FF * ( One - Half * FF + Two * FF**2 ) * Fifth

    ELSE

      ! Levermore:
      EddingtonFactorDerivative &
        = Two * FF / SQRT( Four - Three * FF**2 )

    END IF

    RETURN
  END FUNCTION EddingtonFactorDerivative


END MODULE MomentEquationsUtilitiesModule
