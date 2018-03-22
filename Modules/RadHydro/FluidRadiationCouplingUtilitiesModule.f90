MODULE FluidRadiationCouplingUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, &
    nE, nNodesE
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, MeshE, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, nPF, uAF, nAF
  USE RadiationFieldsModule, ONLY: &
    uPR, nPR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNodes
  PUBLIC :: InitializeNodesX
  PUBLIC :: InitializeWeights
  PUBLIC :: InitializeFluidFields
  PUBLIC :: FinalizeFluidFields
  PUBLIC :: InitializeRadiationFields
  PUBLIC :: FinalizeRadiationFields
  PUBLIC :: MapForward_FluidField
  PUBLIC :: MapBackward_FluidField
  PUBLIC :: MapForward_RadiationField
  PUBLIC :: MapBackward_RadiationField
  PUBLIC :: FermiDirac
  PUBLIC :: dFermiDiracdT
  PUBLIC :: dFermiDiracdY
  PUBLIC :: ENORM

CONTAINS


  SUBROUTINE InitializeNodes( E_N )

    REAL(DP), DIMENSION(:), INTENT(out) :: E_N

    INTEGER :: iE, iNodeE, iNode

    iNode = 0
    DO iE = 1, nE
      DO iNodeE = 1, nNodesE

        iNode = iNode + 1

        E_N(iNode) = NodeCoordinate( MeshE, iE, iNodeE )

      END DO
    END DO

  END SUBROUTINE InitializeNodes


  SUBROUTINE InitializeNodesX( X_N )

    REAL(DP), DIMENSION(:,:), INTENT(out) :: X_N

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX_G

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                X_N(iNodeX_G,1) = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                X_N(iNodeX_G,2) = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
                X_N(iNodeX_G,3) = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeNodesX


  SUBROUTINE InitializeWeights( W2_N, W3_N )

    REAL(DP), DIMENSION(:), INTENT(out) :: W2_N, W3_N

    INTEGER  :: iE, iNodeE, iNode
    REAL(DP) :: E_N

    ASSOCIATE &
      ( dE => MeshE % Width(1:nE) )

    iNode = 0
    DO iE = 1, nE
      DO iNodeE = 1, nNodesE

        iNode = iNode + 1

        E_N = NodeCoordinate( MeshE, iE, iNodeE )

        W2_N(iNode) = dE(iE) * WeightsE(iNodeE) * E_N**2
        W3_N(iNode) = dE(iE) * WeightsE(iNodeE) * E_N**3

      END DO
    END DO

    END ASSOCIATE ! dE

  END SUBROUTINE InitializeWeights


  SUBROUTINE InitializeFluidFields( uPF_N, uAF_N )

    REAL(DP), DIMENSION(:,:), INTENT(out) :: uPF_N, uAF_N

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G
    INTEGER :: iPF, iAF

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                DO iPF = 1, nPF
                  uPF_N(iPF,iNodeX_G) &
                    = uPF(iNodeX,iX1,iX2,iX3,iPF)
                END DO

                DO iAF = 1, nAF
                  uAF_N(iAF,iNodeX_G) &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF)
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFluidFields


  SUBROUTINE FinalizeFluidFields( uPF_N, uAF_N )

    REAL(DP), DIMENSION(:,:), INTENT(in) :: uPF_N, uAF_N

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G
    INTEGER :: iPF, iAF

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                DO iPF = 1, nPF
                  uPF(iNodeX,iX1,iX2,iX3,iPF) &
                    = uPF_N(iPF,iNodeX_G)
                END DO

                DO iAF = 1, nAF
                  uAF(iNodeX,iX1,iX2,iX3,iAF) &
                    = uAF_N(iAF,iNodeX_G)
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE FinalizeFluidFields


  SUBROUTINE InitializeRadiationFields( uPR_N )

    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: uPR_N

    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNode, iNodeE_G
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G
    INTEGER :: iPR

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeE_G = 0
                DO iE = 1, nE
                  DO iNodeE = 1, nNodesE

                    iNodeE_G = iNodeE_G + 1

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    DO iPR = 1, nPR

                      uPR_N(iNodeE_G,iPR,iNodeX_G) &
                        = uPR(iNode,iE,iX1,iX2,iX3,iPR,1)

                    END DO

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeRadiationFields


  SUBROUTINE FinalizeRadiationFields( uPR_N )

    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: uPR_N

    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNode, iNodeE_G
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G
    INTEGER :: iPR

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeE_G = 0
                DO iE = 1, nE
                  DO iNodeE = 1, nNodesE

                    iNodeE_G = iNodeE_G + 1

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    DO iPR = 1, nPR

                      uPR(iNode,iE,iX1,iX2,iX3,iPR,1) &
                        = uPR_N(iNodeE_G,iPR,iNodeX_G)

                    END DO

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE FinalizeRadiationFields


  SUBROUTINE MapForward_FluidField( FF, FF_N )

    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: FF
    REAL(DP), DIMENSION(:),       INTENT(out) :: FF_N

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                FF_N( iNodeX_G ) = FF( iNodeX, iX1, iX2, iX3 )
                    
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapForward_FluidField


  SUBROUTINE MapBackward_FluidField( FF, FF_N )

    REAL(DP), DIMENSION(:,:,:,:), INTENT(out) :: FF
    REAL(DP), DIMENSION(:),       INTENT(in)  :: FF_N

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                FF(iNodeX,iX1,iX2,iX3) = FF_N(iNodeX_G)

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapBackward_FluidField


  SUBROUTINE MapForward_RadiationField( RF, RF_N )

    REAL(DP), DIMENSION(:,:,:,:,:), INTENT(in)  :: RF
    REAL(DP), DIMENSION(:,:),       INTENT(out) :: RF_N

    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNode, iNodeE_G
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeE_G = 0
                DO iE = 1, nE
                  DO iNodeE = 1, nNodesE

                    iNodeE_G = iNodeE_G + 1

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    RF_N(iNodeE_G,iNodeX_G) &
                      = RF(iNode,iE,iX1,iX2,iX3)

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapForward_RadiationField


  SUBROUTINE MapBackward_RadiationField( RF, RF_N )

    REAL(DP), DIMENSION(:,:,:,:,:), INTENT(out) :: RF
    REAL(DP), DIMENSION(:,:),       INTENT(in)  :: RF_N

    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNode, iNodeE_G
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX, iNodeX_G

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX_G = iNodeX_G + 1

                iNodeE_G = 0
                DO iE = 1, nE
                  DO iNodeE = 1, nNodesE

                    iNodeE_G = iNodeE_G + 1

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    RF(iNode,iE,iX1,iX2,iX3) &
                      = RF_N(iNodeE_G,iNodeX_G)

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapBackward_RadiationField


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    REAL(DP) :: Exponent

    Exponent = ( E - Mu ) / kT
    Exponent = MAX( Exponent, - LOG( 1.0d100 ) )
    Exponent = MIN( Exponent, + LOG( 1.0d100 ) )

    FermiDirac &
      = 1.0_DP / ( EXP( Exponent ) + 1.0_DP )

    RETURN
  END FUNCTION FermiDirac


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdT( E, Mu, kT, dMudT, T )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T

    REAL(DP) :: Exponent, FD

    Exponent = ( E - Mu ) / kT
    Exponent = MAX( Exponent, - LOG( 1.0d100 ) )
    Exponent = MIN( Exponent, + LOG( 1.0d100 ) )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdY( E, Mu, kT, dMudY, T )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudY, T

    REAL(DP) :: Exponent, FD

    Exponent = ( E - Mu ) / kT
    Exponent = MAX( Exponent, - LOG( 1.0d100 ) )
    Exponent = MIN( Exponent, + LOG( 1.0d100 ) )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY


  PURE REAL(DP) FUNCTION ENORM( X )

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


END MODULE FluidRadiationCouplingUtilitiesModule
