MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Three, &
    Pi, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE, &
    nDOF
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable1D3D
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSineWaveDiffusion1D
  PUBLIC :: InitializeSineWaveDamping1D

CONTAINS


  SUBROUTINE InitializeSineWaveDiffusion1D

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: Kappa, X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    Kappa = 1.0d2

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable1D3D(2,iNode)
                X1      = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                  = Three * SQRT( FourPi ) &
                    * ( SIN( Third * Pi * X1 ) + One )

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                  = - ( Pi * SQRT( FourPi ) / ( Three * Kappa ) ) &
                    * COS( Third * Pi * X1 )

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                  = Zero

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                  = Zero

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeSineWaveDiffusion1D


  SUBROUTINE InitializeSineWaveDamping1D

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable1D3D(2,iNode)
                X1      = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                  = One + Half * SIN( TwoPi * X1 )

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                  = uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS)

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                  = Zero

                uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                  = Zero

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeSineWaveDamping1D


END MODULE InitializationModule
