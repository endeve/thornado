MODULE ApplicationBoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    xL, xR, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyApplicationBoundaryConditions_Radiation_X1

CONTAINS


  SUBROUTINE ApplyApplicationBoundaryConditions_Radiation_X1( Time )

    REAL(DP), INTENT(in) :: Time

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'GaussianSphericalWave1D' )

        CALL ApplyBC_Radiation_X1_GaussianSphericalWave1D( Time )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A38)') &
          '', 'Application Boundary Condition Missing'
        WRITE(*,*)
        WRITE(*,'(A7,A)') &
          '', TRIM( ProgramName )
        WRITE(*,*)
        STOP

    END SELECT

  END SUBROUTINE ApplyApplicationBoundaryConditions_Radiation_X1


  SUBROUTINE ApplyBC_Radiation_X1_GaussianSphericalWave1D( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    ! -- Inner Boundary:

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uCR(iNode,iE,0,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( xL(1) - Time )**2 ) / xL(1)**2

                    uCR(iNode,iE,0,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,0,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G2,iS) &
                      = 0.0_DP

                    uCR(iNode,iE,0,iX2,iX3,iCR_G3,iS) &
                      = 0.0_DP

                    ! -- Outer Boundary:

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( xR(1) - Time )**2 ) / xR(1)**2

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G2,iS) &
                      = 0.0_DP

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G3,iS) &
                      = 0.0_DP

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyBC_Radiation_X1_GaussianSphericalWave1D


END MODULE ApplicationBoundaryConditionsModule
