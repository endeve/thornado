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


  SUBROUTINE ApplyApplicationBoundaryConditions_Radiation_X1 &
               ( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'GaussianSphericalWave1D' )

        CALL ApplyBC_Radiation_X1_GaussianSphericalWave1D &
               ( Time, LimiterBC_Option )

      CASE ( 'HomogeneousSphere1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

      CASE ( 'GaussianSphericalDiffusion1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

      CASE ( 'CoolingProblem1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

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


  SUBROUTINE ApplyBC_Radiation_X1_GaussianSphericalWave1D &
               ( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    LOGICAL  :: LimiterBC
    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    REAL(DP) :: X1

    LimiterBC = .FALSE.
    IF( PRESENT( LimiterBC_Option ) ) &
      LimiterBC = LimiterBC_Option

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    ! -- Inner Boundary:

                    IF( LimiterBC )THEN

                      X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

                    ELSE

                      X1 = xL(1)

                    END IF

                    uCR(iNode,iE,0,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( X1 - Time )**2 ) / X1**2

                    uCR(iNode,iE,0,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,0,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G2,iS) &
                      = 0.0_DP

                    uCR(iNode,iE,0,iX2,iX3,iCR_G3,iS) &
                      = 0.0_DP

                    ! -- Outer Boundary:

                    IF( LimiterBC )THEN

                      X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

                    ELSE

                      X1 = xR(1)

                    END IF

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( X1 - Time )**2 ) / X1**2

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


  SUBROUTINE ApplyBC_Radiation_X1_OutflowSphericalSymmetry

    INTEGER :: iS, iX1, iX2, iX3, iE
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    INTEGER :: jNodeX1, jNode

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    ! -- Inner Boundary: Reflecting

                    jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1
                    jNode   = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                    uCR(iNode,iE,0,iX2,iX3,iCR_N,iS) &
                      = uCR(jNode,iE,1,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G1,iS) &
                      = - uCR(jNode,iE,1,iX2,iX3,iCR_G1,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G2,iS) &
                      = + uCR(jNode,iE,1,iX2,iX3,iCR_G2,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G3,iS) &
                      = + uCR(jNode,iE,1,iX2,iX3,iCR_G3,iS)

                    ! -- Outer Boundary: Homogeneous

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_G1,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G2,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_G2,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G3,iS) &
                      = + uCR(iNode,iE,nX(1),iX2,iX3,iCR_G3,iS)

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyBC_Radiation_X1_OutflowSphericalSymmetry


END MODULE ApplicationBoundaryConditionsModule
