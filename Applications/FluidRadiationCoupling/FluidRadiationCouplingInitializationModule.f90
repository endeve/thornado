MODULE FluidRadiationCouplingInitializationModule

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    MeV, &
    Kelvin
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumber, &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_T, iAF_Ye, iAF_E, iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRelaxation

CONTAINS


  SUBROUTINE InitializeRelaxation( D, T, Ye, E_0 )

    REAL(DP), INTENT(in) :: D, T, Ye, E_0

    REAL(DP) :: E, kT, Mnu
    INTEGER  :: iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER  :: iNodeX, iNode

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A4,ES10.4E2,A2,A4,ES10.4E2,A2,A5,ES10.4E2,A2,A6,ES10.4E2)') &
      '', 'D = ', D / ( Gram / Centimeter**3 ), &
      '', 'T = ', T / Kelvin, &
      '', 'Ye = ', Ye, &
      '', 'E_0 = ', E_0 / MeV
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D
                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_T)  = T
                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Ye

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

        END DO
      END DO
    END DO

    CALL ApplyEquationOfState

    ! --- Initialize Radiation Fields ---

    ASSOCIATE( kB => BoltzmannConstant )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  kT  = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

                  Mnu = uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                        + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                        - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = 4.0_DP * Pi &
                          * EXP( - 0.5_DP * ( E - E_0 )**2 / ( kB * T )**2 )

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

    END ASSOCIATE ! kB

  END SUBROUTINE InitializeRelaxation


END MODULE FluidRadiationCouplingInitializationModule
