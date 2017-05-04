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
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved
  USE RadiationFieldsModule, ONLY: &
    uCR, nCR, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    Conserved
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRelaxation
  PUBLIC :: InitializeRelaxationNES

CONTAINS


  SUBROUTINE InitializeRelaxation( Density, Temperature, ElectronFraction )

    REAL(DP), INTENT(in) :: Density
    REAL(DP), INTENT(in) :: Temperature
    REAL(DP), INTENT(in) :: ElectronFraction

    REAL(DP) :: E, kT, Mnu
    INTEGER  :: iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER  :: iNodeX, iNode

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A13,ES10.4E2,A2,A10,ES10.4E2,A2,A5,ES10.4E2)') &
      '', 'D [g/cm^3] = ', Density / ( Gram / Centimeter**3 ), &
      '', 'T [MeV] = ', Temperature / MeV, &
      '', 'Ye = ', ElectronFraction
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = Density
                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_T)  = Temperature
                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = ElectronFraction

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

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
                      = 1.0d-8

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


  SUBROUTINE InitializeRelaxationNES &
               ( Density, Temperature, ElectronFraction )

    REAL(DP), INTENT(in) :: Density
    REAL(DP), INTENT(in) :: Temperature
    REAL(DP), INTENT(in) :: ElectronFraction

    INTEGER  :: iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER  :: iNodeX, iNode
    REAL(DP) :: E
    REAL(DP), PARAMETER :: E_0 = 1.0d2 * MeV
    REAL(DP), PARAMETER :: L_0 = 3.0d1 * MeV
    REAL(DP), PARAMETER :: D_0 = 0.75_DP

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A13,ES10.4E2,A2,A10,ES10.4E2,A2,A5,ES10.4E2)') &
      '', 'D [g/cm^3] = ', Density / ( Gram / Centimeter**3 ), &
      '', 'T [MeV] = ', Temperature / MeV, &
      '', 'Ye = ', ElectronFraction
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = Density
                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_T)  = Temperature
                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = ElectronFraction

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

    ! --- Initialize Radiation Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = D_0 * EXP( - ( ( E - E_0 ) / L_0 )**2 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,1) &
                      = 0.0_DP

                    uCR(iNode,iE,iX1,iX2,iX3,1:nCR,1) &
                      = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,1) )

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeRelaxationNES


END MODULE FluidRadiationCouplingInitializationModule
