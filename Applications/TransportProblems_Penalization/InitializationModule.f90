MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Pi, FourPi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV, &
    Kilometer, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm, iAF_Cs
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved
  USE RadiationFieldsModule, ONLY: &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive
  USE MomentEquationsUtilitiesModule, ONLY: &
    Conserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTransportProblem1D

CONTAINS


  SUBROUTINE InitializeTransportProblem1D

    INTEGER             :: iX1, iX2, iX3, iR, iE
    INTEGER             :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER             :: iNodeX, iNode
    REAL(DP)            :: X1, E, Mnu, kT
    REAL(DP), PARAMETER :: MinD = 1.0d08 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: MaxD = 4.0d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: X1_D = 2.0d01 * Kilometer
    REAL(DP), PARAMETER :: H1_D = 1.0d01 * Kilometer
    REAL(DP), PARAMETER :: MinT = 5.0d09 * Kelvin
    REAL(DP), PARAMETER :: MaxT = 2.6d11 * Kelvin
    REAL(DP), PARAMETER :: X1_T = 2.5d01 * Kilometer
    REAL(DP), PARAMETER :: H1_T = 2.0d01 * Kilometer
    REAL(DP), PARAMETER :: MinY = 3.0d-1
    REAL(DP), PARAMETER :: MaxY = 4.6d-1
    REAL(DP), PARAMETER :: X1_Y = 4.5d01 * Kilometer
    REAL(DP), PARAMETER :: H1_Y = 1.0d01 * Kilometer

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! --- Initialize Fluid Fields with Profile ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = 0.5_DP &
                      * ( MaxD*(1.0_DP-TANH((X1-X1_D)/H1_D)) &
                          + MinD*(1.0_DP-TANH((X1_D-X1)/H1_D)) )

                uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                  = 0.5_DP &
                      * ( MaxT*(1.0_DP-TANH((X1-X1_T)/H1_T)) &
                          + MinT*(1.0_DP-TANH((X1_T-X1)/H1_T)) )

                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                  = 0.5_DP &
                      * ( MinY*(1.0_DP-TANH((X1-X1_Y)/H1_Y)) &
                          + MaxY*(1.0_DP-TANH((X1_Y-X1)/H1_Y)) )

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
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
                   uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
                   uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

          uAF(:,iX1,iX2,iX3,iAF_Cs) &
            = SQRT( uAF(:,iX1,iX2,iX3,iAF_Gm) &
                    * uAF(:,iX1,iX2,iX3,iAF_P) &
                    / uPF(:,iX1,iX2,iX3,iPF_D) )

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

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  kT  = BoltzmannConstant &
                        * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

                  Mnu = uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                        + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                        - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = MAX( 1.0d-100, &
                             FourPi / ( EXP( (E-Mnu)/kT ) + 1.0_DP ) )

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

  END SUBROUTINE InitializeTransportProblem1D


END MODULE InitializationModule
