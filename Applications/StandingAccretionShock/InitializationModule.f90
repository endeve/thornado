MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_S, iAF_Ye, iAF_E, iAF_Gm, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeStandingAccretionShock

CONTAINS


  SUBROUTINE InitializeStandingAccretionShock &
               ( mDot, Mass, rShock, Gamma, Mach )

    REAL(DP), INTENT(in) :: mDot, Mass, rShock, Gamma, Mach

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= rShock )THEN

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = 1.0_DP

                ELSE

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  &
                    = ( mDot / FourPi ) / SQRT( 2.0_DP * Mass ) / X1**1.5_DP 
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                    = - SQRT( 2.0_DP * Mass / X1 )
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                    = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) &
                    = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  &
                    = 1.0_DP

                END IF

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ),  &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ),  &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeStandingAccretionShock


  SUBROUTINE ComputeSettlingSpeed_Bisection( r, alpha, gamma, mass )

    REAL(DP), INTENT(in) :: r, alpha, gamma, mass

  END SUBROUTINE ComputeSettlingSpeed_Bisection


  REAL(DP) FUNCTION SettlingSpeedFun( u, r, alpha, gamma, mass )

    REAL(DP), INTENT(in) :: u, r, alpha, gamma, mass

    SettlingSpeedFun &
      = r * u**2 &
        + alpha * r**(3.0_DP-2.0*gamma) * u**(1.0_DP-gamma) &
        - 2.0_DP * mass

    RETURN
  END FUNCTION SettlingSpeedFun


END MODULE InitializationModule
