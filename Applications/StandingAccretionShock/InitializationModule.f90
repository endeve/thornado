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
    REAL(DP) :: X1, Alpha, D_prime, V1_prime, P_prime

    Alpha = 4 * Gamma/((Gamma + 1)*(Gamma - 1))*((Gamma - 1)/(Gamma + 1))**Gamma

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

                  CALL ComputeSettlingSpeed_Bisection(X1, Alpha, Gamma, Mass, uPF(iNode,iX1,iX2,iX3,iPF_V1))                 
                   
                  uPF(iNode,iX1,iX2,iX3,iPF_D)  &
                    = (mDot/FourPi) * (-uPF(iNode,iX1,iX2,iX3,iPF_V1))**(-1.0_DP) * X1**(-2.0_DP)
                  
                  V1_prime &
                    = (gamma - 1_DP)/(gamma + 1.0_DP) * SQRT(2.0_DP * Mass/rShock)
                  D_prime  &
                    = (mDot/FourPi) * (1.0_DP/V1_prime) * rShock**(-2.0_DP)
                  P_prime  &
                    = 2/(gamma + 1.0_DP) * (mDot/FourPi) * SQRT(2 * Mass) * rShock**(-5.0_DP/2.0_DP)
                  
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  &
                    = P_prime * (uPF(iNode,iX1,iX2,iX3,iPF_D)/D_prime)**gamma

                ELSE

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  &
                    = ( mDot / FourPi ) / SQRT( 2.0_DP * Mass ) / (X1**1.5_DP) 
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                    = - SQRT( 2.0_DP * Mass / X1 )
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                    = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) &
                    = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  &
                    = uPF(iNode,iX1,iX2,iX3,iPF_D)/gamma * (-uPF(iNode,iX1,iX2,iX3,iPF_V1)/Mach)**2.0_DP

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


  SUBROUTINE  ComputeSettlingSpeed_Bisection( r, alpha, gamma, mass, V1 )

    REAL(DP), INTENT(in) :: r, alpha, gamma, mass
    
    LOGICAL :: Converged
    INTEGER :: Iter
    REAL(DP) :: a, b, c, ab, F_a, F_b, F_c, F_0
    INTEGER, PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol_ab = 1.0d-8
    REAL(DP), PARAMETER :: Tol_F = 1.0d-8

    REAL(DP), INTENT(out) :: V1

    a = 1.0d-6
    F_a = SettlingSpeedFun(a, r, alpha, gamma, mass)

    b = 1.0_DP
    F_b = SettlingSpeedFun(b, r, alpha, gamma, mass)

    F_0 = F_a
    ab = b - a

    Converged = .FALSE.
    Iter = 0
   
    DO WHILE ( .NOT. Converged)

      Iter = Iter + 1

      ab = 0.5_DP * ab
      c = a + ab
      
      F_c = SettlingSpeedFun(c, r, alpha, gamma, mass)

      IF( F_a * F_c < 0.0_DP ) THEN
     
        b   = c
        F_b = F_c

      ELSE
 
        a   = c
        F_a = F_c

      END IF

      IF( ab < Tol_ab .AND. ABS( F_a ) / F_0 < Tol_F) Converged = .TRUE.

    END DO

    V1 = -a

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
