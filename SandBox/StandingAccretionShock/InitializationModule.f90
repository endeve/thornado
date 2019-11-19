MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B1, iX_E1, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields, ApplyPerturbations

CONTAINS

  SUBROUTINE InitializeFields( mDot, Mass, rShock, Gamma, Mach )

    REAL(DP), INTENT(in) :: mDot, Mass, rShock, Gamma, Mach

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, Alpha, Speed, D_prime, V1_prime, P_prime

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A10,ES12.6E2)') &
      '', 'mDot = ', mDot
    WRITE(*,'(A4,A10,ES12.6E2)') &
      '', 'Mass = ', Mass
    WRITE(*,'(A4,A10,ES12.6E2)') &
      '', 'rShock = ', rShock
    WRITE(*,'(A4,A10,ES12.6E2)') &
      '', 'Gamma = ', Gamma
    WRITE(*,'(A4,A10,ES12.6E2)') &
      '', 'Mach = ', Mach
    WRITE(*,*)

    Alpha = 4.0_DP * Gamma / ( (Gamma + 1.0_DP) * (Gamma - 1.0_DP) ) &
              * ( (Gamma - 1.0_DP) / (Gamma + 1.0_DP) )**Gamma

    ! Loop over elements
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          ! Loop over all nodes in an element
          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) ! Physical coordinate

            IF( X1 <= rShock )THEN

              CALL ComputeSettlingSpeed_Bisection &
                     ( X1, Alpha, Gamma, Mass, Speed )

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = (mDot/FourPi) * Speed**(-1.0_DP) * X1**(-2.0_DP)

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = - Speed

              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP

              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP

              ! --- Post Shock Values from RH Conditions ---

              V1_prime &
                = (Gamma - 1.0_DP)/(Gamma + 1.0_DP) &
                   * SQRT( 2.0_DP * Mass / rShock )

              D_prime  &
                = (mDot/FourPi) * (1.0_DP/V1_prime) * rShock**(-2.0_DP)

              P_prime  &
                = 2.0_DP /(Gamma + 1.0_DP) * (mDot/FourPi) &
                    * SQRT(2.0_DP * Mass) * rShock**(-2.5_DP)

              uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                = P_prime &
                    * ( uPF(iNodeX,iX1,iX2,iX3,iPF_D) / D_prime )**Gamma &
                      / ( Gamma-1.0_DP )

            ELSE

              Speed &
                = SQRT &
                    (1.0_DP/(1.0_DP + 2.0_DP/((Gamma-1.0_DP)*Mach**2))/X1)

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = mDot / (FourPi * X1**2 * Speed )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = -Speed
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uPF(iNodeX, iX1, iX2, iX3, iPF_D) / Gamma &
                    * (Speed / Mach )**2.0_DP /(Gamma-1.0_DP)

            END IF

          END DO

          CALL ComputeConserved_Euler_NonRelativistic &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields


  SUBROUTINE ApplyPerturbations(ShellIn, ShellOut, Order, PerturbParam)

    INTEGER,  INTENT(in) :: Order
    REAL(DP), INTENT(in) :: ShellIn, ShellOut, PerturbParam

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF ( X1 >= ShellIn .AND. X1 <= ShellOut ) THEN

              SELECT CASE( Order )

              CASE ( 0 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = (1.0_DP + PerturbParam) &
                      * uPF(iNodeX, iX1, iX2, iX3, iPF_D)

              CASE ( 1 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = (1.0_DP + PerturbParam * COS(X2)) &
                      * uPF(iNodeX, iX1, iX2, iX3, iPF_D)

              END SELECT

            END IF

          END DO

          CALL ComputeConserved_Euler_NonRelativistic &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO
      END DO
    END DO

  END SUBROUTINE


  SUBROUTINE ComputeSettlingSpeed_Bisection( r, alpha, gamma, mass, V1 )

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

      IF (ab < Tol_ab .AND. ABS( F_a ) / F_0 < Tol_F) Converged = .TRUE.

    END DO

    V1 = a

  END SUBROUTINE ComputeSettlingSpeed_Bisection


  REAL(DP) FUNCTION SettlingSpeedFun( u, r, alpha, gamma, mass )

    REAL(DP), INTENT(in) :: u, r, alpha, gamma, mass

    SettlingSpeedFun &
      = r * u**2 &
        + alpha * r**(3.0_DP-2.0_DP*gamma) * u**(1.0_DP-gamma) &
        - 2.0_DP * mass
    RETURN
  END FUNCTION SettlingSpeedFun


END MODULE InitializationModule
