MODULE DiscretizationModule_Dummy

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Explicit
  PUBLIC :: ComputeIncrement_Implicit

CONTAINS


  SUBROUTINE ComputeIncrement_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX  (1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U_R (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)

                dU_R(:,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                  = 0.0_DP * U_R(:,iZ1,iZ2+1,iZ3,iZ4,iCR,iS) &
                      - 0.0_DP * U_R(:,iZ1,iZ2-1,iZ3,iZ4,iCR,iS)

              END DO
            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE ComputeIncrement_Explicit


  SUBROUTINE ComputeIncrement_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX  (1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U_F (1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U_R (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCF, iCR, iS

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)

                dU_R(:,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                  = 0.0_DP

              END DO
            END DO
          END DO
        END DO

      END DO
    END DO


    DO iCF = 1, nCF

      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)

            dU_F(:,iZ2,iZ3,iZ4,iCF) = 0.0_DP

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputeIncrement_Implicit


END MODULE DiscretizationModule_Dummy
