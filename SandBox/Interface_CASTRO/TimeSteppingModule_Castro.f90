MODULE TimeSteppingModule_Castro

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, nDOF, &
    iX_B0, iX_B1, iX_E0, iX_E0, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE DiscretizationModule_Dummy, ONLY: &
    ComputeIncrement_Implicit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Update_IMEX_PC2

CONTAINS


  SUBROUTINE Update_IMEX_PC2( dt, U_F, U_R )

    use GeometryFieldsModuleE, only : uGE
    use GeometryFieldsModule, only : uGF

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(inout) :: &
      U_F(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U_R(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER  :: &
      iX_SW(3), iZ_SW(4)
    REAL(DP) :: &
      dU_F &
        (1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP) :: &
      U0_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies), &
      dU_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)

    ! --------------------------------------------------------------
    ! --- Precdictor-Corrector IMEX Scheme from McClarren et al. ---
    ! --- Journal of Computational Physics, 227, 7561-7586, 2008 ---
    ! --------------------------------------------------------------

    CALL AddFields_Radiation &
           ( iZ_B1, iZ_E1, One, Zero, U_R, U_R, U0_R )


    ! --- Predictor Step ---


    ! --- Include One Layer of Spatial Ghost Cells in Update

    iX_SW = [ 1, 1, 1 ]
    iZ_SW = [ 0, 1, 1, 1 ]

    ! --- Explicit Step (Radiation Only) ---

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
             uGE, uGF, &
             U_R, dU_R &
                    (1:nDOF, &
                     iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                     iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                     iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                     iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                     1:nCR, 1:nSpecies) )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U0_R, dU_R, U_R )

    ! --- Implicit Step ---

    CALL ComputeIncrement_Implicit &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
             Half * dt, &
             uGE, uGF, &
             U_F, dU_F &
                    (1:nDOFX, &
                     iX_B0(1)-iX_SW(1):iX_E0(1)+iX_SW(1), &
                     iX_B0(2)-iX_SW(2):iX_E0(2)+iX_SW(2), &
                     iX_B0(3)-iX_SW(3):iX_E0(3)+iX_SW(3), &
                     1:nCF), &
             U_R, dU_R &
                    (1:nDOF, &
                     iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                     iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                     iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                     iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                     1:nCR, 1:nSpecies) )

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, Half * dt, U_F, dU_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U_R, dU_R, U_R )


    ! --- Corrector Step ---


    iX_SW = [ 0, 0, 0 ]
    iZ_SW = [ 0, 0, 0, 0 ]

    ! --- Explicit Step (Radiation Only) ---

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
             uGE, uGF, &
             U_R, dU_R &
                    (1:nDOF, &
                     iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                     iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                     iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                     iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                     1:nCR, 1:nSpecies) )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, dt, U0_R, dU_R, U_R )

    ! --- Implicit Step ---

    CALL ComputeIncrement_Implicit &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
             dt, &
             uGE, uGF, &
             U_F, dU_F &
                    (1:nDOFX, &
                     iX_B0(1)-iX_SW(1):iX_E0(1)+iX_SW(1), &
                     iX_B0(2)-iX_SW(2):iX_E0(2)+iX_SW(2), &
                     iX_B0(3)-iX_SW(3):iX_E0(3)+iX_SW(3), &
                     1:nCF), &
             U_R, dU_R &
                    (1:nDOF, &
                     iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                     iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                     iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                     iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                     1:nCR, 1:nSpecies) )

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, dt, U_F, dU_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, dt, U_R, dU_R, U_R )

  END SUBROUTINE Update_IMEX_PC2


  SUBROUTINE AddFields_Fluid( iX_B, iX_E, alpha, beta, A, B )

    ! --- A := alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      B(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFF

    DO iFF = 1, nCR
      DO iX3 = iX_B(3), iX_E(3)
        DO iX2 = iX_B(2), iX_E(2)
          DO iX1 = iX_B(1), iX_E(1)

            A(:,iX1,iX2,iX3,iFF) &
              = alpha * A(:,iX1,iX2,iX3,iFF) &
                  + beta * B(:,iX1,iX2,iX3,iFF)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE AddFields_Fluid


  SUBROUTINE AddFields_Radiation( iZ_B, iZ_E, alpha, beta, A, B, C )

    INTEGER,  INTENT(in)    :: iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      B(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(out)   :: &
      C(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iRF, iS

    DO iS = 1, nSpecies
      DO iRF = 1, nCR
        DO iZ4 = iZ_B(4), iZ_E(4)
          DO iZ3 = iZ_B(3), iZ_E(3)
            DO iZ2 = iZ_B(2), iZ_E(2)
              DO iZ1 = iZ_B(1), iZ_E(1)

                C(:,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
                  = alpha * A(:,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
                      + beta * B(:,iZ1,iZ2,iZ3,iZ4,iRF,iS)

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE AddFields_Radiation


END MODULE TimeSteppingModule_Castro
