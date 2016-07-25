MODULE MomentEquationsLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDimsE, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    MinModB
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Radiation, &
    MapModalToNodal_Radiation
  USE MeshModule, ONLY: &
    MeshX
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: BetaTVD = 2.00_DP
  REAL(DP), PARAMETER :: BetaTVB = 50.0_DP
  REAL(DP), PARAMETER :: Tol_TVD = 1.0d-2

  PUBLIC :: InitializeLimiters_M1_DG
  PUBLIC :: ApplySlopeLimiter_M1_DG

CONTAINS


  SUBROUTINE InitializeLimiters_M1_DG

  END SUBROUTINE InitializeLimiters_M1_DG


  SUBROUTINE ApplySlopeLimiter_M1_DG

    LOGICAL :: LimitPolynomial
    INTEGER :: iS, iE, iX1, iX2, iX3, iCR, iOS
    REAL(DP), DIMENSION(nCR) :: uCR_A, uCR_A_P_X1, uCR_A_N_X1
    REAL(DP), DIMENSION(nDOF,nCR) :: uCR_M, uCR_M_P_X1, uCR_M_N_X1
    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: uCR_X1, uCR_X1_T

    IF( nDOF == 1 ) RETURN

    ALLOCATE &
      ( uCR_X1  (1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR), &
        uCR_X1_T(1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR) )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)) )

    iOS = 1 ! --- Offset to Access Position Space Slopes
    IF( nDimsE == 0 ) iOS = 0

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Limiting Using Modal Representation ---

              DO iCR = 1, nCR

                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1-1,iX2,iX3,iCR,iS), uCR_M_P_X1(:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1,  iX2,iX3,iCR,iS), uCR_M     (:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1+1,iX2,iX3,iCR,iS), uCR_M_N_X1(:,iCR) )

              END DO

              ! --- Cell-Averaged Moments ---

              uCR_A      = uCR_M     (1,1:nCR)
              uCR_A_P_X1 = uCR_M_P_X1(1,1:nCR)
              uCR_A_N_X1 = uCR_M_N_X1(1,1:nCR)

              ! --- Slope From Modal Representation ---

              uCR_X1(iE,iX1,iX2,iX3,1:nCR) = uCR_M(iOS+2,1:nCR) ! X1-Dimension

              uCR_X1_T(iE,iX1,iX2,iX3,1:nCR) &
                = MinModB &
                    ( uCR_X1(iE,iX1,iX2,iX3,1:nCR), &
                      BetaTVD * ( uCR_A - uCR_A_P_X1 ), &
                      BetaTVD * ( uCR_A_N_X1 - uCR_A ), &
                      dX1(iX1), BetaTVB )

            END DO
          END DO
        END DO
      END DO

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              LimitPolynomial = .FALSE.
              LimitPolynomial &
                = ANY( ABS( uCR_X1(iE,iX1,iX2,iX3,1:nCR) &
                              - uCR_X1_T(iE,iX1,iX2,iX3,1:nCR) ) > Tol_TVD )

              IF( LimitPolynomial )THEN

                DO iCR = 1, nCR

                  CALL MapNodalToModal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                END DO

                ! --- Cell-Averaged Moments ---

                uCR_A = uCR_M(1,1:nCR)

                uCR_M(:,1:nCR) = 0.0_DP
                uCR_M(1,1:nCR) &     ! -- Cell-Average
                  = uCR_A(1:nCR)
                uCR_M(iOS+2,1:nCR) & ! -- Slope X1-Direction
                  = uCR_X1_T(iE,iX1,iX2,iX3,1:nCR)

                ! --- Back to Nodal Representation ---

                DO iCR = 1, nCR

                  CALL MapModalToNodal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                END DO

                PRINT*, "iE, iX1, iX2, iX3 = ", iE, iX1, iX2, iX3
                PRINT*, "  uCR_X1   = ", uCR_X1  (iE,iX1,iX2,iX3,1:nCR)
                PRINT*, "  uCR_X1_T = ", uCR_X1_T(iE,iX1,iX2,iX3,1:nCR)
                PRINT*, "  ANY = ", &
                  ( ABS( uCR_X1(iE,iX1,iX2,iX3,1:nCR) &
                         - uCR_X1_T(iE,iX1,iX2,iX3,1:nCR) ) > Tol_TVD )
                PRINT*

              END IF

            END DO
          END DO
        END DO
      END DO

    END DO

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( uCR_X1, uCR_X1_T )

  END SUBROUTINE ApplySlopeLimiter_M1_DG


END MODULE MomentEquationsLimiterModule_DG
