MODULE InitializationModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nNodesX, &
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
    uAF, iAF_P, iAF_Gm
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeConserved_GR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_GR
  PUBLIC :: InitializeFields_RiemannProblem

CONTAINS


  SUBROUTINE InitializeFields_GR

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! Loop over elements
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

           ! Loop over all nodes in an element
           DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) ! Physical coordinate

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = One + Half * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = One * 1.0d-1
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = Zero
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = Zero
            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = 0.01_DP * One
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = ( Gamma_IDEAL - 1.0_DP ) &
                  * uPF(iNodeX,iX1,iX2,iX3,iPF_E)            

           END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_GR


  SUBROUTINE InitializeFields_RiemannProblem &
    ( D_L, V_L, P_L, D_R, V_R, P_R, X_D_Option )

    REAL(DP), INTENT(in) :: D_L, P_L, D_R, P_R
    REAL(DP), INTENT(in) :: V_L(3), V_R(3)
    REAL(DP), INTENT(in), OPTIONAL :: X_D_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X_D, X1

    X_D = 0.5_DP
    IF( PRESENT( X_D_Option ) ) &
      X_D = X_D_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2)') &
      '', 'X_D = ', X_D
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'D_L = ', D_L, '', 'D_R = ', D_R
    WRITE(*,'(A7,A6,3ES10.3E2,A4,A6,3ES10.3E2)') &
      '', 'V_L = ', V_L, '', 'V_R = ', V_R
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'P_L = ', P_L, '', 'P_R = ', P_R

    ! Loop over elements
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          ! Loop over all nodes in an element
          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) ! Physical coordinate

            IF( X1 <= X_D )THEN

              ! -- Left State --

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_L
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_L(1)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_L(2)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_L(3)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_L / (Gamma_IDEAL-One)
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_L

            ELSE

              ! -- Right State --

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_R
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_R(1)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_R(2)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_R(3)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_R / (Gamma_IDEAL-One)
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_R

            END IF

!!$            ! --- Smooth solution ---
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
!!$              = 1.0_DP + 0.1_DP * SIN( TwoPi * X1 )
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6

!!$            ! --- Contact discontinuity solution (Top-hat) ---
!!$
!!$            IF( ( X1 <= 0.25d0 ) .OR. ( X1 >= 0.75d0 ) )THEN
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0d0
!!$            ELSE
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0d0
!!$            END IF
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6

          END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


END MODULE InitializationModule_GR
