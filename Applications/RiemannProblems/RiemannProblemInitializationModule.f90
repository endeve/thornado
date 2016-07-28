MODULE RiemannProblemInitializationModule

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, &
    uAF, iAF_P, iAF_Gm, iAF_Ye, iAF_T
  USE EquationOfStateModule, ONLY: &
    Auxiliary_Fluid, &
    InternalEnergy_Auxiliary, &
    ComputeTemperatureFromPressure
  USE EulerEquationsUtilitiesModule, ONLY: &
    Conserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRiemannProblem1D
  PUBLIC :: InitializeRiemannProblem1D_NuclearEOS
  PUBLIC :: InitializeInteractingBlastWaves1D  

CONTAINS


  SUBROUTINE InitializeRiemannProblem1D( D_L, V_L, P_L, D_R, V_R, P_R )

    REAL(DP),               INTENT(in) :: D_L, P_L, D_R, P_R
    REAL(DP), DIMENSION(3), INTENT(in) :: V_L, V_R

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'D_L = ', D_L, '', 'D_R = ', D_R
    WRITE(*,'(A7,A6,3ES10.3E2,A4,A6,3ES10.3E2)') &
      '', 'V_L = ', V_L, '', 'V_R = ', V_R
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'P_L = ', P_L, '', 'P_R = ', P_R

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= 0.5_DP )THEN

                  ! -- Left State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_L
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_L(1)
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_L(2)
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_L(3)
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_L

                ELSE

                  ! -- Right State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_R
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_R(1)
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_R(2)
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_R(3)
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_R

                END IF

                uPF(iNode,iX1,iX2,iX3,iPF_E) &
                  = InternalEnergy_Auxiliary &
                      ( uPF(iNode,iX1,iX2,iX3,:), uAF(iNode,iX1,iX2,iX3,:) )

                uAF(iNode,iX1,iX2,iX3,:) &
                  = Auxiliary_Fluid( uPF(iNode,iX1,iX2,iX3,:) )

                uCF(iNode,iX1,iX2,iX3,:) &
                  = Conserved( uPF(iNode,iX1,iX2,iX3,:) )

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeRiemannProblem1D

  SUBROUTINE InitializeRiemannProblem1D_NuclearEOS( D_L, V_L, P_L, Ye_L, D_R, V_R, P_R, Ye_R )

  REAL(DP),               INTENT(in) :: D_L, P_L, Ye_L, D_R, P_R, Ye_R
  REAL(DP), DIMENSION(3), INTENT(in) :: V_L, V_R

  INTEGER  :: iX1, iX2, iX3
  INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
  REAL(DP) :: X1

  WRITE(*,*)
  WRITE(*,'(A2,A6,A)') &
  '', 'INFO: ', TRIM( ProgramName )
  WRITE(*,*)
  WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
  '', 'D_L = ', D_L, '', 'D_R = ', D_R
  WRITE(*,'(A7,A6,3ES10.3E2,A4,A6,3ES10.3E2)') &
  '', 'V_L = ', V_L, '', 'V_R = ', V_R
  WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
  '', 'P_L = ', P_L, '', 'P_R = ', P_R
  WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
  '', 'Ye_L = ', Ye_L, '', 'Ye_R = ', Ye_R


  DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
      DO iX1 = 1, nX(1)

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

              iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              IF( X1 <= 0.5_DP )THEN

                ! -- Left State --

                uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_L
                uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_L(1)
                uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_L(2)
                uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_L(3)
                uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_L
                uAF(iNode,iX1,iX2,iX3,iAF_Ye) = Ye_L

              ELSE

                ! -- Right State --

                uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_R
                uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_R(1)
                uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_R(2)
                uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_R(3)
                uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_R
                uAF(iNode,iX1,iX2,iX3,iAF_Ye) = Ye_R

              END IF

              CALL ComputeTemperatureFromPressure &
                ( uPF(iNode:iNode,iX1,iX2,iX3,iPF_D), uAF(iNode:iNode,iX1,iX2,iX3,iAF_P), &
                  uAF(iNode:iNode,iX1,iX2,iX3,iAF_Ye), uAF(iNode:iNode,iX1,iX2,iX3,iAF_T) )

              uPF(iNode,iX1,iX2,iX3,iPF_E) &
                = InternalEnergy_Auxiliary &
                  ( uPF(iNode,iX1,iX2,iX3,:), uAF(iNode,iX1,iX2,iX3,:) )

              uAF(iNode,iX1,iX2,iX3,:) &
                = Auxiliary_Fluid( uPF(iNode,iX1,iX2,iX3,:) )

              uCF(iNode,iX1,iX2,iX3,:) &
                = Conserved( uPF(iNode,iX1,iX2,iX3,:) )

            END DO
          END DO
        END DO

      END DO
    END DO
  END DO

  END SUBROUTINE InitializeRiemannProblem1D_NuclearEOS

  SUBROUTINE InitializeInteractingBlastWaves1D

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

                uPF(iNode,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP

                IF( X1 < 0.1_DP )THEN

                  uAF(iNode,iX1,iX2,iX3,iAF_P) = 1.0d+3

                ELSE IF( X1 > 0.9_DP )THEN

                  uAF(iNode,iX1,iX2,iX3,iAF_P) = 1.0d+2

                ELSE

                  uAF(iNode,iX1,iX2,iX3,iAF_P) = 1.0d-1

                END IF

                uPF(iNode,iX1,iX2,iX3,iPF_E) &
                  = InternalEnergy_Auxiliary &
                      ( uPF(iNode,iX1,iX2,iX3,:), uAF(iNode,iX1,iX2,iX3,:) )

                uAF(iNode,iX1,iX2,iX3,:) &
                  = Auxiliary_Fluid( uPF(iNode,iX1,iX2,iX3,:) )

                uCF(iNode,iX1,iX2,iX3,:) &
                  = Conserved( uPF(iNode,iX1,iX2,iX3,:) )

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeInteractingBlastWaves1D




END MODULE RiemannProblemInitializationModule
