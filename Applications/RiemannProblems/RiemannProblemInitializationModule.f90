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
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_E, iAF_Gm, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    Conserved, &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRiemannProblem1D
  PUBLIC :: InitializeSedov
  PUBLIC :: InitializeInteractingBlastWaves1D
  PUBLIC :: InitializeShockEntropyWaveInteraction1D

CONTAINS


  SUBROUTINE InitializeRiemannProblem1D &
               ( D_L, V_L, P_L, D_R, V_R, P_R, X_D_Option )

    REAL(DP),               INTENT(in) :: D_L, P_L, D_R, P_R
    REAL(DP), DIMENSION(3), INTENT(in) :: V_L, V_R
    REAL(DP),               INTENT(in), OPTIONAL :: X_D_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
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

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= X_D )THEN

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

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uPF(:,iX1,iX2,iX3,iPF_E),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P),  &
                   uAF(:,iX1,iX2,iX3,iAF_T),  uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uAF(:,iX1,iX2,iX3,iAF_Gm), &
                   uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeRiemannProblem1D


  SUBROUTINE InitializeSedov( E_0, nDetonationCells_Option )

    REAL(DP), INTENT(in) :: E_0
    INTEGER,  INTENT(in), OPTIONAL :: nDetonationCells_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    INTEGER  :: nDetonationCells
    REAL(DP) :: X1, R_0, Gm

    nDetonationCells = 1
    IF( PRESENT( nDetonationCells_Option ) ) &
      nDetonationCells = nDetonationCells_Option

    R_0 = REAL( nDetonationCells ) * MeshX(1) % Width(1)

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A,ES10.3E2,A2,I2.2,A2,ES10.3E2)') &
      '', 'E_0, nDetonationCells, R_0 = ', E_0, ', ', nDetonationCells, ', ', R_0

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

                Gm = 1.4_DP

                IF( X1 < R_0 )THEN

                  uAF(iNode, iX1, iX2, iX3, iAF_P) &
                    = ( Gm - 1.0_DP ) * 3.0_DP * E_0 &
                      / ( 4.0_DP * Pi * R_0**3 )

                ELSE

                  uAF(iNode, iX1, iX2, iX3, iAF_P) &
                    = 1.0d-5

                END IF

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uPF(:,iX1,iX2,iX3,iPF_E),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P),  &
                   uAF(:,iX1,iX2,iX3,iAF_T),  uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uAF(:,iX1,iX2,iX3,iAF_Gm), &
                   uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeSedov


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

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uPF(:,iX1,iX2,iX3,iPF_E),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P),  &
                   uAF(:,iX1,iX2,iX3,iAF_T),  uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uAF(:,iX1,iX2,iX3,iAF_Gm), &
                   uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeInteractingBlastWaves1D


  SUBROUTINE InitializeShockEntropyWaveInteraction1D

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X_D, Amplitude, Wavenumber, X1

    X_D        = - 4.0_DP
    Amplitude  = + 0.2_DP
    Wavenumber = + 5.0_DP

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2)') &
      '', 'X_D = ', X_D
    WRITE(*,*)
    WRITE(*,'(A7,A13,ES10.3E2)') '', 'Amplitude  = ', Amplitude
    WRITE(*,'(A7,A13,ES10.3E2)') '', 'Wavenumber = ', Wavenumber
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= X_D )THEN

                  ! -- Left State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = 3.857143_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = 2.629369_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = 10.333333_DP

                ELSE

                  ! -- Right State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D) &
                    = 1.0_DP + Amplitude * SIN( WaveNumber * X1 )
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = 1.0_DP

                END IF

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uPF(:,iX1,iX2,iX3,iPF_E),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P),  &
                   uAF(:,iX1,iX2,iX3,iAF_T),  uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uAF(:,iX1,iX2,iX3,iAF_Gm), &
                   uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeShockEntropyWaveInteraction1D


END MODULE RiemannProblemInitializationModule
