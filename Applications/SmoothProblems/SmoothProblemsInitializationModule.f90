MODULE SmoothProblemsInitializationModule

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
    uAF, iAF_P, iAF_T, iAF_S, iAF_Ye, iAF_E, iAF_Gm, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeLinearWaves1D

CONTAINS


  SUBROUTINE InitializeLinearWaves1D( WaveType, Amplitude_Option )

    CHARACTER(LEN=*) :: WaveType
    REAL(DP), INTENT(in), OPTIONAL :: Amplitude_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1, Amplitude

    Amplitude = 1.0d-6
    IF( PRESENT( Amplitude_Option ) ) &
      Amplitude = Amplitude_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A)') &
      '', TRIM( WaveType )
    WRITE(*,'(A7,A12,ES10.3E2)') &
      '', 'Amplitude = ', Amplitude

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
                uAF(iNode,iX1,iX2,iX3,iAF_P)  = 3.0_DP / 5.0_DP

                SELECT CASE ( TRIM( WaveType ) )

                  CASE ( 'AcousticWaves' )

                    uPF(iNode,iX1,iX2,iX3,iPF_D) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_D) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )
                    uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )
                    uPF(iNode,iX1,iX2,iX3,iAF_P) &
                      = uPF(iNode,iX1,iX2,iX3,iAF_P) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )

                  CASE ( 'EntropyWaves_Velocity' )

                    uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                          + 0.5_DP
                    uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )
                    uPF(iNode,iX1,iX2,iX3,iPF_V3) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )

                  CASE ( 'EntropyWaves_Density' )

                    uPF(iNode,iX1,iX2,iX3,iPF_D) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_D) &
                          + 0.5_DP * SIN( 2.0_DP * Pi * X1 )
                    uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                          + 1.0_DP

                  CASE DEFAULT

                    WRITE(*,'(A7,A18,A)') &
                      '', 'Invalid Wave Type ', TRIM( WaveType )
                    STOP

                END SELECT

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

  END SUBROUTINE InitializeLinearWaves1D


END MODULE SmoothProblemsInitializationModule
