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
    uCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, &
    uAF, iAF_P, iAF_Gm
  USE EquationOfStateModule, ONLY: &
    Auxiliary_Fluid, &
    InternalEnergy_Auxiliary
  USE EulerEquationsUtilitiesModule, ONLY: &
    Conserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeLinearWaves1D

CONTAINS


  SUBROUTINE InitializeLinearWaves1D( WaveType )

    CHARACTER(LEN=*) :: WaveType

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1, Amplitude = 1.0d-6

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A)') &
      '', TRIM( WaveType )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNode,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                uAF(iNode,iX1,iX2,iX3,iAF_P)  = 3.0_DP / 5.0_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V1) = 0.5_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP

                SELECT CASE ( TRIM( WaveType ) )

                  CASE ( 'AcousticWaves' )

                  CASE ( 'EntropyWaves' )

                    uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )
                    uPF(iNode,iX1,iX2,iX3,iPF_V3) &
                      = uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                          + Amplitude * SIN( 2.0_DP * Pi * X1 )

                  CASE DEFAULT

                    WRITE(*,'(A7,A18,A)') &
                      '', 'Invalid Wave Type ', TRIM( WaveType )
                    STOP

                END SELECT

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

  END SUBROUTINE InitializeLinearWaves1D


END MODULE SmoothProblemsInitializationModule
