MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX, nDOFE
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE ReferenceElementModule_Beta, ONLY: &
    NodeNumberTable, &
    OuterProduct1D3D
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE MomentEquationsUtilitiesModule_Beta, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields( Direction_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option

    CHARACTER :: Direction
    INTEGER   :: iE, iX1, iX2, iX3, iS, iK
    INTEGER   :: iNodeE
    INTEGER   :: iNodeX1
    INTEGER   :: iNodeX2
    INTEGER   :: iNode
    REAL(DP)  :: X1, X2
    REAL(DP)  :: Gm_dd_11(nDOF)
    REAL(DP)  :: Gm_dd_22(nDOF)
    REAL(DP)  :: Gm_dd_33(nDOF)
    REAL(DP)  :: Ones(nDOFE)
    REAL(DP)  :: wTime

    Direction = 'X'
    IF( PRESENT( Direction_Option ) ) &
      Direction = TRIM( Direction_Option )

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'Direction = ', TRIM( Direction )
    WRITE(*,*)

    wTime = MPI_WTIME( )

    Ones = 1.0_DP

    DO iS = 1, nSpecies
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)

            Gm_dd_11 &
              = OuterProduct1D3D &
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), nDOFX )

            Gm_dd_22 &
              = OuterProduct1D3D &
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), nDOFX )

            Gm_dd_33 &
              = OuterProduct1D3D &
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), nDOFX )

            DO iE = iE_B0, iE_E0

              SELECT CASE ( TRIM( Direction ) )
              CASE ( 'X' )

                DO iNode = 1, nDOF

                  iNodeE  = NodeNumberTable(1,iNode)
                  iNodeX1 = NodeNumberTable(2,iNode)

                  iK = 1 !( iE - 1 ) * nDOFE + iNodeE

                  X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP + 0.45_DP * SIN( ( TwoPi * iK ) * X1 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                    = 0.0_DP

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                    = 0.0_DP

                END DO

              CASE ( 'Y' )

                DO iNode = 1, nDOF

                  iNodeX2 = NodeNumberTable(3,iNode)

                  X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP + 0.45_DP * SIN( TwoPi * X2 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = 0.0_DP

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                    = uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                    = 0.0_DP

                END DO

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A6,A,A)') &
                  '', 'Invalid Direction: ', TRIM( Direction )
                WRITE(*,*)
                STOP

              END SELECT

              CALL ComputeConserved &
                     ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                       uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                       uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                       uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                       uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                       uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                       uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                       uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                       Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:) )

            END DO
          END DO
        END DO
      END DO
    END DO

    wTime = MPI_WTIME( ) - wTime

    WRITE(*,*)
    WRITE(*,'(A4,A,ES10.4E2)') &
      '', 'InitializeFields: ', wTime
    WRITE(*,*)

  END SUBROUTINE InitializeFields


END MODULE InitializationModule
