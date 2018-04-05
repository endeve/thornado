MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Third, One, Three, &
    Pi, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX, nDOFE
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE ReferenceElementModule, ONLY: &
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

    REAL(DP) :: wTime

    wTime = MPI_WTIME( )

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE ( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming &
               ( Direction_Option = Direction_Option )

      CASE ( 'SineWaveDamping' )

        CALL InitializeFields_SineWaveDamping

      CASE ( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion

      CASE ( 'LineSource' )

        CALL InitializeFields_LineSource

    END SELECT

    wTime = MPI_WTIME( ) - wTime

    WRITE(*,*)
    WRITE(*,'(A4,A,ES10.4E2)') &
      '', 'InitializeFields: ', wTime
    WRITE(*,*)

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( Direction_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option

    CHARACTER(8) :: Direction
    INTEGER      :: iE, iX1, iX2, iX3, iS, iK
    INTEGER      :: iNodeE
    INTEGER      :: iNodeX1
    INTEGER      :: iNodeX2
    INTEGER      :: iNode
    REAL(DP)     :: X1, X2, X, L
    REAL(DP)     :: Gm_dd_11(nDOF)
    REAL(DP)     :: Gm_dd_22(nDOF)
    REAL(DP)     :: Gm_dd_33(nDOF)
    REAL(DP)     :: Ones(nDOFE)

    Direction = 'X'
    IF( PRESENT( Direction_Option ) ) &
      Direction = TRIM( Direction_Option )
    
    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'Direction = ', TRIM( Direction )
    WRITE(*,*)

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

              CASE ( 'XY' )

                DO iNode = 1, nDOF

                  iNodeX1 = NodeNumberTable(2,iNode)
                  iNodeX2 = NodeNumberTable(3,iNode)

                  X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                  X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                  X  = SQRT( 2.0_DP ) * X1 &
                         + SQRT( 2.0_DP ) * X2
                  L  = SQRT( 2.0_DP )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP + 0.45_DP * SIN( TwoPi * X / L )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = SQRT( 2.0_DP ) / 2.0_DP &
                        * uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                    = SQRT( 2.0_DP ) / 2.0_DP &
                        * uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

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

  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_SineWaveDamping

    INTEGER      :: iE, iX1, iX2, iX3, iS, iK
    INTEGER      :: iNodeX1
    INTEGER      :: iNode
    REAL(DP)     :: X1
    REAL(DP)     :: Gm_dd_11(nDOF)
    REAL(DP)     :: Gm_dd_22(nDOF)
    REAL(DP)     :: Gm_dd_33(nDOF)
    REAL(DP)     :: Ones(nDOFE)

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

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                  = 0.5_DP + 0.49_DP * SIN( TwoPi * X1 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                  = uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                  = 0.0_DP

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                  = 0.0_DP

              END DO

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

  END SUBROUTINE InitializeFields_SineWaveDamping


  SUBROUTINE InitializeFields_SineWaveDiffusion

    INTEGER      :: iE, iX1, iX2, iX3, iS, iK
    INTEGER      :: iNodeX1
    INTEGER      :: iNode
    REAL(DP)     :: X1
    REAL(DP)     :: Gm_dd_11(nDOF)
    REAL(DP)     :: Gm_dd_22(nDOF)
    REAL(DP)     :: Gm_dd_33(nDOF)
    REAL(DP)     :: Ones(nDOFE)

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

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                  = 0.49_DP * SIN( Third * Pi * X1 ) + 0.5_DP

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                  = - ( 0.49_DP * Pi / 9.0d2 ) * COS( Third * Pi * X1 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                  = 0.0_DP

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                  = 0.0_DP

              END DO

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

  END SUBROUTINE InitializeFields_SineWaveDiffusion


  SUBROUTINE InitializeFields_LineSource

    INTEGER   :: iE, iX1, iX2, iX3, iS
    INTEGER   :: iNodeX1
    INTEGER   :: iNodeX2
    INTEGER   :: iNode
    REAL(DP)  :: X1, X2, R
    REAL(DP)  :: Gm_dd_11(nDOF)
    REAL(DP)  :: Gm_dd_22(nDOF)
    REAL(DP)  :: Gm_dd_33(nDOF)
    REAL(DP)  :: Ones(nDOFE)

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

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)
                iNodeX2 = NodeNumberTable(3,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                R  = SQRT( X1**2 + X2**2 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                  = MAX( EXP( - 0.5_DP * ( R / 0.03_DP )**2 ), 1.0d-4 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                  = 0.0_DP

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                  = 0.0_DP

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                  = 0.0_DP

              END DO

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

  END SUBROUTINE InitializeFields_LineSource


END MODULE InitializationModule
