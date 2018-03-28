MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, &
    One, Three, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Kelvin, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX, nDOFE, &
    nZ, nNodes
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersX, &
    NodeNumberTable, &
    OuterProduct1D3D
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive
  USE MomentEquationsUtilitiesModule_Beta, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeFields
  PUBLIC :: ComputeError

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

      CASE ( 'PackedBeam' )

        CALL InitializeFields_PackedBeam

      CASE ( 'LineSource' )

        CALL InitializeFields_LineSource

      CASE ( 'HomogeneousSphere' )

        CALL InitializeFields_HomogeneousSphere

      CASE ( 'DeleptonizationWave' )

        CALL InitializeFields_DeleptonizationWave

    END SELECT

    wTime = MPI_WTIME( ) - wTime

!!$    WRITE(*,*)
!!$    WRITE(*,'(A4,A,ES10.4E2)') &
!!$      '', 'InitializeFields: ', wTime
!!$    WRITE(*,*)

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( Direction_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option

    CHARACTER(8) :: Direction
    INTEGER      :: iE, iX1, iX2, iX3, iS
    INTEGER      :: iNodeX1
    INTEGER      :: iNodeX2
    INTEGER      :: iNodeX3
    INTEGER      :: iNode
    REAL(DP)     :: X1, X2, X3, X, L
    REAL(DP)     :: Gm_dd_11(nDOF)
    REAL(DP)     :: Gm_dd_22(nDOF)
    REAL(DP)     :: Gm_dd_33(nDOF)
    REAL(DP)     :: Ones(nDOFE)

    Direction = 'X'
    IF( PRESENT( Direction_Option ) ) &
      Direction = TRIM( Direction_Option )
    
    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'Direction = ', TRIM( Direction )

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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

            DO iE = iE_B0, iE_E0

              SELECT CASE ( TRIM( Direction ) )
              CASE ( 'X' )

                DO iNode = 1, nDOF

                  iNodeX1 = NodeNumberTable(2,iNode)

                  X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )

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
                    = 0.50_DP + 0.49_DP * SIN( TwoPi * X2 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = 0.0_DP

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                    = uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                    = 0.0_DP

                END DO

              CASE ( 'Z' )

                DO iNode = 1, nDOF

                  iNodeX3 = NodeNumberTable(4,iNode)

                  X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP + 0.49_DP * SIN( TwoPi * X3 )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = 0.0_DP

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                    = 0.0_DP

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                    = uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

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
                    = 0.50_DP + 0.49_DP * SIN( TwoPi * X / L )

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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

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


  SUBROUTINE InitializeFields_PackedBeam

    INTEGER   :: iE, iX1, iX2, iX3, iS
    INTEGER   :: iNodeX1
    INTEGER   :: iNode
    REAL(DP)  :: X1, Mu_MP
    REAL(DP)  :: Delta
    REAL(DP)  :: Gm_dd_11(nDOF)
    REAL(DP)  :: Gm_dd_22(nDOF)
    REAL(DP)  :: Gm_dd_33(nDOF)
    REAL(DP)  :: Ones(nDOFE)

    Ones = 1.0_DP

    Mu_MP = 0.0_DP
    Delta = 1.0d-8

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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

            DO iE = iE_B0, iE_E0

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                IF( X1 .LE. ZERO )THEN

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = 0.50_DP * ( 1.0_DP + Delta )

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = 0.25_DP * ( 1.0_DP - Delta )

                ELSE

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                    = Delta

                  uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                    = 0.0d-0

                END IF

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

  END SUBROUTINE InitializeFields_PackedBeam


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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

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


  SUBROUTINE InitializeFields_HomogeneousSphere

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
                  ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

            DO iE = iE_B0, iE_E0

              DO iNode = 1, nDOF

                uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                  = 1.0d-6

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

  END SUBROUTINE InitializeFields_HomogeneousSphere


  SUBROUTINE ComputeError( Time, SigmaA, SigmaS )

    REAL(DP), INTENT(in) :: Time
    REAL(DP), INTENT(in) :: SigmaA
    REAL(DP), INTENT(in) :: SigmaS

    LOGICAL  :: ReportError
    REAL(DP) :: Error_One(nCR)
    REAL(DP) :: Error_Inf(nCR)

    SELECT CASE( TRIM( ProgramName ) )

      CASE ( 'SineWaveStreaming' )

        ReportError = .TRUE.

        CALL ComputeError_SineWaveStreaming &
               ( Time, Error_One, Error_Inf )

      CASE ( 'SineWaveDamping' )

        ReportError = .TRUE.

        CALL ComputeError_SineWaveDamping &
               ( Time, SigmaA, Error_One, Error_Inf )

      CASE ( 'SineWaveDiffusion' )

        ReportError = .TRUE.

        CALL ComputeError_SineWaveDiffusion &
               ( Time, Error_One, Error_Inf )

      CASE DEFAULT

        ReportError = .FALSE.

    END SELECT

    IF( ReportError )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') '', '--- Error Analysis ---'
      WRITE(*,*)
      WRITE(*,'(A4,A10,3I5.4)') '', 'nX = ', nZ(2:4)
      WRITE(*,'(A4,A10,1I5.4)') '', 'nE = ', nZ(1)
      WRITE(*,'(A4,A10,1I5.4)') '', 'nNodes = ', nNodes
      WRITE(*,*)
      WRITE(*,'(A4,A12,A4,A12)') '', 'L_1', '', 'L_Inf'
      WRITE(*,*)
      WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
        '', 'N: ',  Error_One(iCR_N ), '', Error_Inf(iCR_N )
      WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
        '', 'G1: ', Error_One(iCR_G1), '', Error_Inf(iCR_G1)
      WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
        '', 'G2: ', Error_One(iCR_G2), '', Error_Inf(iCR_G2)
      WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
        '', 'G3: ', Error_One(iCR_G3), '', Error_Inf(iCR_G3)
      WRITE(*,*)

    END IF

  END SUBROUTINE ComputeError


  SUBROUTINE ComputeError_SineWaveStreaming( Time, Error_One, Error_Inf )

    REAL(DP), INTENT(in)  :: Time
    REAL(DP), INTENT(out) :: Error_One(nCR)
    REAL(DP), INTENT(out) :: Error_Inf(nCR)

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: X1, N_A, G1_A, G2_A, G3_A

    Error_One = Zero
    Error_Inf = Zero

    DO iS = 1, nSpecies
      DO iX3 = 1, nZ(4)
        DO iX2 = 1, nZ(3)
          DO iX1 = 1, nZ(2)
            DO iE  = 1, nZ(1)

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                N_A  = 0.50_DP + 0.49_DP * SIN( TwoPi * ( X1 - Time ) )
                G1_A = N_A
                G2_A = Zero
                G3_A = Zero

                ! --- L1 Error ---

                Error_One(iCR_N) &
                  = Error_One(iCR_N) &
                    + ABS( N_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) )

                Error_One(iCR_G1) &
                  = Error_One(iCR_G1) &
                    + ABS( G1_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) )

                Error_One(iCR_G2) &
                  = Error_One(iCR_G2) &
                    + ABS( G2_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) )

                Error_One(iCR_G3) &
                  = Error_One(iCR_G3) &
                    + ABS( G3_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) )

                ! --- Infinity Error ---

                Error_Inf(iCR_N) &
                  = MAX( Error_Inf(iCR_N), &
                         ABS( N_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) ) )

                Error_Inf(iCR_G1) &
                  = MAX( Error_Inf(iCR_G1), &
                         ABS( G1_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) ) )

                Error_Inf(iCR_G2) &
                  = MAX( Error_Inf(iCR_G2), &
                         ABS( G2_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) ) )

                Error_Inf(iCR_G3) &
                  = MAX( Error_Inf(iCR_G3), &
                         ABS( G3_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) ) )

              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

    Error_One = Error_One / DBLE( nDOF * PRODUCT( nZ ) * nSpecies )

  END SUBROUTINE ComputeError_SineWaveStreaming


  SUBROUTINE ComputeError_SineWaveDamping( Time, SigmaA, Error_One, Error_Inf )

    REAL(DP), INTENT(in)  :: Time
    REAL(DP), INTENT(in)  :: SigmaA
    REAL(DP), INTENT(out) :: Error_One(nCR)
    REAL(DP), INTENT(out) :: Error_Inf(nCR)

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: X1, N_A, G1_A, G2_A, G3_A

    Error_One = Zero
    Error_Inf = Zero

    DO iS = 1, nSpecies
      DO iX3 = 1, nZ(4)
        DO iX2 = 1, nZ(3)
          DO iX1 = 1, nZ(2)
            DO iE  = 1, nZ(1)

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                N_A  = ( 0.50_DP + 0.49_DP * SIN( TwoPi * ( X1 - Time ) ) ) &
                         * EXP( - SigmaA * Time )
                G1_A = N_A
                G2_A = Zero
                G3_A = Zero

                ! --- L1 Error (Relative) ---

                Error_One(iCR_N) &
                  = Error_One(iCR_N) &
                    + ABS( N_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) ) / N_A

                Error_One(iCR_G1) &
                  = Error_One(iCR_G1) &
                    + ABS( G1_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) ) / N_A

                Error_One(iCR_G2) &
                  = Error_One(iCR_G2) &
                    + ABS( G2_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) ) / N_A

                Error_One(iCR_G3) &
                  = Error_One(iCR_G3) &
                    + ABS( G3_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) ) / N_A

                ! --- Infinity Error (Relative) ---

                Error_Inf(iCR_N) &
                  = MAX( Error_Inf(iCR_N), &
                         ABS(N_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS))/N_A )

                Error_Inf(iCR_G1) &
                  = MAX( Error_Inf(iCR_G1), &
                         ABS(G1_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS))/N_A )

                Error_Inf(iCR_G2) &
                  = MAX( Error_Inf(iCR_G2), &
                         ABS(G2_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS))/N_A )

                Error_Inf(iCR_G3) &
                  = MAX( Error_Inf(iCR_G3), &
                         ABS(G3_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS))/N_A )

              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

    Error_One = Error_One / DBLE( nDOF * PRODUCT( nZ ) * nSpecies )

  END SUBROUTINE ComputeError_SineWaveDamping


  SUBROUTINE ComputeError_SineWaveDiffusion( Time, Error_One, Error_Inf )

    REAL(DP), INTENT(in)  :: Time
    REAL(DP), INTENT(out) :: Error_One(nCR)
    REAL(DP), INTENT(out) :: Error_Inf(nCR)

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: X1, N_A, G1_A, G2_A, G3_A

    Error_One = Zero
    Error_Inf = Zero

    DO iS = 1, nSpecies
      DO iX3 = 1, nZ(4)
        DO iX2 = 1, nZ(3)
          DO iX1 = 1, nZ(2)
            DO iE  = 1, nZ(1)

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable(2,iNode)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                N_A  = ( 0.50_DP + 0.49_DP * SIN( Third * Pi * X1 ) &
                                     * EXP( - Pi**2 * Time / 2.7d3 ) )
                G1_A = - ( 0.49_DP * Pi / 9.0d2 ) * COS( Third * Pi * X1 ) &
                           * EXP( - Pi**2 * Time / 2.7d3 )
                G2_A = Zero
                G3_A = Zero

                ! --- L1 Error ---

                Error_One(iCR_N) &
                  = Error_One(iCR_N) &
                    + ABS( N_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) )

                Error_One(iCR_G1) &
                  = Error_One(iCR_G1) &
                    + ABS( G1_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) )

                Error_One(iCR_G2) &
                  = Error_One(iCR_G2) &
                    + ABS( G2_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) )

                Error_One(iCR_G3) &
                  = Error_One(iCR_G3) &
                    + ABS( G3_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) )

                ! --- Infinity Error ---

                Error_Inf(iCR_N) &
                  = MAX( Error_Inf(iCR_N), &
                         ABS(N_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS)) )

                Error_Inf(iCR_G1) &
                  = MAX( Error_Inf(iCR_G1), &
                         ABS(G1_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS)) )

                Error_Inf(iCR_G2) &
                  = MAX( Error_Inf(iCR_G2), &
                         ABS(G2_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS)) )

                Error_Inf(iCR_G3) &
                  = MAX( Error_Inf(iCR_G3), &
                         ABS(G3_A-uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS)) )

              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

    Error_One = Error_One / DBLE( nDOF * PRODUCT( nZ ) * nSpecies )

  END SUBROUTINE ComputeError_SineWaveDiffusion


  SUBROUTINE InitializeFields_DeleptonizationWave

    CALL InitializeFluidFields_DeleptonizationWave

    CALL InitializeRadiationFields_DeleptonizationWave

  END SUBROUTINE InitializeFields_DeleptonizationWave


  SUBROUTINE InitializeFluidFields_DeleptonizationWave

    ! --- Density Profile ---
    REAL(DP), PARAMETER :: MinD = 1.0d08 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: MaxD = 4.0d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: R_D  = 2.0d01 * Kilometer
    REAL(DP), PARAMETER :: H_D  = 1.0d01 * Kilometer
    ! --- Temperature Profile ---
    REAL(DP), PARAMETER :: MinT = 5.0d09 * Kelvin
    REAL(DP), PARAMETER :: MaxT = 2.6d11 * Kelvin
    REAL(DP), PARAMETER :: R_T  = 2.5d01 * Kilometer
    REAL(DP), PARAMETER :: H_T  = 2.0d01 * Kilometer
    ! --- Electron Fraction Profile ---
    REAL(DP), PARAMETER :: MinY = 3.0d-1
    REAL(DP), PARAMETER :: MaxY = 4.6d-1
    REAL(DP), PARAMETER :: R_Y  = 4.5d01 * Kilometer
    REAL(DP), PARAMETER :: H_Y  = 1.0d01 * Kilometer

    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, R

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)
            iNodeX3 = NodeNumberTableX(3,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            R = SQRT( X1**2 + X2**2 + X3**2 )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = Half * ( MaxD * ( One - TANH( (R-R_D)/H_D ) ) &
                         + MinD * ( One - TANH( (R_D-R)/H_D ) ) )

            uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
              = Half * ( MaxT * ( One - TANH( (R-R_T)/H_T ) ) &
                         + MinT * ( One - TANH( (R_T-R)/H_T ) ) )

            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = Half * ( MinY * ( One - TANH( (R-R_Y)/H_Y ) ) &
                         + MaxY * ( One - TANH( (R_Y-R)/H_Y ) ) )

          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
                   uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
                   uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFluidFields_DeleptonizationWave


  SUBROUTINE InitializeRadiationFields_DeleptonizationWave

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    REAL(DP) :: kT(nDOF)
    REAL(DP) :: Mnu(nDOF), E
    REAL(DP) :: Gm_dd_11(nDOF)
    REAL(DP) :: Gm_dd_22(nDOF)
    REAL(DP) :: Gm_dd_33(nDOF)

    DO iS = 1, nSpecies

      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)

            Gm_dd_11 &
              = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)

            Gm_dd_22 &
              = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)

            Gm_dd_33 &
              = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)

            kT = BoltzmannConstant &
                 * uAF(NodeNumbersX,iX1,iX2,iX3,iAF_T)

            Mnu = uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn)

            DO iE = iE_B0, iE_E0

              DO iNode = 1, nDOF

                iNodeE = NodeNumberTable(1,iNode)

                E = NodeCoordinate( MeshE, iE, iNodeE )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
                  = MAX( One / ( EXP( (E-Mnu(iNode))/kT(iNode) ) + One ), &
                         1.0d-32 )

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
                  = Zero

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
                  = Zero

                uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
                  = Zero

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
                       Gm_dd_11, Gm_dd_22, Gm_dd_33 )

            END DO

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeRadiationFields_DeleptonizationWave


END MODULE InitializationModule
