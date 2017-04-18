MODULE InitializationModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, swX
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX, &
    Interpolate1D_Linear
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, iPF_D, iPF_V1, iPF_E, iPF_Ne, nPF, &
    uAF, iAF_T, iAF_Ye, iAF_E, iAF_P, iAF_S, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitationalCollapse

CONTAINS


  SUBROUTINE InitializeGravitationalCollapse( ProgenitorFile )

    CHARACTER(LEN=*), INTENT(in) :: ProgenitorFile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1
    TYPE(ProgenitorType1D) :: P1D

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    CALL ReadProgenitor1D( ProgenitorFile, P1D )

    ! --- Initialize Fluid Fields ---

    ASSOCIATE &
      ( R1D => P1D % Radius, &
        D1D => P1D % MassDensity, &
        V1D => P1D % RadialVelocity, &
        T1D => P1D % Temperature, &
        Y1D => P1D % ElectronFraction )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1) + swX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                  = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

                uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                  = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                  = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

    END ASSOCIATE ! R1D, etc

  END SUBROUTINE InitializeGravitationalCollapse


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

    END IF

    RETURN
  END FUNCTION Interpolate1D


END MODULE InitializationModule
