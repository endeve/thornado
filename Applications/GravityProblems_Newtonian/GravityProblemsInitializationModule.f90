MODULE GravityProblemsInitializationModule

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX, &
    Interpolate1D_Linear
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Gm, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeHomogeneousSphere
  PUBLIC :: InitializeHydrostaticPolytrope
  PUBLIC :: InitializeHomologousCollapse
  PUBLIC :: InitializeEvrardsCollapse

CONTAINS


  SUBROUTINE InitializeHomogeneousSphere( SphereRadius_Option )

    REAL(DP), INTENT(in), OPTIONAL :: SphereRadius_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: SphereRadius, X1

    SphereRadius = 1.0_DP
    IF( PRESENT( SphereRadius_Option ) ) &
      SphereRadius = SphereRadius_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A16,ES10.3E2)') &
      '', 'Sphere Radius = ', SphereRadius
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= SphereRadius )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 0.0_DP

                END IF

              END DO
            END DO
          END DO

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeHomogeneousSphere


  SUBROUTINE InitializeHydrostaticPolytrope

    INTEGER             :: iX1, iX2, iX3
    INTEGER             :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP)            :: X1
    REAL(DP), PARAMETER :: Gamma = 2.0_DP
    REAL(DP), PARAMETER :: Rho_M = 1.0d-6
    REAL(DP), PARAMETER :: Rho_C = 1.0_DP
    REAL(DP), PARAMETER :: Kappa = 2.0_DP / Pi
    REAL(DP), PARAMETER :: Alpha = Pi

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A8,ES10.3E2)') '', 'Rho_C = ', Rho_C
    WRITE(*,'(A4,A8,ES10.3E2)') '', 'Kappa = ', Kappa
    WRITE(*,'(A4,A8,ES10.3E2)') '', 'Alpha = ', Alpha
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
                  = MAX( Rho_C * SIN( Alpha * X1 ) / ( Alpha * X1 ), Rho_M )
                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                  = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                  = Kappa * uPF(iNodeX,iX1,iX2,iX3,iPF_D)**Gamma

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeHydrostaticPolytrope


  SUBROUTINE InitializeHomologousCollapse

    INTEGER             :: iX1, iX2, iX3, nR
    INTEGER             :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP)            :: X1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Radius, Density, Mass, Pressure

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    nR = 3331
    ALLOCATE( Radius(nR), Density(nR), Mass(nR), Pressure(nR) )

    CALL ReadLaneEmdenProfile &
           ( '../LaneEmden_n_3.dat', nR, Radius, Density, Mass, Pressure )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
                  = Interpolate1D( Radius, Density, nR, X1 )
                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                  = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                  = Interpolate1D( Radius, Pressure, nR, X1 )

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeHomologousCollapse


  SUBROUTINE InitializeEvrardsCollapse

    INTEGER             :: iX1, iX2, iX3
    INTEGER             :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP)            :: X1
    REAL(DP), PARAMETER :: Radius = 1.0_DP
    REAL(DP), PARAMETER :: Mass   = 1.0_DP

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A9,ES10.3E2)') '', 'Radius = ', Radius
    WRITE(*,'(A4,A9,ES10.3E2)') '', 'Mass   = ', Mass
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

!!$                IF( X1 <= Radius )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = Mass / ( 2.0_DP * Pi * Radius**2 * X1 )

!!$                ELSE
!!$
!!$                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
!!$                    = 1.0d-6
!!$
!!$                END IF

                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                  = 0.0_DP
                uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                  = 0.1_DP * uPF(iNodeX,iX1,iX2,iX3,iPF_D) / 3.0_DP

              END DO
            END DO
          END DO

          CALL ComputeInternalEnergyDensityFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeEvrardsCollapse


  SUBROUTINE ReadLaneEmdenProfile( FileName, nR, R, Den, Mass, Prs )

    CHARACTER(*)            :: FileName
    INTEGER                 :: nR
    REAL(DP), DIMENSION(nR) :: R, Den, Mass, Prs

    INTEGER :: i

    open(unit=10, file=Trim(Filename))

    do i = 1, nR
      read(10,*) R(i), Den(i), Mass(i), Prs(i)
    enddo

    close(10)

  END SUBROUTINE ReadLaneEmdenProfile


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


END MODULE GravityProblemsInitializationModule
