MODULE GravitySolutionModule_Newtonian_PointMass_Beta

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGravitationalPotential

CONTAINS


  SUBROUTINE ComputeGravitationalPotential &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ComputeGravitationalPotential_CARTESIAN &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE ( 'SPHERICAL' )

        CALL ComputeGravitationalPotential_SPHERICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE ( 'CYLINDRICAL' )

        CALL ComputeGravitationalPotential_CYLINDRICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

  END SUBROUTINE ComputeGravitationalPotential


  SUBROUTINE ComputeGravitationalPotential_CARTESIAN &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

  END SUBROUTINE ComputeGravitationalPotential_CARTESIAN


  SUBROUTINE ComputeGravitationalPotential_SPHERICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

  END SUBROUTINE ComputeGravitationalPotential_SPHERICAL


  SUBROUTINE ComputeGravitationalPotential_CYLINDRICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

  END SUBROUTINE ComputeGravitationalPotential_CYLINDRICAL


END MODULE GravitySolutionModule_Newtonian_PointMass_Beta
