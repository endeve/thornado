MODULE GravitySolutionModule_Newtonian_PointMass

  USE KindModule, ONLY: &
    DP, Zero, One
  USE UnitsModule, ONLY: &
    GravitationalConstant
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesLX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_L2G
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    iGF_Phi_N, &
    nGF

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
    
    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: XC, dX, xL_q, xG_q
    REAL(DP) :: G_L(nDOFX,nGF)

    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          XC = MeshX(1) % Center(iX1)
          dX = MeshX(1) % Width (iX1)

          DO iNodeX = 1, nDOFX

            xL_q = NodesLX_q(1,iNodeX)
         
            xG_q = XC + dX * xL_q

            ! --- Compute Newtonian Potential 
               
            G_L(iNodeX,iGF_Phi_N) = - GravitationalConstant * Mass / xG_q

          END DO
          
          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_Phi_N), 1, Zero, G(:,iX1,iX2,iX3,iGF_Phi_N), 1 )
          
          END DO
        END DO
      END DO

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


END MODULE GravitySolutionModule_Newtonian_PointMass
