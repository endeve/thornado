MODULE XCFC_UtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    Two, &
    Pi
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Psi
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting, &
    ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6
  PUBLIC :: ApplyBoundaryConditions_Geometry_XCFC
  PUBLIC :: ComputeGravitationalMass

  ! --- GS: Gravity/Geometry Sources ---

  INTEGER, PARAMETER, PUBLIC :: iGS_E  = 1
  INTEGER, PARAMETER, PUBLIC :: iGS_S1 = 2
  INTEGER, PARAMETER, PUBLIC :: iGS_S2 = 3
  INTEGER, PARAMETER, PUBLIC :: iGS_S3 = 4
  INTEGER, PARAMETER, PUBLIC :: iGS_S  = 5
  INTEGER, PARAMETER, PUBLIC :: iGS_Mg = 6
  INTEGER, PARAMETER, PUBLIC :: nGS    = 6

CONTAINS


  SUBROUTINE MultiplyWithPsi6( iX_B1, iX_E1, G, U, Power )

    INTEGER,  INTENT(in)    :: iX_B1(3), iX_E1(3), Power
    REAL(DP), INTENT(in)    :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iFF, nFF

    nFF = SIZE( U, DIM = 5 )

    DO iFF = 1       , nFF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1       , nDOFX

      U(iNX,iX1,iX2,iX3,iFF) &
        = U(iNX,iX1,iX2,iX3,iFF) * G(iNX,iX1,iX2,iX3,iGF_Psi)**( 6 * Power )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MultiplyWithPsi6


  SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    CALL ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

  END SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC


  SUBROUTINE ComputeGravitationalMass &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, GS, GravitationalMass )

    INTEGER,  INTENT(in)   :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)   :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)   :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: &
      GravitationalMass

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: d3X

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    ! --- Mg has been pre-multiplied with \alpha * \psi^6

    ! --- Assuming 1D spherical symmetry ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      d3X = Two / Pi * dX1(iX1) * dX2(iX2) * dX3(iX3)

      GravitationalMass &
        = GravitationalMass + d3X                  &
            * SUM( WeightsX_q * GS(:,iX1,iX2,iX3,iGS_Mg) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGravitationalMass


END MODULE XCFC_UtilitiesModule
