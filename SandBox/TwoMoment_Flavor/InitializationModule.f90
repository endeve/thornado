MODULE InitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Decoherence


CONTAINS


  SUBROUTINE InitializeFields_Decoherence

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1
    REAL(DP) :: R

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        ! Interpolate Chimera Data to R and set variables

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0d14 * Gram / Centimeter**3

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) = 1.0d10 * Kelvin

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 3.0d-1

      END DO

      CALL ApplyEquationOfState_TABLE &
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

  END SUBROUTINE InitializeFields_Decoherence


END MODULE InitializationModule
