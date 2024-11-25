MODULE InitializationModule
  
  USE KindModule, ONLY: &
    DP, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ScalarFieldsModule, ONLY: &
    iSF_u, iSF_v, &
    uSF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS

  SUBROUTINE InitializeFields

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1
    REAL(DP) :: x

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
      
      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        x = NodeCoordinate(MeshX(1),iX1,iNodeX1)

        uSF(iNodeX,iX1,iX2,iX3,iSF_u) = 0.0_dp
        uSF(iNodeX,iX1,iX2,iX3,iSF_v) = SIN(x)

!-- Analytic solution given by u(x,t) = 1/c*SIN(x)*SIN(c*t) & v(x,t)=SIN(x + ct)
        
      END DO
      
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields



END MODULE InitializationModule
