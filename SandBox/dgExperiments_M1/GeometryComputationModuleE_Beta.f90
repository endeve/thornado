MODULE GeometryComputationModuleE_Beta

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFE
  USE ReferenceElementModuleE, ONLY: &
    NodesE_L
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    LE_L2G
  USE MeshModule, ONLY: &
    MeshE
  USE GeometryFieldsModuleE, ONLY: &
    iGE_Ep0, &
    iGE_Ep1, &
    iGE_Ep2, &
    iGE_Ep3, &
    nGE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryE

CONTAINS


  SUBROUTINE ComputeGeometryE( iE_B0, iE_E0, iE_B1, iE_E1, G )

    INTEGER, INTENT(in)     :: &
      iE_B0, iE_E0, iE_B1, iE_E1
    REAL(DP), INTENT(inout) :: &
      G(1:,iE_B1:,1:)

    INTEGER  :: iE, iNodeE
    REAL(DP) :: EC, dE, eG_q
    REAL(DP) :: G_L(nDOFE,nGE)

    DO iE = iE_B1, iE_E1

      EC = MeshE % Center(iE)
      dE = MeshE % Width (iE)

      ! --- Compute Geometry Fields in Lobatto Points ---

      DO iNodeE = 1, nDOFE

        ! --- Global Coordinates (Lobatto Points) ---

        eG_q = EC + dE * NodesE_L(iNodeE)

        ! --- Set Geometry in Lobatto Points ---

        G_L(iNodeE,iGE_Ep0) = eG_q**0
        G_L(iNodeE,iGE_Ep1) = eG_q**1
        G_L(iNodeE,iGE_Ep2) = eG_q**2
        G_L(iNodeE,iGE_Ep3) = eG_q**3

      END DO

      ! --- Interpolate from Lobatto to Gaussian Points ---

      CALL DGEMV &
             ( 'N', nDOFE, nDOFE, One, LE_L2G, nDOFE, &
               G_L(:,iGE_Ep0), 1, Zero, G(:,iE,iGE_Ep0), 1 )
      CALL DGEMV &
             ( 'N', nDOFE, nDOFE, One, LE_L2G, nDOFE, &
               G_L(:,iGE_Ep1), 1, Zero, G(:,iE,iGE_Ep1), 1 )
      CALL DGEMV &
             ( 'N', nDOFE, nDOFE, One, LE_L2G, nDOFE, &
               G_L(:,iGE_Ep2), 1, Zero, G(:,iE,iGE_Ep2), 1 )
      CALL DGEMV &
             ( 'N', nDOFE, nDOFE, One, LE_L2G, nDOFE, &
               G_L(:,iGE_Ep3), 1, Zero, G(:,iE,iGE_Ep3), 1 )

    END DO

  END SUBROUTINE ComputeGeometryE


END MODULE GeometryComputationModuleE_Beta
