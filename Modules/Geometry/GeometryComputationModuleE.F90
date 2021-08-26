MODULE GeometryComputationModuleE

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFE
  USE ReferenceElementModuleE, ONLY: &
    NodesE
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
      G(1:nDOFE,iE_B1:iE_E1,1:nGE)

    INTEGER  :: iE, iNodeE
    REAL(DP) :: EC, dE, eG_q

    DO iE = iE_B1, iE_E1

      EC = MeshE % Center(iE)
      dE = MeshE % Width (iE)

      ! --- Compute Geometry Fields in Gaussian Points ---

      DO iNodeE = 1, nDOFE

        eG_q = EC + dE * NodesE(iNodeE)

        G(iNodeE,iE,iGE_Ep0) = eG_q**0
        G(iNodeE,iE,iGE_Ep1) = eG_q**1
        G(iNodeE,iE,iGE_Ep2) = eG_q**2
        G(iNodeE,iE,iGE_Ep3) = eG_q**3

      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( G )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( G )
#endif

  END SUBROUTINE ComputeGeometryE


END MODULE GeometryComputationModuleE
