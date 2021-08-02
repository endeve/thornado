MODULE ReferenceElementModuleE_Lagrange

  USE KindModule, ONLY: &
    DP, Half, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nNodesE, nDOFE
  USE ReferenceElementModuleE, ONLY: &
    NodesE, NodesE_L
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:)  , ALLOCATABLE, PUBLIC :: LE_Dn
  REAL(DP), DIMENSION(:)  , ALLOCATABLE, PUBLIC :: LE_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LE_L2G

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: InterpMat_E

  PUBLIC :: InitializeReferenceElementE_Lagrange
  PUBLIC :: FinalizeReferenceElementE_Lagrange

CONTAINS


  SUBROUTINE InitializeReferenceElementE_Lagrange

    INTEGER :: iNodeE, jNodeE

    ALLOCATE( LE_Dn(nDOFE), LE_Up(nDOFE) )

    DO iNodeE = 1, nDOFE

      LE_Dn(iNodeE) = LagrangeP( - Half, iNodeE, NodesE, nNodesE )

      LE_Up(iNodeE) = LagrangeP( + Half, iNodeE, NodesE, nNodesE )

    END DO

    ALLOCATE( LE_L2G(nDOFE,nDOFE) )

    DO jNodeE = 1, nDOFE
    DO iNodeE = 1, nDOFE

      LE_L2G(iNodeE,jNodeE) &
        = LagrangeP( NodesE(iNodeE), jNodeE, NodesE_L, nNodesE )

    END DO
    END DO

    ALLOCATE( InterpMat_E(nDOFE+2,nDOFE) )

    InterpMat_E = Zero
    DO iNodeE = 1, nDOFE

      InterpMat_E(iNodeE ,iNodeE) = One
      InterpMat_E(nDOFE+1,iNodeE) = LE_Dn(iNodeE)
      InterpMat_E(nDOFE+2,iNodeE) = LE_Up(iNodeE)

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: InterpMat_E )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( InterpMat_E )
#endif

  END SUBROUTINE InitializeReferenceElementE_Lagrange


  SUBROUTINE FinalizeReferenceElementE_Lagrange

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: InterpMat_E )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( InterpMat_E )
#endif

    DEALLOCATE( LE_Dn, LE_Up, LE_L2G )
    DEALLOCATE( InterpMat_E )

  END SUBROUTINE FinalizeReferenceElementE_Lagrange


END MODULE ReferenceElementModuleE_Lagrange
