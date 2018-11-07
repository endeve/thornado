MODULE TwoMoment_MeshRefinementModule

  USE KindModule, ONLY: &
    DP, Zero, Half
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, nDimsX, &
    nNodesE, nDOFE, nDOF
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    NodeNumberTableX3D, &
    WeightsX_q
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_TwoMoment
  PUBLIC :: FinalizeMeshRefinement_TwoMoment
  PUBLIC :: Refine_TwoMoment
  PUBLIC :: Coarsen_TwoMoment

  INTEGER :: nFine, nFineX(3)
  REAL(DP), ALLOCATABLE :: ProjectionMatrix(:,:,:)

CONTAINS


  SUBROUTINE InitializeMeshRefinement_TwoMoment

    INTEGER :: iFine, iFineX1, iFineX2, iFineX3
    INTEGER :: i, k, iN1, iN2, iN3
    REAL(DP), ALLOCATABLE :: xiX1(:), xiX2(:), xiX3(:)

    nFineX = 1
    DO i = 1, nDimsX
      nFineX(i) = 2 ! --- Refinement Factor of 2 Assumed
    END DO
    nFine = PRODUCT( nFineX )

    ALLOCATE( ProjectionMatrix(nDOFX,nDOFX,nFine) )

    ALLOCATE( xiX1(nNodesX(1)) )
    ALLOCATE( xiX2(nNodesX(2)) )
    ALLOCATE( xiX3(nNodesX(3)) )

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      IF( nFineX(1) > 1 )THEN
        xiX1 = Half * ( NodesX1 + (-1)**iFineX1 * Half )
      ELSE
        xiX1 = Zero
      END IF

      IF( nFineX(2) > 1 )THEN
        xiX2 = Half * ( NodesX2 + (-1)**iFineX2 * Half )
      ELSE
        xiX2 = Zero
      END IF

      IF( nFineX(3) > 1 )THEN
        xiX3 = Half * ( NodesX3 + (-1)**iFineX3 * Half )
      ELSE
        xiX3 = Zero
      END IF

      ProjectionMatrix(:,:,iFine) = 0.0_DP
      DO k = 1, nDOFX
      DO i = 1, nDOFX

        DO iN3 = 1, nNodesX(3)
        DO iN2 = 1, nNodesX(2)
        DO iN1 = 1, nNodesX(1)

          ProjectionMatrix(i,k,iFine) &
            = ProjectionMatrix(i,k,iFine) &
                + WeightsX1(iN1) * WeightsX2(iN2) * WeightsX3(iN3) &
                  * L_X1(IndLX_Q(1,i)) % P( NodesX1(iN1) ) &
                  * L_X2(IndLX_Q(2,i)) % P( NodesX2(iN2) ) &
                  * L_X3(IndLX_Q(3,i)) % P( NodesX3(iN3) ) &
                  * L_X1(IndLX_Q(1,k)) % P( xiX1   (iN1) ) &
                  * L_X2(IndLX_Q(2,k)) % P( xiX2   (iN2) ) &
                  * L_X3(IndLX_Q(3,k)) % P( xiX3   (iN3) )

        END DO
        END DO
        END DO

      END DO
      END DO

    END DO
    END DO
    END DO

    DEALLOCATE( xiX1 )
    DEALLOCATE( xiX2 )
    DEALLOCATE( xiX3 )

  END SUBROUTINE InitializeMeshRefinement_TwoMoment


  SUBROUTINE FinalizeMeshRefinement_TwoMoment

    DEALLOCATE( ProjectionMatrix )

  END SUBROUTINE FinalizeMeshRefinement_TwoMoment


  SUBROUTINE Refine_TwoMoment( nE, nX, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nE, nX(3)
    REAL(DP), INTENT(in)  :: U_Crs(1:nDOF,1:nE)
    REAL(DP), INTENT(out) :: U_Fin(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3))

    INTEGER :: iFine, iFineX1, iFineX2, iFineX3, iE
    REAL(DP) :: U_Crs_P(nDOFX,nDOFE*nE)
    REAL(DP) :: U_Fin_P(nDOFX,nDOFE*nE)

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      U_Crs_P = Pack_TwoMoment( nE, U_Crs )

      U_Fin_P(:,:) = MATMUL( ProjectionMatrix(:,:,iFine), U_Crs_P(:,:) )

      DO iE = 1, nDOFE * nE

        U_Fin_P(:,iE) = U_Fin_P(:,iE) / WeightsX_q(:)

      END DO

      U_Fin(:,:,iFineX1,iFineX2,iFineX3) = Unpack_TwoMoment( nE, U_Fin_P )

    END DO
    END DO
    END DO

  END SUBROUTINE Refine_TwoMoment


  SUBROUTINE Coarsen_TwoMoment

  END SUBROUTINE Coarsen_TwoMoment


  FUNCTION Pack_TwoMoment( nE, U ) RESULT( U_P )

    INTEGER,  INTENT(in) :: nE
    REAL(DP), INTENT(in) :: U(nDOF,nE)
    REAL(DP)             :: U_P(nDOFX,nDOFE*nE)

    INTEGER :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeE, iE, iNode

    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      iNodeX = NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3)

      DO iE     = 1, nE
      DO iNodeE = 1, nNodesE

        iNode = NodeNumberTable4D(iNodeE,iNodeX1,iNodeX2,iNodeX3)

        U_P(iNodeX,(iE-1)*nNodesE+iNodeE) = U(iNode,iE)

      END DO
      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION Pack_TwoMoment


  FUNCTION Unpack_TwoMoment( nE, U_P ) RESULT( U )

    INTEGER,  INTENT(in) :: nE
    REAL(DP), INTENT(in) :: U_P(nDOFX,nDOFE*nE)
    REAL(DP)             :: U(nDOF,nE)

    INTEGER :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeE, iE, iNode

    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      iNodeX = NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3)

      DO iE     = 1, nE
      DO iNodeE = 1, nNodesE

        iNode = NodeNumberTable4D(iNodeE,iNodeX1,iNodeX2,iNodeX3)

        U(iNode,iE) = U_P(iNodeX,(iE-1)*nNodesE+iNodeE)

      END DO
      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION Unpack_TwoMoment


END MODULE TwoMoment_MeshRefinementModule
