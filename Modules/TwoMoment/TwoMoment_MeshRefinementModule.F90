MODULE TwoMoment_MeshRefinementModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_TwoMoment
  PUBLIC :: FinalizeMeshRefinement_TwoMoment
  PUBLIC :: RefineX_TwoMoment
  PUBLIC :: CoarsenX_TwoMoment

  INTEGER  :: nFine, nFineX(3)
  REAL(DP), PUBLIC :: VolumeRatio
  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix  (:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix_T(:,:) ! --- Transpose ---

CONTAINS


  SUBROUTINE InitializeMeshRefinement_TwoMoment

    CHARACTER(2)  :: MeshString
    CHARACTER(32) :: MatrixName
    INTEGER :: iFine, iFineX1, iFineX2, iFineX3
    INTEGER :: i, j, k, iN1, iN2, iN3
    REAL(DP), ALLOCATABLE :: xiX1(:), xiX2(:), xiX3(:)

    nFineX      = 1
    VolumeRatio = One
    DO i = 1, nDimsX
      ! --- Refinement Factor of 2 Assumed:
      nFineX(i)   = 2
      VolumeRatio = VolumeRatio / Two
    END DO
    nFine = PRODUCT( nFineX )

    ALLOCATE( ProjectionMatrix  (nDOFX*nFine,nDOFX      ) )
    ALLOCATE( ProjectionMatrix_T(nDOFX      ,nDOFX*nFine) )

    ALLOCATE( xiX1(nNodesX(1)) )
    ALLOCATE( xiX2(nNodesX(2)) )
    ALLOCATE( xiX3(nNodesX(3)) )

    ProjectionMatrix = 0.0_DP

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

      DO k = 1, nDOFX
      DO i = 1, nDOFX

        j = ( iFine - 1 ) * nDOFX + i

        DO iN3 = 1, nNodesX(3)
        DO iN2 = 1, nNodesX(2)
        DO iN1 = 1, nNodesX(1)

          ProjectionMatrix(j,k) &
            = ProjectionMatrix(j,k) &
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

    ProjectionMatrix_T &
      = TRANSPOSE( ProjectionMatrix )

    DO k = 1, nDOFX
      DO iFine = 1, nFine
        DO i = 1, nDOFX

          j = ( iFine - 1 ) * nDOFX + i

          ProjectionMatrix  (j,k) = ProjectionMatrix  (j,k) / WeightsX_Q(i)
          ProjectionMatrix_T(k,j) = ProjectionMatrix_T(k,j) / WeightsX_Q(k)

        END DO
      END DO
    END DO

    DEALLOCATE( xiX1 )
    DEALLOCATE( xiX2 )
    DEALLOCATE( xiX3 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ProjectionMatrix, ProjectionMatrix_T )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  ProjectionMatrix, ProjectionMatrix_T )
#endif

  END SUBROUTINE InitializeMeshRefinement_TwoMoment


  SUBROUTINE FinalizeMeshRefinement_TwoMoment

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ProjectionMatrix, ProjectionMatrix_T )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       ProjectionMatrix, ProjectionMatrix_T )
#endif

    DEALLOCATE( ProjectionMatrix )
    DEALLOCATE( ProjectionMatrix_T )

  END SUBROUTINE FinalizeMeshRefinement_TwoMoment


  SUBROUTINE RefineX_TwoMoment( nX, nVar, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nX(3), nVar
    REAL(DP), INTENT(in)  :: U_Crs(nDOFX,nX(1),nX(2),nX(3),nVar)
    REAL(DP), INTENT(out) :: U_Fin(nDOFX,nFine,nX(1),nX(2),nX(3),nVar)

    INTEGER :: nP_X

    nP_X = PRODUCT( nX ) * nVar

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    U_Crs, ProjectionMatrix ) &
    !$OMP MAP( alloc: U_Fin )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     U_Crs, Projectionmatrix ) &
    !$ACC CREATE(     U_Fin )
#endif

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX*nFine, nP_X, nDOFX, &
             One, ProjectionMatrix, nDOFX*nFine,  &
             U_Crs, nDOFX, &
             Zero, U_Fin, nDOFX*nFine )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_Fin ) &
    !$OMP MAP( release: U_Crs, ProjectionMatrix )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Fin ) &
    !$ACC DELETE(       U_Crs, ProjectionMatrix )
#endif

  END SUBROUTINE RefineX_TwoMoment


  SUBROUTINE CoarsenX_TwoMoment( nX, nVar, U_Fin, U_Crs )

    INTEGER,  INTENT(in)  :: nX(3), nVar
    REAL(DP), INTENT(in)  :: U_Fin(nDOFX,nFine,nX(1),nX(2),nX(3),nVar)
    REAL(DP), INTENT(out) :: U_Crs(nDOFX,nX(1),nX(2),nX(3),nVar)

    INTEGER :: nP_X

    nP_X = PRODUCT( nX ) * nVar

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    U_Fin, ProjectionMatrix_T ) &
    !$OMP MAP( alloc: U_Crs )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     U_Fin, Projectionmatrix_T ) &
    !$ACC CREATE(     U_Crs )
#endif

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX*nFine, &
             VolumeRatio, ProjectionMatrix_T, nDOFX,  &
             U_Fin, nDOFX*nFine, &
             Zero, U_Crs, nDOFX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_Crs ) &
    !$OMP MAP( release: U_Fin, ProjectionMatrix_T )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Crs ) &
    !$ACC DELETE(       U_Fin, ProjectionMatrix_T )
#endif

  END SUBROUTINE CoarsenX_TwoMoment


END MODULE TwoMoment_MeshRefinementModule
