MODULE TwoMoment_MeshRefinementModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two
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
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_TwoMoment
  PUBLIC :: FinalizeMeshRefinement_TwoMoment
  PUBLIC :: Refine_TwoMoment
  PUBLIC :: RefineX_TwoMoment
  PUBLIC :: RefineX_TwoMoment_Vector
  PUBLIC :: Coarsen_TwoMoment
  PUBLIC :: CoarsenX_TwoMoment
  PUBLIC :: CoarsenX_TwoMoment_Vector

  INTEGER  :: nFine, nFineX(3)
  REAL(DP), PUBLIC :: VolumeRatio
  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix  (:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix_T(:,:,:) ! --- Transpose ---

CONTAINS


  SUBROUTINE InitializeMeshRefinement_TwoMoment

    CHARACTER(2)  :: MeshString
    CHARACTER(32) :: MatrixName
    INTEGER :: iFine, iFineX1, iFineX2, iFineX3
    INTEGER :: i, k, iN1, iN2, iN3
    REAL(DP), ALLOCATABLE :: xiX1(:), xiX2(:), xiX3(:)

    nFineX      = 1
    VolumeRatio = One
    DO i = 1, nDimsX
      ! --- Refinement Factor of 2 Assumed:
      nFineX(i)   = 2
      VolumeRatio = VolumeRatio / Two
    END DO
    nFine = PRODUCT( nFineX )

    ALLOCATE( ProjectionMatrix  (nDOFX,nDOFX,nFine) )
    ALLOCATE( ProjectionMatrix_T(nDOFX,nDOFX,nFine) )

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

      ProjectionMatrix_T(:,:,iFine) &
        = TRANSPOSE( ProjectionMatrix(:,:,iFine) )

    END DO
    END DO
    END DO

    DEALLOCATE( xiX1 )
    DEALLOCATE( xiX2 )
    DEALLOCATE( xiX3 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA
    !$OMP MAP( to: ProjectionMatrix, ProjectionMatrix_T, VolumeRatio )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  ProjectionMatrix, ProjectionMatrix_T, VolumeRatio )
#endif

  END SUBROUTINE InitializeMeshRefinement_TwoMoment


  SUBROUTINE FinalizeMeshRefinement_TwoMoment

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA
    !$OMP MAP( release: ProjectionMatrix, ProjectionMatrix_T, VolumeRatio )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       ProjectionMatrix, ProjectionMatrix_T, VolumeRatio )
#endif

    DEALLOCATE( ProjectionMatrix )
    DEALLOCATE( ProjectionMatrix_T )

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


  FUNCTION RefineX_TwoMoment( iFine, U_Crs ) RESULT( U_Fin )

    INTEGER,  INTENT(in)  :: iFine
    REAL(DP), INTENT(in)  :: U_Crs(1:nDOFX)
    REAL(DP)              :: U_Fin(1:nDOFX)

    U_Fin = MATMUL( ProjectionMatrix(:,:,iFine), U_Crs ) / WeightsX_q

  END FUNCTION RefineX_TwoMoment


  SUBROUTINE RefineX_TwoMoment_Vector( nX, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nX(3)
    REAL(DP), INTENT(in)  :: U_Crs(nDOFX,nX(1),nX(2),nX(3))
    REAL(DP), INTENT(out) :: U_Fin(nDOFX,nX(1),nX(2),nX(3),nFine)

    INTEGER :: iNodeX, iX1, iX2, iX3
    INTEGER :: nP_X, iFine

    nP_X = PRODUCT( nX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA
    !$OMP MAP( to:    nX, U_Crs, ProjectionMatrix, WeightsX_Q ) &
    !$OMP MAP( alloc: U_Fin )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     nX, U_Crs, Projectionmatrix, WeightsX_Q ) &
    !$ACC CREATE(     U_Fin )
#endif

    DO iFine = 1, nFine
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX, nP_X, nDOFX, &
               One, ProjectionMatrix(:,:,iFine), nDOFX,  &
               U_Crs, nDOFX, &
               Zero, U_Fin(:,:,:,:,iFine), nDOFX )
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iFine = 1, nFine
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iNodeX = 1, nDOFX
              U_Fin(iNodeX,iX1,iX2,iX3,iFine) = U_Fin(iNodeX,iX1,iX2,iX3,iFine) / WeightsX_Q(iNodeX)
            END DO
          END DO
        END DO
      END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA
    !$OMP MAP( from:    U_Fin ) &
    !$OMP MAP( release: nX, U_Crs, ProjectionMatrix, WeightsX_Q )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Fin ) &
    !$ACC DELETE(       nX, U_Crs, ProjectionMatrix, WeightsX_Q )
#endif

  END SUBROUTINE RefineX_TwoMoment_Vector


  SUBROUTINE Coarsen_TwoMoment( nE, nX, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nE, nX(3)
    REAL(DP), INTENT(out) :: U_Crs(1:nDOF,1:nE)
    REAL(DP), INTENT(in)  :: U_Fin(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3))

    INTEGER :: iFine, iFineX1, iFineX2, iFineX3, iE
    REAL(DP) :: U_Crs_P(nDOFX,nDOFE*nE)
    REAL(DP) :: U_Fin_P(nDOFX,nDOFE*nE)

    U_Crs = 0.0_DP

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      U_Fin_P = Pack_TwoMoment( nE, U_Fin(:,:,iFineX1,iFineX2,iFineX3) )

      U_Crs_P(:,:) = MATMUL( ProjectionMatrix_T(:,:,iFine), U_Fin_P(:,:) )

      DO iE = 1, nDOFE * nE

        U_Crs_P(:,iE) = U_Crs_P(:,iE) / WeightsX_q(:)

      END DO

      U_Crs(:,:) = U_Crs(:,:) + VolumeRatio * Unpack_TwoMoment( nE, U_Crs_P )

    END DO
    END DO
    END DO

  END SUBROUTINE Coarsen_TwoMoment


  FUNCTION CoarsenX_TwoMoment( iFine, U_Fin ) RESULT( U_Crs )

    INTEGER,  INTENT(in)    :: iFine
    REAL(DP), INTENT(in)    :: U_Fin(1:nDOFX)
    REAL(DP)                :: U_Crs(1:nDOFX)

    U_Crs = VolumeRatio * MATMUL( ProjectionMatrix_T(:,:,iFine), U_Fin ) / WeightsX_q

  END FUNCTION CoarsenX_TwoMoment


  SUBROUTINE CoarsenX_TwoMoment_Vector( nX, U_Fin, U_Crs )

    INTEGER,  INTENT(in)  :: nX(3)
    REAL(DP), INTENT(in)  :: U_Fin(nDOFX,nX(1),nX(2),nX(3),nFine)
    REAL(DP), INTENT(out) :: U_Crs(nDOFX,nX(1),nX(2),nX(3))

    INTEGER :: iNodeX, iX1, iX2, iX3
    INTEGER :: nP_X, iFine

    nP_X = PRODUCT( nX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA
    !$OMP MAP( to:    nX, U_Fin, ProjectionMatrix_T, WeightsX_Q ) &
    !$OMP MAP( alloc: U_Crs )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     nX, U_Fin, Projectionmatrix_T, WeightsX_Q ) &
    !$ACC CREATE(     U_Crs )
#endif

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, &
             One, ProjectionMatrix_T(:,:,1), nDOFX,  &
             U_Fin(:,:,:,:,1), nDOFX, &
             Zero, U_Crs, nDOFX )
    DO iFine = 2, nFine
      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX, nP_X, nDOFX, &
               One, ProjectionMatrix_T(:,:,iFine), nDOFX,  &
               U_Fin(:,:,:,:,iFine), nDOFX, &
               One, U_Crs, nDOFX )
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX = 1, nDOFX
            U_Crs(iNodeX,iX1,iX2,iX3) = VolumeRatio * U_Crs(iNodeX,iX1,iX2,iX3) / WeightsX_Q(iNodeX)
          END DO
        END DO
      END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA
    !$OMP MAP( from:    U_Crs ) &
    !$OMP MAP( release: nX, U_Fin, ProjectionMatrix_T, WeightsX_Q )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Crs ) &
    !$ACC DELETE(       nX, U_Fin, ProjectionMatrix_T, WeightsX_Q )
#endif

  END SUBROUTINE CoarsenX_TwoMoment_Vector


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