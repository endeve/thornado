MODULE Euler_MeshRefinementModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Half
  USE ProgramHeaderModule, ONLY: &
    nDimsX, &
    nNodesX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, &
    NodesX2, &
    NodesX3, &
    WeightsX1, &
    WeightsX2, &
    WeightsX3, &
    WeightsX_q
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndLX_Q, &
    L_X1, &
    L_X2, &
    L_X3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_Euler
  PUBLIC :: FinalizeMeshRefinement_Euler
  PUBLIC :: Refine_Euler
  PUBLIC :: Coarsen_Euler

  INTEGER  :: nFine, nFineX(3)
  REAL(DP) :: VolumeRatio
  REAL(DP), ALLOCATABLE :: ProjectionMatrix  (:,:,:)
  REAL(DP), ALLOCATABLE :: ProjectionMatrix_T(:,:,:) ! --- Transpose ---


CONTAINS


  SUBROUTINE InitializeMeshRefinement_Euler

    INTEGER :: iDim
    INTEGER :: iFine, iFineX1, iFineX2, iFineX3
    INTEGER :: i, k, iN1, iN2, iN3
    REAL(DP), ALLOCATABLE :: xiX1(:), xiX2(:), xiX3(:)

    nFineX      = 1
    VolumeRatio = One
    DO iDim = 1, nDimsX
      ! --- Refinement Factor of 2 Assumed ---
      nFineX(iDim) = 2
      VolumeRatio  = Half * VolumeRatio
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

      ProjectionMatrix(:,:,iFine) = Zero
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

open(100,file='/Users/dunhamsj/Desktop/GaussLegendreWeights.txt')
write(100,'(A)') '{'
do in1=1,ndofx
if(.not.in1.eq.ndofx)then
write(100,'(2x,SP,ES24.16E3,A)') weightsx_q(in1),','
else
write(100,'(2x,SP,ES24.16E3)') weightsx_q(in1)
endif
enddo
write(100,'(A)') '};'
close(100)

open(100,file='/Users/dunhamsj/Desktop/ProjectionMatrix.txt')
write(100,'(A)') '{'
do ifine=1,nfine
write(100,'(A)') '  {'
do i=1,ndofx
write(100,'(A)') '    {'
do k=1,ndofx
if(.not.k.eq.ndofx)then
write(100,'(6x,SP,ES24.16E3,A)')ProjectionMatrix(i,k,iFine), ','
else
write(100,'(6x,SP,ES24.16E3)')ProjectionMatrix(i,k,iFine)
endif
enddo
if(.not.i.eq.ndofx)then
write(100,'(A)') '    },'
else
write(100,'(A)') '    }'
endif
enddo
if(.not.ifine.eq.4)then
write(100,'(A)') '  },'
else
write(100,'(A)') '  }'
endif
enddo
write(100,'(A)') '};'
close(100)

open(100,file='/Users/dunhamsj/Desktop/ProjectionMatrix_T.txt')
write(100,'(A)') '{'
do ifine=1,nfine
write(100,'(A)') '  {'
do i=1,ndofx
write(100,'(A)') '    {'
do k=1,ndofx
if(.not.k.eq.ndofx)then
write(100,'(6x,SP,ES24.16E3,A)')ProjectionMatrix_T(i,k,iFine), ','
else
write(100,'(6x,SP,ES24.16E3)')ProjectionMatrix_T(i,k,iFine)
endif
enddo
if(.not.i.eq.ndofx)then
write(100,'(A)') '    },'
else
write(100,'(A)') '    }'
endif
enddo
if(.not.ifine.eq.4)then
write(100,'(A)') '  },'
else
write(100,'(A)') '  }'
endif
enddo
write(100,'(A)') '};'
close(100)

print*,'  SHAPE( ProjMatrix ): ', shape(ProjectionMatrix)
print*,'SHAPE( ProjMatrix_T ): ', shape(ProjectionMatrix_T)

    DEALLOCATE( xiX1 )
    DEALLOCATE( xiX2 )
    DEALLOCATE( xiX3 )

  END SUBROUTINE InitializeMeshRefinement_Euler


  SUBROUTINE FinalizeMeshRefinement_Euler

    DEALLOCATE( ProjectionMatrix )
    DEALLOCATE( ProjectionMatrix_T )

  END SUBROUTINE FinalizeMeshRefinement_Euler


  SUBROUTINE Refine_Euler( nX, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nX(3)
    REAL(DP), INTENT(in)  :: U_Crs(1:nDOFX)
    REAL(DP), INTENT(out) :: U_Fin(1:nDOFX,1:nX(1),1:nX(2),1:nX(3))

    INTEGER :: iFine, iFineX1, iFineX2, iFineX3

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      U_Fin(:,iFineX1,iFineX2,iFineX3) &
        = MATMUL( ProjectionMatrix(:,:,iFine), U_Crs ) / WeightsX_q

    END DO
    END DO
    END DO

  END SUBROUTINE Refine_Euler


  SUBROUTINE Coarsen_Euler( nX, U_Crs, U_Fin )

    INTEGER,  INTENT(in)  :: nX(3)
    REAL(DP), INTENT(out) :: U_Crs(1:nDOFX)
    REAL(DP), INTENT(in)  :: U_Fin(1:nDOFX,1:nX(1),1:nX(2),1:nX(3))

    INTEGER  :: iFine, iFineX1, iFineX2, iFineX3
    REAL(DP) :: U_Crs_iFine(1:nDOFX)

    U_Crs = Zero

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      U_Crs_iFine = MATMUL( ProjectionMatrix_T(:,:,iFine), &
                            U_Fin(:,iFineX1,iFineX2,iFineX3) ) / WeightsX_q

      U_Crs = U_Crs + VolumeRatio * U_Crs_iFine

    END DO
    END DO
    END DO

  END SUBROUTINE Coarsen_Euler


END MODULE Euler_MeshRefinementModule
