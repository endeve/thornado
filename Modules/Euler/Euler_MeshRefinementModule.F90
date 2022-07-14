MODULE Euler_MeshRefinementModule

#if defined( THORNADO_USE_AMREX ) && defined( REFINE_MESH )

  USE amrex_DGInterfaceModule, ONLY: &
    amrex_InitializeMeshRefinement_DG, &
    amrex_FinalizeMeshRefinement_DG

#endif

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Half
  USE ProgramHeaderModule, ONLY: &
    nDimsX, &
    nNodesX, &
    nDOFX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
     LX_X1_Dn, &
     LX_X1_Up, &
     LX_X2_Dn, &
     LX_X2_Up, &
     LX_X3_Dn, &
     LX_X3_Up
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, &
    NodesX2, &
    NodesX3, &
    WeightsX1, &
    WeightsX2, &
    WeightsX3, &
    WeightsX_q, &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    NodeNumberTableX_X1, &
    NodeNumberTableX_X2, &
    NodeNumberTableX_X3
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    IndLX_Q, &
    L_X1, &
    L_X2, &
    L_X3
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_Euler
  PUBLIC :: FinalizeMeshRefinement_Euler
  PUBLIC :: Refine_Euler
  PUBLIC :: Coarsen_Euler

  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X1_Refined(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X2_Refined(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X3_Refined(:,:,:)

  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X1_Refined_C(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X2_Refined_C(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X3_Refined_C(:)

  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X1_Up_1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X1_Dn_1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X2_Up_1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X2_Dn_1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X3_Up_1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: LX_X3_Dn_1D(:)

  REAL(DP), ALLOCATABLE :: xiX1(:)
  REAL(DP), ALLOCATABLE :: xiX2(:)
  REAL(DP), ALLOCATABLE :: xiX3(:)

  REAL(DP), ALLOCATABLE :: ProjectionMatrix  (:,:,:)
  REAL(DP), ALLOCATABLE :: ProjectionMatrix_c(:)
  REAL(DP), ALLOCATABLE :: ProjectionMatrix_T(:,:,:) ! --- Transpose ---

  INTEGER  :: nFine, nFineX(3)
  REAL(DP) :: VolumeRatio


CONTAINS


  SUBROUTINE InitializeMeshRefinement_Euler

    INTEGER :: iDim
    INTEGER :: iFine, iFineX1, iFineX2, iFineX3
    INTEGER :: i, k, iN1, iN2, iN3, kk, &
               iNX_X_Crse, iNX_X1_Crse, iNX_X2_Crse, iNX_X3_Crse, &
               iNX_X_Fine, iNX_X1_Fine, iNX_X2_Fine, iNX_X3_Fine
    REAL(DP), ALLOCATABLE :: xiX1(:), xiX2(:), xiX3(:)


    nFineX      = 1
    VolumeRatio = One
    DO iDim = 1, nDimsX
      ! --- Refinement Factor of 2 Assumed ---
      nFineX(iDim) = 2
      VolumeRatio  = Half * VolumeRatio
    END DO
    nFine = PRODUCT( nFineX )

    ALLOCATE( LX_X1_Refined(nDOFX_X1,nFineX(2)*nFineX(3),nDOFX_X1) )
    ALLOCATE( LX_X2_Refined(nDOFX_X2,nFineX(1)*nFineX(3),nDOFX_X2) )
    ALLOCATE( LX_X3_Refined(nDOFX_X3,nFineX(1)*nFineX(2),nDOFX_X3) )

    ALLOCATE( LX_X1_Refined_C(nDOFX_X1*nFineX(2)*nFineX(3)*nDOFX_X1) )
    ALLOCATE( LX_X2_Refined_C(nDOFX_X2*nFineX(1)*nFineX(3)*nDOFX_X2) )
    ALLOCATE( LX_X3_Refined_C(nDOFX_X3*nFineX(1)*nFineX(2)*nDOFX_X3) )

    ALLOCATE( LX_X1_Up_1D(nNodesX(1)) )
    ALLOCATE( LX_X1_Dn_1D(nNodesX(1)) )
    ALLOCATE( LX_X2_Up_1D(nNodesX(2)) )
    ALLOCATE( LX_X2_Dn_1D(nNodesX(2)) )
    ALLOCATE( LX_X3_Up_1D(nNodesX(3)) )
    ALLOCATE( LX_X3_Dn_1D(nNodesX(3)) )

    ALLOCATE( xiX1(nNodesX(1)) )
    ALLOCATE( xiX2(nNodesX(2)) )
    ALLOCATE( xiX3(nNodesX(3)) )

    ALLOCATE( ProjectionMatrix  (nDOFX,nDOFX,nFine) )
    ALLOCATE( ProjectionMatrix_c(nDOFX*nDOFX*nFine) )
    ALLOCATE( ProjectionMatrix_T(nDOFX,nDOFX,nFine) )

    DO i = 1, nNodesX(1)

      LX_X1_Up_1D(i) = L_X1(i) % P( +Half )
      LX_X1_Dn_1D(i) = L_X1(i) % P( -Half )

    END DO

    DO i = 1, nNodesX(2)

      LX_X2_Up_1D(i) = L_X2(i) % P( +Half )
      LX_X2_Dn_1D(i) = L_X2(i) % P( -Half )

    END DO

    DO i = 1, nNodesX(3)

      LX_X3_Up_1D(i) = L_X3(i) % P( +Half )
      LX_X3_Dn_1D(i) = L_X3(i) % P( -Half )

    END DO

    kk = 0

    iFine = 0
    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFine = iFine + 1

      IF( nFineX(1) .GT. 1 )THEN
        xiX1 = Half * ( NodesX1 + (-1)**iFineX1 * Half )
      ELSE
        xiX1 = Zero
      END IF

      IF( nFineX(2) .GT. 1 )THEN
        xiX2 = Half * ( NodesX2 + (-1)**iFineX2 * Half )
      ELSE
        xiX2 = Zero
      END IF

      IF( nFineX(3) .GT. 1 )THEN
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

        kk = kk + 1
        ProjectionMatrix_c(kk) = ProjectionMatrix(i,k,iFine)

      END DO
      END DO

      ProjectionMatrix_T(:,:,iFine) &
        = TRANSPOSE( ProjectionMatrix(:,:,iFine) )

    END DO
    END DO
    END DO

    kk = 0
    DO iNX_X1_Crse = 1, nDOFX_X1

      iNX_X2_Crse = NodeNumberTableX_X1(1,iNX_X1_Crse)
      iNX_X3_Crse = NodeNumberTableX_X1(2,iNX_X1_Crse)

      iFine = 0

      DO iFineX3 = 1, nFineX(3)
      DO iFineX2 = 1, nFineX(2)

        IF( nFineX(2) .GT. 1 )THEN
          xiX2 = Half * ( NodesX2 + (-1)**iFineX2 * Half )
        ELSE
          xiX2 = Zero
        END IF

        IF( nFineX(3) .GT. 1 )THEN
          xiX3 = Half * ( NodesX3 + (-1)**iFineX3 * Half )
        ELSE
          xiX3 = Zero
        END IF

        iFine = iFine + 1

        DO iNX_X1_Fine = 1, nDOFX_X1

          iNX_X2_Fine = NodeNumberTableX_X1(1,iNX_X1_Fine)
          iNX_X3_Fine = NodeNumberTableX_X1(2,iNX_X1_Fine)

          LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine) = One
          IF( nDimsX .GT. 1 ) &
            LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine) &
              = LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine) &
                  * Lagrange( xiX2(iNX_X2_Fine), iNX_X2_Crse, NodesX2 )
          IF( nDimsX .GT. 2 ) &
            LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine) &
              = LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine) &
                  * Lagrange( xiX3(iNX_X3_Fine), iNX_X3_Crse, NodesX3 )

          kk = kk + 1
          LX_X1_Refined_C(kk) &
            = LX_X1_Refined(iNX_X1_Crse,iFine,iNX_X1_Fine)

        END DO ! iNX_X1_Fine

      END DO ! iFineX2
      END DO ! iFineX3

    END DO ! iNX_X1_Crse

    IF( nDimsX .GT. 1 )THEN

      kk = 0
      DO iNX_X2_Crse = 1, nDOFX_X2

        iNX_X1_Crse = NodeNumberTableX_X2(1,iNX_X2_Crse)
        iNX_X3_Crse = NodeNumberTableX_X2(2,iNX_X2_Crse)

        iFine = 0

        DO iFineX3 = 1, nFineX(3)
        DO iFineX1 = 1, nFineX(1)

          IF( nFineX(1) .GT. 1 )THEN
            xiX1 = Half * ( NodesX1 + (-1)**iFineX1 * Half )
          ELSE
            xiX1 = Zero
          END IF

          IF( nFineX(3) .GT. 1 )THEN
            xiX3 = Half * ( NodesX3 + (-1)**iFineX3 * Half )
          ELSE
            xiX3 = Zero
          END IF

          iFine = iFine + 1

          DO iNX_X2_Fine = 1, nDOFX_X2

            iNX_X1_Fine = NodeNumberTableX_X2(1,iNX_X2_Fine)
            iNX_X3_Fine = NodeNumberTableX_X2(2,iNX_X2_Fine)

            LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine) = One
            IF( nDimsX .GT. 1 ) &
              LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine) &
                = LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine) &
                    * Lagrange( xiX1(iNX_X1_Fine), iNX_X1_Crse, NodesX1 )
            IF( nDimsX .GT. 2 ) &
              LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine) &
                = LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine) &
                    * Lagrange( xiX3(iNX_X3_Fine), iNX_X3_Crse, NodesX3 )

            kk = kk + 1
            LX_X2_Refined_C(kk) &
              = LX_X2_Refined(iNX_X2_Crse,iFine,iNX_X2_Fine)

          END DO ! iNX_X2_Fine

        END DO ! iFineX1
        END DO ! iFineX3

      END DO ! iNX_X2_Crse

    END IF ! nDimsX .GT. 1

    IF( nDimsX .GT. 2)THEN

      kk = 0
      DO iNX_X3_Crse = 1, nDOFX_X3

        iNX_X1_Crse = NodeNumberTableX_X3(1,iNX_X3_Crse)
        iNX_X2_Crse = NodeNumberTableX_X3(2,iNX_X3_Crse)

        iFine = 0

        DO iFineX2 = 1, nFineX(2)
        DO iFineX1 = 1, nFineX(1)

          IF( nFineX(1) .GT. 1 )THEN
            xiX1 = Half * ( NodesX1 + (-1)**iFineX1 * Half )
          ELSE
            xiX1 = Zero
          END IF

          IF( nFineX(2) .GT. 1 )THEN
            xiX2 = Half * ( NodesX2 + (-1)**iFineX2 * Half )
          ELSE
            xiX2 = Zero
          END IF

          iFine = iFine + 1

          DO iNX_X3_Fine = 1, nDOFX_X3

            iNX_X1_Fine = NodeNumberTableX_X3(1,iNX_X3_Fine)
            iNX_X2_Fine = NodeNumberTableX_X3(2,iNX_X3_Fine)

            LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine) = One
            IF( nDimsX .GT. 1 ) &
              LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine) &
                = LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine) &
                    * Lagrange( xiX1(iNX_X1_Fine), iNX_X1_Crse, NodesX1 )
            IF( nDimsX .GT. 2 ) &
              LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine) &
                = LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine) &
                    * Lagrange( xiX2(iNX_X2_Fine), iNX_X2_Crse, NodesX2 )

            kk = kk + 1
            LX_X3_Refined_C(kk) &
              = LX_X3_Refined(iNX_X3_Crse,iFine,iNX_X3_Fine)

          END DO ! iNX_X3_Fine

        END DO ! iFineX1
        END DO ! iFineX2

      END DO ! iNX_X3_Crse

    END IF ! nDimsX .GT. 2

#if defined( THORNADO_USE_AMREX ) && defined( REFINE_MESH )

    CALL amrex_InitializeMeshRefinement_DG &
           ( nNodesX, ProjectionMatrix_c, WeightsX1, WeightsX2, WeightsX3, &
             LX_X1_Refined_C, LX_X2_Refined_C, LX_X3_Refined_C, &
             LX_X1_Up_1D, LX_X1_Dn_1D, &
             LX_X2_Up_1D, LX_X2_Dn_1D, &
             LX_X3_Up_1D, LX_X3_Dn_1D, iGF_SqrtGm )

#endif

  END SUBROUTINE InitializeMeshRefinement_Euler


  SUBROUTINE FinalizeMeshRefinement_Euler

#ifdef THORNADO_USE_AMREX

    CALL amrex_FinalizeMeshRefinement_DG

#endif

    DEALLOCATE( ProjectionMatrix_T )
    DEALLOCATE( ProjectionMatrix )

    DEALLOCATE( xiX3 )
    DEALLOCATE( xiX2 )
    DEALLOCATE( xiX1 )

    DEALLOCATE( LX_X3_Dn_1D )
    DEALLOCATE( LX_X3_Up_1D )
    DEALLOCATE( LX_X2_Dn_1D )
    DEALLOCATE( LX_X2_Up_1D )
    DEALLOCATE( LX_X1_Dn_1D )
    DEALLOCATE( LX_X1_Up_1D )

    DEALLOCATE( LX_X3_Refined_C )
    DEALLOCATE( LX_X2_Refined_C )
    DEALLOCATE( LX_X1_Refined_C )

    DEALLOCATE( LX_X3_Refined )
    DEALLOCATE( LX_X2_Refined )
    DEALLOCATE( LX_X1_Refined )

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


  REAL(DP) FUNCTION Lagrange( x, i, xn ) RESULT( L )

    REAL(DP), INTENT(in) :: x
    INTEGER , INTENT(in) :: i
    REAL(DP), INTENT(in) :: xn(:)

    INTEGER :: j

    L = One
    DO j = 1, SIZE( xn )

      IF( j .NE. i ) L = L * ( x - xn(j) ) / ( xn(i) - xn(j) )

    END DO

    RETURN
  END FUNCTION Lagrange

END MODULE Euler_MeshRefinementModule
