PROGRAM MeshRefinementTest_TwoMoment

  USE KindModule, ONLY: &
    DP, One
  USE ProgramHeaderModule, ONLY: &
    nDOF, nE, nX, nDimsX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE TwoMoment_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_TwoMoment, &
    FinalizeMeshRefinement_TwoMoment, &
    Refine_TwoMoment, &
    Coarsen_TwoMoment

  IMPLICIT NONE

  INTEGER :: i, iNode, iX1, iX2, iX3, iE
  INTEGER :: iFine, nFine
  INTEGER :: iFineX1, iFineX2, iFineX3, nFineX(3)
  REAL(DP), ALLOCATABLE :: U_0(:,:,:,:,:) ! --- Coarse Data
  REAL(DP), ALLOCATABLE :: U_1(:,:,:,:,:,:) ! --- Fine Data

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'MeshRefinementTest_TwoMoment', &
           nX_Option &
             = [ 04, 04, 01 ], &
           swX_Option &
             = [ 01, 01, 01 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 10, &
           eL_Option &
             = 0.0d0, &
           eR_Option &
             = 1.0d2, &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  ! --- Energy Space Reference Element ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Refine / Coarsen ---

  CALL InitializeMeshRefinement_TwoMoment

  nFineX = 1
  DO i = 1, nDimsX
    nFineX(i) = 2
  END DO
  nFine = PRODUCT( nFineX )

  ALLOCATE( U_0(nDOF,nE,nX(1),nX(2),nX(3)) )
  ALLOCATE( U_1(nDOF,nE,nX(1),nX(2),nX(3),nFine) )

  ! --- Initialize Data on Coarse Level ---

  DO iX3 = 1, nX(3)
  DO iX2 = 1, nX(2)
  DO iX1 = 1, nX(1)
  DO iE  = 1, nE

    SELECT CASE ( nDimsX )

      CASE ( 1 )

        DO iNode = 1, nDOF

          U_0(iNode,iE,iX1,iX2,iX3) = One

        END DO

      CASE ( 2 )

        DO iNode = 1, nDOF

          U_0(iNode,iE,iX1,iX2,iX3) = One

        END DO

      CASE( 3 )

        DO iNode = 1, nDOF

          U_0(iNode,iE,iX1,iX2,iX3) = One

        END DO

    END SELECT

  END DO
  END DO
  END DO
  END DO

  PRINT*, "  SHAPE( U_0 ) = ", SHAPE( U_0 )
  PRINT*, "  SHAPE( U_1 ) = ", SHAPE( U_1 )
  PRINT*, "  nFine  = ", nFine
  PRINT*, "  nFineX = ", nFineX

  ! --- Refine ---

  iFine = 0
  DO iFineX3 = 1, nFineX(3)
  DO iFineX2 = 1, nFineX(2)
  DO iFineX1 = 1, nFineX(1)

    iFine = iFine + 1

    DO iX3 = 1, MAX( nX(3) / 2, 1 )
    DO iX2 = 1, MAX( nX(2) / 2, 1 )
    DO iX1 = 1, MAX( nX(1) / 2, 1 )

      CALL Refine_TwoMoment &
             ( nE, nFineX, U_0(:,:,iX1,iX2,iX3), &
               U_1(:,:,(iX1-1)*nFineX(1)+1:iX1*nFineX(1), &
                       (iX2-1)*nFineX(2)+1:iX2*nFineX(2), &
                       (iX3-1)*nFineX(3)+1:iX3*nFineX(3),iFine) )

    END DO
    END DO
    END DO

  END DO
  END DO
  END DO

  PRINT*, "After Refinement: "
  PRINT*, "  MIN/MAX U_0 = ", MINVAL( U_0 ), MAXVAL( U_0 )
  PRINT*, "  MIN/MAX U_1 = ", MINVAL( U_1 ), MAXVAL( U_1 )

  ! --- Coarsen ---

  iFine = 0
  DO iFineX3 = 1, nFineX(3)
  DO iFineX2 = 1, nFineX(2)
  DO iFineX1 = 1, nFineX(1)

    iFine = iFine + 1

    DO iX3 = 1, MAX( nX(3) / 2, 1 )
    DO iX2 = 1, MAX( nX(2) / 2, 1 )
    DO iX1 = 1, MAX( nX(1) / 2, 1 )

      CALL Coarsen_TwoMoment &
             ( nE, nFineX, U_0(:,:,iX1,iX2,iX3), &
               U_1(:,:,(iX1-1)*nFineX(1)+1:iX1*nFineX(1), &
                       (iX2-1)*nFineX(2)+1:iX2*nFineX(2), &
                       (iX3-1)*nFineX(3)+1:iX3*nFineX(3),iFine) )

    END DO
    END DO
    END DO

  END DO
  END DO
  END DO

  PRINT*, "After Coarsening: "
  PRINT*, "  MIN/MAX U_0 = ", MINVAL( U_0 ), MAXVAL( U_0 )
  PRINT*, "  MIN/MAX U_1 = ", MINVAL( U_1 ), MAXVAL( U_1 )

  CALL FinalizeMeshRefinement_TwoMoment

  DEALLOCATE( U_0, U_1 )

  ! --- Finalize ---

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM MeshRefinementTest_TwoMoment
