MODULE TwoMoment_MeshRefinementModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX, nNodes
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_SpatialMetric
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    MatrixMatrixMultiplyBatched, &
    LinearSolveBatched

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_TwoMoment
  PUBLIC :: FinalizeMeshRefinement_TwoMoment
  PUBLIC :: RefineX_TwoMoment
  PUBLIC :: RefineX_TwoMoment_SIMPLE
  PUBLIC :: RefineX_TwoMoment_CURVILINEAR
  PUBLIC :: CoarsenX_TwoMoment
  PUBLIC :: CoarsenX_TwoMoment_SIMPLE
  PUBLIC :: CoarsenX_TwoMoment_CURVILINEAR

  INTEGER  :: nFine, nFineX(3)
  REAL(DP), PUBLIC :: VolumeRatio

  INTEGER :: nQuadX(3), nQ, nQuad

  LOGICAL :: UseSimpleMeshRefinement

  REAL(DP), ALLOCATABLE, PUBLIC :: RestrictionMatrixSimple (:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProlongationMatrixSimple(:,:,:)

  REAL(DP), ALLOCATABLE, PUBLIC :: RestrictionTensor0 (:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProlongationTensor0(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: MassTensor0        (:,:,:,:,:)

  INTERFACE RefineX_TwoMoment
    MODULE PROCEDURE RefineX_TwoMoment_SIMPLE
    MODULE PROCEDURE RefineX_TwoMoment_CURVILINEAR
  END INTERFACE RefineX_TwoMoment

  INTERFACE CoarsenX_TwoMoment
    MODULE PROCEDURE CoarsenX_TwoMoment_SIMPLE
    MODULE PROCEDURE CoarsenX_TwoMoment_CURVILINEAR
  END INTERFACE CoarsenX_TwoMoment

CONTAINS


  SUBROUTINE InitializeMeshRefinement_TwoMoment &
    ( UseSimpleMeshRefinement_Option, &
      Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: UseSimpleMeshRefinement_Option
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    INTEGER :: i, j, k
    INTEGER :: iFineX, iFineX1, iFineX2, iFineX3
    INTEGER :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX

    INTEGER :: iQuadX1, iQuadX2, iQuadX3
    REAL(DP), ALLOCATABLE :: QuadX1(:), QuadX2(:), QuadX3(:)
    REAL(DP), ALLOCATABLE :: QuadX1_Crse(:), QuadX2_Crse(:), QuadX3_Crse(:)
    REAL(DP), ALLOCATABLE :: WeightsX1(:), WeightsX2(:), WeightsX3(:)

    REAL(DP) :: FineXC_Crse(3), dFineX_Crse(3)

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    SELECT CASE ( TRIM( CoordinateSystem ) )
    CASE ( 'CARTESIAN' )
      nQuad = nNodes
      UseSimpleMeshRefinement = .TRUE.
    CASE ( 'CYLINDRICAL' )
      nQuad = nNodes
      UseSimpleMeshRefinement = .FALSE.
    CASE ( 'SPHERICAL' )
      nQuad = nNodes + 1
      UseSimpleMeshRefinement = .FALSE.
    CASE DEFAULT
      nQuad = nNodes
      UseSimpleMeshRefinement = .FALSE.
    END SELECT

    IF ( PRESENT( UseSimpleMeshRefinement_Option ) ) THEN
      UseSimpleMeshRefinement = UseSimpleMeshRefinement_Option
    END IF

    nQuadX      = 1
    nFineX      = 1
    VolumeRatio = One
    DO i = 1, nDimsX

      ! --- Number of quadrature points for integral
      nQuadX(i)  = nQuad

      ! --- Refinement Factor of 2 Assumed:
      nFineX(i)   = 2
      VolumeRatio = VolumeRatio / Two

    END DO
    nFine = PRODUCT( nFineX )
    nQ    = PRODUCT( nQuadX )

    IF ( UseSimpleMeshRefinement ) THEN
      CALL CreateMeshRefinement_TwoMoment_SIMPLE
    ELSE
      CALL CreateMeshRefinement_TwoMoment_CURVILINEAR
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: MeshRefinement_TwoMoment:'
      WRITE(*,'(A)') '  -------------------------------'
      WRITE(*,*)
      WRITE(*,'(A5,A32,L1)') &
        '', 'UseSimpleMeshRefinement: ', UseSimpleMeshRefinement
      WRITE(*,*)
      WRITE(*,'(A5,A36)') &
        '', 'Quadrature Points'
      WRITE(*,*)
      WRITE(*,'(A7,A9,I2.2)') &
        '', 'nQuad = ', nQuad
      WRITE(*,*)
      DO i = 1, nDimsX
        WRITE(*,'(A9,A4,I1,A10,I1,A4,I2.2)') &
          '', 'i = ', i, ', nQuadX(', i, ') = ', nQuadX(i)
      END DO
      WRITE(*,*)
      WRITE(*,'(A7,A9,I2.2)') &
        '', 'nFine = ', nFine
      WRITE(*,*)
      DO i = 1, nDimsX
        WRITE(*,'(A9,A4,I1,A10,I1,A4,I2.2)') &
          '', 'i = ', i, ', nFineX(', i, ') = ', nFineX(i)
      END DO
      WRITE(*,*)

    END IF

  END SUBROUTINE InitializeMeshRefinement_TwoMoment


  SUBROUTINE CreateMeshRefinement_TwoMoment_SIMPLE

    INTEGER :: i, j, k
    INTEGER :: iFineX, iFineX1, iFineX2, iFineX3
    INTEGER :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX

    INTEGER :: iQuad, iQuadX1, iQuadX2, iQuadX3
    REAL(DP), ALLOCATABLE :: QuadX1(:), QuadX2(:), QuadX3(:)
    REAL(DP), ALLOCATABLE :: QuadX1_Crse(:), QuadX2_Crse(:), QuadX3_Crse(:)
    REAL(DP), ALLOCATABLE :: WeightsX1(:), WeightsX2(:), WeightsX3(:), WeightsX_Q(:)

    REAL(DP) :: FineXC_Crse(3), dFineX_Crse(3)

    ALLOCATE( RestrictionMatrixSimple (nDOFX,nDOFX,nFine) )
    ALLOCATE( ProlongationMatrixSimple(nDOFX,nFine,nDOFX) )

    ProlongationMatrixSimple = 0.0_DP
    RestrictionMatrixSimple  = 0.0_DP

    ALLOCATE( QuadX1(nQuadX(1)) )
    ALLOCATE( QuadX2(nQuadX(2)) )
    ALLOCATE( QuadX3(nQuadX(3)) )

    ALLOCATE( WeightsX1(nQuadX(1)) )
    ALLOCATE( WeightsX2(nQuadX(2)) )
    ALLOCATE( WeightsX3(nQuadX(3)) )

    CALL GetQuadrature( nQuadX(1), QuadX1, WeightsX1 )
    CALL GetQuadrature( nQuadX(2), QuadX2, WeightsX2 )
    CALL GetQuadrature( nQuadX(3), QuadX3, WeightsX3 )

    ALLOCATE( WeightsX_Q(nQ) )

    DO iQuadX3 = 1, nQuadX(3)
    DO iQuadX2 = 1, nQuadX(2)
    DO iQuadX1 = 1, nQuadX(1)

      iQuad = (iQuadX3-1)*nQuadX(1)*nQuadX(2) &
            + (iQuadX2-1)*nQuadX(1) &
            +  iQuadX1

      WeightsX_Q(iQuad) &
        = WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3)

    END DO
    END DO
    END DO

    ALLOCATE( QuadX1_Crse(nQuadX(1)) )
    ALLOCATE( QuadX2_Crse(nQuadX(2)) )
    ALLOCATE( QuadX3_Crse(nQuadX(3)) )

    ! --- Fine element width in coarse element reference coordinates ---

    dFineX_Crse = One / DBLE( nFineX )

    ! --- Loop over fine elements comprising coarse element

    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFineX = (iFineX3-1)*nFineX(1)*nFineX(2) &
             + (iFineX2-1)*nFineX(1) &
             +  iFineX1

      ! --- Fine element center in coarse element reference coordinates ---

      FineXC_Crse(1) = - Half + dFineX_Crse(1) * ( DBLE(iFineX1) - Half )
      FineXC_Crse(2) = - Half + dFineX_Crse(2) * ( DBLE(iFineX2) - Half )
      FineXC_Crse(3) = - Half + dFineX_Crse(3) * ( DBLE(iFineX3) - Half )

      ! --- Fine element quadrature points in coarse element reference coordinates ---

      QuadX1_Crse = FineXC_Crse(1) + dFineX_Crse(1) * QuadX1
      QuadX2_Crse = FineXC_Crse(2) + dFineX_Crse(2) * QuadX2
      QuadX3_Crse = FineXC_Crse(3) + dFineX_Crse(3) * QuadX3

      DO jNodeX = 1, nDOFX
      DO iNodeX = 1, nDOFX

        DO iQuadX3 = 1, nQuadX(3)
        DO iQuadX2 = 1, nQuadX(2)
        DO iQuadX1 = 1, nQuadX(1)

          ProlongationMatrixSimple(iNodeX,iFineX,jNodeX) &
            = ProlongationMatrixSimple(iNodeX,iFineX,jNodeX) &
                + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                  * L_X1(IndLX_Q(1,iNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,iNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,iNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,jNodeX)) % P( QuadX1_Crse(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,jNodeX)) % P( QuadX2_Crse(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,jNodeX)) % P( QuadX3_Crse(iQuadX3) )

        END DO
        END DO
        END DO

        RestrictionMatrixSimple(jNodeX,iNodeX,iFineX) &
          = ProlongationMatrixSimple(iNodeX,iFineX,jNodeX)

        ProlongationMatrixSimple(iNodeX,iFineX,jNodeX) &
          = ProlongationMatrixSimple(iNodeX,iFineX,jNodeX) &
              / WeightsX_Q(iNodeX)

        RestrictionMatrixSimple(jNodeX,iNodeX,iFineX) &
          = RestrictionMatrixSimple(jNodeX,iNodeX,iFineX) &
              / WeightsX_Q(jNodeX)

      END DO
      END DO

    END DO
    END DO
    END DO

    DEALLOCATE( QuadX1 )
    DEALLOCATE( QuadX2 )
    DEALLOCATE( QuadX3 )

    DEALLOCATE( WeightsX1 )
    DEALLOCATE( WeightsX2 )
    DEALLOCATE( WeightsX3 )

    DEALLOCATE( WeightsX_Q )

    DEALLOCATE( QuadX1_Crse )
    DEALLOCATE( QuadX2_Crse )
    DEALLOCATE( QuadX3_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: RestrictionMatrixSimple, ProlongationMatrixSimple )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  RestrictionMatrixSimple, ProlongationMatrixSimple )
#endif

  END SUBROUTINE CreateMeshRefinement_TwoMoment_SIMPLE


  SUBROUTINE CreateMeshRefinement_TwoMoment_CURVILINEAR

    INTEGER :: i, j, k
    INTEGER :: iFineX, iFineX1, iFineX2, iFineX3
    INTEGER :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX

    INTEGER :: iQuadX1, iQuadX2, iQuadX3
    REAL(DP), ALLOCATABLE :: QuadX1(:), QuadX2(:), QuadX3(:)
    REAL(DP), ALLOCATABLE :: QuadX1_Crse(:), QuadX2_Crse(:), QuadX3_Crse(:)
    REAL(DP), ALLOCATABLE :: WeightsX1(:), WeightsX2(:), WeightsX3(:)

    REAL(DP) :: FineXC_Crse(3), dFineX_Crse(3)

    ALLOCATE( RestrictionTensor0 (nDOFX,nDOFX,nDOFX,nDOFX,nDOFX,nFine) )
    ALLOCATE( ProlongationTensor0(nDOFX,nDOFX,nDOFX,nDOFX,nDOFX,nFine) )
    ALLOCATE( MassTensor0        (nDOFX,nDOFX,nDOFX,nDOFX,nDOFX) )

    ALLOCATE( QuadX1(nQuadX(1)) )
    ALLOCATE( QuadX2(nQuadX(2)) )
    ALLOCATE( QuadX3(nQuadX(3)) )

    ALLOCATE( WeightsX1(nQuadX(1)) )
    ALLOCATE( WeightsX2(nQuadX(2)) )
    ALLOCATE( WeightsX3(nQuadX(3)) )

    CALL GetQuadrature( nQuadX(1), QuadX1, WeightsX1 )
    CALL GetQuadrature( nQuadX(2), QuadX2, WeightsX2 )
    CALL GetQuadrature( nQuadX(3), QuadX3, WeightsX3 )

    ALLOCATE( QuadX1_Crse(nQuadX(1)) )
    ALLOCATE( QuadX2_Crse(nQuadX(2)) )
    ALLOCATE( QuadX3_Crse(nQuadX(3)) )

    MassTensor0         = 0.0_DP
    ProlongationTensor0 = 0.0_DP
    RestrictionTensor0  = 0.0_DP

    ! --- Fine element width in coarse element reference coordinates ---

    dFineX_Crse = One / DBLE( nFineX )

    ! --- Loop over fine elements comprising coarse element

    DO iFineX3 = 1, nFineX(3)
    DO iFineX2 = 1, nFineX(2)
    DO iFineX1 = 1, nFineX(1)

      iFineX = (iFineX3-1)*nFineX(1)*nFineX(2) &
             + (iFineX2-1)*nFineX(1) &
             +  iFineX1

      ! --- Fine element center in coarse element reference coordinates ---

      FineXC_Crse(1) = - Half + dFineX_Crse(1) * ( DBLE(iFineX1) - Half )
      FineXC_Crse(2) = - Half + dFineX_Crse(2) * ( DBLE(iFineX2) - Half )
      FineXC_Crse(3) = - Half + dFineX_Crse(3) * ( DBLE(iFineX3) - Half )

      ! --- Fine element quadrature points in coarse element reference coordinates ---

      QuadX1_Crse = FineXC_Crse(1) + dFineX_Crse(1) * QuadX1
      QuadX2_Crse = FineXC_Crse(2) + dFineX_Crse(2) * QuadX2
      QuadX3_Crse = FineXC_Crse(3) + dFineX_Crse(3) * QuadX3

      DO jNodeX = 1, nDOFX
      DO iNodeX = 1, nDOFX

        DO mNodeX = 1, nDOFX
        DO lNodeX = 1, nDOFX
        DO kNodeX = 1, nDOFX

          DO iQuadX3 = 1, nQuadX(3)
          DO iQuadX2 = 1, nQuadX(2)
          DO iQuadX1 = 1, nQuadX(1)

            ProlongationTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) &
              = ProlongationTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) &
                  + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                    * L_X1(IndLX_Q(1,iNodeX)) % P( QuadX1     (iQuadX1) ) &
                    * L_X2(IndLX_Q(2,iNodeX)) % P( QuadX2     (iQuadX2) ) &
                    * L_X3(IndLX_Q(3,iNodeX)) % P( QuadX3     (iQuadX3) ) &
                    * L_X1(IndLX_Q(1,jNodeX)) % P( QuadX1_Crse(iQuadX1) ) &
                    * L_X2(IndLX_Q(2,jNodeX)) % P( QuadX2_Crse(iQuadX2) ) &
                    * L_X3(IndLX_Q(3,jNodeX)) % P( QuadX3_Crse(iQuadX3) ) &
                    * L_X1(IndLX_Q(1,kNodeX)) % P( QuadX1     (iQuadX1) ) &
                    * L_X2(IndLX_Q(2,kNodeX)) % P( QuadX2     (iQuadX2) ) &
                    * L_X3(IndLX_Q(3,kNodeX)) % P( QuadX3     (iQuadX3) ) &
                    * L_X1(IndLX_Q(1,lNodeX)) % P( QuadX1     (iQuadX1) ) &
                    * L_X2(IndLX_Q(2,lNodeX)) % P( QuadX2     (iQuadX2) ) &
                    * L_X3(IndLX_Q(3,lNodeX)) % P( QuadX3     (iQuadX3) ) &
                    * L_X1(IndLX_Q(1,mNodeX)) % P( QuadX1     (iQuadX1) ) &
                    * L_X2(IndLX_Q(2,mNodeX)) % P( QuadX2     (iQuadX2) ) &
                    * L_X3(IndLX_Q(3,mNodeX)) % P( QuadX3     (iQuadX3) )

            RestrictionTensor0(kNodeX,lNodeX,mNodeX,jNodeX,iNodeX,iFineX) &
              = ProlongationTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END DO
      END DO

    END DO
    END DO
    END DO

    DO jNodeX = 1, nDOFX
    DO iNodeX = 1, nDOFX

      DO mNodeX = 1, nDOFX
      DO lNodeX = 1, nDOFX
      DO kNodeX = 1, nDOFX

        DO iQuadX3 = 1, nQuadX(3)
        DO iQuadX2 = 1, nQuadX(2)
        DO iQuadX1 = 1, nQuadX(1)

          MassTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX) &
            = MassTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX) &
                + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                  * L_X1(IndLX_Q(1,iNodeX)) % P( QuadX1(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,iNodeX)) % P( QuadX2(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,iNodeX)) % P( QuadX3(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,jNodeX)) % P( QuadX1(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,jNodeX)) % P( QuadX2(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,jNodeX)) % P( QuadX3(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,kNodeX)) % P( QuadX1(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,kNodeX)) % P( QuadX2(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,kNodeX)) % P( QuadX3(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,lNodeX)) % P( QuadX1(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,lNodeX)) % P( QuadX2(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,lNodeX)) % P( QuadX3(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,mNodeX)) % P( QuadX1(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,mNodeX)) % P( QuadX2(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,mNodeX)) % P( QuadX3(iQuadX3) )

        END DO
        END DO
        END DO

      END DO
      END DO
      END DO

    END DO
    END DO

    DEALLOCATE( QuadX1 )
    DEALLOCATE( QuadX2 )
    DEALLOCATE( QuadX3 )

    DEALLOCATE( WeightsX1 )
    DEALLOCATE( WeightsX2 )
    DEALLOCATE( WeightsX3 )

    DEALLOCATE( QuadX1_Crse )
    DEALLOCATE( QuadX2_Crse )
    DEALLOCATE( QuadX3_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: RestrictionTensor0, ProlongationTensor0, MassTensor0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  RestrictionTensor0, ProlongationTensor0, MassTensor0 )
#endif

  END SUBROUTINE CreateMeshRefinement_TwoMoment_CURVILINEAR


  SUBROUTINE FinalizeMeshRefinement_TwoMoment

    IF ( UseSimpleMeshRefinement ) THEN
      CALL DestroyMeshRefinement_TwoMoment_SIMPLE
    ELSE
      CALL DestroyMeshRefinement_TwoMoment_CURVILINEAR
    END IF

  END SUBROUTINE FinalizeMeshRefinement_TwoMoment


  SUBROUTINE DestroyMeshRefinement_TwoMoment_SIMPLE

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: RestrictionMatrixSimple, ProlongationMatrixSimple )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       RestrictionMatrixSimple, ProlongationMatrixSimple )
#endif

    DEALLOCATE( ProlongationMatrixSimple )
    DEALLOCATE( RestrictionMatrixSimple )


  END SUBROUTINE DestroyMeshRefinement_TwoMoment_SIMPLE


  SUBROUTINE DestroyMeshRefinement_TwoMoment_CURVILINEAR

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: RestrictionTensor0, ProlongationTensor0, MassTensor0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       RestrictionTensor0, ProlongationTensor0, MassTensor0 )
#endif

    DEALLOCATE( MassTensor0 )
    DEALLOCATE( ProlongationTensor0 )
    DEALLOCATE( RestrictionTensor0 )

  END SUBROUTINE DestroyMeshRefinement_TwoMoment_CURVILINEAR


  SUBROUTINE RefineX_TwoMoment_SIMPLE( nX_Crse, nVar, U_Crse, U_Fine )

    INTEGER,  INTENT(in)  :: nX_Crse(3), nVar
    REAL(DP), INTENT(in)  :: U_Crse(nDOFX,      nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP), INTENT(out) :: U_Fine(nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    INTEGER  :: nCrse

    nCrse = PRODUCT( nX_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    U_Crse ) &
    !$OMP MAP( alloc: U_Fine )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     U_Crse ) &
    !$ACC CREATE(     U_Fine )
#endif

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX*nFine, nCrse*nVar, nDOFX, &
             One, ProlongationMatrixSimple, nDOFX*nFine,  &
             U_Crse, nDOFX, &
             Zero, U_Fine, nDOFX*nFine )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_Fine ) &
    !$OMP MAP( release: U_Crse )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Fine ) &
    !$ACC DELETE(       U_Crse )
#endif

  END SUBROUTINE RefineX_TwoMoment_SIMPLE


  SUBROUTINE CoarsenX_TwoMoment_SIMPLE( nX_Crse, nVar, U_Fine, U_Crse )

    INTEGER,  INTENT(in)  :: nX_Crse(3), nVar
    REAL(DP), INTENT(in)  :: U_Fine(nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP), INTENT(out) :: U_Crse(nDOFX,      nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    INTEGER  :: nCrse

    nCrse = PRODUCT( nX_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    U_Fine ) &
    !$OMP MAP( alloc: U_Crse )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     U_Fine ) &
    !$ACC CREATE(     U_Crse )
#endif

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCrse*nVar, nDOFX*nFine, &
             VolumeRatio, RestrictionMatrixSimple, nDOFX,  &
             U_Fine, nDOFX*nFine, &
             Zero, U_Crse, nDOFX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_Crse ) &
    !$OMP MAP( release: U_Fine )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_Crse ) &
    !$ACC DELETE(       U_Fine )
#endif

  END SUBROUTINE CoarsenX_TwoMoment_SIMPLE


  SUBROUTINE RefineX_TwoMoment_CURVILINEAR( nX_Crse, nX_Fine, nVar, G_Crse, U_Crse, G_Fine, U_Fine )

    INTEGER,  INTENT(in)  :: nX_Crse(3), nX_Fine(3), nVar
    REAL(DP), INTENT(in)  :: G_Crse(nDOFX,           nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF)
    REAL(DP), INTENT(in)  :: U_Crse(nDOFX,      nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP), INTENT(in)  :: G_Fine(nDOFX,           nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF)
    REAL(DP), INTENT(out) :: U_Fine(nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: ProlongationMatrix(nDOFX,nFine,nDOFX,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP) :: MassMatrix        (nDOFX,nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: IPIV              (nDOFX,      nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: INFO              (            nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: RHS               (nDOFX,nVar, nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP) :: SUM_M, SUM_P

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iFineX1, iFineX2, iFineX3, iFineX
    INTEGER  :: jX1, jX2, jX3
    INTEGER  :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX
    INTEGER  :: iVar
    INTEGER  :: nCrse

    nCrse = PRODUCT( nX_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine ) &
    !$OMP MAP( alloc: U_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO, RHS )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine ) &
    !$ACC CREATE(     U_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO, RHS )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_M, SUM_P, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_M, SUM_P, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_M, SUM_P, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#endif
    DO jX3 = 1, nX_Crse(3)
    DO jX2 = 1, nX_Crse(2)
    DO jX1 = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO jNodeX = 1, nDOFX
        DO iNodeX = 1, nDOFX

          iFineX3 = MOD( (iFineX-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1
          iFineX2 = MOD( (iFineX-1) / ( nFineX(1)             ), nFineX(2) ) + 1
          iFineX1 = MOD( (iFineX-1)                            , nFineX(1) ) + 1

          iX1 = ( jX1 - 1 ) * nFineX(1) + iFineX1
          iX2 = ( jX2 - 1 ) * nFineX(2) + iFineX2
          iX3 = ( jX3 - 1 ) * nFineX(3) + iFineX3

          SUM_M = Zero
          SUM_P = Zero
          DO mNodeX = 1, nDOFX
          DO lNodeX = 1, nDOFX
          DO kNodeX = 1, nDOFX

            CALL ComputeGeometryX_SpatialMetric &
                   ( G_Fine(kNodeX,iX1,iX2,iX3,iGF_h_1), &
                     G_Fine(lNodeX,iX1,iX2,iX3,iGF_h_2), &
                     G_Fine(mNodeX,iX1,iX2,iX3,iGF_h_3), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

            SUM_M = SUM_M + MassTensor0        (kNodeX,lNodeX,mNodeX,iNodeX,jNodeX       ) * SqrtGm
            SUM_P = SUM_P + ProlongationTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) * SqrtGm

          END DO
          END DO
          END DO

          MassMatrix        (iNodeX,jNodeX,iFineX,jX1,jX2,jX3) = SUM_M
          ProlongationMatrix(iNodeX,iFineX,jNodeX,jX1,jX2,jX3) = SUM_P

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiplyBatched &
           ( 'N', 'N', nDOFX*nFine, nVar, nDOFX, &
             One, ProlongationMatrix, nDOFX*nFine, nDOFX*nFine*nDOFX, &
             U_Crse, nDOFX, nDOFX*nVar, &
             Zero, U_Fine, nDOFX*nFine, nDOFX*nFine*nVar, &
             nCrse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO jX3 = 1, nX_Crse(3)
    DO jX2 = 1, nX_Crse(2)
    DO jX1 = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO iVar = 1, nVar
        DO iNodeX = 1, nDOFX

          RHS(iNodeX,iVar,iFineX,jX1,jX2,jX3) = U_Fine(iNodeX,iFineX,iVar,jX1,jX2,jX3)

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

    CALL LinearSolveBatched &
           ( 'N', nDOFX, nVar, &
             MassMatrix, nDOFX, IPIV, &
             RHS, nDOFX, INFO, nFine*nCrse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO jX3 = 1, nX_Crse(3)
    DO jX2 = 1, nX_Crse(2)
    DO jX1 = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO iVar = 1, nVar
        DO iNodeX = 1, nDOFX

          U_Fine(iNodeX,iFineX,iVar,jX1,jX2,jX3) = RHS(iNodeX,iVar,iFineX,jX1,jX2,jX3)

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( from:    U_Fine ) &
    !$OMP MAP( release: nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine, &
    !$OMP               ProlongationMatrix, MassMatrix, IPIV, INFO, RHS )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYOUT(      U_Fine ) &
    !$ACC DELETE(       nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine, &
    !$ACC               ProlongationMatrix, MassMatrix, IPIV, INFO, RHS )
#endif

  END SUBROUTINE RefineX_TwoMoment_CURVILINEAR


  SUBROUTINE CoarsenX_TwoMoment_CURVILINEAR( nX_Fine, nX_Crse, nVar, G_Fine, U_Fine, G_Crse, U_Crse )

    INTEGER,  INTENT(in)  :: nX_Fine(3), nX_Crse(3), nVar
    REAL(DP), INTENT(in)  :: G_Fine(nDOFX,           nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF)
    REAL(DP), INTENT(in)  :: U_Fine(nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP), INTENT(in)  :: G_Crse(nDOFX,           nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF)
    REAL(DP), INTENT(out) :: U_Crse(nDOFX,      nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: RestrictionMatrix(nDOFX,nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP) :: MassMatrix       (nDOFX,nDOFX,      nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: IPIV             (nDOFX,            nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: INFO             (                  nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP) :: SUM_M, SUM_R

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iFineX1, iFineX2, iFineX3, iFineX
    INTEGER  :: jX1, jX2, jX3
    INTEGER  :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX
    INTEGER  :: nCrse

    nCrse = PRODUCT( nX_Crse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    nX_Fine, nX_Crse, G_Fine, U_Fine, G_Crse ) &
    !$OMP MAP( alloc: U_Crse, RestrictionMatrix, MassMatrix, IPIV, INFO )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     nX_Fine, nX_Crse, G_Fine, U_Fine, G_Crse ) &
    !$ACC CREATE(     U_Crse, RestrictionMatrix, MassMatrix, IPIV, INFO )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_R, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_R, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1, iX2, iX3, SUM_R, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#endif
    DO jX3 = 1, nX_Crse(3)
    DO jX2 = 1, nX_Crse(2)
    DO jX1 = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO jNodeX = 1, nDOFX
        DO iNodeX = 1, nDOFX

          iFineX3 = MOD( (iFineX-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1
          iFineX2 = MOD( (iFineX-1) / ( nFineX(1)             ), nFineX(2) ) + 1
          iFineX1 = MOD( (iFineX-1)                            , nFineX(1) ) + 1

          iX1 = ( jX1 - 1 ) * nFineX(1) + iFineX1
          iX2 = ( jX2 - 1 ) * nFineX(2) + iFineX2
          iX3 = ( jX3 - 1 ) * nFineX(3) + iFineX3

          SUM_R = Zero
          DO mNodeX = 1, nDOFX
          DO lNodeX = 1, nDOFX
          DO kNodeX = 1, nDOFX

            CALL ComputeGeometryX_SpatialMetric &
                   ( G_Fine(kNodeX,iX1,iX2,iX3,iGF_h_1), &
                     G_Fine(lNodeX,iX1,iX2,iX3,iGF_h_2), &
                     G_Fine(mNodeX,iX1,iX2,iX3,iGF_h_3), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

            SUM_R = SUM_R + RestrictionTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) * SqrtGm

          END DO
          END DO
          END DO

          RestrictionMatrix(iNodeX,jNodeX,iFineX,jX1,jX2,jX3) = SUM_R

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( SUM_M, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( SUM_M, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( SUM_M, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#endif
    DO jX3 = 1, nX_Crse(3)
    DO jX2 = 1, nX_Crse(2)
    DO jX1 = 1, nX_Crse(1)

      DO jNodeX = 1, nDOFX
      DO iNodeX = 1, nDOFX

        SUM_M = Zero
        DO mNodeX = 1, nDOFX
        DO lNodeX = 1, nDOFX
        DO kNodeX = 1, nDOFX

          CALL ComputeGeometryX_SpatialMetric &
                 ( G_Crse(kNodeX,jX1,jX2,jX3,iGF_h_1), &
                   G_Crse(lNodeX,jX1,jX2,jX3,iGF_h_2), &
                   G_Crse(mNodeX,jX1,jX2,jX3,iGF_h_3), &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

          SUM_M = SUM_M + MassTensor0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX) * SqrtGm

        END DO
        END DO
        END DO

        MassMatrix(iNodeX,jNodeX,jX1,jX2,jX3) = SUM_M

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiplyBatched &
           ( 'N', 'N', nDOFX, nVar, nDOFX*nFine, &
             VolumeRatio, RestrictionMatrix, nDOFX, nDOFX*nDOFX*nFine, &
             U_Fine, nDOFX*nFine, nDOFX*nFine*nVar, &
             Zero, U_Crse, nDOFX, nDOFX*nVar, &
             nCrse )

    CALL LinearSolveBatched &
           ( 'N', nDOFX, nVar, &
             MassMatrix, nDOFX, IPIV, &
             U_Crse, nDOFX, INFO, nCrse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( from:    U_Crse ) &
    !$OMP MAP( release: nX_Fine, nX_Crse, G_Fine, U_Fine, G_Crse, RestrictionMatrix, MassMatrix, IPIV, INFO )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYOUT(      U_Crse ) &
    !$ACC DELETE(       nX_Fine, nX_Crse, G_Fine, U_Fine, G_Crse, RestrictionMatrix, MassMatrix, IPIV, INFO )
#endif

  END SUBROUTINE CoarsenX_TwoMoment_CURVILINEAR


END MODULE TwoMoment_MeshRefinementModule
