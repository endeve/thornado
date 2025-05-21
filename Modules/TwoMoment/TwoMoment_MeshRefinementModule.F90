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
    MatrixMatrixMultiplyBatched, &
    LinearSolveBatched

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMeshRefinement_TwoMoment
  PUBLIC :: FinalizeMeshRefinement_TwoMoment
  PUBLIC :: CreateMeshRefinement_TwoMoment
  PUBLIC :: DestroyMeshRefinement_TwoMoment
  PUBLIC :: RefineX_TwoMoment
  PUBLIC :: CoarsenX_TwoMoment

  INTEGER  :: nFine, nFineX(3)
  REAL(DP), PUBLIC :: VolumeRatio

  INTEGER :: nQuadX(3), nQ, nQuad

  REAL(DP), ALLOCATABLE, PUBLIC :: RestrictionMatrix0 (:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProlongationMatrix0(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: MassMatrix0        (:,:,:,:,:)

CONTAINS


  SUBROUTINE InitializeMeshRefinement_TwoMoment &
    ( Verbose_Option )

    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    INTEGER :: iFineX, iFineX1, iFineX2, iFineX3

    INTEGER :: i, j, k
    INTEGER :: jCrseNodeX

    INTEGER :: iFineNodeX
    INTEGER :: jFineNodeX
    INTEGER :: kFineNodeX
    INTEGER :: lFineNodeX
    INTEGER :: mFineNodeX

    INTEGER :: iQuad, iQuadX1, iQuadX2, iQuadX3
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
    CASE ( 'CYLINDRICAL' )
      nQuad = nNodes
    CASE ( 'SPHERICAL' )
      nQuad = nNodes + 1
    CASE DEFAULT
      nQuad = nNodes
    END SELECT

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

    ALLOCATE( RestrictionMatrix0 (nDOFX,nDOFX,nDOFX,nDOFX,nDOFX,nFine) )
    ALLOCATE( ProlongationMatrix0(nDOFX,nDOFX,nDOFX,nDOFX,nDOFX,nFine) )
    ALLOCATE( MassMatrix0        (nDOFX,nDOFX,nDOFX,nDOFX,nDOFX) )

    ALLOCATE( QuadX1_Crse(nQuadX(1)) )
    ALLOCATE( QuadX2_Crse(nQuadX(2)) )
    ALLOCATE( QuadX3_Crse(nQuadX(3)) )

    ALLOCATE( QuadX1(nQuadX(1)) )
    ALLOCATE( QuadX2(nQuadX(2)) )
    ALLOCATE( QuadX3(nQuadX(3)) )

    ALLOCATE( WeightsX1(nQuadX(1)) )
    ALLOCATE( WeightsX2(nQuadX(2)) )
    ALLOCATE( WeightsX3(nQuadX(3)) )

    CALL GetQuadrature( nQuadX(1), QuadX1, WeightsX1 )
    CALL GetQuadrature( nQuadX(2), QuadX2, WeightsX2 )
    CALL GetQuadrature( nQuadX(3), QuadX3, WeightsX3 )

    ! --- Fine element width in coarse element reference coordinates ---

    dFineX_Crse = One / DBLE( nFineX )

    ! --- Loop over fine elements comprising coarse element

    MassMatrix0         = 0.0_DP
    ProlongationMatrix0 = 0.0_DP
    RestrictionMatrix0  = 0.0_DP

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

      DO jCrseNodeX = 1, nDOFX
      DO iFineNodeX = 1, nDOFX

      DO mFineNodeX = 1, nDOFX
      DO lFineNodeX = 1, nDOFX
      DO kFineNodeX = 1, nDOFX

        DO iQuadX3 = 1, nQuadX(3)
        DO iQuadX2 = 1, nQuadX(2)
        DO iQuadX1 = 1, nQuadX(1)

          ProlongationMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,iFineNodeX,jCrseNodeX,iFineX) &
            = ProlongationMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,iFineNodeX,jCrseNodeX,iFineX) &
                + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                  * L_X1(IndLX_Q(1,iFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,iFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,iFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,jCrseNodeX)) % P( QuadX1_Crse(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,jCrseNodeX)) % P( QuadX2_Crse(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,jCrseNodeX)) % P( QuadX3_Crse(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,kFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,kFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,kFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,lFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,lFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,lFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,mFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,mFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,mFineNodeX)) % P( QuadX3     (iQuadX3) )

          RestrictionMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,jCrseNodeX,iFineNodeX,iFineX) &
            = RestrictionMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,jCrseNodeX,iFineNodeX,iFineX) &
                + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                  * L_X1(IndLX_Q(1,iFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,iFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,iFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,jCrseNodeX)) % P( QuadX1_Crse(iQuadX1) ) &
                  * L_X2(IndLX_Q(2,jCrseNodeX)) % P( QuadX2_Crse(iQuadX2) ) &
                  * L_X3(IndLX_Q(3,jCrseNodeX)) % P( QuadX3_Crse(iQuadX3) ) &
                  * L_X1(IndLX_Q(1,kFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,kFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,kFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,lFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,lFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,lFineNodeX)) % P( QuadX3     (iQuadX3) ) &
                  * L_X1(IndLX_Q(1,mFineNodeX)) % P( QuadX1     (iQuadX1) ) &
                  * L_X2(IndLX_Q(2,mFineNodeX)) % P( QuadX2     (iQuadX2) ) &
                  * L_X3(IndLX_Q(3,mFineNodeX)) % P( QuadX3     (iQuadX3) )

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

    DO jFineNodeX = 1, nDOFX
    DO iFineNodeX = 1, nDOFX

    DO mFineNodeX = 1, nDOFX
    DO lFineNodeX = 1, nDOFX
    DO kFineNodeX = 1, nDOFX

      DO iQuadX3 = 1, nQuadX(3)
      DO iQuadX2 = 1, nQuadX(2)
      DO iQuadX1 = 1, nQuadX(1)

        MassMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,iFineNodeX,jFineNodeX) &
          = MassMatrix0(kFineNodeX,lFineNodeX,mFineNodeX,iFineNodeX,jFineNodeX) &
              + WeightsX1(iQuadX1) * WeightsX2(iQuadX2) * WeightsX3(iQuadX3) &
                * L_X1(IndLX_Q(1,iFineNodeX)) % P( QuadX1(iQuadX1) ) &
                * L_X2(IndLX_Q(2,iFineNodeX)) % P( QuadX2(iQuadX2) ) &
                * L_X3(IndLX_Q(3,iFineNodeX)) % P( QuadX3(iQuadX3) ) &
                * L_X1(IndLX_Q(1,jFineNodeX)) % P( QuadX1(iQuadX1) ) &
                * L_X2(IndLX_Q(2,jFineNodeX)) % P( QuadX2(iQuadX2) ) &
                * L_X3(IndLX_Q(3,jFineNodeX)) % P( QuadX3(iQuadX3) ) &
                * L_X1(IndLX_Q(1,kFineNodeX)) % P( QuadX1(iQuadX1) ) &
                * L_X2(IndLX_Q(2,kFineNodeX)) % P( QuadX2(iQuadX2) ) &
                * L_X3(IndLX_Q(3,kFineNodeX)) % P( QuadX3(iQuadX3) ) &
                * L_X1(IndLX_Q(1,lFineNodeX)) % P( QuadX1(iQuadX1) ) &
                * L_X2(IndLX_Q(2,lFineNodeX)) % P( QuadX2(iQuadX2) ) &
                * L_X3(IndLX_Q(3,lFineNodeX)) % P( QuadX3(iQuadX3) ) &
                * L_X1(IndLX_Q(1,mFineNodeX)) % P( QuadX1(iQuadX1) ) &
                * L_X2(IndLX_Q(2,mFineNodeX)) % P( QuadX2(iQuadX2) ) &
                * L_X3(IndLX_Q(3,mFineNodeX)) % P( QuadX3(iQuadX3) )

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    END DO
    END DO

    DEALLOCATE( QuadX1_Crse )
    DEALLOCATE( QuadX2_Crse )
    DEALLOCATE( QuadX3_Crse )

    DEALLOCATE( QuadX1 )
    DEALLOCATE( QuadX2 )
    DEALLOCATE( QuadX3 )

    DEALLOCATE( WeightsX1 )
    DEALLOCATE( WeightsX2 )
    DEALLOCATE( WeightsX3 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: RestrictionMatrix0, ProlongationMatrix0, MassMatrix0 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  RestrictionMatrix0, ProlongationMatrix0, MassMatrix0 )
#endif

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: MeshRefinement_TwoMoment:'
      WRITE(*,'(A)') '  -------------------------------'
      WRITE(*,*)
      WRITE(*,'(A5,A36)') '', 'Quadrature Points'
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


  SUBROUTINE CreateMeshRefinement_TwoMoment


  END SUBROUTINE CreateMeshRefinement_TwoMoment


  SUBROUTINE DestroyMeshRefinement_TwoMoment


  END SUBROUTINE DestroyMeshRefinement_TwoMoment


  SUBROUTINE FinalizeMeshRefinement_TwoMoment

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: RestrictionMatrix0, ProlongationMatrix0, MassMatrix0 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       RestrictionMatrix0, ProlongationMatrix0, MassMatrix0 )
#endif

    DEALLOCATE( MassMatrix0 )
    DEALLOCATE( ProlongationMatrix0 )
    DEALLOCATE( RestrictionMatrix0 )

  END SUBROUTINE FinalizeMeshRefinement_TwoMoment


  SUBROUTINE RefineX_TwoMoment( nX_Crse, nX_Fine, nVar, G_Crse, U_Crse, G_Fine, U_Fine )

    INTEGER,  INTENT(in)  :: nX_Crse(3), nX_Fine(3), nVar
    REAL(DP), INTENT(in)  :: G_Crse(nDOFX,           nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF)
    REAL(DP), INTENT(in)  :: U_Crse(nDOFX,nVar,      nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP), INTENT(in)  :: G_Fine(nDOFX,           nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF)
    REAL(DP), INTENT(out) :: U_Fine(nDOFX,nVar,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: ProlongationMatrix(nDOFX,nFine,nDOFX,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    REAL(DP) :: MassMatrix        (nDOFX,nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: IPIV              (nDOFX,      nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))
    INTEGER  :: INFO              (            nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: RHS               (nDOFX,nFine,nVar ,nX_Crse(1),nX_Crse(2),nX_Crse(3))

    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP) :: SUM_M, SUM_P

    INTEGER  :: iX1_Fine, iX2_Fine, iX3_Fine
    INTEGER  :: iFineX1, iFineX2, iFineX3, iFineX
    INTEGER  :: jX1_Crse, jX2_Crse, jX3_Crse
    INTEGER  :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX
    INTEGER  :: iVar
    INTEGER  :: nP_X, nCrse

    nCrse = PRODUCT( nX_Crse )
    nP_X = nCrse * nVar

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine ) &
    !$OMP MAP( alloc: U_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine ) &
    !$ACC CREATE(     U_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_M, SUM_P, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_M, SUM_P, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_M, SUM_P, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#endif
    DO jX3_Crse = 1, nX_Crse(3)
    DO jX2_Crse = 1, nX_Crse(2)
    DO jX1_Crse = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO jNodeX = 1, nDOFX
        DO iNodeX = 1, nDOFX

          iFineX3 = MOD( (iFineX-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1
          iFineX2 = MOD( (iFineX-1) / ( nFineX(1)             ), nFineX(2) ) + 1
          iFineX1 = MOD( (iFineX-1)                            , nFineX(1) ) + 1

          iX1_Fine = ( jX1_Crse - 1 ) * nFineX(1) + iFineX1
          iX2_Fine = ( jX2_Crse - 1 ) * nFineX(2) + iFineX2
          iX3_Fine = ( jX3_Crse - 1 ) * nFineX(3) + iFineX3

          SUM_M = Zero
          SUM_P = Zero
          DO mNodeX = 1, nDOFX
          DO lNodeX = 1, nDOFX
          DO kNodeX = 1, nDOFX

            CALL ComputeGeometryX_SpatialMetric &
                   ( G_Fine(kNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_1), &
                     G_Fine(lNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_2), &
                     G_Fine(mNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_3), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

            SUM_M = SUM_M + MassMatrix0        (kNodeX,lNodeX,mNodeX,iNodeX,jNodeX       ) * SqrtGm
            SUM_P = SUM_P + ProlongationMatrix0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) * SqrtGm

          END DO
          END DO
          END DO

          MassMatrix        (iNodeX,jNodeX,iFineX,jX1_Crse,jX2_Crse,jX3_Crse) = SUM_M
          ProlongationMatrix(iNodeX,iFineX,jNodeX,jX1_Crse,jX2_Crse,jX3_Crse) = SUM_P

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
             Zero, RHS, nDOFX*nFine, nDOFX*nFine*nVar, &
             nCrse )

    DO jX3_Crse = 1, nX_Crse(3)
    DO jX2_Crse = 1, nX_Crse(2)
    DO jX1_Crse = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO iVar = 1, nVar
        DO iNodeX = 1, nDOFX

          U_Fine(iNodeX,iVar,iFineX,jX1_Crse,jX2_Crse,jX3_Crse) = RHS(iNodeX,iFineX,iVar,jX1_Crse,jX2_Crse,jX3_Crse)

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

    CALL LinearSolveBatched &
           ( 'N', nDOFX, nVar, &
             MassMatrix, nDOFX, IPIV, &
             U_Fine, nDOFX, INFO, nFine*nCrse )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( from:    U_Fine ) &
    !$OMP MAP( release: nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYOUT(      U_Fine ) &
    !$ACC DELETE(       nX_Crse, nX_Fine, G_Crse, U_Crse, G_Fine, ProlongationMatrix, MassMatrix, IPIV, INFO )
#endif

  END SUBROUTINE RefineX_TwoMoment


  SUBROUTINE CoarsenX_TwoMoment( nX_Fine, nX_Crse, nVar, G_Fine, U_Fine, G_Crse, U_Crse )

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

    INTEGER  :: iX1_Fine, iX2_Fine, iX3_Fine
    INTEGER  :: iFineX1, iFineX2, iFineX3, iFineX
    INTEGER  :: jX1_Crse, jX2_Crse, jX3_Crse
    INTEGER  :: iNodeX, jNodeX, kNodeX, lNodeX, mNodeX
    INTEGER  :: nP_X, nCrse

    nCrse = PRODUCT( nX_Crse )
    nP_X = nCrse * nVar

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
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_R, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_R, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iFineX1, iFineX2, iFineX3, iX1_Fine, iX2_Fine, iX3_Fine, SUM_R, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#endif
    DO jX3_Crse = 1, nX_Crse(3)
    DO jX2_Crse = 1, nX_Crse(2)
    DO jX1_Crse = 1, nX_Crse(1)

      DO iFineX = 1, nFine

        DO jNodeX = 1, nDOFX
        DO iNodeX = 1, nDOFX

          iFineX3 = MOD( (iFineX-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1
          iFineX2 = MOD( (iFineX-1) / ( nFineX(1)             ), nFineX(2) ) + 1
          iFineX1 = MOD( (iFineX-1)                            , nFineX(1) ) + 1

          iX1_Fine = ( jX1_Crse - 1 ) * nFineX(1) + iFineX1
          iX2_Fine = ( jX2_Crse - 1 ) * nFineX(2) + iFineX2
          iX3_Fine = ( jX3_Crse - 1 ) * nFineX(3) + iFineX3

          SUM_R = Zero
          DO mNodeX = 1, nDOFX
          DO lNodeX = 1, nDOFX
          DO kNodeX = 1, nDOFX

            CALL ComputeGeometryX_SpatialMetric &
                   ( G_Fine(kNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_1), &
                     G_Fine(lNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_2), &
                     G_Fine(mNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_h_3), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

            SUM_R = SUM_R + RestrictionMatrix0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX,iFineX) * SqrtGm

          END DO
          END DO
          END DO

          RestrictionMatrix(iNodeX,jNodeX,iFineX,jX1_Crse,jX2_Crse,jX3_Crse) = SUM_R

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
    DO jX3_Crse = 1, nX_Crse(3)
    DO jX2_Crse = 1, nX_Crse(2)
    DO jX1_Crse = 1, nX_Crse(1)

      DO jNodeX = 1, nDOFX
      DO iNodeX = 1, nDOFX

        SUM_M = Zero
        DO mNodeX = 1, nDOFX
        DO lNodeX = 1, nDOFX
        DO kNodeX = 1, nDOFX

          CALL ComputeGeometryX_SpatialMetric &
                 ( G_Crse(kNodeX,jX1_Crse,jX2_Crse,jX3_Crse,iGF_h_1), &
                   G_Crse(lNodeX,jX1_Crse,jX2_Crse,jX3_Crse,iGF_h_2), &
                   G_Crse(mNodeX,jX1_Crse,jX2_Crse,jX3_Crse,iGF_h_3), &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

          SUM_M = SUM_M + MassMatrix0(kNodeX,lNodeX,mNodeX,iNodeX,jNodeX) * SqrtGm

        END DO
        END DO
        END DO

        MassMatrix(iNodeX,jNodeX,jX1_Crse,jX2_Crse,jX3_Crse) = SUM_M

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

  END SUBROUTINE CoarsenX_TwoMoment


END MODULE TwoMoment_MeshRefinementModule
