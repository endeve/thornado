MODULE SubcellReconstructionModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, iX_B0, iX_E0, iX_B1, iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, iGF_SqrtGm, iGF_h_1, iGF_h_2, iGF_h_3
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_SpatialMetric
  USE MeshModule, ONLY: &
    MeshX
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSubcellReconstruction
  PUBLIC :: FinalizeSubcellReconstruction
  PUBLIC :: CreateSubcellReconstruction
  PUBLIC :: DestroySubcellReconstruction

  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix0(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ReconstructionMatrix(:,:,:,:,:)

CONTAINS


  SUBROUTINE InitializeSubcellReconstruction

    INTEGER  :: iS1, iS2, iS3, iS
    INTEGER  :: jS1, jS2, jS3, jS
    INTEGER  :: kS1, kS2, kS3, kS
    INTEGER  :: lS1, lS2, lS3, lS
    INTEGER  :: mS1, mS2, mS3, mS
    INTEGER  :: iN1, iN2, iN3
    INTEGER  :: INFO
    INTEGER  :: IPIV(nDOFX)
    REAL(DP) :: WORK(nDOFX)
    REAL(DP) :: EtaLo(3), dEta(3)
    REAL(DP) :: Eta_X1(nNodesX(1))
    REAL(DP) :: Eta_X2(nNodesX(2))
    REAL(DP) :: Eta_X3(nNodesX(3))

    ALLOCATE( ProjectionMatrix0(nDOFX,nDOFX,nDOFX,nDOFX,nDOFX) )

    ProjectionMatrix0 = Zero

    ! --- Width of Subcells in Reference Element Coordinates ---

    dEta = One / DBLE( nNodesX )

    iS = 0
    DO iS3 = 1, nNodesX(3)
    DO iS2 = 1, nNodesX(2)
    DO iS1 = 1, nNodesX(1)

      iS = iS + 1

      ! --- Lower Coordinate of Subcells in Reference Element Coordinates ---

      EtaLo = - Half + dEta * DBLE( [ iS1, iS2, iS3 ] - 1 )

      ! --- Node Coordinates of Subcells in Reference Element Coordinates ---

      Eta_X1 = EtaLo(1) + dEta(1) * ( Half + NodesX1 )
      Eta_X2 = EtaLo(2) + dEta(2) * ( Half + NodesX2 )
      Eta_X3 = EtaLo(3) + dEta(3) * ( Half + NodesX3 )

      jS = 0
      DO jS3 = 1, nNodesX(3)
      DO jS2 = 1, nNodesX(2)
      DO jS1 = 1, nNodesX(1)

        jS = jS + 1

        kS = 0
        DO kS3 = 1, nNodesX(3)
        DO kS2 = 1, nNodesX(2)
        DO kS1 = 1, nNodesX(1)

        kS = kS + 1

        lS = 0
        DO lS3 = 1, nNodesX(3)
        DO lS2 = 1, nNodesX(2)
        DO lS1 = 1, nNodesX(1)

        lS = lS + 1

        mS = 0
        DO mS3 = 1, nNodesX(3)
        DO mS2 = 1, nNodesX(2)
        DO mS1 = 1, nNodesX(1)

        mS = mS + 1

        DO iN3 = 1, nNodesX(3)
        DO iN2 = 1, nNodesX(2)
        DO iN1 = 1, nNodesX(1)

          ProjectionMatrix0(kS,jS,lS,mS,iS) &
            = ProjectionMatrix0(kS,jS,lS,mS,iS) &
                + WeightsX1(iN1) * WeightsX2(iN2) * WeightsX3(iN3) &
                    * L_X1(IndLX_Q(1,jS)) % P( Eta_X1(iN1) ) &
                    * L_X2(IndLX_Q(2,jS)) % P( Eta_X2(iN2) ) &
                    * L_X3(IndLX_Q(3,jS)) % P( Eta_X3(iN3) ) &
                    * L_X1(IndLX_Q(1,kS)) % P( Eta_X1(iN1) ) &
                    * L_X2(IndLX_Q(2,kS)) % P( Eta_X2(iN2) ) &
                    * L_X3(IndLX_Q(3,kS)) % P( Eta_X3(iN3) ) &
                    * L_X1(IndLX_Q(1,lS)) % P( Eta_X1(iN1) ) &
                    * L_X2(IndLX_Q(2,lS)) % P( Eta_X2(iN2) ) &
                    * L_X3(IndLX_Q(3,lS)) % P( Eta_X3(iN3) ) &
                    * L_X1(IndLX_Q(1,mS)) % P( Eta_X1(iN1) ) &
                    * L_X2(IndLX_Q(2,mS)) % P( Eta_X2(iN2) ) &
                    * L_X3(IndLX_Q(3,mS)) % P( Eta_X3(iN3) )

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
        END DO

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ProjectionMatrix0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ProjectionMatrix0 )
#endif

  END SUBROUTINE InitializeSubcellReconstruction


  SUBROUTINE CreateSubcellReconstruction

    INTEGER  :: INFO
    INTEGER  :: IPIV(nDOFX)
    REAL(DP) :: WORK(nDOFX)

    INTEGER  :: iX1, iX2, iX3, iDOFX, jDOFX, kDOFX, lDOFX, mDOFX, iS
    REAL(DP) :: SUM1, SqrtGm, Gm_dd_11, Gm_dd_22, Gm_dd_33

    ALLOCATE( ProjectionMatrix    (nDOFX,nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3)) )
    ALLOCATE( ReconstructionMatrix(nDOFX,nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3)) )

    !CALL MatrixMatrixMultiply &
    !       ( 'T', 'N', nDOFX, nDOFX*nDOFX, nX_P, One, ProjectionMatrix0, nDOFX,
    !         SqrtGm, nDOFX, &
    !         Zero, ProjectionMatrix_T, nDOFX*nDOFX )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)


      DO iS = 1, nDOFX
      DO jDOFX = 1, nDOFX

        SUM1 = Zero
        DO kDOFX = 1, nDOFX
        DO lDOFX = 1, nDOFX
        DO mDOFX = 1, nDOFX

          CALL ComputeGeometryX_SpatialMetric &
                 ( uGF(kDOFX,iX1,iX2,iX3,iGF_h_1), &
                   uGF(lDOFX,iX1,iX2,iX3,iGF_h_2), &
                   uGF(mDOFX,iX1,iX2,iX3,iGF_h_3), &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )

          SUM1 = SUM1 + ProjectionMatrix0(jDOFX,kDOFX,lDOFX,mDOFX,iS) * SqrtGm

        END DO
        END DO
        END DO

        ProjectionMatrix(iS,jDOFX,iX1,iX2,iX3) = SUM1

      END DO
      END DO

      !CALL MatrixVectorMultiply &
      !       ( 'T', nDOFX, nDOFX*nDOFX, One, ProjectionMatrix0, &
      !         nDOFX, uGF(1,iX1,iX2,iX3,iGF_SqrtGm), 1, &
      !         Zero, ProjectionMatrix_T(1,1,iX1,iX2,iX3), 1 )

      !DO iS = 1, nDOFX

      !  CALL MatrixVectorMultiply &
      !         ( 'N', nDOFX, nDOFX, One, ProjectionMatrix0(1,1,iS), &
      !           nDOFX, SqrtGm, 1, Zero, ProjectionMatrix_T(1,iS), 1 )

      !END DO

      ReconstructionMatrix(:,:,iX1,iX2,iX3) = ProjectionMatrix(:,:,iX1,iX2,iX3)

      CALL DGETRF( nDOFX, nDOFX, ReconstructionMatrix(1,1,iX1,iX2,iX3), nDOFX, IPIV, INFO )

      IF( INFO /= 0 )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') '', 'InitializeSubcellReconstruction:'
        WRITE(*,*)
        WRITE(*,'(A6,A,I4.4)') '', 'DGETRF returned INFO = ', INFO
        WRITE(*,*)
        STOP
      END IF

      CALL DGETRI( nDOFX, ReconstructionMatrix(1,1,iX1,iX2,iX3), nDOFX, IPIV, WORK, nDOFX, INFO )

      IF( INFO /= 0 )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') '', 'InitializeSubcellReconstruction:'
        WRITE(*,*)
        WRITE(*,'(A6,A,I4.4)') '', 'DGETRI returned INFO = ', INFO
        WRITE(*,*)
        STOP
      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ProjectionMatrix, ReconstructionMatrix )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ProjectionMatrix, ReconstructionMatrix )
#endif

  END SUBROUTINE CreateSubcellReconstruction


  SUBROUTINE DestroySubcellReconstruction

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ProjectionMatrix, ReconstructionMatrix )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ProjectionMatrix, ReconstructionMatrix )
#endif

    DEALLOCATE( ProjectionMatrix )
    DEALLOCATE( ReconstructionMatrix )

  END SUBROUTINE DestroySubcellReconstruction


  SUBROUTINE FinalizeSubcellReconstruction

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ProjectionMatrix0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ProjectionMatrix0 )
#endif

    DEALLOCATE( ProjectionMatrix0 )

  END SUBROUTINE FinalizeSubcellReconstruction


END MODULE SubcellReconstructionModule
