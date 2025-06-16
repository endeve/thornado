MODULE CellMergingModule

    USE KindModule, ONLY: &
      DP, Zero, Half, One, Two, Pi
    USE ProgramHeaderModule, ONLY: & ! Work on replacing input arguments with public variables from ProgramHeader
      iX_B0, iX_E0, iX_B1, iX_E1, nDOFX, nNodesX, nDims
    USE GeometryFieldsModule, ONLY: &
      uGF, iGF_h_2, iGF_h_3, iGF_SqrtGm, nGF ! In spherical coordinates: iGF_h_2 => r, iGF_h_3 => r*sin\theta
    USE FluidFieldsModule, ONLY: &
      nCF
    USE QuadratureModule, ONLY: &
      GetQuadrature
    USE ReferenceElementModuleX, ONLY: &
      NodeNumberTableX, NodeNumberTableX3D
    USE MeshModule, ONLY: &
      MeshX, &
      NodeCoordinate
    USE PolynomialBasisModuleX_Lagrange, ONLY: &
      L_X1, L_X2, L_X3
    USE TimersModule_Euler, ONLY: &
      TimersStart_Euler, TimersStop_Euler, &
      Timer_CellMerging, Timer_CM_UpdateCoefficient

    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: SqrtGm_Merge

    TYPE, PUBLIC :: MergedMeshType
      INTEGER                                 :: NCellsPerMerge
      INTEGER                                 :: NCells
      INTEGER,  DIMENSION(:),     ALLOCATABLE :: FineCellMarker
      INTEGER,  DIMENSION(:),     ALLOCATABLE :: MergeCellMarker
      REAL(DP), DIMENSION(:),     ALLOCATABLE :: MergeWidth
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: MergedBasisCoeff
    END type MergedMeshType
    
    TYPE(MergedMeshType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: MergedMeshX2 ! X2 merging
    TYPE(MergedMeshType), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: MergedMeshX3 ! X3 merging

    PUBLIC :: Initialize_CellMerging
    PUBLIC :: Finalize_CellMerging
    PUBLIC :: Determine_MergedCellsX2
    PUBLIC :: Determine_MergedCellsX3
    PUBLIC :: Determine_BasisCoeffX2
    PUBLIC :: Determine_BasisCoeffX3
    PUBLIC :: Determine_MergedGeometry2D
    PUBLIC :: Determine_MergedGeometry3D
    PUBLIC :: MergeAndRestrict2D
    PUBLIC :: MergeAndRestrict3D
    PUBLIC :: MergeAndRestrict
    PUBLIC :: MergeAndRestrictGeometry ! Delete this subroutine

CONTAINS

    SUBROUTINE Initialize_CellMerging ( nX, nN, Min_NCellsPerMerge_Option )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN
      INTEGER, INTENT(in), OPTIONAL :: Min_NCellsPerMerge_Option

      INTEGER :: iX1, iX2
      INTEGER :: Min_NCellsPerMerge

      IF( PRESENT( Min_NCellsPerMerge_Option ) )THEN

        Min_NCellsPerMerge = Min_NCellsPerMerge_Option

      ELSE

        Min_NCellsPerMerge = 1

      ENDIF

      ! --- Create MergedMesh for X2 direction ---
      
      ALLOCATE( MergedMeshX2(1:nX(1)) )

      MergedMeshX2 % NCellsPerMerge = 1
      MergedMeshX2 % NCells         = nX(2)

      DO iX1 = 1, nX(1)

        ALLOCATE( MergedMeshX2(iX1) % MergeWidth     (1:nX(2)) )
        ALLOCATE( MergedMeshX2(iX1) % FineCellMarker (1:nX(2)) )
        ALLOCATE( MergedMeshX2(iX1) % MergeCellMarker(1:nX(2)) )

      END DO

      CALL Determine_MergedCellsX2( nX, Min_NCellsPerMerge )

      DO iX1 = 1, nX(1)

        ALLOCATE( MergedMeshX2(iX1) % &
                  MergedBasisCoeff(1:nN, &
                                   1:nN, &
                                   1:MergedMeshX2(iX1) % NCellsPerMerge))

      END DO

      CALL Determine_BasisCoeffX2( nX, nN )

      ! --- Create MergedMesh for X3 direction ---
      IF( nDims .GT. 2 )THEN

      ALLOCATE( MergedMeshX3(1:nX(1),1:nX(2)) )

      MergedMeshX3 % NCellsPerMerge = 1
      MergedMeshX3 % NCells         = nX(3)

      DO iX1 = 1, nX(1)
      DO iX2 = 1, nX(2)

        ALLOCATE( MergedMeshX3(iX1,iX2) % MergeWidth     (1:nX(3)) )
        ALLOCATE( MergedMeshX3(iX1,iX2) % FineCellMarker (1:nX(3)) )
        ALLOCATE( MergedMeshX3(iX1,iX2) % MergeCellMarker(1:nX(3)) )

      END DO
      END DO

      CALL Determine_MergedCellsX3( nX, Min_NCellsPerMerge )

      DO iX1 = 1, nX(1)
      DO iX2 = 1, nX(2)

        ALLOCATE( MergedMeshX3(iX1,iX2) % &
                  MergedBasisCoeff(1:nN, &
                                   1:nN, &
                                   1:MergedMeshX3(iX1,iX2) % NCellsPerMerge))

      END DO
      END DO

      CALL Determine_BasisCoeffX3( nX, nN )

      END IF

      ! --- Compute SqrtGm on merged grid ---

      ALLOCATE( SqrtGm_Merge(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) )

      IF( nDims .LT. 3 )THEN

        CALL Determine_MergedGeometry2D( nN )

      ELSE

        CALL Determine_MergedGeometry3D( nN )

      END IF

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging ( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1

      DEALLOCATE( SqrtGm_Merge )

      IF( nDims .GT. 2)THEN

      DO iX1 = 1, nX(1)
      DO iX2 = 1, nX(2)

        DEALLOCATE( MergedMeshX3(iX1,iX2) % MergedBasisCoeff )
        DEALLOCATE( MergedMeshX3(iX1,iX2) % MergeCellMarker  )
        DEALLOCATE( MergedMeshX3(iX1,iX2) % FineCellMarker   )
        DEALLOCATE( MergedMeshX3(iX1,iX2) % MergeWidth       )

      END DO
      END DO

      DEALLOCATE( MergedMeshX3 )

      END IF

      DO iX1 = 1, nX(1)

        DEALLOCATE( MergedMeshX2(iX1) % MergedBasisCoeff )
        DEALLOCATE( MergedMeshX2(iX1) % MergeCellMarker  )
        DEALLOCATE( MergedMeshX2(iX1) % FineCellMarker   )
        DEALLOCATE( MergedMeshX2(iX1) % MergeWidth       )

      END DO

      DEALLOCATE( MergedMeshX2 )

    END SUBROUTINE Finalize_CellMerging

    SUBROUTINE Determine_MergedCellsX2 ( nX, MinCells )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: MinCells

      REAL(DP) :: r(1:nX(1))
      INTEGER  :: n, iX1, iX2

      ! --- Merging is assumed to happen in the X2 direction
      !     such that the total number of cells in the X2 direction
      !     is divisible by 2. Equidistant meshes are also assumed ---

      ASSOCIATE( dX1 => MeshX(1) % Width, &
                 dX2 => MeshX(2) % Width )

      r = MeshX(1) % Center(1:nX(1)) + 0.5_DP * dX1(1:nX(1)) ! It's possible to replace r with uGF?

      DO iX1 = 1, nX(1)

        ASSOCIATE( NCellsX2          => MergedMeshX2(iX1) % NCells, &
                   NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker )

        n = 0

        ! --- Determine how many cells per merge in the X2 direction
        !     and the number of merged cells in the X2 direction ---
        IF ( r(iX1) * dX2(1) .LE. dX1(iX1) / REAL(2**(MinCells-1),DP) ) THEN

          DO WHILE ( REAL(2**n,DP) * r(iX1) * dX2(1) .LE. dX1(iX1) / REAL(2**(MinCells-1),DP) )

            n = n + 1

          END DO

          ! NCellsX2         = MAX(nX(2) / 2**n, 2**MinCells)
          NCellsX2         = nX(2) / 2**n
          ! print *, nX(2), n, 2**n, NCellsX2
          NCellsPerMergeX2 = nX(2) / NCellsX2

        ELSE

          NCellsPerMergeX2 = 1
          NCellsX2         = nX(2)

        END IF

        ! print *, 'NCellsPerMergeX2 = ', NCellsPerMergeX2, 'NCellsX2 = ', NCellsX2

        ! --- Assumes an equidsitant mesh ---
        ! --- Mark the index of the first fine cell in the merged cell ---
        DO iX2 = 1, nX(2), NCellsPerMergeX2

          dX2M(iX2:iX2+NCellsPerMergeX2-1) &
            = dX2(iX2) * REAL(NCellsPerMergeX2,DP)

          MergeCellMarkerX2(iX2:iX2+NCellsPerMergeX2-1) = iX2

        END DO

        ! --- Mark the cells which are merged according to fine cell # within coarse cell ---
        DO iX2 = 1, nX(2)

          FineCellMarkerX2(iX2) = MOD(iX2 - 1, NCellsPerMergeX2) + 1

        END DO

        END ASSOCIATE ! NCellsX2, NCellsPerMergeX2, etc.

      END DO

      END ASSOCIATE ! dX1, dX2

      ! STOP

    END SUBROUTINE Determine_MergedCellsX2

    SUBROUTINE Determine_MergedCellsX3 ( nX, MinCells )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: MinCells

      REAL(DP) :: r(1:nX(1)), thetaL(1:nX(2)), thetaR(1:nX(2)), sintheta
      INTEGER  :: n, iX1, iX2, iX3

      ! --- Merging is assumed to happen in the X2 direction
      !     such that the total number of cells in the X2 direction
      !     is divisible by 2. Equidistant meshes are also assumed ---

      ASSOCIATE( dX1 => MeshX(1) % Width, &
                 dX2 => MeshX(2) % Width, &
                 dX3 => MeshX(3) % Width )

      r      = MeshX(1) % Center(1:nX(1)) + 0.5_DP * dX1(1:nX(1)) ! use uGF for r?
      thetaL = MeshX(2) % Center(1:nX(2)) - 0.5_DP * dX2(1:nX(2)) ! use uGF for r*sin(theta)?
      thetaR = MeshX(2) % Center(1:nX(2)) + 0.5_DP * dX2(1:nX(2)) ! use uGF for r*sin(theta)?

      DO iX1 = 1, nX(1)

        ASSOCIATE ( NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                    MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker )

      DO iX2 = 1, nX(2)

        ASSOCIATE( NCellsX3          => MergedMeshX3(iX1,iX2) % NCells, &
                   NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                   dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                   MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                   FineCellMarkerX3  => MergedMeshX3(iX1,iX2) % FineCellMarker )

        n        = 0
        sintheta = MAX( SIN(thetaL(MergeCellMarkerX2(iX2))), &
                        SIN(thetaR(MergeCellMarkerX2(iX2)+NCellsPerMergeX2-1)) )

        ! --- Determine how many cells per merge in the X3 direction
        !     and the number of merged cells in the X3 direction ---
        IF ( r(iX1) * sintheta * dX3(1) .LE. &
             dX1(iX1) / REAL(2**(MinCells-1),DP) ) THEN

          DO WHILE ( REAL(2**n,DP) * r(iX1) * sintheta * dX3(1) .LE. &
                     dX1(iX1) / REAL(2**(MinCells-1),DP) )

            n = n + 1

          END DO

          ! NCellsX3         = MAX(nX(3) / 2**n, 2**MinCells)
          NCellsX3         = nX(3) / 2**n
          ! print *, nX(2), n, 2**n, NCellsX3
          NCellsPerMergeX3 = nX(3) / NCellsX3

        ELSE

          NCellsPerMergeX3 = 1
          NCellsX3         = nX(3)

        END IF

        ! print *, 'NCellsPerMergeX3 = ', NCellsPerMergeX3, 'NCellsX3 = ', NCellsX3

        ! --- Assumes an equidsitant mesh ---
        ! --- Mark the index of the first fine cell in the merged cell ---
        DO iX3 = 1, nX(3), NCellsPerMergeX3

          dX3M(iX3:iX3+NCellsPerMergeX3-1) &
            = dX3(iX3) * REAL(NCellsPerMergeX3,DP)

          MergeCellMarkerX3(iX3:iX3+NCellsPerMergeX3-1) = iX3

        END DO

        ! --- Mark the cells which are merged according to fine cell # within coarse cell ---
        DO iX3 = 1, nX(3)

          FineCellMarkerX3(iX3) = MOD(iX3 - 1, NCellsPerMergeX3) + 1

        END DO

        END ASSOCIATE ! NCellsX3, NCellsPerMergeX3, etc.

      END DO

        END ASSOCIATE ! NCellsPerMergeX2, etc.

      END DO

      END ASSOCIATE ! dX1, dX2, dX3

      ! STOP

    END SUBROUTINE Determine_MergedCellsX3

    SUBROUTINE Determine_BasisCoeffX2 ( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER  :: iX1, iCell, iFine, iMerge, iXQ
      REAL(DP) :: xQ_F(nN), wQ_F(nN), xQ_M(nN)
      REAL(DP) :: a, b, dX2

      ! --- An equidistant mesh in the X2 direction is assumed ---

      ! --- Initialize Basis Coefficients ---
      DO iX1 = 1 ,nX(1)

        DO iCell  = 1, MergedMeshX2(iX1) % NCellsPerMerge
        DO iFine  = 1, nN
        DO iMerge = 1, nN

          MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
            = Zero

        END DO
        END DO
        END DO

      END DO

      CALL GetQuadrature( nN, xQ_F, wQ_F )

      ! --- Compute integral over the fine cell using
      ! --- the reference interval [-1/2,1/2]
      DO iX1 = 1 ,nX(1)

        ASSOCIATE( NCellsPerMergeX2 => MergedMeshX2(iX1) % NCellsPerMerge, &
                   BasisCoeffX2     => MergedMeshX2(iX1) % MergedBasisCoeff )

        DO iCell  = 1, NCellsPerMergeX2

          ! --- Get interval endpoints for fine cell w.r.t. ref. interval ---
          a = -Half + &
              (REAL(iCell,DP) - One) / REAL(NCellsPerMergeX2,DP)
          b = a + One / REAL(NCellsPerMergeX2,DP)
          dX2 = MeshX(2) % Width(1)

        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN

            ! --- Map ref. quadrature points to the fine cell ---
            xQ_M(iXQ) = (b - a) * xQ_F(iXQ) + (a + b) / Two

            BasisCoeffX2(iMerge,iFine,iCell) &
              = BasisCoeffX2(iMerge,iFine,iCell) + &
                L_X2(iMerge) % P(xQ_M(iXQ)) * L_X2(iFine) % P(xQ_F(iXQ)) * &
                wQ_F(iXQ) * dX2

                ! print *, L_X2(iMerge) % P(xQ_M(iXQ)), L_X2(iFine) % P(xQ_F(iXQ))

          END DO
          
          BasisCoeffX2(iMerge,iFine,iCell) &
              = BasisCoeffX2(iMerge,iFine,iCell) / (wQ_F(iFine) * dX2)

        END DO
        END DO
        END DO

        END ASSOCIATE ! NCellsPerMergeX2, BasisCoeffX2

      END DO
      ! STOP

    END SUBROUTINE Determine_BasisCoeffX2

    SUBROUTINE Determine_BasisCoeffX3 ( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER  :: iX1, iX2, iCell, iFine, iMerge, iXQ
      REAL(DP) :: xQ_F(nN), wQ_F(nN), xQ_M(nN)
      REAL(DP) :: a, b, dX3

      ! --- An equidistant mesh in the X3 direction is assumed ---

      ! --- Initialize Basis Coefficients ---
      DO iX1 = 1 ,nX(1)
      DO iX2 = 1 ,nX(2)

        DO iCell  = 1, MergedMeshX3(iX1,iX2) % NCellsPerMerge
        DO iFine  = 1, nN
        DO iMerge = 1, nN

          MergedMeshX3(iX1,iX2) % MergedBasisCoeff(iMerge,iFine,iCell) &
            = Zero

        END DO
        END DO
        END DO

      END DO
      END DO

      CALL GetQuadrature( nN, xQ_F, wQ_F )

      ! --- Compute integral over the fine cell using
      ! --- the reference interval [-1/2,1/2]
      DO iX1 = 1 ,nX(1)
      DO iX2 = 1 ,nX(2)

        ASSOCIATE( NCellsPerMergeX3 => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                   BasisCoeffX3     => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

        DO iCell  = 1, NCellsPerMergeX3

          ! --- Get interval endpoints for fine cell w.r.t. ref. interval ---
          a = -Half + &
              (REAL(iCell,DP) - One) / REAL(NCellsPerMergeX3,DP)
          b = a + One / REAL(NCellsPerMergeX3,DP)
          dX3 = MeshX(3) % Width(1)

        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN

            ! --- Map ref. quadrature points to the fine cell ---
            xQ_M(iXQ) = (b - a) * xQ_F(iXQ) + (a + b) / Two

            BasisCoeffX3(iMerge,iFine,iCell) &
              = BasisCoeffX3(iMerge,iFine,iCell) + &
                L_X3(iMerge) % P(xQ_M(iXQ)) * L_X3(iFine) % P(xQ_F(iXQ)) * &
                wQ_F(iXQ) * dX3

                ! print *, L_X3(iMerge) % P(xQ_M(iXQ)), L_X3(iFine) % P(xQ_F(iXQ))

          END DO
          
          BasisCoeffX3(iMerge,iFine,iCell) &
              = BasisCoeffX3(iMerge,iFine,iCell) / (wQ_F(iFine) * dX3)

        END DO
        END DO
        END DO

        END ASSOCIATE ! NCellsPerMergeX3, BasisCoeffX3

      END DO
      END DO
      ! STOP

    END SUBROUTINE Determine_BasisCoeffX3

    SUBROUTINE Determine_MergedGeometry2D ( nN ) ! Does this work for 3D? Double check the math.

      INTEGER, INTENT(in)     :: nN

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iX1, iX2, iX3
      INTEGER  :: iGCellX2, iGFineX2
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: geom_sum

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width )

      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

      DO iNodeX = 1,nDOFX

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        ! iNodeX3   = NodeNumberTableX(3,iNodeX)

        geom_sum = Zero
        DO iGCellX2 = 1,NCellsPerMergeX2
        DO iGFineX2 = 1,nN

          geom_sum = &
            geom_sum + &
            BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
            wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
            uGF(MOD(iGFineX2-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
                iX3,iGF_SqrtGm)

        END DO
        END DO

        SqrtGm_Merge(iNodeX,iX1,iX2,iX3) = geom_sum / (wQ(iNodeX2) * dX2M(iX2))

      END DO

        END ASSOCIATE ! NCellsPerMergeX2, dX2M, etc.

      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

    END SUBROUTINE Determine_MergedGeometry2D

    SUBROUTINE Determine_MergedGeometry3D ( nN )

      INTEGER, INTENT(in)     :: nN

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iX1, iX2, iX3
      INTEGER  :: iGCellX2, iGCellX3, iGFineX2, iGFineX3
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: geom_sum

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width, &
                 dX3 => MeshX(3) % Width )

      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1    ) % NCellsPerMerge, &
                   NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1    ) % MergeWidth, &
                   dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1    ) % MergeCellMarker, &
                   MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1    ) % MergedBasisCoeff, &
                   BasisCoeffX3      => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

      DO iNodeX = 1,nDOFX

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)

        geom_sum = Zero

        DO iGCellX2 = 1,NCellsPerMergeX2
        DO iGCellX3 = 1,NCellsPerMergeX3
        DO iGFineX2 = 1,nN
        DO iGFineX3 = 1,nN

          geom_sum = &
            geom_sum + &
            BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
            BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
            wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
            wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
            uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
                iX1, &
                MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
                MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
                iGF_SqrtGm)

        END DO
        END DO
        END DO
        END DO

        SqrtGm_Merge(iNodeX,iX1,iX2,iX3) = geom_sum / &
                                          (wQ(iNodeX2) * dX2M(iX2) * &
                                           wQ(iNodeX3) * dX3M(iX3) )

      END DO

        END ASSOCIATE ! NCellsPerMergeX2, etc.

      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2, dX3

    END SUBROUTINE Determine_MergedGeometry3D

    SUBROUTINE MergeAndRestrict2D ( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iX1, iX2, iX3, iCF
      INTEGER  :: iCell, iFine, iMerge
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: coeff_sum, geom_sum
      REAL(DP) :: uCF(1:nDOFX, &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(1):iX_E0(1), &
                      iX_B0(3):iX_E0(3), &
                      1:nCF)

      CALL TimersStart_Euler( Timer_CellMerging )

      ! --- Permute conserved quantities ---
      
#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)

        DO iNodeX = 1,nDOFX

          uCF(iNodeX,iX2,iX1,iX3,iCF) = U(iNodeX,iX1,iX2,iX3,iCF) * &
                                        uGF(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

        END DO

      END DO
      END DO
      END DO
      END DO

      ! --- Merge and Restrict conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width )

      CALL TimersStart_Euler( Timer_CM_UpdateCoefficient )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeX1, iNodeX2, iNodeX3 ) &
    !$OMP REDUCTION( +: coeff_sum )
#endif      
      DO iCF    = 1, nCF
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iNodeX = 1,nDOFX

        IF ( MergedMeshX2(iX1) % NCellsPerMerge .EQ. 1 ) CYCLE

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)
        
        DO iCell  = 1,NCellsPerMergeX2
        DO iFine  = 1,nN
        DO iMerge = 1,nN

          coeff_sum = &
            coeff_sum + &
            BasisCoeffX2(iMerge,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            BasisCoeffX2(iMerge, &
                         iNodeX2, &
                         FineCellMarkerX2(iX2)) * &
            uCF(MOD(iFine-1,nN)*nN+iNodeX1, & ! Can use NodeNumberTableX3D(iNodeX1,iFine,iNodeX3) instead?
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX1,iX3,iCF) / &
            (wQ(iMerge ) * dX2M(iX2) * &
             SqrtGm_Merge(MOD(iMerge-1,nN)*nN+iNodeX1,iX1,iX2,iX3))

        END DO
        END DO
        END DO

        U(iNodeX,iX1,iX2,iX3,iCF) = coeff_sum

        END ASSOCIATE ! NCellsPerMergeX2, dX2M, etc.

      END DO
      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

      CALL TimersStop_Euler( Timer_CM_UpdateCoefficient)
      CALL TimersStop_Euler( Timer_CellMerging )

    END SUBROUTINE MergeAndRestrict2D

    SUBROUTINE MergeAndRestrict3D ( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iX1, iX2, iX3, iCF
      INTEGER  :: iCellX2, iFineX2, iMergeX2
      INTEGER  :: iCellX3, iFineX3, iMergeX3
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: coeff_sum, geom_sum
      REAL(DP) :: uCF(1:nDOFX, &
                      iX_B0(3):iX_E0(3), &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(1):iX_E0(1), &
                      1:nCF)

      CALL TimersStart_Euler( Timer_CellMerging )

      ! --- Permute conserved quantities ---

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX3 = iX_B0(3), iX_E0(3)

        DO iNodeX = 1,nDOFX

          uCF(iNodeX,iX3,iX2,iX1,iCF) = U(iNodeX,iX1,iX2,iX3,iCF) * &
                                        uGF(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

        END DO

      END DO
      END DO
      END DO
      END DO

      ! --- Merge and Restrict conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width, &
                 dX3 => MeshX(3) % Width )

      CALL TimersStart_Euler( Timer_CM_UpdateCoefficient )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeX1, iNodeX2, iNodeX3 ) &
    !$OMP REDUCTION( +: coeff_sum )
#endif
      DO iCF    = 1, nCF
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iNodeX = 1,nDOFX

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff, &
                   NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                   dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                   MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                   FineCellMarkerX3  => MergedMeshX3(iX1,iX2) % FineCellMarker, &
                   BasisCoeffX3      => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

        IF ( (NCellsPerMergeX2 .EQ. 1) .AND. (NCellsPerMergeX3 .EQ. 1) ) CYCLE

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)

        DO iCellX2  = 1,NCellsPerMergeX2
        DO iCellX3  = 1,NCellsPerMergeX3
        DO iFineX2  = 1,nN
        DO iFineX3  = 1,nN
        DO iMergeX2 = 1,nN
        DO iMergeX3 = 1,nN

          coeff_sum = &
            coeff_sum + &
            BasisCoeffX2(iMergeX2,iFineX2,iCellX2) * &
            BasisCoeffX3(iMergeX3,iFineX3,iCellX3) * &
            wQ(iFineX2) * dX2(MergeCellMarkerX2(iX2) + iCellX2 - 1) * &
            wQ(iFineX3) * dX3(MergeCellMarkerX3(iX3) + iCellX3 - 1) * &
            BasisCoeffX2(iMergeX2, &
                         iNodeX2, &
                         FineCellMarkerX2(iX2)) * &
            BasisCoeffX3(iMergeX3, &
                         iNodeX3, &
                         FineCellMarkerX3(iX3)) * &
            uCF(NodeNumberTableX3D(iNodeX1,iFineX2,iFineX3), &
                MergeCellMarkerX3(iX3) + iCellX3 - 1, &
                MergeCellMarkerX2(iX2) + iCellX2 - 1, &
                iX1, &
                iCF) / &
            (wQ(iMergeX2) * dX2M(iX2) * wQ(iMergeX3) * dX3M(iX3) * &
             SqrtGm_Merge(NodeNumberTableX3D(iNodeX1,iMergeX2,iMergeX3),iX1,iX2,iX3))

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

        U(iNodeX,iX1,iX2,iX3,iCF) = coeff_sum

        END ASSOCIATE ! NCellsPerMergeX2, etc.

      END DO
      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2, dX3

      CALL TimersStop_Euler( Timer_CM_UpdateCoefficient)
      CALL TimersStop_Euler( Timer_CellMerging )

    END SUBROUTINE MergeAndRestrict3D

    SUBROUTINE MergeAndRestrict ( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      IF( nDims .LT. 3 )THEN

        CALL MergeAndRestrict2D( nN, U )

      ELSE

        CALL MergeAndRestrict3D( nN, U )

      END IF

    END SUBROUTINE MergeAndRestrict

    SUBROUTINE MergeAndRestrictGeometry ( nN, G )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: G(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nGF)

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iX1, iX2, iX3, iGF
      INTEGER  :: iCell, iFine, iMerge
      INTEGER  :: iGCell, iGFine, iGMerge
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: coeff_sum, geom_sum, dx
      REAL(DP) :: uGF_P(1:nDOFX, &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(1):iX_E0(1), &
                      iX_B0(3):iX_E0(3), &
                      1:nGF)

      ! --- Permute conserved quantities ---

      DO iGF = 1, nGF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)

        DO iNodeX = 1,nDOFX

          ! uCF(iNodeX,iX2,iX1,iX3,iCF) = U(iNodeX,iX1,iX2,iX3,iCF)
          uGF_P(iNodeX,iX2,iX1,iX3,iGF) = G(iNodeX,iX1,iX2,iX3,iGF)
          ! uCF(iNodeX,iX2,iX1,iX3,iCF) = U(iNodeX,iX1,iX2,iX3,iCF) * uGF(iNodeX,iX1,iX2,iX3,iGF_h_2)

        END DO

      END DO
      END DO
      END DO
      END DO

      ! --- Merge and Restrict conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )

      ! dx = MeshX(2) % Width(1)
      ASSOCIATE( dX2 => MeshX(2) % Width )

      DO iGF    = 1, nGF
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iNodeX = 1,nDOFX

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        
        DO iCell  = 1,MergedMeshX2(iX1) % NCellsPerMerge
        DO iFine  = 1,nN
        DO iMerge = 1,nN

          coeff_sum = &
            coeff_sum + &
            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) * &
            wQ(iFine) * dX2(MergedMeshX2(iX1) % MergeCellMarker(iX2) + iCell - 1) * &
            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge, &
                                                 iNodeX2, &
                                                 MergedMeshX2(iX1) % FineCellMarker(iX2)) * &
            uGF_P(MOD(iFine-1,nN)*nN+iNodeX1, & ! Can use NodeNumberTableX3D(iNodeX1,iFine,iNodeX3) instead?
                  MergedMeshX2(iX1) % MergeCellMarker(iX2) + iCell - 1, &
                  iX1,iX3,iGF) / &
              (wQ(iMerge ) * MergedMeshX2(iX1) % MergeWidth(iX2))

        END DO
        END DO
        END DO

        G(iNodeX,iX1,iX2,iX3,iGF) = coeff_sum

      END DO
      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

    END SUBROUTINE MergeAndRestrictGeometry

END MODULE CellMergingModule