MODULE CellMergingModule

    USE KindModule, ONLY: &
      DP, Zero, Half, One, Two, Pi
    USE ProgramHeaderModule, ONLY: & ! Work on replacing input arguments with public variables from ProgramHeader
      iX_B0, iX_E0, iX_B1, iX_E1, nDOFX, nNodesX, nDimsX, nDims
    USE GeometryFieldsModule, ONLY: &
      uGF, iGF_h_1, iGF_h_2, iGF_h_3, iGF_SqrtGm, nGF, & ! In spherical coordinates: iGF_h_2 => r, iGF_h_3 => r*sin\theta
      iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
    USE FluidFieldsModule, ONLY: &
      nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, &
      iCF_E, iCF_Ne
    USE QuadratureModule, ONLY: &
      GetQuadrature
    USE ReferenceElementModuleX, ONLY: &
      NodesX1, NodesX2, NodesX3, NodesX_q, &
      NodeNumberTableX, NodeNumberTableX3D
    USE MeshModule, ONLY: &
      MeshX, &
      NodeCoordinate
    USE PolynomialBasisModuleX_Lagrange, ONLY: &
      L_X1, L_X2, L_X3
    USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL, ONLY: &
      ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL
    USE TimersModule_Euler, ONLY: &
      TimersStart_Euler, TimersStop_Euler, &
      Timer_CellMerging, Timer_CM_UpdateCoefficient

    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Gm_Merge
    ! REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: U_Merge

    TYPE, PUBLIC :: MergedMeshType
      INTEGER                                 :: NCellsPerMerge
      INTEGER                                 :: NCells
      INTEGER                                 :: nPT
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
    PUBLIC :: ComputeMergeGeometryValues2D
    PUBLIC :: ComputeMergeGeometryValues3D
    PUBLIC :: ComputeMergeFluidValues2D
    PUBLIC :: ComputeMergeFluidValues3D
    PUBLIC :: ProlongFluidValues2D
    PUBLIC :: ProlongFluidValues3D
    PUBLIC :: MergeAndRestrict2D
    PUBLIC :: MergeAndRestrict3D
    PUBLIC :: MergeAndRestrict
    PUBLIC :: MergeAndProlong2D ! These functions will perform a merge step and then a prolong step
    PUBLIC :: MergeAndProlong3D ! which is in opposition to MergeAndRestric, which does both steps in a single step
    PUBLIC :: MergeAndProlong   ! This will allow for the positivity limiter to be applied between merging and prolonging
    PUBLIC :: GetPLRefPoints_MergedCell
    PUBLIC :: GetPLInterpMat

    INTEGER :: nPPX2(5), nPPX3(7)
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: InterpMat_XRef
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: XRef_p

CONTAINS

    SUBROUTINE Initialize_CellMerging ( nX, nN, Min_NCellsPerMerge_Option )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN
      INTEGER, INTENT(in), OPTIONAL :: Min_NCellsPerMerge_Option

      INTEGER :: iX1, iX2, i, Min_NCellsPerMerge

      IF( PRESENT( Min_NCellsPerMerge_Option ) )THEN

        Min_NCellsPerMerge = Min_NCellsPerMerge_Option

      ELSE

        Min_NCellsPerMerge = 1

      ENDIF

      ! --- Create MergedMesh for X2 direction ---
      
      ALLOCATE( MergedMeshX2(1:nX(1)) )

      MergedMeshX2 % NCellsPerMerge = 1
      MergedMeshX2 % NCells         = nX(2)

      nPPX2    = 0
      nPPX2(1) = PRODUCT(nNodesX(1:2))
      DO i = 1, 2
        IF( nNodesX(i) > 1 )THEN
          nPPX2(2*i:2*i+1) = PRODUCT( nNodesX(1:2), MASK = [1,2] .NE. i )
        END IF
      END DO
      MergedMeshX2 % nPT = SUM(nPPX2)

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

      nPPX3    = 0
      nPPX3(1) = PRODUCT(nNodesX)
      DO i = 1, 3
        IF( nNodesX(i) > 1 )THEN
          nPPX3(2*i:2*i+1) = PRODUCT( nNodesX, MASK = [1,2,3] .NE. i )
        END IF
      END DO
      MergedMeshX3 % nPT = SUM(nPPX3)

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

      ALLOCATE( Gm_Merge(1:nDOFX, &
                         iX_B1(1):iX_E1(1), &
                         iX_B1(2):iX_E1(2), &
                         iX_B1(3):iX_E1(3), &
                         1:nGF) )

      IF( nDims .LT. 3 )THEN

        CALL Determine_MergedGeometry2D( nN )

      ELSE

        CALL Determine_MergedGeometry3D( nN )

      END IF

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging ( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1, iX2

      DEALLOCATE( Gm_Merge )

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
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   nPTX2             => MergedMeshX2(iX1) % nPT )

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
        
        ! --- Get reference coordinates of points used in Positivity Limiter ---
        IF( NCellsPerMergeX2 .NE. 1 )THEN
          nPTX2 = nPTX2 & ! Merged Cell interior and boundary points
                + nPPX2(1) * NCellsPerMergeX2 & ! all Fine cell interior points
                + SUM(nPPX2(2:3)) * NCellsPerMergeX2 & ! all Fine Cell X1 boundary points
                + nPPX2(5) + nPPX2(4) * NCellsPerMergeX2 ! all Fine Cell X2 boundary points
        END IF

        END ASSOCIATE ! NCellsX2, NCellsPerMergeX2, etc.

      END DO

      END ASSOCIATE ! dX1, dX2

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
                   FineCellMarkerX3  => MergedMeshX3(iX1,iX2) % FineCellMarker, &
                   nPTX3             => MergedMeshX3(iX1,iX2) % nPT )

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

        ! --- Get reference coordinates of points used in Positivity Limiter ---
        IF( MergedMeshX3(iX1,iX2) % NCellsPerMerge .NE. 1 )THEN
          nPTX3 = nPTX3 & ! Merged Cell interior and boundary points
                + nPPX3(1) * NCellsPerMergeX3 &
                           * NCellsPerMergeX2 & ! all Fine cell interior points
                + SUM(nPPX3(2:3)) * NCellsPerMergeX3 &
                                  * NCellsPerMergeX2 & ! all Fine Cell X1 boundary points
                + nPPX3(4) * NCellsPerMergeX3 &
                           * NCellsPerMergeX2 &
                + nPPX3(5) * NCellsPerMergeX3 &! all Fine Cell X2 boundary points
                + nPPX3(6) * NCellsPerMergeX3 &
                           * NCellsPerMergeX2 & 
                + nPPX3(7) * NCellsPerMergeX2 ! all Fine Cell X3 boundary points
        END IF

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
      REAL(DP) :: geom_sum, Gm_h_Merge(iGF_h_1:iGF_h_3)

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iNodeX1, iNodeX2 ) &
    !$OMP REDUCTION( +: geom_sum )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1,nDOFX

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)

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

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) = geom_sum / (wQ(iNodeX2) * dX2M(iX2))

        END ASSOCIATE ! NCellsPerMergeX2, dX2M, etc.

      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Gm_h_Merge )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1,nDOFX

        CALL ComputeMergeGeometryValues2D( nN, iNodeX, iX1, iX2, iX3, Gm_h_Merge)

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_1) = Gm_h_Merge(iGF_h_1)
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_2) = Gm_h_Merge(iGF_h_2)
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_3) = Gm_h_Merge(iGF_h_3)

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) = Gm_h_Merge(iGF_h_1)**2
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) = Gm_h_Merge(iGF_h_2)**2
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) = Gm_h_Merge(iGF_h_3)**2

      END DO
      END DO
      END DO
      END DO

    END SUBROUTINE Determine_MergedGeometry2D

    SUBROUTINE Determine_MergedGeometry3D ( nN )

      INTEGER, INTENT(in)     :: nN

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iX1, iX2, iX3
      INTEGER  :: iGCellX2, iGCellX3, iGFineX2, iGFineX3
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: geom_sum, Gm_h_Merge(iGF_h_1:iGF_h_3)

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX2 => MeshX(2) % Width, &
                 dX3 => MeshX(3) % Width )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iNodeX1, iNodeX2, iNodeX3 ) &
    !$OMP REDUCTION( +: geom_sum )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1,nDOFX

        ASSOCIATE( NCellsPerMergeX2  => MergedMeshX2(iX1    ) % NCellsPerMerge, &
                   NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1    ) % MergeWidth, &
                   dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1    ) % MergeCellMarker, &
                   MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1    ) % MergedBasisCoeff, &
                   BasisCoeffX3      => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

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

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) = geom_sum / &
                                          (wQ(iNodeX2) * dX2M(iX2) * &
                                           wQ(iNodeX3) * dX3M(iX3) )

        END ASSOCIATE ! NCellsPerMergeX2, etc.

      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2, dX3

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( Gm_h_Merge )
#endif
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1,nDOFX

        CALL ComputeMergeGeometryValues3D( nN, iNodeX, iX1, iX2, iX3, Gm_h_Merge)

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_1) = Gm_h_Merge(iGF_h_1)
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_2) = Gm_h_Merge(iGF_h_2)
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_h_3) = Gm_h_Merge(iGF_h_3)

        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) = Gm_h_Merge(iGF_h_1)**2
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) = Gm_h_Merge(iGF_h_2)**2
        Gm_Merge(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) = Gm_h_Merge(iGF_h_3)**2

      END DO
      END DO
      END DO
      END DO

    END SUBROUTINE Determine_MergedGeometry3D

    SUBROUTINE ComputeMergeGeometryValues2D &
      ( nN, iNodeX, iX1, iX2, iX3, Gm_h_Merge )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(out)  :: Gm_h_Merge(iGF_h_1:iGF_h_3)

      INTEGER  :: iNodeX1, iNodeX2
      INTEGER  :: iCell, iFine
      REAL(DP) :: xQ(nN), wQ(nN)

      ! --- Merge conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )      

        ASSOCIATE( dX2               => MeshX(2) % Width, &
                   NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)

        Gm_h_Merge = Zero

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP REDUCTION( +: Gm_h_Merge )
#endif
        DO iCell  = 1,NCellsPerMergeX2
        DO iFine  = 1,nN

          Gm_h_Merge(iGF_h_1) = &
            Gm_h_Merge(iGF_h_1) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_h_1) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          Gm_h_Merge(iGF_h_2) = &
            Gm_h_Merge(iGF_h_2) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_h_2) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          Gm_h_Merge(iGF_h_3) = &
            Gm_h_Merge(iGF_h_3) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_h_3) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

        END DO
        END DO

        END ASSOCIATE ! dX2, NCellsPerMergeX2, dX2M, etc.

    END SUBROUTINE ComputeMergeGeometryValues2D

    SUBROUTINE ComputeMergeGeometryValues3D &
      ( nN, iNodeX, iX1, iX2, iX3, Gm_h_Merge )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(out)  :: Gm_h_Merge(iGF_h_1:iGF_h_3)

      INTEGER  :: iNodeX1, iNodeX2, iNodeX3
      INTEGER  :: iGCellX2, iGCellX3, iGFineX2, iGFineX3
      REAL(DP) :: xQ(nN), wQ(nN)

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE(  dX2               => MeshX(2) % Width, &
                  dX3               => MeshX(3) % Width, &
                  NCellsPerMergeX2  => MergedMeshX2(iX1    ) % NCellsPerMerge, &
                  NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                  dX2M              => MergedMeshX2(iX1    ) % MergeWidth, &
                  dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                  MergeCellMarkerX2 => MergedMeshX2(iX1    ) % MergeCellMarker, &
                  MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                  BasisCoeffX2      => MergedMeshX2(iX1    ) % MergedBasisCoeff, &
                  BasisCoeffX3      => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

      iNodeX1   = NodeNumberTableX(1,iNodeX)
      iNodeX2   = NodeNumberTableX(2,iNodeX)
      iNodeX3   = NodeNumberTableX(3,iNodeX)

      Gm_h_Merge = Zero

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP REDUCTION( +: Gm_h_Merge )
#endif
      DO iGCellX2 = 1,NCellsPerMergeX2
      DO iGCellX3 = 1,NCellsPerMergeX3
      DO iGFineX2 = 1,nN
      DO iGFineX3 = 1,nN

        Gm_h_Merge(iGF_h_1) = &
          Gm_h_Merge(iGF_h_1) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_h_1) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        Gm_h_Merge(iGF_h_2) = &
          Gm_h_Merge(iGF_h_2) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_h_2) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        Gm_h_Merge(iGF_h_3) = &
          Gm_h_Merge(iGF_h_3) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_h_3) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2, dX3

    END SUBROUTINE ComputeMergeGeometryValues3D

    SUBROUTINE ComputeMergeFluidValues2D &
      ( nN, iNodeX, iX1, iX2, iX3, U, uCF )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(in )  :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)
      REAL(DP), INTENT(out)  :: uCF(iCF_D:iCF_Ne)

      INTEGER  :: iNodeX1, iNodeX2
      INTEGER  :: iCell, iFine
      REAL(DP) :: xQ(nN), wQ(nN)

      ! --- Merge conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )      

        ASSOCIATE( dX2               => MeshX(2) % Width, &
                   NCellsPerMergeX2  => MergedMeshX2(iX1) % NCellsPerMerge, &
                   dX2M              => MergedMeshX2(iX1) % MergeWidth, &
                   MergeCellMarkerX2 => MergedMeshX2(iX1) % MergeCellMarker, &
                   FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)

        uCF = Zero

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP REDUCTION( +: uCF )
#endif
        DO iCell  = 1,NCellsPerMergeX2
        DO iFine  = 1,nN

          uCF(iCF_D) = &
            uCF(iCF_D) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_D) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          uCF(iCF_S1) = &
            uCF(iCF_S1) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_S1) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          uCF(iCF_S2) = &
            uCF(iCF_S2) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_S2) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          uCF(iCF_S3) = &
            uCF(iCF_S3) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_S3) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          uCF(iCF_E) = &
            uCF(iCF_E) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_E) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

          uCF(iCF_Ne) = &
            uCF(iCF_Ne) + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            U(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iCF_Ne) * &
            uGF(MOD(iFine-1,nN)*nN+iNodeX1, &
                iX1, &
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX3,iGF_SqrtGm) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

        END DO
        END DO

        END ASSOCIATE ! dX2, NCellsPerMergeX2, dX2M, etc.

    END SUBROUTINE ComputeMergeFluidValues2D

    SUBROUTINE ComputeMergeFluidValues3D &
      ( nN, iNodeX, iX1, iX2, iX3, U, uCF )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(in )  :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)
      REAL(DP), INTENT(out)  :: uCF(iCF_D:iCF_Ne)

      INTEGER  :: iNodeX1, iNodeX2, iNodeX3
      INTEGER  :: iGCellX2, iGCellX3, iGFineX2, iGFineX3
      REAL(DP) :: xQ(nN), wQ(nN)

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE(  dX2               => MeshX(2) % Width, &
                  dX3               => MeshX(3) % Width, &
                  NCellsPerMergeX2  => MergedMeshX2(iX1    ) % NCellsPerMerge, &
                  NCellsPerMergeX3  => MergedMeshX3(iX1,iX2) % NCellsPerMerge, &
                  dX2M              => MergedMeshX2(iX1    ) % MergeWidth, &
                  dX3M              => MergedMeshX3(iX1,iX2) % MergeWidth, &
                  MergeCellMarkerX2 => MergedMeshX2(iX1    ) % MergeCellMarker, &
                  MergeCellMarkerX3 => MergedMeshX3(iX1,iX2) % MergeCellMarker, &
                  BasisCoeffX2      => MergedMeshX2(iX1    ) % MergedBasisCoeff, &
                  BasisCoeffX3      => MergedMeshX3(iX1,iX2) % MergedBasisCoeff )

      iNodeX1   = NodeNumberTableX(1,iNodeX)
      iNodeX2   = NodeNumberTableX(2,iNodeX)
      iNodeX3   = NodeNumberTableX(3,iNodeX)

      uCF = Zero

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP REDUCTION( +: uCF )
#endif
      DO iGCellX2 = 1,NCellsPerMergeX2
      DO iGCellX3 = 1,NCellsPerMergeX3
      DO iGFineX2 = 1,nN
      DO iGFineX3 = 1,nN

        uCF(iCF_D) = &
          uCF(iCF_D) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_D) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        uCF(iCF_S1) = &
          uCF(iCF_S1) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_S1) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        uCF(iCF_S2) = &
          uCF(iCF_S2) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_S2) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        uCF(iCF_S3) = &
          uCF(iCF_S3) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_S3) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        uCF(iCF_E) = &
          uCF(iCF_E) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_E) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

        uCF(iCF_Ne) = &
          uCF(iCF_Ne) + &
          BasisCoeffX2(iNodeX2,iGFineX2,iGCellX2) * &
          BasisCoeffX3(iNodeX3,iGFineX3,iGCellX3) * &
          wQ(iGFineX2) * dX2(MergeCellMarkerX2(iX2) + iGCellX2 - 1) * &
          wQ(iGFineX3) * dX3(MergeCellMarkerX3(iX3) + iGCellX3 - 1) * &
          U(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iCF_Ne) * &
          uGF(NodeNumberTableX3D(iNodeX1,iGFineX2,iGFineX3), &
              iX1, &
              MergeCellMarkerX2(iX2) + iGCellX2 - 1, &
              MergeCellMarkerX3(iX3) + iGCellX3 - 1, &
              iGF_SqrtGm) / &
          (wQ(iNodeX2) * dX2M(iX2) * wQ(iNodeX3) * dX3M(iX3) * &
           Gm_Merge(NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3),iX1,iX2,iX3,iGF_SqrtGm))

      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2, dX3

    END SUBROUTINE ComputeMergeFluidValues3D

    SUBROUTINE ProlongFluidValues2D( nN, iNodeX, iX1, iX2, iX3, uCF_Merge, uCF )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(in )  :: uCF_Merge(1:nDOFX, &
                                          iCF_D:iCF_Ne)
      REAL(DP), INTENT(out)  :: uCF(iCF_D:iCF_Ne)

      INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iCF, iMerge
      REAL(DP) :: xQ(nN), wQ(nN)

      ! This subroutine should use the merged values associated with the cell
      ! located at iX1, iX2, iX3 to compute the fluid value at iNodeX on that cell
      ! It will be called by MergeAndProlong2D

      CALL GetQuadrature( nN, xQ, wQ )

      uCF = Zero

      iNodeX1   = NodeNumberTableX(1,iNodeX)
      iNodeX2   = NodeNumberTableX(2,iNodeX)
      iNodeX3   = NodeNumberTableX(3,iNodeX)

      ASSOCIATE(FineCellMarkerX2 => MergedMeshX2(iX1) % FineCellMarker, &
                BasisCoeffX2     => MergedMeshX2(iX1) % MergedBasisCoeff, &
                dX2M             => MergedMeshX2(iX1) % MergeWidth )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP REDUCTION( +: uCF )
#endif
      DO iCF = iCF_D, iCF_Ne
      DO iMerge = 1,nN

        ! uCF(iCF) = uCF(iCF) + uCF_Merge(NodeNumberTableX3D(iNodeX1,iMerge,iNodeX3),iCF) * &
        !            BasisCoeffX2(iMerge, &
        !                         iNodeX2, &
        !                         FineCellMarkerX2(iX2))
        uCF(iCF) = uCF(iCF) + uCF_Merge(NodeNumberTableX3D(iNodeX1,iMerge,iNodeX3),iCF) * &
                   wQ(iMerge ) * dX2M(iX2) * &
                   Gm_Merge(MOD(iMerge-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm) * &
                   BasisCoeffX2(iMerge, &
                                iNodeX2, &
                                FineCellMarkerX2(iX2)) / &
                  (wQ(iMerge ) * dX2M(iX2) * &
                   Gm_Merge(MOD(iMerge-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

      END DO
      END DO

      END ASSOCIATE

    END SUBROUTINE ProlongFluidValues2D

    SUBROUTINE ProlongFluidValues3D( nN, iNodeX, iX1, iX2, iX3, uCF_Merge, uCF )

      INTEGER,  INTENT(in )  :: nN, iNodeX, iX1, iX2, iX3
      REAL(DP), INTENT(in )  :: uCF_Merge(1:nDOFX, &
                                          iCF_D:iCF_Ne)
      REAL(DP), INTENT(out)  :: uCF(iCF_D:iCF_Ne)

    END SUBROUTINE ProlongFluidValues3D

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
             Gm_Merge(MOD(iMerge-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

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
             Gm_Merge(NodeNumberTableX3D(iNodeX1,iMergeX2,iMergeX3),iX1,iX2,iX3,iGF_SqrtGm))

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

    SUBROUTINE MergeAndProlong2D( nN, U )

!       INTEGER, INTENT(in)     :: nN
!       REAL(DP), INTENT(inout) :: U(1:nDOFX, &
!                                    iX_B1(1):iX_E1(1), &
!                                    iX_B1(2):iX_E1(2), &
!                                    iX_B1(3):iX_E1(3), &
!                                    1:nCF)

!       INTEGER  :: iX1, iX2, iX3, iNodeX
!       REAL(DP) :: uCF_Merge(1:nDOFX, iCF_D:iCF_Ne)
!       REAL(DP) :: uCF(iCF_D:iCF_Ne)

!       ! This subroutine will call ComputeMergeFluidValues2D to compute all the
!       ! merged fluid values in a single merged cell.
!       ! Once these values are computed for a single merged cell, call ProlongFluidValues2D
!       ! to prolong the values to all fine cells within the merged cell.
!       ! Later, a call to the positivity limiter will be added between ComputeMergeFluidValues2D
!       ! and ProlongFluidValues2D

! #if defined( THORNADO_OMP )
!     !$OMP PARALLEL DO COLLAPSE(3)
! #endif
!       DO iX1    = iX_B0(1), iX_E0(1)
!       DO iX2    = iX_B0(2), iX_E0(2)
!       DO iX3    = iX_B0(3), iX_E0(3)

!         IF ( MergedMeshX2(iX1) % NCellsPerMerge .EQ. 1 ) CYCLE

!         uCF_Merge = Zero

!         DO iNodeX = 1,nDOFX

!           CALL ComputeMergeFluidValues2D &
!             ( nN, iNodeX, iX1, iX2, iX3, U, uCF_Merge(iNodeX,iCF_D:iCF_Ne) )

!         END DO

!         ! Add a call to the positivity limiter to act on uCF_Merge

!         DO iNodeX = 1,nDOFX

!           uCF = Zero

!           CALL ProlongFluidValues2D &
!             ( nN, iNodeX, iX1, iX2, iX3, uCF_Merge, uCF )

!           U(iNodeX,iX1,iX2,iX3,iCF_D:iCF_Ne) = uCF

!         END DO

!       END DO
!       END DO
!       END DO

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
      REAL(DP) :: uCF_Merge(1:nDOFX, &
                            iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3), &
                            1:nCF)

      CALL TimersStart_Euler( Timer_CellMerging )

      uCF_Merge = U

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

      ! --- Merge conserved quantities ---
      
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
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)
        
        DO iCell  = 1,NCellsPerMergeX2
        DO iFine  = 1,nN

          coeff_sum = &
            coeff_sum + &
            BasisCoeffX2(iNodeX2,iFine,iCell) * &
            wQ(iFine) * dX2(MergeCellMarkerX2(iX2) + iCell - 1) * &
            uCF(MOD(iFine-1,nN)*nN+iNodeX1, & ! Can use NodeNumberTableX3D(iNodeX1,iFine,iNodeX3) instead?
                MergeCellMarkerX2(iX2) + iCell - 1, &
                iX1,iX3,iCF) / &
            (wQ(iNodeX2 ) * dX2M(iX2) * &
             Gm_Merge(MOD(iNodeX2-1,nN)*nN+iNodeX1,iX1,iX2,iX3,iGF_SqrtGm))

        END DO
        END DO

        uCF_Merge(iNodeX,iX1,iX2,iX3,iCF) = coeff_sum

        END ASSOCIATE ! NCellsPerMergeX2, dX2M, etc.

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL TimersStop_Euler( Timer_CM_UpdateCoefficient)
      CALL TimersStop_Euler( Timer_CellMerging )      

      ! --- Apply positivity limiter to conserved merged quantities ---

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2), MergedMeshX2(iX1) % NCellsPerMerge

        IF(MergedMeshX2(iX1) % nPT .GT. SUM(nPPX2))THEN

          ALLOCATE( XRef_p(1:3,MergedMeshX2(iX1) % nPT) )
          ALLOCATE( InterpMat_XRef(MergedMeshX2(iX1) % nPT, nDOFX) )

          CALL GetPLRefPoints_MergedCell(iX1,iX2,iX3, &
                                         MergedMeshX2(iX1) % nPT,XRef_p)

          CALL GetPLInterpMat(iX1,iX2,iX3, &
                              MergedMeshX2(iX1) % nPT,XRef_p,InterpMat_XRef)
                              
          CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
            ( [iX1,iX2,iX3], [iX1,iX2,iX3], iX_B1, iX_E1, Gm_Merge, uCF_Merge, &
              nPT_Option = MergedMeshX2(iX1) % nPT, &
              XRef_p_Option = XRef_p, &
              InterpMat_XRef_Option = InterpMat_XRef )

          DEALLOCATE(InterpMat_XRef)
          DEALLOCATE(XRef_p)

        ELSE

          CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
            ( [iX1,iX2,iX3], [iX1,iX2,iX3], iX_B1, iX_E1, Gm_Merge, uCF_Merge )

        END IF

      END DO
      END DO
      END DO

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)

        uCF_Merge(:,iX1,iX2,iX3,:) &
          = uCF_Merge(:,iX1,MergedMeshX2(iX1) % MergeCellMarker(iX2),iX3,:)

      END DO
      END DO
      END DO

      ! --- Prolong the conserved merged quantities ---

      CALL TimersStart_Euler( Timer_CellMerging )
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

        ASSOCIATE( FineCellMarkerX2  => MergedMeshX2(iX1) % FineCellMarker, &
                   BasisCoeffX2      => MergedMeshX2(iX1) % MergedBasisCoeff )

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)
        
        DO iMerge = 1,nN

          coeff_sum = &
            coeff_sum + &
            BasisCoeffX2(iMerge, &
                         iNodeX2, &
                         FineCellMarkerX2(iX2)) * &
            uCF_Merge(NodeNumberTableX3D(iNodeX1,iMerge,iNodeX3), &
                      iX1,iX2,iX3,iCF)
        END DO

        U(iNodeX,iX1,iX2,iX3,iCF) = coeff_sum

        END ASSOCIATE ! FineCellMarkerX2, BasisCoeffX2

      END DO
      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

      CALL TimersStop_Euler( Timer_CM_UpdateCoefficient)
      CALL TimersStop_Euler( Timer_CellMerging )

    END SUBROUTINE MergeAndProlong2D

    SUBROUTINE MergeAndProlong3D( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

    END SUBROUTINE MergeAndProlong3D

    SUBROUTINE MergeAndProlong( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      IF( nDims .LT. 3 )THEN

        CALL MergeAndProlong2D( nN, U )

      ELSE

        CALL MergeAndProlong3D( nN, U )

      END IF

    END SUBROUTINE MergeAndProlong

    SUBROUTINE GetPLRefPoints_MergedCell( iX1, iX2, iX3, nRefPT, XRef_p)

      INTEGER,  INTENT(in)  :: iX1, iX2, iX3, nRefPT
      REAL(DP), INTENT(out) :: XRef_p(1:3,nRefPT)

      INTEGER  :: iRefPT, iNodeX, iCellX1, iCellX2, iCellX3
      INTEGER  :: iCellX2_End, iCellX3_End
      INTEGER  :: iNodeX_X1, iNodeX_X2, iNodeX_X3
      REAL(DP) :: dX2, dX3, CenterX2, CenterX3
      
      iRefPT = 1

      ! ==========================
      ! Gauss nodes of merged cell
      ! ==========================

      DO iNodeX = 1, nDOFX

        XRef_p(1:3, iRefPT) = NodesX_q(1:3,iNodeX)
        iRefPT = iRefPT + 1

      END DO

      ! =========================
      ! Gauss nodes of fine cells
      ! =========================

      iCellX2_End = MergedMeshX2(iX1) % NCellsPerMerge

      IF( nDimsX .LT. 3)THEN
        iCellX3_End = 1
      ELSE
        iCellX3_End = MergedMeshX3(iX1,iX2) % NCellsPerMerge
      END IF

      DO iCellX2 = 1, iCellX2_End
      DO iCellX3 = 1, iCellX3_End

        dX2 = One / REAL(iCellX2_End,DP)
        dX3 = One / REAL(iCellX3_End,DP)

        CenterX2 = -Half + REAL(iCellX2-1,DP) * dX2 + Half * dX2
        CenterX3 = -Half + REAL(iCellX3-1,DP) * dX3 + Half * dX3

        DO iNodeX = 1, nDOFX

          XRef_p(1, iRefPT) = NodesX_q(1,iNodeX)
          XRef_p(2, iRefPT) = NodesX_q(2,iNodeX) * dX2 + CenterX2
          XRef_p(3, iRefPT) = NodesX_q(3,iNodeX) * dX3 + CenterX3

          iRefPT = iRefPT + 1

        END DO

      END DO
      END DO

      IF( nDOFX .EQ. 1 ) RETURN

      ! =============================
      ! Boundary nodes of merged cell
      ! =============================

      ! Left/Right Boundary X1
      IF( nDimsX .GT. 0 )THEN
        DO iNodeX_X2 = 1, nNodesX(2)
        DO iNodeX_X3 = 1, nNodesX(3)

          XRef_p(1, iRefPT) = -Half
          XRef_p(2, iRefPT) = NodesX2(iNodeX_X2)
          XRef_p(3, iRefPT) = NodesX3(iNodeX_X3)

          iRefPT = iRefPT + 1

        END DO
        END DO

        DO iNodeX_X2 = 1, nNodesX(2)
        DO iNodeX_X3 = 1, nNodesX(3)

          XRef_p(1, iRefPT) = Half
          XRef_p(2, iRefPT) = NodesX2(iNodeX_X2)
          XRef_p(3, iRefPT) = NodesX3(iNodeX_X3)

          iRefPT = iRefPT + 1

        END DO
        END DO
      END IF

      ! Left/Right Boundary X2
      IF( nDimsX .GT. 1 )THEN
        DO iNodeX_X1 = 1, nNodesX(1)
        DO iNodeX_X3 = 1, nNodesX(3)

          XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
          XRef_p(2, iRefPT) = -Half
          XRef_p(3, iRefPT) = NodesX3(iNodeX_X3)

          iRefPT = iRefPT + 1

        END DO
        END DO

        DO iNodeX_X1 = 1, nNodesX(1)
        DO iNodeX_X3 = 1, nNodesX(3)

          XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
          XRef_p(2, iRefPT) = Half
          XRef_p(3, iRefPT) = NodesX3(iNodeX_X3)

          iRefPT = iRefPT + 1

        END DO
        END DO
      END IF

      ! Left/Right Boundary X3
      IF( nDimsX .GT. 2 )THEN
        DO iNodeX_X1 = 1, nNodesX(1)
        DO iNodeX_X2 = 1, nNodesX(2)

          XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
          XRef_p(2, iRefPT) = NodesX2(iNodeX_X2)
          XRef_p(3, iRefPT) = -Half

          iRefPT = iRefPT + 1

        END DO
        END DO

        DO iNodeX_X1 = 1, nNodesX(1)
        DO iNodeX_X2 = 1, nNodesX(2)

          XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
          XRef_p(2, iRefPT) = NodesX2(iNodeX_X2)
          XRef_p(3, iRefPT) = Half

          iRefPT = iRefPT + 1

        END DO
        END DO
      END IF

      ! ============================
      ! Boundary nodes of fine cells
      ! ============================

      ! Left/Right Boundary X1
      IF( nDimsX .GT. 0 )THEN
        DO iCellX2 = 1, iCellX2_End
        DO iCellX3 = 1, iCellX3_End

          dX2 = One / REAL(iCellX2_End,DP)
          dX3 = One / REAL(iCellX3_End,DP)

          CenterX2 = -Half + REAL(iCellX2-1,DP) * dX2 + Half * dX2
          CenterX3 = -Half + REAL(iCellX3-1,DP) * dX3 + Half * dX3

          DO iNodeX_X2 = 1, nNodesX(2)
          DO iNodeX_X3 = 1, nNodesX(3)

            XRef_p(1, iRefPT) = -Half
            XRef_p(2, iRefPT) = NodesX2(iNodeX_X2) * dX2 + CenterX2
            XRef_p(3, iRefPT) = NodesX3(iNodeX_X3) * dX3 + CenterX3

            iRefPT = iRefPT + 1

          END DO
          END DO

          DO iNodeX_X2 = 1, nNodesX(2)
          DO iNodeX_X3 = 1, nNodesX(3)

            XRef_p(1, iRefPT) = Half
            XRef_p(2, iRefPT) = NodesX2(iNodeX_X2) * dX2 + CenterX2
            XRef_p(3, iRefPT) = NodesX3(iNodeX_X3) * dX3 + CenterX3

            iRefPT = iRefPT + 1

          END DO
          END DO

        END DO
        END DO
      END IF

      ! Left/Right Boundary X2
      IF( nDimsX .GT. 1 )THEN
        DO iCellX2 = 1, iCellX2_End
        DO iCellX3 = 1, iCellX3_End

          dX3      = One / REAL(iCellX3_End,DP)
          CenterX3 = -Half + REAL(iCellX3-1,DP) * dX3 + Half * dX3

          DO iNodeX_X1 = 1, nNodesX(1)
          DO iNodeX_X3 = 1, nNodesX(3)

            XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
            XRef_p(2, iRefPT) = -Half
            XRef_p(3, iRefPT) = NodesX3(iNodeX_X3) * dX3 + CenterX3

            iRefPT = iRefPT + 1

          END DO
          END DO

          IF( iCellX2 .EQ. iCellX2_End )THEN
            DO iNodeX_X1 = 1, nNodesX(1)
            DO iNodeX_X3 = 1, nNodesX(3)

              XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
              XRef_p(2, iRefPT) = Half
              XRef_p(3, iRefPT) = NodesX3(iNodeX_X3) * dX3 + CenterX3

              iRefPT = iRefPT + 1

            END DO
            END DO
          END IF

        END DO
        END DO
      END IF

      ! Left/Right Boundary X3
      IF( nDimsX .GT. 2 )THEN
        DO iCellX3 = 1, iCellX3_End
        DO iCellX2 = 1, iCellX2_End

          dX2      = One / REAL(iCellX2_End,DP)
          CenterX2 = -Half + REAL(iCellX2-1,DP) * dX2 + Half * dX2

          DO iNodeX_X1 = 1, nNodesX(1)
          DO iNodeX_X2 = 1, nNodesX(2)

            XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
            XRef_p(2, iRefPT) = NodesX2(iNodeX_X2) * dX2 + CenterX2
            XRef_p(3, iRefPT) = -Half

            iRefPT = iRefPT + 1

          END DO
          END DO

          IF( iCellX3 .EQ. iCellX3_End )THEN
            DO iNodeX_X1 = 1, nNodesX(1)
            DO iNodeX_X2 = 1, nNodesX(2)

              XRef_p(1, iRefPT) = NodesX1(iNodeX_X1)
              XRef_p(2, iRefPT) = NodesX2(iNodeX_X2) * dX2 + CenterX2
              XRef_p(3, iRefPT) = Half

              iRefPT = iRefPT + 1

            END DO
            END DO
          END IF

        END DO
        END DO
      END IF

      ! print *, 'RefPt filled in', iRefPT - 1, 'Size of vector', nRefPT
  
    END SUBROUTINE GetPLRefPoints_MergedCell

    
    SUBROUTINE GetPLInterpMat( iX1, iX2, iX3, nRefPT, XRef_p, InterpMat_XRef)

      INTEGER, INTENT(in)   :: iX1, iX2, iX3, nRefPT
      REAL(DP), INTENT(in)  :: XRef_p(1:3,nRefPT)
      REAL(DP), INTENT(out) :: InterpMat_XRef(nRefPT,nDOFX)

      INTEGER  :: iRefPT, iNodeX
      INTEGER  :: iNodeX1, iNodeX2, iNodeX3

      DO iRefPT = 1, nRefPT
      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        InterpMat_XRef(iRefPT,iNodeX) = &
          L_X1(iNodeX1) % P( XRef_p(1,iRefPT) ) &
          * L_X2(iNodeX2) % P( XRef_p(2,iRefPT) ) &
          * L_X3(iNodeX3) % P( XRef_p(3,iRefPT) )


      END DO
      END DO

    END SUBROUTINE GetPLInterpMat

END MODULE CellMergingModule