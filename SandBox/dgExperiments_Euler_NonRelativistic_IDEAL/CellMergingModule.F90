MODULE CellMergingModule

    USE KindModule, ONLY: &
      DP, Zero, Half, One, Two, Pi
    USE ProgramHeaderModule, ONLY: &
      iX_B0, iX_E0, iX_B1, iX_E1, nDOFX, nNodesX
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

    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: SqrtGm_Merge

    TYPE, PUBLIC :: MergedMeshType
      INTEGER                                 :: NCellsPerMerge
      INTEGER                                 :: NCells
      INTEGER,  DIMENSION(:),     ALLOCATABLE :: FineCellMarker
      INTEGER,  DIMENSION(:),     ALLOCATABLE :: MergeCellMarker
      REAL(DP), DIMENSION(:),     ALLOCATABLE :: MergeWidth
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: MergedBasisCoeff
    END type MergedMeshType
    
    TYPE(MergedMeshType), DIMENSION(:), ALLOCATABLE, PUBLIC :: MergedMeshX2 ! X2 merging

    PUBLIC :: Initialize_CellMerging
    PUBLIC :: Finalize_CellMerging
    PUBLIC :: Determine_MergedCells
    PUBLIC :: Determine_BasisCoeff
    PUBLIC :: Determine_MergedGeometry
    PUBLIC :: MergeAndRestrict
    PUBLIC :: MergeAndRestrictGeometry

CONTAINS

    SUBROUTINE Initialize_CellMerging( nX, nN, Min_NCellsPerMerge_Option )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN
      INTEGER, INTENT(in), OPTIONAL :: Min_NCellsPerMerge_Option

      INTEGER :: iX1
      INTEGER :: Min_NCellsPerMerge

      ! --- Create MergedMesh for X2 direction ---
      ALLOCATE( MergedMeshX2(1:nX(1)) )
      MergedMeshX2 % NCellsPerMerge = 1
      MergedMeshX2 % NCells = nX(2)
      DO iX1 = 1, nX(1)
        ALLOCATE( MergedMeshX2(iX1) % MergeWidth     (1:nX(2)) )
        ALLOCATE( MergedMeshX2(iX1) % FineCellMarker (1:nX(2)) )
        ALLOCATE( MergedMeshX2(iX1) % MergeCellMarker(1:nX(2)) )
      END DO

      IF( PRESENT( Min_NCellsPerMerge_Option ) )THEN
        Min_NCellsPerMerge = Min_NCellsPerMerge_Option
      ELSE
        Min_NCellsPerMerge = 1
      ENDIF

      CALL Determine_MergedCells( nX, Min_NCellsPerMerge )

      DO iX1 = 1, nX(1)
        ALLOCATE( MergedMeshX2(iX1) % &
                  MergedBasisCoeff(1:nN, &
                                   1:nN, &
                                   1:MergedMeshX2(iX1) % NCellsPerMerge))
      END DO

      CALL Determine_BasisCoeff( nX, nN )

      ALLOCATE( SqrtGm_Merge(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) )

      CALL Determine_MergedGeometry( nN )

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1

      DEALLOCATE( SqrtGm_Merge )

      DO iX1 = 1, nX(1)
        DEALLOCATE( MergedMeshX2(iX1) % MergedBasisCoeff )
        DEALLOCATE( MergedMeshX2(iX1) % MergeCellMarker  )
        DEALLOCATE( MergedMeshX2(iX1) % FineCellMarker   )
        DEALLOCATE( MergedMeshX2(iX1) % MergeWidth       )
      END DO

      DEALLOCATE( MergedMeshX2 )

    END SUBROUTINE Finalize_CellMerging

    SUBROUTINE Determine_MergedCells( nX, MinCells )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: MinCells

      REAL(DP) :: r(1:nX(1))
      INTEGER  :: n, CellNumber, iX1, iX2
      INTEGER  :: iNodeX, iNodeX1, iNodeX2
      REAL(DP) :: X1,X2

      ! --- Merging is assumed to happen in the X2 direction
      !     such that the total number of cells in the X2 direction
      !     is divisible by 2. Equidistant meshes are also assumed ---

      r = MeshX(1) % Center(1:nX(1)) + 0.5_DP * MeshX(1) % Width(1:nX(1)) ! It's possible to replace r with uGF?

      DO iX1 = 1, nX(1)

        n = 0

        ! --- Determine how many cells per merge in the X2 direction
        !     and the number of merged cells in the X2 direction ---
        IF ( r(iX1) * MeshX(2) % Width(1) .LE. MeshX(1) % Width(iX1) ) THEN

          DO WHILE (2**n * r(iX1) * MeshX(2) % Width(1) &
                    .LE. &
                    MeshX(1) % Width(iX1))

            n = n + 1

          END DO

          MergedMeshX2(iX1) % NCells         = MAX(nX(2) / 2**n, 2**MinCells)
          MergedMeshX2(iX1) % NCellsPerMerge = nX(2) / MergedMeshX2(iX1) % NCells

        ELSE

          MergedMeshX2(iX1) % NCellsPerMerge = 1
          MergedMeshX2(iX1) % NCells         = nX(2)

        END IF

        ! print *, 'NCells = ', MergedMeshX2(iX1) % NCells, 'NCellsPerMerge = ', MergedMeshX2(iX1) % NCellsPerMerge

        ! --- Replace fine cell width with merged cell width ---
        ! --- Assumes an equidsitant mesh ---
        ! --- Mark the index of the first fine cell in the merged cell ---
        DO iX2 = 1, nX(2), MergedMeshX2(iX1) % NCellsPerMerge

          MergedMeshX2(iX1) % MergeWidth(iX2:iX2+&
                                          MergedMeshX2(iX1)%NCellsPerMerge-1)&
            = MeshX(2) % Width(iX2) * REAL(MergedMeshX2(iX1)%NCellsPerMerge,DP)

          MergedMeshX2(iX1) % MergeCellMarker(iX2:iX2+MergedMeshX2(iX1) % NCellsPerMerge - 1) &
            = iX2

        END DO

        ! --- Mark the cells which are merged according to fine cell # within coarse cell ---
        DO iX2 = 1, nX(2)

          MergedMeshX2(iX1) % FineCellMarker(iX2) &
            = MOD(iX2 - 1, MergedMeshX2(iX1) % NCellsPerMerge) + 1

        END DO

      END DO
      ! STOP

    END SUBROUTINE Determine_MergedCells

    SUBROUTINE Determine_BasisCoeff ( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER  :: iX1, iCell, iFine, iMerge, iXQ
      REAL(DP) :: xQ_F(nN), wQ_F(nN), xQ_M(nN)
      REAL(DP) :: a, b, dx

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

        DO iCell  = 1, MergedMeshX2(iX1) % NCellsPerMerge

          ! --- Get interval endpoints for fine cell w.r.t. ref. interval ---
          a = -Half + &
              (REAL(iCell,DP) - One) / &
              REAL(MergedMeshX2(iX1) % NCellsPerMerge,DP)
          b = a + One / REAL(MergedMeshX2(iX1) % NCellsPerMerge,DP)
          dx = MeshX(2) % Width(1)

        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN ! There is a scaling issue

            ! --- Map ref. quadrature points to the fine cell ---
            xQ_M(iXQ) = (b - a) * xQ_F(iXQ) + (a + b) / Two

            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
              = MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) + &
                L_X2(iMerge) % P(xQ_M(iXQ)) * L_X2(iFine) % P(xQ_F(iXQ)) * &
                wQ_F(iXQ) * dx

                ! print *, L_X2(iMerge) % P(xQ_M(iXQ)), L_X2(iFine) % P(xQ_F(iXQ))

          END DO
          
          MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
              = MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) / &
              (wQ_F(iFine) * dx)

        END DO
        END DO
        END DO

      END DO
      ! STOP

    END SUBROUTINE Determine_BasisCoeff

    SUBROUTINE Determine_MergedGeometry ( nN )

      INTEGER, INTENT(in)     :: nN

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iX1, iX2, iX3
      INTEGER  :: iGCell, iGFine
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: geom_sum

      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX1 => MeshX(1) % Width, &
                 dX2 => MeshX(2) % Width )

      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iNodeX = 1,nDOFX

        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        ! iNodeX3   = NodeNumberTableX(3,iNodeX)

        geom_sum = Zero
        DO iGCell = 1,MergedMeshX2(iX1) % NCellsPerMerge
        DO iGFine = 1,nN
        geom_sum = &
          geom_sum + &
          MergedMeshX2(iX1) % MergedBasisCoeff(iNodeX2,iGFine,iGCell) * &
          wQ(iGFine) * dX2(MergedMeshX2(iX1) % MergeCellMarker(iX2) + iGCell - 1) * &
          uGF(MOD(iGFine-1,nN)*nN+iNodeX1, &
              iX1, &
              MergedMeshX2(iX1) % MergeCellMarker(iX2) + iGCell - 1, &
              iX3,iGF_SqrtGm)
        END DO
        END DO

        SqrtGm_Merge(iNodeX,iX1,iX2,iX3) = geom_sum / (wQ(iNodeX2) * MergedMeshX2(iX1) % MergeWidth(iX2))

      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

    END SUBROUTINE Determine_MergedGeometry

    SUBROUTINE MergeAndRestrict ( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iX1, iX2, iX3, iCF
      INTEGER  :: iCell, iFine, iMerge
      INTEGER  :: iGCell, iGFine, iGMerge
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: coeff_sum, geom_sum
      REAL(DP) :: uCF(1:nDOFX, &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(1):iX_E0(1), &
                      iX_B0(3):iX_E0(3), &
                      1:nCF)

      ! --- Permute conserved quantities ---

      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)

        DO iNodeX = 1,nDOFX

          uCF(iNodeX,iX2,iX1,iX3,iCF) = U(iNodeX,iX1,iX2,iX3,iCF) * uGF(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

        END DO

      END DO
      END DO
      END DO
      END DO

      ! --- Merge and Restrict conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )

      ASSOCIATE( dX1 => MeshX(1) % Width, &
                 dX2 => MeshX(2) % Width )

      DO iCF    = 1, nCF
      DO iX3    = iX_B0(3), iX_E0(3)
      DO iX1    = iX_B0(1), iX_E0(1)
      DO iX2    = iX_B0(2), iX_E0(2)
      DO iNodeX = 1,nDOFX

        coeff_sum = Zero
        iNodeX1   = NodeNumberTableX(1,iNodeX)
        iNodeX2   = NodeNumberTableX(2,iNodeX)
        iNodeX3   = NodeNumberTableX(3,iNodeX)
        
        DO iCell  = 1,MergedMeshX2(iX1) % NCellsPerMerge
        DO iFine  = 1,nN
        DO iMerge = 1,nN

          ! geom_sum = Zero
          ! DO iGCell = 1,MergedMeshX2(iX1) % NCellsPerMerge
          ! DO iGFine = 1,nN
          ! geom_sum = &
          !   geom_sum + &
          !   MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iGFine,iGCell) * &
          !   wQ(iGFine) * dX2(MergedMeshX2(iX1) % MergeCellMarker(iX2) + iGCell - 1) * &
          !   uGF(MOD(iGFine-1,nN)*nN+iNodeX1, & ! Can use NodeNumberTableX3D(iNodeX1,iFine,iNodeX3) instead?
          !       iX1, &
          !       MergedMeshX2(iX1) % MergeCellMarker(iX2) + iGCell - 1, &
          !       iX3,iGF_SqrtGm) !/ &
          ! END DO
          ! END DO

          ! geom_sum = geom_sum / (wQ(iMerge ) * MergedMeshX2(iX1) % MergeWidth(iX2))

          coeff_sum = &
            coeff_sum + &
            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) * &
            wQ(iFine) * dX2(MergedMeshX2(iX1) % MergeCellMarker(iX2) + iCell - 1) * &
            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge, &
                                                 iNodeX2, &
                                                 MergedMeshX2(iX1) % FineCellMarker(iX2)) * &
            uCF(MOD(iFine-1,nN)*nN+iNodeX1, & ! Can use NodeNumberTableX3D(iNodeX1,iFine,iNodeX3) instead?
                MergedMeshX2(iX1) % MergeCellMarker(iX2) + iCell - 1, &
                iX1,iX3,iCF) / &
            (wQ(iMerge ) * MergedMeshX2(iX1) % MergeWidth(iX2) * & !* geom_sum)
             SqrtGm_Merge(MOD(iMerge-1,nN)*nN+iNodeX1,iX1,iX2,iX3))

            ! print *, 'Current_method = ', MOD(iFine-1,nN)*nN+iNodeX1
            ! print *, 'New_method     = ', NodeNumberTableX3D(iNodeX1,iFine,iNodeX3)

        END DO
        END DO
        END DO

        U(iNodeX,iX1,iX2,iX3,iCF) = coeff_sum

      END DO
      END DO
      END DO
      END DO
      END DO

      END ASSOCIATE ! dX2

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