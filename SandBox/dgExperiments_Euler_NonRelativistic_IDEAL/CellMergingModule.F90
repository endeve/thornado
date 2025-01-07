MODULE CellMergingModule

    USE KindModule, ONLY: &
        DP
    USE ProgramHeaderModule, ONLY: &
        iX_B0, iX_E0
    USE GeometryFieldsModule, ONLY: &
        uGF, iGF_h_2, iGF_h_3
    USE MeshModule, ONLY: &
        MeshX
    USE LagrangePolynomialsModule, ONLY: &
        LagrangeP

    TYPE, PUBLIC :: MergedMeshType
      INTEGER :: NCellsPerMerge
      INTEGER :: NCells
      INTEGER, DIMENSION(:), ALLOCATABLE :: CellMarker
    END type MergedMeshType
    
    TYPE(MergedMeshType), DIMENSION(:), ALLOCATABLE, PUBLIC :: MergedMeshX2 ! X2 merging

    PUBLIC :: Initialize_CellMerging
    PUBLIC :: Finalize_CellMerging
    PUBLIC :: Determine_MergedCells
    PUBLIC :: MergetoCoarseGrid
    PUBLIC :: RestricttoFineGrid

    ! REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: r

CONTAINS

    SUBROUTINE Initialize_CellMerging( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1

      ! --- Create MergedMesh for X2 direction ---
      ALLOCATE( MergedMeshX2(1:nX(1)) )
      MergedMeshX2 % NCellsPerMerge = 1
      MergedMeshX2 % NCells = nX(2)
      DO iX1 = 1, nX(1)
        ALLOCATE( MergedMeshX2(iX1) % CellMarker(1:nX(2)) )
      END DO

      CALL Determine_MergedCells( nX )

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1

      DO iX1 = 1, nX(1)
        DEALLOCATE( MergedMeshX2(iX1) % CellMarker )
      END DO
      DEALLOCATE( MergedMeshX2 )

      ! NULLIFY( r )

    END SUBROUTINE Finalize_CellMerging

    SUBROUTINE Determine_MergedCells( nX )

      INTEGER, INTENT(in) :: nX(1:3)

        REAL(DP) :: r(1:nX(1))
        INTEGER  :: n, iX1, iX2, CellNumber

        ! --- Merging is assumed to happen in the X2 direction
        !     such that the total number of cells in the X2 direction
        !     is divisible by 2. Equidistant meshes are also assumed ---

        r = MeshX(1) % Center(1:nX(1)) + 0.5_DP * MeshX(1) % Width(1:nX(1))

        DO iX1 = 1, nX(1)

            n = 1

            ! --- Determine how many cells per merge in the X2 direction
            !     and the number of merged cells in the X2 direction ---
            IF ( r(iX1) * MeshX(2) % Width(1) .LE. MeshX(1) % Width(iX1) ) THEN

                DO WHILE (2**n * r(iX1) * MeshX(2) % Width(1) .LE. MeshX(1) % Width(iX1))

                    n = n + 1

                END DO

                MergedMeshX2(iX1) % NCellsPerMerge = 2**n
                MergedMeshX2(iX1) % NCells         = nX(2) / 2**n

            ELSE

                MergedMeshX2(iX1) % NCellsPerMerge = 1
                MergedMeshX2(iX1) % NCells         = nX(2)

            END IF

            ! --- Mark the cells which are merged ---
            CellNumber = 1
            DO iX2 = 1, nX(2), MergedMeshX2(iX1) % NCellsPerMerge

                MergedMeshX2(iX1) % CellMarker(iX2:iX2+MergedMeshX2(iX1) % NCellsPerMerge - 1) &
                 = CellNumber
                
                CellNumber = CellNumber + 1

            END DO

        END DO

        ! Checking to make sure this subroutine does what it should
        ! DO iX1 = 1, 3
        ! DO iX2 = 1, nX(2)

        !   print *, 'MergedMesh(',iX1,') = ', MergedMeshX2(iX1) % CellMarker(iX2)

        ! END DO
        ! END DO

    END SUBROUTINE Determine_MergedCells

    SUBROUTINE MergetoCoarseGrid

    END SUBROUTINE MergetoCoarseGrid

    SUBROUTINE RestricttoFineGrid

    END SUBROUTINE RestricttoFineGrid

END MODULE CellMergingModule