MODULE CellMergingModule

    USE KindModule, ONLY: &
      DP, Zero, Half, One, Two
    USE ProgramHeaderModule, ONLY: &
      iX_B0, iX_E0, nDOFX, nNodesX
    USE GeometryFieldsModule, ONLY: &
      uGF, iGF_h_2, iGF_h_3 ! In spherical coordinates: iGF_h_2 => r, iGF_h_3 => r*sin\theta
    USE QuadratureModule, ONLY: &
      GetQuadrature
    USE MeshModule, ONLY: &
      MeshX
    USE PolynomialBasisModuleX_Lagrange, ONLY: &
      L_X1, L_X2, L_X3

    TYPE, PUBLIC :: MergedMeshType
      INTEGER :: NCellsPerMerge
      INTEGER :: NCells
      INTEGER, DIMENSION(:), ALLOCATABLE :: CellMarker
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: MergedBasisCoeff
    END type MergedMeshType
    
    TYPE(MergedMeshType), DIMENSION(:), ALLOCATABLE, PUBLIC :: MergedMeshX2 ! X2 merging

    PUBLIC :: Initialize_CellMerging
    PUBLIC :: Finalize_CellMerging
    PUBLIC :: Determine_MergedCells
    PUBLIC :: Determine_BasisCoeff
    PUBLIC :: MergetoCoarseGrid
    PUBLIC :: RestricttoFineGrid

CONTAINS

    SUBROUTINE Initialize_CellMerging( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER :: iX1

      ! --- Create MergedMesh for X2 direction ---
      ALLOCATE( MergedMeshX2(1:nX(1)) )
      MergedMeshX2 % NCellsPerMerge = 1
      MergedMeshX2 % NCells = nX(2)
      DO iX1 = 1, nX(1)
        ALLOCATE( MergedMeshX2(iX1) % CellMarker(1:nX(2)) )
      END DO

      CALL Determine_MergedCells( nX )

      DO iX1 = 1, nX(1)
        ALLOCATE( MergedMeshX2(iX1) % MergedBasisCoeff(1:nN,1:nN,1:MergedMeshX2(iX1) % NCellsPerMerge) )
      END DO

      CALL Determine_BasisCoeff( nX, nN )

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging( nX )

      INTEGER, INTENT(in) :: nX(1:3)

      INTEGER :: iX1

      DO iX1 = 1, nX(1)
        DEALLOCATE( MergedMeshX2(iX1) % CellMarker )
      END DO
      DEALLOCATE( MergedMeshX2 )

    END SUBROUTINE Finalize_CellMerging

    SUBROUTINE Determine_MergedCells( nX )

      INTEGER, INTENT(in) :: nX(1:3)

        REAL(DP) :: r(1:nX(1))
        INTEGER  :: n, iX1, iX2, CellNumber

        ! --- Merging is assumed to happen in the X2 direction
        !     such that the total number of cells in the X2 direction
        !     is divisible by 2. Equidistant meshes are also assumed ---

        r = MeshX(1) % Center(1:nX(1)) + 0.5_DP * MeshX(1) % Width(1:nX(1)) ! It's possible to replace r with uGF?

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

        ! DO iX1 = 1, nX(1)
        !   print *, uGF(1,iX1,1,1,iGF_h_2), uGF(2,iX1,1,1,iGF_h_2), uGF(3,iX1,1,1,iGF_h_2),uGF(4,iX1,1,1,iGF_h_2), r(iX1)
        ! END DO

        print *, L_X1(1) % P(-0.5_DP)
        print *, L_X1(2) % P(-0.5_DP)
        write(*,*)
        print *, L_X2(1) % P(-0.5_DP)
        print *, L_X2(2) % P(-0.5_DP)
        write(*,*)
        print *, L_X3(1) % P(-2.0_DP) ! input parameter to L_X# is nNodesX(#)
        write(*,*) nNodesX
        write(*,*) nDOFX

    END SUBROUTINE Determine_MergedCells

    SUBROUTINE Determine_BasisCoeff ( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER  :: iX1, iCell, iFine, iMerge, iXQ
      REAL(DP) :: xQ_F(nN), wQ_F(nN), xQ_M(nN)
      REAL(DP) :: a, b

      ! Initialize Basis Coefficients
      DO iX1 = 1 ,nX(1)

        DO iCell  = 1, MergedMeshX2(iX1) % NCellsPerMerge
        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN

            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
              = Zero

          END DO

        END DO
        END DO
        END DO

      END DO

      CALL GetQuadrature( nN, xQ_F, wQ_F )

      DO iX1 = 1 ,nX(1)

        DO iCell  = 1, MergedMeshX2(iX1) % NCellsPerMerge

          a = -Half + (REAL(iCell,DP) - One) / REAL(MergedMeshX2(iX1) % NCellsPerMerge,DP)
          b = a + One / REAL(MergedMeshX2(iX1) % NCellsPerMerge,DP)

        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN

            xQ_M(iXQ) = (b - a) * xQ_F(iXQ) + (a + b) / Two

            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
              = MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) + &
                L_X2(iMerge) % P(xQ_M(iXQ)) * L_X2(iFine) % P(xQ_F(iXQ)) * wQ_F(iXQ)

          END DO

        END DO
        END DO
        END DO

      END DO

    END SUBROUTINE Determine_BasisCoeff

    SUBROUTINE MergetoCoarseGrid

      print *, MeshX(1) % Nodes ! These are the Gauss nodes on [-1/2, 1/2]. Map these to evaluate integrals.

      ! The mapping from [-1/2,1/2] to [a,b] is given by f(x) = (b-a)*x + (a+b)/2 for x in [-1/2,1/2].
      ! Need to figure out how to get Lagrange polynomials
      ! Can use the reference polynomials on [-1/2,1/2]; for the merged polynomial, chop respective portion and remap xQ?
      ! i.e. if four cells are merged, cut reference interval into fourths to evaluate merged polynomial

      ! Need to add a structure to contain the evaluated integrals. Do I make one large matrix and fill it entry by entry?
      ! Or do I build submatrices and patch them together?
      ! Need to also make sure the matrix vector multiplication "lines" up correctly, i.e., do I need to permute U?

      ! Most likely I should combine the Merge and Restrict as a single operation

    END SUBROUTINE MergetoCoarseGrid

    SUBROUTINE RestricttoFineGrid

    END SUBROUTINE RestricttoFineGrid

END MODULE CellMergingModule