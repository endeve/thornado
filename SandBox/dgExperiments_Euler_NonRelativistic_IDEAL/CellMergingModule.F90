MODULE CellMergingModule

    USE KindModule, ONLY: &
      DP, Zero, Half, One, Two
    USE ProgramHeaderModule, ONLY: &
      iX_B0, iX_E0, iX_B1, iX_E1, nDOFX, nNodesX
    USE GeometryFieldsModule, ONLY: &
      uGF, iGF_h_2, iGF_h_3 ! In spherical coordinates: iGF_h_2 => r, iGF_h_3 => r*sin\theta
    USE FluidFieldsModule, ONLY: &
      nCF
    USE QuadratureModule, ONLY: &
      GetQuadrature
    USE ReferenceElementModuleX, ONLY: &
      NodeNumberTableX
    USE MeshModule, ONLY: &
      MeshX, &
      NodeCoordinate
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
    PUBLIC :: MergeAndRestrict

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
        ALLOCATE( MergedMeshX2(iX1) % &
                  MergedBasisCoeff(1:nN, &
                                   1:nN, &
                                   1:MergedMeshX2(iX1) % NCellsPerMerge) )
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
      INTEGER  :: iNodeX, iNodeX1, iNodeX2
      REAL(DP) :: X1,X2

      ! --- Merging is assumed to happen in the X2 direction
      !     such that the total number of cells in the X2 direction
      !     is divisible by 2. Equidistant meshes are also assumed ---

      r = MeshX(1) % Center(1:nX(1)) + 0.5_DP * MeshX(1) % Width(1:nX(1)) ! It's possible to replace r with uGF?

      DO iX1 = 1, nX(1)

        n = 1

        ! --- Determine how many cells per merge in the X2 direction
        !     and the number of merged cells in the X2 direction ---
        IF ( r(iX1) * MeshX(2) % Width(1) .LE. MeshX(1) % Width(iX1) ) THEN

          DO WHILE (2**n * r(iX1) * MeshX(2) % Width(1) &
                    .LE. &
                    MeshX(1) % Width(iX1))

            n = n + 1

          END DO

          MergedMeshX2(iX1) % NCellsPerMerge = 2**n
          MergedMeshX2(iX1) % NCells         = nX(2) / 2**n

        ELSE

          MergedMeshX2(iX1) % NCellsPerMerge = 1
          MergedMeshX2(iX1) % NCells         = nX(2)

        END IF

        ! --- Mark the cells which are merged according to coarse cell #---
        
        ! CellNumber = 1
        ! DO iX2 = 1, nX(2), MergedMeshX2(iX1) % NCellsPerMerge

        !   MergedMeshX2(iX1) % CellMarker(iX2:iX2+&
        !                                   MergedMeshX2(iX1)%NCellsPerMerge-1)&
        !     = CellNumber
          
        !   CellNumber = CellNumber + 1

        ! END DO

        ! --- Mark the cells which are emerged according to fine cell # within coarse cell ---
        DO iX2 = 1, nX(2)

          MergedMeshX2(iX1) % CellMarker(iX2) &
            = MOD(iX2 - 1, MergedMeshX2(iX1) % NCellsPerMerge) + 1

        END DO

      END DO

      ! DO iX1 = 1, nX(1)
      ! DO iX2 = 1, nX(2)

      !   print *, MergedMeshX2(iX1) % CellMarker(iX2)

      ! END DO

      !   write(*,*)

      ! END DO

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), 1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), 1, iNodeX2 )
        print *, "(",X1,",",X2,")"

      END DO

    END SUBROUTINE Determine_MergedCells

    SUBROUTINE Determine_BasisCoeff ( nX, nN )

      INTEGER, INTENT(in) :: nX(1:3)
      INTEGER, INTENT(in) :: nN

      INTEGER  :: iX1, iCell, iFine, iMerge, iXQ
      REAL(DP) :: xQ_F(nN), wQ_F(nN), xQ_M(nN)
      REAL(DP) :: a, b

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

        DO iFine  = 1, nN
        DO iMerge = 1, nN

          DO iXQ = 1, nN

            ! --- Map ref. quadrature points to the fine cell ---
            xQ_M(iXQ) = (b - a) * xQ_F(iXQ) + (a + b) / Two

            MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) &
              = MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) + &
                L_X2(iMerge) % P(xQ_M(iXQ)) * L_X2(iFine) % P(xQ_F(iXQ)) * &
                wQ_F(iXQ)

          END DO

        END DO
        END DO
        END DO

      END DO

    END SUBROUTINE Determine_BasisCoeff

    SUBROUTINE MergeAndRestrict ( nN, U )

      INTEGER, INTENT(in)     :: nN
      REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                   iX_B1(1):iX_E1(1), &
                                   iX_B1(2):iX_E1(2), &
                                   iX_B1(3):iX_E1(3), &
                                   1:nCF)

      INTEGER  :: iNodeX, iX1, iX2, iX3, iCF
      INTEGER  :: iCell, iFine, iMerge
      INTEGER  :: NodeX(nDOFX) ! will this be necessary?
      REAL(DP) :: xQ(nN), wQ(nN)
      REAL(DP) :: uCF(1:nDOFX, &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(1):iX_E0(1), &
                      iX_B0(3):iX_E0(3), &
                      1:nCF)

      ! The mapping from [-1/2,1/2] to [a,b] is given by f(x) = (b-a)*x + (a+b)/2 for x in [-1/2,1/2].
      ! Need to figure out how to get Lagrange polynomials
      ! Can use the reference polynomials on [-1/2,1/2]; for the merged polynomial, chop respective portion and remap xQ?
      ! i.e. if four cells are merged, cut reference interval into fourths to evaluate merged polynomial

      ! Need to add a structure to contain the evaluated integrals. Do I make one large matrix and fill it entry by entry?
      ! Or do I build submatrices and patch them together?
      ! Need to also make sure the matrix vector multiplication "lines" up correctly, i.e., do I need to permute U?

      ! Most likely I should combine the Merge and Restrict as a single operation

      ! --- Permute node ordering in 2-dimensions ---
      IF( nDOFX .EQ. 4 )THEN

        NodeX(1) = 1
        NodeX(2) = 3
        NodeX(3) = 2
        NodeX(4) = 4

      ELSEIF( nDOFX .EQ. 9)THEN

        NodeX(1) = 1
        NodeX(2) = 4
        NodeX(3) = 7
        NodeX(4) = 2
        NodeX(5) = 1
        NodeX(6) = 8
        NodeX(7) = 3
        NodeX(8) = 6
        NodeX(9) = 9

      END IF

      ! --- Permute conserved quantities ---
      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iX2 = iX_B0(2), iX_E0(2)

        DO iNodeX = 1,nDOFX

          uCF(iNodeX,iX2,iX1,iX3,iCF) = U(NodeX(iNodeX),iX1,iX2,iX3,iCF)

        END DO

      END DO
      END DO
      END DO
      END DO

      ! --- Merge and Restrict conserved quantities ---
      
      CALL GetQuadrature( nN, xQ, wQ )

      ! DO iCF = 1, nCF
      ! DO iX3 = iX_B0(3), iX_E0(3)
      ! DO iX1 = iX_B0(1), iX_E0(1)
      ! DO iX2 = iX_B0(2), iX_E0(2)
      ! DO iNodeX = 1,nDOFX

      !     DO iCell  = 1,MergedMeshX2(iX1) % NCellsPerMerge
      !     DO iFine  = 1,nN
      !     DO iMerge = 1,nN

      !       MergedMeshX2(iX1) % MergedBasisCoeff(iMerge,iFine,iCell) * &
      !       MergedMeshX2(iX1) % MergedBasisCoeff(iMerge, &
      !                                            ?, &
      !                                            MergedMeshX2(iX1) % CellMarker(iX2)) * &
      !       uCF(?,iX2,iX1,iX3,iCF) / wQ(iMerge)

      !     END DO
      !     END DO
      !     END DO

      ! END DO
      ! END DO
      ! END DO
      ! END DO
      ! END DO

    END SUBROUTINE MergeAndRestrict

END MODULE CellMergingModule