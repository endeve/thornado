MODULE SubcellReconstructionModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSubcellReconstruction
  PUBLIC :: FinalizeSubcellReconstruction

  REAL(DP), ALLOCATABLE, PUBLIC :: ProjectionMatrix(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: ReconstructionMatrix(:,:)

CONTAINS


  SUBROUTINE InitializeSubcellReconstruction

    INTEGER  :: iS1, iS2, iS3, iS
    INTEGER  :: jS1, jS2, jS3, jS
    INTEGER  :: iN1, iN2, iN3
    INTEGER  :: INFO
    INTEGER  :: IPIV(nDOFX)
    REAL(DP) :: WORK(nDOFX)
    REAL(DP) :: EtaLo(3), dEta(3)
    REAL(DP) :: Eta_X1(nNodesX(1))
    REAL(DP) :: Eta_X2(nNodesX(2))
    REAL(DP) :: Eta_X3(nNodesX(3))

    ALLOCATE( ProjectionMatrix(nDOFX,nDOFX) )

    ProjectionMatrix = Zero

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

        DO iN3 = 1, nNodesX(3)
        DO iN2 = 1, nNodesX(2)
        DO iN1 = 1, nNodesX(1)

          ProjectionMatrix(iS,jS) &
            = ProjectionMatrix(iS,jS) &
                + WeightsX1(iN1) * WeightsX2(iN2) * WeightsX3(iN3) &
                    * L_X1(IndLX_Q(1,jS)) % P( Eta_X1(iN1) ) &
                    * L_X2(IndLX_Q(2,jS)) % P( Eta_X2(iN2) ) &
                    * L_X3(IndLX_Q(3,jS)) % P( Eta_X3(iN3) )

        END DO
        END DO
        END DO

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    ALLOCATE( ReconstructionMatrix(nDOFX,nDOFX) )

    ! --- Compute Reconstruction Matrix by Inverting ProjectionMatrix ---

    ReconstructionMatrix = ProjectionMatrix

    CALL DGETRF( nDOFX, nDOFX, ReconstructionMatrix, nDOFX, IPIV, INFO )

    IF( INFO /= 0 )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'InitializeSubcellReconstruction:'
      WRITE(*,*)
      WRITE(*,'(A6,A,I4.4)') '', 'DGETRF returned INFO = ', INFO
      WRITE(*,*)
      STOP
    END IF

    CALL DGETRI( nDOFX, ReconstructionMatrix, nDOFX, IPIV, WORK, nDOFX, INFO )

    IF( INFO /= 0 )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'InitializeSubcellReconstruction:'
      WRITE(*,*)
      WRITE(*,'(A6,A,I4.4)') '', 'DGETRI returned INFO = ', INFO
      WRITE(*,*)
      STOP
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ProjectionMatrix, ReconstructionMatrix )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ProjectionMatrix, ReconstructionMatrix )
#endif

  END SUBROUTINE InitializeSubcellReconstruction


  SUBROUTINE FinalizeSubcellReconstruction

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ProjectionMatrix, ReconstructionMatrix )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ProjectionMatrix, ReconstructionMatrix )
#endif

    DEALLOCATE( ProjectionMatrix )
    DEALLOCATE( ReconstructionMatrix )

  END SUBROUTINE FinalizeSubcellReconstruction


END MODULE SubcellReconstructionModule
