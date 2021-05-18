MODULE SubcellReconstructionModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, iX_B0, iX_E0
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    IndLX_Q, L_X1, L_X2, L_X3
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE MeshModule, ONLY: &
    MeshX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSubcellReconstruction
  PUBLIC :: UpdateSubcellReconstruction
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

    INTEGER  :: iX1
    REAL(DP) :: r_l, r_r, r_q1, r_q2

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

!!$ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IF( TRIM( CoordinateSystem ) .eq. 'SPHERICAL' )THEN

    iS = 0
    DO iS3 = 1, nNodesX(3)
    DO iS2 = 1, nNodesX(2)
    DO iS1 = 1, nNodesX(1)

      iS = iS + 1

      jS = 0
      DO jS3 = 1, nNodesX(3)
      DO jS2 = 1, nNodesX(2)
      DO jS1 = 1, nNodesX(1)

        jS = jS + 1

        iX1 = 1
        ! 1D spherical coordinate
        ! use the first element in the block
        IF( iS == 1 )THEN  ! first subgrid cell
          r_l  = MeshX(1) % Center(iX1) - Half*MeshX(1) % Width(iX1)
          r_r  = MeshX(1) % Center(iX1)
        ELSE ! second subgrid cell
          r_l  = MeshX(1) % Center(iX1)
          r_r  = MeshX(1) % Center(iX1) + Half*MeshX(1) % Width(iX1)
        END IF

        IF( jS == 1 )THEN
          r_q1 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(iX1)
          r_q2 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(iX1+1)
        ELSE
          r_q1 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(iX1+1)
          r_q2 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(iX1)
        END IF

        !!$WRITE(*,'(2I4,4ES12.3)') iS, jS, r_l, r_q1, r_q2, r_r

        ProjectionMatrix(iS,jS) &
        =  3.0e0 * ( r_r**4 - r_l**4 ) &
           / ( 4.0e0 * ( r_r**3 - r_l**3 ) ) &
          - r_q2

        ProjectionMatrix(iS,jS) &
        = ProjectionMatrix(iS,jS) &
          / ( r_q1 - r_q2 )

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    END IF

!!$ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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


  SUBROUTINE UpdateSubcellReconstruction( iX )

    INTEGER, INTENT(in) :: iX

    INTEGER  :: iS1, iS2, iS3, iS
    INTEGER  :: jS1, jS2, jS3, jS
    INTEGER  :: iN1, iN2, iN3
    INTEGER  :: INFO
    INTEGER  :: IPIV(nDOFX)
    REAL(DP) :: WORK(nDOFX)

    INTEGER  :: iX1, nNodes
    REAL(DP) :: r_l, r_r, r_q1, r_q2


    IF( TRIM( CoordinateSystem ) .eq. 'SPHERICAL')THEN

      iX1 = iX
      ! Use the mirror physical cell for ghost cell at origin radius
      IF( MeshX(1) % Center( iX ) < 0.0 ) &
        iX1 = iX_B0(1) - iX

    iS = 0
    DO iS3 = 1, nNodesX(3)
    DO iS2 = 1, nNodesX(2)
    DO iS1 = 1, nNodesX(1)

      iS = iS + 1

      jS = 0
      DO jS3 = 1, nNodesX(3)
      DO jS2 = 1, nNodesX(2)
      DO jS1 = 1, nNodesX(1)

        jS = jS + 1

        nNodes = 2 !! ONLY for 2, NEED FIX
        ! 1D spherical coordinate
        IF( iS == 1 )THEN  ! first subgrid cell
          r_l  = MeshX(1) % Center(iX1) - Half*MeshX(1) % Width(iX1)
          r_r  = MeshX(1) % Center(iX1)
        ELSE               ! second subgrid cell
          r_l  = MeshX(1) % Center(iX1)
          r_r  = MeshX(1) % Center(iX1) + Half*MeshX(1) % Width(iX1)
        END IF

        IF( jS == 1 )THEN
          r_q1 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(1)
          r_q2 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(2)
        ELSE
          r_q1 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(2)
          r_q2 = MeshX(1) % Center(iX1) + MeshX(1) % Width(iX1) * MeshX(1) % Nodes(1)
        END IF

        !!$WRITE(*,'(2I4,4ES12.3)') iS, jS, r_l, r_q1, r_q2, r_r

        ProjectionMatrix(iS,jS) &
        =  3.0e0 * ( r_r**4 - r_l**4 ) &
           / ( 4.0e0 * ( r_r**3 - r_l**3 ) ) &
          - r_q2

        ProjectionMatrix(iS,jS) &
        = ProjectionMatrix(iS,jS) &
          / ( r_q1 - r_q2 )

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

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

    END IF ! SPHERICAL

  END SUBROUTINE UpdateSubcellReconstruction


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
