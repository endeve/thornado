MODULE dgDiscretizationModule_Collisions

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, nDOF
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER  :: nE_G, nX_G
  REAL(DP) :: N0, SigmaA, SigmaS, SigmaT
  REAL(DP) :: wTime
  REAL(DP), ALLOCATABLE ::  U_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dU_N(:,:,:,:)

  PUBLIC :: InitializeCollisions
  PUBLIC :: FinalizeCollisions
  PUBLIC :: ComputeIncrement_M1_DG_Implicit
  PUBLIC :: ComputeCorrection_M1_DG_Implicit

CONTAINS


  SUBROUTINE ComputeIncrement_M1_DG_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: iCR, iS, iX_G
    REAL(DP) :: GammaA, GammaT

!!$    PRINT*, "ComputeIncrement_M1_DG_Implicit"

    ! --- Map Data for Collision Update ---

    wTime = MPI_WTIME( )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 U_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    wTime = MPI_WTIME( ) - wTime

!!$    PRINT*
!!$    PRINT*, " Forward: ", wTime

    ! --- Implicit Update ---

    wTime = MPI_WTIME( )

    GammaA = dt * SigmaA
    GammaT = dt * SigmaT

    !$OMP PARALLEL DO PRIVATE ( iX_G, iS )
    DO iX_G = 1, nX_G
      DO iS = 1, nSpecies

        ! --- Number Density ---

        U_N(:,iCR_N, iS,iX_G) &
          = ( GammaA * N0 + U_N(:,iCR_N, iS,iX_G) ) / ( One + GammaA )

        ! --- Number Flux (1) ---

        U_N(:,iCR_G1,iS,iX_G) &
          = U_N(:,iCR_G1,iS,iX_G) / ( One + GammaT )

        ! --- Number Flux (2) ---

        U_N(:,iCR_G2,iS,iX_G) &
          = U_N(:,iCR_G2,iS,iX_G) / ( One + GammaT )

        ! --- Number Flux (3) ---

        dU_N(:,iCR_G3,iS,iX_G) &
          = U_N(:,iCR_G3,iS,iX_G) / ( One + GammaT )

        ! --- Increments ---

        dU_N(:,iCR_N, iS,iX_G) &
          = SigmaA * ( N0 - U_N(:,iCR_N, iS,iX_G) )

        dU_N(:,iCR_G1,iS,iX_G) &
          = - SigmaT * U_N(:,iCR_G1,iS,iX_G)

        dU_N(:,iCR_G2,iS,iX_G) &
          = - SigmaT * U_N(:,iCR_G2,iS,iX_G)

        dU_N(:,iCR_G3,iS,iX_G) &
          = - SigmaT * U_N(:,iCR_G3,iS,iX_G)

      END DO
    END DO
    !$OMP END PARALLEL DO

    wTime = MPI_WTIME( ) - wTime

!!$    PRINT*
!!$    PRINT*, " Update: ", wTime

    ! --- Map Increment Back ---

    wTime = MPI_WTIME( )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                      iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dU_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

!!$    print*,"MAXVAL(dU) = ", MAXVAL(ABS(dU))

    wTime = MPI_WTIME( ) - wTime

!!$    PRINT*
!!$    PRINT*, "Backward: ", wTime

  END SUBROUTINE ComputeIncrement_M1_DG_Implicit


  SUBROUTINE ComputeCorrection_M1_DG_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt2, GE, GX, U, dU )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt2
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: iCR, iS, iX_G
    REAL(DP) :: GammaA2, GammaT2

    ! --- Map Data for Collision Update ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 U_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    GammaA2 = dt2 * SigmaA**2
    GammaT2 = dt2 * SigmaT**2

    !$OMP PARALLEL DO PRIVATE ( iX_G, iS )
    DO iX_G = 1, nX_G
      DO iS = 1, nSpecies

        ! --- Number Density ---

        U_N(:,iCR_N, iS,iX_G) &
          = ( GammaA2 * N0 + U_N(:,iCR_N, iS,iX_G) ) / ( One + GammaA2 )

        ! --- Number Flux Density (1) ---

        U_N(:,iCR_G1,iS,iX_G) &
          = U_N(:,iCR_G1,iS,iX_G) / ( One + GammaT2 )

        ! --- Number Flux Density (2) ---

        U_N(:,iCR_G2,iS,iX_G) &
          = U_N(:,iCR_G2,iS,iX_G) / ( One + GammaT2 )

        ! --- Number Flux Density (3) ---

        U_N(:,iCR_G3,iS,iX_G) &
          = U_N(:,iCR_G3,iS,iX_G) / ( One + GammaT2 )

        ! --- Corrections ---

        dU_N(:,iCR_N, iS,iX_G) &
          = - SigmaA**2 * ( N0 - U_N(:,iCR_N, iS,iX_G) )

        dU_N(:,iCR_G1,iS,iX_G) &
          = SigmaT**2 * U_N(:,iCR_G1,iS,iX_G)

        dU_N(:,iCR_G2,iS,iX_G) &
          = SigmaT**2 * U_N(:,iCR_G2,iS,iX_G)

        dU_N(:,iCR_G3,iS,iX_G) &
          = SigmaT**2 * U_N(:,iCR_G3,iS,iX_G)

      END DO
    END DO
    !$OMP END PARALLEL DO

    ! --- Map Correction Back ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                      iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dU_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

!!$    print*,"MAXVAL(dC) = ", MAXVAL(ABS(dU))

  END SUBROUTINE ComputeCorrection_M1_DG_Implicit


  SUBROUTINE InitializeCollisions( N0_Option, SigmaA_Option, SigmaS_Option )

    REAL(DP), INTENT(in), OPTIONAL :: N0_Option
    REAL(DP), INTENT(in), OPTIONAL :: SigmaA_Option
    REAL(DP), INTENT(in), OPTIONAL :: SigmaS_Option

    N0 = Zero
    IF( PRESENT( N0_Option ) ) &
      N0 = N0_Option

    SigmaA = Zero
    IF( PRESENT( SigmaA_Option ) ) &
      SigmaA = SigmaA_Option

    SigmaS = Zero
    IF( PRESENT( SigmaS_Option ) ) &
      SigmaS = SigmaS_Option

    SigmaT = SigmaA + SigmaS

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializeCollisions'
    WRITE(*,*)
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_A = ', SigmaA
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_S = ', SigmaS
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_T = ', SigmaT
    WRITE(*,*)

    nE_G = nNodesZ(1) * nZ(1)
    nX_G = PRODUCT( nNodesZ(2:4) * nZ(2:4) )

    ALLOCATE(  U_N(nE_G,nCR,nSpecies,nX_G) )
    ALLOCATE( dU_N(nE_G,nCR,nSpecies,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( U_N, dU_N )

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE MapForward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      RF(:,:,:,:,:)
    REAL(DP), INTENT(out) :: &
      RF_N(:,:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4, iNode
    INTEGER :: iE_G, iX_G

    iX_G = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
      DO iZ3 = iZ_B(3), iZ_E(3)
        DO iZ2 = iZ_B(2), iZ_E(2)
          DO iNodeZ4 = 1, nNodesZ(4)
            DO iNodeZ3 = 1, nNodesZ(3)
              DO iNodeZ2 = 1, nNodesZ(2)

                iX_G = iX_G + 1

                iE_G = 0
                DO iZ1 = iZ_B(1), iZ_E(1)
                  DO iNodeZ1 = 1, nNodesZ(1)

                    iE_G = iE_G + 1

                    iNode = NodeNumberTable4D &
                              ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )

                    RF_N(iE_G,iX_G) &
                      = RF(iNode,iZ1,iZ2,iZ3,iZ4)

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapForward_R


  SUBROUTINE MapBackward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(out) :: &
      RF(:,:,:,:,:)
    REAL(DP), INTENT(in)  :: &
      RF_N(:,:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4, iNode
    INTEGER :: iE_G, iX_G

    iX_G = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
      DO iZ3 = iZ_B(3), iZ_E(3)
        DO iZ2 = iZ_B(2), iZ_E(2)
          DO iNodeZ4 = 1, nNodesZ(4)
            DO iNodeZ3 = 1, nNodesZ(3)
              DO iNodeZ2 = 1, nNodesZ(2)

                iX_G = iX_G + 1

                iE_G = 0
                DO iZ1 = iZ_B(1), iZ_E(1)
                  DO iNodeZ1 = 1, nNodesZ(1)

                    iE_G = iE_G + 1

                    iNode = NodeNumberTable4D &
                              ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )

                    RF(iNode,iZ1,iZ2,iZ3,iZ4) &
                      = RF_N(iE_G,iX_G)

                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapBackward_R


END MODULE dgDiscretizationModule_Collisions
