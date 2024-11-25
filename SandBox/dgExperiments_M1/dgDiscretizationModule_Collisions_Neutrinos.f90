MODULE dgDiscretizationModule_Collisions_Neutrinos

  USE KindModule, ONLY: &
    DP, One
  USE ProgramHeaderModule, ONLY: &
    nNodesZ
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE NeutrinoOpacitiesModule, ONLY: &
    f_EQ, opEC, opES

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeCollisions
  PUBLIC :: FinalizeCollisions
  PUBLIC :: ComputeIncrement_M1_DG_Implicit
  PUBLIC :: ComputeCorrection_M1_DG_Implicit

  INTEGER  :: nE_G, nX_G
  REAL(DP) :: wTime
  REAL(DP), ALLOCATABLE ::  R_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dR_N(:,:,:,:)

CONTAINS


  SUBROUTINE ComputeIncrement_M1_DG_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_R, dU_R )

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
      U_R (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iCR, iS, iX_G

    wTime = MPI_WTIME( )

    CALL InitializeCollisions( iZ_B0, iZ_E0 )

    wTime = MPI_WTIME( ) - wTime

!    WRITE(*,'(A4,A32,ES10.4E2)') '', 'InitializeCollisions: ', wTime

    ! --- Map Radiation Data for Collision Update ---

    wTime = MPI_WTIME( )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U_R(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 R_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    wTime = MPI_WTIME( ) - wTime

!    WRITE(*,'(A4,A32,ES10.4E2)') '', 'MapForward_R: ', wTime

    ! --- Implicit Update ---

    DO iX_G = 1, nX_G
      DO iS = 1, nSpecies

        ! --- Number Density ---

        R_N(:,iCR_N,iS,iX_G) &
          = ( dt * opEC(:,iS,iX_G) * f_EQ(:,iS,iX_G) &
              + R_N(:,iCR_N,iS,iX_G) ) / ( One + dt * opEC(:,iS,iX_G) )

        ! --- Number Flux (1) ---

        R_N(:,iCR_G1,iS,iX_G) &
          = R_N(:,iCR_G1,iS,iX_G) &
              / ( One + dt * ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) )

        ! --- Number Flux (2) ---

        R_N(:,iCR_G2,iS,iX_G) &
          = R_N(:,iCR_G2,iS,iX_G) &
              / ( One + dt * ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) )

        ! --- Number Flux (3) ---

        R_N(:,iCR_G3,iS,iX_G) &
          = R_N(:,iCR_G3,iS,iX_G) &
              / ( One + dt * ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) )

        ! --- Increments ---

        dR_N(:,iCR_N,iS,iX_G) &
          = opEC(:,iS,iX_G) * ( f_EQ(:,iS,iX_G) - R_N(:,iCR_N,iS,iX_G) )

        dR_N(:,iCR_G1,iS,iX_G) &
          = - ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) * R_N(:,iCR_G1,iS,iX_G)

        dR_N(:,iCR_G2,iS,iX_G) &
          = - ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) * R_N(:,iCR_G2,iS,iX_G)

        dR_N(:,iCR_G3,iS,iX_G) &
          = - ( opEC(:,iS,iX_G) + opES(:,iS,iX_G) ) * R_N(:,iCR_G3,iS,iX_G)

      END DO
    END DO

    ! --- Map Increment Back ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU_R(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                        iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dR_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    CALL FinalizeCollisions

  END SUBROUTINE ComputeIncrement_M1_DG_Implicit


  SUBROUTINE ComputeCorrection_M1_DG_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, alpha_EX, alpha_IM, GE, GX, U, dU )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt, &
      alpha_EX, &
      alpha_IM
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'WARNING: ComputeCorrection_M1_DG_Implicit not implemented'
    WRITE(*,*)

  END SUBROUTINE ComputeCorrection_M1_DG_Implicit


  SUBROUTINE InitializeCollisions( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    INTEGER :: nZ(4)

    nZ = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    ALLOCATE(  R_N(nE_G,nCR,nSpecies,nX_G) )
    ALLOCATE( dR_N(nE_G,nCR,nSpecies,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( R_N, dR_N )

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


END MODULE dgDiscretizationModule_Collisions_Neutrinos
