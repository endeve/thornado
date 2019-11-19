MODULE dgDiscretizationModule_Collisions

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nZ, nNodesZ, nDOF
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable, &
    NodeNumberTable4D, &
    Weights_q
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER  :: nE_G, nX_G
  REAL(DP) :: wTime
  REAL(DP) :: N0, SigmaA0, SigmaS0, SigmaT0
  REAL(DP), ALLOCATABLE :: SigmaA(:,:)
  REAL(DP), ALLOCATABLE :: SigmaS(:,:)
  REAL(DP), ALLOCATABLE :: SigmaT(:,:)
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
      dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iCR, iS, iX_G

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

    ! --- Implicit Update ---

    DO iX_G = 1, nX_G
      DO iS = 1, nSpecies

        ! --- Number Density ---

        U_N(:,iCR_N, iS,iX_G) &
          = ( dt * SigmaA(:,iX_G) * N0 + U_N(:,iCR_N,iS,iX_G) ) &
            / ( One + dt * SigmaA(:,iX_G) )

        ! --- Number Flux (1) ---

        U_N(:,iCR_G1,iS,iX_G) &
          = U_N(:,iCR_G1,iS,iX_G) / ( One + dt * SigmaT(:,iX_G) )

        ! --- Number Flux (2) ---

        U_N(:,iCR_G2,iS,iX_G) &
          = U_N(:,iCR_G2,iS,iX_G) / ( One + dt * SigmaT(:,iX_G) )

        ! --- Number Flux (3) ---

        U_N(:,iCR_G3,iS,iX_G) &
          = U_N(:,iCR_G3,iS,iX_G) / ( One + dt * SigmaT(:,iX_G) )

        ! --- Increments ---

        dU_N(:,iCR_N, iS,iX_G) &
          = SigmaA(:,iX_G) * ( N0 - U_N(:,iCR_N, iS,iX_G) )

        dU_N(:,iCR_G1,iS,iX_G) &
          = - SigmaT(:,iX_G) * U_N(:,iCR_G1,iS,iX_G)

        dU_N(:,iCR_G2,iS,iX_G) &
          = - SigmaT(:,iX_G) * U_N(:,iCR_G2,iS,iX_G)

        dU_N(:,iCR_G3,iS,iX_G) &
          = - SigmaT(:,iX_G) * U_N(:,iCR_G3,iS,iX_G)

      END DO
    END DO

    ! --- Map Increment Back ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                      iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dU_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

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
      dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER  :: iCR, iS, iX_G

    ! --- Map Data for Collision Update ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 U_N(1:nE_G,iCR,iS,1:nX_G) )

        IF( alpha_EX > Zero )THEN

          CALL MapForward_R &
                 ( iZ_B0, iZ_E0, &
                   dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                        iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                   dU_N(1:nE_G,iCR,iS,1:nX_G) )

        END IF

      END DO
    END DO

    DO iX_G = 1, nX_G
      DO iS = 1, nSpecies

        ! --- Number Density ---

        U_N(:,iCR_N,iS,iX_G) &
          = U_N(:,iCR_N,iS,iX_G) &
              + alpha_EX * dt**2 * SigmaA(:,iX_G) * dU_N(:,iCR_N,iS,iX_G)

        U_N(:,iCR_N,iS,iX_G) &
          = ( alpha_IM * ( dt * SigmaA(:,iX_G) )**2 * N0 &
              + U_N(:,iCR_N,iS,iX_G) ) &
              / ( One + alpha_IM * ( dt * SigmaA(:,iX_G) )**2 )

        ! --- Number Flux Density (1) ---

        U_N(:,iCR_G1,iS,iX_G) &
          = U_N(:,iCR_G1,iS,iX_G) &
              + alpha_EX * dt**2 * SigmaT(:,iX_G) * dU_N(:,iCR_G1,iS,iX_G)

        U_N(:,iCR_G1,iS,iX_G) &
          = U_N(:,iCR_G1,iS,iX_G) &
              / ( One + alpha_IM * ( dt * SigmaT(:,iX_G) )**2 )

        ! --- Number Flux Density (2) ---

        U_N(:,iCR_G2,iS,iX_G) &
          = U_N(:,iCR_G2,iS,iX_G) &
              + alpha_EX * dt**2 * SigmaT(:,iX_G) * dU_N(:,iCR_G2,iS,iX_G)

        U_N(:,iCR_G2,iS,iX_G) &
          = U_N(:,iCR_G2,iS,iX_G) &
              / ( One + alpha_IM * ( dt * SigmaT(:,iX_G) )**2 )

        ! --- Number Flux Density (3) ---

        U_N(:,iCR_G3,iS,iX_G) &
          = U_N(:,iCR_G3,iS,iX_G) &
              + alpha_EX * dt**2 * SigmaT(:,iX_G) * dU_N(:,iCR_G3,iS,iX_G)

        U_N(:,iCR_G3,iS,iX_G) &
          = U_N(:,iCR_G3,iS,iX_G) &
              / ( One + alpha_IM * ( dt * SigmaT(:,iX_G) )**2 )

        ! --- Corrections ---

        dU_N(:,iCR_N,iS,iX_G) &
          = SigmaA(:,iX_G) &
              * ( alpha_EX * dU_N(:,iCR_N,iS,iX_G) &
                  + alpha_IM * SigmaA(:,iX_G) &
                      * ( N0 - U_N(:,iCR_N, iS,iX_G) ) )

        dU_N(:,iCR_G1,iS,iX_G) &
          = SigmaT(:,iX_G) &
              * ( alpha_EX * dU_N(:,iCR_G1,iS,iX_G) &
                  - alpha_IM * SigmaT(:,iX_G) * U_N(:,iCR_G1,iS,iX_G) )

        dU_N(:,iCR_G2,iS,iX_G) &
          = SigmaT(:,iX_G) &
              * ( alpha_EX * dU_N(:,iCR_G2,iS,iX_G) &
                  - alpha_IM * SigmaT(:,iX_G) * U_N(:,iCR_G2,iS,iX_G) )

        dU_N(:,iCR_G3,iS,iX_G) &
          = SigmaT(:,iX_G) &
              * ( alpha_EX * dU_N(:,iCR_G3,iS,iX_G) &
                  - alpha_IM * SigmaT(:,iX_G) * U_N(:,iCR_G3,iS,iX_G) )

      END DO
    END DO

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

  END SUBROUTINE ComputeCorrection_M1_DG_Implicit


  SUBROUTINE InitializeCollisions &
    ( N0_Option, SigmaA0_Option, SigmaS0_Option, Radius_Option )

    REAL(DP), INTENT(in), OPTIONAL :: N0_Option
    REAL(DP), INTENT(in), OPTIONAL :: SigmaA0_Option
    REAL(DP), INTENT(in), OPTIONAL :: SigmaS0_Option
    REAL(DP), INTENT(in), OPTIONAL :: Radius_Option

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iNode
    REAL(DP) :: X1, X2, X3, R, Radius
    REAL(DP), ALLOCATABLE :: SigmaTmp1(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SigmaTmp2(:,:,:,:,:)

    IF( PRESENT( N0_Option ) )THEN
      N0 = N0_Option
    ELSE
      N0 = Zero
    END IF

    IF( PRESENT( SigmaA0_Option ) )THEN
      SigmaA0 = SigmaA0_Option
    ELSE
      SigmaA0 = Zero
    END IF

    IF( PRESENT( SigmaS0_Option ) )THEN
      SigmaS0 = SigmaS0_Option
    ELSE
      SigmaS0 = Zero
    END IF

    IF( PRESENT( Radius_Option ) )THEN
      Radius = Radius_Option
    ELSE
      Radius = HUGE( One )
    END IF

    SigmaT0 = SigmaA0 + SigmaS0

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializeCollisions'
    WRITE(*,*)
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_A = ', SigmaA0
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_S = ', SigmaS0
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Sigma_T = ', SigmaT0
    WRITE(*,'(A6,A,ES10.4E2)') '', 'Radius  = ', Radius

    nE_G = nNodesZ(1) * nZ(1)
    nX_G = PRODUCT( nNodesZ(2:4) * nZ(2:4) )

    ALLOCATE( SigmaA(nE_G,nX_G) )
    ALLOCATE( SigmaS(nE_G,nX_G) )
    ALLOCATE( SigmaT(nE_G,nX_G) )

    ALLOCATE(  U_N(nE_G,nCR,nSpecies,nX_G) )
    ALLOCATE( dU_N(nE_G,nCR,nSpecies,nX_G) )

    SigmaA = SigmaA0
    SigmaS = SigmaS0
    SigmaT = SigmaT0

    ALLOCATE( SigmaTmp1(nDOF,nZ(1),nZ(2),nZ(3),nZ(4)) )
    ALLOCATE( SigmaTmp2(nDOF,nZ(1),nZ(2),nZ(3),nZ(4)) )

    DO iZ4 = 1, nZ(4)
    DO iZ3 = 1, nZ(3)
    DO iZ2 = 1, nZ(2)
    DO iZ1 = 1, nZ(1)

      IF( TRIM( ProgramName ) .EQ. 'HomogeneousSphere' )THEN

        X1 = MeshX(1) % Center(iZ2)
        X2 = MeshX(2) % Center(iZ3)
        X3 = MeshX(3) % Center(iZ4)

        R = SQRT( X1**2 + X2**2 + X3**2 )

        IF( R < Radius )THEN

          SigmaTmp1(:,iZ1,iZ2,iZ3,iZ4) &
            = SigmaA0 / ( (R/0.85_DP)**40 + One )
          SigmaTmp2(:,iZ1,iZ2,iZ3,iZ4) &
            = SigmaS0

        ELSE

          SigmaTmp1(:,iZ1,iZ2,iZ3,iZ4) &
            = SigmaA0 / ( (R/0.85_DP)**40 + One )
          SigmaTmp2(:,iZ1,iZ2,iZ3,iZ4) &
            = Zero

        END IF

      ELSEIF( TRIM( ProgramName ) .EQ. 'HomogeneousSphere_Spherical' )THEN

        R = MeshX(1) % Center(iZ2)

        ! --- Smoothen Absorption Opacity ---

        SigmaTmp1(:,iZ1,iZ2,iZ3,iZ4) &
          = SigmaA0 / ( ( R / Radius )**80 + One )

        SigmaTmp2(:,iZ1,iZ2,iZ3,iZ4) &
          = Zero

      ELSE

        SigmaTmp1(:,iZ1,iZ2,iZ3,iZ4) = SigmaA0
        SigmaTmp2(:,iZ1,iZ2,iZ3,iZ4) = SigmaS0

      END IF

    END DO
    END DO
    END DO
    END DO

    CALL MapForward_R( [1,1,1,1], nZ, SigmaTmp1, SigmaA )
    CALL MapForward_R( [1,1,1,1], nZ, SigmaTmp2, SigmaS )

    SigmaT = SigmaA + SigmaS

    DEALLOCATE( SigmaTmp1 )
    DEALLOCATE( SigmaTmp2 )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( SigmaA, SigmaS, SigmaT )

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

        iNode = NodeNumberTable4D( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )

        RF_N(iE_G,iX_G) = RF(iNode,iZ1,iZ2,iZ3,iZ4)

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

        iNode = NodeNumberTable4D( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )

        RF(iNode,iZ1,iZ2,iZ3,iZ4) = RF_N(iE_G,iX_G)

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
