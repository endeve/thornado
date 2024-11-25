MODULE dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, &
    nNodesE, nDOFE, &
    nDOF
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, &
    nDOF_X2, &
    Weights_q, &
    Weights_X1, &
    Weights_X2, &
    NodeNumberTable, &
    OuterProduct1D3D
  USE ReferenceElementModule_Lagrange, ONLY: &
    dLdX1_q, &
    dLdX2_q, &
    L_X1_Dn, &
    L_X1_Up, &
    L_X2_Dn, &
    L_X2_Up
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
    iGE_Ep0, &
    iGE_Ep2
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE ClosureModule_M1, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE BoundaryConditionsModule_Beta, ONLY: &
    ApplyBoundaryConditions_Radiation
  USE MomentEquationsUtilitiesModule_Beta, ONLY: &
    ComputePrimitive, &
    Flux_X1, &
    Flux_X2, &
    NumericalFlux_LLF

  USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  REAL(DP), PARAMETER :: Ones(16) = One

  PUBLIC :: ComputeIncrement_M1_DG_Explicit

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer
  REAL(DP) :: Timer_RHS
  REAL(DP) :: Timer_VOL, dT_VOL
  REAL(DP) :: Timer_AD1, dT_AD1
  REAL(DP) :: Timer_SUR, dT_SUR
  REAL(DP) :: Timer_INT, dT_INT
  REAL(DP) :: Timer_INT_G, dT_INT_G
  REAL(DP) :: Timer_LFT, dT_LFT
  REAL(DP) :: Timer_RGT, dT_RGT
  REAL(DP) :: Timer_FLX, dT_FLX
  REAL(DP) :: Timer_AD2, dT_AD2
  REAL(DP) :: Timer_AD3, dT_AD3
  REAL(DP) :: Timer_Div_X2
  REAL(DP) :: Timer_INV

CONTAINS


  SUBROUTINE ComputeIncrement_M1_DG_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS
    REAL(DP) :: dZ(4), Tau(nDOF)

    dU  = Zero

    CALL ApplyBoundaryConditions_Radiation &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

!!$    CALL ComputeIncrement_Divergence_X1 &
!!$           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

!!    CALL ComputeIncrement_Divergence_X1_GPU &
!!          ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL ComputeIncrement_Divergence_X1_GPU_TEST &
          ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )



    ! --- Multiply Inverse Mass Matrix ---

    CALL Timer_Start( Timer_INV )

    !$OMP PARALLEL DO PRIVATE &
    !$OMP& ( iZ1, iZ2, iZ3, iZ4, iCR, iS, dZ, Tau )
    DO iZ4 = iZ_B0(4), iZ_E0(4)
      dZ(4) = MeshX(3) % Width(iZ4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
        dZ(3) = MeshX(2) % Width(iZ3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
          dZ(2) = MeshX(1) % Width(iZ2)
          DO iZ1 = iZ_B0(1), iZ_E0(1)
            dZ(1) = MeshE % Width(iZ1)

            Tau(1:nDOF) &
              = OuterProduct1D3D &
                  ( GE(:,iZ1,iGE_Ep2), nDOFE, &
                    GX(:,iZ2,iZ3,iZ4,iGF_SqrtGm), nDOFX )

            DO iS = 1, nSpecies
              DO iCR = 1, nCR

                dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                  = dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                      / ( Weights_q(:) * Tau(:) * PRODUCT( dZ ) )

              END DO
            END DO

          END DO
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

    CALL Timer_Stop( Timer_INV )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Inverse Mass: ', Timer_INV
      WRITE(*,*)

    END IF

  END SUBROUTINE ComputeIncrement_M1_DG_Explicit
 
  SUBROUTINE ComputeIncrement_Divergence_X1_GPU &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: nK, iK, nF, iF, nZ(4), iNode
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iCR
    
    INTEGER  :: size_range(1:4)
    INTEGER  :: Ones_size(1:4)

    REAL(DP) :: wTime, start_t, end_t
    REAL(DP) :: primitive_t, ff_t, ef_t, flux_x1_q_t
    REAL(DP) :: copy_prev_t, num_flux_t, f_dgemm_t
    REAL(DP) :: l_flux_t, r_flux_t, s_num_flux_t
    REAL(DP) :: s_dgemm_t, update_t, total_t
    
    REAL(DP) :: Ones_q(nDOF), Ones_RL(nDOF_X1), dZ(4)
    REAL(DP) :: Tau(1:nDOF), Tau_X1(1:nDOF_X1)
    REAL(DP) :: P_q(1:nDOF,1:nPR), FF_q(1:nDOF), EF_q(1:nDOF)
    REAL(DP), ALLOCATABLE :: Flux_X1_q(:,:,:), dU_q(:,:,:)
    REAL(DP), ALLOCATABLE :: U_P(:,:,:), U_K(:,:,:)
    REAL(DP), ALLOCATABLE :: U_L(:,:,:), Flux_L(:,:,:)
    REAL(DP), ALLOCATABLE :: U_R(:,:,:), Flux_R(:,:,:)
    REAL(DP), ALLOCATABLE :: NumericalFlux(:,:,:), dU_P(:,:,:), dU_K(:,:,:)
    REAL(DP) :: P_L(1:nDOF_X1,1:nPR), FF_L(1:nDOF_X1), EF_L(1:nDOF_X1)
    REAL(DP) :: P_R(1:nDOF_X1,1:nPR), FF_R(1:nDOF_X1), EF_R(1:nDOF_X1)

    nZ = iZ_E0 - iZ_B0 + 1
    nK = PRODUCT( nZ )
    nF = PRODUCT( nZ + [0,1,0,0] )

    Ones_q = One

    Ones_RL = One

    Ones_size = One

    ALLOCATE( Flux_X1_q(1:nDOF,1:nK,1:nCR) )
    ALLOCATE( dU_q     (1:nDOF,1:nK,1:nCR) )

    ALLOCATE( U_P(1:nDOF,1:nF,1:nCR) )
    ALLOCATE( U_K(1:nDOF,1:nF,1:nCR) )

    ALLOCATE( U_L(1:nDOF_X1,1:nF,1:nCR) )
    ALLOCATE( U_R(1:nDOF_X1,1:nF,1:nCR) )

    ALLOCATE( Flux_L(1:nDOF_X1,1:nF,1:nPR) )
    ALLOCATE( Flux_R(1:nDOF_X1,1:nF,1:nPR) )

    ALLOCATE( NumericalFlux(1:nDOF_X1,1:nF,1:nCR) )

    ALLOCATE( dU_P(1:nDOF,1:nF,1:nCR) )
    ALLOCATE( dU_K(1:nDOF,1:nF,1:nCR) )

    ! Init timing summation variables
    primitive_t  = 0.0d0
    ff_t         = 0.0d0
    ef_t         = 0.0d0
    flux_x1_q_t  = 0.0d0
    copy_prev_t  = 0.0d0
    num_flux_t   = 0.0d0
    f_dgemm_t    = 0.0d0
    l_flux_t     = 0.0d0
    r_flux_t     = 0.0d0
    s_num_flux_t = 0.0d0
    s_dgemm_t    = 0.0d0
    update_t     = 0.0d0

    total_t = MPI_Wtime()

    DO iS = 1, nSpecies

      iK = 0
      iF = 0

      wTime = MPI_WTIME( )

      ! Added to accomodate ranges and make iK dependent on iZ
      size_range = iZ_E0 - iZ_B0 + Ones_size

      !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(iZ4, iZ3, iZ2, iZ1, dZ, &
      !$OMP& FF_q, EF_q, Tau, Tau_X1, P_q) FIRSTPRIVATE(iK, iF)
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dZ(1) = MeshE    % Width(iZ1)
              dZ(2) = MeshX(1) % Width(iZ2)
              dZ(3) = MeshX(2) % Width(iZ3)
              dZ(4) = MeshX(3) % Width(iZ4)

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, Ones_q(1:nDOFX), nDOFX )

              Tau_X1(1:nDOF_X1) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, Ones_q(1:nDOFX_X1), nDOFX_X1 )

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                !iK = iK + 1
                iK = (iZ1 - iZ_B0(1) + 1) + (iZ2 - iZ_B0(2)) * size_range(1) +&
                        (iZ3 - iZ_B0(3)) * size_range(2) * size_range(1) + &
                        (iZ4 - iZ_B0(4)) * size_range(3) * size_range(2) * size_range(1)

                !PRINT *, "iK: ", iK
                !PRINT *, "iZ1: ", iZ1, "iZ_B0(1): ", iZ_B0(1), "iZ_E0(1): ", iZ_E0(1)
                !PRINT *, "iZ2: ", iZ2, "iZ_B0(2): ", iZ_B0(2), "iZ_E0(2): ", iZ_E0(2)
                !PRINT *, "iZ3: ", iZ3, "iZ_B0(3): ", iZ_B0(3), "iZ_E0(3): ", iZ_E0(3)
                !PRINT *, "iZ4: ", iZ4, "iZ_B0(4): ", iZ_B0(4), "iZ_E0(4): ", iZ_E0(4)
                !PRINT *, "Sizes: ", size_range

                ! Timing
                start_t = MPI_Wtime()


                CALL ComputePrimitive &
                       ( U(:,iZ1,iZ2,iZ3,iZ4,iCR_N, iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                         P_q(:,iPR_D ), P_q(:,iPR_I1), &
                         P_q(:,iPR_I2), P_q(:,iPR_I3), &
                         Ones_q(:), Ones_q(:), Ones_q(:) )

                ! Timing
                primitive_t = primitive_t + MPI_Wtime() - start_t

                start_t = MPI_Wtime()


                FF_q = FluxFactor &
                         ( P_q(:,iPR_D ), P_q(:,iPR_I1), &
                           P_q(:,iPR_I2), P_q(:,iPR_I3), &
                           Ones_q(:), Ones_q(:), Ones_q(:) )

                ! Timing
                ff_t = ff_t + MPI_Wtime() - start_t

                start_t = MPI_Wtime()

                
                EF_q = EddingtonFactor( P_q(:,iPR_D), FF_q )

                ! Timing
                ef_t = ef_t + MPI_Wtime() - start_t

                start_t = MPI_Wtime()

                DO iNode = 1, nDOF

                  Flux_X1_q(iNode,iK,1:nCR) &
                    = Flux_X1 &
                        ( P_q(iNode,iPR_D ), P_q(iNode,iPR_I1), &
                          P_q(iNode,iPR_I2), P_q(iNode,iPR_I3), &
                          FF_q(iNode), EF_q(iNode), One, One, One )

                END DO

                DO iCR = 1, nCR

                  Flux_X1_q(:,iK,iCR) &
                    = dZ(1) * dZ(3) * dZ(4) * Weights_q(:) &
                        * Tau(:) * Flux_X1_q(:,iK,iCR)

                END DO

                ! Timing
                flux_x1_q_t = flux_x1_q_t + MPI_Wtime() - start_t

              END IF

              ! Timing
              start_t = MPI_Wtime()

              !iF = iF + 1
              iF = (iZ1 - iZ_B0(1) + 1) + (iZ2 - iZ_B0(2)) * size_range(1) +&
                        (iZ3 - iZ_B0(3)) * (size_range(2) + 1) * size_range(1) + &
                        (iZ4 - iZ_B0(4)) * size_range(3) * (size_range(2) + 1) * size_range(1)

              !PRINT *, "iF: ", iF

              U_P(:,iF,1:nCR) = U(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS)
              U_K(:,iF,1:nCR) = U(:,iZ1,iZ2,  iZ3,iZ4,1:nCR,iS)

              ! Timing
              copy_prev_t = copy_prev_t + MPI_Wtime() - start_t

              start_t = MPI_Wtime()

              DO iCR = 1, nCR

                NumericalFlux(:,iF,iCR) &
                  = dZ(1) * dZ(3) * dZ(4) * Weights_X1(:) * Tau_X1(:)

              END DO

              ! Timing
              num_flux_t = num_flux_t + MPI_Wtime() - start_t

            END DO
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO

      wTime = MPI_WTIME( ) - wTime

      PRINT *, "Prep Time: ", wTime

      !PRINT *, "Printing prep time proportions:"
      PRINT *, "primitive: ", primitive_t, "ff: ", ff_t
      PRINT *, "ef: ", ef_t
      PRINT *, "flux_x1_q: ", flux_x1_q_t, "copy left/right: ", copy_prev_t
      PRINT *, "numerical flux: ", num_flux_t
      PRINT *, ""

      wTime = MPI_WTIME( )

      start_t = MPI_Wtime()

      DO iCR = 1, nCR

        CALL DGEMM &
               ( 'T', 'N', nDOF, nK, nDOF, One, dLdX1_q(:,:), nDOF, &
                 Flux_X1_q(:,:,iCR), nDOF, Zero, dU_q(:,:,iCR), nDOF )

        CALL DGEMM &
               ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Up(:,:), nDOF_X1, &
                 U_P(:,:,iCR), nDOF, Zero, U_L(:,:,iCR), nDOF_X1 )

        CALL DGEMM &
               ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Dn(:,:), nDOF_X1, &
                 U_K(:,:,iCR), nDOF, Zero, U_R(:,:,iCR), nDOF_X1 )

      END DO

      ! Timing
      f_dgemm_t = f_dgemm_t + MPI_Wtime() - start_t

      DO iF = 1,nF

        ! --- Left State Primitive, etc. ---

        start_t = MPI_Wtime()

        CALL ComputePrimitive &
               ( U_L(:,iF,iCR_N ), U_L(:,iF,iCR_G1), &
                 U_L(:,iF,iCR_G2), U_L(:,iF,iCR_G3), &
                 P_L(:,   iPR_D ), P_L(:,   iPR_I1), &
                 P_L(:,   iPR_I2), P_L(:,   iPR_I3), &
                 Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        FF_L = FluxFactor &
                 ( P_L(:,iPR_D ), P_L(:,iPR_I1), &
                   P_L(:,iPR_I2), P_L(:,iPR_I3), &
                   Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        EF_L = EddingtonFactor( P_L(:,iPR_D), FF_L )

        DO iNode = 1, nDOF_X1

          Flux_L(iNode,iF,1:nCR) &
            = Flux_X1 &
                ( P_L(iNode,iPR_D ), P_L(iNode,iPR_I1), &
                  P_L(iNode,iPR_I2), P_L(iNode,iPR_I3), &
                  FF_L(iNode), EF_L(iNode), One, One, One )

        END DO

        ! Timing
        l_flux_t = l_flux_t + MPI_Wtime() - start_t

        start_t = MPI_Wtime()

        ! --- Right State Primitive, etc. ---

        CALL ComputePrimitive &
               ( U_R(:,iF,iCR_N ), U_R(:,iF,iCR_G1), &
                 U_R(:,iF,iCR_G2), U_R(:,iF,iCR_G3), &
                 P_R(:,   iPR_D ), P_R(:,   iPR_I1), &
                 P_R(:,   iPR_I2), P_R(:,   iPR_I3), &
                 Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        FF_R = FluxFactor &
                 ( P_R(:,iPR_D ), P_R(:,iPR_I1), &
                   P_R(:,iPR_I2), P_R(:,iPR_I3), &
                   Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        EF_R = EddingtonFactor( P_R(:,iPR_D), FF_R )

        DO iNode = 1, nDOF_X1

          Flux_R(iNode,iF,1:nCR) &
            = Flux_X1 &
                ( P_R(iNode,iPR_D ), P_R(iNode,iPR_I1), &
                  P_R(iNode,iPR_I2), P_R(iNode,iPR_I3), &
                  FF_R(iNode), EF_R(iNode), One, One, One )

        END DO

        ! Timing
        r_flux_t = r_flux_t + MPI_Wtime() - start_t

        start_t = MPI_Wtime()

        ! --- Numerical Flux ---

        DO iCR = 1, nCR

          NumericalFlux(:,iF,iCR) &
            = NumericalFlux_LLF &
                ( U_L   (:,iF,iCR), U_R   (:,iF,iCR), &
                  Flux_L(:,iF,iCR), Flux_R(:,iF,iCR), &
                  Ones_RL(:) ) * NumericalFlux(:,iF,iCR)

        END DO

        ! Timing
        s_num_flux_t = s_num_flux_t + MPI_Wtime() - start_t

      END DO ! --- iF

      ! Timing
      start_t = MPI_Wtime()

      DO iCR = 1, nCR

        CALL DGEMM &
               ( 'T', 'N', nDOF, nF, nDOF_X1, One, L_X1_Dn, nDOF_X1, &
                 NumericalFlux(:,:,iCR), nDOF_X1, Zero, dU_K(:,:,iCR), nDOF )

        CALL DGEMM &
               ( 'T', 'N', nDOF, nF, nDOF_X1, One, L_X1_Up, nDOF_X1, &
                 NumericalFlux(:,:,iCR), nDOF_X1, Zero, dU_P(:,:,iCR), nDOF )

      END DO

      ! Timing
      s_dgemm_t = s_dgemm_t + MPI_Wtime() - start_t

      start_t = MPI_Wtime()

      iK = 0
      iF = 0

      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                ! --- Volume Term ---

                iK = iK + 1

                dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                      + dU_q(:,iK,1:nCR)

              END IF

              IF = iF + 1

              ! --- Flux Terms ---

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                ! --- Contribution to this Element ---

                dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                      + dU_K(:,iF,1:nCR)

              END IF

              IF( iZ2 > iZ_B0(2) )THEN

                ! --- Contribution to this Element ---

                dU(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS) &
                      - dU_P(:,iF,1:nCR)

              END IF

            END DO
          END DO
        END DO
      END DO

      ! Timing
      update_t = update_t + MPI_Wtime() - start_t

      !PRINT *, "Printing calculation times:"
      PRINT *, "first dgemm calls: ", f_dgemm_t, "left flux: ", l_flux_t
      PRINT *, "right flux: ", r_flux_t
      PRINT *, "second num flux update: ", s_num_flux_t
      PRINT *, "second dgemm calls: ", s_dgemm_t
      PRINT *, "update of result: ", update_t
      PRINT *, ""

      wTime = MPI_WTIME( ) - wTime

      PRINT *,  "Comp Time: ", wTime


      PRINT *, "Total time: ", MPI_Wtime() - total_t
      PRINT *, ""
      PRINT *, ""
    
    END DO ! --- nSpecies

    DEALLOCATE( Flux_X1_q, dU_q, U_P, U_K )

    DEALLOCATE( U_L, U_R, Flux_L, Flux_R, NumericalFlux, dU_P, dU_K )

  END SUBROUTINE ComputeIncrement_Divergence_X1_GPU

  SUBROUTINE ComputeIncrement_Divergence_X1_GPU_TEST &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: nK, iK, nF, iF, nZ(4), iNode
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iCR
    
    INTEGER  :: size_range(1:4)
    INTEGER  :: Ones_size(1:4)

    REAL(DP) :: wTime, start_t, end_t
    REAL(DP) :: primitive_t, ff_t, ef_t, flux_x1_q_t
    REAL(DP) :: copy_prev_t, num_flux_t, f_dgemm_t
    REAL(DP) :: l_flux_t, r_flux_t, s_num_flux_t
    REAL(DP) :: s_dgemm_t, update_t, total_t
    
    REAL(DP) :: Ones_q(nDOF), Ones_RL(nDOF_X1), dZ(4)
    REAL(DP) :: Tau(1:nDOF), Tau_X1(1:nDOF_X1)
    REAL(DP) :: P_q(1:nDOF,1:nPR), FF_q(1:nDOF), EF_q(1:nDOF)
    REAL(DP), ALLOCATABLE :: Flux_X1_q(:,:,:), dU_q(:,:,:)
    REAL(DP), ALLOCATABLE :: U_P(:,:,:), U_K(:,:,:)
    REAL(DP), ALLOCATABLE :: U_L(:,:,:), Flux_L(:,:,:)
    REAL(DP), ALLOCATABLE :: U_R(:,:,:), Flux_R(:,:,:)
    REAL(DP), ALLOCATABLE :: NumericalFlux(:,:,:), dU_P(:,:,:), dU_K(:,:,:)
    REAL(DP) :: P_L(1:nDOF_X1,1:nPR), FF_L(1:nDOF_X1), EF_L(1:nDOF_X1)
    REAL(DP) :: P_R(1:nDOF_X1,1:nPR), FF_R(1:nDOF_X1), EF_R(1:nDOF_X1)

    nZ = iZ_E0 - iZ_B0 + 1
    nK = PRODUCT( nZ )
    nF = PRODUCT( nZ + [0,1,0,0] )

    Ones_q = One

    Ones_RL = One

    Ones_size = One

    ALLOCATE( Flux_X1_q(1:nDOF,1:nK,1:nCR) )
    ALLOCATE( dU_q     (1:nDOF,1:nK,1:nCR) )

    ALLOCATE( U_P(1:nDOF,1:nF,1:nCR) )
    ALLOCATE( U_K(1:nDOF,1:nF,1:nCR) )

    ALLOCATE( U_L(1:nDOF_X1,1:nF,1:nCR) )
    ALLOCATE( U_R(1:nDOF_X1,1:nF,1:nCR) )

    ALLOCATE( Flux_L(1:nDOF_X1,1:nF,1:nPR) )
    ALLOCATE( Flux_R(1:nDOF_X1,1:nF,1:nPR) )

    ALLOCATE( NumericalFlux(1:nDOF_X1,1:nF,1:nCR) )

    ALLOCATE( dU_P(1:nDOF,1:nF,1:nCR) )
    ALLOCATE( dU_K(1:nDOF,1:nF,1:nCR) )

    ! Init timing summation variables
    primitive_t  = 0.0d0
    ff_t         = 0.0d0
    ef_t         = 0.0d0
    flux_x1_q_t  = 0.0d0
    copy_prev_t  = 0.0d0
    num_flux_t   = 0.0d0
    f_dgemm_t    = 0.0d0
    l_flux_t     = 0.0d0
    r_flux_t     = 0.0d0
    s_num_flux_t = 0.0d0
    s_dgemm_t    = 0.0d0
    update_t     = 0.0d0

    total_t = MPI_Wtime()

    DO iS = 1, nSpecies

      iK = 0
      iF = 0

      wTime = MPI_WTIME( )

      ! Added to accomodate ranges and make iK dependent on iZ
      size_range = iZ_E0 - iZ_B0 + Ones_size

      !$OMP TARGET MAP(alloc:) MAP(tofrom: )

      !$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) PRIVATE(iZ4, iZ3, iZ2, iZ1, dZ, &
      !$OMP& FF_q, EF_q, Tau, Tau_X1, P_q) FIRSTPRIVATE(iK, iF)
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dZ(1) = MeshE    % Width(iZ1)
              dZ(2) = MeshX(1) % Width(iZ2)
              dZ(3) = MeshX(2) % Width(iZ3)
              dZ(4) = MeshX(3) % Width(iZ4)

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, Ones_q(1:nDOFX), nDOFX )

              Tau_X1(1:nDOF_X1) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, Ones_q(1:nDOFX_X1), nDOFX_X1 )

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                !iK = iK + 1
                iK = (iZ1 - iZ_B0(1) + 1) + (iZ2 - iZ_B0(2)) * size_range(1) +&
                        (iZ3 - iZ_B0(3)) * size_range(2) * size_range(1) + &
                        (iZ4 - iZ_B0(4)) * size_range(3) * size_range(2) * size_range(1)

                !PRINT *, "iK: ", iK
                !PRINT *, "iZ1: ", iZ1, "iZ_B0(1): ", iZ_B0(1), "iZ_E0(1): ", iZ_E0(1)
                !PRINT *, "iZ2: ", iZ2, "iZ_B0(2): ", iZ_B0(2), "iZ_E0(2): ", iZ_E0(2)
                !PRINT *, "iZ3: ", iZ3, "iZ_B0(3): ", iZ_B0(3), "iZ_E0(3): ", iZ_E0(3)
                !PRINT *, "iZ4: ", iZ4, "iZ_B0(4): ", iZ_B0(4), "iZ_E0(4): ", iZ_E0(4)
                !PRINT *, "Sizes: ", size_range

                ! Timing
                !start_t = MPI_Wtime()


                CALL ComputePrimitive &
                       ( U(:,iZ1,iZ2,iZ3,iZ4,iCR_N, iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                         U(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                         P_q(:,iPR_D ), P_q(:,iPR_I1), &
                         P_q(:,iPR_I2), P_q(:,iPR_I3), &
                         Ones_q(:), Ones_q(:), Ones_q(:) )

                ! Timing
                !primitive_t = primitive_t + MPI_Wtime() - start_t

                !start_t = MPI_Wtime()


                FF_q = FluxFactor &
                         ( P_q(:,iPR_D ), P_q(:,iPR_I1), &
                           P_q(:,iPR_I2), P_q(:,iPR_I3), &
                           Ones_q(:), Ones_q(:), Ones_q(:) )

                ! Timing
                !ff_t = ff_t + MPI_Wtime() - start_t

                !start_t = MPI_Wtime()

                
                EF_q = EddingtonFactor( P_q(:,iPR_D), FF_q )

                ! Timing
                !ef_t = ef_t + MPI_Wtime() - start_t

                !start_t = MPI_Wtime()

                DO iNode = 1, nDOF

                  Flux_X1_q(iNode,iK,1:nCR) &
                    = Flux_X1 &
                        ( P_q(iNode,iPR_D ), P_q(iNode,iPR_I1), &
                          P_q(iNode,iPR_I2), P_q(iNode,iPR_I3), &
                          FF_q(iNode), EF_q(iNode), One, One, One )

                END DO

                DO iCR = 1, nCR

                  Flux_X1_q(:,iK,iCR) &
                    = dZ(1) * dZ(3) * dZ(4) * Weights_q(:) &
                        * Tau(:) * Flux_X1_q(:,iK,iCR)

                END DO

                ! Timing
                !flux_x1_q_t = flux_x1_q_t + MPI_Wtime() - start_t

              END IF

              ! Timing
              !start_t = MPI_Wtime()

              !iF = iF + 1
              iF = (iZ1 - iZ_B0(1) + 1) + (iZ2 - iZ_B0(2)) * size_range(1) +&
                        (iZ3 - iZ_B0(3)) * (size_range(2) + 1) * size_range(1) + &
                        (iZ4 - iZ_B0(4)) * size_range(3) * (size_range(2) + 1) * size_range(1)

              !PRINT *, "iF: ", iF

              U_P(:,iF,1:nCR) = U(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS)
              U_K(:,iF,1:nCR) = U(:,iZ1,iZ2,  iZ3,iZ4,1:nCR,iS)

              ! Timing
              !copy_prev_t = copy_prev_t + MPI_Wtime() - start_t

              !start_t = MPI_Wtime()

              DO iCR = 1, nCR

                NumericalFlux(:,iF,iCR) &
                  = dZ(1) * dZ(3) * dZ(4) * Weights_X1(:) * Tau_X1(:)

              END DO

              ! Timing
              !num_flux_t = num_flux_t + MPI_Wtime() - start_t

            END DO
          END DO
        END DO
      END DO
      !$OMP END TARGET

      wTime = MPI_WTIME( ) - wTime

      PRINT *, "Prep Time: ", wTime

      !PRINT *, "Printing prep time proportions:"
      PRINT *, "primitive: ", primitive_t, "ff: ", ff_t
      PRINT *, "ef: ", ef_t
      PRINT *, "flux_x1_q: ", flux_x1_q_t, "copy left/right: ", copy_prev_t
      PRINT *, "numerical flux: ", num_flux_t
      PRINT *, ""

      wTime = MPI_WTIME( )

      start_t = MPI_Wtime()

      DO iCR = 1, nCR

        CALL DGEMM &
               ( 'T', 'N', nDOF, nK, nDOF, One, dLdX1_q(:,:), nDOF, &
                 Flux_X1_q(:,:,iCR), nDOF, Zero, dU_q(:,:,iCR), nDOF )

        CALL DGEMM &
               ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Up(:,:), nDOF_X1, &
                 U_P(:,:,iCR), nDOF, Zero, U_L(:,:,iCR), nDOF_X1 )

        CALL DGEMM &
               ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Dn(:,:), nDOF_X1, &
                 U_K(:,:,iCR), nDOF, Zero, U_R(:,:,iCR), nDOF_X1 )

      END DO

      ! Timing
      f_dgemm_t = f_dgemm_t + MPI_Wtime() - start_t

      DO iF = 1,nF

        ! --- Left State Primitive, etc. ---

        start_t = MPI_Wtime()

        CALL ComputePrimitive &
               ( U_L(:,iF,iCR_N ), U_L(:,iF,iCR_G1), &
                 U_L(:,iF,iCR_G2), U_L(:,iF,iCR_G3), &
                 P_L(:,   iPR_D ), P_L(:,   iPR_I1), &
                 P_L(:,   iPR_I2), P_L(:,   iPR_I3), &
                 Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        FF_L = FluxFactor &
                 ( P_L(:,iPR_D ), P_L(:,iPR_I1), &
                   P_L(:,iPR_I2), P_L(:,iPR_I3), &
                   Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        EF_L = EddingtonFactor( P_L(:,iPR_D), FF_L )

        DO iNode = 1, nDOF_X1

          Flux_L(iNode,iF,1:nCR) &
            = Flux_X1 &
                ( P_L(iNode,iPR_D ), P_L(iNode,iPR_I1), &
                  P_L(iNode,iPR_I2), P_L(iNode,iPR_I3), &
                  FF_L(iNode), EF_L(iNode), One, One, One )

        END DO

        ! Timing
        l_flux_t = l_flux_t + MPI_Wtime() - start_t

        start_t = MPI_Wtime()

        ! --- Right State Primitive, etc. ---

        CALL ComputePrimitive &
               ( U_R(:,iF,iCR_N ), U_R(:,iF,iCR_G1), &
                 U_R(:,iF,iCR_G2), U_R(:,iF,iCR_G3), &
                 P_R(:,   iPR_D ), P_R(:,   iPR_I1), &
                 P_R(:,   iPR_I2), P_R(:,   iPR_I3), &
                 Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        FF_R = FluxFactor &
                 ( P_R(:,iPR_D ), P_R(:,iPR_I1), &
                   P_R(:,iPR_I2), P_R(:,iPR_I3), &
                   Ones_RL(:), Ones_RL(:), Ones_RL(:) )

        EF_R = EddingtonFactor( P_R(:,iPR_D), FF_R )

        DO iNode = 1, nDOF_X1

          Flux_R(iNode,iF,1:nCR) &
            = Flux_X1 &
                ( P_R(iNode,iPR_D ), P_R(iNode,iPR_I1), &
                  P_R(iNode,iPR_I2), P_R(iNode,iPR_I3), &
                  FF_R(iNode), EF_R(iNode), One, One, One )

        END DO

        ! Timing
        r_flux_t = r_flux_t + MPI_Wtime() - start_t

        start_t = MPI_Wtime()

        ! --- Numerical Flux ---

        DO iCR = 1, nCR

          NumericalFlux(:,iF,iCR) &
            = NumericalFlux_LLF &
                ( U_L   (:,iF,iCR), U_R   (:,iF,iCR), &
                  Flux_L(:,iF,iCR), Flux_R(:,iF,iCR), &
                  Ones_RL(:) ) * NumericalFlux(:,iF,iCR)

        END DO

        ! Timing
        s_num_flux_t = s_num_flux_t + MPI_Wtime() - start_t

      END DO ! --- iF

      ! Timing
      start_t = MPI_Wtime()

      DO iCR = 1, nCR

        CALL DGEMM &
               ( 'T', 'N', nDOF, nF, nDOF_X1, One, L_X1_Dn, nDOF_X1, &
                 NumericalFlux(:,:,iCR), nDOF_X1, Zero, dU_K(:,:,iCR), nDOF )

        CALL DGEMM &
               ( 'T', 'N', nDOF, nF, nDOF_X1, One, L_X1_Up, nDOF_X1, &
                 NumericalFlux(:,:,iCR), nDOF_X1, Zero, dU_P(:,:,iCR), nDOF )

      END DO

      ! Timing
      s_dgemm_t = s_dgemm_t + MPI_Wtime() - start_t

      start_t = MPI_Wtime()

      iK = 0
      iF = 0

      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                ! --- Volume Term ---

                iK = iK + 1

                dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                      + dU_q(:,iK,1:nCR)

              END IF

              IF = iF + 1

              ! --- Flux Terms ---

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                ! --- Contribution to this Element ---

                dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                      + dU_K(:,iF,1:nCR)

              END IF

              IF( iZ2 > iZ_B0(2) )THEN

                ! --- Contribution to this Element ---

                dU(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS) &
                  = dU(:,iZ1,iZ2-1,iZ3,iZ4,1:nCR,iS) &
                      - dU_P(:,iF,1:nCR)

              END IF

            END DO
          END DO
        END DO
      END DO

      ! Timing
      update_t = update_t + MPI_Wtime() - start_t

      !PRINT *, "Printing calculation times:"
      PRINT *, "first dgemm calls: ", f_dgemm_t, "left flux: ", l_flux_t
      PRINT *, "right flux: ", r_flux_t
      PRINT *, "second num flux update: ", s_num_flux_t
      PRINT *, "second dgemm calls: ", s_dgemm_t
      PRINT *, "update of result: ", update_t
      PRINT *, ""

      wTime = MPI_WTIME( ) - wTime

      PRINT *,  "Comp Time: ", wTime


      PRINT *, "Total time: ", MPI_Wtime() - total_t
      PRINT *, ""
      PRINT *, ""
    
    END DO ! --- nSpecies

    DEALLOCATE( Flux_X1_q, dU_q, U_P, U_K )

    DEALLOCATE( U_L, U_R, Flux_L, Flux_R, NumericalFlux, dU_P, dU_K )

  END SUBROUTINE ComputeIncrement_Divergence_X1_GPU_TEST



  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNode
    INTEGER  :: iGF, iCR
    REAL(DP) :: dZ(4)
    REAL(DP) :: FF, EF
    REAL(DP) :: GX_P(nDOFX,   nGF)
    REAL(DP) :: GX_K(nDOFX,   nGF)
    REAL(DP) :: GX_F(nDOFX_X1,nGF)
    REAL(DP) :: G_K(nDOF,nGF)
    REAL(DP) :: G_F(nDOF_X1,nGF)
    REAL(DP), DIMENSION(nDOF_X1)     :: absLambda_L
    REAL(DP), DIMENSION(nDOF_X1)     :: absLambda_R
    REAL(DP), DIMENSION(nDOF_X1)     :: alpha
    REAL(DP), DIMENSION(nDOF)        :: Tau
    REAL(DP), DIMENSION(nDOF_X1)     :: Tau_X1
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: uCR_L, uCR_R
    REAL(DP), DIMENSION(nDOF_X1,nPR) :: uPR_L, uPR_R
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: Flux_X1_L
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: Flux_X1_R
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: NumericalFlux
    REAL(DP), DIMENSION(nDOF   ,nCR) :: uCR_P, uCR_K
    REAL(DP), DIMENSION(nDOF   ,nPR) :: uPR_K
    REAL(DP), DIMENSION(nDOF   ,nCR) :: Flux_X1_q

    Timer_VOL   = 0.0_DP
    Timer_AD1   = 0.0_DP
    Timer_SUR   = 0.0_DP
    Timer_INT   = 0.0_DP
    Timer_INT_G = 0.0_DP
    Timer_LFT   = 0.0_DP
    Timer_RGT   = 0.0_DP
    Timer_FLX   = 0.0_DP
    Timer_AD2   = 0.0_DP
    Timer_AD3   = 0.0_DP

    CALL Timer_Start( Timer_RHS )

    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)

        dZ(4) = MeshX(3) % Width(iZ4)

        DO iZ3 = iZ_B0(3), iZ_E0(3)

          dZ(3) = MeshX(2) % Width(iZ3)

          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1

            ! --- Geometry Fields in Element Nodes ---

            DO iGF = 1, nGF

              GX_P(:,iGF) = GX(:,iZ2-1,iZ3,iZ4,iGF) ! --- Previous Element
              GX_K(:,iGF) = GX(:,iZ2,  iZ3,iZ4,iGF) ! --- This     Element

              G_K(1:nDOF,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_K(1:nDOFX,iGF), nDOFX )

            END DO

            ! --- Interpolate Geometry Fields on Shared Face ---

            CALL Timer_Start( dT_INT_G )

            ! --- Face States (Average of Left and Right States) ---

            ! --- Scale Factors ---

            DO iGF = iGF_h_1, iGF_h_3

              CALL DGEMV &
                     ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                       GX_P(:,iGF), 1, Zero, GX_F(:,iGF), 1 )
              CALL DGEMV &
                     ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                       GX_K(:,iGF), 1, Half, GX_F(:,iGF), 1 )

              GX_F(1:nDOFX_X1,iGF) &
                = MAX( GX_F(1:nDOFX_X1,iGF), SqrtTiny )

              G_F(1:nDOF_X1,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_F(1:nDOFX_X1,iGF), nDOFX_X1 )

            END DO

            CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

            ! --- Lapse Function ---

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     GX_P(:,iGF_Alpha), 1, Zero, GX_F(:,iGF_Alpha), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     GX_K(:,iGF_Alpha), 1, Half, GX_F(:,iGF_Alpha), 1 )

            GX_F(1:nDOFX_X1,iGF_Alpha) &
              = MAX( GX_F(1:nDOFX_X1,iGF_Alpha), SqrtTiny )

            G_F(1:nDOF_X1,iGF_Alpha) &
              = OuterProduct1D3D &
                  ( Ones(1:nDOFE), nDOFE, &
                    GX_F(1:nDOFX_X1,iGF_Alpha), nDOFX_X1 )

            CALL Timer_Stop( dT_INT_G )

            CALL Timer_Add( Timer_INT_G, dT_INT_G )

            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dZ(1) = MeshE % Width(iZ1)

              ! --- Volume Jacobian in Energy-Position Element ---

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, G_K(:,iGF_SqrtGm), nDOFX )

              Tau_X1(1:nDOF_X1) &
                = OuterProduct1D3D &
                    ( GE(:,iZ1,iGE_Ep2), nDOFE, G_F(:,iGF_SqrtGm), nDOFX_X1 )

              DO iCR = 1, nCR

                uCR_P(:,iCR) = U(:,iZ1,iZ2-1,iZ3,iZ4,iCR,iS)
                uCR_K(:,iCR) = U(:,iZ1,iZ2,  iZ3,iZ4,iCR,iS)

              END DO

              !--------------------
              ! --- Volume Term ---
              !--------------------

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                CALL Timer_Start( dT_VOL )

                CALL ComputePrimitive &
                       ( uCR_K(:,iCR_N ), uCR_K(:,iCR_G1), &
                         uCR_K(:,iCR_G2), uCR_K(:,iCR_G3), &
                         uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                DO iNode = 1, nDOF

                  FF = FluxFactor &
                         ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                           uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                           G_K(iNode,iGF_Gm_dd_11), &
                           G_K(iNode,iGF_Gm_dd_22), &
                           G_K(iNode,iGF_Gm_dd_33) )


                  EF = EddingtonFactor( uPR_K(iNode,iPR_D), FF )

                  Flux_X1_q(iNode,1:nCR) &
                    = Flux_X1 &
                        ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                          uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                          FF, EF, &
                          G_K(iNode,iGF_Gm_dd_11), &
                          G_K(iNode,iGF_Gm_dd_22), &
                          G_K(iNode,iGF_Gm_dd_33) )

                END DO

                CALL Timer_Start( dT_AD1 )

                DO iCR = 1, nCR

                  Flux_X1_q(:,iCR) &
                    = dZ(1) * dZ(3) * dZ(4) * Weights_q(:) &
                        * G_K(:,iGF_Alpha) * Tau(:) * Flux_X1_q(:,iCR)

                  CALL DGEMV &
                         ( 'T', nDOF, nDOF, One, dLdX1_q, nDOF, &
                           Flux_X1_q(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS), 1 )

                END DO

                CALL Timer_Stop( dT_AD1 )

                CALL Timer_Add( Timer_AD1, dT_AD1 )

                CALL Timer_Stop( dT_VOL )

                CALL Timer_Add( Timer_VOL, dT_VOL )

              END IF

              !---------------------
              ! --- Surface Term ---
              !---------------------

              CALL Timer_Start( dT_SUR )

              CALL Timer_Start( dT_INT )

              ! --- Interpolate Radiation Fields ---

              DO iCR = 1, nCR

                ! --- Interpolate Left State ---

                CALL DGEMV &
                       ( 'N', nDOF_X1, nDOF, One, L_X1_Up, nDOF_X1, &
                         uCR_P(:,iCR), 1, Zero, uCR_L(:,iCR), 1 )

                ! --- Interpolate Right State ---

                CALL DGEMV &
                       ( 'N', nDOF_X1, nDOF, One, L_X1_Dn, nDOF_X1, &
                         uCR_K(:,iCR), 1, Zero, uCR_R(:,iCR), 1 )

              END DO

              CALL Timer_Stop( dT_INT )

              CALL Timer_Add( Timer_INT, dT_INT )

              ! --- Left State Primitive, etc. ---

              CALL ComputePrimitive &
                     ( uCR_L(:,iCR_N ), uCR_L(:,iCR_G1), &
                       uCR_L(:,iCR_G2), uCR_L(:,iCR_G3), &
                       uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                       uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              CALL Timer_Start( dT_LFT )

              DO iNode = 1, nDOF_X1

                FF = FluxFactor &
                       ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                         uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11), &
                         G_F(iNode,iGF_Gm_dd_22), &
                         G_F(iNode,iGF_Gm_dd_33) )

                EF = EddingtonFactor( uPR_L(iNode,iPR_D), FF )

                Flux_X1_L(iNode,1:nCR) &
                  = Flux_X1 &
                      ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                        uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                        FF, EF, &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

                absLambda_L(iNode) = 1.0_DP

              END DO

              CALL Timer_Stop( dT_LFT )

              CALL Timer_Add( Timer_LFT, dT_LFT )

              ! --- Right State Primitive, etc. ---

              CALL ComputePrimitive &
                     ( uCR_R(:,iCR_N ), uCR_R(:,iCR_G1), &
                       uCR_R(:,iCR_G2), uCR_R(:,iCR_G3), &
                       uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                       uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              CALL Timer_Start( dT_RGT )

              DO iNode = 1, nDOF_X1

                FF = FluxFactor &
                       ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                         uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11), &
                         G_F(iNode,iGF_Gm_dd_22), &
                         G_F(iNode,iGF_Gm_dd_33) )

                EF = EddingtonFactor( uPR_R(iNode,iPR_D), FF )

                Flux_X1_R(iNode,1:nCR) &
                  = Flux_X1 &
                      ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                        uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                        FF, EF, &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

                absLambda_R(iNode) = 1.0_DP

              END DO

              CALL Timer_Stop( dT_RGT )

              CALL Timer_Add( Timer_RGT, dT_RGT )

              ! --- Numerical Flux ---

              CALL Timer_Start( dT_FLX )

              alpha = MAX( absLambda_L, absLambda_R )

              DO iCR = 1, nCR

                NumericalFlux(:,iCR) &
                  = NumericalFlux_LLF &
                      ( uCR_L    (:,iCR), &
                        uCR_R    (:,iCR), &
                        Flux_X1_L(:,iCR), &
                        Flux_X1_R(:,iCR), alpha(:) )

                NumericalFlux(:,iCR) &
                  = dZ(1) * dZ(3) * dZ(4) * Weights_X1(:) &
                      * G_F(:,iGF_Alpha) * Tau_X1(:) * NumericalFlux(:,iCR)

              END DO

              CALL Timer_Stop( dT_FLX )

              CALL Timer_Add( Timer_FLX, dT_FLX )

              ! --- Contribution to this Element ---

              CALL Timer_Start( dT_AD2 )

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X1, nDOF, + One, L_X1_Dn, &
                           nDOF_X1, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2  ,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              CALL Timer_Stop( dT_AD2 )

              CALL Timer_Add( Timer_AD2, dT_AD2 )

              ! --- Contribution to Previous Element ---

              CALL Timer_Start( dT_AD3 )

              IF( iZ2 > iZ_B0(2) )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X1, nDOF, - One, L_X1_Up, &
                           nDOF_X1, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2-1,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              CALL Timer_Stop( dT_AD3 )

              CALL Timer_Add( Timer_AD3, dT_AD3 )

              CALL Timer_Stop( dT_SUR )

              CALL Timer_Add( Timer_SUR, dT_SUR )

            END DO ! iZ1
          END DO ! iZ2
        END DO ! iZ3
      END DO ! iZ4
    END DO ! iS

    CALL Timer_Stop( Timer_RHS )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeRHS: ', Timer_RHS
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Volume Term: ', Timer_VOL
      WRITE(*,*)
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Add 1: ', Timer_AD1
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Surface Term: ', Timer_SUR
      WRITE(*,*)
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Interp: ', Timer_INT
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Int (G): ', Timer_INT_G
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Left: ', Timer_LFT
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Right: ', Timer_RGT
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Flux: ', Timer_FLX
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Add 2: ', Timer_AD2
      WRITE(*,'(A6,A10,ES10.4E2)') &
        '', 'Add 3: ', Timer_AD3
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Inverse Mass: ', Timer_INV
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Sum: ', Timer_VOL + Timer_SUR + Timer_INV

    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    IF( .NOT. DisplayTimers ) RETURN

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE dgDiscretizationModule
