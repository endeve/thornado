MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Pi, TwoPi, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Phi_N, &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_BoundaryConditionsModule, ONLY: &
    Euler_ApplyBoundaryConditions
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive, &
    AlphaPlus, &
    AlphaMinus, &
    AlphaMiddle, &
    Flux_X1, &
    Flux_X2, &
    StressTensor_Diagonal, &
    NumericalFlux_HLL, &
    NumericalFlux_X1_HLLC, &
    NumericalFlux_X2_HLLC, &
    NumericalFlux_X3_HLLC
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_ComputeIncrement_DG_Explicit

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer_RHS
  REAL(DP) :: Timer_RHS_1, dT_RHS_1
  REAL(DP) :: Timer_RHS_2, dT_RHS_2
  REAL(DP) :: Timer_RHS_3, dT_RHS_3
  REAL(DP) :: Timer_INT_F, dT_INT_F
  REAL(DP) :: Timer_INT_G, dT_INT_G
  REAL(DP) :: Timer_FLX_N, dT_FLX_N
  REAL(DP) :: Timer_Div_X2
  REAL(DP) :: Timer_Geo

CONTAINS


  SUBROUTINE Euler_ComputeIncrement_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)          :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    LOGICAL  :: SuppressBC

    dU = Zero

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL Euler_ApplyBoundaryConditions &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    ! --- Multiply Inverse Mass Matrix ---

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)

            dX1 = MeshX(1) % Width(iX1)
            dX2 = MeshX(2) % Width(iX2)
            dX3 = MeshX(3) % Width(iX3)

            dU(:,iX1,iX2,iX3,iCF) &
              = dU(:,iX1,iX2,iX3,iCF) &
                  / ( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                        * dX1 * dX2 * dX3 )

          END DO
        END DO
      END DO
    END DO

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

  END SUBROUTINE Euler_ComputeIncrement_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X1
    REAL(DP) :: dX2, dX3
    REAL(DP) :: AlphaPls, AlphaMns, AlphaMdl
    REAL(DP), DIMENSION(nDOFX_X1)     :: P_L, Cs_L
    REAL(DP), DIMENSION(nDOFX_X1)     :: P_R, Cs_R
    REAL(DP), DIMENSION(nDOFX)        :: P_K
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(nDOFX_X1,nGF) :: G_F
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_L
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: NumericalFlux
    REAL(DP), DIMENSION(nDOFX,   nCF) :: uCF_P, uCF_K
    REAL(DP), DIMENSION(nDOFX,   nGF) :: G_P, G_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: Flux_X1_q
    REAL(DP), DIMENSION(nDOFX_X1,nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(nDOFX,   nPF) :: uPF_K

    IF( iX_E0(1) == iX_B0(1) ) RETURN

    Timer_RHS_1 = 0.0_DP
    Timer_RHS_2 = 0.0_DP
    Timer_RHS_3 = 0.0_DP
    Timer_INT_F = 0.0_DP
    Timer_INT_G = 0.0_DP
    Timer_FLX_N = 0.0_DP

    CALL Timer_Start( Timer_RHS )

    !$OMP PARALLEL DO PRIVATE &
    !$OMP& ( iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X1, dX2, dX3, &
    !$OMP&   uCF_P, uCF_K, uCF_L, uCF_R, uPF_K, uPF_L, uPF_R, P_K, &
    !$OMP&   P_L, P_R, Cs_L, Cs_R, G_P, G_K, G_F, Flux_X1_q, &
    !$OMP&   Flux_X1_L, Flux_X1_R, NumericalFlux )
    DO iX3 = iX_B0(3), iX_E0(3)

      dX3 = MeshX(3) % Width(iX3)

      DO iX2 = iX_B0(2), iX_E0(2)

        dX2 = MeshX(2) % Width(iX2)

        DO iX1 = iX_B0(1), iX_E0(1) + 1

          DO iCF = 1, nCF

            uCF_P(:,iCF) = U(:,iX1-1,iX2,iX3,iCF)
            uCF_K(:,iCF) = U(:,iX1,  iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            G_P(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
            G_K(:,iGF) = G(:,iX1,  iX2,iX3,iGF)

          END DO

          !--------------------
          ! --- Volume Term ---
          !--------------------

          IF( iX1 < iX_E0(1) + 1 )THEN

            CALL ComputePrimitive &
                   ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                     uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                     uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                     uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                     G_K(:,iGF_Gm_dd_11), &
                     G_K(:,iGF_Gm_dd_22), &
                     G_K(:,iGF_Gm_dd_33) )

            CALL ComputePressureFromPrimitive &
                   ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

            DO iNodeX = 1, nDOFX

              Flux_X1_q(iNodeX,1:nCF) &
                = Flux_X1 &
                    ( uPF_K(iNodeX,iPF_D ), &
                      uPF_K(iNodeX,iPF_V1), &
                      uPF_K(iNodeX,iPF_V2), &
                      uPF_K(iNodeX,iPF_V3), &
                      uPF_K(iNodeX,iPF_E ), &
                      uPF_K(iNodeX,iPF_Ne), &
                      P_K  (iNodeX), &
                      G_K(iNodeX,iGF_Gm_dd_11), &
                      G_K(iNodeX,iGF_Gm_dd_22), &
                      G_K(iNodeX,iGF_Gm_dd_33) )

            END DO

            CALL Timer_Start( dT_RHS_1 )

            DO iCF = 1, nCF

              Flux_X1_q(:,iCF) &
                = dX2 * dX3 * WeightsX_q(:) * G_K(:,iGF_SqrtGm) &
                    * Flux_X1_q(:,iCF)

              CALL DGEMV &
                     ( 'T', nDOFX, nDOFX, One, dLXdX1_q, nDOFX, &
                       Flux_X1_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

            CALL Timer_Stop( dT_RHS_1 )

            CALL Timer_Add( Timer_RHS_1, dT_RHS_1 )

          END IF

          !---------------------
          ! --- Surface Term ---
          !---------------------

          ! --- Interpolate Fluid Fields ---

          CALL Timer_Start( dT_INT_F )

          DO iCF = 1, nCF

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                     uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                     uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

          END DO

          CALL Timer_Stop( dT_INT_F )

          CALL Timer_Add( Timer_INT_F, dT_INT_F )

          ! --- Interpolate Geometry Fields ---

          CALL Timer_Start( dT_INT_G )

          ! --- Face States (Average of Left and Right States) ---

          DO iGF = iGF_h_1, iGF_h_3

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

            G_F(1:nDOFX_X1,iGF) &
              = MAX( G_F(1:nDOFX_X1,iGF), SqrtTiny )

          END DO

          CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

          CALL Timer_Stop( dT_INT_G )

          CALL Timer_Add( Timer_INT_G, dT_INT_G )

          ! --- Left State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
                   uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
                   uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
                   uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L(:) )

          DO iNodeX_X1 = 1, nDOFX_X1

            Flux_X1_L(iNodeX_X1,1:nCF) &
              = Flux_X1 &
                  ( uPF_L(iNodeX_X1,iPF_D ), &
                    uPF_L(iNodeX_X1,iPF_V1), &
                    uPF_L(iNodeX_X1,iPF_V2), &
                    uPF_L(iNodeX_X1,iPF_V3), &
                    uPF_L(iNodeX_X1,iPF_E ), &
                    uPF_L(iNodeX_X1,iPF_Ne), &
                    P_L  (iNodeX_X1), &
                    G_F(iNodeX_X1,iGF_Gm_dd_11), &
                    G_F(iNodeX_X1,iGF_Gm_dd_22), &
                    G_F(iNodeX_X1,iGF_Gm_dd_33) )

          END DO

          ! --- Right State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                   uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                   uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                   uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R(:) )

          DO iNodeX_X1 = 1, nDOFX_X1

            Flux_X1_R(iNodeX_X1,1:nCF) &
              = Flux_X1 &
                  ( uPF_R(iNodeX_X1,iPF_D ), &
                    uPF_R(iNodeX_X1,iPF_V1), &
                    uPF_R(iNodeX_X1,iPF_V2), &
                    uPF_R(iNodeX_X1,iPF_V3), &
                    uPF_R(iNodeX_X1,iPF_E ), &
                    uPF_R(iNodeX_X1,iPF_Ne), &
                    P_R(iNodeX_X1), &
                    G_F(iNodeX_X1,iGF_Gm_dd_11), &
                    G_F(iNodeX_X1,iGF_Gm_dd_22), &
                    G_F(iNodeX_X1,iGF_Gm_dd_33) )

          END DO

          ! --- Numerical Flux ---

          CALL Timer_Start( dT_FLX_N )

          DO iNodeX_X1 = 1, nDOFX_X1

            AlphaPls &
              = AlphaPlus &
                  ( uPF_L(iNodeX_X1,iPF_V1), Cs_L(iNodeX_X1), &
                    uPF_R(iNodeX_X1,iPF_V1), Cs_R(iNodeX_X1), &
                    G_F(iNodeX_X1,iGF_Gm_dd_11) )

            AlphaMns &
              = AlphaMinus &
                  ( uPF_L(iNodeX_X1,iPF_V1), Cs_L(iNodeX_X1), &
                    uPF_R(iNodeX_X1,iPF_V1), Cs_R(iNodeX_X1), &
                    G_F(iNodeX_X1,iGF_Gm_dd_11) )

            AlphaMdl &
              = AlphaMiddle &
                  ( uCF_L(iNodeX_X1,iCF_D ), Flux_X1_L(iNodeX_X1,iCF_D ), &
                    uCF_L(iNodeX_X1,iCF_S1), Flux_X1_L(iNodeX_X1,iCF_S1), &
                    uCF_R(iNodeX_X1,iCF_D ), Flux_X1_R(iNodeX_X1,iCF_D ), &
                    uCF_R(iNodeX_X1,iCF_S1), Flux_X1_R(iNodeX_X1,iCF_S1), &
                    AlphaPls, AlphaMns, G_F(iNodeX_X1,iGF_Gm_dd_11) )

            NumericalFlux(iNodeX_X1,:) &
!              = NumericalFlux_X1_HLLC &
              = NumericalFlux_HLL &
                  ( uCF_L(iNodeX_X1,:), Flux_X1_L(iNodeX_X1,:), &
                    uCF_R(iNodeX_X1,:), Flux_X1_R(iNodeX_X1,:), &
                    AlphaPls, AlphaMns, AlphaMdl, G_F(iNodeX_X1,iGF_Gm_dd_11) )

          END DO

          DO iCF = 1, nCF

            NumericalFlux(:,iCF) &
              = dX2 * dX3 * WeightsX_X1(:) * G_F(:,iGF_SqrtGm) &
                  * NumericalFlux(:,iCF)

          END DO

          CALL Timer_Stop( dT_FLX_N )

          CALL Timer_Add( Timer_FLX_N, dT_FLX_N )

          ! --- Contribution to This Element ---

          CALL Timer_Start( dT_RHS_2 )

          IF( iX1 < iX_E0(1) + 1 )THEN

            DO iCF = 1, nCF

              CALL DGEMV &
                     ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                       NumericalFlux(:,iCF), 1, One, &
                       dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          CALL Timer_Stop( dT_RHS_2 )

          CALL Timer_Add( Timer_RHS_2, dT_RHS_2 )

          ! --- Contribution to Previous Element ---

          CALL Timer_Start( dT_RHS_3 )

          IF( iX1 > iX_B0(1) )THEN

            DO iCF = 1, nCF

              CALL DGEMV &
                     ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                       NumericalFlux(:,iCF), 1, One, &
                       dU(:,iX1-1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          CALL Timer_Stop( dT_RHS_3 )

          CALL Timer_Add( Timer_RHS_3, dT_RHS_3 )

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

    CALL Timer_Stop( Timer_RHS )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A24,ES10.4E2)') &
        '', 'ComputeRHS_Euler: ', Timer_RHS
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 1: ', Timer_RHS_1
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 2: ', Timer_RHS_2
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 3: ', Timer_RHS_3
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT F: ', Timer_INT_F
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT G: ', Timer_INT_G
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'FLX N: ', Timer_FLX_N
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'Sum: ', Timer_RHS_1+Timer_RHS_2+Timer_RHS_3+Timer_INT_F &
        + Timer_INT_G + Timer_FLX_N

    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X2
    REAL(DP) :: dX1, dX3
    REAL(DP) :: AlphaPls
    REAL(DP) :: AlphaMns
    REAL(DP) :: AlphaMdl
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: P_L(nDOFX_X2), Cs_L(nDOFX_X2)
    REAL(DP) :: P_R(nDOFX_X2), Cs_R(nDOFX_X2)
    REAL(DP) :: uCF_L(nDOFX_X2,nCF), uCF_R(nDOFX_X2,nCF)
    REAL(DP) :: uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uPF_L(nDOFX_X2,nPF)
    REAL(DP) :: uPF_R(nDOFX_X2,nPF)
    REAL(DP) :: G_F(nDOFX_X2,nGF)
    REAL(DP) :: G_P(nDOFX,nGF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: Flux_X2_q(nDOFX,nCF)
    REAL(DP) :: Flux_X2_L(nDOFX_X2,nCF)
    REAL(DP) :: Flux_X2_R(nDOFX_X2,nCF)
    REAL(DP) :: NumericalFlux(nDOFX_X2,nCF)

    IF( iX_E0(2) == iX_B0(2) ) RETURN

    CALL Timer_Start( Timer_Div_X2 )

    DO iX3 = iX_B0(3), iX_E0(3)

      dX3 = MeshX(3) % Width(iX3)

      DO iX2 = iX_B0(2), iX_E0(2) + 1

        DO iX1 = iX_B0(1), iX_E0(1)

          dX1 = MeshX(1) % Width(iX1)

          DO iCF = 1, nCF

            uCF_P(:,iCF) = U(:,iX1,iX2-1,iX3,iCF)
            uCF_K(:,iCF) = U(:,iX1,iX2,  iX3,iCF)

          END DO

          DO iGF = 1, nGF

            G_P(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
            G_K(:,iGF) = G(:,iX1,iX2,  iX3,iGF)

          END DO

          !--------------------
          ! --- Volume Term ---
          !--------------------

          IF( iX2 < iX_E0(2) + 1 )THEN

            CALL ComputePrimitive &
                   ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                     uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                     uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                     uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                     G_K(:,iGF_Gm_dd_11), &
                     G_K(:,iGF_Gm_dd_22), &
                     G_K(:,iGF_Gm_dd_33) )

            CALL ComputePressureFromPrimitive &
                   ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

            DO iNodeX = 1, nDOFX

              Flux_X2_q(iNodeX,1:nCF) &
                = Flux_X2 &
                    ( uPF_K(iNodeX,iPF_D ), &
                      uPF_K(iNodeX,iPF_V1), &
                      uPF_K(iNodeX,iPF_V2), &
                      uPF_K(iNodeX,iPF_V3), &
                      uPF_K(iNodeX,iPF_E ), &
                      uPF_K(iNodeX,iPF_Ne), &
                      P_K  (iNodeX), &
                      G_K(iNodeX,iGF_Gm_dd_11), &
                      G_K(iNodeX,iGF_Gm_dd_22), &
                      G_K(iNodeX,iGF_Gm_dd_33) )

            END DO

            DO iCF = 1, nCF

              Flux_X2_q(:,iCF) &
                = dX1 * dX3 * WeightsX_q(:) * G_K(:,iGF_SqrtGm) &
                    * Flux_X2_q(:,iCF)

              CALL DGEMV &
                     ( 'T', nDOFX, nDOFX, One, dLXdX2_q, nDOFX, &
                       Flux_X2_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          !---------------------
          ! --- Surface Term ---
          !---------------------

          ! --- Interpolate Fluid Fields ---

          DO iCF = 1, nCF

            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                     uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                     uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

          END DO

          ! --- Interpolate Geometry Fields ---

          ! --- Face States (Average of Left and Right States) ---

          DO iGF = iGF_h_1, iGF_h_3

            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                     G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                     G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

            G_F(1:nDOFX_X2,iGF) &
              = MAX( G_F(1:nDOFX_X2,iGF), SqrtTiny )

          END DO

          CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

          ! --- Left State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
                   uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
                   uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
                   uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L(:) )

          DO iNodeX_X2 = 1, nDOFX_X2

            Flux_X2_L(iNodeX_X2,1:nCF) &
              = Flux_X2 &
                  ( uPF_L(iNodeX_X2,iPF_D ), &
                    uPF_L(iNodeX_X2,iPF_V1), &
                    uPF_L(iNodeX_X2,iPF_V2), &
                    uPF_L(iNodeX_X2,iPF_V3), &
                    uPF_L(iNodeX_X2,iPF_E ), &
                    uPF_L(iNodeX_X2,iPF_Ne), &
                    P_L  (iNodeX_X2), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_33) )

          END DO

          ! --- Right State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                   uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                   uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                   uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R(:) )

          DO iNodeX_X2 = 1, nDOFX_X2

            Flux_X2_R(iNodeX_X2,1:nCF) &
              = Flux_X2 &
                  ( uPF_R(iNodeX_X2,iPF_D ), &
                    uPF_R(iNodeX_X2,iPF_V1), &
                    uPF_R(iNodeX_X2,iPF_V2), &
                    uPF_R(iNodeX_X2,iPF_V3), &
                    uPF_R(iNodeX_X2,iPF_E ), &
                    uPF_R(iNodeX_X2,iPF_Ne), &
                    P_R  (iNodeX_X2), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                    G_F  (iNodeX_X2,iGF_Gm_dd_33) )

          END DO

          ! --- Numerical Flux ---

          DO iNodeX_X2 = 1, nDOFX_X2

            AlphaPls &
              = AlphaPlus &
                  ( uPF_L(iNodeX_X2,iPF_V2), Cs_L(iNodeX_X2), &
                    uPF_R(iNodeX_X2,iPF_V2), Cs_R(iNodeX_X2), &
                    G_F(iNodeX_X2,iGF_Gm_dd_22) )

            AlphaMns &
              = AlphaMinus &
                  ( uPF_L(iNodeX_X2,iPF_V2), Cs_L(iNodeX_X2), &
                    uPF_R(iNodeX_X2,iPF_V2), Cs_R(iNodeX_X2), &
                    G_F(iNodeX_X2,iGF_Gm_dd_22) )

            AlphaMdl &
              = AlphaMiddle &
                  ( uCF_L(iNodeX_X2,iCF_D ), Flux_X2_L(iNodeX_X2,iCF_D ), &
                    uCF_L(iNodeX_X2,iCF_S2), Flux_X2_L(iNodeX_X2,iCF_S2), &
                    uCF_R(iNodeX_X2,iCF_D ), Flux_X2_R(iNodeX_X2,iCF_D ), &
                    uCF_R(iNodeX_X2,iCF_S2), Flux_X2_R(iNodeX_X2,iCF_S2), &
                    AlphaPls, AlphaMns, G_F(iNodeX_X2,iGF_Gm_dd_22) )

            NumericalFlux(iNodeX_X2,:) &
!              = NumericalFlux_X2_HLLC &
              = NumericalFlux_HLL &
                  ( uCF_L(iNodeX_X2,:), Flux_X2_L(iNodeX_X2,:), &
                    uCF_R(iNodeX_X2,:), Flux_X2_R(iNodeX_X2,:), &
                    AlphaPls, AlphaMns, AlphaMdl, G_F(iNodeX_X2,iGF_Gm_dd_22) )

          END DO

          DO iCF = 1, nCF

            NumericalFlux(:,iCF) &
              = dX1 * dX3 * WeightsX_X2(:) * G_F(:,iGF_SqrtGm) &
                  * NumericalFlux(:,iCF)

          END DO

          ! --- Contribution to This Element ---

          IF( iX2 < iX_E0(2) + 1 )THEN

            DO iCF = 1, nCF

              CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Dn, nDOFX_X2, &
                          NumericalFlux(:,iCF), 1, One, &
                          dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          ! --- Contribution to Previous Element ---

          IF( iX2 > iX_B0(2) )THEN

            DO iCF = 1, nCF

              CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Up, nDOFX_X2, &
                          NumericalFlux(:,iCF), 1, &
                          One, dU(:,iX1,iX2-1,iX3,iCF), 1 )

            END DO

          END IF

        END DO
      END DO
    END DO

    CALL Timer_Stop( Timer_Div_X2 )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A32,ES10.4E2)') &
        '', 'ComputeIncrement_Divergence_X2: ', Timer_Div_X2

    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    IF( iX_E0(3) == iX_B0(3) ) RETURN

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh2dX1(nDOFX), dh3dX1(nDOFX), dh3dX2(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF)
    REAL(DP) :: G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)
    REAL(DP) :: G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF)

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    CALL Timer_Start( Timer_Geo )

    DO iX3 = iX_B0(3), iX_E0(3)

      DO iX2 = iX_B0(2), iX_E0(2)

        dX2 = MeshX(2) % Width(iX2)

        DO iX1 = iX_B0(1), iX_E0(1)

          dX1 = MeshX(1) % Width(iX1)

!          print*,"iX1, iX2, iX3 = ", iX1, iX2, iX3

          DO iCF = 1, nCF

            uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
            G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
            G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)
            G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
            G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)

          END DO

          CALL ComputePrimitive &
                 ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   G_K(:,iGF_Gm_dd_11), &
                   G_K(:,iGF_Gm_dd_22), &
                   G_K(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

          DO iNodeX = 1, nDOFX

            Stress(iNodeX,1:3) &
              = StressTensor_Diagonal &
                  ( uPF_K(iNodeX,iPF_D ), &
                    uPF_K(iNodeX,iPF_V1), &
                    uPF_K(iNodeX,iPF_V2), &
                    uPF_K(iNodeX,iPF_V3), &
                    P_K  (iNodeX), &
                    G_K(iNodeX,iGF_Gm_dd_11), &
                    G_K(iNodeX,iGF_Gm_dd_22), &
                    G_K(iNodeX,iGF_Gm_dd_33) )

          END DO

          ! --- Scale Factor Derivatives wrt X1 ---

          ! --- Face States (Average of Left and Right States) ---

          DO iGF = iGF_h_2, iGF_h_3

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

            G_X1_Dn(1:nDOFX_X1,iGF) &
              = MAX( G_X1_Dn(1:nDOFX_X1,iGF), SqrtTiny )

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

            G_X1_Up(1:nDOFX_X1,iGF) &
              = MAX( G_X1_Up(1:nDOFX_X1,iGF), SqrtTiny )

          END DO

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

          dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

          dh3dx1 = dh3dx1 / ( WeightsX_q(:) * dX1 )

          ! --- Scale Factor Derivatives wrt X2 ---

          ! --- Face States (Average of Left and Right States) ---

          DO iGF = iGF_h_3, iGF_h_3

            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                     G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                     G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

            G_X2_Dn(1:nDOFX_X2,iGF) &
              = MAX( G_X2_Dn(1:nDOFX_X2,iGF), SqrtTiny )

            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                     G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                     G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

            G_X2_Up(1:nDOFX_X2,iGF) &
              = MAX( G_X2_Up(1:nDOFX_X2,iGF), SqrtTiny )

          END DO

          CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                      WeightsX_X2(:) * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )
          CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                      WeightsX_X2(:) * G_X2_Dn(:,iGF_h_3), 1,  One, dh3dX2, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX2, 1 )

          dh3dx2 = dh3dx2 / ( WeightsX_q(:) * dX2 )

          ! --- Compute Increments ---

          dU(:,iX1,iX2,iX3,iCF_S1) &
            = dU(:,iX1,iX2,iX3,iCF_S1) &
                + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)  &
                + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3)

          dU(:,iX1,iX2,iX3,iCF_S2) &
            = dU(:,iX1,iX2,iX3,iCF_S2) &
                + ( Stress(:,3) * dh3dX2(:) ) / G_K(:,iGF_h_3)

        END DO
      END DO
    END DO

    CALL Timer_Stop( Timer_Geo )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES10.4E2)') &
        '', 'ComputeIncrement_Geometry: ', Timer_Geo

    END IF

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Gravity &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: Phi_P_X1(nDOFX)
    REAL(DP) :: Phi_K   (nDOFX)
    REAL(DP) :: Phi_N_X1(nDOFX)
    REAL(DP) :: dPhidX1(nDOFX)
    REAL(DP) :: Phi_X1_Dn(nDOFX_X1)
    REAL(DP) :: Phi_X1_Up(nDOFX_X1)
    REAL(DP) :: uCF_K(nDOFX,nCF)

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          dX1 = MeshX(1) % Width(iX1)
          dX2 = MeshX(2) % Width(iX2)
          dX3 = MeshX(3) % Width(iX3)

          DO iCF = 1, nCF

            uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

          END DO

          Phi_P_X1(:) = G(:,iX1-1,iX2,iX3,iGF_Phi_N)
          Phi_K   (:) = G(:,iX1,  iX2,iX3,iGF_Phi_N)
          Phi_N_X1(:) = G(:,iX1+1,iX2,iX3,iGF_Phi_N)

          ! --- Derivative of Potential wrt X1 ---

          ! --- Face States (Average of Left and Right States) ---

          ! --- Face at X1_L ---

          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                   Phi_P_X1(:), 1, Zero, Phi_X1_Dn(:), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                   Phi_K   (:), 1, Half, Phi_X1_Dn(:), 1 )

          ! --- Face at X1_H ---

          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                   Phi_K   (:), 1, Zero, Phi_X1_Up(:), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                   Phi_N_X1(:), 1, Half, Phi_X1_Up(:), 1 )

          ! --- dPhidX1 ---

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * Phi_X1_Up(:), 1, Zero, dPhidX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * Phi_X1_Dn(:), 1,  One, dPhidX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * Phi_K    (:), 1,  One, dPhidX1, 1 )

          dPhidX1 = dPhidX1 / ( WeightsX_q(:) * dX1 )

          ! --- Increments ---

          dU(:,iX1,iX2,iX3,iCF_S1) &
            = dU(:,iX1,iX2,iX3,iCF_S1) &
                - uCF_K(:,iCF_D) * dPhidX1(:)

          dU(:,iX1,iX2,iX3,iCF_E) &
            = dU(:,iX1,iX2,iX3,iCF_E) &
                - uCF_K(:,iCF_S1) * dPhidX1(:)

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeIncrement_Gravity


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


END MODULE Euler_dgDiscretizationModule
