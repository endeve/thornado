MODULE dgDiscretizationModule_Euler

  USE KindModule, ONLY: &
    DP, Zero, Half, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    WeightsX_q, &
    WeightsX_X1, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    LX_X1_Dn, &
    LX_X1_Up
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    rhsCF
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputePrimitive, &
    Eigenvalues, &
    AlphaPlus, &
    AlphaMinus, &
    AlphaMiddle, &
    Flux_X1, &
    NumericalFlux_X1_HLLC
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeRHS_Euler

  LOGICAL, PARAMETER :: DisplayTimers = .TRUE.
  REAL(DP) :: Timer_RHS
  REAL(DP) :: Timer_RHS_1, dT_RHS_1
  REAL(DP) :: Timer_RHS_2, dT_RHS_2
  REAL(DP) :: Timer_RHS_3, dT_RHS_3
  REAL(DP) :: Timer_INT_F, dT_INT_F
  REAL(DP) :: Timer_INT_G, dT_INT_G
  REAL(DP) :: Timer_FLX_N, dT_FLX_N

CONTAINS


  SUBROUTINE ComputeRHS_Euler

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X1, iNodeX1
    REAL(DP) :: AlphaPls, AlphaMns, AlphaMdl
    REAL(DP) :: ErrorL1, ErrorIn, Error, X1
    REAL(DP), DIMENSION(nDOFX_X1)     :: Pressure_L, SoundSpeed_L, Lambda_L
    REAL(DP), DIMENSION(nDOFX_X1)     :: Pressure_R, SoundSpeed_R, Lambda_R
    REAL(DP), DIMENSION(nDOFX)        :: Pressure_K
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(nDOFX_X1,nGF) :: uGF_L, uGF_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_L
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: NumericalFlux
    REAL(DP), DIMENSION(nDOFX,   nCF) :: uCF_P, uCF_K
    REAL(DP), DIMENSION(nDOFX,   nGF) :: uGF_P, uGF_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: Flux_X1_q
    REAL(DP), DIMENSION(nDOFX_X1,nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(nDOFX,   nPF) :: uPF_K

    Timer_RHS_1 = 0.0_DP
    Timer_RHS_2 = 0.0_DP
    Timer_RHS_3 = 0.0_DP
    Timer_INT_F = 0.0_DP
    Timer_INT_G = 0.0_DP
    Timer_FLX_N = 0.0_DP

    CALL Timer_Start( Timer_RHS )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1) + 1

          DO iCF = 1, nCF

            uCF_P(:,iCF) = uCF(:,iX1-1,iX2,iX3,iCF)
            uCF_K(:,iCF) = uCF(:,iX1,  iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            uGF_P(:,iGF) = uGF(:,iX1-1,iX2,iX3,iGF)
            uGF_K(:,iGF) = uGF(:,iX1,  iX2,iX3,iGF)

          END DO

          !--------------------
          ! --- Volume Term ---
          !--------------------

          CALL ComputePrimitive &
                 ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   Pressure_K )

          DO iNodeX = 1, nDOFX

            Flux_X1_q(iNodeX,1:nCF) &
              = Flux_X1 &
                  ( uPF_K(iNodeX,iPF_D ), &
                    uPF_K(iNodeX,iPF_V1), &
                    uPF_K(iNodeX,iPF_V2), &
                    uPF_K(iNodeX,iPF_V3), &
                    uPF_K(iNodeX,iPF_E ), &
                    uPF_K(iNodeX,iPF_Ne), &
                    Pressure_K(iNodeX), &
                    uGF_K(iNodeX,iGF_Gm_dd_11), &
                    uGF_K(iNodeX,iGF_Gm_dd_22), &
                    uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO

          CALL Timer_Start( dT_RHS_1 )

          DO iCF = 1, nCF

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = MATMUL &
                  ( TRANSPOSE( dLXdX1_q(1:nDOFX,1:nDOFX) ), &
                    WeightsX_q(1:nDOFX) * Flux_X1_q(1:nDOFX,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_1 )

          CALL Timer_Add( Timer_RHS_1, dT_RHS_1 )

          !------------------------
          ! --- Divergence Term ---
          !------------------------

          ! --- Interpolate Fluid Fields ---

          CALL Timer_Start( dT_INT_F )

          DO iCF = 1, nCF

            ! --- Left States ---

            uCF_L(1:nDOFX_X1,iCF) &
              = MATMUL( LX_X1_Up(:,:), uCF_P(:,iCF) )

            ! --- Right States ---

            uCF_R(1:nDOFX_X1,iCF) &
              = MATMUL( LX_X1_Dn(:,:), uCF_K(:,iCF) )

          END DO

          CALL Timer_Stop( dT_INT_F )

          CALL Timer_Add( Timer_INT_F, dT_INT_F )

          ! --- Interpolate Geometry Fields ---

          CALL Timer_Start( dT_INT_G )

          ! --- Left States ---

          uGF_L(1:nDOFX_X1,iGF_Gm_dd_11) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_11) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_11) ) )
          uGF_L(1:nDOFX_X1,iGF_Gm_dd_22) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_22) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_22) ) )
          uGF_L(1:nDOFX_X1,iGF_Gm_dd_33) &
            = Half * ( MATMUL( LX_X1_Up(:,:), uGF_P(:,iGF_Gm_dd_33) ) &
                       + MATMUL( LX_X1_Dn(:,:), uGF_K(:,iGF_Gm_dd_33) ) )

          ! --- Right States (Equal to Left States) ---

          uGF_R(1:nDOFX_X1,iGF_Gm_dd_11) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_11)
          uGF_R(1:nDOFX_X1,iGF_Gm_dd_22) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_22)
          uGF_R(1:nDOFX_X1,iGF_Gm_dd_33) &
            = uGF_L(1:nDOFX_X1,iGF_Gm_dd_33)

          CALL Timer_Stop( dT_INT_G )

          CALL Timer_Add( Timer_INT_G, dT_INT_G )

          ! --- Left State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
                   uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
                   uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
                   uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
                   uGF_L(:,iGF_Gm_dd_11), &
                   uGF_L(:,iGF_Gm_dd_22), &
                   uGF_L(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), &
                   Pressure_L(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), &
                   SoundSpeed_L(:) )

          DO iNodeX_X1 = 1, nDOFX_X1

            Flux_X1_L(iNodeX_X1,1:nCF) &
              = Flux_X1 &
                  ( uPF_L(iNodeX_X1,iPF_D ), &
                    uPF_L(iNodeX_X1,iPF_V1), &
                    uPF_L(iNodeX_X1,iPF_V2), &
                    uPF_L(iNodeX_X1,iPF_V3), &
                    uPF_L(iNodeX_X1,iPF_E ), &
                    uPF_L(iNodeX_X1,iPF_Ne), &
                    Pressure_L(iNodeX_X1), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_11), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_22), &
                    uGF_L(iNodeX_X1,iGF_Gm_dd_33) )

          END DO

          ! --- Right State Primitive, etc. ---

          CALL ComputePrimitive &
                 ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                   uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                   uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                   uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                   uGF_R(:,iGF_Gm_dd_11), &
                   uGF_R(:,iGF_Gm_dd_22), &
                   uGF_R(:,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), &
                   Pressure_R(:) )

          CALL ComputeSoundSpeedFromPrimitive &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), &
                   SoundSpeed_R(:) )

          DO iNodeX_X1 = 1, nDOFX_X1

            Flux_X1_R(iNodeX_X1,1:nCF) &
              = Flux_X1 &
                  ( uPF_R(iNodeX_X1,iPF_D ), &
                    uPF_R(iNodeX_X1,iPF_V1), &
                    uPF_R(iNodeX_X1,iPF_V2), &
                    uPF_R(iNodeX_X1,iPF_V3), &
                    uPF_R(iNodeX_X1,iPF_E ), &
                    uPF_R(iNodeX_X1,iPF_Ne), &
                    Pressure_R(iNodeX_X1), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_11), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_22), &
                    uGF_R(iNodeX_X1,iGF_Gm_dd_33) )

          END DO

          ! --- Numerical Flux ---

          CALL Timer_Start( dT_FLX_N )

          DO iNodeX_X1 = 1, nDOFX_X1

            AlphaPls &
              = AlphaPlus &
                  ( uPF_L(iNodeX_X1,iPF_V1), SoundSpeed_L(iNodeX_X1), &
                    uPF_R(iNodeX_X1,iPF_V1), SoundSpeed_R(iNodeX_X1) )

            AlphaMns &
              = AlphaMinus &
                  ( uPF_L(iNodeX_X1,iPF_V1), SoundSpeed_L(iNodeX_X1), &
                    uPF_R(iNodeX_X1,iPF_V1), SoundSpeed_R(iNodeX_X1) )

            AlphaMdl &
              = AlphaMiddle &
                  ( uCF_L(iNodeX_X1,iCF_D ), Flux_X1_L(iNodeX_X1,iCF_D ), &
                    uCF_L(iNodeX_X1,iCF_S1), Flux_X1_L(iNodeX_X1,iCF_S1), &
                    uCF_R(iNodeX_X1,iCF_D ), Flux_X1_R(iNodeX_X1,iCF_D ), &
                    uCF_R(iNodeX_X1,iCF_S1), Flux_X1_R(iNodeX_X1,iCF_S1), &
                    AlphaPls, AlphaMns )

            NumericalFlux(iNodeX_X1,:) &
              = NumericalFlux_X1_HLLC &
                  ( uCF_L(iNodeX_X1,:), Flux_X1_L(iNodeX_X1,:), &
                    uCF_R(iNodeX_X1,:), Flux_X1_R(iNodeX_X1,:), &
                    AlphaPls, AlphaMns, AlphaMdl )

          END DO

          CALL Timer_Stop( dT_FLX_N )

          CALL Timer_Add( Timer_FLX_N, dT_FLX_N )

          ! --- Contribution to This Element ---

          CALL Timer_Start( dT_RHS_2 )

          DO iCF = 1, nCF

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = rhsCF(:,iX1,iX2,iX3,iCF) &
                  + MATMUL &
                      ( TRANSPOSE( LX_X1_Dn(:,:) ), &
                        WeightsX_X1(:) * NumericalFlux(:,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_2 )

          CALL Timer_Add( Timer_RHS_2, dT_RHS_2 )

          ! --- Contribution to Previous Element ---

          CALL Timer_Start( dT_RHS_3 )

          DO iCF = 1, nCF

            rhsCF(:,iX1-1,iX2,iX3,iCF) &
              = rhsCF(:,iX1-1,iX2,iX3,iCF) &
                  - MATMUL &
                      ( TRANSPOSE( LX_X1_Up(:,:) ), &
                        WeightsX_X1(:) * NumericalFlux(:,iCF) )

          END DO

          CALL Timer_Stop( dT_RHS_3 )

          CALL Timer_Add( Timer_RHS_3, dT_RHS_3 )

        END DO
      END DO
    END DO

    ! --- Multiply Inverse Mass Matrix ---

    DO iCF = 1, nCF
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = rhsCF(:,iX1,iX2,iX3,iCF) &
                  / ( WeightsX_q(:) * MeshX(1) % Width(iX1) )

          END DO
        END DO
      END DO
    END DO

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

    ! --- Compute Error ---

    ErrorL1 = 0.0_DP
    ErrorIn = 0.0_DP

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            Error &
              = ABS( - Pi * COS( TwoPi * X1 ) &
                     - rhsCF(iNodeX,iX1,iX2,iX3,1) )

            ErrorL1 = ErrorL1 + Error
            ErrorIn = MAX( ErrorIn, Error )

          END DO

        END DO
      END DO
    END DO

    ErrorL1 = ErrorL1 / REAL( nDOFX*nX(1)*nX(2)*nX(3) )

    WRITE(*,*)
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorL1: ', ErrorL1
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorIn: ', ErrorIn
    WRITE(*,*)

  END SUBROUTINE ComputeRHS_Euler


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE dgDiscretizationModule_Euler
