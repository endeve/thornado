MODULE ScalarWave_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, WeightsX_X1, &
    nDOFX_X2, WeightsX_X2, &
    nDOFX_X3, WeightsX_X3, &
    WeightsX_q,            &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
 USE ScalarFieldsModule, ONLY: &
    nSF, iSF_u, iSF_v, &
    Flux_X1, &
    NumericalFlux_X1_Upwind, &
    Source_X1
 USE ScalarWave_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_ScalarWave

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_ScalarWave_DG_Explicit


CONTAINS


  SUBROUTINE ComputeIncrement_ScalarWave_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)         :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)           :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    INTEGER  :: iX1, iX2, iX3, iSF
    REAL(DP) :: dX1, dX2, dX3

    dU = Zero

    CALL ApplyBoundaryConditions_ScalarWave &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    CALL ComputeIncrement_SourceTerms_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    ! --- Multiply Inverse Mass Matrix ---

    DO iSF = 1, nSF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        dX1 = MeshX(1) % Width(iX1)
        dX2 = MeshX(2) % Width(iX2)
        dX3 = MeshX(3) % Width(iX3)

        dU(:,iX1,iX2,iX3,iSF) &
          = dU(:,iX1,iX2,iX3,iSF) &
              / ( WeightsX_q * dX1 * dX2 * dX3 )

      END DO
      END DO
      END DO
    END DO

  END SUBROUTINE ComputeIncrement_ScalarWave_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nSF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nSF)

    INTEGER  :: iX1, iX2, iX3, iSF, iNodeX, iNodeX_X1
    REAL(DP) :: dX2, dX3
    REAL(DP) :: uSF_L(nDOFX_X1,nSF), uSF_R(nDOFX_X1,nSF)
    REAL(DP) :: uSF_P(nDOFX,nSF), uSF_K(nDOFX,nSF)
    REAL(DP) :: Flux_X1_L(nDOFX_X1,nSF), Flux_X1_R(nDOFX_X1,nSF)
    REAL(DP) :: Flux_X1_q(nDOFX,nSF)
    REAL(DP) :: NumericalFlux(nDOFX_X1,nSF)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    !$OMP PARALLEL DO PRIVATE &
    !$OMP& ( iX1, iX2, iX3, iSF, iNodeX, iNodeX_X1, dX2, dX3, &
    !$OMP&   uSF_P, uSF_K, uSF_L, uSF_R, &
    !$OMP&   Flux_X1_q, Flux_X1_L, Flux_X1_R, NumericalFlux )
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1) + 1

      dX3 = MeshX(3) % Width(iX3)
      dX2 = MeshX(2) % Width(iX2)

      DO iSF = 1, nSF

        uSF_P(:,iSF) = U(:,iX1-1,iX2,iX3,iSF)
        uSF_K(:,iSF) = U(:,iX1,  iX2,iX3,iSF)

      END DO

      !--------------------
      ! --- Volume Term ---
      !--------------------

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iNodeX = 1, nDOFX

          Flux_X1_q(iNodeX,:) &
            = Flux_X1 &
                ( uSF_K(iNodeX,iSF_U ),       &
                  uSF_K(iNodeX,iSF_V) )

        END DO

        DO iSF = 1, nSF

          Flux_X1_q(:,iSF) &
            = dX2 * dX3 * WeightsX_q * Flux_X1_q(:,iSF)

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX1_q, nDOFX, &
                   Flux_X1_q(:,iSF), 1, One, dU(:,iX1,iX2,iX3,iSF), 1 )
        END DO

      END IF

      !---------------------
      ! --- Surface Term ---
      !---------------------

      ! --- Interpolate Scalar Fields ---

      DO iSF = 1, nSF

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 uSF_P(:,iSF), 1, Zero, uSF_L(:,iSF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 uSF_K(:,iSF), 1, Zero, uSF_R(:,iSF), 1 )

      END DO

      ! --- Left State ---

      DO iNodeX_X1 = 1, nDOFX_X1

        Flux_X1_L(iNodeX_X1,:) &
          = Flux_X1 &
              ( uSF_L(iNodeX_X1,iSF_U ),       &
                uSF_L(iNodeX_X1,iSF_V) )

      END DO

      ! --- Right State ---


      DO iNodeX_X1 = 1, nDOFX_X1

        Flux_X1_R(iNodeX_X1,:) &
          = Flux_X1 &
              ( uSF_R(iNodeX_X1,iSF_U ),       &
                uSF_R(iNodeX_X1,iSF_V) )

      END DO

      ! --- Numerical Flux ---

      DO iNodeX_X1 = 1, nDOFX_X1

        NumericalFlux(iNodeX_X1,:) &
          = NumericalFlux_X1_Upwind &
              ( uSF_L    (iNodeX_X1,:), &
                uSF_R    (iNodeX_X1,:), &
                Flux_X1_L(iNodeX_X1,:), &
                Flux_X1_R(iNodeX_X1,:), &
                [1.0_dp, -1.0_dp] )

      END DO

      DO iSF = 1, nSF

        NumericalFlux(:,iSF) &
          = dX2 * dX3 * WeightsX_X1 * NumericalFlux(:,iSF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iSF = 1, nSF

          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                   NumericalFlux(:,iSF), 1, One, dU(:,iX1,iX2,iX3,iSF), 1 )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX1 .GT. iX_B0(1) )THEN

        DO iSF = 1, nSF

          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                   NumericalFlux(:,iSF), 1, One, dU(:,iX1-1,iX2,iX3,iSF), 1 )
        END DO

      END IF

    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_SourceTerms_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nSF)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nSF)

    INTEGER  :: iX1, iX2, iX3, iSF, iNodeX, iNodeX_X1
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: uSF_K(nDOFX,nSF)
    REAL(DP) :: Source_X1_q(nDOFx, nSF)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    !$OMP PARALLEL DO PRIVATE &
    !$OMP& ( iX1, iX2, iX3, iSF, iNodeX, &
    !$OMP&   iNodeX_X1,dX1, dX2, dX3, uSF_K )
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1) + 1

      dX3 = MeshX(3) % Width(iX3)
      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

      DO iSF = 1, nSF

        uSF_K(:,iSF) = U(:,iX1,iX2,iX3,iSF)

      END DO

      !--------------------
      ! --- Source Term ---
      !--------------------

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iNodeX = 1, nDOFX

          Source_X1_q(iNodeX,:) &
            = Source_X1 &
                ( uSF_K(iNodeX,iSF_U ),       &
                  uSF_K(iNodeX,iSF_V) )

        END DO

        DO iSF = 1, nSF

          Source_X1_q(:,iSF) &
            = dX1 * dX2 * dX3 * WeightsX_q * Source_X1_q(:,iSF)
     
          dU(:,iX1,iX2,iX3,iSF) = Source_X1_q(:,iSF) + dU(:,iX1,iX2,iX3,iSF)

        END DO

      END IF
    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

   END SUBROUTINE ComputeIncrement_SourceTerms_X1


END MODULE ScalarWave_dgDiscretizationModule
