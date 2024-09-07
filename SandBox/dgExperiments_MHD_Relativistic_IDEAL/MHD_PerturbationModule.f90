MODULE MHD_PerturbationModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Pi, &
    FourPi
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nPM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nAM, &
    iAM_P, &
    iAM_T, &
    iAM_Ye, &
    iAM_S, &
    iAM_E, &
    iAM_Gm, &
    iAM_Cs, &
    nDM, &
    iDM_Sh_X1, &
    iDM_Sh_X2, &
    iDM_Sh_X3, &
    iDM_Div
   USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic, &
    ComputeFromConserved_MHD_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRandPerturbations, &
            ApplyRandPerturbations, &
            FinalizeRandPerturbations

  REAL(DP), PUBLIC :: Rand_Amplitude
  REAL(DP), PUBLIC, ALLOCATABLE :: Random_r(:,:,:,:), Random_z(:,:,:,:), &
                                   Random_theta(:,:,:,:)

CONTAINS


  SUBROUTINE InitializeRandPerturbations &
               ( iX_B0, iX_E0, nDOFX, Rand_Amplitude_Option )

    INTEGER,  INTENT(in) :: iX_B0(3), iX_E0(3), nDOFX
    REAL(DP), INTENT(in), OPTIONAL :: Rand_Amplitude_Option

    INTEGER :: iNX, iX1, iX2, iX3
    REAL(DP) :: Rand_r, Rand_z, Rand_theta

    IF( PRESENT( Rand_Amplitude_Option ) )THEN
      Rand_Amplitude = Rand_Amplitude_Option
    ELSE
      Rand_Amplitude = Zero
    END IF

    ALLOCATE(Random_r(nDOFX, iX_E0(1) - iX_B0(1) + 1, &
                             iX_E0(2) - iX_B0(2) + 1, &
                             iX_E0(3) - iX_B0(3) + 1))
    ALLOCATE(Random_z(nDOFX, iX_E0(1) - iX_B0(1) + 1, &
                             iX_E0(2) - iX_B0(2) + 1, &
                             iX_E0(3) - iX_B0(3) + 1))
    ALLOCATE(Random_theta(nDOFX, iX_E0(1) - iX_B0(1) + 1, &
                                 iX_E0(2) - iX_B0(2) + 1, &
                                 iX_E0(3) - iX_B0(3) + 1))

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      CALL RANDOM_SEED()

      CALL RANDOM_NUMBER( Rand_r )

      CALL RANDOM_SEED()

      CALL RANDOM_NUMBER( Rand_z )

      CALL RANDOM_SEED()

      CALL RANDOM_NUMBER( Rand_theta )

      Random_r(iNX,iX1,iX2,iX3)     = Two * Rand_r - One
      Random_z(iNX,iX1,iX2,iX3)     = Two * Rand_z - One
      Random_theta(iNX,iX1,iX2,iX3) = Two * Rand_theta - One

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRandPerturbations


  SUBROUTINE ApplyRandPerturbations( iX_B0, iX_E0, iX_B1, iX_E1, nDOFX, &
                                     G, U, EvolveOnlyMagnetic )

    ! --- Still in progress. ---

    INTEGER, INTENT(in) :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), nDOFX
    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2

    REAL(DP) :: &
      P(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nPM)
    REAL(DP) :: &
      A(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nAM)

    REAL(DP) :: X1, X2
    REAL(DP) :: kz

    ! --- Applying the random radial velocity       ---
    ! --- perturbations from Rembiasz et al. (2016) ---

    CALL ComputeFromConserved_MHD_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A, EvolveOnlyMagnetic )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1, iNX)
      iNodeX2 = NodeNumberTableX(2, iNX)

      X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

      kz = Two * Pi / ( Half * Kilometer )

      P(iNX,iX1,iX2,iX3,iPM_V1) &
        = ( 0.1_DP * Rand_Amplitude * Random_r(iNX,iX1,iX2,iX3) &
            + 0.2d-5 * SIN( kz * X2 ) ) &
          * X1 * P(iNX,iX1,iX2,iX3,iPM_V3)

      P(iNX,iX1,iX2,iX3,iPM_V2) &
        = Rand_Amplitude * Random_z(iNX,iX1,iX2,iX3) &
          * X1 * P(iNX,iX1,iX2,iX3,iPM_V3)

      P(iNX,iX1,iX2,iX3,iPM_V3) &
        = ( One + Rand_Amplitude * Random_theta(iNX,iX1,iX2,iX3) ) &
          * P(iNX,iX1,iX2,iX3,iPM_V3)

    END DO

     CALL ComputeConserved_MHD_Relativistic &
             ( P(:,iX1,iX2,iX3,iPM_D ), P(:,iX1,iX2,iX3,iPM_V1 ), &
               P(:,iX1,iX2,iX3,iPM_V2), P(:,iX1,iX2,iX3,iPM_V3 ), &
               P(:,iX1,iX2,iX3,iPM_E ), P(:,iX1,iX2,iX3,iPM_Ne ), &
               P(:,iX1,iX2,iX3,iPM_B1), P(:,iX1,iX2,iX3,iPM_B2 ), &
               P(:,iX1,iX2,iX3,iPM_B3), P(:,iX1,iX2,iX3,iPM_Chi), &
               U(:,iX1,iX2,iX3,iCM_D ), U(:,iX1,iX2,iX3,iCM_S1 ), &
               U(:,iX1,iX2,iX3,iCM_S2), U(:,iX1,iX2,iX3,iCM_S3 ), &
               U(:,iX1,iX2,iX3,iCM_E ), U(:,iX1,iX2,iX3,iCM_Ne ), &
               U(:,iX1,iX2,iX3,iCM_B1), U(:,iX1,iX2,iX3,iCM_B2 ), &
               U(:,iX1,iX2,iX3,iCM_B3), U(:,iX1,iX2,iX3,iCM_Chi), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(:,iX1,iX2,iX3,iGF_Alpha   ), &
               G(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               G(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               G(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               A(:,iX1,iX2,iX3,iAM_P       ), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyRandPerturbations


  SUBROUTINE FinalizeRandPerturbations

    DEALLOCATE( Random_r )
    DEALLOCATE( Random_z )
    DEALLOCATE( Random_theta )

  END SUBROUTINE FinalizeRandPerturbations


END MODULE MHD_PerturbationModule
