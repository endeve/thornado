MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Third, Three, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFX, nDOFE, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    Euler_ComputeConserved_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeConserved_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields


CONTAINS


  SUBROUTINE InitializeFields( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( V_0 )

      CASE( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion( V_0 )

    END SELECT


  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeZ2
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1

#ifndef MOMENT_CLOSURE_MINERBO

    WRITE(*,*)
    WRITE(*,'(A8,A)') &
      '', 'Must use Minerbo closure for this application'
    STOP

#endif

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        iNodeZ2 = NodeNumberTable(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I1,iS) &
          = uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_SineWaveDiffusion( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeZ2
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        iNodeZ2 = NodeNumberTable(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 0.49_DP * SIN( Third * Pi * X1 ) + 0.5_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I1,iS) &
          = - ( 0.49_DP * Pi / ( 9.0_DP * 1.0d2 ) ) * COS( Third * Pi * X1 )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ3,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SineWaveDiffusion


END MODULE InitializationModule
