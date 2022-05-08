MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Three, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFX, nDOFE, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    xL, xR
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, nPR, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeConserved_TwoMoment
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    uOP, iOP_Sigma

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields
  PUBLIC :: ComputeError


CONTAINS


  SUBROUTINE InitializeFields( V_0, LengthScale, Direction, Spectrum )

    REAL(DP),      INTENT(in) :: V_0(3)
    REAL(DP),      INTENT(in) :: LengthScale
    CHARACTER(2),  INTENT(in) :: Direction
    CHARACTER(32), INTENT(in) :: Spectrum

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming &
               ( V_0, Direction )

      CASE( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion( V_0 )

      CASE( 'SphericalDiffusion' )

        CALL InitializeFields_SphericalDiffusion

      CASE( 'IsotropicRadiation' )

        CALL InitializeFields_IsotropicRadiation( V_0 )

      CASE( 'StreamingDopplerShift' )

        CALL InitializeFields_StreamingDopplerShift &
               ( V_0, Direction, Spectrum )

      CASE( 'TransparentTurbulence' )

        CALL InitializeFields_TransparentTurbulence &
               ( V_0, Direction )

      CASE( 'TransparentShock' )

        CALL InitializeFields_TransparentShock &
               ( V_0, LengthScale, Direction )

      CASE( 'TransparentVortex' )

        CALL InitializeFields_TransparentVortex &
               ( V_0, Direction )

      CASE( 'RadiatingSphere' )

        CALL InitializeFields_RadiatingSphere

      CASE( 'GaussianDiffusion' )

        CALL InitializeFields_GaussianDiffusion( V_0 )

      CASE( 'ExpandingAtmosphere' )

        CALL InitializeFields_ExpandingAtmosphere( V_0 )

      CASE( 'HomogeneousSphere1D' )

        CALL InitializeFields_HomogeneousSphere1D

      CASE( 'HomogeneousSphere2D' )

        CALL InitializeFields_HomogeneousSphere2D( V_0 )

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uCF, uPF, uCR, uPR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uCF, uPF, uCR, uPR )
#endif

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( V_0, Direction )

    REAL(DP),     INTENT(in) :: V_0(3)
    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeZ2, iNodeZ3, iNodeZ4
    REAL(DP) :: X1, X2, X3

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

      CALL ComputeConserved_Euler_NonRelativistic &
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

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)
        iNodeZ4 = NodeNumberTableZ(4,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3 )
        X3 = NodeCoordinate( MeshX(3), iZ4, iNodeZ4 )

        SELECT CASE( TRIM( Direction ) )
        CASE( 'X' )

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
            = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        CASE( 'Y' )

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
            = 0.50_DP + 0.49_DP * SIN( TwoPi * X2 )
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        CASE( 'Z' )

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
            = 0.50_DP + 0.49_DP * SIN( TwoPi * X3 )
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

        CASE( 'XY' )

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A8,A)')    '', 'InitializeFields_SineWaveStreaming'
          WRITE(*,'(A8,A,A2)') '', 'Invalid Direction: ', TRIM( Direction )
          WRITE(*,*)
          STOP

        END SELECT

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

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
    REAL(DP) :: X1, Sigma

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

      CALL ComputeConserved_Euler_NonRelativistic &
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

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        Sigma = uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS)

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 0.49_DP * SIN( Third * Pi * X1 ) + 0.5_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = - ( 0.49_DP * Pi / ( 9.0_DP * Sigma ) ) * COS( Third * Pi * X1 )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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


  SUBROUTINE InitializeFields_SphericalDiffusion

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeZ2
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: R, Sigma, t_0 = 1.0d0

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

      CALL ComputeConserved_Euler_NonRelativistic &
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

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

        R = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        Sigma = uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS)

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = ( Sigma / t_0 )**( 1.5_DP ) &
              * EXP( - 3.0_DP * Sigma * R**2 / ( 4.0_DP *t_0 ) )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) * R / ( 2.0_DP * t_0 )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE InitializeFields_SphericalDiffusion


  SUBROUTINE InitializeFields_IsotropicRadiation( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    REAL(DP), PARAMETER :: L_X = 0.25_DP

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeX1
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1_N, X1_E

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1_N = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X1_E = MeshX(1) % Center(iX1) &
                 + Half * MeshX(1) % Width(iX1)

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        IF( X1_E .LE. 0.25_DP )THEN
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * ( X1_N - 0.0_DP ) / L_X
        ELSEIF( X1_E .GT. 0.25_DP .AND. X1_E .LE. 0.75_DP )THEN
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * ( 0.5_DP - X1_N ) / L_X
        ELSE
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * ( X1_N - 1.0_DP ) / L_X
        END IF
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

      CALL ComputeConserved_Euler_NonRelativistic &
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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 0.5_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE InitializeFields_IsotropicRadiation


  SUBROUTINE InitializeFields_StreamingDopplerShift &
    ( V_0, Direction, Spectrum )

    REAL(DP),      INTENT(in) :: V_0(3)
    CHARACTER(2),  INTENT(in) :: Direction
    CHARACTER(32), INTENT(in) :: Spectrum

    REAL(DP), PARAMETER :: X_0 = 2.0_DP
    REAL(DP), PARAMETER :: X_1 = 3.5_DP
    REAL(DP), PARAMETER :: X_2 = 6.5_DP
    REAL(DP), PARAMETER :: X_3 = 8.0_DP
    REAL(DP), PARAMETER :: L_X = 6.0_DP

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP

        SELECT CASE( TRIM( Direction ) )
        CASE( 'X' )

          IF( X1 .LT. X_0 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.0_DP
          ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1)
          ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSE
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.0_DP
          END IF
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)

        CASE( 'Y' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
          IF( X2 .LT. X_0 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP
          ELSEIF( X2 .GE. X_0 .AND. X2 .LT. X_1 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
          ELSEIF( X2 .GE. X_1 .AND. X2 .LT. X_2 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = V_0(2)
          ELSEIF( X2 .GE. X_2 .AND. X2 .LT. X_3 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
          ELSE
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP
          END IF
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)

        CASE( 'Z' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
          IF( X3 .LT. X_0 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP
          ELSEIF( X3 .GE. X_0 .AND. X3 .LT. X_1 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
          ELSEIF( X3 .GE. X_1 .AND. X3 .LT. X_2 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = V_0(3)
          ELSEIF( X3 .GE. X_2 .AND. X3 .LT. X_3 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
          ELSE
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP
          END IF

        END SELECT

        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 1.0d-40
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = Zero
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )
      
         END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL SetInnerBoundary_StreamingDopplerShift( Direction, Spectrum )

  END SUBROUTINE InitializeFields_StreamingDopplerShift


  SUBROUTINE SetInnerBoundary_StreamingDopplerShift( Direction, Spectrum )

    CHARACTER(2),  INTENT(in) :: Direction
    CHARACTER(32), INTENT(in) :: Spectrum

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeE
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: loX(3), hiX(3)
    INTEGER  :: loZ(4), hiZ(4)
    REAL(DP) :: E

    ! --- Fluid Fields ---

    loX = iX_B0
    hiX = iX_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loX(1) = iX_B1(1)
      hiX(1) = iX_B1(1)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loX(2) = iX_B1(2)
      hiX(2) = iX_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loX(3) = iX_B1(3)
      hiX(3) = iX_B1(3)
    END IF

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1), hiX(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    loZ = iZ_B0
    hiZ = iZ_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loZ(2) = iZ_B1(2)
      hiZ(2) = iZ_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loZ(3) = iZ_B1(3)
      hiZ(3) = iZ_B1(3)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loZ(4) = iZ_B1(4)
      hiZ(4) = iZ_B1(4)
    END IF
    
    DO iS  = 1, nSpecies
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2), hiZ(2)
    DO iZ1 = loZ(1), hiZ(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        SELECT CASE( TRIM( Spectrum ) )

          CASE( 'Fermi-Dirac' )

            uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
              = One / ( EXP( E / Three - Three ) + One )

          CASE( 'Bose-Einstein' )

            uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
              = One / ( EXP( E ) - One )

          CASE DEFAULT

            uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
              = One / ( EXP( E / Three - Three ) + One )

          END SELECT

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS)

        END IF

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )
      
         END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE SetInnerBoundary_StreamingDopplerShift


  SUBROUTINE InitializeFields_TransparentTurbulence( V_0, Direction )

    REAL(DP),     INTENT(in) :: V_0(3)
    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, A_0

    WRITE(*,*)
    WRITE(*,'(A6,A13,3ES9.2E2)') '', 'V_0 = ', V_0
    WRITE(*,'(A6,A13,A2)'      ) '', 'Direction = ', TRIM( Direction )
    WRITE(*,*)

    ! --- Fluid Fields ---

    A_0 = 2.5_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One

        SELECT CASE( TRIM( Direction ) )
        CASE( 'X' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * A_0 * SIN( TwoPi * X1 ) * EXP( - ( TwoPi * X1 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = V_0(2) * A_0 * SIN( TwoPi * X1 ) * EXP( - ( TwoPi * X1 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
            = V_0(3) * A_0 * SIN( TwoPi * X1 ) * EXP( - ( TwoPi * X1 )**2 )

        CASE( 'Y' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * A_0 * SIN( TwoPi * X2 ) * EXP( - ( TwoPi * X2 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = V_0(2) * A_0 * SIN( TwoPi * X2 ) * EXP( - ( TwoPi * X2 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
            = V_0(3) * A_0 * SIN( TwoPi * X2 ) * EXP( - ( TwoPi * X2 )**2 )

        CASE( 'Z' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * A_0 * SIN( TwoPi * X3 ) * EXP( - ( TwoPi * X3 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = V_0(2) * A_0 * SIN( TwoPi * X3 ) * EXP( - ( TwoPi * X3 )**2 )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
            = V_0(3) * A_0 * SIN( TwoPi * X3 ) * EXP( - ( TwoPi * X3 )**2 )

        END SELECT

        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0d-0


        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 1.0d-8
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

    CALL SetInnerBoundary_TransparentTurbulence( Direction )

  END SUBROUTINE InitializeFields_TransparentTurbulence


  SUBROUTINE SetInnerBoundary_TransparentTurbulence( Direction )

    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeE
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: loX(3), hiX(3)
    INTEGER  :: loZ(4), hiZ(4)
    REAL(DP) :: E

    ! --- Fluid Fields ---

    loX = iX_B0
    hiX = iX_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loX(1) = iX_B1(1)
      hiX(1) = iX_B1(1)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loX(2) = iX_B1(2)
      hiX(2) = iX_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loX(3) = iX_B1(3)
      hiX(3) = iX_B1(3)
    END IF

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1), hiX(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    loZ = iZ_B0
    hiZ = iZ_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loZ(2) = iZ_B1(2)
      hiZ(2) = iZ_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loZ(3) = iZ_B1(3)
      hiZ(3) = iZ_B1(3)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loZ(4) = iZ_B1(4)
      hiZ(4) = iZ_B1(4)
    END IF

    DO iS  = 1, nSpecies
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2), hiZ(2)
    DO iZ1 = loZ(1), hiZ(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = One / ( EXP( E / Three - Three ) + One )

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

        END IF

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE SetInnerBoundary_TransparentTurbulence


  SUBROUTINE InitializeFields_TransparentShock( V_0, ShockWidth, Direction )

    REAL(DP),     INTENT(in) :: V_0(3)
    REAL(DP),     INTENT(in) :: ShockWidth
    CHARACTER(2), INTENT(in) :: Direction

    REAL(DP), PARAMETER :: X_Shock = 1.0d0

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3

    WRITE(*,*)
    WRITE(*,'(A6,A13,3ES9.2E2)') '', 'V_0 = ', V_0
    WRITE(*,'(A6,A13,1ES9.2E2)') '', 'ShockWidth = ', ShockWidth
    WRITE(*,'(A6,A13,A2)'      ) '', 'Direction = ', TRIM( Direction )
    WRITE(*,*)

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One

        SELECT CASE( TRIM( Direction ) )
        CASE( 'X' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = Half * V_0(1) * ( One + TANH( (X1-X_Shock)/ShockWidth ) )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        CASE( 'Y' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = Half * V_0(2) * ( One + TANH( (X2-X_Shock)/ShockWidth ) )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        CASE( 'Z' )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
            = Half * V_0(3) * ( One + TANH( (X3-X_Shock)/ShockWidth ) )

        END SELECT

        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0d-0


        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 1.0d-8
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

    CALL SetInnerBoundary_TransparentShock( Direction )

  END SUBROUTINE InitializeFields_TransparentShock


  SUBROUTINE SetInnerBoundary_TransparentShock( Direction )

    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeE
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: loX(3), hiX(3)
    INTEGER  :: loZ(4), hiZ(4)
    REAL(DP) :: E

    ! --- Fluid Fields ---

    loX = iX_B0
    hiX = iX_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loX(1) = iX_B1(1)
      hiX(1) = iX_B1(1)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loX(2) = iX_B1(2)
      hiX(2) = iX_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loX(3) = iX_B1(3)
      hiX(3) = iX_B1(3)
    END IF

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1), hiX(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    loZ = iZ_B0
    hiZ = iZ_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loZ(2) = iZ_B1(2)
      hiZ(2) = iZ_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loZ(3) = iZ_B1(3)
      hiZ(3) = iZ_B1(3)
    ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
      loZ(4) = iZ_B1(4)
      hiZ(4) = iZ_B1(4)
    END IF

    DO iS  = 1, nSpecies
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2), hiZ(2)
    DO iZ1 = loZ(1), hiZ(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = One / ( EXP( E / Three - Three ) + One )

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

        END IF

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE SetInnerBoundary_TransparentShock


  SUBROUTINE InitializeFields_TransparentVortex( V_0, Direction )

    REAL(DP),     INTENT(in) :: V_0(3)
    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, Beta, R

    Beta = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        R  = SQRT( X1**2 + X2**2 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
          = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = - X2 * Beta * EXP( Half * ( One - R**2 ) )
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = + X1 * Beta * EXP( Half * ( One - R**2 ) )
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) &
          = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) &
          = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 1.0d-8
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

    CALL SetInnerBoundary_TransparentVortex( Direction )

  END SUBROUTINE InitializeFields_TransparentVortex


  SUBROUTINE SetInnerBoundary_TransparentVortex( Direction )

    CHARACTER(2), INTENT(in) :: Direction

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeE
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: loX(3), hiX(3)
    INTEGER  :: loZ(4), hiZ(4)
    REAL(DP) :: E, Mu_0

    ! --- Fluid Fields ---

    loX = iX_B0
    hiX = iX_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loX(1) = iX_B1(1)
      hiX(1) = iX_B1(1)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loX(2) = iX_B1(2)
      hiX(2) = iX_B1(2)
    END IF

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1), hiX(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    loZ = iZ_B0
    hiZ = iZ_E0
    IF(     TRIM( Direction ) .EQ. 'X' )THEN
      loZ(2) = iZ_B1(2)
      hiZ(2) = iZ_B1(2)
    ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN
      loZ(3) = iZ_B1(3)
      hiZ(3) = iZ_B1(3)
    END IF

    DO iS  = 1, nSpecies
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2), hiZ(2)
    DO iZ1 = loZ(1), hiZ(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        Mu_0 = 0.9_DP

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.5_DP * ( One + Mu_0 ) * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
            = 0.0_DP
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
            = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
            = 0.0_DP

        END IF

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE SetInnerBoundary_TransparentVortex


  SUBROUTINE InitializeFields_RadiatingSphere

    CHARACTER(32) :: Profile = 'Shock'
    INTEGER       :: iNodeX, iX1, iX2, iX3
    INTEGER       :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER       :: iNodeX1
    REAL(DP)      :: X1, Theta, V_Max = 0.30_DP

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = One

        SELECT CASE( TRIM( Profile ) )

        CASE( 'Shock' )

          IF( X1 <= 135.0_DP )THEN

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = Zero

          ELSEIF( X1 > 135.0_DP .AND. X1 <= 150.0_DP )THEN

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = - V_Max * ( X1 - 135.0_DP ) / 15.0_DP

          ELSEIF( X1 > 150.0_DP )THEN

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = - V_Max * ( 150.0_DP / X1 )**2

          END IF

        CASE( 'Collapse' )

          Theta = Half * ( One + TANH( ( X1 - 2.0d2 ) / 3.0d1 ) )

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = - ( One - Theta ) * V_Max * ( X1 - 1.0d1 ) / ( 2.0d2 - 1.0d1 ) &
              - Theta * V_Max * ( 2.0d2 / X1 )**2

        END SELECT

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0d-0

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 1.0d-40
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

    CALL SetInnerBoundary_RadiatingSphere

  END SUBROUTINE InitializeFields_RadiatingSphere


  SUBROUTINE SetInnerBoundary_RadiatingSphere

    INTEGER  :: iNodeZ, iNodeX, iNodeE
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: E

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B1(2), iZ_B1(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = One / ( EXP( E / Three - Three ) + One )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = 0.999_DP * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE SetInnerBoundary_RadiatingSphere


  SUBROUTINE InitializeFields_GaussianDiffusion( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeZ2, iNodeZ3
    REAL(DP) :: X1, X2, t_0, D_min, D_0, X1_0, X2_0

    D_min = 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One
    X2_0  = One

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

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3 )

        D_0 = One / ( Three * uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = EXP( - ( (X1-X1_0)**2 + (X2-X2_0)**2 ) / ( 4.0_DP * t_0 * D_0 ) )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = Half * (X1-X1_0) * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) / t_0
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = Half * (X2-X2_0) * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) / t_0
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        IF( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) < D_min )THEN

          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = D_min
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0d-00
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0d-00
          uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0d-00

        END IF

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE InitializeFields_GaussianDiffusion


  SUBROUTINE InitializeFields_ExpandingAtmosphere( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    REAL(DP), PARAMETER :: Rmin = 01.0_DP
    REAL(DP), PARAMETER :: Rmax = 11.0_DP

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeE
    REAL(DP) :: R

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP

        IF( R .LT. Rmin )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP

        ELSEIF( R .GE. Rmin .AND. R .LE. Rmax )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = V_0(1) * ( R - Rmin ) / ( Rmax - Rmin )

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)

        END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 1.0d-10
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = Zero
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )
      
         END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_ExpandingAtmosphere


  SUBROUTINE InitializeFields_HomogeneousSphere1D

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = 1.0d-8
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE InitializeFields_HomogeneousSphere1D


  SUBROUTINE InitializeFields_HomogeneousSphere2D( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    REAL(DP), PARAMETER :: R_Sh = 2.0_DP
    REAL(DP), PARAMETER :: L_R  = 1.0d-2

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, R, absV, ShapeFunction

    absV = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        R = SQRT( X1**2 + X2**2 + X3**2 )

        ShapeFunction &
          = Half * ( One + tanh( ( R - R_Sh ) / L_R ) ) * sqrt( R_Sh / R )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = - absV * ShapeFunction * X1 / R
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = - absV * ShapeFunction * X2 / R
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = - absV * ShapeFunction * X3 / R
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

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

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 1.0d-8
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
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

  END SUBROUTINE InitializeFields_HomogeneousSphere2D


  SUBROUTINE ComputeError( t )

    REAL(DP), INTENT(in) :: t

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL ComputeError_SineWaveStreaming( t )

      CASE DEFAULT

    END SELECT

  END SUBROUTINE ComputeError


  SUBROUTINE ComputeError_SineWaveStreaming( t )

    REAL(DP), INTENT(in) :: t

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeZ, iNodeX1
    REAL(DP) :: X1, D_A, I1_A, I2_A, I3_A
    REAL(DP) :: MaxError(nPR,nSpecies)

    MaxError = Zero

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeZ = 1, nDOFZ

        iNodeX1 = NodeNumberTableZ(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        D_A  = 0.50_DP + 0.49_DP * SIN( TwoPi * ( X1 - t ) )
        I1_A = D_A
        I2_A = Zero
        I3_A = Zero

        MaxError(iPR_D, iS) &
          = MAX( ABS( D_A  - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) ), &
                 MaxError(iPR_D ,iS) )

        MaxError(iPR_I1,iS) &
          = MAX( ABS( I1_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) ), &
                 MaxError(iPR_I1,iS) )

        MaxError(iPR_I2,iS) &
          = MAX( ABS( I2_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) ), &
                 MaxError(iPR_I2,iS) )

        MaxError(iPR_I3,iS) &
          = MAX( ABS( I3_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) ), &
                 MaxError(iPR_I3,iS) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,*)
    WRITE(*,'(A2,A)') '', 'INFO: SineWaveStreaming Error'
    WRITE(*,*)
    WRITE(*,'(A4,A2,4A12)') '', 'Sp', 'N', 'G1', 'G2', 'G3'
    DO iS = 1, nSpecies
      WRITE(*,'(A4,I2.2,4ES12.4E2)') '', iS, MaxError(:,iS)
    END DO

  END SUBROUTINE ComputeError_SineWaveStreaming


END MODULE InitializationModule
