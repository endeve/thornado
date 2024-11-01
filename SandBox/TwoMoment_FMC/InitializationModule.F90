MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Two, Three, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFX, nDOFE, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    xL, xR
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    uPM, iPM_J, iPM_H1, iPM_H2, iPM_H3, nPM, &
    uCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, nCM
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    ComputeConserved_TwoMoment_FMC
  USE TwoMoment_OpacityModule_FMC, ONLY: &
    uOP, iOP_Sigma

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields( V_0 , LengthScale, Direction, Spectrum )

    REAL(DP),      INTENT(in) :: V_0(3)
    REAL(DP),      INTENT(in) :: LengthScale
    CHARACTER(2),  INTENT(in) :: Direction
    CHARACTER(32), INTENT(in) :: Spectrum


    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( V_0 )

      CASE( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion( V_0 )

      CASE( 'StreamingDopplerShift')

        CALL InitializeFields_StreamingDopplerShift( V_0, Direction, Spectrum )

      CASE ( 'TransparentShock' )

        CALL InitializeFields_TransparentShock &
               ( V_0, LengthScale, Direction )

      CASE( 'TransparentVortex' )

        CALL InitializeFields_TransparentVortex &
               ( V_0, Direction )

      CASE( 'GaussianDiffusion1D' )

        CALL InitializeFields_GaussianDiffusion1D( V_0 )

      CASE( 'DiffusionMovingMed')

        CALL InitializeFields_DiffusionMovingMed ( V_0 )

      CASE( 'ShadowCasting2D_Cartesian' )

        CALL InitializeFields_ShadowCasting2D_Cartesian

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uCF, uPF, uCM, uPM )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uCF, uPF, uCM, uPM )
#endif

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeZ2, iNodeZ3, iNodeZ4
    Real(DP) :: X1, X2, X3, W, vMagSq

    vMagSq = V_0(1)**2 + V_0(2)**2 + V_0(3)**2
    W = One / SQRT( One - vMagSq )

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

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE,nDOFX ) + 1

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)
        iNodeZ4 = NodeNumberTableZ(4,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2)
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3)
        X3 = NodeCoordinate( MeshX(3), iZ4, iNodeZ4)

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) &
          = 0.50_DP + 0.49_DP * SIN( TwoPi * X1)
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = W * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS)
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
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

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

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

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) &
          = 0.49_DP * SIN( Third * Pi * X1 ) + 0.5_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = - ( 0.49_DP * Pi / ( 9.0_DP * Sigma ) ) * COS( Third * Pi * X1 )
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  SUBROUTINE InitializeFields_StreamingDopplerShift( V_0 , Direction, Spectrum )

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
              ! = 2 * V_0(1) / 3 * ( X1 - X_0 )
              = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1)
          ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              ! = - 2 * V_0(1) / 3 * ( X1 - X_3 )
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

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---
    
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) &
          = 1.0d-40
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = Zero
        
        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

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

            uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) &
              = E / ( EXP( E / Three - Three ) + One )

          CASE( 'Bose-Einstein' )

            uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) &
              = One / ( EXP( E ) - One )

          CASE DEFAULT

            uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) &
              = E / ( EXP( E / Three - Three ) + One )

          END SELECT

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = ( One - 1.0d-3 ) * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.999_DP * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.999_DP * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)

        END IF

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  SUBROUTINE InitializeFields_TransparentShock( V_0, ShockWidth, Direction)

    REAL(DP),     INTENT(in) :: V_0(3)
    REAL(DP),     INTENT(in) :: ShockWidth
    CHARACTER(2), INTENT(in) :: Direction

    REAL(DP), PARAMETER :: X_Shock = 1.0d0
    REAL(DP), PARAMETER :: delta   = 0.00000001_DP ! Set as 0.0_DP if using analytic IC, otherwise use 0.001_DP

    INTEGER  :: iNodeX, iX1, iX2, iX3, i, iXi
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeXi
    INTEGER  :: iNodeE
    REAL(DP) :: X1, X2, X3, Xi, X1_C
    REAL(DP) :: E, sSq, vMagSq, W
    REAL(DP) :: SUM_K, SUM_KPF
    REAL(DP) :: uPF_smooth(nDOFX, &
                           iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))

    WRITE(*,*)
    WRITE(*,'(A6,A13,3ES9.2E2)') '', 'V_0 = ', V_0
    WRITE(*,'(A6,A13,1ES9.2E2)') '', 'ShockWidth = ', ShockWidth
    WRITE(*,'(A6,A13,A2)'      ) '', 'Direction = ', TRIM( Direction )
    WRITE(*,*)

    ! --- Fluid Fields ---

    ! --- Unsmoothed velocity ---

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

      END DO

    END DO
    END DO
    END DO

    ! ! --- Smooth velocity ---

    ! ASSOCIATE( dX1 => MeshX(1) % Width )

    ! DO iX3 = iX_B0(3), iX_E0(3)
    ! DO iX2 = iX_B0(2), iX_E0(2)
    ! DO iX1 = iX_B0(1), iX_E0(1)
  
    !     DO iNodeX = 1, nDOFX

    !       iNodeX1 = NodeNumberTableX(1,iNodeX)
    !       X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
    !       X1_C = MeshX(1) % Center(iX1)

    !       SUM_K = Zero
    !       SUM_KPF = Zero

    !       DO iXi = iX_B0(1), iX_E0(1)
    !       DO i = 1, nDOFX

    !         iNodeXi = NodeNumberTableX(1,i)

    !         Xi = NodeCoordinate( MeshX(1), iXi, iNodeXi )

    !         IF( ABS(Xi-X1_C)<3.0_DP * dX1(iX1) / 2.0_DP )THEN

    !           SUM_K = SUM_K + &
    !                   EXP(-((X1-Xi)**2)/(Two*(dX1(iX1)/4.0_DP)**2))

    !           SUM_KPF = SUM_KPF + &
    !                     EXP(-((X1-Xi)**2)/(Two*(dX1(iX1)/4.0_DP)**2)) * uPF(iNodeXi,iXi,iX2,iX3,iPF_V1)

    !         END IF

    !       END DO
    !       END DO

    !       uPF_smooth(iNodeX,iX1,iX2,iX3) = SUM_KPF / SUM_K

    !     END DO

    ! END DO
    ! END DO
    ! END DO

    ! END ASSOCIATE !dX1

    ! ! --- Replace velocity by smoothed velocity ---

    ! DO iX3 = iX_B0(3), iX_E0(3)
    ! DO iX2 = iX_B0(2), iX_E0(2)
    ! DO iX1 = iX_B0(1), iX_E0(1)
  
    !     DO iNodeX = 1, nDOFX

    !       uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = uPF_smooth(iNodeX,iX1,iX2,iX3)

    !     END DO

    ! END DO
    ! END DO
    ! END DO

    ! --- Two-Moment Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        sSq = (One + uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)) / &
              (One - uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1))

        vMagSq = uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)**2 + &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2)**2 + &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3)**2

        W = One / SQRT( One - vMagSq )

        ! uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = 1.0d-8
        ! uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) = 0.0_DP
        ! uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) = 0.0_DP
        ! uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) = 0.0_DP
        ! Replace J by steady state with a delta parameter

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = 1.0d-8
        ! uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = sSq * E / ( EXP( SQRT(sSQ) * E / Three - Three ) + One )
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) = (One - delta) * W * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

    REAL(DP), PARAMETER :: delta   = 0.00000001_DP ! Set as 0.0_DP if using analytic IC, otherwise use 0.001_DP

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

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

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

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) &
          ! = One / ( EXP( E / Three - Three ) + One )
          = E / ( EXP( E / Three - Three ) + One )

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = (One - delta) * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.999_DP * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.999_DP * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS)

        END IF

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = 1.0d-8
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) = 0.0_DP
        
        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

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

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) &
          = E * 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )

        IF(     TRIM( Direction ) .EQ. 'X' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.5_DP * ( One + Mu_0 ) * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
            = 0.0_DP
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
            = 0.5_DP * ( One + Mu_0 ) * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS)
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
            = 0.0_DP

        END IF

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  SUBROUTINE InitializeFields_GaussianDiffusion1D( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeZ2, iNodeZ3
    REAL(DP) :: X1, t_0, J_min, J_0, X1_0

    J_min = 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One

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

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

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

        J_0 = One / ( Three * uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) )

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) &
          = EXP( - ( (X1-X1_0)**2 ) / ( 4.0_DP * t_0 * J_0 ) )
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = Half * (X1-X1_0) * uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) / t_0
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = 0.0_DP

        IF( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) < J_min )THEN

          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = J_min
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) = 0.0d-00
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) = 0.0d-00
          uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) = 0.0d-00

        END IF

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  END SUBROUTINE InitializeFields_GaussianDiffusion1D

  SUBROUTINE InitializeFields_DiffusionMovingMed( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeZ2
    REAL(DP) :: X1, WSq, vMagSq

    vMagSq = V_0(1)**2 + V_0(2)**2 + V_0(3)**2
    WSq = One / ( One - vMagSq )

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

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J, iS) &
          = 3.0_DP * EXP(-9.0_DP * X1**2) / ( 4.0_DP * WSq - One )
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = Zero

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  END SUBROUTINE InitializeFields_DiffusionMovingMed

  SUBROUTINE InitializeFields_ShadowCasting2D_Cartesian

    INTEGER :: iNodeX, iX1, iX2, iX3
    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---
    
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS) = 1.0d-10
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) = Zero
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) = Zero
        
        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
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

  END SUBROUTINE InitializeFields_ShadowCasting2D_Cartesian

END MODULE InitializationModule
