MODULE InitializationModule

  USE KindModule, ONLY: &
    DP,       &
    Zero,     &
    Half,     &
    One,      &
    Two,      &
    Three,    &
    Four,     &
    Pi,       &
    TwoPi,    &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX,     &
    nDimsX,      &
    nDOFX,       &
    iX_B0,       &
    iX_B1,       &
    iX_E0,       &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    uGF,          &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_Psi
  USE MagnetofluidFieldsModule, ONLY: &
    nPM,     &
    uPM,     &
    iPM_D,   &
    iPM_V1,  &
    iPM_V2,  &
    iPM_V3,  &
    iPM_E,   &
    iPM_Ne,  &
    iPM_B1,  &
    iPM_B2,  &
    iPM_B3,  &
    iPM_Chi, &
    uCM,     &
    iCM_D,   &
    iCM_S1,  &
    iCM_S2,  &
    iCM_S3,  &
    iCM_E,   &
    iCM_Ne,  &
    iCM_B1,  &
    iCM_B2,  &
    iCM_B3,  &
    iCM_Chi, &
    uAM,     &
    iAM_P
  USE EquationOfStateModule, ONLY : &
      ComputePressureFromPrimitive
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Second, &
    Erg, &
    Gauss
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: HDFERR

  PUBLIC :: InitializeFields_Relativistic_MHD


CONTAINS


  SUBROUTINE InitializeFields_Relativistic_MHD &
               ( AdvectionProfile_Option, SmoothProfile_Option, &
                 ConstantDensity_Option, &
                 Angle_Option, RiemannProblemName_Option, &
                 MMBlastWaveB0_Option, MMBlastWavePhi_Option, &
                 OTScaleFactor_Option, EvolveOnlyMagnetic_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: EvolveOnlyMagnetic_Option
    LOGICAL,          INTENT(in), OPTIONAL :: SmoothProfile_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ConstantDensity_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Angle_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MMBlastWaveB0_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MMBlastWavePhi_Option
    REAL(DP),         INTENT(in), OPTIONAL :: OTScaleFactor_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'HydroSineWaveX1'
    CHARACTER(LEN=64) :: RiemannProblemName = 'IsolatedContact'
    LOGICAL           :: SmoothProfile = .TRUE.
    LOGICAL           :: ConstantDensity = .TRUE.
    REAL(DP)          :: Angle = Pi / Four
    REAL(DP)          :: MMBlastWaveB0 = 0.5_DP
    REAL(DP)          :: MMBlastWavePhi = 0.0_DP
    REAL(DP)          :: OTScaleFactor = 100

    LOGICAL           :: EvolveOnlyMagnetic = .FALSE.

    uPM(:,:,:,:,iPM_Ne) = Zero

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( SmoothProfile_Option ) ) &
      SmoothProfile = SmoothProfile_Option

    IF( PRESENT( ConstantDensity_Option ) ) &
      ConstantDensity = ConstantDensity_Option

    IF( PRESENT( Angle_Option ) ) &
      Angle = Angle_Option

    IF( PRESENT( RiemannProblemName_Option ) ) &
      RiemannProblemName = TRIM( RiemannProblemName_Option )

    IF( PRESENT( MMBlastWaveB0_Option ) ) &
      MMBlastWaveB0 = MMBlastWaveB0_Option

    IF( PRESENT( MMBlastWavePhi_Option ) ) &
      MMBlastWavePhi = MMBlastWavePhi_Option

    IF( PRESENT( OTScaleFactor_Option ) ) &
      OTScaleFactor = OTScaleFactor_Option

    IF( PRESENT( EvolveOnlyMagnetic_Option ) ) &
      EvolveOnlyMagnetic = EvolveOnlyMagnetic_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection1D' )

        CALL InitializeFields_Advection1D &
               ( TRIM( AdvectionProfile ), EvolveOnlyMagnetic )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D &
               ( TRIM( AdvectionProfile ), Angle, EvolveOnlyMagnetic )

      CASE( 'Advection3D' )

        CALL InitializeFields_Advection3D &
               ( TRIM( AdvectionProfile ), EvolveOnlyMagnetic )

      CASE( 'Cleaning1D' )

        CALL InitializeFields_Cleaning1D &
               ( SmoothProfile, EvolveOnlyMagnetic )

      CASE( 'Cleaning2D' )

        CALL InitializeFields_Cleaning2D &
               ( ConstantDensity, EvolveOnlyMagnetic )

      CASE( 'Riemann1D' )

        CALL InitializeFields_Riemann1D &
               ( TRIM( RiemannProblemName ), EvolveOnlyMagnetic )

      CASE( 'MMBlastWave2D' )

        CALL InitializeFields_MMBlastWave2D &
               ( MMBlastWaveB0, MMBlastWavePhi, EvolveOnlyMagnetic )

      CASE( 'OrszagTang2D' )

        CALL InitializeFields_OrszagTang2D( OTScaleFactor, EvolveOnlyMagnetic )

      CASE( 'ShearingDisk' )

        CALL InitializeFields_ShearingDisk( EvolveOnlyMagnetic )

      CASE( 'MagneticKH' )

        CALL InitializeFields_MagneticKH( EvolveOnlyMagnetic )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic_MHD


  SUBROUTINE InitializeFields_Advection1D( AdvectionProfile, EvolveOnlyMagnetic )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile
    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1
    REAL(DP) :: Eta, h, P, VA, W, k, L, VdotB, Kappa, &
                V1_Transport, V2_Transport, V3_Transport

    REAL(DP) :: V1, V2, V3, CB1, CB2, CB3

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          ! Sine wave advection problem to test hydro portion of code.

          CASE( 'HydroSineWaveX1' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One + 0.1_DP * SIN ( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE( 'MagneticSineWaveX1' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP * SIN( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * COS( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          ! Variant of Top Hat advection problem from Evans and
          ! Hawley (1988).

          CASE( 'TopHat' )

          ! Circularly polarized Alfven wave with the exact solution
          ! from Del Zanna et al. (2007) and Mattia and Mignone (2022).

          CASE( 'CPAlfvenX1' )

            Eta = One
            k   = One
            h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
            VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                        * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

            W = One / SQRT( One - VA**2 * Eta**2 )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = Zero
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = -VA * Eta * COS( k * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = -VA * Eta * SIN( k * X1 )
            uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            VdotB = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) * Eta * COS( k * X1 ) &
                      + uPM(iNodeX,iX1,iX2,iX3,iPM_V3) * Eta * SIN( k * X1 )

            V1_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V2_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V3_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V3) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
              = W * VdotB * V1_Transport + ( One / W )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) &
              = W * VdotB * V2_Transport + Eta * COS( k * X1 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) &
              = W * VdotB * V3_Transport + Eta * SIN( k * X1 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  HydroSineWaveX1'
            WRITE(*,'(A)') '  MagneticSineWaveX1'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,'(A)') '  CPAlfvenX1'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection1D


  SUBROUTINE InitializeFields_Advection2D( AdvectionProfile, Angle, EvolveOnlyMagnetic )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile
    LOGICAL,  INTENT(in) :: EvolveOnlyMagnetic
    REAL(DP), INTENT(in)         :: Angle

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: Eta, h, P, VA, W, kappa, L, k, VdotB, &
                V1_Transport, V2_Transport, V3_Transport
    REAL(DP) :: V1, V2, V3, VSq, CB1, CB2, CB3

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          ! Sine wave advection problems to test hydro portion of code.

         CASE( 'HydroSineWaveX2' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D) &
              = One + 0.1_DP * SIN ( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            CALL ComputePressureFromPrimitive &
                   ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
                     uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
                     0.0_DP, P )

          CASE( 'HydroSineWaveX1X2' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D) &
              = One + 0.1_DP * SIN ( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

         CASE( 'MagneticSineWaveX2' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP * SIN( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * COS( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE( 'MagneticSineWaveX1X2' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          ! Circularly polarized Alfven wave based
          ! on the exact solution from, e.g.,
          ! Del Zanna et al. (2007) and Mattia and Mignone (2022).

          CASE( 'CPAlfvenX2' )

            Eta = One
            k   = One
            h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
            VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                        * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

            W = One / SQRT( One - VA**2 * Eta**2 )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = -VA * Eta * COS( k * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = -VA * Eta * SIN( k * X2 )
            uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            VdotB = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) * Eta * COS( k * X2 ) &
                      + uPM(iNodeX,iX1,iX2,iX3,iPM_V3) * Eta * SIN( k * X2 )

            V1_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V2_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V3_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V3) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
              = W * VdotB * V1_Transport + Eta * COS( k * X2 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) &
              = W * VdotB * V2_Transport + ( One / W )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) &
              = W * VdotB * V3_Transport + Eta * SIN( k * X2 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          ! See oblique case conditions from Londrillo
          ! and Del Zanna (2004) (with sin and cos switched).

          CASE( 'CPAlfvenOblique' )

            L = X1 * COS( Angle ) + X2 * SIN( Angle )

            Eta = One
            k   = One
            h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
            VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                        * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

            W = One / SQRT( One - VA**2 * Eta**2 )

            V1 =  VA * Eta * SIN( Angle ) * COS( TwoPi * k * L )
            V2 = -VA * Eta * COS( Angle ) * COS( TwoPi * k * L )
            V3 = -VA * Eta * SIN( TwoPi * k * L )

            CB1 = COS( Angle ) + V1 / ( -VA )
            CB2 = SIN( Angle ) + V2 / ( -VA )
            CB3 = V3 / ( -VA )

            VSq = V1**2 + V2**2 + V3**2

            VdotB = V1 * CB1 + V2 * CB2 + V3 * CB3

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = V1
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = V2
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = V3
            uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
              = W * VdotB * V1 + CB1 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) &
              = W * VdotB * V2 + CB2 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) &
              = W * VdotB * V3 + CB3 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          ! Loop advection problem from Section 5.5 of Mosta et al. (2014).

          CASE( 'LoopAdvection' )

            V1 = One / Two
            V2 = Zero
            V3 = Zero

            W = One / SQRT( One - V1**2  - V2**2 - V3**2 )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = V1
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = V2
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = V3

            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = Three
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
                = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            IF( SQRT( X1**2 + X2**2 ) .LE. 0.3_DP )THEN

              CB1 = -1.0d-3 * ( X2 / SQRT( X1**2 + X2**2) )
              CB2 =  1.0d-3 * ( X1 / SQRT( X1**2 + X2**2) )
              CB3 = Zero

            ELSE

              CB1 = Zero
              CB2 = Zero
              CB3 = Zero

            END IF

            VdotB = V1 * CB1 + V2 * CB2 + V3 * CB3

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  HydroSineWaveX2'
            WRITE(*,'(A)') '  HydroSineWaveX1X2'
            WRITE(*,'(A)') '  MagneticSineWaveX2'
            WRITE(*,'(A)') '  MagneticSineWaveX1X2'
            WRITE(*,'(A)') '  CPAlfvenX2'
            WRITE(*,'(A)') '  CPAlfvenOblique'
            WRITE(*,'(A)') '  LoopAdvection'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_Advection3D( AdvectionProfile, EvolveOnlyMagnetic )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile
    LOGICAL,  INTENT(in) :: EvolveOnlyMagnetic

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3
    REAL(DP) :: Eta, h, P, VA, W, k, VdotB, &
                V1_Transport, V2_Transport, V3_Transport

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

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

        SELECT CASE( TRIM( AdvectionProfile ) )

          ! Sine wave advection problems to test hydro portion of code.

         CASE( 'HydroSineWaveX3' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D) &
              = One + 0.1_DP * SIN ( TwoPi * X3 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.1_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            CALL ComputePressureFromPrimitive &
                   ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
                     uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
                     0.0_DP, P )

          CASE( 'MagneticSineWaveX3' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.1_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP * SIN( TwoPi * X3 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP * COS( TwoPi * X3 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE( 'CPAlfvenX3' )

            Eta = One
            k   = One
            h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
            VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                        * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

            W = One / SQRT( One - VA**2 * Eta**2 )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = -VA * Eta * COS( k * X3 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = -VA * Eta * SIN( k * X3 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            VdotB = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) * Eta * COS( k * X3 ) &
                      + uPM(iNodeX,iX1,iX2,iX3,iPM_V2) * Eta * SIN( k * X3 )

            V1_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V2_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            V3_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V3) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) )

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
              = W * VdotB * V1_Transport + Eta * COS( k * X3 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) &
              = W * VdotB * V2_Transport + Eta * SIN( k * X3 ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) &
              = W * VdotB * V3_Transport + ( One / W )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  HydroSineWaveX3'
            WRITE(*,'(A)') '  MagneticSineWaveX3'
            WRITE(*,'(A)') '  CPAlfvenX3'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection3D


  SUBROUTINE InitializeFields_Cleaning1D( SmoothProfile, EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: SmoothProfile
    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    ! 1D divergence cleaning test from Section 5.1 of Derigs et al. (2018)
    ! with option to use only the smooth part of the initial condition.

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1
    REAL(DP) :: V1, V2, V3, W, CB1, CB2, CB3, VdotB

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        V1 = Zero
        V2 = Zero
        V3 = Zero

        W = One / SQRT( One - V1**2 - V2**2 - V3**2 )

        uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
        uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = V1
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = V2
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = V3
        uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
        uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
          = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

        IF( SmoothProfile )THEN

          IF( ( X1 >= -One ) .AND. ( X1 <= -0.6_DP ) )THEN
            CB1 = Zero
          ELSE IF( ( X1 >= 0.6_DP ) .AND. ( X1 <= One ) )THEN
            CB1 = Zero
          ELSE
            CB1 = EXP( -( X1 / 0.11_DP )**2 / Two )
          END IF

        ELSE

          IF( ( X1 > -0.8_DP ) .AND. ( X1 <= -0.6_DP ) )THEN
            CB1 = -Two * ( X1 + 0.8_DP )
          ELSE IF( ( X1 > -0.6_DP ) .AND. ( X1 <= 0.6_DP ) )THEN
            CB1 = EXP( -( X1 / 0.11_DP )**2 / Two )
          ELSE IF( ( X1 > 0.6_DP ) .AND. ( X1 < One ) )THEN
            CB1 = Half
          ELSE
            CB1 = Zero
          END IF

        END IF

        CB2 = Zero
        CB3 = Zero

        VdotB = V1 * CB1 + V2 * CB2 + V3 * CB3

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Cleaning1D


  SUBROUTINE InitializeFields_Cleaning2D( ConstantDensity, EvolveOnlyMagnetic )

    ! 2D divergence cleaning test from Section 5.2 of
    ! Derigs et al. (2018) with option to use
    ! constant density.

    LOGICAL, INTENT(in) :: ConstantDensity
    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    REAL(DP) :: D, V1, V2, V3, W
    REAL(DP) :: CB1, CB2, CB3, VdotB
    REAL(DP) :: R, R0

    R0 = One / SQRT( 8.0_DP )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        R = SQRT( X1**2 + X2**2 )

        IF( ConstantDensity )THEN

          D = One

        ELSE

          IF( X1 .LE. 0.5_DP )THEN
            D = One
          ELSE
            D = Two
          END IF

        END IF

        V1 = Zero
        V2 = Zero
        V3 = Zero

        W = One / SQRT( One - V1**2 - V2**2 - V3**2 )

        uPM(iNodeX,iX1,iX2,iX3,iPM_D )  = D
        uPM(iNodeX,iX1,iX2,iX3,iPM_V1 ) = V1
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2 ) = V2
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3 ) = V3
        uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 6.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
          = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

        IF( R .LE. R0 )THEN

          CB1 = ( One / SQRT( FourPi ) ) &
                * ( ( R / R0 )**8 - Two * ( R / R0 )**4 + One )

        ELSE

          CB1 = Zero

        END IF

        CB2 = Zero
        CB3 = One / SQRT( FourPi )

        VdotB = V1 * CB1 + V2 * CB2 + V3 * CB3

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = W * VdotB * V1 + CB1 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = W * VdotB * V2 + CB2 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = W * VdotB * V3 + CB3 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Cleaning2D


  SUBROUTINE InitializeFields_Riemann1D &
               ( RiemannProblemName, EvolveOnlyMagnetic )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName
    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    INTEGER iX1, iX2, iX3
    INTEGER iNodeX, iNodeX1
    REAL(DP) :: X1, XD
    REAL(DP) :: CB1_L, CB2_L, CB3_L, &
                CB1_R, CB2_R, CB3_R
    REAL(DP) :: VSq_L, VSq_R, W_L, W_R, &
                VdotB_L, VdotB_R

    REAL(DP) :: LeftState(nPM), RightState(nPM)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', '1D Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'HydroIsolatedContact' )

        ! --- Pure hydro version of isolated contact ---
        ! --- from Mattia & Mignone, 2022, MRNAS,    ---
        ! --- 510, 481-499, Table 1                  ---

        XD = Half

        LeftState(iPM_D  ) = 5.9718209694880811e0_DP
        LeftState(iPM_V1 ) = 0.01_DP
        LeftState(iPM_V2 ) = 0.0_DP
        LeftState(iPM_V3 ) = 0.0_DP
        LeftState(iPM_E  ) = One / ( Gamma_IDEAL - One )

        LeftState(iPM_B1 ) = 0.0_DP
        LeftState(iPM_B2 ) = 0.0_DP
        LeftState(iPM_B3 ) = 0.0_DP
        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D  ) = One
        RightState(iPM_V1 ) = 0.01_DP
        RightState(iPM_V2 ) = 0.0_DP
        RightState(iPM_V3 ) = 0.0_DP
        RightState(iPM_E  ) = One / ( Gamma_IDEAL - One )

        RightState(iPM_B1 ) = 0.0_DP
        RightState(iPM_B2 ) = 0.0_DP
        RightState(iPM_B3 ) = 0.0_DP
        RightState(iPM_Chi) = 0.0_DP

     CASE( 'IsolatedContact' )

        ! --- Isolated contact from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                                ---

        XD = Half

        LeftState(iPM_D ) = 10.0_DP
        LeftState(iPM_V1) = 0.0_DP
        LeftState(iPM_V2) = 0.7_DP
        LeftState(iPM_V3) = 0.2_DP
        LeftState(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CB1_L = 5.0_DP
        CB2_L = One
        CB3_L = Half

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) = 1.0_DP
        RightState(iPM_V1) = 0.0_DP
        RightState(iPM_V2) = 0.7_DP
        RightState(iPM_V3) = 0.2_DP
        RightState(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CB1_R = 5.0_DP
        CB2_R = One
        CB3_R = Half

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'RotationalWave' )

        ! --- Rotational wave from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                               ---

        XD = Half

        LeftState(iPM_D ) =  1.0_DP
        LeftState(iPM_V1) =  0.4_DP
        LeftState(iPM_V2) = -0.3_DP
        LeftState(iPM_V3) =  0.5_DP
        LeftState(iPM_E ) =  1.0_DP / ( Gamma_IDEAL - One )

        CB1_L =  2.4_DP
        CB2_L =  One
        CB3_L = -1.6_DP

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) =  1.0_DP
        RightState(iPM_V1) =  0.377237_DP
        RightState(iPM_V2) = -0.482389_DP
        RightState(iPM_V3) =  0.424190_DP
        RightState(iPM_E ) =  1.0_DP / ( Gamma_IDEAL - One )

        CB1_R =  2.4_DP
        CB2_R = -0.1_DP
        CB3_R = -2.178213_DP

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'HydroShockTube1' )

        ! --- Pure hydro version of ST1 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                                         ---

        XD = Half

        LeftState(iPM_D  ) = 1.0_DP
        LeftState(iPM_V1 ) = 0.0_DP
        LeftState(iPM_V2 ) = 0.0_DP
        LeftState(iPM_V3 ) = 0.0_DP
        LeftState(iPM_E  ) = 1.0_DP / ( Gamma_IDEAL - One )

        LeftState(iPM_B1 ) = 0.0_DP
        LeftState(iPM_B2 ) = 0.0_DP
        LeftState(iPM_B3 ) = 0.0_DP
        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D  ) = 0.125_DP
        RightState(iPM_V1 ) = 0.0_DP
        RightState(iPM_V2 ) = 0.0_DP
        RightState(iPM_V3 ) = 0.0_DP
        RightState(iPM_E  ) = 0.1_DP / ( Gamma_IDEAL - One )

        RightState(iPM_B1 ) = 0.0_DP
        RightState(iPM_B2 ) = 0.0_DP
        RightState(iPM_B3 ) = 0.0_DP
        RightState(iPM_Chi) = 0.0_DP

      CASE( 'ShockTube1' )

        ! --- ST1 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                   ---

        XD = Half

        LeftState(iPM_D ) =  1.0_DP
        LeftState(iPM_V1) =  0.0_DP
        LeftState(iPM_V2) =  0.0_DP
        LeftState(iPM_V3) =  0.0_DP
        LeftState(iPM_E ) =  1.0_DP / ( Gamma_IDEAL - One )

        CB1_L =  0.5_DP
        CB2_L =  One
        CB3_L =  0.0_DP

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) =  0.125_DP
        RightState(iPM_V1) =  0.0_DP
        RightState(iPM_V2) =  0.0_DP
        RightState(iPM_V3) =  0.0_DP
        RightState(iPM_E ) =  0.1_DP / ( Gamma_IDEAL - One )

        CB1_R =  0.5_DP
        CB2_R = -One
        CB3_R =  0.0_DP

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'HydroShockTube2' )

        ! --- Pure hydro version of ST2 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                                         ---

        XD = Half

        LeftState(iPM_D  ) = 1.08_DP
        LeftState(iPM_V1 ) = 0.4_DP
        LeftState(iPM_V2 ) = 0.3_DP
        LeftState(iPM_V3 ) = 0.2_DP
        LeftState(iPM_E  ) = 0.95_DP / ( Gamma_IDEAL - One )

        LeftState(iPM_B1 ) = 0.0_DP
        LeftState(iPM_B2 ) = 0.0_DP
        LeftState(iPM_B3 ) = 0.0_DP
        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D  ) = One
        RightState(iPM_V1 ) = -0.45_DP
        RightState(iPM_V2 ) = -0.2_DP
        RightState(iPM_V3 ) =  0.2_DP
        RightState(iPM_E  ) = One / ( Gamma_IDEAL - One )

        RightState(iPM_B1 ) = 0.0_DP
        RightState(iPM_B2 ) = 0.0_DP
        RightState(iPM_B3 ) = 0.0_DP
        RightState(iPM_Chi) = 0.0_DP

      CASE( 'ShockTube2' )

        ! --- ST2 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                   ---

        XD = Half

        LeftState(iPM_D ) = 1.08_DP
        LeftState(iPM_V1) = 0.4_DP
        LeftState(iPM_V2) = 0.3_DP
        LeftState(iPM_V3) = 0.2_DP
        LeftState(iPM_E ) = 0.95_DP / ( Gamma_IDEAL - One )

        CB1_L = 2.0_DP
        CB2_L = 0.3_DP
        CB3_L = 0.3_DP

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) = One
        RightState(iPM_V1) = -0.45_DP
        RightState(iPM_V2) = -0.2_DP
        RightState(iPM_V3) =  0.2_DP
        RightState(iPM_E ) = One / ( Gamma_IDEAL - One )

        CB1_R =  2.0_DP
        CB2_R = -0.7_DP
        CB3_R =  0.5_DP

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'HydroShockTube3' )

        ! --- Pure hydro version of ST3 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                                         ---

        XD = Half

        LeftState(iPM_D ) = One
        LeftState(iPM_V1) = 0.999_DP
        LeftState(iPM_V2) = 0.0_DP
        LeftState(iPM_V3) = 0.0_DP
        LeftState(iPM_E ) = 0.1_DP / ( Gamma_IDEAL - One )

        LeftState(iPM_B1) = 0.0_DP
        LeftState(iPM_B2) = 0.0_DP
        LeftState(iPM_B3) = 0.0_DP
        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) = One
        RightState(iPM_V1) = -0.999_DP
        RightState(iPM_V2) =  0.0_DP
        RightState(iPM_V3) =  0.0_DP
        RightState(iPM_E ) =  0.1_DP / ( Gamma_IDEAL - One )

        RightState(iPM_B1) = 0.0_DP
        RightState(iPM_B2) = 0.0_DP
        RightState(iPM_B3) = 0.0_DP
        RightState(iPM_Chi) = 0.0_DP

      CASE( 'ShockTube3' )

        ! --- ST3 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                   ---

        XD = Half

        LeftState(iPM_D ) = One
        LeftState(iPM_V1) = 0.999_DP
        LeftState(iPM_V2) = 0.0_DP
        LeftState(iPM_V3) = 0.0_DP
        LeftState(iPM_E ) = 0.1_DP / ( Gamma_IDEAL - One )

        CB1_L = 10.0_DP
        CB2_L = 7.0_DP
        CB3_L = 7.0_DP

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) = One
        RightState(iPM_V1) = -0.999_DP
        RightState(iPM_V2) =  0.0_DP
        RightState(iPM_V3) =  0.0_DP
        RightState(iPM_E ) =  0.1_DP / ( Gamma_IDEAL - One )

        CB1_R =  10.0_DP
        CB2_R = -7.0_DP
        CB3_R = -7.0_DP ! Typo (7.0_DP instead of -7.0_DP) in Mattia & Mignone (2022).

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'ShockTube4' )

        ! --- ST1 from Mattia & Mignone, 2022, MRNAS, ---
        ! --- 510, 481-499, Table 1                   ---

        XD = Half

        LeftState(iPM_D ) =  One
        LeftState(iPM_V1) =  0.0_DP
        LeftState(iPM_V2) =  0.3_DP
        LeftState(iPM_V3) =  0.4_DP
        LeftState(iPM_E ) =  5.0_DP / ( Gamma_IDEAL - One )

        CB1_L =  One
        CB2_L =  6.0_DP
        CB3_L =  2.0_DP

        VSq_L = LeftState(iPM_V1)**2 + LeftState(iPM_V2)**2 + LeftState(iPM_V3)**2

        W_L = One / SQRT( One - VSq_L )

        VdotB_L = LeftState(iPM_V1) * CB1_L &
                    + LeftState(iPM_V2) * CB2_L &
                    + LeftState(iPM_V3) * CB3_L

        LeftState(iPM_B1) = W_L * VdotB_L * LeftState(iPM_V1) + CB1_L / W_L
        LeftState(iPM_B2) = W_L * VdotB_L * LeftState(iPM_V2) + CB2_L / W_L
        LeftState(iPM_B3) = W_L * VdotB_L * LeftState(iPM_V3) + CB3_L / W_L

        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D ) =  0.9_DP
        RightState(iPM_V1) =  0.0_DP
        RightState(iPM_V2) =  0.0_DP
        RightState(iPM_V3) =  0.0_DP
        RightState(iPM_E ) =  5.3_DP / ( Gamma_IDEAL - One )

        CB1_R = One
        CB2_R = 5.0_DP
        CB3_R = 2.0_DP

        VSq_R = RightState(iPM_V1)**2 + RightState(iPM_V2)**2 + RightState(iPM_V3)**2

        W_R = One / SQRT( One - VSq_R )

        VdotB_R = RightState(iPM_V1) * CB1_R &
                    + RightState(iPM_V2) * CB2_R &
                    + RightState(iPM_V3) * CB3_R

        RightState(iPM_B1) = W_R * VdotB_R * RightState(iPM_V1) + CB1_R / W_R
        RightState(iPM_B2) = W_R * VdotB_R * RightState(iPM_V2) + CB2_R / W_R
        RightState(iPM_B3) = W_R * VdotB_R * RightState(iPM_V3) + CB3_R / W_R

        RightState(iPM_Chi) = 0.0_DP

      CASE( 'HydroMBProblem1' )

        ! --- Problem 1 from Mignone & Bodo, 2005, MRNAS, ---
        ! --- 364, 126                                    ---

        XD = 0.5_DP

        LeftState(iPM_D  ) = 1.0_DP
        LeftState(iPM_V1 ) = 0.9_DP
        LeftState(iPM_V2 ) = 0.0_DP
        LeftState(iPM_V3 ) = 0.0_DP
        LeftState(iPM_E  ) = 1.0_DP / ( Gamma_IDEAL - One )
        LeftState(iPM_B1 ) = 0.0_DP
        LeftState(iPM_B2 ) = 0.0_DP
        LeftState(iPM_B3 ) = 0.0_DP
        LeftState(iPM_Chi) = 0.0_DP

        RightState(iPM_D  ) = 1.0_DP
        RightState(iPM_V1 ) = 0.0_DP
        RightState(iPM_V2 ) = 0.0_DP
        RightState(iPM_V3 ) = 0.0_DP
        RightState(iPM_E  ) = 10.0_DP / ( Gamma_IDEAL - One )
        RightState(iPM_B1 ) = 0.0_DP
        RightState(iPM_B2 ) = 0.0_DP
        RightState(iPM_B3 ) = 0.0_DP
        RightState(iPM_Chi) = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', RiemannProblemName
        WRITE(*,'(A)') 'Valid choices:'
         WRITE(*,'(A)') &
          "  'HydroIsolatedContact' - &
          Pure hydro version of isolated contact problem from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'IsolatedContact' - &
          Isolated contact problem from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'RotationalWave' - &
          Rotational wave problem from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'HydroShockTube1' - &
          Pure hydro version of 1st shock tube problem (ST1) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'ShockTube1' - &
          1st shock tube problem (ST1) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'HydroShockTube2' - &
          Pure hydro version of 2nd shock tube problem (ST2) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'ShockTube2' - &
          2nd shock tube problem (ST2) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'HydroShockTube3' - &
          Pure hydro version of the 3rd shock tube problem (ST3) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'ShockTube3' - &
          3rd shock tube problem (ST3) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'ShockTube4' - &
          4th shock tube problem (ST4) from Mattia & Mignone, 2022, MRNAS, 510, 481-499, Table 1"
        WRITE(*,'(A)') &
          "  'HydroMBProblem1' - &
          1st hydro shock tube problem from Mignone & Bodo, 2005, MRNAS, 364, 126"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Right State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_D   = ', RightState(iPM_D  )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V1  = ', RightState(iPM_V1 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V2  = ', RightState(iPM_V2 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V3  = ', RightState(iPM_V3 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_E   = ', RightState(iPM_E  )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B1  = ', RightState(iPM_B1 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B2  = ', RightState(iPM_B2 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B3  = ', RightState(iPM_B3 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_Chi = ', RightState(iPM_Chi)
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Left State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_D   = ', LeftState(iPM_D  )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V1  = ', LeftState(iPM_V1 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V2  = ', LeftState(iPM_V2 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_V3  = ', LeftState(iPM_V3 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_E   = ', LeftState(iPM_E  )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B1  = ', LeftState(iPM_B1 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B2  = ', LeftState(iPM_B2 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_B3  = ', LeftState(iPM_B3 )
    WRITE(*,'(8x,A,ES24.16E3)') 'PM_Chi = ', LeftState(iPM_Chi)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. XD )THEN

          uPM(iNodeX,iX1,iX2,iX3,iPM_D  ) = LeftState(iPM_D  )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V1 ) = LeftState(iPM_V1 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V2 ) = LeftState(iPM_V2 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V3 ) = LeftState(iPM_V3 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_E  ) = LeftState(iPM_E  )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = LeftState(iPM_B1 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = LeftState(iPM_B2 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = LeftState(iPM_B3 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = LeftState(iPM_Chi)

        ELSE

          uPM(iNodeX,iX1,iX2,iX3,iPM_D  ) = RightState(iPM_D  )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V1 ) = RightState(iPM_V1 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V2 ) = RightState(iPM_V2 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_V3 ) = RightState(iPM_V3 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_E  ) = RightState(iPM_E  )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = RightState(iPM_B1 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = RightState(iPM_B2 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = RightState(iPM_B3 )
          uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = RightState(iPM_Chi)

        END IF

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_E ), &
               uPM(:,iX1,iX2,iX3,iPM_Ne), uAM(:,iX1,iX2,iX3,iAM_P) )

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Riemann1D


  SUBROUTINE InitializeFields_MMBlastWave2D &
               ( MMBlastWaveB0, MMBlastWavePhi, EvolveOnlyMagnetic )

    LOGICAL,  INTENT(in) :: EvolveOnlyMagnetic
    REAL(DP), INTENT(in) :: MMBlastWaveB0, MMBlastWavePhi

    INTEGER iX1, iX2, iX3
    INTEGER iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: B0, Phi

    ! --- 2D cylindrical blast wave in Cartesian coordinates ---
    ! --- from Section 4.4 of Mattia & Mignone, 2022, MRNAS, ---
    ! --- 510, 481-499                                       ---

    B0 = MMBlastWaveB0
    Phi = MMBlastWavePhi

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        IF( SQRT( X1**2 + X2**2 ) .LE. 0.8_DP ) THEN

          uPM(iNodeX,iX1,iX2,iX3,iPM_D   ) = 0.01_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_V1  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_V2  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_V3  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_E   ) = One / ( Gamma_IDEAL - One )

        ELSE

          uPM(iNodeX,iX1,iX2,iX3,iPM_D   ) = 1.0d-4
          uPM(iNodeX,iX1,iX2,iX3,iPM_V1  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_V2  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_V3  ) = 0.0_DP
          uPM(iNodeX,iX1,iX2,iX3,iPM_E   ) = 5.0d-3 / ( Gamma_IDEAL - One )

        END IF

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = B0 * COS( Phi )
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = B0 * SIN( Phi )
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_E ), &
               uPM(:,iX1,iX2,iX3,iPM_Ne), uAM(:,iX1,iX2,iX3,iAM_P) )

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_MMBlastWave2D


  SUBROUTINE InitializeFields_OrszagTang2D( OTScaleFactor, EvolveOnlyMagnetic )

    LOGICAL,  INTENT(in) :: EvolveOnlyMagnetic
    REAL(DP), INTENT(in) :: OTScaleFactor

    INTEGER iX1, iX2, iX3
    INTEGER iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: V1, V2, V3, VSq, W
    REAL(DP) :: CB1, CB2, CB3, VdotB

    ! --- 2D Orszag-Tang vortex based on the approach taken  ---
    ! --- in Section 4.6 of Anninos et al., 2017, ApJ Supp., ---
    ! --- 231:17                                             ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        V1 = -SIN( TwoPi * X2 ) / OTScaleFactor
        V2 =  SIN( TwoPi * X1 ) / OTScaleFactor
        V3 = Zero

        CB1 = -SIN( TwoPi  * X2 ) / ( SQRT( FourPi ) * OTScaleFactor )
        CB2 =  SIN( FourPi * X1 ) / ( SQRT( FourPi ) * OTScaleFactor )
        CB3 =  Zero

        uPM(iNodeX,iX1,iX2,iX3,iPM_D   ) = Gamma_IDEAL**2 / ( FourPi )
        uPM(iNodeX,iX1,iX2,iX3,iPM_V1  ) = V1
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2  ) = V2
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3  ) = V3
        uPM(iNodeX,iX1,iX2,iX3,iPM_E   ) = Gamma_IDEAL &
                                           / ( FourPi * ( Gamma_IDEAL - One ) &
                                               * OTScaleFactor**2 )
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = Zero

        VSq = V1**2 + V2**2 + V3**2

        VdotB = V1 * CB1 &
                + V2 * CB2 &
                + V3 * CB3

        W = One / SQRT( One - VSq )

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = W * VdotB * V1 + CB1 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = W * VdotB * V2 + CB2 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = W * VdotB * V3 + CB3 / W

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_E ), &
               uPM(:,iX1,iX2,iX3,iPM_Ne), uAM(:,iX1,iX2,iX3,iAM_P) )

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_OrszagTang2D


  SUBROUTINE InitializeFields_ShearingDisk( EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    CHARACTER(256) :: FileName

    INTEGER(HID_T) :: FILE_ID
    INTEGER        :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3
    REAL(DP) :: V1, V2, V3, VSq, W, CB1, CB2, CB3, VdotB
    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), V3Arr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    FileName = "/home/jbuffal/thornado_MHD_3D/Workflow/MHD/ShearingDisk/GR_LR_diffrot.h5"

    ! --- Populate arrays ---

    CALL H5OPEN_F( HDFERR )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    nX = 10000

    ALLOCATE( PressureArr(nX), DensityArr(nX), V3Arr(nX), AlphaArr(nX), &
              PsiArr(nX), X1Arr(nX) )

    CALL ReadDataset1DHDF( PsiArr,      '/psi',   FILE_ID )
    CALL ReadDataset1DHDF( AlphaArr,    '/alpha', FILE_ID )
    CALL ReadDataset1DHDF( X1Arr,       '/r',     FILE_ID )
    CALL ReadDataset1DHDF( PressureArr, '/pres',  FILE_ID )
    CALL ReadDataset1DHDF( DensityArr,  '/rho',   FILE_ID )
    CALL ReadDataset1DHDF( V3Arr,       '/V3',    FILE_ID )

    X1Arr       = X1Arr       * Centimeter
    DensityArr  = DensityArr  * ( Gram / Centimeter**3 )
    PressureArr = PressureArr * ( Erg  / Centimeter**3 )
    V3Arr       = V3Arr       * ( One  / Second )

    ! --- Map to 3D domain ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        ! --- Geometry Fields ---

        uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha) &
          = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Psi) &
          = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

        uGF(iNodeX,iX1,iX2,iX3,iGF_h_1) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_2) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_3) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1

        CALL ComputeGeometryX_FromScaleFactors( uGF(:,iX1,iX2,iX3,:) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        uPM(iNodeX,iX1,iX2,iX3,iPM_D) &
          = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

        V1 = Zero
        V2 = Zero
        V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

        VSq = uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) * V1**2 &
              + uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) * V2**2 &
              + uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) * V3**2

        W = One / SQRT( One - VSq )

        uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = V1
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = V2
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = V3

        uPM(iNodeX,iX1,iX2,iX3,iPM_E) &
          = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
            / ( Gamma_IDEAL - One )

        CB1 = Zero
        CB2 = 2.0 * 1.0d13 * Gauss
        CB3 = Zero

        VdotB = uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                + uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                + uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W

        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = Zero

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_E ), &
               uPM(:,iX1,iX2,iX3,iPM_Ne), uAM(:,iX1,iX2,iX3,iAM_P) )

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

    DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, V3Arr, PressureArr )

  END SUBROUTINE InitializeFields_ShearingDisk


  SUBROUTINE InitializeFields_MagneticKH( EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    INTEGER iX1, iX2, iX3
    INTEGER iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: V1, V2, V3, VSq, W
    REAL(DP) :: CB1, CB2, CB3, VdotB

    REAL(DP) :: A0, VSh, a, sigma
    REAL(DP) :: rhoL, rhoH

    A0 = 1.0d-1
    VSh = Half
    a = 1.0d-2
    sigma = 0.1

    rhoL = 1.0d-2
    rhoH = One

    ! --- 2D magnetic Kelvin-Helmholtz setup from Mattia and Mignone (2022) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        V1 = DSIGN( One, X2 ) * VSh &
             * TANH( ( Two * X2 - DSIGN( One, X2 ) ) / ( Two * a ) )
        V2 = DSIGN( One, X2 ) * A0 * VSh &
             * SIN( Two * Pi * X1 ) &
             * EXP( - ( ( Two * X2 - DSIGN( One, X2 ) ) / ( Two * sigma ) )**2 )
        V3 = Zero

        CB1 = 1.0d-3
        CB2 = Zero
        CB3 = Zero

        uPM(iNodeX,iX1,iX2,iX3,iPM_D   ) = Half * ( rhoL + rhoH ) &
                                           + Half * ( rhoH - rhoL ) * ( V1 / VSh )
        uPM(iNodeX,iX1,iX2,iX3,iPM_V1  ) = V1
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2  ) = V2
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3  ) = V3
        uPM(iNodeX,iX1,iX2,iX3,iPM_E   ) = One / ( Gamma_IDEAL - One )
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = Zero

        VSq = V1**2 + V2**2 + V3**2

        VdotB = V1 * CB1 &
                + V2 * CB2 &
                + V3 * CB3

        W = One / SQRT( One - VSq )

        uPM(iNodeX,iX1,iX2,iX3,iPM_B1 ) = W * VdotB * V1 + CB1 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B2 ) = W * VdotB * V2 + CB2 / W
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3 ) = W * VdotB * V3 + CB3 / W

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_E ), &
               uPM(:,iX1,iX2,iX3,iPM_Ne), uAM(:,iX1,iX2,iX3,iAM_P) )

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P), &
               EvolveOnlyMagnetic )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_MagneticKH


  SUBROUTINE ReadDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(1)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset1DHDF


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    !PRINT*, 'i: ', i

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

      !PRINT*, 'x(1): ', x(1)
      !PRINT*, 'x(2): ', x(2)
      !PRINT*, 'y(1): ', y(1)
      !PRINT*, 'y(2): ', y(2)

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

      !PRINT*, 'x(n-1): ', x(n-1)
      !PRINT*, 'x(n): ',   x(n)
      !PRINT*, 'y(n-1): ', y(n-1)
      !PRINT*, 'y(n): ',   y(n)


    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

      !PRINT*, 'x(i): ', x(i)
      !PRINT*, 'x(i+1): ', x(i+1)
      !PRINT*, 'y(i): ', y(i)
      !PRINT*, 'y(i+1): ', y(i+1)

    END IF

    RETURN

  END FUNCTION Interpolate1D


END MODULE InitializationModule
