MODULE MF_InitializationModule_Relativistic_IDEAL

  ! --- AMReX Modules ---

  USE amrex_box_module,        ONLY: &
    amrex_box
  USE amrex_geometry_module,   ONLY: &
    amrex_geometry
  USE amrex_multifab_module,   ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,   ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module,  ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,     ONLY: &
    nDOFX,   &
    nX,      &
    nNodesX, &
    swX,     &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX1
  USE MeshModule,              ONLY: &
    MeshType,    &
    CreateMesh,  &
    DestroyMesh, &
    NodeCoordinate
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule,    ONLY: &
    nGF,          &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi
  USE MagnetofluidFieldsModule,       ONLY: &
    nCM,     &
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
    nPM,     &
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
    nAM,     &
    iAM_P
  USE MHD_UtilitiesModule,   ONLY: &
    ComputeConserved_MHD
  USE EquationOfStateModule,   ONLY: &
    ComputePressureFromPrimitive
  USE UnitsModule,             ONLY: &
    Kilometer,    &
    Second,       &
    SolarMass,    &
    Gram,         &
    Centimeter,   &
    Erg,          &
    SpeedOfLight, &
    GravitationalConstant, &
    Millisecond, &
    Gauss
  USE UtilitiesModule,         ONLY: &
    NodeNumberX, &
    Locate, &
    Interpolate1D_Linear

  ! --- Local Modules ---

  USE MF_KindModule,           ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    Pi, &
    TwoPi, &
    FourPi
  USE InputParsingModule,      ONLY: &
    nLevels,            &
    xL,                 &
    xR,                 &
    Gamma_IDEAL,        &
    UseTiling,          &
    t_end,              &
    EvolveOnlyMagnetic
  USE MF_UtilitiesModule,      ONLY: &
    amrex2thornado_X_Global

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: HDFERR

  PUBLIC :: MF_InitializeFields_Relativistic_IDEAL


CONTAINS


  SUBROUTINE MF_InitializeFields_Relativistic_IDEAL &
    ( ProgramName, MF_uGF, MF_uCM, GEOM )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )

    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection( MF_uGF, MF_uCM )

      CASE( 'Cleaning1D' )

        CALL InitializeFields_Cleaning1D( MF_uGF, MF_uCM )

      CASE( 'Cleaning2D' )

        CALL InitializeFields_Cleaning2D( MF_uGF, MF_uCM )

      CASE( 'OrszagTang2D' )

        CALL InitializeFields_OrszagTang2D( MF_uGF, MF_uCM )

      CASE( 'ShearingDisk' )

        CALL InitializeFields_ShearingDisk( MF_uGF, MF_uCM )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Advection'
          WRITE(*,'(6x,A)')     'Cleaning1D'
          WRITE(*,'(6x,A)')     'Cleaning2D'
          WRITE(*,'(6x,A)')     'OrszagTang2D'
          WRITE(*,'(6x,A)')     'MagnetizedKH_3D'
          WRITE(*,'(6x,A)')     'ShearingDisk'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields_Relativistic_IDEAL


  SUBROUTINE InitializeFields_Advection( MF_uGF, MF_uCM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2, iNX3
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCM_K(nDOFX,nCM)
    REAL(DP)       :: uPM_K(nDOFX,nPM)
    REAL(DP)       :: uAM_K(nDOFX,nAM)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: AdvectionProfile
    REAL(DP) :: W, CB1, CB2, CB3, V1, V2, V3, VdotB
    REAL(DP) :: Eta, k, h, VA

    AdvectionProfile = 'HydroSineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( AdvectionProfile ) .NE. 'HydroSineWaveX1' &
        .AND. TRIM( AdvectionProfile ) .NE. 'HydroSineWaveX2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'HydroSineWaveX3' &
        .AND. TRIM( AdvectionProfile ) .NE. 'HydroSineWaveX1X2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'MagneticSineWaveX1' &
        .AND. TRIM( AdvectionProfile ) .NE. 'MagneticSineWaveX2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'MagneticSineWaveX1X2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'MagneticSineWaveX3' &
        .AND. TRIM( AdvectionProfile ) .NE. 'CPAlfvenX1' &
        .AND. TRIM( AdvectionProfile ) .NE. 'CPAlfvenX2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'CPAlfvenX3' &
        .AND. TRIM( AdvectionProfile ) .NE. 'CPAlfvenOblique' &
        .AND. TRIM( AdvectionProfile ) .NE. 'LoopAdvection2D' &
        .AND. TRIM( AdvectionProfile ) .NE. 'LoopAdvection3D' )THEN

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for AdvectionProfile: ', &
          TRIM( AdvectionProfile )
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') '  HydroSineWaveX1'
        WRITE(*,'(A)') '  HydroSineWaveX2'
        WRITE(*,'(A)') '  HydroSineWaveX3'
        WRITE(*,'(A)') '  HydroSineWaveX1X2'
        WRITE(*,'(A)') '  MagneticSineWaveX1'
        WRITE(*,'(A)') '  MagneticSineWaveX2'
        WRITE(*,'(A)') '  MagneticSineWaveX1X2'
        WRITE(*,'(A)') '  MagneticSineWaveX3'
        WRITE(*,'(A)') '  CPAlfvenX1'
        WRITE(*,'(A)') '  CPAlfvenX2'
        WRITE(*,'(A)') '  CPAlfvenX3'
        WRITE(*,'(A)') '  CPAlfvenOblique'
        WRITE(*,'(A)') '  LoopAdvection2D'
        WRITE(*,'(A)') '  LoopAdvection3D'
        WRITE(*,*)
        WRITE(*,'(A)') 'Stopping...'

      END IF

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(4x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )

    END IF

    uGF_K = Zero
    uCM_K = Zero
    uPM_K = Zero
    uAM_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCM )
        hi_F = UBOUND( uCM )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)
            iNX3 = NodeNumberTableX(3,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

            IF     ( TRIM( AdvectionProfile ) .EQ. 'HydroSineWaveX1' )THEN

              uPM_K(iNX,iPM_D )  = One + 0.1_DP * SIN( TwoPi * X1 )
              uPM_K(iNX,iPM_V1)  = 0.1_DP
              uPM_K(iNX,iPM_V2)  = Zero
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = One / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = Zero
              uPM_K(iNX,iPM_B2)  = Zero
              uPM_K(iNX,iPM_B3)  = Zero
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'HydroSineWaveX2' )THEN

              uPM_K(iNX,iPM_D )  = One + 0.1_DP * SIN( TwoPi * X2 )
              uPM_K(iNX,iPM_V1)  = Zero
              uPM_K(iNX,iPM_V2)  = 0.1_DP
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = One / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = Zero
              uPM_K(iNX,iPM_B2)  = Zero
              uPM_K(iNX,iPM_B3)  = Zero
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'HydroSineWaveX3' )THEN

              uPM_K(iNX,iPM_D )  = One + 0.1_DP * SIN( TwoPi * X3 )
              uPM_K(iNX,iPM_V1)  = Zero
              uPM_K(iNX,iPM_V2)  = Zero
              uPM_K(iNX,iPM_V3)  = 0.1_DP
              uPM_K(iNX,iPM_E )  = One / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = Zero
              uPM_K(iNX,iPM_B2)  = Zero
              uPM_K(iNX,iPM_B3)  = Zero
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'HydroSineWaveX1X2' )THEN

              uPM_K(iNX,iPM_D ) &
                = One + 0.1_DP * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
              uPM_K(iNX,iPM_V1)  = 0.1_DP / SQRT( Two )
              uPM_K(iNX,iPM_V2)  = 0.1_DP / SQRT( Two )
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = One / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = Zero
              uPM_K(iNX,iPM_B2)  = Zero
              uPM_K(iNX,iPM_B3)  = Zero
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'MagneticSineWaveX1' )THEN

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1)  = 0.1_DP
              uPM_K(iNX,iPM_V2)  = Zero
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = 1.0d-4 / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = 1.0d-4
              uPM_K(iNX,iPM_B2)  = 1.0d-4 * SIN( TwoPi * X1 )
              uPM_K(iNX,iPM_B3)  = 1.0d-4 * COS( TwoPi * X1 )
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'MagneticSineWaveX2' )THEN

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1)  = 0.1_DP
              uPM_K(iNX,iPM_V2)  = Zero
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = 1.0d-4 / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = 1.0d-4 * SIN( TwoPi * X2 )
              uPM_K(iNX,iPM_B2)  = 1.0d-4
              uPM_K(iNX,iPM_B3)  = 1.0d-4 * COS( TwoPi * X2 )
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'MagneticSineWaveX1X2' )THEN

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1)  = 0.1_DP / SQRT( Two )
              uPM_K(iNX,iPM_V2)  = 0.1_DP / SQRT( Two )
              uPM_K(iNX,iPM_V3)  = Zero
              uPM_K(iNX,iPM_E )  = 1.0d-4 / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = 1.0d-4 / SQRT( Two )
              uPM_K(iNX,iPM_B2)  = 1.0d-4 / SQRT( Two )
              uPM_K(iNX,iPM_B3)  = 1.0d-4 * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'MagneticSineWaveX3' )THEN

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1)  = Zero
              uPM_K(iNX,iPM_V2)  = Zero
              uPM_K(iNX,iPM_V3)  = 0.1_DP
              uPM_K(iNX,iPM_E )  = 1.0d-4 / ( Gamma_IDEAL - One )
              uPM_K(iNX,iPM_B1)  = 1.0d-4 * SIN( TwoPi * X3 )
              uPM_K(iNX,iPM_B2)  = 1.0d-4 * COS( TwoPi * X3 )
              uPM_K(iNX,iPM_B3)  = 1.0d-4
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'CPAlfvenX1' )THEN

              ! Circularly polarized Alfven wave with the exact solution
              ! from Del Zanna et al. (2007) and Mattia and Mignone (2022).

              Eta = One
              k   = One
              h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
              VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                          * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

              W = One / SQRT( One - VA**2 * Eta**2 )

              CB1 = One
              CB2 = Eta * COS( k * X1 )
              CB3 = Eta * SIN( k * X1 )

              V1 = Zero
              V2 = -VA * CB2
              V3 = -VA * CB3

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1) = V1
              uPM_K(iNX,iPM_V2) = V2
              uPM_K(iNX,iPM_V3) = V3
              uAM_K(iNX,iAM_P )  = One
              uPM_K(iNX,iPM_E )  &
                = uAM_K(iNX,iAM_P) / ( Gamma_IDEAL - One )

              VdotB = -VA * Eta**2

              uPM_K(iNX,iPM_B1) &
                = W * VdotB * V1 + CB1 / W
              uPM_K(iNX,iPM_B2) &
                = W * VdotB * V2 + CB2 / W
              uPM_K(iNX,iPM_B3) &
                = W * VdotB * V3 + CB3 / W
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'CPAlfvenX2' )THEN

              ! Circularly polarized Alfven wave with the exact solution
              ! from Del Zanna et al. (2007) and Mattia and Mignone (2022).

              Eta = One
              k   = One
              h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
              VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                          * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

              W = One / SQRT( One - VA**2 * Eta**2 )

              CB1 = Eta * COS( k * X2 )
              CB2 = One
              CB3 = Eta * SIN( k * X2 )

              V1 = -VA * CB1
              V2 = Zero
              V3 = -VA * CB3

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1) = V1
              uPM_K(iNX,iPM_V2) = V2
              uPM_K(iNX,iPM_V3) = V3
              uAM_K(iNX,iAM_P )  = One
              uPM_K(iNX,iPM_E )  &
                = uAM_K(iNX,iAM_P) / ( Gamma_IDEAL - One )

              VdotB = -VA * Eta**2

              uPM_K(iNX,iPM_B1) &
                = W * VdotB * V1 + CB1 / W
              uPM_K(iNX,iPM_B2) &
                = W * VdotB * V2 + CB2 / W
              uPM_K(iNX,iPM_B3) &
                = W * VdotB * V3 + CB3 / W
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'CPAlfvenX3' )THEN

              ! Circularly polarized Alfven wave with the exact solution
              ! from Del Zanna et al. (2007) and Mattia and Mignone (2022).

              Eta = One
              k   = One
              h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
              VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                          * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

              W = One / SQRT( One - VA**2 * Eta**2 )

              CB1 = Eta * COS( k * X3 )
              CB2 = Eta * SIN( k * X3 )
              CB3 = One

              V1 = -VA * CB1
              V2 = -VA * CB2
              V3 = Zero

              uPM_K(iNX,iPM_D )  = One
              uPM_K(iNX,iPM_V1) = V1
              uPM_K(iNX,iPM_V2) = V2
              uPM_K(iNX,iPM_V3) = V3
              uAM_K(iNX,iAM_P )  = One
              uPM_K(iNX,iPM_E )  &
                = uAM_K(iNX,iAM_P) / ( Gamma_IDEAL - One )

              VdotB = -VA * Eta**2

              uPM_K(iNX,iPM_B1) &
                = W * VdotB * V1 + CB1 / W
              uPM_K(iNX,iPM_B2) &
                = W * VdotB * V2 + CB2 / W
              uPM_K(iNX,iPM_B3) &
                = W * VdotB * V3 + CB3 / W
              uPM_K(iNX,iPM_Chi) = Zero

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'LoopAdvection2D' )THEN

              ! V3 = 0 variant of the 2D field loop advection
              ! problem from Section 5.5 of Mosta et al. (2014).

              uPM_K(iNX,iPM_D ) = One
              uPM_K(iNX,iPM_V1) = Half
              uPM_K(iNX,iPM_V2) = One / 24.0_DP
              uPM_K(iNX,iPM_V3) = Zero
              uPM_K(iNX,iPM_E ) = Three / ( Gamma_IDEAL - One )

              W = One / SQRT( One - ( uPM_K(iNX,iPM_V1)**2 &
                                      + uPM_K(iNX,iPM_V2)**2 &
                                      + uPM_K(iNX,iPM_V3)**2 ) )

              IF( SQRT( X1**2 + X2**2 ) < 0.3_DP )THEN

                CB1 = -1.0d-3 * X2 / SQRT( X1**2 + X2**2 )
                CB2 =  1.0d-3 * X1 / SQRT( X1**2 + X2**2 )
                CB3 =  Zero

              ELSE

                CB1 = Zero
                CB2 = Zero
                CB3 = Zero

              END IF

              VdotB = uPM_K(iNX,iPM_V1) * CB1 &
                      + uPM_K(iNX,iPM_V2) * CB2 &
                      + uPM_K(iNX,iPM_V3) * CB3

              uPM_K(iNX,iPM_B1) = W * VdotB * uPM_K(iNX,iPM_V1) + CB1 / W
              uPM_K(iNX,iPM_B2) = W * VdotB * uPM_K(iNX,iPM_V2) + CB2 / W
              uPM_K(iNX,iPM_B3) = W * VdotB * uPM_K(iNX,iPM_V3) + CB3 / W

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'LoopAdvection3D' )THEN

              ! 3D field loop advection problem from Section 5.5 of
              ! Mosta et al. (2014).

              uPM_K(iNX,iPM_D ) = One
              uPM_K(iNX,iPM_V1) = 0.2_DP * SQRT( Two )
              uPM_K(iNX,iPM_V2) = 0.2_DP
              uPM_K(iNX,iPM_V3) = 0.1_DP
              uPM_K(iNX,iPM_E ) = Three / ( Gamma_IDEAL - One )

              W = One / SQRT( One - ( uPM_K(iNX,iPM_V1)**2 &
                                      + uPM_K(iNX,iPM_V2)**2 &
                                      + uPM_K(iNX,iPM_V3)**2 ) )

              IF( SQRT( X1**2 + X2**2 ) .LE. 0.3_DP )THEN

                CB1 = -( SQRT( Two ) / Two ) * 1.0d-3 * X2 / SQRT( X1**2 + X2**2 )
                CB2 =  1.0d-3 * X1 / SQRT( X1**2 + X2**2 )
                CB3 =  ( SQRT( Two ) / Two ) * 1.0d-3 * X2 / SQRT( X1**2 + X2**2 )

              ELSE

                CB1 = Zero
                CB2 = Zero
                CB3 = Zero

              END IF

              VdotB = uPM_K(iNX,iPM_V1) * CB1 &
                      + uPM_K(iNX,iPM_V2) * CB2 &
                      + uPM_K(iNX,iPM_V3) * CB3

              uPM_K(iNX,iPM_B1) = W * VdotB * uPM_K(iNX,iPM_V1) + CB1 / W
              uPM_K(iNX,iPM_B2) = W * VdotB * uPM_K(iNX,iPM_V2) + CB2 / W
              uPM_K(iNX,iPM_B3) = W * VdotB * uPM_K(iNX,iPM_V3) + CB3 / W

            END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPM_K(:,iPM_D), uPM_K(:,iPM_E), uPM_K(:,iPM_Ne), &
                   uAM_K(:,iAM_P) )

          CALL ComputeConserved_MHD &
                 ( uPM_K(:,iPM_D ), uPM_K(:,iPM_V1), uPM_K(:,iPM_V2), &
                   uPM_K(:,iPM_V3), uPM_K(:,iPM_E ), uPM_K(:,iPM_Ne), &
                   uPM_K(:,iPM_B1), uPM_K(:,iPM_B2), uPM_K(:,iPM_B3), &
                   uPM_K(:,iPM_Chi), &
                   uCM_K(:,iCM_D  ), uCM_K(:,iCM_S1), uCM_K(:,iCM_S2), &
                   uCM_K(:,iCM_S3 ), uCM_K(:,iCM_E ), uCM_K(:,iCM_Ne), &
                   uCM_K(:,iCM_B1 ), uCM_K(:,iCM_B2), uCM_K(:,iCM_B3), &
                   uCM_K(:,iCM_Chi), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uGF_K(:,iGF_Alpha   ), &
                   uGF_K(:,iGF_Beta_1  ), &
                   uGF_K(:,iGF_Beta_2  ), &
                   uGF_K(:,iGF_Beta_3  ), &
                   uAM_K(:,iAM_P), &
                   EvolveOnlyMagnetic )

          uCM(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCM_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_OrszagTang2D( MF_uGF, MF_uCM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2, iNX3
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCM_K(nDOFX,nCM)
    REAL(DP)       :: uPM_K(nDOFX,nPM)
    REAL(DP)       :: uAM_K(nDOFX,nAM)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    REAL(DP) :: OTScaleFactor
    REAL(DP) :: CB1, CB2, CB3, V1, V2, V3, VSq, W, VdotB

    OTScaleFactor = 1.0d2
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'OTScaleFactor', OTScaleFactor )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCM_K = Zero
    uPM_K = Zero
    uAM_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCM )
        hi_F = UBOUND( uCM )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)
            iNX3 = NodeNumberTableX(3,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

            V1 = -SIN( TwoPi * X2 ) / OTScaleFactor
            V2 =  SIN( TwoPi * X1 ) / OTScaleFactor
            V3 = Zero

            CB1 = -SIN( TwoPi  * X2 ) / ( SQRT( FourPi ) * OTScaleFactor )
            CB2 =  SIN( FourPi * X1 ) / ( SQRT( FourPi ) * OTScaleFactor )
            CB3 = Zero

            uPM_K(iNX,iPM_D )  = Gamma_IDEAL**2 / FourPi
            uPM_K(iNX,iPM_V1)  = V1
            uPM_K(iNX,iPM_V2)  = V2
            uPM_K(iNX,iPM_V3)  = V3
            uPM_K(iNX,iPM_E )  = Gamma_IDEAL &
                                 / ( FourPi * ( Gamma_IDEAL - One ) &
                                     * OTScaleFactor**2 )
            uPM_K(iNX,iPM_Chi) = Zero

            VSq = V1**2 + V2**2 + V3**2

            VdotB = V1 * CB1 &
                    + V2 * CB2 &
                    + V3 * CB3

            W = One / SQRT( One - VSq )

            uPM_K(iNX,iPM_B1) = W * VdotB * V1 + CB1 / W
            uPM_K(iNX,iPM_B2) = W * VdotB * V2 + CB2 / W
            uPM_K(iNX,iPM_B3) = W * VdotB * V3 + CB3 / W

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPM_K(:,iPM_D), uPM_K(:,iPM_E), uPM_K(:,iPM_Ne), &
                   uAM_K(:,iAM_P) )

          CALL ComputeConserved_MHD &
                 ( uPM_K(:,iPM_D ), uPM_K(:,iPM_V1), uPM_K(:,iPM_V2), &
                   uPM_K(:,iPM_V3), uPM_K(:,iPM_E ), uPM_K(:,iPM_Ne), &
                   uPM_K(:,iPM_B1), uPM_K(:,iPM_B2), uPM_K(:,iPM_B3), &
                   uPM_K(:,iPM_Chi), &
                   uCM_K(:,iCM_D  ), uCM_K(:,iCM_S1), uCM_K(:,iCM_S2), &
                   uCM_K(:,iCM_S3 ), uCM_K(:,iCM_E ), uCM_K(:,iCM_Ne), &
                   uCM_K(:,iCM_B1 ), uCM_K(:,iCM_B2), uCM_K(:,iCM_B3), &
                   uCM_K(:,iCM_Chi), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uGF_K(:,iGF_Alpha   ), &
                   uGF_K(:,iGF_Beta_1  ), &
                   uGF_K(:,iGF_Beta_2  ), &
                   uGF_K(:,iGF_Beta_3  ), &
                   uAM_K(:,iAM_P), &
                   EvolveOnlyMagnetic )

          uCM(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCM_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_OrszagTang2D


  SUBROUTINE InitializeFields_Cleaning1D( MF_uGF, MF_uCM )

    ! 1D divergence cleaning test from Section 5.1 of
    ! Derigs et al. (2018) with option to use
    ! constant density.

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCM_K(nDOFX,nCM)
    REAL(DP)       :: uPM_K(nDOFX,nPM)
    REAL(DP)       :: uAM_K(nDOFX,nAM)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    LOGICAL  :: SmoothProfile
    REAL(DP) :: CB1, CB2, CB3, D, V1, V2, V3, VSq, W, VdotB

    SmoothProfile = .TRUE.
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'SmoothProfile', SmoothProfile )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCM_K = Zero
    uPM_K = Zero
    uAM_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCM )
        hi_F = UBOUND( uCM )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

            D  = One
            V1 = Zero
            V2 = Zero
            V3 = Zero

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

            uPM_K(iNX,iPM_D )  = D
            uPM_K(iNX,iPM_V1)  = V1
            uPM_K(iNX,iPM_V2)  = V2
            uPM_K(iNX,iPM_V3)  = V3
            uPM_K(iNX,iPM_E )  = One / ( Gamma_IDEAL - One )
            uPM_K(iNX,iPM_Chi) = Zero

            VSq = V1**2 + V2**2 + V3**2

            VdotB = V1 * CB1 &
                    + V2 * CB2 &
                    + V3 * CB3

            W = One / SQRT( One - VSq )

            uPM_K(iNX,iPM_B1) = W * VdotB * V1 + CB1 / W
            uPM_K(iNX,iPM_B2) = W * VdotB * V2 + CB2 / W
            uPM_K(iNX,iPM_B3) = W * VdotB * V3 + CB3 / W

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPM_K(:,iPM_D), uPM_K(:,iPM_E), uPM_K(:,iPM_Ne), &
                   uAM_K(:,iAM_P) )

          CALL ComputeConserved_MHD &
                 ( uPM_K(:,iPM_D ), uPM_K(:,iPM_V1), uPM_K(:,iPM_V2), &
                   uPM_K(:,iPM_V3), uPM_K(:,iPM_E ), uPM_K(:,iPM_Ne), &
                   uPM_K(:,iPM_B1), uPM_K(:,iPM_B2), uPM_K(:,iPM_B3), &
                   uPM_K(:,iPM_Chi), &
                   uCM_K(:,iCM_D  ), uCM_K(:,iCM_S1), uCM_K(:,iCM_S2), &
                   uCM_K(:,iCM_S3 ), uCM_K(:,iCM_E ), uCM_K(:,iCM_Ne), &
                   uCM_K(:,iCM_B1 ), uCM_K(:,iCM_B2), uCM_K(:,iCM_B3), &
                   uCM_K(:,iCM_Chi), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uGF_K(:,iGF_Alpha   ), &
                   uGF_K(:,iGF_Beta_1  ), &
                   uGF_K(:,iGF_Beta_2  ), &
                   uGF_K(:,iGF_Beta_3  ), &
                   uAM_K(:,iAM_P), &
                   EvolveOnlyMagnetic )

          uCM(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCM_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Cleaning1D


  SUBROUTINE InitializeFields_Cleaning2D( MF_uGF, MF_uCM )

    ! 2D divergence cleaning test from Section 5.2 of
    ! Derigs et al. (2018) with option to use
    ! constant density.

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2, iNX3
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCM_K(nDOFX,nCM)
    REAL(DP)       :: uPM_K(nDOFX,nPM)
    REAL(DP)       :: uAM_K(nDOFX,nAM)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    LOGICAL  :: ConstantDensity
    REAL(DP) :: CB1, CB2, CB3, D, V1, V2, V3, VSq, W, VdotB
    REAL(DP) :: R, R0

    R0 = One / SQRT( 8.0_DP )

    ConstantDensity = .TRUE.
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'ConstantDensity', ConstantDensity )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCM_K = Zero
    uPM_K = Zero
    uAM_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCM )
        hi_F = UBOUND( uCM )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)
            iNX3 = NodeNumberTableX(3,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

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

            IF( R .LE. R0 )THEN

              CB1 = ( One / SQRT( FourPi ) ) &
                    * ( ( R / R0 )**8 - Two * ( R / R0 )**4 + One )

            ELSE

              CB1 = Zero

            END IF

            CB2 = Zero
            CB3 = One / SQRT( FourPi )

            uPM_K(iNX,iPM_D )  = D
            uPM_K(iNX,iPM_V1)  = V1
            uPM_K(iNX,iPM_V2)  = V2
            uPM_K(iNX,iPM_V3)  = V3
            uPM_K(iNX,iPM_E )  = 6.0_DP / ( Gamma_IDEAL - One )
            uPM_K(iNX,iPM_Chi) = Zero

            VSq = V1**2 + V2**2 + V3**2

            VdotB = V1 * CB1 &
                    + V2 * CB2 &
                    + V3 * CB3

            W = One / SQRT( One - VSq )

            uPM_K(iNX,iPM_B1) = W * VdotB * V1 + CB1 / W
            uPM_K(iNX,iPM_B2) = W * VdotB * V2 + CB2 / W
            uPM_K(iNX,iPM_B3) = W * VdotB * V3 + CB3 / W

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPM_K(:,iPM_D), uPM_K(:,iPM_E), uPM_K(:,iPM_Ne), &
                   uAM_K(:,iAM_P) )

          CALL ComputeConserved_MHD &
                 ( uPM_K(:,iPM_D ), uPM_K(:,iPM_V1), uPM_K(:,iPM_V2), &
                   uPM_K(:,iPM_V3), uPM_K(:,iPM_E ), uPM_K(:,iPM_Ne), &
                   uPM_K(:,iPM_B1), uPM_K(:,iPM_B2), uPM_K(:,iPM_B3), &
                   uPM_K(:,iPM_Chi), &
                   uCM_K(:,iCM_D  ), uCM_K(:,iCM_S1), uCM_K(:,iCM_S2), &
                   uCM_K(:,iCM_S3 ), uCM_K(:,iCM_E ), uCM_K(:,iCM_Ne), &
                   uCM_K(:,iCM_B1 ), uCM_K(:,iCM_B2), uCM_K(:,iCM_B3), &
                   uCM_K(:,iCM_Chi), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uGF_K(:,iGF_Alpha   ), &
                   uGF_K(:,iGF_Beta_1  ), &
                   uGF_K(:,iGF_Beta_2  ), &
                   uGF_K(:,iGF_Beta_3  ), &
                   uAM_K(:,iAM_P), &
                   EvolveOnlyMagnetic )

          uCM(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCM_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Cleaning2D


  SUBROUTINE InitializeFields_ShearingDisk( MF_uGF, MF_uCM )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2, iNX3
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCM_K(nDOFX,nCM)
    REAL(DP)       :: uPM_K(nDOFX,nPM)
    REAL(DP)       :: uAM_K(nDOFX,nAM)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    CHARACTER(256) :: FileName

    INTEGER        :: nX_Data
    INTEGER(HID_T) :: FILE_ID

    REAL(DP) :: CB1, CB2, CB3, V1, V2, V3, VSq, W, VdotB
    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), V3Arr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    FileName &
      = "/home/jbuffal/thornado_MHD_3D/Workflow/MHD/ShearingDisk/GR_LR_diffrot.h5"

    ! --- Populate arrays ---

    CALL H5OPEN_F( HDFERR )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    nX_Data = 10000

    ALLOCATE( PressureArr(nX_Data), DensityArr(nX_Data), V3Arr(nX_Data), AlphaArr(nX_Data), &
              PsiArr(nX_Data), X1Arr(nX_Data) )

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

    uGF_K = Zero
    uCM_K = Zero
    uPM_K = Zero
    uAM_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCM )
        hi_F = UBOUND( uCM )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)
            iNX3 = NodeNumberTableX(3,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

            uGF_K(iNX,iGF_Alpha) &
               = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

            uGF_K(iNX,iGF_Psi) &
               = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

            uGF_K(iNX,iGF_h_1) &
              = uGF_K(iNX,iGF_Psi)**2
            uGF_K(iNX,iGF_h_2) &
              = uGF_K(iNX,iGF_Psi)**2
            uGF_K(iNX,iGF_h_3) &
              = uGF_K(iNX,iGF_Psi)**2 * X1

            CALL ComputeGeometryX_FromScaleFactors( uGF_K(:,:) )

            uGF_K(iNX,iGF_Beta_1) = Zero
            uGF_K(iNX,iGF_Beta_2) = Zero
            uGF_K(iNX,iGF_Beta_3) = Zero

            ! --- Fluid Fields ---

            uPM_K(iNX,iPM_D) &
              = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

            V1 = Zero
            V2 = Zero
            V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

            VSq = uGF_K(iNX,iGF_Gm_dd_11) * V1**2 &
                  + uGF_K(iNX,iGF_Gm_dd_22) * V2**2 &
                  + uGF_K(iNX,iGF_Gm_dd_33) * V3**2

            W = One / SQRT( One - VSq )

            uPM_K(iNX,iPM_V1) = V1
            uPM_K(iNX,iPM_V2) = V2
            uPM_K(iNX,iPM_V3) = V3

            uPM_K(iNX,iPM_E) &
              = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
                / ( Gamma_IDEAL - One )

            CB1 = Zero
            CB2 = 2.0 * 1.0d13 * Gauss
            CB3 = Zero

            VdotB = uGF_K(iNX,iGF_Gm_dd_11) * V1 * CB1 &
                    + uGF_K(iNX,iGF_Gm_dd_22) * V2 * CB2 &
                    + uGF_K(iNX,iGF_Gm_dd_33) * V3 * CB3

            uPM_K(iNX,iPM_B1) = W * VdotB * V1 + CB1 / W
            uPM_K(iNX,iPM_B2) = W * VdotB * V2 + CB2 / W
            uPM_K(iNX,iPM_B3) = W * VdotB * V3 + CB3 / W

            uPM_K(iNX,iPM_Chi) = Zero

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPM_K(:,iPM_D), uPM_K(:,iPM_E), uPM_K(:,iPM_Ne), &
                   uAM_K(:,iAM_P) )

          CALL ComputeConserved_MHD &
                 ( uPM_K(:,iPM_D ), uPM_K(:,iPM_V1), uPM_K(:,iPM_V2), &
                   uPM_K(:,iPM_V3), uPM_K(:,iPM_E ), uPM_K(:,iPM_Ne), &
                   uPM_K(:,iPM_B1), uPM_K(:,iPM_B2), uPM_K(:,iPM_B3), &
                   uPM_K(:,iPM_Chi), &
                   uCM_K(:,iCM_D  ), uCM_K(:,iCM_S1), uCM_K(:,iCM_S2), &
                   uCM_K(:,iCM_S3 ), uCM_K(:,iCM_E ), uCM_K(:,iCM_Ne), &
                   uCM_K(:,iCM_B1 ), uCM_K(:,iCM_B2), uCM_K(:,iCM_B3), &
                   uCM_K(:,iCM_Chi), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uGF_K(:,iGF_Alpha   ), &
                   uGF_K(:,iGF_Beta_1  ), &
                   uGF_K(:,iGF_Beta_2  ), &
                   uGF_K(:,iGF_Beta_3  ), &
                   uAM_K(:,iAM_P), &
                   EvolveOnlyMagnetic )

          uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)) &
            = RESHAPE( uGF_K, [ hi_G(4) - lo_G(4) + 1 ] )

          uCM(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCM_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_ShearingDisk


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


END MODULE MF_InitializationModule_Relativistic_IDEAL
