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
  USE GeometryFieldsModule,    ONLY: &
    nGF,          &
    iGF_Alpha,    &
    iGF_Psi,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
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
    Millisecond
  USE UtilitiesModule,         ONLY: &
    NodeNumberX

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

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_Relativistic_IDEAL


CONTAINS


  SUBROUTINE MF_InitializeFields_Relativistic_IDEAL &
    ( ProgramName, MF_uGF, MF_uCM, GEOM )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )

    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection( MF_uGF, MF_uCM )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Advection'
          WRITE(*,'(6x,A)')     'MagnetizedKH_3D'
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
        .AND. TRIM( AdvectionProfile ) .NE. 'CPAlfvenOblique' &
        .AND. TRIM( AdvectionProfile ) .NE. 'LoopAdvection' )THEN

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
        WRITE(*,'(A)') '  CPAlfvenOblique'
        WRITE(*,'(A)') '  LoopAdvection'
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


END MODULE MF_InitializationModule_Relativistic_IDEAL
