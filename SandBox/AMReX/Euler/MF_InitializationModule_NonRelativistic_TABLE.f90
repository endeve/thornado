MODULE MF_InitializationModule_NonRelativistic_TABLE

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

  USE UnitsModule,             ONLY: &
    Gram,       &
    Kilogram,   &
    Centimeter, &
    Meter,      &
    Kilometer,  &
    Erg,        &
    Joule,      &
    Second,     &
    Kelvin,     &
    SpeedOfLight
  USE ProgramHeaderModule,     ONLY: &
    nDOFX,   &
    nX,      &
    nNodesX, &
    swX,     &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule,              ONLY: &
    MeshType,    &
    CreateMesh,  &
    DestroyMesh, &
    NodeCoordinate
  USE GeometryFieldsModule,    ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule,       ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iPF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iCF_Ne, &
    nAF,    &
    iAF_P,  &
    iAF_E,  &
    iAF_T,  &
    iAF_Ye
  USE Euler_UtilitiesModule,   ONLY: &
    ComputeConserved_Euler
  USE EquationOfStateModule,   ONLY: &
    ComputePressureFromPrimitive,   &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE UtilitiesModule,         ONLY: &
    Locate
  USE Euler_ErrorModule,       ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE MF_KindModule,           ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Pi, &
    TwoPi, &
    FourPi
  USE InputParsingModule,      ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    Gamma_IDEAL, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_NonRelativistic_TABLE


CONTAINS


  SUBROUTINE MF_InitializeFields_NonRelativistic_TABLE &
    ( ProgramName, MF_uGF, MF_uCF, GEOM )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)


    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )

    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection_TABLE' )

        CALL InitializeFields_Advection_TABLE( MF_uGF, MF_uCF )

      CASE( 'Sod_TABLE' )

        CALL InitializeFields_Sod_TABLE( MF_uGF, MF_uCF )

      CASE( 'SphericalSod_TABLE' )

        CALL InitializeFields_SphericalSod_TABLE( MF_uGF, MF_uCF )

      CASE( 'CylindricalSod_TABLE' )

        CALL InitializeFields_CylindricalSod_TABLE( MF_uGF, MF_uCF )

      CASE( 'Jet_TABLE' )

        CALL InitializeFields_Jet_TABLE( MF_uGF, MF_uCF )

      CASE( 'Implosion_TABLE' )

        CALL InitializeFields_Implosion_TABLE( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Advection_TABLE'
          WRITE(*,'(6x,A)')     'Sod_TABLE'
          WRITE(*,'(6x,A)')     'SphericalSod_TABLE'
          WRITE(*,'(6x,A)')     'CylindricalSod_TABLE'
          WRITE(*,'(6x,A)')     'Jet_TABLE'
          WRITE(*,'(6x,A)')     'Implosion_TABLE'

        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

  END SUBROUTINE MF_InitializeFields_NonRelativistic_TABLE


  SUBROUTINE InitializeFields_Advection_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---
    REAL(DP), PARAMETER :: D_0 = 1.0e12_DP * Gram / Centimeter**3
    REAL(DP), PARAMETER :: Amp = 1.0e11_DP * Gram / Centimeter**3
    REAL(DP), PARAMETER :: L   = 1.0e02_DP * Kilometer

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            uPF_K(iNodeX,iPF_D ) &
              = D_0 + Amp * SIN( TwoPi * X1 / L )
            uPF_K(iNodeX,iPF_V1) &
              = 0.1_DP * SpeedOfLight
            uPF_K(iNodeX,iPF_V2) &
              = 0.0_DP * Kilometer / Second
            uPF_K(iNodeX,iPF_V3) &
              = 0.0_DP * Kilometer / Second
            uAF_K(iNodeX,iAF_P ) &
              = 1.0e-2_DP * D_0 * SpeedOfLight**2
            uAF_K(iNodeX,iAF_Ye) &
              = 0.3_DP

          END DO

          CALL ComputeTemperatureFromPressure &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                   uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                   uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                   uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

          CALL ComputeConserved_Euler &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uAF_K(:,iAF_P) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Advection_TABLE


  SUBROUTINE InitializeFields_Sod_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 .LE. Zero ) THEN

              uPF_K(iNodeX,iPF_D ) = 1.0e12_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero      * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_DP

            ELSE

              uPF_K(iNodeX,iPF_D ) = 1.25e11_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e31_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3_DP

            END IF

          END DO


          CALL ComputeTemperatureFromPressure &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                   uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                   uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                   uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

          CALL ComputeConserved_Euler &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uAF_K(:,iAF_P) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Sod_TABLE


  SUBROUTINE InitializeFields_SphericalSod_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 .LE. 5.0_DP * Kilometer ) THEN

              uPF_K(iNodeX,iPF_D ) = 1.0e12_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero      * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_DP

            ELSE

              uPF_K(iNodeX,iPF_D ) = 1.25e11_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e31_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_DP

            END IF

          END DO

          CALL ComputeTemperatureFromPressure &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                   uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                   uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                   uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

          CALL ComputeConserved_Euler &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uAF_K(:,iAF_P) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_SphericalSod_TABLE


  SUBROUTINE InitializeFields_CylindricalSod_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)       :: X1, X2
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF( SQRT( X1**2 + X2**2 ) <= 5.0_DP * Kilometer ) THEN

              uPF_K(iNodeX,iPF_D ) = 1.0e12_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero      * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero      * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_DP

            ELSE

              uPF_K(iNodeX,iPF_D ) = 1.25e11_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e31_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_DP

            END IF

          END DO


          CALL ComputeTemperatureFromPressure &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                   uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                   uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                   uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

          CALL ComputeConserved_Euler &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uAF_K(:,iAF_P) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_CylindricalSod_TABLE


  SUBROUTINE InitializeFields_Jet_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)       :: X1, X2
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF( X1 .LE. Half * Kilometer .AND. X2 .LE. Half * Kilometer )THEN

              !SW
              uPF_K(iNodeX,iPF_D ) = 0.80e12_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3e0_DP

            ELSE IF( X1 .LE. Half * Kilometer &
                       .AND. X2 .GT. Half * Kilometer )THEN
              !NW
              uPF_K(iNodeX,iPF_D ) = 1.0e12_DP  * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = 7.275e4_DP * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3e0_DP

            ELSE IF( X1 .GT. Half * Kilometer &
                       .AND. X2 .GT. Half * Kilometer )THEN

              !NE
              uPF_K(iNodeX,iPF_D ) = 0.5313e12_DP * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero         * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = Zero         * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero         * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 0.4e32_DP    * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3e0_DP

            ELSE

              !SE
              uPF_K(iNodeX,iPF_D ) = 1.0e12_DP  * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = Zero       * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = 7.275e4_DP * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = Zero       * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_DP  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3e0_DP

            END IF

          END DO

        CALL ComputeTemperatureFromPressure &
               ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                 uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                 uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                 uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

        CALL ComputeConserved_Euler &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Jet_TABLE


  SUBROUTINE InitializeFields_Implosion_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)       :: X1, X2
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---
    REAL(DP), PARAMETER :: D_0  = 1.25e13_DP * ( Gram / Centimeter**3 )
    REAL(DP), PARAMETER :: P_0  = 1.0e32_DP  * ( Erg / Centimeter**3 )
    REAL(DP), PARAMETER :: Ye_0 = 1.35e-1_DP
    REAL(DP), PARAMETER :: D_1  = 1.0e14_DP  * ( Gram / Centimeter**3 )
    REAL(DP), PARAMETER :: P_1  = 1.0e33_DP  * ( Erg / Centimeter**3 )
    REAL(DP), PARAMETER :: Ye_1 = 1.5e-1_DP

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF( X1 + X2 .LT. 0.15_DP * Kilometer )THEN

               uPF_K(iNodeX,iPF_D) &
                 = D_0
               uAF_K(iNodeX,iAF_P) &
                 = P_0
               uAF_K(iNodeX,iAF_Ye) &
                 = Ye_0

             ELSE

               uPF_K(iNodeX,iPF_D) &
                 = D_1
               uAF_K(iNodeX,iAF_P) &
                 = P_1
               uAF_K(iNodeX,iAF_Ye) &
                 = Ye_1

             ENDIF

             uPF_K(iNodeX,iPF_V1) &
               = Zero
             uPF_K(iNodeX,iPF_V2) &
               = Zero
             uPF_K(iNodeX,iPF_V3) &
               = Zero

          END DO

          CALL ComputeTemperatureFromPressure &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_P), &
                   uAF_K(:,iAF_Ye), uAF_K(:,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(:,iPF_D ), uAF_K(:,iAF_T ), &
                   uAF_K(:,iAF_Ye), uPF_K(:,iPF_E ), &
                   uAF_K(:,iAF_E ), uPF_K(:,iPF_Ne) )

          CALL ComputeConserved_Euler &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), &
                   uGF_K(:,iGF_Gm_dd_22), &
                   uGF_K(:,iGF_Gm_dd_33), &
                   uAF_K(:,iAF_P) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Implosion_TABLE


END MODULE MF_InitializationModule_NonRelativistic_TABLE
