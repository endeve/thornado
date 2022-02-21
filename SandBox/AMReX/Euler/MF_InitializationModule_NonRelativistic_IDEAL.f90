MODULE MF_InitializationModule_NonRelativistic_IDEAL

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
  USE UtilitiesModule,         ONLY: &
    NodeNumberX
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
    iAF_P
  USE Euler_UtilitiesModule,   ONLY: &
    ComputeConserved_Euler
  USE EquationOfStateModule,   ONLY: &
    ComputePressureFromPrimitive
  USE UnitsModule,             ONLY: &
    GravitationalConstant, &
    SolarMass,             &
    Second,                &
    Kilometer,             &
    Centimeter,            &
    Erg,                   &
    Gram,                  &
    SpeedOfLight,          &
    Millisecond
  USE Euler_ErrorModule,       ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    Five, &
    Pi, &
    TwoPi, &
    FourPi
  USE InputParsingModule,      ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    Gamma_IDEAL, &
    UseTiling,   &
    t_end
  USE MF_AccretionShockUtilitiesModule, ONLY: &
    WriteNodal1DIC_SAS, &
    FileName_Nodal1DIC_SAS, &
    AccretionShockDiagnosticsFileName

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_NonRelativistic_IDEAL


CONTAINS


  SUBROUTINE MF_InitializeFields_NonRelativistic_IDEAL &
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

      CASE ( 'IsentropicVortex' )

        CALL InitializeFields_IsentropicVortex( MF_uGF, MF_uCF )

      CASE ( 'Sod' )

        CALL InitializeFields_Sod( MF_uGF, MF_uCF )

      CASE ( 'SphericalSod' )

        CALL InitializeFields_SphericalSod( MF_uGF, MF_uCF )

      CASE( 'TopHatAdvection' )

        CALL InitializeFields_TopHatAdvection( MF_uGF, MF_uCF )

      CASE( 'Implosion' )

        CALL InitializeFields_Implosion( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock' )

        CALL InitializeFields_StandingAccretionShock( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock_Units' )

        CALL InitializeFields_StandingAccretionShock_Units( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock_Units_ConstantEntropy' )

        CALL InitializeFields_StandingAccretionShock_Units_ConstantEntropy &
               ( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'IsentropicVortex'
          WRITE(*,'(6x,A)')     'Sod'
          WRITE(*,'(6x,A)')     'SphericalSod'
          WRITE(*,'(6x,A)')     'TopHatAdvection'
          WRITE(*,'(6x,A)')     'Implosion'
          WRITE(*,'(6x,A)')     'StandingAccretionShock'
          WRITE(*,'(6x,A)')     'StandingAccretionShock_Units'
          WRITE(*,'(6x,A)')     'StandingAccretionShock_Units_ConstantEntropy'

        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

  END SUBROUTINE MF_InitializeFields_NonRelativistic_IDEAL


  SUBROUTINE InitializeFields_IsentropicVortex( MF_uGF, MF_uCF )

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
    REAL(DP), PARAMETER :: Beta = 5.0_DP
    REAL(DP)            :: R

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

            R = SQRT( X1**2 + X2**2 )

            uPF_K(iNodeX,iPF_D) &
              = ( One - ( Gamma_IDEAL - One ) * Beta**2 &
                        / ( 8.0_DP * Gamma_IDEAL * Pi**2 ) * EXP( One - R**2 ) &
                )**( One / ( Gamma_IDEAL - One ) )

            uPF_K(iNodeX,iPF_V1) &
              = One - X2 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

            uPF_K(iNodeX,iPF_V2) &
              = One + X1 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

            uPF_K(iNodeX,iPF_V3) &
              = Zero

            uPF_K(iNodeX,iPF_E) &
              = uPF_K(iNodeX,iPF_D)**Gamma_IDEAL / ( Gamma_IDEAL - One )

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_IsentropicVortex


  SUBROUTINE InitializeFields_Sod( MF_uGF, MF_uCF )

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

            IF( X1 .LE. Half ) THEN

              uPF_K(iNodeX,iPF_D)  = One
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = One / (Gamma_IDEAL - One)

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_DP
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = 0.1_DP / (Gamma_IDEAL - One)

            END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_Sod


  SUBROUTINE InitializeFields_SphericalSod( MF_uGF, MF_uCF )

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

            IF( X1 <= One ) THEN

              uPF_K(iNodeX,iPF_D)  = One
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = One / ( Gamma_IDEAL - One )

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_DP
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = 0.1_DP / ( Gamma_IDEAL - One )

            END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_SphericalSod


  SUBROUTINE InitializeFields_TopHatAdvection( MF_uGF, MF_uCF )

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

            IF( ABS( X1 - Half ) .LE. 0.25_DP )THEN

              uPF_K(iNodeX,iPF_D) = Two

            ELSE

              uPF_K(iNodeX,iPF_D) = One

            END IF

          END DO

          uPF_K(:,iPF_V1) = One
          uPF_K(:,iPF_V2) = Zero
          uPF_K(:,iPF_V3) = Zero
          uPF_K(:,iPF_E)  = 0.01_DP

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_TopHatAdvection


  SUBROUTINE InitializeFields_Implosion( MF_uGF, MF_uCF )

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
    REAL(DP), PARAMETER :: Beta = 5.0_DP
    REAL(DP), PARAMETER :: D_0  = 0.125_DP
    REAL(DP), PARAMETER :: E_0  = 0.350_DP
    REAL(DP), PARAMETER :: D_1  = 1.000_DP
    REAL(DP), PARAMETER :: E_1  = 2.500_DP

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

            IF( X1 + X2 .LT. 0.15_DP )THEN

              uPF_K(iNodeX,iPF_D) = D_0
              uPF_K(iNodeX,iPF_E) = E_0

            ELSE

              uPF_K(iNodeX,iPF_D) = D_1
              uPF_K(iNodeX,iPF_E) = E_1

            END IF

            uPF_K(iNodeX,iPF_V1) = Zero
            uPF_K(iNodeX,iPF_V2) = Zero
            uPF_K(iNodeX,iPF_V3) = Zero

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_Implosion


  SUBROUTINE InitializeFields_StandingAccretionShock( MF_uGF, MF_uCF )

    ! --- Initialize standing accretion shock as outlined ---
    ! --- in Blondin et al., (2003), ApJ, 584, 971.       ---

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
    REAL(DP) :: Alpha, Speed, D_Prime, V1_Prime, P_Prime
    REAL(DP) :: mDot, Mass, rShock, Mach

    LOGICAL  :: Perturb
    INTEGER  :: PerturbOrder
    REAL(DP) :: ShellIn, ShellOut, PerturbAmplitude

    TYPE(amrex_parmparse) :: PP

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'mDot',             mDot )
      CALL PP % get( 'Mass',             Mass )
      CALL PP % get( 'rShock',           rShock )
      CALL PP % get( 'Mach',             Mach )
      CALL PP % get( 'Perturb',          Perturb )
      CALL PP % get( 'ShellIn',          ShellIn )
      CALL PP % get( 'ShellOut',         ShellOut )
      CALL PP % get( 'PerturbOrder',     PerturbOrder )
      CALL PP % get( 'PerturbAmplitude', PerturbAmplitude )
    CALL amrex_parmparse_destroy( PP )

    mDot = mDot * Pi

    Alpha = 4.0_DP * Gamma_IDEAL &
              / ( ( Gamma_IDEAL + One ) &
                  * ( Gamma_IDEAL - One ) ) &
              * ( ( Gamma_IDEAL - One ) &
                  / ( Gamma_IDEAL + One ) )**Gamma_IDEAL

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
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

        DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
        DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF( X1 .LE. rShock ) THEN

              CALL ComputeSettlingSpeed_Bisection &
                     ( X1, Alpha, Mass, rShock, Speed )

              uPF_K(iNodeX,iPF_D)  &
                = ( mDot / FourPi ) &
                    * Speed**(-1) * X1**(-2)

              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = Zero
              uPF_K(iNodeX,iPF_V3) &
                = Zero

              V1_prime &
                = ( Gamma_IDEAL - One ) &
                    / ( Gamma_IDEAL + One ) &
                    * SQRT( Two * GravitationalConstant &
                              * Mass / rShock )

              D_prime &
                = ( mDot / FourPi ) * ( One / V1_prime ) &
                    * rShock**(-2)

              P_prime &
                = Two / ( Gamma_IDEAL + One ) &
                    * ( mDot / FourPi ) &
                    * SQRT( Two * GravitationalConstant &
                              * Mass ) * rShock**(-2.5_DP)

              uPF_K(iNodeX,iPF_E)  &
                = P_Prime &
                    * ( uPF_K(iNodeX,iPF_D) / D_prime )**Gamma_IDEAL &
                      / ( Gamma_IDEAL - One )

            ELSE

              Speed &
                = SQRT( One / ( One &
                        + Two &
                        / ( ( Gamma_IDEAL - One) * Mach**2 ) ) / X1 )

              uPF_K(iNodeX,iPF_D)  &
                = mDot / ( FourPi * X1**2 * Speed )
              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = Zero
              uPF_K(iNodeX,iPF_V3) &
                = Zero
              uPF_K(iNodeX,iPF_E)  &
                = uPF_K(iNodeX,iPF_D) / Gamma_IDEAL &
                    * ( Speed / Mach )**2 / ( Gamma_IDEAL - One )

            END IF

            IF( Perturb .AND. ( X1 .GE. ShellIn ) &
                        .AND. ( X1 .LE. ShellOut ) ) THEN

              SELECT CASE( PerturbOrder )

              CASE( 0 )

                uPF_K(iNodeX, iPF_D) &
                  = ( One + PerturbAmplitude ) &
                      * uPF_K(iNodeX,iPF_D)

              CASE( 1 )

                uPF_K(iNodeX, iPF_D) &
                  = ( One + PerturbAmplitude * COS(X2) ) &
                      * uPF_K(iNodeX,iPF_D)

               END SELECT

             END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_StandingAccretionShock


  SUBROUTINE InitializeFields_StandingAccretionShock_Units( MF_uGF, MF_uCF )

    ! --- Initialize standing accretion shock (with physical units) ---
    ! --- as outlined in Blondin et al., (2003), ApJ, 584, 971.     ---

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
    TYPE(amrex_parmparse) :: PP

    ! --- Problem-dependent Parameters ---
    REAL(DP) :: Alpha, Speed, D_Prime, V1_Prime, P_Prime
    REAL(DP) :: mDot, Mass, rShock, Mach

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'mDot',             mDot )
      CALL PP % get( 'Mass',             Mass )
      CALL PP % get( 'rShock',           rShock )
      CALL PP % get( 'Mach',             Mach )
    CALL amrex_parmparse_destroy( PP )

    Mass     = Mass      * SolarMass
    mDot     = mDot      * SolarMass / Second
    rShock   = rShock    * Kilometer

    Alpha = Four * Gamma_IDEAL &
              / ( ( Gamma_IDEAL + One ) &
                  * ( Gamma_IDEAL - One ) ) &
              * ( ( Gamma_IDEAL - One ) &
                  / ( Gamma_IDEAL + One ) )**Gamma_IDEAL

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
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

        DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
        DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF( X1 .LE. rShock ) THEN

              CALL ComputeSettlingSpeed_Bisection &
                ( X1, Alpha, Mass, rShock, Speed )

              uPF_K(iNodeX,iPF_D)  &
                = ( mDot / FourPi ) &
                    * Speed**(-1) * X1**(-2)

              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = Zero
              uPF_K(iNodeX,iPF_V3) &
                = Zero

              V1_prime &
                = ( Gamma_IDEAL - One ) &
                    / ( Gamma_IDEAL + One ) &
                    * SQRT( Two * GravitationalConstant &
                            * Mass / rShock )

              D_prime &
                = ( mDot / FourPi ) * ( One / V1_prime ) &
                    * rShock**(-2)

              P_prime &
                = Two / ( Gamma_IDEAL + One ) &
                    * ( mDot / FourPi ) &
                    * SQRT( Two * GravitationalConstant &
                            * Mass ) * rShock**(-2.5_DP)

              uPF_K(iNodeX,iPF_E)  &
                = P_Prime &
                    * ( uPF_K(iNodeX,iPF_D) / D_prime )**Gamma_IDEAL &
                      / ( Gamma_IDEAL - One )

            ELSE

              Speed &
                = SQRT( Two * GravitationalConstant * Mass &
                        / ( One + Two / ( ( Gamma_IDEAL - One) &
                                          * Mach**2 ) ) / X1 )

              uPF_K(iNodeX,iPF_D)  &
                = mDot / ( FourPi * X1**2 * Speed )
              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = Zero
              uPF_K(iNodeX,iPF_V3) &
                = Zero
              uPF_K(iNodeX,iPF_E)  &
                = uPF_K(iNodeX,iPF_D) / Gamma_IDEAL &
                    * ( Speed / Mach )**2 / ( Gamma_IDEAL - One )

            END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

  END SUBROUTINE InitializeFields_StandingAccretionShock_Units


  SUBROUTINE InitializeFields_StandingAccretionShock_Units_ConstantEntropy( MF_uGF, MF_uCF )

    ! --- Initialize standing accretion shock (with physical units) ---
    ! --- as outlined in Blondin et al., (2003), ApJ, 584, 971.     ---

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
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
    TYPE(amrex_parmparse) :: PP

    ! --- Problem-dependent Parameters ---
    REAL(DP) :: mDot, Mass, rShock, Entropy

    INTEGER  :: iX_B1(3), iX_E1(3)
    REAL(DP) :: BernoulliConstant

    REAL(DP), ALLOCATABLE :: R(:,:), D(:,:), V(:,:), P(:,:)

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'mDot',             mDot )
      CALL PP % get( 'Mass',             Mass )
      CALL PP % get( 'rShock',           rShock )
      CALL PP % get( 'Entropy',          Entropy )
    CALL amrex_parmparse_destroy( PP )

    Mass     = Mass      * SolarMass
    mDot     = mDot      * SolarMass / Second
    rShock   = rShock    * Kilometer
    Entropy  = Entropy

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
               xL(iDim), xR(iDim) )

    END DO

    iX_B1 = [1,1,1] - swX
    iX_E1 = nX      + swX

    ALLOCATE( R(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( D(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( V(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( P(1:nNodesX(1),iX_B1(1):iX_E1(1)) )

    DO iX1 = iX_B1(1), iX_E1(1)
      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R(iNodeX1,iX1) &
          = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

      END DO
    END DO

    BernoulliConstant = 0.0_DP

    CALL ComputeInitialConditions_NewtonRaphson &
           ( R, mDot, Mass, rShock, Entropy, &
             BernoulliConstant, D, V, P )

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

        DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
        DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            uPF_K(iNodeX,iPF_D ) = D(iNodeX1,iX1)
            uPF_K(iNodeX,iPF_V1) = -V(iNodeX1,iX1)
            uPF_K(iNodeX,iPF_E ) = P(iNodeX1,iX1) / ( Gamma_IDEAL - One )

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

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

    DEALLOCATE(R)
    DEALLOCATE(D)
    DEALLOCATE(V)
    DEALLOCATE(P)

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_StandingAccretionShock_Units_ConstantEntropy


  SUBROUTINE ComputeSettlingSpeed_Bisection( r, Alpha, Mass, rShock, V1 )

    REAL(DP), INTENT(in) :: r, Alpha, Mass, rShock

    LOGICAL             :: Converged
    INTEGER             :: Iter
    REAL(DP)            :: a, b, c, ab, F_a, F_b, F_c, F_0
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol_ab  = 1.0e-8_DP
    REAL(DP), PARAMETER :: Tol_F   = 1.0e-8_DP

    REAL(DP), INTENT(out) :: V1

    a = 1.0e-6_DP
    F_a = SettlingSpeedFun( a, r, Alpha, Mass, rShock )

    b = 1.0e-1_DP
    F_b = SettlingSpeedFun( b, r, Alpha, Mass, rShock )

    F_0 = F_a
    ab = b - a

    Converged = .FALSE.
    Iter = 0

    DO WHILE ( .NOT. Converged)

      Iter = Iter + 1

      ab = Half * ab

      c = a + ab

      F_c = SettlingSpeedFun( c, r, Alpha, Mass, rShock )

      IF( F_a * F_c < Zero ) THEN

        b   = c

        F_b = F_c

      ELSE

        a   = c

        F_a = F_c

      END IF

      IF (ab < Tol_ab .AND. ABS( F_a ) / F_0 < Tol_F) Converged = .TRUE.

    END DO

    V1 = a

  END SUBROUTINE ComputeSettlingSpeed_Bisection


  REAL(DP) FUNCTION SettlingSpeedFun( u, r, Alpha, Mass, rShock )

    REAL(DP), INTENT(in) :: u, r, Alpha, Mass, rShock

    SettlingSpeedFun &
      = r * u**2 &
        + Alpha * ( Two * GravitationalConstant &
                      * Mass)**( ( Gamma_Ideal + One ) / Two ) &
        * rShock**( ( Three * Gamma_IDEAL - Five ) / Two ) &
        * r**( Three - Two * Gamma_Ideal ) &
        * u**(One - Gamma_IDEAL) &
        - Two * GravitationalConstant * Mass

    RETURN
  END FUNCTION SettlingSpeedFun


  SUBROUTINE ComputeInitialConditions_NewtonRaphson &
               ( R, mDot, Mass, rShock, EntropyConstant, &
                 BernoulliConstant, D, V, P )

    REAL(DP), INTENT(in)    :: mDot, Mass, rShock
    REAL(DP), INTENT(inout) :: BernoulliConstant
    REAL(DP), INTENT(inout) :: EntropyConstant
    REAL(DP), INTENT(in)    :: R(1:,0:)
    REAL(DP), INTENT(out)   :: D(1:,0:)
    REAL(DP), INTENT(out)   :: V(1:,0:)
    REAL(DP), INTENT(out)   :: P(1:,0:)
    INTEGER, PARAMETER :: iD = 1, iV = 2, iP = 3

    LOGICAL                  :: Converged
    LOGICAL                  :: PostShockGuessed
    INTEGER                  :: iX1
    INTEGER                  :: iX_B1(3), iX_E1(3)
    INTEGER                  :: iNodeX, iNodeX1
    INTEGER                  :: Iter, INFO
    INTEGER,  DIMENSION(3)   :: IPIV
    REAL(DP), PARAMETER      :: Tol = 1.0d-10
    REAL(DP)                 :: D_Previous, V_Previous, P_Previous
    REAL(DP)                 :: E, S
    REAL(DP)                 :: dEdD, dEdP, dSdD, dSdP
    REAL(DP), DIMENSION(3)   :: FVEC, UVEC, dUVEC
    REAL(DP), DIMENSION(3,3) :: FJAC

    iX_B1 = [1,1,1] - swX
    iX_E1 = nX      + swX

    PostShockGuessed = .FALSE.

    DO iX1 = iX_E1(1), iX_B1(1), -1
      DO iNodeX = nDOFX, 1, -1

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        ! --- Initial guess ---

        IF( ( iX1 .EQ. iX_E1(1) ) .AND. ( iNodeX1 .EQ. nNodesX(1) ) )THEN

          UVEC(iD) = mDot / ( FourPi * R(iNodeX1,iX1)**2 &
                          * SQRT( Two * GravitationalConstant &
                                    * Mass / R(iNodeX1,iX1) ) )

          UVEC(iV) = SQRT( Two * GravitationalConstant * Mass &
                             / R(iNodeX1,iX1) )

          UVEC(iP) = EXP( EntropyConstant ) * UVEC(iD)**Gamma_IDEAL

        ELSE IF( ( R(iNodeX1,iX1) .LE. rShock ) .AND. &
                 ( PostShockGuessed .EQV. .FALSE. ) )THEN

          ! Grab values from just above shock.

          D_Previous = UVEC(iD)
          V_Previous = UVEC(iV)
          P_Previous = UVEC(iP)

          ! Jump condition values for the fluid immediately below the shock.

          UVEC(iD) =  ( ( Gamma_IDEAL + One ) / ( Gamma_IDEAL - One ) ) &
                      * ( mDot / FourPi ) &
                      * SQRT( One / ( Two * GravitationalConstant &
                                           * Mass * R(iNodeX1,iX1)**3 ) )

          UVEC(iV) = ( ( Gamma_IDEAL - One ) / ( Gamma_IDEAL + One ) ) &
                     * SQRT( Two * GravitationalConstant * Mass / R(iNodeX1,iX1) )

          UVEC(iP) = ( Two / ( Gamma_IDEAL + One ) * ( mDot / FourPi ) ) &
                     * SQRT( Two * GravitationalConstant * Mass / R(iNodeX1,iX1)**5 )

          CALL ComputeLeftState &
                 ( D_Previous, V_Previous, P_Previous, &
                   UVEC(iD), UVEC(iV), UVEC(iP) )

          ! Set new constants for fluid below shock.

          EntropyConstant = LOG( UVEC(iP) / UVEC(iD)**Gamma_IDEAL )

          E = UVEC(iP) / ( ( Gamma_IDEAL - One ) * UVEC(iD) )

          BernoulliConstant = Half * UVEC(iV)**2 + E + UVEC(iP) / UVEC(iD) &
                                - ( GravitationalConstant * Mass / R(iNodeX1,iX1) )

          PostShockGuessed = .TRUE.

        END IF

        Converged = .FALSE.
        Iter      = 0
        DO WHILE( .NOT. Converged )

          Iter = Iter + 1

          ! Compute specific internal energy and derivatives.

          E = UVEC(iP) / ( ( Gamma_IDEAL - One ) * UVEC(iD) )

          dEdD = -( E / UVEC(iD) )
          dEdP = E / UVEC(iP)

          ! Compute entropy and derivatives.

          S = LOG( UVEC(iP) / ( UVEC(iD)**Gamma_IDEAL) )

          dSdD = -Gamma_IDEAL / UVEC(iD)
          dSdP = One / UVEC(iP)

          CALL ComputeFVEC_InitialConditions &
                 ( R(iNodeX1,iX1), mDot, Mass, BernoulliConstant, &
                   EntropyConstant, UVEC(iD), UVEC(iV), UVEC(iP), S, E, FVEC )

          CALL ComputeFJAC_InitialConditions &
                 ( R(iNodeX1,iX1), UVEC(iD), UVEC(iV), UVEC(iP), E, &
                   dEdD, dEdP, dSdD, dSdP, FJAC )

          dUVEC = FVEC

          CALL DGESV( 3, 1, FJAC, 3, IPIV, dUVEC, 3, INFO )

          !WRITE(*,*)
          !WRITE(*,'(A6,A12,I2.2)') '', 'Iteration: ', Iter
          !WRITE(*,'(A6,A12,4ES20.10E3)') '', '|FVEC| = ', &
          !  ABS( FVEC )
          !WRITE(*,'(A6,A12,4ES20.10E3)') '', '|dU/U| = ', &
          !  ABS(  dUVEC / ( UVEC + TINY(1.0_DP) ) )
          !WRITE(*,*)

          UVEC = UVEC - dUVEC

          IF( ALL( ABS( dUVEC / ( UVEC + TINY(1.0_DP) ) ) < Tol ) ) &
            Converged = .TRUE.

        END DO

        D(iNodeX1,iX1) = UVEC(iD)
        V(iNodeX1,iX1) = UVEC(iV)
        P(iNodeX1,iX1) = UVEC(iP)

      END DO
    END DO

  END SUBROUTINE ComputeInitialConditions_NewtonRaphson


  SUBROUTINE ComputeFVEC_InitialConditions &
               ( r, mDot, Mass, BernoulliConstant, EntropyConstant, &
                 D, V, P, S, E, FVEC )

    REAL(DP), INTENT(in)  :: mDot, Mass, r
    REAL(DP), INTENT(in)  :: BernoulliConstant
    REAL(DP), INTENT(in)  :: EntropyConstant
    REAL(DP), INTENT(in)  :: D, V, P, S, E
    REAL(DP), INTENT(out) :: FVEC(3)

    FVEC(1) = FourPi * r**2 * D * V - mDot
    FVEC(2) = Half * V**2 + E + P / D - Mass / r - BernoulliConstant
    FVEC(3) = S - EntropyConstant

  END SUBROUTINE ComputeFVEC_InitialConditions


  SUBROUTINE ComputeFJAC_InitialConditions &
               ( r, D, V, P, E, &
                 dEdD, dEdP, &
                 dSdD, dSdP, &
                 FJAC )

    REAL(DP), INTENT(in)  :: r, D, V, E, P
    REAL(DP), INTENT(in)  :: dEdD, dEdP
    REAL(DP), INTENT(in)  :: dSdD, dSdP
    REAL(DP), INTENT(out) :: FJAC(3,3)

    FJAC(1,1) = FourPi * r**2 * V
    FJAC(1,2) = FourPi * r**2 * D
    FJAC(1,3) = 0.0_DP

    FJAC(2,1) = ( dEdD - ( P / D**2 ) )
    FJAC(2,2) = V
    FJAC(2,3) = ( dEdP + ( 1.0_DP / D ) )

    FJAC(3,1) = dSdD
    FJAC(3,2) = 0.0_DP
    FJAC(3,3) = dSdP

  END SUBROUTINE ComputeFJAC_InitialConditions


  SUBROUTINE ComputeLeftState( D_R, V_R, P_R, D_L, V_L, P_L )

    LOGICAL                  :: Converged
    INTEGER                  :: Iter, INFO
    INTEGER,  DIMENSION(3)   :: IPIV
    REAL(DP), PARAMETER      :: Tol = 1.0d-10
    REAL(DP), INTENT(inout)  :: D_L, V_L, P_L
    REAL(DP) :: E_L
    REAL(DP), INTENT(inout)  :: D_R, V_R, P_R
    REAL(DP) :: E_R
    REAL(DP) :: dEdD_L, dEdP_L
    REAL(DP), DIMENSION(3)   :: FVEC, UVEC, dUVEC
    REAL(DP), DIMENSION(3,3) :: FJAC

    ! --- Right State ---

    E_R = P_R / ( ( Gamma_IDEAL - One ) * D_R )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Right State:'
    WRITE(*,*)
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_R  =', D_R
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_R  =', V_R
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_R  =', E_R
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_R  =', P_R

    WRITE(*,*)

    E_L = P_L / ( ( Gamma_IDEAL - One ) * D_L )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Left State (Initial Guess):'
    WRITE(*,*)
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_L  =', D_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_L  =', V_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_L  =', E_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_L  =', P_L

    Converged = .FALSE.
    Iter      = 0
    DO WHILE( .NOT. Converged )

      Iter = Iter + 1

      E_L = P_L / ( ( Gamma_IDEAL- One ) * D_L )
      dEdD_L = -( E_L / D_L )
      dEdP_L = E_L / P_L

      CALL ComputeFVEC_LeftState &
             ( D_R, V_R, P_R, E_R, &
               D_L, V_L, P_L, E_L, &
               FVEC )

      CALL ComputeFJAC_LeftState &
             ( D_R, V_R, P_R, E_R, &
               D_L, V_L, P_L, E_L, &
               dEdD_L, dEdP_L, FJAC )

      UVEC  = [ D_L, V_L, P_L ]
      dUVEC = FVEC

      CALL DGESV( 3, 1, FJAC, 3, IPIV, dUVEC, 3, INFO )

      !WRITE(*,*)
      !WRITE(*,'(A6,A12,I2.2)') '', 'Iteration: ', Iter
      !WRITE(*,'(A6,A12,4ES20.10E3)') '', '|FVEC| = ', &
      !  ABS( FVEC )
      !WRITE(*,'(A6,A12,4ES20.10E3)') '', '|dU/U| = ', &
      !  ABS(  dUVEC / ( UVEC + TINY(1.0_DP) ) )
      !WRITE(*,*)

      D_L = UVEC(1) - dUVEC(1)
      V_L = UVEC(2) - dUVEC(2)
      P_L = UVEC(3) - dUVEC(3)

      IF( ALL( ABS(  dUVEC / ( UVEC + TINY(1.0_DP) ) ) < Tol ) ) &
        Converged = .TRUE.

    END DO

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Left State:'
    WRITE(*,*)
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_L  =', D_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_L  =', V_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_L  =', E_L
    WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_L  =', P_L
    WRITE(*,*)

  END SUBROUTINE ComputeLeftState


  SUBROUTINE ComputeFVEC_LeftState &
               ( D_R, V_R, P_R, E_R, D_L, V_L, P_L, E_L, FVEC )

    REAL(DP),               INTENT(in)  :: D_R, V_R, P_R, E_R
    REAL(DP),               INTENT(in)  :: D_L, V_L, P_L, E_L
    REAL(DP), DIMENSION(3), INTENT(out) :: FVEC

    FVEC(1) = D_L * V_L - D_R * V_R
    FVEC(2) = D_L * V_L**2 - D_R * V_R**2 + P_L - P_R
    FVEC(3) = E_L - E_R + Half * ( V_L**2 - V_R**2 ) + P_L / D_L - P_R / D_R

  END SUBROUTINE ComputeFVEC_LeftState


  SUBROUTINE ComputeFJAC_LeftState &
               ( D_R, V_R, P_R, E_R, D_L, V_L, P_L, E_L, &
                 dEdD_L, dEdP_L, FJAC )

    REAL(DP),                 INTENT(in)  :: D_R, V_R, P_R, E_R
    REAL(DP),                 INTENT(in)  :: D_L, V_L, P_L, E_L
    REAL(DP),                 INTENT(in)  :: dEdD_L, dEdP_L
    REAL(DP), DIMENSION(3,3), INTENT(out) :: FJAC

    FJAC(1,1) = V_L
    FJAC(1,2) = D_L
    FJAC(1,3) = 0.0_DP

    FJAC(2,1) = V_L**2
    FJAC(2,2) = 2.0_DP * D_L * V_L
    FJAC(2,3) = One

    FJAC(3,1) = dEdD_L - ( P_L / D_L**2 )
    FJAC(3,2) = V_L
    FJAC(3,3) = dEdP_L + ( One / D_L )

  END SUBROUTINE ComputeFJAC_LeftState


END MODULE MF_InitializationModule_NonRelativistic_IDEAL
