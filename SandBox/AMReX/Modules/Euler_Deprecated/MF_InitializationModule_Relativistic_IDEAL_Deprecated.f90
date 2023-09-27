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
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE
  USE Euler_UtilitiesModule,   ONLY: &
    ComputeConserved_Euler
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
  USE Euler_ErrorModule,       ONLY: &
    DescribeError_Euler

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
    t_end
  USE MF_UtilitiesModule,      ONLY: &
    amrex2thornado_X_Global
  USE MF_AccretionShockUtilitiesModule, ONLY: &
    WriteNodal1DIC_SAS, &
    FileName_Nodal1DIC_SAS, &
    AccretionShockDiagnosticsFileName

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_Relativistic_IDEAL


CONTAINS


  SUBROUTINE MF_InitializeFields_Relativistic_IDEAL &
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

      CASE( 'Advection' )

        CALL InitializeFields_Advection( MF_uGF, MF_uCF )

      CASE( 'RiemannProblem1D' )

        CALL InitializeFields_RiemannProblem1D( MF_uGF, MF_uCF )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D( MF_uGF, MF_uCF )

      CASE( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz' )

        CALL InitializeFields_KelvinHelmholtz( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz3D' )

        CALL InitializeFields_KelvinHelmholtz3D( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock_Relativistic' )

        CALL InitializeFields_StandingAccretionShock_Relativistic &
               ( MF_uGF, MF_uCF, GEOM )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Advection1D'
          WRITE(*,'(6x,A)')     'Advection2D'
          WRITE(*,'(6x,A)')     'Advection3D'
          WRITE(*,'(6x,A)')     'RiemannProblem1D'
          WRITE(*,'(6x,A)')     'RiemannProblem2D'
          WRITE(*,'(6x,A)')     'RiemannProblemSpherical'
          WRITE(*,'(6x,A)')     'KelvinHelmholtz'
          WRITE(*,'(6x,A)')     'KelvinHelmholtz3D'
          WRITE(*,'(6x,A)')     'StandingAccretionShock_Relativistic'
        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

  END SUBROUTINE MF_InitializeFields_Relativistic_IDEAL


  SUBROUTINE InitializeFields_Advection( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2
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
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: AdvectionProfile

    AdvectionProfile = 'SineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( AdvectionProfile ) .NE. 'SineWaveX1' &
        .AND. TRIM( AdvectionProfile ) .NE. 'SineWaveX2' &
        .AND. TRIM( AdvectionProfile ) .NE. 'SineWaveX1X2' )THEN

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for AdvectionProfile: ', &
          TRIM( AdvectionProfile )
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') '  SineWaveX1'
        WRITE(*,'(A)') '  SineWaveX2'
        WRITE(*,'(A)') '  SineWaveX1X2'
        WRITE(*,*)
        WRITE(*,'(A)') 'Stopping...'

      END IF

      CALL DescribeError_Euler( 99 )

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(4x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )

    END IF

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            IF     ( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

              uPF_K(iNX,iPF_D ) = One + 0.1_DP * SIN( TwoPi * X1 )
              uPF_K(iNX,iPF_V1) = 0.1_DP
              uPF_K(iNX,iPF_V2) = Zero
              uPF_K(iNX,iPF_V3) = Zero
              uPF_K(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX2' )THEN

              uPF_K(iNX,iPF_D ) = One + 0.1_DP * SIN( TwoPi * X2 )
              uPF_K(iNX,iPF_V1) = Zero
              uPF_K(iNX,iPF_V2) = 0.1_DP
              uPF_K(iNX,iPF_V3) = Zero
              uPF_K(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1X2' )THEN

              uPF_K(iNX,iPF_D ) &
                = One + 0.1_DP * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
              uPF_K(iNX,iPF_V1) = 0.1_DP * COS( Pi / Four )
              uPF_K(iNX,iPF_V2) = 0.1_DP * SIN( Pi / Four )
              uPF_K(iNX,iPF_V3) = Zero
              uPF_K(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )

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

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_RiemannProblem1D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    INTEGER        :: iDim
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_parmparse)         :: PP
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-Specific Parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(DP)                      :: XD, Vs
    REAL(DP)                      :: LeftState(nPF), RightState(nPF)

    RiemannProblemName = 'Sod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') &
        '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        XD = Half

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.1_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        XD = Half

        Vs = 0.01_DP

        RightState(iPF_D)  = 1.0_DP
        RightState(iPF_V1) = -0.9_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E)  = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs,                 &
                 RightState(iPF_D ), &
                 RightState(iPF_V1), &
                 RightState(iPF_E ) * ( Gamma_IDEAL - One ), &
                 LeftState (iPF_D ), &
                 LeftState (iPF_V1), &
                 LeftState (iPF_E ) )

        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP

      CASE( 'IsolatedContact' )

        Vs = 0.01_DP
        XD = Half

        LeftState(iPF_D ) = 5.9718209694880811e0_DP
        LeftState(iPF_V1) = Vs
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = Vs
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem1' )

        XD = Half

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.9_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 10.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem4' )

        XD = Half

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0e3_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0e-2_DP / ( Gamma_IDEAL - One )

      CASE( 'PerturbedShockTube' )

        XD = Half

        LeftState(iPF_D ) = 5.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 50.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.0_DP ! --- Dummy ---
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 5.0_DP / ( Gamma_IDEAL - One )

      CASE( 'ShockReflection' )

        XD = One

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.99999_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        ! --- All of these are dummies ---
        RightState(iPF_D ) = 0.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.0_DP

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(A,A)') &
            'Invalid choice for RiemannProblemName: ', &
            TRIM( RiemannProblemName )
          WRITE(*,'(A)') 'Valid choices:'
          WRITE(*,'(A)') &
            "  'Sod' - &
            Sod's shock tube"
          WRITE(*,'(A)') &
            "  'MBProblem1' - &
            Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
          WRITE(*,'(A)') &
            "  'MBProblem4' - &
            Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
          WRITE(*,'(A)') &
            "  'PerturbedShockTube' - &
            Del Zanna & Bucciantini (2002) AA, 390, 1177, &
            Sinusoidal density perturbation"
          WRITE(*,'(A)') &
            "  'ShockReflection' - &
            Del Zanna & Bucciantini (2002) AA, 390, 1177, &
            Planar shock reflection"
          WRITE(*,*)
          WRITE(*,'(A)') 'Stopping...'

        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

    IF( amrex_parallel_ioprocessor() )THEN

      IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

        WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
        WRITE(*,*)

      END IF

      WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
      WRITE(*,*)
      WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Right State:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', RightState(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', RightState(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', RightState(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', RightState(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', RightState(iPF_E )
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Left State:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', LeftState(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', LeftState(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', LeftState(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', LeftState(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', LeftState(iPF_E )
      WRITE(*,*)

    END IF

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

            IF( X1 .LE. XD ) THEN

              uPF_K(iNX,iPF_D)  = LeftState(iPF_D )
              uPF_K(iNX,iPF_V1) = LeftState(iPF_V1)
              uPF_K(iNX,iPF_V2) = LeftState(iPF_V2)
              uPF_K(iNX,iPF_V3) = LeftState(iPF_V3)
              uPF_K(iNX,iPF_E)  = LeftState(iPF_E )

            ELSE

              uPF_K(iNX,iPF_D)  = RightState(iPF_D )
              uPF_K(iNX,iPF_V1) = RightState(iPF_V1)
              uPF_K(iNX,iPF_V2) = RightState(iPF_V2)
              uPF_K(iNX,iPF_V3) = RightState(iPF_V3)
              uPF_K(iNX,iPF_E)  = RightState(iPF_E )

              IF( TRIM( RiemannProblemName ) .EQ. 'PerturbedShockTube' ) &
                uPF_K(iNX,iPF_D) &
                  = 2.0_DP + 0.3_DP * SIN( 50.0_DP * X1 )

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

  END SUBROUTINE InitializeFields_RiemannProblem1D


  SUBROUTINE InitializeFields_RiemannProblem2D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2
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
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-specific parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(DP)                      :: X1D, X2D, Vs, V2
    REAL(DP)                      :: NE(nPF), NW(nPF), SE(nPF), SW(nPF)

    RiemannProblemName = 'DzB2002'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') &
        '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'DzB2002' )

        X1D = Half
        X2D = Half

        NE(iPF_D ) = 0.1_DP
        NE(iPF_V1) = 0.0_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        NW(iPF_D ) = 0.1_DP
        NW(iPF_V1) = 0.99_DP
        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP
        NW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SW(iPF_D ) = 0.5_DP
        SW(iPF_V1) = 0.0_DP
        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP
        SW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SE(iPF_D ) = 0.1_DP
        SE(iPF_V1) = 0.0_DP
        SE(iPF_V2) = 0.99_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        X1D = Half
        X2D = Half

        Vs  = 0.01_DP

        NE(iPF_D ) = 1.0_DP
        NE(iPF_V1) = -0.9_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 NE(iPF_D ), &
                 NE(iPF_V1), &
                 NE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 NW(iPF_D ), &
                 NW(iPF_V1), &
                 NW(iPF_E ) )

        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP

        SE(iPF_D ) = 1.0_DP
        SE(iPF_V1) = -0.9_DP
        SE(iPF_V2) = 0.0_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 SE(iPF_D ), &
                 SE(iPF_V1), &
                 SE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 SW(iPF_D ), &
                 SW(iPF_V1), &
                 SW(iPF_E ) )

        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(A,A)') &
            'Invalid choice for RiemannProblemName: ', &
            TRIM( RiemannProblemName )
          WRITE(*,'(A)') 'Valid choices:'
          WRITE(*,'(A)') &
            "  'DzB2002' - &
            Del Zanna & Bucciantini (2002) AA, 390, 1177, Figure 6"
          WRITE(*,'(A)') '  IsolatedShock'
          WRITE(*,*)
          WRITE(*,'(A)') 'Stopping...'

        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

    IF( amrex_parallel_ioprocessor() )THEN

      IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

        WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
        WRITE(*,*)

      END IF

      WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
      WRITE(*,*)
      WRITE(*,'(6x,A,F8.6)') 'X1D = ', X1D
      WRITE(*,'(6x,A,F8.6)') 'X2D = ', X2D
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'NE:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NE(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NE(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NE(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NE(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NE(iPF_E )
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'NW:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NW(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NW(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NW(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NW(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NW(iPF_E )
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'SW:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SW(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SW(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SW(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SW(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SW(iPF_E )
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'SE:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SE(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SE(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SE(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SE(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SE(iPF_E )
      WRITE(*,*)

    END IF

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            ! --- NE ---
            IF     ( X1 .GT. X1D .AND. X2 .GT. X2D )THEN

              uPF_K(iNX,iPF_D ) = NE(iPF_D )
              uPF_K(iNX,iPF_V1) = NE(iPF_V1)
              uPF_K(iNX,iPF_V2) = NE(iPF_V2)
              uPF_K(iNX,iPF_V3) = NE(iPF_V3)
              uPF_K(iNX,iPF_E ) = NE(iPF_E )

            ! --- NW ---
            ELSE IF( X1 .LE. X1D .AND. X2 .GT. X2D )THEN

              uPF_K(iNX,iPF_D ) = NW(iPF_D )
              uPF_K(iNX,iPF_V1) = NW(iPF_V1)
              uPF_K(iNX,iPF_V2) = NW(iPF_V2)
              uPF_K(iNX,iPF_V3) = NW(iPF_V3)
              uPF_K(iNX,iPF_E ) = NW(iPF_E )

            ! --- SW ---
            ELSE IF( X1 .LE. X1D .AND. X2 .LE. X2D )THEN

              uPF_K(iNX,iPF_D ) = SW(iPF_D )
              uPF_K(iNX,iPF_V1) = SW(iPF_V1)
              uPF_K(iNX,iPF_V2) = SW(iPF_V2)
              uPF_K(iNX,iPF_V3) = SW(iPF_V3)
              uPF_K(iNX,iPF_E ) = SW(iPF_E )

            ! --- SE ---
            ELSE

              uPF_K(iNX,iPF_D ) = SE(iPF_D )
              uPF_K(iNX,iPF_V1) = SE(iPF_V1)
              uPF_K(iNX,iPF_V2) = SE(iPF_V2)
              uPF_K(iNX,iPF_V3) = SE(iPF_V3)
              uPF_K(iNX,iPF_E ) = SE(iPF_E )

            END IF

            IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

              ! --- Perturb velocity in X2-direction ---
              CALL RANDOM_NUMBER( V2 )
              uPF_K(iNX,iPF_V2) = 1.0e-13_DP * ( Two * V2 - One )

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

  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_RiemannProblemSpherical( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1
    REAL(DP)       :: X1
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    INTEGER        :: iDim
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---

    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_parmparse)         :: PP
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-Specific Parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(DP)                      :: XD
    REAL(DP)                      :: LeftState(nPF), RightState(nPF)

    RiemannProblemName = 'SphericalSod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') &
        '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'SphericalSod' )

        XD = One

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.1_DP / ( Gamma_IDEAL - One )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(A,A)') &
            'Invalid choice for RiemannProblemName: ', &
            TRIM( RiemannProblemName )
          WRITE(*,'(A)') 'Valid choices:'
          WRITE(*,'(A)') &
            "  'SphericalSod' - &
            Spherical Sod shock tube"
          WRITE(*,*)
          WRITE(*,'(A)') 'Stopping...'

        END IF

        CALL DescribeError_Euler( 99 )

    END SELECT

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
      WRITE(*,*)
      WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Right State:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', RightState(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', RightState(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', RightState(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', RightState(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', RightState(iPF_E )
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Left State:'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', LeftState(iPF_D )
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', LeftState(iPF_V1)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', LeftState(iPF_V2)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', LeftState(iPF_V3)
      WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', LeftState(iPF_E )
      WRITE(*,*)

    END IF

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

            IF( X1 .LE. XD ) THEN

              uPF_K(iNX,iPF_D)  = LeftState(iPF_D )
              uPF_K(iNX,iPF_V1) = LeftState(iPF_V1)
              uPF_K(iNX,iPF_V2) = LeftState(iPF_V2)
              uPF_K(iNX,iPF_V3) = LeftState(iPF_V3)
              uPF_K(iNX,iPF_E)  = LeftState(iPF_E )

            ELSE

              uPF_K(iNX,iPF_D)  = RightState(iPF_D )
              uPF_K(iNX,iPF_V1) = RightState(iPF_V1)
              uPF_K(iNX,iPF_V2) = RightState(iPF_V2)
              uPF_K(iNX,iPF_V3) = RightState(iPF_V3)
              uPF_K(iNX,iPF_E)  = RightState(iPF_E )

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

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Radice & Rezzolla, (2012), AA, 547, A26 ---
  SUBROUTINE InitializeFields_KelvinHelmholtz( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2
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

    REAL(DP) :: a      = 0.01_DP
    REAL(DP) :: Vshear = Half
    REAL(DP) :: A0     = 0.1_DP ! --- Perturbation amplitude ---
    REAL(DP) :: sigma  = 0.1_DP
    REAL(DP) :: rho0   = 0.505_DP
    REAL(DP) :: rho1   = 0.495_DP

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            ! --- V1 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_V1) &
                = +Vshear * TANH( ( X2 - Half ) / a )

            ELSE

              ! --- Paper has a typo here, the minus sign is required ---
              uPF_K(iNX,iPF_V1) &
                = -Vshear * TANH( ( X2 + Half ) / a )

            END IF

            ! --- V2 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_V2) &
                =  A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma**2 ) )

            ELSE

              uPF_K(iNX,iPF_V2) &
                = -A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 + Half )**2 / sigma**2 ) )

            END IF

            ! --- rho ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - Half ) / a )

            ELSE

              uPF_K(iNX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + Half ) / a )

            END IF

            uPF_K(iNX,iPF_V3) = Zero
            uPF_K(iNX,iPF_E)  = One / ( Gamma_IDEAL - One )

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

  END SUBROUTINE InitializeFields_KelvinHelmholtz


  ! --- Relativistic 3D Kelvin-Helmholtz instability a la
  !     Radice & Rezzolla, (2012), AA, 547, A26 ---
  SUBROUTINE InitializeFields_KelvinHelmholtz3D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2
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

    REAL(DP) :: a      = 0.01_DP
    REAL(DP) :: Vshear = Half
    REAL(DP) :: A0     = 0.1_DP ! --- Perturbation amplitude ---
    REAL(DP) :: sigma  = 0.1_DP
    REAL(DP) :: rho0   = 0.505_DP
    REAL(DP) :: rho1   = 0.495_DP
    REAL(DP) :: Vz

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

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            ! --- V1 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_V1) &
                = +Vshear * TANH( ( X2 - Half ) / a )

            ELSE

              ! --- Paper has a typo here, the minus sign is required ---
              uPF_K(iNX,iPF_V1) &
                = -Vshear * TANH( ( X2 + Half ) / a )

            END IF

            ! --- V2 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_V2) &
                =  A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma**2 ) )

            ELSE

              uPF_K(iNX,iPF_V2) &
                = -A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 + Half )**2 / sigma**2 ) )

            END IF

            ! --- rho ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - Half ) / a )

            ELSE

              uPF_K(iNX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + Half ) / a )

            END IF

            CALL RANDOM_NUMBER( Vz )

            uPF_K(iNX,iPF_V3) = 0.01_DP * Vz

            uPF_K(iNX,iPF_E)  = One / ( Gamma_IDEAL - One )

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

  END SUBROUTINE InitializeFields_KelvinHelmholtz3D


  SUBROUTINE InitializeFields_StandingAccretionShock_Relativistic &
    ( MF_uGF, MF_uCF, GEOM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNX, iNX1, iNX2
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
    TYPE(amrex_parmparse)         :: PP

    ! --- Problem-dependent Parameters ---

    REAL(DP) :: MassPNS, ShockRadius, AccretionRate, PolytropicConstant
    LOGICAL  :: ApplyPerturbation
    INTEGER  :: PerturbationOrder
    REAL(DP) :: PerturbationAmplitude
    REAL(DP) :: rPerturbationInner
    REAL(DP) :: rPerturbationOuter

    INTEGER  :: iX1_1, iX1_2, iNX1_1, iNX1_2
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)
    REAL(DP) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_1, P_2
    REAL(DP) :: D0, V0, P0
    REAL(DP) :: Ka, Kb, Mdot, AdvectionTime
    REAL(DP), ALLOCATABLE :: G    (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D    (:,:)
    REAL(DP), ALLOCATABLE :: V    (:,:)
    REAL(DP), ALLOCATABLE :: P    (:,:)
    REAL(DP), ALLOCATABLE :: Alpha(:,:)
    REAL(DP), ALLOCATABLE :: Psi  (:,:)
    LOGICAL               :: InitializeFromFile, ResetEndTime
    INTEGER, PARAMETER    :: nX_LeastSquares = 5

    ApplyPerturbation                 = .FALSE.
    PerturbationOrder                 = 0
    PerturbationAmplitude             = Zero
    rPerturbationInner                = Zero
    rPerturbationOuter                = Zero
    InitializeFromFile                = .FALSE.
    ResetEndTime                      = .FALSE.
    WriteNodal1DIC_SAS                = .FALSE.
    FileName_Nodal1DIC_SAS            = 'Nodal1DIC_SAS.dat'
    AccretionShockDiagnosticsFileName = 'AccretionShockDiagnostics.dat'
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'Mass', &
                        MassPNS )
      CALL PP % get  ( 'AccretionRate', &
                        AccretionRate )
      CALL PP % get  ( 'ShockRadius', &
                        ShockRadius )
      CALL PP % query( 'ApplyPerturbation', &
                        ApplyPerturbation )
      CALL PP % query( 'PerturbationOrder', &
                        PerturbationOrder )
      CALL PP % query( 'PerturbationAmplitude', &
                        PerturbationAmplitude )
      CALL PP % query( 'rPerturbationInner', &
                        rPerturbationInner )
      CALL PP % query( 'rPerturbationOuter', &
                        rPerturbationOuter )
      CALL PP % query( 'InitializeFromFile', &
                        InitializeFromFile )
      CALL PP % query( 'ResetEndTime', &
                        ResetEndTime )
      CALL PP % query( 'WriteNodal1DIC_SAS', &
                        WriteNodal1DIC_SAS )
      CALL PP % query( 'FileName_Nodal1DIC_SAS', &
                        FileName_Nodal1DIC_SAS )
      CALL PP % query( 'AccretionShockDiagnosticsFileName', &
                        AccretionShockDiagnosticsFileName )
    CALL amrex_parmparse_destroy( PP )

    MassPNS            = MassPNS            * SolarMass
    AccretionRate      = AccretionRate      * ( SolarMass / Second )
    ShockRadius        = ShockRadius        * Kilometer
    PolytropicConstant = 2.0e14_DP &
                           * ( Erg / Centimeter**3 &
                           / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) )
    rPerturbationInner = rPerturbationInner * Kilometer
    rPerturbationOuter = rPerturbationOuter * Kilometer

    Mdot = AccretionRate
    Ka   = PolytropicConstant

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'InitializeFromFile:              ', &
        InitializeFromFile

      WRITE(*,'(6x,A,A)') &
        'FileName_Nodal1DIC_SAS:          ', &
        TRIM( FileName_Nodal1DIC_SAS )

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Shock radius:                    ', &
        ShockRadius / Kilometer, &
        ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'PNS Mass:                        ', &
        MassPNS / SolarMass, &
        ' Msun'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Accretion Rate:                  ', &
        AccretionRate / ( SolarMass / Second ), &
        ' Msun/s'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Polytropic Constant (pre-shock): ', &
        Ka / ( ( Erg / Centimeter**3 ) &
          / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) ), &
        ' erg/cm^3 / (g/cm^3)^( Gamma )'

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'Apply Perturbation:           ', &
        ApplyPerturbation

      WRITE(*,'(6x,A,I1)') &
        'Perturbation order:           ', &
        PerturbationOrder

      WRITE(*,'(6x,A,ES9.2E3)') &
        'Perturbation amplitude:       ', &
         PerturbationAmplitude

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Inner radius of perturbation: ', &
        rPerturbationInner / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Outer radius of perturbation: ', &
        rPerturbationOuter / Kilometer, ' km'

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'Reset end-time:  ', &
        ResetEndTime

    END IF

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO

    iX_B1 = [1,1,1] - swX
    iX_E1 = nX      + swX

    ALLOCATE( G    (1:nDOFX     ,iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1:nGF) )
    ALLOCATE( D    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( V    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( P    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Alpha(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Psi  (1:nNodesX(1),iX_B1(1):iX_E1(1)) )

    IF( InitializeFromFile )THEN

      CALL ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )

    ELSE

      ! --- Make local copies of Lapse and Conformal Factor ---

      CALL amrex2thornado_X_Global( GEOM, MF_uGF, nGF, G )

      DO iX1 = iX_B1(1), iX_E1(1)

        Alpha(:,iX1) = G(1:nNodesX(1),iX1,1,1,iGF_Alpha)
        Psi  (:,iX1) = G(1:nNodesX(1),iX1,1,1,iGF_Psi)

      END DO

      ! --- Quantities with _1 are pre-shock, those with _2 are post-shock ---

      CALL LocateFirstUnShockedElement &
             ( iX_B1, iX_E1, ShockRadius, MeshX, &
               iX1_1, iX1_2, iNX1_1, iNX1_2, X1_1, X1_2 )

      ! --- Pre-shock Fields ---

      X1 = NodeCoordinate( MeshX(1), iX_E1(1), nNodesX(1) )

      ! --- Use Newtonian values as initial guesses ---

      V0 = -SQRT( Two * GravitationalConstant * MassPNS / X1 )
      D0 = -Mdot / ( FourPi * X1**2 * V0 )
      P0 = Ka * D0**( Gamma_IDEAL )

      DO iX1 = iX_E1(1), iX1_1, -1

        DO iNX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Ka, Mdot, &
                   Alpha(iNX1,iX1), Psi(iNX1,iX1), D0, V0, P0, &
                   D(iNX1,iX1), V(iNX1,iX1), P(iNX1,iX1) )

          D0 = D(iNX1,iX1)
          V0 = V(iNX1,iX1)
          P0 = P(iNX1,iX1)

        END DO

      END DO

      ! --- Apply Jump Conditions ---

      D_1 = D(iNX1_1,iX1_1)
      V_1 = V(iNX1_1,iX1_1)
      P_1 = P(iNX1_1,iX1_1)

      CALL ApplyJumpConditions_SAS &
             ( Psi(iNX1_1,iX1), D_1, V_1, P_1, D_2, V_2, P_2 )

      Kb = P_2 / D_2**( Gamma_IDEAL )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(6x,A)') 'Jump Conditions'
        WRITE(*,'(6x,A)') '---------------'
        WRITE(*,*)
        WRITE(*,'(8x,A)') 'Pre-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_1
        WRITE(*,'(10x,A,I2.2)')       'iNX1     = ', iNX1_1
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_1 / Kilometer, '  km'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
          D_1 / ( Gram / Centimeter**3 ), '  g/cm^3'
        WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
          V_1 / ( Kilometer / Second ), ' km/s'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
          P_1 / ( Erg / Centimeter**3 ), '  erg/cm^3'
        WRITE(*,*)
        WRITE(*,'(8x,A)') 'Post-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_2
        WRITE(*,'(10x,A,I2.2)')       'iNX1     = ', iNX1_2
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_2 / Kilometer, '  km'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
          D_2 / ( Gram / Centimeter**3 ), '  g/cm^3'
        WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
          V_2 / ( Kilometer / Second ), ' km/s'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
          P_2 / ( Erg / Centimeter**3 ), '  erg/cm^3'
        WRITE(*,*)

      END IF

      ! --- Post-shock Fields ---

      AdvectionTime = Zero

      D0 = D_2
      V0 = V_2
      P0 = P_2

      DO iX1 = iX1_2, iX_B1(1), -1

        DO iNX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          IF( X1 .GT. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Kb, Mdot, &
                   Alpha(iNX1,iX1), Psi(iNX1,iX1), D0, V0, P0, &
                   D(iNX1,iX1), V(iNX1,iX1), P(iNX1,iX1) )

          D0 = D(iNX1,iX1)
          V0 = V(iNX1,iX1)
          P0 = P(iNX1,iX1)

          AdvectionTime &
            = AdvectionTime &
                + WeightsX1(iNX1) * MeshX(1) % Width(iX1) / ABS( V(iNX1,iX1) )

        END DO

      END DO

      IF( ResetEndTime ) &
        t_end = 4.0_DP * AdvectionTime

    END IF

    ! --- Map to 3D domain ---

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

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iX_B(1) = iX_B0(1)
        iX_E(1) = iX_E0(1)

        IF( BX % lo(1) .EQ. 1     ) iX_B(1) = iX_B1(1)
        IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = iX_E1(1)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B (1), iX_E (1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            IF( ApplyPerturbation )THEN

              X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
              X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

              IF( X1 .GE. rPerturbationInner &
                    .AND. X1 .LE. rPerturbationOuter )THEN

                IF( PerturbationOrder .EQ. 0 ) &
                  uPF_K(iNX,iPF_D) &
                    = D(iNX1,iX1) &
                        * ( One + PerturbationAmplitude )

                IF( PerturbationOrder .EQ. 1 ) &
                  uPF_K(iNX,iPF_D) &
                    = D(iNX1,iX1) &
                        * ( One + PerturbationAmplitude * COS( X2 ) )

              ELSE

                uPF_K(iNX,iPF_D) = D(iNX1,iX1)

              END IF

            ELSE

              uPF_K(iNX,iPF_D) = D(iNX1,iX1)

            END IF

            uPF_K(iNX,iPF_V1) = V(iNX1,iX1)
            uPF_K(iNX,iPF_V2) = Zero
            uPF_K(iNX,iPF_V3) = Zero
            uPF_K(iNX,iPF_E ) = P(iNX1,iX1) / ( Gamma_IDEAL - One )

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

    END DO

    DEALLOCATE( Psi     )
    DEALLOCATE( Alpha   )
    DEALLOCATE( P )
    DEALLOCATE( V )
    DEALLOCATE( D )
    DEALLOCATE( G )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    CALL ComputeExtrapolationExponents( MF_uCF, GEOM, nX_LeastSquares )

    IF( amrex_parallel_ioprocessor() )THEN

      IF( .NOT. InitializeFromFile ) &
        WRITE(*,'(6x,A,ES13.6E3,A)') &
          'Advection time:  ', &
          AdvectionTime / Millisecond, ' ms'

      WRITE(*,'(6x,A,I2.2)') &
        'nX_LeastSquares: ', &
        nX_LeastSquares
      WRITE(*,'(6x,A,F8.6)') &
        'ExpD:            ', &
        ExpD
      WRITE(*,'(6x,A,F8.6)') &
        'ExpE:            ', &
        ExpE

    END IF

  END SUBROUTINE InitializeFields_StandingAccretionShock_Relativistic


  ! --- Auxiliary utilities for SAS problem ---


  SUBROUTINE ComputeExtrapolationExponents( MF_uCF, GEOM, nX_LeastSquares )

    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    INTEGER             , INTENT(in) :: nX_LeastSquares

    REAL(DP) :: lnR   (nNodesX(1),1-swX(1):nX(1)+swX(1))
    REAL(DP) :: lnD   (nNodesX(1),1-swX(1):nX(1)+swX(1))
    REAL(DP) :: lnE   (nNodesX(1),1-swX(1):nX(1)+swX(1))
    REAL(DP) :: U(nDOFX          ,1-swX(1):nX(1)+swX(1), &
                                  1-swX(2):nX(2)+swX(2), &
                                  1-swX(3):nX(3)+swX(3),nCF)
    REAL(DP) :: lnR_LS(nNodesX(1),nX_LeastSquares)
    REAL(DP) :: lnD_LS(nNodesX(1),nX_LeastSquares)
    REAL(DP) :: lnE_LS(nNodesX(1),nX_LeastSquares)

    INTEGER  :: iX1, iNX1, iDim
    REAL(DP) :: n

    TYPE(MeshType) :: MeshX(3)

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO

    ! --- Make local copies of X1, D, and tau ---

    CALL amrex2thornado_X_Global( GEOM, MF_uCF, nCF, U )

    DO iX1 = 1-swX(1), nX(1)+swX(1)

      DO iNX1 = 1, nNodesX(1)

        lnR(iNX1,iX1) = LOG( NodeCoordinate( MeshX(1), iX1, iNX1 ) )
        lnD(iNX1,iX1) = U(iNX1,iX1,1,1,iCF_D)
        lnE(iNX1,iX1) = U(iNX1,iX1,1,1,iCF_E)

      END DO

    END DO

    lnD = LOG( lnD )
    lnE = LOG( lnE )

    lnR_LS = lnR(:,1:nX_LeastSquares)
    lnD_LS = lnD(:,1:nX_LeastSquares)
    lnE_LS = lnE(:,1:nX_LeastSquares)

    n = DBLE( nNodesX(1) ) * DBLE( nX_LeastSquares )

    ! --- Expression for exponents from:
    !     https://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html ---

    ExpD = -( n * SUM( lnR_LS * lnD_LS ) - SUM( lnR_LS ) * SUM( lnD_LS ) ) &
             / ( n * SUM( lnR_LS**2 ) - SUM( lnR_LS )**2 )

    ExpE = -( n * SUM( lnR_LS * lnE_LS ) - SUM( lnR_LS ) * SUM( lnE_LS ) ) &
             / ( n * SUM( lnR_LS**2 ) - SUM( lnR_LS )**2 )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE ComputeExtrapolationExponents


  SUBROUTINE ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )

    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(out) :: D(1:,iX_B1(1):), V(1:,iX_B1(1):), P(1:,iX_B1(1):)

    CHARACTER(LEN=16) :: FMT
    INTEGER           :: iX1

    OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) )

    READ(101,*) FMT

    DO iX1 = iX_B1(1), iX_E1(1)

      READ(101,TRIM(FMT)) D(:,iX1)
      READ(101,TRIM(FMT)) V(:,iX1)
      READ(101,TRIM(FMT)) P(:,iX1)

    END DO

    CLOSE( 101 )

  END SUBROUTINE ReadFluidFieldsFromFile


  SUBROUTINE NewtonRaphson_SAS &
    ( X1, MassPNS, K, Mdot, Alpha, Psi, D0, V0, P0, D, V, P )

    REAL(DP), INTENT(in)  :: X1, MassPNS, K, &
                             Mdot, Alpha, Psi, D0, V0, P0
    REAL(DP), INTENT(out) :: D ,V ,P

    REAL(DP) :: W
    REAL(DP) :: Jac(3,3), invJac(3,3)
    REAL(DP) :: f(3), uO(3), uN(3), du(3)

    LOGICAL             :: CONVERGED
    INTEGER             :: ITER
    REAL(DP), PARAMETER :: Tolu = 1.0e-16_DP
    REAL(DP), PARAMETER :: Tolf = 1.0e-16_DP
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolu ) /  LOG( Two ) )

    uO(1) = One
    uO(2) = One
    uO(3) = One

    CONVERGED = .FALSE.
    ITER      = 0
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      W = LorentzFactor( Psi, V0 * uO(2) )

      f(1) &
        = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
            * uO(1) * W * uO(2) + One
      f(2) &
        = Alpha * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                      * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1) ) * W &
            - One
      f(3) &
        = P0 * D0**( -Gamma_IDEAL ) / K * uO(3) * uO(1)**( -Gamma_IDEAL ) - One

      Jac(1,1) = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
                   * W * uO(2)
      Jac(1,2) = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
                   * uO(1) * W**3
      Jac(1,3) = Zero
      Jac(2,1) = -Alpha * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1)**2 * W
      Jac(2,2) = Alpha * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1) ) &
                   * W**3 * Psi**4 * V0**2 / SpeedOfLight**2 * uO(2)
      Jac(2,3) = Alpha * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) / uO(1) * W
      Jac(3,1) = -P0 * D0**( -Gamma_IDEAL ) / K &
                   * Gamma_IDEAL * uO(3) * uO(1)**( -Gamma_IDEAL - One )
      Jac(3,2) = Zero
      Jac(3,3) = P0 * D0**( -Gamma_IDEAL ) / K &
                   * uO(1)**( -Gamma_IDEAL )

      InvJac = Inv3x3( Jac )

      uN = uO - MATMUL( InvJac, f )

      du = uN - uO

      IF( MAXVAL( ABS( du / uO ) ) .LT. Tolu ) CONVERGED = .TRUE.

      uO = uN

    END DO

    D = uN(1) * D0
    V = uN(2) * V0
    P = uN(3) * P0

  END SUBROUTINE NewtonRaphson_SAS


  SUBROUTINE ApplyJumpConditions_SAS( Psi, D_1, V_1, P_1, D_2, V_2, P_2 )

    REAL(DP), INTENT(in)  :: Psi, D_1, V_1, P_1
    REAL(DP), INTENT(out) ::      D_2, V_2, P_2

    REAL(DP) :: C1, C2, C3, A, B, C, D, E
    REAL(DP) :: W_1, h_1
    REAL(DP) :: dx, xa, xb, xc, fa, fb, fc, W_2

    INTEGER             :: ITER
    INTEGER,  PARAMETER :: MAX_ITER = 1000
    REAL(DP), PARAMETER :: TolChi = 1.0e-16_DP

    LOGICAL :: CONVERGED

    W_1 = LorentzFactor( Psi, V_1 )

    h_1 = SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) * P_1 / D_1

    C1 = D_1 * W_1 * V_1
    C2 = D_1 * h_1 * W_1**2 * V_1**2 / SpeedOfLight**2 + Psi**( -4 ) * P_1
    C3 = D_1 * h_1 * W_1**2 * V_1

    A = SpeedOfLight**( -4 ) * ( C3 / C1 )**2 - One
    B = -Two * SpeedOfLight**( 3 ) * C2 * C3 / C1**2 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One ) * Psi**4
    C = SpeedOfLight**( -2 ) * ( C2 / C1 )**2 &
          * ( Gamma_IDEAL / ( Gamma_IDEAL - One ) )**2 * Psi**8 &
          + Two * SpeedOfLight**( -4 ) * ( C3 / C1 )**2 &
          / ( Gamma_IDEAL - One ) * Psi**4 + Psi**4
    D = -Two * SpeedOfLight**( -3 ) * C2 * C3 / C1**2 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One )**2 * Psi**8
    E = SpeedOfLight**( -4 ) * ( C3 / C1 )**2 &
          / ( Gamma_IDEAL - One )**2 * Psi**8

    ! --- Solve with bisection ---

    ! Add 1 km/s to exclude smooth solution
    xa = ( V_1 + One * Kilometer / Second ) / SpeedOfLight

    xb = 1.0e-10_DP * xa

    fa = A * xa**( -2 ) + B * xa**( -1 ) &
             + C + D * xa + E * xa**2
    fb = A * xb**( -2 ) + B * xb**( -1 ) &
             + C + D * xb + E * xb**2

    dx = xb - xa

    ITER      = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      dx = Half * dx

      xc = xa + dx

      fc = A * xc**( -2 ) + B * xc**( -1 ) + C + D * xc + E * xc**2

      IF( fa * fc .LT. Zero )THEN

        xb = xc
        fb = fc

      ELSE IF( fa * fc .GT. Zero )THEN

        xa = xc
        fa = fc

      ELSE

        CONVERGED = .TRUE.

      END IF

      IF( ABS( dx / xa ) .LT. TolChi ) CONVERGED = .TRUE.

    END DO

    V_2 = xc * SpeedOfLight

    W_2 = LorentzFactor( Psi, V_2 )

    D_2 = SpeedOfLight**( -1 ) * ABS( C1 ) * SQRT( xc**( -2 ) - Psi**4 )
    P_2 = ( C3 - D_2 * SpeedOfLight**2 * W_2**2 * V_2 ) &
            / ( Gamma_IDEAL / ( Gamma_IDEAL - One ) * W_2**2 * V_2 )

  END SUBROUTINE ApplyJumpConditions_SAS


  ! --- From: http://fortranwiki.org/fortran/show/Matrix+inversion ---
  FUNCTION Inv3x3( A ) RESULT( invA )

    ! --- Performs a direct calculation of the inverse of a 3×3 matrix ---

    REAL(DP), INTENT(in) :: A   (3,3)
    REAL(DP)             :: invA(3,3)
    REAL(DP)             :: InvDet

    ! --- Calculate the inverse of the determinant of the matrix ---

    InvDet = One / ( A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)     &
                       - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                       + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) )

    ! --- Calculate the inverse of the matrix ---

    invA(1,1) = +InvDet * ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    invA(2,1) = -InvDet * ( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    invA(3,1) = +InvDet * ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    invA(1,2) = -InvDet * ( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    invA(2,2) = +InvDet * ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    invA(3,2) = -InvDet * ( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    invA(1,3) = +InvDet * ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    invA(2,3) = -InvDet * ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    invA(3,3) = +InvDet * ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

    RETURN
  END FUNCTION Inv3x3


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


  SUBROUTINE LocateFirstUnShockedElement &
    ( iX_B1, iX_E1, ShockRadius, MeshX, &
      iX1_1, iX1_2, iNX1_1, iNX1_2, X1_1, X1_2 )

    INTEGER,        INTENT(in)  :: iX_B1(3), iX_E1(3)
    REAL(DP),       INTENT(in)  :: ShockRadius
    TYPE(MeshType), INTENT(in)  :: MeshX(3)
    INTEGER,        INTENT(out) :: iX1_1, iX1_2, iNX1_1, iNX1_2
    REAL(DP),       INTENT(out) :: X1_1, X1_2

    REAL(DP) :: X1, dX1
    INTEGER  :: iX1, iNX1
    LOGICAL  :: FirstPreShockElement = .FALSE.

    X1 = Zero

    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNX1 = 1, nNodesX(1)

        dX1 = NodeCoordinate( MeshX(1), iX1, iNX1 ) - X1
        X1  = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

          iX1_1  = iX1
          iNX1_1 = iNX1
          X1_1   = X1
          X1_2   = X1 - dX1

          IF( iNX1_1 .EQ. 1 )THEN

            iX1_2  = iX1_1 - 1
            iNX1_2 = nNodesX(1)

          ELSE

            iX1_2  = iX1_1
            iNX1_2 = iNX1_1 - 1

          END IF

          FirstPreShockElement = .TRUE.

        END IF

      END DO

    END DO

  END SUBROUTINE LocateFirstUnShockedElement


  ! --- Auxiliary utilities for computine left state ---


  SUBROUTINE ComputeLeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    CALL ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    ! --- Return energy-density instead of pressure ---
    PL = PL / ( Gamma_IDEAL - One )

  END SUBROUTINE ComputeLeftState


  SUBROUTINE ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    REAL(DP), PARAMETER :: EPS = 1.0e-15_DP

    REAL(DP), PARAMETER :: ToldV = EPS
    REAL(DP), PARAMETER :: TolF  = EPS
    INTEGER,  PARAMETER :: nMaxIter = 1000

    INTEGER :: ITERATION
    REAL(DP) :: D, V, P, F
    REAL(DP) :: Vmin, Vmax, Fmin, Fmax, VV, FF

    IF( VR .LT. Zero )THEN

      Vmin = VR   + EPS
      Vmax = +One - EPS

    ELSE

      Vmin = -One + EPS
      Vmax = VR   - EPS

    END IF

    D = Density ( Vs, DR, VR, Vmin )
    P = Pressure( Vs, DR, VR, PR, D, Vmin )
    Fmin = PostShockVelocity( Vs, DR, VR, PR, D, Vmin, P )

    D = Density ( Vs, DR, VR, Vmax )
    P = Pressure( Vs, DR, VR, PR, D, Vmax )
    Fmax = PostShockVelocity( Vs, DR, VR, PR, D, Vmax, P )

    IF( .NOT. Fmin * Fmax .LT. Zero )THEN

      WRITE(*,*) 'Root not bracketed. Stopping...'
      WRITE(*,*) 'Fmin = ', Fmin
      WRITE(*,*) 'Fmax = ', Fmax

      CALL DescribeError_Euler( 10 )

    END IF

    IF( Fmin .GT. Zero )THEN

      VV = Vmax
      FF = Fmax

      Vmax = Vmin
      Vmin = VV

      Fmax = Fmin
      Fmin = FF

    END IF

    ITERATION = 0
    DO WHILE( ITERATION .LT. nMaxIter )

      ITERATION = ITERATION + 1

      V = ( Vmin + Vmax ) / Two

      D = Density ( Vs, DR, VR, V )
      P = Pressure( Vs, DR, VR, PR, D, V )

      F = PostShockVelocity( Vs, DR, VR, PR, D, V, P )

      IF( ABS( V - Vmin ) / MAX( ABS( Vmax ), ABS( Vmin ) ) .LT. ToldV ) EXIT

      IF( F .GT. Zero )THEN

        Vmax = V
        Fmax = F

     ELSE

        Vmin = V
        Fmin = F

     END IF

    END DO

!!$    WRITE(*,*) 'Converged at iteration ', ITERATION
!!$    WRITE(*,*) '|F|:  ' , ABS( F )
!!$    WRITE(*,*) 'dV/V: ', ABS( V - Vmax ) / ABS( Vmax )

    VL = V
    DL = Density ( Vs, DR, VR, VL )
    PL = Pressure( Vs, DR, VR, PR, DL, VL )

  END SUBROUTINE ApplyJumpConditions_LeftState


  REAL(DP) FUNCTION Density( Vs, DR, VR, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, VL

    REAL(DP) :: WR, WL

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    Density = DR * ( WR * ( VR - Vs ) ) / ( WL * ( VL - Vs ) )

    RETURN
  END FUNCTION Density


  REAL(DP) FUNCTION Pressure( Vs, DR, VR, PR, DL, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    Pressure = ( PR * ( One + tau * WR**2 * VR * ( VR - Vs ) ) &
                 - DL * WL**2 * VL**2 + DR * WR**2 * VR**2 &
                 + Vs * ( DL * WL**2 * VL - DR * WR**2 * VR ) ) &
               / ( One + tau * WL**2 * VL * ( VL - Vs ) )

    RETURN
  END FUNCTION Pressure


  REAL(DP) FUNCTION PostShockVelocity( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL, PL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    PostShockVelocity &
      = ( DL + tau * PL ) * WL**2 * ( VL - Vs ) &
          - ( DR + tau * PR ) * WR**2 * ( VR - Vs ) + Vs * ( PL - PR )

    RETURN
  END FUNCTION PostShockVelocity


END MODULE MF_InitializationModule_Relativistic_IDEAL
