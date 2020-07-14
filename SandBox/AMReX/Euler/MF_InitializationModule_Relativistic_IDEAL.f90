MODULE MF_InitializationModule_Relativistic_IDEAL

  ! --- AMReX Modules ---

  USE amrex_fort_module,       ONLY: &
    AR => amrex_real
  USE amrex_box_module,        ONLY: &
    amrex_box
  USE amrex_multifab_module,   ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,   ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_max
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
    NodeNumberTableX
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
    GravitationalConstant
  USE UtilitiesModule,         ONLY: &
    NodeNumberX
  USE Euler_ErrorModule,       ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE InputParsingModule,      ONLY: &
    nLevels,            &
    xL,                 &
    xR,                 &
    Gamma_IDEAL,        &
    InitializeFromFile, &
    NodalDataFileNameBase

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_Relativistic_IDEAL

  REAL(AR), PARAMETER :: Zero     = 0.0_AR
  REAL(AR), PARAMETER :: SqrtTiny = SQRT( TINY( 1.0_AR ) )
  REAL(AR), PARAMETER :: Half     = 0.5_AR
  REAL(AR), PARAMETER :: One      = 1.0_AR
  REAL(AR), PARAMETER :: Two      = 2.0_AR
  REAL(AR), PARAMETER :: Three    = 3.0_AR
  REAL(AR), PARAMETER :: Four     = 4.0_AR
  REAL(AR), PARAMETER :: Pi       = ACOS( -1.0_AR )
  REAL(AR), PARAMETER :: TwoPi    = 2.0_AR * Pi
  REAL(AR), PARAMETER :: FourPi   = 4.0_AR * Pi

  LOGICAL, PUBLIC :: WriteNodalData_SAS


CONTAINS


  SUBROUTINE MF_InitializeFields_Relativistic_IDEAL &
    ( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )

    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection1D' )

        CALL InitializeFields_Advection1D( MF_uGF, MF_uCF )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D( MF_uGF, MF_uCF )

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
               ( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Advection1D'
          WRITE(*,'(6x,A)')     'Advection2D'
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


  SUBROUTINE InitializeFields_Advection1D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(AR)       :: X1
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    INTEGER        :: iDim
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---
    CHARACTER(LEN=:), ALLOCATABLE :: AdvectionProfile

    AdvectionProfile = 'SineWave'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( AdvectionProfile ) .NE. 'SineWave' )THEN

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for AdvectionProfile: ', &
          TRIM( AdvectionProfile )
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') '  SineWave'
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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            IF     ( TRIM( AdvectionProfile ) .EQ. 'SineWave' )THEN

              uPF_K(iNodeX,iPF_D ) = One + 0.1_AR * SIN( TwoPi * X1 )
              uPF_K(iNodeX,iPF_V1) = 0.1_AR
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )

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

  END SUBROUTINE InitializeFields_Advection1D


  SUBROUTINE InitializeFields_Advection2D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(AR)       :: X1, X2
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            IF     ( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

              uPF_K(iNodeX,iPF_D ) = One + 0.1_AR * SIN( TwoPi * X1 )
              uPF_K(iNodeX,iPF_V1) = 0.1_AR
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX2' )THEN

              uPF_K(iNodeX,iPF_D ) = One + 0.1_AR * SIN( TwoPi * X2 )
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = 0.1_AR
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )

            ELSE IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1X2' )THEN

              uPF_K(iNodeX,iPF_D ) &
                = One + 0.1_AR * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
              uPF_K(iNodeX,iPF_V1) = 0.1_AR * COS( Pi / Four )
              uPF_K(iNodeX,iPF_V2) = 0.1_AR * SIN( Pi / Four )
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )

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

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_RiemannProblem1D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(AR)       :: X1
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    INTEGER        :: iDim
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_parmparse)         :: PP
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-Specific Parameters ---
    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(AR)                      :: XD, Vs
    REAL(AR)                      :: LeftState(nPF), RightState(nPF)

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

        LeftState(iPF_D ) = 1.0_AR
        LeftState(iPF_V1) = 0.0_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_AR
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 0.1_AR / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        XD = Half

        Vs = 0.01_AR

        RightState(iPF_D)  = 1.0_AR
        RightState(iPF_V1) = -0.9_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E)  = 1.0_AR / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs,                 &
                 RightState(iPF_D ), &
                 RightState(iPF_V1), &
                 RightState(iPF_E ) * ( Gamma_IDEAL - One ), &
                 LeftState (iPF_D ), &
                 LeftState (iPF_V1), &
                 LeftState (iPF_E ) )

        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR

      CASE( 'IsolatedContact' )

        Vs = 0.01_AR
        XD = Half

        LeftState(iPF_D ) = 5.9718209694880811e0_AR
        LeftState(iPF_V1) = Vs
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_AR
        RightState(iPF_V1) = Vs
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

      CASE( 'MBProblem1' )

        XD = Half

        LeftState(iPF_D ) = 1.0_AR
        LeftState(iPF_V1) = 0.9_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_AR
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 10.0_AR / ( Gamma_IDEAL - One )

      CASE( 'MBProblem4' )

        XD = Half

        LeftState(iPF_D ) = 1.0_AR
        LeftState(iPF_V1) = 0.0_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 1.0e3_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_AR
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 1.0e-2_AR / ( Gamma_IDEAL - One )

      CASE( 'PerturbedShockTube' )

        XD = Half

        LeftState(iPF_D ) = 5.0_AR
        LeftState(iPF_V1) = 0.0_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 50.0_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.0_AR ! --- Dummy ---
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 5.0_AR / ( Gamma_IDEAL - One )

      CASE( 'ShockReflection' )

        XD = One

        LeftState(iPF_D ) = 1.0_AR
        LeftState(iPF_V1) = 0.99999_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 0.01_AR / ( Gamma_IDEAL - One )

        ! --- All of these are dummies ---
        RightState(iPF_D ) = 0.0_AR
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 0.0_AR

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            IF( X1 .LE. XD ) THEN

              uPF_K(iNodeX,iPF_D)  = LeftState(iPF_D )
              uPF_K(iNodeX,iPF_V1) = LeftState(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = LeftState(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = LeftState(iPF_V3)
              uPF_K(iNodeX,iPF_E)  = LeftState(iPF_E )

            ELSE

              uPF_K(iNodeX,iPF_D)  = RightState(iPF_D )
              uPF_K(iNodeX,iPF_V1) = RightState(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = RightState(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = RightState(iPF_V3)
              uPF_K(iNodeX,iPF_E)  = RightState(iPF_E )

              IF( TRIM( RiemannProblemName ) .EQ. 'PerturbedShockTube' ) &
                uPF_k(iNodeX,iPF_D) &
                  = 2.0_AR + 0.3_AR * SIN( 50.0_AR * X1 )

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
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(AR)       :: X1, X2
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-specific parameters ---
    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(AR)                      :: X1D, X2D, Vs, V2
    REAL(AR)                      :: NE(nPF), NW(nPF), SE(nPF), SW(nPF)

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

        NE(iPF_D ) = 0.1_AR
        NE(iPF_V1) = 0.0_AR
        NE(iPF_V2) = 0.0_AR
        NE(iPF_V3) = 0.0_AR
        NE(iPF_E ) = 0.01_AR / ( Gamma_IDEAL - One )

        NW(iPF_D ) = 0.1_AR
        NW(iPF_V1) = 0.99_AR
        NW(iPF_V2) = 0.0_AR
        NW(iPF_V3) = 0.0_AR
        NW(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        SW(iPF_D ) = 0.5_AR
        SW(iPF_V1) = 0.0_AR
        SW(iPF_V2) = 0.0_AR
        SW(iPF_V3) = 0.0_AR
        SW(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        SE(iPF_D ) = 0.1_AR
        SE(iPF_V1) = 0.0_AR
        SE(iPF_V2) = 0.99_AR
        SE(iPF_V3) = 0.0_AR
        SE(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        X1D = Half
        X2D = Half

        Vs  = 0.01_AR

        NE(iPF_D ) = 1.0_AR
        NE(iPF_V1) = -0.9_AR
        NE(iPF_V2) = 0.0_AR
        NE(iPF_V3) = 0.0_AR
        NE(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 NE(iPF_D ), &
                 NE(iPF_V1), &
                 NE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 NW(iPF_D ), &
                 NW(iPF_V1), &
                 NW(iPF_E ) )

        NW(iPF_V2) = 0.0_AR
        NW(iPF_V3) = 0.0_AR

        SE(iPF_D ) = 1.0_AR
        SE(iPF_V1) = -0.9_AR
        SE(iPF_V2) = 0.0_AR
        SE(iPF_V3) = 0.0_AR
        SE(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 SE(iPF_D ), &
                 SE(iPF_V1), &
                 SE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 SW(iPF_D ), &
                 SW(iPF_V1), &
                 SW(iPF_E ) )

        SW(iPF_V2) = 0.0_AR
        SW(iPF_V3) = 0.0_AR

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            ! --- NE ---
            IF     ( X1 .GT. X1D .AND. X2 .GT. X2D )THEN

              uPF_K(iNodeX,iPF_D ) = NE(iPF_D )
              uPF_K(iNodeX,iPF_V1) = NE(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = NE(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = NE(iPF_V3)
              uPF_K(iNodeX,iPF_E ) = NE(iPF_E )

            ! --- NW ---
            ELSE IF( X1 .LE. X1D .AND. X2 .GT. X2D )THEN

              uPF_K(iNodeX,iPF_D ) = NW(iPF_D )
              uPF_K(iNodeX,iPF_V1) = NW(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = NW(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = NW(iPF_V3)
              uPF_K(iNodeX,iPF_E ) = NW(iPF_E )

            ! --- SW ---
            ELSE IF( X1 .LE. X1D .AND. X2 .LE. X2D )THEN

              uPF_K(iNodeX,iPF_D ) = SW(iPF_D )
              uPF_K(iNodeX,iPF_V1) = SW(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = SW(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = SW(iPF_V3)
              uPF_K(iNodeX,iPF_E ) = SW(iPF_E )

            ! --- SE ---
            ELSE

              uPF_K(iNodeX,iPF_D ) = SE(iPF_D )
              uPF_K(iNodeX,iPF_V1) = SE(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = SE(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = SE(iPF_V3)
              uPF_K(iNodeX,iPF_E ) = SE(iPF_E )

            END IF

            IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

              ! --- Perturb velocity in X2-direction ---
              CALL RANDOM_NUMBER( V2 )
              uPF_K(iNodeX,iPF_V2) = 1.0e-13_AR * ( Two * V2 - One )

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
    INTEGER        :: iNodeX, iNodeX1
    REAL(AR)       :: X1
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    INTEGER        :: iDim
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_parmparse)         :: PP
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-Specific Parameters ---
    CHARACTER(LEN=:), ALLOCATABLE :: RiemannProblemName
    REAL(AR)                      :: XD
    REAL(AR)                      :: LeftState(nPF), RightState(nPF)

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

        LeftState(iPF_D ) = 1.0_AR
        LeftState(iPF_V1) = 0.0_AR
        LeftState(iPF_V2) = 0.0_AR
        LeftState(iPF_V3) = 0.0_AR
        LeftState(iPF_E ) = 1.0_AR / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_AR
        RightState(iPF_V1) = 0.0_AR
        RightState(iPF_V2) = 0.0_AR
        RightState(iPF_V3) = 0.0_AR
        RightState(iPF_E ) = 0.1_AR / ( Gamma_IDEAL - One )

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            IF( X1 .LE. XD ) THEN

              uPF_K(iNodeX,iPF_D)  = LeftState(iPF_D )
              uPF_K(iNodeX,iPF_V1) = LeftState(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = LeftState(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = LeftState(iPF_V3)
              uPF_K(iNodeX,iPF_E)  = LeftState(iPF_E )

            ELSE

              uPF_K(iNodeX,iPF_D)  = RightState(iPF_D )
              uPF_K(iNodeX,iPF_V1) = RightState(iPF_V1)
              uPF_K(iNodeX,iPF_V2) = RightState(iPF_V2)
              uPF_K(iNodeX,iPF_V3) = RightState(iPF_V3)
              uPF_K(iNodeX,iPF_E)  = RightState(iPF_E )

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
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(AR)       :: X1, X2
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---
    REAL(AR) :: a      = 0.01_AR
    REAL(AR) :: Vshear = Half
    REAL(AR) :: A0     = 0.1_AR ! --- Perturbation amplitude ---
    REAL(AR) :: sigma  = 0.1_AR
    REAL(AR) :: rho0   = 0.505_AR
    REAL(AR) :: rho1   = 0.495_AR

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            ! --- V1 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_V1) &
                = +Vshear * TANH( ( X2 - Half ) / a )

            ELSE

              ! --- Paper has a typo here, the minus sign is required ---
              uPF_K(iNodeX,iPF_V1) &
                = -Vshear * TANH( ( X2 + Half ) / a )

            END IF

            ! --- V2 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_V2) &
                =  A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma ) )

            ELSE

              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 + Half )**2 / sigma ) )

            END IF

            ! --- rho ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - Half ) / a )

            ELSE

              uPF_K(iNodeX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + Half ) / a )

            END IF

            uPF_K(iNodeX,iPF_V3) = Zero
            uPF_K(iNodeX,iPF_E)  = One / ( Gamma_IDEAL - One )

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
    INTEGER        :: iNodeX, iNodeX1, iNodeX2
    REAL(AR)       :: X1, X2
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent Parameters ---
    REAL(AR) :: a      = 0.01_AR
    REAL(AR) :: Vshear = Half
    REAL(AR) :: A0     = 0.1_AR ! --- Perturbation amplitude ---
    REAL(AR) :: sigma  = 0.1_AR
    REAL(AR) :: rho0   = 0.505_AR
    REAL(AR) :: rho1   = 0.495_AR
    REAL(AR) :: Vz

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

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

            ! --- V1 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_V1) &
                = +Vshear * TANH( ( X2 - Half ) / a )

            ELSE

              ! --- Paper has a typo here, the minus sign is required ---
              uPF_K(iNodeX,iPF_V1) &
                = -Vshear * TANH( ( X2 + Half ) / a )

            END IF

            ! --- V2 ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_V2) &
                =  A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma ) )

            ELSE

              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( TwoPi * X1 ) &
                    * EXP( -( ( X2 + Half )**2 / sigma ) )

            END IF

            ! --- rho ---
            IF( X2 .GT. Zero )THEN

              uPF_K(iNodeX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - Half ) / a )

            ELSE

              uPF_K(iNodeX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + Half ) / a )

            END IF

            CALL RANDOM_NUMBER( Vz )

            uPF_K(iNodeX,iPF_V3) = 0.01_AR * Vz

            uPF_K(iNodeX,iPF_E)  = One / ( Gamma_IDEAL - One )

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
    ( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2
    REAL(AR)       :: uGF_K(nDOFX,nGF)
    REAL(AR)       :: uCF_K(nDOFX,nCF)
    REAL(AR)       :: uPF_K(nDOFX,nPF)
    REAL(AR)       :: uAF_K(nDOFX,nAF)
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    TYPE(amrex_parmparse)         :: PP

    ! --- Problem-dependent Parameters ---
    REAL(AR) :: MassPNS, ShockRadius, AccretionRate, PolytropicConstant
    LOGICAL  :: ApplyPerturbation
    INTEGER  :: PerturbationOrder
    REAL(AR) :: PerturbationAmplitude
    REAL(AR) :: rPerturbationInner
    REAL(AR) :: rPerturbationOuter

    INTEGER  :: iX1_1, iX1_2, iNodeX1_1, iNodeX1_2
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)
    REAL(AR) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_1, P_2, C1, C2, C3
    REAL(AR) :: D0, V0, P0
    REAL(AR) :: W, dX1, Ka, Kb, Mdot
    REAL(AR), ALLOCATABLE :: D     (:,:)
    REAL(AR), ALLOCATABLE :: V     (:,:)
    REAL(AR), ALLOCATABLE :: P     (:,:)
    REAL(AR), ALLOCATABLE :: Alpha (:,:)
    REAL(AR), ALLOCATABLE :: Psi   (:,:)
    REAL(AR), ALLOCATABLE :: Alpha0(:,:,:)
    REAL(AR), ALLOCATABLE :: Psi0  (:,:,:)
    LOGICAL  :: FirstPreShockElement = .FALSE.

    ! --- Quantities with (1) are pre-shock, those with (2) are post-shock ---

    ApplyPerturbation     = .FALSE.
    PerturbationOrder     = 0
    PerturbationAmplitude = Zero
    rPerturbationInner    = Zero
    rPerturbationOuter    = Zero
    WriteNodalData_SAS    = .FALSE.
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'Mass'                 , MassPNS               )
      CALL PP % get  ( 'AccretionRate'        , AccretionRate         )
      CALL PP % get  ( 'ShockRadius'          , ShockRadius           )
      CALL PP % get  ( 'PolytropicConstant'   , PolytropicConstant    )
      CALL PP % query( 'ApplyPerturbation'    , ApplyPerturbation     )
      CALL PP % query( 'PerturbationOrder'    , PerturbationOrder     )
      CALL PP % query( 'PerturbationAmplitude', PerturbationAmplitude )
      CALL PP % query( 'rPerturbationInner'   , rPerturbationInner    )
      CALL PP % query( 'rPerturbationOuter'   , rPerturbationOuter    )
      CALL PP % query( 'WriteNodalData_SAS'   , WriteNodalData_SAS    )
    CALL amrex_parmparse_destroy( PP )

    MassPNS            = MassPNS            * SolarMass
    AccretionRate      = AccretionRate      * ( SolarMass / Second )
    ShockRadius        = ShockRadius        * Kilometer
    PolytropicConstant = PolytropicConstant &
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

    ALLOCATE( D     (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( V     (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( P     (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Alpha (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Psi   (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Alpha0(1:nNodesX(1),iX_B1(1):iX_E1(1),0:nLevels-1) )
    ALLOCATE( Psi0  (1:nNodesX(1),iX_B1(1):iX_E1(1),0:nLevels-1) )

    IF( InitializeFromFile )THEN

      CALL ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )

    ELSE

      ! --- Make local copies of Lapse and Conformal Factor ---

      Alpha0 = Zero
      Psi0   = Zero

      DO iLevel = 0, nLevels-1

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

        DO WHILE( MFI % next() )

          uGF => MF_uGF(iLevel) % DataPtr( MFI )

          BX = MFI % tilebox()

          lo_G = LBOUND( uGF )
          hi_G = UBOUND( uGF )

          iX_B = BX % lo
          iX_E = BX % hi

          ! -- Get physical ghost cells right ---

          IF( BX % lo(1) .EQ. 1     ) iX_B(1) = iX_B(1) - swX(1)
          IF( BX % lo(2) .EQ. 1     ) iX_B(2) = iX_B(2) - swX(2)
          IF( BX % lo(3) .EQ. 1     ) iX_B(3) = iX_B(3) - swX(3)
          IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = iX_E(1) + swX(1)
          IF( BX % hi(2) .EQ. nX(2) ) iX_E(2) = iX_E(2) + swX(2)
          IF( BX % hi(3) .EQ. nX(3) ) iX_E(3) = iX_E(3) + swX(3)

          DO iX3 = iX_B(3), iX_E(3)
          DO iX2 = iX_B(2), iX_E(2)
          DO iX1 = iX_B(1), iX_E(1)

            uGF_K &
              = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

            DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX(iNodeX1,iNodeX2,iNodeX3)

              Alpha0(iNodeX1,iX1,iLevel) = uGF_K(iNodeX,iGF_Alpha)
              Psi0  (iNodeX1,iX1,iLevel) = uGF_K(iNodeX,iGF_Psi)

            END DO
            END DO
            END DO

          END DO
          END DO
          END DO

        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO

      ! --- Combine data from different grids ---

      DO iX1 = iX_B1(1), iX_E1(1)

        DO iNodeX1 = 1, nNodesX(1)

          CALL amrex_parallel_reduce_max( Alpha0(iNodeX1,iX1,:), nLevels )
          CALL amrex_parallel_reduce_max( Psi0  (iNodeX1,iX1,:), nLevels )

          Alpha(iNodeX1,iX1) = Alpha0(iNodeX1,iX1,0)
          Psi  (iNodeX1,iX1) = Psi0  (iNodeX1,iX1,0)

        END DO

      END DO

      ! --- Locate first element with un-shocked fluid ---

      X1 = Zero

      DO iX1 = iX_B1(1), iX_E1(1)

        DO iNodeX1 = 1, nNodesX(1)

          dX1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) - X1
          X1  = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

            iX1_1     = iX1
            iNodeX1_1 = iNodeX1
            X1_1      = X1
            X1_2      = X1 - dX1

            IF( iNodeX1_1 .EQ. 1 )THEN

              iX1_2     = iX1_1 - 1
              iNodeX1_2 = nNodesX(1)

            ELSE

              iX1_2     = iX1_1
              iNodeX1_2 = iNodeX1_1 - 1

            END IF

            FirstPreShockElement = .TRUE.

          END IF

        END DO

      END DO

      ! --- Pre-shock Fields ---

      X1 = NodeCoordinate( MeshX(1), iX_E1(1), nNodesX(1) )

      ! --- Use Newtonian values as initial guesses ---

      V0 = -SQRT( Two * GravitationalConstant * MassPNS / X1 )
      D0 = -Mdot / ( FourPi * X1**2 * V0 )
      P0 = Ka * D0**( Gamma_IDEAL )

      DO iX1 = iX_E1(1), iX1_1, -1

        DO iNodeX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Ka, Mdot, &
                   Alpha(iNodeX1,iX1), Psi(iNodeX1,iX1), D0, V0, P0, &
                   D(iNodeX1,iX1), V(iNodeX1,iX1), P(iNodeX1,iX1) )

          D0 = D(iNodeX1,iX1)
          V0 = V(iNodeX1,iX1)
          P0 = P(iNodeX1,iX1)

        END DO

      END DO

      ! --- Apply Jump Conditions ---

      D_1 = D(iNodeX1_1,iX1_1)
      V_1 = V(iNodeX1_1,iX1_1)
      P_1 = P(iNodeX1_1,iX1_1)

      W = LorentzFactor( Psi(iNodeX1_1,iX1), V_1 )

      C1 = D_1 * W * V_1
      C2 = D_1 &
             * ( SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P_1 / D_1  ) * W**2 * V_1**2 / SpeedOfLight**2 &
             + Psi(iNodeX1_1,iX1)**( -4 ) * P_1
      C3 = D_1 &
             * ( SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P_1 / D_1  ) * W**2 * V_1

      CALL ApplyJumpConditions_SAS &
             ( Psi(iNodeX1_1,iX1), V_1, C1, C2, C3, D_2, V_2, P_2 )

      Kb = P_2 / D_2**( Gamma_IDEAL )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(6x,A)') 'Shock location:'
        WRITE(*,'(8x,A)') 'Pre-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_1
        WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_1
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_1 / Kilometer, ' km'
        WRITE(*,'(8x,A)') 'Post-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_2
        WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_2
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_2 / Kilometer, ' km'
        WRITE(*,*)
        WRITE(*,'(6x,A,ES13.6E3)') &
          'Compression Ratio LOG10(D_2/D_1) = ', LOG( D_2 / D_1 ) / LOG( 1.0d1 )
        WRITE(*,*)

      END IF

      ! --- Post-shock Fields ---

      D0 = D_2
      V0 = V_2
      P0 = P_2

      DO iX1 = iX1_2, iX_B1(1), -1

        DO iNodeX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .GT. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Kb, Mdot, &
                   Alpha(iNodeX1,iX1), Psi(iNodeX1,iX1), D0, V0, P0, &
                   D(iNodeX1,iX1), V(iNodeX1,iX1), P(iNodeX1,iX1) )

          D0 = D(iNodeX1,iX1)
          V0 = V(iNodeX1,iX1)
          P0 = P(iNodeX1,iX1)

        END DO

      END DO

    END IF

    ! --- Map to 3D domain ---

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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

        iX_E(1) = iX_E0(1)
        iX_B(1) = iX_B0(1)

        IF( BX % lo(1) .EQ. 1     ) iX_B(1) = iX_B1(1)
        IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = iX_E1(1)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B (1), iX_E (1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            IF( ApplyPerturbation )THEN

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
              X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

              IF( X1 .GE. rPerturbationInner &
                    .AND. X1 .LE. rPerturbationOuter )THEN

                IF( PerturbationOrder .EQ. 0 ) &
                  uPF_K(iNodeX,iPF_D) &
                    = D(iNodeX1,iX1) &
                        * ( One + PerturbationAmplitude )

                IF( PerturbationOrder .EQ. 1 ) &
                  uPF_K(iNodeX,iPF_D) &
                    = D(iNodeX1,iX1) &
                        * ( One + PerturbationAmplitude * COS( X2 ) )

              ELSE

                uPF_K(iNodeX,iPF_D) = D(iNodeX1,iX1)

              END IF

            ELSE

              uPF_K(iNodeX,iPF_D) = D(iNodeX1,iX1)

            END IF

            uPF_K(iNodeX,iPF_V1) = V(iNodeX1,iX1)
            uPF_K(iNodeX,iPF_V2) = Zero
            uPF_K(iNodeX,iPF_V3) = Zero
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

    END DO

    DEALLOCATE( Psi0   )
    DEALLOCATE( Alpha0 )
    DEALLOCATE( Psi    )
    DEALLOCATE( Alpha  )
    DEALLOCATE( P )
    DEALLOCATE( V )
    DEALLOCATE( D )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_StandingAccretionShock_Relativistic


  ! --- Auxiliary functions/subroutines for SAS problem ---


  SUBROUTINE ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )

    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3)
    REAL(AR), INTENT(out) :: D(1:,iX_B1(1):), V(1:,iX_B1(1):), P(1:,iX_B1(1):)

    CHARACTER(LEN=16) :: FMT
    INTEGER           :: iX1

    OPEN( UNIT = 101, FILE = TRIM( NodalDataFileNameBase ) // '_D.dat' )
    OPEN( UNIT = 102, FILE = TRIM( NodalDataFileNameBase ) // '_V.dat' )
    OPEN( UNIT = 103, FILE = TRIM( NodalDataFileNameBase ) // '_P.dat' )

    READ(101,*) FMT
    READ(102,*) FMT
    READ(103,*) FMT

    DO iX1 = iX_B1(1), iX_E1(1)

      READ(101,TRIM(FMT)) D(:,iX1)
      READ(102,TRIM(FMT)) V(:,iX1)
      READ(103,TRIM(FMT)) P(:,iX1)

    END DO

    CLOSE( 103 )
    CLOSE( 102 )
    CLOSE( 101 )

  END SUBROUTINE ReadFluidFieldsFromFile


  SUBROUTINE NewtonRaphson_SAS &
    ( X1, MassPNS, K, Mdot, Alpha, Psi, D0, V0, P0, D, V, P )

    REAL(AR), INTENT(in)  :: X1, MassPNS, K, &
                             Mdot, Alpha, Psi, D0, V0, P0
    REAL(AR), INTENT(out) :: D ,V ,P

    REAL(AR) :: W
    REAL(AR) :: Jac(3,3), invJac(3,3)
    REAL(AR) :: f(3), uO(3), uN(3), du(3)

    LOGICAL             :: CONVERGED
    INTEGER             :: ITER
    REAL(AR), PARAMETER :: Tolu = 1.0e-16_AR
    REAL(AR), PARAMETER :: Tolf = 1.0e-16_AR
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


  SUBROUTINE ApplyJumpConditions_SAS( Psi, V_1, C1, C2, C3, D_2, V_2, P_2 )

    REAL(AR), INTENT(in)  :: Psi, V_1, C1, C2, C3
    REAL(AR), INTENT(out) :: D_2, V_2, P_2

    REAL(AR) :: A, B, C, D, E
    REAL(AR) :: dx, xa, xb, xc, fa, fb, fc, W

    INTEGER             :: ITER
    INTEGER,  PARAMETER :: MAX_ITER = 1000
    REAL(AR), PARAMETER :: TolChi = 1.0e-16_AR

    LOGICAL :: CONVERGED

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

    xb = 1.0e-10_AR * xa

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

    W = LorentzFactor( Psi, V_2 )

    D_2 = SpeedOfLight**( -1 ) * ABS( C1 ) * SQRT( xc**( -2 ) - Psi**4 )
    P_2 = ( C3 - D_2 * SpeedOfLight**2 * W**2 * V_2 ) &
            / ( Gamma_IDEAL / ( Gamma_IDEAL - One ) * W**2 * V_2 )

  END SUBROUTINE ApplyJumpConditions_SAS


  ! --- From: http://fortranwiki.org/fortran/show/Matrix+inversion ---
  FUNCTION Inv3x3( A ) RESULT( invA )

    ! --- Performs a direct calculation of the inverse of a 3×3 matrix ---

    REAL(AR), INTENT(in) :: A   (3,3)
    REAL(AR)             :: invA(3,3)
    REAL(AR)             :: InvDet

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


  REAL(AR) FUNCTION LorentzFactor( Psi, V )

    REAL(AR), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


  ! --- Auxiliary functions/subroutines for computine left state ---


  SUBROUTINE ComputeLeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(AR), INTENT(in)  :: Vs, DR, VR, PR
    REAL(AR), INTENT(out) ::     DL, VL, PL

    CALL ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    ! --- Return energy-density instead of pressure ---
    PL = PL / ( Gamma_IDEAL - One )

  END SUBROUTINE ComputeLeftState


  SUBROUTINE ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(AR), INTENT(in)  :: Vs, DR, VR, PR
    REAL(AR), INTENT(out) ::     DL, VL, PL

    REAL(AR), PARAMETER :: EPS = 1.0e-15_AR

    REAL(AR), PARAMETER :: ToldV = EPS
    REAL(AR), PARAMETER :: TolF  = EPS
    INTEGER,  PARAMETER :: nMaxIter = 1000

    INTEGER :: ITERATION
    REAL(AR) :: D, V, P, F
    REAL(AR) :: Vmin, Vmax, Fmin, Fmax, VV, FF

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


  REAL(AR) FUNCTION Density( Vs, DR, VR, VL )

    REAL(AR), INTENT(in) :: Vs, DR, VR, VL

    REAL(AR) :: WR, WL

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    Density = DR * ( WR * ( VR - Vs ) ) / ( WL * ( VL - Vs ) )

    RETURN
  END FUNCTION Density


  REAL(AR) FUNCTION Pressure( Vs, DR, VR, PR, DL, VL )

    REAL(AR), INTENT(in) :: Vs, DR, VR, PR, DL, VL

    REAL(AR) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    Pressure = ( PR * ( One + tau * WR**2 * VR * ( VR - Vs ) ) &
                 - DL * WL**2 * VL**2 + DR * WR**2 * VR**2 &
                 + Vs * ( DL * WL**2 * VL - DR * WR**2 * VR ) ) &
               / ( One + tau * WL**2 * VL * ( VL - Vs ) )

    RETURN
  END FUNCTION Pressure


  REAL(AR) FUNCTION PostShockVelocity( Vs, DR, VR, PR, DL, VL, PL )

    REAL(AR), INTENT(in) :: Vs, DR, VR, PR, DL, VL, PL

    REAL(AR) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    PostShockVelocity &
      = ( DL + tau * PL ) * WL**2 * ( VL - Vs ) &
          - ( DR + tau * PR ) * WR**2 * ( VR - Vs ) + Vs * ( PL - PR )

    RETURN
  END FUNCTION PostShockVelocity


END MODULE MF_InitializationModule_Relativistic_IDEAL
