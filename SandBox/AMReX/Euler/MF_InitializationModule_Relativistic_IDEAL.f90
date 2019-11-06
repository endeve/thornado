MODULE MF_InitializationModule_Relativistic_IDEAL

  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    AR => amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,  ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
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
  USE UnitsModule, ONLY: &
    Meter, &
    Kilometer, &
    Kilogram, &
    Second, &
    Joule
  USE UtilitiesModule, ONLY: &
    Locate

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_Relativistic_IDEAL

  REAL(AR), PARAMETER :: Zero   = 0.0_AR
  REAL(AR), PARAMETER :: Half   = 0.5_AR
  REAL(AR), PARAMETER :: One    = 1.0_AR
  REAL(AR), PARAMETER :: Two    = 2.0_AR
  REAL(AR), PARAMETER :: Three  = 3.0_AR
  REAL(AR), PARAMETER :: Pi     = ACOS( -1.0_AR )
  REAL(AR), PARAMETER :: TwoPi  = 2.0_AR * Pi
  REAL(AR), PARAMETER :: FourPi = 4.0_AR * Pi


CONTAINS


  SUBROUTINE MF_InitializeFields_Relativistic_IDEAL &
    ( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )
    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Sod_Relativistic' )

        CALL InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz_Relativistic' )

        CALL InitializeFields_KelvinHelmholtz_Relativistic( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz_Relativistic_3D' )

        CALL InitializeFields_KelvinHelmholtz_Relativistic_3D( MF_uGF, MF_uCF )

      CASE( 'RiemannProblem_2D_Relativistic' )

        CALL InitializeFields_RiemannProblem_2D_Relativistic( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock_Relativistic' )

        CALL InitializeFields_StandingAccretionShock_Relativistic &
               ( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'Sod_Relativistic'
          WRITE(*,'(6x,A)')     'KelvinHelmholtz_Relativistic'
          WRITE(*,'(6x,A)')     'KelvinHelmholtz_Relativistic_3D'
          WRITE(*,'(6x,A)')     'RiemannProblem_2D_Relativistic'
          WRITE(*,'(6x,A)')     'StandingAccretionShock_Relativistic'
          STOP 'MF_InitializationModule.f90'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields_Relativistic_IDEAL


  SUBROUTINE InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2, X3
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

          X1 = MeshX(1) % Center(iX1)

          DO iNodeX = 1, nDOFX

            IF( X1 .LE. Half ) THEN

              uPF_K(iNodeX,iPF_D)  = One
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = One / (Gamma_IDEAL - One)

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_AR
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = 0.1_AR / (Gamma_IDEAL - One)

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

  END SUBROUTINE InitializeFields_Sod_Relativistic


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Radice & Rezzolla, (2012), AA, 547, A26 ---
  SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2, X3
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
                =  A0 * Vshear * SIN( 2.0_AR * Pi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma ) )
            ELSE
              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( 2.0_AR * Pi * X1 ) &
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

  END SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic


  ! --- Relativistic 3D Kelvin-Helmholtz instability a la
  !     Radice & Rezzolla, (2012), AA, 547, A26 ---
  SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic_3D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2, X3
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
                =  A0 * Vshear * SIN( 2.0_AR * Pi * X1 ) &
                    * EXP( -( ( X2 - Half )**2 / sigma ) )
            ELSE
              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( 2.0_AR * Pi * X1 ) &
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

  END SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic_3D


  ! --- Relativistic 2D Riemann problem from
  !     Del Zanna & Bucciantini, (2002), A&A, 390, 1177 ---
  SUBROUTINE InitializeFields_RiemannProblem_2D_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2, X3
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

            ! --- NE ---
            IF     ( X1 .GT. Half .AND. X2 .GT. Half )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_AR
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_E ) = 0.01_AR / ( Gamma_IDEAL - One )
            ! --- NW ---
            ELSE IF( X1 .LE. Half .AND. X2 .GT. Half )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_AR
              uPF_K(iNodeX,iPF_V1) = 0.99_AR
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )
            ! --- SW ---
            ELSE IF( X1 .LE. Half .AND. X2 .LE. Half )THEN
              uPF_K(iNodeX,iPF_D ) = Half
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )

            ! --- SE ---
            ELSE
              uPF_K(iNodeX,iPF_D ) = 0.1_AR
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = 0.99_AR
              uPF_K(iNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )
            END IF

          END DO

          uPF_K(:,iPF_V3) = Zero

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

  END SUBROUTINE InitializeFields_RiemannProblem_2D_Relativistic


  SUBROUTINE InitializeFields_StandingAccretionShock_Relativistic &
    ( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(AR)       :: X1, X2, X3
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
    INTEGER, PARAMETER    :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER               :: iL, nLines
    REAL(AR), ALLOCATABLE :: FluidFieldData(:,:), FluidFieldParameters(:)
    LOGICAL               :: ApplyPerturbation
    INTEGER               :: PerturbationOrder
    REAL(AR)              :: rPerturbationInner, rPerturbationOuter, &
                               PerturbationAmplitude

    TYPE(amrex_parmparse) :: PP

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'ApplyPerturbation',     ApplyPerturbation )
      CALL PP % get( 'PerturbationOrder',     PerturbationOrder )
      CALL PP % get( 'rPerturbationInner',    rPerturbationInner )
      CALL PP % get( 'rPerturbationOuter',    rPerturbationOuter )
      CALL PP % get( 'PerturbationAmplitude', PerturbationAmplitude )
    CALL amrex_parmparse_destroy( PP )
    rPerturbationInner = rPerturbationInner * Kilometer
    rPerturbationOuter = rPerturbationOuter * Kilometer

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
               xL(iDim), xR(iDim) )

    END DO

    CALL ReadParameters &
      ( '../Euler_Relativistic_IDEAL/StandingAccretionShock_Parameters.dat', &
             FluidFieldParameters )
    CALL ReadData &
      ( '../Euler_Relativistic_IDEAL/StandingAccretionShock_Data.dat', &
        nLines, FluidFieldData )

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
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            ! --- Get lower index of input array
            !     (FluidFieldData) corresponding to physical coordinate (X1) ---
            iL = Locate( X1, FluidFieldData(:,i_r), nLines )

            ! --- Interpolate to the physical point X1 ---

            uPF_K(iNodeX,iPF_D) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_D, i_r, iL, X1, FluidFieldData )

            uPF_K(iNodeX,iPF_V1) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_V1, i_r, iL, X1, FluidFieldData )

            uPF_K(iNodeX,iPF_V2) = Zero

            uPF_K(iNodeX,iPF_V3) = Zero

            uPF_K(iNodeX,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            uPF_K(iNodeX,iPF_Ne) = Zero

            ! --- Apply perturbations ---
            IF( ApplyPerturbation )THEN

              IF     ( PerturbationOrder .EQ. 0 )THEN

                IF( X1 .GE. rPerturbationInner &
                      .AND. X1 .LE. rPerturbationOuter )THEN

                  uPF_K(iNodeX,iPF_D) &
                    = uPF_K(iNodeX,iPF_D) &
                        * ( One + PerturbationAmplitude )

                END IF

              ELSE IF( PerturbationOrder .EQ. 1 )THEN

                IF( nDimsX .EQ. 1 )THEN
                  WRITE(*,'(A)') 'Cannot have l = 1 mode perturbation in 1D.'
                  WRITE(*,'(A)') 'Stopping...'
                  STOP
                END IF

                IF( X1 .GE. rPerturbationInner &
                      .AND. X1 .LE. rPerturbationOuter )THEN

                  uPF_K(iNodeX,iPF_D) &
                    = uPF_K(iNodeX,iPF_D) &
                        * ( One + PerturbationAmplitude * COS( X2 ) )
                END IF

              ELSE

                WRITE(*,'(A)') 'Fatal Error'
                WRITE(*,'(A)') '-----------'
                WRITE(*,'(A,I1)') 'Invalid value of PerturbationOrder: ', &
                             PerturbationOrder
                WRITE(*,'(A)') 'Valid values: 0, 1'
                WRITE(*,'(A)') 'Stopping...'
                STOP

              END IF

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

  END SUBROUTINE InitializeFields_StandingAccretionShock_Relativistic


  ! --- Auxiliary functions/subroutines for relativistic SAS problem ---

  REAL(AR) FUNCTION InterpolateInitialConditionsOntoGrid &
    ( iVar, i_r, iL, X, FluidFieldData )  RESULT( yInterp )

    INTEGER,  INTENT(in) :: iVar, i_r, iL
    REAL(AR), INTENT(in) :: X
    REAL(AR), INTENT(in) :: FluidFieldData(:,:)

    REAL(AR)             :: Xa, Xb, Ya, Yb, m

    Xa = FluidFieldData(iL,i_r)
    Xb = FLuidFieldData(iL+1,i_r)
    Ya = FluidFieldData(iL,iVar)
    Yb = FluidFieldData(iL+1,iVar)

    m = ( Yb - Ya ) / ( Xb - Xa )

    yInterp = m * ( X - Xa ) + Ya

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid


  SUBROUTINE ReadParameters( FILEIN, FluidFieldParameters )

    CHARACTER(LEN=*), INTENT(in)               :: FILEIN
    REAL(AR),         INTENT(out), ALLOCATABLE :: FluidFieldParameters(:)

    INTEGER :: i, nParams

    ! --- Get number of parameters ---
    nParams = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nParams = nParams + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in parameters ---
    ALLOCATE( FluidFieldParameters(nParams) )

    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nParams
       READ( 100, '(ES23.16E2)' ) FluidFieldParameters(i)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldParameters(1) = FluidFieldParameters(1) * Kilogram
    FluidFieldParameters(2) = FluidFieldParameters(2)
    FluidFieldParameters(3) = FluidFieldParameters(3) * Meter
    FluidFieldParameters(4) = FluidFieldParameters(4) * Meter
    FluidFieldParameters(5) = FluidFieldParameters(5) * Meter
    FluidFieldParameters(6) = FluidFieldParameters(6) * Meter
    FluidFieldParameters(7) = FluidFieldParameters(7) * Kilogram / Second

  END SUBROUTINE ReadParameters


  SUBROUTINE ReadData( FILEIN, nLines, FluidFieldData )

    CHARACTER(LEN=*), INTENT(in)               :: FILEIN
    INTEGER,          INTENT(inout)            :: nLines
    REAL(AR),         INTENT(out), ALLOCATABLE :: FluidFieldData(:,:)

    INTEGER :: i

    ! --- Get number of lines in data file ---
    nLines = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nLines = nLines + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in data ---
    ALLOCATE( FluidFieldData( 1:nLines, 4 ) )

    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nLines
       READ( 100, '(ES22.16E2,1x,ES22.16E2,1x,ES23.16E2,1x,ES22.16E2)' ) &
         FluidFieldData(i,:)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldData(:,1) = FluidFieldData(:,1) * Meter
    FluidFieldData(:,2) = FluidFieldData(:,2) * Kilogram / Meter**3
    FluidFieldData(:,3) = FluidFieldData(:,3) * Meter / Second
    FluidFieldData(:,4) = FluidFieldData(:,4) * Joule / Meter**3

  END SUBROUTINE ReadData


END MODULE MF_InitializationModule_Relativistic_IDEAL
