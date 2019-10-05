MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,  ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP, Half, One, Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Erg, Second, Kelvin
  USE ProgramHeaderModule,     ONLY: &
    nDOFX, nX, nNodesX, swX, nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule,              ONLY: &
    MeshType,                        &
    CreateMesh, DestroyMesh,         &
    NodeCoordinate
  USE GeometryFieldsModule,    ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule,       ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iPF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_Ne, &
    nAF, iAF_P, iAF_E, iAF_T, iAF_Ye
  USE Euler_UtilitiesModule,   ONLY: &
    ComputeConserved_Euler
  USE EquationOfStateModule,   ONLY: &
    ComputePressureFromPrimitive, &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE UnitsModule, ONLY: &
    Meter, Kilometer, Kilogram, Second, Joule
  USE UtilitiesModule, ONLY: &
    Locate

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, xL, xR, Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields


CONTAINS


  SUBROUTINE MF_InitializeFields( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )
    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      ! --- Non-relativistic ---

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

      ! --- Table ---

      CASE( 'Sod_TABLE' )

        CALL InitializeFields_Sod_TABLE( MF_uGF, MF_uCF )

      ! --- Relativistic ---

      CASE( 'Sod_Relativistic' )

        CALL InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz_Relativistic' )

        CALL InitializeFields_KelvinHelmholtz_Relativistic( MF_uGF, MF_uCF )

      CASE( 'RiemannProblem_2D_Relativistic' )

        CALL InitializeFields_RiemannProblem_2D_Relativistic( MF_uGF, MF_uCF )

      CASE( 'StandingAccretionShock_Relativistic' )

        CALL InitializeFields_StandingAccretionShock_Relativistic &
               ( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A,A)') '', 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(A4,A)')   '', 'Valid Options:'
          WRITE(*,'(A6,A)')   '', 'Non-Relativistic (IDEAL):'
          WRITE(*,'(A8,A)')   '', 'IsentropicVortex'
          WRITE(*,'(A8,A)')   '', 'Sod'
          WRITE(*,'(A8,A)')   '', 'SphericalSod'
          WRITE(*,'(A8,A)')   '', 'TopHatAdvection'
          WRITE(*,'(A8,A)')   '', 'Implosion'
          WRITE(*,'(A8,A)')   '', 'StandingAccretionShock'
          WRITE(*,'(A6,A)')   '', 'Non-Relativistic (TABLE):'
          WRITE(*,'(A8,A)')   '', 'Sod_TABLE'
          WRITE(*,'(A6,A)')   '', 'Relativistic:'
          WRITE(*,'(A8,A)')   '', 'Sod_Relativistic'
          WRITE(*,'(A8,A)')   '', 'KelvinHelmholtz_Relativistic'
          WRITE(*,'(A8,A)')   '', 'RiemannProblem_2D_Relativistic'
          WRITE(*,'(A8,A)')   '', 'StandingAccretionShock_Relativistic'
          STOP 'MF_InitializationModule.f90'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields


  SUBROUTINE InitializeFields_IsentropicVortex( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    REAL(amrex_real), PARAMETER :: Beta  = 5.0_amrex_real

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: iNodeX1, iNodeX2
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1, X2, R
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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
              = 0.0_DP

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( X1 .LE. 0.5_amrex_real ) THEN

              uPF_K(iNodeX,iPF_D)  = 1.0_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 1.0_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 0.1_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( X1 <= 1.0_DP ) THEN

              uPF_K(iNodeX,iPF_D)  = 1.0_DP
              uPF_K(iNodeX,iPF_V1) = 0.0_DP
              uPF_K(iNodeX,iPF_V2) = 0.0_DP
              uPF_K(iNodeX,iPF_V3) = 0.0_DP
              uPF_K(iNodeX,iPF_E)  = 1.0_DP / ( Gamma_IDEAL - 1.0_DP )

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_DP
              uPF_K(iNodeX,iPF_V1) = 0.0_DP
              uPF_K(iNodeX,iPF_V2) = 0.0_DP
              uPF_K(iNodeX,iPF_V3) = 0.0_DP
              uPF_K(iNodeX,iPF_E)  = 0.1_DP / ( Gamma_IDEAL - 1.0_DP )

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: iNodeX1
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( ABS( X1 - 0.5_DP ) .LE. 0.25_DP )THEN
              uPF_K(iNodeX,iPF_D) = 2.0_DP
            ELSE
              uPF_K(iNodeX,iPF_D) = 1.0_DP
            END IF

          END DO

          uPF_K(:,iPF_V1) = 1.0_DP
          uPF_K(:,iPF_V2) = 0.0_DP
          uPF_K(:,iPF_V3) = 0.0_DP
          uPF_K(:,iPF_E)  = 1.0d-2

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    REAL(amrex_real), PARAMETER :: Beta  = 5.0_amrex_real

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: iNodeX1, iNodeX2
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1, X2
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(amrex_real), PARAMETER :: D_0 = 0.125_amrex_real
    REAL(amrex_real), PARAMETER :: E_0 = 0.350_amrex_real
    REAL(amrex_real), PARAMETER :: D_1 = 1.000_amrex_real
    REAL(amrex_real), PARAMETER :: E_1 = 2.500_amrex_real

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( X1 + X2 .LT. 0.15_amrex_real )THEN

              uPF_K(iNodeX,iPF_D) &
                = D_0
              uPF_K(iNodeX,iPF_E) &
                = E_0

            ELSE

              uPF_K(iNodeX,iPF_D) &
                = D_1
              uPF_K(iNodeX,iPF_E) &
                = E_1

            END IF

            uPF_K(iNodeX,iPF_V1) &
              = 0.0_amrex_real
            uPF_K(iNodeX,iPF_V2) &
              = 0.0_amrex_real
            uPF_K(iNodeX,iPF_V3) &
              = 0.0_amrex_real

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    ! --- AMReX variables ---
    INTEGER               :: iDim, iLevel
    INTEGER               :: lo_G(4), hi_G(4)
    INTEGER               :: lo_F(4), hi_F(4)
    REAL(amrex_real)      :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)      :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)      :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)      :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(MeshType)        :: MeshX(3)
    TYPE(amrex_parmparse) :: PP
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- thornado variables ---
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real)   :: X1, X2, Alpha, Speed, D_Prime, V1_Prime, P_Prime

    REAL(amrex_real) :: mDot, Mass, rShock, Mach
    INTEGER          :: PerturbOrder
    REAL(amrex_real) :: ShellIn, ShellOut, PerturbAmplitude
    LOGICAL          :: Perturb

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

    Alpha = 4.0_amrex_real * Gamma_IDEAL &
              / ( ( Gamma_IDEAL + 1.0_amrex_real ) &
                  * ( Gamma_IDEAL - 1.0_amrex_real ) ) &
              * ( ( Gamma_IDEAL - 1.0_amrex_real ) &
                  / ( Gamma_IDEAL + 1.0_amrex_real ) )**Gamma_IDEAL

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

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
                     ( X1, Alpha, Gamma_IDEAL, Mass, Speed )

              uPF_K(iNodeX,iPF_D)  &
                = ( mDot / FourPi ) &
                    * Speed**(-1) * X1**(-2)

              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) &
                = 0.0_amrex_real

              V1_prime &
                = ( Gamma_IDEAL - 1.0_amrex_real ) &
                    / ( Gamma_IDEAL + 1.0_amrex_real ) &
                    * SQRT( 2.0_amrex_real * Mass / rShock )

              D_prime &
                = ( mDot / FourPi ) * ( 1.0_amrex_real / V1_prime ) &
                    * rShock**(-2)

              P_prime &
                = 2.0_amrex_real / ( Gamma_IDEAL + 1.0_amrex_real ) &
                    * ( mDot / FourPi ) &
                    * SQRT( 2.0_amrex_real * Mass ) * rShock**(-2.5_amrex_real)

              uPF_K(iNodeX,iPF_E)  &
                = P_Prime &
                    * ( uPF_K(iNodeX,iPF_D) / D_prime )**Gamma_IDEAL &
                      / ( Gamma_IDEAL - 1.0_amrex_real )

            ELSE

              Speed &
                = SQRT( 1.0_amrex_real / ( 1.0_amrex_real &
                    + 2.0_amrex_real &
                    / ( ( Gamma_IDEAL - 1.0_amrex_real) * Mach**2 ) ) / X1 )

              uPF_K(iNodeX,iPF_D)  &
                = mDot / ( FourPi * X1**2 * Speed )
              uPF_K(iNodeX,iPF_V1) &
                = - Speed
              uPF_K(iNodeX,iPF_V2) &
                = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) &
                = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  &
                = uPF_K(iNodeX,iPF_D) / Gamma_IDEAL &
                    * ( Speed / Mach )**2 / ( Gamma_IDEAL - 1.0_amrex_real )

            END IF

            IF( Perturb .AND. ( X1 .GE. ShellIn ) &
                        .AND. ( X1 .LE. ShellOut ) ) THEN

              SELECT CASE( PerturbOrder )

              CASE( 0 )

                uPF_K(iNodeX, iPF_D) &
                  = ( 1.0_amrex_real + PerturbAmplitude ) &
                      * uPF_K(iNodeX,iPF_D)

              CASE( 1 )

                uPF_K(iNodeX, iPF_D) &
                  = ( 1.0_amrex_real + PerturbAmplitude * COS(X2) ) &
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


  SUBROUTINE ComputeSettlingSpeed_Bisection( r, Alpha, Gamma, Mass, V1 )

    REAL(amrex_real), INTENT(in) :: r, Alpha, Gamma, Mass

    LOGICAL :: Converged
    INTEGER :: Iter
    REAL(amrex_real) :: a, b, c, ab, F_a, F_b, F_c, F_0
    INTEGER, PARAMETER :: MaxIter = 128
    REAL(amrex_real), PARAMETER :: Tol_ab = 1.0d-8
    REAL(amrex_real), PARAMETER :: Tol_F = 1.0d-8

    REAL(amrex_real), INTENT(out) :: V1

    a = 1.0d-6
    F_a = SettlingSpeedFun(a, r, Alpha, Gamma, Mass)

    b = 1.0_amrex_real
    F_b = SettlingSpeedFun(b, r, Alpha, Gamma, Mass)

    F_0 = F_a
    ab = b - a

    Converged = .FALSE.
    Iter = 0

    DO WHILE ( .NOT. Converged)

      Iter = Iter + 1

      ab = 0.5_amrex_real * ab
      c = a + ab

      F_c = SettlingSpeedFun(c, r, Alpha, Gamma, Mass)

      IF( F_a * F_c < 0.0_amrex_real ) THEN

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


  REAL(amrex_real) FUNCTION SettlingSpeedFun( u, r, Alpha, Gamma, Mass )

    REAL(amrex_real), INTENT(in) :: u, r, Alpha, Gamma, Mass

    SettlingSpeedFun &
      = r * u**2 &
        + Alpha * r**(3.0_amrex_real-2.0_amrex_real*Gamma) &
        * u**(1.0_amrex_real-Gamma) &
        - 2.0_amrex_real * Mass
    RETURN
  END FUNCTION SettlingSpeedFun


  SUBROUTINE InitializeFields_Sod_TABLE( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real
    uAF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( X1 .LE. 0.0_amrex_real ) THEN

              uPF_K(iNodeX,iPF_D ) = 1.00e12_amrex_real * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real     * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real     * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real     * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e32_amrex_real  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.4_amrex_real

            ELSE

              uPF_K(iNodeX,iPF_D ) = 1.25e11_amrex_real * Gram / Centimeter**3
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real     * Kilometer / Second
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real     * Kilometer / Second
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real     * Kilometer / Second
              uAF_K(iNodeX,iAF_P ) = 1.0e31_amrex_real  * Erg / Centimeter**3
              uAF_K(iNodeX,iAF_Ye) = 0.3_amrex_real

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


  SUBROUTINE InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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

            IF( X1 .LE. 0.5_amrex_real ) THEN

              uPF_K(iNodeX,iPF_D)  = 1.0_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 1.0_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 0.1_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER          :: iX1, iX2, iX3
    INTEGER          :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real) :: X1, X2
    REAL(amrex_real) :: a, Vshear, A0, sigma, rho0, rho1

    INTEGER            :: iDim, iLevel
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    ! --- Problem-dependent parameters ---
    a      = 0.01_amrex_real
    Vshear = 0.5_amrex_real

    A0    = 0.1_amrex_real
    sigma = 0.1_amrex_real

    rho0 = 0.505_amrex_real
    rho1 = 0.495_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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
            IF( X2 .GT. 0.0_amrex_real )THEN
              uPF_K(iNodeX,iPF_V1) &
                = +Vshear * TANH( ( X2 - 0.5_amrex_real ) / a )
            ELSE
              ! --- Paper has a typo here, the minus sign is required ---
              uPF_K(iNodeX,iPF_V1) &
                = -Vshear * TANH( ( X2 + 0.5_amrex_real ) / a )
            END IF

            ! --- V2 ---
            IF( X2 .GT. 0.0_amrex_real )THEN
              uPF_K(iNodeX,iPF_V2) &
                =  A0 * Vshear * SIN( 2.0_amrex_real * Pi * X1 ) &
                    * EXP( -( ( X2 - 0.5_amrex_real )**2 / sigma ) )
            ELSE
              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( 2.0_amrex_real * Pi * X1 ) &
                    * EXP( -( ( X2 + 0.5_amrex_real )**2 / sigma ) )
            END IF

            ! --- rho ---
            IF( X2 .GT. 0.0_amrex_real )THEN
              uPF_K(iNodeX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - 0.5_amrex_real ) / a )
            ELSE
              uPF_K(iNodeX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + 0.5_amrex_real ) / a )
            END IF

            uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
            uPF_K(iNodeX,iPF_E)  = 1.0_amrex_real / ( Gamma_IDEAL - One )

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


  ! --- Relativistic 2D Riemann problem from
  !     Del Zanna & Bucciantini, (2002), A&A, 390, 1177 ---
  SUBROUTINE InitializeFields_RiemannProblem_2D_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER          :: iX1, iX2, iX3
    INTEGER          :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real) :: X1, X2

    INTEGER            :: iDim, iLevel
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

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
            IF     ( X1 .GT. 0.5_amrex_real .AND. X2 .GT. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 0.01_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            ! --- NW ---
            ELSE IF( X1 .LE. 0.5_amrex_real .AND. X2 .GT. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.99_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            ! --- SW ---
            ELSE IF( X1 .LE. 0.5_amrex_real .AND. X2 .LE. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.5_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )

            ! --- SE ---
            ELSE
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.99_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            END IF

          END DO

          uPF_K(:,iPF_V3) = 0.0_amrex_real

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

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER          :: iX1, iX2, iX3
    INTEGER          :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real) :: X1, X2

    INTEGER            :: iDim, iLevel
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-specific parameters ---
    INTEGER, PARAMETER    :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER               :: iL, nLines
    REAL(DP), ALLOCATABLE :: FluidFieldData(:,:), FluidFieldParameters(:)
    TYPE(amrex_parmparse) :: PP
    LOGICAL               :: ApplyPerturbation
    INTEGER               :: PerturbationOrder
    REAL(amrex_real)      :: rPerturbationInner, rPerturbationOuter, &
                             PerturbationAmplitude

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'ApplyPerturbation',     ApplyPerturbation )
      CALL PP % get( 'PerturbationOrder',     PerturbationOrder )
      CALL PP % get( 'rPerturbationInner',    rPerturbationInner )
      CALL PP % get( 'rPerturbationOuter',    rPerturbationOuter )
      CALL PP % get( 'PerturbationAmplitude', PerturbationAmplitude )
    CALL amrex_parmparse_destroy( PP )
    rPerturbationInner = rPerturbationInner * Kilometer
    rPerturbationOuter = rPerturbationOuter * Kilometer

    uGF_K = 0.0_amrex_real
    uPF_K = 0.0_amrex_real
    uCF_K = 0.0_amrex_real

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(1), &
               xL(iDim), xR(iDim) )

    END DO

    CALL ReadParameters &
           ( '../Euler_Relativistic/StandingAccretionShock_Parameters.dat', &
             FluidFieldParameters )
    CALL ReadData &
           ( '../Euler_Relativistic/StandingAccretionShock_Data.dat', &
             nLines, FluidFieldData )

    DO iLevel = 0, nLevels

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

            uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real

            uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real

            uPF_K(iNodeX,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            uPF_K(iNodeX,iPF_Ne) = 0.0_amrex_real

            ! --- Apply perturbations ---
            IF( ApplyPerturbation )THEN

              IF     ( PerturbationOrder .EQ. 0 )THEN

                IF( X1 .GE. rPerturbationInner &
                      .AND. X1 .LE. rPerturbationOuter )THEN

                  uPF_K(iNodeX,iPF_D) &
                    = uPF_K(iNodeX,iPF_D) &
                        * ( 1.0e0_amrex_real &
                              + PerturbationAmplitude )

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
                        * ( 1.0e0_amrex_real &
                              + PerturbationAmplitude * COS( X2 ) )
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

  REAL(amrex_real) FUNCTION InterpolateInitialConditionsOntoGrid &
    ( iVar, i_r, iL, X, FluidFieldData )  RESULT( yInterp )

    INTEGER,  INTENT(in)         :: iVar, i_r, iL
    REAL(amrex_real), INTENT(in) :: X
    REAL(amrex_real), INTENT(in) :: FluidFieldData(:,:)
    REAL(amrex_real)             :: X1, X2, Y1, Y2, m

    X1 = FluidFieldData(iL,i_r)
    X2 = FLuidFieldData(iL+1,i_r)
    Y1 = FluidFieldData(iL,iVar)
    Y2 = FluidFieldData(iL+1,iVar)

    m = ( Y2 - Y1 ) / ( X2 - X1 )

    yInterp = m * ( X - X1 ) + Y1

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid


  SUBROUTINE ReadParameters( FILEIN, FluidFieldParameters )

    CHARACTER(LEN=*), INTENT(in)               :: FILEIN
    REAL(amrex_real), INTENT(out), ALLOCATABLE :: FluidFieldParameters(:)

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
    REAL(amrex_real), INTENT(out), ALLOCATABLE :: FluidFieldData(:,:)

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


END MODULE MF_InitializationModule
