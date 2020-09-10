MODULE MF_InitializationModule_NonRelativistic_IDEAL

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
    NodeNumberTableX
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
    Kilometer
  USE Euler_ErrorModule,       ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE InputParsingModule,      ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields_NonRelativistic_IDEAL

  REAL(AR), PARAMETER :: Zero   = 0.0_AR
  REAL(AR), PARAMETER :: Half   = 0.5_AR
  REAL(AR), PARAMETER :: One    = 1.0_AR
  REAL(AR), PARAMETER :: Two    = 2.0_AR
  REAL(AR), PARAMETER :: Three  = 3.0_AR
  REAL(AR), PARAMETER :: Four   = 4.0_AR
  REAL(AR), PARAMETER :: Five   = 5.0_AR
  REAL(AR), PARAMETER :: Pi     = ACOS( -1.0_AR )
  REAL(AR), PARAMETER :: TwoPi  = 2.0_AR * Pi
  REAL(AR), PARAMETER :: FourPi = 4.0_AR * Pi


CONTAINS


  SUBROUTINE MF_InitializeFields_NonRelativistic_IDEAL &
    ( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

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

      CASE( 'SASI' )

        CALL InitializeFields_SASI( MF_uGF, MF_uCF )

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
          WRITE(*,'(6x,A)')     'SASI'

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
    REAL(AR), PARAMETER :: Beta = 5.0_AR
    REAL(AR)            :: R

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

            R = SQRT( X1**2 + X2**2 )

            uPF_K(iNodeX,iPF_D) &
              = ( One - ( Gamma_IDEAL - One ) * Beta**2 &
                        / ( 8.0_AR * Gamma_IDEAL * Pi**2 ) * EXP( One - R**2 ) &
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
    REAL(AR)       :: X1
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

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

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

  END SUBROUTINE InitializeFields_Sod


  SUBROUTINE InitializeFields_SphericalSod( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3
    INTEGER        :: iNodeX, iNodeX1
    REAL(AR)       :: X1
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

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 <= One ) THEN

              uPF_K(iNodeX,iPF_D)  = One
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = One / ( Gamma_IDEAL - One )

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_AR
              uPF_K(iNodeX,iPF_V1) = Zero
              uPF_K(iNodeX,iPF_V2) = Zero
              uPF_K(iNodeX,iPF_V3) = Zero
              uPF_K(iNodeX,iPF_E)  = 0.1_AR / ( Gamma_IDEAL - One )

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
    REAL(AR)       :: X1
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

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( ABS( X1 - Half ) .LE. 0.25_AR )THEN

              uPF_K(iNodeX,iPF_D) = Two

            ELSE

              uPF_K(iNodeX,iPF_D) = One

            END IF

          END DO

          uPF_K(:,iPF_V1) = One
          uPF_K(:,iPF_V2) = Zero
          uPF_K(:,iPF_V3) = Zero
          uPF_K(:,iPF_E)  = 0.01_AR

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
    REAL(AR), PARAMETER :: Beta = 5.0_AR
    REAL(AR), PARAMETER :: D_0  = 0.125_AR
    REAL(AR), PARAMETER :: E_0  = 0.350_AR
    REAL(AR), PARAMETER :: D_1  = 1.000_AR
    REAL(AR), PARAMETER :: E_1  = 2.500_AR

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

            IF( X1 + X2 .LT. 0.15_AR )THEN

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
    REAL(AR) :: Alpha, Speed, D_Prime, V1_Prime, P_Prime
    REAL(AR) :: mDot, Mass, rShock, Mach

    LOGICAL  :: Perturb
    INTEGER  :: PerturbOrder
    REAL(AR) :: ShellIn, ShellOut, PerturbAmplitude

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

    Alpha = 4.0_AR * Gamma_IDEAL &
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
                              * Mass ) * rShock**(-2.5_AR)

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
    TYPE(amrex_parmparse) :: PP

    ! --- Problem-dependent Parameters ---
    REAL(AR) :: Alpha, Speed, D_Prime, V1_Prime, P_Prime
    REAL(AR) :: mDot, Mass, rShock, Mach

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
                            * Mass ) * rShock**(-2.5_AR)

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
    TYPE(amrex_parmparse) :: PP

    ! --- Problem-dependent Parameters ---
    REAL(AR) :: mDot, Mass, rShock, Entropy

    INTEGER  :: iX_B1(3), iX_E1(3)
    REAL(AR) :: BernoulliConstant

    REAL(AR), ALLOCATABLE :: R(:,:), D(:,:), V(:,:), P(:,:)

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

    BernoulliConstant = 0.0_AR

    CALL ComputeInitialConditions_NewtonRaphson &
           ( R, mDot, Mass, rShock, Entropy, &
             BernoulliConstant, D, V, P )

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


  SUBROUTINE InitializeFields_SASI &
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
    INTEGER               :: iX1_1, iX1_2, iNodeX1_1, iNodeX1_2
    REAL(AR)              :: X1_1, X1_2, dX1
    REAL(AR)              :: alpha, PolytropicConstant, V0
    LOGICAL               :: FirstPreShockElement = .FALSE.
    INTEGER               :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(AR), ALLOCATABLE :: D(:,:), V(:,:), P(:,:)
    REAL(AR)              :: MassPNS, ShockRadius, AccretionRate, MachNumber
    LOGICAL               :: ApplyPerturbation
    INTEGER               :: PerturbationOrder
    REAL(AR)              :: PerturbationAmplitude, &
                             rPerturbationInner, rPerturbationOuter

    ApplyPerturbation     = .FALSE.
    PerturbationOrder     = 0
    PerturbationAmplitude = Zero
    rPerturbationInner    = Zero
    rPerturbationOuter    = Zero
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'Mass'                 , MassPNS               )
      CALL PP % get  ( 'AccretionRate'        , AccretionRate         )
      CALL PP % get  ( 'ShockRadius'          , ShockRadius           )
      CALL PP % get  ( 'MachNumber'           , MachNumber            )
      CALL PP % query( 'ApplyPerturbation'    , ApplyPerturbation     )
      CALL PP % query( 'PerturbationOrder'    , PerturbationOrder     )
      CALL PP % query( 'PerturbationAmplitude', PerturbationAmplitude )
      CALL PP % query( 'rPerturbationInner'   , rPerturbationInner    )
      CALL PP % query( 'rPerturbationOuter'   , rPerturbationOuter    )
    CALL amrex_parmparse_destroy( PP )

    MassPNS            = MassPNS            * SolarMass
    AccretionRate      = AccretionRate      * SolarMass / Second
    ShockRadius        = ShockRadius        * Kilometer
    rPerturbationInner = rPerturbationInner * Kilometer
    rPerturbationOuter = rPerturbationOuter * Kilometer

    alpha = Four * Gamma_IDEAL &
              / ( ( Gamma_IDEAL - One ) * ( Gamma_IDEAL + One ) ) &
              * ( ( Gamma_IDEAL - One ) / ( Gamma_IDEAL + One ) &
                )**Gamma_IDEAL

    PolytropicConstant &
      = Two / ( Gamma_IDEAL + One ) &
          * ( ( Gamma_IDEAL - One ) / ( Gamma_IDEAL + One ) )**( Gamma_IDEAL ) &
          * ( AccretionRate / FourPi )**( One - Gamma_IDEAL ) &
          * ( Two * GravitationalConstant * MassPNS &
            )**( ( Gamma_IDEAL + One ) / Two ) &
          * ShockRadius**( ( Three * Gamma_IDEAL - Five ) / Two )

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Shock radius:   ', ShockRadius / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'PNS Mass:       ', MassPNS / SolarMass, ' Msun'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Accretion Rate: ', AccretionRate / ( SolarMass / Second ), &
        ' Msun/s'

      WRITE(*,'(6x,A,ES9.2E3)') &
        'Mach number:    ', MachNumber

      WRITE(*,'(6x,A,ES9.2E3)') &
        'alpha:          ', alpha

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'Apply Perturbation: ', ApplyPerturbation

      WRITE(*,'(6x,A,I1)') &
        'Perturbation order: ', PerturbationOrder

      WRITE(*,'(6x,A,ES9.2E3)') &
        'Perturbation amplitude: ', PerturbationAmplitude

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Inner radius of perturbation: ', rPerturbationInner / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Outer radius of perturbation: ', rPerturbationOuter / Kilometer, ' km'

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

    ALLOCATE( D(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( V(1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( P(1:nNodesX(1),iX_B1(1):iX_E1(1)) )

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
          X1_2      = X1  - dX1

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

    ! --- Compute fields, pre-shock ---

    DO iX1 = iX_E1(1), iX1_1, -1

      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        V(iNodeX1,iX1) = -SQRT( Two * GravitationalConstant * MassPNS / X1 )

        D(iNodeX1,iX1) = -AccretionRate / ( FourPi * X1**2 * V(iNodeX1,iX1) )

        P(iNodeX1,iX1) = D(iNodeX1,iX1) * V(iNodeX1,iX1)**2 &
                           / ( Gamma_IDEAL * MachNumber**2 )

      END DO

    END DO

    ! --- Apply jump conditions (Blondin et al., (2003), ApJ, 584, 971) ---

    V(iNodeX1_2,iX1_2) &
      = - ( Gamma_IDEAL - One ) / ( Gamma_IDEAL + One ) &
          * ( Two * GravitationalConstant * MassPNS )**( Half ) &
          * ShockRadius**( -Half )

    D(iNodeX1_2,iX1_2) &
      = ( Gamma_IDEAL + One ) / ( Gamma_IDEAL - One ) &
          * ( AccretionRate / FourPi ) &
          * ( Two * GravitationalConstant * MassPNS )**( -Half ) &
          * ShockRadius**( -Three / Two )

    P(iNodeX1_2,iX1_2) &
      = Two / ( Gamma_IDEAL + One ) &
          * ( AccretionRate / FourPi ) &
          * ( Two * GravitationalConstant * MassPNS )**( Half ) &
          * ShockRadius**( -Five / Two )

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
        'Compression Ratio LOG10(D_2/D_1) = ', &
        LOG( D(iNodeX1_2,iX1_2) / D(iNodeX1_1,iX1_1) ) / LOG( 10.0_AR )
      WRITE(*,*)

    END IF

    ! --- Compute fields, post-shock ---

    V0 = V(iNodeX1_2,iX1_2)

    DO iX1 = iX1_2, iX_B1(1), -1

      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .GT. ShockRadius ) CYCLE

        V0 = ABS( V0 )

        CALL NewtonRaphson_PostShockVelocity &
               ( X1, alpha, MassPNS, ShockRadius, V0 )

        V(iNodeX1,iX1) = -V0

        D(iNodeX1,iX1) &
          = ( ( Gamma_IDEAL - One ) &
              / ( Two * Gamma_IDEAL * PolytropicConstant ) &
              * ( Two * GravitationalConstant * MassPNS / X1 &
                    - V(iNodeX1,iX1)**2 ) &
            )**( One / ( Gamma_IDEAL - One ) )

        P(iNodeX1,iX1) &
          = PolytropicConstant * D(iNodeX1,iX1)**( Gamma_IDEAL )

      END DO

    END DO

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

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

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
          END DO
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

    DEALLOCATE( P )
    DEALLOCATE( V )
    DEALLOCATE( D )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_SASI


  SUBROUTINE NewtonRaphson_PostShockVelocity &
    ( X1, alpha, MassPNS, ShockRadius, V )

    REAL(AR), INTENT(in)    :: X1, alpha, MassPNS, ShockRadius
    REAL(AR), INTENT(inout) :: V

    REAL(AR) :: a0, a1, a2, f, df, dV
    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION

    INTEGER,  PARAMETER :: MAX_ITER  = 10
    REAL(AR), PARAMETER :: TOLERANCE = 1.0e-15_AR

    CONVERGED = .FALSE.
    ITERATION = 0
    DO WHILE( .NOT. CONVERGED .AND. ITERATION .LT. MAX_ITER )

      ITERATION = ITERATION + 1

      a0 = -Two * GravitationalConstant * MassPNS

      a1 = alpha * ( -a0 )**( ( Gamma_IDEAL + One ) / Two ) &
             * ShockRadius**( ( Three * Gamma_IDEAL - Five ) / Two ) &
             * X1**( Three - Two * Gamma_IDEAL )

      a2 = X1

      f  = a0 + a1 * V**( One - Gamma_IDEAL ) + a2 * V**2

      df = ( One - Gamma_IDEAL ) * a1 * V**( -Gamma_IDEAL ) + Two * a2 * V

      dV = -f / df
      V  = V + dV

      IF( ABS( dV / V ) .LT. TOLERANCE ) &
        CONVERGED = .TRUE.

    END DO

  END SUBROUTINE NewtonRaphson_PostShockVelocity


  SUBROUTINE ComputeSettlingSpeed_Bisection( r, Alpha, Mass, rShock, V1 )

    REAL(AR), INTENT(in) :: r, Alpha, Mass, rShock

    LOGICAL             :: Converged
    INTEGER             :: Iter
    REAL(AR)            :: a, b, c, ab, F_a, F_b, F_c, F_0
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(AR), PARAMETER :: Tol_ab  = 1.0e-8_AR
    REAL(AR), PARAMETER :: Tol_F   = 1.0e-8_AR

    REAL(AR), INTENT(out) :: V1

    a = 1.0e-6_AR
    F_a = SettlingSpeedFun( a, r, Alpha, Mass, rShock )

    b = 1.0e-1_AR
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


  REAL(AR) FUNCTION SettlingSpeedFun( u, r, Alpha, Mass, rShock )

    REAL(AR), INTENT(in) :: u, r, Alpha, Mass, rShock

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

    REAL(AR), INTENT(in)    :: mDot, Mass, rShock
    REAL(AR), INTENT(inout) :: BernoulliConstant
    REAL(AR), INTENT(inout) :: EntropyConstant
    REAL(AR), INTENT(in)    :: R(1:,0:)
    REAL(AR), INTENT(out)   :: D(1:,0:)
    REAL(AR), INTENT(out)   :: V(1:,0:)
    REAL(AR), INTENT(out)   :: P(1:,0:)
    INTEGER, PARAMETER :: iD = 1, iV = 2, iP = 3

    LOGICAL                  :: Converged
    LOGICAL                  :: PostShockGuessed
    INTEGER                  :: iX1
    INTEGER                  :: iX_B1(3), iX_E1(3)
    INTEGER                  :: iNodeX, iNodeX1
    INTEGER                  :: Iter, INFO
    INTEGER,  DIMENSION(3)   :: IPIV
    REAL(AR), PARAMETER      :: Tol = 1.0d-10
    REAL(AR)                 :: D_Previous, V_Previous, P_Previous
    REAL(AR)                 :: E, S
    REAL(AR)                 :: dEdD, dEdP, dSdD, dSdP
    REAL(AR), DIMENSION(3)   :: FVEC, UVEC, dUVEC
    REAL(AR), DIMENSION(3,3) :: FJAC

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
          !  ABS(  dUVEC / ( UVEC + TINY(1.0_AR) ) )
          !WRITE(*,*)

          UVEC = UVEC - dUVEC

          IF( ALL( ABS( dUVEC / ( UVEC + TINY(1.0_AR) ) ) < Tol ) ) &
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

    REAL(AR), INTENT(in)  :: mDot, Mass, r
    REAL(AR), INTENT(in)  :: BernoulliConstant
    REAL(AR), INTENT(in)  :: EntropyConstant
    REAL(AR), INTENT(in)  :: D, V, P, S, E
    REAL(AR), INTENT(out) :: FVEC(3)

    FVEC(1) = FourPi * r**2 * D * V - mDot
    FVEC(2) = Half * V**2 + E + P / D - Mass / r - BernoulliConstant
    FVEC(3) = S - EntropyConstant

  END SUBROUTINE ComputeFVEC_InitialConditions


  SUBROUTINE ComputeFJAC_InitialConditions &
               ( r, D, V, P, E, &
                 dEdD, dEdP, &
                 dSdD, dSdP, &
                 FJAC )

    REAL(AR), INTENT(in)  :: r, D, V, E, P
    REAL(AR), INTENT(in)  :: dEdD, dEdP
    REAL(AR), INTENT(in)  :: dSdD, dSdP
    REAL(AR), INTENT(out) :: FJAC(3,3)

    FJAC(1,1) = FourPi * r**2 * V
    FJAC(1,2) = FourPi * r**2 * D
    FJAC(1,3) = 0.0_AR

    FJAC(2,1) = ( dEdD - ( P / D**2 ) )
    FJAC(2,2) = V
    FJAC(2,3) = ( dEdP + ( 1.0_AR / D ) )

    FJAC(3,1) = dSdD
    FJAC(3,2) = 0.0_AR
    FJAC(3,3) = dSdP

  END SUBROUTINE ComputeFJAC_InitialConditions


  SUBROUTINE ComputeLeftState( D_R, V_R, P_R, D_L, V_L, P_L )

    LOGICAL                  :: Converged
    INTEGER                  :: Iter, INFO
    INTEGER,  DIMENSION(3)   :: IPIV
    REAL(AR), PARAMETER      :: Tol = 1.0d-10
    REAL(AR), INTENT(inout)  :: D_L, V_L, P_L
    REAL(AR) :: E_L
    REAL(AR), INTENT(inout)  :: D_R, V_R, P_R
    REAL(AR) :: E_R
    REAL(AR) :: dEdD_L, dEdP_L
    REAL(AR), DIMENSION(3)   :: FVEC, UVEC, dUVEC
    REAL(AR), DIMENSION(3,3) :: FJAC

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
      !  ABS(  dUVEC / ( UVEC + TINY(1.0_AR) ) )
      !WRITE(*,*)

      D_L = UVEC(1) - dUVEC(1)
      V_L = UVEC(2) - dUVEC(2)
      P_L = UVEC(3) - dUVEC(3)

      IF( ALL( ABS(  dUVEC / ( UVEC + TINY(1.0_AR) ) ) < Tol ) ) &
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

    REAL(AR),               INTENT(in)  :: D_R, V_R, P_R, E_R
    REAL(AR),               INTENT(in)  :: D_L, V_L, P_L, E_L
    REAL(AR), DIMENSION(3), INTENT(out) :: FVEC

    FVEC(1) = D_L * V_L - D_R * V_R
    FVEC(2) = D_L * V_L**2 - D_R * V_R**2 + P_L - P_R
    FVEC(3) = E_L - E_R + Half * ( V_L**2 - V_R**2 ) + P_L / D_L - P_R / D_R

  END SUBROUTINE ComputeFVEC_LeftState


  SUBROUTINE ComputeFJAC_LeftState &
               ( D_R, V_R, P_R, E_R, D_L, V_L, P_L, E_L, &
                 dEdD_L, dEdP_L, FJAC )

    REAL(AR),                 INTENT(in)  :: D_R, V_R, P_R, E_R
    REAL(AR),                 INTENT(in)  :: D_L, V_L, P_L, E_L
    REAL(AR),                 INTENT(in)  :: dEdD_L, dEdP_L
    REAL(AR), DIMENSION(3,3), INTENT(out) :: FJAC

    FJAC(1,1) = V_L
    FJAC(1,2) = D_L
    FJAC(1,3) = 0.0_AR

    FJAC(2,1) = V_L**2
    FJAC(2,2) = 2.0_AR * D_L * V_L
    FJAC(2,3) = One

    FJAC(3,1) = dEdD_L - ( P_L / D_L**2 )
    FJAC(3,2) = V_L
    FJAC(3,3) = dEdP_L + ( One / D_L )

  END SUBROUTINE ComputeFJAC_LeftState


END MODULE MF_InitializationModule_NonRelativistic_IDEAL
