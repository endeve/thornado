MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_min

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    swX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL, &
    Gamma_IDEAL
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi, &
    nGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE MagnetofluidFieldsModule, ONLY: &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nCM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nPM
  USE MHD_UtilitiesModule, ONLY: &
    ComputeConserved_MHD
  USE UnitsModule, ONLY: &
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
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Pi, &
    TwoPi
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X, &
    amrex2thornado_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    UseTiling, &
    EvolveOnlyMagnetic, &
    ApplyPerturbations, &
    Rand_Amplitude
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF  

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: HDFERR

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'Advection1D' )

        CALL InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCM )

      CASE( 'RiemannProblem1D' )

        CALL InitializeFields_RiemannProblem1D( iLevel, MF_uGF, MF_uCM )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D( iLevel, MF_uGF, MF_uCM )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCM )

      CASE( 'KelvinHelmholtz2D' )

        CALL InitializeFields_KelvinHelmholtz2D( iLevel, MF_uGF, MF_uCM )

      CASE( 'Advection3D' )

        CALL InitializeFields_Advection3D( iLevel, MF_uGF, MF_uCM )

      CASE( 'ShearingDisk' )

        CALL InitializeFields_ShearingDisk( iLevel, MF_uGF, MF_uCM )

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 201, Message_Option &
                        = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P, Ne
    REAL(DP) :: Amp
    REAL(DP) :: X1_0
    REAL(DP) :: sigmaX1

    uPM = Zero

    AdvectionProfile = 'SineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( AdvectionProfile ) )

      CASE( 'SineWaveX1' )

        D_0 = 1.0_DP
        Amp = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') 'D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') 'Amp: ', Amp
          WRITE(*,'(8x,A,F5.3)') ' V1: ', V1
          WRITE(*,'(8x,A,F5.3)') ' V2: ', V2
          WRITE(*,'(8x,A,F5.3)') ' V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '  P: ', P
          WRITE(*,*)

        END IF

      CASE( 'Gaussian' )

        D_0     = 1.0_DP
        X1_0    = 0.5_DP
        sigmaX1 = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') '    D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') '   X1_0: ', X1_0
          WRITE(*,'(8x,A,F5.3)') 'sigmaX1: ', sigmaX1
          WRITE(*,'(8x,A,F5.3)') '     V1: ', V1
          WRITE(*,'(8x,A,F5.3)') '     V2: ', V2
          WRITE(*,'(8x,A,F5.3)') '     V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '      P: ', P
          WRITE(*,*)

        END IF

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 202, Message_Option &
                        = 'Invalid AdvectionProfile: ' &
                            // TRIM( AdvectionProfile ) )

    END SELECT

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

          uPM(iNX,iPM_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Ne

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uPM(iNX,iPM_D) &
            = D_0 * EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) )
          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Ne

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection1D


  SUBROUTINE InitializeFields_RiemannProblem1D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1, X1_D

    CHARACTER(:), ALLOCATABLE :: RiemannProblemName

    REAL(DP) :: uPM_L(nPM), uPM_R(nPM)

    RiemannProblemName = 'Sod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        X1_D = 0.5_DP

        uPM_L(iPM_D ) = 1.0_DP
        uPM_L(iPM_V1) = 0.0_DP
        uPM_L(iPM_V2) = 0.0_DP
        uPM_L(iPM_V3) = 0.0_DP
        uPM_L(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_L(iPM_Ne) = 0.0_DP

        uPM_R(iPM_D ) = 0.125_DP
        uPM_R(iPM_V1) = 0.0_DP
        uPM_R(iPM_V2) = 0.0_DP
        uPM_R(iPM_V3) = 0.0_DP
        uPM_R(iPM_E ) = 0.1_DP / ( Gamma_IDEAL - One )
        uPM_R(iPM_Ne) = 0.0_DP

      CASE( 'SphericalSod' )

        X1_D = 1.0_DP

        uPM_L(iPM_D ) = 1.0_DP
        uPM_L(iPM_V1) = 0.0_DP
        uPM_L(iPM_V2) = 0.0_DP
        uPM_L(iPM_V3) = 0.0_DP
        uPM_L(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_L(iPM_Ne) = 0.0_DP

        uPM_R(iPM_D ) = 0.125_DP
        uPM_R(iPM_V1) = 0.0_DP
        uPM_R(iPM_V2) = 0.0_DP
        uPM_R(iPM_V3) = 0.0_DP
        uPM_R(iPM_E ) = 0.1_DP / ( Gamma_IDEAL - One )
        uPM_R(iPM_Ne) = 0.0_DP

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 203, Message_Option &
                        = 'Invalid RiemannProblemName: ' &
                            // TRIM( RiemannProblemName ) )

    END SELECT

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A,A)') 'RiemannProblemName: ', TRIM( RiemannProblemName )
      WRITE(*,*)
      WRITE(*,'(8x,A,F5.3)') 'X1_D: ', X1_D
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'Left State'
      WRITE(*,'(8x,A)')      '----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_L(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_L(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_L(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_L(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_L(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_L(iPM_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'Right State'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_R(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_R(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_R(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_R(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_R(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_R(iPM_Ne)
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( X1 .LT. X1_D )THEN

          uPM(iNX,iPM_D ) = uPM_L(iPM_D )
          uPM(iNX,iPM_V1) = uPM_L(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_L(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_L(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_L(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_L(iPM_Ne)

        ELSE

          uPM(iNX,iPM_D ) = uPM_R(iPM_D )
          uPM(iNX,iPM_V1) = uPM_R(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_R(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_R(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_R(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_R(iPM_Ne)

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_RiemannProblem1D


  SUBROUTINE InitializeFields_RiemannProblem2D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2, X1_D, X2_D

    CHARACTER(:), ALLOCATABLE :: RiemannProblemName

    REAL(DP) :: uPM_NW(nPM), uPM_NE(nPM), uPM_SW(nPM), uPM_SE(nPM)

    RiemannProblemName = 'Sod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        X1_D = 0.5_DP
        X2_D = 0.5_DP

        uPM_SW(iPM_D ) = 1.0_DP
        uPM_SW(iPM_V1) = 0.0_DP
        uPM_SW(iPM_V2) = 0.0_DP
        uPM_SW(iPM_V3) = 0.0_DP
        uPM_SW(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_SW(iPM_Ne) = 0.0_DP

        uPM_SE(iPM_D ) = 0.125_DP
        uPM_SE(iPM_V1) = 0.0_DP
        uPM_SE(iPM_V2) = 0.0_DP
        uPM_SE(iPM_V3) = 0.0_DP
        uPM_SE(iPM_E ) = 0.1_DP / ( Gamma_IDEAL - One )
        uPM_SE(iPM_Ne) = 0.0_DP

        uPM_NE = uPM_SE
        uPM_NW = uPM_SE

      CASE( 'dZB2002' )

        X1_D = 0.5_DP
        X2_D = 0.5_DP

        uPM_NE(iPM_D ) = 0.1_DP
        uPM_NE(iPM_V1) = 0.0_DP
        uPM_NE(iPM_V2) = 0.0_DP
        uPM_NE(iPM_V3) = 0.0_DP
        uPM_NE(iPM_E ) = 0.01_DP / ( Gamma_IDEAL - One )
        uPM_NE(iPM_Ne) = 0.0_DP

        uPM_NW(iPM_D ) = 0.1_DP
        uPM_NW(iPM_V1) = 0.99_DP
        uPM_NW(iPM_V2) = 0.0_DP
        uPM_NW(iPM_V3) = 0.0_DP
        uPM_NW(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_NW(iPM_Ne) = 0.0_DP

        uPM_SW(iPM_D ) = 0.5_DP
        uPM_SW(iPM_V1) = 0.0_DP
        uPM_SW(iPM_V2) = 0.0_DP
        uPM_SW(iPM_V3) = 0.0_DP
        uPM_SW(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_SW(iPM_Ne) = 0.0_DP

        uPM_SE(iPM_D ) = 0.1_DP
        uPM_SE(iPM_V1) = 0.0_DP
        uPM_SE(iPM_V2) = 0.99_DP
        uPM_SE(iPM_V3) = 0.0_DP
        uPM_SE(iPM_E ) = 1.0_DP / ( Gamma_IDEAL - One )
        uPM_SE(iPM_Ne) = 0.0_DP

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 203, Message_Option &
                        = 'Invalid RiemannProblemName: ' &
                            // TRIM( RiemannProblemName ) )

    END SELECT

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A,A)') 'RiemannProblemName: ', TRIM( RiemannProblemName )
      WRITE(*,*)
      WRITE(*,'(8x,A,F5.3)') 'X1_D: ', X1_D
      WRITE(*,'(8x,A,F5.3)') 'X2_D: ', X2_D
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'NW Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_NW(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_NW(iPM_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_NW(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_NW(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_NW(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_NW(iPM_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'NE Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_NE(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_NE(iPM_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_NE(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_NE(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_NE(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_NE(iPM_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'SW Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_SW(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_SW(iPM_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_SW(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_SW(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_SW(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_SW(iPM_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'SE Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_D ): ', uPM_SE(iPM_D )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V1): ', uPM_SE(iPM_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V2): ', uPM_SE(iPM_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_V3): ', uPM_SE(iPM_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_E ): ', uPM_SE(iPM_E )
      WRITE(*,'(8x,A,F5.3)') 'uPM(iPM_Ne): ', uPM_SE(iPM_Ne)
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2   = NodeCoordinate( MeshX(2), iX2, iNX2 )

        IF( X1 .LT. X1_D .AND. X2 .GT. X2_D )THEN

          uPM(iNX,iPM_D ) = uPM_NW(iPM_D )
          uPM(iNX,iPM_V1) = uPM_NW(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_NW(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_NW(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_NW(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_NW(iPM_Ne)

        ELSE IF( X1 .GT. X1_D .AND. X2 .GT. X2_D )THEN

          uPM(iNX,iPM_D ) = uPM_NE(iPM_D )
          uPM(iNX,iPM_V1) = uPM_NE(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_NE(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_NE(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_NE(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_NE(iPM_Ne)

        ELSE IF( X1 .LT. X1_D .AND. X2 .LT. X2_D )THEN

          uPM(iNX,iPM_D ) = uPM_SW(iPM_D )
          uPM(iNX,iPM_V1) = uPM_SW(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_SW(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_SW(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_SW(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_SW(iPM_Ne)

        ELSE

          uPM(iNX,iPM_D ) = uPM_SE(iPM_D )
          uPM(iNX,iPM_V1) = uPM_SE(iPM_V1)
          uPM(iNX,iPM_V2) = uPM_SE(iPM_V2)
          uPM(iNX,iPM_V3) = uPM_SE(iPM_V3)
          uPM(iNX,iPM_E ) = uPM_SE(iPM_E )
          uPM(iNX,iPM_Ne) = uPM_SE(iPM_Ne)

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P, Ne
    REAL(DP) :: Amp
    REAL(DP) :: X1_0, X2_0
    REAL(DP) :: sigmaX1, sigmaX2

    AdvectionProfile = 'SineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( AdvectionProfile ) )

      CASE( 'SineWaveX1' )

        D_0 = 1.0_DP
        Amp = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') 'D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') 'Amp: ', Amp
          WRITE(*,'(8x,A,F5.3)') ' V1: ', V1
          WRITE(*,'(8x,A,F5.3)') ' V2: ', V2
          WRITE(*,'(8x,A,F5.3)') ' V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '  P: ', P
          WRITE(*,*)

        END IF

      CASE( 'Gaussian' )

        D_0     = 1.0_DP
        X1_0    = 0.5_DP
        X2_0    = 0.5_DP
        sigmaX1 = 0.1_DP
        sigmaX2 = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') '    D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') '   X1_0: ', X1_0
          WRITE(*,'(8x,A,F5.3)') '   X2_0: ', X2_0
          WRITE(*,'(8x,A,F5.3)') 'sigmaX1: ', sigmaX1
          WRITE(*,'(8x,A,F5.3)') 'sigmaX2: ', sigmaX2
          WRITE(*,'(8x,A,F5.3)') '     V1: ', V1
          WRITE(*,'(8x,A,F5.3)') '     V2: ', V2
          WRITE(*,'(8x,A,F5.3)') '     V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '      P: ', P
          WRITE(*,*)

        END IF

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 204, Message_Option &
                        = 'Invalid AdvectionProfile: ' &
                            // TRIM( AdvectionProfile ) )

    END SELECT

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        iNX2 = NodeNumberTableX(2,iNX)
        X2   = NodeCoordinate( MeshX(2), iX2, iNX2 )

        IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

          uPM(iNX,iPM_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Ne

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uPM(iNX,iPM_D) &
            = D_0 * EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) ) &
                  * EXP( -( X2 - X2_0 )**2 / ( Two * sigmaX2**2 ) )
          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Ne

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection2D


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Radice & Rezzolla, (2012), AA, 547, A26 ---
  SUBROUTINE InitializeFields_KelvinHelmholtz2D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific Parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2

    REAL(DP) :: a      = 0.01_DP
    REAL(DP) :: Vshear = Half
    REAL(DP) :: A0     = 0.1_DP ! --- Perturbation amplitude ---
    REAL(DP) :: sigma  = 0.1_DP
    REAL(DP) :: rho0   = 0.505_DP
    REAL(DP) :: rho1   = 0.495_DP

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A)')   'Kelvin--Helmholtz'
      WRITE(*,'(6x,A,A)') '-----------------'
      WRITE(*,*)
      WRITE(*,'(8x,A,F5.3)') '     a: ', a
      WRITE(*,'(8x,A,F5.3)') 'Vshear: ', Vshear
      WRITE(*,'(8x,A,F5.3)') ' sigma: ', sigma
      WRITE(*,'(8x,A,F5.3)') '  rho0: ', rho0
      WRITE(*,'(8x,A,F5.3)') '  rho1: ', rho1
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        ! --- V1 ---
        IF( X2 .GT. Zero )THEN

          uPM(iNX,iPM_V1) &
            = +Vshear * TANH( ( X2 - Half ) / a )

        ELSE

          ! --- Paper has a typo here, the minus sign is required ---
          uPM(iNX,iPM_V1) &
            = -Vshear * TANH( ( X2 + Half ) / a )

        END IF

        ! --- V2 ---
        IF( X2 .GT. Zero )THEN

          uPM(iNX,iPM_V2) &
            =  A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 - Half )**2 / sigma**2 ) )

        ELSE

          uPM(iNX,iPM_V2) &
            = -A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 + Half )**2 / sigma**2 ) )

        END IF

        ! --- rho ---
        IF( X2 .GT. Zero )THEN

          uPM(iNX,iPM_D) &
            = rho0 + rho1 * TANH( ( X2 - Half ) / a )

        ELSE

          uPM(iNX,iPM_D) &
            = rho0 - rho1 * TANH( ( X2 + Half ) / a )

        END IF

        uPM(iNX,iPM_V3) = Zero
        uPM(iNX,iPM_E ) = One / ( Gamma_IDEAL - One )
        uPM(iNX,iPM_Ne) = Zero

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_KelvinHelmholtz2D


  SUBROUTINE InitializeFields_Advection3D( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2, iNX3
    REAL(DP) :: X1, X2, X3

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P, Ne
    REAL(DP) :: Amp
    REAL(DP) :: X1_0, X2_0, X3_0
    REAL(DP) :: sigmaX1, sigmaX2, sigmaX3

    AdvectionProfile = 'SineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( AdvectionProfile ) )

      CASE( 'SineWaveX1' )

        D_0 = 1.0_DP
        Amp = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') 'D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') 'Amp: ', Amp
          WRITE(*,'(8x,A,F5.3)') ' V1: ', V1
          WRITE(*,'(8x,A,F5.3)') ' V2: ', V2
          WRITE(*,'(8x,A,F5.3)') ' V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '  P: ', P
          WRITE(*,*)

        END IF

      CASE( 'Gaussian' )

        D_0     = 1.0_DP
        X1_0    = 0.5_DP
        X2_0    = 0.5_DP
        X3_0    = 0.5_DP
        sigmaX1 = 0.1_DP
        sigmaX2 = 0.1_DP
        sigmaX3 = 0.1_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,'(6x,A,A)') '------------------ '
          WRITE(*,*)
          WRITE(*,'(8x,A,F5.3)') '    D_0: ', D_0
          WRITE(*,'(8x,A,F5.3)') '   X1_0: ', X1_0
          WRITE(*,'(8x,A,F5.3)') '   X2_0: ', X2_0
          WRITE(*,'(8x,A,F5.3)') '   X3_0: ', X3_0
          WRITE(*,'(8x,A,F5.3)') 'sigmaX1: ', sigmaX1
          WRITE(*,'(8x,A,F5.3)') 'sigmaX2: ', sigmaX2
          WRITE(*,'(8x,A,F5.3)') 'sigmaX3: ', sigmaX3
          WRITE(*,'(8x,A,F5.3)') '     V1: ', V1
          WRITE(*,'(8x,A,F5.3)') '     V2: ', V2
          WRITE(*,'(8x,A,F5.3)') '     V3: ', V3
          WRITE(*,'(8x,A,F5.3)') '      P: ', P
          WRITE(*,*)

        END IF

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 205, Message_Option &
                        = 'Invalid AdvectionProfile: ' &
                            // TRIM( AdvectionProfile ) )

    END SELECT

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        iNX2 = NodeNumberTableX(2,iNX)
        X2   = NodeCoordinate( MeshX(2), iX2, iNX2 )

        iNX3 = NodeNumberTableX(3,iNX)
        X3   = NodeCoordinate( MeshX(3), iX3, iNX3 )

        IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

          uPM(iNX,iPM_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Zero

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uPM(iNX,iPM_D) &
            = D_0 * EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) ) &
                  * EXP( -( X2 - X2_0 )**2 / ( Two * sigmaX2**2 ) ) &
                  * EXP( -( X3 - X3_0 )**2 / ( Two * sigmaX3**2 ) )

          uPM(iNX,iPM_V1) = V1
          uPM(iNX,iPM_V2) = V2
          uPM(iNX,iPM_V3) = V3
          uPM(iNX,iPM_E ) = P / ( Gamma_IDEAL - One )
          uPM(iNX,iPM_Ne) = Zero

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection3D


  SUBROUTINE InitializeFields_ShearingDisk( iLevel, MF_uGF, MF_uCM )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPM(nDOFX,nPM)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific Parameters ---

    INTEGER  :: iNX1, iNX2, iNX3
    REAL(DP) :: X1, X2, X3

    CHARACTER(256) :: FileName

    INTEGER        :: nX_Data
    INTEGER(HID_T) :: FILE_ID

    REAL(DP) :: CB1, CB2, CB3, V1, V2, V3, VSq, W, VdotB
    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), V3Arr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    REAL(DP) :: kz
    REAL(DP) :: Rand_r, Rand_z, Rand_theta
    REAL(DP) :: Random_r, Random_z, Random_theta

    uPM = Zero

    FileName &
      = "/home/jbuffal/thornado_MHD_3D/Workflow/MHD/ShearingDisk/GR_LR_diffrot.h5"
    nX_Data = 10000

    ! --- Populate arrays ---

    CALL H5OPEN_F( HDFERR )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    ALLOCATE( PressureArr(nX_Data), DensityArr(nX_Data), V3Arr(nX_Data), &
              AlphaArr(nX_Data), PsiArr(nX_Data), X1Arr(nX_Data) )

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

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A)')   'Shearing Disk'
      WRITE(*,'(6x,A,A)') '-----------------'
      WRITE(*,*)
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)
        iNX3 = NodeNumberTableX(3,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

        G(iNX,iX1,iX2,iX3,iGF_Alpha) &
           = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

        G(iNX,iX1,iX2,iX3,iGF_Psi) &
             = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

        G(iNX,iX1,iX2,iX3,iGF_h_1) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2
        G(iNX,iX1,iX2,iX3,iGF_h_2) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2
        G(iNX,iX1,iX2,iX3,iGF_h_3) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2 * X1

        CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

        G(iNX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        G(iNX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        G(iNX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        uPM(iNX,iPM_D) &
          = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

        V1 = Zero
        V2 = Zero
        V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

        VSq = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * V1**2 &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * V2**2 &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * V3**2

        W = One / SQRT( One - VSq )

        uPM(iNX,iPM_V1) = V1
        uPM(iNX,iPM_V2) = V2
        uPM(iNX,iPM_V3) = V3

        IF( ApplyPerturbations )THEN

          PRINT*, 'Applying perturbations.'

          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER( Rand_r )

          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER( Rand_z )

          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER( Rand_theta )

          Random_r       =  Two * Rand_r - One
          Random_z       =  Two * Rand_z - One
          Random_theta   =  Two * Rand_theta - One

          kz = Two * Pi / ( Half * Kilometer )

          uPM(iNX,iPM_V1) = ( 0.1_DP * Rand_Amplitude * Random_r &
                              + 0.2d-5 * SIN( kz * X2 ) ) &
                            * X1 * uPM(iNX,iPM_V3) 
          uPM(iNX,iPM_V2) = Rand_Amplitude * Random_z &
                            * X1 * uPM(iNX,iPM_V3) 
          uPM(iNX,iPM_V3) = ( One + Rand_Amplitude * Random_theta ) &
                            * uPM(iNX,iPM_V3)

        END IF

        uPM(iNX,iPM_E) &
          = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
            / ( Gamma_IDEAL - One )

        CB1 = Zero
        CB2 = 2.0 * 1.0d13 * Gauss
        CB3 = Zero

        VdotB = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

        uPM(iNX,iPM_B1) = W * VdotB * V1 + CB1 / W
        uPM(iNX,iPM_B2) = W * VdotB * V2 + CB2 / W
        uPM(iNX,iPM_B3) = W * VdotB * V3 + CB3 / W

        uPM(iNX,iPM_Chi) = Zero

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPM(iNX,iPM_D ), uPM(iNX,iPM_E), &
                 uPM(iNX,iPM_Ne), Pressure )

        CALL ComputeConserved_MHD &
               ( uPM(iNX,iPM_D  ), &
                 uPM(iNX,iPM_V1 ), &
                 uPM(iNX,iPM_V2 ), &
                 uPM(iNX,iPM_V3 ), &
                 uPM(iNX,iPM_E  ), &
                 uPM(iNX,iPM_Ne ), &
                 uPM(iNX,iPM_B1 ), &
                 uPM(iNX,iPM_B2 ), &
                 uPM(iNX,iPM_B3 ), &
                 uPM(iNX,iPM_Chi), &
                 U  (iNX,iX1,iX2,iX3,iCM_D  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_S3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_E  ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Ne ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B1 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B2 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_B3 ), &
                 U  (iNX,iX1,iX2,iX3,iCM_Chi), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 G  (iNX,iX1,iX2,iX3,iGF_Alpha), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_1), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_2), &
                 G  (iNX,iX1,iX2,iX3,iGF_Beta_3), &
                 Pressure, &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( Zero, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, LBOUND( uCM ), iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

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


END MODULE MF_InitializationModule
