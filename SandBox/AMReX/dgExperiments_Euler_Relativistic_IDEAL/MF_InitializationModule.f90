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
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nAF, &
    iAF_P
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    epsMin_Euler_GR

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    TwoPi
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X, &
    amrex2thornado_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    swX, &
    Gamma_IDEAL, &
    ProgramName, &
    UseTiling
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'Advection1D' )

        CALL InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCF )

      CASE( 'RiemannProblem1D' )

        CALL InitializeFields_RiemannProblem1D( iLevel, MF_uGF, MF_uCF )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D( iLevel, MF_uGF, MF_uCF )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz2D' )

        CALL InitializeFields_KelvinHelmholtz2D( iLevel, MF_uGF, MF_uCF )

      CASE( 'Advection3D' )

        CALL InitializeFields_Advection3D( iLevel, MF_uGF, MF_uCF )

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 201, Message_Option &
                        = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P
    REAL(DP) :: Amp
    REAL(DP) :: X1_0
    REAL(DP) :: sigmaX1

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
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D) &
            = D_0 + EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) )

          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection1D


  SUBROUTINE InitializeFields_RiemannProblem1D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1, X1_D

    CHARACTER(:), ALLOCATABLE :: RiemannProblemName

    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uAF_L(nAF), uAF_R(nAF)

    RiemannProblemName = 'Sod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        X1_D = 0.5_DP

        uAF_L(iAF_P ) = 1.0_DP

        uPF_L(iPF_D ) = 1.0_DP
        uPF_L(iPF_V1) = 0.0_DP
        uPF_L(iPF_V2) = 0.0_DP
        uPF_L(iPF_V3) = 0.0_DP
        uPF_L(iPF_E ) = uAF_L(iAF_P) / ( Gamma_IDEAL - One )
        uPF_L(iPF_Ne) = 0.0_DP

        uAF_R(iAF_P ) = 0.1_DP

        uPF_R(iPF_D ) = 0.125_DP
        uPF_R(iPF_V1) = 0.0_DP
        uPF_R(iPF_V2) = 0.0_DP
        uPF_R(iPF_V3) = 0.0_DP
        uPF_R(iPF_E ) = uAF_R(iAF_P) / ( Gamma_IDEAL - One )
        uPF_R(iPF_Ne) = 0.0_DP

      CASE( 'SphericalSod' )

        X1_D = 1.0_DP

        uAF_L(iAF_P ) = 1.0_DP

        uPF_L(iPF_D ) = 1.0_DP
        uPF_L(iPF_V1) = 0.0_DP
        uPF_L(iPF_V2) = 0.0_DP
        uPF_L(iPF_V3) = 0.0_DP
        uPF_L(iPF_E ) = uAF_L(iAF_P) / ( Gamma_IDEAL - One )
        uPF_L(iPF_Ne) = 0.0_DP

        uAF_R(iAF_P ) = 0.1_DP

        uPF_R(iPF_D ) = 0.125_DP
        uPF_R(iPF_V1) = 0.0_DP
        uPF_R(iPF_V2) = 0.0_DP
        uPF_R(iPF_V3) = 0.0_DP
        uPF_R(iPF_E ) = uAF_R(iAF_P) / ( Gamma_IDEAL - One )
        uPF_R(iPF_Ne) = 0.0_DP

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
      WRITE(*,'(8x,A)')     'Left State'
      WRITE(*,'(8x,A)')     '----------'
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_L(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_L(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_L(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_L(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_L(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_L(iPF_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')     'Right State'
      WRITE(*,'(8x,A)')     '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_R(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_R(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_R(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_R(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_R(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_R(iPF_Ne)
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uPF(iNX,iPF_D ) = uPF_L(iPF_D )
          uPF(iNX,iPF_V1) = uPF_L(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_L(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_L(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_L(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_L(iPF_Ne)

          uAF(iNX,iAF_P) = uAF_L(iAF_P)

        ELSE

          uPF(iNX,iPF_D ) = uPF_R(iPF_D )
          uPF(iNX,iPF_V1) = uPF_R(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_R(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_R(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_R(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_R(iPF_Ne)

          uAF(iNX,iAF_P) = uAF_R(iAF_P)

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_RiemannProblem1D


  SUBROUTINE InitializeFields_RiemannProblem2D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2, X1_D, X2_D

    CHARACTER(:), ALLOCATABLE :: RiemannProblemName

    REAL(DP) :: uPF_NW(nPF), uPF_NE(nPF), uPF_SW(nPF), uPF_SE(nPF)
    REAL(DP) :: uAF_NW(nAF), uAF_NE(nAF), uAF_SW(nAF), uAF_SE(nAF)

    RiemannProblemName = 'Sod'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RiemannProblemName', RiemannProblemName )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        X1_D = 0.5_DP
        X2_D = 0.5_DP

        uAF_SW(iAF_P ) = 1.0_DP
        uPF_SW(iPF_D ) = 1.0_DP
        uPF_SW(iPF_V1) = 0.0_DP
        uPF_SW(iPF_V2) = 0.0_DP
        uPF_SW(iPF_V3) = 0.0_DP
        uPF_SW(iPF_E ) = uAF_SW(iAF_P) / ( Gamma_IDEAL - One )
        uPF_SW(iPF_Ne) = 0.0_DP

        uAF_SE(iAF_P ) = 0.1_DP
        uPF_SE(iPF_D ) = 0.125_DP
        uPF_SE(iPF_V1) = 0.0_DP
        uPF_SE(iPF_V2) = 0.0_DP
        uPF_SE(iPF_V3) = 0.0_DP
        uPF_SE(iPF_E ) = uAF_SE(iAF_P) / ( Gamma_IDEAL - One )
        uPF_SE(iPF_Ne) = 0.0_DP

        uPF_NE = uPF_SE
        uPF_NW = uPF_SE
        uAF_NE = uAF_SE
        uAF_NW = uAF_SE

      CASE( 'dZB2002' )

        X1_D = 0.5_DP
        X2_D = 0.5_DP

        uPF_NE(iPF_D ) = 0.1_DP
        uPF_NE(iPF_V1) = 0.0_DP
        uPF_NE(iPF_V2) = 0.0_DP
        uPF_NE(iPF_V3) = 0.0_DP
        uAF_NE(iAF_P ) = 0.01_DP
        uPF_NE(iPF_E ) = uAF_NE(iAF_P) / ( Gamma_IDEAL - One )
        uPF_NE(iPF_Ne) = 0.0_DP

        uPF_NW(iPF_D ) = 0.1_DP
        uPF_NW(iPF_V1) = 0.99_DP
        uPF_NW(iPF_V2) = 0.0_DP
        uPF_NW(iPF_V3) = 0.0_DP
        uAF_NW(iAF_P ) = 1.0_DP
        uPF_NW(iPF_E ) = uAF_NW(iAF_P) / ( Gamma_IDEAL - One )
        uPF_NW(iPF_Ne) = 0.0_DP

        uPF_SW(iPF_D ) = 0.5_DP
        uPF_SW(iPF_V1) = 0.0_DP
        uPF_SW(iPF_V2) = 0.0_DP
        uPF_SW(iPF_V3) = 0.0_DP
        uAF_SW(iAF_P ) = 1.0_DP
        uPF_SW(iPF_E ) = uAF_SW(iAF_P) / ( Gamma_IDEAL - One )
        uPF_SW(iPF_Ne) = 0.0_DP

        uPF_SE(iPF_D ) = 0.1_DP
        uPF_SE(iPF_V1) = 0.0_DP
        uPF_SE(iPF_V2) = 0.99_DP
        uPF_SE(iPF_V3) = 0.0_DP
        uAF_SE(iAF_P ) = 1.0_DP
        uPF_SE(iPF_E ) = uAF_SE(iAF_P) / ( Gamma_IDEAL - One )
        uPF_SE(iPF_Ne) = 0.0_DP

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
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_NW(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_NW(iPF_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_NW(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_NW(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_NW(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_NW(iPF_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'NE Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_NE(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_NE(iPF_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_NE(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_NE(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_NE(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_NE(iPF_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'SW Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_SW(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_SW(iPF_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_SW(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_SW(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_SW(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_SW(iPF_Ne)
      WRITE(*,*)
      WRITE(*,'(8x,A)')      'SE Quadrant'
      WRITE(*,'(8x,A)')      '-----------'
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_D ): ', uPF_SE(iPF_D )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V1): ', uPF_SE(iPF_V1)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V2): ', uPF_SE(iPF_V2)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_V3): ', uPF_SE(iPF_V3)
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_E ): ', uPF_SE(iPF_E )
      WRITE(*,'(8x,A,F5.3)') 'uPF(iPF_Ne): ', uPF_SE(iPF_Ne)
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uPF(iNX,iPF_D ) = uPF_NW(iPF_D )
          uPF(iNX,iPF_V1) = uPF_NW(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_NW(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_NW(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_NW(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_NW(iPF_Ne)

          uAF(iNX,iAF_P ) = uAF_NW(iAF_P)

        ELSE IF( X1 .GT. X1_D .AND. X2 .GT. X2_D )THEN

          uPF(iNX,iPF_D ) = uPF_NE(iPF_D )
          uPF(iNX,iPF_V1) = uPF_NE(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_NE(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_NE(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_NE(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_NE(iPF_Ne)

          uAF(iNX,iAF_P ) = uAF_NE(iAF_P)

        ELSE IF( X1 .LT. X1_D .AND. X2 .LT. X2_D )THEN

          uPF(iNX,iPF_D ) = uPF_SW(iPF_D )
          uPF(iNX,iPF_V1) = uPF_SW(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_SW(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_SW(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_SW(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_SW(iPF_Ne)

          uAF(iNX,iAF_P ) = uAF_SW(iAF_P)

        ELSE

          uPF(iNX,iPF_D ) = uPF_SE(iPF_D )
          uPF(iNX,iPF_V1) = uPF_SE(iPF_V1)
          uPF(iNX,iPF_V2) = uPF_SE(iPF_V2)
          uPF(iNX,iPF_V3) = uPF_SE(iPF_V3)
          uPF(iNX,iPF_E ) = uPF_SE(iPF_E )
          uPF(iNX,iPF_Ne) = uPF_SE(iPF_Ne)

          uAF(iNX,iAF_P ) = uAF_SE(iAF_P)

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P
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
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D) &
            = D_0 + EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) ) &
                  * EXP( -( X2 - X2_0 )**2 / ( Two * sigmaX2**2 ) )

          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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
  SUBROUTINE InitializeFields_KelvinHelmholtz2D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

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

    epsMin_Euler_GR = HUGE( One )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uPF(iNX,iPF_V1) &
            = +Vshear * TANH( ( X2 - Half ) / a )

        ELSE

          ! --- Paper has a typo here, the minus sign is required ---
          uPF(iNX,iPF_V1) &
            = -Vshear * TANH( ( X2 + Half ) / a )

        END IF

        ! --- V2 ---
        IF( X2 .GT. Zero )THEN

          uPF(iNX,iPF_V2) &
            =  A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 - Half )**2 / sigma**2 ) )

        ELSE

          uPF(iNX,iPF_V2) &
            = -A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 + Half )**2 / sigma**2 ) )

        END IF

        ! --- rho ---
        IF( X2 .GT. Zero )THEN

          uPF(iNX,iPF_D) &
            = rho0 + rho1 * TANH( ( X2 - Half ) / a )

        ELSE

          uPF(iNX,iPF_D) &
            = rho0 - rho1 * TANH( ( X2 + Half ) / a )

        END IF

        uAF(iNX,iAF_P ) = One
        uPF(iNX,iPF_V3) = Zero
        uPF(iNX,iPF_E)  = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
        uPF(iNX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

        epsMin_Euler_GR &
          = MIN( epsMin_Euler_GR, uPF(iNX,iPF_E) / uPF(iNX,iPF_D) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_parallel_reduce_min( epsMin_Euler_GR )

    epsMin_Euler_GR = 1.0e-12_DP * epsMin_Euler_GR

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_KelvinHelmholtz2D


  SUBROUTINE InitializeFields_Advection3D( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: uAF(nDOFX,nAF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2, iNX3
    REAL(DP) :: X1, X2, X3

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    REAL(DP) :: D_0, V1, V2, V3, P
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
      uCF => MF_uCF % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          uAF(iNX,iAF_P ) = P

          uPF(iNX,iPF_D) &
            = D_0 + EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) ) &
                  * EXP( -( X2 - X2_0 )**2 / ( Two * sigmaX2**2 ) ) &
                  * EXP( -( X3 - X3_0 )**2 / ( Two * sigmaX3**2 ) )

          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = uAF(iNX,iAF_P) / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection3D


END MODULE MF_InitializationModule
