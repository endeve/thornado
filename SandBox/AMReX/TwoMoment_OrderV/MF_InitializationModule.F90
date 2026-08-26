MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_error_module, ONLY: &
   amrex_abort

  ! --- thornado Modules ---

  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    nX, &
    nE, &
    swX, &
    nNodesX, &
    nNodesE, &
    nDOFE, &
    nDOFZ, &
    iZ_B0, &
    iZ_E0, &
    iZ_B1, &
    iZ_E1, &
    iE_B0, &
    iE_E0
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    uPR, &
    uCR, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nCR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    nPR, &
    nCR, &
    nSpecies
  USE FluidFieldsModule, ONLY: &
    uPF, &
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
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Phi_N, &
    iGF_Psi, &
    iGF_SqrtGm
  USE TwoMoment_OpacityModule, ONLY: &
   uOP, &
   iOP_Sigma, &
   iOP_Chi, &
   iOP_D0
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    NodesX1, &
    NodesX2
  USE MeshModule, ONLY: &
    MeshType, &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    SolarMass, &
    SpeedOfLight, &
    Erg, &
    Gram, &
    Centimeter, &
    MeV, &
    BoltzmannConstant, &
    GravitationalConstant
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Pi, &
    Half, &
    TwoPi, &
    Three, &
    SqrtTiny
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    xR,        &
    xL,        &
    nNodes
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_base_module, ONLY: &
  amrex_problo

  USE amrex_amrcore_module, ONLY: amrex_geom
  USE MF_MeshModule       , ONLY: CreateMesh_MF, DestroyMesh_MF
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    uCR_iBC
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    uCF_iBC


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR, MF_uCF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'ALT_SineWaveStreaming' )

        CALL InitializeFields_ALT_SineWaveStreaming &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'GaussianDiffusion' )

        CALL InitializeFields_GaussianDiffusion &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'ALT_GaussianDiffusion' )

        CALL InitializeFields_ALT_GaussianDiffusion &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'ALT_GaussianDiffusion2D' )
        CALL InitializeFields_ALT_GaussianDiffusion2D &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'ShadowCasting2D_Cartesian' )
        CALL InitializeFields_ShadowCasting2D &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'TransparentShock' )

        CALL InitializeFields_TransparentShock &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'StreamingDopplerShift' )

        CALL InitializeFields_StreamingDopplerShift &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'TransparentVortex' )

        CALL InitializeFields_ALT_TransparentVortex &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'TransparentVortex_Spherical' )

        CALL InitializeFields_ALT_TransparentVortex_Spherical &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

     CASE( 'ALT_RadiatingSphere' )

       CALL InitializeFields_ALT_RadiatingSphere &
              ( iLevel, MF_uGF, MF_uCR, MF_uCF )
    
    CASE( 'ExpandingAtmosphere' )

        CALL InitializeFields_ExpandingAtmosphere &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 301, Message_Option &
                        = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


SUBROUTINE InitializeFields_ALT_SineWaveStreaming &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER        :: iDim, iE
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3, X_2D, L
    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iZ_B(4), iZ_E(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    CHARACTER(:), ALLOCATABLE :: Direction
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    Direction = 'X'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query ( 'Direction', &
                         Direction )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )
      
      loX = BX % lo
      hiX = BX % hi
      
      iZ_B = [iZ_B0(1), loX(1), loX(2), loX(3)]
      iZ_E = [iZ_E0(1), hiX(1), hiX(2), hiX(3)]

      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_B0(1):iZ_E0(1), &
                      loX(1):hiX(1), &
                      loX(2):hiX(2), &
                      loX(3):hiX(3), &
                      nCR, nSpecies) )
      
      uCR_K = Zero

      
      CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), loX(1), loX(2), loX(3)], &  
                            [iZ_E0(1), hiX(1), hiX(2), hiX(3)], &
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = V_0(1)
          uPF_K(iNodeX,iPF_V2) = V_0(2)
          uPF_K(iNodeX,iPF_V3) = V_0(3)
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
          iNodeE = MOD( (iNodeZ-1), nDOFE ) + 1

          iNodeZ2 = NodeNumberTable(2,iNodeZ)
          iNodeZ3 = NodeNumberTable(3,iNodeZ)
          iNodeZ4 = NodeNumberTable(4,iNodeZ)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )
          X3 = NodeCoordinate( MeshX(3), iX3, iNodeZ4 )

          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)  ! Loop over energy bounds

            SELECT CASE( TRIM( Direction ) )
              CASE( 'X' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Y' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X2 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Z' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X3 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

              CASE( 'XY' )

                X_2D = SQRT( 2.0_DP ) * X1 &
                         + SQRT( 2.0_DP ) * X2

                L  = SQRT( 2.0_DP )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X_2D / L )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A8,A)')    &
                  '', 'InitializeFields_SineWaveStreaming'
                WRITE(*,'(A8,A,A2)') &
                  '', 'Invalid Direction: ', TRIM( Direction )
                WRITE(*,*)
                STOP

            END SELECT

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO 
          END DO 

        END DO 

      END DO ! iX1 = loX(1), hiX(1)
      END DO ! iX2 = loX(2), hiX(2)
      END DO ! iX3 = loX(3), hiX(3)
      
      ! HERE
      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), loX(1), loX(2), loX(3)], & 
                            [iZ_E0(1), hiX(1), hiX(2), hiX(3)], & 
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DEALLOCATE( uCR_K )

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ALT_SineWaveStreaming


  SUBROUTINE InitializeFields_SineWaveStreaming &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER        :: iDim, iE
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3, X_2D, L
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    CHARACTER(:), ALLOCATABLE :: Direction
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    Direction = 'X'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query ( 'Direction', &
                         Direction )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = V_0(1)
          uPF_K(iNodeX,iPF_V2) = V_0(2)
          uPF_K(iNodeX,iPF_V3) = V_0(3)
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNodeZ4 )

            SELECT CASE( TRIM( Direction ) )
              CASE( 'X' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Y' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X2 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Z' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X3 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )


              CASE( 'XY' )

                X_2D = SQRT( 2.0_DP ) * X1 &
                         + SQRT( 2.0_DP ) * X2

                L  = SQRT( 2.0_DP )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  =0.50_DP + 0.49_DP * SIN( TwoPi * X_2D / L )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A8,A)')    &
                  '', 'InitializeFields_SineWaveStreaming'
                WRITE(*,'(A8,A,A2)') &
                  '', 'Invalid Direction: ', TRIM( Direction )
                WRITE(*,*)
                STOP

          END SELECT

          CALL ComputeConserved_TwoMoment &
                 ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                   uPF_K(iNodeX,iPF_V1), &
                   uPF_K(iNodeX,iPF_V2), &
                   uPF_K(iNodeX,iPF_V3), &
                   uGF_K(iNodeX,iGF_Gm_dd_11), &
                   uGF_K(iNodeX,iGF_Gm_dd_22), &
                   uGF_K(iNodeX,iGF_Gm_dd_33) )


          END DO 
          END DO 

        END DO 

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_SineWaveStreaming


SUBROUTINE InitializeFields_ALT_GaussianDiffusion &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies, loX(3), hiX(3)
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3, t_0, D_min, D_0, X1_0, X2_0
    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)
    REAL(DP)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(DP)       :: uGF_K( nDOFX, nGF )
    REAL(DP)       :: uPF_K( nDOFX, nPF )
    REAL(DP)       :: uCF_K( nDOFX, nCF )
    REAL(DP)       :: uAF_K( nDOFX, nAF )

    ! --- AMReX ---
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iZ_B(4), iZ_E(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: Three, EN, Sigma
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    D_min = 1.0d-06 !Was 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One
    X2_0  = One

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get( 'Sigma', &
                         Sigma )                         
    CALL amrex_parmparse_destroy( PP )
    Three = 3.0_DP


    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )
      
      loX = BX % lo
      hiX = BX % hi
      
      iZ_B = [iZ_B0(1), loX(1), loX(2), loX(3)]
      iZ_E = [iZ_E0(1), hiX(1), hiX(2), hiX(3)]

      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_B0(1):iZ_E0(1), &
                      loX(1):hiX(1), &
                      loX(2):hiX(2), &
                      loX(3):hiX(3), &
                      nCR, nSpecies) )
      
      uCR_K = Zero

      IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
      ALLOCATE( uOP(1:nDOFZ, &
                    iZ_B0(1):iZ_E0(1), &
                    loX(1):hiX(1), &
                    loX(2):hiX(2), &
                    loX(3):hiX(3), &
                    1:3, 1:nSpecies) )

      !uOP(:,:,:,:,:,iOP_D0,   :) = Zero
      !uOP(:,:,:,:,:,iOP_Chi,  :) = Zero
      uOP(:,:,:,:,:,iOP_Sigma,:) = Sigma

      
      CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), loX(1), loX(2), loX(3)], &  
                            [iZ_E0(1), hiX(1), hiX(2), hiX(3)], &
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = V_0(1)
          uPF_K(iNodeX,iPF_V2) = V_0(2)
          uPF_K(iNodeX,iPF_V3) = V_0(3)
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
          iNodeE = MOD( (iNodeZ-1), nDOFE ) + 1

          iNodeZ2 = NodeNumberTable(2,iNodeZ)
          iNodeZ3 = NodeNumberTable(3,iNodeZ)
          iNodeZ4 = NodeNumberTable(4,iNodeZ)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )
          X3 = NodeCoordinate( MeshX(3), iX3, iNodeZ4 )

          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)  ! Loop over energy bounds


              D_0 = One / (Three * Sigma) !( Three * uOP(iNodeZ,iZ1,iX1,iX2,iX3,iOP_Sigma,iS) )

              EN = EXP( - ( (X1-X1_0)*(X1-X1_0) ) / ( 4.0_DP * t_0 * D_0 ) )

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  EN

              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = Half * (X1-X1_0) * uPR_K( iNodeZ, iZ1, iPR_D, iS )/t_0

              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP

              IF( uPR_K( iNodeZ, iZ1, iPR_D, iS ) < D_min ) THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) = D_min
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) = 0.0d-00

              END IF

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO 
          END DO 

        END DO 

      END DO ! iX1 = loX(1), hiX(1)
      END DO ! iX2 = loX(2), hiX(2)
      END DO ! iX3 = loX(3), hiX(3)
      
      ! HERE
      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), loX(1), loX(2), loX(3)], & 
                            [iZ_E0(1), hiX(1), hiX(2), hiX(3)], & 
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DEALLOCATE( uCR_K )
      IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ALT_GaussianDiffusion


  SUBROUTINE InitializeFields_ALT_GaussianDiffusion2D &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---
    INTEGER        :: iX1, iX2, iX3, iZ1, iS, iNodeZ, loX(3), hiX(3)
    INTEGER        :: iNodeX, iNodeE, iNodeZ2, iNodeZ3, iNodeZ4
    REAL(DP)       :: X1, X2, t_0, D_min, D_0, X1_0, X2_0
    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)
    REAL(DP)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(DP)       :: uGF_K( nDOFX, nGF )
    REAL(DP)       :: uPF_K( nDOFX, nPF )
    REAL(DP)       :: uCF_K( nDOFX, nCF )
    REAL(DP)       :: uAF_K( nDOFX, nAF )

    ! --- AMReX ---
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iZ_B(4), iZ_E(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: Three, EN, R2, Sigma
    TYPE(amrex_parmparse)         :: PP
    REAL(DP)    , ALLOCATABLE     :: V_0(:)

    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    D_min = 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One
    X2_0  = One
    Three = 3.0_DP

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0',   V_0   )
    CALL amrex_parmparse_destroy( PP )

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get(    'Sigma', Sigma )
    CALL amrex_parmparse_destroy( PP )

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF ) ; hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCR ) ; hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF ) ; hi_F = UBOUND( uCF )

      loX  = BX % lo
      hiX  = BX % hi

      iZ_B = [ iZ_B0(1), loX(1), loX(2), loX(3) ]
      iZ_E = [ iZ_E0(1), hiX(1), hiX(2), hiX(3) ]

      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_B0(1):iZ_E0(1), &
                      loX(1):hiX(1), &
                      loX(2):hiX(2), &
                      loX(3):hiX(3), &
                      nCR, nSpecies) )

      uCR_K = Zero

      IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
      ALLOCATE( uOP(1:nDOFZ, &
                    iZ_B0(1):iZ_E0(1), &
                    loX(1):hiX(1), &
                    loX(2):hiX(2), &
                    loX(3):hiX(3), &
                    1:3, 1:nSpecies) )

      uOP(:,:,:,:,:,iOP_Sigma,:) = Sigma

      CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                             [iZ_B0(1), loX(1), loX(2), loX(3)], &
                             [iZ_E0(1), hiX(1), hiX(2), hiX(3)], &
                             lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = V_0(1)
          uPF_K(iNodeX,iPF_V2) = V_0(2)
          uPF_K(iNodeX,iPF_V3) = V_0(3)
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        ! --- Radiation Fields (2D Gaussian in X1,X2) ---

        DO iNodeZ = 1, nDOFZ

          iNodeX  = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
          iNodeE  = MOD( (iNodeZ-1), nDOFE ) + 1

          iNodeZ2 = NodeNumberTable(2,iNodeZ)
          iNodeZ3 = NodeNumberTable(3,iNodeZ)
          iNodeZ4 = NodeNumberTable(4,iNodeZ)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )

          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)

            D_0 = One / ( Three * Sigma )

            R2  = (X1-X1_0)*(X1-X1_0) + (X2-X2_0)*(X2-X2_0)
            EN  = EXP( - R2 / ( 4.0_DP * t_0 * D_0 ) )

            uPR_K( iNodeZ, iZ1, iPR_D , iS ) = EN
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = Half * (X1-X1_0) * uPR_K( iNodeZ, iZ1, iPR_D, iS ) / t_0
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = Half * (X2-X2_0) * uPR_K( iNodeZ, iZ1, iPR_D, iS ) / t_0
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) = 0.0_DP

            IF( uPR_K( iNodeZ, iZ1, iPR_D, iS ) < D_min ) THEN
              uPR_K( iNodeZ, iZ1, iPR_D , iS ) = D_min
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) = 0.0_DP
            END IF

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO
          END DO

        END DO

      END DO
      END DO
      END DO

      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                             [iZ_B0(1), loX(1), loX(2), loX(3)], &
                             [iZ_E0(1), hiX(1), hiX(2), hiX(3)], &
                             lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DEALLOCATE( uCR_K )
      IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ALT_GaussianDiffusion2D
  


  SUBROUTINE InitializeFields_GaussianDiffusion &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeZ3, iNodeE
    REAL(DP)       :: X1, X2, X3, t_0, D_min, D_0, X1_0, X2_0
    REAL(DP)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(DP)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(DP)       :: uGF_K( nDOFX, nGF )
    REAL(DP)       :: uPF_K( nDOFX, nPF )
    REAL(DP)       :: uCF_K( nDOFX, nCF )
    REAL(DP)       :: uAF_K( nDOFX, nAF )

    ! --- AMReX ---
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: Three, EN, Sigma
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    D_min = 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One
    X2_0  = One

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get( 'Sigma', &
                         Sigma )                         
    CALL amrex_parmparse_destroy( PP )
    Three = 3.0_DP


      CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF % DataPtr( MFI )
        uCR => MF_uCR % DataPtr( MFI )
        uCF => MF_uCF % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_DP
            uPF_K(iNodeX,iPF_V1) = V_0(1)
            uPF_K(iNodeX,iPF_V2) = V_0(2)
            uPF_K(iNodeX,iPF_V3) = V_0(3)
            uPF_K(iNodeX,iPF_E ) = 0.1_DP
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)
              iNodeZ3 = NodeNumberTable(3,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

              D_0 = One / ( Three * uOP(iNodeZ,iZ1,iX1,iX2,iX3,iOP_Sigma,iS) )

              EN = EXP( - ( (X1-X1_0)*(X1-X1_0) ) / ( 4.0_DP * t_0 * D_0 ) )

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  EN

              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = Half * (X1-X1_0) * uPR_K( iNodeZ, iZ1, iPR_D, iS )/t_0

              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP

              IF( uPR_K( iNodeZ, iZ1, iPR_D, iS ) < D_min ) THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) = D_min
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0d-00

              END IF

          CALL ComputeConserved_TwoMoment &
                 ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                   uPF_K(iNodeX,iPF_V1), &
                   uPF_K(iNodeX,iPF_V2), &
                   uPF_K(iNodeX,iPF_V3), &
                   uGF_K(iNodeX,iGF_Gm_dd_11), &
                   uGF_K(iNodeX,iGF_Gm_dd_22), &
                   uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO 
          END DO 

        END DO 

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_GaussianDiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE InitializeFields_ShadowCasting2D &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeZ3, iNodeE
    REAL(DP)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(DP)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(DP)       :: uGF_K( nDOFX, nGF )
    REAL(DP)       :: uPF_K( nDOFX, nPF )
    REAL(DP)       :: uCF_K( nDOFX, nCF )
    REAL(DP)       :: uAF_K( nDOFX, nAF )

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero
    uPR_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = 0.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
               (  uPF_K(:,iPF_D ), &
                  uPF_K(:,iPF_V1), &
                  uPF_K(:,iPF_V2), &
                  uPF_K(:,iPF_V3), &
                  uPF_K(:,iPF_E ), &
                  uPF_K(:,iPF_Ne), &
                  uCF_K(:,iCF_D ), &
                  uCF_K(:,iCF_S1), &
                  uCF_K(:,iCF_S2), &
                  uCF_K(:,iCF_S3), &
                  uCF_K(:,iCF_E ), &
                  uCF_K(:,iCF_Ne), &
                  uGF_K(:,iGF_Gm_dd_11), &
                  uGF_K(:,iGF_Gm_dd_22), &
                  uGF_K(:,iGF_Gm_dd_33) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            uPR_K( iNodeZ, iZ1, iPR_D , iS ) = 1.0d-10
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO
          END DO

        END DO

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ShadowCasting2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  SUBROUTINE InitializeFields_TransparentShock &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: ShockWidth, E

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    ShockWidth = 0.01_DP

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX )

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          uPF_K(iNodeX,iPF_V1) &
            = Half * V_0(1) * ( One + TANH( (X1-1.0_DP)/ShockWidth ) )
        
        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-8
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

!!!! Applying boundary conditions

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
    
    ! --- BC for Radiation Fields ---

    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = One / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.999_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO
    
!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_TransparentShock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE InitializeFields_ALT_RadiatingSphere &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    INTEGER  :: iX1, iX2, iX3, iZ1, iS, iNodeZ
    INTEGER  :: iNodeX, iNodeX1, iNodeE
    REAL(DP) :: X1, Theta, E
    REAL(DP) :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP) :: uGF_K(nDOFX,nGF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iZ_B(4), iZ_E(4), loZ1
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(amrex_parmparse)     :: PP
    CHARACTER(:), ALLOCATABLE :: Profile
    REAL(DP)                  :: V_Max

    Profile = 'Shock'
    V_Max   = 0.20_DP
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'RadiatingSphere_Profile', Profile )
      CALL PP % query( 'RadiatingSphere_V_Max'  , V_Max   )
    CALL amrex_parmparse_destroy( PP )

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(6x,A,A)')       'Profile = ', TRIM( Profile )
      WRITE(*,'(6x,A,ES9.2E2)') 'V_Max   = ', V_Max
      WRITE(*,*)
    END IF

    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uPR_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF ); hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCR ); hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF ); hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      loZ1 = loX(1)
      IF( loX(1) .EQ. 0 ) loZ1 = loX(1) - swX(1)

      iZ_B = [ iZ_B0(1), loZ1  , loX(2), loX(3) ]
      iZ_E = [ iZ_E0(1), hiX(1), hiX(2), hiX(3) ]

      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_B0(1):iZ_E0(1), &
                      loZ1:hiX(1), &
                      loX(2):hiX(2), &
                      loX(3):hiX(3), &
                      nCR, nSpecies) )

      uCR_K = Zero

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          X1      = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          uPF_K(iNodeX,iPF_D ) = One

          SELECT CASE( TRIM( Profile ) )

          CASE( 'Shock' )

            IF( X1 .LE. 135.0_DP )THEN
              uPF_K(iNodeX,iPF_V1) = Zero
            ELSE IF( X1 .LE. 150.0_DP )THEN
              uPF_K(iNodeX,iPF_V1) = - V_Max * ( X1 - 135.0_DP ) / 15.0_DP
            ELSE
              uPF_K(iNodeX,iPF_V1) = - V_Max * ( 150.0_DP / X1 )**2
            END IF

          CASE( 'Collapse' )

            Theta = Half * ( One + TANH( ( X1 - 2.0d2 ) / 3.0d1 ) )
            uPF_K(iNodeX,iPF_V1) &
              = - ( One - Theta ) * V_Max * ( X1 - 1.0d1 ) / ( 2.0d2 - 1.0d1 ) &
                - Theta * V_Max * ( 2.0d2 / X1 )**2

          END SELECT

          uPF_K(iNodeX,iPF_V2) = Zero
          uPF_K(iNodeX,iPF_V3) = Zero
          uPF_K(iNodeX,iPF_E ) = 1.0d-1
          uPF_K(iNodeX,iPF_Ne) = Zero

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), uGF_K(:,iGF_Gm_dd_22), uGF_K(:,iGF_Gm_dd_33) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)

            uPR_K(iNodeZ,iZ1,iPR_D ,iS) = 1.0d-40
            uPR_K(iNodeZ,iZ1,iPR_I1,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I2,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I3,iS) = Zero

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO
          END DO

        END DO

      END DO
      END DO
      END DO

      IF( loX(1) .EQ. 0 )THEN

        DO iX3 = loX(3), hiX(3)
        DO iX2 = loX(2), hiX(2)
        DO iX1 = loX(1)-swX(1), loX(1)-1

          uGF_K = RESHAPE( uGF(loX(1),iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX
            uPF_K(iNodeX,iPF_D ) = One
            uPF_K(iNodeX,iPF_V1) = Zero
            uPF_K(iNodeX,iPF_V2) = Zero
            uPF_K(iNodeX,iPF_V3) = Zero
            uPF_K(iNodeX,iPF_E ) = 1.0d-1
            uPF_K(iNodeX,iPF_Ne) = Zero
          END DO

          CALL ComputeConserved_Euler_NonRelativistic &
                 ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                   uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                   uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                   uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                   uGF_K(:,iGF_Gm_dd_11), uGF_K(:,iGF_Gm_dd_22), uGF_K(:,iGF_Gm_dd_33) )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        END DO
        END DO
        END DO

      END IF

      IF( loX(1) .EQ. 0 )THEN

        DO iX3 = loX(3), hiX(3)
        DO iX2 = loX(2), hiX(2)
        DO iX1 = loX(1)-swX(1), loX(1)-1

          uGF_K = RESHAPE( uGF(loX(1),iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeZ = 1, nDOFZ

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            DO iS  = 1, nSpecies
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              E = NodeCoordinate( MeshE, iZ1, iNodeE )

              uPR_K(iNodeZ,iZ1,iPR_D ,iS) &
                = One / ( EXP( E / Three - Three ) + One )
              uPR_K(iNodeZ,iZ1,iPR_I1,iS) &
                = 0.999_DP * uPR_K(iNodeZ,iZ1,iPR_D,iS)
              uPR_K(iNodeZ,iZ1,iPR_I2,iS) = Zero
              uPR_K(iNodeZ,iZ1,iPR_I3,iS) = Zero

              CALL ComputeConserved_TwoMoment &
                     ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                       uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                       uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                       uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                       uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                       uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                       uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                       uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                       uPF_K(iNodeX,iPF_V1), uPF_K(iNodeX,iPF_V2), &
                       uPF_K(iNodeX,iPF_V3), &
                       uGF_K(iNodeX,iGF_Gm_dd_11), &
                       uGF_K(iNodeX,iGF_Gm_dd_22), &
                       uGF_K(iNodeX,iGF_Gm_dd_33) )

            END DO
            END DO

          END DO

        END DO
        END DO
        END DO

      END IF

      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                             iZ_B, iZ_E, lo_C, iZ_B, iZ_E, uCR, uCR_K )

      DEALLOCATE( uCR_K )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ALT_RadiatingSphere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE InitializeFields_StreamingDopplerShift &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---

    REAL(DP), PARAMETER :: X_0 = 2.0_DP
    REAL(DP), PARAMETER :: X_1 = 3.5_DP
    REAL(DP), PARAMETER :: X_2 = 6.5_DP
    REAL(DP), PARAMETER :: X_3 = 8.0_DP
    REAL(DP), PARAMETER :: L_X = 6.0_DP

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    REAL(DP)       :: E

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )


    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX )

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          IF( X1 .LT. X_0 )THEN
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
          ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
            uPF_K(iNodeX,iPF_V1) = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
            uPF_K(iNodeX,iPF_V1) = V_0(1)
          ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
            uPF_K(iNodeX,iPF_V1)= V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSE
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
          END IF
        
        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
        END DO

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-40
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

!!!! Applying boundary conditions

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
    
    ! --- BC for Radiation Fields ---

    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = One / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.999_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO
    
!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_StreamingDopplerShift



SUBROUTINE InitializeFields_TransparentVortex &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )
    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3

    REAL(DP)       :: uCR_K(nDOFZ,iZ_B0(1):iZ_E0(1),1:nX(1),1:nX(2),1:nX(3),nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    REAL(DP)       :: Beta, R
    REAL(DP)       :: E, Mu_0
    ! --- AMReX ---
    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iZ_B(4), iZ_E(4) 
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    ! --- Problem-dependent parameters ---
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )
    Beta = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )
    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero
    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )
    DO WHILE( MFI % next() )
      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      BX = MFI % tilebox()
      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )
      loX = BX % lo
      hiX = BX % hi
      
      iZ_B = [iZ_B0(1), loX(1), loX(2), loX(3)]
      iZ_E = [iZ_E0(1), hiX(1), hiX(2), hiX(3)]
      CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), 1, 1, 1], &
                            [iZ_E0(1), nX(1), nX(2), nX(3)], &
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )
      
      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)
        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )
        DO iNodeX = 1, nDOFX
          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)
          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
          R  = SQRT( X1**2 + X2**2 )
          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) &
            = - X2 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V2) &
            = + X1 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP
        END DO
        
        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))
        
        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
              = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
        
        DO iNodeZ = 1, nDOFZ
          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)
            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1
            E = NodeCoordinate( MeshE, iZ1, iNodeE )
            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-8
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP
            
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies
        END DO ! iNodeZ = 1, nDOFZ
        
      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)
      
      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                            [iZ_B0(1), 1, 1, 1], &
                            [iZ_E0(1), nX(1), nX(2), nX(3)], &
                            lo_C, iZ_B, iZ_E, uCR, uCR_K )
      
!!!! Applying boundary conditions

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1
      DO iNodeX = 1, nDOFX
        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero
        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))
      END DO
      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
    END DO
    END DO
    END DO

!!!!Looping in y direction

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2)-sWX(2), loX(2)-1
    DO iX1 = loX(1), hiX(1)
      DO iNodeX = 1, nDOFX
        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero
        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))
      END DO
      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
    END DO
    END DO
    END DO
    
    
    ! --- BC for Radiation Fields ---
!! Looping in X direction
    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)
    
    iZ_B = [loZ(1), loZ(2)-sWX(1), loZ(3), loZ(4)]
    iZ_E = [hiZ(1), loZ(2)-1, hiZ(3), hiZ(4)]
    CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                          [iZ_B0(1), 1, 1, 1], &
                          [iZ_E0(1), nX(1), nX(2), nX(3)], &
                          lo_C, iZ_B, iZ_E, uCR, uCR_K )
    
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1
      DO iNodeZ = 1, nDOFZ
        DO iS  = 1, nSpecies
          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1
            E = NodeCoordinate( MeshE, iZ1, iNodeE )
            Mu_0 = 0.9_DP
            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.5_DP * ( One + Mu_0 ) * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO
      
    END DO
    END DO
    END DO
    
    CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                          [iZ_B0(1), 1, 1, 1], &
                          [iZ_E0(1), nX(1), nX(2), nX(3)], &
                          lo_C, iZ_B, iZ_E, uCR, uCR_K )
    
!!! Loop for y direction Radiation

    iZ_B = [loZ(1), loZ(2), loZ(3)-sWX(2), loZ(4)]
    iZ_E = [hiZ(1), hiZ(2), loZ(3)-1, hiZ(4)]
    CALL amrex2thornado_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                          [iZ_B0(1), 1, 1, 1], &
                          [iZ_E0(1), nX(1), nX(2), nX(3)], &
                          lo_C, iZ_B, iZ_E, uCR, uCR_K )
    
    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3)-sWX(2), loZ(3)-1
    DO iZ2 = loZ(2), hiZ(2)
      DO iNodeZ = 1, nDOFZ
        DO iS  = 1, nSpecies
          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1
            E = NodeCoordinate( MeshE, iZ1, iNodeE )
            Mu_0 = 0.9_DP
            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.5_DP * ( One + Mu_0 ) * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO
      
    END DO
    END DO
    END DO
    
    CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                          [iZ_B0(1), 1, 1, 1], &
                          [iZ_E0(1), nX(1), nX(2), nX(3)], &
                          lo_C, iZ_B, iZ_E, uCR, uCR_K )
!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_TransparentVortex



SUBROUTINE InitializeFields_ExpandingAtmosphere &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    REAL(DP), PARAMETER :: Rmin  = 1.0_DP
    REAL(DP), PARAMETER :: Rmax  = 11.0_DP
    REAL(DP), PARAMETER :: Dnu_0 = 1.0d-10

    INTEGER  :: iX1, iX2, iX3, iZ1, iS
    INTEGER  :: iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: R, dX1

    REAL(DP) :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP) :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP) :: uGF_K(nDOFX,nGF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uCF_K(nDOFX,nCF)

    INTEGER  :: lo_C(4), hi_C(4)
    INTEGER  :: lo_G(4), hi_G(4)
    INTEGER  :: lo_F(4), hi_F(4)
    INTEGER  :: loX(3), hiX(3), iX_B(3), iX_E(3)
    INTEGER  :: iX_Domain_E1

    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(amrex_parmparse) :: PP
    REAL(DP), ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', V_0 )
    CALL amrex_parmparse_destroy( PP )

    dX1          = ( xR(1) - xL(1) ) / ( DBLE( nX(1) ) * Two**iLevel )
    iX_Domain_E1 = nX(1) * 2**iLevel - 1

    uGF_K = Zero
    uPF_K = Zero
    uCF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF ); hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCR ); hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF ); hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      iX_B = loX
      iX_E = hiX

      IF( loX(1) .EQ. 0            ) iX_B(1) = loX(1) - swX(1)
      IF( hiX(1) .EQ. iX_Domain_E1 ) iX_E(1) = hiX(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)

        uGF_K = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)

          R = xL(1) + ( DBLE( iX1 ) + Half ) * dX1 &
                + NodesX1(iNodeX1) * dX1

          uPF_K(iNodeX,iPF_D ) = One

          IF( R .LT. Rmin )THEN

            uPF_K(iNodeX,iPF_V1) = Zero

          ELSE IF( R .LE. Rmax )THEN

            uPF_K(iNodeX,iPF_V1) = V_0(1) * ( R - Rmin ) / ( Rmax - Rmin )

          ELSE

            uPF_K(iNodeX,iPF_V1) = V_0(1)

          END IF

          uPF_K(iNodeX,iPF_V2) = Zero
          uPF_K(iNodeX,iPF_V3) = Zero
          uPF_K(iNodeX,iPF_E ) = 1.0d-1
          uPF_K(iNodeX,iPF_Ne) = Zero

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        uCR_K = Zero

        DO iNodeZ = 1, nDOFZ

          iNodeX = MOD( ( iNodeZ - 1 ) / nDOFE, nDOFX ) + 1

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            uPR_K(iNodeZ,iZ1,iPR_D ,iS) = Dnu_0
            uPR_K(iNodeZ,iZ1,iPR_I1,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I2,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I3,iS) = Zero

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO
          END DO

        END DO

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ExpandingAtmosphere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE InitializeFields_ALT_TransparentVortex_Spherical &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    REAL(DP), PARAMETER :: X_0 = 6.0_DP / SQRT( 2.0_DP )
    REAL(DP), PARAMETER :: Y_0 = 6.0_DP / SQRT( 2.0_DP )
    REAL(DP), PARAMETER :: D_1 = 1.0_DP
    REAL(DP), PARAMETER :: D_2 = 0.4_DP

    REAL(DP), PARAMETER :: D_Background = 1.0d-8

    INTEGER  :: iX1, iX2, iX3, iZ1, iS, iCR, iNodeZ, iNodeE
    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iDim
    REAL(DP) :: R, Theta, X, Y, R_0, AbsV, V_X, V_Y, Beta
    REAL(DP) :: E, D_FD, I1_FD
    REAL(DP) :: h_1, h_2, h_3

    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)
    REAL(DP) :: uPR_K(nDOFZ,iZ_B0(1):iZ_E0(1),nPR,nSpecies)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uGm_K(nDOFX,3)

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: loX(3), hiX(3)
    INTEGER                       :: iX_L(3), iX_H(3)
    INTEGER                       :: iZ_L(4), iZ_H(4)
    INTEGER                       :: iX_B0_D(3), iX_E0_D(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(MeshType) :: MeshX_L(3)

    TYPE(amrex_parmparse) :: PP
    REAL(DP), ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', V_0 )
    CALL amrex_parmparse_destroy( PP )

    Beta = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )

    iX_B0_D = amrex_geom(iLevel) % domain % lo
    iX_E0_D = amrex_geom(iLevel) % domain % hi

    CALL CreateMesh_MF( iLevel, MeshX_L )

    uGm_K = One

    IF( .NOT. ALLOCATED( uCR_iBC ) )THEN

      ALLOCATE( uCF_iBC(nDOFX,nCF) )
      ALLOCATE( uCR_iBC(nDOFZ,iZ_B0(1):iZ_E0(1),nCR,nSpecies) )

      uPF_K = Zero
      uCF_K = Zero

      uPF_K(:,iPF_D) = One
      uPF_K(:,iPF_E) = 1.0d-1

      CALL ComputeConserved_Euler_NonRelativistic                 &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
               uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
               uGm_K(:,1), uGm_K(:,1), uGm_K(:,1) )

      uCF_iBC = uCF_K

      DO iNodeZ = 1, nDOFZ

        iNodeE = MOD( iNodeZ - 1, nDOFE ) + 1

        DO iS  = 1, nSpecies
        DO iZ1 = iZ_B0(1), iZ_E0(1)

          E = NodeCoordinate( MeshE, iZ1, iNodeE )

          D_FD  = One / ( EXP( E / Three - Three ) + One )
          I1_FD = 0.5_DP * D_FD

          CALL ComputeConserved_TwoMoment          &
                 ( D_FD, I1_FD, Zero, Zero,        &
                   uCR_iBC(iNodeZ,iZ1,iCR_N ,iS),  &
                   uCR_iBC(iNodeZ,iZ1,iCR_G1,iS),  &
                   uCR_iBC(iNodeZ,iZ1,iCR_G2,iS),  &
                   uCR_iBC(iNodeZ,iZ1,iCR_G3,iS),  &
                   Zero, Zero, Zero,               &
                   One, One, One )

        END DO
        END DO

      END DO

    END IF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(6x,A)')         'TransparentVortex_Spherical'
      WRITE(*,'(6x,A,ES10.3E2)') 'v_max        = ', Beta
      WRITE(*,'(6x,A,2ES10.3E2)')'(x_0,y_0)    = ', X_0, Y_0
      WRITE(*,'(6x,A,ES10.3E2)') 'D_Background = ', D_Background
      WRITE(*,*)
    END IF

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = .FALSE. )

    DO WHILE( MFI % next() )

      BX  = MFI % tilebox()
      loX = BX % lo
      hiX = BX % hi

      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      lo_C = LBOUND( uCR ); hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF ); hi_F = UBOUND( uCF )

      iX_L = loX
      iX_H = hiX

      DO iDim = 1, 3

        IF( loX(iDim) .EQ. iX_B0_D(iDim) ) &
          iX_L(iDim) = loX(iDim) - swX(iDim)

        IF( hiX(iDim) .EQ. iX_E0_D(iDim) ) &
          iX_H(iDim) = hiX(iDim) + swX(iDim)

      END DO

      iZ_L = [ iZ_B0(1), iX_L(1), iX_L(2), iX_L(3) ]
      iZ_H = [ iZ_E0(1), iX_H(1), iX_H(2), iX_H(3) ]

      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_L(1):iZ_H(1), iZ_L(2):iZ_H(2), &
                      iZ_L(3):iZ_H(3), iZ_L(4):iZ_H(4), &
                      nCR, nSpecies) )

      uCR_K = Zero

      DO iX3 = iX_L(3), iX_H(3)
      DO iX2 = iX_L(2), iX_H(2)
      DO iX1 = iX_L(1), iX_H(1)

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          R     = NodeCoordinate( MeshX_L(1), iX1, iNodeX1 )
          Theta = NodeCoordinate( MeshX_L(2), iX2, iNodeX2 )

          h_1 = One
          h_2 = R
          h_3 = R * SIN( Theta )

          uGm_K(iNodeX,1) = h_1**2
          uGm_K(iNodeX,2) = h_2**2
          uGm_K(iNodeX,3) = MAX( h_3**2, SqrtTiny )

          X = R * SIN( Theta )
          Y = R * COS( Theta )

          R_0 = SQRT( ( X - X_0 )**2 + ( Y - Y_0 )**2 )

          IF( R_0 .LT. SqrtTiny )THEN

            V_X = Zero
            V_Y = Zero

          ELSE

            AbsV = Beta * EXP( - ( R_0 - D_1 )**2 / D_2**2 )

            V_X = - AbsV * ( Y - Y_0 ) / R_0
            V_Y = + AbsV * ( X - X_0 ) / R_0

          END IF

          uPF_K(iNodeX,iPF_D ) = One
          uPF_K(iNodeX,iPF_V1) &
            = ( SIN( Theta ) * V_X + COS( Theta ) * V_Y ) / h_1
          uPF_K(iNodeX,iPF_V2) &
            = ( COS( Theta ) * V_X - SIN( Theta ) * V_Y ) / h_2
          uPF_K(iNodeX,iPF_V3) = Zero
          uPF_K(iNodeX,iPF_E ) = 1.0d-1
          uPF_K(iNodeX,iPF_Ne) = Zero

        END DO

        CALL ComputeConserved_Euler_NonRelativistic                 &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGm_K(:,1), uGm_K(:,2), uGm_K(:,3) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        uPR_K = Zero

        DO iNodeZ = 1, nDOFZ

          iNodeX = MOD( ( iNodeZ - 1 ) / nDOFE, nDOFX ) + 1

          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)

            uPR_K(iNodeZ,iZ1,iPR_D ,iS) = D_Background
            uPR_K(iNodeZ,iZ1,iPR_I1,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I2,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I3,iS) = Zero

            CALL ComputeConserved_TwoMoment                                    &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS),                  &
                     uPF_K(iNodeX,iPF_V1), uPF_K(iNodeX,iPF_V2),               &
                     uPF_K(iNodeX,iPF_V3),                                     &
                     uGm_K(iNodeX,1), uGm_K(iNodeX,2), uGm_K(iNodeX,3) )

          END DO
          END DO

        END DO

      END DO
      END DO
      END DO

      IF( loX(1) .EQ. iX_B0_D(1) )THEN

        DO iX3 = iX_L(3), iX_H(3)
        DO iX2 = iX_L(2), iX_H(2)
        DO iX1 = loX(1)-swX(1), loX(1)-1

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_iBC, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iS  = 1, nSpecies
          DO iCR = 1, nCR

            uCR_K(:,:,iX1,iX2,iX3,iCR,iS) = uCR_iBC(:,:,iCR,iS)

          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                             iZ_L, iZ_H, lo_C, iZ_L, iZ_H, uCR, uCR_K )

      DEALLOCATE( uCR_K )

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX_L )

    DEALLOCATE( V_0 )

  END SUBROUTINE InitializeFields_ALT_TransparentVortex_Spherical


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE InitializeFields_ALT_TransparentVortex &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )
 
    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
 
    REAL(DP), PARAMETER :: Mu_0         = 0.9_DP
    REAL(DP), PARAMETER :: D_Background = 1.0d-8
 
    INTEGER  :: iX1, iX2, iX3, iZ1, iS, iCR, iNodeZ, iNodeE
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, R, Beta, E, D_Beam
 
    REAL(DP), ALLOCATABLE :: uCR_K(:,:,:,:,:,:,:)
    REAL(DP) :: uPR_K(nDOFZ,iZ_B0(1):iZ_E0(1),nPR,nSpecies)
    REAL(DP) :: uGF_K(nDOFX,nGF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uGm_K(nDOFX)
 
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: loX(3), hiX(3)
    INTEGER                       :: gX_B(3), gX_E(3)
    INTEGER                       :: iZ_B(4), iZ_E(4), iZ_L(4), iZ_H(4)
    INTEGER                       :: iX_B0_D(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
 
    TYPE(amrex_parmparse) :: PP
    REAL(DP), ALLOCATABLE :: V_0(:)
 
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', V_0 )
    CALL amrex_parmparse_destroy( PP )
 
    Beta = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )
 
    iX_B0_D = amrex_geom(iLevel) % domain % lo
  
    IF( .NOT. ALLOCATED( uCR_iBC ) )THEN
 
      ALLOCATE( uCF_iBC(nDOFX,nCF) )
      ALLOCATE( uCR_iBC(nDOFZ,iZ_B0(1):iZ_E0(1),nCR,nSpecies) )
 
      uPF_K = Zero
      uCF_K = Zero
      uGm_K = One
 
      uPF_K(:,iPF_D ) = One
      uPF_K(:,iPF_E ) = 1.0d-1
 
      CALL ComputeConserved_Euler_NonRelativistic                 &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
               uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
               uGm_K, uGm_K, uGm_K )
 
      uCF_iBC = uCF_K
 
      DO iNodeZ = 1, nDOFZ
 
        iNodeE = MOD( iNodeZ - 1, nDOFE ) + 1
 
        DO iS  = 1, nSpecies
        DO iZ1 = iZ_B0(1), iZ_E0(1)
 
          E = NodeCoordinate( MeshE, iZ1, iNodeE )
 
          D_Beam = Half * ( One - Mu_0 ) &
                     / ( EXP( E / Three - Three ) + One )
 
          CALL ComputeConserved_TwoMoment                       &
                 ( D_Beam,                                      &
                   Half * ( One + Mu_0 ) * D_Beam,              &
                   Zero, Zero,                                  &
                   uCR_iBC(iNodeZ,iZ1,iCR_N ,iS),               &
                   uCR_iBC(iNodeZ,iZ1,iCR_G1,iS),               &
                   uCR_iBC(iNodeZ,iZ1,iCR_G2,iS),               &
                   uCR_iBC(iNodeZ,iZ1,iCR_G3,iS),               &
                   Zero, Zero, Zero,                            &
                   One, One, One )
 
        END DO
        END DO
 
      END DO
 
    END IF
 
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uPR_K = Zero
 
    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = .FALSE. )
 
    DO WHILE( MFI % next() )
 
      BX  = MFI % tilebox()
      loX = BX % lo
      hiX = BX % hi
 
      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
 
      lo_G = LBOUND( uGF ); hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCR ); hi_C = UBOUND( uCR )
      lo_F = LBOUND( uCF ); hi_F = UBOUND( uCF )
  
      iZ_L = [ iZ_B0(1), loX(1)-swX(1), loX(2)-swX(2), loX(3)-swX(3) ]
      iZ_H = [ iZ_E0(1), hiX(1)+swX(1), hiX(2)+swX(2), hiX(3)+swX(3) ]
 
      ALLOCATE( uCR_K(nDOFZ, &
                      iZ_L(1):iZ_H(1), iZ_L(2):iZ_H(2), &
                      iZ_L(3):iZ_H(3), iZ_L(4):iZ_H(4), &
                      nCR, nSpecies) )
 
      uCR_K = Zero
 
      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)
 
        uGF_K = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )
 
        DO iNodeX = 1, nDOFX
 
          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)
 
          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
 
          R = SQRT( X1**2 + X2**2 )
 
          uPF_K(iNodeX,iPF_D ) = One
          uPF_K(iNodeX,iPF_V1) = - X2 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V2) = + X1 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V3) = Zero
          uPF_K(iNodeX,iPF_E ) = 1.0d-1
          uPF_K(iNodeX,iPF_Ne) = Zero
 
        END DO
 
        CALL ComputeConserved_Euler_NonRelativistic                 &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), uGF_K(:,iGF_Gm_dd_22),      &
                 uGF_K(:,iGF_Gm_dd_33) )
 
        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
 
        uPR_K = Zero
 
        DO iNodeZ = 1, nDOFZ
 
          iNodeX = MOD( ( iNodeZ - 1 ) / nDOFE, nDOFX ) + 1
 
          DO iS  = 1, nSpecies
          DO iZ1 = iZ_B0(1), iZ_E0(1)
 
            uPR_K(iNodeZ,iZ1,iPR_D ,iS) = D_Background
            uPR_K(iNodeZ,iZ1,iPR_I1,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I2,iS) = Zero
            uPR_K(iNodeZ,iZ1,iPR_I3,iS) = Zero
 
            CALL ComputeConserved_TwoMoment                                    &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N ,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G1,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G2,iS),                  &
                     uCR_K(iNodeZ,iZ1,iX1,iX2,iX3,iCR_G3,iS),                  &
                     uPF_K(iNodeX,iPF_V1), uPF_K(iNodeX,iPF_V2),               &
                     uPF_K(iNodeX,iPF_V3),                                     &
                     uGF_K(iNodeX,iGF_Gm_dd_11), uGF_K(iNodeX,iGF_Gm_dd_22),   &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
 
          END DO
          END DO
 
        END DO
 
      END DO
      END DO
      END DO
 
      iZ_B = [ iZ_B0(1), loX(1), loX(2), loX(3) ]
      iZ_E = [ iZ_E0(1), hiX(1), hiX(2), hiX(3) ]
 
      CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                             iZ_L, iZ_H, lo_C, iZ_B, iZ_E, uCR, uCR_K )
  
      IF( loX(1) .EQ. iX_B0_D(1) )THEN
 
        gX_B = [ loX(1)-swX(1), loX(2)-swX(2), loX(3)-swX(3) ]
        gX_E = [ loX(1)-1     , hiX(2)+swX(2), hiX(3)+swX(3) ]
 
        DO iX3 = gX_B(3), gX_E(3)
        DO iX2 = gX_B(2), gX_E(2)
        DO iX1 = gX_B(1), gX_E(1)
 
          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_iBC, [ hi_F(4) - lo_F(4) + 1 ] )
 
          DO iS  = 1, nSpecies
          DO iCR = 1, nCR
 
            uCR_K(:,:,iX1,iX2,iX3,iCR,iS) = uCR_iBC(:,:,iCR,iS)
 
          END DO
          END DO
 
        END DO
        END DO
        END DO
 
        iZ_B = [ iZ_B0(1), gX_B(1), gX_B(2), gX_B(3) ]
        iZ_E = [ iZ_E0(1), gX_E(1), gX_E(2), gX_E(3) ]
 
        CALL thornado2amrex_Z( nCR, nSpecies, nE, iZ_B0(1), iZ_E0(1), &
                               iZ_L, iZ_H, lo_C, iZ_B, iZ_E, uCR, uCR_K )
 
      END IF
 
      DEALLOCATE( uCR_K )
 
    END DO
 
    CALL amrex_mfiter_destroy( MFI )
 
    DEALLOCATE( V_0 )
 
  END SUBROUTINE InitializeFields_ALT_TransparentVortex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE MF_InitializationModule
