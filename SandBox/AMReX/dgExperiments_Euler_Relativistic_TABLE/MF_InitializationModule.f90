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
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS
  USE UnitsModule, ONLY: &
    Gram, &
    Erg, &
    Centimeter, &
    Kilometer, &
    Millisecond
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_E, &
    nAF
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
    t_end, &
    dt_wrt, &
    xL, &
    xR, &
    nX, &
    swX, &
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

  REAL(DP), PARAMETER :: LengthUnits  = Kilometer
  REAL(DP), PARAMETER :: DensityUnits = Gram / Centimeter**3

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_parmparse) :: PP

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'Advection1D' )

        CALL InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCF )

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 201, Message_Option &
                        = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'dt_wrt', dt_wrt )
    CALL amrex_parmparse_destroy( PP )

    dt_wrt = dt_wrt * t_end

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields_Advection1D( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
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
    REAL(DP) :: L, X1_0
    REAL(DP) :: sigmaX1

    AdvectionProfile = 'SineWave'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
      CALL PP % get  ( 't_end'           , t_end            )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( AdvectionProfile ) )

      CASE( 'SineWave' )

        D_0 = 1.0e12_DP * DensityUnits
        Amp = 0.1_DP * D_0

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0e-2_DP * D_0

        L = xR(1) - xL(1)

        t_end &
          = t_end * L / V1

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )
          WRITE(*,*)
          WRITE(*,'(8x,A,ES14.6E3,A)') '  D_0: ', D_0 / DensityUnits, ' gram/cm^3'
          WRITE(*,'(8x,A,ES14.6E3,A)') '  Amp: ', Amp / DensityUnits, ' gram/cm^3'
          WRITE(*,'(8x,A,ES14.6E3)')   '   V1: ', V1
          WRITE(*,'(8x,A,ES14.6E3)')   '   V2: ', V2
          WRITE(*,'(8x,A,ES14.6E3)')   '   V3: ', V3
          WRITE(*,'(8x,A,ES14.6E3,A)') '    P: ', P / D_0, ' D_0'
          WRITE(*,'(8x,A,ES14.6E3,A)') 't_end: ', t_end / Millisecond, ' ms'
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
      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( TRIM( AdvectionProfile ) .EQ. 'SineWave' )THEN

          uPF(iNX,iPF_D ) = D_0 + Amp * SIN( TwoPi * X1 / L )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_Ne) = Zero

          uAF(iNX,iAF_P ) = P
          uAF(iNX,iAF_Ye) = 0.3_DP

        END IF

        CALL ComputeTemperatureFromPressure &
               ( uPF(iNX,iPF_D ), &
                 uAF(iNX,iAF_P ), &
                 uAF(iNX,iAF_Ye), &
                 uAF(iNX,iAF_T ) )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNX,iPF_D ), &
                 uAF(iNX,iAF_T ), &
                 uAF(iNX,iAF_Ye), &
                 uPF(iNX,iPF_E ), &
                 uAF(iNX,iAF_E ), &
                 uPF(iNX,iPF_Ne) )

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


END MODULE MF_InitializationModule
