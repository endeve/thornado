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
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    nGF
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
    nPF
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

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    CALL InitializeFields( iLevel, MF_uGF, MF_uCF )

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1

    CHARACTER(:), ALLOCATABLE :: Profile

    REAL(DP) :: D_0, V1, V2, V3, P, Ne
    REAL(DP) :: Amp
    REAL(DP) :: X1_0
    REAL(DP) :: sigmaX1

    Profile = 'SineWaveX1'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'Profile', Profile )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( Profile ) )

      CASE( 'SineWaveX1' )

        D_0 = 1.0_DP
        Amp = 0.5_DP

        V1 = 0.1_DP
        V2 = 0.0_DP
        V3 = 0.0_DP
        P  = 1.0_DP
        Ne = 0.0_DP

        IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(6x,A,A)') 'Profile: ', TRIM( Profile )
          WRITE(*,'(6x,A,A)') '-------- '
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

          WRITE(*,'(6x,A,A)') 'Profile: ', TRIM( Profile )
          WRITE(*,'(6x,A,A)') '-------- '
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
                        = 'Invalid Profile: ' &
                            // TRIM( Profile ) )

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

        IF( TRIM( Profile ) .EQ. 'SineWaveX1' )THEN

          uPF(iNX,iPF_D ) = D_0 + Amp * SIN( TwoPi * X1 )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = P / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Ne

        ELSE IF( TRIM( Profile ) .EQ. 'Gaussian' )THEN

          uPF(iNX,iPF_D) &
            = D_0 + EXP( -( X1 - X1_0 )**2 / ( Two * sigmaX1**2 ) )
          uPF(iNX,iPF_V1) = V1
          uPF(iNX,iPF_V2) = V2
          uPF(iNX,iPF_V3) = V3
          uPF(iNX,iPF_E ) = P / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Ne

        END IF

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(iNX,iPF_D ), uPF(iNX,iPF_E), &
                 uPF(iNX,iPF_Ne), Pressure )

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
                 Pressure )

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

  END SUBROUTINE InitializeFields


END MODULE MF_InitializationModule
