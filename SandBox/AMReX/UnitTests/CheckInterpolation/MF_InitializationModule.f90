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
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
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
    ( iLevel, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    CALL InitializeFields( iLevel, MF_uCF )

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields( iLevel, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1

    REAL(DP) :: D_0, S1, S2, S3, tau
    REAL(DP) :: Amp

    D_0 = One
    Amp = Half

    S1  = Zero
    S2  = Zero
    S3  = Zero
    tau = One

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A)') 'Fluid:'
      WRITE(*,'(6x,A)') '------'
      WRITE(*,*)
      WRITE(*,'(8x,A,F5.3)') 'D_0: ', D_0
      WRITE(*,'(8x,A,F5.3)') 'Amp: ', Amp
      WRITE(*,'(8x,A,F5.3)') ' S1: ', S1
      WRITE(*,'(8x,A,F5.3)') ' S2: ', S2
      WRITE(*,'(8x,A,F5.3)') ' S3: ', S3
      WRITE(*,'(8x,A,F5.3)') 'tau: ', tau
      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uCF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        U(iNX,iX1,iX2,iX3,iCF_D ) = D_0 + Amp * SIN( TwoPi * X1 )
        U(iNX,iX1,iX2,iX3,iCF_S1) = S1
        U(iNX,iX1,iX2,iX3,iCF_S2) = S2
        U(iNX,iX1,iX2,iX3,iCF_S3) = S3
        U(iNX,iX1,iX2,iX3,iCF_E ) = tau
        U(iNX,iX1,iX2,iX3,iCF_Ne) = Zero

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

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields


END MODULE MF_InitializationModule
