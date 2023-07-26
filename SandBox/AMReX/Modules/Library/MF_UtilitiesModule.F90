!> Module for operations on MultiFabs
MODULE MF_UtilitiesModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_multifab, &
    amrex_imultifab
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_myproc
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ, &
    nNodesX, &
    nNodesE, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate
  USE UnitsModule, ONLY: &
    MeV
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
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
    iPF_D, &
    iPF_E, &
    iPF_Ne, &
    unitsPF, &
    nAF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Cs, &
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh, &
    unitsAF
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nPR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeAuxiliary_Fluid_TABLE, &
    ApplyEquationOfState_TABLE
#ifndef THORNADO_NOTRANSPORT
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment
#endif

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nX, &
    nE, &
    swX, &
    swE, &
    nSpecies, &
    StepNo
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
#ifndef THORNADO_NOTRANSPORT
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
#endif
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_Allocate_X, &
    Timer_AMReX_Allocate_Z, &
    Timer_AMReX_PermuteData_X
  USE MF_EdgeMapModule, ONLY: &
    ConstructEdgeMap, &
    EdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: ShowVariableFromMultiFab_Single
  PUBLIC :: MultiplyWithMetric
  PUBLIC :: amrex2thornado_X
  PUBLIC :: thornado2amrex_X
  PUBLIC :: amrex2thornado_Z
  PUBLIC :: thornado2amrex_Z
  PUBLIC :: amrex2thornado_Integrated
  PUBLIC :: thornado2amrex_Integrated
  PUBLIC :: amrex2amrex_permute_Z
  PUBLIC :: amrex_permute2amrex_Z
  PUBLIC :: MF_amrex2amrex_permute_Z_Level
  PUBLIC :: MF_amrex_permute2amrex_Z_Level
  PUBLIC :: thornado2amrex_X_F
  PUBLIC :: amrex2thornado_X_F
  PUBLIC :: WriteNodalDataToFile
  PUBLIC :: WriteEulerToFile
  PUBLIC :: AllocateArray_X
  PUBLIC :: DeallocateArray_X
  PUBLIC :: AllocateArray_Z
  PUBLIC :: DeallocateArray_Z
  PUBLIC :: AllocateArray_Integrated
  PUBLIC :: DeallocateArray_Integrated
  PUBLIC :: PrintBoxArray

  INTERFACE MultiplyWithMetric
    MODULE PROCEDURE MultiplyWithMetric_uGF
    MODULE PROCEDURE MultiplyWithMetric_uCF
  END INTERFACE MultiplyWithMetric

  INTERFACE ShowVariableFromMultiFab
    MODULE PROCEDURE ShowVariableFromMultiFab_Single
    MODULE PROCEDURE ShowVariableFromMultiFab_Vector
  END INTERFACE ShowVariableFromMultiFab

  INTERFACE

    SUBROUTINE print_boxarray( BA ) BIND(c)
        IMPORT
        IMPLICIT NONE
        TYPE(c_ptr), VALUE :: BA
    END SUBROUTINE print_boxarray

  END INTERFACE

CONTAINS


  SUBROUTINE ShowVariableFromMultiFab_Single &
    ( iLevel, MF, iField, iMF_FineMask_Option, &
      swXX_Option, WriteToFile_Option, FileNameBase_Option )

    INTEGER              , INTENT(in) :: iLevel, iField
    TYPE(amrex_multifab) , INTENT(in) :: MF
    TYPE(amrex_imultifab), INTENT(in), OPTIONAL :: iMF_FineMask_Option
    INTEGER              , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL              , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)         , INTENT(in), OPTIONAL :: FileNameBase_Option

    INTEGER                       :: iX1, iX2, iX3, iNX
    INTEGER                       :: lo(4), hi(4), iX_B1(3), iX_E1(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: F       (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    INTEGER                       :: swXX(3)
    INTEGER                       :: iFileNo
    LOGICAL                       :: WriteToFile
    CHARACTER(128)                :: FMT
    CHARACTER(128)                :: FileNameBase, FileName

    TYPE(MeshType) :: MeshX(3)

    REAL(DP) :: NodesX1(nNodesX(1))
    REAL(DP) :: NodesX2(nNodesX(2))
    REAL(DP) :: NodesX3(nNodesX(3))

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    WRITE(FMT,'(A,I2.2,A,I2.2,A,I2.2,A,I3.3,A)') &
      '(I2.2,3I8.6,SP3ES25.16E3,SP', &
      nNodesX(1),    'ES25.16E3,SP', &
      nNodesX(2),    'ES25.16E3,SP', &
      nNodesX(3),    'ES25.16E3,SP', &
      nDOFX     ,    'ES25.16E3)'

    CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

    CALL CreateMesh_MF( iLevel, MeshX )

    DO WHILE( MFI % next() )

      IF( WriteToFile )THEN

        iFileNo = 100 + amrex_parallel_myproc()

        FileNameBase = 'NodalData'
        IF( PRESENT( FileNameBase_Option ) ) &
          FileNameBase = TRIM( FileNameBase_Option )

        WRITE(FileName,'(A,A,I3.3,A,I8.8,A)') &
          TRIM( FileNameBase ), '_proc', &
          amrex_parallel_myproc(), '_', StepNo(0), '.dat'

        OPEN( iFileNo, FILE = TRIM( FileName ), POSITION = 'APPEND' )

      END IF

      IF( PRESENT( iMF_FineMask_Option ) ) &
        FineMask => iMF_FineMask_Option % DataPtr( MFI )

      F => MF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo = LBOUND( F ); hi = UBOUND( F )

      iX_B1 = BX % lo - swXX
      iX_E1 = BX % hi + swXX

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        IF( PRESENT( iMF_FineMask_Option ) )THEN

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

        END IF

        DO iNX = 1, nNodesX(1)
          NodesX1(iNX) = NodeCoordinate( MeshX(1), iX1, iNX )
        END DO
        DO iNX = 1, nNodesX(2)
          NodesX2(iNX) = NodeCoordinate( MeshX(2), iX2, iNX )
        END DO
        DO iNX = 1, nNodesX(3)
          NodesX3(iNX) = NodeCoordinate( MeshX(3), iX3, iNX )
        END DO

        IF( WriteToFile )THEN

          WRITE(iFileNo,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            MeshX(1) % Width(iX1), &
            MeshX(2) % Width(iX2), &
            MeshX(3) % Width(iX3), &
            NodesX1, NodesX2, NodesX3, &
            F(iX1,iX2,iX3,1+nDOFX*(iField-1):nDOFX*iField)

        ELSE

          WRITE(*,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            F(iX1,iX2,iX3,1+nDOFX*(iField-1):nDOFX*iField)

        END IF

      END DO
      END DO
      END DO

      IF( WriteToFile )THEN

        WRITE(iFileNo,*)

        CLOSE( iFileNo )

      ELSE IF( ANY( FineMask(:,:,:,1) .EQ. 0 ) )THEN

        WRITE(*,*)

      END IF

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

    IF( .NOT. WriteToFile ) WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab_Single


  SUBROUTINE ShowVariableFromMultiFab_Vector &
    ( MF, iField, swXX_Option, WriteToFile_Option, FileNameBase_Option )

    INTEGER             , INTENT(in) :: iField
    TYPE(amrex_multifab), INTENT(in) :: MF(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL             , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)        , INTENT(in), OPTIONAL :: FileNameBase_Option

    INTEGER :: iLevel

    INTEGER        :: swXX(3)
    LOGICAL        :: WriteToFile
    CHARACTER(128) :: FileNameBase

    TYPE(amrex_imultifab) :: iMF_FineMask

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    FileNameBase = ''
    IF( PRESENT( FileNameBase_Option ) ) &
      FileNameBase = TRIM( FileNameBase_Option )

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF % BA, MF % DM )

      CALL ShowVariableFromMultiFab_Single &
             ( iLevel, MF(iLevel), iField, &
               iMF_FineMask_Option = iMF_FineMask, &
               swXX_Option = swXX, &
               WriteToFile_Option = WriteToFile, &
               FileNameBase_Option = TRIM( FileNameBase ) )

      CALL DestroyFineMask( iMF_FineMask )

    END DO

  END SUBROUTINE ShowVariableFromMultiFab_Vector


  SUBROUTINE MultiplyWithMetric_uGF &
    ( iLevel, MF_uGF, nFd, Power, swXX_Option )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    INTEGER             , INTENT(in)    :: iLevel, nFd, Power
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    INTEGER            :: iX1, iX2, iX3, iNX, iFd, swXX(3)
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(DP)           :: G_K(nDOFX,nFd)

    REAL(DP)                      :: SqrtGm(nDOFX)
    REAL(DP), CONTIGUOUS, POINTER :: G(:,:,:,:)

    swXX = swX
    IF( PRESENT( swXX_Option ) ) &
      swXX = swXX_Option

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( lo_G, hi_G, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP          BX, MFI, G_K, SqrtGm, G )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      G => MF_uGF(iLevel) % DataPtr( MFI )

      lo_G = LBOUND( G ); hi_G = UBOUND( G )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swXX
      iX_E1 = iX_E0 + swXX

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        G_K(1:nDOFX,1:nFd) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nFd ] )

        DO iNX = 1, nDOFX

          SqrtGm(iNX) = SQRT( G_K(iNX,iGF_SqrtGm) )

        END DO

        DO iFd = 1, nFd
        DO iNX = 1, nDOFX

          G_K(iNX,iFd) = G_K(iNX,iFd) * SqrtGm(iNX)**( Power )

        END DO
        END DO

      END DO
      END DO
      END DO

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

  END SUBROUTINE MultiplyWithMetric_uGF


  SUBROUTINE MultiplyWithMetric_uCF &
    ( iLevel, MF_SqrtGm, MF, nFd, Power, swXX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_SqrtGm
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    INTEGER             , INTENT(in)    :: iLevel, nFd, Power
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    INTEGER            :: iX1, iX2, iX3, iNX, iFd, swXX(3)
    INTEGER            :: lo_F(4), hi_F(4)
    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(DP)           :: F_K(nDOFX,nFd)

    REAL(DP), CONTIGUOUS, POINTER :: SqrtGm  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: F       (:,:,:,:)

    swXX = swX
    IF( PRESENT( swXX_Option ) ) &
      swXX = swXX_Option

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( lo_F, hi_F, iX_B0, iX_E0, iX_B1, iX_E1, &
    !$OMP          BX, MFI, F_K, SqrtGm, F )
#endif

    CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      SqrtGm => MF_SqrtGm  % DataPtr( MFI )
      F      => MF(iLevel) % DataPtr( MFI )

      lo_F = LBOUND( F ); hi_F = UBOUND( F )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swXX
      iX_E1 = iX_E0 + swXX

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

        F_K(1:nDOFX,1:nFd) &
          = RESHAPE( F(iX1,iX2,iX3,lo_F(4):hi_F(4)), [ nDOFX, nFd ] )

        DO iFd = 1, nFd
        DO iNX = 1, nDOFX

          F_K(iNX,iFd) = F_K(iNX,iFd) * SqrtGm(iX1,iX2,iX3,iNX)**( Power )

        END DO
        END DO

      END DO
      END DO
      END DO

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

  END SUBROUTINE MultiplyWithMetric_uCF


  SUBROUTINE amrex2thornado_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX( Timer_AMReX_PermuteData_X )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX( Timer_AMReX_PermuteData_X )

  END SUBROUTINE amrex2thornado_X


  SUBROUTINE thornado2amrex_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX( Timer_AMReX_PermuteData_X )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd) &
        = Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd)
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX( Timer_AMReX_PermuteData_X )

  END SUBROUTINE thornado2amrex_X


  SUBROUTINE amrex2thornado_Z &
    ( nFields, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iFd, iD, iNodeZ

    DO iS  = 1      , nS
    DO iFd = 1      , nFields
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iZ1    = iE_B0, iE_E0 ! always want iZ1 to not include ghost cells
      DO iNodeZ = 1    , nDOFZ

        iD = ( iS - 1 ) * nFields * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iFd - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        Data_thornado(iNodeZ,iZ1,iZ2,iZ3,iZ4,iFd,iS) &
          = Data_amrex(iZ2,iZ3,iZ4,iD)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE amrex2thornado_Z


  SUBROUTINE thornado2amrex_Z &
    ( nFields, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iFd, iNodeZ, iD

    DO iS  = 1      , nS
    DO iFd = 1      , nFields
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iZ1    = iE_B0, iE_E0
      DO iNodeZ = 1    , nDOFZ

        iD = ( iS - 1 ) * nFields * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
               + ( iFd - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
               + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        Data_amrex(iZ2,iZ3,iZ4,iD) &
          = Data_thornado(iNodeZ,iZ1,iZ2,iZ3,iZ4,iFd,iS)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE thornado2amrex_Z

  SUBROUTINE amrex2thornado_Integrated &
    ( nFields, nS, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields, nS
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iX_B(3), iX_E(3), iLo_MF(4)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iX1, iX2, iX3, iS, iFd, iD, iNodeX

    DO iS  = 1      , nS
    DO iFd = 1      , nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iNodeX = 1    , nDOFX

        iD = ( iS - 1 ) * nFields * nDOFX &
                + ( iFd - 1 ) * nDOFX &
                + iNodeX

        Data_thornado(iNodeX,iX1,iX2,iX3,iFd,iS) &
          = Data_amrex(iX1,iX2,iX3,iD)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE amrex2thornado_Integrated

  SUBROUTINE thornado2amrex_Integrated &
    ( nFields, nS, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields, nS
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iX_B(3), iX_E(3), iLo_MF(4)
    REAL(DP), INTENT(out)  :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iX1, iX2, iX3, iS, iFd, iD, iNodeX

    DO iS  = 1      , nS
    DO iFd = 1      , nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iNodeX = 1    , nDOFX

        iD = ( iS - 1 ) * nFields * nDOFX &
                + ( iFd - 1 ) * nDOFX &
                + iNodeX

        Data_amrex(iX1,iX2,iX3,iD) &
          = Data_thornado(iNodeX,iX1,iX2,iX3,iFd,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE thornado2amrex_Integrated

  SUBROUTINE thornado2amrex_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX( Timer_AMReX_PermuteData_X )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd) &
        = Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX( Timer_AMReX_PermuteData_X )

  END SUBROUTINE thornado2amrex_X_F


  SUBROUTINE amrex2thornado_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX( Timer_AMReX_PermuteData_X )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX( Timer_AMReX_PermuteData_X )

  END SUBROUTINE amrex2thornado_X_F


  SUBROUTINE WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uCR, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

#ifndef THORNADO_NOTRANSPORT

    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3),     &
                          iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4),     &
                          iZ_B (4), iZ_E (4), iE_B, iE_E,             &
                          iX_B (3), iX_E (3), iEL, iER, iLo_MF(4),    &
                          iLevel, nCompGF, nCompCF, nCompCR,          &
                          iX1, iX2, iX3, iCF, iGF, iCR, i,            &
                          iZ1, iZ2, iZ3, iZ4, iS, iE,                 &
                          iNodeZ,                                     &
                          iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeE
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(EdgeMap)      :: Edge_Map

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR (:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR (:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR (:,:,:,:,:,:,:)

    REAL(DP) :: D, I1, I2, I3, N, G1, G2, G3

    CALL AllocateArray_Z &
           ( [ 1             , &
               1 - swE       , &
               1 - swX(1)    , &
               1 - swX(2)    , &
               1 - swX(3)    , &
               1             , &
               1        ]    , &
             [ nDOFZ         , &
               nE + swE      , &
               nX(1) + swX(1), &
               nX(2) + swX(2), &
               nX(3) + swX(3), &
               nCR           , &
               nSpecies ]    , &
             CR )

    CALL AllocateArray_Z &
           ( [ 1             , &
               1 - swE       , &
               1 - swX(1)    , &
               1 - swX(2)    , &
               1 - swX(3)    , &
               1             , &
               1        ]    , &
             [ nDOFZ         , &
               nE + swE      , &
               nX(1) + swX(1), &
               nX(2) + swX(2), &
               nX(3) + swX(3), &
               nPR           , &
               nSpecies ]    , &
             PR )

    CALL AllocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nGF ], &
             G )

    CALL AllocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nCF ], &
             CF )

    CALL AllocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nPF ], &
             PF )

    G  = 0.0_DP
    CF = 0.0_DP
    PF = 0.0_DP
    CR = 0.0_DP
    PR = 0.0_DP

    ! --- Convert AMReX MultiFabs to thornado arrays ---

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF     => MF_uGF(iLevel) % DataPtr( MFI )
        nCompGF =  MF_uGF(iLevel) % nComp()

        uCF     => MF_uCF(iLevel) % DataPtr( MFI )
        nCompCF =  MF_uCF(iLevel) % nComp()

        uCR     => MF_uCR(iLevel) % DataPtr( MFI )
        nCompCR =  MF_uCr(iLevel) % nComp()

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iX_B = iX_B0
        iX_E = iX_E0

        i=1

        DO WHILE (i<=4)

          IF (i==1) THEN

            iZ_B0(i)=iE_B0
            iZ_E0(i)=iE_E0
          ELSE

            iZ_B0(i)=iX_B0(i-1)
            iZ_E0(i)=iX_E0(i-1)
          END IF
          i = i + 1
        END DO

        IF( iX_B0(1) .EQ. 1     ) iX_B(1) = 1     - swX(1)
        IF( iX_B0(2) .EQ. 1     ) iX_B(2) = 1     - swX(2)
        IF( iX_B0(3) .EQ. 1     ) iX_B(3) = 1     - swX(3)
        IF( iX_E0(1) .EQ. nX(1) ) iX_E(1) = nX(1) + swX(1)
        IF( iX_E0(2) .EQ. nX(2) ) iX_E(2) = nX(2) + swX(2)
        IF( iX_E0(3) .EQ. nX(3) ) iX_E(3) = nX(3) + swX(3)
        IF( iE_B0 .EQ. 1     ) iE_B = 1     - swE
        IF( iE_E0 .EQ. nE ) iE_E = nE + swE


        i=1

        DO WHILE (i<=4)

          IF (i==1) THEN

            iZ_B(i)=iE_B
            iZ_E(i)=iE_E
            iZ_B1(i)=iE_B1
            iZ_E1(i)=iE_E1
          ELSE

            iZ_B(i)=iX_B(i-1)
            iZ_E(i)=iX_E(i-1)
            iZ_B1(i)=iX_B1(i-1)
            iZ_E1(i)=iX_E1(i-1)
          END IF
          i = i + 1
        END DO

        CALL amrex2thornado_X( nGF, iX_B, iX_E, iLo_MF, iX_B, iX_E, uGF, &
                               G(1:nDOFX,iX_B(1):iX_E(1), &
                                         iX_B(2):iX_E(2), &
                                         iX_B(3):iX_E(3),1:nGF) )

        CALL amrex2thornado_X( nCF, iX_B, iX_E, iLo_MF, iX_B, iX_E, uCF, &
                               CF(1:nDOFX,iX_B(1):iX_E(1), &
                                          iX_B(2):iX_E(2), &
                                          iX_B(3):iX_E(3),1:nCF) )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, &
                 CR(1:nDOFZ,iZ_B1(1):iZ_E1(1), &
                            iZ_B1(2):iZ_E1(2), &
                            iZ_B1(3):iZ_E1(3), &
                            iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )



        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, CR, Edge_Map )

      END DO

    END DO

    DO iX3    = 1-swX(3), nX(3)+swX(3)
    DO iX2    = 1-swX(2), nX(2)+swX(2)
    DO iX1    = 1-swX(1), nX(1)+swX(1)
    DO iNodeX = 1       , nDOFX

      DO iGF = 1, nGF

        CALL amrex_parallel_reduce_sum( G(iNodeX,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CALL amrex_parallel_reduce_sum( CF(iNodeX,iX1,iX2,iX3,iCF) )

      END DO

    END DO
    END DO
    END DO
    END DO

    DO iS     = 1, nSpecies
    DO iZ4    = 1-swX(3), nX(3)+swX(3)
    DO iZ3    = 1-swX(2), nX(2)+swX(2)
    DO iZ2    = 1-swX(1), nX(1)+swX(1)
    DO iZ1    = 1-swE, nE+swE
    DO iNodeZ = 1       , nDOFZ


      DO iCR = 1, nCR

        CALL amrex_parallel_reduce_sum( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
      DO iX1 = 1, nX(1)

        CALL ComputePrimitive_Euler_Relativistic &
               ( CF(:,iX1,iX2,iX3,iCF_D ),       &
                 CF(:,iX1,iX2,iX3,iCF_S1),       &
                 CF(:,iX1,iX2,iX3,iCF_S2),       &
                 CF(:,iX1,iX2,iX3,iCF_S3),       &
                 CF(:,iX1,iX2,iX3,iCF_E ),       &
                 CF(:,iX1,iX2,iX3,iCF_Ne),       &
                 PF(:,iX1,iX2,iX3,iPF_D ),       &
                 PF(:,iX1,iX2,iX3,iPF_V1),       &
                 PF(:,iX1,iX2,iX3,iPF_V2),       &
                 PF(:,iX1,iX2,iX3,iPF_V3),       &
                 PF(:,iX1,iX2,iX3,iPF_E ),       &
                 PF(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )


      END DO
      END DO
      END DO

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      DO iS = 1, nSpecies
      DO iZ4 = 1, nX(3)
      DO iZ3 = 1, nX(2)
      DO iZ2 = 1, nX(1)
      DO iZ1 = 1, nE

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE
          CALL ComputePrimitive_TwoMoment &
               ( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G2, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G3, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_D, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I1, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I2, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I3, iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4, iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4, iPF_V3), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 G(iNodeX ,iZ2,iZ3,iZ4,iGF_Alpha  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_1  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_2  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_3  ) )

        END DO
        END DO

      END DO
      END DO
      END DO
      END DO
      END DO

      OPEN( UNIT = 100, FILE = TRIM( FileNameBase ) // 'r.dat'      )
      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // 'D.dat'      )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // 'I1.dat'  )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // 'I2.dat'    )
      OPEN( UNIT = 104, FILE = TRIM( FileNameBase ) // 'I3.dat' )
      OPEN( UNIT = 105, FILE = TRIM( FileNameBase ) // 'N.dat'      )
      OPEN( UNIT = 106, FILE = TRIM( FileNameBase ) // 'G1.dat'      )
      OPEN( UNIT = 107, FILE = TRIM( FileNameBase ) // 'G2.dat'      )
      OPEN( UNIT = 108, FILE = TRIM( FileNameBase ) // 'G3.dat'      )
      OPEN( UNIT = 109, FILE = TRIM( FileNameBase ) // 'D_middle.dat'      )
      OPEN( UNIT = 110, FILE = TRIM( FileNameBase ) // 'D_spatial.dat'      )


      ! --- Hacked to work only for 1D and 2D problems ---
      DO iS = 1, nSpecies
      DO iX3     = 1, nX(3)
      DO iNodeX3 = 1, nNodesX(3)

        DO iX2     = 1, nX(2)
        DO iNodeX2 = 1, nNodesX(2)

        DO iX1     = 1, nX(1)
        DO iNodeX1 = 1, nNodesX(1)
          IF     ( iNodeX1 .EQ. 1 )THEN

            iEL = 1
            iER = nNodesE

          ELSE IF( iNodeX1 .EQ. 2 )THEN

            iEL = nNodesE + 1
            iER = 2 * nNodesE

          ELSE IF( iNodeX1 .EQ. 3 )THEN

            iEL = 2 * nNodesE + 1
            iER = 3 * nNodesE

          END IF

          DO iE = 1, nE
          DO iNodeE = iEL, iER
            D = PR(iNodeE,iE,iX1,iX2,iX3,iPR_D,iS)
            I1 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I1,iS)
            I2 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I2,iS)
            I3 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I3,iS)
            N = CR(iNodeE,iE,iX1,iX2,iX3,iCR_N,iS)
            G1 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G1,iS)
            G2 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G2,iS)
            G3 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G3,iS)

            WRITE(101,'(ES24.16E3,1x)',ADVANCE='NO') &
              D
            WRITE(102,'(ES24.16E3,1x)',ADVANCE='NO') &
              I1
            WRITE(103,'(ES24.16E3,1x)',ADVANCE='NO') &
              I2
            WRITE(104,'(ES24.16E3,1x)',ADVANCE='NO') &
              I3
            WRITE(105,'(ES24.16E3,1x)',ADVANCE='NO') &
              N
            WRITE(106,'(ES24.16E3,1x)',ADVANCE='NO') &
              G1
            WRITE(107,'(ES24.16E3,1x)',ADVANCE='NO') &
              G2
            WRITE(108,'(ES24.16E3,1x)',ADVANCE='NO') &
              G3
            IF (iX1 .EQ. nX(1)/2 .AND. iNodeX1 .EQ. 1) THEN
              WRITE(109,'(ES24.16E3,1x)',ADVANCE='NO') &
                D
            END IF
            IF (iE .EQ. 1 .AND. iNodeE .EQ. 1) THEN
              WRITE(110,'(ES24.16E3,1x)',ADVANCE='NO') &
                D
            END IF
          END DO
          END DO
        WRITE(101,*)
        WRITE(102,*)
        WRITE(103,*)
        WRITE(104,*)
        WRITE(105,*)
        WRITE(106,*)
        WRITE(107,*)
        WRITE(108,*)
        END DO
        END DO
        END DO
        END DO


      END DO
      END DO
      END DO

      CLOSE( 110 )
      CLOSE( 109 )
      CLOSE( 108 )
      CLOSE( 107 )
      CLOSE( 106 )
      CLOSE( 105 )
      CLOSE( 104 )
      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )
      CLOSE( 100 )

    END IF

    CALL DeallocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nPF ], &
             PF )

    CALL DeallocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nCF ], &
             CF )

    CALL DeallocateArray_X &
           ( [ 1    , 1    -swX(1), 1    -swX(2), 1    -swX(3), 1   ], &
             [ nDOFX, nX(1)+swX(1), nX(2)+swX(2), nX(3)+swX(3), nGF ], &
             G )

    CALL DeallocateArray_Z &
           ( [ 1             , &
               1 - swE       , &
               1 - swX(1)    , &
               1 - swX(2)    , &
               1 - swX(3)    , &
               1             , &
               1        ]    , &
             [ nDOFZ         , &
               nE + swE      , &
               nX(1) + swX(1), &
               nX(2) + swX(2), &
               nX(3) + swX(3), &
               nPR           , &
               nSpecies ]    , &
             PR )

    CALL DeallocateArray_Z &
           ( [ 1             , &
               1 - swE       , &
               1 - swX(1)    , &
               1 - swX(2)    , &
               1 - swX(3)    , &
               1             , &
               1        ]    , &
             [ nDOFZ         , &
               nE + swE      , &
               nX(1) + swX(1), &
               nX(2) + swX(2), &
               nX(3) + swX(3), &
               nCR           , &
               nSpecies ]    , &
             CR )

    print*, "Writing Nodal Values"

#endif

  END SUBROUTINE WriteNodalDataToFile


  SUBROUTINE WriteEulerToFile( MF_uCF, MF_uGF, n )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    INTEGER, INTENT(in) :: n

#ifndef THORNADO_NOTRANSPORT

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)


    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: AF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF(:,:,:,:,:)

    REAL(DP) :: Mu(1:nDOFX), num

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iX1, iX2, iX3, l

    character(len=64) :: nm
    CHARACTER(LEN=16) :: FMT

    WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 AF )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )


        DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          CALL ComputePrimitive_Euler_Relativistic &
               ( U(:,iX1,iX2,iX3,iCF_D ),       &
                 U(:,iX1,iX2,iX3,iCF_S1),       &
                 U(:,iX1,iX2,iX3,iCF_S2),       &
                 U(:,iX1,iX2,iX3,iCF_S3),       &
                 U(:,iX1,iX2,iX3,iCF_E ),       &
                 U(:,iX1,iX2,iX3,iCF_Ne),       &
                 PF(:,iX1,iX2,iX3,iPF_D ),       &
                 PF(:,iX1,iX2,iX3,iPF_V1),       &
                 PF(:,iX1,iX2,iX3,iPF_V2),       &
                 PF(:,iX1,iX2,iX3,iPF_V3),       &
                 PF(:,iX1,iX2,iX3,iPF_E ),       &
                 PF(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

           CALL ComputeAuxiliary_Fluid_TABLE &
               ( PF(:,iX1,iX2,iX3,iPF_D ), &
                 PF(:,iX1,iX2,iX3,iPF_E ), &
                 PF(:,iX1,iX2,iX3,iPF_Ne), &
                 AF(:,iX1,iX2,iX3,iAF_P ), &
                 AF(:,iX1,iX2,iX3,iAF_T ), &
                 AF(:,iX1,iX2,iX3,iAF_Ye), &
                 AF(:,iX1,iX2,iX3,iAF_S ), &
                 AF(:,iX1,iX2,iX3,iAF_E ), &
                 AF(:,iX1,iX2,iX3,iAF_Gm), &
                 AF(:,iX1,iX2,iX3,iAF_Cs) )



          CALL ApplyEquationOfState_TABLE &
               ( PF(:,iX1,iX2,iX3,iPF_D ), &
                 AF(:,iX1,iX2,iX3,iAF_T ), &
                 AF(:,iX1,iX2,iX3,iAF_Ye), &
                 AF(:,iX1,iX2,iX3,iAF_P ), &
                 AF(:,iX1,iX2,iX3,iAF_S ), &
                 AF(:,iX1,iX2,iX3,iAF_E ), &
                 AF(:,iX1,iX2,iX3,iAF_Me), &
                 AF(:,iX1,iX2,iX3,iAF_Mp), &
                 AF(:,iX1,iX2,iX3,iAF_Mn), &
                 AF(:,iX1,iX2,iX3,iAF_Xp), &
                 AF(:,iX1,iX2,iX3,iAF_Xn), &
                 AF(:,iX1,iX2,iX3,iAF_Xa), &
                 AF(:,iX1,iX2,iX3,iAF_Xh), &
                 AF(:,iX1,iX2,iX3,iAF_Gm) )


        END DO
        END DO
        END DO
        WRITE(nm,*) n

        nm = ADJUSTL(nm)

        OPEN( UNIT = 101, FILE = 'DF.dat'      )
        OPEN( UNIT = 102, FILE = 'V1.dat'      )
        OPEN( UNIT = 103, FILE = 'V2.dat'      )
        OPEN( UNIT = 104, FILE = 'V3.dat'      )
        OPEN( UNIT = 105, FILE = 'Ye.dat'      )
        OPEN( UNIT = 106, FILE = 'Mu'//trim(nm) //'.dat'      )
        OPEN( UNIT = 107, FILE = 'E.dat'      )
        OPEN( UNIT = 108, FILE = 'T'//trim(nm) //'.dat'      )


        DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          Mu(:) = AF(:,iX1,iX2,iX3,iAF_Me) + AF(:,iX1,iX2,iX3,iAF_Mp) - AF(:,iX1,iX2,iX3,iAF_Mn)


          WRITE(101,FMT) &
            PF(:,iX1,iX2,iX3,iPF_D) &
              / unitsPF(iPF_D)

          WRITE(102,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V1) &
              / unitsPF(iPF_V1)

          WRITE(103,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V2) &
              / unitsPF(iPF_V2)

          WRITE(104,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V3) &
              / unitsPF(iPF_V3)

          WRITE(105,FMT) &
            AF(:,iX1,iX2,iX3,iAF_Ye) &
              / unitsAF(iAF_Ye)

          WRITE(106,FMT) &
            Mu(:) &
              / unitsAF(iAF_Me)

          WRITE(107,FMT) &
            PF(:,iX1,iX2,iX3,iPF_E) &
              / unitsPF(iPF_E)

          WRITE(108,FMT) &
            AF(:,iX1,iX2,iX3,iAF_T) &
              / MeV
        END DO
        END DO
        END DO

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 AF )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

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

    END DO

#endif

  END SUBROUTINE WriteEulerToFile


  SUBROUTINE AllocateArray_X( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(5), iHi(5)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:)

    CALL TimersStart_AMReX( Timer_AMReX_Allocate_X )

    ALLOCATE( A(iLo(1):iHi(1), &
                iLo(2):iHi(2), &
                iLo(3):iHi(3), &
                iLo(4):iHi(4), &
                iLo(5):iHi(5)) )

    CALL TimersStop_AMReX( Timer_AMReX_Allocate_X )

  END SUBROUTINE AllocateArray_X


  SUBROUTINE DeallocateArray_X( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(5), iHi(5)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:)

    CALL TimersStart_AMReX( Timer_AMReX_Allocate_X )

    DEALLOCATE( A )

    CALL TimersStop_AMReX( Timer_AMReX_Allocate_X )

  END SUBROUTINE DeallocateArray_X


  SUBROUTINE AllocateArray_Z( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(7), iHi(7)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:,:,:)

    CALL TimersStart_AMReX( Timer_AMReX_Allocate_Z )

    ALLOCATE( A(iLo(1):iHi(1), &
                iLo(2):iHi(2), &
                iLo(3):iHi(3), &
                iLo(4):iHi(4), &
                iLo(5):iHi(5), &
                iLo(6):iHi(6), &
                iLo(7):iHi(7)) )

    CALL TimersStop_AMReX( Timer_AMReX_Allocate_Z )

  END SUBROUTINE AllocateArray_Z


  SUBROUTINE DeallocateArray_Z( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(7), iHi(7)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:,:,:)

    CALL TimersStart_AMReX( Timer_AMReX_Allocate_Z )

    DEALLOCATE( A )

    CALL TimersStop_AMReX( Timer_AMReX_Allocate_Z )

  END SUBROUTINE DeallocateArray_Z

  SUBROUTINE AllocateArray_Integrated( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(6), iHi(6)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:,:)


    ALLOCATE( A(iLo(1):iHi(1), &
                iLo(2):iHi(2), &
                iLo(3):iHi(3), &
                iLo(4):iHi(4), &
                iLo(5):iHi(5), &
                iLo(6):iHi(6)) )


  END SUBROUTINE AllocateArray_Integrated


  SUBROUTINE DeallocateArray_Integrated( iLo, iHi, A )

    INTEGER ,              INTENT(in)    :: iLo(6), iHi(6)
    REAL(DP), ALLOCATABLE, INTENT(inout) :: A(:,:,:,:,:,:)


    DEALLOCATE( A )


  END SUBROUTINE DeallocateArray_Integrated

  SUBROUTINE amrex2amrex_permute_Z &
    ( nFields, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_amrex_permute )

    INTEGER,  INTENT(in)  :: nFields, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      Data_amrex           (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_amrex_permute   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iFd, iD, iNodeZ, iD_permute, iNodeX, iNodeE

    iD_permute = 0

    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iS  = 1      , nS
      DO iZ1    = iE_B0, iE_E0 ! always want iZ1 to not include ghost cells
      DO iNodeE = 1    , nDOFE
      DO iNodeX = 1    , nDOFX
      DO iFd = 1      , nFields


        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

        iD      = ( iS - 1 ) * nFields * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iFd - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        iD_permute = iD_permute + 1

        Data_amrex_permute(iZ2,iZ3,iZ4,iD_permute) &
          = Data_amrex(iZ2,iZ3,iZ4,iD)

      END DO
      END DO
      END DO
      END DO
      END DO
iD_permute = 0
    END DO
    END DO
    END DO


  END SUBROUTINE amrex2amrex_permute_Z

  SUBROUTINE amrex_permute2amrex_Z &
    ( nFields, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_amrex_permute )

    INTEGER,  INTENT(in)  :: nFields, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      Data_amrex_permute    (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_amrex            (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iFd, iD, iNodeZ, iD_permute, iNodeX, iNodeE

    iD_permute = 0

    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iS  = 1      , nS
      DO iZ1    = iE_B0, iE_E0 ! always want iZ1 to not include ghost cells
      DO iNodeE = 1    , nDOFE
      DO iNodeX = 1    , nDOFX
      DO iFd = 1      , nFields


        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

        iD      = ( iS - 1 ) * nFields * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iFd - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        iD_permute = iD_permute + 1
        Data_amrex(iZ2,iZ3,iZ4,iD) &
          = Data_amrex_permute(iZ2,iZ3,iZ4,iD_permute)

      END DO
      END DO
      END DO
      END DO
      END DO
iD_permute = 0
    END DO
    END DO
    END DO


  END SUBROUTINE amrex_permute2amrex_Z



  SUBROUTINE MF_amrex2amrex_permute_Z_Level &
    ( iLevel, nFields, MF_uGF, MF_uCR, MF_uCR_permute )

    INTEGER,  INTENT(in)  :: iLevel, nFields

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout)   :: MF_uCR_permute

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX


    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR_permute(:,:,:,:)

    INTEGER :: iLo_MF(4)
    INTEGER :: iZ_E0(4), iZ_E1(4), iZ_B0(4), iZ_B1(4)
    INTEGER :: iX_E0(3), iX_E1(3), iX_B0(3), iX_B1(3)


    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX



      iZ_B0(1)=iE_B0
      iZ_E0(1)=iE_E0

      iZ_B0(2:4)=iX_B0(1:3)
      iZ_E0(2:4)=iX_E0(1:3)

      iZ_B1(1)=iE_B1
      iZ_E1(1)=iE_E1

      iZ_B1(2:4)=iX_B1(1:3)
      iZ_E1(2:4)=iX_E1(1:3)

      uGF  => MF_uGF % DataPtr( MFI )
      uCR  => MF_uCR % DataPtr( MFI )
      uCR_permute => MF_uCR_permute % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

      CALL amrex2amrex_permute_Z &
           ( nFields, nSpecies, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
             iZ_B1, iZ_E1, uCR, uCR_permute )

    END DO

    CALL amrex_mfiter_destroy( MFI )






  END SUBROUTINE MF_amrex2amrex_permute_Z_Level

  SUBROUTINE MF_amrex_permute2amrex_Z_Level &
    ( iLevel, nFields, MF_uGF, MF_uCR, MF_uCR_permute )

    INTEGER,  INTENT(in)  :: iLevel, nFields

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR_permute

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX


    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR_permute(:,:,:,:)

    INTEGER :: iLo_MF(4)

    INTEGER :: iZ_E0(4), iZ_E1(4), iZ_B0(4), iZ_B1(4)
    INTEGER :: iX_E0(3), iX_E1(3), iX_B0(3), iX_B1(3)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )
    DO WHILE( MFI % next() )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX



      iZ_B0(1)=iE_B0
      iZ_E0(1)=iE_E0

      iZ_B0(2:4)=iX_B0(1:3)
      iZ_E0(2:4)=iX_E0(1:3)

      iZ_B1(1)=iE_B1
      iZ_E1(1)=iE_E1

      iZ_B1(2:4)=iX_B1(1:3)
      iZ_E1(2:4)=iX_E1(1:3)


      uGF  => MF_uGF  % DataPtr( MFI )
      uCR  => MF_uCR  % DataPtr( MFI )
      uCR_permute => MF_uCR_permute % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

      CALL amrex_permute2amrex_Z &
           ( nFields, nSpecies, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
             iZ_B1, iZ_E1, uCR, uCR_permute )

    END DO

    CALL amrex_mfiter_destroy( MFI )






  END SUBROUTINE MF_amrex_permute2amrex_Z_Level


  SUBROUTINE PrintBoxArray( BA )

    TYPE(amrex_boxarray), INTENT(in) :: BA

    CALL print_boxarray( BA % p )

  END SUBROUTINE PrintBoxArray


END MODULE MF_UtilitiesModule
