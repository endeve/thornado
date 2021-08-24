MODULE InputOutputModuleAMReX

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_plotfile_module, ONLY: &
    amrex_write_plotfile
  USE amrex_string_module, ONLY: &
    amrex_string, &
    amrex_string_build
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build, &
    amrex_geometry_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_amrcore, &
    amrex_max_level, &
    amrex_ref_ratio
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_amr_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    ShortNamesGF, &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    ShortNamesCF, &
    nCF, &
    ShortNamesPF, &
    nPF, &
    ShortNamesAF, &
    nAF, &
    ShortNamesDF, &
    nDF
  USE UnitsModule, ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    MaxGridSizeX, &
    dt, &
    StepNo, &
    swX, &
    t_new, &
    PlotFileBaseName, &
    nX, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteFieldsAMReX_PlotFile

CONTAINS


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uGF, &
      MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, &
      MF_uAF_Option, MF_uDF_Option )

    REAL(DP),             INTENT(in) :: Time
    INTEGER,              INTENT(in) :: StepNo(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uGF_Option(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uCF_Option(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uPF_Option(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uAF_Option(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uDF_Option(0:amrex_max_level)

    CHARACTER(08)                   :: NumberString
    CHARACTER(64)                   :: PlotFileName
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C
    LOGICAL                         :: WriteFF_P
    LOGICAL                         :: WriteFF_A
    LOGICAL                         :: WriteFF_D
    INTEGER                         :: iComp, iOS, iLevel, nF
    TYPE(amrex_multifab)            :: MF_plt(0:amrex_max_level)
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    nF = 3

    WriteGF = .FALSE.
    IF( PRESENT( MF_uGF_Option ) )THEN

      WriteGF = .TRUE.
      nF = nF + nGF

    END IF

    WriteFF_C = .FALSE.
    IF( PRESENT( MF_uCF_Option ) )THEN

      WriteFF_C = .TRUE.
      nF = nF + nCF

    END IF

    WriteFF_P = .FALSE.
    IF( PRESENT( MF_uPF_Option ) )THEN

      WriteFF_P = .TRUE.
      nF = nF + nPF

    END IF

    WriteFF_A = .FALSE.
    IF( PRESENT( MF_uAF_Option ) )THEN

      WriteFF_A = .TRUE.
      nF = nF + nAF

    END IF

    WriteFF_D = .FALSE.
    IF( PRESENT( MF_uDF_Option ) )THEN

      WriteFF_D = .TRUE.
      nF = nF + nDF

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(6x,A,I8.8)') 'Writing PlotFile: ', StepNo(0)
      WRITE(*,*)

    END IF

    WRITE(NumberString,'(I8.8)') StepNo(0)

    PlotFileName = TRIM( PlotFileBaseName ) // '_' // NumberString

    ALLOCATE( VarNames(nF) )

    CALL amrex_string_build( VarNames( 1 ), 'X1_C' )
    CALL amrex_string_build( VarNames( 2 ), 'X2_C' )
    CALL amrex_string_build( VarNames( 3 ), 'X3_C' )

    iOS = 3

    IF( WriteGF )THEN

      DO iComp = 1, nGF

        CALL amrex_string_build &
               ( VarNames( iComp + iOS ), TRIM( ShortNamesGF(iComp) ) )

      END DO

      iOS = iOS + nGF

    END IF

    IF( WriteFF_C )THEN

      DO iComp = 1, nCF

        CALL amrex_string_build &
               ( VarNames( iComp + iOS ), TRIM( ShortNamesCF(iComp) ) )

      END DO

      iOS = iOS + nCF

    END IF

    IF( WriteFF_P )THEN

      DO iComp = 1, nPF

        CALL amrex_string_build &
               ( VarNames( iComp + iOS ), TRIM( ShortNamesPF(iComp) ) )

      END DO

      iOS = iOS + nPF

    END IF

    IF( WriteFF_A )THEN

      DO iComp = 1, nAF

        CALL amrex_string_build &
               ( VarNames( iComp + iOS ), TRIM( ShortNamesAF(iComp) ) )

      END DO

      iOS = iOS + nAF

    END IF

    IF( WriteFF_D )THEN

      DO iComp = 1, nDF

        CALL amrex_string_build &
               ( VarNames( iComp + iOS ), TRIM( ShortNamesDF(iComp) ) )

      END DO

      iOS = iOS + nDF

    END IF

    DO iLevel = 0, amrex_max_level

      CALL amrex_multifab_build &
             ( MF_plt(iLevel), MF_uGF(iLevel) % BA, &
                               MF_uGF(iLevel) % DM, &
               nF, 0 )
      CALL MF_plt(iLevel) % setVal( Zero )

      CALL WriteMesh( iLevel, MF_uGF(iLevel), MF_plt(iLevel) )

      iOS = 3

      IF( WriteGF )THEN

        CALL ComputeCellAverage_MF &
               ( nGF, MF_uGF(iLevel), MF_uGF_Option(iLevel), &
                 iOS, 'GF', MF_plt(iLevel) )

        iOS = iOS + nGF

      END IF

      IF( WriteFF_C )THEN

        CALL ComputeCellAverage_MF &
               ( nCF, MF_uGF(iLevel), MF_uCF_Option(iLevel), &
                 iOS, 'CF', MF_plt(iLevel) )

        iOS = iOS + nCF

      END IF

      IF( WriteFF_P )THEN

        CALL ComputeCellAverage_MF &
               ( nPF, MF_uGF(iLevel), MF_uPF_Option(iLevel), &
                 iOS, 'PF', MF_plt(iLevel) )

        iOS = iOS + nPF

      END IF

      IF( WriteFF_A )THEN

        CALL ComputeCellAverage_MF &
               ( nAF, MF_uGF(iLevel), MF_uAF_Option(iLevel), &
                 iOS, 'AF', MF_plt(iLevel) )

        iOS = iOS + nAF

      END IF

      IF( WriteFF_D )THEN

        CALL ComputeCellAverage_MF &
               ( nDF, MF_uGF(iLevel), MF_uDF_Option(iLevel), &
                 iOS, 'DF', MF_plt(iLevel) )

        iOS = iOS + nDF

      END IF

    END DO ! End of loop over levels

    CALL amrex_write_plotfile &
           ( PlotFileName, amrex_max_level+1, MF_plt, VarNames, &
             amrex_geom, Time / UnitsDisplay % TimeUnit, StepNo, &
             amrex_ref_ratio )

    DO iLevel = 0, amrex_max_level

      CALL amrex_multifab_destroy ( MF_plt(iLevel) )

    END DO

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE ComputeCellAverage_MF &
    ( nComp, MF_uGF, MF, iOS, Field, MF_plt )

    INTEGER,              INTENT(in)    :: nComp, iOS
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(in)    :: MF
    CHARACTER(2),         INTENT(in)    :: Field
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER                       :: iX1, iX2, iX3, iComp
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_U(4), hi_U(4)
    REAL(DP)                      :: G_K(nDOFX,nGF)
    REAL(DP)                      :: U_K(nDOFX,nComp)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: G    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      G     => MF_uGF % DataPtr( MFI )
      U     => MF     % DataPtr( MFI )
      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        G_K(1:nDOFX,1:nGF) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        U_K(1:nDOFX,1:nComp) &
          = RESHAPE( U(iX1,iX2,iX3,lo_U(4):hi_U(4)), [ nDOFX, nComp ] )

        DO iComp = 1, nComp

          U_plt(iX1,iX2,iX3,iComp+iOS) &
            = SUM( WeightsX_q * U_K(:,iComp) * G_K(:,iGF_SqrtGm) ) &
                / SUM( WeightsX_q * G_K(:,iGF_SqrtGm) )

        END DO

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeCellAverage_MF


  SUBROUTINE WriteMesh( iLevel, MF_uGF, MF_plt )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER            :: iX1, iX2, iX3
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP) :: X1, X2, X3
    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        X1 = MeshX(1) % Center(iX1)
        X2 = MeshX(2) % Center(iX2)
        X3 = MeshX(3) % Center(iX3)

        U_plt(iX1,iX2,iX3,1) = X1
        U_plt(iX1,iX2,iX3,2) = X2
        U_plt(iX1,iX2,iX3,3) = X3

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE WriteMesh


END MODULE InputOutputModuleAMReX
