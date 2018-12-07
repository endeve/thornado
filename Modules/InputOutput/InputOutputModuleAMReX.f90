MODULE InputOutputModuleAMReX

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP
  USE ProgramHeaderModule,     ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE GeometryFieldsModule,    ONLY: &
    nGF, ShortNamesGF
  USE FluidFieldsModule,       ONLY: &
    nCF, ShortNamesCF, &
    nPF, ShortNamesPF, &
    nAF, ShortNamesAF

  IMPLICIT NONE
  PRIVATE

  CHARACTER(8) :: BaseFileName = 'thornado'
  INTEGER      :: PlotFileNumber = 0

  PUBLIC :: WriteFieldsAMReX_PlotFile

CONTAINS


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, GEOM, MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, MF_uAF_Option )

    REAL(DP),             INTENT(in)           :: Time
    TYPE(amrex_geometry), INTENT(in)           :: GEOM
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uAF_Option

    CHARACTER(08)                   :: NumberString
    CHARACTER(32)                   :: PlotFileName
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C, WriteFF_P, WriteFF_A
    INTEGER                         :: iComp, nLevels, nF = 0
    TYPE(amrex_multifab)            :: MF_PF
    TYPE(amrex_boxarray)            :: BA
    TYPE(amrex_distromap)           :: DM
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    INTEGER :: iTemp

    WriteGF   = .FALSE.
    IF( PRESENT( MF_uGF_Option ) )THEN
      WriteGF   = .TRUE.
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

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A18,I9.8)') '', 'Writing PlotFile: ', PlotFileNumber

    END IF

    WRITE(NumberString,'(I8.8)') PlotFileNumber

    nLevels = 1

    PlotFileName = TRIM( BaseFileName ) // '_' // NumberString

    ALLOCATE( VarNames(nF) )

    iTemp = 0
    IF( WriteGF )THEN
      DO iComp = 1, nGF
        CALL amrex_string_build &
               ( VarNames(iComp + iTemp), TRIM( ShortNamesGF(iComp) ) )
      END DO
      iTemp = iTemp + nGF
      BA = MF_uGF_Option % BA
      DM = MF_uGF_Option % DM
    END IF

    IF( WriteFF_C )THEN
      DO iComp = 1, nCF
        CALL amrex_string_build &
               ( VarNames(iComp + iTemp), TRIM( ShortNamesCF(iComp) ) )
      END DO
      iTemp = iTemp + nCF
      BA = MF_uCF_Option % BA
      DM = MF_uCF_Option % DM
    END IF

    IF( WriteFF_P )THEN
      DO iComp = 1, nPF
        CALL amrex_string_build &
               ( VarNames(iComp + iTemp), TRIM( ShortNamesPF(iComp) ) )
      END DO
      iTemp = iTemp + nPF
      BA = MF_uPF_Option % BA
      DM = MF_uPF_Option % DM
    END IF

    IF( WriteFF_A )THEN
      DO iComp = 1, nAF
        CALL amrex_string_build &
               ( VarNames(iComp + iTemp), TRIM( ShortNamesAF(iComp) ) )
      END DO
      iTemp = iTemp + nAF
      BA = MF_uAF_Option % BA
      DM = MF_uAF_Option % DM
    END IF

    CALL amrex_multifab_build &
           ( MF_PF, BA, DM, nF, 0 )

    iTemp = 0
    IF( WriteGF   ) &
      CALL MF_ComputeCellAverage( nGF, MF_uGF_Option, MF_PF, iTemp )
    IF( WriteFF_C ) &
      CALL MF_ComputeCellAverage( nCF, MF_uCF_Option, MF_PF, iTemp )
    IF( WriteFF_P ) &
      CALL MF_ComputeCellAverage( nPF, MF_uPF_Option, MF_PF, iTemp )
    IF( WriteFF_A ) &
      CALL MF_ComputeCellAverage( nAF, MF_uAF_Option, MF_PF, iTemp )

    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels, [ MF_PF ], VarNames, &
             [ GEOM ], Time, [ PlotFileNumber ], [1] )

    CALL amrex_multifab_destroy( MF_PF )

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeCellAverage( nComp, MF, MF_A, iTemp )

    INTEGER,              INTENT(in   ) :: nComp
    INTEGER,              INTENT(inout) :: iTemp
    TYPE(amrex_multifab), INTENT(in   ) :: MF
    TYPE(amrex_multifab), INTENT(inout) :: MF_A

    INTEGER            :: iX1, iX2, iX3, iComp
    INTEGER            :: lo(4), hi(4)
    REAL(amrex_real)   :: u_K(nDOFX,nComp)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: u  (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: u_A(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      u   => MF   % DataPtr( MFI )
      u_A => MF_A % DataPtr( MFI )

      BX = MFI % tilebox()

      lo = LBOUND( u ); hi = UBOUND( u )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        u_K(1:nDOFX,1:nComp) &
          = RESHAPE( u(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nComp ] )

        ! --- Compute cell-average ---
        DO iComp = 1, nComp

          u_A(iX1,iX2,iX3,iComp + iTemp) &
            = DOT_PRODUCT( WeightsX_q(:), u_K(:,iComp) )

        END DO

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )
    iTemp = iTemp + nComp

  END SUBROUTINE MF_ComputeCellAverage


END MODULE InputOutputModuleAMReX
