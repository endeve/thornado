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
    INTEGER                         :: iComp, nLevels
    TYPE(amrex_multifab)            :: MF_PF
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    WriteGF   = .FALSE.
    IF( PRESENT( MF_uGF_Option ) ) WriteGF   = .TRUE.

    WriteFF_C = .FALSE.
    IF( PRESENT( MF_uCF_Option ) ) WriteFF_C = .TRUE.

    WriteFF_P = .FALSE.
    IF( PRESENT( MF_uPF_Option ) ) WriteFF_P = .TRUE.

    WriteFF_A = .FALSE.
    IF( PRESENT( MF_uAF_Option ) ) WriteFF_A = .TRUE.


    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A18,I9.8)') '', 'Writing PlotFile: ', PlotFileNumber

    END IF

    WRITE(NumberString,'(I8.8)') PlotFileNumber

    nLevels = 1

    ! --- Geometry fields ---
    IF( WriteGF )THEN

      PlotFileName = TRIM( BaseFileName ) // '_GeometryFields_' // NumberString

      ALLOCATE( VarNames(nGF) )

      DO iComp = 1, nGF

        CALL amrex_string_build &
               ( VarNames(iComp), TRIM( ShortNamesGF(iComp) ) )

      END DO

      CALL amrex_multifab_build &
             ( MF_PF, MF_uGF_Option % BA, MF_uGF_Option % DM, nGF, 0 )

      CALL MF_ComputeCellAverage( nGF, MF_uGF_Option, MF_PF )

      CALL amrex_write_plotfile &
             ( PlotFileName, nLevels, [ MF_PF ], VarNames, &
               [ GEOM ], Time, [ PlotFileNumber ], [1] )

      CALL amrex_multifab_destroy( MF_PF )

      DEALLOCATE( VarNames )

    END IF

    ! --- Conserved fluid fields ---
    IF( WriteFF_C )THEN

      PlotFileName = TRIM( BaseFileName ) // '_FluidFields_C_' // NumberString

      ALLOCATE( VarNames(nCF) )

      DO iComp = 1, nCF

         CALL amrex_string_build &
                ( VarNames(iComp), TRIM( ShortNamesCF(iComp) ) )

      END DO

      CALL amrex_multifab_build &
             ( MF_PF, MF_uCF_Option % BA, MF_uCF_Option % DM, nCF, 0 )

      CALL MF_ComputeCellAverage( nCF, MF_uCF_Option, MF_PF )

      CALL amrex_write_plotfile &
             ( PlotFileName, nLevels, [ MF_PF ], VarNames, &
               [ GEOM ], Time, [ PlotFileNumber ], [1] )

      CALL amrex_multifab_destroy( MF_PF )

      DEALLOCATE( VarNames )

    END IF

    ! --- Primitive fluid fields ---
    IF( WriteFF_P )THEN

      PlotFileName = TRIM( BaseFileName ) // '_FluidFields_P_' // NumberString

      ALLOCATE( VarNames(nPF) )

      DO iComp = 1, nPF

         CALL amrex_string_build &
                ( VarNames(iComp), TRIM( ShortNamesPF(iComp) ) )

      END DO

      CALL amrex_multifab_build &
             ( MF_PF, MF_uPF_Option % BA, MF_uPF_Option % DM, nPF, 0 )

      CALL MF_ComputeCellAverage( nPF, MF_uPF_Option, MF_PF )

      CALL amrex_write_plotfile &
             ( PlotFileName, nLevels, [ MF_PF ], VarNames, &
               [ GEOM ], Time, [ PlotFileNumber ], [1] )

      CALL amrex_multifab_destroy( MF_PF )

      DEALLOCATE( VarNames )

    END IF

    ! --- Auxiliary fluid fields ---
    IF( WriteFF_A )THEN

      PlotFileName = TRIM( BaseFileName ) // '_FluidFields_A_' // NumberString

      ALLOCATE( VarNames(nAF) )

      DO iComp = 1, nAF

         CALL amrex_string_build &
                ( VarNames(iComp), TRIM( ShortNamesAF(iComp) ) )

      END DO

      CALL amrex_multifab_build &
             ( MF_PF, MF_uAF_Option % BA, MF_uAF_Option % DM, nAF, 0 )

      CALL MF_ComputeCellAverage( nAF, MF_uAF_Option, MF_PF )

      CALL amrex_write_plotfile &
             ( PlotFileName, nLevels, [ MF_PF ], VarNames, &
               [ GEOM ], Time, [ PlotFileNumber ], [1] )

      CALL amrex_multifab_destroy( MF_PF )

      DEALLOCATE( VarNames )

    END IF

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeCellAverage( nComp, MF, MF_A )

    INTEGER,              INTENT(in   ) :: nComp
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

          u_A(iX1,iX2,iX3,iComp) &
            = DOT_PRODUCT( WeightsX_q(:), u_K(:,iComp) )

        END DO

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeCellAverage


END MODULE InputOutputModuleAMReX
