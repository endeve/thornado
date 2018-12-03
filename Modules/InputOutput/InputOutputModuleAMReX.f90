MODULE InputOutputModuleAMReX

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    nGF

  IMPLICIT NONE
  PRIVATE

  INTEGER :: PlotFileNumber = 0

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

    LOGICAL :: WriteGF
    LOGICAL :: WriteFF_C
    INTEGER :: MF_nComp
    TYPE(amrex_multifab) :: MF_PF

    IF( PRESENT( MF_uGF_Option ) )THEN
      WriteGF  = .TRUE.
    ELSE
      WriteGF  = .FALSE.
    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A18,I9.8)') '', 'Writing PlotFile: ', PlotFileNumber

    END IF

    IF( WriteGF )THEN

      CALL amrex_multifab_build &
             ( MF_PF, MF_uGF_Option % BA, MF_uGF_Option % DM, nGF, 0 )

      CALL MF_ComputeElementAverage( nGF, MF_uGF_Option, MF_PF )

      CALL amrex_multifab_destroy( MF_PF )

    END IF

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeElementAverage( nComp, MF, MF_A )

    INTEGER,              INTENT(in   ) :: nComp
    TYPE(amrex_multifab), INTENT(in   ) :: MF
    TYPE(amrex_multifab), INTENT(inout) :: MF_A

  END SUBROUTINE MF_ComputeElementAverage


END MODULE InputOutputModuleAMReX
