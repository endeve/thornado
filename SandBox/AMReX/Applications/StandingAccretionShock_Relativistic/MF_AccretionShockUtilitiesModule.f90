MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDimsX
  USE FluidFieldsModule, ONLY: &
    iPF_D, &
    iPF_V1, &
    iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE

  ! --- Local Modules ---

  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab
  USE InputParsingModule, ONLY: &
    nLevels
  USE FillPatchModule, ONLY: &
    FillPatch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteNodal1DICToFile_SAS

  LOGICAL,          PUBLIC              :: WriteNodal1DIC_SAS
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: FileName_Nodal1DIC_SAS

CONTAINS


  SUBROUTINE WriteNodal1DICToFile_SAS

    TYPE(amrex_parmparse) :: PP

    INTEGER :: iLevel

    CHARACTER(256) :: FileName, FMT

    IF( nDimsX .GT. 1 ) RETURN

    DO iLevel = 0, nLevels - 1

      CALL FillPatch( iLevel, MF_uGF )
      CALL FillPatch( iLevel, MF_uGF, MF_uCF )

    END DO

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
             swXX_Option = swX )

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'FileName_Nodal1DIC_SAS', &
                        FileName_Nodal1DIC_SAS )
    CALL amrex_parmparse_destroy( PP )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_D'
    CALL ShowVariableFromMultiFab &
      ( MF_uPF, iPF_D, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_V1'
    CALL ShowVariableFromMultiFab &
      ( MF_uPF, iPF_V1, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_AF_P'
    CALL ShowVariableFromMultiFab &
      ( MF_uAF, iAF_P, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    IF( amrex_parallel_ioprocessor() )THEN

      OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) // '_BC.dat' )

      WRITE(101,'(ES24.16E3)') ExpD
      WRITE(101,'(ES24.16E3)') ExpE

    END IF

  END SUBROUTINE WriteNodal1DICToFile_SAS


END MODULE MF_AccretionShockUtilitiesModule
