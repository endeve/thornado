MODULE MF_EquationOfStateModule_MHD

  ! --- AMReX Modules ---

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeEquationOfState_MF
  PUBLIC :: FinalizeEquationOfState_MF

  CHARACTER(:), ALLOCATABLE, PUBLIC :: EosTableName

CONTAINS


  SUBROUTINE InitializeEquationOfState_MF

    TYPE(amrex_parmparse) :: PP

    CHARACTER(:), ALLOCATABLE :: EquationOfState
    REAL(DP)                  :: Gamma_IDEAL

    Gamma_IDEAL  = -HUGE( One )
    EosTableName = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % get  ( 'EquationOfState', &
                        EquationOfState )
      CALL PP % query( 'Gamma_IDEAL', &
                        Gamma_IDEAL )
      CALL PP % query( 'EosTableName', &
                        EosTableName )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN

      CALL InitializeEquationOfState &
             ( EquationOfState_Option          = TRIM( EquationOfState ), &
               EquationOfStateTableName_Option = TRIM( EosTableName ), &
               Verbose_Option = amrex_parallel_ioprocessor() )

    ELSE IF( TRIM( EquationOfState ) .EQ. 'IDEAL' )THEN

      CALL InitializeEquationOfState &
               ( EquationOfState_Option = TRIM( EquationOfState ), &
                 Gamma_IDEAL_Option     = Gamma_IDEAL, &
                 Verbose_Option = amrex_parallel_ioprocessor() )

    ELSE

      CALL DescribeError_MF &
        ( 106, Message_Option = TRIM( EquationOfState ) )

    END IF

  END SUBROUTINE InitializeEquationOfState_MF


  SUBROUTINE FinalizeEquationOfState_MF

    CALL FinalizeEquationOfState

  END SUBROUTINE FinalizeEquationOfState_MF


END MODULE MF_EquationOfStateModule_MHD
