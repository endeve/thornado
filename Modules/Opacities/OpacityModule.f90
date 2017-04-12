MODULE OpacityModule

  USE KindModule, ONLY: &
    DP
  USE OpacityModule_IDEAL, ONLY: &
    InitializeOpacities_IDEAL, &
    FinalizeOpacities_IDEAL, &
    ComputeAbsorptionOpacity_IDEAL, &
    ComputeScatteringOpacity_ES_IDEAL
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE, &
    ComputeAbsorptionOpacity_TABLE, &
    ComputeScatteringOpacity_ES_TABLE

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5) :: &
    Opacity

  ! ---
  ! --- Interfaces for Various Opacity Functions and Subroutines ---
  ! ---

  INTERFACE
    SUBROUTINE ComputeOpacity_A( E, D, T, Y, Opacity )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y
      REAL(DP), DIMENSION(:,:), INTENT(out) :: Opacity
    END SUBROUTINE ComputeOpacity_A
  END INTERFACE

  ! ---
  ! --- Declaration of Opacity Functions and Subroutines ---
  ! ---

  PROCEDURE (ComputeOpacity_A), POINTER, PUBLIC :: &
    ComputeAbsorptionOpacity    => NULL(), &
    ComputeScatteringOpacity_ES => NULL() ! --- Elastic Scattering

  PUBLIC :: InitializeOpacities
  PUBLIC :: FinalizeOpacities

CONTAINS


  SUBROUTINE InitializeOpacities &
               ( Opacity_Option, OpacityTableName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Opacity_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Option

    Opacity = 'IDEAL'
    IF( PRESENT( Opacity_Option ) ) &
      Opacity = TRIM( Opacity_Option )

    WRITE(*,*)
    WRITE(*,'(A5,A11,A)') &
      '', 'Opacities: ', TRIM( Opacity )
    WRITE(*,'(A5,A11)') &
      '', '---------- '

    SELECT CASE ( TRIM( Opacity ) )
      CASE( 'IDEAL' )

        CALL InitializeOpacities_IDEAL

        ComputeAbsorptionOpacity &
          => ComputeAbsorptionOpacity_IDEAL
        ComputeScatteringOpacity_ES &
          => ComputeScatteringOpacity_ES_IDEAL

      CASE( 'TABLE' )

        CALL InitializeOpacities_TABLE &
               ( OpacityTableName_Option &
                   = OpacityTableName_Option )

        ComputeAbsorptionOpacity &
          => ComputeAbsorptionOpacity_TABLE
        ComputeScatteringOpacity_ES &
          => ComputeScatteringOpacity_ES_TABLE

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A17,A5)') &
          '', 'Invalid Opacity: ', TRIM( Opacity )
        STOP

    END SELECT

  END SUBROUTINE InitializeOpacities


  SUBROUTINE FinalizeOpacities

    SELECT CASE ( TRIM( Opacity ) )
      CASE( 'IDEAL' )

        CALL FinalizeOpacities_IDEAL

      CASE( 'TABLE' )

        CALL FinalizeOpacities_TABLE

    END SELECT

    NULLIFY &
      ( ComputeAbsorptionOpacity, &
        ComputeScatteringOpacity_ES )

  END SUBROUTINE FinalizeOpacities


END MODULE OpacityModule
