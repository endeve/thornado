MODULE OpacityModule

  ! --- weaklib modules --------------------------

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType

  ! ----------------------------------------------

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5) :: &
    Opacity &
      = 'IDEAL'
  CHARACTER(256) :: &
    OpacityTableName &
      = 'OpacityTable.h5'
  TYPE(OpacityTableType) :: &
    OPACITIES

  PUBLIC :: InitializeOpacities
  PUBLIC :: FinalizeOpacities

CONTAINS


  SUBROUTINE InitializeOpacities &
               ( Opacity_Option, OpacityTableName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Opacity_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Option

    IF( PRESENT( Opacity_Option ) )THEN
      Opacity = TRIM( Opacity_Option )
    END IF

    IF( PRESENT( OpacityTableName_Option ) )THEN
      OpacityTableName = TRIM( OpacityTableName_Option )
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A9,A)') &
      '', 'Opacity: ', TRIM( Opacity )
    WRITE(*,'(A5,A9)') &
      '', '-------- '

    SELECT CASE ( TRIM( Opacity ) )
      CASE( 'IDEAL' )

      CASE( 'TABLE' )

        WRITE(*,*)
        WRITE(*,'(A7,A12,A)') &
          '', 'Table Name: ', TRIM( OpacityTableName )

        CALL InitializeHDF( )

        CALL ReadOpacityTableHDF &
               ( OPACITIES, TRIM( OpacityTableName ) )

        CALL FinalizeHDF( )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A17,A5)') &
          '', 'Invalid Opacity: ', TRIM( Opacity )
        STOP

    END SELECT

    STOP

  END SUBROUTINE InitializeOpacities


  SUBROUTINE FinalizeOpacities

  END SUBROUTINE FinalizeOpacities


END MODULE OpacityModule
