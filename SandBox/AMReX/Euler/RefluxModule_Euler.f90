MODULE RefluxModule_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- thornado Modules ---

  USE FluidFieldsModule, ONLY: &
    nCF
  USE MeshModule, ONLY: &
    MeshX

  ! --- Local Modules ---

  USE MF_FieldsModule, ONLY: &
    FluxRegister
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE AverageDownModule, ONLY: &
    AverageDownTo
  USE InputParsingModule, ONLY: &
    nLevels, &
    do_reflux

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Reflux_Euler_MF

  INTERFACE Reflux_Euler_MF
    MODULE PROCEDURE Reflux_Euler_MF_SingleLevel
    MODULE PROCEDURE Reflux_Euler_MF_MultipleLevels
  END INTERFACE Reflux_Euler_MF


CONTAINS


  SUBROUTINE Reflux_Euler_MF_MultipleLevels( MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF(0:nLevels-1)

    INTEGER :: iLevel

    IF( .NOT. do_reflux ) RETURN

    DO iLevel = 0, nLevels-1

      IF( iLevel .GT. 0 ) &
        CALL Reflux_Euler_MF_SingleLevel( iLevel, MF )

    END DO

  END SUBROUTINE Reflux_Euler_MF_MultipleLevels


  SUBROUTINE Reflux_Euler_MF_SingleLevel( FineLevel, MF )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:nLevels-1)

    CALL CreateMesh_MF( FineLevel-1, MeshX )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL FluxRegister( FineLevel ) &
           % reflux_dg( MF(FineLevel-1), nCF, dX1, dX2, dX3 )

    END ASSOCIATE

    CALL DestroyMesh_MF( MeshX )

    CALL AverageDownTo( FineLevel-1, MF )

  END SUBROUTINE Reflux_Euler_MF_SingleLevel


END MODULE RefluxModule_Euler
