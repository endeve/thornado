MODULE RefluxModule_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- thornado Modules ---

  USE GeometryFieldsModule, ONLY: &
    nGF
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
  USE MF_UtilitiesModule, ONLY: &
    MultiplyWithMetric
  USE InputParsingModule, ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Reflux_Euler_MF

  INTERFACE Reflux_Euler_MF
    MODULE PROCEDURE Reflux_Euler_MF_SingleLevel
    MODULE PROCEDURE Reflux_Euler_MF_MultipleLevels
  END INTERFACE Reflux_Euler_MF

CONTAINS


  SUBROUTINE Reflux_Euler_MF_MultipleLevels( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:nLevels-1)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      IF( iLevel .GT. 0 ) &
        CALL Reflux_Euler_MF_SingleLevel( iLevel, MF_uGF, MF )

    END DO

  END SUBROUTINE Reflux_Euler_MF_MultipleLevels


  SUBROUTINE Reflux_Euler_MF_SingleLevel( FineLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:nLevels-1)

    CALL CreateMesh_MF( FineLevel-1, MeshX )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL FluxRegister( FineLevel ) &
           % reflux_dg( MF_uGF(FineLevel-1), MF(FineLevel-1), &
                        nCF, dX1, dX2, dX3 )

    END ASSOCIATE

    CALL DestroyMesh_MF( MeshX )

    ! --- MF_uGF must be LAST ---
    CALL MultiplyWithMetric( MF_uGF(FineLevel), MF    (FineLevel), nCF, +1 )
    CALL MultiplyWithMetric( MF_uGF(FineLevel), MF_uGF(FineLevel), nGF, +1 )

    CALL AverageDownTo( FineLevel-1, MF_uGF )
    CALL AverageDownTo( FineLevel-1, MF     )

    ! --- MF_uGF must be FIRST ---
    CALL MultiplyWithMetric( MF_uGF(FineLevel  ), MF_uGF(FineLevel  ), nGF, -1 )
    CALL MultiplyWithMetric( MF_uGF(FineLevel  ), MF    (FineLevel  ), nCF, -1 )
    CALL MultiplyWithMetric( MF_uGF(FineLevel-1), MF_uGF(FineLevel-1), nGF, -1 )
    CALL MultiplyWithMetric( MF_uGF(FineLevel-1), MF    (FineLevel-1), nCF, -1 )

  END SUBROUTINE Reflux_Euler_MF_SingleLevel


END MODULE RefluxModule_Euler
