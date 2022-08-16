MODULE RefluxModule_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
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
    AverageDown
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Reflux_Euler_MF

  INTERFACE Reflux_Euler_MF
    MODULE PROCEDURE Reflux_Euler_MF_SingleLevel
    MODULE PROCEDURE Reflux_Euler_MF_MultipleLevels
  END INTERFACE Reflux_Euler_MF

CONTAINS


  SUBROUTINE Reflux_Euler_MF_MultipleLevels( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      IF( iLevel .GT. 0 ) &
        CALL Reflux_Euler_MF_SingleLevel( iLevel, MF_uGF, MF )

    END DO

    CALL AverageDown( MF_uGF, MF )

  END SUBROUTINE Reflux_Euler_MF_MultipleLevels


  SUBROUTINE Reflux_Euler_MF_SingleLevel( FineLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iErr

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL Reflux_Euler_MF_SingleLevel, FineLevel: ', FineLevel

    END IF

    CALL CreateMesh_MF( FineLevel-1, MeshX )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL FluxRegister( FineLevel ) &
           % reflux_dg( MF_uGF(FineLevel-1), MF(FineLevel-1), &
                        nCF, dX1, dX2, dX3 )

    END ASSOCIATE

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE Reflux_Euler_MF_SingleLevel


END MODULE RefluxModule_Euler
