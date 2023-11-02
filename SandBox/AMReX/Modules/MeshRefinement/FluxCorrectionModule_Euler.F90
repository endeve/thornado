MODULE FluxCorrectionModule_Euler

  USE Euler_MeshRefinementModule, ONLY: &
    pNodeNumberTableX_X1_c, &
    pNodeNumberTableX_X2_c, &
    pNodeNumberTableX_X3_c, &
    pWeightsX_q_c, &
    pLX_X1_Up_c, pLX_X1_Dn_c, &
    pLX_X2_Up_c, pLX_X2_Dn_c, &
    pLX_X3_Up_c, pLX_X3_Dn_c

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
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF
  USE ReferenceElementModuleX_Lagrange, ONLY: &
     LX_X1_Dn, &
     LX_X1_Up, &
     LX_X2_Dn, &
     LX_X2_Up, &
     LX_X3_Dn, &
     LX_X3_Up
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    NodeNumberTableX, &
    NodeNumberTableX_X1, &
    NodeNumberTableX_X2, &
    NodeNumberTableX_X3, &
    WeightsX_q

  ! --- Local Modules ---

  USE MF_FieldsModule_Euler, ONLY: &
    FluxRegister_Euler
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyFluxCorrection_Euler_MF

  INTERFACE ApplyFluxCorrection_Euler_MF
    MODULE PROCEDURE ApplyFluxCorrection_Euler_MF_SingleLevel
    MODULE PROCEDURE ApplyFluxCorrection_Euler_MF_MultipleLevels
  END INTERFACE ApplyFluxCorrection_Euler_MF

CONTAINS


  SUBROUTINE ApplyFluxCorrection_Euler_MF_MultipleLevels( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      IF( iLevel .GT. 0 ) &
        CALL ApplyFluxCorrection_Euler_MF_SingleLevel( iLevel, MF_uGF, MF )

    END DO

    CALL AverageDown( MF_uGF, MF )

  END SUBROUTINE ApplyFluxCorrection_Euler_MF_MultipleLevels


  SUBROUTINE ApplyFluxCorrection_Euler_MF_SingleLevel( FineLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iErr

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL ApplyFluxCorrection_Euler_MF_SingleLevel, FineLevel: ', FineLevel

    END IF

#if defined( THORNADO_USE_MESHREFINEMENT )

    CALL CreateMesh_MF( FineLevel-1, MeshX )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    CALL FluxRegister_Euler( FineLevel ) &
           % reflux_dg( MF_uGF(FineLevel-1), MF(FineLevel-1), &
                        nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, nCF, iGF_SqrtGm, &
                        pNodeNumberTableX_X1_c, &
                        pNodeNumberTableX_X2_c, &
                        pNodeNumberTableX_X3_c, &
                        pWeightsX_q_c, &
                        pLX_X1_Up_c, pLX_X1_Dn_c, &
                        pLX_X2_Up_c, pLX_X2_Dn_c, &
                        pLX_X3_Up_c, pLX_X3_Dn_c, &
                        MINVAL( dX1 ), MINVAL( dX2 ), MINVAL( dX3 ) )

    END ASSOCIATE

    CALL DestroyMesh_MF( MeshX )

#endif

  END SUBROUTINE ApplyFluxCorrection_Euler_MF_SingleLevel


END MODULE FluxCorrectionModule_Euler
