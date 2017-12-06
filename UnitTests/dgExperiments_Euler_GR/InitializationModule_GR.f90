MODULE InitializationModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P, iAF_Gm
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeConserved_GR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_GR

CONTAINS


  SUBROUTINE InitializeFields_GR

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! Loop over elements
    DO iX3 = 0, nX(3) + 1
      DO iX2 = 0, nX(2) + 1
        DO iX1 = 0, nX(1) + 1

           ! Loop over all nodes in an element
           DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) ! Physical coordinate

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = One + Half * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = One * 1.0d-1
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = Zero
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = Zero
            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = 0.01_DP * One
            uAF( iNodeX , iX1 , iX2 , iX3 , iAF_P )       &
              = ( 4.0_DP / 3.0_DP - 1.0_DP )              &
                * uPF( iNodeX , iX1 , iX2 , iX3 , iPF_E )            

           END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_GR


END MODULE InitializationModule_GR
