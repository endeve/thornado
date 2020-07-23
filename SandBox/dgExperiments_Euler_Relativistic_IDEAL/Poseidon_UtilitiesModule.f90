MODULE Poseidon_UtilitiesModule

  USE KindModule, ONLY: &
    DP, Three
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeSourceTerms_Poseidon


CONTAINS


  SUBROUTINE ComputeSourceTerms_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, U_Poseidon )

     INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
     REAL(DP), INTENT(in)  :: G         (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
     REAL(DP), INTENT(in)  :: U         (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
     REAL(DP), INTENT(out) :: U_Poseidon(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

     REAL(DP) :: uPF(nDOFX,nPF), Pressure(nDOFX)
     INTEGER  :: iX1, iX2, iX3

     DO iX3 = iX_B0(3), iX_E0(3)
     DO iX2 = iX_B0(2), iX_E0(2)
     DO iX1 = iX_B0(1), iX_E0(1)

       ! --- Compute trace of stress tensor ---

       CALL ComputePrimitive_Euler_Relativistic &
              ( U  (:,iX1,iX2,iX3,iCF_D ), &
                U  (:,iX1,iX2,iX3,iCF_S1), &
                U  (:,iX1,iX2,iX3,iCF_S2), &
                U  (:,iX1,iX2,iX3,iCF_S3), &
                U  (:,iX1,iX2,iX3,iCF_E ), &
                U  (:,iX1,iX2,iX3,iCF_Ne), &
                uPF(:,iPF_D ), &
                uPF(:,iPF_V1), &
                uPF(:,iPF_V2), &
                uPF(:,iPF_V3), &
                uPF(:,iPF_E ), &
                uPF(:,iPF_Ne), &
                G  (:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                G  (:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                G  (:,iX1,iX2,iX3,iGF_Gm_dd_33) )

       CALL ComputePressureFromPrimitive_IDEAL &
              ( uPF(:,iPF_D), uPF(:,iPF_E), uPF(:,iPF_Ne), Pressure )

       U_Poseidon(:,iX1,iX2,iX3,1) &
         = U(:,iX1,iX2,iX3,iCF_E) + U(:,iX1,iX2,iX3,iCF_D)

       U_Poseidon(:,iX1,iX2,iX3,2) &
         = U(:,iX1,iX2,iX3,iCF_S1) * uPF(:,iPF_V1) &
             + U(:,iX1,iX2,iX3,iCF_S2) * uPF(:,iPF_V2) &
             + U(:,iX1,iX2,iX3,iCF_S3) * uPF(:,iPF_V3) &
             + Three * Pressure

       U_Poseidon(:,iX1,iX2,iX3,3) &
         = U(:,iX1,iX2,iX3,iCF_S1) / G(:,iX1,iX2,iX3,iGF_Gm_dd_11)

       U_Poseidon(:,iX1,iX2,iX3,4) &
         = U(:,iX1,iX2,iX3,iCF_S2) / G(:,iX1,iX2,iX3,iGF_Gm_dd_22)

       U_Poseidon(:,iX1,iX2,iX3,5) &
         = U(:,iX1,iX2,iX3,iCF_S3) / G(:,iX1,iX2,iX3,iGF_Gm_dd_33)

     END DO
     END DO
     END DO

  END SUBROUTINE ComputeSourceTerms_Poseidon


END MODULE Poseidon_UtilitiesModule

