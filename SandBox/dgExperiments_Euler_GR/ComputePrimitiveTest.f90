PROGRAM ComputePrimitiveTest

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputePrimitive_GR

  IMPLICIT NONE

  REAL(DP) :: U(1,nCF), P(1,nPF), G(1,nGF), Pressure(1)

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = 4.0_DP / 3.0_DP )

  U(:,iCF_D)  =  2.3052949117970522E-017_DP
  U(:,ICF_S1) = -3.5282989837921351E-018_DP
  U(:,iCF_S2) =  0.0_DP
  U(:,iCF_S3) =  0.0_DP
  U(:,iCF_E)  =  2.6254050919391916E-019_DP
  U(:,iCF_Ne) =  0.0_DP

  G(:,iGF_Gm_dd_11) = 1.0226525967671742E+000_DP
  G(:,iGF_Gm_dd_22) = 3.5000285124356529E+010_DP
  G(:,iGF_Gm_dd_33) = 3.5000285124356529E+010_DP

  Pressure(:) = 0.0_DP

  PRINT*, "P(in)  = ", Pressure

  CALL ComputePrimitive_GR &
         ( U(:,iCF_D), U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), U(:,iCF_E), U(:,iCF_Ne), &
           P(:,iPF_D), P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), P(:,iPF_E), P(:,iPF_Ne), &
           Pressure(:), G(:,iGF_Gm_dd_11), G(:,iGF_Gm_dd_22), G(:,iGF_Gm_dd_33) )

  PRINT*, "P(out) = ", Pressure

  CALL FinalizeEquationOfState

END PROGRAM ComputePrimitiveTest
