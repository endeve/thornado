MODULE Poseidon_UtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    One, &
    Two, &
    Three
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputePressureFromPrimitive_TABLE
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeSourceTerms_Poseidon


CONTAINS


  SUBROUTINE ComputeSourceTerms_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Sources )

     INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
     REAL(DP), INTENT(in)  :: G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
     REAL(DP), INTENT(in)  :: U      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
     REAL(DP), INTENT(out) :: Sources(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

     REAL(DP) :: uGF(nGF), uPF(nPF), &
                 Pressure, LorentzFactor, Enthalpy, BetaDotV
     INTEGER  :: iNX, iX1, iX2, iX3, iGF

     INTEGER :: iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                             iX_B0(2):iX_E0(2), &
                             iX_B0(3):iX_E0(3))

     DO iX3 = iX_B0(3), iX_E0(3)
     DO iX2 = iX_B0(2), iX_E0(2)
     DO iX1 = iX_B0(1), iX_E0(1)
     DO iNX = 1, nDOFX

       iErr(iNX,iX1,iX2,iX3) = 0

       DO iGF = 1, nGF

         uGF(iGF) = G(iNX,iX1,iX2,iX3,iGF)

       END DO

       ! --- Compute trace of stress tensor ---

       CALL ComputePrimitive_Euler_Relativistic &
              ( U   (iNX,iX1,iX2,iX3,iCF_D ), &
                U   (iNX,iX1,iX2,iX3,iCF_S1), &
                U   (iNX,iX1,iX2,iX3,iCF_S2), &
                U   (iNX,iX1,iX2,iX3,iCF_S3), &
                U   (iNX,iX1,iX2,iX3,iCF_E ), &
                U   (iNX,iX1,iX2,iX3,iCF_Ne), &
                uPF (iPF_D ), &
                uPF (iPF_V1), &
                uPF (iPF_V2), &
                uPF (iPF_V3), &
                uPF (iPF_E ), &
                uPF (iPF_Ne), &
                uGF (iGF_Gm_dd_11), &
                uGF (iGF_Gm_dd_22), &
                uGF (iGF_Gm_dd_33), &
                iErr(iNX,iX1,iX2,iX3) )

#ifdef MICROPHYSICS_WEAKLIB

       CALL ComputePressureFromPrimitive_TABLE &
              ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

#else

       CALL ComputePressureFromPrimitive_IDEAL &
              ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

#endif

       Sources(iNX,iX1,iX2,iX3,1) &
         = U(iNX,iX1,iX2,iX3,iCF_E) + U(iNX,iX1,iX2,iX3,iCF_D)

       Sources(iNX,iX1,iX2,iX3,2) &
         =   U(iNX,iX1,iX2,iX3,iCF_S1) * uPF(iPF_V1) &
           + U(iNX,iX1,iX2,iX3,iCF_S2) * uPF(iPF_V2) &
           + U(iNX,iX1,iX2,iX3,iCF_S3) * uPF(iPF_V3) &
           + Three * Pressure

       Sources(iNX,iX1,iX2,iX3,3) &
         = U(iNX,iX1,iX2,iX3,iCF_S1) / uGF(iGF_Gm_dd_11)

       Sources(iNX,iX1,iX2,iX3,4) &
         = U(iNX,iX1,iX2,iX3,iCF_S2) / uGF(iGF_Gm_dd_22)

       Sources(iNX,iX1,iX2,iX3,5) &
         = U(iNX,iX1,iX2,iX3,iCF_S3) / uGF(iGF_Gm_dd_33)

       LorentzFactor &
         = One / SQRT( One                              &
             - ( uGF(iGF_Gm_dd_11) * uPF(iPF_V1)**2 &
               + uGF(iGF_Gm_dd_22) * uPF(iPF_V2)**2 &
               + uGF(iGF_Gm_dd_33) * uPF(iPF_V3)**2 ) )

       BetaDotV =   uGF(iGF_Gm_dd_11) * uGF(iGF_Beta_1) * uPF(iPF_V1) &
                  + uGF(iGF_Gm_dd_22) * uGF(iGF_Beta_2) * uPF(iPF_V2) &
                  + uGF(iGF_Gm_dd_33) * uGF(iGF_Beta_3) * uPF(iPF_V3)

       Enthalpy = uPF(iPF_D) + uPF(iPF_E) + Pressure

       Sources(iNX,iX1,iX2,iX3,6) &
         = Enthalpy * ( Two * LorentzFactor**2             &
             * ( One - BetaDotV / uGF(iGF_Alpha) ) - One ) &
             + Two * Pressure

     END DO
     END DO
     END DO
     END DO

     IF( ANY( iErr .NE. 0 ) )THEN

       WRITE(*,*) 'ERROR'
       WRITE(*,*) 'Module: Poseidon_UtilitiesModule'
       WRITE(*,*) 'Subroutine: ComputeSourceTerms_Poseidon'

       DO iX3 = iX_B0(3), iX_E0(3)
       DO iX2 = iX_B0(2), iX_E0(2)
       DO iX1 = iX_B0(1), iX_E0(1)
       DO iNX = 1, nDOFX

         IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

           WRITE(*,'(2x,A,4I5.4)') 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3
           CALL DescribeError_Euler( iErr(iNX,iX1,iX2,iX3) )

         END IF

       END DO
       END DO
       END DO
       END DO

     END IF

  END SUBROUTINE ComputeSourceTerms_Poseidon


END MODULE Poseidon_UtilitiesModule

