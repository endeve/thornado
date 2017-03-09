MODULE EulerEquationsUtilitiesModule_GR

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, nCF, &
    iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    iAF_P
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy    

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeConserved
  PUBLIC :: Flux_X1

CONTAINS


  SUBROUTINE ComputeConserved( uPF, uAF, uGF, uCF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uPF
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uAF
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uGF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uCF

    INTEGER  :: i
    REAL(DP) :: Gm_11, Gm_22, Gm_33, vSq, W, h

    DO i = 1, SIZE( uCF, DIM = 1 )

      Gm_11 = uGF(i,iGF_Gm_dd_11)
      Gm_22 = uGF(i,iGF_Gm_dd_22)
      Gm_33 = uGF(i,iGF_Gm_dd_33)

      vSq = Gm_11 * uPF(i,iPF_V1)**2 &
            + Gm_22 * uPF(i,iPF_V2)**2 &
            + Gm_33 * uPF(i,iPF_V3)**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = 1.0_DP + ( uPF(i,iPF_E) + uAF(i,iAF_P) ) / uPF(i,iPF_D)

      uCF(i,iCF_D)  = W * uPF(i,iPF_D)
      uCF(i,iCF_S1) = h * W**2 * Gm_11 * uPF(i,iPF_D) * uPF(i,iPF_V1)
      uCF(i,iCF_S2) = h * W**2 * Gm_22 * uPF(i,iPF_D) * uPF(i,iPF_V2)
      uCF(i,iCF_S3) = h * W**2 * Gm_33 * uPF(i,iPF_D) * uPF(i,iPF_V3)
      uCF(i,iCF_E)  = h * W**2 * uPF(i,iPF_D) - uAF(i,iAF_P) &
                        - W * uPF(i,iPF_D)
      uCF(i,iCF_Ne) = W * uPF(i,iPF_Ne)

    END DO

  END SUBROUTINE ComputeConserved


  SUBROUTINE ComputePrimitive( uCF, uGF, uPF, uAF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uCF
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uGF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uPF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uAF

    LOGICAL  :: Converged
    INTEGER  :: i
    REAL(DP) :: Gm11, Gm22, Gm33, SSq, P

    DO i = 1, SIZE( uPF, DIM = 1 )

      Gm11 = uGF(i,iGF_Gm_dd_11)
      Gm22 = uGF(i,iGF_Gm_dd_22)
      Gm33 = uGF(i,iGF_Gm_dd_33)

      SSq = Gm11 * uCF(i,iCF_S1)**2 &
            + Gm22 * uCF(i,iCF_S2)**2 &
            + Gm33 * uCF(i,iCF_S3)**2

      P = uAF(i,iAF_P)
      Converged = .FALSE.
      DO WHILE ( .NOT. Converged )

      END DO

    END DO

  END SUBROUTINE ComputePrimitive


  PURE FUNCTION Flux_X1 &
    ( D, V1, V2, V3, E, P, Ne, Alpha, Beta1, Gm11, Gm22, Gm33 )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, P, Ne
    REAL(DP), INTENT(in)       :: Alpha, Beta1, Gm11, Gm22, Gm33
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1

    REAL(DP) :: vSq, W, h

    vSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = 1.0_DP / SQRT( 1.0_DP - vSq )

    h   = 1.0_DP + ( E + P ) / D

    Flux_X1(iCF_D)  &
      = W * D * ( Alpha * V1 - Beta1 )

    Flux_X1(iCF_S1) &
      = D * h * W**2 * Gm11 * V1 * ( Alpha * V1 - Beta1 ) + Alpha * P

    Flux_X1(iCF_S2) &
      = D * h * W**2 * Gm22 * V2 * ( Alpha * V1 - Beta1 )

    Flux_X1(iCF_S3) &
      = D * h * W**2 * Gm33 * V3 * ( Alpha * V1 - Beta1 )

    Flux_X1(iCF_E)  &
      = ( D * h * W**2 - D * W ) * ( Alpha * V1 - Beta1 ) + Beta1 * P

    Flux_X1(iCF_Ne) &
      = ( Alpha * V1 - Beta1 ) * W * Ne

    RETURN
  END FUNCTION Flux_X1


  SUBROUTINE ComputeFunJac_P( D, SSq, E, P, Fun_P, Jac_P )

    REAL(DP), INTENT(in)  :: D, SSq, E, P
    REAL(DP), INTENT(out) :: Fun_P, Jac_P

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS
    REAL(DP), DIMENSION(1) :: Pbar

    HSq = ( E + P + D )**2

    RHO = D * SQRT( HSq - SSq ) / SQRT( HSq )
    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - D ) / D

    CALL ComputePressureFromSpecificInternalEnergy &
           ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    Fun_P = P - Pbar(1)

    dRHO = D * SSq / ( SQRT( HSq - SSq ) * HSq )
    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )

    Jac_P = 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS )

  END SUBROUTINE ComputeFunJac_P


END MODULE EulerEquationsUtilitiesModule_GR
