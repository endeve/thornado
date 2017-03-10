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
  PUBLIC :: ComputePrimitive
  PUBLIC :: Flux_X1

CONTAINS


  SUBROUTINE ComputeConserved( uPF, uAF, uGF, uCF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uPF
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uAF
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uGF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uCF

    INTEGER  :: i
    REAL(DP) :: Gm11, Gm22, Gm33, vSq, W, h

    DO i = 1, SIZE( uCF, DIM = 1 )

      Gm11 = uGF(i,iGF_Gm_dd_11)
      Gm22 = uGF(i,iGF_Gm_dd_22)
      Gm33 = uGF(i,iGF_Gm_dd_33)

      vSq = Gm11 * uPF(i,iPF_V1)**2 &
            + Gm22 * uPF(i,iPF_V2)**2 &
            + Gm33 * uPF(i,iPF_V3)**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = 1.0_DP + ( uPF(i,iPF_E) + uAF(i,iAF_P) ) / uPF(i,iPF_D)

      uCF(i,iCF_D)  = W * uPF(i,iPF_D)
      uCF(i,iCF_S1) = h * W**2 * uPF(i,iPF_D) * Gm11 * uPF(i,iPF_V1)
      uCF(i,iCF_S2) = h * W**2 * uPF(i,iPF_D) * Gm22 * uPF(i,iPF_V2)
      uCF(i,iCF_S3) = h * W**2 * uPF(i,iPF_D) * Gm33 * uPF(i,iPF_V3)
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
    INTEGER  :: i, nIter
    REAL(DP) :: Gm11, Gm22, Gm33, SSq, vSq, W, h
    REAL(DP) :: Pold, Pnew, FunP, JacP
    REAL(DP), PARAMETER :: TolP = 1.0d-12

    DO i = 1, SIZE( uPF, DIM = 1 )

      Gm11 = uGF(i,iGF_Gm_dd_11)
      Gm22 = uGF(i,iGF_Gm_dd_22)
      Gm33 = uGF(i,iGF_Gm_dd_33)

      SSq = Gm11 * uCF(i,iCF_S1)**2 &
            + Gm22 * uCF(i,iCF_S2)**2 &
            + Gm33 * uCF(i,iCF_S3)**2

      ! --- Find Pressure with Newton's Method ---

      Pold = uAF(i,iAF_P) ! --- Initial Guess

      Converged = .FALSE.
      nIter     = 0
      DO WHILE ( .NOT. Converged )

        nIter = nIter + 1

        CALL ComputeFunJacP &
               ( uCF(i,iCF_D), SSq, uCF(i,iCF_E), Pold, FunP, JacP )

        Pnew = Pold - FunP / JacP

        IF( ABS( Pnew / Pold - 1.0_DP ) < TolP ) Converged = .TRUE.

        IF( nIter > 10 )THEN
          PRINT*, "ComputePrimitive"
          PRINT*, "  nIter = ", nIter
          PRINT*, "  Pold, Pnew, dP = ", &
            Pold, Pnew, ABS( Pnew / Pold - 1.0_DP )
          PRINT*, "  FunP = ", FunP
          PRINT*, "  Converged = ", Converged
          IF( nIter > 100) STOP
        END IF

        Pold = Pnew

      END DO

      uAF(i,iAF_P) = Pnew

      vSq = SSq / ( uCF(i,iCF_E) + Pnew + uCF(i,iCF_D) )**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = ( uCF(i,iCF_E) + Pnew + uCF(i,iCF_D) ) / ( W * uCF(i,iCF_D) )

      ! --- Recover Primitive Variables ---

      uPF(i,iPF_D)  = uCF(i,iCF_D) / W

      uPF(i,iPF_V1) = uCF(i,iCF_S1) / ( uCF(i,iCF_D) * h * W * Gm11 )

      uPF(i,iPF_V2) = uCF(i,iCF_S2) / ( uCF(i,iCF_D) * h * W * Gm22 )

      uPF(i,iPF_V3) = uCF(i,iCF_S3) / ( uCF(i,iCF_D) * h * W * Gm22 )

      uPF(i,iPF_E)  = uCF(i,iCF_D) * ( h - 1.0_DP ) / W - Pnew

      uPF(i,iPF_Ne) = uCF(i,iCF_Ne) / W

    END DO

  END SUBROUTINE ComputePrimitive


  FUNCTION Flux_X1 &
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
      = W * Ne * ( Alpha * V1 - Beta1 )

    RETURN
  END FUNCTION Flux_X1


  SUBROUTINE ComputeFunJacP( D, SSq, E, P, FunP, JacP )

    REAL(DP), INTENT(in)  :: D, SSq, E, P
    REAL(DP), INTENT(out) :: FunP, JacP

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS
    REAL(DP), DIMENSION(1) :: Pbar

    HSq = ( E + P + D )**2

    RHO = D * SQRT( HSq - SSq ) / SQRT( HSq )
    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - D ) / D

    CALL ComputePressureFromSpecificInternalEnergy &
           ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    FunP = P - Pbar(1)

    dRHO = D * SSq / ( SQRT( HSq - SSq ) * HSq )
    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )

    JacP = 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS )

  END SUBROUTINE ComputeFunJacP


END MODULE EulerEquationsUtilitiesModule_GR
