MODULE EulerEquationsUtilitiesModule_Beta_GR

  USE ProgramHeaderModule, ONLY: &
    nX
  USE KindModule, ONLY: &
    DP, Zero, Half
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne

  IMPLICIT NONE
  PRIVATE :: ComputeFunJacP

  PUBLIC :: ComputePrimitive_GR
  PUBLIC :: ComputeConserved_GR
  PUBLIC :: Eigenvalues_GR
  PUBLIC :: AlphaC_GR
  PUBLIC :: ComputeSoundSpeed_GR
  PUBLIC :: Flux_X1_GR
  PUBLIC :: NumericalFlux_X1_HLLC_GR

CONTAINS


  SUBROUTINE ComputePrimitive_GR( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                                  PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                                  GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
                                  AF_P, K )

    REAL(DP), DIMENSION(:), INTENT(in)  :: CF_D, CF_S1, CF_S2, CF_S3, &
                                           CF_E, CF_Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: PF_D, PF_V1, PF_V2, PF_V3, &
                                           PF_E, PF_Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33
    REAL(DP), DIMENSION( 1 : K ), INTENT(out) :: AF_P

    LOGICAL :: Converged
    INTEGER  :: K, i, nIter, nNodes

    REAL(DP), DIMENSION( 1 : K ) :: SSq, Pold, vSq, W, h, Pnew

    REAL(DP) :: FunP, JacP
    REAL(DP), PARAMETER :: TolP = 1.0d-12

    SSq = GF_Gm_dd_11 * CF_S1**2 &
        + GF_Gm_dd_22 * CF_S2**2 &
        + GF_Gm_dd_33 * CF_S3**2

    ! --- Find Pressure with Newton's Method ---

    Pold = AF_P ! -- Initial guess
    nNodes = SIZE( Pold )
    
    ! Loop through all the nodes
    DO i = 1, nNodes
       
      Converged = .FALSE.
      nIter     = 0

      DO WHILE ( .NOT. Converged )

        nIter = nIter + 1

        CALL ComputeFunJacP( CF_D(i), SSq(i), CF_E(i), Pold(i), FunP, JacP )

        Pnew(i) = Pold(i) - FunP / JacP

        ! Check if Newton's method has converged
        IF( ABS( Pnew(i) / Pold(i) - 1.0_DP ) <= TolP ) Converged = .TRUE.

        ! For de-bugging
        IF( nIter == 10)THEN
          WRITE(*,*) 'No convergence, |ERROR|:', &
                      ABS( Pnew(i) / Pold(i) - 1.0_DP )
          WRITE(*,*) 'Pold:                   ', Pold(i)
          WRITE(*,*) 'Pnew:                   ', Pnew(i)
          STOP
        END IF

        Pold(i) = Pnew(i)

      END DO

    END DO

    AF_P = Pnew

    vSq = SSq / ( CF_E + Pnew + CF_D )**2

    W = 1.0_DP / SQRT( 1.0_DP - vSq )

    h = ( CF_E + Pnew + CF_D ) / ( W * CF_D )

    ! --- Recover Primitive Variables ---

    PF_D  = CF_D / W

    PF_V1 = CF_S1 / ( CF_D * h * W * GF_Gm_dd_11 )

    PF_V2 = CF_S2 / ( CF_D * h * W * GF_Gm_dd_22 )

    PF_V3 = CF_S3 / ( CF_D * h * W * GF_Gm_dd_33 )

    PF_E  = CF_D * ( h - 1.0_DP ) / W - Pnew

    PF_Ne = CF_Ne / W

  END SUBROUTINE ComputePrimitive_GR


  SUBROUTINE ComputeConserved_GR( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                                  CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                                  GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
                                  AF_P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: PF_D, PF_V1, PF_V2, PF_V3, &
                                           PF_E, PF_Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: CF_D, CF_S1, CF_S2, CF_S3, &
                                           CF_E, CF_Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33
    REAL(DP), DIMENSION(:), INTENT(in)  :: AF_P

    REAL(DP), DIMENSION( 1 : nX( 1 ) ) :: vSq, W, h

    vSq = GF_Gm_dd_11 * PF_V1**2 &
        + GF_Gm_dd_22 * PF_V2**2 &
        + GF_Gm_dd_33 * PF_V3**2

    W = 1.0_DP / SQRT( 1.0_DP - vSq )
    h = 1.0_DP + ( PF_E + AF_P ) / PF_D
    
    CF_D   = W * PF_D
    CF_S1  = h * W**2 * PF_D * GF_Gm_dd_11 * PF_V1
    CF_S2  = h * W**2 * PF_D * GF_Gm_dd_22 * PF_V2
    CF_S3  = h * W**2 * PF_D * GF_Gm_dd_33 * PF_V3
    CF_E   = h * W**2 * PF_D - AF_P - W * PF_D
    CF_Ne  = W * PF_Ne

  END SUBROUTINE ComputeConserved_GR
  

  SUBROUTINE ComputeFunJacP( D, SSq, E, P, FunP, JacP )

    REAL(DP), INTENT(in)  :: D, SSq, E, P
    REAL(DP), INTENT(out) :: FunP, JacP

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS
    REAL(DP), DIMENSION(1) :: Pbar

    HSq = ( E + P + D )**2

    RHO = D * SQRT( HSq - SSq ) / SQRT( HSq )

    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - D ) / D

!    CALL ComputePressureFromSpecificInternalEnergy &
!         ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )
!    FunP = P - Pbar(1)
    
    Pbar(1) = ( 4.0_DP / 3.0_DP - 1.0_DP ) * RHO * EPS

    FunP = P - Pbar(1)

    dRHO = D * SSq / ( SQRT( HSq - SSq ) * HSq )
    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )
 
    JacP = 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS )

  END SUBROUTINE ComputeFunJacP

  
  SUBROUTINE Eigenvalues_GR( V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                Alpha, Beta1, Cs )

    ! Alpha is the lapse function
    ! Vi is the contravariant component V^i
    ! Beta1 is the contravariant component Beta^1

    REAL(DP), INTENT(in) :: V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                            Alpha, Beta1, Cs
    REAL(DP) :: VSq
    INTEGER :: i
    REAL(DP), DIMENSION(1:nX(1), 1:5) :: Eigvals_GR

    VSq = Gm_dd_11 * V1**2 + Gm_dd_22 * V2**2 + Gm_dd_33 * V3**2

    DO i = 1, nX(1)

       Eigvals_GR(i, 1:5) = &
          [ Alpha / ( 1.0_DP - VSq * Cs**2 ) * ( V1 * ( 1.0_DP - Cs**2 ) - Cs &
          * SQRT( ( 1.0_DP - VSq ) * ( 1.0_DP / Gm_dd_11                      &
          * ( 1.0_DP - VSq * Cs**2 ) - V1**2 * ( 1.0_DP - Cs**2 )))) - Beta1, &

          Alpha * V1 - Beta1, &

          Alpha / ( 1.0_DP - VSq * Cs**2 ) * ( V1 * ( 1.0_DP - Cs**2 ) + Cs   &
          * SQRT( ( 1.0_DP - VSq ) * ( 1.0_DP / Gm_dd_11                      &
          * ( 1.0_DP - VSq * Cs**2 ) - V1**2 * ( 1.0_DP - Cs**2 )))) - Beta1, &

          Alpha * V1 - Beta1, &

          Alpha * V1 - Beta1 ]
    END DO

  END SUBROUTINE Eigenvalues_GR


  SUBROUTINE ComputeSoundSpeed_GR( p, e, rho, Gamma )

    REAL(DP), DIMENSION(:), INTENT(in) :: p, e, rho, Gamma
    REAL(DP), DIMENSION( 1 : nX( 1 ) ) :: Cs

    Cs = SQRT ( Gamma * p / ( rho + e + p ) )

  END SUBROUTINE ComputeSoundSpeed_GR


  FUNCTION Flux_X1_GR &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Alpha, Beta1 )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, P, Ne
    REAL(DP), INTENT(in)       :: Alpha, Beta1, Gm11, Gm22, Gm33
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1_GR

    REAL(DP) :: vSq, W, h

    vSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = 1.0_DP / SQRT( 1.0_DP - vSq )

    h   = 1.0_DP + ( E + P ) / D

    Flux_X1_GR(iCF_D)  &
      = W * D * ( Alpha * V1 - Beta1 )

    Flux_X1_GR(iCF_S1) &
      = D * h * W**2 * Gm11 * V1 * ( Alpha * V1 - Beta1 ) + Alpha * P

    Flux_X1_GR(iCF_S2) &
      = D * h * W**2 * Gm22 * V2 * ( Alpha * V1 - Beta1 )

    Flux_X1_GR(iCF_S3) &
      = D * h * W**2 * Gm33 * V3 * ( Alpha * V1 - Beta1 )

    Flux_X1_GR(iCF_E)  &
      = ( D * h * W**2 - D * W ) * ( Alpha * V1 - Beta1 ) + Beta1 * P

    Flux_X1_GR(iCF_Ne) &
      = W * Ne * ( Alpha * V1 - Beta1 )

    RETURN
  END FUNCTION Flux_X1_GR


  REAL(DP) FUNCTION AlphaC_GR( U_L, U_R, F_L, F_R, aP, aM, &
                               V1_L, V1_R, Beta_u_1, Gm_dd_11 )

    ! --- Middle Wavespeed as suggested by Mignone and Bodo (2005) ---

    REAL(DP), DIMENSION(1:nCF), INTENT(in) :: U_L, U_R, F_L, F_R
    REAL(DP),                   INTENT(in) :: aP, aM, V1_L, V1_R, &
                                              Beta_u_1, Gm_dd_11
    REAL(DP), DIMENSION(1:nCF)             :: U, F, U_LL, U_RR, F_LL, F_RR
    REAL(DP)                               :: A, B, C, eps

    eps = TINY( 1.0_DP )

    ! --- Make sure that tau -> E for conserved variables and fluxes

    U_LL = U_L
    U_RR = U_R

    F_LL = F_L
    F_RR = F_R
    
    ! --- E = tau + D
    U_LL( iCF_E ) = U_L( iCF_E ) + U_L( iCF_D )
    U_RR( iCF_E ) = U_R( iCF_E ) + U_R( iCF_D )

    ! --- F_E = F_tau + D * ( V1 - Beta1 )
    F_LL( iCF_E ) = F_L( iCF_E ) + U_L( iCF_D ) * ( V1_L - Beta_u_1 )
    F_RR( iCF_E ) = F_R( iCF_E ) + U_R( iCF_D ) * ( V1_R - Beta_u_1 )

    ! --- Calculate the HLL conserved variable vector and flux
    ! --- Mignone & Bodo (2005) (Note the sign change on aM which is due
    ! --- to it being read in as positive but the formulae assuming
    ! --- it is negative)

    U = aP * U_RR + aM * U_LL + F_LL - F_RR
    F = aP * F_LL + aM * F_RR - aP * aM * (U_RR - U_LL )

    A = Gm_dd_11**2 * ( F( iCF_E ) + Beta_u_1 * U( iCF_E ) )
    B = -Gm_dd_11 * ( U( iCF_E ) + F( iCF_S1 ) + Beta_u_1 * U( iCF_S1 ) )
    C = U( iCF_S1 )

    ! --- Accounting for special cases of the solution to a
    ! --- quadratic equation when A = 0

    IF     ( ( ABS( A ) .LT. eps ) .AND. ( ABS( B ) .LT. eps ) &
            .AND. ( ABS( C ) .LT. eps ) )THEN
      WRITE(*,*) 'AlphaC is undefined'
      AlphaC_GR = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( B ) .LT. eps ) )THEN
      WRITE(*,*) 'AlphaC is undefined'
      WRITE(*,*) 'C:', C
      AlphaC_GR = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( C ) .LT. eps ) )THEN
      AlphaC_GR = 0.0_DP
    ELSE IF( ABS( A ) .LT. eps )THEN
      AlphaC_GR = -C / B
    ELSE
      AlphaC_GR = ( -B - SQRT( B**2 - 4.0_DP * A * C ) ) / ( 2.0_DP * A )
    END IF

    RETURN
    
  END FUNCTION AlphaC_GR


  PURE FUNCTION NumericalFlux_X1_HLLC_GR &
      ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, alpha_C, nF, &
        V1_L, V1_R, p_L, p_R, Beta_u_1, Gm_dd_11 )

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: alpha, alpha_P, alpha_M,      &
                                              alpha_C, V1_L, V1_R,          &
                                              p_L, p_R, Beta_u_1, Gm_dd_11
    REAL(DP)                               :: p, D, S1, S2, S3, E, Ne
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_HLLC_GR

    IF( alpha_M .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_L

    ELSEIF( alpha_P .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_R

    ELSE

      ! --- Note the sign change on alpha_M which is due to it being
      ! --- read in as positive but the formulae assuming it is negative

      IF( alpha_C .GE. 0.0_DP )THEN    

        ! -- UL_star

        p  = ( Gm_dd_11 * alpha_C * ( u_L( iCF_E ) + u_L( iCF_D ) ) &
             * ( -alpha_M + Beta_u_1 ) - u_L( iCF_S1 ) * ( alpha_C  &
             - alpha_M - V1_L + Beta_u_1 ) + p_L ) / ( 1.0_DP       &
             - Gm_dd_11 * alpha_C * ( -alpha_M + Beta_u_1 ) )

        D  = u_L( iCF_D  ) *  ( -alpha_M - V1_L    + Beta_u_1 ) &
                           /  ( -alpha_M - alpha_C + Beta_u_1 )

        S1 = u_L( iCF_S1 ) *  ( -alpha_M - V1_L    + Beta_u_1 ) &
                           /  ( -alpha_M - alpha_C + Beta_u_1 ) &
             + ( p - p_L ) /  ( -alpha_M - alpha_C + Beta_u_1 )

        S2 = u_L( iCF_S2 ) *  ( -alpha_M - V1_L    + Beta_u_1 ) &
                           /  ( -alpha_M - alpha_C + Beta_u_1 )

        S3 = u_L( iCF_S3 ) *  ( -alpha_M - V1_L    + Beta_u_1 ) &
                           /  ( -alpha_M - alpha_C + Beta_u_1 )

        E  = ( ( u_L( iCF_E ) + u_L( iCF_D ) ) * ( -alpha_M + Beta_u_1 )    &
             + 1.0_DP / Gm_dd_11 * S1 - 1.0_DP / Gm_dd_11 * u_L( iCF_S1 ) ) &
             / ( -alpha_M + Beta_u_1 )

        Ne = u_L( iCF_Ne ) *  ( -alpha_M - V1_L    + Beta_u_1 ) &
                           /  ( -alpha_M - alpha_C + Beta_u_1 )

      ELSE

        ! -- UR_star

        p  = ( Gm_dd_11 * alpha_C * ( u_R( iCF_E ) + u_R( iCF_D ) ) &
             * ( alpha_P + Beta_u_1 ) - u_R( iCF_S1 ) * ( alpha_C   &
             + alpha_P - V1_R + Beta_u_1 ) + p_R ) / ( 1.0_DP       &
             - Gm_dd_11 * alpha_C * ( alpha_P + Beta_u_1 ) )

        D  = u_R( iCF_D  ) *  ( alpha_P - V1_R    + Beta_u_1 ) &
                           /  ( alpha_P - alpha_C + Beta_u_1 )

        S1 = u_R( iCF_S1 ) *  ( alpha_P - V1_R    + Beta_u_1 ) &
                           /  ( alpha_P - alpha_C + Beta_u_1 ) &
             + ( p - p_R ) /  ( alpha_P - alpha_C + Beta_u_1 )

        S2 = u_R( iCF_S2 ) *  ( alpha_P - V1_R    + Beta_u_1 ) &
                           /  ( alpha_P - alpha_C + Beta_u_1 )

        S3 = u_R( iCF_S3 ) *  ( alpha_P - V1_R    + Beta_u_1 ) &
                           /  ( alpha_P - alpha_C + Beta_u_1 )

        E  = ( ( u_R( iCF_E ) + u_R( iCF_D ) ) * ( alpha_P + Beta_u_1 )     &
             + 1.0_DP / Gm_dd_11 * S1 - 1.0_DP / Gm_dd_11 * u_R( iCF_S1 ) ) &
             / ( alpha_P + Beta_u_1 )

        Ne  = u_R( iCF_Ne ) *  ( alpha_P - V1_R    + Beta_u_1 ) &
                            /  ( alpha_P - alpha_C + Beta_u_1 )

      END IF

      NumericalFlux_X1_HLLC_GR( iCF_D  ) &
        = D  * ( alpha_C - Beta_u_1 )
      NumericalFlux_X1_HLLC_GR( iCF_S1 ) &
        = S1 * ( alpha_C - Beta_u_1 ) + p
      NumericalFlux_X1_HLLC_GR( iCF_S2 ) &
        = S2 * ( alpha_C - Beta_u_1 )
      NumericalFlux_X1_HLLC_GR( iCF_S3 ) &
        = S3 * ( alpha_C - Beta_u_1 )
      NumericalFlux_X1_HLLC_GR( iCF_E  ) &
        = 1.0_DP / Gm_dd_11 * S1 - Beta_u_1 * ( E - D ) - alpha_C * D
      NumericalFlux_X1_HLLC_GR( iCF_Ne ) &
        = Ne * ( alpha_C - Beta_u_1 )

    END IF

    RETURN

  END FUNCTION NumericalFlux_X1_HLLC_GR


END MODULE EulerEquationsUtilitiesModule_Beta_GR
