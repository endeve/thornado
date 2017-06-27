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
  PUBLIC :: Eigenvalues
  PUBLIC :: ComputeSoundSpeed
  PUBLIC :: Flux_X1
  PUBLIC :: AlphaC

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

      SSq =   Gm11 * uCF(i,iCF_S1)**2 &
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

        IF( ABS( Pnew / Pold - 1.0_DP ) <= TolP ) Converged = .TRUE.

        IF( 1 == 1 )THEN
          IF( nIter .GT. 10 )THEN
          WRITE(*,'(A6,1x,I2,3x,A10,E27.20)') 'nIter:', &
                nIter,'|ERROR| =',ABS( Pnew / Pold - 1.0_DP )
          END IF
        END IF
        
        IF( nIter == 100)THEN
          WRITE(*,*) 'No convergence, |ERROR|:', ABS( Pnew / Pold - 1.0_DP )
          WRITE(*,*) 'Pold:                   ', Pold
          WRITE(*,*) 'Pnew:                   ', Pnew
          STOP
        END IF

        IF( ISNAN( Pnew ) )THEN
          WRITE(*,*) 'nIter:              ', nIter
          WRITE(*,*) 'Pold:               ', Pold
          WRITE(*,*) 'Pnew:               ', Pnew
          WRITE(*,*) 'D:                  ', uCF(i,iCF_D)
          WRITE(*,*) 'tau+D-sqrt(D^2+S^2):', &
            uCF( i, iCF_E )+uCF( i, iCF_D )-SQRT( uCF( i, iCF_D )**2 + SSq )
          STOP
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

      uPF(i,iPF_V3) = uCF(i,iCF_S3) / ( uCF(i,iCF_D) * h * W * Gm33 )

      uPF(i,iPF_E)  = uCF(i,iCF_D) * ( h - 1.0_DP ) / W - Pnew

      uPF(i,iPF_Ne) = uCF(i,iCF_Ne) / W

    END DO

  END SUBROUTINE ComputePrimitive


  PURE FUNCTION Eigenvalues( V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                              Gm_uu_11, Alpha, Beta1, Cs )

    ! Alpha is the lapse function
    ! Vi is the contravariant component V^i
    ! Beta1 is the contravariant component Beta^1

    REAL(DP), INTENT(in) :: V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                              Gm_uu_11, Alpha, Beta1, Cs
    REAL(DP) :: VSq
    REAL(DP), DIMENSION(1:5) :: Eigenvalues

    Vsq = Gm_dd_11 * V1**2 + Gm_dd_22 * V2**2 + Gm_dd_33 * V3**2

    Eigenvalues(1:5) = &
      [ Alpha / ( 1.0_DP - VSq * Cs**2 ) * ( V1 * ( 1.0_DP - Cs**2 ) - Cs &
        * SQRT( ( 1.0_DP - VSq ) * ( Gm_uu_11 * ( 1.0_DP - VSq * Cs**2 )  &
        - V1**2 * ( 1.0_DP - Cs**2 )))) - Beta1, &

        Alpha * V1 - Beta1, &

        Alpha / ( 1.0_DP - VSq * Cs**2 ) * ( V1 * ( 1.0_DP - Cs**2 ) + Cs &
        * SQRT( ( 1.0_DP - VSq ) * ( Gm_uu_11 * ( 1.0_DP - VSq * Cs**2 )  &
        - V1**2 * ( 1.0_DP - Cs**2 )))) - Beta1, &

        Alpha * V1 - Beta1, &

        Alpha * V1 - Beta1 ]

    RETURN
  END FUNCTION Eigenvalues


  PURE FUNCTION ComputeSoundSpeed( p, e, rho, Gamma ) RESULT( Cs )

    REAL(DP), INTENT(in) :: p, e, rho, Gamma
    REAL(DP) :: Cs

    Cs = SQRT ( Gamma * p / ( rho + e + p ) )

  END FUNCTION ComputeSoundSpeed


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


  REAL(DP) FUNCTION AlphaC( U_L, U_R, F_L, F_R, aP, aM, &
                            V1_L, V1_R, Beta_u_1, Gm_dd_11 )

    ! --- Middle Wavespeed as Suggested by Mignone and Bodo (2005) ---

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
      AlphaC = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( B ) .LT. eps ) )THEN
      WRITE(*,*) 'AlphaC is undefined'
      WRITE(*,*) 'C:', C
      AlphaC = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( C ) .LT. eps ) )THEN
      AlphaC = 0.0_DP
    ELSE IF( ABS( A ) .LT. eps )THEN
      AlphaC = 0.0_DP
    ELSE
      AlphaC = ( -B - SQRT( B**2 - 4.0_DP * A * C ) ) / ( 2.0_DP * A )
    END IF

    RETURN
    
  END FUNCTION AlphaC


END MODULE EulerEquationsUtilitiesModule_GR

