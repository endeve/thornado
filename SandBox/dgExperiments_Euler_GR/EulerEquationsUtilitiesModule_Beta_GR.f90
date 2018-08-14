MODULE EulerEquationsUtilitiesModule_Beta_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_Cs
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR
  
  IMPLICIT NONE
  PRIVATE :: ComputeFunJacP

  PUBLIC :: ComputeFromConserved
  PUBLIC :: ComputePrimitive_GR
  PUBLIC :: ComputeConserved_GR
  PUBLIC :: Eigenvalues_GR
  PUBLIC :: Flux_X1_GR
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: AlphaC_GR
  PUBLIC :: NumericalFlux_X1_LLF_GR
  PUBLIC :: NumericalFlux_X1_HLL_GR
  PUBLIC :: NumericalFlux_X1_HLLC_GR

  LOGICAL, PARAMETER :: DEBUG = .FALSE.
  INTEGER, PARAMETER :: MAX_IT = 100

CONTAINS


  SUBROUTINE ComputeFromConserved( iX_B0, iX_E0, G, U, P, A )

    INTEGER, INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout)  :: &
      P(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
      A(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    INTEGER :: iX1, iX2, iX3

    ! --- Update primitive variables, pressure, and sound speed
    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputePrimitive_GR &
                 ( U(1:,iX1,iX2,iX3,iCF_D),        &
                   U(1:,iX1,iX2,iX3,iCF_S1),       &
                   U(1:,iX1,iX2,iX3,iCF_S2),       &
                   U(1:,iX1,iX2,iX3,iCF_S3),       &
                   U(1:,iX1,iX2,iX3,iCF_E),        &
                   U(1:,iX1,iX2,iX3,iCF_Ne),       &
                   P(1:,iX1,iX2,iX3,iPF_D),        &
                   P(1:,iX1,iX2,iX3,iPF_V1),       &
                   P(1:,iX1,iX2,iX3,iPF_V2),       &
                   P(1:,iX1,iX2,iX3,iPF_V3),       &
                   P(1:,iX1,iX2,iX3,iPF_E),        &
                   P(1:,iX1,iX2,iX3,iPF_Ne),       &
                   A(1:,iX1,iX2,iX3,iAF_P),        &
                   G(1:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   G(1:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   G(1:,iX1,iX2,iX3,iGF_Gm_dd_33) )

          CALL ComputeSoundSpeedFromPrimitive_GR &
                 ( P(1:,iX1,iX2,iX3,iPF_D), P(1:,iX1,iX2,iX3,iPF_E), &
                     P(1:,iX1,iX2,iX3,iPF_Ne), A(1:,iX1,iX2,iX3,iAF_Cs) )

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeFromConserved


  SUBROUTINE ComputePrimitive_GR &
              ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                AF_P, GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: CF_D, CF_S1, CF_S2, CF_S3, &
                                           CF_E, CF_Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: PF_D, PF_V1, PF_V2, PF_V3, &
                                           PF_E, PF_Ne, AF_P
    REAL(DP), DIMENSION(:), INTENT(in)  :: GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33

    LOGICAL :: CONVERGED

    INTEGER :: i, ITERATION, nNodes

    REAL(DP) :: SSq, Pold, vSq, W, h, Pnew, q(SIZE(CF_D))

    REAL(DP) :: FunP, JacP
    REAL(DP), PARAMETER :: TolP = 1.0d-8, TolFunP = 1.0d-12

    q = CF_E(:) + CF_D(:) - SQRT( CF_D(:)**2 &
                                + CF_S1(:)**2 / GF_Gm_dd_11(:)  &
                                + CF_S2(:)**2 / GF_Gm_dd_22(:)  &
                                + CF_S3(:)**2 / GF_Gm_dd_33(:) )

    nNodes = SIZE( CF_D )

    ! --- Loop through all the nodes ---
    DO i = 1, nNodes

      !IF( DEBUG )THEN
         IF( q(i) .LT. Zero )THEN
           WRITE(*,'(A6,ES18.10E3)') 'q(i): ', q(i)
!           WRITE(*,'(A6,ES18.10E3)') 'Gm11: ', GF_Gm_dd_11(:)
!           WRITE(*,'(A6,ES18.10E3)') 'Gm22: ', GF_Gm_dd_22(:)
!           WRITE(*,'(A6,ES18.10E3)') 'Gm33: ', GF_Gm_dd_33(:)
         END IF
      !END IF
    
      CONVERGED = .FALSE.
      ITERATION = 0

      SSq = CF_S1(i)**2 / GF_Gm_dd_11(i) &
            + CF_S2(i)**2 / GF_Gm_dd_22(i) &
            + CF_S3(i)**2 / GF_Gm_dd_33(i)

      ! --- Find Pressure with Newton's Method ---

      ! --- Approximation for pressure assuming h^2~=1 ---
!      Pold = MAX( SQRT( SSq + CF_D(i)**2 ) - CF_D(i) - CF_E(i), SqrtTiny )
      Pold = MAX( -q(i), SqrtTiny )

      DO WHILE ( .NOT. CONVERGED )

        ITERATION = ITERATION + 1

        CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pold, FunP, JacP )

        Pnew = Pold - FunP / JacP

        ! --- Check if Newton's method has converged ---
        IF( ABS( FunP / JacP ) / ABS( Pnew ) .LT. TolP )THEN
!        IF( ABS( Pnew - Pold ) / ABS( Pnew ) .LT. TolP )THEN
          CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pnew, FunP, JacP )
          IF( ABS( FunP ) .LT. TolFunP ) THEN
            CONVERGED = .TRUE.
          ELSE
            WRITE(*,'(A)') 'No convergence...'
            WRITE(*,'(A7,ES24.16E3)') 'Pold:  ', Pold
            WRITE(*,'(A7,ES24.16E3)') 'Pnew:  ', Pnew
            WRITE(*,'(A7,ES24.16E3)') 'FunP:  ', FunP
            WRITE(*,'(A)') 'Stopping...'
            STOP
          END IF
        END IF

        ! --- STOP after MAX_IT iterations ---
        IF( ITERATION .GE. MAX_IT )THEN
           
          WRITE(*,*) 'PF_D:   ', PF_D(:)
          WRITE(*,*) 'PF_V1:  ', PF_V1(:)
          WRITE(*,*) 'PF_V2:  ', PF_V2(:)
          WRITE(*,*) 'PF_V3:  ', PF_V3(:)
          WRITE(*,*) 'PF_E:   ', PF_E(:)

          WRITE(*,*) 'CF_D:   ', CF_D(:)
          WRITE(*,*) 'CF_S1:  ', CF_S1(:)
          WRITE(*,*) 'CF_S2:  ', CF_S2(:)
          WRITE(*,*) 'CF_S3:  ', CF_S3(:)
          WRITE(*,*) 'CF_E:   ', CF_E(:)

          WRITE(*,*) 'Gm11:   ', GF_Gm_dd_11(:)
          WRITE(*,*) 'Gm22:   ', GF_Gm_dd_22(:)
          WRITE(*,*) 'Gm33:   ', GF_Gm_dd_33(:)
           
          WRITE(*,'(A,ES24.16E3)') 'Pold:  ', Pold
          WRITE(*,'(A,ES24.16E3)') 'Pnew:  ', Pnew
          WRITE(*,'(A,ES24.16E3)') 'FunP:  ', FunP
          WRITE(*,'(A,ES24.16E3)') &
                  '|Pnew-Pold|/|Pnew|: ', ABS( Pnew - Pold ) / ABS( Pnew )
          IF( ITERATION == MAX_IT + 3 )THEN
            WRITE(*,'(A)') &
                   'Max allowed iterations reached, no convergence. Stopping...'
            STOP
          END IF
        END IF

        Pold = Pnew

      END DO

      AF_P(i) = Pnew

      IF( AF_P(i) .LT. Zero ) WRITE(*,'(A3,ES24.10E3)') 'P: ', AF_P(i)
      vSq = SSq / ( CF_E(i) + AF_P(i) + CF_D(i) )**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = ( CF_E(i) + AF_P(i) + CF_D(i) ) / ( W * CF_D(i) )

      !IF( DEBUG )THEN
        IF( h .LT. 1.0_DP ) WRITE(*,'(A6,ES18.10E3)') 'h:    ', h
      !END IF

      ! --- Recover Primitive Variables ---

      PF_D(i)  = CF_D(i) / W
      
      PF_V1(i) = CF_S1(i) / ( CF_D(i) * h * W * GF_Gm_dd_11(i) )

      PF_V2(i) = CF_S2(i) / ( CF_D(i) * h * W * GF_Gm_dd_22(i) )

      PF_V3(i) = CF_S3(i) / ( CF_D(i) * h * W * GF_Gm_dd_33(i) )

      PF_E(i)  = CF_D(i) * ( h - 1.0_DP ) / W - AF_P(i)
      
      PF_Ne(i) = CF_Ne(i) / W

   END DO

  END SUBROUTINE ComputePrimitive_GR


  SUBROUTINE ComputeFunJacP( D, E, SSq, P, FunP, JacP )

    REAL(DP), INTENT(in)  :: D, E, SSq, P
    REAL(DP), INTENT(out) :: FunP, JacP

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS
    REAL(DP), DIMENSION(1) :: Pbar

    HSq = ( E + P + D )**2

    RHO = D * SQRT( HSq - SSq ) / SQRT( HSq )

    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - D ) / D

    !EPS = MAX( EPS, SqrtTiny )
    CALL ComputePressureFromSpecificInternalEnergy &
         ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    FunP = P - Pbar(1)
!    FunP = 1.0_DP - Pbar(1) / P
    dRHO = D * SSq / ( SQRT( HSq - SSq ) * HSq )
    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )

    JacP = 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS )
!    JacP = Pbar(1) / P * ( 1.0_DP / P - dRHO / RHO - dEPS / EPS )

  END SUBROUTINE ComputeFunJacP

  
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

    REAL(DP), DIMENSION( 1 : SIZE(PF_D) ) :: vSq, W, h

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
  

  PURE FUNCTION Eigenvalues_GR &
    ( V_1, V_2, V_3, V_i, Cs, Gm_11, Gm_22, Gm_33, Gm_ii, Alpha, Beta_i )

    ! Alpha is the lapse function
    ! V_i is the contravariant component V^i
    ! Beta_1 is the contravariant component Beta^1

    REAL(DP)             :: Eigenvalues_GR(1:nCF)
    REAL(DP), INTENT(in) :: V_1, V_2, V_3, V_i, Cs
    REAL(DP), INTENT(in) :: Gm_11, Gm_22, Gm_33, Gm_ii, Alpha, Beta_i

    REAL(DP) :: VSq

    VSq = Gm_11 * V_1**2 + Gm_22 * V_2**2 + Gm_33 * V_3**2

    Eigenvalues_GR(1) &
      = Alpha / ( One - VSq * Cs**2 ) * ( V_i * ( One - Cs**2 ) &
        - Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gm_ii &
           - V_i**2 * ( One - Cs**2 ) ) ) ) - Beta_i

    Eigenvalues_GR(2) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(3) &
      = Alpha / ( One - VSq * Cs**2 ) * ( V_i * ( One - Cs**2 ) &
        + Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gm_ii &
           - V_i**2 * ( One - Cs**2 ) ) ) ) - Beta_i

    Eigenvalues_GR(4) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(5) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(6) &
      = Alpha * V_i - Beta_i

    RETURN
  END FUNCTION Eigenvalues_GR


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
      = D * W * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_S1) &
      = D * h * W**2 * Gm11 * V1 * ( V1 - Beta1 / Alpha ) + P

    Flux_X1_GR(iCF_S2) &
      = D * h * W**2 * Gm22 * V2 * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_S3) &
      = D * h * W**2 * Gm33 * V3 * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_E)  &
      = D * W * ( h * W - One ) * ( V1 - Beta1 / Alpha ) + Beta1 / Alpha * P

    Flux_X1_GR(iCF_Ne) &
      = Ne * W * ( V1 - Beta1 / Alpha )

    RETURN
  END FUNCTION Flux_X1_GR


  PURE FUNCTION StressTensor_Diagonal( S_1, S_2, S_3, V_1, V_2, V_3, P )

    REAL(DP)             :: StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: S_1, S_2, S_3, V_1, V_2, V_3, P

    StressTensor_Diagonal(1) = S_1 * V_1 + P
    StressTensor_Diagonal(2) = S_2 * V_2 + P
    StressTensor_Diagonal(3) = S_3 * V_3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal


  REAL(DP) FUNCTION AlphaC_GR &
             ( U_L, F_L, U_R, F_R, aP, aM, Gm_dd_11, Beta_u_1, LapseFunction )

    ! --- Middle Wavespeed as suggested by Mignone and Bodo (2005) ---

    REAL(DP), INTENT(in) :: U_L(nCF), F_L(nCF), U_R(nCF), F_R(nCF), &
                            aP, aM, Gm_dd_11, Beta_u_1, LapseFunction
    REAL(DP)             :: U_S1, F_S1, U_E, F_E, A, B, C, eps

    eps = SqrtTiny

    ! --- Note the sign change on aM which is due
    !     to it being read in as positive but Mignone assuming
    !     it is negative ---

    U_S1 = aP * U_R(iCF_S1) + aM * U_L(iCF_S1) + F_L(iCF_S1) - F_R(iCF_S1)

    U_E  = aP * ( U_R(iCF_E) + U_R(iCF_D) ) + aM * ( U_L(iCF_E) + U_L(iCF_D) ) &
          + ( F_L(iCF_E) + F_L(iCF_D) ) - ( F_R(iCF_E) + F_R(iCF_D) )

    F_S1 =  aP * F_L(iCF_S1) + aM * F_R(iCF_S1) &
          - aP * aM * ( U_R(iCF_S1) - U_L(iCF_S1 ) )

    F_E  =  aP * ( F_L(iCF_E) + F_L(iCF_D) ) &
          + aM * ( F_R(iCF_E) + F_R(iCF_D) ) &
          - aP * aM * ( ( U_R(iCF_E) + U_R(iCF_D) ) &
          - ( U_L(iCF_E) + U_L(iCF_D) ) )

    ! --- A, B, and C from quadratic equation ---
    A = Gm_dd_11**2 * ( F_E + Beta_u_1 * U_E )
    B = -Gm_dd_11 * ( LapseFunction * U_E + F_S1 + Beta_u_1 * U_S1 )
    C = LapseFunction * U_S1

    ! --- Accounting for special cases of the solution to a
    !     quadratic equation when A = 0 ---

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


  PURE FUNCTION NumericalFlux_X1_LLF_GR &
      ( u_L, u_R, Flux_L, Flux_R, alpha_P, alpha_M, alpha_C, nF, &
        V1_L, V1_R, p_L, p_R, Beta_u_1, Gm_dd_11 )    

    ! --- Local Lax-Friedrichs Flux ---

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: alpha_P, alpha_M,      &
                                              alpha_C, V1_L, V1_R,          &
                                              p_L, p_R, Beta_u_1, Gm_dd_11
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_LLF_GR(1:nF)
    REAL(DP) :: alpha

    alpha    = MAX( alpha_m, alpha_p )

    NumericalFlux_X1_LLF_GR &
      = 0.5_DP * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_X1_LLF_GR


  PURE FUNCTION NumericalFlux_X1_HLL_GR &
      ( u_L, u_R, Flux_L, Flux_R, alpha_P, alpha_M, alpha_C, nF, &
        V1_L, V1_R, p_L, p_R, Beta_u_1, Gm_dd_11 )

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: alpha_P, alpha_M,    &
                                              alpha_C, V1_L, V1_R, &
                                              p_L, p_R, Beta_u_1, Gm_dd_11
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_HLL_GR(1:nF)

    NumericalFlux_X1_HLL_GR &
      = ( alpha_P * flux_L + alpha_M * flux_R &
            - alpha_P * alpha_M * ( u_R - u_L ) ) / ( alpha_P + alpha_M )

  END FUNCTION NumericalFlux_X1_HLL_GR


  PURE FUNCTION NumericalFlux_X1_HLLC_GR &
      ( u_L, u_R, Flux_L, Flux_R, alpha_P, alpha_M, alpha_C, nF, &
        V1_L, V1_R, p_L, p_R, Beta_u_1, Gm_dd_11 )

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: alpha_P, alpha_M,      &
                                              alpha_C, V1_L, V1_R,          &
                                              p_L, p_R, Beta_u_1, Gm_dd_11
    REAL(DP)                               :: p, D, S1, S2, S3, E, Ne
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_HLLC_GR

    IF( alpha_M .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_L

    ELSEIF( alpha_P .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_R

    ELSE

      ! --- From Mignone & Bodo (2005)
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
