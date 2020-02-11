MODULE TwoMoment_UtilitiesModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, Five
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor_Relativistic, &
    EddingtonFactor, &
    HeatFluxFactor


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeEddingtonTensorComponents_dd
  PUBLIC :: ComputeEddingtonTensorComponents_ud
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3

CONTAINS

  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, nIterations_Option )

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    REAL(DP) :: I_d_1, I_d_2, I_d_3, A_d_1, A_d_2, A_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)
    REAL(DP) :: W, DT    

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3
    
    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N / W
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )


    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      CALL ComputeEddingtonTensorComponents_dd &
             ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
               alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3   )

      A_d_1 = V_u_1 * k_dd_11 + V_u_2 * k_dd_12 + V_u_3 * k_dd_13
      A_d_2 = V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_23
      A_d_3 = V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33

      DET = W**2 - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 )

      LMAT(1,1) = W
      LMAT(2,1) = - A_d_1
      LMAT(3,1) = - A_d_2
      LMAT(4,1) = - A_d_3

      LMAT(1,2) = - V_u_1
      LMAT(2,2) = W - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 ) / W
      LMAT(3,2) = V_u_1 * A_d_2 / W
      LMAT(4,2) = V_u_1 * A_d_3 / W

      LMAT(1,3) = - V_u_2
      LMAT(2,3) = V_u_2 * A_d_1 / W
      LMAT(3,3) = W - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 ) / W
      LMAT(4,3) = V_u_2 * A_d_3 / W

      LMAT(1,4) = - V_u_3
      LMAT(2,4) = V_u_3 * A_d_1 / W
      LMAT(3,4) = V_u_3 * A_d_2 / W
      LMAT(4,4) = W - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 ) / W

      LMAT = LMAT / DET

      CALL DGEMV( 'N', 4, 4, One, LMAT, 4, CVEC, 1, Zero, GVEC(:,mk), 1 )

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        BVEC = - FVEC(:,mk)

        AMAT(:,1:mk-1) &
          = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 4, mk-1, 1, AMAT(:,1:mk-1), 4, BVEC, 4, &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVEC(1:mk-1)
        Alpha(mk)     = One - SUM( Alpha(1:mk-1) )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF
      !this may be the problem converting abck and forth from upper to lower and back
      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = ( 1.0_DP - B_d_1 * V_u_1 / alp ) * I_d_1 / Gm_dd_11  &
            - I_d_2 * B_d_1 * V_u_2 / ( alp *Gm_dd_11 ) - I_d_3 * B_d_1 * V_u_3 / ( Gm_dd_11 * alp )
      I_d_2 = UVEC(3); I_u_2 = ( 1.0_DP - B_d_2 * V_u_2 / alp ) * I_d_2 / Gm_dd_22  &
            - I_d_1 * B_d_2 * V_u_1 / ( alp *Gm_dd_22 ) - I_d_3 * B_d_2 * V_u_3 / ( Gm_dd_22 * alp )
      I_d_3 = UVEC(4); I_u_3 = ( 1.0_DP - B_d_3 * V_u_3 / alp ) * I_d_3 / Gm_dd_33  &
            - I_d_1 * B_d_3 * V_u_1 / ( alp *Gm_dd_33 ) - I_d_2 * B_d_3 * V_u_2 / ( Gm_dd_33 * alp )

    END DO

    IF( PRESENT( nIterations_Option ) )THEN

      nIterations_Option = k

    END IF

    IF( k == MaxIterations )THEN

      PRINT*
      PRINT*, "ComputePrimitive_TwoMoment"
      PRINT*
      PRINT*, "  N     = ", N
      PRINT*, "  G_d_1 = ", G_d_1
      PRINT*, "  G_d_2 = ", G_d_2
      PRINT*, "  G_d_3 = ", G_d_3
      PRINT*
      PRINT*, "  V_u_1 = ", V_u_1
      PRINT*, "  V_u_2 = ", V_u_2
      PRINT*, "  V_u_3 = ", V_u_3
      PRINT*

      PRINT*, "  Converged with k = ", k

      PRINT*
      PRINT*, "  FVECm = ", FVECm
      PRINT*

      PRINT*
      PRINT*, "  D     = ", D
      PRINT*, "  I_u_1 = ", I_u_1
      PRINT*, "  I_u_2 = ", I_u_2
      PRINT*, "  I_u_3 = ", I_u_3
      PRINT*

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment

  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_u_1, I_u_2, I_u_3, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: W  
    REAL(DP) :: DT 
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3 ! --- Index Up

    CALL ComputeEddingtonTensorComponents_dd &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
             alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )


   

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

 
    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    ! --- Conserved Number Density ---

    N = W * D + V_u_1 * I_d_1 &
              + V_u_2 * I_d_2 &
              + V_u_3 * I_d_3

    ! --- Conserved Number Flux Density (1) ---

    G_d_1 = W * I_d_1 &
                 + (   V_u_1 * k_dd_11 &
                     + V_u_2 * k_dd_12 &
                     + V_u_3 * k_dd_13 ) * D

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = W * I_d_2 &
                 + (   V_u_1 * k_dd_12 &
                     + V_u_2 * k_dd_22 &
                     + V_u_3 * k_dd_23 ) * D

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = W * I_d_3 &
                 + (   V_u_1 * k_dd_13 &
                     + V_u_2 * k_dd_23 &
                     + V_u_3 * k_dd_33 ) * D


  END SUBROUTINE ComputeConserved_TwoMoment

  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )

    EF = EddingtonFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    u_d_1 = B_d_1 * W / alp + Gm_dd_11 * ( V_u_1 - B_u_1 / alp )
    u_d_2 = B_d_2 * W / alp + Gm_dd_22 * ( V_u_2 - B_u_2 / alp )
    u_d_3 = B_d_3 * W / alp + Gm_dd_33 * ( V_u_3 - B_u_3 / alp )
 
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    ! --- Diagonal Eddington Tensor Components ---

    k_dd_11 = a * ( Gm_dd_11 + u_d_1 * u_d_1 ) + b * h_d_1 * h_d_1
    k_dd_22 = a * ( Gm_dd_22 + u_d_2 * u_d_2 ) + b * h_d_2 * h_d_2
    k_dd_33 = a * ( Gm_dd_33 + u_d_3 * u_d_3 ) + b * h_d_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_dd_12 = a * u_d_1 * u_d_2 + b * h_d_1 * h_d_2
    k_dd_13 = a * u_d_1 * u_d_3 + b * h_d_1 * h_d_3
    k_dd_23 = a * u_d_2 * u_d_3 + b * h_d_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_dd

  SUBROUTINE ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33

    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_u_1, h_u_2, h_u_3, h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, u_u_1, u_u_2, u_u_3, W

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )

    EF = EddingtonFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    u_d_1 = B_d_1 * W / alp + Gm_dd_11 * ( V_u_1 - B_u_1 / alp )
    u_d_2 = B_d_2 * W / alp + Gm_dd_22 * ( V_u_2 - B_u_2 / alp )
    u_d_3 = B_d_3 * W / alp + Gm_dd_33 * ( V_u_3 - B_u_3 / alp )

    u_u_1 = u_d_1 / Gm_dd_11 
    u_u_2 = u_d_2 / Gm_dd_22
    u_u_3 = u_d_3 / Gm_dd_33
 
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    ! --- Diagonal Eddington Tensor Components ---

    k_ud_11 = a * ( 1.0_DP + u_u_1 * u_d_1 ) + b * h_u_1 * h_d_1
    k_ud_22 = a * ( 1.0_DP + u_u_2 * u_d_2 ) + b * h_u_2 * h_d_2
    k_ud_33 = a * ( 1.0_DP + u_u_3 * u_d_3 ) + b * h_u_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_ud_12 = a * u_u_1 * u_d_2 + b * h_u_1 * h_d_2
    k_ud_13 = a * u_u_1 * u_d_3 + b * h_u_1 * h_d_3
    k_ud_23 = a * u_u_2 * u_d_3 + b * h_u_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_ud

  FUNCTION Flux_X1( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X1(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_11, k_ud_22, k_ud_33, k_ud_12, k_ud_13, k_ud_23
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )

    Flux_X1(1) = alp * I_u_1 + ( alp * V_u_1 - B_u_1 ) * W * D
    Flux_X1(2) = alp * k_ud_11 + (alp * V_u_1 - B_u_1 ) * W * I_d_1
    Flux_X1(3) = alp * k_ud_12 + (alp * V_u_1 - B_u_1 ) * W * I_d_2
    Flux_X1(4) = alp * k_ud_13 + (alp * V_u_1 - B_u_1 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X1

  FUNCTION Flux_X2( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X2(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_11, k_ud_22, k_ud_33, k_ud_12, k_ud_13, k_ud_23
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )

    Flux_X2(1) = alp * I_u_2 + ( alp * V_u_2 - B_u_2 ) * W * D
    Flux_X2(2) = alp * k_ud_12 + (alp * V_u_2 - B_u_2 ) * W * I_d_1
    Flux_X2(3) = alp * k_ud_22 + (alp * V_u_2 - B_u_2 ) * W * I_d_2
    Flux_X2(4) = alp * k_ud_23 + (alp * V_u_2 - B_u_2 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X2

  FUNCTION Flux_X3( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X3(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_11, k_ud_22, k_ud_33, k_ud_12, k_ud_13, k_ud_23
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3  )

    Flux_X3(1) = alp * I_u_3 + ( alp * V_u_3 - B_u_3 ) * W * D
    Flux_X3(2) = alp * k_ud_13 + (alp * V_u_3 - B_u_3 ) * W * I_d_1
    Flux_X3(3) = alp * k_ud_23 + (alp * V_u_3 - B_u_3 ) * W * I_d_2
    Flux_X3(4) = alp * k_ud_33 + (alp * V_u_3 - B_u_3 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X3



END MODULE TwoMoment_UtilitiesModule_Relativistic
