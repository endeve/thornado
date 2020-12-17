MODULE TwoMoment_UtilitiesModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, Five
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE
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
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
    USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor_Relativistic, &
    EddingtonFactor, &
    HeatFluxFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeFromConserved_TwoMoment
  PUBLIC :: ComputeEddingtonTensorComponents_dd
  PUBLIC :: ComputeEddingtonTensorComponents_uu
  PUBLIC :: ComputeEddingtonTensorComponents_ud
  PUBLIC :: ComputeHeatFluxTensorComponents_ddd_Lagrangian 
  PUBLIC :: ComputeHeatFluxTensorComponents_uud_Lagrangian
  PUBLIC :: ComputeHeatFluxTensorComponents_udd_Lagrangian
 ! PUBLIC :: ComputeStressEnergyComponents_Eulerian
  !PUBLIC :: ComputeHeatFluxTensorComponents_udd_Eulerian
  !PUBLIC :: ComputeXYZ
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: Flux_E
  PUBLIC :: Source_E
  !PUBLIC :: Flux_G
  PUBLIC :: NumericalFlux_LLF

CONTAINS

  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23, alp, B_u_1, B_u_2, B_u_3, nIterations_Option )

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23
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
    REAL(DP) :: k_dd_ij(1:3,1:3)
    REAL(DP) :: h_dd_11, h_dd_22, h_dd_33, h_dd_12, h_dd_13, h_dd_23, u_d_1, u_d_2, u_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)
    REAL(DP) :: W, DT, DTG
    
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

!    u_d_1 = B_d_1 * W / alp + Gm_dd_11 * ( V_u_1 - B_u_1 / alp )
!    u_d_2 = B_d_2 * W / alp + Gm_dd_22 * ( V_u_2 - B_u_2 / alp )
!    u_d_3 = B_d_3 * W / alp + Gm_dd_33 * ( V_u_3 - B_u_3 / alp )
!
!    h_dd_11 = Gm_dd_11 + u_d_1 * u_d_1 
!    h_dd_22 = Gm_dd_22 + u_d_2 * u_d_2 
!    h_dd_33 = Gm_dd_33 + u_d_3 * u_d_3 
!    h_dd_12 = u_d_1 * u_d_2 
!    h_dd_13 = u_d_1 * u_d_2 
!    h_dd_23 = u_d_1 * u_d_2 
!
!    DTG = ( -3.0_DP * W**2 + h_dd_11 * V_u_2**2 + 2.0_DP * h_dd_12 * V_u_1 * V_u_2 &
!        + 2.0_DP * h_dd_13 * V_u_1 * V_u_3 + 2.0_DP * h_dd_22 * V_u_2 * V_u_2 &
!        + 2.0_DP * h_dd_23 * V_u_2 * V_u_3 + 2.0_DP * h_dd_33 * V_u_3 * V_u_3  ) 
!
    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---
    D     = N
    I_u_1 = G_d_1
    I_u_2 = G_d_2
    I_u_3 = G_d_3
!    D     = ( -3.0_DP * N + 3.0_DP * V_u_1 * G_d_1 + 3.0_DP * V_u_2 * G_d_2 + 3.0_DP * V_u_3 * G_d_3 ) / DTG
!    I_d_1 = ( h_dd_11 * V_u_1 + h_dd_12 * V_u_2 + h_dd_13 * V_u_3 ) * N  &
!          + ( h_dd_22*V_u_2**2 + 2*h_dd_23*V_u_2*V_u_3 + V_u_1*h_dd_12*V_u_2 &
!          + h_dd_33*V_u_3**2 + V_u_1*h_dd_13*V_u_3 - 3*W**2 ) * G_d_1 &
!          -(V_u_2**2*h_dd_12 + V_u_1*V_u_2*h_dd_11 + V_u_2*V_u_3*h_dd_13) * G_d_2  &
!          -(V_u_3**2*h_dd_13 + V_u_1*V_u_3*h_dd_11 + V_u_2*V_u_3*h_dd_12) * G_d_3
!    I_d_1 = I_d_1 / ( W * DTG )
!    I_d_2 = (V_u_3*h_dd_23 + V_u_1*h_dd_12 + V_u_2*h_dd_22) * N  &
!          -(V_u_1**2*h_dd_12 + V_u_1*V_u_3*h_dd_23 + V_u_1*V_u_2*h_dd_22) * G_d_1 &
!          + (h_dd_11*V_u_1**2 + 2*h_dd_13*V_u_1*V_u_3 + V_u_2*h_dd_12*V_u_1 &
!          + h_dd_33*V_u_3**2 + V_u_2*h_dd_23*V_u_3 - 3*W**2) * G_d_2 &
!          -(V_u_3**2*h_dd_23 + V_u_1*V_u_3*h_dd_12 + V_u_2*V_u_3*h_dd_22) * G_d_3
!    I_d_2 = I_d_2 / ( W * DTG )
!    I_d_3 = (V_u_2*h_dd_23 + V_u_1*h_dd_13 + V_u_3*h_dd_33) * N &
!          -(V_u_1**2*h_dd_13 + V_u_1*V_u_2*h_dd_23 + V_u_1*V_u_3*h_dd_33) *G_d_1 &
!          -(V_u_2**2*h_dd_23 + V_u_1*V_u_2*h_dd_13 + V_u_2*V_u_3*h_dd_33) * G_d_2 &
!          + (h_dd_11*V_u_1**2 + 2*h_dd_12*V_u_1*V_u_2 + V_u_3*h_dd_13*V_u_1 &
!          + h_dd_22*V_u_2**2 + V_u_3*h_dd_23*V_u_2 - 3*W**2) * G_d_3 
!    I_d_3 = I_d_3 / ( W * DTG ) 
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
               alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij   )

      A_d_1 = V_u_1 * k_dd_ij(1,1) + V_u_2 * k_dd_ij(1,2) + V_u_3 * k_dd_ij(1,3)
      A_d_2 = V_u_1 * k_dd_ij(1,2) + V_u_2 * k_dd_ij(2,2) + V_u_3 * k_dd_ij(2,3)
      A_d_3 = V_u_1 * k_dd_ij(1,3) + V_u_2 * k_dd_ij(2,3) + V_u_3 * k_dd_ij(3,3)

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
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23, alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    REAL(DP) :: k_dd_ij(1:3,1:3)
    REAL(DP) :: W  
    REAL(DP) :: DT 
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3 ! --- Index Up

    CALL ComputeEddingtonTensorComponents_dd &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij  )



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
                 + (   V_u_1 * k_dd_ij(1,1) &
                     + V_u_2 * k_dd_ij(1,2) &
                     + V_u_3 * k_dd_ij(1,3) ) * D

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = W * I_d_2 &
                 + (   V_u_1 * k_dd_ij(1,2) &
                     + V_u_2 * k_dd_ij(2,2) &
                     + V_u_3 * k_dd_ij(2,3) ) * D

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = W * I_d_3 &
                 + (   V_u_1 * k_dd_ij(1,3) &
                     + V_u_2 * k_dd_ij(2,3) &
                     + V_u_3 * k_dd_ij(3,3) ) * D


  END SUBROUTINE ComputeConserved_TwoMoment


  SUBROUTINE ComputeFromConserved_TwoMoment &
     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF, CR, PR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      CF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(in)  :: &
      CR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nPR,1:nSpecies)

    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX
    REAL(DP) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nPF)

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( CF(iNodeX,iZ2,iZ3,iZ4,iCF_D ), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S1), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S2), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S3), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_E ), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_Ne), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_D ), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_E ), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_Ne), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        CALL ComputePrimitive_TwoMoment &
               ( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 Zero, Zero, Zero, &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha ), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_TwoMoment


  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij  )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      k_dd_ij(1:3,1:3)

    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    EF = EddingtonFactor( D, FF )
    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3


    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3
 
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
   
    k_dd_ij(1,1) = k_dd_11
 
    k_dd_ij(1,2) = k_dd_12
 
    k_dd_ij(1,3) = k_dd_13
 
    k_dd_ij(2,1) = k_dd_12
 
    k_dd_ij(2,2) = k_dd_22
 
    k_dd_ij(2,3) = k_dd_23
 
    k_dd_ij(3,1) = k_dd_13
 
    k_dd_ij(3,2) = k_dd_23
 
    k_dd_ij(3,3) = k_dd_33

 
  END SUBROUTINE ComputeEddingtonTensorComponents_dd

  SUBROUTINE ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


	    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
	    REAL(DP), INTENT(out) :: &
                k_uu_munu(0:3,0:3)


	    REAL(DP) :: FF, EF, a, b, Gm_uu_11, Gm_uu_22, Gm_uu_33
	    REAL(DP) :: h_u_1, h_u_2, h_u_3
	    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, V_0
	    REAL(DP) :: u_u_1, u_u_2, u_u_3, W
            REAL(DP) :: k_uu_11, k_uu_12, k_uu_13, k_uu_22, k_uu_23, k_uu_33, k_uu_10, k_uu_20, k_uu_30, k_uu_00

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    EF = EddingtonFactor( D, FF )
    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   

    Gm_uu_11 = 1.0_DP / Gm_dd_11 
    Gm_uu_22 = 1.0_DP / Gm_dd_22 
    Gm_uu_33 = 1.0_DP / Gm_dd_33 

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_2 - B_u_2 / alp ) 


    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )


    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )
   
    ! --- Diagonal Euuington Tensor Components ---

    k_uu_11 = a * ( Gm_uu_11 + u_u_1 * u_u_1 ) + b * h_u_1 * h_u_1
    k_uu_22 = a * ( Gm_uu_22 + u_u_2 * u_u_2 ) + b * h_u_2 * h_u_2
    k_uu_33 = a * ( Gm_uu_33 + u_u_3 * u_u_3 ) + b * h_u_3 * h_u_3
    ! --- Off-Diagonal Euuington Tensor Components ---
    k_uu_12 = a * u_u_1 * u_u_2 + b * h_u_1 * h_u_2
    k_uu_13 = a * u_u_1 * u_u_3 + b * h_u_1 * h_u_3
    k_uu_23 = a * u_u_2 * u_u_3 + b * h_u_2 * h_u_3


    k_uu_10 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_11 + V_d_2 * k_uu_12 + V_d_3 * k_uu_13 ) / alp

    k_uu_20 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_12 + V_d_2 * k_uu_22 + V_d_3 * k_uu_23 ) / alp

    k_uu_30 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_13 + V_d_2 * k_uu_23 + V_d_3 * k_uu_33 ) / alp
    
    k_uu_00 =  ( 1.0_DP/( 1.0_DP - V_0 / alp )**2 ) *( ( V_d_1**2 * k_uu_11 + V_d_2**2 * k_uu_22 + V_d_3**2 * k_uu_33 ) &
               + 2.0_DP * ( V_d_1 * V_d_2 * k_uu_12 + V_d_1 * V_d_3 * k_uu_13 + V_d_2 * V_d_3 * k_uu_23 ))/alp**2

    k_uu_munu(0,0) = k_uu_00
 
    k_uu_munu(0,1) = k_uu_10
 
    k_uu_munu(0,2) = k_uu_20
 
    k_uu_munu(0,3) = k_uu_30
 
    k_uu_munu(1,0) = k_uu_10
 
    k_uu_munu(1,1) = k_uu_11
 
    k_uu_munu(1,2) = k_uu_12
 
    k_uu_munu(1,3) = k_uu_13
 
    k_uu_munu(2,0) = k_uu_20
 
    k_uu_munu(2,1) = k_uu_12
 
    k_uu_munu(2,2) = k_uu_22
 
    k_uu_munu(2,3) = k_uu_23
 
    k_uu_munu(3,0) = k_uu_30
 
    k_uu_munu(3,1) = k_uu_13
 
    k_uu_munu(3,2) = k_uu_23
 
    k_uu_munu(3,3) = k_uu_33



  END SUBROUTINE ComputeEddingtonTensorComponents_uu



  SUBROUTINE ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )
   REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
               k_ud_munu(0:3,0:3)
    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_u_1, h_u_2, h_u_3, h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, V_0
    REAL(DP) :: u_d_1, u_d_2, u_d_3, u_u_1, u_u_2, u_u_3, W, x, y
    REAL(DP) :: k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, k_ud_21, k_ud_31, k_ud_32, & 
                k_ud_10, k_ud_20, k_ud_30, k_ud_01, k_ud_02, k_ud_03,  k_ud_00 

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )

    EF = EddingtonFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 
 
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

    k_ud_21 = a * u_u_2 * u_d_1 + b * h_u_2 * h_d_1
    k_ud_31 = a * u_u_3 * u_d_1 + b * h_u_3 * h_d_1
    k_ud_32 = a * u_u_3 * u_d_2 + b * h_u_3 * h_d_2

    x = ( 1.0_DP / alp ) * ( 1.0_DP / ( 1.0_DP - V_0 / alp ) ) 
    y = ( 1.0_DP / alp ) 
 
    k_ud_01 = x * ( V_d_1 * k_ud_11 + V_d_2 * k_ud_21 + V_d_3 * k_ud_31 )
    k_ud_02 = x * ( V_d_1 * k_ud_12 + V_d_2 * k_ud_22 + V_d_3 * k_ud_32 )
    k_ud_03 = x * ( V_d_1 * k_ud_13 + V_d_2 * k_ud_23 + V_d_3 * k_ud_33 )

    k_ud_10 = y * ( V_u_1 * k_ud_11 + V_u_2 * k_ud_12 + V_u_3 * k_ud_13 )
    k_ud_20 = y * ( V_u_1 * k_ud_21 + V_u_2 * k_ud_22 + V_u_3 * k_ud_23 )
    k_ud_30 = y * ( V_u_1 * k_ud_31 + V_u_2 * k_ud_32 + V_u_3 * k_ud_33 )

    k_ud_00 = x * y * ( V_u_1 * V_d_1 * k_ud_11 + V_u_1 * V_d_2 * k_ud_21 + V_u_1 * V_d_3 * k_ud_31 &
            + V_u_2 * V_d_1 * k_ud_12 + V_u_2 * V_d_2 * k_ud_22 + V_u_2 * V_d_3 * k_ud_32 &
            + V_u_3 * V_d_1 * k_ud_13 + V_u_3 * V_d_2 * k_ud_23 + V_u_3 * V_d_3 * k_ud_33 )
            
    k_ud_munu(0,0) = k_ud_00
 
    k_ud_munu(0,1) = k_ud_01
 
    k_ud_munu(0,2) = k_ud_02
 
    k_ud_munu(0,3) = k_ud_03
 
    k_ud_munu(1,0) = k_ud_10
 
    k_ud_munu(1,1) = k_ud_11
 
    k_ud_munu(1,2) = k_ud_12
 
    k_ud_munu(1,3) = k_ud_13
 
    k_ud_munu(2,0) = k_ud_20
 
    k_ud_munu(2,1) = k_ud_21
 
    k_ud_munu(2,2) = k_ud_22
 
    k_ud_munu(2,3) = k_ud_23
 
    k_ud_munu(3,0) = k_ud_30
 
    k_ud_munu(3,1) = k_ud_31
 
    k_ud_munu(3,2) = k_ud_32
 
    k_ud_munu(3,3) = k_ud_33
 

  END SUBROUTINE ComputeEddingtonTensorComponents_ud

  SUBROUTINE ComputeHeatFluxTensorComponents_ddd_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_ddd_ijk )
 
    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3

    REAL(DP), INTENT(out) :: &
      l_ddd_ijk(1:3,1:3,1:3)

    REAL(DP) :: FF, HF, a, b, DT
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W
    REAL(DP) :: h_dd_11, h_dd_22, h_dd_33, h_dd_12, h_dd_13, h_dd_23, h_dd_21, h_dd_31, h_dd_32
    REAL(DP) :: &
      l_ddd_111, l_ddd_112, l_ddd_113, l_ddd_121, l_ddd_122, & 
      l_ddd_123, l_ddd_131, l_ddd_132, l_ddd_133, l_ddd_211, & 
      l_ddd_212, l_ddd_213, l_ddd_221, l_ddd_222, l_ddd_223, & 
      l_ddd_231, l_ddd_232, l_ddd_233, l_ddd_311, l_ddd_312, & 
      l_ddd_313, l_ddd_321, l_ddd_322, l_ddd_323, l_ddd_331, & 
      l_ddd_332, l_ddd_333
   

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3


    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_dd_11 = Gm_dd_11 + u_d_1 * u_d_1
    h_dd_22 = Gm_dd_22 + u_d_2 * u_d_2
    h_dd_33 = Gm_dd_33 + u_d_3 * u_d_3
    h_dd_12 = u_d_1 * u_d_2
    h_dd_13 = u_d_1 * u_d_3
    h_dd_23 = u_d_2 * u_d_3
    h_dd_21 = u_d_1 * u_d_2
    h_dd_31 = u_d_1 * u_d_3
    h_dd_32 = u_d_2 * u_d_3
 
    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )


     l_ddd_111 & 
      = a * ( h_d_1 * h_dd_11 + h_d_1 * h_dd_11 + h_d_1 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_1

    l_ddd_112 & 
      = a * ( h_d_1 * h_dd_12 + h_d_1 * h_dd_12 + h_d_2 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_2

    l_ddd_113 & 
      = a * ( h_d_1 * h_dd_13 + h_d_1 * h_dd_13 + h_d_3 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_3

    l_ddd_121 & 
      = a * ( h_d_1 * h_dd_21 + h_d_2 * h_dd_11 + h_d_1 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_1

    l_ddd_122 & 
      = a * ( h_d_1 * h_dd_22 + h_d_2 * h_dd_12 + h_d_2 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_2

    l_ddd_123 & 
      = a * ( h_d_1 * h_dd_23 + h_d_2 * h_dd_13 + h_d_3 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_3

    l_ddd_131 & 
      = a * ( h_d_1 * h_dd_31 + h_d_3 * h_dd_11 + h_d_1 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_1

    l_ddd_132 & 
      = a * ( h_d_1 * h_dd_32 + h_d_3 * h_dd_12 + h_d_2 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_2

    l_ddd_133 & 
      = a * ( h_d_1 * h_dd_33 + h_d_3 * h_dd_13 + h_d_3 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_3

    l_ddd_211 & 
      = a * ( h_d_2 * h_dd_11 + h_d_1 * h_dd_21 + h_d_1 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_1

    l_ddd_212 & 
      = a * ( h_d_2 * h_dd_12 + h_d_1 * h_dd_22 + h_d_2 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_2

    l_ddd_213 & 
      = a * ( h_d_2 * h_dd_13 + h_d_1 * h_dd_23 + h_d_3 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_3

    l_ddd_221 & 
      = a * ( h_d_2 * h_dd_21 + h_d_2 * h_dd_21 + h_d_1 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_1

    l_ddd_222 & 
      = a * ( h_d_2 * h_dd_22 + h_d_2 * h_dd_22 + h_d_2 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_2

    l_ddd_223 & 
      = a * ( h_d_2 * h_dd_23 + h_d_2 * h_dd_23 + h_d_3 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_3

    l_ddd_231 & 
      = a * ( h_d_2 * h_dd_31 + h_d_3 * h_dd_21 + h_d_1 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_1

    l_ddd_232 & 
      = a * ( h_d_2 * h_dd_32 + h_d_3 * h_dd_22 + h_d_2 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_2

    l_ddd_233 & 
      = a * ( h_d_2 * h_dd_33 + h_d_3 * h_dd_23 + h_d_3 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_3

    l_ddd_311 & 
      = a * ( h_d_3 * h_dd_11 + h_d_1 * h_dd_31 + h_d_1 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_1

    l_ddd_312 & 
      = a * ( h_d_3 * h_dd_12 + h_d_1 * h_dd_32 + h_d_2 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_2

    l_ddd_313 & 
      = a * ( h_d_3 * h_dd_13 + h_d_1 * h_dd_33 + h_d_3 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_3

    l_ddd_321 & 
      = a * ( h_d_3 * h_dd_21 + h_d_2 * h_dd_31 + h_d_1 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_1

    l_ddd_322 & 
      = a * ( h_d_3 * h_dd_22 + h_d_2 * h_dd_32 + h_d_2 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_2

    l_ddd_323 & 
      = a * ( h_d_3 * h_dd_23 + h_d_2 * h_dd_33 + h_d_3 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_3

    l_ddd_331 & 
      = a * ( h_d_3 * h_dd_31 + h_d_3 * h_dd_31 + h_d_1 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_1

    l_ddd_332 & 
      = a * ( h_d_3 * h_dd_32 + h_d_3 * h_dd_32 + h_d_2 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_2

    l_ddd_333 & 
      = a * ( h_d_3 * h_dd_33 + h_d_3 * h_dd_33 + h_d_3 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_3

    l_ddd_ijk(1,1,1) = l_ddd_111
 
    l_ddd_ijk(1,1,2) = l_ddd_112
 
    l_ddd_ijk(1,1,3) = l_ddd_113
 
    l_ddd_ijk(1,2,1) = l_ddd_121
 
    l_ddd_ijk(1,2,2) = l_ddd_122
 
    l_ddd_ijk(1,2,3) = l_ddd_123
 
    l_ddd_ijk(1,3,1) = l_ddd_131
 
    l_ddd_ijk(1,3,2) = l_ddd_132
 
    l_ddd_ijk(1,3,3) = l_ddd_133
 
    l_ddd_ijk(2,1,1) = l_ddd_211
 
    l_ddd_ijk(2,1,2) = l_ddd_212
 
    l_ddd_ijk(2,1,3) = l_ddd_213
 
    l_ddd_ijk(2,2,1) = l_ddd_221
 
    l_ddd_ijk(2,2,2) = l_ddd_222
 
    l_ddd_ijk(2,2,3) = l_ddd_223
 
    l_ddd_ijk(2,3,1) = l_ddd_231
 
    l_ddd_ijk(2,3,2) = l_ddd_232
 
    l_ddd_ijk(2,3,3) = l_ddd_233
 
    l_ddd_ijk(3,1,1) = l_ddd_311
 
    l_ddd_ijk(3,1,2) = l_ddd_312
 
    l_ddd_ijk(3,1,3) = l_ddd_313
 
    l_ddd_ijk(3,2,1) = l_ddd_321
 
    l_ddd_ijk(3,2,2) = l_ddd_322
 
    l_ddd_ijk(3,2,3) = l_ddd_323
 
    l_ddd_ijk(3,3,1) = l_ddd_331
 
    l_ddd_ijk(3,3,2) = l_ddd_332
 
    l_ddd_ijk(3,3,3) = l_ddd_333

  END SUBROUTINE ComputeHeatFluxTensorComponents_ddd_Lagrangian

  SUBROUTINE ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )
    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
     REAL(DP), INTENT(out)  :: &  
      l_uud_munurho(0:3,0:3,0:3)

    REAL(DP) :: FF, HF, a, b, DT, V_0, x, y
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W, u_u_1, u_u_2, u_u_3
    REAL(DP) :: h_uu_11, h_uu_22, h_uu_33, h_uu_12, h_uu_13, h_uu_23, h_uu_21, h_uu_31, h_uu_32
    REAL(DP) :: h_ud_11, h_ud_22, h_ud_33, h_ud_12, h_ud_21, h_ud_13, h_ud_31, h_ud_23, h_ud_32
    REAL(DP) :: Gm_uu_11, Gm_uu_22, Gm_uu_33  
 
    REAL(DP) :: &
      l_uud_000, l_uud_001, l_uud_002, l_uud_003, & 
      l_uud_010, l_uud_011, l_uud_012, l_uud_013, & 
      l_uud_020, l_uud_021, l_uud_022, l_uud_023, & 
      l_uud_030, l_uud_031, l_uud_032, l_uud_033, & 
      l_uud_100, l_uud_101, l_uud_102, l_uud_103, & 
      l_uud_110, l_uud_111, l_uud_112, l_uud_113, & 
      l_uud_120, l_uud_121, l_uud_122, l_uud_123, & 
      l_uud_130, l_uud_131, l_uud_132, l_uud_133, & 
      l_uud_200, l_uud_201, l_uud_202, l_uud_203, & 
      l_uud_210, l_uud_211, l_uud_212, l_uud_213, & 
      l_uud_220, l_uud_221, l_uud_222, l_uud_223, & 
      l_uud_230, l_uud_231, l_uud_232, l_uud_233, & 
      l_uud_300, l_uud_301, l_uud_302, l_uud_303, & 
      l_uud_310, l_uud_311, l_uud_312, l_uud_313, & 
      l_uud_320, l_uud_321, l_uud_322, l_uud_323, & 
      l_uud_330, l_uud_331, l_uud_332, l_uud_333

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    Gm_uu_11 = 1.0_DP /  Gm_dd_11
    Gm_uu_22 = 1.0_DP /  Gm_dd_22
    Gm_uu_33 = 1.0_DP /  Gm_dd_33

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )


    h_uu_11 = Gm_uu_11 + u_u_1 * u_u_1
    h_uu_22 = Gm_uu_22 + u_u_2 * u_u_2
    h_uu_33 = Gm_uu_33 + u_u_3 * u_u_3
    h_uu_12 = u_u_1 * u_u_2
    h_uu_13 = u_u_1 * u_u_3
    h_uu_23 = u_u_2 * u_u_3
    h_uu_21 = u_u_1 * u_u_2
    h_uu_31 = u_u_1 * u_u_3
    h_uu_32 = u_u_2 * u_u_3


    h_ud_11 = 1.0_DP + u_u_1 * u_d_1
    h_ud_22 = 1.0_DP + u_u_2 * u_d_2
    h_ud_33 = 1.0_DP + u_u_3 * u_d_3
    h_ud_12 = u_u_1 * u_d_2
    h_ud_13 = u_u_1 * u_d_3
    h_ud_23 = u_u_2 * u_d_3
    h_ud_21 = u_u_2 * u_d_1
    h_ud_31 = u_u_3 * u_d_1
    h_ud_32 = u_u_3 * u_d_2
 




    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    ! --- Diagonal Heat Flux Tensor Components ---
    l_uud_111 & 
      = a * ( h_u_1 * h_ud_11 + h_u_1 * h_ud_11 + h_d_1 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_1

    l_uud_112 & 
      = a * ( h_u_1 * h_ud_12 + h_u_1 * h_ud_12 + h_d_2 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_2

    l_uud_113 & 
      = a * ( h_u_1 * h_ud_13 + h_u_1 * h_ud_13 + h_d_3 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_3

    l_uud_121 & 
      = a * ( h_u_1 * h_ud_21 + h_u_2 * h_ud_11 + h_d_1 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_1

    l_uud_122 & 
      = a * ( h_u_1 * h_ud_22 + h_u_2 * h_ud_12 + h_d_2 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_2

    l_uud_123 & 
      = a * ( h_u_1 * h_ud_23 + h_u_2 * h_ud_13 + h_d_3 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_3

    l_uud_131 & 
      = a * ( h_u_1 * h_ud_31 + h_u_3 * h_ud_11 + h_d_1 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_1

    l_uud_132 & 
      = a * ( h_u_1 * h_ud_32 + h_u_3 * h_ud_12 + h_d_2 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_2

    l_uud_133 & 
      = a * ( h_u_1 * h_ud_33 + h_u_3 * h_ud_13 + h_d_3 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_3

    l_uud_211 & 
      = a * ( h_u_2 * h_ud_11 + h_u_1 * h_ud_21 + h_d_1 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_1

    l_uud_212 & 
      = a * ( h_u_2 * h_ud_12 + h_u_1 * h_ud_22 + h_d_2 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_2

    l_uud_213 & 
      = a * ( h_u_2 * h_ud_13 + h_u_1 * h_ud_23 + h_d_3 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_3

    l_uud_221 & 
      = a * ( h_u_2 * h_ud_21 + h_u_2 * h_ud_21 + h_d_1 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_1

    l_uud_222 & 
      = a * ( h_u_2 * h_ud_22 + h_u_2 * h_ud_22 + h_d_2 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_2

    l_uud_223 & 
      = a * ( h_u_2 * h_ud_23 + h_u_2 * h_ud_23 + h_d_3 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_3

    l_uud_231 & 
      = a * ( h_u_2 * h_ud_31 + h_u_3 * h_ud_21 + h_d_1 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_1

    l_uud_232 & 
      = a * ( h_u_2 * h_ud_32 + h_u_3 * h_ud_22 + h_d_2 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_2

    l_uud_233 & 
      = a * ( h_u_2 * h_ud_33 + h_u_3 * h_ud_23 + h_d_3 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_3

    l_uud_311 & 
      = a * ( h_u_3 * h_ud_11 + h_u_1 * h_ud_31 + h_d_1 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_1

    l_uud_312 & 
      = a * ( h_u_3 * h_ud_12 + h_u_1 * h_ud_32 + h_d_2 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_2

    l_uud_313 & 
      = a * ( h_u_3 * h_ud_13 + h_u_1 * h_ud_33 + h_d_3 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_3

    l_uud_321 & 
      = a * ( h_u_3 * h_ud_21 + h_u_2 * h_ud_31 + h_d_1 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_1

    l_uud_322 & 
      = a * ( h_u_3 * h_ud_22 + h_u_2 * h_ud_32 + h_d_2 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_2

    l_uud_323 & 
      = a * ( h_u_3 * h_ud_23 + h_u_2 * h_ud_33 + h_d_3 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_3

    l_uud_331 & 
      = a * ( h_u_3 * h_ud_31 + h_u_3 * h_ud_31 + h_d_1 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_1

    l_uud_332 & 
      = a * ( h_u_3 * h_ud_32 + h_u_3 * h_ud_32 + h_d_2 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_2

    l_uud_333 & 
      = a * ( h_u_3 * h_ud_33 + h_u_3 * h_ud_33 + h_d_3 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    x = ( 1.0_DP / alp ) * ( 1.0_DP / ( 1.0_DP - V_0 / alp ) ) 
    y = ( 1.0_DP / alp ) 

    l_uud_011 = x * (V_d_1 * l_uud_111 + V_d_2 * l_uud_211 + V_d_3 * l_uud_311)
 
    l_uud_012 = x * (V_d_1 * l_uud_112 + V_d_2 * l_uud_212 + V_d_3 * l_uud_312) 
 
    l_uud_013 = x * (V_d_1 * l_uud_113 + V_d_2 * l_uud_213 + V_d_3 * l_uud_313) 
 
    l_uud_021 = x * (V_d_1 * l_uud_121 + V_d_2 * l_uud_221 + V_d_3 * l_uud_321) 
 
    l_uud_022 = x * (V_d_1 * l_uud_122 + V_d_2 * l_uud_222 + V_d_3 * l_uud_322) 
 
    l_uud_023 = x * (V_d_1 * l_uud_123 + V_d_2 * l_uud_223 + V_d_3 * l_uud_323) 
 
    l_uud_031 = x * (V_d_1 * l_uud_131 + V_d_2 * l_uud_231 + V_d_3 * l_uud_331) 
 
    l_uud_032 = x * (V_d_1 * l_uud_132 + V_d_2 * l_uud_232 + V_d_3 * l_uud_332) 

    l_uud_033 = x * (V_d_1 * l_uud_133 + V_d_2 * l_uud_233 + V_d_3 * l_uud_333)
 
    l_uud_101 = x * (V_d_1 * l_uud_111 + V_d_2 * l_uud_121 + V_d_3 * l_uud_131) 
 
    l_uud_102 = x * (V_d_1 * l_uud_112 + V_d_2 * l_uud_122 + V_d_3 * l_uud_132) 
 
    l_uud_103 = x * (V_d_1 * l_uud_113 + V_d_2 * l_uud_123 + V_d_3 * l_uud_133) 
 
    l_uud_201 = x * (V_d_1 * l_uud_211 + V_d_2 * l_uud_221 + V_d_3 * l_uud_231) 
 
    l_uud_202 = x * (V_d_1 * l_uud_212 + V_d_2 * l_uud_222 + V_d_3 * l_uud_232) 
 
    l_uud_203 = x * (V_d_1 * l_uud_213 + V_d_2 * l_uud_223 + V_d_3 * l_uud_233) 
 
    l_uud_301 = x * (V_d_1 * l_uud_311 + V_d_2 * l_uud_321 + V_d_3 * l_uud_331) 
 
    l_uud_302 = x * (V_d_1 * l_uud_312 + V_d_2 * l_uud_322 + V_d_3 * l_uud_332) 
 
    l_uud_303 = x * (V_d_1 * l_uud_313 + V_d_2 * l_uud_323 + V_d_3 * l_uud_333)
 
    l_uud_110 = y * (V_u_1 * l_uud_111 + V_u_2 * l_uud_112 + V_u_3 * l_uud_113) 
 
    l_uud_120 = y * (V_u_1 * l_uud_121 + V_u_2 * l_uud_122 + V_u_3 * l_uud_123) 
 
    l_uud_130 = y * (V_u_1 * l_uud_131 + V_u_2 * l_uud_132 + V_u_3 * l_uud_133)
 
    l_uud_210 = y * (V_u_1 * l_uud_211 + V_u_2 * l_uud_212 + V_u_3 * l_uud_213)
 
    l_uud_220 = y * (V_u_1 * l_uud_221 + V_u_2 * l_uud_222 + V_u_3 * l_uud_223)
 
    l_uud_230 = y * (V_u_1 * l_uud_231 + V_u_2 * l_uud_232 + V_u_3 * l_uud_233)
 
    l_uud_310 = y * (V_u_1 * l_uud_311 + V_u_2 * l_uud_312 + V_u_3 * l_uud_313)
 
    l_uud_320 = y * (V_u_1 * l_uud_321 + V_u_2 * l_uud_322 + V_u_3 * l_uud_323)
 
    l_uud_330 = y * (V_u_1 * l_uud_331 + V_u_2 * l_uud_332 + V_u_3 * l_uud_333)

    l_uud_001 = x**2 * (  V_d_1 * V_d_1 * l_uud_111 + V_d_1 * V_d_2 * l_uud_121 + V_d_1 * V_d_3 * l_uud_131 & 
              + V_d_2 * V_d_1 * l_uud_211 + V_d_2 * V_d_2 * l_uud_221 + V_d_2 * V_d_3 * l_uud_231 & 
              + V_d_3 * V_d_1 * l_uud_311 + V_d_3 * V_d_2 * l_uud_321 + V_d_3 * V_d_3 * l_uud_331 ) 
             
 
    l_uud_002 = x**2 * (  V_d_1 * V_d_1 * l_uud_112 + V_d_1 * V_d_2 * l_uud_122 + V_d_1 * V_d_3 * l_uud_132 & 
              + V_d_2 * V_d_1 * l_uud_212 + V_d_2 * V_d_2 * l_uud_222 + V_d_2 * V_d_3 * l_uud_232 & 
              + V_d_3 * V_d_1 * l_uud_312 + V_d_3 * V_d_2 * l_uud_322 + V_d_3 * V_d_3 * l_uud_332 ) 
             
 
    l_uud_003 = x**2 * (  V_d_1 * V_d_1 * l_uud_113 + V_d_1 * V_d_2 * l_uud_123 + V_d_1 * V_d_3 * l_uud_133 & 
              + V_d_2 * V_d_1 * l_uud_213 + V_d_2 * V_d_2 * l_uud_223 + V_d_2 * V_d_3 * l_uud_233 & 
              + V_d_3 * V_d_1 * l_uud_313 + V_d_3 * V_d_2 * l_uud_323 + V_d_3 * V_d_3 * l_uud_333 ) 

    l_uud_100 = x * y * (  V_d_1 * V_u_1 * l_uud_111 + V_d_1 * V_u_2 * l_uud_112 + V_d_1 * V_u_3 * l_uud_113 & 
              + V_d_2 * V_u_1 * l_uud_121 + V_d_2 * V_u_2 * l_uud_122 + V_d_2 * V_u_3 * l_uud_123 & 
              + V_d_3 * V_u_1 * l_uud_131 + V_d_3 * V_u_2 * l_uud_132 + V_d_3 * V_u_3 * l_uud_133 ) 
             
 
    l_uud_200 = x * y * (  V_d_1 * V_u_1 * l_uud_211 + V_d_1 * V_u_2 * l_uud_212 + V_d_1 * V_u_3 * l_uud_213 & 
              + V_d_2 * V_u_1 * l_uud_221 + V_d_2 * V_u_2 * l_uud_222 + V_d_2 * V_u_3 * l_uud_223 & 
              + V_d_3 * V_u_1 * l_uud_231 + V_d_3 * V_u_2 * l_uud_232 + V_d_3 * V_u_3 * l_uud_233 ) 
             
    l_uud_300 = x * y * (  V_d_1 * V_u_1 * l_uud_311 + V_d_1 * V_u_2 * l_uud_312 + V_d_1 * V_u_3 * l_uud_313 & 
              + V_d_2 * V_u_1 * l_uud_321 + V_d_2 * V_u_2 * l_uud_322 + V_d_2 * V_u_3 * l_uud_323 & 
              + V_d_3 * V_u_1 * l_uud_331 + V_d_3 * V_u_2 * l_uud_332 + V_d_3 * V_u_3 * l_uud_333 ) 


    l_uud_010 = x * y * (  V_d_1 * V_u_1 * l_uud_111 + V_d_1 * V_u_2 * l_uud_112 + V_d_1 * V_u_3 * l_uud_113 & 
              + V_d_2 * V_u_1 * l_uud_211 + V_d_2 * V_u_2 * l_uud_212 + V_d_2 * V_u_3 * l_uud_213 & 
              + V_d_3 * V_u_1 * l_uud_311 + V_d_3 * V_u_2 * l_uud_312 + V_d_3 * V_u_3 * l_uud_313 ) 
             
 
    l_uud_020 = x * y * (  V_d_1 * V_u_1 * l_uud_121 + V_d_1 * V_u_2 * l_uud_122 + V_d_1 * V_u_3 * l_uud_123 & 
              + V_d_2 * V_u_1 * l_uud_221 + V_d_2 * V_u_2 * l_uud_222 + V_d_2 * V_u_3 * l_uud_223 & 
              + V_d_3 * V_u_1 * l_uud_321 + V_d_3 * V_u_2 * l_uud_322 + V_d_3 * V_u_3 * l_uud_323 ) 
             
 
    l_uud_030 = x * y * (  V_d_1 * V_u_1 * l_uud_131 + V_d_1 * V_u_2 * l_uud_132 + V_d_1 * V_u_3 * l_uud_133 & 
              + V_d_2 * V_u_1 * l_uud_231 + V_d_2 * V_u_2 * l_uud_232 + V_d_2 * V_u_3 * l_uud_233 & 
              + V_d_3 * V_u_1 * l_uud_331 + V_d_3 * V_u_2 * l_uud_332 + V_d_3 * V_u_3 * l_uud_333 ) 

    l_uud_000 &
      = x**2 * y * (   V_d_1 * V_d_1 * V_u_1 * l_uud_111 &
                     + V_d_1 * V_d_1 * V_u_2 * l_uud_112 &
                     + V_d_1 * V_d_1 * V_u_3 * l_uud_113 & 
                     + V_d_1 * V_d_2 * V_u_1 * l_uud_121 &
                     + V_d_1 * V_d_2 * V_u_2 * l_uud_122 &
                     + V_d_1 * V_d_2 * V_u_3 * l_uud_123 & 
                     + V_d_1 * V_d_3 * V_u_1 * l_uud_131 &
                     + V_d_1 * V_d_3 * V_u_2 * l_uud_132 &
                     + V_d_1 * V_d_3 * V_u_3 * l_uud_133 & 
                     + V_d_2 * V_d_1 * V_u_1 * l_uud_211 &
                     + V_d_2 * V_d_1 * V_u_2 * l_uud_212 &
                     + V_d_2 * V_d_1 * V_u_3 * l_uud_213 & 
                     + V_d_2 * V_d_2 * V_u_1 * l_uud_221 &
                     + V_d_2 * V_d_2 * V_u_2 * l_uud_222 &
                     + V_d_2 * V_d_2 * V_u_3 * l_uud_223 & 
                     + V_d_2 * V_d_3 * V_u_1 * l_uud_231 &
                     + V_d_2 * V_d_3 * V_u_2 * l_uud_232 &
                     + V_d_2 * V_d_3 * V_u_3 * l_uud_233 & 
                     + V_d_3 * V_d_1 * V_u_1 * l_uud_311 &
                     + V_d_3 * V_d_1 * V_u_2 * l_uud_312 &
                     + V_d_3 * V_d_1 * V_u_3 * l_uud_313 & 
                     + V_d_3 * V_d_2 * V_u_1 * l_uud_321 &
                     + V_d_3 * V_d_2 * V_u_2 * l_uud_322 &
                     + V_d_3 * V_d_2 * V_u_3 * l_uud_323 & 
                     + V_d_3 * V_d_3 * V_u_1 * l_uud_331 &
                     + V_d_3 * V_d_3 * V_u_2 * l_uud_332 &
                     + V_d_3 * V_d_3 * V_u_3 * l_uud_333 )

    l_uud_munurho(0,0,0) = l_uud_000
 
    l_uud_munurho(0,0,1) = l_uud_001
 
    l_uud_munurho(0,0,2) = l_uud_002
 
    l_uud_munurho(0,0,3) = l_uud_003
 
    l_uud_munurho(0,1,0) = l_uud_010
 
    l_uud_munurho(0,1,1) = l_uud_011
 
    l_uud_munurho(0,1,2) = l_uud_012
 
    l_uud_munurho(0,1,3) = l_uud_013
 
    l_uud_munurho(0,2,0) = l_uud_020
 
    l_uud_munurho(0,2,1) = l_uud_021
 
    l_uud_munurho(0,2,2) = l_uud_022
 
    l_uud_munurho(0,2,3) = l_uud_023
 
    l_uud_munurho(0,3,0) = l_uud_030
 
    l_uud_munurho(0,3,1) = l_uud_031
 
    l_uud_munurho(0,3,2) = l_uud_032
 
    l_uud_munurho(0,3,3) = l_uud_033
 
    l_uud_munurho(1,0,0) = l_uud_100
 
    l_uud_munurho(1,0,1) = l_uud_101
 
    l_uud_munurho(1,0,2) = l_uud_102
 
    l_uud_munurho(1,0,3) = l_uud_103
 
    l_uud_munurho(1,1,0) = l_uud_110
 
    l_uud_munurho(1,1,1) = l_uud_111
 
    l_uud_munurho(1,1,2) = l_uud_112
 
    l_uud_munurho(1,1,3) = l_uud_113
 
    l_uud_munurho(1,2,0) = l_uud_120
 
    l_uud_munurho(1,2,1) = l_uud_121
 
    l_uud_munurho(1,2,2) = l_uud_122
 
    l_uud_munurho(1,2,3) = l_uud_123
 
    l_uud_munurho(1,3,0) = l_uud_130
 
    l_uud_munurho(1,3,1) = l_uud_131
 
    l_uud_munurho(1,3,2) = l_uud_132
 
    l_uud_munurho(1,3,3) = l_uud_133
 
    l_uud_munurho(2,0,0) = l_uud_200
 
    l_uud_munurho(2,0,1) = l_uud_201
 
    l_uud_munurho(2,0,2) = l_uud_202
 
    l_uud_munurho(2,0,3) = l_uud_203
 
    l_uud_munurho(2,1,0) = l_uud_210
 
    l_uud_munurho(2,1,1) = l_uud_211
 
    l_uud_munurho(2,1,2) = l_uud_212
 
    l_uud_munurho(2,1,3) = l_uud_213
 
    l_uud_munurho(2,2,0) = l_uud_220
 
    l_uud_munurho(2,2,1) = l_uud_221
 
    l_uud_munurho(2,2,2) = l_uud_222
 
    l_uud_munurho(2,2,3) = l_uud_223
 
    l_uud_munurho(2,3,0) = l_uud_230
 
    l_uud_munurho(2,3,1) = l_uud_231
 
    l_uud_munurho(2,3,2) = l_uud_232
 
    l_uud_munurho(2,3,3) = l_uud_233
 
    l_uud_munurho(3,0,0) = l_uud_300
 
    l_uud_munurho(3,0,1) = l_uud_301
 
    l_uud_munurho(3,0,2) = l_uud_302
 
    l_uud_munurho(3,0,3) = l_uud_303
 
    l_uud_munurho(3,1,0) = l_uud_310
 
    l_uud_munurho(3,1,1) = l_uud_311
 
    l_uud_munurho(3,1,2) = l_uud_312
 
    l_uud_munurho(3,1,3) = l_uud_313
 
    l_uud_munurho(3,2,0) = l_uud_320
 
    l_uud_munurho(3,2,1) = l_uud_321
 
    l_uud_munurho(3,2,2) = l_uud_322
 
    l_uud_munurho(3,2,3) = l_uud_323
 
    l_uud_munurho(3,3,0) = l_uud_330
 
    l_uud_munurho(3,3,1) = l_uud_331
 
    l_uud_munurho(3,3,2) = l_uud_332
 
    l_uud_munurho(3,3,3) = l_uud_333



  END SUBROUTINE ComputeHeatFluxTensorComponents_uud_Lagrangian





  SUBROUTINE ComputeHeatFluxTensorComponents_udd_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_udd_ijk )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      l_udd_ijk(1:3,1:3,1:3)      

    REAL(DP) :: &
      l_udd_111, l_udd_112, l_udd_113, l_udd_121, l_udd_122, & 
      l_udd_123, l_udd_131, l_udd_132, l_udd_133, l_udd_211, & 
      l_udd_212, l_udd_213, l_udd_221, l_udd_222, l_udd_223, & 
      l_udd_231, l_udd_232, l_udd_233, l_udd_311, l_udd_312, & 
      l_udd_313, l_udd_321, l_udd_322, l_udd_323, l_udd_331, & 
      l_udd_332, l_udd_333 
    REAL(DP) :: FF, HF, a, b, DT
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W, u_u_1, u_u_2, u_u_3
    REAL(DP) :: h_dd_11, h_dd_22, h_dd_33, h_dd_12, h_dd_13, h_dd_23, h_dd_21, h_dd_31, h_dd_32
    REAL(DP) :: h_ud_11, h_ud_22, h_ud_33, h_ud_12, h_ud_21, h_ud_13, h_ud_31, h_ud_23, h_ud_32
   

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_dd_11 = Gm_dd_11 + u_d_1 * u_d_1
    h_dd_22 = Gm_dd_22 + u_d_2 * u_d_2
    h_dd_33 = Gm_dd_33 + u_d_3 * u_d_3
    h_dd_12 = u_d_1 * u_d_2
    h_dd_13 = u_d_1 * u_d_3
    h_dd_23 = u_d_2 * u_d_3
    h_dd_21 = u_d_1 * u_d_2
    h_dd_31 = u_d_1 * u_d_3
    h_dd_32 = u_d_2 * u_d_3


    h_ud_11 = 1.0_DP + u_u_1 * u_d_1
    h_ud_22 = 1.0_DP + u_u_2 * u_d_2
    h_ud_33 = 1.0_DP + u_u_3 * u_d_3
    h_ud_12 = u_u_1 * u_d_2
    h_ud_13 = u_u_1 * u_d_3
    h_ud_23 = u_u_2 * u_d_3
    h_ud_21 = u_u_2 * u_d_1
    h_ud_31 = u_u_3 * u_d_1
    h_ud_32 = u_u_3 * u_d_2
 




    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    ! --- Diagonal Heat Flux Tensor Components ---

    l_udd_111 & 
      = a * ( h_u_1 * h_dd_11 + h_d_1 * h_ud_11 + h_d_1 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_1

    l_udd_112 & 
      = a * ( h_u_1 * h_dd_12 + h_d_1 * h_ud_12 + h_d_2 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_2

    l_udd_113 & 
      = a * ( h_u_1 * h_dd_13 + h_d_1 * h_ud_13 + h_d_3 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_3

    l_udd_121 & 
      = a * ( h_u_1 * h_dd_21 + h_d_2 * h_ud_11 + h_d_1 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_1

    l_udd_122 & 
      = a * ( h_u_1 * h_dd_22 + h_d_2 * h_ud_12 + h_d_2 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_2

    l_udd_123 & 
      = a * ( h_u_1 * h_dd_23 + h_d_2 * h_ud_13 + h_d_3 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_3

    l_udd_131 & 
      = a * ( h_u_1 * h_dd_31 + h_d_3 * h_ud_11 + h_d_1 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_1

    l_udd_132 & 
      = a * ( h_u_1 * h_dd_32 + h_d_3 * h_ud_12 + h_d_2 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_2

    l_udd_133 & 
      = a * ( h_u_1 * h_dd_33 + h_d_3 * h_ud_13 + h_d_3 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_3

    l_udd_211 & 
      = a * ( h_u_2 * h_dd_11 + h_d_1 * h_ud_21 + h_d_1 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_1

    l_udd_212 & 
      = a * ( h_u_2 * h_dd_12 + h_d_1 * h_ud_22 + h_d_2 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_2

    l_udd_213 & 
      = a * ( h_u_2 * h_dd_13 + h_d_1 * h_ud_23 + h_d_3 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_3

    l_udd_221 & 
      = a * ( h_u_2 * h_dd_21 + h_d_2 * h_ud_21 + h_d_1 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_1

    l_udd_222 & 
      = a * ( h_u_2 * h_dd_22 + h_d_2 * h_ud_22 + h_d_2 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_2

    l_udd_223 & 
      = a * ( h_u_2 * h_dd_23 + h_d_2 * h_ud_23 + h_d_3 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_3

    l_udd_231 & 
      = a * ( h_u_2 * h_dd_31 + h_d_3 * h_ud_21 + h_d_1 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_1

    l_udd_232 & 
      = a * ( h_u_2 * h_dd_32 + h_d_3 * h_ud_22 + h_d_2 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_2

    l_udd_233 & 
      = a * ( h_u_2 * h_dd_33 + h_d_3 * h_ud_23 + h_d_3 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_3

    l_udd_311 & 
      = a * ( h_u_3 * h_dd_11 + h_d_1 * h_ud_31 + h_d_1 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_1

    l_udd_312 & 
      = a * ( h_u_3 * h_dd_12 + h_d_1 * h_ud_32 + h_d_2 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_2

    l_udd_313 & 
      = a * ( h_u_3 * h_dd_13 + h_d_1 * h_ud_33 + h_d_3 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_3

    l_udd_321 & 
      = a * ( h_u_3 * h_dd_21 + h_d_2 * h_ud_31 + h_d_1 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_1

    l_udd_322 & 
      = a * ( h_u_3 * h_dd_22 + h_d_2 * h_ud_32 + h_d_2 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_2

    l_udd_323 & 
      = a * ( h_u_3 * h_dd_23 + h_d_2 * h_ud_33 + h_d_3 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_3

    l_udd_331 & 
      = a * ( h_u_3 * h_dd_31 + h_d_3 * h_ud_31 + h_d_1 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_1

    l_udd_332 & 
      = a * ( h_u_3 * h_dd_32 + h_d_3 * h_ud_32 + h_d_2 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_2

    l_udd_333 & 
      = a * ( h_u_3 * h_dd_33 + h_d_3 * h_ud_33 + h_d_3 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_3

    l_udd_ijk(1,1,1) = l_udd_111
 
    l_udd_ijk(1,1,2) = l_udd_112
 
    l_udd_ijk(1,1,3) = l_udd_113
 
    l_udd_ijk(1,2,1) = l_udd_121
 
    l_udd_ijk(1,2,2) = l_udd_122
 
    l_udd_ijk(1,2,3) = l_udd_123
 
    l_udd_ijk(1,3,1) = l_udd_131
 
    l_udd_ijk(1,3,2) = l_udd_132
 
    l_udd_ijk(1,3,3) = l_udd_133
 
    l_udd_ijk(2,1,1) = l_udd_211
 
    l_udd_ijk(2,1,2) = l_udd_212
 
    l_udd_ijk(2,1,3) = l_udd_213
 
    l_udd_ijk(2,2,1) = l_udd_221
 
    l_udd_ijk(2,2,2) = l_udd_222
 
    l_udd_ijk(2,2,3) = l_udd_223
 
    l_udd_ijk(2,3,1) = l_udd_231
 
    l_udd_ijk(2,3,2) = l_udd_232
 
    l_udd_ijk(2,3,3) = l_udd_233
 
    l_udd_ijk(3,1,1) = l_udd_311
 
    l_udd_ijk(3,1,2) = l_udd_312
 
    l_udd_ijk(3,1,3) = l_udd_313
 
    l_udd_ijk(3,2,1) = l_udd_321
 
    l_udd_ijk(3,2,2) = l_udd_322
 
    l_udd_ijk(3,2,3) = l_udd_323
 
    l_udd_ijk(3,3,1) = l_udd_331
 
    l_udd_ijk(3,3,2) = l_udd_332
 
    l_udd_ijk(3,3,3) = l_udd_333
 

  END SUBROUTINE ComputeHeatFluxTensorComponents_udd_Lagrangian





!  SUBROUTINE ComputeStressEnergyComponents_Eulerian &
!                        ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3, &
!                         alp, B_u_1, B_u_2, B_u_3, EP, F_u_1, F_u_2, F_u_3, &
!                         S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32 )
!
!    REAL(DP), INTENT(in)  :: &
!      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
!    REAL(DP), INTENT(out) :: &
!      EP, F_u_1, F_u_2, F_u_3
!    REAL(DP), INTENT(out) :: &
!      S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32
!
!    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
!    REAL(DP) :: k_ud_11, k_ud_22, k_ud_33, k_ud_12, k_ud_13, k_ud_23
!    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3
!    REAL(DP) :: V_d_1, V_d_2, V_d_3
!    REAL(DP) :: n_u_1, n_u_2, n_u_3
!
!
!    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
!               + Gm_dd_22 * V_u_2 * V_u_2 &  
!               + Gm_dd_33 * V_u_3 * V_u_3) )
!
!    B_d_1 = Gm_dd_11 * B_u_1
!    B_d_2 = Gm_dd_22 * B_u_2
!    B_d_3 = Gm_dd_33 * B_u_3
!
!    V_d_1 = Gm_dd_11 * V_u_1
!    V_d_2 = Gm_dd_22 * V_u_2
!    V_d_3 = Gm_dd_33 * V_u_3
!
!    n_u_1 = -B_u_1 / alp    
!    n_u_2 = -B_u_2 / alp    
!    n_u_3 = -B_u_3 / alp    
!
!
!    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )
!
!    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
!          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
!    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
!          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
!    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
!          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )
!    CALL ComputeEddingtonTensorComponents_dd &
!           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!             k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
!             alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3   )
!
!    CALL ComputeEddingtonTensorComponents_ud &
!    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33 )
!
!    EP = ( W**2 * D + 2.0_DP * W * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 ) )
!    EP = EP + ( V_u_1 * V_u_1 * k_dd_11 * D + V_u_2 * V_u_2 * k_dd_22 * D + V_u_3 * V_u_3 * k_dd_33 * D )
!    EP = EP + 2.0_DP * ( V_u_1 * V_u_2 * k_dd_12 * D + V_u_1 * V_u_3 * k_dd_13 * D + V_u_2 * V_u_3 * k_dd_23 * D )
!
!    F_u_1 = ( W * I_u_1 ) 
!    F_u_1 = F_u_1 + V_u_1 * ( n_u_1 * W * I_d_1 + k_ud_11 * D ) + V_u_2 * ( n_u_1 * W * I_d_2 + k_ud_12 * D ) + V_u_3 * ( n_u_1 * W * I_d_3 + k_ud_13 * D )
!    F_u_1 = F_u_1 + ( W**2 * V_u_1 * D )
!    F_u_1 = F_u_1 +  n_u_1 * ( V_u_1**2 * k_dd_11 * D + V_u_2**2 * k_dd_22 * D + V_u_3**2 * k_dd_33 * D ) 
!    F_u_1 = F_u_1 + 2.0_DP * n_u_1 * ( V_u_1 * V_u_2 * k_dd_12 * D + V_u_1 * V_u_3 * k_dd_13 * D + V_u_2 * V_u_3 * k_dd_23 * D )
!    F_u_1 = F_u_1 +  W * V_u_1 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    F_u_2 = ( W * I_u_2 ) 
!    F_u_2 = F_u_2 +  V_u_1 * ( n_u_2 * W * I_d_1 + k_ud_12 * D ) + V_u_2 * ( n_u_2 * W * I_d_2 + k_ud_22 * D ) + V_u_3 * ( n_u_2 * W * I_d_3 + k_ud_23 * D )
!    F_u_2 = F_u_2 + ( W**2 * V_u_2 * D )
!    F_u_2 = F_u_2 +  n_u_2 * ( V_u_1**2 * k_dd_11 * D + V_u_2**2 * k_dd_22 * D + V_u_3**2 * k_dd_33 * D ) 
!    F_u_2 = F_u_2 + 2.0_DP * n_u_2 * ( V_u_1 * V_u_2 * k_dd_12 * D + V_u_1 * V_u_3 * k_dd_13 * D + V_u_2 * V_u_3 * k_dd_23 * D )
!    F_u_2 = F_u_2 + W * V_u_2 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    F_u_3 = ( W * I_u_3 ) 
!    F_u_3 = F_u_3 +  V_u_1 * ( n_u_3 * W * I_d_1 + k_ud_13 * D ) + V_u_2 * ( n_u_3 * W * I_d_2 + k_ud_23 * D ) + V_u_3 * ( n_u_3 * W * I_d_3 + k_ud_33 * D )
!    F_u_3 = F_u_3 + ( W**2 * V_u_3 * D )
!    F_u_3 = F_u_3 +  n_u_3 * ( V_u_1**2 * k_dd_11 * D + V_u_2**2 * k_dd_22 * D + V_u_3**2 * k_dd_33 * D ) 
!    F_u_3 = F_u_3 + 2.0_DP *  n_u_3 * ( V_u_1 * V_u_2 * k_dd_12 * D + V_u_1 * V_u_3 * k_dd_13 * D + V_u_2 * V_u_3 * k_dd_23 * D )
!    F_u_3 = F_u_3 + W * V_u_3 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!
!
!    S_ud_11 =  ( k_ud_11 * D + W * ( I_u_1 * V_d_1 + V_u_1 * I_d_1 ) + V_d_1 * V_u_1 * W**2 * D )
!    S_ud_11 = S_ud_11 -  n_u_1 * ( V_u_1 * k_dd_11 * D + V_u_2 * k_dd_12 * D + V_u_3 * k_dd_13 * D )
!    S_ud_11 = S_ud_11 - V_d_1 * n_u_1 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_22 =  ( k_ud_22 * D + W * ( I_u_2 * V_d_2 + V_u_2 * I_d_2 ) + V_d_2 * V_u_2 * W**2 * D )
!    S_ud_22 = S_ud_22 -  n_u_2 * ( V_u_1 * k_dd_12 * D + V_u_2 * k_dd_22 * D + V_u_3 * k_dd_23 * D )
!    S_ud_22 = S_ud_22 - V_d_2 * n_u_2 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_33 = ( k_ud_33 * D + W * ( I_u_3 * V_d_3 + V_u_3 * I_d_3 ) + V_d_3 * V_u_3 * W**2 * D )
!    S_ud_33 = S_ud_33 - n_u_3 * ( V_u_1 * k_dd_13 * D + V_u_2 * k_dd_23 * D + V_u_3 * k_dd_33 * D )
!    S_ud_33 = S_ud_33 -  V_d_3 * n_u_3 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_12 = ( k_ud_12 * D + W * ( I_u_1 * V_d_2 + V_u_1 * I_d_2 ) + V_d_1 * V_u_2 * W**2 * D )
!    S_ud_12 = S_ud_12 - n_u_1 * ( V_u_1 * k_dd_12 * D + V_u_2 * k_dd_22 * D + V_u_3 * k_dd_23 * D )
!    S_ud_12 = S_ud_12 - V_d_2 * n_u_1 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_23 = ( k_ud_23 * D + W * ( I_u_2 * V_d_3 + V_u_2 * I_d_3 ) + V_d_2 * V_u_3 * W**2 * D )
!    S_ud_23 = S_ud_23 - n_u_2 * ( V_u_1 * k_dd_13 * D + V_u_2 * k_dd_23 * D + V_u_3 * k_dd_33 * D )
!    S_ud_23 = S_ud_23 - V_d_3 * n_u_2 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_13 = ( k_ud_13 * D + W * ( I_u_1 * V_d_3 + V_u_1 * I_d_3 ) + V_d_1 * V_u_3 * W**2 * D )
!    S_ud_13 = S_ud_13 - n_u_1 * ( V_u_1 * k_dd_13 * D + V_u_2 * k_dd_23 * D + V_u_3 * k_dd_33 * D )
!    S_ud_13 = S_ud_13 - V_d_3 * n_u_1 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_21 = ( k_ud_12 * D + W * ( I_u_2 * V_d_1 + V_u_2 * I_d_1 ) + V_d_2 * V_u_1 * W**2 * D )
!    S_ud_21 = S_ud_21 - n_u_2 * ( V_u_1 * k_dd_11 * D + V_u_2 * k_dd_12 * D + V_u_3 * k_dd_13 * D )
!    S_ud_21 = S_ud_21 - V_d_1 * n_u_2 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_32 =  ( k_ud_23 * D + W * ( I_u_3 * V_d_2 + V_u_3 * I_d_2 ) + V_d_3 * V_u_2 * W**2 * D )
!    S_ud_32 = S_ud_32 -  n_u_3 * ( V_u_1 * k_dd_12 * D + V_u_2 * k_dd_22 * D + V_u_3 * k_dd_23 * D )
!    S_ud_32 = S_ud_32 - V_d_2 * n_u_3 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!    S_ud_31 = ( k_ud_13 * D + W * ( I_u_3 * V_d_1 + V_u_3 * I_d_1 ) + V_d_3 * V_u_1 * W**2 * D )
!    S_ud_31 = S_ud_31 - n_u_3 * ( V_u_1 * k_dd_11 * D + V_u_2 * k_dd_12 * D + V_u_3 * k_dd_13 * D )
!    S_ud_31 = S_ud_31 - V_d_1 * n_u_3 * ( V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3 )
!
!
!  END SUBROUTINE ComputeStressEnergyComponents_Eulerian
!
!
!  SUBROUTINE ComputeHeatFluxTensorComponents_udd_Eulerian &
!                                          ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                                           alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, &
!                                           W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!                                           W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!                                           W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!                                           W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!                                           W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!                                           W_udd_332, W_udd_333) 
! 
!    REAL(DP), INTENT(in)  :: &
!      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
!
!    REAL(DP), INTENT(out)  :: &
!      W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!      W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!      W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!      W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!      W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!      W_udd_332, W_udd_333
!
!
!    REAL(DP) :: &
!      l_udd_111, l_udd_112, l_udd_113, l_udd_121, l_udd_122, & 
!      l_udd_123, l_udd_131, l_udd_132, l_udd_133, l_udd_211, & 
!      l_udd_212, l_udd_213, l_udd_221, l_udd_222, l_udd_223, & 
!      l_udd_231, l_udd_232, l_udd_233, l_udd_311, l_udd_312, & 
!      l_udd_313, l_udd_321, l_udd_322, l_udd_323, l_udd_331, & 
!      l_udd_332, l_udd_333 
!
!    REAL(DP) :: &
!      l_ddd_111, l_ddd_112, l_ddd_113, l_ddd_121, l_ddd_122, & 
!      l_ddd_123, l_ddd_131, l_ddd_132, l_ddd_133, l_ddd_211, & 
!      l_ddd_212, l_ddd_213, l_ddd_221, l_ddd_222, l_ddd_223, & 
!      l_ddd_231, l_ddd_232, l_ddd_233, l_ddd_311, l_ddd_312, & 
!      l_ddd_313, l_ddd_321, l_ddd_322, l_ddd_323, l_ddd_331, & 
!      l_ddd_332, l_ddd_333 
!
!    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
!                k_dd_21, k_dd_31, k_dd_32
!    REAL(DP) :: k_ud_11, k_ud_22, k_ud_33, k_ud_12, k_ud_13, k_ud_23, &
!                k_ud_21, k_ud_31, k_ud_32
!
!    REAL(DP) :: n_u_1, n_u_2, n_u_3
!    REAL(DP) :: V_d_1, V_d_2, V_d_3
!    REAL(DP) :: I_d_1, I_d_2, I_d_3
!    REAL(DP) :: B_d_1, B_d_2, B_d_3
!    REAL(DP) :: W, DT
!
!    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
!               + Gm_dd_22 * V_u_2 * V_u_2 &  
!               + Gm_dd_33 * V_u_3 * V_u_3) )
!
!    B_d_1 = Gm_dd_11 * B_u_1
!    B_d_2 = Gm_dd_22 * B_u_2
!    B_d_3 = Gm_dd_33 * B_u_3
!
!    V_d_1 = Gm_dd_11 * V_u_1
!    V_d_2 = Gm_dd_22 * V_u_2
!    V_d_3 = Gm_dd_33 * V_u_3
!
!    n_u_1 = -B_u_1 / alp    
!    n_u_2 = -B_u_2 / alp    
!    n_u_3 = -B_u_3 / alp    
!
!
!    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )
!
!    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
!          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
!    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
!          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
!    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
!          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )
!
!
!
!
!    CALL ComputeEddingtonTensorComponents_ud &
!    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33 )
!
!    CALL ComputeEddingtonTensorComponents_dd &
!           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!             k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
!             alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3   )
!
!    k_dd_21 = k_dd_12
!    k_dd_31 = k_dd_13
!    k_dd_32 = k_dd_23
!
!    k_ud_21 = k_ud_12
!    k_ud_31 = k_ud_13
!    k_ud_32 = k_ud_23
!
!   CALL ComputeHeatFluxTensorComponents_udd_Lagrangian &
!    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3,         &
!      l_udd_111, l_udd_112, l_udd_113, l_udd_121, l_udd_122, & 
!      l_udd_123, l_udd_131, l_udd_132, l_udd_133, l_udd_211, & 
!      l_udd_212, l_udd_213, l_udd_221, l_udd_222, l_udd_223, & 
!      l_udd_231, l_udd_232, l_udd_233, l_udd_311, l_udd_312, & 
!      l_udd_313, l_udd_321, l_udd_322, l_udd_323, l_udd_331, & 
!      l_udd_332, l_udd_333 )
!
!    CALL ComputeHeatFluxTensorComponents_ddd_Lagrangian &
!      ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!        alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3,         &
!        l_ddd_111, l_ddd_112, l_ddd_113, l_ddd_121, l_ddd_122, & 
!        l_ddd_123, l_ddd_131, l_ddd_132, l_ddd_133, l_ddd_211, & 
!        l_ddd_212, l_ddd_213, l_ddd_221, l_ddd_222, l_ddd_223, & 
!        l_ddd_231, l_ddd_232, l_ddd_233, l_ddd_311, l_ddd_312, & 
!        l_ddd_313, l_ddd_321, l_ddd_322, l_ddd_323, l_ddd_331, & 
!        l_ddd_332, l_ddd_333 )
!
!    n_u_1 = -B_u_1 / alp    
!    n_u_2 = -B_u_2 / alp    
!    n_u_3 = -B_u_3 / alp    
!
!    V_d_1 = Gm_dd_11 * V_u_1
!    V_d_2 = Gm_dd_22 * V_u_2
!    V_d_3 = Gm_dd_33 * V_u_3
!
!    !pulled E**2 out here not sure if its rihgt
!   
!    W_udd_111 = ( W**3 * V_u_1 * V_d_1 * V_d_1 * D & 
!      + W**2 * V_d_1 * V_d_1 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_1 * I_d_1 & 
!      + W**2 * V_u_1 * V_d_1 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_1 * k_dd_11 * D & 
!      + W * V_d_1 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_111 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_111 + V_u_2 * l_ddd_211 + V_u_3 * l_ddd_311 ) * D )
! 
!    W_udd_112 = ( W**3 * V_u_1 * V_d_1 * V_d_2 * D & 
!      + W**2 * V_d_1 * V_d_2 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_1 * I_d_1 & 
!      + W**2 * V_u_1 * V_d_1 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_1 * k_dd_12 * D & 
!      + W * V_d_1 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_112 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_112 + V_u_2 * l_ddd_212 + V_u_3 * l_ddd_312 ) * D ) 
! 
!    W_udd_113 = ( W**3 * V_u_1 * V_d_1 * V_d_3 * D & 
!      + W**2 * V_d_1 * V_d_3 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_1 * I_d_1 & 
!      + W**2 * V_u_1 * V_d_1 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_1 * k_dd_13 * D & 
!      + W * V_d_1 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_113 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_113 + V_u_2 * l_ddd_213 + V_u_3 * l_ddd_313 ) * D ) 
! 
!    W_udd_121 = ( W**3 * V_u_1 * V_d_2 * V_d_1 * D & 
!      + W**2 * V_d_2 * V_d_1 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_1 * I_d_2 & 
!      + W**2 * V_u_1 * V_d_2 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_1 * k_dd_21 * D & 
!      + W * V_d_2 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_121 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_121 + V_u_2 * l_ddd_221 + V_u_3 * l_ddd_321 ) * D ) 
! 
!    W_udd_122 =  ( W**3 * V_u_1 * V_d_2 * V_d_2 * D & 
!      + W**2 * V_d_2 * V_d_2 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_1 * I_d_2 & 
!      + W**2 * V_u_1 * V_d_2 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_1 * k_dd_22 * D & 
!      + W * V_d_2 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_122 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_122 + V_u_2 * l_ddd_222 + V_u_3 * l_ddd_322 ) * D ) 
! 
!    W_udd_123 =  ( W**3 * V_u_1 * V_d_2 * V_d_3 * D & 
!      + W**2 * V_d_2 * V_d_3 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_1 * I_d_2 & 
!      + W**2 * V_u_1 * V_d_2 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_1 * k_dd_23 * D & 
!      + W * V_d_2 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_123 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_123 + V_u_2 * l_ddd_223 + V_u_3 * l_ddd_323 ) * D ) 
! 
!    W_udd_131 =  ( W**3 * V_u_1 * V_d_3 * V_d_1 * D & 
!      + W**2 * V_d_3 * V_d_1 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_1 * I_d_3 & 
!      + W**2 * V_u_1 * V_d_3 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_1 * k_dd_31 * D & 
!      + W * V_d_3 * ( k_ud_11 * D - n_u_1 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_131 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_131 + V_u_2 * l_ddd_231 + V_u_3 * l_ddd_331 ) * D ) 
! 
!    W_udd_132 =  ( W**3 * V_u_1 * V_d_3 * V_d_2 * D & 
!      + W**2 * V_d_3 * V_d_2 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_1 * I_d_3 & 
!      + W**2 * V_u_1 * V_d_3 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_1 * k_dd_32 * D & 
!      + W * V_d_3 * ( k_ud_12 * D - n_u_1 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_132 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_132 + V_u_2 * l_ddd_232 + V_u_3 * l_ddd_332 ) * D ) 
! 
!    W_udd_133 =  ( W**3 * V_u_1 * V_d_3 * V_d_3 * D & 
!      + W**2 * V_d_3 * V_d_3 * ( I_u_1 - n_u_1 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_1 * I_d_3 & 
!      + W**2 * V_u_1 * V_d_3 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_1 * k_dd_33 * D & 
!      + W * V_d_3 * ( k_ud_13 * D - n_u_1 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_133 * D & 
!      + n_u_1 * ( V_u_1 * l_ddd_133 + V_u_2 * l_ddd_233 + V_u_3 * l_ddd_333 ) * D ) 
! 
!    W_udd_211 =  ( W**3 * V_u_2 * V_d_1 * V_d_1 * D & 
!      + W**2 * V_d_1 * V_d_1 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_2 * I_d_1 & 
!      + W**2 * V_u_2 * V_d_1 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_2 * k_dd_11 * D & 
!      + W * V_d_1 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_211 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_111 + V_u_2 * l_ddd_211 + V_u_3 * l_ddd_311 ) * D ) 
! 
!    W_udd_212 =  ( W**3 * V_u_2 * V_d_1 * V_d_2 * D & 
!      + W**2 * V_d_1 * V_d_2 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_2 * I_d_1 & 
!      + W**2 * V_u_2 * V_d_1 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_2 * k_dd_12 * D & 
!      + W * V_d_1 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_212 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_112 + V_u_2 * l_ddd_212 + V_u_3 * l_ddd_312 ) * D ) 
! 
!    W_udd_213 =  ( W**3 * V_u_2 * V_d_1 * V_d_3 * D & 
!      + W**2 * V_d_1 * V_d_3 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_2 * I_d_1 & 
!      + W**2 * V_u_2 * V_d_1 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_2 * k_dd_13 * D & 
!      + W * V_d_1 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_213 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_113 + V_u_2 * l_ddd_213 + V_u_3 * l_ddd_313 ) * D ) 
! 
!    W_udd_221 =  ( W**3 * V_u_2 * V_d_2 * V_d_1 * D & 
!      + W**2 * V_d_2 * V_d_1 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_2 * I_d_2 & 
!      + W**2 * V_u_2 * V_d_2 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_2 * k_dd_21 * D & 
!      + W * V_d_2 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_221 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_121 + V_u_2 * l_ddd_221 + V_u_3 * l_ddd_321 ) * D ) 
! 
!    W_udd_222 = ( W**3 * V_u_2 * V_d_2 * V_d_2 * D & 
!      + W**2 * V_d_2 * V_d_2 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_2 * I_d_2 & 
!      + W**2 * V_u_2 * V_d_2 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_2 * k_dd_22 * D & 
!      + W * V_d_2 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_222 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_122 + V_u_2 * l_ddd_222 + V_u_3 * l_ddd_322 ) * D ) 
! 
!    W_udd_223 = ( W**3 * V_u_2 * V_d_2 * V_d_3 * D & 
!      + W**2 * V_d_2 * V_d_3 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_2 * I_d_2 & 
!      + W**2 * V_u_2 * V_d_2 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_2 * k_dd_23 * D & 
!      + W * V_d_2 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_223 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_123 + V_u_2 * l_ddd_223 + V_u_3 * l_ddd_323 ) * D ) 
! 
!    W_udd_231 = ( W**3 * V_u_2 * V_d_3 * V_d_1 * D & 
!      + W**2 * V_d_3 * V_d_1 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_2 * I_d_3 & 
!      + W**2 * V_u_2 * V_d_3 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_2 * k_dd_31 * D & 
!      + W * V_d_3 * ( k_ud_21 * D - n_u_2 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_231 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_131 + V_u_2 * l_ddd_231 + V_u_3 * l_ddd_331 ) * D ) 
! 
!    W_udd_232 =  ( W**3 * V_u_2 * V_d_3 * V_d_2 * D & 
!      + W**2 * V_d_3 * V_d_2 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_2 * I_d_3 & 
!      + W**2 * V_u_2 * V_d_3 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_2 * k_dd_32 * D & 
!      + W * V_d_3 * ( k_ud_22 * D - n_u_2 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_232 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_132 + V_u_2 * l_ddd_232 + V_u_3 * l_ddd_332 ) * D ) 
! 
!    W_udd_233 =  ( W**3 * V_u_2 * V_d_3 * V_d_3 * D & 
!      + W**2 * V_d_3 * V_d_3 * ( I_u_2 - n_u_2 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_2 * I_d_3 & 
!      + W**2 * V_u_2 * V_d_3 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_2 * k_dd_33 * D & 
!      + W * V_d_3 * ( k_ud_23 * D - n_u_2 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_233 * D & 
!      + n_u_2 * ( V_u_1 * l_ddd_133 + V_u_2 * l_ddd_233 + V_u_3 * l_ddd_333 ) * D ) 
! 
!    W_udd_311 =  ( W**3 * V_u_3 * V_d_1 * V_d_1 * D & 
!      + W**2 * V_d_1 * V_d_1 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_3 * I_d_1 & 
!      + W**2 * V_u_3 * V_d_1 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_3 * k_dd_11 * D & 
!      + W * V_d_1 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_311 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_111 + V_u_2 * l_ddd_211 + V_u_3 * l_ddd_311 ) * D ) 
! 
!    W_udd_312 = ( W**3 * V_u_3 * V_d_1 * V_d_2 * D & 
!      + W**2 * V_d_1 * V_d_2 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_3 * I_d_1 & 
!      + W**2 * V_u_3 * V_d_1 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_3 * k_dd_12 * D & 
!      + W * V_d_1 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_312 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_112 + V_u_2 * l_ddd_212 + V_u_3 * l_ddd_312 ) * D ) 
! 
!    W_udd_313 =  ( W**3 * V_u_3 * V_d_1 * V_d_3 * D & 
!      + W**2 * V_d_1 * V_d_3 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_3 * I_d_1 & 
!      + W**2 * V_u_3 * V_d_1 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + W * V_u_3 * k_dd_13 * D & 
!      + W * V_d_1 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_313 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_113 + V_u_2 * l_ddd_213 + V_u_3 * l_ddd_313 ) * D ) 
! 
!    W_udd_321 = ( W**3 * V_u_3 * V_d_2 * V_d_1 * D & 
!      + W**2 * V_d_2 * V_d_1 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_3 * I_d_2 & 
!      + W**2 * V_u_3 * V_d_2 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_3 * k_dd_21 * D & 
!      + W * V_d_2 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_321 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_121 + V_u_2 * l_ddd_221 + V_u_3 * l_ddd_321 ) * D ) 
! 
!    W_udd_322 = ( W**3 * V_u_3 * V_d_2 * V_d_2 * D & 
!      + W**2 * V_d_2 * V_d_2 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_3 * I_d_2 & 
!      + W**2 * V_u_3 * V_d_2 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_3 * k_dd_22 * D & 
!      + W * V_d_2 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_322 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_122 + V_u_2 * l_ddd_222 + V_u_3 * l_ddd_322 ) * D ) 
! 
!    W_udd_323 =  ( W**3 * V_u_3 * V_d_2 * V_d_3 * D & 
!      + W**2 * V_d_2 * V_d_3 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_3 * I_d_2 & 
!      + W**2 * V_u_3 * V_d_2 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + W * V_u_3 * k_dd_23 * D & 
!      + W * V_d_2 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_323 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_123 + V_u_2 * l_ddd_223 + V_u_3 * l_ddd_323 ) * D ) 
! 
!    W_udd_331 = ( W**3 * V_u_3 * V_d_3 * V_d_1 * D & 
!      + W**2 * V_d_3 * V_d_1 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_1 * V_u_3 * I_d_3 & 
!      + W**2 * V_u_3 * V_d_3 * I_d_1 & 
!      + W * V_d_1 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_3 * k_dd_31 * D & 
!      + W * V_d_3 * ( k_ud_31 * D - n_u_3 * ( V_u_1 * k_dd_11 + V_u_2 * k_dd_21 + V_u_3 * k_dd_31 ) * D ) & 
!      + l_udd_331 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_131 + V_u_2 * l_ddd_231 + V_u_3 * l_ddd_331 ) * D ) 
! 
!    W_udd_332 =  ( W**3 * V_u_3 * V_d_3 * V_d_2 * D & 
!      + W**2 * V_d_3 * V_d_2 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_2 * V_u_3 * I_d_3 & 
!      + W**2 * V_u_3 * V_d_3 * I_d_2 & 
!      + W * V_d_2 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_3 * k_dd_32 * D & 
!      + W * V_d_3 * ( k_ud_32 * D - n_u_3 * ( V_u_1 * k_dd_12 + V_u_2 * k_dd_22 + V_u_3 * k_dd_32 ) * D ) & 
!      + l_udd_332 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_132 + V_u_2 * l_ddd_232 + V_u_3 * l_ddd_332 ) * D ) 
! 
!    W_udd_333 =  ( W**3 * V_u_3 * V_d_3 * V_d_3 * D & 
!      + W**2 * V_d_3 * V_d_3 * ( I_u_3 - n_u_3 * (V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3)) & 
!      + W**2 * V_d_3 * V_u_3 * I_d_3 & 
!      + W**2 * V_u_3 * V_d_3 * I_d_3 & 
!      + W * V_d_3 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + W * V_u_3 * k_dd_33 * D & 
!      + W * V_d_3 * ( k_ud_33 * D - n_u_3 * ( V_u_1 * k_dd_13 + V_u_2 * k_dd_23 + V_u_3 * k_dd_33 ) * D ) & 
!      + l_udd_333 * D & 
!      + n_u_3 * ( V_u_1 * l_ddd_133 + V_u_2 * l_ddd_233 + V_u_3 * l_ddd_333 ) * D ) 
! 
!
!  END SUBROUTINE ComputeHeatFluxTensorComponents_udd_Eulerian
!
!  SUBROUTINE ComputeXYZ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                         alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, &
!                         X, Y_d_1, Y_d_2, Y_d_3, Z_ud_11, Z_ud_12, Z_ud_13, Z_ud_21, Z_ud_22, & 
!                         Z_ud_23, Z_ud_31, Z_ud_32, Z_ud_33 )
!
!    REAL(DP), INTENT(in)  :: &
!      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
!
!    REAL(DP), INTENT(out)  :: &
!      X, Y_d_1, Y_d_2, Y_d_3, Z_ud_11, Z_ud_12, Z_ud_13, Z_ud_21, Z_ud_22, &
!      Z_ud_23, Z_ud_31, Z_ud_32, Z_ud_33 
!
!    REAL(DP) &
!      EP, F_u_1, F_u_2, F_u_3, F_d_1, F_d_2, F_d_3
!    REAL(DP) :: &
!      S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32
!    REAL(DP)  :: &
!      S_dd_11, S_dd_22, S_dd_33, S_dd_12, S_dd_13, S_dd_23, S_dd_21, S_dd_31, S_dd_32
!    REAL(DP)  :: &
!      W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!      W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!      W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!      W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!      W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!      W_udd_332, W_udd_333
!    REAL(DP)  :: &
!      W_ddd_111, W_ddd_112, W_ddd_113, W_ddd_121, W_ddd_122, & 
!      W_ddd_123, W_ddd_131, W_ddd_132, W_ddd_133, W_ddd_211, & 
!      W_ddd_212, W_ddd_213, W_ddd_221, W_ddd_222, W_ddd_223, & 
!      W_ddd_231, W_ddd_232, W_ddd_233, W_ddd_311, W_ddd_312, & 
!      W_ddd_313, W_ddd_321, W_ddd_322, W_ddd_323, W_ddd_331, & 
!      W_ddd_332, W_ddd_333
!
!    REAL(DP) :: W
!
!    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
!               + Gm_dd_22 * V_u_2 * V_u_2 &  
!               + Gm_dd_33 * V_u_3 * V_u_3) )
!
!    CALL ComputeStressEnergyComponents_Eulerian &
!                        ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3, &
!                         alp, B_u_1, B_u_2, B_u_3,EP, F_u_1, F_u_2, F_u_3, &
!                         S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32 )
!
!    CALL ComputeHeatFluxTensorComponents_udd_Eulerian &
!                                          ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                                           alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, &
!                                           W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!                                           W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!                                           W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!                                           W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!                                           W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!                                           W_udd_332, W_udd_333)
!
!    F_d_1 = Gm_dd_11 * F_u_1
!
!    F_d_2 = Gm_dd_22 * F_u_2
!
!    F_d_3 = Gm_dd_33 * F_u_3
!
!    S_dd_11 =  Gm_dd_11 * S_ud_11 
! 
!    S_dd_12 =  Gm_dd_11 * S_ud_12 
! 
!    S_dd_13 =  Gm_dd_11 * S_ud_13 
! 
!    S_dd_21 =  Gm_dd_22 * S_ud_21 
! 
!    S_dd_22 =  Gm_dd_22 * S_ud_22 
! 
!    S_dd_23 =  Gm_dd_22 * S_ud_23 
! 
!    S_dd_31 =  Gm_dd_33 * S_ud_31 
! 
!    S_dd_32 =  Gm_dd_33 * S_ud_32 
! 
!    S_dd_33 =  Gm_dd_33 * S_ud_33 
!
!    W_ddd_111 =  Gm_dd_11 * W_udd_111 
! 
!    W_ddd_112 =  Gm_dd_11 * W_udd_112 
! 
!    W_ddd_113 =  Gm_dd_11 * W_udd_113 
! 
!    W_ddd_121 =  Gm_dd_11 * W_udd_121 
! 
!    W_ddd_122 =  Gm_dd_11 * W_udd_122 
! 
!    W_ddd_123 =  Gm_dd_11 * W_udd_123 
! 
!    W_ddd_131 =  Gm_dd_11 * W_udd_131 
! 
!    W_ddd_132 =  Gm_dd_11 * W_udd_132 
! 
!    W_ddd_133 =  Gm_dd_11 * W_udd_133 
! 
!    W_ddd_211 =  Gm_dd_22 * W_udd_211 
! 
!    W_ddd_212 =  Gm_dd_22 * W_udd_212 
! 
!    W_ddd_213 =  Gm_dd_22 * W_udd_213 
! 
!    W_ddd_221 =  Gm_dd_22 * W_udd_221 
! 
!    W_ddd_222 =  Gm_dd_22 * W_udd_222 
! 
!    W_ddd_223 =  Gm_dd_22 * W_udd_223 
! 
!    W_ddd_231 =  Gm_dd_22 * W_udd_231 
! 
!    W_ddd_232 =  Gm_dd_22 * W_udd_232 
! 
!    W_ddd_233 =  Gm_dd_22 * W_udd_233 
! 
!    W_ddd_311 =  Gm_dd_33 * W_udd_311 
! 
!    W_ddd_312 =  Gm_dd_33 * W_udd_312 
! 
!    W_ddd_313 =  Gm_dd_33 * W_udd_313 
! 
!    W_ddd_321 =  Gm_dd_33 * W_udd_321 
! 
!    W_ddd_322 =  Gm_dd_33 * W_udd_322 
! 
!    W_ddd_323 =  Gm_dd_33 * W_udd_323 
! 
!    W_ddd_331 =  Gm_dd_33 * W_udd_331 
! 
!    W_ddd_332 =  Gm_dd_33 * W_udd_332 
! 
!    W_ddd_333 =  Gm_dd_33 * W_udd_333 
!
!!pulled out Energy here not sure if that was correct
!    X = EP + ( V_u_1 * F_d_1 + V_u_2 * F_d_2 + V_u_3 * F_d_3 ) &
!      +  ( V_u_1 * V_u_1 * S_dd_11 + V_u_1 * V_u_2 * S_dd_12 + V_u_1 * V_u_3 * S_dd_13 & 
!      + V_u_2 * V_u_1 * S_dd_21 + V_u_2 * V_u_2 * S_dd_22 + V_u_2 * V_u_3 * S_dd_23 & 
!      + V_u_3 * V_u_1 * S_dd_31 + V_u_3 * V_u_2 * S_dd_32 + V_u_3 * V_u_3 * S_dd_33 ) &
!      +  W * ( V_u_1 * V_u_1 * V_u_1 * W_ddd_111 + V_u_1 * V_u_1 * V_u_2 * W_ddd_112 + V_u_1 * V_u_1 * V_u_3 * W_ddd_113 & 
!      + V_u_1 * V_u_2 * V_u_1 * W_ddd_121 + V_u_1 * V_u_2 * V_u_2 * W_ddd_122 + V_u_1 * V_u_2 * V_u_3 * W_ddd_123 & 
!      + V_u_1 * V_u_3 * V_u_1 * W_ddd_131 + V_u_1 * V_u_3 * V_u_2 * W_ddd_132 + V_u_1 * V_u_3 * V_u_3 * W_ddd_133 & 
!      + V_u_2 * V_u_1 * V_u_1 * W_ddd_211 + V_u_2 * V_u_1 * V_u_2 * W_ddd_212 + V_u_2 * V_u_1 * V_u_3 * W_ddd_213 & 
!      + V_u_2 * V_u_2 * V_u_1 * W_ddd_221 + V_u_2 * V_u_2 * V_u_2 * W_ddd_222 + V_u_2 * V_u_2 * V_u_3 * W_ddd_223 & 
!      + V_u_2 * V_u_3 * V_u_1 * W_ddd_231 + V_u_2 * V_u_3 * V_u_2 * W_ddd_232 + V_u_2 * V_u_3 * V_u_3 * W_ddd_233 & 
!      + V_u_3 * V_u_1 * V_u_1 * W_ddd_311 + V_u_3 * V_u_1 * V_u_2 * W_ddd_312 + V_u_3 * V_u_1 * V_u_3 * W_ddd_313 & 
!      + V_u_3 * V_u_2 * V_u_1 * W_ddd_321 + V_u_3 * V_u_2 * V_u_2 * W_ddd_322 + V_u_3 * V_u_2 * V_u_3 * W_ddd_323 & 
!      + V_u_3 * V_u_3 * V_u_1 * W_ddd_331 + V_u_3 * V_u_3 * V_u_2 * W_ddd_332 + V_u_3 * V_u_3 * V_u_3 * W_ddd_333 )  
!    X = X / W
!
!    Y_d_1 = F_d_1 + ( V_u_1 * S_dd_11 + V_u_2 * S_dd_12 + V_u_3 * S_dd_13 ) &
!          + W * ( V_u_1 * V_u_1 * W_ddd_111 + V_u_1 * V_u_2 * W_ddd_112 + V_u_1 * V_u_3 * W_ddd_113 & 
!          + V_u_2 * V_u_1 * W_ddd_121 + V_u_2 * V_u_2 * W_ddd_122 + V_u_2 * V_u_3 * W_ddd_123 & 
!          + V_u_3 * V_u_1 * W_ddd_131 + V_u_3 * V_u_2 * W_ddd_132 + V_u_3 * V_u_3 * W_ddd_133 )
!    Y_d_1 = Y_d_1 / W
!
!    Y_d_2 = F_d_2 + ( V_u_1 * S_dd_21 + V_u_2 * S_dd_22 + V_u_3 * S_dd_23 ) &
!          + W * ( V_u_1 * V_u_1 * W_ddd_211 + V_u_1 * V_u_2 * W_ddd_212 + V_u_1 * V_u_3 * W_ddd_213 & 
!          + V_u_2 * V_u_1 * W_ddd_221 + V_u_2 * V_u_2 * W_ddd_222 + V_u_2 * V_u_3 * W_ddd_223 & 
!          + V_u_3 * V_u_1 * W_ddd_231 + V_u_3 * V_u_2 * W_ddd_232 + V_u_3 * V_u_3 * W_ddd_233 )
!    Y_d_2 = Y_d_2 / W
!
!    Y_d_3 = F_d_3 + ( V_u_1 * S_dd_31 + V_u_2 * S_dd_32 + V_u_3 * S_dd_33 ) &
!          + W * ( V_u_1 * V_u_1 * W_ddd_311 + V_u_1 * V_u_2 * W_ddd_312 + V_u_1 * V_u_3 * W_ddd_313 & 
!          + V_u_2 * V_u_1 * W_ddd_321 + V_u_2 * V_u_2 * W_ddd_322 + V_u_2 * V_u_3 * W_ddd_323 & 
!          + V_u_3 * V_u_1 * W_ddd_331 + V_u_3 * V_u_2 * W_ddd_332 + V_u_3 * V_u_3 * W_ddd_333 )
!    Y_d_3 = Y_d_3 / W
!
!    Z_ud_11 = S_ud_11 + W * ( V_u_1 * W_udd_111 + V_u_2 * W_udd_112 + V_u_3 * W_udd_113 )
!    Z_ud_11 = Z_ud_11 / W
!
!    Z_ud_12 = S_ud_12 + W * ( V_u_1 * W_udd_121 + V_u_2 * W_udd_122 + V_u_3 * W_udd_123 )
!    Z_ud_12 = Z_ud_12 / W
!
!    Z_ud_13 = S_ud_13 + W * ( V_u_1 * W_udd_131 + V_u_2 * W_udd_132 + V_u_3 * W_udd_133 )
!    Z_ud_13 = Z_ud_13 / W
!
!    Z_ud_21 = S_ud_21 + W * ( V_u_1 * W_udd_211 + V_u_2 * W_udd_212 + V_u_3 * W_udd_213 )
!    Z_ud_21 = Z_ud_21 / W
!
!    Z_ud_22 = S_ud_22 + W * ( V_u_1 * W_udd_221 + V_u_2 * W_udd_222 + V_u_3 * W_udd_223 )
!    Z_ud_22 = Z_ud_22 / W
!
!    Z_ud_23 = S_ud_23 + W * ( V_u_1 * W_udd_231 + V_u_2 * W_udd_232 + V_u_3 * W_udd_233 )
!    Z_ud_23 = Z_ud_23 / W
!
!    Z_ud_31 = S_ud_31 + W * ( V_u_1 * W_udd_311 + V_u_2 * W_udd_312 + V_u_3 * W_udd_313 )
!    Z_ud_31 = Z_ud_31 / W
!
!    Z_ud_32 = S_ud_32 + W * ( V_u_1 * W_udd_321 + V_u_2 * W_udd_322 + V_u_3 * W_udd_323 )
!    Z_ud_32 = Z_ud_32 / W
!
!    Z_ud_33 = S_ud_33 + W * ( V_u_1 * W_udd_331 + V_u_2 * W_udd_332 + V_u_3 * W_udd_333 )
!    Z_ud_33 = Z_ud_33 / W
!
!
!
!  END SUBROUTINE 
!
!



  FUNCTION Flux_X1( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X1(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
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
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )

    Flux_X1(1) = alp * I_u_1 + ( alp * V_u_1 - B_u_1 ) * W * D
    Flux_X1(2) = alp * k_ud_ij(1,1) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_1
    Flux_X1(3) = alp * k_ud_ij(1,2) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_2
    Flux_X1(4) = alp * k_ud_ij(1,3) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X1

  FUNCTION Flux_X2( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X2(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
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
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )

    Flux_X2(1) = alp * I_u_2 + ( alp * V_u_2 - B_u_2 ) * W * D
    Flux_X2(2) = alp * k_ud_ij(1,2) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_1
    Flux_X2(3) = alp * k_ud_ij(2,2) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_2
    Flux_X2(4) = alp * k_ud_ij(2,3) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X2

  FUNCTION Flux_X3( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X3(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
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
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )


    Flux_X3(1) = alp * I_u_3 + ( alp * V_u_3 - B_u_3 ) * W * D
    Flux_X3(2) = alp * k_ud_ij(1,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_1
    Flux_X3(3) = alp * k_ud_ij(2,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_2
    Flux_X3(4) = alp * k_ud_ij(3,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X3

  FUNCTION Flux_E( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3, U, dU_dX0, dU_dX1, dU_dX2, dU_dX3)

! negative sign has been incorperated into Flux_E

    REAL(DP)             :: Flux_E(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), INTENT(in) :: U(0:3), dU_dX0(0:3), dU_dX1(0:3), dU_dX2(0:3), dU_dX3(0:3)

    REAL(DP) :: V_0, B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: I(0:3), k_uu_munu(0:3,0:3), k_ud_munu(0:3,0:3), l_uud_munurho(0:3,0:3,0:3), dU_dX(0:3,0:3)
    INTEGER :: mu, nu

    CALL ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )

    CALL ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )

    B_d_1 = Gm_dd_11 * B_u_1 
    B_d_2 = Gm_dd_22 * B_u_2 
    B_d_3 = Gm_dd_33 * B_u_3 

    V_d_1 = Gm_dd_11 * V_u_1 
    V_d_2 = Gm_dd_22 * V_u_2 
    V_d_3 = Gm_dd_33 * V_u_3
 
    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3   

    I(1) = I_u_1
    I(2) = I_u_2
    I(3) = I_u_3
    I(0) = 1.0_DP / ( 1.0_DP - V_0 / alp ) * ( V_d_1 * I(1) + V_d_2 * I(2) + V_d_3 * I(3) ) / alp



    dU_dX(0,0) = dU_dX0(0) 
    dU_dX(0,1) = dU_dX0(1) 
    dU_dX(0,2) = dU_dX0(2) 
    dU_dX(0,3) = dU_dX0(3) 

    dU_dX(1,0) = dU_dX1(0) 
    dU_dX(1,1) = dU_dX1(1) 
    dU_dX(1,2) = dU_dX1(2) 
    dU_dX(1,3) = dU_dX1(3) 

    dU_dX(2,0) = dU_dX2(0) 
    dU_dX(2,1) = dU_dX2(1) 
    dU_dX(2,2) = dU_dX2(2) 
    dU_dX(2,3) = dU_dX2(3)
 
    dU_dX(3,0) = dU_dX3(0) 
    dU_dX(3,1) = dU_dX3(1) 
    dU_dX(3,2) = dU_dX3(2) 
    dU_dX(3,3) = dU_dX3(3) 


    Flux_E = 0.0_DP
   
    DO mu = 0,3
    DO nu = 0,3

      Flux_E(1) = Flux_E(1) + ( I(nu) * U(mu) + k_uu_munu(mu,nu) * D ) * dU_dX(mu,nu) 
      Flux_E(2) = Flux_E(2) + ( k_ud_munu(nu,1) * D * U(mu) + l_uud_munurho(mu,nu,1) * D ) * dU_dX(mu,nu)
      Flux_E(3) = Flux_E(3) + ( k_ud_munu(nu,2) * D * U(mu) + l_uud_munurho(mu,nu,2) * D ) * dU_dX(mu,nu)
      Flux_E(4) = Flux_E(4) + ( k_ud_munu(nu,3) * D * U(mu) + l_uud_munurho(mu,nu,3) * D ) * dU_dX(mu,nu)

    END DO
    END DO


    Flux_E = - Flux_E

    RETURN

  END FUNCTION Flux_E

  FUNCTION Source_E( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3, U_u, U_d, dU_dX0, dU_dX1, dU_dX2, dU_dX3)

! negative sign has been incorperated into Flux_E

    REAL(DP)             :: Source_E(3)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), INTENT(in) :: U_u(0:3), U_d(0:3), dU_dX0(0:3), dU_dX1(0:3), dU_dX2(0:3), dU_dX3(0:3)

    REAL(DP) :: V_0, B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: I(0:3), k_uu_munu(0:3,0:3), k_ud_munu(0:3,0:3), l_uud_munurho(0:3,0:3,0:3), dU_dX(0:3,0:3)
    INTEGER :: mu, nu, j

    CALL ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )

    CALL ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )



    B_d_1 = Gm_dd_11 * B_u_1 
    B_d_2 = Gm_dd_22 * B_u_2 
    B_d_3 = Gm_dd_33 * B_u_3 

    V_d_1 = Gm_dd_11 * V_u_1 
    V_d_2 = Gm_dd_22 * V_u_2 
    V_d_3 = Gm_dd_33 * V_u_3
 
    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3   

    I(1) = I_u_1
    I(2) = I_u_2
    I(3) = I_u_3
    I(0) = 1.0_DP / ( 1.0_DP - V_0 / alp ) * ( V_d_1 * I(1) + V_d_2 * I(2) + V_d_3 * I(3) ) / alp

    dU_dX(0,0) = dU_dX0(0) 
    dU_dX(0,1) = dU_dX0(1) 
    dU_dX(0,2) = dU_dX0(2) 
    dU_dX(0,3) = dU_dX0(3) 

    dU_dX(1,0) = dU_dX1(0) 
    dU_dX(1,1) = dU_dX1(1) 
    dU_dX(1,2) = dU_dX1(2) 
    dU_dX(1,3) = dU_dX1(3) 

    dU_dX(2,0) = dU_dX2(0) 
    dU_dX(2,1) = dU_dX2(1) 
    dU_dX(2,2) = dU_dX2(2) 
    dU_dX(2,3) = dU_dX2(3)
 
    dU_dX(3,0) = dU_dX3(0) 
    dU_dX(3,1) = dU_dX3(1) 
    dU_dX(3,2) = dU_dX3(2) 
    dU_dX(3,3) = dU_dX3(3) 

    Source_E = 0.0_DP



    DO j = 1,3
    DO mu = 0,3
    DO nu = 0,3

      Source_E(j) = Source_E(j) - ( I(nu) * U_d(j) * U_u(mu) + l_uud_munurho(mu,nu,j) * D &
                    + k_ud_munu(nu,j) * D * U_u(mu) + k_uu_munu(nu,mu) * D * U_d(j) ) * dU_dX(mu,nu)

    END DO
    END DO
    END DO

    DO j = 1,3
    DO mu = 0,3

      Source_E(j) = Source_E(j) + ( D * U_u(mu) + I(mu) ) * dU_dX(mu,j)

    END DO
    END DO

    RETURN

  END FUNCTION Source_E

!  FUNCTION Flux_E( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                      alp, B_u_1, B_u_2, B_u_3, dW_dX1, dW_dX2, dW_dX3, &
!                      dWV_u_1_dX1, dWV_u_2_dX1, dWV_u_3_dX1, &
!                      dWV_u_1_dX2, dWV_u_2_dX2, dWV_u_3_dX2, &
!                      dWV_u_1_dX3, dWV_u_2_dX3, dWV_u_3_dX3 )
!! negative sign has been incorperated into Flux_E
!
!    REAL(DP)             :: Flux_E(4)
!    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
!    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
!    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
!    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
!    REAL(DP), INTENT(in) :: dW_dX1, dW_dX2, dW_dX3
!    REAL(DP), INTENT(in) :: dWV_u_1_dX1, dWV_u_2_dX1, dWV_u_3_dX1, &
!                            dWV_u_1_dX2, dWV_u_2_dX2, dWV_u_3_dX2, &
!                            dWV_u_1_dX3, dWV_u_2_dX3, dWV_u_3_dX3 
!
!
!
!    REAL(DP) :: &
!      EP, F_u_1, F_u_2, F_u_3
!    REAL(DP) :: &
!      S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32
!    REAL(DP) :: &
!      X, Y_d_1, Y_d_2, Y_d_3, Z_ud_11, Z_ud_12, Z_ud_13, Z_ud_21, Z_ud_22, &
!      Z_ud_23, Z_ud_31, Z_ud_32, Z_ud_33 
!    REAL(DP) :: &
!      W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!      W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!      W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!      W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!      W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!      W_udd_332, W_udd_333
!
!    REAL(DP) :: &
!      YT_d_1, YT_d_2, YT_d_3, ZT_ud_11, ZT_ud_12, ZT_ud_13, ZT_ud_21, ZT_ud_22, &
!      ZT_ud_23, ZT_ud_31, ZT_ud_32, ZT_ud_33 
!
!    REAL(DP) :: &
!      WT_udd_111, WT_udd_112, WT_udd_113, WT_udd_121, WT_udd_122, & 
!      WT_udd_123, WT_udd_131, WT_udd_132, WT_udd_133, WT_udd_211, & 
!      WT_udd_212, WT_udd_213, WT_udd_221, WT_udd_222, WT_udd_223, & 
!      WT_udd_231, WT_udd_232, WT_udd_233, WT_udd_311, WT_udd_312, & 
!      WT_udd_313, WT_udd_321, WT_udd_322, WT_udd_323, WT_udd_331, & 
!      WT_udd_332, WT_udd_333
!
!    REAL(DP) :: W, V_d_1, V_d_2, V_d_3
!    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
!               + Gm_dd_22 * V_u_2 * V_u_2 &  
!               + Gm_dd_33 * V_u_3 * V_u_3) )
!
!    V_d_1 = Gm_dd_11 * V_u_1
!    V_d_2 = Gm_dd_22 * V_u_2
!    V_d_3 = Gm_dd_33 * V_u_3
!
!
!!have XYZ have W as variables so I dont have to calculate them mroe than once
!    CALL ComputeXYZ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                         alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, &
!                         X, Y_d_1, Y_d_2, Y_d_3, Z_ud_11, Z_ud_12, Z_ud_13, Z_ud_21, Z_ud_22, & 
!                         Z_ud_23, Z_ud_31, Z_ud_32, Z_ud_33 )
!
!    CALL ComputeStressEnergyComponents_Eulerian &
!                    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3, &
!                      alp, B_u_1, B_u_2, B_u_3, EP, F_u_1, F_u_2, F_u_3, &
!                      S_ud_11, S_ud_22, S_ud_33, S_ud_12, S_ud_13, S_ud_23, S_ud_21, S_ud_31, S_ud_32 )
!
!    CALL ComputeHeatFluxTensorComponents_udd_Eulerian &
!                    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, &
!                      W_udd_111, W_udd_112, W_udd_113, W_udd_121, W_udd_122, & 
!                      W_udd_123, W_udd_131, W_udd_132, W_udd_133, W_udd_211, & 
!                      W_udd_212, W_udd_213, W_udd_221, W_udd_222, W_udd_223, & 
!                      W_udd_231, W_udd_232, W_udd_233, W_udd_311, W_udd_312, & 
!                      W_udd_313, W_udd_321, W_udd_322, W_udd_323, W_udd_331, & 
!                      W_udd_332, W_udd_333) 
!
!
!    YT_d_1 = Y_d_1 - W * V_d_1 * EP
!    YT_d_2 = Y_d_2 - W * V_d_2 * EP
!    YT_d_3 = Y_d_3 - W * V_d_3 * EP
!
!    ZT_ud_11 = Z_ud_11 - W * V_d_1 * F_u_1
!    ZT_ud_12 = Z_ud_12 - W * V_d_2 * F_u_1
!    ZT_ud_13 = Z_ud_13 - W * V_d_3 * F_u_1
!    ZT_ud_21 = Z_ud_21 - W * V_d_1 * F_u_2
!    ZT_ud_22 = Z_ud_22 - W * V_d_2 * F_u_2
!    ZT_ud_23 = Z_ud_23 - W * V_d_3 * F_u_2
!    ZT_ud_31 = Z_ud_31 - W * V_d_1 * F_u_3
!    ZT_ud_32 = Z_ud_32 - W * V_d_2 * F_u_3
!    ZT_ud_33 = Z_ud_33 - W * V_d_3 * F_u_3
!
!    WT_udd_111 = W_udd_111 - W * V_d_1 * S_ud_11
!    WT_udd_112 = W_udd_112 - W * V_d_1 * S_ud_12
!    WT_udd_113 = W_udd_113 - W * V_d_1 * S_ud_13
!    WT_udd_121 = W_udd_121 - W * V_d_2 * S_ud_11
!    WT_udd_122 = W_udd_122 - W * V_d_2 * S_ud_12
!    WT_udd_123 = W_udd_123 - W * V_d_2 * S_ud_13
!    WT_udd_131 = W_udd_131 - W * V_d_3 * S_ud_11
!    WT_udd_132 = W_udd_132 - W * V_d_3 * S_ud_12
!    WT_udd_133 = W_udd_133 - W * V_d_3 * S_ud_13
!    WT_udd_211 = W_udd_211 - W * V_d_1 * S_ud_21
!    WT_udd_212 = W_udd_212 - W * V_d_1 * S_ud_22
!    WT_udd_213 = W_udd_213 - W * V_d_1 * S_ud_23
!    WT_udd_221 = W_udd_221 - W * V_d_2 * S_ud_21
!    WT_udd_222 = W_udd_222 - W * V_d_2 * S_ud_22
!    WT_udd_223 = W_udd_223 - W * V_d_2 * S_ud_23
!    WT_udd_231 = W_udd_231 - W * V_d_3 * S_ud_21
!    WT_udd_232 = W_udd_232 - W * V_d_3 * S_ud_22
!    WT_udd_233 = W_udd_233 - W * V_d_3 * S_ud_23
!    WT_udd_311 = W_udd_311 - W * V_d_1 * S_ud_31
!    WT_udd_312 = W_udd_312 - W * V_d_1 * S_ud_32
!    WT_udd_313 = W_udd_313 - W * V_d_1 * S_ud_33
!    WT_udd_321 = W_udd_321 - W * V_d_2 * S_ud_31
!    WT_udd_322 = W_udd_322 - W * V_d_2 * S_ud_32
!    WT_udd_323 = W_udd_323 - W * V_d_2 * S_ud_33
!    WT_udd_331 = W_udd_331 - W * V_d_3 * S_ud_31
!    WT_udd_332 = W_udd_332 - W * V_d_3 * S_ud_32
!    WT_udd_333 = W_udd_333 - W * V_d_3 * S_ud_33
!
!    Flux_E(1) = - ( F_u_1 * dW_dX1 + F_u_2 * dW_dX2 + F_u_3 * dW_dX3 ) &
!               + ( S_ud_11 * dWV_u_1_dX1 + S_ud_12 * dWV_u_2_dX1 + S_ud_13 * dWV_u_3_dX1 & 
!               + S_ud_21 * dWV_u_1_dX2 + S_ud_22 * dWV_u_2_dX2 + S_ud_23 * dWV_u_3_dX2 & 
!               + S_ud_31 * dWV_u_1_dX3 + S_ud_32 * dWV_u_2_dX3 + S_ud_33 * dWV_u_3_dX3 ) 
!    Flux_E(1) = -Flux_E(1)
!
!
!    Flux_E(2) = - ( alp * ZT_ud_11 - YT_d_1 * B_u_1 ) * dW_dX1 &
!                - ( alp * ZT_ud_21 - YT_d_1 * B_u_2 ) * dW_dX2 &
!                - ( alp * ZT_ud_31 - YT_d_1 * B_u_3 ) * dW_dX3 &
!                + alp * WT_udd_111 * dWV_u_1_dX1 + alp * WT_udd_112 * dWV_u_2_dX1& 
!                + alp * WT_udd_113 * dWV_u_3_dX1 + alp * WT_udd_211 * dWV_u_1_dX2& 
!                + alp * WT_udd_212 * dWV_u_2_dX2 + alp * WT_udd_213 * dWV_u_3_dX2& 
!                + alp * WT_udd_311 * dWV_u_1_dX3 + alp * WT_udd_312 * dWV_u_2_dX3& 
!                + alp * WT_udd_313 * dWV_u_3_dX3
!    Flux_E(2) = -Flux_E(2) / alp
!
!
!
!    Flux_E(3) = - ( alp * ZT_ud_12 - YT_d_2 * B_u_1 ) * dW_dX1 &
!                - ( alp * ZT_ud_22 - YT_d_2 * B_u_2 ) * dW_dX2 &
!                - ( alp * ZT_ud_32 - YT_d_2 * B_u_3 ) * dW_dX3 &
!                + alp * WT_udd_121 * dWV_u_1_dX1 + alp * WT_udd_122 * dWV_u_2_dX1& 
!                + alp * WT_udd_123 * dWV_u_3_dX1 + alp * WT_udd_221 * dWV_u_1_dX2& 
!                + alp * WT_udd_222 * dWV_u_2_dX2 + alp * WT_udd_223 * dWV_u_3_dX2& 
!                + alp * WT_udd_321 * dWV_u_1_dX3 + alp * WT_udd_322 * dWV_u_2_dX3& 
!                + alp * WT_udd_323 * dWV_u_3_dX3
!    Flux_E(3) = -Flux_E(3) / alp
!
!    Flux_E(4) = - ( alp * ZT_ud_13 - YT_d_3 * B_u_1 ) * dW_dX1 &
!                - ( alp * ZT_ud_23 - YT_d_3 * B_u_2 ) * dW_dX2 &
!                - ( alp * ZT_ud_33 - YT_d_3 * B_u_3 ) * dW_dX3 &
!                + alp * WT_udd_131 * dWV_u_1_dX1 + alp * WT_udd_132 * dWV_u_2_dX1& 
!                + alp * WT_udd_133 * dWV_u_3_dX1 + alp * WT_udd_231 * dWV_u_1_dX2& 
!                + alp * WT_udd_232 * dWV_u_2_dX2 + alp * WT_udd_233 * dWV_u_3_dX2& 
!                + alp * WT_udd_331 * dWV_u_1_dX3 + alp * WT_udd_332 * dWV_u_2_dX3& 
!                + alp * WT_udd_333 * dWV_u_3_dX3
!
!    Flux_E(4) = -Flux_E(4) / alp
!
!    RETURN
!
!  END FUNCTION Flux_E
!

  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF

!  FUNCTION Flux_G( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
!                      alp, B_u_1, B_u_2, B_u_3 )
!
!    REAL(DP)             :: Flux_G(3)
!    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
!    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
!    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
!    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
!
!    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3, vI
!
!    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
!               + Gm_dd_22 * V_u_2 * V_u_2 &  
!               + Gm_dd_33 * V_u_3 * V_u_3) )
!
!
!
!    B_d_1 = Gm_dd_11 * B_u_1
!    B_d_2 = Gm_dd_22 * B_u_2
!    B_d_3 = Gm_dd_33 * B_u_3
!
!    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )
!
!    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
!          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
!    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
!          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
!    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
!          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )
!
!    vI = V_u_1 * I_d_1 + V_u_2 * I_d_2 + V_u_3 * I_d_3
!
!
!   Flux_G(1) = I_u_1 + W * D * V_u_1 - ( B_u_1 / alp ) * vI
!   Flux_G(2) = I_u_2 + W * D * V_u_2 - ( B_u_2 / alp ) * vI
!   Flux_G(3) = I_u_3 + W * D * V_u_3 - ( B_u_3 / alp ) * vI
!
!   END FUNCTION Flux_G

END MODULE TwoMoment_UtilitiesModule_Relativistic!
