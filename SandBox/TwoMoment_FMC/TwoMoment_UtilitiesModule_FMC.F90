MODULE TwoMoment_UtilitiesModule_FMC

  USE KindModule, ONLY: &
    DP, Half, Zero, One, Three, Five
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_V1, iPF_V2, iPF_V3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor_Relativistic, EddingtonFactor, HeatFluxFactor
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, &
    nPM, iPM_J, iPM_H1, iPM_H2, iPM_H3, &
    nAM, &
    nGM

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeConserved_TwoMoment_FMC
  PUBLIC :: ComputeFromConserved_TwoMoment_FMC
  PUBLIC :: ComputePrimitive_TwoMoment_Richardson_FMC
  PUBLIC :: ComputeHatMomentsFromConserved
  PUBLIC :: ComputeHatMomentsFromPrimitive
  PUBLIC :: EddingtonTensorComponents_dd
  PUBLIC :: HeatFluxTensorComponents_uuu
  PUBLIC :: ComputeHeatFluxTensorComponents_ddd_Lagrangian
  PUBLIC :: ComputeHeatFluxTensorComponents_uud_Lagrangian
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: Flux_E
  PUBLIC :: NumericalFlux_LLF

  CONTAINS

  SUBROUTINE ComputeConserved_TwoMoment_FMC &
    ( J, H_d_1, H_d_2, H_d_3, E, F_d_1, F_d_2, F_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in)  :: J, H_d_1, H_d_2, H_d_3 ! --- Index down
    REAL(DP), INTENT(out) :: E, F_d_1, F_d_2, F_d_3 ! --- Index down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index up 
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: k_dd(3,3), vMag, W, vFhat
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3 ! --- Index down

    k_dd = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33)

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )

    CALL ComputeHatMomentsFromPrimitive &
      ( J, H_d_1, H_d_2, H_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Eulerian Energy Density ---
        E = W*E_hat+vFhat

    ! --- Eulerian Momentum Density (1) ---
        F_d_1 = F_hat_d_1 + W*V_u_1*E_hat

    ! --- Eulerian Momentum Density (2) ---
        F_d_2 = F_hat_d_2 + W*V_u_2*E_hat

    ! --- Eulerian Momentum Density (3) ---
        F_d_3 = F_hat_d_3 + W*V_u_3*E_hat

  END SUBROUTINE ComputeConserved_TwoMoment_FMC

  SUBROUTINE ComputeFromConserved_TwoMoment_FMC &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CM, PM, AM, GM )
    
    ! --- Input/Output variables ---
    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nPF)
    REAL(DP), INTENT(in)  :: &
      CM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      AM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      GM(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGM,1:nSpecies)

    ! --- Local variables ---
    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        CALL ComputePrimitive_TwoMoment_Richardson_FMC &
               ( CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_TwoMoment_FMC

  SUBROUTINE ComputePrimitive_TwoMoment_Richardson_FMC &
    ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---
    INTEGER, PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: TOL = 1.0d-08

    ! --- Local variables ---
    INTEGER :: iteration, i, k
    REAL(DP) :: vMag, vH, vK, lambda, W
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP) :: Mvec(4), Uvec(4), Fvec(4), Gvec(4), dX(4), k_dd(3,3)

    LOGICAL :: CONVERGED

    ! --- Initial guess ---
    CALL ComputeHatMomentsFromConserved &
      ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    Uvec = [ E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3 ]
    Mvec(:) = Uvec(:)

    ! --- Richardson update ---
    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
                 + V_u_2 * Gm_dd_22 * V_u_2 &
                 + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    lambda = One / ( W + W * vMag )  

    iteration = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. iteration < MaxIterations )

      iteration = iteration + 1
      
      k_dd = EddingtonTensorComponents_dd &
        ( Mvec(1), Mvec(2), Mvec(3), Mvec(4), V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )
      
      vH = V_u_1 * Mvec(2) + V_u_2 * Mvec(3) + V_u_3 * Mvec(4)
      Fvec(1) = Uvec(1) - vH
      Gvec(1) = (One - lambda * W) * Mvec(1) + lambda * Fvec(1)


      DO i = 1, 3

        vK = V_u_1 * k_dd(i,1) + V_u_2 * k_dd(i,2) + V_u_3 * k_dd(i,3)
        Fvec(i+1) = Uvec(i+1) - vK * Mvec(1)
        Gvec(i+1) = (One - lambda * W) * Mvec(i+1) + lambda * Fvec(i+1)

      END DO

      !dX = ABS ( Gvec - Mvec )
      CONVERGED = SQRT( SUM( (Fvec-W*Mvec)**2 ) ) < TOL

      ! print *, "k: ", iteration
      ! print *, "||F(M)||: ", SQRT( SUM( (Fvec-W*Mvec)**2 ) )
      Mvec(:) = Gvec(:)

    END DO

    ! --- Converged solution ---
    J = Mvec(1)
    H_d_1 = Mvec(2)
    H_d_2 = Mvec(3)
    H_d_3 = Mvec(4)

    !print *, iteration

    IF( PRESENT( nIterations_Option ) ) THEN

      nIterations_Option = iteration

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_Richardson_FMC

  SUBROUTINE ComputeHatMomentsFromConserved &
    ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! -- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: W, vMag, vF
    
    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
                 + V_u_2 * Gm_dd_22 * V_u_2 &
                 + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    vF = V_u_1 * F_d_1 + V_u_2 * F_d_2 + V_u_3 * F_d_3

    E_hat = W * ( E - vF )
    F_hat_d_1 = F_d_1 - W * V_u_1 * E_hat
    F_hat_d_2 = F_d_2 - W * V_u_2 * E_hat
    F_hat_d_3 = F_d_3 - W * V_u_3 * E_hat

  END SUBROUTINE ComputeHatMomentsFromConserved

  SUBROUTINE ComputeHatMomentsFromPrimitive &
    ( J, H_d_1, H_d_2, H_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
      V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33)
    
    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(out) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: k_dd(3,3)
    REAL(DP) :: vMag, W, vH, vK

    k_dd = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
                 + V_u_2 * Gm_dd_22 * V_u_2 &
                 + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    vH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3

    E_hat = W*J + vH
    F_hat_d_1 = W*H_d_1 + V_u_1*k_dd(1,1)*J + V_u_2*k_dd(1,2)*J + V_u_3*k_dd(1,3)*J
    F_hat_d_2 = W*H_d_2 + V_u_1*k_dd(2,1)*J + V_u_2*k_dd(2,2)*J + V_u_3*k_dd(2,3)*J
    F_hat_d_3 = W*H_d_3 + V_u_1*k_dd(3,1)*J + V_u_2*k_dd(3,2)*J + V_u_3*k_dd(3,3)*J

  END SUBROUTINE ComputeHatMomentsFromPrimitive

  FUNCTION EddingtonTensorComponents_dd &
    ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: EddingtonTensorComponents_dd(3,3)

    ! --- Local variables ---
    INTEGER :: i, k
    REAL(DP) :: FF, EF, a, b, W, vMag
    REAL(DP) :: h_hat_d(3), u(3), Gm_dd(3,3) 

    FF = FluxFactor_Relativistic &
      ( J, H_d_1, H_d_2, H_d_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33, &
        -One, Zero, Zero, Zero, &
        V_u_1, V_u_2, V_u_3 )

    EF = EddingtonFactor( J, FF )

    a = Half * ( One - EF )
    b= Half * ( Three * EF - One )

    h_hat_d(1) = H_d_1 / ( FF * J )
    h_hat_d(2) = H_d_2 / ( FF * J )
    h_hat_d(3) = H_d_3 / ( FF * J )

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )

    u(1) = W * V_u_1
    u(2) = W * V_u_2
    u(3) = W * V_u_3

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    DO k = 1, 3
      DO i = 1, 3

        EddingtonTensorComponents_dd(i,k) &
          = a * ( Gm_dd(i,k) + u(i)*u(k) ) + b * h_hat_d(i) * h_hat_d(k)

      END DO
    END DO

    RETURN
  END FUNCTION EddingtonTensorComponents_dd

  SUBROUTINE HeatFluxTensorComponents_uuu &
    ( J, H_u_1, H_u_2, H_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      V_u_1, V_u_2, V_u_3, l_uuu_munurho )

    ! Loops to shorten code and get rid of symmetric variables?

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: l_uuu_munurho(0:3,0:3,0:3)

    ! --- Local variables ---
    REAL(DP) :: FF, HF
    REAL(DP) :: W, vMag, H_u_0
    REAL(DP) :: h_hat_u_1, h_hat_u_2, h_hat_u_3
    REAL(DP) :: h_uu_11, h_uu_12, h_uu_13, h_uu_22, h_uu_23, h_uu_33
    REAL(DP) :: a, b
    REAL(DP) :: l_uuu_111, l_uuu_112, l_uuu_113, &
                l_uuu_121, l_uuu_122, l_uuu_123, &
                l_uuu_131, l_uuu_132, l_uuu_133, &
                l_uuu_211, l_uuu_212, l_uuu_213, &
                l_uuu_221, l_uuu_222, l_uuu_223, &
                l_uuu_231, l_uuu_232, l_uuu_233, &
                l_uuu_311, l_uuu_312, l_uuu_313, &
                l_uuu_321, l_uuu_322, l_uuu_323, &
                l_uuu_331, l_uuu_332, l_uuu_333
    REAL(DP) :: l_uuu_011, l_uuu_012, l_uuu_013, &
                l_uuu_021, l_uuu_022, l_uuu_023, &
                l_uuu_031, l_uuu_032, l_uuu_033, &
                l_uuu_101, l_uuu_102, l_uuu_103, &
                l_uuu_201, l_uuu_202, l_uuu_203, &
                l_uuu_301, l_uuu_302, l_uuu_303, &
                l_uuu_110, l_uuu_120, l_uuu_130, &
                l_uuu_210, l_uuu_220, l_uuu_230, &
                l_uuu_310, l_uuu_320, l_uuu_330, &
                l_uuu_001, l_uuu_002, l_uuu_003, &
                l_uuu_100, l_uuu_200, l_uuu_300, &
                l_uuu_010, l_uuu_020, l_uuu_030, &
                l_uuu_000


    FF = FluxFactor_Relativistic &
      ( J, H_u_1, H_u_2, H_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33, &
          -One, Zero, Zero, Zero, &
          V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor(J, FF)

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )

    h_hat_u_1 = H_u_1 / ( FF * J )
    h_hat_u_2 = H_u_2 / ( FF * J )
    h_hat_u_3 = H_u_3 / ( FF * J )

    h_uu_11 = Gm_dd_11 + W**2 * V_u_1 * V_u_1
    h_uu_12 = W**2 * V_u_1 * V_u_2
    h_uu_13 = W**2 * V_u_1 * V_u_3
    h_uu_22 = Gm_dd_22 + W**2 * V_u_2 * V_u_2
    h_uu_23 = W**2 * V_u_2 * V_u_3
    h_uu_33 = Gm_dd_33 + W**2 * V_u_3 * V_u_3

    a = Half * (FF - HF)
    b = Half * (Five * HF - Three * FF)

    ! --- ijk components ---

    l_uuu_111 = a * (h_hat_u_1 * h_uu_11 + h_hat_u_1 * h_uu_11 + h_hat_u_1 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_1)
    l_uuu_112 = a * (h_hat_u_1 * h_uu_12 + h_hat_u_1 * h_uu_12 + h_hat_u_2 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_2)
    l_uuu_113 = a * (h_hat_u_1 * h_uu_13 + h_hat_u_1 * h_uu_13 + h_hat_u_3 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_3)
    l_uuu_121 = l_uuu_112
    l_uuu_122 = a * (h_hat_u_1 * h_uu_22 + h_hat_u_2 * h_uu_12 + h_hat_u_2 * h_uu_12) &
                + b * (h_hat_u_1 * h_hat_u_2 * h_hat_u_2)
    l_uuu_123 = a * (h_hat_u_1 * h_uu_23 + h_hat_u_2 * h_uu_13 + h_hat_u_3 * h_uu_12) &
                + b * (h_hat_u_1 * h_hat_u_2 * h_hat_u_3)
    l_uuu_131 = l_uuu_113
    l_uuu_132 = l_uuu_123
    l_uuu_133 = a * (h_hat_u_1 * h_uu_33 + h_hat_u_3 * h_uu_13 + h_hat_u_3 * h_uu_13) &
                + b * (h_hat_u_1 * h_hat_u_3 * h_hat_u_3)
    l_uuu_211 = l_uuu_112
    l_uuu_212 = l_uuu_122
    l_uuu_213 = l_uuu_123
    l_uuu_221 = l_uuu_122
    l_uuu_222 = a * (h_hat_u_2 * h_uu_22 + h_hat_u_2 * h_uu_22 + h_hat_u_2 * h_uu_22) &
                + b * (h_hat_u_2 * h_hat_u_2 * h_hat_u_2)
    l_uuu_223 = a * (h_hat_u_2 * h_uu_23 + h_hat_u_2 * h_uu_23 + h_hat_u_3 * h_uu_22) &
                + b * (h_hat_u_2 * h_hat_u_2 * h_hat_u_3)
    l_uuu_231 = l_uuu_123
    l_uuu_232 = l_uuu_223
    l_uuu_233 = a * (h_hat_u_2 * h_uu_33 + h_hat_u_3 * h_uu_23 + h_hat_u_3 * h_uu_23) &
                + b * (h_hat_u_2 * h_hat_u_3 * h_hat_u_3)
    l_uuu_311 = l_uuu_113
    l_uuu_312 = l_uuu_123
    l_uuu_313 = l_uuu_133
    l_uuu_321 = l_uuu_123
    l_uuu_322 = l_uuu_223
    l_uuu_323 = l_uuu_233
    l_uuu_331 = l_uuu_133
    l_uuu_332 = l_uuu_233
    l_uuu_333 = a * (h_hat_u_3 * h_uu_33 + h_hat_u_3 * h_uu_33 + h_hat_u_3 * h_uu_33) &
                + b * (h_hat_u_3 * h_hat_u_3 * h_hat_u_3)

    ! --- 0th components ---

    l_uuu_011 = V_u_1 * l_uuu_111 + V_u_2 * l_uuu_211 + V_u_3 * l_uuu_311
    l_uuu_012 = V_u_1 * l_uuu_112 + V_u_2 * l_uuu_212 + V_u_3 * l_uuu_312
    l_uuu_013 = V_u_1 * l_uuu_113 + V_u_2 * l_uuu_213 + V_u_3 * l_uuu_313
    l_uuu_021 = l_uuu_012
    l_uuu_022 = V_u_1 * l_uuu_122 + V_u_2 * l_uuu_222 + V_u_3 * l_uuu_322
    l_uuu_023 = V_u_1 * l_uuu_123 + V_u_2 * l_uuu_223 + V_u_3 * l_uuu_323
    l_uuu_031 = l_uuu_013
    l_uuu_032 = l_uuu_023
    l_uuu_033 = V_u_1 * l_uuu_133 + V_u_2 * l_uuu_233 + V_u_3 * l_uuu_333
    l_uuu_101 = l_uuu_011
    l_uuu_102 = l_uuu_012
    l_uuu_103 = l_uuu_013
    l_uuu_201 = l_uuu_012
    l_uuu_202 = l_uuu_022
    l_uuu_203 = l_uuu_023
    l_uuu_301 = l_uuu_013
    l_uuu_302 = l_uuu_023
    l_uuu_303 = l_uuu_033
    l_uuu_110 = l_uuu_011
    l_uuu_120 = l_uuu_012
    l_uuu_130 = l_uuu_013
    l_uuu_210 = l_uuu_012
    l_uuu_220 = l_uuu_022
    l_uuu_230 = l_uuu_023
    l_uuu_310 = l_uuu_013
    l_uuu_320 = l_uuu_023
    l_uuu_330 = l_uuu_033
    l_uuu_001 = V_u_1 * l_uuu_011 + V_u_2 * l_uuu_021 + V_u_3 * l_uuu_031
    l_uuu_002 = V_u_1 * l_uuu_012 + V_u_2 * l_uuu_022 + V_u_3 * l_uuu_032
    l_uuu_003 = V_u_1 * l_uuu_013 + V_u_2 * l_uuu_023 + V_u_3 * l_uuu_033
    l_uuu_100 = l_uuu_001
    l_uuu_200 = l_uuu_002
    l_uuu_300 = l_uuu_003
    l_uuu_010 = l_uuu_001
    l_uuu_020 = l_uuu_002
    l_uuu_030 = l_uuu_003
    l_uuu_000 = V_u_1 * l_uuu_001 + V_u_2 * l_uuu_002 + V_u_3 * l_uuu_003

    ! --- Array assignment ---

    l_uuu_munurho(0,0,0) = l_uuu_000
    l_uuu_munurho(0,0,1) = l_uuu_001
    l_uuu_munurho(0,0,2) = l_uuu_002
    l_uuu_munurho(0,0,3) = l_uuu_003
    l_uuu_munurho(0,1,0) = l_uuu_010
    l_uuu_munurho(0,1,1) = l_uuu_011
    l_uuu_munurho(0,1,2) = l_uuu_012
    l_uuu_munurho(0,1,3) = l_uuu_013
    l_uuu_munurho(0,2,0) = l_uuu_020
    l_uuu_munurho(0,2,1) = l_uuu_021
    l_uuu_munurho(0,2,2) = l_uuu_022
    l_uuu_munurho(0,2,3) = l_uuu_023
    l_uuu_munurho(0,3,0) = l_uuu_030
    l_uuu_munurho(0,3,1) = l_uuu_031
    l_uuu_munurho(0,3,2) = l_uuu_032
    l_uuu_munurho(0,3,3) = l_uuu_033
    l_uuu_munurho(1,0,0) = l_uuu_100
    l_uuu_munurho(1,0,1) = l_uuu_101
    l_uuu_munurho(1,0,2) = l_uuu_102
    l_uuu_munurho(1,0,3) = l_uuu_103
    l_uuu_munurho(1,1,0) = l_uuu_110
    l_uuu_munurho(1,1,1) = l_uuu_111
    l_uuu_munurho(1,1,2) = l_uuu_112
    l_uuu_munurho(1,1,3) = l_uuu_113
    l_uuu_munurho(1,2,0) = l_uuu_120
    l_uuu_munurho(1,2,1) = l_uuu_121
    l_uuu_munurho(1,2,2) = l_uuu_122
    l_uuu_munurho(1,2,3) = l_uuu_123
    l_uuu_munurho(1,3,0) = l_uuu_130
    l_uuu_munurho(1,3,1) = l_uuu_131
    l_uuu_munurho(1,3,2) = l_uuu_132
    l_uuu_munurho(1,3,3) = l_uuu_133
    l_uuu_munurho(2,0,0) = l_uuu_200
    l_uuu_munurho(2,0,1) = l_uuu_201
    l_uuu_munurho(2,0,2) = l_uuu_202
    l_uuu_munurho(2,0,3) = l_uuu_203
    l_uuu_munurho(2,1,0) = l_uuu_210
    l_uuu_munurho(2,1,1) = l_uuu_211
    l_uuu_munurho(2,1,2) = l_uuu_212
    l_uuu_munurho(2,1,3) = l_uuu_213
    l_uuu_munurho(2,2,0) = l_uuu_220
    l_uuu_munurho(2,2,1) = l_uuu_221
    l_uuu_munurho(2,2,2) = l_uuu_222
    l_uuu_munurho(2,2,3) = l_uuu_223
    l_uuu_munurho(2,3,0) = l_uuu_230
    l_uuu_munurho(2,3,1) = l_uuu_231
    l_uuu_munurho(2,3,2) = l_uuu_232
    l_uuu_munurho(2,3,3) = l_uuu_233
    l_uuu_munurho(3,0,0) = l_uuu_300
    l_uuu_munurho(3,0,1) = l_uuu_301
    l_uuu_munurho(3,0,2) = l_uuu_302
    l_uuu_munurho(3,0,3) = l_uuu_303
    l_uuu_munurho(3,1,0) = l_uuu_310
    l_uuu_munurho(3,1,1) = l_uuu_311
    l_uuu_munurho(3,1,2) = l_uuu_312
    l_uuu_munurho(3,1,3) = l_uuu_313
    l_uuu_munurho(3,2,0) = l_uuu_320
    l_uuu_munurho(3,2,1) = l_uuu_321
    l_uuu_munurho(3,2,2) = l_uuu_322
    l_uuu_munurho(3,2,3) = l_uuu_323
    l_uuu_munurho(3,3,0) = l_uuu_330
    l_uuu_munurho(3,3,1) = l_uuu_331
    l_uuu_munurho(3,3,2) = l_uuu_332
    l_uuu_munurho(3,3,3) = l_uuu_333
  
  END SUBROUTINE HeatFluxTensorComponents_uuu

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

    Gm_uu_11 = 1.0_DP /  Gm_dd_11 - B_u_1**2 / alp**2
    Gm_uu_22 = 1.0_DP /  Gm_dd_22 - B_u_2**2 / alp**2
    Gm_uu_33 = 1.0_DP /  Gm_dd_33 - B_u_3**2 / alp**2

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

  FUNCTION Flux_X1 &
    ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X1(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMag, W, vH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_1

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    vH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_1 = V_u_1 * k_ud(1,1) + V_u_2 * k_ud(1,2) + V_u_3 * k_ud(1,3)

    Flux_X1(1) = W * H_u_1 + W * V_u_1 * (W * J + vH) + vK_u_1 * J
    Flux_X1(2) = k_ud(1,1) + W * (H_u_1 * V_u_1 + V_u_1 * H_u_1) + W * W * V_u_1 * V_u_1 * J
    Flux_X1(3) = k_ud(1,2) + W * (H_u_1 * V_u_2 + V_u_1 * H_u_2) + W * W * V_u_1 * V_u_2 * J
    Flux_X1(4) = k_ud(1,3) + W * (H_u_1 * V_u_3 + V_u_1 * H_u_3) + W * W * V_u_1 * V_u_3 * J
    RETURN

  END FUNCTION Flux_X1

  FUNCTION Flux_X2 ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X2(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMag, W, vH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_2

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    vH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_2 = V_u_1 * k_ud(2,1) + V_u_2 * k_ud(2,2) + V_u_3 * k_ud(2,3)

    Flux_X2(1) = W * H_u_2 + W * V_u_2 * (W * J + vH) + vK_u_2 * J
    Flux_X2(2) = k_ud(2,1) + W * (H_u_2 * V_u_1 + V_u_2 * H_u_1) + W * W * V_u_2 * V_u_1 * J
    Flux_X2(3) = k_ud(2,2) + W * (H_u_2 * V_u_2 + V_u_2 * H_u_2) + W * W * V_u_2 * V_u_2 * J
    Flux_X2(4) = k_ud(2,3) + W * (H_u_2 * V_u_3 + V_u_2 * H_u_3) + W * W * V_u_2 * V_u_3 * J
    RETURN

  END FUNCTION Flux_X2

  FUNCTION Flux_X3 ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X3(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMag, W, vH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_3

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    vH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_3 = V_u_1 * k_ud(3,1) + V_u_2 * k_ud(3,2) + V_u_3 * k_ud(3,3)

    Flux_X3(1) = W * H_u_3 + W * V_u_3 * (W * J + vH) + vK_u_3 * J
    Flux_X3(2) = k_ud(3,1) + W * (H_u_3 * V_u_1 + V_u_3 * H_u_1) + W * W * V_u_3 * V_u_1 * J
    Flux_X3(3) = k_ud(3,2) + W * (H_u_3 * V_u_2 + V_u_3 * H_u_2) + W * W * V_u_3 * V_u_2 * J
    Flux_X3(4) = k_ud(3,3) + W * (H_u_3 * V_u_3 + V_u_3 * H_u_3) + W * W * V_u_3 * V_u_3 * J
    RETURN

  END FUNCTION Flux_X3

  FUNCTION Flux_E( J, H_d_1, H_d_2, H_d_3, &
    V_u_1, V_u_2, V_u_3, Jacobian_U, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_E(4)
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3, Jacobian_U(0:3,0:3)
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMag, W, u_u(0:3)
    REAL(DP) :: vH, H_u(0:3)
    REAL(DP) :: k_uu(0:3,0:3), l_uuu(0:3,0:3,0:3), Q_uuu(0:3,0:3,0:3)
    INTEGER :: mu, nu, rho

    vMag = SQRT( V_u_1 * Gm_dd_11 * V_u_1 &
      + V_u_2 * Gm_dd_22 * V_u_2 &
      + V_u_3 * Gm_dd_33 * V_u_3 )
    W = One / SQRT ( One - vMag**2 )
    u_u(0) = W
    u_u(1) = W * V_u_1
    u_u(2) = W * V_u_2
    u_u(3) = W * V_u_3

    vH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3
    H_u(0) = vH
    H_u(1) = H_d_1
    H_u(2) = H_d_2
    H_u(3) = H_d_3

    k_uu(1:3,1:3) = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )
    DO nu = 1,3

      k_uu(0,nu) = V_u_1 * k_uu(nu,1) + V_u_2 * k_uu(nu,2) + V_u_3 * k_uu(nu,3)
      k_uu(nu,0) = k_uu(0,nu)

    END DO
    k_uu(0,0) = V_u_1 * k_uu(0,1) + V_u_2 * k_uu(0,2) + V_u_3 * k_uu(0,3)

    CALL HeatFluxTensorComponents_uuu &
      ( J, H_d_1, H_d_2, H_d_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
        V_u_1, V_u_2, V_u_3, l_uuu )

    Flux_E = 0.0_DP
    DO rho = 0,3
      DO nu = 0,3
        DO mu = 0,3

          Q_uuu(mu,nu,rho) = J*u_u(mu)*u_u(nu)*u_u(rho)+ &
                             H_u(mu)*u_u(nu)*u_U(rho)+ &
                             H_u(nu)*u_u(mu)*u_u(rho)+ &
                             H_u(rho)*u_u(mu)*u_u(nu)+ &
                             k_uu(mu,nu)*J*u_u(rho)+ &
                             k_uu(mu,rho)*J*u_u(nu)+ &
                             k_uu(nu,rho)*J*u_u(mu)+ &
                             l_uuu(mu,nu,rho)

        END DO
        Flux_E(1) = Flux_E(1) + Q_uuu(0,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(2) = Flux_E(2) + Q_uuu(1,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(3) = Flux_E(3) + Q_uuu(2,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(4) = Flux_E(4) + Q_uuu(3,nu,rho) * Jacobian_U(nu,rho)
      END DO
    END DO
    Flux_E = -Flux_E
    RETURN

  END FUNCTION Flux_E

  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


END MODULE TwoMoment_UtilitiesModule_FMC
