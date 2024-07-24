!> Perform computations related to the 3+1, CFA Euler equations.
!> Find the equations in Rezzolla & Zanotti, Relativistic Hydrodynamics, 2013,
!> Equation 7.234.
MODULE MHD_UtilitiesModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Half, &
    One, &
    Two, &
    Four, &
    Fourth
  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX, &
    nDimsX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    dLXdX3_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up 
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nPM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nAM, &
    iAM_P, &
    iAM_T, &
    iAM_Ye, &
    iAM_S, &
    iAM_E, &
    iAM_Gm, &
    iAM_Cs, &
    nDM, &
    iDM_Sh_X1, &
    iDM_Sh_X2, &
    iDM_Sh_X3, &
    iDM_Div
  USE MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputeAuxiliary_Fluid, &
    ComputePressure, &
    ComputePressureFromPrimitive, &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSpecificInternalEnergy
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    Min_Y, &
    Max_Y
  USE UnitsModule, ONLY: &
    AtomicMassUnit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_MHD_Relativistic
  PUBLIC :: ComputeConserved_MHD_Relativistic
  PUBLIC :: ComputeMagneticDivergence_MHD_Relativistic
  PUBLIC :: ComputeWeakMagneticDivergence_MHD_Relativistic
  PUBLIC :: ComputeFromConserved_MHD_Relativistic
  PUBLIC :: ComputeTimeStep_MHD_Relativistic
  PUBLIC :: Eigenvalues_MHD_Relativistic
  PUBLIC :: Flux_X1_MHD_Relativistic
  PUBLIC :: Flux_X2_MHD_Relativistic
  PUBLIC :: Flux_X3_MHD_Relativistic
  PUBLIC :: NumericalFlux_HLL_MHD_Relativistic

  REAL(DP), POINTER, CONTIGUOUS :: &
    Gm_dd_11_K(:), Gm_dd_22_K(:), Gm_dd_33_K(:), SqrtGm_K(:), &
    Gm_dd_11_F(:), Gm_dd_22_F(:), Gm_dd_33_F(:), SqrtGm_F(:), &
    Beta_1_K  (:), Beta_2_K  (:), Beta_3_K  (:), Alpha_K (:), &
    Beta_1_F  (:), Beta_2_F  (:), Beta_3_F  (:), Alpha_F (:)

  REAL(DP), POINTER, CONTIGUOUS :: &
    uD_K(:), uS1_K(:), uS2_K(:), uS3_K(:), uE_K(:), uNe_K(:), &
    uB1_K(:), uB2_K(:), uB3_K(:), uChi_K(:), &
    uD_L(:), uS1_L(:), uS2_L(:), uS3_L(:), uE_L(:), uNe_L(:), &
    uB1_L(:), uB2_L(:), uB3_L(:), uChi_L(:), &
    uD_R(:), uS1_R(:), uS2_R(:), uS3_R(:), uE_R(:), uNe_R(:), &
    uB1_R(:), uB2_R(:), uB3_R(:), uChi_R(:)
 
  REAL(DP), ALLOCATABLE :: &
    pD_K(:), pV1_K(:), pV2_K(:), pV3_K(:), pE_K(:), pNe_K(:), &
    pB1_K(:), pB2_K(:), pB3_K(:), pChi_K(:), &
    pD_L(:), pV1_L(:), pV2_L(:), pV3_L(:), pE_L(:), pNe_L(:), &
    pB1_L(:), pB2_L(:), pB3_L(:), pChi_L(:), & 
    pD_R(:), pV1_R(:), pV2_R(:), pV3_R(:), pE_R(:), pNe_R(:), &
    pB1_R(:), pB2_R(:), pB3_R(:), pChi_R(:)

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableX_F(:,:), IndexTableX_V(:,:)

  INTEGER :: nX   (3), nX_K,  nNodesX_K
  INTEGER :: nX_X1(3), nX1_X, nNodesX_X1
  INTEGER :: nX_X2(3), nX2_X, nNodesX_X2
  INTEGER :: nX_X3(3), nX3_X, nNodesX_X3

  INTERFACE ComputePrimitive_MHD_Relativistic
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector
  END INTERFACE ComputePrimitive_MHD_Relativistic

  INTERFACE ComputeConserved_MHD_Relativistic
    MODULE PROCEDURE ComputeConserved_Scalar
    MODULE PROCEDURE ComputeConserved_Vector
  END INTERFACE ComputeConserved_MHD_Relativistic


CONTAINS


  !> Compute the primitive variables from the conserved variables,
  !> a la Kastaun, Kalinani, and Ciolfi (2021), Phys. Rev. D., 103, 023018
  !> @todo Incorporate magnetic fields and move to separate modules.
  SUBROUTINE ComputePrimitive_Scalar &
    ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi, &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi, &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3, &
      EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    REAL(DP), INTENT(in)    :: &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi
    REAL(DP), INTENT(in)    :: &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3
    REAL(DP), INTENT(out)   :: &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi

    REAL(DP) :: B0u, B1, B2, B3
    REAL(DP) :: D, D_bar, tau, Ne, Ne_bar
    REAL(DP) :: Sd1, Sd2, Sd3
    REAL(DP) :: q, q_bar
    REAL(DP) :: rd1, rd2, rd3, ru1, ru2, ru3, r
    REAL(DP) :: bd1, bd2, bd3, bu1, bu2, bu3, b
    REAL(DP) :: bSq, r_barSq, rb
    REAL(DP) :: h0, r0, z0, v0
    REAL(DP) :: mu
    REAL(DP) :: MinEps, MinP
    REAL(DP) :: x, W, eps, Vsq

    REAL(DP) :: AM_P, Cs

    IF( .NOT. EvolveOnlyMagnetic )THEN

     !PRINT*, 'Magnetofluid coupling on for conserved-to-primitive.'

      B1 = CM_B1
      B2 = CM_B2
      B3 = CM_B3

    ELSE

     !PRINT*, 'Magnetofluid coupling off for conserved-to-primitive.'

      B1 = 0.0_DP
      B2 = 0.0_DP
      B3 = 0.0_DP

    END IF

    ! --- Eqs. 18-20                                                    ---
    ! --- Note: Without magnetic fields, barred and unbarred quantities ---
    ! --- are equivalent.                                               ---

    D     = CM_D
    D_bar = D

   !PRINT*
   !PRINT*, 'Computing primitive variables.'
   !PRINT*, 'CM_D: ',  CM_D
   !PRINT*, 'CM_S1: ', CM_S1
   !PRINT*, 'CM_S2: ', CM_S2
   !PRINT*, 'CM_S3: ', CM_S3
   !PRINT*, 'CM_B1: ', CM_B1
   !PRINT*, 'CM_B2: ', CM_B2
   !PRINT*, 'CM_B3: ', CM_B3
   !PRINT*, 'CM_E: ',  CM_E
   !PRINT*, 'CM_Ne: ', CM_Ne

   !PRINT*, 'D_bar: ', D

    Ne     = CM_Ne
    Ne_bar = Ne

   !PRINT*, 'Ne_bar: ', Ne_bar

    tau = CM_E

   !*, 'tau: ', tau

    Sd1 = CM_S1
    Sd2 = CM_S2
    Sd3 = CM_S3

   !PRINT*, 'Sd1: ', Sd1
   !PRINT*, 'Sd2: ', Sd2
   !PRINT*, 'Sd3: ', Sd3

    ! --- Eqs. 22-24 ---

    q = tau / D_bar

   !PRINT*, 'q: ', q

    rd1 = Sd1 / D_bar
    rd2 = Sd2 / D_bar
    rd3 = Sd3 / D_bar

   !PRINT*, 'rd1: ', rd1
   !PRINT*, 'rd2: ', rd2
   !PRINT*, 'rd3: ', rd3

    ru1 = rd1 / GF_Gm11 
    ru2 = rd2 / GF_Gm22
    ru3 = rd3 / GF_Gm33

   !PRINT*, 'ru1: ', ru1
   !PRINT*, 'ru2: ', ru2
   !PRINT*, 'ru3: ', ru3

    r = SQRT( ru1 * rd1 + ru2 * rd2 + ru3 * rd3 )

   !PRINT*, 'r: ', r

    bu1 = B1 / SQRT( D_bar )
    bu2 = B2 / SQRT( D_bar )
    bu3 = B3 / SQRT( D_bar )

   !PRINT*, 'bu1: ', bu1
   !PRINT*, 'bu2: ', bu2
   !PRINT*, 'bu3: ', bu3

    bd1 = GF_Gm11 * bu1
    bd2 = GF_Gm22 * bu2
    bd3 = GF_Gm33 * bu3

   !PRINT*, 'bd1: ', bd1
   !PRINT*, 'bd2: ', bd2
   !PRINT*, 'bd3: ', bd3

    ! --- Squared magnitude of b ---

    bSq = ( bu1 * bd1 + bu2 * bd2 + bu3 * bd3 )

   !PRINT*, 'bSq: ', bSq

    ! --- Contraction of r with b ---

    rb = ( bu1 * rd1 + bu2 * rd2 + bu3 * rd3 )

   !PRINT*, 'rb: ', rb

#ifdef MICROPHYSICS_WEAKLIB

    !*, 'MinD: ', MinD
    !*, 'MinT: ', MinT
    !*, 'MinY: ', MinY

    CALL ComputePressure &
           ( Min_D, Min_T, Min_Y, Min_P )

    !*, 'MinP: ', MinP

    CALL ComputeSpecificInternalEnergy &
           ( Min_D, Min_T, Min_Y, MinEps )

    !*, 'MinEps: ', MinEps
 

    ! --- Not valid. Need better method of getting         ---
    ! --- minimum enthalpy (and deal with ideal gas case). ---

    h0 = One + MinEps + Min_P / Min_D  

    !*, 'h0: ', /0

#else

    h0 = One

#endif

    ! --- Upper velocity limit (Eqs. 32 and 33) --- 

    z0 = r / h0

   !PRINT*, 'z0: ', z0

    v0 = z0 / SQRT( One + z0**2 )

   !PRINT*, 'v0: ', v0

    CALL SolveMu_Bisection &
           ( D_bar, Ne_bar, q, r, &
             bSq, rb, h0, v0, mu )
 
   !PRINT*, 'mu: ', mu

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

   !PRINT*, 'x: ', x

    ! --- Eqs. 38 and 39                              ---
    ! --- Note: Eq. 39 rewritten using Eqs. 25 and 30 ---
    ! --- to avoid B = 0 degeneracy.                  ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

    q_bar = q - Half * bSq - Half * mu**2 * x**2 * ( bSq * r**2 - rb**2 )

   !PRINT*, 'q_bar: ', q_bar

    ! --- Eq. 68 ---

    PM_V1 = mu * x * ( ru1 + mu * rb * bu1 )
    PM_V2 = mu * x * ( ru2 + mu * rb * bu2 )
    PM_V3 = mu * x * ( ru3 + mu * rb * bu3 )

   !PRINT*, 'PM_V1: ', PM_V1
   !PRINT*, 'PM_V2: ', PM_V2
   !PRINT*, 'PM_V3: ', PM_V3 

    VSq = PM_V1**2 * GF_Gm11 + PM_V2**2 * GF_Gm22 + PM_V3**2 * GF_Gm33

    ! --- Eq. 40 ---

    W = One / SQRT( One - VSq )

   !PRINT*, 'W: ', W

    ! --- Eq. 41 ---

    PM_D = D_bar / W

   !PRINT*, 'PM_D: ', PM_D

    PM_Ne = CM_Ne / W

    ! --- Eq. 42 ---

    ! See line 203 of RePrimAnd's con2prim_imhd.cc (these should be equivalent).

    eps = W * (q_bar - mu * r_barSq * ( One - mu * W / ( One + W ) ) )

    !PRINT*, 'eps: ', eps

    PM_E = PM_D * eps

   !PRINT*, 'PM_E: ', PM_E

    CALL ComputePressureFromPrimitive &
           ( PM_D, PM_E, PM_Ne, AM_P )
 
   !PRINT*, 'mu: ', One / ( W * ( One + eps + AM_P / PM_D ) )

    B0u = ( W / GF_alpha ) &
            * ( GF_Gm11 * PM_V1 * CM_B1 &
                + GF_Gm22 * PM_V2 * CM_B2 &
                + GF_Gm33 * PM_V3 * CM_B3 )
                      
   !PRINT* , 'B0u: ', B0u
                            
    PM_B1 = ( CM_B1 / W ) + GF_Alpha * B0u &
                            * ( PM_V1 - ( GF_Beta1 / GF_Alpha ) ) 
   !PRINT*, 'PM_B1: ', PM_B1
 
    PM_B2 = ( CM_B2 / W ) + GF_Alpha * B0u &
                            * ( PM_V2 - ( GF_Beta2 / GF_Alpha ) )

   !PRINT*, 'PM_B2: ', PM_B2
 
    PM_B3 = ( CM_B3 / W ) + GF_Alpha * B0u &
                            * ( PM_V3 - ( GF_Beta3 / GF_Alpha ) )

   !PRINT*, 'PM_B3: ', PM_B3
 
    PM_Chi = CM_Chi

    CALL ComputeSoundSpeedFromPrimitive &
           ( PM_D, PM_E, PM_Ne, Cs )
 
   !PRINT*, 'Cs: ', Cs

   !PRINT*, 'Master Function Derivative Bound: ', One - VSq * Cs**2

   !PRINT*

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
    ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi, &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi, &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3, &
      EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    REAL(DP), INTENT(in)    :: &
      CM_D(:), CM_S1(:), CM_S2(:), CM_S3(:), CM_E(:), CM_Ne(:), &
      CM_B1(:), CM_B2(:), CM_B3(:), CM_Chi(:)
    REAL(DP), INTENT(in)    :: &
      GF_Gm11(:), GF_Gm22(:), GF_Gm33(:), &
      GF_Alpha(:), GF_Beta1(:), GF_Beta2(:), GF_Beta3(:)
    REAL(DP), INTENT(out)   :: &
      PM_D(:), PM_V1(:), PM_V2(:), PM_V3(:), PM_E(:), PM_Ne(:), &
      PM_B1(:), PM_B2(:), PM_B3(:), PM_Chi(:)

    INTEGER :: iNX

    DO iNX = 1, SIZE( CM_D )

      CALL ComputePrimitive_Scalar &
             ( CM_D   (iNX), &
               CM_S1  (iNX), &
               CM_S2  (iNX), &
               CM_S3  (iNX), &
               CM_E   (iNX), &
               CM_Ne  (iNX), &
               CM_B1  (iNX), &
               CM_B2  (iNX), &
               CM_B3  (iNX), &
               CM_Chi (iNX), &
               PM_D   (iNX), &
               PM_V1  (iNX), &
               PM_V2  (iNX), &
               PM_V3  (iNX), &
               PM_E   (iNX), &
               PM_Ne  (iNX), &
               PM_B1  (iNX), &
               PM_B2  (iNX), &
               PM_B3  (iNX), &
               PM_Chi (iNX), &
               GF_Gm11(iNX), &
               GF_Gm22(iNX), &
               GF_Gm33(iNX), &
               GF_Alpha(iNX), &
               GF_Beta1(iNX), &
               GF_Beta2(iNX), &
               GF_Beta3(iNX), &
               EvolveOnlyMagnetic )

    END DO

  END SUBROUTINE ComputePrimitive_Vector


  !> Compute conserved variables from primitive variables.
  SUBROUTINE ComputeConserved_Scalar &
    ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi, &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi, &
      GF_Gm11, GF_Gm22, GF_Gm33,  &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3, &
      AM_P, EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    REAL(DP), INTENT(in)  :: PM_D, PM_V1, PM_V2, PM_V3, &
                             PM_E, PM_Ne, PM_B1, PM_B2, &
                             PM_B3, PM_Chi, AM_P, &
                             GF_Gm11, GF_Gm22, GF_Gm33, &
                             GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3
    REAL(DP), INTENT(out) :: CM_D, CM_S1, CM_S2, CM_S3, &
                             CM_E, CM_Ne, CM_B1, CM_B2, &
                             CM_B3, CM_Chi

    REAL(DP) :: VSq, W, B0u, B0d, BSq, h, hStar, p, pStar

   !PRINT*
   !PRINT*, 'Computing conserved variables.'

   !PRINT*, 'PM_D : ', PM_D
   !PRINT*, 'PM_V1: ', PM_V1
   !PRINT*, 'PM_V2: ', PM_V2
   !PRINT*, 'PM_V3: ', PM_V3
   !PRINT*, 'PM_E : ', PM_E
   !PRINT*, 'PM_Ne: ', PM_Ne
   !PRINT*, 'PM_B1: ', PM_B1
   !PRINT*, 'PM_B2: ', PM_B2
   !PRINT*, 'PM_B3: ', PM_B3

    VSq = GF_Gm11 * PM_V1**2 + GF_Gm22 * PM_V2**2 + GF_Gm33 * PM_V3**2
    W = One / SQRT( One - VSq )

   !PRINT*, 'VSq: ', VSq 
   !PRINT*, 'W: ', W

    B0u = ( GF_Gm11 * PM_V1 * PM_B1 &
             + GF_Gm22 * PM_V2 * PM_B2 &
             + GF_Gm33 * PM_V3 * PM_B3 ) &
         / ( GF_Alpha - GF_Gm11 * PM_V1 * GF_Beta1 &
             - GF_Gm22 * PM_V2 * GF_Beta2 &
             - GF_Gm33 * PM_V3 * GF_Beta3 )

    B0d = B0u * ( - GF_Alpha**2 + GF_Gm11 * GF_Beta1**2 &
                                + GF_Gm22 * GF_Beta2**2 &
                                + GF_Gm33 * GF_Beta3**2 ) &
          + ( GF_Gm11 * GF_Beta1 * PM_B1 &
              + GF_Gm22 * GF_Beta2 * PM_B2 &
              + GF_Gm33 * GF_Beta3 * PM_B3 )
                 
   !PRINT*, 'B0u: ', B0u
   !PRINT*, 'B0d: ', B0d

    BSq = B0d * B0u &
          + B0u * ( GF_Gm11 * GF_Beta1 * PM_B1 &
                    + GF_Gm22 * GF_Beta2 * PM_B2 &
                    + GF_Gm33 * GF_Beta3 * PM_B3 ) &
          + ( GF_Gm11 * PM_B1**2 &
              + GF_Gm22 * PM_B2**2 &
              + GF_Gm33 * PM_B3**2 )

   !PRINT*, 'BSq: ', BSq

    h = One + ( PM_E + AM_P ) / PM_D
    p = AM_P
    hStar = h + BSq / PM_D
    pStar = p + BSq / 2.0_DP

   !PRINT*, 'hStar: ', hStar
   !PRINT*, 'pStar: ', pStar

    CM_D   = W * PM_D

    IF( .NOT. EvolveOnlyMagnetic )THEN

     !PRINT*, 'Magnetofluid coupling on for primitive-to-conserved.'

      CM_S1  = hStar * W**2 * PM_D * GF_Gm11 * PM_V1 &
               - GF_Alpha * B0u**2 * ( GF_Gm11 * GF_Beta1 ) &
               - GF_Alpha * B0u * ( GF_Gm11 * PM_B1 )
      CM_S2  = hStar * W**2 * PM_D * GF_Gm22 * PM_V2 &
               - GF_Alpha * B0u**2 * ( GF_Gm22 * GF_Beta2 ) &
               - GF_Alpha * B0u * ( GF_Gm22 * PM_B2 )
      CM_S3  = hStar * W**2 * PM_D * GF_Gm33 * PM_V3 &
               - GF_Alpha * B0u**2 * ( GF_Gm33 * GF_Beta3 ) &
               - GF_Alpha * B0u * ( GF_Gm33 * PM_B3 )
      CM_E   = hStar * W**2 * PM_D - pStar - ( GF_Alpha * B0u )**2 - W * PM_D

    ELSE

     !PRINT*, 'Magnetofluid coupling off for primitive-to-conserved.'

      CM_S1  = h * W**2 * PM_D * GF_Gm11 * PM_V1
      CM_S2  = h * W**2 * PM_D * GF_Gm22 * PM_V2
      CM_S3  = h * W**2 * PM_D * GF_Gm33 * PM_V3
      CM_E   = h * W**2 * PM_D - p - W * PM_D

    END IF

    CM_Ne  = W * PM_Ne
    CM_B1  = -W * GF_Alpha * B0u * ( PM_V1 - ( GF_Beta1 / GF_Alpha  ) ) &
             + W * PM_B1 
    CM_B2  = -W * GF_Alpha * B0u * ( PM_V2 - ( GF_Beta2 / GF_Alpha  ) ) &
             + W * PM_B2
    CM_B3  = -W * GF_Alpha * B0u * ( PM_V3 - ( GF_Beta3 / GF_Alpha  ) ) &
             + W * PM_B3
    CM_Chi = PM_Chi

   !PRINT*, 'CM_D : ', CM_D
   !PRINT*, 'CM_S1: ', CM_S1
   !PRINT*, 'CM_S2: ', CM_S2
   !PRINT*, 'CM_S3: ', CM_S3
   !PRINT*, 'CM_E : ', CM_E
   !PRINT*, 'CM_Ne: ', CM_Ne
   !PRINT*, 'CM_B1: ', CM_B1
   !PRINT*, 'CM_B2: ', CM_B2
   !PRINT*, 'CM_B3: ', CM_B3
   !PRINT*

  END SUBROUTINE ComputeConserved_Scalar


  SUBROUTINE ComputeConserved_Vector &
    ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi, &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi, &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3, &
      AM_P, EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    REAL(DP), INTENT(in)  :: PM_D(:), PM_V1(:), PM_V2(:), PM_V3(:), &
                             PM_E(:), PM_Ne(:), PM_B1(:), PM_B2(:), &
                             PM_B3(:), PM_Chi(:), AM_P(:), &
                             GF_Gm11(:), GF_Gm22(:), GF_Gm33(:), &
                             GF_Alpha(:), GF_Beta1(:), GF_Beta2(:), GF_Beta3(:)
    REAL(DP), INTENT(out) :: CM_D(:), CM_S1(:), CM_S2(:), CM_S3(:), &
                             CM_E(:), CM_Ne(:), CM_B1(:), CM_B2(:), &
                             CM_B3(:), CM_Chi(:)

    INTEGER :: iNX

    DO iNX = 1, SIZE( PM_D )

      CALL ComputeConserved_Scalar &
             ( PM_D (iNX),  &
               PM_V1(iNX),  &
               PM_V2(iNX),  &
               PM_V3(iNX),  &
               PM_E (iNX),  &
               PM_Ne(iNX),  &
               PM_B1(iNX),  &
               PM_B2(iNX),  &
               PM_B3(iNX),  &
               PM_Chi(iNX), &
               CM_D (iNX),  &
               CM_S1(iNX),  &
               CM_S2(iNX),  &
               CM_S3(iNX),  &
               CM_E (iNX),  &
               CM_Ne(iNX),  &
               CM_B1(iNX),  &
               CM_B2(iNX),  &
               CM_B3(iNX),  &
               CM_Chi(iNX), &
               GF_Gm11 (iNX),  &
               GF_Gm22 (iNX),  &
               GF_Gm33 (iNX),  &
               GF_Alpha(iNX),  &
               GF_Beta1(iNX),  &
               GF_Beta2(iNX),  &
               GF_Beta3(iNX),  &
               AM_P (iNX), &
               EvolveOnlyMagnetic )

    END DO

  END SUBROUTINE ComputeConserved_Vector


  SUBROUTINE ComputeMagneticDivergence_MHD_Relativistic &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    REAL(DP), INTENT(in)  :: t
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ApplyBoundaryConditions_MHD &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL InitializeIncrement &
      ( iX_B0, iX_E0, iX_B1, iX_E1 )

    CALL MagneticDivergence_X1_Relativistic &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL MagneticDivergence_X2_Relativistic &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL MagneticDivergence_X3_Relativistic &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

  END SUBROUTINE ComputeMagneticDivergence_MHD_Relativistic


  SUBROUTINE ComputeWeakMagneticDivergence_MHD_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      WeakDiv(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):)

    CALL InitializeIncrement &
      ( iX_B0, iX_E0, iX_B1, iX_E1 )

    CALL WeakMagneticDivergence_X1_Relativistic &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

    CALL WeakMagneticDivergence_X2_Relativistic &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

  END SUBROUTINE ComputeWeakMagneticDivergence_MHD_Relativistic


  SUBROUTINE MagneticDivergence_X1_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(inout) :: &
      D(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)-1:iX_E0(1)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)  :iX_E0(1)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCM_K(nDOFX, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)-1:iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_K(1,1,iX_B0(2)        :iX_E0(2), &
                iX_B0(3)        :iX_E0(3), &
                iX_B0(1)        :iX_E0(1))

    REAL(DP) :: &
      Div(nDOFX_X1,1,iX_B0(2)    :iX_E0(2),   &
                     iX_B0(3)    :iX_E0(3),   &
                     iX_B0(1)-1  :iX_E0(1)+1)

    REAL(DP) :: &
      Vol(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3) )

    REAL(DP) :: &
      Ones(nDOFX_X1)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(2) ; iXP_E0(1) = iX_E0(2)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(1) ; iXP_E0(3) = iX_E0(1)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_MagneticDivergence &
           ( iXP_B0, iXP_E0, nDOFX_X1, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    DO iGF = 1           , nGF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX2,iX3,iX1,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iCM = 1           , nCM
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX2,iX3,iX1,iCM) = U(iNX,iX1,iX2,iX3,iCM    )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1)        , iX_E0(1)
    DO iX3 = iX_B0(3)        , iX_E0(3)
    DO iX2 = iX_B0(2)        , iX_E0(2)

      uDM_K(1,1,iX2,iX3,iX1) = 0.0_DP

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )

    ! --- Compute metric and metric determinant on faces ---

    DO iX1   = iX_B0(1), iX_E0(1) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iNX_X = 1       , nDOFX_X1

      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_2   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_3   ), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

    END DO

    DO iNX_X = 1, nNodesX_X1

      iNX = IndexTableX_F(1,iNX_X)
      iX2 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX1 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX2,iX3,iX1,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX2,iX3,iX1,iCM)

      END DO

      Div(iNX,1,iX2,iX3,iX1) = Half * (uCM_L_nCM(iCM_B1) + uCM_R_nCM(iCM_B1)) &
                               * SqrtGm_F(iNX_X) * dX2(iX2) * dX3(iX3) &
                               * WeightsX_X1(iNX)

    END DO ! iNX_X

    Ones = One

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', 1, nX_K, nDOFX_X1, - One, Ones, nDOFX_X1, &
             Div(1,1,iX_B0(2),iX_B0(3),iX_B0(1) ), &
             nDOFX_X1, Zero, uDM_K, 1 )

    ! --- Contribution from Right Face

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', 1, nX_K, nDOFX_X1, + One, Ones, nDOFX_X1, &
             Div(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1 ), &
             nDOFX_X1, One, uDM_K, 1 )

    Vol = Zero

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       ,     nDOFX

      Vol(iX1,iX2,iX3) = Vol(iX1,iX2,iX3) + WeightsX_q(iNX) &
                                            * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                                            * dX1(iX1) * dX2(iX2) * dX3(iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       , nDOFX

      D(iNX,iX1,iX2,iX3,iDM_Div) = uDM_K(1,1,iX2,iX3,iX1) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1, sWX(1)
    DO iNX = 1, nDOFX

      D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_Div) = 0.0_DP
      D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_Div) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_MagneticDivergence

    END ASSOCIATE

  END SUBROUTINE MagneticDivergence_X1_Relativistic


  SUBROUTINE MagneticDivergence_X2_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(inout) :: &
      D(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)-1:iX_E0(2)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)  :iX_E0(2)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCM_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)-1:iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_K(1,1,iX_B0(1)      :iX_E0(1), &
                iX_B0(3)      :iX_E0(3), &
                iX_B0(2)      :iX_E0(2))

    REAL(DP) :: &
      Div(nDOFX_X1,1,iX_B0(1)    :iX_E0(1),   &
                     iX_B0(3)    :iX_E0(3),   &
                     iX_B0(2)-1  :iX_E0(2)+1)

    REAL(DP) :: &
      Vol(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3) )

    REAL(DP) :: &
      Ones(nDOFX_X2)

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(2) ; iXP_E0(3) = iX_E0(2)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_MagneticDivergence &
           ( iXP_B0, iXP_E0, nDOFX_X2, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    DO iGF = 1           , nGF
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX1,iX3,iX2,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iCM = 1           , nCM
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX1,iX3,iX2,iCM) = U(iNX,iX1,iX2,iX3,iCM    )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B0(2)        , iX_E0(2)
    DO iX3 = iX_B0(3)        , iX_E0(3)
    DO iX1 = iX_B0(1)        , iX_E0(1)

      uDM_K(1,1,iX1,iX3,iX2) = 0.0_DP

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )

    ! --- Compute metric and metric determinant on faces ---

    ! --- Compute metric and metric determinant on faces ---

    DO iX2   = iX_B0(2), iX_E0(2) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iNX_X = 1       , nDOFX_X2

      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_2   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_3   ), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

    END DO

    DO iNX_X = 1, nNodesX_X2

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX2 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX1,iX3,iX2,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX1,iX3,iX2,iCM)

      END DO

      Div(iNX,1,iX1,iX3,iX2) = Half * (uCM_L_nCM(iCM_B2) + uCM_R_nCM(iCM_B2)) &
                               * SqrtGm_F(iNX_X) * dX1(iX1) * dX3(iX3) &
                               * WeightsX_X2(iNX)

    END DO ! iNX_X

    Ones = One

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', 1, nX_K, nDOFX_X2, - One, Ones, nDOFX_X2, &
             Div(1,1,iX_B0(1),iX_B0(3),iX_B0(2) ), &
             nDOFX_X2, Zero, uDM_K, 1 )

    ! --- Contribution from Right Face

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', 1, nX_K, nDOFX_X2, + One, Ones, nDOFX_X2, &
             Div(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1 ), &
             nDOFX_X2, One, uDM_K, 1 )

    Vol = Zero

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       ,     nDOFX

      Vol(iX1,iX2,iX3) = Vol(iX1,iX2,iX3) + WeightsX_q(iNX) &
                                            * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                                            * dX1(iX1) * dX2(iX2) * dX3(iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       , nDOFX

      D(iNX,iX1,iX2,iX3,iDM_Div) = uDM_K(1,1,iX1,iX3,iX2) / Vol(iX1,iX2,iX3) + D(iNX,iX1,iX2,iX3,iDM_Div)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = 1, sWX(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      D(iNX,iX1,iX_B0(2)-iX2,iX3,iDM_Div) = 0.0_DP
      D(iNX,iX1,iX_E0(2)+iX2,iX3,iDM_Div) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_MagneticDivergence

    END ASSOCIATE

  END SUBROUTINE MagneticDivergence_X2_Relativistic


  SUBROUTINE MagneticDivergence_X3_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(inout) :: &
      D(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)

  END SUBROUTINE MagneticDivergence_X3_Relativistic


  SUBROUTINE WeakMagneticDivergence_X1_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(out) :: &
      WeakDiv(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)
    REAL(DP) :: Div_K

    REAL(DP) :: WeakDiv_X1(nDOFX,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1):iX_E0(1))
    REAL(DP) :: Div_q(nDOFX,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1):iX_E0(1))

    REAL(DP) :: &
      Vol(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3) )

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)-1:iX_E0(1)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)  :iX_E0(1)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCM_K(nDOFX, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)-1:iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      Div(nDOFX_X1,iX_B0(2)    :iX_E0(2),   &
                   iX_B0(3)    :iX_E0(3),   &
                   iX_B0(1)    :iX_E0(1)+1)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(2) ; iXP_E0(1) = iX_E0(2)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(1) ; iXP_E0(3) = iX_E0(1)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_MagneticDivergence &
           ( iXP_B0, iXP_E0, nDOFX_X1, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    DO iGF = 1           , nGF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX2,iX3,iX1,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iCM = 1           , nCM
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX2,iX3,iX1,iCM) = U(iNX,iX1,iX2,iX3,iCM    )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1)        , iX_E0(1)
    DO iX3 = iX_B0(3)        , iX_E0(3)
    DO iX2 = iX_B0(2)        , iX_E0(2)
    DO iNX = 1               , nDOFX

      WeakDiv_X1(iNX,iX2,iX3,iX1) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )

    ! --- Compute metric and metric determinant on faces ---

    DO iX1   = iX_B0(1), iX_E0(1) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iNX_X = 1       , nDOFX_X1

      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_2   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_3   ), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

    END DO

    DO iNX_X = 1, nNodesX_X1

      iNX = IndexTableX_F(1,iNX_X)
      iX2 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX1 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX2,iX3,iX1,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX2,iX3,iX1,iCM)

      END DO

      Div(iNX,iX2,iX3,iX1) = Half * (uCM_L_nCM(iCM_B1) + uCM_R_nCM(iCM_B1)) &
                             * SqrtGm_F(iNX_X) * dX2(iX2) * dX3(iX3) &
                             * WeightsX_X1(iNX)

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             Div(1,iX_B0(2),iX_B0(3),iX_B0(1) ), &
             nDOFX_X1, Zero, WeakDiv_X1, nDOFX )

    ! --- Contribution from Right Face

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             Div(1,iX_B0(2),iX_B0(3),iX_B0(1)+1 ), &
             nDOFX_X1, One, WeakDiv_X1, nDOFX )

    DO iNX_K = 1, nNodesX_K

      iNX = IndexTableX_V(1,iNX_K)
      iX2 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX1 = IndexTableX_V(4,iNX_K)

      Div_K = uCM_K(iNX,iX2,iX3,iX1,iCM_B1)

      Div_q(iNX,iX2,iX3,iX1) &
        = Div_K &
            * SqrtGm_K(iNX_K) &
            * dX2(iX2) * dX3(iX3) * WeightsX_q(iNX)

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX, -One, dLXdX1_q, nDOFX, &
             Div_q, nDOFX, One, WeakDiv_X1, nDOFX )

    Vol = Zero

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       ,     nDOFX

      Vol(iX1,iX2,iX3) = Vol(iX1,iX2,iX3) + WeightsX_q(iNX) &
                                            * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                                            * dX1(iX1) * dX2(iX2) * dX3(iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      WeakDiv(iNX,iX1,iX2,iX3) &
        = WeakDiv_X1(iNX,iX2,iX3,iX1) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1, sWX(1)
    DO iNX = 1, nDOFX

      WeakDiv(iNX,iX_B0(1)-iX1,iX2,iX3) = 0.0_DP
      WeakDiv(iNX,iX_E0(1)+iX1,iX2,iX3) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_MagneticDivergence

    END ASSOCIATE

  END SUBROUTINE WeakMagneticDivergence_X1_Relativistic


  SUBROUTINE WeakMagneticDivergence_X2_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, WeakDiv )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(inout) :: &
      WeakDiv(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)
    REAL(DP) :: Div_K

    REAL(DP) :: WeakDiv_X2(nDOFX,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2):iX_E0(2))
    REAL(DP) :: Div_q(nDOFX,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2):iX_E0(2))

    REAL(DP) :: &
      Vol(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(1):iX_E0(3) )

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)-1:iX_E0(2)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(1)  :iX_E0(1),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(2)  :iX_E0(2)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCM_K(nDOFX, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)-1:iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(1)  :iX_E0(1),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(2)  :iX_E0(2)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      Div(nDOFX_X1,iX_B0(1)    :iX_E0(1),   &
                   iX_B0(3)    :iX_E0(3),   &
                   iX_B0(2)-1  :iX_E0(2)+1)

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(1) ; iXP_E0(1) = iX_E0(1)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(2) ; iXP_E0(3) = iX_E0(2)

    ASSOCIATE( dX1 => MeshX(1) % Width, dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_MagneticDivergence &
           ( iXP_B0, iXP_E0, nDOFX_X2, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    DO iGF = 1           , nGF
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX1,iX3,iX2,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iCM = 1           , nCM
    DO iX2 = iX_B0(2) - 1, iX_E0(2) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX1 = iX_B0(1)    , iX_E0(1)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX1,iX3,iX2,iCM) = U(iNX,iX1,iX2,iX3,iCM    )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B0(2)        , iX_E0(2)
    DO iX3 = iX_B0(3)        , iX_E0(3)
    DO iX1 = iX_B0(1)        , iX_E0(1)
    DO iNX = 1               , nDOFX

      WeakDiv_X2(iNX,iX1,iX3,iX2) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_1), nDOFX_X2 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_2), nDOFX_X2 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             G_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iGF_h_3), nDOFX_X2 )

    ! --- Compute metric and metric determinant on faces ---

    ! --- Compute metric and metric determinant on faces ---

    DO iX2   = iX_B0(2), iX_E0(2) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX1   = iX_B0(1), iX_E0(1)
    DO iNX_X = 1       , nDOFX_X2

      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX1,iX3,iX2,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX1,iX3,iX2,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX1,iX3,iX2,iGF_h_1   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_2   ) &
                 * G_F(iNX_X,iX1,iX3,iX2,iGF_h_3   ), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, nX2_X, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               uCM_K(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(1),iX_B0(3),iX_B0(2)  ,iCM), nDOFX_X2 )

    END DO

    DO iNX_X = 1, nNodesX_X2

      iNX = IndexTableX_F(1,iNX_X)
      iX1 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX2 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX1,iX3,iX2,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX1,iX3,iX2,iCM)

      END DO

      Div(iNX,iX1,iX3,iX2) = Half * (uCM_L_nCM(iCM_B2) + uCM_R_nCM(iCM_B2)) &
                             * SqrtGm_F(iNX_X) * dX1(iX1) * dX3(iX3) &
                             * WeightsX_X2(iNX)

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             Div(1,iX_B0(1),iX_B0(3),iX_B0(2) ), &
             nDOFX_X2, Zero, WeakDiv_X2, nDOFX )

    ! --- Contribution from Right Face

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             Div(1,iX_B0(1),iX_B0(3),iX_B0(2)+1 ), &
             nDOFX_X2, One, WeakDiv_X2, nDOFX )

    DO iNX_K = 1, nNodesX_K

      iNX = IndexTableX_V(1,iNX_K)
      iX1 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX2 = IndexTableX_V(4,iNX_K)

      Div_K = uCM_K(iNX,iX1,iX3,iX2,iCM_B2)

      Div_q(iNX,iX1,iX3,iX2) &
        = Div_K &
            * SqrtGm_K(iNX_K) &
            * dX1(iX1) * dX3(iX3) * WeightsX_q(iNX)

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K, nDOFX, -One, dLXdX2_q, nDOFX, &
             Div_q, nDOFX, One, WeakDiv_X2, nDOFX )

    Vol = Zero

    DO iX3 = iX_B0(3),     iX_E0(3)
    DO iX2 = iX_B0(2),     iX_E0(2)
    DO iX1 = iX_B0(1),     iX_E0(1)
    DO iNX = 1       ,     nDOFX

      Vol(iX1,iX2,iX3) = Vol(iX1,iX2,iX3) + WeightsX_q(iNX) &
                                            * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                                            * dX1(iX1) * dX2(iX2) * dX3(iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      WeakDiv(iNX,iX1,iX2,iX3) &
        = WeakDiv(iNX,iX1,iX2,iX3) &
            + WeakDiv_X2(iNX,iX1,iX3,iX2) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = 1, sWX(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      WeakDiv(iNX,iX1,iX_B0(2)-iX2,iX3) = 0.0_DP
      WeakDiv(iNX,iX1,ix_B0(2)+iX2,iX3) = 0.0_DP

    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_MagneticDivergence

    END ASSOCIATE

  END SUBROUTINE WeakMagneticDivergence_X2_Relativistic


  SUBROUTINE InitializeIncrement &
               ( iX_B0, iX_E0, iX_B1, iX_E1 )

    INTEGER,  INTENT(in)           :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    nX        = iX_E0 - iX_B0 + 1 ! Number of Elements per Dimension
    nX_K      = PRODUCT( nX )     ! Number of Elements in Position Space
    nNodesX_K = nDOFX * nX_K      ! Number of Nodes in Elements

    nX_X1 = nX + [1,0,0] ! Number of X1 Faces per Dimension
    nX_X2 = nX + [0,1,0] ! Number of X2 Faces per Dimension
    nX_X3 = nX + [0,0,1] ! Number of X3 Faces per Dimension

    nX1_X = PRODUCT( nX_X1 ) ! Number of X1 Faces
    nX2_X = PRODUCT( nX_X2 ) ! Number of X2 Faces
    nX3_X = PRODUCT( nX_X3 ) ! Number of X3 Faces

    nNodesX_X1 = nDOFX_X1 * nX1_X ! Number of Nodes on X1 Faces
    nNodesX_X2 = nDOFX_X2 * nX2_X ! Number of Nodes on X2 Faces
    nNodesX_X3 = nDOFX_X3 * nX3_X ! Number of Nodes on X3 Faces

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE InitializeIncrement
  

  SUBROUTINE InitializeIncrement_MagneticDivergence &
    ( iXP_B0, iXP_E0, nDOFX_X, G_K, G_F, uCM_K, uCM_L, uCM_R)

    INTEGER, INTENT(in) :: iXP_B0(3), iXP_E0(3) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X              ! nDOFX_X1, ...

    ! --- Geometry Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      G_K (1:nDOFX, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)-1:iXP_E0(3)+1, &
           1:nGF), &
      G_F (1:nDOFX_X, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)  :iXP_E0(3)+1, &
           1:nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCM_K(1:nDOFX, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)-1:iXP_E0(3)+1, &
            1:nCM), &
      uCM_L(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM), &
      uCM_R(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM)

    INTEGER :: iNX, iX1, iX2, iX3, iCM
    INTEGER :: nXP(3), nXP_X(3), nX_X, nNodesX_X, iX_F, iX_V

    nXP       = iXP_E0 - iXP_B0 + 1
    nXP_X     = nXP + [0,0,1]
    nX_X      = PRODUCT( nXP_X )
    nNodesX_X = nDOFX_X * nX_X

    Gm_dd_11_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_SqrtGm  )

    Gm_dd_11_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_SqrtGm  )

    uB1_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B1)
    uB2_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B2)
    uB3_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B3)
    uChi_K(1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_Chi)

    uB1_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B1)
    uB2_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B2)
    uB3_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B3)
    uChi_L(1:nNodesX_X) => uCM_L(:,:,:,:,iCM_Chi)

    uB1_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B1)
    uB2_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B2)
    uB3_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B3)
    uChi_R(1:nNodesX_X) => uCM_R(:,:,:,:,iCM_Chi)

    ALLOCATE( IndexTableX_F(4,nNodesX_X) )
    ALLOCATE( IndexTableX_V(4,nNodesX_K) )

    DO iX3 = iXP_B0(3), iXP_E0(3) + 1
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX_X

      iX_F = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX_X &
               + ( iX2 - iXP_B0(2) ) * nDOFX_X * nXP_X(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX_X * nXP_X(1) * nXP_X(2)

      IndexTableX_F(:,iX_F) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iXP_B0(3), iXP_E0(3)
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX

      iX_V = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX &
               + ( iX2 - iXP_B0(2) ) * nDOFX * nXP(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX * nXP(1) * nXP(2)

      IndexTableX_V(:,iX_V) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeIncrement_MagneticDivergence


  SUBROUTINE FinalizeIncrement_MagneticDivergence

    DEALLOCATE( IndexTableX_F, IndexTableX_V )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY(   Beta_1_K,   Beta_2_K,   Beta_3_K,  Alpha_K )
    NULLIFY(   Beta_1_F,   Beta_2_F,   Beta_3_F,  Alpha_F )
    NULLIFY( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K )
    NULLIFY( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L )
    NULLIFY( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R )

  END SUBROUTINE FinalizeIncrement_MagneticDivergence


  !> Compute primitive variables, pressure, and sound-speed from conserved
  !> variables for a data block.
  SUBROUTINE ComputeFromConserved_MHD_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A, &
      EvolveOnlyMagnetic )

    LOGICAL, INTENT(in) :: EvolveOnlyMagnetic

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iAM

    ! --- Update primitive variables, pressure, and sound speed ---

    DO iAM = 1, nAM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      A(iNX,iX1,iX2,iX3,iAM) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      CALL ComputePrimitive_MHD_Relativistic &
             ( U   (iNX,iX1,iX2,iX3,iCM_D ),        &
               U   (iNX,iX1,iX2,iX3,iCM_S1),        &
               U   (iNX,iX1,iX2,iX3,iCM_S2),        &
               U   (iNX,iX1,iX2,iX3,iCM_S3),        &
               U   (iNX,iX1,iX2,iX3,iCM_E ),        &
               U   (iNX,iX1,iX2,iX3,iCM_Ne),        &
               U   (iNX,iX1,iX2,iX3,iCM_B1),        &
               U   (iNX,iX1,iX2,iX3,iCM_B2),        &
               U   (iNX,iX1,iX2,iX3,iCM_B3),        &
               U   (iNX,iX1,iX2,iX3,iCM_Chi),       &
               P   (iNX,iX1,iX2,iX3,iPM_D ),        &
               P   (iNX,iX1,iX2,iX3,iPM_V1),        &
               P   (iNX,iX1,iX2,iX3,iPM_V2),        &
               P   (iNX,iX1,iX2,iX3,iPM_V3),        &
               P   (iNX,iX1,iX2,iX3,iPM_E ),        &
               P   (iNX,iX1,iX2,iX3,iPM_Ne),        &
               P   (iNX,iX1,iX2,iX3,iPM_B1),        &
               P   (iNX,iX1,iX2,iX3,iPM_B2),        &
               P   (iNX,iX1,iX2,iX3,iPM_B3),        &
               P   (iNX,iX1,iX2,iX3,iPM_Chi),       &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_33),  &
               G   (iNX,iX1,iX2,iX3,iGF_Alpha),     &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_1),    &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_2),    &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_3), &
               EvolveOnlyMagnetic )

      CALL ComputeAuxiliary_Fluid &
             ( P(iNX,iX1,iX2,iX3,iPM_D ), &
               P(iNX,iX1,iX2,iX3,iPM_E ), &
               P(iNX,iX1,iX2,iX3,iPM_Ne), &
               A(iNX,iX1,iX2,iX3,iAM_P ), &
               A(iNX,iX1,iX2,iX3,iAM_T ), &
               A(iNX,iX1,iX2,iX3,iAM_Ye), &
               A(iNX,iX1,iX2,iX3,iAM_S ), &
               A(iNX,iX1,iX2,iX3,iAM_E ), &
               A(iNX,iX1,iX2,iX3,iAM_Gm), &
               A(iNX,iX1,iX2,iX3,iAM_Cs) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_MHD_Relativistic


  !> Loop over all the elements in the spatial domain and compute the minimum
  !> required time-step for numerical stability.
  SUBROUTINE ComputeTimeStep_MHD_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep, &
      UseDivergenceCleaning, EvolveOnlyMagnetic )

    LOGICAL,  INTENT(in)  :: &
      UseDivergenceCleaning, &
      EvolveOnlyMagnetic
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNX, iDimX
    REAL(DP) :: dX(3), dt
    REAL(DP) :: P(nPM), Cs, EigVals(2)

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    TimeStep = HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

     !PRINT*, 'In cell: ', iX1, iX2, iX3
     !PRINT*, 'In node: ', iNX

      dX(1) = dX1(iX1)
      dX(2) = dX2(iX2)
      dX(3) = dX3(iX3)

     !PRINT*, 'Computing the primitive variables for the eigenvalue / timestep calculation.'

      CALL ComputePrimitive_MHD_Relativistic &
             ( U   (iNX,iX1,iX2,iX3,iCM_D ),  &
               U   (iNX,iX1,iX2,iX3,iCM_S1),  &
               U   (iNX,iX1,iX2,iX3,iCM_S2),  &
               U   (iNX,iX1,iX2,iX3,iCM_S3),  &
               U   (iNX,iX1,iX2,iX3,iCM_E ),  &
               U   (iNX,iX1,iX2,iX3,iCM_Ne),  &
               U   (iNX,iX1,iX2,iX3,iCM_B1),  &
               U   (iNX,iX1,iX2,iX3,iCM_B2),  &
               U   (iNX,iX1,iX2,iX3,iCM_B3),  &
               U   (iNX,iX1,iX2,iX3,iCM_Chi), &
               P   (iPM_D ), P(iPM_V1), P(iPM_V2), &
               P   (iPM_V3), P(iPM_E ), P(iPM_Ne), &
               P   (iPM_B1), P(iPM_B2), P(iPM_B3), &
               P   (iPM_Chi),                      &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G   (iNX,iX1,iX2,iX3,iGF_Alpha),    &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_1),   &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_2),   &
               G   (iNX,iX1,iX2,iX3,iGF_Beta_3),   &
               EvolveOnlyMagnetic )

     !PRINT*, 'Computing sound speed.'

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(iPM_D), P(iPM_E), P(iPM_Ne), Cs )
 
     !PRINT*, 'Sound speed: ', Cs

      DO iDimX = 1, nDimsX

       !PRINT*, 'Computing the ', iDimX, ' eigenvalue.'

        EigVals &
          = Eigenvalues_MHD_Relativistic &
              ( P(iPM_V1+(iDimX-1)), Cs, &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11+(iDimX-1)), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_1+(iDimX-1)), &
                P(iPM_D), P(iPM_V1), P(iPM_V2), P(iPM_V3), &
                P(iPM_E), P(iPM_Ne), &
                P(iPM_B1), P(iPM_B2), P(iPM_B3), P(iPM_Chi), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                G(iNX,iX1,iX2,iX3,iGF_Alpha), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_1), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_2), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_3), &
                UseDivergenceCleaning )

       !PRINT*, 'The eigenvalues are: ', EigVals

        dt = dX(iDimX) / MAX( SqrtTiny, MAXVAL( ABS( EigVals ) ) )

        TimeStep = MIN( TimeStep, dt )

      END DO

    END DO
    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

   !PRINT*, 'The timestep is: ', TimeStep

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_MHD_Relativistic


  FUNCTION Eigenvalues_MHD_Relativistic &
    ( Vi, Cs, Gmii, Shifti, D, V1, V2, V3, E, Ne, &
      B1, B2, B3, Chi, Gm11, Gm22, Gm33, &
      Lapse, Shift1, Shift2, Shift3, &
      UseDivergenceCleaning )

    LOGICAL,  INTENT(in) :: UseDivergenceCleaning
    REAL(DP), INTENT(in) :: Vi, Cs, Gmii, Shifti, D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, P, h, B0u, B0d, BSq, &
                Ca, aSq, Eigenvalues_MHD_Relativistic(2)

   !PRINT*
   !PRINT*, 'Input to Eigvenalues_MHD_Relativistic'
   !PRINT*, '-------------------------------------'
   !PRINT*, 'Vi: ', Vi
   !PRINT*, 'Cs: ', Cs
   !PRINT*, 'Gmii: ', Gmii
   !PRINT*, 'D: ', D
   !PRINT*, 'V1: ', V1
   !PRINT*, 'V2: ', V2
   !PRINT*, 'V3: ', V3
   !PRINT*, 'E:  ', E
   !PRINT*, 'Ne: ', Ne
   !PRINT*, 'B1: ', B1
   !PRINT*, 'B2: ', B2
   !PRINT*, 'B3: ', B3
   !PRINT*, 'Chi: ', Chi
   !PRINT*, 'Gm11: ', Gm11
   !PRINT*, 'Gm22: ', Gm22
   !PRINT*, 'Gm33: ', Gm33
   !PRINT*, 'Lapse: ', Lapse
   !PRINT*, 'Shift1: ', Shift1
   !PRINT*, 'Shift2: ', Shift2
   !PRINT*, 'Shift3: ', Shift3
   !PRINT*, '-------------------------------------' 
   !PRINT*

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = One / SQRT( One - VSq )

   !PRINT*
   !PRINT*, 'Inside Eigenvalues_MHD_Relativistic'
   !PRINT*, '-----------------------------------'
   !PRINT*, 'Lorentz factor: ', W

    CALL ComputePressureFromPrimitive &
           ( D, E, Ne, P )
 
   !PRINT*, 'Pressure: ', P

    h = One + E + P / D

   !PRINT*, 'Non-magnetic enthalpy: ', h

    B0u = ( Gm11 * V1 * B1 &
            + Gm22 * V2 * B2 &
            + Gm33 * V3 * B3 ) &
          / ( Lapse &
              - Gm11 * V1 * Shift1 &
              - Gm22 * V2 * Shift2 &
              - Gm33 * V3 * Shift3 )

    B0d =  - ( Lapse / W ) &
             * ( Gm11 * V1 * B1 &
                 + Gm22 * V2 * B2 &
                 + Gm33 * V3 * B3 )

    BSq = B0d * B0u &
          + B0u * ( Gm11 * Shift1 * B1 &
                    + Gm22 * Shift2 * B2 &
                    + Gm33 * Shift3 * B3 ) &
          + ( Gm11 * B1**2 &
              + Gm22 * B2**2 &
              + Gm33 * B3**2 )

    Ca = SQRT( BSq / ( D * h + BSq ) )

   !PRINT*, 'The Alfven speed is: ', Ca

    aSq = Cs**2 + Ca**2 - Cs**2 * Ca**2

   !PRINT*, 'aSq is: ', aSq
   !PRINT*, '-----------------------------------'
   !PRINT*

    IF( UseDivergenceCleaning )THEN

      Eigenvalues_MHD_Relativistic(1) = One
      Eigenvalues_MHD_Relativistic(2) = -One

    ELSE
 
      ! Estimate of max/min fast magnetosonic
      ! eigenvalues from Del Zanna et al. (2007)

      Eigenvalues_MHD_Relativistic(1) &
        = ( ( One - aSq ) * Vi &
            + SQRT( aSq * ( One - VSq ) &
                    * ( ( One - VSq * aSq ) * ( One / Gmii ) &
                    - ( One - aSq ) * Vi**2 ) ) ) &
          / ( One - VSq * aSq )

      Eigenvalues_MHD_Relativistic(1) &
        = Lapse * Eigenvalues_MHD_Relativistic(1) - Shifti

      Eigenvalues_MHD_Relativistic(2) &
         = ( ( One - aSq ) * Vi &
            - SQRT( aSq * ( One - VSq ) &
                    * ( ( One - VSq * aSq ) * ( One / Gmii ) &
                    - ( One - aSq ) * Vi**2 ) ) ) &
          / ( One - VSq * aSq )

      Eigenvalues_MHD_Relativistic(2) &
        = Lapse * Eigenvalues_MHD_Relativistic(2) - Shifti

    END IF

    RETURN
  END FUNCTION Eigenvalues_MHD_Relativistic


  !> Compute the physical flux in the X1-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X1_MHD_Relativistic &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, &
      P, Gm11, Gm22, Gm33, Lapse, & 
      Shift1, Shift2, Shift3, &
      UseDivergenceCleaning )

    LOGICAL,  INTENT(in) :: UseDivergenceCleaning
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, P, &
                            Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, Pstar, h, hStar, B0u, B0d, BSq, Flux_X1_MHD_Relativistic(nCM)


  !PRINT*
  !PRINT*, 'Input to Flux_X1_MHD_Relativistic Routine'
  !PRINT*, '-----------------------------------------'
  !PRINT*, 'D: ', D
  !PRINT*, 'V1: ', V1
  !PRINT*, 'V2: ', V2
  !PRINT*, 'V3: ', V3
  !PRINT*, 'E: ', E
  !PRINT*, 'P: ', P
  !PRINT*, 'Ne: ', Ne
  !PRINT*, 'B1: ', B1
  !PRINT*, 'B2: ', B2
  !PRINT*, 'B3: ', B3
  !PRINT*, 'Chi: ', Chi
  !PRINT*, 'Gm11: ', Gm11
  !PRINT*, 'Gm22: ', Gm22
  !PRINT*, 'Gm33: ', Gm33
  !PRINT*, 'Lapse: ', Lapse
  !PRINT*, 'Shift1: ', Shift1
  !PRINT*, 'Shift2: ', Shift2
  !PRINT*, 'Shift3: ', Shift3
  !PRINT*, '-----------------------------------------'
  !PRINT*

    VSq   = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2 
    W   = One / SQRT( One - VSq )

    B0u = ( Gm11 * V1 * B1 &
            + Gm22 * V2 * B2 &
            + Gm33 * V3 * B3 ) &
          / ( Lapse &
              - Gm11 * V1 * Shift1 &
              - Gm22 * V2 * Shift2 &
              - Gm33 * V3 * Shift3 )

    B0d =  - ( Lapse / W ) &
             * ( Gm11 * V1 * B1 &
                 + Gm22 * V2 * B2 &
                 + Gm33 * V3 * B3 )

    BSq = B0d * B0u &
          + B0u * ( Gm11 * Shift1 * B1 &
                    + Gm22 * Shift2 * B2 &
                    + Gm33 * Shift3 * B3 ) &
          + ( Gm11 * B1**2 &
              + Gm22 * B2**2 &
              + Gm33 * B3**2 )

    hStar = One + ( E + P ) / D + BSq / D
    pStar = P + BSq / 2.0_DP

    Flux_X1_MHD_Relativistic(iCM_D)  &
      = D * W * ( V1 - Shift1 / Lapse )

   !PRINT*, 'iCM_D Flux: ', Flux_X1_MHD_Relativistic(iCM_D)

    Flux_X1_MHD_Relativistic(iCM_S1) &
      = D * hStar * W**2 * Gm11 * V1  * ( V1 - Shift1 / Lapse ) + pStar &
        - Gm11 * ( B1**2 + B0u * B1 * Shift1 &
                   + ( Lapse * B0u )**2 * Shift1**2 )

   !PRINT*, 'iCM_S1 Flux: ', Flux_X1_MHD_Relativistic(iCM_S1)

    Flux_X1_MHD_Relativistic(iCM_S2) &
      = D * hStar * W**2 * Gm22 * V2  * ( V1 - Shift1 / Lapse ) &
        - Gm22 * ( B1 * B2 + B0u * B1 * Shift2 &
                   + ( Lapse * B0u )**2 * Shift1 * Shift2 )

   !PRINT*, 'iCM_S2 Flux: ', Flux_X1_MHD_Relativistic(iCM_S2)

    Flux_X1_MHD_Relativistic(iCM_S3) &
      = D * hStar * W**2 * Gm33 * V3  * ( V1 - Shift1 / Lapse ) &
        - Gm33 * ( B1 * B3 + B0u * B1 * Shift3 &
                   + ( Lapse * B0u )**2 * Shift1 * Shift3 )

   !PRINT*, 'iCM_S3 Flux: ', Flux_X1_MHD_Relativistic(iCM_S3)

    Flux_X1_MHD_Relativistic(iCM_E)  &
      = D * W * ( hStar * W - One ) * ( V1 - Shift1 / Lapse ) &
        + Shift1 / Lapse * pStar &
        - Lapse * B0u * B1 + Two * Lapse * (B0u)**2 * Shift1

   !PRINT*, 'iCM_E Flux: ', Flux_X1_MHD_Relativistic(iCM_E)

    Flux_X1_MHD_Relativistic(iCM_Ne) &
      = Ne * W * ( V1 - Shift1 / Lapse )

   !PRINT*, 'iCM_Ne Flux: ', Flux_X1_MHD_Relativistic(iCM_Ne)

    IF( UseDivergenceCleaning )THEN

      Flux_X1_MHD_Relativistic(iCM_B1) &
        = Lapse * W * B0u * ( V1 - ( Shift1 / Lapse ) ) * ( Shift1 / Lapse ) &
          - W * ( Shift1 / Lapse ) * B1 &
          + ( Chi / Gm11 )

      Flux_X1_MHD_Relativistic(iCM_B2) &
        = Lapse * W * B0u * ( V1 - ( Shift1 / Lapse ) ) * ( Shift2 / Lapse ) &
          + W * V1 * B2 - W * V2 * B1 &
          - W * ( Shift1 / Lapse ) * B2

      Flux_X1_MHD_Relativistic(iCM_B3) &
        = Lapse * W * B0u * ( V1 - ( Shift1 / Lapse ) ) * ( Shift3 / Lapse ) &
          + W * V1 * B3 - W * V3 * B1 &
          - W * ( Shift1 / Lapse ) * B3

      Flux_X1_MHD_Relativistic(iCM_Chi) &
        = - Lapse * W * B0u * ( V1 - ( Shift1 / Lapse ) ) + W * B1 - Chi * ( Shift1 / Lapse )

    ELSE

      Flux_X1_MHD_Relativistic(iCM_B1) &
        = 0.0_DP

     !PRINT*, 'iCM_B1 Flux: ', Flux_X1_MHD_Relativistic(iCM_B1)

      Flux_X1_MHD_Relativistic(iCM_B2) &
        = W * ( V1 - Shift1 / Lapse ) * B2 - W * ( V2 - Shift2 / Lapse ) * B1

      Flux_X1_MHD_Relativistic(iCM_B3) &
        = W * ( V1 - Shift1 / Lapse ) * B3 - W * ( V3 - Shift3 / Lapse ) * B1

      Flux_X1_MHD_Relativistic(iCM_Chi) &
        = 0.0_DP

    END IF

    RETURN
  END FUNCTION Flux_X1_MHD_Relativistic


  !> Compute the physical flux in the X2-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X2_MHD_Relativistic &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, &
      P, Gm11, Gm22, Gm33, Lapse, & 
      Shift1, Shift2, Shift3, &
      UseDivergenceCleaning )

    LOGICAL,  INTENT(in) :: UseDivergenceCleaning
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, P, &
                            Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, Pstar, h, hStar, B0u, B0d, BSq, Flux_X2_MHD_Relativistic(nCM)


  !PRINT*
  !PRINT*, 'Input to Flux_X1_MHD_Relativistic Routine'
  !PRINT*, '-----------------------------------------'
  !PRINT*, 'D: ', D
  !PRINT*, 'V1: ', V1
  !PRINT*, 'V2: ', V2
  !PRINT*, 'V3: ', V3
  !PRINT*, 'E: ', E
  !PRINT*, 'P: ', P
  !PRINT*, 'Ne: ', Ne
  !PRINT*, 'B1: ', B1
  !PRINT*, 'B2: ', B2
  !PRINT*, 'B3: ', B3
  !PRINT*, 'Chi: ', Chi
  !PRINT*, 'Gm11: ', Gm11
  !PRINT*, 'Gm22: ', Gm22
  !PRINT*, 'Gm33: ', Gm33
  !PRINT*, 'Lapse: ', Lapse
  !PRINT*, 'Shift1: ', Shift1
  !PRINT*, 'Shift2: ', Shift2
  !PRINT*, 'Shift3: ', Shift3
  !PRINT*, '-----------------------------------------'
  !PRINT*

    VSq   = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2 
    W   = One / SQRT( One - VSq )

    B0u = ( Gm11 * V1 * B1 &
            + Gm22 * V2 * B2 &
            + Gm33 * V3 * B3 ) &
          / ( Lapse &
              - Gm11 * V1 * Shift1 &
              - Gm22 * V2 * Shift2 &
              - Gm33 * V3 * Shift3 )

    B0d =  - ( Lapse / W ) &
             * ( Gm11 * V1 * B1 &
                 + Gm22 * V2 * B2 &
                 + Gm33 * V3 * B3 )

    BSq = B0d * B0u &
          + B0u * ( Gm11 * Shift1 * B1 &
                    + Gm22 * Shift2 * B2 &
                    + Gm33 * Shift3 * B3 ) &
          + ( Gm11 * B1**2 &
              + Gm22 * B2**2 &
              + Gm33 * B3**2 )

    hStar = One + ( E + P ) / D + BSq / D
    pStar = P + BSq / 2.0_DP

    Flux_X2_MHD_Relativistic(iCM_D)  &
      = D * W * ( V2 - Shift2 / Lapse )

   !PRINT*, 'iCM_D Flux: ', Flux_X2_MHD_Relativistic(iCM_D)

    Flux_X2_MHD_Relativistic(iCM_S1) &
      = D * hStar * W**2 * Gm11 * V1  * ( V2 - Shift2 / Lapse ) &
        - Gm11 * ( B2 * B1 + B0u * B2 * Shift1 &
                   + ( Lapse * B0u )**2 * Shift2 * Shift1 )

   !PRINT*, 'iCM_S1 Flux: ', Flux_X2_MHD_Relativistic(iCM_S1)

    Flux_X2_MHD_Relativistic(iCM_S2) &
      = D * hStar * W**2 * Gm22 * V2  * ( V2 - Shift2 / Lapse ) + pStar &
        - Gm22 * ( B2**2 + B0u * B2 * Shift2 &
                   + ( Lapse * B0u )**2 * Shift2**2 )

   !PRINT*, 'iCM_S2 Flux: ', Flux_X2_MHD_Relativistic(iCM_S2)

    Flux_X2_MHD_Relativistic(iCM_S3) &
      = D * hStar * W**2 * Gm33 * V3  * ( V2 - Shift2 / Lapse ) &
        - Gm33 * ( B2 * B3 + B0u * B2 * Shift3 &
                   + ( Lapse * B0u )**2 * Shift2 * Shift3 )

   !PRINT*, 'iCM_S3 Flux: ', Flux_X2_MHD_Relativistic(iCM_S3)

    Flux_X2_MHD_Relativistic(iCM_E)  &
      = D * W * ( hStar * W - One ) * ( V2 - Shift2 / Lapse ) &
        + Shift2 / Lapse * pStar &
        - Lapse * B0u * B2 + Two * Lapse * (B0u)**2 * Shift2

   !PRINT*, 'iCM_E Flux: ', Flux_X2_MHD_Relativistic(iCM_E)

    Flux_X2_MHD_Relativistic(iCM_Ne) &
      = Ne * W * ( V2 - Shift2 / Lapse )

   !PRINT*, 'iCM_Ne Flux: ', Flux_X2_MHD_Relativistic(iCM_Ne)

    IF( UseDivergenceCleaning )THEN

      Flux_X2_MHD_Relativistic(iCM_B1) &
        = Lapse * W * B0u * ( V2 - ( Shift2 / Lapse ) ) * ( Shift1 / Lapse ) &
          + W * V2 * B1 - W * V1 * B2 &
          - W * ( Shift2 / Lapse ) * B1

      Flux_X2_MHD_Relativistic(iCM_B2) &
        = Lapse * W * B0u * ( V2 - ( Shift2 / Lapse ) ) * ( Shift2 / Lapse ) &
          - W * ( Shift2 / Lapse ) * B2 &
          + ( Chi / Gm22 )

      Flux_X2_MHD_Relativistic(iCM_B3) &
        = Lapse * W * B0u * ( V2 - ( Shift2 / Lapse ) ) * ( Shift3 / Lapse ) &
          + W * V2 * B3 - W * V3 * B2 &
          - W * ( Shift2 / Lapse ) * B3

      Flux_X2_MHD_Relativistic(iCM_Chi) &
        = - Lapse * W * B0u * ( V2 - ( Shift2 / Lapse ) )  + W * B2 - Chi * ( Shift2 / Lapse )

    ELSE

      Flux_X2_MHD_Relativistic(iCM_B1) &
        = W * ( V2 - Shift2 / Lapse ) * B1 - W * ( V1 - Shift1 / Lapse ) * B2

     !PRINT*, 'iCM_B2 Flux: ', Flux_X2_MHD_Relativistic(iCM_B1)

      Flux_X2_MHD_Relativistic(iCM_B2) &
        = 0.0_DP

      Flux_X2_MHD_Relativistic(iCM_B3) &
        = W * ( V2 - Shift2 / Lapse ) * B3 - W * ( V3 - Shift3 / Lapse ) * B2

      Flux_X2_MHD_Relativistic(iCM_Chi) &
        = 0.0_DP

    END IF

    RETURN
  END FUNCTION Flux_X2_MHD_Relativistic


  !> Compute the physical flux in the X2-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X3_MHD_Relativistic &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, &
      P, Gm11, Gm22, Gm33, Lapse, & 
      Shift1, Shift2, Shift3, &
      UseDivergenceCleaning )

    LOGICAL,  INTENT(in) :: UseDivergenceCleaning
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, P, &
                            Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, Pstar, h, hStar, B0u, B0d, BSq, Flux_X3_MHD_Relativistic(nCM)


  !PRINT*
  !PRINT*, 'Input to Flux_X1_MHD_Relativistic Routine'
  !PRINT*, '-----------------------------------------'
  !PRINT*, 'D: ', D
  !PRINT*, 'V1: ', V1
  !PRINT*, 'V2: ', V2
  !PRINT*, 'V3: ', V3
  !PRINT*, 'E: ', E
  !PRINT*, 'P: ', P
  !PRINT*, 'Ne: ', Ne
  !PRINT*, 'B1: ', B1
  !PRINT*, 'B2: ', B2
  !PRINT*, 'B3: ', B3
  !PRINT*, 'Chi: ', Chi
  !PRINT*, 'Gm11: ', Gm11
  !PRINT*, 'Gm22: ', Gm22
  !PRINT*, 'Gm33: ', Gm33
  !PRINT*, 'Lapse: ', Lapse
  !PRINT*, 'Shift1: ', Shift1
  !PRINT*, 'Shift2: ', Shift2
  !PRINT*, 'Shift3: ', Shift3
  !PRINT*, '-----------------------------------------'
  !PRINT*

    VSq   = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2 
    W   = One / SQRT( One - VSq )

    B0u = ( Gm11 * V1 * B1 &
            + Gm22 * V2 * B2 &
            + Gm33 * V3 * B3 ) &
          / ( Lapse &
              - Gm11 * V1 * Shift1 &
              - Gm22 * V2 * Shift2 &
              - Gm33 * V3 * Shift3 )

    B0d =  - ( Lapse / W ) &
             * ( Gm11 * V1 * B1 &
                 + Gm22 * V2 * B2 &
                 + Gm33 * V3 * B3 )

    BSq = B0d * B0u &
          + B0u * ( Gm11 * Shift1 * B1 &
                    + Gm22 * Shift2 * B2 &
                    + Gm33 * Shift3 * B3 ) &
          + ( Gm11 * B1**2 &
              + Gm22 * B2**2 &
              + Gm33 * B3**2 )

    hStar = One + ( E + P ) / D + BSq / D
    pStar = P + BSq / 2.0_DP

    Flux_X3_MHD_Relativistic(iCM_D)  &
      = D * W * ( V3 - Shift3 / Lapse )

   !PRINT*, 'iCM_D Flux: ', Flux_X2_MHD_Relativistic(iCM_D)

    Flux_X3_MHD_Relativistic(iCM_S1) &
      = D * hStar * W**2 * Gm11 * V1  * ( V3 - Shift3 / Lapse ) &
        - Gm11 * ( B3 * B1 + B0u * B3 * Shift1 &
                   + ( Lapse * B0u )**2 * Shift3 * Shift1 )

   !PRINT*, 'iCM_S1 Flux: ', Flux_X2_MHD_Relativistic(iCM_S1)

    Flux_X3_MHD_Relativistic(iCM_S2) &
      = D * hStar * W**2 * Gm22 * V2  * ( V3 - Shift3 / Lapse ) &
        - Gm22 * ( B3 * B2 + B0u * B3 * Shift2 &
                   + ( Lapse * B0u )**2 * Shift3 * Shift2 )

   !PRINT*, 'iCM_S2 Flux: ', Flux_X2_MHD_Relativistic(iCM_S2)

    Flux_X3_MHD_Relativistic(iCM_S3) &
      = D * hStar * W**2 * Gm33 * V3  * ( V3 - Shift3 / Lapse ) + pStar &
        - Gm33 * ( B3**2 + B0u * B3 * Shift3 &
                   + ( Lapse * B0u )**2 * Shift3**2 )

   !PRINT*, 'iCM_S3 Flux: ', Flux_X2_MHD_Relativistic(iCM_S3)

    Flux_X3_MHD_Relativistic(iCM_E)  &
      = D * W * ( hStar * W - One ) * ( V3 - Shift3 / Lapse ) &
        + Shift3 / Lapse * pStar &
        - Lapse * B0u * B3 + Two * Lapse * (B0u)**2 * Shift3

   !PRINT*, 'iCM_E Flux: ', Flux_X2_MHD_Relativistic(iCM_E)

    Flux_X3_MHD_Relativistic(iCM_Ne) &
      = Ne * W * ( V3 - Shift3 / Lapse )

   !PRINT*, 'iCM_Ne Flux: ', Flux_X2_MHD_Relativistic(iCM_Ne)

    IF( UseDivergenceCleaning )THEN

      Flux_X3_MHD_Relativistic(iCM_B1) &
        = Lapse * W * B0u * ( V3 - ( Shift3 / Lapse ) ) * ( Shift1 / Lapse ) &
          + W * V3 * B1 - W * V1 * B3 &
          - W * ( Shift3 / Lapse ) * B1

      Flux_X3_MHD_Relativistic(iCM_B2) &
        = Lapse * W * B0u * ( V3 - ( Shift3 / Lapse ) ) * ( Shift2 / Lapse ) &
          + W * V3 * B2 - W * V2 * B3 &
          - W * ( Shift3 / Lapse ) * B2 &
          + ( Chi / Gm33 )

      Flux_X3_MHD_Relativistic(iCM_B3) &
        = Lapse * W * B0u * ( V3 - ( Shift3 / Lapse ) ) * ( Shift3 / Lapse ) &
          - W * ( Shift3 / Lapse ) * B3

      Flux_X3_MHD_Relativistic(iCM_Chi) &
        = - Lapse * W * B0u * ( V3 - ( Shift3 / Lapse ) )  + W * B3 - Chi * ( Shift3 / Lapse )

    ELSE

      Flux_X3_MHD_Relativistic(iCM_B1) &
        = W * ( V3 - Shift3 / Lapse ) * B1 - W * ( V1 - Shift1 / Lapse ) * B3

     !PRINT*, 'iCM_B2 Flux: ', Flux_X2_MHD_Relativistic(iCM_B1)

      Flux_X3_MHD_Relativistic(iCM_B2) &
        = W * ( V3 - Shift3 / Lapse ) * B2 - W * ( V2 - Shift2 / Lapse ) * B3

      Flux_X3_MHD_Relativistic(iCM_B3) &
        = 0.0_DP

      Flux_X3_MHD_Relativistic(iCM_Chi) &
        = 0.0_DP

    END IF

    RETURN
  END FUNCTION Flux_X3_MHD_Relativistic


  !> Compute the Harten-Lax-van-Leer numerical flux at a given element
  !> interface, in a given dimension.
  FUNCTION NumericalFlux_HLL_MHD_Relativistic &
    ( uL, uR, fL, fR, aP, aM )

    REAL(DP), INTENT(in) :: uL(nCM), uR(nCM), fL(nCM), fR(nCM), aP, aM

    REAL(DP) :: NumericalFlux_HLL_MHD_Relativistic(nCM)

    NumericalFlux_HLL_MHD_Relativistic &
      = ( aP * fL + aM * fR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL_MHD_Relativistic


  ! --- Auxiliary utilities for ComputePrimitive ---

  SUBROUTINE ComputeFunMu &
               ( D_bar, Ne_bar, q, r, bSq, rb, v0, mu, f )

    REAL(DP), INTENT(in) :: D_bar, Ne_bar
    REAL(DP), INTENT(in) :: q
    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: v0
    REAL(DP), INTENT(in) :: mu

    REAL(DP), INTENT(out) :: f

    REAL(DP) :: x, r_barSq, q_bar
    REAL(DP) :: v, W, rho, rho_e, eps
    REAL(DP) :: P, a, h
    REAL(DP) :: nu, nua, nub

   !PRINT*, 'Inside function evaluation.'

   !PRINT*, 'mu: ', mu

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

   !PRINT*, 'x: ', x

    ! --- Eqs. 38 and 39                              ---
    ! --- Note: Eq. 39 rewritten using Eqs. 25 and 30 ---
    ! --- to avoid B = 0 degeneracy.                  ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

   !PRINT*, 'r_barSq: ', r_barSq

   !PRINT*, 'SQRT( r_barSq ): ', SQRT( r_barSq )

    ! See lines 178 and 146 of RePrimAnd's con2prim_imhd.cc file (think this is correct).

    q_bar = q - Half * bSq - Half * mu**2 * x**2 * ( bSq * r**2 - rb**2 )

    !PRINT*, 'q_bar: ', q_bar

    ! --- Eq. 40 ---

   !PRINT*, 'v0: ', v0

    v = MIN( mu * SQRT( r_barSq ), v0 )

   !PRINT*, 'v: ', v

    W = One / SQRT( One - v**2 )

   !PRINT*, 'W: ', W

    ! --- Eq. 41 ---

    rho = D_bar / W

   !PRINT*, 'rho: ', rho

    ! Use Line 229 of RePrimAnd's con2prim_imhd.cc to
    ! limit rho to physical values?????

    rho_e = Ne_bar / W

   !PRINT*, 'rho_e: ', rho_e

    ! --- Eq. 42 ---

    ! See again line 203 of RePrimAnd's con2prim_imhd.cc

    eps = W * ( q_bar - mu * r_barSq * ( One - mu * W / ( One + W ) ) )

   !PRINT*, 'eps:   ', eps

    ! --- Eq. 43 ---

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rho, eps, AtomicMassUnit * rho_e / rho, P )

    a = P / ( rho * ( One + eps ) )

   !PRINT*, 'a: ', a

    ! Lines 240-244 of RePrimAnd's con2prim_imhd.cc

    h = (One + eps) * (One + a)

   !PRINT*, 'h: ', h

    ! --- Eqs. 46-48 ---

    nua = h / W

   !PRINT*, 'nua: ', nua

    nub = ( One + a ) * ( One + q_bar - mu * r_barSq )

   !PRINT*, 'nub: ', nub

    nu = MAX( nua, nub )

   !PRINT*, 'nu: ', nu

    ! --- Eq. 44 --- 

   !PRINT*, 'mu: ', mu
   !PRINT*, '1 / ( nu + r_barSq * mu ): ', One / ( nu + r_barSq * mu )

    f = mu - One / ( nu + r_barSq * mu )

  END SUBROUTINE ComputeFunMu


  SUBROUTINE ComputeAuxiliaryFunMu &
               ( r, bSq, rb, h0, mu, f )

    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0
    REAL(DP), INTENT(in) :: mu

    REAL(DP), INTENT(out) :: f

    REAL(DP) :: x, r_barSq

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

    !*, 'x: ', x

    ! --- Eq. 38 ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

    !*, 'r_barSq: ', r_barSq

    ! --- Eq. 49 --- 

    f = mu * SQRT( h0**2 + r_barSq ) - One

  END SUBROUTINE ComputeAuxiliaryFunMu


  SUBROUTINE ComputeAuxiliaryFunDerivativeMu &
               ( r, bSq, rb, h0, mu, df )

    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0
    REAL(DP), INTENT(in) :: mu

    REAL(DP), INTENT(out) :: df

    REAL(DP) :: x, dx, r_barSq, dr_barSq

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

    dx = -bSq / ( One + mu * bSq )**2 

    ! --- Eq. 38 ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

    dr_barSq = Two * r**2 * x * dx + ( x + mu * dx + x**2 + Two * mu * x * dx ) * rb**2

    df = SQRT( h0**2 + r_barSq ) + ( mu / Two ) * ( dr_barSq / SQRT( h0**2 + r_barSq ) ) 

  END SUBROUTINE 


  SUBROUTINE SolveMuBound_NewtonRaphson &
               ( r, bSq, rb, h0, mu )

    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0
    REAL(DP), INTENT(inout) :: mu

    REAL(DP) :: muO, muN, dmu
    REAL(DP) :: f, df

    LOGICAL :: CONVERGED
    INTEGER :: ITER
    REAL(DP), PARAMETER :: Tolmu = 1.0d-08
    REAL(DP), PARAMETER :: Tolf = 1.0d-08
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolmu ) / LOG( Two ) )

   !PRINT*, 'Solving for the upper bound with the Newton-Raphson method.'

    muO = One / h0

   !PRINT*, 'Initial guess: ', One / h0

   !PRINT*, 'MAX_ITER: ', MAX_ITER

    CONVERGED = .FALSE.
    ITER = 0

    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

     !PRINT*
     !PRINT*, 'ITER: ', ITER
     !PRINT*

     !PRINT*, 'muO: ', muO

     !PRINT*, 'Computing the auxiliary function for mu = : ', muO

      CALL ComputeAuxiliaryFunMu( r, bSq, rb, h0, muO, f )

     !PRINT*, 'f: ', f

     !PRINT*, 'Computing the auxiliary function derivative for mu = ', muO

      CALL ComputeAuxiliaryFunDerivativeMu( r, bSq, rb, h0, muO, df )

     !PRINT*, 'df: ', df

      muN = muO - f / df 

     !PRINT*, 'muN: ', muN

      dmu = muN - muO

     !PRINT*, 'dmu: ', dmu

      IF( ABS ( dmu / muO ) .LT. Tolmu ) CONVERGED = .TRUE.

      muO = muN

    END DO

    mu = muO  

  END SUBROUTINE

  
  SUBROUTINE SolveMu_Bisection &
               ( D_bar, Ne_bar, q, r, &
                 bSq, rb, h0, v0, mu )

    REAL(DP), INTENT(in) :: D_bar, Ne_bar, q
    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0, v0

    REAL(DP), INTENT(out) :: mu

    LOGICAL             :: CONVERGED
    INTEGER             :: ITERATION
    REAL(DP)            :: mua, mub, muc, dmu
    REAL(DP)            :: fa, fb, fc
    REAL(DP), PARAMETER :: dmu_min = 1.0d-08
    INTEGER,  PARAMETER :: MAX_IT = 4 - INT( LOG( dmu_min) / LOG( Two ) )

    IF( r < h0 ) THEN

     !PRINT*, 'r < h0'

      mua = Zero
      mub = One / h0

    ELSE

     !PRINT*, 'r >= h0'

      mua = Zero

      CALL SolveMuBound_NewtonRaphson( r, bSq, rb, h0, mub )

    END IF

   !PRINT*

   !PRINT*, 'Computing function for lower bound: ', mua 

   !PRINT*, '------------------------------------------------------------'

    CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, mua, fa )

   !PRINT*, 'fa: ', fa

   !PRINT*

   !PRINT*, 'Computing function for upper bound: ', mub

   !PRINT*, '------------------------------------------------------------'    

    CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, mub, fb )

   !PRINT*, 'fb: ', fb

    ! --- Check that sign of FunZ changes across bounds ---

    IF( .NOT. fa * fb .LT. 0 ) THEN
      PRINT*, 'Cannot perform bisection in primitive recovery!'
      STOP
    END IF

    dmu = mub - mua

   !PRINT*, 'dmu: ', dmu

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

     !PRINT*, 'ITERATION: ', ITERATION

      ! --- Compute midpoint, muc ---

      dmu = Half * dmu

     !PRINT*, 'dmu: ', dmu

      ! --- Bisection ---

      muc = mua + dmu

      ! --- Compute f(muc) for midpoint muc ---

     !PRINT*, 'Computing function for midpoint: ', muc

      CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, muc, fc )

     !PRINT*, 'fc: ', fc

      ! --- Change muc to mua or mub, depending on sign of fc ---

      IF( fa * fc .LT. Zero )THEN

        mub = muc
        fb = fc

      ELSE IF( fa * fc .GT. Zero )THEN

        mua = muc
        fa = fc

      ELSE

        CONVERGED = .TRUE.

      END IF

      IF( ABS( dmu ) / MAX( ABS( muc ), SqrtTiny ) .LE. dmu_min ) &
        CONVERGED = .TRUE.

    END DO

    mu = muc

   !PRINT*, 'mu, dmu: ', mu, dmu

  END SUBROUTINE SolveMu_Bisection 

  ! --- Algorithm 4.1 of Alefeld and Porta (1995) ---

!  SUBROUTINE SolveMu_TOM748 &
!               ( D_bar, Ne_bar, q, ru1, ru2, ru3, r, &
!                 bu1, bu2, bu3, bSq, rb, h0, v0, mu )
!
!    REAL(DP), INTENT(in) :: D_bar, Ne_bar, q
!    REAL(DP), INTENT(in) :: ru1, ru2, ru3
!    REAL(DP), INTENT(in) :: bu1, bu2, bu3
!    REAL(DP), INTENT(in) :: r, bSq, rb
!    REAL(DP), INTENT(in) :: h0, v0
!
!    REAL(DP), INTENT(out) :: mu
!
!    REAL(DP) :: mua, mub
!
!    mu = Zero
!
!    IF( r < h0 ) THEN
!
!      mua = Zero
!      mub = One / h0
!
!    ELSE
!
!      mua = Zero
!
!      CALL SolveMuBound_NewtonRaphson( r, bSq, rb, h0, mub )
!
!    END IF
!
!  END SUBROUTINE SolveMu_TOM748

END MODULE MHD_UtilitiesModule_Relativistic
