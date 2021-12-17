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
    nDOFX, &
    nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
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
    iAM_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputeAuxiliary_Fluid, &
    ComputePressure, &
    ComputePressureFromPrimitive, &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSpecificInternalEnergy
  USE EquationOfStateModule_TABLE, ONLY: &
    MinD, &
    MaxD, &
    MinT, &
    MaxT, &
    MinY, &
    MaxY
  USE UnitsModule, ONLY: &
    AtomicMassUnit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_MHD_Relativistic
  PUBLIC :: ComputeConserved_MHD_Relativistic
  PUBLIC :: ComputeFromConserved_MHD_Relativistic
  PUBLIC :: ComputeTimeStep_MHD_Relativistic
  PUBLIC :: Eigenvalues_MHD_Relativistic
  PUBLIC :: Flux_X1_MHD_Relativistic
  PUBLIC :: NumericalFlux_HLL_MHD_Relativistic

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
      iErr )

    REAL(DP), INTENT(in)    :: &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi
    REAL(DP), INTENT(in)    :: &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3
    REAL(DP), INTENT(out)   :: &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi
    INTEGER,  INTENT(inout), OPTIONAL :: &
      iErr

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

    B1 = CM_B1
    B2 = CM_B2
    B3 = CM_B3

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

   !*, 'D_bar: ', D

    Ne     = CM_Ne
    Ne_bar = Ne

   !*, 'Ne_bar: ', Ne_bar

    tau = CM_E

   !*, 'tau: ', tau

    Sd1 = CM_S1
    Sd2 = CM_S2
    Sd3 = CM_S3

   !*, 'Sd1: ', Sd1
   !*, 'Sd2: ', Sd2
   !*, 'Sd3: ', Sd3

    ! --- Eqs. 22-24 ---

    q = tau / D_bar

    !*, 'q: ', q

    rd1 = Sd1 / D_bar
    rd2 = Sd2 / D_bar
    rd3 = Sd3 / D_bar

    !*, 'rd1: ', rd1
    !*, 'rd2: ', rd2
    !*, 'rd3: ', rd3

    ru1 = rd1 / GF_Gm11 
    ru2 = rd2 / GF_Gm22
    ru3 = rd3 / GF_Gm33

    !*, 'ru1: ', ru1
    !*, 'ru2: ', ru2
    !*, 'ru3: ', ru3

    r = SQRT( ru1 * rd1 + ru2 * rd2 + ru3 * rd3 )

    !*, 'r: ', r

    bu1 = B1 / SQRT( D_bar )
    bu2 = B2 / SQRT( D_bar )
    bu3 = B3 / SQRT( D_bar )

    bd1 = GF_Gm11 * bu1
    bd2 = GF_Gm22 * bu2
    bd3 = GF_Gm33 * bu3

    ! --- Squared magnitude of b ---

    bSq = ( bu1 * bd1 + bu2 * bd2 + bu3 * bd3 )

   !*, 'bSq: ', bSq

    ! --- Contraction of r with b ---

    rb = ( bu1 * rd1 + bu2 * rd2 + bu3 * rd3 )

   !*, 'rb: ', rb

#ifdef MICROPHYSICS_WEAKLIB

    !*, 'MinD: ', MinD
    !*, 'MinT: ', MinT
    !*, 'MinY: ', MinY

    CALL ComputePressure &
           ( MinD, MinT, MinY, MinP )

    !*, 'MinP: ', MinP

    CALL ComputeSpecificInternalEnergy &
           ( MinD, MinT, MinY, MinEps )

    !*, 'MinEps: ', MinEps

    ! --- Not valid. Need better method of getting         ---
    ! --- minimum enthalpy (and deal with ideal gas case). ---

    h0 = One + MinEps + MinP / MinD  

    !*, 'h0: ', h0

#else

    h0 = One

#endif

    ! --- Upper velocity limit (Eqs. 32 and 33) --- 

    z0 = r / h0

   !*, 'z0: ', z0

    v0 = z0 / SQRT( One + z0**2 )

   !*, 'v0: ', v0

    CALL SolveMu_Bisection &
           ( D_bar, Ne_bar, q, ru1, ru2, ru3, r, &
             bu1, bu2, bu3, bSq, rb, h0, v0, mu )

!    CALL SolveMu_TOM748 &
!           ( D_bar, Ne_bar, q, ru1, ru2, ru3, r, &
!             bu1, bu2, bu3, bSq, rb, h0, v0, mu )

   !*, 'mu: ', mu

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

    ! --- Eqs. 38 and 39                              ---
    ! --- Note: Eq. 39 rewritten using Eqs. 25 and 30 ---
    ! --- to avoid B = 0 degeneracy.                  ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

    q_bar = q - Half * bSq - Half * mu**2 &
                             * ( bSq * r_barSq - rb**2 )

    ! --- Eq. 68 ---

    PM_V1 = mu * x * ( ru1 + mu * rb * bu1 )
    PM_V2 = mu * x * ( ru2 + mu * rb * bu2 )
    PM_V3 = mu * x * ( ru3 + mu * rb * bu3 )

    !PRINT*, 'PM_V1: ', PM_V1
 
    VSq = PM_V1**2 * GF_Gm11 + PM_V2**2 * GF_Gm22 + PM_V3**2 * GF_Gm33

    ! --- Eq. 40 ---

    W = One / SQRT( One - VSq )

    ! --- Eq. 41 ---

    PM_D = D_bar / W

    !PRINT*, 'PM_D: ', PM_D

    PM_Ne = CM_Ne / W

    ! --- Eq. 42 ---

    eps = W * ( q_bar - mu * r_barSq ) &
          + VSq * ( W**2 / ( One + W ) )

    PM_E = PM_D * eps

    !PRINT*, 'PM_E: ', PM_E

    B0u = ( W / GF_alpha ) &
            * ( GF_Gm11 * PM_V1 * CM_B1 &
                + GF_Gm22 * PM_V2 * CM_B2 &
                + GF_Gm33 * PM_V3 * CM_B3 )
 
   !PRINT*, 'PM_D: ',  PM_D
   !PRINT*, 'PM_V1: ', PM_V1
   !PRINT*, 'PM_V2: ', PM_V2
   !PRINT*, 'PM_V3: ', PM_V3
   !PRINT*, 'PM_E: ',  PM_E
   !PRINT*, 'PM_Ne: ', PM_Ne
                      
   !PRINT* , 'B0u: ', B0u
                            
    PM_B1 = ( CM_B1 / W ) + GF_Alpha * B0u &
                            * ( PM_V1 - ( GF_Beta1 / GF_Alpha  ) ) 
    PM_B2 = ( CM_B2 / W ) + GF_Alpha * B0u &
                            * ( PM_V2 - ( GF_Beta2 / GF_Alpha  ) )
    PM_B3 = ( CM_B3 / W ) + GF_Alpha * B0u &
                            * ( PM_V3 - ( GF_Beta3 / GF_Alpha  ) )

   !PRINT*, 'PM_B1: ', PM_B1
   !PRINT*, 'PM_B2: ', PM_B2
   !PRINT*, 'PM_B3: ', PM_B3
   !PRINT*

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
    ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi, &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi, &
      GF_Gm11, GF_Gm22, GF_Gm33, &
      GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3, &
      iErr )

    REAL(DP), INTENT(in)    :: &
      CM_D(:), CM_S1(:), CM_S2(:), CM_S3(:), CM_E(:), CM_Ne(:), &
      CM_B1(:), CM_B2(:), CM_B3(:), CM_Chi(:)
    REAL(DP), INTENT(in)    :: &
      GF_Gm11(:), GF_Gm22(:), GF_Gm33(:), &
      GF_Alpha(:), GF_Beta1(:), GF_Beta2(:), GF_Beta3(:)
    REAL(DP), INTENT(out)   :: &
      PM_D(:), PM_V1(:), PM_V2(:), PM_V3(:), PM_E(:), PM_Ne(:), &
      PM_B1(:), PM_B2(:), PM_B3(:), PM_Chi(:)
    INTEGER,  INTENT(inout), OPTIONAL :: &
      iErr(:)

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
               iErr   (iNX) )

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
      AM_P )

    REAL(DP), INTENT(in)  :: PM_D, PM_V1, PM_V2, PM_V3, &
                             PM_E, PM_Ne, PM_B1, PM_B2, &
                             PM_B3, PM_Chi, AM_P, &
                             GF_Gm11, GF_Gm22, GF_Gm33, &
                             GF_Alpha, GF_Beta1, GF_Beta2, GF_Beta3
    REAL(DP), INTENT(out) :: CM_D, CM_S1, CM_S2, CM_S3, &
                             CM_E, CM_Ne, CM_B1, CM_B2, &
                             CM_B3, CM_Chi

    REAL(DP) :: VSq, W, B0u, B0d, BSq, hStar, pStar 

    VSq = GF_Gm11 * PM_V1**2 + GF_Gm22 * PM_V2**2 + GF_Gm33 * PM_V3**2
    W = One / SQRT( One - VSq )

    B0u = ( GF_Gm11 * PM_V1 * PM_B1 &
             + GF_Gm22 * PM_V2 * PM_B2 &
             + GF_Gm33 * PM_V3 * PM_B3 ) &
         / ( GF_Alpha - GF_Gm11 * PM_V1 * GF_Beta1 &
             - GF_Gm22 * PM_V2 * GF_Beta2 &
             - GF_Gm33 * PM_V3 * GF_Beta3 )

    B0d =  - ( GF_Alpha / W ) &
             * ( GF_Gm11 * PM_V1 * PM_B1 &
                 + GF_Gm22 * PM_V2 * PM_B2 &
                 + GF_Gm33 * PM_V3 * PM_B3 )

   !PRINT*
   !PRINT*, 'Computing conserved variables.'

   !PRINT*, 'PM_D: ', PM_D
   !PRINT*, 'PM_V1: ', PM_V1
   !PRINT*, 'PM_V2: ', PM_V2
   !PRINT*, 'PM_V3: ', PM_V3
   !PRINT*, 'PM_E: ', PM_E
   !PRINT*, 'PM_Ne: ', PM_Ne
   !PRINT*, 'PM_B1: ', PM_B1
   !PRINT*, 'PM_B2: ', PM_B2
   !PRINT*, 'PM_B3: ', PM_B3
 
   !PRINT*, 'B0u: ', B0u

    BSq = B0d * B0u &
          + B0u * ( GF_Gm11 * GF_Beta1 * PM_B1 &
                    + GF_Gm22 * GF_Beta2 * PM_B2 &
                    + GF_Gm33 * GF_Beta3 * PM_B3 ) &
          + ( GF_Gm11 * PM_B1**2 &
              + GF_Gm22 * PM_B2**2 &
              + GF_Gm33 * PM_B3**2 )

    hStar = One + ( PM_E + AM_P ) / PM_D + BSq / PM_D
    pStar = AM_P + BSq / 2.0_DP

    CM_D   = W * PM_D
    CM_S1  = hStar * W**2 * PM_D * GF_Gm11 * PM_V1 &
             - GF_Alpha * B0u * ( GF_Gm11 * PM_B1 )
    CM_S2  = hStar * W**2 * PM_D * GF_Gm22 * PM_V2 &
             - GF_Alpha * B0u * ( GF_Gm22 * PM_B2 )
    CM_S3  = hStar * W**2 * PM_D * GF_Gm33 * PM_V3 &
             - GF_Alpha * B0u * ( GF_Gm33 * PM_B3 ) 
    CM_E   = hStar * W**2 * PM_D - pStar - ( GF_Alpha * B0u )**2 - W * PM_D
    CM_Ne  = W * PM_Ne
    CM_B1  = -W * GF_Alpha * B0u * ( PM_V1 - ( GF_Beta1 / GF_Alpha  ) ) &
             + W * PM_B1 
    CM_B2  = -W * GF_Alpha * B0u * ( PM_V2 - ( GF_Beta2 / GF_Alpha  ) ) &
             + W * PM_B2
    CM_B3  = -W * GF_Alpha * B0u * ( PM_V3 - ( GF_Beta3 / GF_Alpha  ) ) &
             + W * PM_B3
    CM_Chi = PM_Chi

   !PRINT*, 'CM_D: ', CM_D
   !PRINT*, 'CM_S1: ', CM_S1
   !PRINT*, 'CM_S2: ', CM_S2
   !PRINT*, 'CM_S3: ', CM_S3
   !PRINT*, 'CM_E: ', CM_E
   !PRINT*, 'CM_Ne: ', CM_Ne
   !PRINT*, 'CM_E: ', CM_E
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
      AM_P )

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
               AM_P (iNX) )

    END DO

  END SUBROUTINE ComputeConserved_Vector


  !> Compute primitive variables, pressure, and sound-speed from conserved
  !> variables for a data block.
  SUBROUTINE ComputeFromConserved_MHD_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iAM
    INTEGER :: iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

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

     !PRINT*, 'In cell: ', iX1

      iErr(iNX,iX1,iX2,iX3) = 0

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
               G   (iNX,iX1,iX2,iX3,iGF_Beta_3),    &
               iErr(iNX,iX1,iX2,iX3) )

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
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep )

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
    REAL(DP) :: P(nPM), VSq, W, Pressure, h, Cs, &
                B0u, B0d, BSq, Ca, aSq, EigVals(2)
    INTEGER  :: iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                             iX_B0(2):iX_E0(2), &
                             iX_B0(3):iX_E0(3))

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    TimeStep = HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iErr(iNX,iX1,iX2,iX3) = 0

      dX(1) = dX1(iX1)
      dX(2) = dX2(iX2)
      dX(3) = dX3(iX3)

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
               iErr(iNX,iX1,iX2,iX3) )

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
                P(iPM_D), P(iPM_V1), P(iPM_V2), P(iPM_V3), &
                P(iPM_E), P(iPM_Ne), &
                P(iPM_B1), P(iPM_B2), P(iPM_B3), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                G(iNX,iX1,iX2,iX3,iGF_Alpha), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_1), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_2), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_3) )

        PRINT*, 'The eigenvalues are: ', EigVals

        dt = dX(iDimX) / MAX( SqrtTiny, MAXVAL( ABS( EigVals ) ) )

        TimeStep = MIN( TimeStep, dt )

      END DO

    END DO
    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

    PRINT*, 'The timestep is: ', TimeStep

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_MHD_Relativistic


  FUNCTION Eigenvalues_MHD_Relativistic &
    ( Vi, Cs, Gmii, D, V1, V2, V3, E, Ne, &
      B1, B2, B3, Gm11, Gm22, Gm33, &
      Lapse, Shift1, Shift2, Shift3 )

    REAL(DP), INTENT(in) :: Vi, Cs, Gmii, D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, P, h, B0u, B0d, BSq, &
                Ca, aSq, Eigenvalues_MHD_Relativistic(2)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = One / SQRT( One - VSq )

    !PRINT*, 'Lorentz factor: ', W

    CALL ComputePressureFromPrimitive &
           ( D, E, Ne, P )
 
    !PRINT*, 'Pressure: ', Pressure

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

    ! Estimate of max/min fast magnetosonic 
    ! eigenvalues from Del Zanna et al. (2007)

    Eigenvalues_MHD_Relativistic(1) &
      = ( ( One - aSq ) * V1 &
          + SQRT( aSq * ( One - VSq ) &
                  * ( ( One - VSq * aSq ) * ( One / Gm11 ) &
                  - ( One - aSq ) * V1**2 ) ) ) &
        / ( One - VSq * aSq )

    !Eigenvalues_MHD_Relativistic(1) = 0.90_DP + 0.005_DP
    !  = Lapse * Eigenvalues_MHD_Relativistic(1) - Shift1

    Eigenvalues_MHD_Relativistic(2) &
       = ( ( One - aSq ) * V1 &
          - SQRT( aSq * ( One - VSq ) & 
                  * ( ( One - VSq * aSq ) * ( One / Gm11 ) &
                  - ( One - aSq ) * V1**2 ) ) ) &
        / ( One - VSq * aSq )

    !Eigenvalues_MHD_Relativistic(2) = 0.90_DP - 0.005_DP
    !  = Lapse * Eigenvalues_MHD_Relativistic(2) - Shift1

    RETURN
  END FUNCTION Eigenvalues_MHD_Relativistic


  !> Compute the physical flux in the X1-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X1_MHD_Relativistic &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, &
      P, Gm11, Gm22, Gm33, Lapse, & 
      Shift1, Shift2, Shift3 )

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, P, &
                            Gm11, Gm22, Gm33, Lapse, &
                            Shift1, Shift2, Shift3

    REAL(DP) :: VSq, W, Pstar, h, hStar, B0u, B0d, BSq, Flux_X1_MHD_Relativistic(nCM)

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

    Flux_X1_MHD_Relativistic(iCM_S1) &
      = D * hStar * W**2 * Gm11 * V1  * ( V1 - Shift1 / Lapse ) + pStar &
        - Gm11 * ( B1**2 + B0u * B1 * Shift1 &
                   + ( Lapse * B0u )**2 * Shift1**2 )

    Flux_X1_MHD_Relativistic(iCM_S2) &
      = D * hStar * W**2 * Gm22 * V2  * ( V1 - Shift1 / Lapse ) &
        - Gm22 * ( B1 * B2 + B0u * B1 * Shift2 &
                   + ( Lapse * B0u )**2 * Shift1 * Shift2 )

    Flux_X1_MHD_Relativistic(iCM_S3) &
      = D * hStar * W**2 * Gm33 * V3  * ( V1 - Shift1 / Lapse ) &
        - Gm33 * ( B1 * B3 + B0u * B1 * Shift3 &
                   + ( Lapse * B0u )**2 * Shift1 * Shift3 )

    Flux_X1_MHD_Relativistic(iCM_E)  &
      = D * W * ( hStar * W - One ) * ( V1 - Shift1 / Lapse ) &
        + Shift1 / Lapse * pStar &
        - Lapse * B0u * B1 + Two * Lapse * (B0u)**2 * Shift1

    Flux_X1_MHD_Relativistic(iCM_Ne) &
      = Ne * W * ( V1 - Shift1 / Lapse )

    Flux_X1_MHD_Relativistic(iCM_B1) &
      = 0.0_DP

    Flux_X1_MHD_Relativistic(iCM_B2) &
      = W * ( V1 - Shift1 / Lapse ) * B2 - W * ( V2 - Shift2 / Lapse ) * B1

    Flux_X1_MHD_Relativistic(iCM_B3) &
      = W * ( V1 - Shift1 / Lapse ) * B3 - W * ( V3 - Shift3 / Lapse ) * B1

    Flux_X1_MHD_Relativistic(iCM_Chi) &
      = 0.0_DP

    RETURN
  END FUNCTION Flux_X1_MHD_Relativistic


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


    !*, 'mu: ', mu

    ! --- Eq. 26 ---

    x = One / ( One + mu * bSq )

    ! --- Eqs. 38 and 39                              ---
    ! --- Note: Eq. 39 rewritten using Eqs. 25 and 30 ---
    ! --- to avoid B = 0 degeneracy.                  ---

    r_barSq = r**2 * x**2 + mu * x * ( One + x ) * rb**2

    !*, 'r_barSq: ', r_barSq

    !*, 'SQRT( r_barSq ): ', SQRT( r_barSq )

    q_bar = q - Half * bSq - Half * mu**2 &
                             * ( bSq * r_barSq - rb**2 )

    !*, 'q_bar: ', q_bar

    ! --- Eq. 40 ---

    !*, 'v0: ', v0

    v = MIN( mu * SQRT( r_barSq ), v0 )

    !*, 'v: ', v

    W = One / SQRT( One - v**2 )

    !*, 'W: ', W

    ! --- Eq. 41 ---

    rho = D_bar / W

    rho_e = Ne_bar / W

    ! --- Eq. 42 ---

    eps = W * ( q_bar - mu * r_barSq ) &
          + v**2 * ( W**2 / ( 1.0_DP + W ) )

    ! --- Eq. 43 ---

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rho, eps, AtomicMassUnit * rho_e / rho, P )

    a = P / ( rho * ( One + eps ) )

    ! --- Eqs. 46-48 ---

    nua = ( One + a ) * ( ( One + eps ) / W )

    nub = ( One + a ) * ( One + q_bar - mu * r_barSq )

    nu = MAX( nua, nub )

    ! --- Eq. 44 --- 

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
    REAL(DP), PARAMETER :: Tolmu = 1.0e-16_DP
    REAL(DP), PARAMETER :: Tolf = 1.0e-16_DP
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolmu ) / LOG( Two ) )

    !*, 'Solving for the upper bound with the Newton-Raphson method.'

    muO = One / h0

    !*, 'Initial guess: ', One / h0

    !*, 'MAX_ITER: ', MAX_ITER

    CONVERGED = .FALSE.
    ITER = 0
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      !*
      !*, 'ITER: ', ITER
      !*

      !*, 'muO: ', muO

      !*, 'Computing the auxiliary function for mu = : ', muO

      CALL ComputeAuxiliaryFunMu( r, bSq, rb, h0, muO, f )

      !*, 'f: ', f

      !*, 'Computing the auxiliary function derivative for mu = ', muO

      CALL ComputeAuxiliaryFunDerivativeMu( r, bSq, rb, h0, muO, df )

      !*, 'df: ', df

      muN = muO - f / df 

      !*, 'muN: ', muN

      dmu = muN - muO

      !*, 'dmu: ', dmu

      IF( ABS ( dmu / muO ) .LT. Tolmu ) CONVERGED = .TRUE.

      muO = muN

    END DO

    mu = muO  

  END SUBROUTINE

  
  SUBROUTINE SolveMu_Bisection &
               ( D_bar, Ne_bar, q, ru1, ru2, ru3, r, &
                 bu1, bu2, bu3, bSq, rb, h0, v0, mu )

    REAL(DP), INTENT(in) :: D_bar, Ne_bar, q
    REAL(DP), INTENT(in) :: ru1, ru2, ru3
    REAL(DP), INTENT(in) :: bu1, bu2, bu3
    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0, v0

    REAL(DP), INTENT(out) :: mu

    LOGICAL             :: CONVERGED
    INTEGER             :: ITERATION
    REAL(DP)            :: mua, mub, muc, dmu
    REAL(DP)            :: fa, fb, fc
    REAL(DP), PARAMETER :: dmu_min = 1.0e-16_DP
    INTEGER,  PARAMETER :: MAX_IT = 4 - INT( LOG( dmu_min) / LOG( Two ) )

    IF( r < h0 ) THEN

      !*, 'r < h0'

      mua = Zero
      mub = One / h0

    ELSE

      !*, 'r >= h0'

      mua = Zero

      CALL SolveMuBound_NewtonRaphson( r, bSq, rb, h0, mub )

    END IF

   !*

   !*, 'Computing function for lower bound: ', mua 

   !*, '------------------------------------------------------------'

    CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, mua, fa )

   !*, 'fa: ', fa

   !*

   !*, 'Computing function for upper bound: ', mub

   !*, '------------------------------------------------------------'    

    CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, mub, fb )

   !*, 'fb: ', fb

    ! --- Check that sign of FunZ changes across bounds ---

    !IF( .NOT. fa * fb .LT. 0 ) &
      !PRINT*, 'Cannot perform bisection!'

    dmu = mub - mua

   !*, 'dmu: ', dmu

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

      !*, 'ITERATION: ', ITERATION

      ! --- Compute midpoint, muc ---

      dmu = Half * dmu

      ! --- Bisection ---

      muc = mua + dmu

      ! --- Compute f(muc) for midpoint muc ---

      !*, 'Computing function for midpoint: ', muc

      CALL ComputeFunMu( D_bar, Ne_bar, q, r, bSq, rb, v0, muc, fc )

      !*, 'fc: ', fc

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

  SUBROUTINE SolveMu_TOM748 &
               ( D_bar, Ne_bar, q, ru1, ru2, ru3, r, &
                 bu1, bu2, bu3, bSq, rb, h0, v0, mu )

    REAL(DP), INTENT(in) :: D_bar, Ne_bar, q
    REAL(DP), INTENT(in) :: ru1, ru2, ru3
    REAL(DP), INTENT(in) :: bu1, bu2, bu3
    REAL(DP), INTENT(in) :: r, bSq, rb
    REAL(DP), INTENT(in) :: h0, v0

    REAL(DP), INTENT(out) :: mu

    REAL(DP) :: mua, mub

    mu = Zero

    IF( r < h0 ) THEN

      mua = Zero
      mub = One / h0

    ELSE

      mua = Zero

      CALL SolveMuBound_NewtonRaphson( r, bSq, rb, h0, mub )

    END IF

  END SUBROUTINE SolveMu_TOM748

END MODULE MHD_UtilitiesModule_Relativistic
