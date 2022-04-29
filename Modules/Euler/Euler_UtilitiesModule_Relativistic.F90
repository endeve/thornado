!> Perform computations related to the 3+1, CFA Euler equations.
!> Find the equations in Rezzolla & Zanotti, Relativistic Hydrodynamics, 2013,
!> Equation 7.234.
MODULE Euler_UtilitiesModule_Relativistic

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
    iGF_Beta_1
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nAF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputeAuxiliary_Fluid, &
    ComputePressureFromSpecificInternalEnergy
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    ComputeSpecificInternalEnergy_TABLE
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE TimersModule_Euler,ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_ComputePrimitive, &
    Timer_Euler_CP_CopyIn, &
    Timer_Euler_CP_CopyOut, &
    Timer_Euler_CP_Permute, &
    Timer_Euler_CP_GetBounds, &
    Timer_Euler_CP_Bisection, &
    Timer_Euler_CP_RecoverPrimitives, &
    Timer_Euler_ComputeTimeStep, &
    Timer_Euler_CTS_ComputeTimeStep, &
    Timer_Euler_CTS_CopyIn, &
    Timer_Euler_CTS_CopyOut, &
    Timer_Euler_ComputeFromConserved, &
    Timer_Euler_CFC_CopyIn, &
    Timer_Euler_CFC_CopyOut, &
    Timer_Euler_CFC_ComputePrimitive, &
    Timer_Euler_CFC_ErrorCheck
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_Euler_Relativistic
  PUBLIC :: ComputeConserved_Euler_Relativistic
  PUBLIC :: ComputeFromConserved_Euler_Relativistic
  PUBLIC :: ComputeTimeStep_Euler_Relativistic
  PUBLIC :: Eigenvalues_Euler_Relativistic
  PUBLIC :: AlphaMiddle_Euler_Relativistic
  PUBLIC :: Flux_X1_Euler_Relativistic
  PUBLIC :: Flux_X2_Euler_Relativistic
  PUBLIC :: Flux_X3_Euler_Relativistic
  PUBLIC :: StressTensor_Diagonal_Euler_Relativistic
  PUBLIC :: NumericalFlux_LLF_Euler_Relativistic
  PUBLIC :: NumericalFlux_HLL_Euler_Relativistic
  PUBLIC :: NumericalFlux_X1_HLLC_Euler_Relativistic
  PUBLIC :: NumericalFlux_X2_HLLC_Euler_Relativistic
  PUBLIC :: NumericalFlux_X3_HLLC_Euler_Relativistic

  INTERFACE ComputePrimitive_Euler_Relativistic
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector_old
  END INTERFACE ComputePrimitive_Euler_Relativistic

  INTERFACE ComputeConserved_Euler_Relativistic
    MODULE PROCEDURE ComputeConserved_Scalar
    MODULE PROCEDURE ComputeConserved_Vector
  END INTERFACE ComputeConserved_Euler_Relativistic

  REAL(DP), PARAMETER :: Offset_Temperature = 1.0e-12_DP
  REAL(DP), PARAMETER :: Offset_Epsilon     = 1.0e-12_DP


CONTAINS


  SUBROUTINE ComputePrimitive_Vector_old &
    ( uD, uS1, uS2, uS3, uE, uNe, &
      pD, pV1, pV2, pV3, pE, pNe, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: uD(:), uS1(:), uS2(:), uS3(:), uE(:), uNe(:)
    REAL(DP), INTENT(out) :: pD(:), pV1(:), pV2(:), pV3(:), pE(:), pNe(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

    INTEGER :: N, ErrorExists
    INTEGER :: iNX, iErr(SIZE(uD))

    N = SIZE(uD)

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_CP_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    uD, uS1, uS2, uS3, uE, uNe, &
    !$OMP             Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    !$OMP MAP( alloc: pD, pV1, pV2, pV3, pE, pNe, iErr )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC             Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    !$ACC CREATE(     pD, pV1, pV2, pV3, pE, pNe, iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_CP_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_CP_RecoverPrimitives )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP REDUCTION( +:ErrorExists )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC REDUCTION( +:ErrorExists ) &
    !$ACC PRESENT( uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC          pD, pV1, pV2, pV3, pE, pNe, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP REDUCTION( +:ErrorExists )
#endif
    DO iNX = 1, N

      iErr(iNX) = 0

      CALL ComputePrimitive_Scalar &
             ( uD(iNX), uS1(iNX), uS2(iNX), uS3(iNX), uE(iNX), uNe(iNX), &
               pD(iNX), pV1(iNX), pV2(iNX), pV3(iNX), pE(iNX), pNe(iNX), &
               Gm_dd_11(iNX), Gm_dd_22(iNX), Gm_dd_33(iNX), iErr(iNX) )

      ErrorExists = ErrorExists + iErr(iNX)

    END DO

    CALL TimersStop_Euler( Timer_Euler_CP_RecoverPrimitives )

    CALL TimersStart_Euler( Timer_Euler_CP_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    pD, pV1, pV2, pV3, pE, pNe, iErr ) &
    !$OMP MAP( release: uD, uS1, uS2, uS3, uE, uNe, &
    !$OMP               Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      pD, pV1, pV2, pV3, pE, pNe, iErr ) &
    !$ACC DELETE(       uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC               Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CP_CopyOut )

    IF( ErrorExists .GT. 0 )THEN

      DO iNX = 1, N

        CALL DescribeError_Euler &
          ( iErr(iNX), &
            Int_Option = [ iNX ], &
            Real_Option = [ uD(iNX), uS1(iNX), uS2(iNX), uS3(iNX), &
                            uE(iNX), uNe(iNX), &
                            Gm_dd_11(iNX), Gm_dd_22(iNX), Gm_dd_33(iNX) ] )

      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

  END SUBROUTINE ComputePrimitive_Vector_old


  SUBROUTINE ComputePrimitive_Vector_new &
    ( uD, uS1, uS2, uS3, uE, uNe, &
      pD, pV1, pV2, pV3, pE, pNe, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: uD(:), uS1(:), uS2(:), uS3(:), uE(:), uNe(:)
    REAL(DP), INTENT(out) :: pD(:), pV1(:), pV2(:), pV3(:), pE(:), pNe(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)

    INTEGER  :: N, iX, ITERATION, ErrorExists
    INTEGER  :: IsPhysical(SIZE(uD)), iErr(SIZE(uD))

    REAL(DP) :: W, eps, DhW, p
    REAL(DP) :: q (SIZE(uD)), r (SIZE(uD)), k (SIZE(uD)), &
                dz(SIZE(uD)), &
                za(SIZE(uD)), zb(SIZE(uD)), zc(SIZE(uD)), &
                fa(SIZE(uD)), fb(SIZE(uD)), fc(SIZE(uD))

    LOGICAL  :: ITERATE(SIZE(uD))

    CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_CP_CopyIn )

    ITERATE = .TRUE.

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    uD, uS1, uS2, uS3, uE, uNe, &
    !$OMP             Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP             ITERATE ) &
    !$OMP MAP( alloc: pD, pV1, pV2, pV3, pE, pNe, &
    !$OMP             q, r, k, &
    !$OMP             dz, &
    !$OMP             za, zb, zc, &
    !$OMP             fa, fb, fc, &
    !$OMP             iErr )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC             Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC             ITERATE ) &
    !$ACC CREATE(     pD, pV1, pV2, pV3, pE, pNe, &
    !$ACC             q, r, k, &
    !$ACC             dz, &
    !$ACC             za, zb, zc, &
    !$ACC             fa, fb, fc, &
    !$ACC             iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_CP_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_CP_GetBounds )

    N = SIZE( uD )

    ! --- Eqs. C2/C23/C25 ---

    CALL GetBoundsForBisection &
           ( N, uD, uS1, uS2, uS3, uE, uNe, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             dz, za, zb, fa, fb, q, r, k, iErr )

    CALL TimersStop_Euler( Timer_Euler_CP_GetBounds )

    CALL TimersStart_Euler( Timer_Euler_CP_Bisection )

    CALL SolveZ_Bisection_Vector &
           ( N, uD, uNe, q, r, k, &
             za, zb, fa, fb, &
             dz, zc, fc, ITERATE )

    CALL TimersStop_Euler( Timer_Euler_CP_Bisection )

    CALL TimersStart_Euler( Timer_Euler_CP_RecoverPrimitives )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( W, eps, DhW, p )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( W, eps, DhW, p ) &
    !$ACC PRESENT( zc, q, r, &
    !$ACC          uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC          pD, pV1, pV2, pV3, pE, pNe, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( W, eps, DhW, p )
#endif
    DO iX = 1, N

      IF( uD(iX) .GT. Min_D )THEN

        ! --- Eqs. C15/C16 ---

        W       = SQRT( One + zc(iX)**2 )
        pD (iX) = uD (iX) / W
        eps     = W * q(iX) - zc(iX) * r(iX) + zc(iX)**2 / ( One + W )

        pNe(iX) = uNe(iX) / W

        CALL ComputePressureFromSpecificInternalEnergy &
               ( pD(iX), eps, AtomicMassUnit * pNe(iX) / pD(iX), p )

        DhW = uD(iX) * ( One + eps + p / pD(iX) ) * W

        ! --- Eq. C26 ---

        pV1(iX) = ( uS1(iX) / Gm_dd_11(iX) ) / DhW
        pV2(iX) = ( uS2(iX) / Gm_dd_22(iX) ) / DhW
        pV3(iX) = ( uS3(iX) / Gm_dd_33(iX) ) / DhW

        pE (iX) = pD(iX) * eps

      END IF

    END DO

    CALL TimersStop_Euler( Timer_Euler_CP_RecoverPrimitives )

    CALL TimersStart_Euler( Timer_Euler_CP_Permute )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP REDUCTION( +:ErrorExists )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC REDUCTION( +:ErrorExists ) &
    !$ACC PRESENT( uD, uE, uNe, pD, pV1, pV2, pV3, pE, pNe, iErr ) &
    !$ACC COPYOUT(     ErrorExists )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP REDUCTION( +:ErrorExists )
#endif
    DO iX = 1, N

      ErrorExists = ErrorExists + iErr(iX)

      IF( uD(iX) .LT. Min_D )THEN

        pD (iX) = 1.01_DP * Min_D
        pV1(iX) = Zero
        pV2(iX) = Zero
        pV3(iX) = Zero
        pE (iX) = MAX( uE(iX), SqrtTiny )
        pNe(iX) = uNe(iX) / uD(iX)

      END IF

    END DO

    CALL TimersStop_Euler( Timer_Euler_CP_Permute )

    CALL TimersStart_Euler( Timer_Euler_CP_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr, &
    !$OMP               pD, pV1, pV2, pV3, pE, pNe ) &
    !$OMP MAP( release: uD, uS1, uS2, uS3, uE, uNe, &
    !$OMP               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP               q, r, k, &
    !$OMP               dz, &
    !$OMP               za, zb, zc, &
    !$OMP               fa, fb, fc, &
    !$OMP               ITERATE )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr, &
    !$ACC               pD, pV1, pV2, pV3, pE, pNe ) &
    !$ACC DELETE(       uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC               q, r, k, &
    !$ACC               dz, &
    !$ACC               za, zb, zc, &
    !$ACC               fa, fb, fc, &
    !$ACC               ITERATE )
#endif

    IF( ErrorExists .GT. 0 )THEN

      DO iX = 1, N

        CALL DescribeError_Euler &
          ( iErr(iX), &
            Int_Option = [ iX ], &
            Real_Option = [ uD(iX), uS1(iX), uS2(iX), uS3(iX), &
                            uE(iX), uNe(iX), &
                            Gm_dd_11(iX), Gm_dd_22(iX), Gm_dd_33(iX) ] )

      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_CP_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

  END SUBROUTINE ComputePrimitive_Vector_new


  !> Compute the primitive variables from the conserved variables,
  !> a la Galeazzi et al., (2013), Phys. Rev. D., 88, 064009, Appendix C
  !> @todo Modify for tabular EOS
  SUBROUTINE ComputePrimitive_Scalar &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm11, GF_Gm22, GF_Gm33, iErr )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne
    REAL(DP), INTENT(in)    :: &
      GF_Gm11, GF_Gm22, GF_Gm33
    REAL(DP), INTENT(out)   :: &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne
    INTEGER,  INTENT(inout), OPTIONAL :: &
      iErr

    REAL(DP) :: S, q, r, k, z0
    REAL(DP) :: W, eps, p, h, DhW

    S = SQRT( CF_S1**2 / GF_Gm11 + CF_S2**2 / GF_Gm22 + CF_S3**2 / GF_Gm33 )

    ! --- Eq. C2 ---

    q = CF_E / CF_D
    r = S    / CF_D
    k = r    / ( One + q )

    ! --- Ensure primitive fields can be recovered ---

    IF( CF_D .LT. Min_D )THEN

      PF_D  = 1.01_DP * Min_D
      PF_V1 = Zero
      PF_V2 = Zero
      PF_V3 = Zero
      PF_E  = MAX( CF_E, SqrtTiny ) ! What to put here
      PF_Ne = CF_Ne / CF_D

      RETURN

    END IF

    IF( q .LT. Zero )THEN

      r = k
      q = Zero

    END IF

    ! --- Solve for primitive ---

    CALL SolveZ_Bisection_Scalar( CF_D, CF_Ne, q, r, k, z0, iErr )

    ! --- Eq. C15 ---

    W     = SQRT( One + z0**2 )
    PF_D  = CF_D  / W
    PF_Ne = CF_Ne / W

    ! --- Eq. C16 ---

    eps = W * q - z0 * r + z0**2 / ( One + W )

    CALL ComputePressureFromSpecificInternalEnergy &
           ( PF_D, eps, AtomicMassUnit * PF_Ne / PF_D, p )

    h = One + eps + p / PF_D

    DhW = CF_D * h * W

    PF_V1 = ( CF_S1 / GF_Gm11 ) / DhW
    PF_V2 = ( CF_S2 / GF_Gm22 ) / DhW
    PF_V3 = ( CF_S3 / GF_Gm33 ) / DhW
    PF_E  = PF_D * eps

  END SUBROUTINE ComputePrimitive_Scalar


  !> Compute conserved variables from primitive variables.
  SUBROUTINE ComputeConserved_Scalar &
    ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      Gm11, Gm22, Gm33, &
      AF_P )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: PF_D, PF_V1, PF_V2, PF_V3, &
                             PF_E, PF_Ne, AF_P, &
                             Gm11, Gm22, Gm33
    REAL(DP), INTENT(out) :: CF_D, CF_S1, CF_S2, CF_S3, &
                             CF_E, CF_Ne

    REAL(DP) :: VSq, W, h

    VSq = Gm11 * PF_V1**2 + Gm22 * PF_V2**2 + Gm33 * PF_V3**2

    W = One / SQRT( One - VSq )
    h = One + ( PF_E + AF_P ) / PF_D

    CF_D  = W * PF_D
    CF_S1 = h * W**2 * PF_D * Gm11 * PF_V1
    CF_S2 = h * W**2 * PF_D * Gm22 * PF_V2
    CF_S3 = h * W**2 * PF_D * Gm33 * PF_V3
    CF_E  = h * W**2 * PF_D - AF_P - W * PF_D
    CF_Ne = W * PF_Ne

  END SUBROUTINE ComputeConserved_Scalar


  SUBROUTINE ComputeConserved_Vector &
    ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      Gm11, Gm22, Gm33, &
      AF_P )

    REAL(DP), INTENT(in)  :: PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), &
                             PF_E(:), PF_Ne(:), AF_P(:), &
                             Gm11(:), Gm22(:), Gm33(:)
    REAL(DP), INTENT(out) :: CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), &
                             CF_E(:), CF_Ne(:)

    INTEGER :: iNX

    DO iNX = 1, SIZE( PF_D )

      CALL ComputeConserved_Scalar &
             ( PF_D (iNX), &
               PF_V1(iNX), &
               PF_V2(iNX), &
               PF_V3(iNX), &
               PF_E (iNX), &
               PF_Ne(iNX), &
               CF_D (iNX), &
               CF_S1(iNX), &
               CF_S2(iNX), &
               CF_S3(iNX), &
               CF_E (iNX), &
               CF_Ne(iNX), &
               Gm11 (iNX), &
               Gm22 (iNX), &
               Gm33 (iNX), &
               AF_P (iNX) )

    END DO

  END SUBROUTINE ComputeConserved_Vector


  !> Compute primitive variables, pressure, and sound-speed from conserved
  !> variables for a data block.
  SUBROUTINE ComputeFromConserved_Euler_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iAF
    INTEGER :: iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    CALL TimersStart_Euler( Timer_Euler_ComputeFromConserved )

    ! --- Update primitive variables, pressure, and sound speed ---

    CALL TimersStart_Euler( Timer_Euler_CFC_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U ) &
    !$OMP MAP( alloc: P, A, iErr )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U ) &
    !$ACC CREATE(     P, A, iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_CFC_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_CFC_ComputePrimitive )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, A )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iAF = 1, nAF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      A(iNX,iX1,iX2,iX3,iAF) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G, U, P, A, iErr )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iErr(iNX,iX1,iX2,iX3) = 0

      CALL ComputePrimitive_Euler_Relativistic &
             ( U   (iNX,iX1,iX2,iX3,iCF_D ),        &
               U   (iNX,iX1,iX2,iX3,iCF_S1),        &
               U   (iNX,iX1,iX2,iX3,iCF_S2),        &
               U   (iNX,iX1,iX2,iX3,iCF_S3),        &
               U   (iNX,iX1,iX2,iX3,iCF_E ),        &
               U   (iNX,iX1,iX2,iX3,iCF_Ne),        &
               P   (iNX,iX1,iX2,iX3,iPF_D ),        &
               P   (iNX,iX1,iX2,iX3,iPF_V1),        &
               P   (iNX,iX1,iX2,iX3,iPF_V2),        &
               P   (iNX,iX1,iX2,iX3,iPF_V3),        &
               P   (iNX,iX1,iX2,iX3,iPF_E ),        &
               P   (iNX,iX1,iX2,iX3,iPF_Ne),        &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_33),  &
               iErr(iNX,iX1,iX2,iX3) )

      CALL ComputeAuxiliary_Fluid &
             ( P(iNX,iX1,iX2,iX3,iPF_D ), &
               P(iNX,iX1,iX2,iX3,iPF_E ), &
               P(iNX,iX1,iX2,iX3,iPF_Ne), &
               A(iNX,iX1,iX2,iX3,iAF_P ), &
               A(iNX,iX1,iX2,iX3,iAF_T ), &
               A(iNX,iX1,iX2,iX3,iAF_Ye), &
               A(iNX,iX1,iX2,iX3,iAF_S ), &
               A(iNX,iX1,iX2,iX3,iAF_E ), &
               A(iNX,iX1,iX2,iX3,iAF_Gm), &
               A(iNX,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_CFC_ComputePrimitive )

    CALL TimersStart_Euler( Timer_Euler_CFC_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    P, A, iErr ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, U )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      P, A, iErr ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, U )
#endif

    CALL TimersStop_Euler( Timer_Euler_CFC_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_CFC_ErrorCheck )

    IF( ANY( iErr .NE. 0 ) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'ERROR: ComputeFromConserved_Euler_Relativistic'

          WRITE(*,*) 'iX1, iX2, iX3 = ', iX1, iX2, iX3

          CALL DescribeError_Euler &
            ( iErr(iNX,iX1,iX2,iX3), &
              Int_Option = [ iNX ], &
              Real_Option = [ U(iNX,iX1,iX2,iX3,iCF_D ), &
                              U(iNX,iX1,iX2,iX3,iCF_S1), &
                              U(iNX,iX1,iX2,iX3,iCF_S2), &
                              U(iNX,iX1,iX2,iX3,iCF_S3), &
                              U(iNX,iX1,iX2,iX3,iCF_E ), &
                              U(iNX,iX1,iX2,iX3,iCF_Ne), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ] )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_CFC_ErrorCheck )

    CALL TimersStop_Euler( Timer_Euler_ComputeFromConserved )

  END SUBROUTINE ComputeFromConserved_Euler_Relativistic


  !> Loop over all the elements in the spatial domain and compute the minimum
  !> required time-step for numerical stability.
  SUBROUTINE ComputeTimeStep_Euler_Relativistic &
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
    REAL(DP) :: P(nPF), Cs, EigVals(nCF)
    INTEGER  :: iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                             iX_B0(2):iX_E0(2), &
                             iX_B0(3):iX_E0(3))

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    CALL TimersStart_Euler( Timer_Euler_ComputeTimeStep )

    CALL TimersStart_Euler( Timer_Euler_CTS_CopyIn )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    G, U, iX_B0, iX_E0, dX1, dX2, dX3 ) &
    !$OMP MAP( alloc: iErr )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     G, U, iX_B0, iX_E0, dX1, dX2, dX3 ) &
    !$ACC CREATE(     iErr )
#endif

    CALL TimersStop_Euler( Timer_Euler_CTS_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_CTS_ComputeTimeStep )

    TimeStep = HUGE( One )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( dX, dt, P, Cs, EigVals ) &
    !$OMP REDUCTION( MIN: TimeStep )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( dX, dt, P, Cs, EigVals ) &
    !$ACC REDUCTION( MIN: TimeStep ) &
    !$ACC PRESENT( G, U, iX_B0, iX_E0, dX1, dX2, dX3, iErr )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( dX, dt, P, Cs, EigVals ) &
    !$OMP REDUCTION( MIN: TimeStep )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iErr(iNX,iX1,iX2,iX3) = 0

      dX(1) = dX1(iX1)
      dX(2) = dX2(iX2)
      dX(3) = dX3(iX3)

      CALL ComputePrimitive_Euler_Relativistic &
             ( U   (iNX,iX1,iX2,iX3,iCF_D ), &
               U   (iNX,iX1,iX2,iX3,iCF_S1), &
               U   (iNX,iX1,iX2,iX3,iCF_S2), &
               U   (iNX,iX1,iX2,iX3,iCF_S3), &
               U   (iNX,iX1,iX2,iX3,iCF_E ), &
               U   (iNX,iX1,iX2,iX3,iCF_Ne), &
               P   (iPF_D ), P(iPF_V1), P(iPF_V2), &
               P   (iPF_V3), P(iPF_E ), P(iPF_Ne), &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G   (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               iErr(iNX,iX1,iX2,iX3) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(iPF_D), P(iPF_E), P(iPF_Ne), Cs )

      DO iDimX = 1, nDimsX

        EigVals &
          = Eigenvalues_Euler_Relativistic &
              ( P(iPF_V1+(iDimX-1)), Cs, &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11+(iDimX-1)), &
                P(iPF_V1), P(iPF_V2), P(iPF_V3), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                G(iNX,iX1,iX2,iX3,iGF_Alpha), &
                G(iNX,iX1,iX2,iX3,iGF_Beta_1+(iDimX-1)) )

        dt = dX(iDimX) / MAX( SqrtTiny, MAXVAL( ABS( EigVals ) ) )

        TimeStep = MIN( TimeStep, dt )

      END DO

    END DO
    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

    CALL TimersStop_Euler( Timer_Euler_CTS_ComputeTimeStep )

    CALL TimersStart_Euler( Timer_Euler_CTS_CopyOut )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    iErr ) &
    !$OMP MAP( release: G, U, iX_B0, iX_E0, dX1, dX2, dX3 )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      iErr ) &
    !$ACC DELETE(       G, U, iX_B0, iX_E0, dX1, dX2, dX3 )
#endif

    CALL TimersStop_Euler( Timer_Euler_CTS_CopyOut )

    END ASSOCIATE ! dX1, etc.

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR: ComputeTimeStep_Euler_Relativistic'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

          WRITE(*,*) 'iX1, iX2, iX3 = ', iX1, iX2, iX3

          CALL DescribeError_Euler &
            ( iErr(iNX,iX1,iX2,iX3), &
              Int_Option = [ iNX ], &
              Real_Option = [ U(iNX,iX1,iX2,iX3,iCF_D ), &
                              U(iNX,iX1,iX2,iX3,iCF_S1), &
                              U(iNX,iX1,iX2,iX3,iCF_S2), &
                              U(iNX,iX1,iX2,iX3,iCF_S3), &
                              U(iNX,iX1,iX2,iX3,iCF_E ), &
                              U(iNX,iX1,iX2,iX3,iCF_Ne), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ] )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_Euler_ComputeTimeStep )

  END SUBROUTINE ComputeTimeStep_Euler_Relativistic


  FUNCTION Eigenvalues_Euler_Relativistic &
    ( Vi, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: Vi, Cs, Gmii, V1, V2, V3, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, Eigenvalues_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    Eigenvalues_Euler_Relativistic(1) &
      = Lapse / ( One - VSq * Cs**2 ) * ( Vi * ( One - Cs**2 ) &
        - Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gmii &
           - Vi**2 * ( One - Cs**2 ) ) ) ) - Shift

    Eigenvalues_Euler_Relativistic(2) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(3) &
      = Lapse / ( One - VSq * Cs**2 ) * ( Vi * ( One - Cs**2 ) &
        + Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gmii &
           - Vi**2 * ( One - Cs**2 ) ) ) ) - Shift

    Eigenvalues_Euler_Relativistic(4) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(5) &
      = Lapse * Vi - Shift

    Eigenvalues_Euler_Relativistic(6) &
      = Lapse * Vi - Shift

    RETURN
  END FUNCTION Eigenvalues_Euler_Relativistic


  !> Estimate the contact wave-speed as suggested by
  !> Mignone & Bodo, (2005), MNRAS, 364, 126.
  !> @param Shift The ith contravariant component of the shift-vector.
  !> @param Gmii The ith covariant component of the spatial three-metric.
  !> @todo Optimize special cases of quadratic formula solutions.
  REAL(DP) FUNCTION AlphaMiddle_Euler_Relativistic &
    ( DL, SL, tauL, F_DL, F_SL, F_tauL, DR, SR, tauR, F_DR, F_SR, F_tauR, &
      Gmii, aP, aM, Lapse, Shift, iErr )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: DL, SL, tauL, F_DL, F_SL, F_tauL, &
                               DR, SR, tauR, F_DR, F_SR, F_tauR, &
                               Gmii, aP, aM, Lapse, Shift
    INTEGER,  INTENT(inout) :: iErr

    REAL(DP) :: EL, F_EL, ER, F_ER, a2, a1, a0
    REAL(DP) :: E_HLL, S_HLL, FE_HLL, FS_HLL

#ifdef HYDRO_RIEMANN_SOLVER_HLL

    AlphaMiddle_Euler_Relativistic = 1.0e1_DP

    RETURN

#endif

    EL   = tauL + DL
    F_EL = F_tauL + F_DL
    ER   = tauR + DR
    F_ER = F_tauR + F_DR

    E_HLL  = aP * ER + aM * EL + Lapse * ( F_EL - F_ER )
    S_HLL  = aP * SR + aM * SL + Lapse * ( F_SL - F_SR )
    FE_HLL = Lapse * ( aP * F_EL + aM * F_ER ) - aM * aP * ( ER - EL )
    FS_HLL = Lapse * ( aP * F_SL + aM * F_SR ) - aM * aP * ( SR - SL )

    ! --- Coefficients in quadratic equation ---

    a2 = Gmii**2 * ( FE_HLL + Shift * E_HLL )
    a1 = -Gmii * ( Lapse * E_HLL + FS_HLL + Shift * S_HLL )
    a0 = Lapse * S_HLL

    ! --- Accounting for special cases of the solution to a
    !     quadratic equation when a2 = 0 ---

    IF     ( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a1 ) .LT. SqrtTiny ) &
            .AND. ( ABS( a0 ) .LT. SqrtTiny ) )THEN

      iErr = 9

    ELSE IF( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a1 ) .LT. SqrtTiny ) )THEN

      iErr = 9

    ELSE IF( ( ABS( a2 ) .LT. SqrtTiny ) .AND. ( ABS( a0 ) .LT. SqrtTiny ) )THEN

      AlphaMiddle_Euler_Relativistic = Zero

    ELSE IF( ABS( a2 ) .LT. SqrtTiny )THEN

      AlphaMiddle_Euler_Relativistic = -a0 / a1

    ELSE IF( ABS( a0 ) .LT. SqrtTiny )THEN

       AlphaMiddle_Euler_Relativistic = Zero

    ELSE

      AlphaMiddle_Euler_Relativistic &
        = ( -a1 - SQRT( MAX( a1**2 - Four * a2 * a0, SqrtTiny ) ) ) &
            / ( Two * a2 )
    END IF

    RETURN
  END FUNCTION AlphaMiddle_Euler_Relativistic


  !> Compute the physical flux in the X1-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X1_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X1_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X1_Euler_Relativistic(iCF_D)  &
      = D * W * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V1 - Shift / Lapse ) + P

    Flux_X1_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V1 - Shift / Lapse )

    Flux_X1_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V1 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X1_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V1 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X1_Euler_Relativistic


  !> Compute the physical flux in the X2-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X2_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X2_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X2_Euler_Relativistic(iCF_D)  &
      = D * W * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V2 - Shift / Lapse ) + P

    Flux_X2_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V2 - Shift / Lapse )

    Flux_X2_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V2 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X2_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V2 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X2_Euler_Relativistic


  !> Compute the physical flux in the X3-direction.
  !> @param Vi The ith contravariant components of the three-velocity.
  !> @param Gmii The ith covariant components of the spatial three-metric.
  !> @param Shift The first contravariant component of the shift-vector.
  FUNCTION Flux_X3_Euler_Relativistic &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, &
                            Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: VSq, W, h, Flux_X3_Euler_Relativistic(nCF)

    VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2
    W   = One / SQRT( One - VSq )
    h   = One + ( E + P ) / D

    Flux_X3_Euler_Relativistic(iCF_D)  &
      = D * W * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V3 - Shift / Lapse )

    Flux_X3_Euler_Relativistic(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V3 - Shift / Lapse ) + P

    Flux_X3_Euler_Relativistic(iCF_E)  &
      = D * W * ( h * W - One ) * ( V3 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X3_Euler_Relativistic(iCF_Ne) &
      = Ne * W * ( V3 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X3_Euler_Relativistic


  !> Compute the diagonal elements of the stress-tensor, needed for the
  !> source-terms in the hydro equations.
  !> @param Si The ith covariant components of the conserved momentum-density.
  !> @param Vi The ith contravavriant components of the three-velocity.
  FUNCTION StressTensor_Diagonal_Euler_Relativistic &
    ( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: StressTensor_Diagonal_Euler_Relativistic(3)

    StressTensor_Diagonal_Euler_Relativistic(1) = S1 * V1 + P
    StressTensor_Diagonal_Euler_Relativistic(2) = S2 * V2 + P
    StressTensor_Diagonal_Euler_Relativistic(3) = S3 * V3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal_Euler_Relativistic


  !> Compute the Local-Lax-Friedrichs numerical flux at a given element
  !> interface, in a given dimension.
  FUNCTION NumericalFlux_LLF_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    ! --- Local Lax-Friedrichs Flux ---

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), aP, aM

    REAL(DP) :: NumericalFlux_LLF_Euler_Relativistic(nCF)

    REAL(DP) :: alpha

    alpha = MAX( aM, aP )

    NumericalFlux_LLF_Euler_Relativistic &
      = Half * ( fL + fR - alpha * ( uR - uL ) )

    RETURN
  END FUNCTION NumericalFlux_LLF_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer numerical flux at a given element
  !> interface, in a given dimension.
  FUNCTION NumericalFlux_HLL_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), aP, aM

    REAL(DP) :: NumericalFlux_HLL_Euler_Relativistic(nCF)

    NumericalFlux_HLL_Euler_Relativistic &
      = ( aP * fL + aM * fR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X1-direction.
  !> @param Shift The first contravariant component of the shift-vector.
  !> @param Gm11 The first covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X1_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio, &
                NumericalFlux_X1_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S1) / Gm11 - Shift / Lapse * UE
        FS = uL(iCF_S1) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm11 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S1) - Lapse * FS ) ) &
             / ( Lapse - Gm11 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        S2 = uL(iCF_S2) * VelocityRatio

        S3 = uL(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S1) / Gm11 - Shift / Lapse * UE
        FS = uR(iCF_S1) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm11 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S1) - Lapse * FS ) ) &
               / ( Lapse - Gm11 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        S2 = uR(iCF_S2) * VelocityRatio

        S3 = uR(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_E)  &
        = S1 / Gm11 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X1_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X1_HLLC_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X2-direction.
  !> @param Shift The second contravariant component of the shift-vector.
  !> @param Gm22 The second covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X2_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm22, vL, vR, pL, pR, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm22, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio
    REAL(DP) :: NumericalFlux_X2_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S2) / Gm22 - Shift / Lapse * UE
        FS = uL(iCF_S2) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm22 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S2) - Lapse * FS ) ) &
             / ( Lapse - Gm22 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio

        S2 = uL(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        S3 = uL(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S2) / Gm22 - Shift / Lapse * UE
        FS = uR(iCF_S2) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm22 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S2) - Lapse * FS ) ) &
               / ( Lapse - Gm22 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio

        S2 = uR(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        S3 = uR(iCF_S3) * VelocityRatio

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_E)  &
        = S2 / Gm22 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X2_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X2_HLLC_Euler_Relativistic


  !> Compute the Harten-Lax-van-Leer-Contact numerical flux at a given element
  !> in the X3-direction.
  !> @param Shift The third contravariant component of the shift-vector.
  !> @param Gm33 The third covariant component of the spatial three-metric.
  FUNCTION NumericalFlux_X3_HLLC_Euler_Relativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm33, vL, vR, pL, pR, Lapse, Shift )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    ! --- Shift is the third contravariant component of the shift-vector
    !     Gm is the third covariant component of the spatial three-metric ---

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm33, vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: p, D, S1, S2, S3, E, Ne, UE, FE, FS, VelocityRatio
    REAL(DP) :: NumericalFlux_X3_HLLC_Euler_Relativistic(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_Relativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_Relativistic = fR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. Zero )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- uL_star ---
        UE = uL(iCF_E) + uL(iCF_D)
        FE = uL(iCF_S3) / Gm33 - Shift / Lapse * UE
        FS = uL(iCF_S3) * ( vL - Shift / Lapse ) + pL

        p  = ( Gm33 * aC * ( -aM * UE - Lapse * FE ) &
               - ( -aM * uL(iCF_S3) - Lapse * FS ) ) &
             / ( Lapse - Gm33 * aC * ( -aM + Shift ) )

        D  = uL(iCF_D)  * VelocityRatio

        S1 = uL(iCF_S1) * VelocityRatio

        S2 = uL(iCF_S2) * VelocityRatio

        S3 = uL(iCF_S3) * VelocityRatio &
               + Lapse * ( p - pL ) / ( -aM - Lapse * aC + Shift )

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pL * vL ) &
               / ( -aM - Lapse * aC + Shift )

        Ne = uL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- uR_star ---
        UE = uR(iCF_E) + uR(iCF_D)
        FE = uR(iCF_S3) / Gm33 - Shift / Lapse * UE
        FS = uR(iCF_S3) * ( vR - Shift / Lapse ) + pR

        p  = ( Gm33 * aC * ( aP * UE - Lapse * FE ) &
               - ( aP * uR(iCF_S3) - Lapse * FS ) ) &
               / ( Lapse - Gm33 * aC * ( aP + Shift ) )

        D  = uR(iCF_D)  * VelocityRatio

        S1 = uR(iCF_S1) * VelocityRatio

        S2 = uR(iCF_S2) * VelocityRatio

        S3 = uR(iCF_S3) * VelocityRatio &
               + Lapse * ( p - pR ) / ( aP - Lapse * aC + Shift )

        E  = UE         * VelocityRatio + Lapse * ( p * aC - pR * vR ) &
               / ( aP - Lapse * aC + Shift )

        Ne = uR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_D)  &
        = D  * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S1) &
        = S1 * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S2) &
        = S2 * ( aC - Shift / Lapse )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_S3) &
        = S3 * ( aC - Shift / Lapse ) + p

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_E)  &
        = S3 / Gm33 - D * aC - Shift / Lapse * ( E - D )

      NumericalFlux_X3_HLLC_Euler_Relativistic(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X3_HLLC_Euler_Relativistic


  ! --- Auxiliary utilities for ComputePrimitive ---


  SUBROUTINE GetBoundsForBisection &
    ( N, uD, uS1, uS2, uS3, uE, uNe, &
      Gm11, Gm22, Gm33, &
      dz, za, zb, fa, fb, q, r, k, iErr )

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: uD(N), uS1(N), uS2(N), &
                               uS3(N), uE(N), uNe(N), &
                               Gm11(N), Gm22(N), Gm33(N)
    REAL(DP), INTENT(out)   :: dz(N), za(N), zb(N), fa(N), fb(N), &
                               q(N), r(N), k(N)
    INTEGER,  INTENT(inout) :: iErr(N)

    INTEGER  :: iX
    REAL(DP) :: S

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( S )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( S ) &
    !$ACC PRESENT( uD, uS1, uS2, uS3, uE, uNe, &
    !$ACC          Gm11, Gm22, Gm33, &
    !$ACC          q, r, k, dz, za, zb, fa, fb, &
    !$ACC          iErr )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( S )
#endif
     DO iX = 1, N

       iErr(iX) = 0

       IF( uD(iX) .GT. Min_D )THEN

         ! --- Eq. C2 ---

         S = SQRT(   uS1(iX)**2 / Gm11(iX) &
                   + uS2(iX)**2 / Gm22(iX) &
                   + uS3(iX)**2 / Gm33(iX) )

         q(iX) = uE(iX) / uD(iX)
         r(iX) = S      / uD(iX)
         k(iX) = r (iX) / ( One + q(iX) )

         IF( q(iX) .LT. Zero )THEN

           r(iX) = k(iX)
           q(iX) = Zero

         END IF

         ! --- Eq. C23 ---

         za(iX) = Half * k(iX) / SQRT( One - Fourth * k(iX)**2 ) - SqrtTiny
         zb(iX) = k(iX)        / SQRT( One - k(iX)**2 )          + SqrtTiny

         dz(iX) = zb(iX) - za(iX)

         ! --- Eq. C25 ---

         CALL ComputeFunZ_Scalar &
                ( uD(iX), uNe(iX), r(iX), za(iX), q(iX), fa(iX) )
         CALL ComputeFunZ_Scalar &
                ( uD(iX), uNe(iX), r(iX), zb(iX), q(iX), fb(iX) )

         IF( .NOT. fa(iX) * fb(iX) .LT. Zero ) iErr(iX) = 8

       END IF ! IF( uD(iX) .GT. Min_D )

     END DO

  END SUBROUTINE GetBoundsForBisection


#ifdef MICROPHYSICS_WEAKLIB

  SUBROUTINE ComputeFunZ_Scalar( D, Ne, r, z, q, FunZ )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: D, Ne, r, z
    REAL(DP), INTENT(inout) :: q
    REAL(DP), INTENT(out)   :: FunZ

    REAL(DP) :: epst, at, ht, Wt, epsh, rhoh, ph, Ye, Min_E, Max_E

    ! --- Eq. C15 ---

    Wt = SQRT( One + z**2 )

    ! --- Eq. C16 ---

    epst = Wt * q - z * r + z**2 / ( One + Wt )

    ! --- Eq. C17 ---

    rhoh = MAX( MIN( Max_D, D / Wt ), Min_D )

    ! --- Eq. C18 ---

    Ye = Ne * AtomicMassUnit / D

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( rhoh, ( One + Offset_Temperature ) * Min_T, Ye, Min_E )
    CALL ComputeSpecificInternalEnergy_TABLE &
           ( rhoh, ( One - Offset_Temperature ) * Max_T, Ye, Max_E )

    Min_E = ( One + Offset_Epsilon ) * Min_E
    Max_E = ( One - Offset_Epsilon ) * Max_E

    epsh = MAX( MIN( Max_E, epst ), Min_E )

    ! --- Eq. C27 ---

    IF( epst .LT. Min_E )THEN

      q = ( One + q ) * ( One + epsh ) / ( One + epst ) - One

      epst = epsh

    ELSE IF( epst .GT. Max_E )THEN

      q = ( One + q ) * ( One + epsh ) / ( One + epst ) - One

      epst = epsh

    END IF

    ! --- Eqs. C19/C20 ---

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rhoh, epsh, Ye, ph )

    at = ph / ( rhoh * ( One + epsh ) )

    ! --- Eq. C21 ---

    ht = ( One + epst ) * ( One + at )

    ! --- Eq. C22 ---

    FunZ = z - r / ht

  END SUBROUTINE ComputeFunZ_Scalar


  SUBROUTINE ComputeFunZ_Vector( N, D, Ne, r, z, q, FunZ, ITERATE )

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: D(N), Ne(N), r(N), z(N)
    REAL(DP), INTENT(inout) :: q(N)
    REAL(DP), INTENT(out)   :: FunZ(N)
    LOGICAL,  INTENT(inout) :: ITERATE(N)

    REAL(DP) :: Wt, Min_E, Max_E, at, ht
    REAL(DP) :: epst(N), rhoh(N), Ye(N), epsh(N), ph(N)

    INTEGER :: iX

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    D, Ne, r, z, q, FunZ, ITERATE ) &
    !$OMP MAP( alloc: epst, rhoh, Ye, epsh, ph )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC        ENTER DATA &
    !$ACC COPYIN(     D, Ne, r, z, q, FunZ, ITERATE ) &
    !$ACC CREATE(     epst, rhoh, Ye, epsh, ph )
#endif

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Wt, Min_E, Max_E )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Wt, Min_E, Max_E ) &
    !$ACC PRESENT( ITERATE, z, epst, q, r, rhoh, D, Ye, Ne, epsh )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Wt, Min_E, Max_E )
#endif
    DO iX = 1, N

      IF( ITERATE(iX) )THEN

        ! --- Eq. C15 ---

        Wt = SQRT( One + z(iX)**2 )

        ! --- Eq. C16 ---

        epst(iX) = Wt * q(iX) - z(iX) * r(iX) + z(iX)**2 / ( One + Wt )

        ! --- Eq. C17 ---

        rhoh(iX) = MAX( MIN( Max_D, D(iX) / Wt ), Min_D )

        ! --- Eq. C18 ---

        Ye(iX) = Ne(iX) * AtomicMassUnit / D(iX)

        CALL ComputeSpecificInternalEnergy_TABLE &
               ( rhoh(iX), ( One + Offset_Temperature ) * Min_T, Ye(iX), Min_E )
        CALL ComputeSpecificInternalEnergy_TABLE &
               ( rhoh(iX), ( One - Offset_Temperature ) * Max_T, Ye(iX), Max_E )

        Min_E = ( One + Offset_Epsilon ) * Min_E
        Max_E = ( One - Offset_Epsilon ) * Max_E

        epsh(iX) = MAX( MIN( Max_E, epst(iX) ), Min_E )

        ! --- Eq. C27 ---

        IF( epst(iX) .LT. Min_E )THEN

          q(iX) = ( One + q(iX) ) &
                    * ( One + epsh(iX) ) / ( One + epst(iX) ) - One

          epst(iX) = epsh(iX)

        ELSE IF( epst(iX) .GT. Max_E )THEN

          q(iX) = ( One + q(iX) ) &
                    * ( One + epsh(iX) ) / ( One + epst(iX) ) - One

          epst(iX) = epsh(iX)

        END IF ! epst(iX) < Min_E || epst(iX) > Max_E

      END IF ! ITERATE(iX)

    END DO ! iX = 1, N

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( ITERATE, rhoh, epsh, Ye, ph )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iX = 1, N

      IF( ITERATE(iX) )THEN

        ! --- Eqs. C19/C20 ---

        CALL ComputePressureFromSpecificInternalEnergy &
               ( rhoh(iX), epsh(iX), Ye(iX), ph(iX) )

      END IF ! ITERATE(iX)

    END DO ! iX = 1, N

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( at, ht )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( at, ht ) &
    !$ACC PRESENT( ITERATE, ph, rhoh, epsh, epst, z, r, FunZ )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( at, ht )
#endif
    DO iX = 1, N

      IF( ITERATE(iX) )THEN

        at = ph(iX) / ( rhoh(iX) * ( One + epsh(iX) ) )

        ! --- Eq. C21 ---

        ht = ( One + epst(iX) ) * ( One + at )

        ! --- Eq. C22 ---

        FunZ(iX) = z(iX) - r(iX) / ht

      END IF ! ITERATE(iX)

    END DO ! iX = 1, N

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    FunZ ) &
    !$OMP MAP( release: D, Ne, r, z, q, ITERATE, &
    !$OMP               epst, rhoh, Ye, epsh, ph )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC        EXIT DATA &
    !$ACC COPYOUT(      FunZ ) &
    !$ACC DELETE(       D, Ne, r, z, q, ITERATE, &
    !$ACC               epst, rhoh, Ye, epsh, ph ) &
    !$ACC
#endif

  END SUBROUTINE ComputeFunZ_Vector

#else

  SUBROUTINE ComputeFunZ_Scalar( D, Ne, r, z, q, FunZ )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ne, r, z, q
    REAL(DP), INTENT(out) :: FunZ

    REAL(DP) :: Wt, rhot, epst, pt, at, Ye, ht

    ! --- Eq. C15 ---

    Wt = SQRT( One + z**2 )
    rhot = D / Wt

    ! --- Eq. C16 ---

    epst = Wt * q - z * r + z**2 / ( One + Wt )

    Ye = Ne * AtomicMassUnit / D

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rhot, epst, Ye, pt )

    ! --- Eq. C20 ---

    at = pt / ( rhot * ( One + epst ) )

    ! --- Eq. C21 ---

    ht = ( One + epst ) * ( One + at )

    ! --- Eq. C22 ---

    FunZ = z - r / ht

  END SUBROUTINE ComputeFunZ_Scalar


  SUBROUTINE ComputeFunZ_Vector( N, D, Ne, r, z, q, FunZ, ITERATE )

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: D(N), Ne(N), r(N), z(N), q(N)
    REAL(DP), INTENT(out)   :: FunZ(N)
    LOGICAL,  INTENT(inout) :: ITERATE(N)

    INTEGER  :: iX
    REAL(DP) :: Wt, rhot, epst, Ye, pt, at, ht

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Wt, rhot, epst, Ye, pt, at, ht )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( D, Ne, r, z, q, FunZ, ITERATE ) &
    !$ACC PRIVATE( Wt, rhot, epst, Ye, pt, at, ht )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Wt, rhot, epst, Ye, pt, at, ht )
#endif
    DO iX = 1, N

      IF( ITERATE(iX) )THEN

        ! --- Eq. C15 ---

        Wt = SQRT( One + z(iX)**2 )
        rhot = D(iX) / Wt

        ! --- Eq. C16 ---

        epst = Wt * q(iX) - z(iX) * r(iX) + z(iX)**2 / ( One + Wt )

        Ye = Ne(iX) * AtomicMassUnit / D(iX)

        CALL ComputePressureFromSpecificInternalEnergy &
               ( rhot, epst, Ye, pt )

        ! --- Eq. C20 ---

        at = pt / ( rhot * ( One + epst ) )

        ! --- Eq. C21 ---

        ht = ( One + epst ) * ( One + at )

        ! --- Eq. C22 ---

        FunZ(iX) = z(iX) - r(iX) / ht

      END IF ! ITERATE(iX)

    END DO ! iX = 1, N

  END SUBROUTINE ComputeFunZ_Vector

#endif

  SUBROUTINE SolveZ_Bisection_Scalar &
    ( CF_D, CF_Ne, q, r, k, z0, iErr, dzMin_Option )

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: CF_D, CF_Ne, r, k
    REAL(DP), INTENT(inout) :: q
    REAL(DP), INTENT(out)   :: z0
    INTEGER,  INTENT(inout) :: iErr
    REAL(DP), INTENT(in), OPTIONAL :: dzMin_Option

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: za, zb, zc, dz, dzMin
    REAL(DP) :: fa, fb, fc
    INTEGER  :: MAX_IT

    dzMin = 1.0e-08_DP
    IF( PRESENT( dzMin_Option ) ) &
      dzMin = dzMin_Option

    MAX_IT = 4 - INT( LOG( dzMin ) / LOG( Two ) )

    ! --- Eq. C23 ---

    za = SQRT( One - Fourth * k**2 )
    zb = SQRT( One - k**2 )
    za = Half * k / za - SqrtTiny
    zb = k        / zb + SqrtTiny

    ! --- Compute FunZ for upper and lower bounds ---

    CALL ComputeFunZ_Scalar( CF_D, CF_Ne, r, za, q, fa )
    CALL ComputeFunZ_Scalar( CF_D, CF_Ne, r, zb, q, fb )

    ! --- Check that sign of FunZ changes across bounds ---

    IF( .NOT. fa * fb .LT. 0 ) &
      iErr = 8

    dz = zb - za

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

      ! --- Compute midpoint, zc ---

      dz = Half * dz

      ! --- Bisection ---

      zc = za + dz

!!$      ! --- Regula Falsi ---
!!$
!!$      zc = ( za * fb - zb * fa  ) / ( fb - fa )

!!$      ! --- Illinois ---
!!$
!!$      IF( fa * fc .GT. Zero )THEN
!!$
!!$        zc = ( fb * za - Half * fa * zb ) / ( fb - Half * fa )
!!$
!!$      ELSE
!!$
!!$        zc = ( Half * fb * za - fa * zb ) / ( Half * fb - fa )
!!$
!!$      ENDIF

      ! --- Compute f(zc) for midpoint zc ---

      CALL ComputeFunZ_Scalar( CF_D, CF_Ne, r, zc, q, fc )

      ! --- Change zc to za or zb, depending on sign of fc ---

      IF( fa * fc .LT. Zero )THEN

        zb = zc
        fb = fc

      ELSE IF( fa * fc .GT. Zero )THEN

        za = zc
        fa = fc

      ELSE

        CONVERGED = .TRUE.

      END IF

      IF( ABS( dz ) / MAX( ABS( zc ), SqrtTiny ) .LE. dzMin ) &
        CONVERGED = .TRUE.

!!$      IF( ITERATION .GT. MAX_IT - 3 )THEN
!!$
!!$        WRITE(*,*) 'Iter   = ', ITERATION
!!$        WRITE(*,*) 'za, zb = ', za, zb
!!$        WRITE(*,*) 'dz     = ', dz
!!$        WRITE(*,*)
!!$
!!$      END IF

    END DO

    z0 = zc

  END SUBROUTINE SolveZ_Bisection_Scalar


  SUBROUTINE SolveZ_Bisection_Vector &
    ( N, uD, uNe, q, r, k, &
      za, zb, fa, fb, &
      dz, zc, fc, ITERATE, dzMin_Option )

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: uD(N), uNe(N), r(N), k(N)
    REAL(DP), INTENT(inout) :: za(N), zb(N), fa(N), fb(N), q(N), dz(N)
    REAL(DP), INTENT(out)   :: zc(N), fc(N)
    LOGICAL,  INTENT(inout) :: ITERATE(N)
    REAL(DP), INTENT(in), OPTIONAL :: dzMin_Option

    REAL(DP) :: dzMin
    INTEGER  :: MAX_IT

    INTEGER :: ITERATION, iX

    dzMin = 1.0e-08
    IF( PRESENT( dzMin_Option ) ) &
      dzMin = dzMin_Option

    MAX_IT = 4 - INT( LOG( dzMin ) / LOG( Two ) )

    ITERATION = 0
    DO WHILE( ANY( ITERATE ) .AND. ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

      ! --- Compute midpoints, zc ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( za, zc, dz, uD, ITERATE )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO
#endif
      DO iX = 1, N

        IF( uD(iX) .GT. Min_D .AND. ITERATE(iX) )THEN

          dz(iX) = Half * dz(iX)
          zc(iX) = za(iX) + dz(iX)

        END IF

      END DO

      CALL ComputeFunZ_Vector( N, uD, uNe, r, zc, q, fc, ITERATE )

      ! --- Change either za or zb to zc, depending on sign of fc ---

#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( uD, uNe, q, r, k, &
      !$ACC          za, zb, fa, fb, &
      !$ACC          dz, zc, fc, ITERATE )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO
#endif
      DO iX = 1, N

        IF( uD(iX) .GT. Min_D .AND. ITERATE(iX) )THEN

          IF( fa(iX) * fc(iX) .LT. Zero )THEN

            zb(iX) = zc(iX)
            fb(iX) = fc(iX)

          ELSE

            za(iX) = zc(iX)
            fa(iX) = fc(iX)

          END IF

          !IF( ABS( dz(iX) ) / MAX( ABS( zc(iX) ), SqrtTiny ) .LE. dzMin ) &
          IF( ABS( dz(iX) ) .LE. ABS( zc(iX) ) * dzMin ) &
            ITERATE(iX) = .FALSE.

!!$          IF( ITERATION .GT. MAX_IT - 3 )THEN
!!$
!!$            WRITE(*,*) 'iX     = ', iX
!!$            WRITE(*,*) 'Iter   = ', ITERATION
!!$            WRITE(*,*) 'za, zb = ', za(iX), zb(iX)
!!$            WRITE(*,*) 'dz     = ', dz(iX)
!!$            WRITE(*,*)
!!$
!!$          END IF

        ELSE

          ITERATE(iX) = .FALSE.

        END IF ! uD(iX) .GT. Min_D .AND. ITERATE(iX)

      END DO ! iX
#if   defined( THORNADO_OMP_OL ) && !defined( THORNADO_EULER_NOGPU )
      !$OMP TARGET UPDATE FROM( ITERATE )
#elif defined( THORNADO_OACC   ) && !defined( THORNADO_EULER_NOGPU )
      !$ACC UPDATE HOST       ( ITERATE )
#endif

    END DO ! WHILE( ANY( ITERATE ) .AND. ITERATION .LT. MAX_IT )

  END SUBROUTINE SolveZ_Bisection_Vector

END MODULE Euler_UtilitiesModule_Relativistic
