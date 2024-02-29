MODULE MF_InitializationModule

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    swX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL, &
    Gamma_IDEAL
  USE UnitsModule, ONLY: &
    Gram, &
    Erg, &
    Centimeter, &
    Kilometer, &
    SolarMass
  USE GeometryFieldsModule, ONLY: &
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
    iGF_Beta_3, &
    iGF_Psi, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    iAF_P
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE MF_KindModule, ONLY: &
    DP, &
    SqrtTiny, &
    Zero, &
    One, &
    Two, &
    Three, &
    TwoPi, &
    FourPi
  USE InputParsingModule, ONLY: &
    xR, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

  REAL(DP) :: PolytropicConstant
  REAL(DP) :: CentralPressure
  REAL(DP) :: GravitationalMass
  REAL(DP) :: NeutronStarRadius
  REAL(DP), PARAMETER :: DeltaR = 1.0e-4_DP * Kilometer
  REAL(DP), PARAMETER :: TOLERANCE = 1.0e-10_DP

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF, MF_uPF, MF_uAF

    LOGICAL :: Verbose

    TYPE(amrex_parmparse) :: PP

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    INTEGER :: iX_B(3), iX_E(3), iX1, iX2, iX3, iNX

    CHARACTER(LEN=32) :: FMT

    ! --- Problem-Dependent Parameters ---

    INTEGER  :: N, iNX1, iNX2, iX_L
    REAL(DP) :: X1, X2

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      RadiusArr, DensityArr, PressureArr, AlphaArr, PsiArr

    REAL(DP) :: CentralDensity

    CALL amrex_parmparse_build( PP, 'TOV' )
      CALL PP % get( 'CentralDensity'     , CentralDensity     )
      CALL PP % get( 'PolytropicConstant' , PolytropicConstant )
    CALL amrex_parmparse_destroy( PP )

    CentralDensity &
      = CentralDensity  * ( Gram / Centimeter**3 )
    PolytropicConstant &
      = PolytropicConstant &
          * ( Erg / Centimeter**3 / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) )

    CentralPressure = PolytropicConstant * CentralDensity**( Gamma_IDEAL )

    N = ( xR(1) + 5.0_DP * Kilometer ) / DeltaR

    ALLOCATE( RadiusArr  (N) )
    ALLOCATE( DensityArr (N) )
    ALLOCATE( PressureArr(N) )
    ALLOCATE( AlphaArr   (N) )
    ALLOCATE( PsiArr     (N) )

    CALL IntegrateTOVEquation &
           ( N, RadiusArr, DensityArr, PressureArr, AlphaArr, PsiArr )

    Verbose = .FALSE.
    IF( amrex_parallel_ioprocessor() .AND. iLevel .EQ. 0 ) Verbose = .TRUE.

    IF( Verbose )THEN

      WRITE(FMT,'(A)') '(A52,ES11.3E3)'

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'Initializing: ', TRIM( ProgramName )
      WRITE(*,'(4x,A)')   '-------------'
      WRITE(*,*)
      WRITE(*,TRIM(FMT)) &
        'Central Density [g/cm^3]: ', &
          CentralDensity / ( Gram / Centimeter**3 )
      WRITE(*,TRIM(FMT)) &
        'PolytropicConstant [erg/cm^3/(g/cm^3)^Gamma]: ', &
        PolytropicConstant &
          / ( Erg / Centimeter**3 / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) )
      WRITE(*,TRIM(FMT)) &
        'Central Pressure [erg/cm^3]: ', &
          CentralPressure / ( Erg / Centimeter**3 )
      WRITE(*,TRIM(FMT)) &
        'Gravitational Mass [Msun]: ', &
           GravitationalMass / SolarMass
      WRITE(*,TRIM(FMT)) &
        'Neutron Star Radius [km]: ', &
          NeutronStarRadius / Kilometer
      WRITE(*,*)

    END IF

    ! --- Map to computational domain ---

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uPF => MF_uPF % DataPtr( MFI )
      uAF => MF_uAF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B = BX % lo
      iX_E = BX % hi

      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) ) &
        iX_E(1) = iX_E(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)
      DO iNX = 1      , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        iX_L = Locate( X1, RadiusArr, N )

        uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
          = Interpolate1D_Linear &
              ( X1, RadiusArr(iX_L), RadiusArr(iX_L+1), &
                    AlphaArr (iX_L), AlphaArr (iX_L+1) )

        uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) &
          = Interpolate1D_Linear &
              ( X1, RadiusArr(iX_L), RadiusArr(iX_L+1), &
                    PsiArr   (iX_L), PsiArr   (iX_L+1) )

        uGF    (iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**2
        uGF    (iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**2 * X1
        uGF    (iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**2 * X1 * SIN( X2 )

        uGF    (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1     -1)+iNX)**2
        uGF    (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2     -1)+iNX)**2
        uGF    (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
          = uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3     -1)+iNX)**2

        uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) &
          =   uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
            * uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
            * uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX)

        uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1 -1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2 -1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3 -1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX) = Zero
        uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX) = Zero

        uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX) &
          = Interpolate1D_Linear &
              ( X1, RadiusArr (iX_L), RadiusArr (iX_L+1), &
                    DensityArr(iX_L), DensityArr(iX_L+1) )

        uPF(iX1,iX2,iX3,nDOFX*(iPF_V1-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V2-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V3-1)+iNX) = Zero

        uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) &
          = Interpolate1D_Linear &
              ( X1, RadiusArr  (iX_L), RadiusArr  (iX_L+1), &
                    PressureArr(iX_L), PressureArr(iX_L+1) ) &
              / ( Gamma_IDEAL - One )

        uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX) = Zero

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_E -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX), &
                 uAF(iX1,iX2,iX3,nDOFX*(iAF_P -1)+iNX) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D       -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V1      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V2      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V3      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_E       -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_D       -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S1      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S2      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S3      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_E       -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne      -1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                 uAF(iX1,iX2,iX3,nDOFX*(iAF_P       -1)+iNX) )

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE MFI % next()

    CALL amrex_mfiter_destroy( MFI )

    DEALLOCATE( PsiArr      )
    DEALLOCATE( AlphaArr    )
    DEALLOCATE( PressureArr )
    DEALLOCATE( DensityArr  )
    DEALLOCATE( RadiusArr   )

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE FUNCTIONS/SUBROUTINES ---


  SUBROUTINE IntegrateTOVEquation &
    ( N, RadiusArr, DensityArr, PressureArr, AlphaArr, PsiArr )

    INTEGER , INTENT(in)  :: N
    REAL(DP), INTENT(out) :: RadiusArr(N), DensityArr(N), PressureArr(N), &
                             AlphaArr(N), PsiArr(N)

    LOGICAL :: CONVERGED
    INTEGER :: ITER, iX1, iXX1
    INTEGER, PARAMETER :: MAX_ITER = 100

    REAL(DP) :: Radius, Pressure, Psi, Phi, Epsi, Ephi
    REAL(DP) :: PressureN, EpsiN, EphiN
    REAL(DP) :: Alpha_Iso, Psi_Iso, dAlpha, dPsi
    REAL(DP) :: Alpha0
    REAL(DP) :: Psi0

    Alpha0 = 0.8_DP
    Psi0   = 1.4_DP

    Radius   = SqrtTiny * Kilometer
    Pressure = CentralPressure
    Psi      = Psi0
    Phi      = Alpha0 * Psi0
    Epsi     = Zero
    Ephi     = Zero
    GravitationalMass = Zero

    ITER = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      DO iX1 = 1, N

        RadiusArr  (iX1) = Radius
        PressureArr(iX1) = Pressure
        DensityArr (iX1) &
          = ( Pressure / PolytropicConstant )**( One / Gamma_IDEAL )

        PressureN &
          = Pressure &
              + DeltaR &
                  * dpdr &
                      ( Radius, DensityArr(iX1), Pressure, &
                        Psi, Phi, Epsi, Ephi )

        EpsiN = Epsi + DeltaR * dEpsidr( Radius, Pressure, Psi )
        EphiN = Ephi + DeltaR * dEphidr( Radius, Pressure, Psi, Phi )

        Radius = Radius + DeltaR
        Psi    = Psi    + DeltaR * dPsidr( Radius, EpsiN )
        Phi    = Phi    + DeltaR * dPhidr( Radius, EphiN )

        PsiArr  (iX1) = Psi
        AlphaArr(iX1) = Phi / Psi

        Pressure = PressureN
        Epsi     = EpsiN
        Ephi     = EphiN

        GravitationalMass &
          = GravitationalMass &
              + FourPi * Radius**2 &
                  * ( DensityArr(iX1) + Pressure / ( Gamma_IDEAL - One ) &
                        + Pressure + Three * Pressure ) &
                  * Phi * Psi**5 * DeltaR

        IF( Pressure .LT. 1.0e-8_DP * CentralPressure )THEN

          NeutronStarRadius = Radius

          Alpha_Iso =   ( One - GravitationalMass / ( Two * Radius ) ) &
                      / ( One + GravitationalMass / ( Two * Radius ) )
          Psi_Iso   = One + GravitationalMass / ( Two * Radius )

          dAlpha = Alpha_Iso - AlphaArr(iX1)
          dPsi   = Psi_Iso   - PsiArr  (iX1)

          IF( MAX( ABS( dAlpha ) / Alpha_Iso, ABS( dPsi ) / Psi_Iso ) &
                .LT. TOLERANCE )THEN

            CONVERGED = .TRUE.

            ! --- For Locate function to work properly ---
            DO iXX1 = iX1, N

              RadiusArr(iXX1) = RadiusArr(iXX1-1) + DeltaR

            END DO

          ELSE

            Alpha0 = Alpha0 + dAlpha
            Psi0   = Psi0   + dPsi

            Radius   = SqrtTiny * Kilometer
            Pressure = CentralPressure
            Psi      = Psi0
            Phi      = Alpha0 * Psi0
            Epsi     = Zero
            Ephi     = Zero
            GravitationalMass = Zero

          END IF

          EXIT

        END IF ! Pressure .LT. 1.0e-8_DP * CentralPressure

      END DO ! iX1 = 1, N

    END DO ! WHILE

  END SUBROUTINE IntegrateTOVEquation


  REAL(DP) FUNCTION dpdr( r, rho, p, Psi, Phi, Epsi, Ephi )

    REAL(DP), INTENT(in) :: r, rho, p, Psi, Phi, Epsi, Ephi

    REAL(DP) :: Lapse, e

    e = p / ( Gamma_IDEAL - One )

    Lapse = Phi / Psi

    dpdr = - ( rho + e + p ) &
             * ( dPhidr( r, Ephi ) - Lapse * dPsidr( r, Epsi ) ) / Phi

    RETURN
  END FUNCTION dpdr


  REAL(DP) FUNCTION dEpsidr( r, p, Psi )

    REAL(DP), INTENT(in) :: r, p, Psi

    REAL(DP) :: rho, e

    rho = ( p / PolytropicConstant )**( One / Gamma_IDEAL )

    e = p / ( Gamma_IDEAL - One )

    dEpsidr = TwoPi * Psi**5 * ( rho + e ) * r**2

    RETURN
  END FUNCTION dEpsidr


  REAL(DP) FUNCTION dEphidr( r, p, Psi, Phi )

    REAL(DP), INTENT(in) :: r, p, Psi, Phi

    REAL(DP) :: rho, e

    rho = ( p / PolytropicConstant )**( One / Gamma_IDEAL )

    e = p / ( Gamma_IDEAL - One )

    dEphidr = TwoPi * Phi * Psi**4 * ( rho + e + 6.0_DP * p ) * r**2

    RETURN
  END FUNCTION dEphidr


  REAL(DP) FUNCTION dPsidr( r, Epsi )

    REAL(DP), INTENT(in) :: r, Epsi

    dPsidr = -Epsi / r**2

    RETURN
  END FUNCTION dPsidr


  REAL(DP) FUNCTION dPhidr( r, Ephi )

    REAL(DP), INTENT(in) :: r, Ephi

    dPhidr = Ephi / r**2

    RETURN
  END FUNCTION dPhidr


END MODULE MF_InitializationModule
