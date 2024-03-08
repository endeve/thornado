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
  USE GeometryComputationModule, ONLY: &
    ConformalFactor, &
    LapseFunction
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
    Three, &
    TwoPi, &
    FourPi, &
    Half
  USE InputParsingModule, ONLY: &
    xR, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

  REAL(DP) :: PolytropicConstant
  REAL(DP) :: CentralPressure
  REAL(DP) :: CentralDensity
  REAL(DP) :: GravitationalMass
  REAL(DP) :: NeutronStarRadius
  REAL(DP), PARAMETER :: DeltaR    = 1.0e-4_DP * Kilometer
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

    CALL amrex_parmparse_build( PP, 'TOV' )
      CALL PP % get( 'CentralDensity'    , CentralDensity     )
      CALL PP % get( 'PolytropicConstant', PolytropicConstant )
    CALL amrex_parmparse_destroy( PP )

    CentralDensity &
      = CentralDensity  * ( Gram / Centimeter**3 )
    PolytropicConstant &
      = PolytropicConstant &
          * ( Erg / Centimeter**3 / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) )

    CentralPressure = PolytropicConstant * CentralDensity**( Gamma_IDEAL )

    N = ( xR(1) + 3.0_DP * Kilometer ) / DeltaR

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
    INTEGER :: ITER
    INTEGER, PARAMETER :: MAX_ITER = 100
    REAL(DP) :: Psi0, Phi0

    Psi0 = 1.4_DP
    Phi0 = 1.2_DP

    ITER = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      CALL IntegrateOutwards &
             ( N, RadiusArr, DensityArr, PressureArr, AlphaArr, PsiArr, &
               CONVERGED, Psi0, Phi0 )

    END DO ! WHILE

  END SUBROUTINE IntegrateTOVEquation


  SUBROUTINE IntegrateOutwards &
    ( N, RadiusArr, DensityArr, PressureArr, AlphaArr, PsiArr, &
      CONVERGED, Psi0, Phi0 )

    INTEGER , INTENT(in)    :: N
    REAL(DP), INTENT(out)   :: RadiusArr(N), DensityArr(N), PressureArr(N), &
                               AlphaArr(N), PsiArr(N)
    LOGICAL , INTENT(inout) :: CONVERGED
    REAL(DP), INTENT(inout) :: Psi0, Phi0

    INTEGER :: i, j

    REAL(DP) :: ra, Epsi_ra, Ephi_ra, PsiPrime_ra, PhiPrime_ra, Mass
    REAL(DP) :: Radius, Density, Pressure
    REAL(DP) :: Psi, Phi, Epsi, Ephi
    REAL(DP) :: PsiR, PhiR, dPsi, dPhi
    REAL(DP) :: PressureN, InternalEnergyDensity, PsiPrime, PhiPrime

    ra          = Zero
    Epsi_ra     = Zero
    Ephi_ra     = Zero
    PsiPrime_ra = Zero
    PhiPrime_ra = Zero
    Mass        = Zero

    Radius   = SqrtTiny
    Density  = CentralDensity
    Pressure = CentralPressure
    Psi      = Psi0
    Phi      = Phi0
    Epsi     = Epsi_ra
    Ephi     = Epsi_ra

    RadiusArr  (1) = Radius
    PressureArr(1) = Pressure
    DensityArr (1) = Density
    PsiArr     (1) = Psi
    AlphaArr   (1) = Phi / Psi

    DO i = 2, N

      Density &
        = ( Pressure / PolytropicConstant )**( One / Gamma_IDEAL )

      InternalEnergyDensity &
        = Pressure / ( Gamma_IDEAL - One )

      PsiPrime = dPsidr( Radius, Epsi, ra, PsiPrime_ra, Epsi_ra )
      PhiPrime = dPhidr( Radius, Ephi, ra, PhiPrime_ra, Ephi_ra )

      PressureN &
        = Pressure &
            + DeltaR &
                * dpdr( Radius, Density, InternalEnergyDensity, Pressure, &
                        Phi, PhiPrime, Psi, PsiPrime )

      Epsi &
        = Epsi &
            + DeltaR * dEpsidr( Radius, Density, InternalEnergyDensity, Psi )
      Ephi &
        = Ephi &
            + DeltaR * dEphidr( Radius, Density, InternalEnergyDensity, &
                                Pressure, Psi, Phi )

      Radius = Radius + DeltaR
      Psi    = Psi + DeltaR * dPsidr( Radius, Epsi, ra, PsiPrime_ra, Epsi_ra )
      Phi    = Phi + DeltaR * dPhidr( Radius, Ephi, ra, PhiPrime_ra, Ephi_ra )

      Pressure = PressureN

      Density &
        = ( Pressure / PolytropicConstant )**( One / Gamma_IDEAL )

      InternalEnergyDensity &
        = Pressure / ( Gamma_IDEAL - One )

      Mass &
        = Mass &
            + FourPi * Radius**2 &
                * ( Density + InternalEnergyDensity + Three * Pressure ) &
                * Phi * Psi**5 * DeltaR

      RadiusArr  (i) = Radius
      DensityArr (i) = Density
      PressureArr(i) = Pressure
      PsiArr     (i) = Psi
      AlphaArr   (i) = Phi / Psi

      IF( Pressure .LT. 1.0e-8_DP * CentralPressure )THEN

        NeutronStarRadius = Radius
        GravitationalMass = Mass

        PsiR = Psi_Iso( NeutronStarRadius, GravitationalMass )
        PhiR = Phi_Iso( NeutronStarRadius, GravitationalMass )

        dPsi = PsiR - Psi
        dPhi = PhiR - Phi

        IF( MAX( ABS( dPsi ) / PsiR, ABS( dPhi ) / PhiR ) .LT. TOLERANCE )THEN

          CONVERGED = .TRUE.

          ! --- Extrapolate beyond neutron star surface
          !     with constant hydro and isotropic metric ---

          DO j = i + 1, N

            RadiusArr  (j) = RadiusArr(j-1) + DeltaR
            PressureArr(j) = Pressure
            DensityArr (j) = Density
            PsiArr     (j) = Psi_Iso( RadiusArr(j), GravitationalMass )
            AlphaArr   (j) = Phi_Iso( RadiusArr(j), GravitationalMass )

          END DO

        ELSE

          Psi0 = Psi0 + dPsi
          Phi0 = Phi0 + dPhi

        END IF

        EXIT

      END IF ! Pressure .LT. 1.0e-8_DP * CentralPressure

    END DO ! i = 2, N

  END SUBROUTINE IntegrateOutwards


  REAL(DP) FUNCTION dpdr( r, rho, e, p, Phi, PhiPrime, Psi, PsiPrime )

    REAL(DP), INTENT(in) :: r, rho, e, p, Phi, PhiPrime, Psi, PsiPrime

    dpdr = - ( rho + e + p ) * ( PhiPrime / Phi - PsiPrime / Psi )

    RETURN
  END FUNCTION dpdr


  REAL(DP) FUNCTION dEpsidr( r, rho, e, Psi )

    REAL(DP), INTENT(in) :: r, rho, e, Psi

    dEpsidr = TwoPi * Psi**5 * ( rho + e ) * r**2

    RETURN
  END FUNCTION dEpsidr


  REAL(DP) FUNCTION dEphidr( r, rho, e, p, Psi, Phi )

    REAL(DP), INTENT(in) :: r, rho, e, p, Psi, Phi

    dEphidr = TwoPi * Phi * Psi**4 * ( rho + e + 6.0_DP * p ) * r**2

    RETURN
  END FUNCTION dEphidr


  REAL(DP) FUNCTION dPsidr( r, Epsi, ra, PsiPrime_ra, Epsi_ra )

    REAL(DP), INTENT(in) :: r, Epsi, ra, PsiPrime_ra, Epsi_ra

    dPsidr = One / r**2 * ( ra**2 * PsiPrime_ra - Epsi + Epsi_ra )

    RETURN
  END FUNCTION dPsidr


  REAL(DP) FUNCTION dPhidr( r, Ephi, ra, PhiPrime_ra, Ephi_ra )

    REAL(DP), INTENT(in) :: r, Ephi, ra, PhiPrime_ra, Ephi_ra

    dPhidr = One / r**2 * ( ra**2 * PhiPrime_ra + Ephi - Ephi_ra )

    RETURN
  END FUNCTION dPhidr


  REAL(DP) FUNCTION Psi_Iso( r, M )

    REAL(DP), INTENT(in) :: r, M

    Psi_Iso = ConformalFactor( r, M )

    RETURN
  END FUNCTION Psi_Iso


  REAL(DP) FUNCTION Phi_Iso( r, M )

    REAL(DP), INTENT(in) :: r, M

    Phi_Iso = Psi_Iso( r, M ) * LapseFunction( r, M )

    RETURN
  END FUNCTION Phi_Iso


END MODULE MF_InitializationModule
