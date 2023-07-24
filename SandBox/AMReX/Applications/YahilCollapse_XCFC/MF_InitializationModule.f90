MODULE MF_InitializationModule

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_mfiter
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse_build, &
    amrex_parmparse_destroy, &
    amrex_parmparse

  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL
  USE UnitsModule, ONLY: &
    Gram, &
    Erg, &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    GravitationalConstant
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
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
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh, &
    iAF_Gm
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    FourPi
  USE InputParsingModule, ONLY: &
    ProgramName, &
    nLevels, &
    xL, &
    xR, &
    nX, &
    swX, &
    UseTiling, &
    Gamma_IDEAL, &
    iOS_CPP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    LOGICAL :: Verbose

    TYPE(amrex_parmparse) :: PP

    ! --- Problem-Dependent Parameters ---

    REAL(DP) :: D0
    REAL(DP) :: CentralDensity
    REAL(DP) :: CentralPressure
    REAL(DP) :: CoreRadius
    REAL(DP) :: CollapseTime
    REAL(DP) :: PolytropicConstant
    REAL(DP) :: dXdr, drhodD, dvdV, dmdM, TotalEnclosedMass

    CALL amrex_parmparse_build( PP, 'YC' )
      CALL PP % get( 'D0'             , D0 )
      CALL PP % get( 'CentralDensity' , CentralDensity )
      CALL PP % get( 'CentralPressure', CentralPressure )
      CALL PP % get( 'CoreRadius'     , CoreRadius )
      CALL PP % get( 'CollapseTime'   , CollapseTime )
    CALL amrex_parmparse_destroy( PP )

    Verbose = .FALSE.
    IF( amrex_parallel_ioprocessor() .AND. iLevel .EQ. 0 ) Verbose = .TRUE.

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'Initializing: ', TRIM( ProgramName )
      WRITE(*,'(4x,A)')   '-------------'
      WRITE(*,*)
      WRITE(*,'(6x,A,ES11.3E3)') '             D0: ', D0
      WRITE(*,'(6x,A,ES11.3E3)') ' CentralDensity: ', CentralDensity
      WRITE(*,'(6x,A,ES11.3E3)') 'CentralPressure: ', CentralPressure
      WRITE(*,'(6x,A,ES11.3E3)') '     CoreRadius: ', CoreRadius
      WRITE(*,'(6x,A,ES11.3E3)') '   CollapseTime: ', CollapseTime
      WRITE(*,*)

    END IF

    CentralDensity  = CentralDensity  * ( Gram / Centimeter**3 )
    CentralPressure = CentralPressure * ( Erg  / Centimeter**3 )
    CoreRadius      = CoreRadius      * ( Kilometer )
    CollapseTime    = CollapseTime    * ( Millisecond )

    PolytropicConstant = CentralPressure / CentralDensity**Gamma_IDEAL

    dXdr   = PolytropicConstant**( -Half ) &
               * GravitationalConstant**( ( Gamma_IDEAL - One ) / Two ) &
               * CollapseTime**( Gamma_IDEAL - Two )
    drhodD = GravitationalConstant**( -1 ) * CollapseTime**( -2 )
    dvdV   = PolytropicConstant**( Half ) &
               * GravitationalConstant**( ( One - Gamma_IDEAL ) / Two ) &
               * CollapseTime**( One - Gamma_IDEAL )
    dmdM   = PolytropicConstant**( Three / Two ) &
               * CollapseTime**( Four - Three * Gamma_IDEAL ) &
               * GravitationalConstant**( ( One - Three * Gamma_IDEAL ) / Two )

    CALL InitializeFields_YahilCollapse_FromScratch &
           ( iLevel, MF_uGF, MF_uCF, &
             dXdr, drhodD, dvdV, dmdM, PolytropicConstant, &
             CoreRadius, D0, CollapseTime, TotalEnclosedMass )

  END SUBROUTINE InitializeFields_MF


  ! --- Auxiliary functions for Yahil collapse problem ---


  SUBROUTINE IntegrateD( dX, X, D, U, M, Numer, Denom )

    REAL(DP), INTENT(in)    :: dX
    REAL(DP), INTENT(inout) :: X(:), D(:), U(:), M(:), Numer(:), Denom(:)

    REAL(DP)            :: dDdX, dMdX, XC, X0, &
                           NumerC, DenomC, NumerPrime, DenomPrime
    INTEGER             :: iX1
    LOGICAL             :: FirstTime
    REAL(DP), PARAMETER :: Threshold = 0.015_DP

    FirstTime = .TRUE.

    DO iX1 = 2, SIZE(X)

      dDdX = Numer(iX1-1) / Denom(iX1-1)
      dMdX = FourPi * X(iX1-1)**2 * D(iX1-1)

      X(iX1) = X(iX1-1) + dX
      D(iX1) = D(iX1-1) + dX * dDdX
      M(iX1) = M(iX1-1) + dX * dMdX

      U(iX1) = ( Four - Three * Gamma_IDEAL ) * M(iX1) &
                 / ( FourPi * X(iX1)**2 * D(iX1) )

      Numer(iX1) = Numerator  ( X(iX1), D(iX1), U(iX1), M(iX1) )
      Denom(iX1) = Denominator( D(iX1), U(iX1) )

      IF( ABS( Denom(iX1) ) .LT. Threshold .AND. FirstTime )THEN

        XC     = X(iX1)
        NumerC = Numer(iX1)
        DenomC = Denom(iX1)

        DenomPrime = ( Denom(iX1) - Denom(iX1-1) ) / dX
        X0 = XC - DenomC / DenomPrime;
        NumerPrime = NumerC / ( DenomC / DenomPrime )

        FirstTime = .FALSE.

      ELSE IF( ABS( Denom(iX1) ) .LT. Threshold )THEN

        Numer(iX1) = NumerC + NumerPrime + ( X(iX1) - XC )
        Denom(iX1) = DenomC + DenomPrime + ( X(iX1) - XC )

      END IF

    END DO

  END SUBROUTINE IntegrateD


  REAL(DP) FUNCTION Numerator( X, D, U, M )

    REAL(DP), INTENT(in) :: X, D, U, M

    Numerator = D * ( -M / X**2 + Two * U**2 / X + ( Gamma_IDEAL - One ) * U &
                  + ( Gamma_IDEAL - One ) * ( Two - Gamma_IDEAL ) * X )

  END FUNCTION Numerator


  REAL(DP) FUNCTION Denominator( D, U )

    REAL(DP), INTENT(in) :: D, U

    Denominator = Gamma_IDEAL * D**( Gamma_IDEAL - One ) - U**2

  END FUNCTION Denominator


  ! --- End of auxiliary functions for Yahil collapse problem ---


  SUBROUTINE InitializeFields_YahilCollapse_FromScratch &
    ( iLevel, MF_uGF, MF_uCF, &
      dXdr, drhodD, dvdV, dmdM, PolytropicConstant, &
      CoreRadius, D0, CollapseTime, TotalEnclosedMass )

    INTEGER, INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF
    REAL(DP), INTENT(in)  :: dXdr, drhodD, dvdV, dmdM, PolytropicConstant, &
                             CoreRadius, D0, CollapseTime
    REAL(DP), INTENT(out) :: TotalEnclosedMass

    INTEGER               :: N, iX1, iX2, iX3, iX_L, iNX, iNX1
    REAL(DP)              :: dr, dX, XX, R
    REAL(DP), ALLOCATABLE :: X(:), D(:), U(:), V(:), M(:), &
                             Numer(:), Denom(:)

    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uAF_K(nDOFX,nAF)
    INTEGER, PARAMETER :: NX = 2048

    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: iX_B(3), iX_E(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    dr = 1.0e-2_DP * Kilometer
    dX = dXdr * dr

    N = 1.1_DP * CoreRadius * dXdr / dX

    ALLOCATE( Numer(N) )
    ALLOCATE( Denom(N) )
    ALLOCATE( X    (N) )
    ALLOCATE( D    (N) )
    ALLOCATE( U    (N) )
    ALLOCATE( V    (N) )
    ALLOCATE( M    (N) )

    X    (1) = 1.0e-5_DP
    D    (1) = D0
    U    (1) = Zero
    M    (1) = Zero
    Numer(1) = Numerator  ( X(1), D(1), U(1), M(1) )
    Denom(1) = Denominator( D(1), U(1) )

    CALL IntegrateD( dX, X, D, U, M, Numer, Denom )

    TotalEnclosedMass = M(N)

    V = ( Gamma_IDEAL - Two ) * X + U

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      iX_B = BX % lo
      iX_E = BX % hi

      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) ) &
        iX_E(1) = iX_E(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)
      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)

        R = NodeCoordinate( MeshX(1), iX1, iNX1 )
        XX = dXdr * R

        iX_L = Locate( XX, X, N )

        uPF_K(iNX,iPF_D ) &
          = drhodD * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                               D(iX_L), D(iX_L+1) )

        uPF_K(iNX,iPF_V1) &
          = dvdV   * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                               V(iX_L), V(iX_L+1) )

        uPF_K(iNX,iPF_V2) = Zero

        uPF_K(iNX,iPF_V3) = Zero

        uPF_K(iNX,iPF_E ) &
          = PolytropicConstant * uPF_K(iNX,iPF_D)**( Gamma_IDEAL ) &
              / ( Gamma_IDEAL - One )

        uPF_K(iNX,iPF_Ne) = Zero

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF_K(iNX,iPF_D ), uPF_K(iNX,iPF_E), &
                 uPF_K(iNX,iPF_Ne), uAF_K(iNX,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(iNX,iPF_D ), uPF_K(iNX,iPF_V1), &
                 uPF_K(iNX,iPF_V2), uPF_K(iNX,iPF_V3), &
                 uPF_K(iNX,iPF_E ), uPF_K(iNX,iPF_Ne), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                 uAF_K(iNX,iAF_P) )

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE MFI % next()

    DEALLOCATE( M )
    DEALLOCATE( V )
    DEALLOCATE( U )
    DEALLOCATE( D )
    DEALLOCATE( X )
    DEALLOCATE( Denom )
    DEALLOCATE( Numer )

  END SUBROUTINE InitializeFields_YahilCollapse_FromScratch


END MODULE MF_InitializationModule
