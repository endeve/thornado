MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
!!$  USE GravitySolutionModule_XCFC_Poseidon, ONLY: &
!!$    SolveGravity_XCFC_Poseidon
  USE GravitySolutionModule_XCFC_Poseidon, ONLY: &
    ComputeConformalFactor_Poseidon, &
    ComputeGeometry_Poseidon
  USE Poseidon_UtilitiesModule, ONLY: &
    MultiplyByPsi6, &
    DivideByPsi6, &
    ComputeMatterSources_Poseidon, &
    ComputePressureTensorTrace_Poseidon
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    uPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    uAF, &
    iAF_P, &
    uDF
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Kilometer, &
    SolarMass, &
    Gram, &
    Centimeter, &
    Erg, &
    Millisecond, &
    GravitationalConstant
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE Euler_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    ApplySlopeLimiter_Euler_Relativistic_IDEAL
  USE Euler_PositivityLimiterModule_Relativistic_IDEAL, ONLY: &
    ApplyPositivityLimiter_Euler_Relativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic

CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
    ( ReadFromFile_Option, FileName_Option, &
      D0_Option, CentralDensity_Option, CentralPressure_Option, &
      CoreRadius_Option, CollapseTime_Option )

    LOGICAL,           INTENT(in), OPTIONAL :: ReadFromFile_Option
    CHARACTER(LEN=64), INTENT(in), OPTIONAL :: FileName_Option
    REAL(DP),          INTENT(in), OPTIONAL :: D0_Option
    REAL(DP),          INTENT(in), OPTIONAL :: CentralDensity_Option
    REAL(DP),          INTENT(in), OPTIONAL :: CentralPressure_Option
    REAL(DP),          INTENT(in), OPTIONAL :: CoreRadius_Option
    REAL(DP),          INTENT(in), OPTIONAL :: CollapseTime_Option

    LOGICAL           :: ReadFromFile    = .TRUE.
    CHARACTER(LEN=64) :: &
      FileName = 'YahilHomologousCollapse_Gm1.30_t1.500E+002ms.dat'
    REAL(DP)          :: D0              = 1.75_DP
    REAL(DP)          :: CentralDensity  = 7.0e9_DP  * ( Gram / Centimeter**3 )
    REAL(DP)          :: CentralPressure = 6.0e27_DP * ( Erg  / Centimeter**3 )
    REAL(DP)          :: CoreRadius      = 1.0e5_DP  * Kilometer
    REAL(DP)          :: CollapseTime    = 1.50e2_DP * Millisecond

    uPF(:,:,:,:,iPF_Ne) = Zero

    IF( PRESENT( ReadFromFile_Option ) ) &
      ReadFromFile = ReadFromFile_Option
    IF( PRESENT( FileName_Option ) ) &
      FileName = TRIM( FileName_Option )
    IF( PRESENT( D0_Option ) ) &
      D0 = D0_Option
    IF( PRESENT( CentralDensity_Option ) ) &
      CentralDensity = CentralDensity_Option
    IF( PRESENT( CentralPressure_Option ) ) &
      CentralPressure = CentralPressure_Option
    IF( PRESENT( CoreRadius_Option ) ) &
      CoreRadius = CoreRadius_Option
    IF( PRESENT( CollapseTime_Option ) ) &
      CollapseTime = CollapseTime_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'YahilCollapse' )

         CALL InitializeFields_YahilCollapse &
                ( ReadFromFile, FileName, &
                  D0, CentralDensity, CentralPressure, &
                  CoreRadius, CollapseTime )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_YahilCollapse &
    ( ReadFromFile, FileName, D0, CentralDensity, CentralPressure, &
      CoreRadius, CollapseTime )

    LOGICAL,           INTENT(in) :: ReadFromFile
    CHARACTER(LEN=64), INTENT(in) :: FileName
    REAL(DP),          INTENT(in) :: D0
    REAL(DP),          INTENT(in) :: CentralDensity
    REAL(DP),          INTENT(in) :: CentralPressure
    REAL(DP),          INTENT(in) :: CoreRadius
    REAL(DP),          INTENT(in) :: CollapseTime

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: PolytropicConstant, dXdr, drhodD, dvdV, dmdM, TotalEnclosedMass

    REAL(DP) :: E (nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: Si(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),3)
    REAL(DP) :: S (nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: Mg(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))

    INTEGER  :: ITER
    REAL(DP) :: dAlpha, dPsi
    LOGICAL  :: CONVERGED

    REAL(DP) :: dAl1(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dCF1(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dAl2(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dCF2(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))

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

    IF( ReadFromFile )THEN

      CALL InitializeFields_YahilCollapse_FromFile &
             ( FileName, dXdr, drhodD, dvdV, &
               PolytropicConstant, TotalEnclosedMass )

    ELSE

      CALL InitializeFields_YahilCollapse_FromScratch &
             ( dXdr, drhodD, dvdV, dmdM, &
               PolytropicConstant, CoreRadius, D0, CollapseTime, &
               TotalEnclosedMass )

    END IF

    ! --- Iterate to incorporate gravity in initial conditions ---

    CONVERGED = .FALSE.
    ITER      = 0

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      dAl1 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Alpha)
      dCF1 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Psi  )

      CALL MultiplyByPsi6( iX_B1, iX_E1, uGF, uCF )

      CALL ComputeMatterSources_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, E, Si, Mg )

      CALL ComputeConformalFactor_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, uGF )

      CALL ComputePressureTensorTrace_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, S )

      CALL ComputeGeometry_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, uGF )

      dAl2 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Alpha)
      dCF2 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Psi  )

      dAlpha = MINVAL( ABS( dAl2 - dAl1 ) / ( Half * ( dAl1 + dAl2 ) ) )
      dPsi   = MINVAL( ABS( dCF2 - dCF1 ) / ( Half * ( dCF1 + dCF2 ) ) )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E1(1)

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                 uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                 uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                 uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                 uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                 uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(:,iX1,iX2,iX3,iAF_P) )

      END DO
      END DO
      END DO

      IF( MAX( dAlpha, dPsi ) .LT. 1.0e-12_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. 20 )THEN

        WRITE(*,*) 'Could not initialize fields. Exiting...'
        STOP

      END IF

    END DO

    CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

    CALL MultiplyByPsi6( iX_B1, iX_E1, uGF, uCF )

    CALL ComputeMatterSources_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, E, Si, Mg )

    CALL ComputeConformalFactor_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, uGF )

    CALL ComputePressureTensorTrace_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, S )

    CALL ComputeGeometry_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, uGF )

    CALL DivideByPsi6( iX_B1, iX_E1, uGF, uCF )

    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'ReadFromFile:        ', ReadFromFile
    WRITE(*,'(6x,A,F5.3)') &
      'Adiabatic Gamma:     ', &
      Gamma_IDEAL
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Polytropic Constant: ', &
      PolytropicConstant / ( Erg / Centimeter**3 &
              / ( Gram / Centimeter**3  )**( Gamma_IDEAL ) ), &
      ' ( erg/cm^3 ) / ( g/cm^3 )^( Gamma )'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Core Radius:         ', &
      CoreRadius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Collapse Time:       ', &
      CollapseTime / Millisecond, ' ms'
    WRITE(*,'(6x,A,ES10.3E3,A)') &
      'Mass:                ', &
      TotalEnclosedMass * dmdM / SolarMass, ' Msun'
    WRITE(*,*)

  END SUBROUTINE InitializeFields_YahilCollapse


  SUBROUTINE InitializeFields_YahilCollapse_FromFile &
    ( FileName, dXdr, drhodD, dvdV, PolytropicConstant, TotalEnclosedMass )

    CHARACTER(LEN=64), INTENT(in)  :: FileName
    REAL(DP),          INTENT(in)  :: dXdr, drhodD, dvdV, PolytropicConstant
    REAL(DP),          INTENT(out) :: TotalEnclosedMass

    INTEGER               :: nLines
    INTEGER               :: iX1, iX2, iX3, iNodeX, iNodeX1, iX_L
    REAL(DP)              :: R, XX
    REAL(DP), ALLOCATABLE :: X(:), D(:), V(:), M(:)

    ! --- https://stackoverflow.com/questions/30692424/
    !     how-to-read-number-of-lines-in-fortran-90-from-a-text-file
    nLines = 0
    OPEN(100,FILE=TRIM(FileName))
    READ(100,*)
    DO
      READ(100,*,END=10)
      nLines = nLines + 1
    END DO
    10 CLOSE(100)

    ALLOCATE( X(nLines) )
    ALLOCATE( D(nLines) )
    ALLOCATE( V(nLines) )
    ALLOCATE( M(nLines) )

    OPEN(100,FILE=TRIM(FileName))
    READ(100,*)

    DO iX1 = 1, nLines

      READ(100,*) X(iX1), D(iX1), V(iX1), M(iX1)

    END DO

    CLOSE(100)

    TotalEnclosedMass = M(nLines)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

     DO iNodeX = 1, nDOFX

       iNodeX1 = NodeNumberTableX(1,iNodeX)

       R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
       XX = dXdr * R

       iX_L = Locate( XX, X, nLines )

       uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
         = drhodD * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                              D(iX_L), D(iX_L+1) )

       uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
         = dvdV   * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                              V(iX_L), V(iX_L+1) )

       uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero

       uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

       uPF(iNodeX,iX1,iX2,iX3,iPF_E ) &
         = PolytropicConstant * uPF(iNodeX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) &
             / ( Gamma_IDEAL - One )

     END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    DEALLOCATE( M )
    DEALLOCATE( V )
    DEALLOCATE( D )
    DEALLOCATE( X )

  END SUBROUTINE InitializeFields_YahilCollapse_FromFile


  SUBROUTINE InitializeFields_YahilCollapse_FromScratch &
    ( dXdr, drhodD, dvdV, dmdM, PolytropicConstant, &
      CoreRadius, D0, CollapseTime, TotalEnclosedMass )

    REAL(DP), INTENT(in)  :: dXdr, drhodD, dvdV, dmdM, PolytropicConstant, &
                             CoreRadius, D0, CollapseTime
    REAL(DP), INTENT(out) :: TotalEnclosedMass

    INTEGER               :: N, iX1, iX2, iX3, iNodeX, iNodeX1, iX_L
    REAL(DP)              :: dr, dX, XX, R
    REAL(DP), ALLOCATABLE :: X(:), D(:), U(:), V(:), M(:), &
                             Numer(:), Denom(:)

    LOGICAL, PARAMETER    :: WriteToFile = .FALSE.
    CHARACTER(LEN=64)     :: FileName
    INTEGER               :: nLines, iLine
    REAL(DP), ALLOCATABLE :: FileX(:), FileD(:), FileV(:), FileM(:)

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

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

     DO iNodeX = 1, nDOFX

       iNodeX1 = NodeNumberTableX(1,iNodeX)

       R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
       XX = dXdr * R

       iX_L = Locate( XX, X, N )

       uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
         = drhodD * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                              D(iX_L), D(iX_L+1) )

       uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
         = dvdV * Interpolate1D_Linear( XX, X(iX_L), X(iX_L+1), &
                                            V(iX_L), V(iX_L+1) )

       uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero

       uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

       uPF(iNodeX,iX1,iX2,iX3,iPF_E ) &
         = PolytropicConstant * uPF(iNodeX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) &
             / ( Gamma_IDEAL - One )

     END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    IF( WriteToFile )THEN

      nLines = N

      ALLOCATE( FileX(nLines) )
      ALLOCATE( FileD(nLines) )
      ALLOCATE( FileV(nLines) )
      ALLOCATE( FileM(nLines) )

      FileX(1)      = X(1)
      FileD(1)      = D(1)
      FileV(1)      = V(1)
      FileM(1)      = M(1)
      FileX(nLines) = X(N)
      FileD(nLines) = D(N)
      FileV(nLines) = V(N)
      FileM(nLines) = M(N)

      IF( nLines .NE. N ) dX = ( X(N) - X(1) ) / DBLE( nLines - 1 )

      DO iLine = 2, nLines-1

        FileX(iLine) = FileX(iLine-1) + dX

        iX_L = Locate( FileX(iLine), X, N )

        FileD(iLine) &
          = Interpolate1D_Linear( FileX(iLine), X(iX_L), X(iX_L+1), &
                                                D(iX_L), D(iX_L+1) )

        FileV(iLine) &
          = Interpolate1D_Linear( FileX(iLine), X(iX_L), X(iX_L+1), &
                                                V(iX_L), V(iX_L+1) )

        FileM(iLine) &
          = Interpolate1D_Linear( FileX(iLine), X(iX_L), X(iX_L+1), &
                                                M(iX_L), M(iX_L+1) )

      END DO

      WRITE( FileName, '(A,F4.2,A,ES10.3E3,A)' ) &
        'YahilHomologousCollapse_Gm', Gamma_IDEAL, '_t', &
                 CollapseTime / Millisecond, 'ms.dat'

      OPEN( 100, FILE = TRIM( FileName ) )

      WRITE( 100, '(A)' ) '# X D V M'

      DO iLine = 1, nLines

        WRITE( 100, '(ES24.16E3,1x,ES24.16E3,1x,ES24.16E3,1x,ES24.16E3)' ) &
          FileX(iLine), FileD(iLine), FileV(iLine), FileM(iLine)

      END DO

      CLOSE( 100 )

      DEALLOCATE( FileM )
      DEALLOCATE( FileV )
      DEALLOCATE( FileD )
      DEALLOCATE( FileX )

    END IF ! WriteToFile

    DEALLOCATE( M )
    DEALLOCATE( V )
    DEALLOCATE( U )
    DEALLOCATE( D )
    DEALLOCATE( X )
    DEALLOCATE( Denom )
    DEALLOCATE( Numer )

  END SUBROUTINE InitializeFields_YahilCollapse_FromScratch


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

END MODULE InitializationModule_Relativistic
