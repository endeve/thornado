MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    nNodesX, &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX, &
    Interpolate1D_Linear
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
    iGF_Alpha,    &
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
    iAF_Gm, &
    uDF
  USE Euler_SlopeLimiterModule_Relativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_Relativistic_TABLE
  USE Euler_PositivityLimiterModule_Relativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_Relativistic_TABLE
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic, &
    ComputeFromConserved_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputeThermodynamicStates_Primitive, &
    ApplyEquationOfState
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic


CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
    ( ProgenitorFileName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL ::  ProgenitorFileName_Option

    CHARACTER(LEN=32) :: ProgenitorFileName

    ProgenitorFileName = '../Progenitors/WH07_15M_Sun.h5'
    IF( PRESENT( ProgenitorFileName_Option ) ) &
       ProgenitorFileName = TRIM( ProgenitorFileName_Option )

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'AdiabaticCollapse_XCFC' )

         CALL InitializeFields_AdiabaticCollapse_XCFC &
                ( TRIM( ProgenitorFileName ) )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_AdiabaticCollapse_XCFC &
    ( ProgenitorFileName )

    CHARACTER(LEN=*), INTENT(in) :: ProgenitorFileName

    INTEGER                :: iX1, iX2, iX3
    INTEGER                :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP)               :: X1
    TYPE(ProgenitorType1D) :: P1D

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

    WRITE(*,*)
    WRITE(*,'(6x,A,A)') &
      'ProgenitorFileName: ', TRIM( ProgenitorFileName )
    WRITE(*,*)

    CALL ReadProgenitor1D( TRIM( ProgenitorFileName ), P1D )

    ! --- Initialize Fluid Fields ---

    ASSOCIATE &
      ( R1D => P1D % Radius, &
        D1D => P1D % MassDensity, &
        V1D => P1D % RadialVelocity, &
        T1D => P1D % Temperature, &
        Y1D => P1D % ElectronFraction )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
          = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ApplyEquationOfState &
                ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xa), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xh), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Gm) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ) )

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! R1D, etc

    ! --- Iterate to incorporate gravity in initial conditions ---

    CONVERGED = .FALSE.
    ITER = 0

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

      IF( MAX( dAlpha, dPsi ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. 10 )THEN

        WRITE(*,*) 'Could not initialize fields. Exiting...'
        STOP

      END IF

    END DO

    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
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

  END SUBROUTINE InitializeFields_AdiabaticCollapse_XCFC


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

    END IF

    RETURN

  END FUNCTION Interpolate1D


END MODULE InitializationModule_Relativistic
