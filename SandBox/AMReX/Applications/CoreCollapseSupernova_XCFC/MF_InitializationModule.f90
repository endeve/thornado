MODULE MF_InitializationModule

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_mfiter
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse_build, &
    amrex_parmparse_destroy, &
    amrex_parmparse

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE EquationOfStateModule, ONLY: &
    ComputeThermodynamicStates_Primitive, &
    ApplyEquationOfState
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE UnitsModule, ONLY: &
    UnitsDisplay
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
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nCR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    nPR, &
    nSpecies
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear, &
    Interpolate2D_BiLinear
  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D
  USE InputParsingModule, ONLY: &
    ProgramName, &
    swX, &
    nE, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF, MF_uPF, MF_uAF, MF_uCR, MF_uPR )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uCF, MF_uPF, MF_uAF, MF_uCR, MF_uPR

    ! --- thornado ---

    INTEGER  :: iX1, iX2, iX3, iZ1, iS
    INTEGER  :: iNX, iNX1, iNZ, iNX_Z, iN_E
    REAL(DP) :: X1, Eq
    REAL(DP) :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP) :: uPR_K(nDOFZ,nE,nPR,nSpecies)

    ! --- AMReX ---

    INTEGER                       :: lo_CR(4), hi_CR(4)
    INTEGER                       :: lo_PR(4), hi_PR(4)
    INTEGER                       :: iX_B(3), iX_E(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    REAL(DP) :: R_custom(1:542), D_custom(1:542), Ye_custom(1:542), T_custom(1:542)
    REAL(DP) :: Dnu_custom(1:542,1:20)
    REAL(DP) :: Dnubar_custom(1:542,1:20)
    REAL(DP) :: Inu_custom(1:542,1:20)
    REAL(DP) :: Inubar_custom(1:542,1:20)
    REAL(DP) :: E_custom(1:20)

    ! --- Problem-dependent parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: ProgenitorFileName
    TYPE(ProgenitorType1D)        :: P1D
    LOGICAL                       :: Verbose
    LOGICAL                       :: Custom_Start
 
    CALL amrex_parmparse_build( PP, 'AC' )
      CALL PP % get( 'ProgenitorFileName', ProgenitorFileName )
    CALL amrex_parmparse_destroy( PP )

    Verbose = .FALSE.
    IF( amrex_parallel_ioprocessor() .AND. iLevel .EQ. 0 ) Verbose = .TRUE.

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'Initializing: ', TRIM( ProgramName )
      WRITE(*,'(4x,A)')   '-------------'
      WRITE(*,*)

    END IF

    CALL ReadProgenitor1D &
           ( TRIM( ProgenitorFileName ), P1D, &
             Verbose_Option = Verbose )

    ASSOCIATE &
      ( R1D => P1D % Radius, &
        D1D => P1D % MassDensity, &
        V1D => P1D % RadialVelocity, &
        T1D => P1D % Temperature, &
        Y1D => P1D % ElectronFraction )

    uCR_K = Zero
    uPR_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )


    Custom_Start = .FALSE.

    IF ( Custom_Start ) THEN

      CALL ReadCustom_Fluid( 1, 542, R_custom, D_custom, Ye_custom, T_custom )

      CALL ReadCustom_Transport( 1, 542, 1, 20, &
             Dnu_custom, Dnubar_custom, Inu_custom, Inubar_custom, E_custom )

    END IF


    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uPR => MF_uPR % DataPtr( MFI )
      uPF => MF_uPF % DataPtr( MFI )
      uAF => MF_uAF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_CR = LBOUND( uCR )
      hi_CR = UBOUND( uCR )

      lo_PR = LBOUND( uPR )
      hi_PR = UBOUND( uPR )

      iX_B = BX % lo
      iX_E = BX % hi

      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) ) &
        iX_E(1) = iX_E(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)

        DO iNX = 1, nDOFX

          iNX1 = NodeNumberTableX(1,iNX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX) &
            = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

          uPF(iX1,iX2,iX3,nDOFX*(iPF_V1-1)+iNX) &
            = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

          uPF(iX1,iX2,iX3,nDOFX*(iPF_V2-1)+iNX) &
            = Zero

          uPF(iX1,iX2,iX3,nDOFX*(iPF_V3-1)+iNX) &
            = Zero

          uAF(iX1,iX2,iX3,nDOFX*(iAF_T -1)+iNX) &
            = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

          uAF(iX1,iX2,iX3,nDOFX*(iAF_Ye-1)+iNX) &
            = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )



          IF ( Custom_Start ) THEN

            uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX) &
              = Interpolate1D( R_custom, D_custom, SIZE( R_custom ), X1 )

            uPF(iX1,iX2,iX3,nDOFX*(iPF_V1-1)+iNX) &
              = Zero

            uPF(iX1,iX2,iX3,nDOFX*(iPF_V2-1)+iNX) &
              = Zero

            uPF(iX1,iX2,iX3,nDOFX*(iPF_V3-1)+iNX) &
              = Zero

            uAF(iX1,iX2,iX3,nDOFX*(iAF_T -1)+iNX) &
              = Interpolate1D( R_custom, T_custom, SIZE( R_custom ), X1 )

            uAF(iX1,iX2,iX3,nDOFX*(iAF_Ye-1)+iNX) &
              = Interpolate1D( R_custom, Ye_custom, SIZE( R_custom ), X1 )

          END IF


          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX), &
                   uAF(iX1,iX2,iX3,nDOFX*(iAF_T -1)+iNX), &
                   uAF(iX1,iX2,iX3,nDOFX*(iAF_Ye-1)+iNX), &
                   uPF(iX1,iX2,iX3,nDOFX*(iPF_E -1)+iNX), &
                   uAF(iX1,iX2,iX3,nDOFX*(iAF_E -1)+iNX), &
                   uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX) )

          CALL ApplyEquationOfState &
                  ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_T -1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Ye-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_P -1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_S -1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_E -1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Me-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Mp-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Mn-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Xp-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Xn-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Xa-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Xh-1)+iNX), &
                    uAF(iX1,iX2,iX3,nDOFX*(iAF_Gm-1)+iNX) )

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

          ! --- Initialize radiation fields ---

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            DO iNZ = 1, nDOFZ

              iNX_Z = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

              uPR_K(iNZ,iZ1,iPR_D ,iS) = 1.0e-40_DP
              uPR_K(iNZ,iZ1,iPR_I1,iS) = Zero
              uPR_K(iNZ,iZ1,iPR_I2,iS) = Zero
              uPR_K(iNZ,iZ1,iPR_I3,iS) = Zero


             IF ( Custom_Start ) THEN

               iN_E = MOD( (iNZ-1)        , nDOFE ) + 1

               Eq = NodeCoordinate( MeshE, iZ1, iN_E )
               IF (iS .EQ. 1) THEN

                 uPR_K(iNZ,iZ1,iPR_D ,iS) &
                   = Interpolate2D( R_custom, E_custom, Dnu_custom, &
                                    SIZE(R_custom), SIZE(E_custom), X1, Eq )

                 uPR_K(iNZ,iZ1,iPR_I1,iS) &
                   = Interpolate2D( R_custom, E_custom, Inu_custom, &
                                    SIZE(R_custom), SIZE(E_custom), X1, Eq )

                 uPR_K(iNZ,iZ1,iPR_I2,iS) = Zero

                 uPR_K(iNZ,iZ1,iPR_I3,iS) = Zero

               ELSE IF (iS .EQ. 2) THEN

                 uPR_K(iNZ,iZ1,iPR_D ,iS) &
                   = Interpolate2D( R_custom, E_custom, Dnubar_custom, &
                                    SIZE(R_custom), SIZE(E_custom), X1, Eq )

                 uPR_K(iNZ,iZ1,iPR_I1,iS) &
                   = Interpolate2D( R_custom, E_custom, Inubar_custom, &
                                    SIZE(R_custom), SIZE(E_custom), X1, Eq )

                 uPR_K(iNZ,iZ1,iPR_I2,iS) = Zero

                 uPR_K(iNZ,iZ1,iPR_I3,iS) = Zero

               END IF
             END IF



              CALL ComputeConserved_TwoMoment &
                     ( uPR_K(iNZ,iZ1,iPR_D ,iS), &
                       uPR_K(iNZ,iZ1,iPR_I1,iS), &
                       uPR_K(iNZ,iZ1,iPR_I2,iS), &
                       uPR_K(iNZ,iZ1,iPR_I3,iS), &
                       uCR_K(iNZ,iZ1,iCR_N ,iS), &
                       uCR_K(iNZ,iZ1,iCR_G1,iS), &
                       uCR_K(iNZ,iZ1,iCR_G2,iS), &
                       uCR_K(iNZ,iZ1,iCR_G3,iS), &
                       uPF(iX1,iX2,iX3,nDOFX*(iPF_V1      -1)+iNX_Z), &
                       uPF(iX1,iX2,iX3,nDOFX*(iPF_V2      -1)+iNX_Z), &
                       uPF(iX1,iX2,iX3,nDOFX*(iPF_V3      -1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX_Z), &
                       0.0_DP, 0.0_DP, 0.0_DP, & ! off-diagonal metric comp.
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha   -1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1  -1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2  -1)+iNX_Z), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3  -1)+iNX_Z) )

            END DO ! iNZ = 1, nDOFZ

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNX = 1, nDOFX

        uPR(iX1,iX2,iX3,lo_PR(4):hi_PR(4)) &
          = RESHAPE( uPR_K, [ hi_PR(4) - lo_PR(4) + 1 ] )

        uCR(iX1,iX2,iX3,lo_CR(4):hi_CR(4)) &
          = RESHAPE( uCR_K, [ hi_CR(4) - lo_CR(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO ! MFI % next()

    CALL amrex_mfiter_destroy( MFI )

    END ASSOCIATE ! P1D variables

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    IF( i .EQ. 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

    ELSE IF( i .EQ. n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

    END IF

    RETURN

  END FUNCTION Interpolate1D

  REAL(DP) FUNCTION Interpolate2D( x, E, U, n, n2, xq, Eq )
 
    INTEGER,                   INTENT(in) :: n, n2
    REAL(DP), DIMENSION(n),    INTENT(in) :: x
    REAL(DP), DIMENSION(n2),   INTENT(in) :: E
    REAL(DP), DIMENSION(n,n2), INTENT(in) :: U
    REAL(DP),                  INTENT(in) :: xq, Eq
 
    INTEGER :: i, j
 
    i = Locate( xq, x, n )
    j = Locate( Eq, E, n2 )
 
    IF( i .EQ. 0 )THEN
 
      ! --- Extrapolate Left ---
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(1), x(2), Eq, E(j), E(j+1), U(1,j), U(1,j+1), U(2,j), U(2,j+1) )
 
    ELSE IF( j .EQ. 0 )THEN
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(i), x(i+1), Eq,  E(1), E(2), U(i,1), U(i,2), U(i+1,1), U(i+1,2) )
 
    ELSE IF( i .EQ. 0 .AND. j .EQ. 0 )THEN
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(1), x(2), Eq,  E(1), E(2), U(1,1), U(1,2), U(2,1), U(2,2) )
 
    ELSE IF( i .EQ. n )THEN
 
      ! --- Extrapolate Right ---
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(n-1), x(n), Eq,  E(j), E(j+1), U(n-1,j), U(n-1,j+1), U(n,j), U(n,j+1) )
 
    ELSE IF( j .EQ. n2 )THEN
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(i), x(i+1), Eq,  E(n2-1), E(n2), U(i,n2-1), U(i,n2), U(i+1,n2-1), U(i+1,n2) )
 
    ELSE IF( i .EQ. n .AND. j .EQ. n2 )THEN
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(n-1), x(n), Eq,  E(n2-1), E(n2), U(n-1,n2-1), U(n-1,n2), U(n,n2-1), U(n,n2) )
 
    ELSE
 
      Interpolate2D &
        = Interpolate2D_BiLinear( xq, x(i), x(i+1), Eq,  E(j), E(j+1), U(i,j), U(i,j+1), U(i+1,j), U(i+1,j+1) )
 
    END IF
 
    RETURN
 
  END FUNCTION Interpolate2D
 
  SUBROUTINE ReadCustom_Fluid( iX_B1, iX_E1, R, D, Ye, T )

    INTEGER,  INTENT(in)  :: iX_B1, iX_E1
    REAL(DP), INTENT(out) :: R(iX_B1:), D(iX_B1:), Ye(iX_B1:), T(iX_B1:)

    INTEGER               :: iX1



    OPEN( UNIT = 101, FILE = "./Custom/fluidvariables.dat", status='old' )

    DO iX1 = iX_B1, iX_E1
      READ(101, *) R(iX1), D(iX1), Ye(iX1), T(iX1)

      R(iX1) = R(iX1) * UnitsDisplay % LengthX1Unit
      D(iX1) = D(iX1) * UnitsDisplay % MassDensityUnit
      T(iX1) = T(iX1) * UnitsDisplay % TemperatureUnit

    END DO

    CLOSE( 101 )

  END SUBROUTINE


  SUBROUTINE ReadCustom_Transport( iX_B1, iX_E1, iE_B1, iE_E1, Dnu, Dnubar, Inu, Inubar, E )

    INTEGER,  INTENT(in)  :: iX_B1, iX_E1, iE_B1, iE_E1
    REAL(DP), INTENT(out) :: Dnu(iX_B1:, iE_B1:)
    REAL(DP), INTENT(out) :: Dnubar(iX_B1:, iE_B1:)
    REAL(DP), INTENT(out) :: Inu(iX_B1:, iE_B1:)
    REAL(DP), INTENT(out) :: Inubar(iX_B1:, iE_B1:)
    REAL(DP), INTENT(out) :: E(iE_B1:)

    INTEGER               :: iX1, iE



    OPEN( UNIT = 101, FILE = "./Custom/NeutrinoD.dat", status='old' )
    OPEN( UNIT = 102, FILE = "./Custom/NeutrinoDbar.dat", status='old' )
    OPEN( UNIT = 103, FILE = "./Custom/NeutrinoI.dat", status='old' )
    OPEN( UNIT = 104, FILE = "./Custom/NeutrinoIbar.dat", status='old' )
    OPEN( UNIT = 105, FILE = "./Custom/E.dat", status='old' )

    DO iX1 = iX_B1, iX_E1

       READ(101, *) Dnu(iX1,:)
       READ(102, *) Dnubar(iX1,:)
       READ(103, *) Inu(iX1,:)
       READ(104, *) Inubar(iX1,:)

    END DO


    DO iE = iE_B1, iE_E1

       READ(105, *) E(iE)

      E(iE) = E(iE) * UnitsDisplay % EnergyUnit
    END DO

    CLOSE( 101 )
    CLOSE( 102 )
    CLOSE( 103 )
    CLOSE( 104 )
    CLOSE( 105 )


  END SUBROUTINE




END MODULE MF_InitializationModule
