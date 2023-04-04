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
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE EquationOfStateModule, ONLY: &
    ComputeThermodynamicStates_Primitive, &
    ApplyEquationOfState
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
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
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nPR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D
  USE InputParsingModule, ONLY: &
    ProgramName, &
    nLevels, &
    xL, &
    xR, &
    nX, &
    swX, &
    nSpecies, &
    nE, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF, MF_uCR )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF, MF_uCF, MF_uCR

    ! --- thornado ---

    INTEGER  :: iX1, iX2, iX3, iZ1, iS
    INTEGER  :: iNX, iNX1, iNZ, iNX_Z
    REAL(DP) :: X1
    REAL(DP) :: uGF_K(nDOFX,nGF)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uAF_K(nDOFX,nAF)
    REAL(DP) :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP) :: uPR_K(nDOFZ,nE,nPR,nSpecies)

    ! --- AMReX ---

    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    INTEGER                       :: lo_R(4), hi_R(4)
    INTEGER                       :: iX_B(3), iX_E(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(amrex_parmparse)         :: PP
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)

    ! --- Problem-dependent parameters ---

    CHARACTER(LEN=:), ALLOCATABLE :: ProgenitorFileName
    TYPE(ProgenitorType1D)        :: P1D
    LOGICAL                       :: Verbose

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

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero
    uCR_K = Zero
    uPR_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      lo_R = LBOUND( uCR )
      hi_R = UBOUND( uCR )

      iX_B = BX % lo
      iX_E = BX % hi

      IF( BX % lo(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) ) &
        iX_B(1) = iX_B(1) - swX(1)

      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) ) &
        iX_E(1) = iX_E(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNX = 1, nDOFX

          iNX1 = NodeNumberTableX(1,iNX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          uPF_K(iNX,iPF_D) &
            = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

          uPF_K(iNX,iPF_V1) &
            = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

          uPF_K(iNX,iPF_V2) &
            = Zero

          uPF_K(iNX,iPF_V3) &
            = Zero

          uAF_K(iNX,iAF_T) &
            = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

          uAF_K(iNX,iAF_Ye) &
            = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF_K(iNX,iPF_D ), &
                   uAF_K(iNX,iAF_T ), &
                   uAF_K(iNX,iAF_Ye), &
                   uPF_K(iNX,iPF_E ), &
                   uAF_K(iNX,iAF_E ), &
                   uPF_K(iNX,iPF_Ne) )

          CALL ApplyEquationOfState &
                  ( uPF_K(iNX,iPF_D ), &
                    uAF_K(iNX,iAF_T ), &
                    uAF_K(iNX,iAF_Ye), &
                    uAF_K(iNX,iAF_P ), &
                    uAF_K(iNX,iAF_S ), &
                    uAF_K(iNX,iAF_E ), &
                    uAF_K(iNX,iAF_Me), &
                    uAF_K(iNX,iAF_Mp), &
                    uAF_K(iNX,iAF_Mn), &
                    uAF_K(iNX,iAF_Xp), &
                    uAF_K(iNX,iAF_Xn), &
                    uAF_K(iNX,iAF_Xa), &
                    uAF_K(iNX,iAF_Xh), &
                    uAF_K(iNX,iAF_Gm) )

          CALL ComputeConserved_Euler_Relativistic &
                 ( uPF_K(iNX,iPF_D ), &
                   uPF_K(iNX,iPF_V1), &
                   uPF_K(iNX,iPF_V2), &
                   uPF_K(iNX,iPF_V3), &
                   uPF_K(iNX,iPF_E ), &
                   uPF_K(iNX,iPF_Ne), &
                   uCF_K(iNX,iCF_D ), &
                   uCF_K(iNX,iCF_S1), &
                   uCF_K(iNX,iCF_S2), &
                   uCF_K(iNX,iCF_S3), &
                   uCF_K(iNX,iCF_E ), &
                   uCF_K(iNX,iCF_Ne), &
                   uGF_K(iNX,iGF_Gm_dd_11), &
                   uGF_K(iNX,iGF_Gm_dd_22), &
                   uGF_K(iNX,iGF_Gm_dd_33), &
                   uAF_K(iNX,iAF_P) )

          ! --- Initialize radiation fields ---

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            DO iNZ = 1, nDOFZ

              iNX_Z = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

              uPR_K(iNZ,iZ1,iPR_D ,iS) = 1.0e-40_DP
              uPR_K(iNZ,iZ1,iPR_I1,iS) = Zero
              uPR_K(iNZ,iZ1,iPR_I2,iS) = Zero
              uPR_K(iNZ,iZ1,iPR_I3,iS) = Zero

              CALL ComputeConserved_TwoMoment &
                     ( uPR_K(iNZ,iZ1,iPR_D ,iS),  &
                       uPR_K(iNZ,iZ1,iPR_I1,iS),  &
                       uPR_K(iNZ,iZ1,iPR_I2,iS),  &
                       uPR_K(iNZ,iZ1,iPR_I3,iS),  &
                       uCR_K(iNZ,iZ1,iCR_N ,iS),  &
                       uCR_K(iNZ,iZ1,iCR_G1,iS),  &
                       uCR_K(iNZ,iZ1,iCR_G2,iS),  &
                       uCR_K(iNZ,iZ1,iCR_G3,iS),  &
                       uPF_K(iNX_Z,iPF_V1),       &
                       uPF_K(iNX_Z,iPF_V2),       &
                       uPF_K(iNX_Z,iPF_V3),       &
                       uGF_K(iNX_Z,iGF_Gm_dd_11), &
                       uGF_K(iNX_Z,iGF_Gm_dd_22), &
                       uGF_K(iNX_Z,iGF_Gm_dd_33), &
                       0.0_DP, 0.0_DP, 0.0_DP,    & ! off-diagonal metric comp.
                       uGF_K(iNX_Z,iGF_Alpha),    &
                       uGF_K(iNX_Z,iGF_Beta_1),   &
                       uGF_K(iNX_Z,iGF_Beta_2),   &
                       uGF_K(iNX_Z,iGF_Beta_3) )


            END DO ! iNZ = 1, nDOFZ

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNX = 1, nDOFX

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        uCR(iX1,iX2,iX3,lo_R(4):hi_R(4)) &
          = RESHAPE( uCR_K, [ hi_R(4) - lo_R(4) + 1 ] )

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


END MODULE MF_InitializationModule
