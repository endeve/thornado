MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_imultifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDimsX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshType
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    Centimeter, &
    Second

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    FourPi
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab
  USE InputParsingModule, ONLY: &
    nLevels, &
    t_new
  USE FillPatchModule, ONLY: &
    FillPatch
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteNodal1DICToFile_SAS
  PUBLIC :: ComputePowerInLegendreModes

  LOGICAL,          PUBLIC              :: WriteNodal1DIC_SAS
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: FileName_Nodal1DIC_SAS

  INTEGER, PARAMETER :: nLeg = 5 ! number of Legendre polynomials

  LOGICAL :: ComputePowerInSitu

CONTAINS


  SUBROUTINE WriteNodal1DICToFile_SAS

    TYPE(amrex_parmparse) :: PP

    INTEGER :: iLevel

    CHARACTER(256) :: FileName

    IF( nDimsX .GT. 1 ) RETURN

    DO iLevel = 0, nLevels - 1

      CALL FillPatch( iLevel, MF_uGF )
      CALL FillPatch( iLevel, MF_uGF, MF_uCF )

    END DO

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
             swXX_Option = swX )

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'FileName_Nodal1DIC_SAS', &
                      FileName_Nodal1DIC_SAS )
    CALL amrex_parmparse_destroy( PP )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_D'
    CALL ShowVariableFromMultiFab &
      ( MF_uPF, iPF_D, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_V1'
    CALL ShowVariableFromMultiFab &
      ( MF_uPF, iPF_V1, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_AF_P'
    CALL ShowVariableFromMultiFab &
      ( MF_uAF, iAF_P, swXX_Option = swX, UseFineMask_Option = .FALSE., &
        WriteToFile_Option = .TRUE., &
        FileNameBase_Option = TRIM( FileName ) )

    IF( amrex_parallel_ioprocessor() )THEN

      OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) // '_BC.dat' )

      WRITE(101,'(ES24.16E3)') ExpD
      WRITE(101,'(ES24.16E3)') ExpE

    END IF

  END SUBROUTINE WriteNodal1DICToFile_SAS


  SUBROUTINE ComputePowerInLegendreModes( MF_uGF, MF_uPF )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uPF(0:nLevels-1)

    TYPE(amrex_multifab) :: PowerDensity      (0:nLevels-1)
    TYPE(amrex_multifab) :: RadialPowerDensity(0:nLevels-1)
    TYPE(amrex_multifab) :: Psi_K             (0:nLevels-1)

    INTEGER :: iLevel

    REAL(DP) :: Power(0:nLeg-1)

    TYPE(amrex_parmparse) :: PP

    ComputePowerInSitu = .FALSE.
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'ComputePower', ComputePowerInSitu )
    CALL amrex_parmparse_destroy( PP )

    IF( .NOT. ComputePowerInSitu ) RETURN

    Power = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( PowerDensity(iLevel), &
               MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, 1, swX )
      CALL PowerDensity      (iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( RadialPowerDensity(iLevel), &
               MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, nLeg, 0 )
      CALL RadialPowerDensity(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( Psi_K(iLevel), &
               MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, 1, 0 )
      CALL Psi_K(iLevel) % SetVal( Zero )

    END DO

    CALL ComputePowerDensity( MF_uGF, MF_uPF, PowerDensity, Psi_K )

    CALL ComputeRadialPowerDensity( PowerDensity, RadialPowerDensity )

    CALL ComputePower( Psi_K, RadialPowerDensity, Power )

    CALL WriteToFile( Power )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( Psi_K             (iLevel) )
      CALL amrex_multifab_destroy( RadialPowerDensity(iLevel) )
      CALL amrex_multifab_destroy( PowerDensity      (iLevel) )

    END DO

  END SUBROUTINE ComputePowerInLegendreModes


  SUBROUTINE ComputePowerDensity( MF_uGF, MF_uPF, PowerDensity, Psi_K )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF      (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF      (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: PowerDensity(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: Psi_K       (0:nLevels-1)

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_imultifab) :: iMF_FineMask

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX1, iX2, iX3
    INTEGER :: iLevel, iLo_G, iHi_G, iLo_F, iHi_F, iLo_P, iHi_P

    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: Psi     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: PD      (:,:,:,:)

    REAL(DP) :: V2_K(2), VolK(2), X2(3)

    TYPE(MeshType) :: MeshX(3)

    iLo_G = 1+nDOFX*(iGF_SqrtGm-1)
    iHi_G = nDOFX*iGF_SqrtGm

    iLo_F = 1+nDOFX*(iPF_V2-1)
    iHi_F = nDOFX*iPF_V2

    iLo_P = 1+nDOFX*(iGF_Psi-1)
    iHi_P = nDOFX*iGF_Psi

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel) )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask         % DataPtr( MFI )
        uGF      => MF_uGF      (iLevel) % DataPtr( MFI )
        Psi      => Psi_K       (iLevel) % DataPtr( MFI )
        uPF      => MF_uPF      (iLevel) % DataPtr( MFI )
        PD       => PowerDensity(iLevel) % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          VolK(1) = SUM( WeightsX_q * uGF(iX1,iX2-1,iX3,iLo_G:iHi_G) )
          VolK(2) = SUM( WeightsX_q * uGF(iX1,iX2+1,iX3,iLo_G:iHi_G) )

          V2_K(1) = SUM( WeightsX_q * uGF(iX1,iX2-1,iX3,iLo_G:iHi_G) &
                                    * uPF(iX1,iX2-1,iX3,iLo_F:iHi_F) ) / VolK(1)
          V2_K(2) = SUM( WeightsX_q * uGF(iX1,iX2+1,iX3,iLo_G:iHi_G) &
                                    * uPF(iX1,iX2+1,iX3,iLo_F:iHi_F) ) / VolK(2)

          X2(1) = MeshX(2) % Center(iX2-1)
          X2(2) = MeshX(2) % Center(iX2+1)
          X2(3) = MeshX(2) % Center(iX2  )

          PD(iX1,iX2,iX3,1) &
            = One / SIN( X2(3) ) &
                * ( V2_K(2) * SIN( X2(2) ) - V2_K(1) * SIN( X2(1) ) ) &
                / ( X2(2) - X2(1) )

          Psi(iX1,iX2,iX3,1) &
            = SUM( WeightsX_q * uGF(iX1,iX2,iX3,iLo_G:iHi_G) &
                              * uGF(iX1,iX2,iX3,iLo_P:iHi_P) ) &
               / SUM( WeightsX_q * uGF(iX1,iX2,iX3,iLo_G:iHi_G) )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel

  END SUBROUTINE ComputePowerDensity


  SUBROUTINE ComputeRadialPowerDensity( PowerDensity, RadialPowerDensity )

    TYPE(amrex_multifab), INTENT(in)    :: PowerDensity      (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: RadialPowerDensity(0:nLevels-1)

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_imultifab) :: iMF_FineMask

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX1, iX2, ell
    INTEGER :: iLevel

    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: PD      (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: RPD     (:,:,:,:)

    REAL(DP) :: X2, dX2

    TYPE(MeshType) :: MeshX(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask &
             ( iLevel, iMF_FineMask, PowerDensity % BA, PowerDensity % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, PowerDensity(iLevel) )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask               % DataPtr( MFI )
        PD       => PowerDensity      (iLevel) % DataPtr( MFI )
        RPD      => RadialPowerDensity(iLevel) % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        DO iX1 = iX_B0(1), iX_E0(1)

          DO iX2 = iX_B0(2), iX_E0(2)

            IF( IsNotLeafElement( FineMask(iX1,iX2,iX_B0(3),1) ) ) CYCLE

            X2  = MeshX(2) % Center(iX2)
            dX2 = MeshX(2) % Width (iX2)

            DO ell = 0, nLeg-1

              RPD(iX1,iX_B0(2),iX_B0(3),ell+1) &
                = RPD(iX1,iX_B0(2),iX_B0(3),ell+1) &
                    + PD(iX1,iX2,iX_B0(3),1) &
                        * LegendrePolynomial( ell, COS( X2 ) ) &
                        * SIN( X2 ) * dX2

            END DO ! ell

          END DO ! iX2

        END DO ! iX1

      END DO ! WHILE

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel

  END SUBROUTINE ComputeRadialPowerDensity


  SUBROUTINE ComputePower( Psi_K, RadialPowerDensity, Power )

    TYPE(amrex_multifab), INTENT(in)    :: Psi_K             (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: RadialPowerDensity(0:nLevels-1)
    REAL(DP)            , INTENT(inout) :: Power(0:nLeg-1)

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_imultifab) :: iMF_FineMask

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX1, ell
    INTEGER :: iLevel

    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: Psi     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: RPD     (:,:,:,:)

    REAL(DP) :: X1, dX1, ra, rb, ShockRadius

    TYPE(MeshType) :: MeshX(3)

    TYPE(amrex_parmparse) :: PP

    ra = 0.8_DP
    rb = 0.9_DP
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'ShellBottom', ra )
      CALL PP % query( 'ShellTop'   , rb )
      CALL PP % get  ( 'ShockRadius', ShockRadius )
    CALL amrex_parmparse_destroy( PP )

    ra = ra * ShockRadius * Kilometer
    rb = rb * ShockRadius * Kilometer

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask &
             ( iLevel, iMF_FineMask, Psi_K % BA, Psi_K % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, Psi_K(iLevel) )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask               % DataPtr( MFI )
        Psi      => Psi_K             (iLevel) % DataPtr( MFI )
        RPD      => RadialPowerDensity(iLevel) % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        DO iX1 = iX_B0(1), iX_E0(1)

          IF( IsNotLeafElement( FineMask(iX1,iX_B0(2),iX_B0(3),1) ) ) CYCLE

            X1  = MeshX(1) % Center(iX1)
            dX1 = MeshX(1) % Width (iX1)

            IF( X1 .GT. ra .AND. X1 .LT. rb )THEN

              DO ell = 0, nLeg-1

                Power(ell) &
                  = Power(ell) &
                      + RPD(iX1,iX_B0(2),iX_B0(3),ell+1)**2 &
                          * Psi(iX1,iX_B0(2),iX_B0(3),1)**6 &
                          * X1**2 * dX1

              END DO ! ell

          END IF

        END DO ! iX1

      END DO ! WHILE

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel

    Power = FourPi * Power

    DO ell = 0, nLeg-1

      CALL amrex_parallel_reduce_sum( Power(ell) )

    END DO

  END SUBROUTINE ComputePower


  SUBROUTINE WriteToFile( Power )

    REAL(DP), INTENT(in) :: Power(0:nLeg-1)

    TYPE(amrex_parmparse) :: PP

    CHARACTER(256) :: FMT
    CHARACTER(:), ALLOCATABLE :: PowerFileName

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'PowerFileName', PowerFileName )
    CALL amrex_parmparse_destroy( PP )

    WRITE(FMT,'(A3,I2.2,A10)') '(SP', nLeg+1, 'ES25.16E3)'

    IF( amrex_parallel_ioprocessor() )THEN

      OPEN( 100, FILE = TRIM( PowerFileName ), POSITION = 'APPEND' )

      WRITE(100,TRIM(FMT)) &
        t_new(0) / Millisecond, Power / ( Centimeter**3 / Second**2 )

      CLOSE( 100 )

    END IF

  END SUBROUTINE WriteToFile


  REAL(DP) FUNCTION LegendrePolynomial( ell, x ) RESULT( Pell )

    INTEGER , INTENT(in) :: ell
    REAL(DP), INTENT(in) :: x

    SELECT CASE( ell )

    CASE( 0 )

      Pell = SQRT( 1.0_DP / 2.0_DP )

    CASE( 1 )

      Pell = SQRT( 3.0_DP / 2.0_DP ) * x

    CASE( 2 )

      Pell = SQRT( 5.0_DP / 2.0_DP ) * ( 3.0_DP * x**2 - 1.0_DP ) / 2.0_DP

    CASE( 3 )

      Pell = SQRT( 7.0_DP / 2.0_DP ) &
               * 1.0_DP / 2.0_DP * ( 5.0_DP * x**3 - 3.0_DP * x )

    CASE( 4 )

      Pell = SQRT( 9.0_DP / 2.0_DP ) &
               * 1.0_DP / 8.0_DP * ( 35.0_DP * x**4 - 30.0_DP * x**2 + 3.0_DP )

    END SELECT

    RETURN
  END FUNCTION LegendrePolynomial


END MODULE MF_AccretionShockUtilitiesModule
