MODULE MF_InitializationModule

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    CoordinateSystem
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
  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Two, &
    Three, &
    FourPi
  USE InputParsingModule, ONLY: &
    UseTiling, &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

  REAL(DP), PUBLIC :: IntE_Min_Euler_PL
  REAL(DP) :: X_D, Eblast, eIntBlast, rho0
  INTEGER  :: nDetCells

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

    INTEGER :: iX_B0(3), iX_E0(3), iX1, iX2, iX3, iNX

    CHARACTER(LEN=32) :: FMT

    ! --- Problem-Dependent Parameters ---

    INTEGER  :: iNX1, iNX2, iNX3, nDetCells
    REAL(DP) :: X1, X2, X3, Radius, dX1, dX2, dX3

    Verbose = .FALSE.
    IF( amrex_parallel_ioprocessor() .AND. iLevel .EQ. 0 ) Verbose = .TRUE.

    IF( iLevel .EQ. 0 )THEN

      nDetCells = -HUGE( 1 )
      X_D       = -HUGE( 1.0_DP )
      CALL amrex_parmparse_build( PP, 'Sedov' )
        CALL PP % get  ( 'Eblast'   , Eblast )
        CALL PP % query( 'nDetCells', nDetCells )
        CALL PP % query( 'X_D'      , X_D )
      CALL amrex_parmparse_destroy( PP )

      IF( X_D .LT. Zero )THEN

        IF( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

          dX1 = MeshX(1) % Width(0)

          X_D = dX1 * DBLE( nDetCells )

        ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

          dX1 = MeshX(1) % Width(0)
          dX2 = MeshX(2) % Width(0)

          X_D = SQRT( dX1**2 + dX2**2 ) * DBLE( nDetCells )

        ELSE

          dX1 = MeshX(1) % Width(0)
          dX2 = MeshX(2) % Width(0)
          dX3 = MeshX(3) % Width(0)

          X_D = SQRT( dX1**2 + dX2**2 + dX3**2 ) * DBLE( nDetCells )

        END IF

      END IF

      eIntBlast = Eblast / ( FourPi / Three * ( X_D / Two**( nMaxLevels-1 ) )**3 )

      ! --- Enforce non-relativistic specific enthalpy ---
      rho0 = 1.0e2_DP * eIntBlast

    END IF

    IF( Verbose )THEN

      WRITE(FMT,'(A)') '(A20,ES11.3E3)'

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'Initializing: ', TRIM( ProgramName )
      WRITE(*,'(4x,A)')   '-------------'
      WRITE(*,*)
      WRITE(*,'(A20,I4.4)') 'nDetCells: ', nDetCells
      WRITE(*,TRIM(FMT)) 'X_D: ', X_D
      WRITE(*,TRIM(FMT)) 'Eblast: ', Eblast
      WRITE(*,TRIM(FMT)) 'eblast: ', eIntBlast
      WRITE(*,TRIM(FMT)) 'rho0: ', rho0
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

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)
        iNX3 = NodeNumberTableX(3,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

        IF( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

          Radius = SQRT( X1**2 + X2**2 + X3**2 )

        ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

          Radius = SQRT( X1**2 + X2**2 )

        ELSE

          Radius = SQRT( X1**2 )

        END IF

        uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX) = rho0
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V1-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V2-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V3-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX) = Zero

        IF( Radius .LT. X_D )THEN

          uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) &
            = eIntBlast

        ELSE

          ! --- Ensure pre-shock pressure is negligible ---
          uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) &
            = 1.0e-10_DP * eIntBlast

        END IF

        IntE_Min_Euler_PL &
          = MIN( IntE_Min_Euler_PL, &
                 1.0e-12_DP * uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) )

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

  END SUBROUTINE InitializeFields_MF


END MODULE MF_InitializationModule
