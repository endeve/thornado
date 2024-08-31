MODULE TimeSteppingModule_SSPRK

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1, &
    nDOFX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE MHD_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    ApplySlopeLimiter_MHD_Relativistic_IDEAL
  USE MHD_DiscretizationModule_Relativistic, ONLY: &
    OffGridFlux_MHD_X1_Inner, &
    OffGridFlux_MHD_X1_Outer, &
    OffGridFlux_MHD_X2_Inner, &
    OffGridFlux_MHD_X2_Outer, &
    OffGridFlux_MHD_X3_Inner, &
    OffGridFlux_MHD_X3_Outer
  USE MHD_TallyModule_Relativistic, ONLY: &
    IncrementOffGridTally_MHD_Relativistic

  IMPLICIT NONE
  PRIVATE

  LOGICAL  :: EvolveOnlyMagnetic
  LOGICAL  :: UseDivergenceCleaning
  LOGICAL  :: UsePowellSource
  INTEGER  :: nStages_SSPRK
  REAL(DP) :: DampingParameter

  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: U_SSPRK
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: D_SSPRK

  PUBLIC :: InitializeMagnetofluid_SSPRK
  PUBLIC :: UpdateMagnetofluid_SSPRK
  PUBLIC :: FinalizeMagnetofluid_SSPRK

  INTERFACE
    SUBROUTINE MagnetofluidIncrement &
      ( t, CFL, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
        SuppressBC_Option, &
        EvolveOnlyMagnetic_Option, &
        UseDivergenceCleaning_Option, &
        DampingParameter_Option, &
        UsePowellSource_Option, &
        SurfaceFlux_X1_Option, &
        SurfaceFlux_X2_Option, &
        SurfaceFlux_X3_Option )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in)     :: t, CFL
      INTEGER,  INTENT(in)     :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(in)    :: &
        G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
        D (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(out)   :: &
        dU(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      LOGICAL,  INTENT(in), OPTIONAL :: &
        SuppressBC_Option
      LOGICAL,  INTENT(in), OPTIONAL :: &
        EvolveOnlyMagnetic_Option
      LOGICAL,  INTENT(in), OPTIONAL :: &
        UseDivergenceCleaning_Option, &
        UsePowellSource_Option
      REAL(DP), INTENT(in), OPTIONAL :: &
        DampingParameter_Option
      REAL(DP), INTENT(out), OPTIONAL :: &
        SurfaceFlux_X1_Option(:,:,:,:,:), &
        SurfaceFlux_X2_Option(:,:,:,:,:), &
        SurfaceFlux_X3_Option(:,:,:,:,:) 
    END SUBROUTINE MagnetofluidIncrement
  END INTERFACE

CONTAINS


  SUBROUTINE InitializeMagnetofluid_SSPRK &
               ( nStages, EvolveOnlyMagnetic_Option, &
                 UseDivergenceCleaning_Option, DampingParameter_Option, &
                 UsePowellSource_Option )

    INTEGER, INTENT(in) :: nStages

    LOGICAL,  INTENT(in), OPTIONAL :: EvolveOnlyMagnetic_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UseDivergenceCleaning_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePowellSource_Option
    REAL(DP), INTENT(in), OPTIONAL :: DampingParameter_Option

    INTEGER :: i

    IF( PRESENT( EvolveOnlyMagnetic_Option ) )THEN
      EvolveOnlyMagnetic = EvolveOnlyMagnetic_Option
    ELSE
      EvolveOnlyMagnetic = .FALSE.
    END IF

    DampingParameter = 0.0_DP
    IF( PRESENT( DampingParameter_Option ) )THEN
      DampingParameter = DampingParameter_Option
    END IF

    UseDivergenceCleaning = .FALSE.
    IF( PRESENT( UseDivergenceCleaning_Option ) ) THEN
      UseDivergenceCleaning = UseDivergenceCleaning_Option
    END IF

    UsePowellSource = .FALSE.
    IF( PRESENT( UsePowellSource_Option ) ) THEN
      UsePowellSource = UsePowellSource_Option
    END IF

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A)') '    INFO: TimeSteppingModule_SSPRK'
    WRITE(*,'(A)') '    ------------------------------'

    WRITE(*,*)
    WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Butcher Table:'
    WRITE(*,'(A5,A)') '', '--------------'
    DO i = 1, nStages
      WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(i), a_SSPRK(i,1:nStages)
    END DO
    WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)

    ALLOCATE( U_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCM) )

    ALLOCATE( D_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCM,1:nStages) )

  END SUBROUTINE InitializeMagnetofluid_SSPRK


  SUBROUTINE FinalizeMagnetofluid_SSPRK

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( U_SSPRK, D_SSPRK )

  END SUBROUTINE FinalizeMagnetofluid_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )
      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_DP
        w_SSPRK(1)   = 1.0_DP

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSPRK(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSPRK(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSPRK(1:3)   = [ 1.0_DP / 6.0_DP, &
                           1.0_DP / 6.0_DP, &
                           2.0_DP / 3.0_DP ]

    END SELECT

    DO iS = 1, nStages
      c_SSPRK(iS) = SUM( a_SSPRK(iS,1:iS-1) )
    END DO

  END SUBROUTINE InitializeSSPRK


  SUBROUTINE AllocateButcherTables_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE UpdateMagnetofluid_SSPRK &
               ( t, dt, CFL, G, U, D, ComputeIncrement_Magnetofluid )

    REAL(DP), INTENT(in) :: &
      t, dt, CFL
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    PROCEDURE (MagnetofluidIncrement) :: &
      ComputeIncrement_Magnetofluid

    INTEGER :: iNX, iX1, iX2, iX3, iCM
    INTEGER :: iS, jS

    REAL(DP) :: dM_OffGrid_MHD(nCM)

    dM_OffGrid_MHD = Zero

    U_SSPRK = Zero ! --- State
    D_SSPRK = Zero ! --- Increment

    DO iS = 1, nStages_SSPRK

      DO iCM = 1, nCM
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        U_SSPRK(iNX,iX1,iX2,iX3,iCM) = U(iNX,iX1,iX2,iX3,iCM)

      END DO
      END DO
      END DO
      END DO
      END DO

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          CALL AddIncrement_Magnetofluid &
                 ( jS, One, U_SSPRK, dt * a_SSPRK(iS,jS), D_SSPRK )

        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( iS .NE. 1)THEN

          CALL ApplySlopeLimiter_MHD_Relativistic_IDEAL &
                 ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U_SSPRK, D )

        END IF

        CALL ComputeIncrement_Magnetofluid &
               ( t, CFL, iX_B0, iX_E0, iX_B1, iX_E1, &
                 G, U_SSPRK, D, D_SSPRK(:,:,:,:,:,iS), &
                 EvolveOnlyMagnetic_Option = EvolveOnlyMagnetic, &
                 UseDivergenceCleaning_Option = UseDivergenceCleaning, &
                 DampingParameter_Option = DampingParameter, &
                 UsePowellSource_Option = UsePowellSource )

        dM_OffGrid_MHD &
          = dM_OffGrid_MHD &
              + dt * w_SSPRK(iS) &
                  * (   OffGridFlux_MHD_X1_Outer &
                      - OffGridFlux_MHD_X1_Inner &
                      + OffGridFlux_MHD_X2_Outer &
                      - OffGridFlux_MHD_X2_Inner &
                      + OffGridFlux_MHD_X3_Outer &
                      - OffGridFlux_MHD_X3_Inner )

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. Zero )THEN

        CALL AddIncrement_Magnetofluid &
               ( iS, One, U, dt * w_SSPRK(iS), D_SSPRK )

      END IF

    END DO

    CALL ApplySlopeLimiter_MHD_Relativistic_IDEAL &
           ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL IncrementOffGridTally_MHD_Relativistic( dM_OffGrid_MHD )

  END SUBROUTINE UpdateMagnetofluid_SSPRK


  SUBROUTINE AddIncrement_Magnetofluid( iS, alpha, U, beta, D )

    INTEGER,  INTENT(in)    :: &
      iS
    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iCM

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      IF( EvolveOnlyMagnetic )THEN

        IF( iCM .LE. 6 )THEN

          U(iNX,iX1,iX2,iX3,iCM) &
            = U(iNX,iX1,iX2,iX3,iCM)

        ELSE

          U(iNX,iX1,iX2,iX3,iCM) &
            = alpha * U(iNX,iX1,iX2,iX3,iCM) &
                + beta * D(iNX,iX1,iX2,iX3,iCM,iS)

        END IF

      ELSE

          U(iNX,iX1,iX2,iX3,iCM) &
            = alpha * U(iNX,iX1,iX2,iX3,iCM) &
                + beta * D(iNX,iX1,iX2,iX3,iCM,iS)

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE AddIncrement_Magnetofluid


END MODULE TimeSteppingModule_SSPRK
