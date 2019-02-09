MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---
  USE amrex_base_module, ONLY: &
    amrex_box,               &
    amrex_geometry,          &
    amrex_boxarray,          &
    amrex_boxarray_build,    &
    amrex_boxarray_destroy,  &
    amrex_distromap,         &
    amrex_distromap_build,   &
    amrex_distromap_destroy, &
    amrex_multifab,          &
    amrex_multifab_build,    &
    amrex_multifab_destroy
  USE amrex_fort_module, ONLY: &
    amrex_real

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nX
  USE FluidFieldsModule,        ONLY: &
    nCF
  USE GeometryFieldsModule,     ONLY: &
    nGF

  ! --- Local Modules ---
  USE MF_SlopeLimiterModule_Euler,      ONLY: &
    MF_ApplySlopeLimiter_Euler
  USE MF_PositivityLimiterModule_Euler, ONLY: &
    MF_ApplyPositivityLimiter_Euler
  USE MF_UtilitiesModule,               ONLY: &
    LinComb


  IMPLICIT NONE
  PRIVATE

  INTEGER :: nStages_SSPRK
  REAL(amrex_real), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(amrex_real), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(amrex_real), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  TYPE(amrex_multifab), DIMENSION(:), ALLOCATABLE :: MF_U
  TYPE(amrex_multifab), DIMENSION(:), ALLOCATABLE :: MF_D

  PUBLIC :: MF_InitializeFluid_SSPRK
  PUBLIC :: MF_UpdateFluid_SSPRK
  PUBLIC :: MF_FinalizeFluid_SSPRK

  INTERFACE
    SUBROUTINE MF_FluidIncrement &
      ( nLevels, GEOM, MF_uGF, MF_uCF, MF_duCF, iS )
      USE amrex_base_module, ONLY: &
        amrex_multifab, &
        amrex_geometry
      INTEGER,              INTENT(in)    :: nLevels, iS
      TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels)
      TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels)
      TYPE(amrex_multifab), INTENT(inout) :: MF_uCF (0:nLevels)
      TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:nLevels)
    END SUBROUTINE MF_FluidIncrement
  END INTERFACE

CONTAINS


  SUBROUTINE MF_InitializeFluid_SSPRK( nLevels, nStages, BA, DM )

    INTEGER,               INTENT(in) :: nLevels
    INTEGER,               INTENT(in) :: nStages
    TYPE(amrex_boxarray),  INTENT(in) :: BA(0:nLevels)
    TYPE(amrex_distromap), INTENT(in) :: DM(0:nLevels)

    INTEGER         :: i, iLevel
    TYPE(amrex_box) :: BX

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    WRITE(*,*)
    WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Butcher Table:'
    WRITE(*,'(A5,A)') '', '--------------'
    DO i = 1, nStages
      WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(i), a_SSPRK(i,1:nStages)
    END DO
    WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)
    WRITE(*,*)

    ALLOCATE( MF_U(0:nLevels) )
    ALLOCATE( MF_D(0:nLevels) )

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    DO iLevel = 0, nLevels
      CALL amrex_multifab_build &
        ( MF_U(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )
      CALL amrex_multifab_build &
        ( MF_D(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF * nStages, swX(1) )
    END DO

  END SUBROUTINE MF_InitializeFluid_SSPRK


  SUBROUTINE MF_FinalizeFluid_SSPRK( nLevels )

    INTEGER, INTENT(in) :: nLevels
    INTEGER :: iLevel

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DO iLevel = 0, nLevels
      CALL amrex_multifab_destroy( MF_U(iLevel) )
      CALL amrex_multifab_destroy( MF_D(iLevel) )
    END DO
    DEALLOCATE( MF_U )
    DEALLOCATE( MF_D )

  END SUBROUTINE MF_FinalizeFluid_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )
      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_amrex_real
        w_SSPRK(1)   = 1.0_amrex_real

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_amrex_real, 0.0_amrex_real ]
        a_SSPRK(2,1:2) = [ 1.0_amrex_real, 0.0_amrex_real ]
        w_SSPRK(1:2)   = [ 0.5_amrex_real, 0.5_amrex_real ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_amrex_real, 0.00_amrex_real, 0.00_amrex_real ]
        a_SSPRK(2,1:3) = [ 1.00_amrex_real, 0.00_amrex_real, 0.00_amrex_real ]
        a_SSPRK(3,1:3) = [ 0.25_amrex_real, 0.25_amrex_real, 0.00_amrex_real ]
        w_SSPRK(1:3)   = [ 1.0_amrex_real / 6.0_amrex_real, &
                           1.0_amrex_real / 6.0_amrex_real, &
                           2.0_amrex_real / 3.0_amrex_real ]

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

    a_SSPRK = 0.0_amrex_real
    c_SSPRK = 0.0_amrex_real
    w_SSPRK = 0.0_amrex_real

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE MF_UpdateFluid_SSPRK &
              ( nLevels, t, dt, MF_uGF, MF_uCF, &
                GEOM, MF_ComputeIncrement_Fluid )

    INTEGER,              INTENT(in)    :: nLevels
    REAL(amrex_real),     INTENT(in)    :: t, dt (0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels)
    PROCEDURE(MF_FluidIncrement)        :: MF_ComputeIncrement_Fluid

    INTEGER :: iLevel
    INTEGER :: iS, jS

    DO iLevel = 0, nLevels
      CALL MF_U(iLevel) % setval( 0.0_amrex_real )
      CALL MF_D(iLevel) % setval( 0.0_amrex_real )
    END DO

    DO iS = 1, nStages_SSPRK

      DO iLevel = 0, nLevels
        CALL MF_U(iLevel) &
               % COPY( MF_uCF(iLevel), 1, 1, MF_uCF(iLevel) % nComp(), swX(1) )
      END DO

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. 0.0_amrex_real ) &
          CALL LinComb( nLevels, 1.0_amrex_real, MF_U, &
                                 dt * a_SSPRK(iS,jS), MF_D, jS )

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. 0.0_amrex_real ) &
          .OR. ( w_SSPRK(iS) .NE. 0.0_amrex_real ) )THEN

        CALL MF_ApplySlopeLimiter_Euler     ( nLevels, MF_uGF, MF_U, GEOM )
        CALL MF_ApplyPositivityLimiter_Euler( nLevels, MF_uGF, MF_U )

        CALL MF_ComputeIncrement_Fluid &
            ( nLevels, GEOM, MF_uGF, MF_U, MF_D, iS )

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. 0.0_amrex_real ) &
        CALL LinComb( nLevels, 1.0_amrex_real, MF_uCF, &
                               dt * w_SSPRK(iS), MF_D, iS )

    END DO

    CALL MF_ApplySlopeLimiter_Euler     ( nLevels, MF_uGF, MF_uCF, GEOM )
    CALL MF_ApplyPositivityLimiter_Euler( nLevels, MF_uGF, MF_uCF )

  END SUBROUTINE MF_UpdateFluid_SSPRK


END MODULE MF_TimeSteppingModule_SSPRK
