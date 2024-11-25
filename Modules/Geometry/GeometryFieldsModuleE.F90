MODULE GeometryFieldsModuleE

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFE

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- 1D Momentum Space (Energy) Geometry Fields ---

  INTEGER, PUBLIC, PARAMETER :: iGE_Ep0 = 1 ! E**0
  INTEGER, PUBLIC, PARAMETER :: iGE_Ep1 = 2 ! E**1
  INTEGER, PUBLIC, PARAMETER :: iGE_Ep2 = 3 ! E**2
  INTEGER, PUBLIC, PARAMETER :: iGE_Ep3 = 4 ! E**3
  INTEGER, PUBLIC, PARAMETER :: nGE     = 4

  CHARACTER(32), DIMENSION(nGE), PUBLIC, PARAMETER :: &
    namesGE = [ 'Energy Power 0                              ', &
                'Energy Power 1                              ', &
                'Energy Power 2                              ', &
                'Energy Power 3                              ' ]

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: uGE

  PUBLIC :: CreateGeometryFieldsE
  PUBLIC :: DestroyGeometryFieldsE

CONTAINS


  SUBROUTINE CreateGeometryFieldsE( nE, swE, Verbose_Option )

    INTEGER, INTENT(in)           :: nE, swE
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iGE

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Geometry Fields (Energy)'
      WRITE(*,*)
      DO iGE = 1, nGE
        WRITE(*,'(A5,A32)') '', TRIM( namesGE(iGE) )
      END DO
    END IF

    ALLOCATE( uGE(1:nDOFE,1-swE:nE+swE,1:nGE) )

    uGE(:,:,iGE_Ep0) = 0.0_DP
    uGE(:,:,iGE_Ep1) = 0.0_DP
    uGE(:,:,iGE_Ep2) = 0.0_DP
    uGE(:,:,iGE_Ep3) = 0.0_DP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uGE )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uGE )
#endif

  END SUBROUTINE CreateGeometryFieldsE


  SUBROUTINE DestroyGeometryFieldsE

    IF (ALLOCATED( uGE )) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: uGE )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( uGE )
#endif

      DEALLOCATE( uGE )

    END IF

  END SUBROUTINE DestroyGeometryFieldsE


END MODULE GeometryFieldsModuleE
