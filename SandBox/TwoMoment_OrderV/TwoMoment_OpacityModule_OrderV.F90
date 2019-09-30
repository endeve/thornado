MODULE TwoMoment_OpacityModule_OrderV

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFZ
  USE RadiationFieldsModule, ONLY: &
    nSpecies

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC :: iOP_D0    = 1
  INTEGER, PUBLIC :: iOP_Chi   = 2
  INTEGER, PUBLIC :: iOP_Sigma = 3
  INTEGER, PUBLIC :: nOP       = 3

  REAL(DP), PUBLIC, ALLOCATABLE :: uOP(:,:,:,:,:,:,:)

  PUBLIC :: SetConstantOpacities
  PUBLIC :: CreateOpacities
  PUBLIC :: DestroyOpacities

CONTAINS


  SUBROUTINE SetConstantOpacities( D0, Chi, Sigma )

    REAL(DP), INTENT(in) :: D0, Chi, Sigma

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Setting Constant Opacities:'
    WRITE(*,*)
    WRITE(*,'(A7,A8,ES10.4E2)') '',    'D0 = ', D0
    WRITE(*,'(A7,A8,ES10.4E2)') '',   'Chi = ', Chi
    WRITE(*,'(A7,A8,ES10.4E2)') '', 'Sigma = ', Sigma

    uOP(:,:,:,:,:,iOP_D0   ,:) = D0
    uOP(:,:,:,:,:,iOP_Chi  ,:) = Chi
    uOP(:,:,:,:,:,iOP_Sigma,:) = Sigma

  END SUBROUTINE SetConstantOpacities


  SUBROUTINE CreateOpacities( nX, swX, nE, swE, Verbose_Option )

    INTEGER, INTENT(in) :: nX(3), swX(3), nE, swE
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A29,I2.2)') &
        '', 'Simple Opacities, nSpecies = ', nSpecies
    END IF

    ALLOCATE &
      ( uOP(1:nDOFZ, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nOP,1:nSpecies) )

  END SUBROUTINE CreateOpacities


  SUBROUTINE DestroyOpacities

    DEALLOCATE( uOP )

  END SUBROUTINE DestroyOpacities


END MODULE TwoMoment_OpacityModule_OrderV
