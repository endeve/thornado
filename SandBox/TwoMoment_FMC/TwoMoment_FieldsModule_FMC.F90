MODULE TwoMoment_FieldsModule_FMC

  USE KindModule, ONLY: &
    DP, One
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFX

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  INTEGER, PUBLIC            :: nSpecies
  INTEGER, PUBLIC, PARAMETER :: iNuE     = 1 ! Electron Neutrino
  INTEGER, PUBLIC, PARAMETER :: iNuE_Bar = 2 ! Electron Antineutrino
  INTEGER, PUBLIC, PARAMETER :: iNuM     = 3 ! Muon Neutrino
  INTEGER, PUBLIC, PARAMETER :: iNuM_Bar = 4 ! Muon Antineutrino
  INTEGER, PUBLIC, PARAMETER :: iNuT     = 5 ! Tau Neutrino
  INTEGER, PUBLIC, PARAMETER :: iNuT_Bar = 6 ! Tau Antineutrino

  REAL(DP), DIMENSION(6), PUBLIC :: LeptonNumber

  ! --- Eulerian (Conserved) Two-Moment Fields ---

  INTEGER, PUBLIC, PARAMETER :: iCM_E  = 1  ! Eulerian Energy Density
  INTEGER, PUBLIC, PARAMETER :: iCM_F1 = 2  ! Eulerian Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iCM_F2 = 3  ! Eulerian Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iCM_F3 = 4  ! Eulerian Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: nCM    = 4  ! n Eulerian Two-Moment Fields

  CHARACTER(32), DIMENSION(nCM), PUBLIC, PARAMETER :: &
    namesCM = [ 'Eulerian Energy Density         ', &
                'Eulerian Momentum Density (1)   ', &
                'Eulerian Momentum Density (2)   ', &
                'Eulerian Momentum Density (3)   ' ]

  CHARACTER(5),  DIMENSION(nCM), PUBLIC, PARAMETER :: &
    ShortNamesCM = [ 'CM_E ', &
                     'CM_F1', &
                     'CM_F2', &
                     'CM_F3' ]

  REAL(DP), DIMENSION(nCM), PUBLIC :: unitsCM

  REAL(DP), ALLOCATABLE, PUBLIC :: uCM(:,:,:,:,:,:,:)

  ! --- Lagrangian (Primitive) Two-Moment Fields ---

  INTEGER, PUBLIC, PARAMETER :: iPM_J  = 1  ! Lagrangian Energy Density
  INTEGER, PUBLIC, PARAMETER :: iPM_H1 = 2  ! Lagrangian Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iPM_H2 = 3  ! Lagrangian Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iPM_H3 = 4  ! Lagrangian Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: nPM    = 4  ! n Lagrangian Two-Moment Fields

  CHARACTER(32), DIMENSION(nPM), PUBLIC, PARAMETER :: &
    namesPM = [ 'Lagrangian Energy Density       ', &
                'Lagrangian Momentum Density (1) ', &
                'Lagrangian Momentum Density (2) ', &
                'Lagrangian Momentum Density (3) ' ]

  CHARACTER(5),  DIMENSION(nPM), PUBLIC, PARAMETER :: &
    ShortNamesPR = [ 'PM_J ', &
                     'PM_H1', &
                     'PM_H2', &
                     'PM_H3' ]

  REAL(DP), DIMENSION(nPM), PUBLIC :: unitsPM

  REAL(DP), ALLOCATABLE, PUBLIC :: uPM(:,:,:,:,:,:,:)

  ! --- Auxiliary Two-Moment Fields ---

  INTEGER, PUBLIC, PARAMETER :: iAM_F = 1  ! Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iAM_K = 2  ! Eddington Factor
  INTEGER, PUBLIC, PARAMETER :: iAM_Q = 3  ! Heat Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iAM_N = 4  ! Eulerian Number Density
  INTEGER, PUBLIC, PARAMETER :: nAM   = 4  ! n Auxiliary Two-Moment Fields

  CHARACTER(32), DIMENSION(nAM), PUBLIC, PARAMETER :: &
    namesAM = [ 'Lagrangian Flux Factor          ', &
                'Lagrangian Eddington Factor     ', &
                'Lagrangian Heat Flux Factor     ', &
                'Eulerian Number Density         ' ]

  REAL(DP), DIMENSION(nAM), PUBLIC :: unitsAM

  REAL(DP), ALLOCATABLE, PUBLIC :: uAM(:,:,:,:,:,:,:)

  ! --- Grey (Energy-Integrated) Two-Moment Variables ---

  INTEGER, PUBLIC, PARAMETER :: iGM_E   = 1  ! Eulerian Energy Density
  INTEGER, PUBLIC, PARAMETER :: iGM_F1  = 2  ! Eulerian Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iGM_F2  = 3  ! Eulerian Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iGM_F3  = 4  ! Eulerian Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: iGM_J   = 5  ! Lagrangian Energy Density
  INTEGER, PUBLIC, PARAMETER :: iGM_H1  = 6  ! Lagrangian Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iGM_H2  = 7  ! Lagrangian Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iGM_H3  = 8  ! Lagrangian Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: iGM_RMS = 9  ! RMS Energy
  INTEGER, PUBLIC, PARAMETER :: iGM_F   = 10 ! Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iGM_K   = 11 ! Eddington Factor
  INTEGER, PUBLIC, PARAMETER :: iGM_Q   = 12 ! Heat Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iGM_N   = 13 ! Eulerian   Number Density
  INTEGER, PUBLIC, PARAMETER :: nGM     = 13 ! n Gray Two-Moment Fields

  CHARACTER(32), DIMENSION(nGM), PUBLIC, PARAMETER :: &
    namesGM = [ 'Eulerian Energy Density         ', &
                'Eulerian Momentum Density (1)   ', &
                'Eulerian Momentum Density (2)   ', &
                'Eulerian Momentum Density (3)   ', &
                'Lagrangian Energy Density       ', &
                'Lagrangian Momentum Density (1) ', &
                'Lagrangian Momentum Density (2) ', &
                'Lagrangian Momentum Density (3) ', &
                'RMS Energy                      ', &
                'Lagrangian Flux Factor          ', &
                'Lagrangian Eddington Factor     ', &
                'Lagrangian Heat Flux Factor     ', &
                'Eulerian Number Density         ' ]

  CHARACTER(6), DIMENSION(nGM), PUBLIC, PARAMETER :: &
    ShortNamesGR = [ 'GM_E  ', &
                     'GM_F1 ', &
                     'GM_F2 ', &
                     'GM_F3 ', &
                     'GM_J  ', &
                     'GM_H1 ', &
                     'GM_H2 ', &
                     'GM_H3 ', &
                     'GM_RMS', &
                     'GM_F  ', &
                     'GM_K  ', &
                     'GM_Q  ', &
                     'GM_N  ' ]

  REAL(DP), DIMENSION(nGM), PUBLIC :: unitsGM

  REAL(DP), ALLOCATABLE, PUBLIC :: uGM(:,:,:,:,:,:)

  PUBLIC :: CreateTwoMomentFields
  PUBLIC :: DestroyTwoMomentFields
  PUBLIC :: SetUnitsTwoMomentFields
  PUBLIC :: DescribeTwoMomentFields_Conserved
  PUBLIC :: DescribeTwoMomentFields_Primitive
  PUBLIC :: SetNumberOfSpecies

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET( LeptonNumber )
#elif defined(THORNADO_OACC)
    !$ACC DECLARE CREATE( LeptonNumber )
#endif

CONTAINS


  SUBROUTINE SetNumberOfSpecies( nS )

    INTEGER, INTENT(in) :: nS

    nSpecies = nS

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A30,I2.2)') &
        '', 'Two-Moment Fields, nSpecies = ', nSpecies
    END IF

    LeptonNumber = [ One, - One, One, - One, One, - One ]

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( LeptonNumber )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( LeptonNumber )
#endif

  END SUBROUTINE SetNumberOfSpecies


  SUBROUTINE CreateTwoMomentFields &
    ( nX, swX, nE, swE, nSpecies_Option, Verbose_Option )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE,    swE
    INTEGER, INTENT(in), OPTIONAL :: nSpecies_Option
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: nS

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    nS = 1
    IF( PRESENT( nSpecies_Option ) ) &
      nS = nSpecies_Option

    CALL SetNumberOfSpecies( nS )

    CALL CreateTwoMomentFields_Conserved ( nX, swX, nE, swE )
    CALL CreateTwoMomentFields_Primitive ( nX, swX, nE, swE )
    CALL CreateTwoMomentFields_Auxiliary ( nX, swX, nE, swE )
    CALL CreateTwoMomentFields_Gray      ( nX, swX )

    CALL SetUnitsTwoMomentFields

  END SUBROUTINE CreateTwoMomentFields


  SUBROUTINE CreateTwoMomentFields_Conserved( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE,    swE

    CALL DescribeTwoMomentFields_Conserved( Verbose )

    ALLOCATE &
      ( uCM(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nCM, 1:nSpecies) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: uCM )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( uCM )
#endif

  END SUBROUTINE CreateTwoMomentFields_Conserved


  SUBROUTINE DescribeTwoMomentFields_Conserved( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iCM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A29)') '', 'Two-Moment Fields (Conserved)'
      WRITE(*,*)
      DO iCM = 1, nCM
        WRITE(*,'(A5,A32)') '', TRIM( namesCM(iCM) )
      END DO

    END IF

  END SUBROUTINE DescribeTwoMomentFields_Conserved


  SUBROUTINE CreateTwoMomentFields_Primitive( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE,    swE

    CALL DescribeTwoMomentFields_Primitive( Verbose )

    ALLOCATE &
      ( uPM(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nPM, 1:nSpecies) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: uPM )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( uPM )
#endif

  END SUBROUTINE CreateTwoMomentFields_Primitive


  SUBROUTINE DescribeTwoMomentFields_Primitive( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iPM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A29)') '', 'Two-Moment Fields (Primitive)'
      WRITE(*,*)
      DO iPM = 1, nPM
        WRITE(*,'(A5,A32)') '', TRIM( namesPM(iPM) )
      END DO

    END IF

  END SUBROUTINE DescribeTwoMomentFields_Primitive


  SUBROUTINE CreateTwoMomentFields_Auxiliary( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE, swE

    INTEGER :: iAM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A29)') '', 'Two-Moment Fields (Auxiliary)'
      WRITE(*,*)
      DO iAM = 1, nAM
        WRITE(*,'(A5,A)') '', TRIM( namesAM(iAM) )
      END DO

    END IF

    ALLOCATE &
      ( uAM(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nAM, 1:nSpecies) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: uAM )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( uAM )
#endif

  END SUBROUTINE CreateTwoMomentFields_Auxiliary


  SUBROUTINE CreateTwoMomentFields_Gray( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iGM

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Two-Moment Fields (Gray)'
      WRITE(*,*)
      DO iGM = 1, nGM
        WRITE(*,'(A5,A)') '', TRIM( namesGM(iGM) )
      END DO
    END IF

    ALLOCATE &
      ( uGM(1:nDOFX, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nGM, 1:nSpecies) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: uGM )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( uGM )
#endif

  END SUBROUTINE CreateTwoMomentFields_Gray


  SUBROUTINE DestroyTwoMomentFields

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uCM, uPM, uAM, uGM, &
    !$OMP               unitsCM, unitsPM, unitsAM, unitsGM )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( uCM, uPM, uPM, uAM, uGM, &
    !$ACC         unitsCM, unitsPM, unitsAM, unitsGM, unitsDM )
#endif

    DEALLOCATE( uCM, uPM, uAM, uGM )

  END SUBROUTINE DestroyTwoMomentFields


  SUBROUTINE SetUnitsTwoMomentFields

    unitsCM = One
    unitsPM = One
    unitsAM = One
    unitsGM = One

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: unitsCM, unitsPM, unitsAM, unitsGM )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( unitsCM, unitsPM, unitsAM, unitsGM )
#endif

  END SUBROUTINE SetUnitsTwoMomentFields


END MODULE TwoMoment_FieldsModule_FMC
