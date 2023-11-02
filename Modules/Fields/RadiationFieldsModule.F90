MODULE RadiationFieldsModule

  USE KindModule, ONLY: &
    DP, One, Zero
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

  ! --- Eulerian (Conserved) Radiation Fields ---

  INTEGER, PUBLIC, PARAMETER :: iCR_N  = 1  ! Eulerian Number Density
  INTEGER, PUBLIC, PARAMETER :: iCR_G1 = 2  ! Eulerian Number Flux Density 1
  INTEGER, PUBLIC, PARAMETER :: iCR_G2 = 3  ! Eulerian Number Flux Density 2
  INTEGER, PUBLIC, PARAMETER :: iCR_G3 = 4  ! Eulerian Number Flux Density 3
  INTEGER, PUBLIC, PARAMETER :: nCR    = 4  ! n Eulerian Radiation Fields

  CHARACTER(32), DIMENSION(nCR), PUBLIC, PARAMETER :: &
    namesCR = [ 'Eulerian Number Density         ', &
                'Eulerian Number Flux Density (1)', &
                'Eulerian Number Flux Density (2)', &
                'Eulerian Number Flux Density (3)' ]

  CHARACTER(5),  DIMENSION(nCR), PUBLIC, PARAMETER :: &
    ShortNamesCR = [ 'CR_N ', &
                     'CR_G1', &
                     'CR_G2', &
                     'CR_G3' ]

  REAL(DP), DIMENSION(nCR), PUBLIC :: unitsCR

  REAL(DP), ALLOCATABLE, PUBLIC :: uCR  (:,:,:,:,:,:,:)

  ! --- Lagrangian (Primitive) Radiation Fields ---

  INTEGER, PUBLIC, PARAMETER :: iPR_D  = 1  ! Lagrangian Number Density
  INTEGER, PUBLIC, PARAMETER :: iPR_I1 = 2  ! Lagrangian Number Flux 1
  INTEGER, PUBLIC, PARAMETER :: iPR_I2 = 3  ! Lagrangian Number Flux 2
  INTEGER, PUBLIC, PARAMETER :: iPR_I3 = 4  ! Lagrangian Number Flux 3
  INTEGER, PUBLIC, PARAMETER :: nPR    = 4  ! n Lagrangian Radiation Fields

  CHARACTER(34), DIMENSION(nPR), PUBLIC, PARAMETER :: &
    namesPR = [ 'Lagrangian Number Density         ', &
                'Lagrangian Number Flux Density (1)', &
                'Lagrangian Number Flux Density (2)', &
                'Lagrangian Number Flux Density (3)' ]

  CHARACTER(5),  DIMENSION(nPR), PUBLIC, PARAMETER :: &
    ShortNamesPR = [ 'PR_D ', &
                     'PR_I1', &
                     'PR_I2', &
                     'PR_I3' ]

  REAL(DP), DIMENSION(nPR), PUBLIC :: unitsPR

  REAL(DP), ALLOCATABLE, PUBLIC :: uPR(:,:,:,:,:,:,:)

  ! --- Auxiliary Radiation Fields ---

  INTEGER, PUBLIC, PARAMETER :: iAR_F = 1  ! Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iAR_K = 2  ! Eddington Factor
  INTEGER, PUBLIC, PARAMETER :: iAR_Q = 3  ! Heat Flux Factor
  INTEGER, PUBLIC, PARAMETER :: nAR   = 3  ! n Auxiliary Radiation Fields

  CHARACTER(32), DIMENSION(nAR), PUBLIC, PARAMETER :: &
    namesAR = [ 'Lagrangian Flux Factor          ', &
                'Lagrangian Eddington Factor     ', &
                'Lagrangian Heat Flux Factor     ' ]

  REAL(DP), DIMENSION(nAR), PUBLIC :: unitsAR

  REAL(DP), ALLOCATABLE, PUBLIC :: uAR(:,:,:,:,:,:,:)

  ! --- Grey (Energy-Integrated) Radiation Variables ---

  INTEGER, PUBLIC, PARAMETER :: iGR_N   = 1  ! Eulerian   Number Density
  INTEGER, PUBLIC, PARAMETER :: iGR_D   = 2  ! Lagrangian Number Density
  INTEGER, PUBLIC, PARAMETER :: iGR_I1  = 3  ! Lagrangian Number Flux 1
  INTEGER, PUBLIC, PARAMETER :: iGR_I2  = 4  ! Lagrangian Number Flux 2
  INTEGER, PUBLIC, PARAMETER :: iGR_I3  = 5  ! Lagrangian Number Flux 3
  INTEGER, PUBLIC, PARAMETER :: iGR_J   = 6  ! Lagrangian Energy Density
  INTEGER, PUBLIC, PARAMETER :: iGR_H1  = 7  ! Lagrangian Energy Flux 1
  INTEGER, PUBLIC, PARAMETER :: iGR_H2  = 8  ! Lagrangian Energy Flux 2
  INTEGER, PUBLIC, PARAMETER :: iGR_H3  = 9  ! Lagrangian Energy Flux 3
  INTEGER, PUBLIC, PARAMETER :: iGR_RMS = 10 ! RMS Energy
  INTEGER, PUBLIC, PARAMETER :: iGR_F   = 11 ! Flux Factor
  INTEGER, PUBLIC, PARAMETER :: iGR_K   = 12 ! Eddington Factor
  INTEGER, PUBLIC, PARAMETER :: iGR_Q   = 13 ! Heat Flux Factor
  INTEGER, PUBLIC, PARAMETER :: nGR     = 13 ! n Gray Radiation Fields

  CHARACTER(34), DIMENSION(nGR), PUBLIC, PARAMETER :: &
    namesGR = [ 'Eulerian Number Density           ', &
                'Lagrangian Number Density         ', &
                'Lagrangian Number Flux Density (1)', &
                'Lagrangian Number Flux Density (2)', &
                'Lagrangian Number Flux Density (3)', &
                'Lagrangian Energy Density         ', &
                'Lagrangian Energy Flux Density (1)', &
                'Lagrangian Energy Flux Density (2)', &
                'Lagrangian Energy Flux Density (3)', &
                'RMS Energy                        ', &
                'Lagrangian Flux Factor            ', &
                'Lagrangian Eddington Factor       ', &
                'Lagrangian Heat Flux Factor       ' ]

  CHARACTER(6), DIMENSION(nGR), PUBLIC, PARAMETER :: &
    ShortNamesGR = [ 'GR_N  ', &
                     'GR_D  ', &
                     'GR_I1 ', &
                     'GR_I2 ', &
                     'GR_I3 ', &
                     'GR_J  ', &
                     'GR_H1 ', &
                     'GR_H2 ', &
                     'GR_H3 ', &
                     'GR_RMS', &
                     'GR_F  ', &
                     'GR_K  ', &
                     'GR_Q  ' ]
  REAL(DP), DIMENSION(nGR), PUBLIC :: unitsGR

  REAL(DP), ALLOCATABLE, PUBLIC :: uGR(:,:,:,:,:,:)

  ! --- Diagnostic Radiation Variables ---

  INTEGER, PUBLIC, PARAMETER :: iDR_iter_outer  = 1
  INTEGER, PUBLIC, PARAMETER :: iDR_iter_inner  = 2
  INTEGER, PUBLIC, PARAMETER :: iDR_PL_Theta_1  = 3
  INTEGER, PUBLIC, PARAMETER :: iDR_PL_Theta_2  = 4
  INTEGER, PUBLIC, PARAMETER :: iDR_PL_dEnergy  = 5
  INTEGER, PUBLIC, PARAMETER :: nDR             = 5

  CHARACTER(32), DIMENSION(nDR), PUBLIC, PARAMETER :: &
    namesDR = [ 'Outer Iterations                ', &
                'Inner Iterations                ', &
                'Positivity Limiter Theta 1      ', &
                'Positivity Limiter Theta 2      ', &
                'Positivity Limiter Energy Change' ]

  REAL(DP), DIMENSION(nDR), PUBLIC :: unitsDR

  REAL(DP), ALLOCATABLE, PUBLIC :: uDR(:,:,:,:)

  PUBLIC :: CreateRadiationFields
  PUBLIC :: DestroyRadiationFields
  PUBLIC :: SetUnitsRadiationFields
  PUBLIC :: DescribeRadiationFields_Conserved
  PUBLIC :: DescribeRadiationFields_Primitive
  PUBLIC :: SetNumberOfSpecies

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET( LeptonNumber )
#elif defined(THORNADO_OACC)
    !$ACC DECLARE CREATE( LeptonNumber )
#endif

CONTAINS


  SUBROUTINE SetNumberOfSpecies( nS, Verbose_Option )

    INTEGER, INTENT(in) :: nS
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    nSpecies = nS

    IF( PRESENT( Verbose_Option ) )THEN
      IF( Verbose_Option )THEN
        WRITE(*,*)
        WRITE(*,'(A5,A29,I2.2)') &
          '', 'Radiation Fields, nSpecies = ', nSpecies
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,'(A5,A29,I2.2)') &
        '', 'Radiation Fields, nSpecies = ', nSpecies
    END IF

    LeptonNumber = [ One, - One, One, - One, One, - One ]

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( LeptonNumber )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( LeptonNumber )
#endif

  END SUBROUTINE SetNumberOfSpecies


  SUBROUTINE CreateRadiationFields &
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

    CALL SetNumberOfSpecies( nS, Verbose_Option = Verbose )

    CALL CreateRadiationFields_Conserved ( nX, swX, nE, swE )
    CALL CreateRadiationFields_Primitive ( nX, swX, nE, swE )
    CALL CreateRadiationFields_Auxiliary ( nX, swX, nE, swE )
    CALL CreateRadiationFields_Gray      ( nX, swX )
    CALL CreateRadiationFields_Diagnostic( nX, swX )

    CALL SetUnitsRadiationFields

  END SUBROUTINE CreateRadiationFields


  SUBROUTINE CreateRadiationFields_Conserved( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE,    swE

    CALL DescribeRadiationFields_Conserved( Verbose )

    ALLOCATE &
      ( uCR(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nCR, 1:nSpecies) )

    uCR = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uCR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uCR )
#endif

  END SUBROUTINE CreateRadiationFields_Conserved


  SUBROUTINE DescribeRadiationFields_Conserved( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iCR

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A28)') '', 'Radiation Fields (Conserved)'
      WRITE(*,*)
      DO iCR = 1, nCR
        WRITE(*,'(A5,A32)') '', TRIM( namesCR(iCR) )
      END DO

    END IF

  END SUBROUTINE DescribeRadiationFields_Conserved


  SUBROUTINE CreateRadiationFields_Primitive( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE,    swE

    CALL DescribeRadiationFields_Primitive( Verbose )

    ALLOCATE &
      ( uPR(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nPR, 1:nSpecies) )

    uPR = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uPR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uPR )
#endif

  END SUBROUTINE CreateRadiationFields_Primitive


  SUBROUTINE DescribeRadiationFields_Primitive( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iPR

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A28)') '', 'Radiation Fields (Primitive)'
      WRITE(*,*)
      DO iPR = 1, nPR
        WRITE(*,'(A5,A34)') '', TRIM( namesPR(iPR) )
      END DO

    END IF

  END SUBROUTINE DescribeRadiationFields_Primitive


  SUBROUTINE CreateRadiationFields_Auxiliary( nX, swX, nE, swE )

    INTEGER, INTENT(in) :: nX(3), swX(3)
    INTEGER, INTENT(in) :: nE, swE

    INTEGER :: iAR

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A28)') '', 'Radiation Fields (Auxiliary)'
      WRITE(*,*)
      DO iAR = 1, nAR
        WRITE(*,'(A5,A)') '', TRIM( namesAR(iAR) )
      END DO

    END IF

    ALLOCATE &
      ( uAR(1:nDOF, &
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nAR, 1:nSpecies) )

    uAR = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uAR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uAR )
#endif

  END SUBROUTINE CreateRadiationFields_Auxiliary


  SUBROUTINE CreateRadiationFields_Gray( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iGR

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A23)') '', 'Radiation Fields (Gray)'
      WRITE(*,*)
      DO iGR = 1, nGR
        WRITE(*,'(A5,A)') '', TRIM( namesGR(iGR) )
      END DO
    END IF

    ALLOCATE &
      ( uGR(1:nDOFX, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nGR, 1:nSpecies) )

    uGR = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uGR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uGR )
#endif

  END SUBROUTINE CreateRadiationFields_Gray


  SUBROUTINE CreateRadiationFields_Diagnostic( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iDR

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A29)') '', 'Radiation Fields (Diagnostic)'
      WRITE(*,*)
      DO iDR = 1, nDR
        WRITE(*,'(A5,A)') '', TRIM( namesDR(iDR) )
      END DO
    END IF

    ALLOCATE( uDR(1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nDR) )

    uDR = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uDR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uDR )
#endif

  END SUBROUTINE CreateRadiationFields_Diagnostic


  SUBROUTINE DestroyRadiationFields

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uCR, uPR, uAR, uGR, uDR, &
    !$OMP               unitsCR, unitsPR, unitsAR, unitsGR, unitsDR )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( uCR, uPR, uPR, uAR, uGR, uDR, &
    !$ACC         unitsCR, unitsPR, unitsAR, unitsGR, unitsDR )
#endif

    DEALLOCATE( uCR, uPR, uAR, uGR, uDR )

  END SUBROUTINE DestroyRadiationFields


  SUBROUTINE SetUnitsRadiationFields

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Centimeter, &
      Second, &
      Erg, MeV

    unitsCR = One
    unitsPR = One
    unitsAR = One
    unitsGR = One
    unitsDR = One

    IF( UnitsActive )THEN

      ! --- Gray Units ---

      unitsGR(iGR_N)   = One / Centimeter**3
      unitsGR(iGR_D)   = One / Centimeter**3
      unitsGR(iGR_I1)  = One / Centimeter**2 / Second
      unitsGR(iGR_I2)  = One / Centimeter**2 / Second
      unitsGR(iGR_I3)  = One / Centimeter**2 / Second
      unitsGR(iGR_J)   = Erg / Centimeter**3
      unitsGR(iGR_H1)  = Erg / Centimeter**2 / Second
      unitsGR(iGR_H2)  = Erg / Centimeter**2 / Second
      unitsGR(iGR_H3)  = Erg / Centimeter**2 / Second
      unitsGR(iGR_RMS) = MeV

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: unitsCR, unitsPR, unitsAR, unitsGR, unitsDR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( unitsCR, unitsPR, unitsAR, unitsGR, unitsDR )
#endif

  END SUBROUTINE SetUnitsRadiationFields


END MODULE RadiationFieldsModule
