MODULE UnitsModule

  USE KindModule, ONLY: &
    DP
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS, &
    GravitationalConstantMKS, &
    BoltzmannConstantMKS, &
    ElectronVoltMKS, &
    PlanckConstantMKS, &
    AvogadroConstantMKS

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC :: UnitsActive = .FALSE.

  REAL(DP), PUBLIC, PARAMETER :: &
    SpeedOfLight          = 1.0_DP, &
    GravitationalConstant = 1.0_DP, &
    BoltzmannConstant     = 1.0_DP

  ! --- Length ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Meter      = 1.0_DP, &
    Centimeter = 1.0e-2_DP * Meter, &
    Kilometer  = 1.0e+3_DP * Meter

  ! --- Time ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Second      = SpeedOfLightMKS / SpeedOfLight * Meter, &
    Millisecond = 1.0e-3_DP * Second, &
    Microsecond = 1.0e-6_DP * Second

  ! --- Mass ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Kilogram  = GravitationalConstantMKS / GravitationalConstant &
                  * Meter**3 / Second**2, &
    Gram      = 1.0e-3_DP * Kilogram, &
    SolarMass = 1.98892e30_DP * Kilogram

  ! --- Other Units of Measure and Constants ---

  REAL(DP), PUBLIC, PARAMETER :: &
    Joule          = Kilogram * ( Meter / Second )**2, &
    Erg            = Gram * ( Centimeter / Second )**2, &
    Bethe          = 1.0e51_DP * Erg, &
    ElectronVolt   = ElectronVoltMKS * Joule, &
    MeV            = 1.0e6_DP * ElectronVolt, &
    Kelvin         = BoltzmannConstantMKS / BoltzmannConstant * Joule, &
    Newton         = Joule / Meter, &
    Dyne           = Erg / Centimeter, &
    PlanckConstant = PlanckConstantMKS * Joule * Second, &
    AtomicMassUnit = Gram / AvogadroConstantMKS


  ! --- Units Displayed During Execution and for IO ---

  CHARACTER(16), PRIVATE, PARAMETER :: &
    DisplayLabel_Null              = '', &
    DisplayLabel_Length_L          = 'km', &
    DisplayLabel_Length_A          = 'rad', &
    DisplayLabel_Time              = 'ms', &
    DisplayLabel_Mass              = 'M_sun', &
    DisplayLabel_MassDensity       = 'g/cm^3', &
    DisplayLabel_ParticleDensity   = '1/cm^3', &
    DisplayLabel_Velocity_L        = 'km/s', &
    DisplayLabel_Velocity_A        = '1/s', &
    DisplayLabel_Momentum_L        = 'g cm/s', &
    DisplayLabel_Momentum_A        = 'g cm^2/s', &
    DisplayLabel_MomentumDensity_L = 'g/cm^2/s', &
    DisplayLabel_MomentumDensity_A = 'g/cm/s', &
    DisplayLabel_Energy            = 'MeV', &
    DisplayLabel_EnergyGlobal      = 'B', &
    DisplayLabel_EnergyDensity     = 'erg/cm^3', &
    DisplayLabel_Pressure          = 'erg/cm^3', &
    DisplayLabel_Temperature       = 'K'

  REAL(DP), PRIVATE, PARAMETER :: &
    DisplayUnit_Length_L          = Kilometer, &
    DisplayUnit_Length_A          = 1.0_DP, &
    DisplayUnit_Time              = Millisecond, &
    DisplayUnit_Mass              = SolarMass, &
    DisplayUnit_MassDensity       = Gram / Centimeter**3, &
    DisplayUnit_ParticleDensity   = 1.0_DP / Centimeter**3, &
    DisplayUnit_Velocity_L        = Kilometer / Second, &
    DisplayUnit_Velocity_A        = 1.0_DP / Second, &
    DisplayUnit_Momentum_L        = Gram * Centimeter / Second, &
    DisplayUnit_Momentum_A        = Gram * Centimeter**2 / Second, &
    DisplayUnit_MomentumDensity_L = Gram / Centimeter**2 / Second, &
    DisplayUnit_MomentumDensity_A = Gram / Centimeter / Second, &
    DisplayUnit_Energy            = MeV, &
    DisplayUnit_EnergyGlobal      = Bethe, &
    DisplayUnit_EnergyDensity     = Erg / Centimeter**3, &
    DisplayUnit_Pressure          = Erg / Centimeter**3, &
    DisplayUnit_Temperature       = Kelvin

  TYPE, PRIVATE :: UnitsDisplayType
    LOGICAL  :: &
      Active = .FALSE.
    CHARACTER(16) :: &
      LengthX1Label          = DisplayLabel_Null, &
      LengthX2Label          = DisplayLabel_Null, &
      LengthX3Label          = DisplayLabel_Null, &
      TimeLabel              = DisplayLabel_Null, &
      MassLabel              = DisplayLabel_Null, &
      MassDensityLabel       = DisplayLabel_Null, &
      ParticleDensityLabel   = DisplayLabel_Null, &
      VelocityX1Label        = DisplayLabel_Null, &
      VelocityX2Label        = DisplayLabel_Null, &
      VelocityX3Label        = DisplayLabel_Null, &
      MomentumX1Label        = DisplayLabel_Null, &
      MomentumX2Label        = DisplayLabel_Null, &
      MomentumX3Label        = DisplayLabel_Null, &
      MomentumDensityX1Label = DisplayLabel_Null, &
      MomentumDensityX2Label = DisplayLabel_Null, &
      MomentumDensityX3Label = DisplayLabel_Null, &
      EnergyLabel            = DisplayLabel_Null, &
      EnergyGlobalLabel      = DisplayLabel_Null, &
      EnergyDensityLabel     = DisplayLabel_Null, &
      PressureLabel          = DisplayLabel_Null, &
      TemperatureLabel       = DisplayLabel_Null
    REAL(DP) :: &
      LengthX1Unit          = 1.0_DP, &
      LengthX2Unit          = 1.0_DP, &
      LengthX3Unit          = 1.0_DP, &
      TimeUnit              = 1.0_DP, &
      MassUnit              = 1.0_DP, &
      MassDensityUnit       = 1.0_DP, &
      ParticleDensityUnit   = 1.0_DP, &
      VelocityX1Unit        = 1.0_DP, &
      VelocityX2Unit        = 1.0_DP, &
      VelocityX3Unit        = 1.0_DP, &
      MomentumX1Unit        = 1.0_DP, &
      MomentumX2Unit        = 1.0_DP, &
      MomentumX3Unit        = 1.0_DP, &
      MomentumDensityX1Unit = 1.0_DP, &
      MomentumDensityX2Unit = 1.0_DP, &
      MomentumDensityX3Unit = 1.0_DP, &
      EnergyUnit            = 1.0_DP, &
      EnergyGlobalUnit      = 1.0_DP, &
      EnergyDensityUnit     = 1.0_DP, &
      PressureUnit          = 1.0_DP, &
      TemperatureUnit       = 1.0_DP
  END type UnitsDisplayType

  TYPE(UnitsDisplayType), PUBLIC :: UnitsDisplay

  PUBLIC :: ActivateUnitsDisplay
  PUBLIC :: DescribeUnitsDisplay

CONTAINS


  SUBROUTINE ActivateUnitsDisplay( CoordinateSystem_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL ::  CoordinateSystem_Option

    CHARACTER(LEN=16)  :: CoordinateSystem

    CoordinateSystem = 'CARTESIAN'
    IF( PRESENT( CoordinateSystem_Option ) ) &
      CoordinateSystem = TRIM( CoordinateSystem_Option )

    UnitsActive = .TRUE.

    UnitsDisplay % Active = .TRUE.

    UnitsDisplay % TimeLabel            = DisplayLabel_Time
    UnitsDisplay % MassLabel            = DisplayLabel_Mass
    UnitsDisplay % MassDensityLabel     = DisplayLabel_MassDensity
    UnitsDisplay % ParticleDensityLabel = DisplayLabel_ParticleDensity

    SELECT CASE( TRIM( CoordinateSystem ) )

      CASE( 'CARTESIAN' )

        UnitsDisplay % LengthX1Label = DisplayLabel_Length_L
        UnitsDisplay % LengthX2Label = DisplayLabel_Length_L
        UnitsDisplay % LengthX3Label = DisplayLabel_Length_L

        UnitsDisplay % VelocityX1Label = DisplayLabel_Velocity_L
        UnitsDisplay % VelocityX2Label = DisplayLabel_Velocity_L
        UnitsDisplay % VelocityX3Label = DisplayLabel_Velocity_L

        UnitsDisplay % MomentumX1Label = DisplayLabel_Momentum_L
        UnitsDisplay % MomentumX2Label = DisplayLabel_Momentum_L
        UnitsDisplay % MomentumX3Label = DisplayLabel_Momentum_L

        UnitsDisplay % MomentumDensityX1Label = DisplayLabel_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Label = DisplayLabel_MomentumDensity_L
        UnitsDisplay % MomentumDensityX3Label = DisplayLabel_MomentumDensity_L

      CASE( 'CYLINDRICAL' )

        UnitsDisplay % LengthX1Label = DisplayLabel_Length_L
        UnitsDisplay % LengthX2Label = DisplayLabel_Length_L
        UnitsDisplay % LengthX3Label = DisplayLabel_Length_A

        UnitsDisplay % VelocityX1Label = DisplayLabel_Velocity_L
        UnitsDisplay % VelocityX2Label = DisplayLabel_Velocity_L
        UnitsDisplay % VelocityX3Label = DisplayLabel_Velocity_A

        UnitsDisplay % MomentumX1Label = DisplayLabel_Momentum_L
        UnitsDisplay % MomentumX2Label = DisplayLabel_Momentum_L
        UnitsDisplay % MomentumX3Label = DisplayLabel_Momentum_A

        UnitsDisplay % MomentumDensityX1Label = DisplayLabel_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Label = DisplayLabel_MomentumDensity_L
        UnitsDisplay % MomentumDensityX3Label = DisplayLabel_MomentumDensity_A

      CASE( 'SPHERICAL' )

        UnitsDisplay % LengthX1Label = DisplayLabel_Length_L
        UnitsDisplay % LengthX2Label = DisplayLabel_Length_A
        UnitsDisplay % LengthX3Label = DisplayLabel_Length_A

        UnitsDisplay % VelocityX1Label = DisplayLabel_Velocity_L
        UnitsDisplay % VelocityX2Label = DisplayLabel_Velocity_A
        UnitsDisplay % VelocityX3Label = DisplayLabel_Velocity_A

        UnitsDisplay % MomentumX1Label = DisplayLabel_Momentum_L
        UnitsDisplay % MomentumX2Label = DisplayLabel_Momentum_A
        UnitsDisplay % MomentumX3Label = DisplayLabel_Momentum_A

        UnitsDisplay % MomentumDensityX1Label = DisplayLabel_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Label = DisplayLabel_MomentumDensity_A
        UnitsDisplay % MomentumDensityX3Label = DisplayLabel_MomentumDensity_A

      CASE DEFAULT

        WRITE(*,*) 'Invalid choice for coodinate system: ', &
                   TRIM( CoordinateSystem )

    END SELECT

    UnitsDisplay % EnergyLabel          = DisplayLabel_Energy
    UnitsDisplay % EnergyGlobalLabel    = DisplayLabel_EnergyGlobal
    UnitsDisplay % EnergyDensityLabel   = DisplayLabel_EnergyDensity
    UnitsDisplay % PressureLabel        = DisplayLabel_Pressure
    UnitsDisplay % TemperatureLabel     = DisplayLabel_Temperature

    UnitsDisplay % TimeUnit            = DisplayUnit_Time
    UnitsDisplay % MassUnit            = DisplayUnit_Mass
    UnitsDisplay % MassDensityUnit     = DisplayUnit_MassDensity
    UnitsDisplay % ParticleDensityUnit = DisplayUnit_ParticleDensity

    SELECT CASE( TRIM( CoordinateSystem ) )

      CASE( 'CARTESIAN' )

        UnitsDisplay % LengthX1Unit = DisplayUnit_Length_L
        UnitsDisplay % LengthX2Unit = DisplayUnit_Length_L
        UnitsDisplay % LengthX3Unit = DisplayUnit_Length_L

        UnitsDisplay % VelocityX1Unit = DisplayUnit_Velocity_L
        UnitsDisplay % VelocityX2Unit = DisplayUnit_Velocity_L
        UnitsDisplay % VelocityX3Unit = DisplayUnit_Velocity_L

        UnitsDisplay % MomentumX1Unit = DisplayUnit_Momentum_L
        UnitsDisplay % MomentumX2Unit = DisplayUnit_Momentum_L
        UnitsDisplay % MomentumX3Unit = DisplayUnit_Momentum_L

        UnitsDisplay % MomentumDensityX1Unit = DisplayUnit_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Unit = DisplayUnit_MomentumDensity_L
        UnitsDisplay % MomentumDensityX3Unit = DisplayUnit_MomentumDensity_L

      CASE( 'CYLINDRICAL' )

        UnitsDisplay % LengthX1Unit = DisplayUnit_Length_L
        UnitsDisplay % LengthX2Unit = DisplayUnit_Length_L
        UnitsDisplay % LengthX3Unit = DisplayUnit_Length_A

        UnitsDisplay % VelocityX1Unit = DisplayUnit_Velocity_L
        UnitsDisplay % VelocityX2Unit = DisplayUnit_Velocity_L
        UnitsDisplay % VelocityX3Unit = DisplayUnit_Velocity_A

        UnitsDisplay % MomentumX1Unit = DisplayUnit_Momentum_L
        UnitsDisplay % MomentumX2Unit = DisplayUnit_Momentum_L
        UnitsDisplay % MomentumX3Unit = DisplayUnit_Momentum_A

        UnitsDisplay % MomentumDensityX1Unit = DisplayUnit_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Unit = DisplayUnit_MomentumDensity_L
        UnitsDisplay % MomentumDensityX3Unit = DisplayUnit_MomentumDensity_A

      CASE( 'SPHERICAL' )

        UnitsDisplay % LengthX1Unit = DisplayUnit_Length_L
        UnitsDisplay % LengthX2Unit = DisplayUnit_Length_A
        UnitsDisplay % LengthX3Unit = DisplayUnit_Length_A

        UnitsDisplay % VelocityX1Unit = DisplayUnit_Velocity_L
        UnitsDisplay % VelocityX2Unit = DisplayUnit_Velocity_A
        UnitsDisplay % VelocityX3Unit = DisplayUnit_Velocity_A

        UnitsDisplay % MomentumX1Unit = DisplayUnit_Momentum_L
        UnitsDisplay % MomentumX2Unit = DisplayUnit_Momentum_A
        UnitsDisplay % MomentumX3Unit = DisplayUnit_Momentum_A

        UnitsDisplay % MomentumDensityX1Unit = DisplayUnit_MomentumDensity_L
        UnitsDisplay % MomentumDensityX2Unit = DisplayUnit_MomentumDensity_A
        UnitsDisplay % MomentumDensityX3Unit = DisplayUnit_MomentumDensity_A

      CASE DEFAULT

        WRITE(*,*) 'Invalid choice for coodinate system: ', &
                   TRIM( CoordinateSystem )

    END SELECT

    UnitsDisplay % EnergyUnit          = DisplayUnit_Energy
    UnitsDisplay % EnergyGlobalUnit    = DisplayUnit_EnergyGlobal
    UnitsDisplay % EnergyDensityUnit   = DisplayUnit_EnergyDensity
    UnitsDisplay % PressureUnit        = DisplayUnit_Pressure
    UnitsDisplay % TemperatureUnit     = DisplayUnit_Temperature

  END SUBROUTINE ActivateUnitsDisplay


  SUBROUTINE DescribeUnitsDisplay

    WRITE(*,*)
    WRITE(*,'(A5,A25,L2)') &
      '', 'Units Activation Status =', UnitsDisplay % Active
    WRITE(*,'(A5,A32)') &
      '', '--------------------------------'
    WRITE(*,*)
    WRITE(*,'(A7,A29,A)') &
      '', 'Length (X1) Units: ', &
      TRIM( UnitsDisplay % LengthX1Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Length (X2) Units: ', &
      TRIM( UnitsDisplay % LengthX2Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Length (X3) Units: ', &
      TRIM( UnitsDisplay % LengthX3Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Time Units: ', &
      TRIM( UnitsDisplay % TimeLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Mass Units: ', &
      TRIM( UnitsDisplay % MassLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Mass Density Units: ', &
      TRIM( UnitsDisplay % MassDensityLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Particle Density Units: ', &
      TRIM( UnitsDisplay % ParticleDensityLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Velocity (X1) Units: ', &
      TRIM( UnitsDisplay % VelocityX1Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Velocity (X2) Units: ', &
      TRIM( UnitsDisplay % VelocityX2Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Velocity (X3) Units: ', &
      TRIM( UnitsDisplay % VelocityX3Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum (X1) Units: ', &
      TRIM( UnitsDisplay % MomentumX1Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum (X2) Units: ', &
      TRIM( UnitsDisplay % MomentumX2Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum (X3) Units: ', &
      TRIM( UnitsDisplay % MomentumX3Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum Density (X1) Units: ', &
      TRIM( UnitsDisplay % MomentumDensityX1Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum Density (X2) Units: ', &
      TRIM( UnitsDisplay % MomentumDensityX2Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Momentum Density (X3) Units: ', &
      TRIM( UnitsDisplay % MomentumDensityX3Label )
    WRITE(*,'(A7,A29,A)') &
      '', 'Energy Units: ', &
      TRIM( UnitsDisplay % EnergyLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Energy Global Units: ', &
      TRIM( UnitsDisplay % EnergyGlobalLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Energy Density Units: ', &
      TRIM( UnitsDisplay % EnergyDensityLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Pressure Units: ', &
      TRIM( UnitsDisplay % PressureLabel )
    WRITE(*,'(A7,A29,A)') &
      '', 'Temperature Units: ', &
      TRIM( UnitsDisplay % TemperatureLabel )
    WRITE(*,*)

  END SUBROUTINE DescribeUnitsDisplay


END MODULE UnitsModule
