MODULE MF_ErrorModule

  ! --- thornado Modules ---

  USE UtilitiesModule, ONLY: &
    thornado_abort

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DescribeError_MF

CONTAINS


  SUBROUTINE DescribeError_MF &
    ( iErr, Message_Option, Int_Option, Real_Option, Char_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option
    INTEGER,          INTENT(in), OPTIONAL :: Int_Option(:)
    REAL(DP),         INTENT(in), OPTIONAL :: Real_Option(:)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Char_Option(:)

    CHARACTER(LEN=128) :: Message

    Message = ''
    IF( PRESENT( Message_Option ) ) &
      Message = TRIM( Message_Option )

    SELECT CASE( iErr )

      CASE( 000 )

        RETURN

      CASE( 101 )

        CALL WriteHeader( 'InputParsingModule', 'InitializeParameters' )
        WRITE(*,'(2x,A)') &
          'iCycleW and dt_wrt cannot both be greater than zero.'
        WRITE(*,*)
        WRITE(*,'(2x,A9,I9.8)')      'iCycleW: ', Int_Option(1)
        WRITE(*,'(2x,A9,ES24.16E3)') 'dt_wrt: ', Real_Option(1)
        WRITE(*,*)

        CALL thornado_abort

      CASE( 102 )

        CALL WriteHeader( 'InputParsingModule', 'InitializeParameters' )
        WRITE(*,'(2x,A)') &
          'iCycleChk and dt_chk cannot both be greater than zero.'
        WRITE(*,*)
        WRITE(*,'(2x,A11,I9.8)')      'iCycleChk: ', Int_Option(1)
        WRITE(*,'(2x,A11,ES24.16E3)') 'dt_chk: ', Real_Option(1)
        WRITE(*,*)

        CALL thornado_abort

      CASE( 103 )

        CALL WriteHeader( 'InputParsingModule', 'InitializeParameters' )
        WRITE(*,'(2x,A,1x,I3.2)') 'Invalid CoordSys:', Int_Option(1)
        WRITE(*,*)
        WRITE(*,'(2x,A)') 'Valid choices'
        WRITE(*,'(2x,A)') '-------------'
        WRITE(*,'(2x,A)') '  0 (CARTESIAN)'
        WRITE(*,'(2x,A)') '  1 (CYLINDRICAL)'
        WRITE(*,'(2x,A)') '  2 (SPHERICAL)'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 104 )

        CALL WriteHeader( 'InputParsingModule', 'InitializeParameters' )
        WRITE(*,'(2x,A)') 'thornado nDimsX different from AMReX amrex_spacedim.'
        WRITE(*,'(2x,A,I2.2)') 'thornado: nDimsX: ', Int_Option(1)
        WRITE(*,'(2x,A,I2.2)') 'amrex:    nDimsX: ', Int_Option(2)
        WRITE(*,'(2x,A)') 'Check DIM parameter in GNUmakefile.'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 105 )

        CALL WriteHeader( 'MF_TimeSteppingModule', 'Initialize_IMEX_RK' )
        WRITE(*,'(2x,A,A)') &
          'Unknown Time Stepping Scheme: ', TRIM( Char_Option(1) )
        WRITE(*,*)
        WRITE(*,'(2x,A)') 'Available Options'
        WRITE(*,'(2x,A)') '-----------------'
        WRITE(*,*)
        WRITE(*,'(2x,A)') '  BackwardEuler'
        WRITE(*,'(2x,A)') '  SSPRK1'
        WRITE(*,'(2x,A)') '  SSPRK2'
        WRITE(*,'(2x,A)') '  SSPRK3'
        WRITE(*,'(2x,A)') '  IMEX_ARS_111'
        WRITE(*,'(2x,A)') '  IMEX_PDARS'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 201 )

        CALL WriteHeader( 'MF_InitializationModule', 'InitializeFields_MF' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  Advection1D'
        WRITE(*,'(2x,A)')   '  RiemannProblem1D'
        WRITE(*,'(2x,A)')   '  RiemannProblem2D'
        WRITE(*,'(2x,A)')   '  Advection2D'
        WRITE(*,'(2x,A)')   '  KelvinHelmholtz2D'
        WRITE(*,'(2x,A)')   '  Advection3D'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 202 )

        CALL WriteHeader( 'MF_InitializationModule', &
                          'InitializeFields_Advection_1D' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  SineWaveX1'
        WRITE(*,'(2x,A)')   '  Gaussian'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 203 )

        CALL WriteHeader( 'MF_InitializationModule', &
                          'InitializeFields_RiemannProblem1D' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  Sod'
        WRITE(*,'(2x,A)')   '  SphericalSod'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 204 )

        CALL WriteHeader( 'MF_InitializationModule', &
                          'InitializeFields_Advection2D' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  SineWaveX1'
        WRITE(*,'(2x,A)')   '  Gaussian'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 205 )

        CALL WriteHeader( 'MF_InitializationModule', &
                          'InitializeFields_Advection3D' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  SineWaveX1'
        WRITE(*,'(2x,A)')   '  Gaussian'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 301 )

        CALL WriteHeader( 'MF_InitializationModule', &
                          'InitializeFields_MF' )
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  SineWaveStreaming'
        WRITE(*,*)

        CALL thornado_abort

      CASE( 901 )

        CALL WriteHeader( 'InitializationModule', &
                          'InitializeProgram' )
        WRITE(*,'(2x,A)') 'Boundary condition exponent < 0'
        WRITE(*,'(2x,A)') TRIM( Message )

        CALL thornado_abort

      CASE( 902 )

        WRITE(*,*)
        CALL WriteHeader( 'MF_GravitySolutionModule_XCFC_Poseidon', &
                          'InitializeMetric_Euler_MF_Poseidon' )
        WRITE(*,'(2x,A)') 'Metric fields did not converge'
        WRITE(*,'(2x,A,ES24.16E3)') 'MAX( LaxLF, MaxCF ): ', Real_Option(1)
        WRITE(*,'(2x,A)') TRIM( Message )

        CALL thornado_abort

      CASE( 903 )

        WRITE(*,*)
        CALL WriteHeader( 'MF_GravitySolutionModule_XCFC_Poseidon', &
                          'InitializeMetric_TwoMoment_MF_Poseidon' )
        WRITE(*,'(2x,A)') 'Metric fields did not converge'
        WRITE(*,'(2x,A,ES24.16E3)') 'MAX( LaxLF, MaxCF ): ', Real_Option(1)
        WRITE(*,'(2x,A)') TRIM( Message )

        CALL thornado_abort

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') TRIM( Message )
          WRITE(*,*)

          CALL thornado_abort

    END SELECT

  END SUBROUTINE DescribeError_MF


  SUBROUTINE WriteHeader( ModuleName, SubroutineName )

    CHARACTER(*), INTENT(in) :: ModuleName, SubroutineName

    WRITE(*,*)
    WRITE(*,'(2x,A,A)') '-----------'
    WRITE(*,'(2x,A,A)') 'FATAL ERROR'
    WRITE(*,'(2x,A,A)') '-----------'
    WRITE(*,*)
    WRITE(*,'(2x,A,A)') &
      'MODULE: ', TRIM( ModuleName )
    WRITE(*,'(2x,A,A)') &
      'SUBROUTINE: ', TRIM( SubroutineName )
    WRITE(*,*)

  END SUBROUTINE WriteHeader


END MODULE MF_ErrorModule
