MODULE MF_GravitySolutionModule_Newtonian_Poseidon

  use mf_utilitiesmodule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX, &
    swX, &
    nX
  USE MeshModule, ONLY: &
    MeshX
  USE UnitsModule, ONLY: &
    Kilometer, &
    Kilogram, &
    Gram, &
    Meter, &
    Erg, &
    Second, &
    GravitationalConstant, &
    Centimeter
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iAF_Cs
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N, &
    iGF_SqrtGm, &
    nGF
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL

  USE MF_UtilitiesModule, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    amrex2thornado_X

  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask

  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q

  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    ComputeTotalBaryonMass

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    Half, &
    Zero, &
    DP, &
    Pi, &
    FourPi
  USE InputParsingModule, ONLY: &
    xR, &
    nLevels, &
    UseTiling, &
    ProgramName
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE AverageDownModule, ONLY: &
    AverageDown

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Initialization, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Interface_Boundary_Conditions, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Interface_Source_Input, ONLY: &
    Poseidon_Input_Sources
  USE Poseidon_Interface_Run, ONLY: &
    Poseidon_Run
  USE Poseidon_Interface_Return_Routines, ONLY: &
    Poseidon_Return_Newtonian_Potential
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_MF_Poseidon
  PUBLIC :: FinalizeGravitySolver_Newtonian_MF_Poseidon
  PUBLIC :: ComputeGravitationalPotential_Newtonian_MF_Poseidon
  PUBLIC :: ComputeEnclosedMass_MF

  INTEGER, PARAMETER :: iGS_D = 1
  INTEGER, PARAMETER :: nGS   = 1

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_MF_Poseidon( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') &
        'INFO: Gravity Solver (Poseidon, Newtonian)'
      WRITE(*,'(4x,A)') &
        '------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Only implemented for 1D spherical symmetry.'
      WRITE(*,*)

    END IF

    CALL Initialize_Poseidon &
           ( Source_NQ                    = nNodesX,          &
             Source_xL                    = [ -Half, +Half ], &
             Source_RQ_xlocs              = MeshX(1) % Nodes, &
             Source_TQ_xlocs              = MeshX(2) % Nodes, &
             Source_PQ_xlocs              = MeshX(3) % Nodes, &
             Source_Units                 = 'G',              &
             Source_Radial_Boundary_Units = 'km',             &
	     Newtonian_Mode_Option        = .TRUE. ,          &
             Verbose_Option               = .TRUE.,           &
             Print_Setup_Option           = .TRUE.,           &
             Print_Results_Option         = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_Newtonian_MF_Poseidon


  SUBROUTINE ComputeGravitationalPotential_Newtonian_MF_Poseidon( MF_uCF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    REAL(DP) :: Phi_r_outer
    INTEGER  :: iLevel

    TYPE(amrex_multifab) :: MF_uGS(0:nLevels-1)
    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1)

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
              ( MF_uMF(iLevel), MF_uCF(iLevel) % BA, &
                                MF_uCF(iLevel) % DM, nDOFX, 0 )

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uCF(iLevel) % BA, &
                               MF_uCF(iLevel) % DM, nGS * nDOFX, 0 )

      CALL MF_uGS(iLevel) % COPY &
             ( MF_uCF(iLevel), nDOFX*(iCF_D-1)+1, 1, nDOFX, 0 )

    END DO

call showvariablefrommultifab(mf_ugs,1)    

    CALL Poseidon_Input_Sources( MF_uGS )

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
  	       ( BC_Location_Input = 'I', &
    	       BC_Type_Input     = 'N', &
    	       BC_Value_Input    = Zero )

    CALL ComputeBoundaryValues( MF_uGF, MF_uGS, Phi_r_outer )

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( BC_Location_Input = 'O', &
    	       BC_Type_Input     = 'D', &
    	       BC_Value_Input    = Phi_r_outer )

   CALL Poseidon_Run()

   CALL Poseidon_Return_Newtonian_Potential( MF_uMF )

   DO iLevel = 0, nLevels-1

     CALL MF_uGF(iLevel) % COPY &
            ( MF_uMF(iLevel), 1, nDOFX*(iGF_Phi_N-1)+1, nDOFX, 0 )

     CALL amrex_multifab_destroy( MF_uGS(iLevel) )

     CALL amrex_multifab_destroy( MF_uMF(iLevel) )

   END DO

   CALL AverageDown( MF_uGF )

   CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )





#endif

  END SUBROUTINE ComputeGravitationalPotential_Newtonian_MF_Poseidon


  SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon


  SUBROUTINE ComputeBoundaryValues( MF_uGF, MF_uGS, Phi_r_outer )

      TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:)
      TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:)
      REAL(DP)            , INTENT(out) :: Phi_r_outer

      REAL(DP)              :: GravitationalMass
      

      TYPE(amrex_parmparse) :: PP
      REAL(DP) :: Phi_r, Phi_r_converted, P_c, rho_c, r_c, SolarRadius, R_star, Domain
      REAL(DP) :: Part1, Part2, Phi_r_2, K_polytrope, R_neutron_star, Outer_Domain, alpha, alpha2

      CALL amrex_parmparse_build( PP, 'thornado' )
        CALL PP % get( 'ProgramName', ProgramName )
      CALL amrex_parmparse_destroy(PP)

      SELECT CASE ( TRIM( ProgramName ) )

        CASE ( 'PoissonSolverTest_Newtonian' )

          CALL amrex_parmparse_build( PP, 'inputs' )
            CALL PP % get( 'rho_c', rho_c )
            CALL PP % get( 'R_star', R_star )
            CALL PP % get( 'SolarRadius', SolarRadius )
            CALL PP % get( 'r_c', r_c )
          CALL amrex_parmparse_destroy( PP )  

          rho_c       = rho_c * ( Gram / Centimeter**3 )
          SolarRadius = SolarRadius * Kilometer
          R_star      = R_star * SolarRadius
          r_c         = r_c * SolarRadius
          Domain      = 2.0_DP * R_star  

          Part1           = -4.0_DP * Pi * rho_c * ( r_c**3 / Domain )
          Part2           = ( R_star / r_c ) - atan( R_star / r_c )
          Phi_r           = Part1 * Part2
          Phi_r_converted = Phi_r / ( Erg / Gram )
          Phi_r_outer     = Phi_r  

          WRITE(*,'(4x,A)') &
            '------------------------------------------'
          WRITE(*,*)
          WRITE(*,'(4x,A)')'Using Outer Boundary Values from Pochik et al 2021 '
          WRITE(*,'(8x,A,ES15.6E3)')'Phi_r in CGS: ', Phi_r_converted
          WRITE(*,'(4x,A)') &
            '------------------------------------------'
          WRITE(*,*)

        CASE ( 'Hydrostatic_Polytrope' )

          CALL amrex_parmparse_build( PP, 'inputs' )
            CALL PP % get( 'rho_c', rho_c )
            CALL PP % get( 'R_neutron_star', R_neutron_star )
            CALL PP % get( 'Outer_Domain', Outer_Domain )
          CALL amrex_parmparse_destroy( PP )  

          rho_c           = rho_c * ( Gram / Centimeter**3 )
          R_neutron_star  = R_neutron_star * Kilometer
          Outer_Domain    = Outer_Domain * Kilometer
          K_polytrope     = ( 2.0_DP / Pi ) * Outer_Domain**2
          P_c             = K_polytrope * rho_c**Gamma_IDEAL
          alpha           = R_neutron_star / Pi
          Phi_r_2         = 4.0_DP * Pi * alpha**2 * rho_c * &
                              ( (- ( alpha * sin( Outer_Domain / alpha ) ) &
                              + Outer_Domain * cos( Outer_Domain / alpha ) ) &
                              / Outer_Domain ) 
          Phi_r_outer     = Phi_r_2
          Phi_r_converted = Phi_r_2 / ( Erg / Gram )  

          WRITE(*,'(4x,A)') &
            '------------------------------------------'
          WRITE(*,*)
          WRITE(*,'(4x,A)')'Using Outer Boundary Values for Polytrope Test '
          WRITE(*,'(8x,A,ES15.6E3)')'Phi_r in CGS: ', Phi_r_converted
          WRITE(*,'(8x,A,ES15.6E3)')'Sound Speed in Kilometer/Second: ', iAF_Cs/( Kilometer / Second)
          WRITE(*,'(4x,A)') &
            '------------------------------------------'

        CASE DEFAULT 

          CALL ComputeEnclosedMass_MF( MF_uGF, MF_uGS, GravitationalMass )

          Phi_r_outer  = - GravitationalMass / xR(1)
          PRINT *, Phi_r_outer

          WRITE(*,'(4x,A)') &
            '------------------------------------------'
          WRITE(*,*)
          WRITE(*,'(4x,A)')'Using Outer Boundary Values using Enclosed Mass '
          WRITE(*,'(8x,A,ES15.6E3)')'Phi_r in CGS: ', Phi_r_outer / ( Erg / Gram )
          WRITE(*,'(4x,A)') &
            '------------------------------------------'

    END SELECT

  END SUBROUTINE ComputeBoundaryValues


  SUBROUTINE ComputeEnclosedMass_MF( MF_uGF, MF_uGS, GravitationalMass )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:)
    REAL(DP)            , INTENT(out) :: GravitationalMass

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: GF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: GS(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLevel

    GravitationalMass = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, GF, GS, uGF, uGS, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 ) &
      !$OMP REDUCTION( +:GravitationalMass )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )



      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 GF )        

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGS ], &
                 GS )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, GF )

        CALL amrex2thornado_X &
               ( nGS, iX_B1, iX_E1, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

        CALL ComputeTotalBaryonMass &
                ( iX_B0, iX_E0, iX_B1, iX_E1, &
                  GF, &
                  GS(1:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),iGS_D), &
                  GravitationalMass, Mask_Option = uFM )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGS ], &
                 GS )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 GF )        

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

    CALL amrex_parallel_reduce_sum( GravitationalMass )

  END SUBROUTINE ComputeEnclosedMass_MF


END MODULE MF_GravitySolutionModule_Newtonian_Poseidon
