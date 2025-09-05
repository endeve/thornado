PROGRAM main

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    iCF_D , &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E , &
    iCF_Ne, &
    iPF_D , &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E , &
    iPF_Ne, &
    iAF_P, &
    ShortNamesPF
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic, &
    ComputePrimitive_Euler_Relativistic

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(amrex_parmparse) :: PP

  INTEGER  :: ITERATION, iErr
  REAL(DP) :: PF_D(2), PF_V1(2), PF_V2(2), PF_V3(2), PF_E(2), PF_Ne(2), AF_P, &
              CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
              GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33

  CALL InitializeProgram

  ! Initialize primitives

  CALL amrex_parmparse_build( PP )
    CALL PP % get( 'PF_D'       , PF_D (1)       )
    CALL PP % get( 'PF_V1'      , PF_V1(1)       )
    CALL PP % get( 'PF_V2'      , PF_V2(1)       )
    CALL PP % get( 'PF_V3'      , PF_V3(1)       )
    CALL PP % get( 'PF_E'       , PF_E (1)       )
    CALL PP % get( 'PF_Ne'      , PF_Ne(1)       )
    CALL PP % get( 'GF_Gm_dd_11', GF_Gm_dd_11 )
    CALL PP % get( 'GF_Gm_dd_22', GF_Gm_dd_22 )
    CALL PP % get( 'GF_Gm_dd_33', GF_Gm_dd_33 )
  CALL amrex_parmparse_destroy( PP )

  AF_P = ( Gamma_IDEAL - One ) * PF_E(1)

  CALL ComputeConserved_Euler_Relativistic &
    ( PF_D(1), PF_V1(1), PF_V2(1), PF_V3(1), PF_E(1), PF_Ne(1), &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
      AF_P )

  iErr = 0
  CALL ComputePrimitive_Euler_Relativistic &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D(2), PF_V1(2), PF_V2(2), PF_V3(2), PF_E(2), PF_Ne(2), &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
      ITERATION_Option = ITERATION, iErr_Option = iErr )

  WRITE(*,'(3A25)') 'Initial', 'Final', 'Difference'
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_D ), PF_D (1), PF_D (2), PF_D (2) - PF_D (1)
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_V1), PF_V1(1), PF_V1(2), PF_V1(2) - PF_V1(1)
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_V2), PF_V2(1), PF_V2(2), PF_V2(2) - PF_V2(1)
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_V3), PF_V3(1), PF_V3(2), PF_V3(2) - PF_V3(1)
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_E ), PF_E (1), PF_E (2), PF_E (2) - PF_E (1)
  WRITE(*,'(A6,3ES25.16E3)') ShortNamesPF(iPF_Ne), PF_Ne(1), PF_Ne(2), PF_Ne(2) - PF_Ne(1)
  WRITE(*,'(A,I4.4)') 'Number of Iterations = ', ITERATION
  WRITE(*,'(A,I2.2)') 'iErr = ', iErr

  CALL FinalizeProgram

END PROGRAM main
