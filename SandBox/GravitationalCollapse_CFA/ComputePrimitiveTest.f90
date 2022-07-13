PROGRAM ComputePrimitiveTest

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE UnitsModule, ONLY: &
    Kilometer, &
    Gram, &
    Second, &
    Erg, &
    Centimeter
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState, &
    ComputeSpecificInternalEnergy
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_T
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER,  PARAMETER :: MaxIter = 10000000
  INTEGER             :: iErr
  REAL(DP)            :: U(nCF), P(nPF), G(nGF), rhot, Ye, MinE
  REAL(DP), PARAMETER :: UnitsD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitsS1   = Gram / Centimeter**2 / Second
  REAL(DP), PARAMETER :: UnitsE    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsNe   = One  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsGm11 = One
  REAL(DP), PARAMETER :: UnitsGm22 = Kilometer**2
  REAL(DP), PARAMETER :: UnitsGm33 = Kilometer**2
  REAL(DP), PARAMETER :: UnitsV1   = Kilometer / Second

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'wl-EOS-SFHo-25-50-100.h5' )

! bad
 U(iCF_D) =    7.0649187518464309E-011
 U(iCF_S1) =   -1.0090020859170580E-011
 U(iCF_S2) =    0.0000000000000000
 U(iCF_S3) =    0.0000000000000000
 U(iCF_E) =    4.6814468682533512E-012
 U(iCF_Ne) =    2.5320751005275095d043
 G(iGF_Gm_dd_11) =    1.3084305249726724
 G(iGF_Gm_dd_22) =    52098906.726469778
 G(iGF_Gm_dd_33) =    52098906.726469778

!! good
! U(iCF_D) =    7.1310685995766948E-011
! U(iCF_S1) =   -9.6043957705295021E-012
! U(iCF_S2) =    0.0000000000000000
! U(iCF_S3) =    0.0000000000000000
! U(iCF_E) =    4.7298870460777394E-012
! U(iCF_Ne) =    2.5552046053893353d+043
! G(iGF_Gm_dd_11) =    1.3072681090882421
! G(iGF_Gm_dd_22) =    52052367.494140975
! G(iGF_Gm_dd_33) =    52052367.494140975

  WRITE(*,*)
  WRITE(*,'(A)') 'Inputs:'
  WRITE(*,'(A,ES24.16E3)') 'CF_D  = ', U(iCF_D ) / UnitsD
  WRITE(*,'(A,ES24.16E3)') 'CF_S1 = ', U(iCF_S1) / UnitsS1
  WRITE(*,'(A,ES24.16E3)') 'CF_S2 = ', U(iCF_S2)
  WRITE(*,'(A,ES24.16E3)') 'CF_S3 = ', U(iCF_S3)
  WRITE(*,'(A,ES24.16E3)') 'CF_E  = ', U(iCF_E ) / UnitsE
  WRITE(*,'(A,ES24.16E3)') 'CF_Ne = ', U(iCF_Ne) / UnitsNe
  WRITE(*,'(A,ES24.16E3)') 'GF_g1 = ', G(iGF_Gm_dd_11) / UnitsGm11
  WRITE(*,'(A,ES24.16E3)') 'GF_g2 = ', G(iGF_Gm_dd_22) / UnitsGm22
  WRITE(*,'(A,ES24.16E3)') 'GF_g3 = ', G(iGF_Gm_dd_33) / UnitsGm33
  WRITE(*,*)

 rhot =    7.0375498702540617d-011
 Ye =   0.44187368273397742_DP
  CALL ComputeSpecificInternalEnergy( rhot, Min_T, Ye, MinE )
print*
print*,'From Table: ', MinE
print*,'From Code:  ', 5.9387860185921022d-002
print*
  iErr = 0

  CALL ComputePrimitive_Euler_Relativistic &
         ( U(iCF_D ), &
           U(iCF_S1), &
           U(iCF_S2), &
           U(iCF_S3), &
           U(iCF_E ), &
           U(iCF_Ne), &
           P(iPF_D ), &
           P(iPF_V1), &
           P(iPF_V2), &
           P(iPF_V3), &
           P(iPF_E ), &
           P(iPF_Ne), &
           G(iGF_Gm_dd_11), &
           G(iGF_Gm_dd_22), &
           G(iGF_Gm_dd_33), &
           iErr )

  IF( iErr .NE. 0 )THEN

    print*,'iErr = ', iErr
    stop 'fail'

  END IF

  WRITE(*,'(A)') 'Outputs:'
  WRITE(*,'(A,ES24.16E3)') 'PF_D  = ', P(iPF_D ) / UnitsD
  WRITE(*,'(A,ES24.16E3)') 'PF_V1 = ', P(iPF_V1) / UnitsV1
  WRITE(*,'(A,ES24.16E3)') 'PF_V2 = ', P(iPF_V2)
  WRITE(*,'(A,ES24.16E3)') 'PF_V3 = ', P(iPF_V3)
  WRITE(*,'(A,ES24.16E3)') 'PF_E  = ', P(iPF_E ) / UnitsE
  WRITE(*,'(A,ES24.16E3)') 'PF_Ne = ', P(iPF_Ne) / UnitsNe

  CALL FinalizeEquationOfState

END PROGRAM ComputePrimitiveTest
