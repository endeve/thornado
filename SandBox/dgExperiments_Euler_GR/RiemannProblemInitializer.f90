MODULE RiemannProblemInitializer

  USE KindModule, ONLY: &
    DP, Pi, TwoPi

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: RiemannProblemChoice

CONTAINS

  SUBROUTINE RiemannProblemChoice &
               ( D_L, V_L, P_L, D_R, V_R, P_R, &
                 xL, xR, x_D, K, t, t_end, CFL, Gamma, bcX, CS, iRP )
    
    REAL(DP), INTENT(out) :: D_L, V_L(3), P_L, D_R, V_R(3), P_R
    REAL(DP), INTENT(out) :: xL(3), xR(3), x_D, t, t_end, CFL, Gamma
    INTEGER,  INTENT(out) :: bcX(3)
    INTEGER,  INTENT(in)  :: iRP
    INTEGER,  INTENT(out) :: K
    CHARACTER( LEN = 11 ) , INTENT(out) :: CS

    SELECT CASE( iRP )

      CASE( 0 )

        WRITE(*,*) 'Sods Shock Tube'

        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0_DP

        D_R = 0.125_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 0.1_DP

        t     = 0.0_DP
        t_end = 0.4_DP
        CFL   = 0.1_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 128
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 1 )

        WRITE(*,*) 'Blast Wave 1 (Del Zanna & Bucciantini (2002))'
         
        D_L = 10.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.33d+1

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-6

        t     = 0.0_DP
        t_end = 4.0d-1
        CFL   = 0.5_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 2 )

        WRITE(*,*) 'Blast Wave 2 (Del Zanna & Bucciantini (2002))'
         
        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+3

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-2

        t     = 0.0_DP
        t_end = 3.5d-1
        CFL   = 0.5_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 3 )

        WRITE(*,*) 'Strong Blast Wave (Qin et al., (2016)))'
         
        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+4

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-6

        t     = 0.0_DP
        t_end = 3.0d-1
        CFL   = 0.5_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 400
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 4 )

        WRITE(*,*) 'Shock Reflection (Del Zanna & Bucciantini (2002))'
         
        D_L = 1.0_DP
        V_L = [  0.99999_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d-2

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-2

        t     = 0.0_DP
        t_end = 7.0d-1
        CFL   = 0.1_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 1.0_DP
        K     = 250
        Gamma = 5.0_DP / 3.0_DP

        bcX = [ 3, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 5 )

        WRITE(*,*) 'Contact Discontinuity 1 (Mignone & Bodo (2005))'
         
        D_L = 1.0_DP
        V_L = [ 0.9_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+0

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d+1

        t     = 0.0_DP
        t_end = 4.0d-1
        CFL   = 0.8_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 100
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 6 )

        WRITE(*,*) 'Contact Discontinuity 2 (Mignone & Bodo (2005))'
         
        D_L = 1.0d+0
        V_L = [ -0.6d+0, 0.0d+0, 0.0d+0 ]
        P_L = 1.0d+1

        D_R = 1.0d+1
        V_R = [ 0.5d+0, 0.0d+0, 0.0d+0 ]
        P_R = 2.0d+1

        t     = 0.0d+0
        t_end = 4.0d+0
        CFL   = 0.8d+0
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5d+0
        K     = 100
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 7 )

        WRITE(*,*) 'Non-zero transverse velocity 1 (Qin et al., (2016))'
         
        D_L = 1.0d+0
        V_L = [ 0.0d+0, 0.9d+0, 0.0d+0 ]
        P_L = 1.0d+3

        D_R = 1.0d+0
        V_R = [ 0.0d+0, 0.9d+0, 0.0d+0 ]
        P_R = 1.0d-2

        t     = 0.0d+0
        t_end = 6.0d-1
        CFL   = 0.15d+0
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5d+0
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

      CASE( 8 )

        WRITE(*,*) 'Smooth problem: advected sine wave'

        ! --- Fluid variables here are dummies, they are set in
        !     InitializationModule_GR.f90 ---
        D_L = 999.0_DP
        V_L = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_L = 999.0_DP

        D_R = 999.0_DP
        V_R = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_R = 999.0_DP

        t     = 0.0_DP
        CFL   = 0.1_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 128
        t_end = 1.0d1! * CFL * xR / ( 1.0d0 * K )
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 1, 0, 0 ]        

        CS = 'CARTESIAN'

      CASE( 9 )

        WRITE(*,*) 'Contact discontinuity, top-hat'

        ! --- Fluid variables here are dummies, they are set in
        !     InitializationModule.f90 ---
        D_L = 999.0_DP
        V_L = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_L = 999.0_DP

        D_R = 999.0_DP
        V_R = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_R = 999.0_DP

        t     = 0.0_DP
        CFL   = 0.1_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 128
        t_end = 1.0d1
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 1, 0, 0 ]

        CS = 'CARTESIAN'
        
      CASE( 10 )

        WRITE(*,*) 'Perturbed shock tube'

        ! --- Fluid variables here are dummies, they are set in
        !     InitializationModule.f90 ---
        D_L = 999.0_DP
        V_L = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_L = 999.0_DP

        D_R = 999.0_DP
        V_R = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_R = 999.0_DP

        t     = 0.0_DP
        CFL   = 0.1_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 128
        t_end = 3.5d-1
        Gamma = 5.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]        

        CS = 'CARTESIAN'

      CASE( 11 )

        WRITE(*,*) 'Acoustics problem in spherical symmetry'

        ! --- Fluid variables here are dummies, they are set in
        !     InitializationModule.f90 ---
        D_L = 999.0_DP
        V_L = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_L = 999.0_DP

        D_R = 999.0_DP
        V_R = [ 999.0_DP, 999.0_DP, 999.0_DP ]
        P_R = 999.0_DP

        t     = 0.0d+0
        t_end = 0.5d+0
        CFL   = 0.5d+0
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, Pi, TwoPi ]
        x_D   = 0.5d+0
        K     = 64
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'SPHERICAL'

      CASE( 12 )

        WRITE(*,*) 'Sods Shock Tube, Spherical'

        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
!        P_L = 1.0_DP
        P_L = 0.1_DP

!        D_R = 0.125_DP
        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 0.1_DP

        t     = 0.0_DP
        t_end = 0.5_DP
        CFL   = 0.1_DP
        xL    = [ 0.0001d0, 0.0d0, 0.0d0 ]
        xR    = [ 2.0d0, Pi, TwoPi ]
        x_D   = 1.0_DP
        K     = 128
        Gamma = 4.0_DP / 3.0_DP

        bcX = [ 2, 0, 0 ]

        CS = 'SPHERICAL'

      CASE( 13 )

        WRITE(*,*) &
              'Stationary Contact Discontinuity (Test 5) (Liska & Wendroff (2003))'
         
        D_L = 1.4_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+0

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d+0

        t     = 0.0_DP
        t_end = 2.0d0
        CFL   = 0.8_DP
        xL    = [ 0.0d0, 0.0d0, 0.0d0 ]
        xR    = [ 1.0d0, 1.0d0, 1.0d0 ]
        x_D   = 0.5_DP
        K     = 128
        Gamma = 1.4_DP

        bcX = [ 2, 0, 0 ]

        CS = 'CARTESIAN'

     CASE DEFAULT

        WRITE(*,*) 'Invalid choice for iRP, exiting...'
        STOP

    END SELECT

  END SUBROUTINE RiemannProblemChoice

END MODULE RiemannProblemInitializer

