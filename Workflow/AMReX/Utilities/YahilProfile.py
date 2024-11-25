#!/usr/bin/env python3
import numpy as np

#+101+##############################################!
#                                                   !
#   LoadYahilProfile                                !
#                                                   !
####################################################!
def LoadYahilProfile( FileName ):

    from os.path import isfile
    import csv
    
    
    if not ( isfile(FileName) ):
        
        print('File does not exist : '+FileName)
        
        return
    

    with open(FileName) as file:
        reader = csv.reader(file,delimiter=' ')
        
        row_count = sum(1 for row in file)
                
        YahilProfile = np.zeros((row_count,4))
        
        file.seek(0)
        next(file)
        i = 0
        for row in reader:
            i += 1
            YahilProfile[i] = row

    return YahilProfile





#+201+##############################################!
#                                                   !
#   GetYahilValues                                  !
#                                                   !
####################################################!
def GetYahilValues( rlocs, kappa, gamma, Time, YahilProfile ):

    X_Col = 0
    D_Col = 1
    V_Col = 2
    
    xlocs = [-1.0,+1.0]

#    GravConstG = 1.0
#    Erg = 8.2611082525066313E-052
#    Gram = 7.4247138240457958E-031
#    Millisecond = 299792.45799999998
#    Second = 1000*Millisecond
#    Centimeter = 1E-2
#    Meter = 100*Centimeter
#    Kilometer = 1000*Meter

    C_MKS = 2.99792458e8
    G_MKS = 6.673e-11
#    Erg = 8.2611082525066313E-052
    Gram = 1.0
    Kilogram = 1000.0*Gram
    
    Second = 1.0
    Millisecond = Second / 1000.0
    Centimeter = 1.0
    Meter = 100.0*Centimeter
    Kilometer = 1000.0*Meter
    
    C          = C_MKS*(Meter/Second)
    GravConstG = G_MKS*(Meter*Meter*Meter)/(Kilogram*Second*Second)
    Erg        = Gram*( Centimeter/Second)**2
    
    New_Time = 150 - Time
    t = New_Time*Millisecond
    
    kappa_wUnits = kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**gamma)

    D_Factor = 1.0/(GravConstG*t*t)

    V_Factor = np.sqrt(kappa_wUnits) \
             * GravConstG**((1.0 - gamma)/2.0) \
             *t**(1.0 - gamma )
    
    X_Factor = kappa_wUnits**(-0.5) \
         * GravConstG**((gamma -1.0)/2.0) \
         * t**(gamma - 2.0)
         
    R_Factor = np.sqrt(kappa_wUnits) \
        *GravConstG**((1.0-gamma)/2.0) \
        *t**(2.0-gamma)

    rlocs = rlocs*Kilometer
    rlen = len(rlocs)
    
    Density = np.zeros((rlen,1))
    Velocity = np.zeros((rlen,1))
    
    for r in range(rlen):
        
        xloc = X_Factor*rlocs[r]
        line = FindLine( xloc, YahilProfile[:,X_Col])
              
        x = MapToXSpace(xloc,YahilProfile[line,X_Col],YahilProfile[line+1,X_Col],xlocs)
        if x > 1:
            x = 1
        
        LagPoly_Values =  LagPoly(x, xlocs )
    
        Tmp_D = ( YahilProfile[line,D_Col]*LagPoly_Values[0] \
                  + YahilProfile[line+1,D_Col]*LagPoly_Values[1] ) \
                  * D_Factor

        Tmp_V = ( YahilProfile[line,V_Col]*LagPoly_Values[0] \
                   + YahilProfile[line+1,V_Col]*LagPoly_Values[1] ) \
                   * V_Factor
                   
        Density[r] = Tmp_D/(Gram/Centimeter**3)
        Velocity[r] = Tmp_V/(Kilometer/Second)
                   
    return Density,Velocity



#+202+##############################################!
#                                                   !
#   GetYahilPotential                               !
#                                                   !
####################################################!
def GetYahilPotential( kappa, gamma, Time, Profile ):

    
    X_Col = 0
    D_Col = 1
    V_Col = 2
    M_Col = 3
    
    xlocs = [-1.0,+1.0]


    C_MKS = 2.99792458e8
    G_MKS = 6.673e-11
    Gram = 1.0
    Kilogram = 1000.0*Gram
    
    Second = 1.0
    Millisecond = Second / 1000.0
    Centimeter = 1.0
    Meter = 100.0*Centimeter
    Kilometer = 1000.0*Meter
    
    C          = C_MKS*(Meter/Second)
    GravConstG = G_MKS*(Meter*Meter*Meter)/(Kilogram*Second*Second)
    Erg        = Gram*( Centimeter/Second)**2
    
    CSquare = C*C
    
    New_Time = 150 - Time
    t = New_Time*Millisecond
    
    kappa_wUnits = kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**gamma)
     
    R_Factor = np.sqrt(kappa_wUnits)            \
        * GravConstG**((1.0-gamma)/2.0)         \
        * (t**(2.0-gamma))
        
    M_Factor = kappa_wUnits**(1.50)             \
        * GravConstG**((1.0 - 3.0*gamma)/2.0)   \
        * (t**(4.0 - 3.0*gamma))


    NumLines = len(Profile)
    
    EnclosedMass = np.zeros((NumLines,1))
    R            = np.zeros((NumLines,1))
    for line in range(NumLines):
        R[line]            = R_Factor*Profile[line,X_Col]
        EnclosedMass[line] = M_Factor*Profile[line,M_Col]
        
        
    Potential = np.zeros((NumLines,1))
    
    i = NumLines-1
    Potential[i]=-GravConstG*EnclosedMass[i]/R[i]
    
    
    
    for i in reversed(range(1,NumLines-1)):
        Potential[i] = Potential[i+1]                   \
                     - GravConstG*EnclosedMass[i]       \
                     * (R[i+1] - R[i])                  \
                     / (R[i]*R[i])
      
    Potential[0] = Potential[1]                     \
                 - 3*GravConstG*EnclosedMass[1]     \
                 /(2*R[1])

    
    return R, Potential
    



#+203+##############################################!
#                                                   !
#   GetYahilPotential                               !
#                                                   !
####################################################!
def GetYahilMetric( rlocs, kappa, gamma, Time, Profile ):

    X_Col = 0
    xlocs = [-1.0,+1.0]
    
    C_MKS = 2.99792458e8
    G_MKS = 6.673e-11

    Gram = 1.0
    Kilogram = 1000.0*Gram
    
    Second = 1.0
    Millisecond = Second / 1000.0
    Centimeter = 1.0
    Meter = 100.0*Centimeter
    Kilometer = 1000.0*Meter
    
    C          = C_MKS*(Meter/Second)
    GravConstG = G_MKS*(Meter*Meter*Meter)/(Kilogram*Second*Second)
    Erg        = Gram*( Centimeter/Second)**2
    
    CSquare    = C*C
    
    
    
    R, Potential = GetYahilPotential( kappa, gamma, Time, Profile )
    
    rlocs = rlocs*Kilometer
    rlen = len(rlocs)
    
    kappa_wUnits = kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**gamma)
    
    t = (150-Time)*Millisecond
    X_Factor = kappa_wUnits**(-0.5) \
             * GravConstG**((gamma -1.0)/2.0) \
             * t**(gamma - 2.0)
    
    
    
    
    Psi   = np.zeros((rlen,1))
    Alpha = np.zeros((rlen,1))
    for r in range(rlen):
    
        xloc = X_Factor*rlocs[r]
        line = FindLine( xloc, Profile[:,X_Col])
        
        x = MapToXSpace(xloc,Profile[line,X_Col],Profile[line+1,X_Col],xlocs)
        if x > 1:
            x = 1
        
        LagPoly_Values =  LagPoly(x, xlocs )
    
        TmpPot = Potential[line]*LagPoly_Values[0]     \
               + Potential[line+1]*LagPoly_Values[1]

        
        Psi[r] = 1.0 - 0.5 * TmpPot/CSquare
        Alpha[r] = (1.0 + 0.5*TmpPot/CSquare)/Psi[r]

#        print(rlocs[r],TmpPot,CSquare)
    
    return Psi, Alpha




#+301+##############################################!
#                                                   !
#   FindLine                                        !
#                                                   !
####################################################!
def FindLine( x, xlocs ):

    xlen = len(xlocs)-1
    up = xlen
    down = 0
    
    while (up - down > 1):
        mid = (up + down)//2
        if ( (xlocs[xlen]>=xlocs[1]) == (x>=xlocs[mid]) ):
            down = mid
        else:
            up = mid
            
    if x == xlocs[0]:
        return 0
    elif x == xlocs[xlen]:
        return xlen-1
    elif down == xlen:
        return xlen-1
    else:
        return down
     
     
            
#+302+##############################################!
#                                                   !
#   MapToXSpace                                     !
#                                                   !
####################################################!
def MapToXSpace(r, rleft, rright, xlocs ):

    return ((xlocs[1]-xlocs[0])/(rright-rleft))*(r-rleft) + xlocs[0]





#+303+##############################################!
#                                                   !
#   LagPoly                                         !
#                                                   !
####################################################!
def LagPoly(x, xlocs ):

    xlen = len(xlocs)
    tmp = [1.0]*xlen
    
    for j in range(xlen):
        for i in range(xlen):
            if i != j:
                tmp[i] = tmp[i] * (x - xlocs[j])/(xlocs[i]-xlocs[j])
            
    return tmp










#+401+##############################################!
#                                                   !
#   CalcPotential                                   !
#                                                   !
####################################################!
def CalcPotential( X1_C, dX1, Density ):

    C_MKS = 2.99792458e8
    G_MKS = 6.673e-11

    Gram = 1.0
    Kilogram = 1000.0*Gram
    
    Second = 1.0
    Millisecond = Second / 1000.0
    Centimeter = 1.0
    Meter = 100.0*Centimeter
    Kilometer = 1000.0*Meter
    
    C          = C_MKS*(Meter/Second)
    GravConstG = G_MKS*(Meter*Meter*Meter)/(Kilogram*Second*Second)
    Erg        = Gram*( Centimeter/Second)**2
    
    CSquare    = C*C

    M = len(Density)
    N = len(dX1)
    
    if (M != N):
        print('M,N:',M,N)
        return -1
        
        
    rlocs = X1_C*Kilometer
    dr    = dX1*Kilometer
    
    EnclosedMass = np.zeros((M+1,1))
    EnclosedMass[0] = 0.0
    EnclosedMass[1] = 4.0/3.0 * np.pi * dr[0]**3 * Density[0]
    for shell in range(1,M):
        RO = rlocs[shell] + dr[shell]/2
        RI = rlocs[shell] - dr[shell]/2
    
        EnclosedMass[shell+1] = EnclosedMass[shell]             \
                             + 4.0/3.0 * np.pi * (RO**3 - RI**3)\
                             * Density[shell]

    
    Potential = np.zeros((M+1,1))
    i = M
    RO = rlocs[i-1] + dr[i-1]/2
    Potential[i]=-GravConstG*EnclosedMass[i]/RO
    
    for i in reversed(range(1,M)):
        RO = rlocs[i] + dr[i]/2
        RI = rlocs[i] - dr[i]/2
        Potential[i] = Potential[i+1]                   \
                     - GravConstG*EnclosedMass[i]       \
                     * (RO - RI)                        \
                     / (RI*RI)
      
    RO = rlocs[0] + dr[0]/2
    Potential[0] = Potential[1]                     \
                 - 3*GravConstG*EnclosedMass[0]     \
                 /(2*RO)
    
    
    
    Psi   = np.zeros((M,1))
    Alpha = np.zeros((M,1))
    for r in range(M):
    
        TmpPot = Potential[r]*0.5     \
               + Potential[r+1]*0.5

        Psi[r] = 1.0 - 0.5 * TmpPot/CSquare
        Alpha[r] = (1.0 + 0.5*TmpPot/CSquare)/Psi[r]
        
        
    return Psi, Alpha
