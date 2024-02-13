#!/usr/bin/env python3


#+101+##############################################!
#                                                   !
#   LoadYahilProfile                                !
#                                                   !
####################################################!
def LoadYahilProfile( FileName ):

    from os.path import isfile
    import numpy as np
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





#+102+##############################################!
#                                                   !
#   GetYahilValues                                  !
#                                                   !
####################################################!
def GetYahilValues( rlocs, kappa, gamma, Time, YahilProfile ):

    import numpy as np
    
    X_Col = 0
    D_Col = 1
    V_Col = 2
    
    xlocs = [-1.0,+1.0]

    GravConstG = 1.0
    Erg = 8.2611082525066313E-052
    Gram = 7.4247138240457958E-031
    Millisecond = 299792.45799999998
    Second = 1000*Millisecond
    Centimeter = 1E-2
    Meter = 100*Centimeter
    Kilometer = 1000*Meter
    
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








#+201+##############################################!
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
     
     
     
     
     
            
#+202+##############################################!
#                                                   !
#   MapToXSpace                                     !
#                                                   !
####################################################!
def MapToXSpace(r, rleft, rright, xlocs ):

    return ((xlocs[1]-xlocs[0])/(rright-rleft))*(r-rleft) + xlocs[0]





#+201+##############################################!
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

