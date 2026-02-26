def _Approx_Fluid_CM_S3(field, data):
    rho = data['boxlib', 'PM_D' ]
    r   = data['boxlib', 'X1_C' ]
    V3  = data['boxlib', 'PM_V3']
    return 1.00e+05 * r * rho * V3

def _Approx_Fluid_CM_S3_Flux_X1(field, data):
    rho = data['boxlib', 'PM_D' ]
    r   = data['boxlib', 'X1_C' ]
    V1  = data['boxlib', 'PM_V1']
    V3  = data['boxlib', 'PM_V3']
    B1  = data['boxlib', 'CM_B1']
    B3  = data['boxlib', 'CM_B3']
    return 1.00e+05 * r * ( 1.00e+05 * rho * V1 * V3 - B1 * B3 )

def _Approx_Fluid_CM_S3_Flux_X2(field, data):
    rho = data['boxlib', 'PM_D' ]
    r   = data['boxlib', 'X1_C' ]
    V2  = data['boxlib', 'PM_V2']
    V3  = data['boxlib', 'PM_V3']
    B2  = data['boxlib', 'CM_B2']
    B3  = data['boxlib', 'CM_B3']
    return 1.00e+05 * r * ( 1.00e+05 * rho * V2 * V3 - B2 * B3 )
