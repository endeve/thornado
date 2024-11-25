#!/usr/local/bin/python3

# --- Import libraries ---
import warnings
warnings.simplefilter( action = 'ignore' , category = FutureWarning )

# --- Get rid of "Too many open files" error ---

import numpy as np
import h5py as h5

def ReadFields( PathToData, Snapshots ):


    GF_root = PathToData + 'thornado_Christoffel_000000'

    nFiles = len( Snapshots )


    DataFileNames_Geometry = np.empty( nFiles, object )


    # Arrays to hold data

    Data_GF = np.empty( nFiles, object )

    GF = np.empty( nFiles, object )

    for i in range( nFiles ):



        DataFileNames_Geometry[i] = GF_root + '.h5'
        Data_GF[i]                = h5.File( DataFileNames_Geometry[i], 'r' )

    # Get the spatial grid
    X  = Data_GF[0][ 'Spatial Grid' ]
    X1 = np.array( X[ 'X1' ] )
    X2 = np.array( X[ 'X2' ] )
    X3 = np.array( X[ 'X3' ] )


    # Geometry fields
    GF_000=np.empty( nFiles, object )
    GF_001=np.empty( nFiles, object )
    GF_002=np.empty( nFiles, object )
    GF_003=np.empty( nFiles, object )
    GF_010=np.empty( nFiles, object )
    GF_011=np.empty( nFiles, object )
    GF_012=np.empty( nFiles, object )
    GF_013=np.empty( nFiles, object )
    GF_020=np.empty( nFiles, object )
    GF_021=np.empty( nFiles, object )
    GF_022=np.empty( nFiles, object )
    GF_023=np.empty( nFiles, object )
    GF_030=np.empty( nFiles, object )
    GF_031=np.empty( nFiles, object )
    GF_032=np.empty( nFiles, object )
    GF_033=np.empty( nFiles, object )
    GF_100=np.empty( nFiles, object )
    GF_101=np.empty( nFiles, object )
    GF_102=np.empty( nFiles, object )
    GF_103=np.empty( nFiles, object )
    GF_110=np.empty( nFiles, object )
    GF_111=np.empty( nFiles, object )
    GF_112=np.empty( nFiles, object )
    GF_113=np.empty( nFiles, object )
    GF_120=np.empty( nFiles, object )
    GF_121=np.empty( nFiles, object )
    GF_122=np.empty( nFiles, object )
    GF_123=np.empty( nFiles, object )
    GF_130=np.empty( nFiles, object )
    GF_131=np.empty( nFiles, object )
    GF_132=np.empty( nFiles, object )
    GF_133=np.empty( nFiles, object )
    GF_200=np.empty( nFiles, object )
    GF_201=np.empty( nFiles, object )
    GF_202=np.empty( nFiles, object )
    GF_203=np.empty( nFiles, object )
    GF_210=np.empty( nFiles, object )
    GF_211=np.empty( nFiles, object )
    GF_212=np.empty( nFiles, object )
    GF_213=np.empty( nFiles, object )
    GF_220=np.empty( nFiles, object )
    GF_221=np.empty( nFiles, object )
    GF_222=np.empty( nFiles, object )
    GF_223=np.empty( nFiles, object )
    GF_230=np.empty( nFiles, object )
    GF_231=np.empty( nFiles, object )
    GF_232=np.empty( nFiles, object )
    GF_233=np.empty( nFiles, object )
    GF_300=np.empty( nFiles, object )
    GF_301=np.empty( nFiles, object )
    GF_302=np.empty( nFiles, object )
    GF_303=np.empty( nFiles, object )
    GF_310=np.empty( nFiles, object )
    GF_311=np.empty( nFiles, object )
    GF_312=np.empty( nFiles, object )
    GF_313=np.empty( nFiles, object )
    GF_320=np.empty( nFiles, object )
    GF_321=np.empty( nFiles, object )
    GF_322=np.empty( nFiles, object )
    GF_323=np.empty( nFiles, object )
    GF_330=np.empty( nFiles, object )
    GF_331=np.empty( nFiles, object )
    GF_332=np.empty( nFiles, object )
    GF_333=np.empty( nFiles, object )

    for i in range( nFiles ):

        # First level groups

        GF[i] = Data_GF[i][ 'Geometry Fields' ]


        GF_000=GF[i]['000'][:][:][:]
        GF_001=GF[i]['001'][:][:][:]
        GF_002=GF[i]['002'][:][:][:]
        GF_003=GF[i]['003'][:][:][:]
        GF_010=GF[i]['010'][:][:][:]
        GF_011=GF[i]['011'][:][:][:]
        GF_012=GF[i]['012'][:][:][:]
        GF_013=GF[i]['013'][:][:][:]
        GF_020=GF[i]['020'][:][:][:]
        GF_021=GF[i]['021'][:][:][:]
        GF_022=GF[i]['022'][:][:][:]
        GF_023=GF[i]['023'][:][:][:]
        GF_030=GF[i]['030'][:][:][:]
        GF_031=GF[i]['031'][:][:][:]
        GF_032=GF[i]['032'][:][:][:]
        GF_033=GF[i]['033'][:][:][:]
        GF_100=GF[i]['100'][:][:][:]
        GF_101=GF[i]['101'][:][:][:]
        GF_102=GF[i]['102'][:][:][:]
        GF_103=GF[i]['103'][:][:][:]
        GF_110=GF[i]['110'][:][:][:]
        GF_111=GF[i]['111'][:][:][:]
        GF_112=GF[i]['112'][:][:][:]
        GF_113=GF[i]['113'][:][:][:]
        GF_120=GF[i]['120'][:][:][:]
        GF_121=GF[i]['121'][:][:][:]
        GF_122=GF[i]['122'][:][:][:]
        GF_123=GF[i]['123'][:][:][:]
        GF_130=GF[i]['130'][:][:][:]
        GF_131=GF[i]['131'][:][:][:]
        GF_132=GF[i]['132'][:][:][:]
        GF_133=GF[i]['133'][:][:][:]
        GF_200=GF[i]['200'][:][:][:]
        GF_201=GF[i]['201'][:][:][:]
        GF_202=GF[i]['202'][:][:][:]
        GF_203=GF[i]['203'][:][:][:]
        GF_210=GF[i]['210'][:][:][:]
        GF_211=GF[i]['211'][:][:][:]
        GF_212=GF[i]['212'][:][:][:]
        GF_213=GF[i]['213'][:][:][:]
        GF_220=GF[i]['220'][:][:][:]
        GF_221=GF[i]['221'][:][:][:]
        GF_222=GF[i]['222'][:][:][:]
        GF_223=GF[i]['223'][:][:][:]
        GF_230=GF[i]['230'][:][:][:]
        GF_231=GF[i]['231'][:][:][:]
        GF_232=GF[i]['232'][:][:][:]
        GF_233=GF[i]['233'][:][:][:]
        GF_300=GF[i]['300'][:][:][:]
        GF_301=GF[i]['301'][:][:][:]
        GF_302=GF[i]['302'][:][:][:]
        GF_303=GF[i]['303'][:][:][:]
        GF_310=GF[i]['310'][:][:][:]
        GF_311=GF[i]['311'][:][:][:]
        GF_312=GF[i]['312'][:][:][:]
        GF_313=GF[i]['313'][:][:][:]
        GF_320=GF[i]['320'][:][:][:]
        GF_321=GF[i]['321'][:][:][:]
        GF_322=GF[i]['322'][:][:][:]
        GF_323=GF[i]['323'][:][:][:]
        GF_330=GF[i]['330'][:][:][:]
        GF_331=GF[i]['331'][:][:][:]
        GF_332=GF[i]['332'][:][:][:]
        GF_333=GF[i]['333'][:][:][:]

    names =                             \
    {                                   \
      'X1'                 : X1       , \
      'X2'                 : X2       , \
      'X3'                 : X3       , \
      'GF_000'             : GF_000   , \
      'GF_001'             :GF_001       ,  \
      'GF_002'             :GF_002       ,  \
      'GF_003'             :GF_003       ,  \
      'GF_010'             :GF_010       ,  \
      'GF_011'             :GF_011       ,  \
      'GF_012'             :GF_012       ,  \
      'GF_013'             :GF_013       ,  \
      'GF_020'             :GF_020       ,  \
      'GF_021'             :GF_021       ,  \
      'GF_022'             :GF_022       ,  \
      'GF_023'             :GF_023       ,  \
      'GF_030'             :GF_030       ,  \
      'GF_031'             :GF_031       ,  \
      'GF_032'             :GF_032       , \
      'GF_033'             :GF_033       ,  \
      'GF_100'             :GF_100       ,  \
      'GF_101'             :GF_101       , \
      'GF_102'             :GF_102       ,  \
      'GF_103'             :GF_103       ,  \
      'GF_110'             :GF_110       , \
      'GF_111'             :GF_111       ,  \
      'GF_112'             :GF_112       ,  \
      'GF_113'             :GF_113       ,  \
      'GF_120'             :GF_120       , \
      'GF_121'             :GF_121       ,  \
      'GF_122'             :GF_122       ,  \
      'GF_123'             :GF_123       ,  \
      'GF_130'             :GF_130       ,  \
      'GF_131'             :GF_131       , \
      'GF_132'             :GF_132       ,  \
      'GF_133'             :GF_133       ,  \
      'GF_200'             :GF_200       ,  \
      'GF_201'             :GF_201       ,  \
      'GF_202'             :GF_202       ,  \
      'GF_203'             :GF_203       , \
      'GF_210'             :GF_210       ,  \
      'GF_211'             :GF_211       ,  \
      'GF_212'             :GF_212       ,  \
      'GF_213'             :GF_213       ,  \
      'GF_220'             :GF_220       ,  \
      'GF_221'             :GF_221       , \
      'GF_222'             :GF_222       ,  \
      'GF_223'             :GF_223       ,  \
      'GF_230'             :GF_230       ,  \
      'GF_231'             :GF_231       ,  \
      'GF_232'             :GF_232       , \
      'GF_233'             :GF_233       ,  \
      'GF_300'             :GF_300       ,  \
      'GF_301'             :GF_301       ,  \
      'GF_302'             :GF_302       , \
      'GF_303'             :GF_303       ,  \
      'GF_310'             :GF_310       ,  \
      'GF_311'             :GF_311       ,  \
      'GF_312'             :GF_312       ,  \
      'GF_313'             :GF_313       ,  \
      'GF_320'             :GF_320       ,  \
      'GF_321'             :GF_321       ,  \
      'GF_322'             :GF_322       ,  \
      'GF_323'             :GF_323       , \
      'GF_330'             :GF_330       ,  \
      'GF_331'             :GF_331       ,  \
      'GF_332'             :GF_332       ,  \
      'GF_333'             :GF_333
    }


    return names
