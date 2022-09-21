f2py3 -c Constants_Module.f90 Utils_Module.f GW_Chimera3D_Module.f90 Meridional_Int_Module.f90 -m Meridional | tee f2py_out.txt
python3 thornadoStrain_Analytical.py | tee strain_out.txt
