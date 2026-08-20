declare -A params

Perturb=('.FALSE.')
PerturbTxt=('NP')
BField=('0.00e+00')
BFieldTxt=('GRHD')
Order=('2')
OrderTxt=('2OS_2OT_SLOn_BTVD1.75e+00')
BC=('44')
BCTxt=('44BCs-FPD-BSP-CH')
IC=('zerorot')
ICTxt=('ZRot')

i=0
for indx_p in "${!Perturb[@]}"
do
	for indx_f in "${!BField[@]}"
	do
		for indx_c in "${!IC[@]}"
		do
			for indx_o in "${!Order[@]}"
			do
				for indx_b in "${!BC[@]}"
				do
					params[$i,0]="NAT_SD_1D_${BFieldTxt[$indx_f]}_${ICTxt[$indx_c]}_${PerturbTxt[$indx_p]}_${BCTxt[$indx_b]}_${OrderTxt[$indx_o]}_Short_DEB"
					params[$i,1]="${IC[$indx_c]}"
					params[$i,2]="${Order[$indx_o]}"
					params[$i,3]="${BC[$indx_b]}"
					params[$i,4]="${Perturb[$indx_p]}"
					params[$i,5]="${BField[$indx_f]}"
					((i++))
				done
			done
		done
	done
done
