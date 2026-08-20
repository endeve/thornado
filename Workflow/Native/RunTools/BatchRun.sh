source ./BatchRun_Params.sh

EXE_DIR="$PWD"
REL_OUTPUT_DIR=../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output

for ((j=0; j<$i; j++))
do
	PROBLEM_DIR=${params["$j,0"]}
	rm -r ../Output/${PROBLEM_DIR}
	mkdir ../Output/${PROBLEM_DIR}
	sed -i "153c\  ProgramName = 'ShearingDisk'" ../ApplicationDriver.F90
	sed -i "706c\      bcX = [ ${params["$j,3"]}, 1, 1 ]" ../ApplicationDriver.F90
	sed -i "760c\  nNodes = ${params["$j,2"]}" ../ApplicationDriver.F90
	sed -i "766c\  nStagesSSPRK = ${params["$j,2"]}" ../ApplicationDriver.F90
	sed -i "702c\      SDICFileName = './GR_LR_${params["$j,1"]}.h5'" ../ApplicationDriver.F90
	sed -i "686c\      ApplyRandomPerturbations = ${params["$j,4"]}" ../ApplicationDriver.F90
	sed -i "689c\      SDInitialField = ${params["$j,5"]} * Gauss" ../ApplicationDriver.F90
	make -j 4
	./ApplicationDriver_${THORNADO_MACHINE} | tee out.txt
	mv ../Output/ShearingDisk_*.h5 ../Output/${PROBLEM_DIR}
	mv ../Output/ShearingDisk_Tally_*.dat ../Output/${PROBLEM_DIR}
	mv out.txt ../Output/${PROBLEM_DIR}
	cp ../ApplicationDriver.F90 ../Output/${PROBLEM_DIR}
	cd $THORNADO_DIR/Workflow/Native/Python
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Primitive' 'Comoving Baryon Density'
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Primitive' 'Three-Velocity (3)'
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Primitive' 'Three-Velocity (1)'
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Conserved' 'Conserved Magnetic Field (1)'
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Conserved' 'Conserved Magnetic Field (2)'
	python3 NativeMovie.py ${REL_OUTPUT_DIR}/${PROBLEM_DIR} 'Conserved' 'Conserved Magnetic Field (3)'
	python3 Energetics.py  ${REL_OUTPUT_DIR}/${PROBLEM_DIR}
	mv Images_${PROBLEM_DIR}_* ${REL_OUTPUT_DIR}/${PROBLEM_DIR}
	mv movie_${PROBLEM_DIR}_*.mp4 ${REL_OUTPUT_DIR}/${PROBLEM_DIR}
    mv ${PROBLEM_DIR}_Energetics.png ${REL_OUTPUT_DIR}/${PROBLEM_DIR}
	cd $EXE_DIR
done
