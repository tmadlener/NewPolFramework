#!/bin/sh

########## INPUTS ##########

JobID=Upsilon2011_full

RequestTrigger=1 ###set1

#following flags decide if the step is executed (1) or not (0):
execute_prepareEvents=0
execute_runMassFit=1
execute_runBoostAngles=1
execute_runReshuffleNch=1


#Real pp 2011 upsilon data tree:
inputTree1=/hadoop/store/user/cferraio/Polarization_2011pp/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Upsi.root

############################

rejectCowboys=false

JobDir=$HOME/work/PhD/tests_Ups_v_Nch_new_pol/DataFiles/${JobID}

mkdir -p ${JobDir}/Figures
mkdir -p ${JobDir}/tmpFiles
mkdir -p ${JobDir}/PDF


cp Makefile ${JobDir}/Makefile

cp prepareEvents.C ${JobDir}/prepareEvents.C
cp prepareEvents.h ${JobDir}/prepareEvents.h
cp runPrepareEvents.cc ${JobDir}/runPrepareEvents.cc
cp rootIncludes.inc ${JobDir}/rootIncludes.inc
cp runMassFit.cc ${JobDir}/runMassFit.cc
cp upsilon_2StepFit.C ${JobDir}/upsilon_2StepFit.C
cp CBFunction.C ${JobDir}/CBFunction.C
cp runBoostAngles.cc ${JobDir}/runBoostAngles.cc
cp BoostAngles.C ${JobDir}/BoostAngles.C
cp runReshuffleNch.cc ${JobDir}/runReshuffleNch.cc
cp ReshuffleNch.C ${JobDir}/ReshuffleNch.C

## already processed / prepared data
selEvents_data=$HOME/cernbox/Chic/NewFitInputFiles/selEvents_data_Ups.root
cp ${selEvents_data} ${JobDir}/tmpFiles

cd ${JobDir}

make

if [ ${execute_prepareEvents} -eq 1 ]
then
  ./runPrepareEvents rejectCowboys=${rejectCowboys} ${inputTree1}=inputTree1 RequestTrigger=${RequestTrigger}
fi

if [ ${execute_runMassFit} -eq 1 ]
then
  ./runMassFit
fi

if [ ${execute_runBoostAngles} -eq 1 ]
then
  ./runBoostAngles
fi

if [ ${execute_runReshuffleNch} -eq 1 ]
then
  ./runReshuffleNch
fi

rm *_C.d
rm *_C.so
rm runMassFit
rm runPrepareEvents
rm runBoostAngles
cd ..
cd ..
