#!/bin/bash
#SBATCH -J ToyMCGeneration
#SBATCH -D /afs/hephy.at/work/t/tmadlener/NewPolMethodTests/Framework/polFit/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/NewPolMethodTests/ToyMCTests/logfiles/runPolGen_%A_%a.out

function cleanup_and_exit() {
  # check if the last action has been succesful (exited with 0), if so remove the exetuable
  # (second argument)
  # else leave it and exit with the passed exitcode (first argument)
  if [ ${1} -eq 0 ]; then
    echo ${2} "exited with 0. Cleaning up executable"
    rm ${2}
  else
    echo ${2} "exited with "${1}". Leaving executable in place"
    echo "-------------------- end: " $(date) " --------------------"
    exit ${1}
  fi
}

echo "-------------------- start: " $(date) " --------------------"
## input arguments
outFileBase=${1}
lthsig=${2}
lthbkg=${3}
lphsig=${4}
lphbkg=${5}
ltpsig=${6}
ltpbkg=${7}
nEvents=${8}
workDir=${9}

exebase=${WORK}/NewPolMethodTests/Framework/polFit/runPolGen

## do the "sandboxing" by giving each executable a unique identifier
iGen=$SLURM_ARRAY_TASK_ID
uniqueJobId="_gen_"${iGen}

mkdir -p ${workDir}
cd ${workDir}

exe="runPolGen"${uniqueJobId}
cp ${exebase} ${exe}

outfile=${outFileBase}${uniqueJobId}.root

./${exe} --outfile ${outfile} \
  --lthsig ${lthsig} --lthbkg ${lthbkg} \
  --lphsig ${lphsig} --lphbkg ${lphbkg} \
  --ltpsig ${ltpsig} --ltpbkg ${ltpbkg} \
  --nEvents ${nEvents}

cleanup_and_exit $? ${exe}

## run polSub on the just created file
polSub=${WORK}/NewPolMethodTests/Framework/polFit/runPolSub

${polSub} --file ${outfile}
echo "-------------------- end: " $(date) " --------------------"
