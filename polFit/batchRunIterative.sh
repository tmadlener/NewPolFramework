#!/bin/bash
#SBATCH -J RunIterativeFit
#SBATCH -D /afs/hephy.at/work/t/tmadlener/NewPolMethodTests/Framework/polFit/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/NewPolMethodTests/ToyMCTests/logfiles/runIterativeFit_%A.out

echo '-------------------- start: '$(date)' --------------------'

dataFileName=${1}
refFileName=${2}
treeName=${3} # same name for both ref and data
outFileName=${4} # make sure that the folder this file points to already exists! (it will not be checked)

cd $(dirname ${outFileName})

exe=${WORK}/NewPolMethodTests/Framework/polFit/run_iterative_fit.py

python ${exe} --tree ${treeName} --nMaxIterations 10 --stopSignificance 0.5 \
       --lthRefStart 0.5 ${dataFileName} ${refFileName} ${outFileName}

# capture exitcode
exitcode=$?

# cleanup
rm pTweight_fitted_iter*.pdf
rm Polarization_fitted_iter*.pdf
rm Lambdas_starting.pdf

echo '-------------------- end: '$(date)' --------------------'
exit ${exitcode}
