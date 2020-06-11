#!/bin/bash
wdir=/work/opa/witoil-dev/medslik-sanifs/meglob-fabio/EXE
source ${wdir}/set_env.sh
JJ=$1

# Invoke medslik
MSHome=/work/opa/witoil-dev/medslik-sanifs/meglob-fabio/EXE
MSCommand=./RUN.sh
# bsub -e ${wdir}/log/%J.err -o ${wdir}/log/%J.out -q s_short -J ${JJ} $MSHome/$MSCommand
bsub -Is -q s_short -J ${JJ} $MSHome/$MSCommand
