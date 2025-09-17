#!/bin/bash

PREFIX='si'

#should be in the work directory of PHonon calculation
echo `pwd`
mkdir -p save
mkdir -p save/${PREFIX}.phsave

DIR="./tmp/_ph0"
echo $DIR

#copy prefix.phsave
cp ${DIR}/${PREFIX}.phsave/* save/${PREFIX}.phsave/

#copy dyn files
cp ${PREFIX}.dyn* save/

# For NQ=1
NQ=1
cp ${DIR}/${PREFIX}.dvscf save/${PREFIX}.dvscf_q${NQ}

# Other Q's
for ((NQ=2; NQ<=3; NQ++ ))
do
   cp ${DIR}/${PREFIX}.q_${NQ}/${PREFIX}.dvscf save/${PREFIX}.dvscf_q${NQ}
done
