#!/bin/bash
date

folder=TestingOriginal #whatever you want the output folder to be named
geomtag=dev2022 # dev2022, dev2022m    (ideal, misaligned)
njobs=1000 # each job is set for 500 events
dir=/star/u/gwilks3/fst/star-fwd-integration
#dir=$(echo "`pwd`" | sed 's:/:\\/:g')

echo "$dir"

#mkdir -p /gpfs01/star/pwg/gwilks3/${folder}/out
#mkdir -p /gpfs01/star/pwg/gwilks3/${folder}/log
#mkdir -p `pwd`/schedinfo

#sed "s/njobs/${njobs}/g;s/geomtag/${geomtag}/g;s/folder/${folder}/g;s/dir/${dir}/g" schedule_starsim.xml > schedule_starsim_${folder}_${geomtag}.xml

logdir=/gpfs01/star/pwg/gwilks3/${folder}/log    
outdir=/gpfs01/star/pwg/gwilks3/${folder}/out 

mkdir -p ${logdir}
mkdir -p ${outdir}
mkdir -p schedinfo

star-submit-template -template starsim.xml -entities dir=$dir,geomtag=$geomtag,njobs=$njobs,logdir=$logdir,outdir=$outdir
