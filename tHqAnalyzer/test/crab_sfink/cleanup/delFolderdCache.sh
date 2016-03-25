#!/bin/bash

DIR=$1

echo `ls $DIR`

for i in `ls $DIR`
do 
  lcg-del -v -b -l -T  srmv2 "srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=$DIR/$i"
  #echo "srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=$DIR/$i"  
done

srmrmdir "srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=$DIR"
