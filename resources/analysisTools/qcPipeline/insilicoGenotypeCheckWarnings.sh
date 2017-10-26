#!/bin/bash

#PBS -l walltime=0:2:00
#PBS -l nodes=1
#PBS -m ae


[ -f ${WARNING_FILE} ] && exit -100

sleep 60