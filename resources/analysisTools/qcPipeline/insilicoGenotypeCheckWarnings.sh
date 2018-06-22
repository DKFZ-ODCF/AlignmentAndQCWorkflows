#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

#PBS -l walltime=0:2:00
#PBS -l nodes=1
#PBS -m ae


[ -f ${WARNING_FILE} ] && exit -100

sleep 60