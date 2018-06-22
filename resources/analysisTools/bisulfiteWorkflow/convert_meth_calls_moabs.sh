#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

# Reads methylation calls from methylCtools and writes to stdout
input_filename=$1
context=$2 # Either 'CG' or 'CH'

# Convert methylCTools meth calls to MOABS compatible meth calls
echo -e "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC\tsnpScore"

if [[ ${context} == 'CG' ]]; then
	# MOABS conversion CG case
	awk 'BEGIN{FS="\t"; OFS="\t";\
	        COUNTER=0;\
	        CHROM="";\
	        START=0;\
	        E=0;\
	        RATIO=0;\
	        TOTALCB=0;\
	        METHCB=0;\
	        STRAND="";\
	        NEXT="G";\
	        TOTALCP=0;\
	        TOTALCM=0;\
	        METHCP=0;\
	        METHCM=0;\
	        SNPSCORE=0;\
	        TOTALC=0;\
	        METHC=0}\
	        {\
	                if($4 == "CG"){\
	                        SNPSCORE = $5;\
	                        COUNTER += 1;\
	                        if($3 == "+"){\
	                                TOTALCB += $6+$7;\
	                                METHCB += $6;\
	                                TOTALCP += $6+$7;\
	                                METHCP += $6;\
	                                START = $2;\
	                        }\
	                        else{\
	                                TOTALCB += $6+$7;\
	                                METHCB += $6;\
	                                TOTALCM += $6+$7;\
	                                METHCM += $6;\
	                                E = $2;\
	                        }\
	                        if(COUNTER == 2){\
	                                if(TOTALCM == 0 && TOTALCP > 0){\
	                                        STRAND="+";\
	                                        RATIO=METHCB/TOTALCB;\
	                                }\
	                                else if(TOTALCP == 0 && TOTALCM > 0){\
	                                        STRAND="-";\
	                                        RATIO=METHCB/TOTALCB;\
	                                }\
	                                else if(TOTALCP == 0 && TOTALCM == 0){\
	                                        STRAND="B";\
	                                        RATIO=0;\
	                                }
	                                else{\
	                                        STRAND="B";\
	                                        RATIO=METHCB/TOTALCB;\
	                                }\
	                                print $1, START, E+1, RATIO, TOTALCB, METHCB, STRAND, NEXT, "+", TOTALCP, METHCP, "-", TOTALCM, METHCM, SNPSCORE;\
	                                COUNTER=0;\
	                                CHROM="";\
	                                START=0;\
	                                E=0;\
	                                RATIO=0;\
	                                TOTALCB=0;\
	                                METHCB=0;\
	                                STRAND="";\
	                                NEXT="G";\
	                                TOTALCP=0;\
	                                TOTALCM=0;\
	                                METHCP=0;\
	                                METHCM=0;\
	                                SNPSCORE=0;\
	                                TOTALC=0;\
	                                METHC=0}\
	                        }\
	        }' ${input_filename}
elif [[ ${context} == 'CH' ]]; then
	# MOABS conversion CH case
	awk 'BEGIN{FS="\t"; OFS="\t"}{\
	                if($4 == "CH"){\
	                        NEXT = "H";\
	                        SNPSCORE = $5;\
	                        TOTALC = $6+$7;\
	                        METHC = $6;\
	                        START = $2-1;\
	                        E = $2;\
	                        STRAND = $3;\
	                        RATIO = 0;\
	                        if(TOTALC > 0){\
	                                RATIO=METHC/TOTALC;\
	                        }\
	                        TOTALCP=0;\
	                        METHCP=0;\
	                        TOTALCM=0;\
	                        METHCM=0;\
	                        if(STRAND == "+"){\
	                                TOTALCP = TOTALC;\
	                                METHCP = METHC;\
	                        }\
	                        else{\
	                                TOTALCM = TOTALC;\
	                                METHCM = METHC;\
	                        }\
	                        print $1, START, E, RATIO, TOTALC, METHC, STRAND, NEXT, "+", TOTALCP, METHCP, "-", TOTALCM, METHCM, SNPSCORE;\
	                }\
	        }' ${input_filename}
fi;
