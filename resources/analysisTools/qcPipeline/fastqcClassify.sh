#!/usr/bin/env bash
inputFile="${1:?No input file given}"
gawk 'BEGIN{FS="\t";OFS="\t";sum=0;n=0;down=1==2}/Median/{flag=1;next}/^>>END/{flag=0}flag{sum+=$4;n++;down=down || $3<=20}END{status="PASS"; avg=sum/n;if (avg<=28){if (avg<=20) {status="FAIL"} else {status="WARN"}}; if (down && status == "PASS"){status="WARN"}; print status}' "$inputFile"


