#!/bin/bash

set -o pipefail

# Total number of junctions
TOTJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | cut -f1-6 | sort | uniq | wc -l)

# number of annotated junctions (column 6=1)
ANNJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($6==1)' | cut -f1-6 | sort | uniq |  wc -l)

# number of junctions in scaffold sequence
SCAFFOLDJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 ~ var)' | cut -f1-6 | sort | uniq |  wc -l)

# number of non-canonical junctions (column 5 > 0)
NCJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($5==0)' | cut -f1-6 | sort | uniq |  wc -l) #total number of non-canonical junctions = 615

# total not annotated, not in scaffold, canonical
TOTREM=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | cut -f1-6 | sort | uniq | wc -l)

# obtain filtered list
cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab

FINALNUM=$(less ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab | wc -l)


echo "The total number of junctions was ${TOTJUNC}, ${ANNJUNC} were already annotated"
echo "There were ${SCAFFOLDJUNC} junctions on scaffold sequence"
echo "There were ${NCJUNC} non-canonical junctions"
echo "Removing annotated, non-canonical, and scaffold junctions left ${TOTREM} junctions"
echo "Out of ${TOTREM}, ${FINALNUM} junctions were supported by at least ${UNIQUE_NUM} uniquely mapped reads"
echo "The final filtered list of ${FINALNUM} junctions is: ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab and can be used for 2nd-pass read mapping"