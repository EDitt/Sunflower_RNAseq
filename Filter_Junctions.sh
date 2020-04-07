#!/bin/bash

set -o pipefail

# Total number of junctions
TOTJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | cut -f1-6 | sort | uniq | wc -l)

# number of annotated junctions (column 6=1)
ANNJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($6==1)' | cut -f1-6 | sort | uniq |  wc -l)

# number of non-canonical junctions (column 5 > 0)
NCJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($5==0)' | cut -f1-6 | sort | uniq |  wc -l) #total number of non-canonical junctions = 615

# obtain filtered list
#cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab

echo "The total number of junctions was ${TOTJUNC}, ${ANNJUNC} were already annotated"


if [[ "$SCAFFOLD_STRING" == "NA" ]]; then
	echo "No scaffold sequence specified"
	if [[ "$REMOVE_NC_JUNC" == "yes" ]]; then
		echo "Removing ${NCJUNC} non-canonical junctions"
		cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($5 > 0 && $6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
	else
		echo "Not removing ${NCJUNC} non-canonical junctions"
		cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
	fi
else
	SCAFFOLDJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 ~ var)' | cut -f1-6 | sort | uniq |  wc -l)
	echo "Removing ${SCAFFOLDJUNC} junctions from scaffold sequence"
	if [[ "$REMOVE_NC_JUNC" == "yes" ]]; then
		echo "Removing ${NCJUNC} non-canonical junctions"
		#TOTREM=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | cut -f1-6 | sort | uniq | wc -l) # total not annotated, not in scaffold, canonical
		#echo "Removing annotated, non-canonical, and scaffold junctions left ${TOTREM} junctions"
		cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
		
	else
		echo "Not removing ${NCJUNC} non-canonical junctions"
		cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $6==0)' | awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | cut -f1-6 | sort | uniq  > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
	fi
fi

FINALNUM=$(less ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab | wc -l)
#echo "Out of ${TOTREM}, ${FINALNUM} junctions were supported by at least ${UNIQUE_NUM} uniquely mapped reads"
echo "After removing junctions supported by less than ${UNIQUE_NUM} uniquely mapped reads..."
echo "The final filtered list of ${FINALNUM} junctions is: ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab and can be used for 2nd-pass read mapping"
