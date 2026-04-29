#!/usr/bin/env bash
set -euo pipefail

SAG_SPECIES_LIST="Streptococcus anginosus|Streptococcus intermedius|Streptococcus constellatus"

# Dominant SAG species from Bracken report
SAG_RECORD=$(awk -F"\t" -v targets="$SAG_SPECIES_LIST" '
    NR > 1 && $3 == "S" && $1 ~ targets {
        print $0
    }
' "$BRACKEN_REPORT" | sort -t $'\t' -nr -k7,7 | head -n 1)

SAG_FRACTION=$(awk -F"\t" '{ print $7 }' <<< "$SAG_RECORD")
SAG_PERCENTAGE=$(awk -v f="${SAG_FRACTION:-0}" 'BEGIN { printf "%.2f", f * 100 }')

SAG_SPECIES=$(awk -F"\t" '{
    name=$1
    gsub(/ /, "_", name)
    print name
}' <<< "$SAG_RECORD")

# Top non-SAG species/genus from Bracken report
TOP_NON_SAG_RECORD=$(awk -F"\t" -v targets="$SAG_SPECIES_LIST" '
    NR > 1 && $3 == "S" && $1 !~ targets {
        print $0
    }
' "$BRACKEN_REPORT" | sort -t $'\t' -nr -k7,7 | head -n 1)

TOP_NON_SAG=$(awk -F"\t" '{ print $1 }' <<< "$TOP_NON_SAG_RECORD")
TOP_NON_SAG_FRACTION=$(awk -F"\t" '{ print $7 }' <<< "$TOP_NON_SAG_RECORD")
TOP_NON_SAG_PERCENTAGE=$(awk -v f="${TOP_NON_SAG_FRACTION:-0}" 'BEGIN { printf "%.2f", f * 100 }')

if [ -z "${SAG_SPECIES:-}" ]; then
    SAG_SPECIES="Unknown"
fi

if [[ "$(echo "$SAG_PERCENTAGE >= $QC_SAG_SPECIES_PERCENTAGE" | bc -l)" == 1 ]] && \
   [[ "$(echo "$TOP_NON_SAG_PERCENTAGE <= $QC_TOP_NON_SAG_PERCENTAGE" | bc -l)" == 1 ]]; then
    TAXONOMY_QC="PASS"
else
    TAXONOMY_QC="FAIL"
fi

export TAXONOMY_QC
export SAG_SPECIES

echo "\"Taxonomy_QC\",\"SAG_Species\",\"SAG_Species_%\",\"Top_Non_SAG_Species\",\"Top_Non_SAG_Species_%\"" > "$TAXONOMY_QC_REPORT"
echo "\"$TAXONOMY_QC\",\"$SAG_SPECIES\",\"$SAG_PERCENTAGE\",\"${TOP_NON_SAG:-}\",\"$TOP_NON_SAG_PERCENTAGE\"" >> "$TAXONOMY_QC_REPORT"