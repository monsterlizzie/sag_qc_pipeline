# Determine overall QC result based on File Validity, Read QC, Assembly QC, Taxonomy QC, and Mapping QC

assign_overall_qc() {
    if [[ "$FILE_VALIDITY" == "null" ]]; then
        OVERALL_QC="FILE VALIDATION FAILURE"
        return
    fi
    
    if [[ "$FILE_VALIDITY" != "PASS" ]]; then
        OVERALL_QC="$FILE_VALIDITY"
        return
    fi
    
    if [[ "$READ_QC" == "null" ]]; then
        OVERALL_QC="PREPROCESS MODULE FAILURE"
        return
    fi 
    
    if [[ "$READ_QC" == "FAIL" ]]; then
        OVERALL_QC="FAIL"
        return
    fi 

    if [[ "$ASSEMBLY_QC" == "null" ]]; then
        OVERALL_QC="ASSEMBLY MODULE FAILURE"
        return
    fi

    if [[ "$TAXONOMY_QC" == "null" ]]; then
        OVERALL_QC="TAXONOMY MODULE FAILURE"
        return
    fi

    if [[ "$TAXONOMY_QC" == "FAIL" ]]; then
        OVERALL_QC="FAIL"
        return
    fi

    # Mapping is only expected if taxonomy passed
    if [[ "$MAPPING_QC" == "null" ]]; then
        OVERALL_QC="MAPPING MODULE FAILURE"
        return
    fi

    if [[ "$READ_QC" == "PASS" ]] && \
       [[ "$ASSEMBLY_QC" == "PASS" ]] && \
       [[ "$TAXONOMY_QC" == "PASS" ]] && \
       [[ "$MAPPING_QC" == "PASS" ]]; then
        OVERALL_QC="PASS"
    else
        OVERALL_QC="FAIL"
    fi
}

assign_overall_qc

echo \"Overall_QC\" > "$OVERALL_QC_REPORT"
echo \""$OVERALL_QC"\" >> "$OVERALL_QC_REPORT"