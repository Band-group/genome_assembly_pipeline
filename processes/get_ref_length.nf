process GET_REF_LENGTH {
    label "process_single"

    input:
        val(ref) // meta.reference
    
    output:
        val(ref_length), emit: ref_length

    script:
        """
        ref_length=\$(awk 'BEGIN { total=0 } !/^>/ { total += length(\$0) } END { print total }' ${ref})
        """
}