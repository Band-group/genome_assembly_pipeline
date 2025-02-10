include { RAGTAG_CORRECT  } from '../processes/ragtag_correct'
include { RAGTAG_SCAFFOLD } from '../processes/ragtag_scaffold'

workflow SCAFFOLD {
    take:
        ch_oriented_ctgs
    
    main:
        RAGTAG_CORRECT ( ch_oriented_ctgs )
        rt_corrected_fa = RAGTAG_CORRECT.out.rt_corrected_fa
        RAGTAG_SCAFFOLD ( rt_corrected_fa )
    // ragtag.py correct 
    // ragtag.py scaffold

    // longstitch correct and scaffold 

    // ragtag.py merge

}