nextflow.enable.dsl=2

// Define a parameter for a dummy samplesheet file (update the path as needed)
params.dummy_samplesheet = '/path/to/dummy_samplesheet.csv'

// You can comment out or remove upstream processes that are not needed for testing CURIOTREKKER
/*
process CELLRANGER { ... }
process CELLBENDER { ... }
process CELLBENDERPOSTPROCESSING { ... }
process CURIOTREKKERSAMPLESETUP { ... }
*/

process CURIOTREKKER {
    publishDir "${id}_results/CurioTrekker"
    label "curiotrekker_level"

    input:
        tuple val(id), path(samplesheet)

    output:
        path "output"

    script:
    """
    module load StdEnv/2023
    module load apptainer/1.3.4

    # Get the absolute path to the samplesheet
    SAMPLESHEET_PATH=\$(readlink -f ${samplesheet})

    # Use the absolute path
    cd "/home/hudsonhu/scratch/FreshCurioTrekker/curiotrekker-v1.1.0/"
    bash nuclei_locater_toplevel.sh \$SAMPLESHEET_PATH
    """
}

workflow {
    // Create a dummy input channel for the CURIOTREKKER process.
    // The tuple contains a dummy id ("dummy_sample") and the dummy samplesheet file.
    Channel
        .of( tuple("dummy_sample", file(params.dummy_samplesheet)) )
        .set { curiotrekker_input }

    CURIOTREKKER(curiotrekker_input)
}
