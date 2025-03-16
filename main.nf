nextflow.enable.dsl=2

params.samplesheet = '/path/to/samplesheet.csv'

process CELLRANGER {
    publishDir "${id}_results/CellRanger"
    label "cellranger_level"
    
    input:
        tuple val(id), val(name), val(transcriptome), val(cDNA_fastq_dir)

    output:
        tuple val(id), path("${id}/outs")

    script:
    """
    ${projectDir}/PipelineFolder/CellRanger/cellranger-9.0.1/cellranger count \
        --id=${id} \
        --transcriptome=${transcriptome} \
        --fastqs=${cDNA_fastq_dir} \
        --create-bam false \
        --tenx-cloud-token-path=${projectDir}/PipelineFolder/CellRanger/cellranger-9.0.1/config/txg/credentials
    """
}

process CELLBENDER {
    publishDir "${id}_results/CellBender"
    label "cellbender_level"
    
    // Takes the output directory from CELLRANGER as input.
    input:
        tuple val(id), path(cellranger_out)

    output:
        tuple val(id), path("cellbender_output")

    script:
    """
    mkdir -p cellbender_output
    apptainer exec --nv ${projectDir}/PipelineFolder/CellBender/cellbender_latest.sif cellbender remove-background \
        --cuda \
        --input ${cellranger_out}/raw_feature_bc_matrix.h5 \
        --model "ambient" \
        --output cellbender_output/output_file.h5 \
        --epochs 150
    """
}

process CELLBENDERPOSTPROCESSING {
    publishDir "${id}_results/CellBender/"

    input:
        tuple val(id), path(cellbender_out)
    
    output:
        tuple val(id), path("sc_out")
    
    script:
    """
    module load StdEnv/2023 apptainer/1.3.4
    export CELLBENDER_OUT=${cellbender_out}
    apptainer exec ${projectDir}/PipelineFolder/Containers/slide_tag_env_latest.sif python ${projectDir}/process_data.py
    """
}

process CURIOTREKKERSAMPLESETUP {
    label "process_default"
    // Expects a tuple with: name, sp_fastq1, sp_fastq2, and the cellbender-processed output
    input:
        tuple val(name), file(sp_fastq1), file(sp_fastq2), file(bead_barcodes), val(chip)
        tuple val(id), path(cellbender_processed)

    output:
        tuple val(id), path("curiotrekker_samplesheet.csv")

    script:
    def current_date = new java.text.SimpleDateFormat("yyyy-MM-dd").format(new Date())
    """
    echo sample,sample_sc,experiment_date,barcode_file,fastq_1,fastq_2,sc_outdir,sc_platform,profile,subsample,cores > curiotrekker_samplesheet.csv
    echo ${name},${name}_sc,${current_date},${bead_barcodes},${sp_fastq1},${sp_fastq2},${cellbender_processed},${chip},singularity,no,16 >> curiotrekker_samplesheet.csv
    """
}

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
    cd ${projectDir}/PipelineFolder/CurioTrekker/curiotrekker-v1.1.0/
    bash nuclei_locater_toplevel.sh \$SAMPLESHEET_PATH
    """
}

workflow {
    /*
      Read the CSV, skip the header, and create a unicast channel of sample tuples.
      Each tuple contains:
      (id, name, fastq_cDNA_dir, fastq_R1_SP, fastq_R2_SP, transcriptome, barcode_file)
    */
    def sampleChannel = Channel.fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
        // Remove any quotes from the id value
        def cleanId = row.id.replaceAll(/['"]/, '')
        tuple(cleanId, row.name, row.fastq_cDNA_dir, row.fastq_R1_SP, row.fastq_R2_SP, row.transcriptome, row.barcode_file, row.chip)
    }
    
    cellranger_in = sampleChannel.map {
        id, name, fastq_cDNA_dir, _fastq_R1_SP, _fastq_R2_SP, transcriptome, _barcode_file, _chip ->
        tuple(id, name, transcriptome, fastq_cDNA_dir)
    }
    cellranger_out = CELLRANGER(cellranger_in)
    cellbender_out = CELLBENDER(cellranger_out)
    cellbender_processed = CELLBENDERPOSTPROCESSING(cellbender_out)
    
    trekkersample_in = sampleChannel.map {
        _id, name, _fastq_cDNA_dir, fastq_R1_SP, fastq_R2_SP, _transcriptome, barcode_file, chip ->
        tuple(name, fastq_R1_SP, fastq_R2_SP, barcode_file, chip)
    }
    trekkersamplesheet = CURIOTREKKERSAMPLESETUP(trekkersample_in, cellbender_processed)
    CURIOTREKKER(trekkersamplesheet)
}