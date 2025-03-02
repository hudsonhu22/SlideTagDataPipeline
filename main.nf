nextflow.enable.dsl=2

params.samplesheet = '/path/to/samplesheet.csv'

process CELLRANGER {
    publishDir "${id}_results/CellRanger"
    label "cellranger_level"
    
    input:
        tuple val(id), val(name), val(transcriptome), val(cDNA_fastq_dir)

    // Use the id without extra quotes to capture the output folder
    output:
        tuple val(id), path("${id}/outs")

    script:
    """
    /home/hudsonhu/scratch/SlideTagNextflow/PipelineFolder/CellRanger/cellranger-9.0.1/cellranger count \\
        --id=${id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=${cDNA_fastq_dir} \\
        --create-bam false \\
        --sample=${name}
    ls -la .
    """
}

process CELLBENDER {
    container 'file:///home/hudsonhu/scratch/SlideTagNextflow/PipelineFolder/CellBender/cellbender_latest.sif'
    containerOptions '--nv'
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
    module load StdEnv/2023 
    module load apptainer/1.3.4
    cellbender remove-background \\
        --cuda \\
        --input ${cellranger_out}/raw_feature_bc_matrix.h5 \\
        --model "ambient" \\
        --output cellbender_output/output_file.h5 \\
        --epochs 150
    """
}

process CELLBENDERPOSTPROCESSING {
    conda '/home/hudsonhu/scratch/SlideTagNextflow/environment.yml'
    publishDir "${id}_results/CellBender/sc_out/"
    
    input:
        tuple val(id), path(cellbender_out)

    output:
        tuple val(id), path("sc_out")

    script:
    """
    import os
    import h5py
    import numpy as np
    import scipy.sparse as sp
    from scipy.io import mmwrite

    h5_file = os.path.join('${cellbender_out}', "output_file_filtered.h5")

    with h5py.File(h5_file, 'r') as f:
        mat_group = f['matrix']
        barcodes = mat_group['barcodes'][:]
        features_group = mat_group['features']
        # I want to make a tsv with id, name, and type
        id = features_group['id'][:]
        name = features_group['name'][:]
        type = features_group['feature_type'][:]
        features = np.column_stack((id, name, type))
        data = mat_group['data'][:]
        indices = mat_group['indices'][:]
        indptr = mat_group['indptr'][:]
        shape = tuple(mat_group['shape'][:])
        expression_matrix = sp.csc_matrix((data, indices, indptr), shape=shape)

    os.makedirs("sc_out", exist_ok=True)
    mmwrite('sc_out/matrix.mtx', expression_matrix)
    np.savetxt('sc_out/barcodes.tsv', barcodes, fmt='%s', delimiter='\t')
    np.savetxt('sc_out/features.tsv', features, fmt='%s', delimiter='\t')
    """
}

// process CURIOTREKKERMODIFY {
//     script:
//     """
//     CURRENT_DIR=\$(pwd)
//     CURRENT_DATE=\$(date +"%Y-%m-%d")

//     # Instead of referencing the local working directory, reference the pipeline's directory:
//     NEW_SCRIPT_DIR="${workflow.projectDir}/PipelineFolder/CurioTrekker/curiotrekker-v1.1.0/"
//     NEW_OUT_DIR="${workflow.projectDir}/results/CurioTrekker/"
//     bash_script="${NEW_SCRIPT_DIR}/nuclei_locater_toplevel.sh"

//     echo "Current Directory: \${CURRENT_DIR}"
//     echo "Current Date: \${CURRENT_DATE}"
//     echo "NEW_SCRIPT_DIR: \${NEW_SCRIPT_DIR}"
//     echo "NEW_OUT_DIR: \${NEW_OUT_DIR}"

//     sed -i "s|SCRIPT_DIR=.*|SCRIPT_DIR=\${NEW_SCRIPT_DIR}|" \${bash_script}
//     sed -i "s|OUT_DIR=.*|OUT_DIR=\${NEW_OUT_DIR}|" \${bash_script}
//     """
// }


process CURIOTREKKERSAMPLESETUP {
    // Expects a tuple with: name, sp_fastq1, sp_fastq2, and the cellbender-processed output
    input:
        tuple val(name), file(sp_fastq1), file(sp_fastq2), file(BEADED_BARCODES)
        tuple val(id), path(cellbender_processed)

    output:
        tuple val(id), path("curiotrekker_samplesheet.csv")

    script:
    """
    CURRENT_DATE=\$(date +"%Y-%m-%d")
    
    echo "sample,sample_sc,experiment_date,barcode_file,fastq_1,fastq_2,sc_outdir,sc_platform,profile,subsample,cores" > sample_sheet_trekker.csv
    echo "${name},${name}_sc,\${CURRENT_DATE},\${BEADED_BARCODES},${sp_fastq1},${sp_fastq2},${cellbender_processed}/output_file_filtered.h5,TrekkerU_C,singularity,no,16" >> sample_sheet_trekker.csv
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
    cd "/home/hudsonhu/scratch/SlideTagNextflow/PipelineFolder/CurioTrekker/curiotrekker-v1.1.0/"
    bash nuclei_locater_toplevel.sh ${samplesheet}
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
        tuple(cleanId, row.name, row.fastq_cDNA_dir, row.fastq_R1_SP, row.fastq_R2_SP, row.transcriptome, row.barcode_file, row.processing)
    }

    
    cellranger_in = sampleChannel.map {
        id, name, fastq_cDNA_dir, _fastq_R1_SP, _fastq_R2_SP, transcriptome, _barcode_file, _processing ->
        tuple(id, name, transcriptome, fastq_cDNA_dir)
    }
    cellranger_out = CELLRANGER(cellranger_in)
    cellbender_out = CELLBENDER(cellranger_out)
    cellbender_processed = CELLBENDERPOSTPROCESSING(cellbender_out)
    
    // CURIOTREKKERMODIFY()
    
    trekkersample_in = sampleChannel.map {
        _id, name, _fastq_cDNA_dir, fastq_R1_SP, fastq_R2_SP, _transcriptome, barcode_file, _processing ->
        tuple(name, fastq_R1_SP, fastq_R2_SP, barcode_file)
    }
    trekkersamplesheet = CURIOTREKKERSAMPLESETUP(trekkersample_in, cellbender_processed)
    CURIOTREKKER(trekkersamplesheet)
}