params.reference = "mouse"

process setupBender {
	publishDir "PipelineFolder/CellBender/"

	output:
	path "CellBender"
	path "cellbender_latest.sif"

	script:
	"""
	git clone https://github.com/broadinstitute/CellBender.git # to clone repo
	module load apptainer # so we can pull docker image
	apptainer pull docker://us.gcr.io/broad-dsde-methods/cellbender:latest
	"""
}

process setupTrekker {
	publishDir "PipelineFolder/CurioTrekker/"

	output:
	path "curiotrekker-v1.1.0"

	script:
	"""
	wget https://curiotrekkerbioinformatics.s3.us-west-1.amazonaws.com/CurioTrekker_v1.1.0/curiotrekker-v1.1.0.tar.gz -O - |\
tar -xzf -
	"""
}

process setupMouseReference {
	publishDir "PipelineFolder/References/Mouse/"

	output:
	path "Mus_musculus_GRCm39"

	script:
	"""
	wget https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Mus_musculus_GRCm39.tar.gz -O - | tar -xzf -
	"""
}

process setupRatReference {
	publishDir "PipelineFolder/References/Rat/"

	output:
	path "mRatBN7"

	script:
	"""
	wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mRatBN7-2-2024-A.tar.gz -O - | tar -xzf -
	"""
}

workflow {
    setupBender()
    setupTrekker()
    
    // if (params.reference == "rat") {
    //     setupRatReference()
    // }
    // else {
    //     setupMouseReference()
    // }
}
