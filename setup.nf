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

process setupContainer {
	publishDir "PipelineFolder/Containers/"

	output:
	path "slide_tag_env_latest.sif"

	script:
	"""
	module load StdEnv/2023 apptainer/1.3.4
	apptainer pull docker://hmdh22/slide_tag_env:latest
	"""
}

workflow {
    setupBender()
    setupTrekker()
    setupContainer()
}
