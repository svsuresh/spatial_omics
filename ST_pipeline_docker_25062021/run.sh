#! /usr/bin/env bash
sudo docker run -it -v ~/references:/home/docker_user/references \
	-v ~/raw_data:/home/docker_user/raw_data \
	-v ~/st_pipeline_docker_05062021/docker_st_results_05062021/tmp:/home/docker_user/results/tmp  \
	-v ~/st_pipeline_docker_05062021/docker_st_results_05062021:/home/docker_user/results \
	st_docker st_pipeline_run.py \
	--output-folder results \
	--ids raw_data/spatial_barcodes_index.txt \
	--ref-map references/star_index \
	--ref-annotation references/GCF_000001635.27_GRCm39_genomic.gtf \
	--expName st_docker_test \
	--htseq-no-ambiguous \
	--verbose \
	--log-file results/st_docker_test_log.txt \
	--temp-folder results/tmp/  \
	--demultiplexing-kmer 5 \
	--threads 20 \
	--no-clean-up \
	--umi-start-position 16 \
	--umi-end-position 26 \
	--demultiplexing-overhang 0 \
	--min-length-qual-trimming 10 E11-2L_S41_L004_R1_001.fastq.gz E11-2L_S41_L004_R2_001.fastq.gz
