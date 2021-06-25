#! /usr/bin/env bash
sudo docker run -it -v ~/references:/home/docker_user/references \
	-v $(readlink -f $1):/home/docker_user/raw_data \
	-v $(readlink -f $2):/home/docker_user/results \
	-v $(readlink -f $2)/tmp:/home/docker_user/results/tmp  \
	st_docker st_pipeline_run.py \
	--output-folder results \
	--ids raw_data/spatial_barcodes_index.txt \
	--ref-map references/star_index \
	--ref-annotation references/GCF_000001635.27_GRCm39_genomic.gtf \
	--expName $5 \
	--htseq-no-ambiguous \
	--verbose \
	--log-file results/$5_logs.txt \
	--temp-folder results/tmp/  \
	--demultiplexing-kmer 5 \
	--threads 20 \
	--no-clean-up \
	--umi-start-position 16 \
	--umi-end-position 26 \
	--demultiplexing-overhang 0 \
	--min-length-qual-trimming 10 raw_data/$3 raw_data/$4
