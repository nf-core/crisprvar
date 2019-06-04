test:
	nextflow run main.nf \
		-profile docker \
		--reads 's3://czb-seqbot/fastqs/190326_M05295_0256_000000000-CC47Y/rawdata/*{R1,R2}*.fastq.gz' \
		--samplesheet test-data/samplesheet.csv -ansi-log false

docker_build:
	docker build -t czbiohub/crispresso .
