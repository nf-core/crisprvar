test:
	nextflow run main.nf \
		-profile docker \
		--reads 'test-data/fastqs/*{R1,R2}*.fastq.gz' \
		--samplesheet test-data/samplesheet.csv -ansi-log false \
		--adapters test-data/TruSeq3-PE-2.fa.txt

docker_build:
	docker build -t czbiohub/crisprvar:dev .

docker_push:
	sudo docker login
	sudo docker push czbiohub/crisprvar:dev
	docker images
