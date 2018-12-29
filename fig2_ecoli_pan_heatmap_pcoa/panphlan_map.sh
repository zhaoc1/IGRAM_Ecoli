for file in `cat samples_panphlan.txt`; do
        sample=$(echo $file | sed 's/\..*//')
        echo $sample
	input="decontam/${sample}.fastq"
	echo $input
	python local/panphlan/panphlan_map.py -c Ecoli -i ${input} --fastx fastq -o "${sample}.csv" --verbose
done
