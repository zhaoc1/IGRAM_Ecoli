for file in *.fa; do
        SAMPLEID="${file%.*}"
        echo $SAMPLEID
	OUTPUT_DIR="Result-$SAMPLEID"
	echo $OUTPUT_DIR
	GENOME_DIR="$SAMPLEID.fa"
        echo $GENOME_DIR
	prokka $GENOME_DIR --outdir $OUTPUT_DIR --prefix $SAMPLEID 
done
