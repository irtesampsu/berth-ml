# conda init
# conda activate irtesam-berth
if [ $# -lt 2 ]; then
    echo "Not enough arguments supplied"
    exit 1
fi



data_home="/datadisk1/ixk5174/long_reads_compare/long_reads_out/berth/"
ref_annotation="${data_home}refSeq.tsstes"


if [ "$1" == "preprocess" ]; then
    python src/preprocess.py --ref "$ref_annotation" --tss "${data_home}tss_merged_features-$2.tsv" --tes "${data_home}tes_merged_features-$2.tsv" --mapping "data/GRCh38_RefSeq2UCSC.txt"
elif [ "$1" == "train" ]; then
    python src/train.py --tss "data/tss_data.csv" --tes "data/tes_data.csv" --oversample -1 > "logs/train-log-cv-current.txt" 
elif [ "$1" == "benchmark" ]; then
    methods=("stringtie" "isoquant")
    if [ $# -lt 2 ]; then
        echo "No method specified"
        exit 1
    fi
    input_file="${data_home}${methods[$2]}.tsstes"
    python src/benchmark.py --ref "$ref_annotation" --input_file "$input_file" --mapping "data/GRCh38_RefSeq2UCSC.txt" --tss_predictions "out/tss/tss_predictions.tsv" --tes_predictions "out/tes/tes_predictions.tsv" > "logs/benchmark-log-${methods[$2]}.txt" --method "${methods[$2]}"
else
    echo "Invalid argument"
    exit 1
fi