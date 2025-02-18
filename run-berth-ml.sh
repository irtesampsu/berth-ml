# conda init
# conda activate irtesam-berth
if [ $# -eq 0 ]; then
    echo "No arguments supplied"
    exit 1
fi

data_home="/datadisk1/ixk5174/long_reads_compare/long_reads_out/berth/"
ref_annotation="${data_home}refSeq.tsstes"

if [ "$1" == "preprocess" ]; then
    python src/preprocess.py --ref "$ref_annotation" --tss "${data_home}tss_merged_features.tsv" --tes "${data_home}tes_merged_features.tsv"
elif [ "$1" == "train" ]; then
    python src/train.py --tss "data/tss_data.csv" --tes "data/tes_data.csv" --oversample 0.4 > "logs/train-log-cv.txt" 
elif [ "$1" == "benchmark" ]; then
    methods=("stringtie" "isoquant")
    if [ $# -lt 2 ]; then
        echo "No method specified"
        exit 1
    fi
    input_file="${data_home}${methods[$2]}.tsstes"
    python src/benchmark.py --ref "$ref_annotation" --input_file "$input_file"
else
    echo "Invalid argument"
    exit 1
fi