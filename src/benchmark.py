import argparse
import pandas as pd
import preprocess
from sklearn.metrics import accuracy_score, classification_report
import config

def main():
    parser = argparse.ArgumentParser(description='Benchmark TSS and TES feature files.')
    parser.add_argument('--ref', type=str, help='Reference annotation file', required=True)
    parser.add_argument('--input_file', type=str, help='List of feature files', required=True)
    args = parser.parse_args()
    cfg = config.config()
    ref_annotation = pd.read_csv(args.ref, sep=' ', header=None, names=["type", "chrom", "position", "+_strand_transcripts" , "-_strand_transcripts"], index_col = False)
    input_annotation = pd.read_csv(args.input_file, sep=' ', header=None, names=["type", "chrom", "position", "+_strand_transcripts" , "-_strand_transcripts"], index_col = False)
    tss_data = input_annotation[input_annotation['type'] == 'TSS'].reset_index(drop=True).copy()
    tes_data = input_annotation[input_annotation['type'] == 'TES'].reset_index(drop=True).copy()
    print(tss_data.shape)
    print(tes_data.shape)
    tss_labels, tes_labels = preprocess.get_labels(ref_annotation,tss_data, tes_data,cfg)
    print(tss_data.head())
    print(tes_data.head())
    tss_data['label'] = tss_labels
    tes_data['label'] = tes_labels
    tss_data['gt'] = 1
    tes_data['gt'] = 1

    print("TSS data")
    print(accuracy_score(tss_data['gt'].values, tss_data['label'].values))
    print("TES data")   
    print(accuracy_score(tes_data['gt'].values, tes_data['label'].values))

if __name__ == '__main__':
    main()