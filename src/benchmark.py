import time
import argparse
import pandas as pd
import numpy as np
import preprocess
from sklearn.metrics import accuracy_score, classification_report
import config
import os

prediction_folder_tss = "out/tss/"
prediction_folder_tes = "out/tes/"
data_folder = "data/"
# method = "stringtie"
output_folder = "/datadisk1/ixk5174/long_reads_compare/long_reads_out/berth"


def filter_unknown(method, output_folder):
    tss_predictions_path = os.path.join(prediction_folder_tss, f'tss_predictions_{method}.tsv')
    tss_predictions_unkown_path = os.path.join(output_folder, f'tss_predictions_unk_{method}.bed')

    tes_predictions_path = os.path.join(prediction_folder_tes, f'tes_predictions_{method}.tsv')
    tes_predictions_unkown_path = os.path.join(output_folder, f'tes_predictions_unk_{method}.bed')

    tss_predictions = pd.read_csv(tss_predictions_path, sep='\t')
    tss_predictions_unk = tss_predictions[tss_predictions['prediction'] == 2 ].copy()
    tss_predictions_unk.loc[:, "end_position"] = tss_predictions_unk["position"] + 50
    tss_predictions_unk.loc[:, "strand"] = np.where(tss_predictions_unk["+_strand_transcripts"] >= tss_predictions_unk["-_strand_transcripts"], '+', '-')
    tss_predictions_unk.loc[:,"name"] = "" + tss_predictions_unk["position"].astype(str) + "__" + tss_predictions_unk["strand"] + "__" + tss_predictions_unk["label"].astype(str)
    tss_predictions_unk = tss_predictions_unk[["chrom","position", "end_position","name", "label", "strand"]]
    tss_predictions_unk.to_csv(tss_predictions_unkown_path, index=False, sep='\t', header=False)

    tes_predictions = pd.read_csv(tes_predictions_path, sep='\t')
    tes_predictions_unk = tes_predictions[tes_predictions['prediction'] == 2].copy()
    tes_predictions_unk.loc[:, "start_position"] = tes_predictions_unk["position"] - 50
    tes_predictions_unk.loc[:, "strand"] = np.where(tes_predictions_unk["+_strand_transcripts"] >= tes_predictions_unk["-_strand_transcripts"], '+', '-')
    tes_predictions_unk.loc[:,"name"] = "" + tes_predictions_unk["position"].astype(str) + "__" + tes_predictions_unk["strand"] + "__" + tes_predictions_unk["label"].astype(str)
    tes_predictions_unk = tes_predictions_unk[["chrom", "start_position", "position", "name", "label", "strand"]]
    tes_predictions_unk.to_csv(tes_predictions_unkown_path, index=False, sep='\t', header=False)

    # return tss_predictions_unk, tes_predictions_unk

def get_neighbourhood_predictions(berth_predictions, benchmark_data, cfg):
    #sort berth_data based on chrom and position
    berth_data = berth_predictions.sort_values(by=['chrom', 'position'])
    benchmark_predictions = []
    for _, row in benchmark_data.iterrows():
        chrom = row['chrom']
        position = row['position']
        berth_chrom = berth_data[berth_data['chrom'] == chrom]
        berth_positions = berth_chrom['position'].values
        idx = preprocess.binary_search(berth_positions, position)

        # assert idx > 0 and idx < len(berth_positions), f"Index {idx}:({position} --> {berth_positions[0]}, {berth_positions[-1]} :: {len(berth_positions)}) out of bounds for binary search" # check if the index is within bounds
        idx = min(idx, len(berth_positions) - 1)
        idx = max(idx, 0)

        window_start = position - cfg.eval_neighbourhood
        window_end = position + cfg.eval_neighbourhood
        left_idx = max(idx-1, 0)
        right_idx = min(idx+1, len(berth_positions)-1)
        neighbor_preds = []
        if berth_positions[idx] >= window_start and berth_positions[idx] <= window_end:
            neighbor_preds.append(berth_chrom.iloc[idx]['prediction'])

        while left_idx >= 0 and berth_positions[left_idx] >= window_start:
            neighbor_preds.append(berth_chrom.iloc[left_idx]['prediction'])
            left_idx -= 1

        while right_idx < len(berth_positions) and berth_positions[right_idx] <= window_end:
            neighbor_preds.append(berth_chrom.iloc[right_idx]['prediction'])
            right_idx += 1

        if len(neighbor_preds) == 0:
            benchmark_predictions.append(2)
        else:
            # assign if any of the neighbors is 1 // TODO: check if this is the correct way to assign the label
            benchmark_predictions.append(1 if 1 in neighbor_preds else 0)

        if idx % 5000 == 0:
            print(f"Processed {idx}")
        
    return benchmark_predictions


def main():
    # start timer 
    start = time.time()

    # ------------------ Initialize Arguments ------------------
    parser = argparse.ArgumentParser(description='Benchmark TSS and TES feature files.')
    parser.add_argument('--ref', type=str, help='Reference annotation file', required=True)
    parser.add_argument('--input_file', type=str, help='List of feature files', required=True)
    parser.add_argument('--mapping', required=True, help='Path to the chromosome mapping file')
    parser.add_argument('--tss_predictions', required=True, help='Path to the TSS predictions file')
    parser.add_argument('--tes_predictions', required=True, help='Path to the TES predictions file')
    parser.add_argument('--method', type=str, help='Method used for benchmarking', required=True)
    # Parse the arguments and set parameters
    args = parser.parse_args()
    cfg = config.config()
    cfg.eval_neighbourhood = 100
    # print(cfg.save_data)

    # ------------------ Read and Map Reference Data ------------------
    ref_annotation = pd.read_csv(args.ref, sep=' ', header=None, names=["type", "chrom", "position", "+_strand_transcripts" , "-_strand_transcripts"], index_col = False)
    chrom_mapping = pd.read_csv(args.mapping, sep='\t', header=None, names=["chrom", "mapped"], index_col = False)
    # convert it to a dictionary
    chrom_mapping = chrom_mapping.set_index('chrom')['mapped'].to_dict()
    ref_annotation = preprocess.map_chromosomes(ref_annotation, chrom_mapping)
    print(f"Time taken to read and map reference data: {time.time() - start} seconds")

    # ------------------ Read Benchmark Data ------------------
    input_annotation = pd.read_csv(args.input_file, sep=' ', header=None, names=["type", "chrom", "position", "+_strand_transcripts" , "-_strand_transcripts"], index_col = False)
    benchmark_tss_data = input_annotation[input_annotation['type'] == 'TSS'].reset_index(drop=True).copy()
    benchmark_tes_data = input_annotation[input_annotation['type'] == 'TES'].reset_index(drop=True).copy()

    # ----------------- Get TSS and TES labels -----------------
    tss_labels, tes_labels = preprocess.get_labels(ref_annotation,benchmark_tss_data, benchmark_tes_data,cfg)
    benchmark_tss_data['label'] = tss_labels
    benchmark_tes_data['label'] = tes_labels
    print("Benchmark Data(TSS):")
    print(benchmark_tss_data.head())
    print("Benchmark Data(TES):")
    print(benchmark_tes_data.head())
    print(f"Time taken to get TSS and TES labels: {time.time() - start} seconds")
    
    # ----------------- Read Predictions -----------------
    tss_predictions = pd.read_csv(args.tss_predictions, sep='\t')
    tes_predictions = pd.read_csv(args.tes_predictions, sep='\t')

    # ----------------- Get Neighbourhood Predictions -----------------
    tss_preds = get_neighbourhood_predictions(tss_predictions, benchmark_tss_data, cfg)
    tes_preds = get_neighbourhood_predictions(tes_predictions, benchmark_tes_data, cfg)
    benchmark_tss_data['prediction'] = tss_preds
    benchmark_tes_data['prediction'] = tes_preds
    # print prediction class frequency distribution
    print("TSS Prediction Class Distribution:")
    print(benchmark_tss_data['prediction'].value_counts())
    print("TES Prediction Class Distribution:")
    print(benchmark_tes_data['prediction'].value_counts())
    
    # ----------------- Write Predictions -----------------
    benchmark_tss_data.to_csv(f'out/tss/tss_predictions_{args.method}.tsv', sep='\t', index=False)
    benchmark_tes_data.to_csv(f'out/tes/tes_predictions_{args.method}.tsv', sep='\t', index=False)
    print(f"Time taken to get neighbourhood predictions: {time.time() - start} seconds")    


    filter_unknown(args.method, output_folder)

    # ----------------- Evaluate Predictions -----------------
    benchmark_tss_data['prediction'] = benchmark_tss_data['prediction'].replace(2, 1)
    benchmark_tes_data['prediction'] = benchmark_tes_data['prediction'].replace(2, 1)
    print("TSS Classification Report")
    print(classification_report(benchmark_tss_data['label'], benchmark_tss_data['prediction']))
    print("TES Classification Report")
    print(classification_report(benchmark_tes_data['label'], benchmark_tes_data['prediction']))


    print(f"Time taken to write predictions: {time.time() - start} seconds")

if __name__ == '__main__':
    main()