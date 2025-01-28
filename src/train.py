import pandas as pd
import argparse
import config
import numpy as np
import copy

def read_tss_tes_features(fname_tss, fname_tes):
    colnames = ["chrom","strand", "position", "weight_sg", "weight_berth", "read_density", "leading_clip_length", "trailing_clip_length", "junction_start_cnt", "junction_end_cnt", "junction_cross_cnt"]
    tss_data = pd.read_csv(fname_tss, sep='\t', header=None, names=colnames, index_col = False)
    tes_data = pd.read_csv(fname_tes, sep='\t', header=None, names=colnames, index_col = False)
    
    return tss_data, tes_data


def get_labels(ref_annotation, tss_data, tes_data, cfg):
    tss_labels = []
    tes_labels = []

    def binary_search(arr, x):
        lo, hi = 0, len(arr) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if arr[mid] == x:
                return mid
            elif arr[mid] < x:
                lo = mid + 1
            else:
                hi = mid - 1
        return lo

    def find_neighborhood(ref_positions, position, window):
        start = position - window
        end = position + window
        idx = binary_search(ref_positions, position)
        idx = min(idx, len(ref_positions) - 1)
        idx = max(idx, 0)
        # print(f"Position: {position}, Start: {start}, End: {end}, Index: {idx}")
        # Check left side of the position
        for i in range(idx, -1, -1):
            if ref_positions[i] < start:
                break
            if start <= ref_positions[i] and ref_positions[i] <= end:
                return True
        
        # Check right side of the position
        for i in range(idx, len(ref_positions)):
            if ref_positions[i] > end:
                break
            if start <= ref_positions[i] and ref_positions[i] <= end:
                return True
        
        return False
    
    ref_tss = ref_annotation[ref_annotation['type'] == 'TSS'].sort_values(by='position')
    ref_tes = ref_annotation[ref_annotation['type'] == 'TES'].sort_values(by='position')

    print(f"Starting TSS: {tss_data.shape[0]}")
    for it, row in tss_data.iterrows():
        chrom = row['chrom']
        position = row['position']
        
        # Filter ref_annotation for the same chromosome
        ref_chrom_data = ref_tss[ref_tss['chrom'] == chrom]['position'].values

        if ref_chrom_data.size > 0:
            ref_positions = ref_chrom_data
            if find_neighborhood(ref_positions, position, cfg.eval_neighbourhood):
                tss_labels.append(1)
            else:
                tss_labels.append(0)
        else:
            tss_labels.append(0)
        
        if it % 5000 == 0:
            print(f"Processed {it} TSS")
    
    print(f"Starting TES: {tes_data.shape[0]}")
    for it, row in tes_data.iterrows():
        chrom = row['chrom']
        position = row['position']
        
        # Filter ref_annotation for the same chromosome
        ref_chrom_data = ref_tes[ref_tes['chrom'] == chrom]['position'].values
        
        if ref_chrom_data.size > 0:
            ref_positions = ref_chrom_data
            if find_neighborhood(ref_positions, position, cfg.eval_neighbourhood):
                tes_labels.append(1)
            else:
                tes_labels.append(0)
        else:
            tes_labels.append(0)
        
        if it % 5000 == 0:
            print(f"Processed {it} TES")
    
    tss_data['label'] = tss_labels
    tes_data['label'] = tes_labels

    print(tss_data.describe())
    print(tes_data.describe())
    
    tss_data.to_csv("tss_data.csv", index=False)
    tes_data.to_csv("tes_data.csv", index=False)
    return tss_labels, tes_labels


def main():
    # two command line arguments to take tss and tes filename as input with --tss and --tes
    parser = argparse.ArgumentParser(description='Process TSS and TES feature files.')
    parser.add_argument('--tss', required=True, help='Path to the TSS feature file')
    parser.add_argument('--tes', required=True, help='Path to the TES feature file')
    # argument for reference annotation with --ref
    parser.add_argument('--ref', required=True, help='Path to the reference annotation file')
    args = parser.parse_args()
    print(args)
    cfg = config.config()

    tss_data, tes_data = read_tss_tes_features(args.tss, args.tes)
    ref_annotation = pd.read_csv(args.ref, sep=' ', header=None, names=["type", "chrom", "position", "+_strand_transcripts" , "-_strand_transcripts"], index_col = False)
    print(ref_annotation.head())
    # tss_data_sample = copy.deepcopy(tss_data.iloc[0:1000])
    # tes_data_sample = copy.deepcopy(tes_data.iloc[0:1000])
    # tss_labels, tes_labels = get_labels(ref_annotation, tss_data_sample, tes_data_sample, cfg)
    tss_labels, tes_labels = get_labels(ref_annotation, tss_data, tes_data, cfg)
    

    return 0

if __name__ == "__main__":
    main()