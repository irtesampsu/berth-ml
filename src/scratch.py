import pandas as pd
import os
from matplotlib import pyplot as plt

def set_home():
    os.chdir("/datadisk1/ixk5174/tools/berth-ml/")

prediction_folder_tss = "out/tss/"
prediction_folder_tes = "out/tes/"
data_folder = "data/"
method = "stringtie"
output_folder = "/datadisk1/ixk5174/long_reads_compare/long_reads_out/berth"

def analyze_data():
    tss_data_path = os.path.join(data_folder, 'tss_data.csv')
    tes_data_path = os.path.join(data_folder, 'tes_data.csv')

    tss_data = pd.read_csv(tss_data_path)
    tes_data = pd.read_csv(tes_data_path)

    print("-----------------------TSS Data-----------------------")
    tss_data = tss_data.sort_values(by='position')
    # print(tss_data.head())
    print(tss_data[tss_data['position'] == 0][['chrom', 'strand', 'position', 'weight_sg', 'weight_berth' ,'label']])
    print(tss_data[tss_data['position'] == 0].shape)
    print(tss_data[(tss_data['weight_berth'] > 0)  & (tss_data['weight_sg'] == 0) & (tss_data['position'] == 0)].shape)
    print(tss_data[(tss_data['weight_berth'] > 0)  & (tss_data['weight_sg'] == 0) & (tss_data['position'] == 0)].describe())
    # print(tss_data.describe())

    print("-----------------------TES Data-----------------------")
    tes_data = tes_data.sort_values(by='position')
    # print(tes_data.head())
    print(tes_data[tes_data['position'] == 0][['chrom', 'strand', 'position', 'weight_sg', 'weight_berth' ,'label']])
    print(tes_data[tes_data['position'] == 0].shape)
    print(tes_data[(tes_data['weight_berth'] > 0 ) & (tes_data['weight_sg'] == 0) & (tes_data['position'] == 0)].shape)
    print(tes_data[(tes_data['weight_berth'] > 0 ) & (tes_data['weight_sg'] == 0) & (tes_data['position'] == 0)].describe())
    # print(tes_data.describe())


def main():
    tss_predictions_path = os.path.join(prediction_folder_tss, 'tss_predictions.tsv')
    tss_predictions_wrong_path = os.path.join(prediction_folder_tss, 'tss_predictions_wrong.bed')
    tss_predictions_correct_path = os.path.join(prediction_folder_tss, 'tss_predictions_correct.bed')

    tes_predictions_path = os.path.join(prediction_folder_tes, 'tes_predictions.tsv')
    tes_predictions_wrong_path = os.path.join(prediction_folder_tes, 'tes_predictions_wrong.bed')
    tes_predictions_correct_path = os.path.join(prediction_folder_tes, 'tes_predictions_correct.bed')

    tss_predictions = pd.read_csv(tss_predictions_path, sep='\t')
    tss_predictions_wrong = tss_predictions[tss_predictions['label'] != tss_predictions['prediction']].copy()
    tss_predictions_wrong.loc[:, "end_position"] = tss_predictions_wrong["position"] + 50
    tss_predictions_wrong = tss_predictions_wrong[["chrom","position", "end_position","name", "label", "strand"]]
    tss_predictions_wrong.to_csv(tss_predictions_wrong_path, index=False, sep='\t', header=False)

    tes_predictions = pd.read_csv(tes_predictions_path, sep='\t')
    tes_predictions_wrong = tes_predictions[tes_predictions['label'] != tes_predictions['prediction']].copy()
    tes_predictions_wrong.loc[:, "start_position"] = tes_predictions_wrong["position"] - 50
    tes_predictions_wrong = tes_predictions_wrong[["chrom", "start_position", "position", "name", "label", "strand"]]
    tes_predictions_wrong.to_csv(tes_predictions_wrong_path, index=False, sep='\t', header=False)

    tss_predictions_correct = tss_predictions[tss_predictions['label'] == tss_predictions['prediction']].copy()
    tss_predictions_correct.loc[:, "end_position"] = tss_predictions_correct["position"] + 50
    tss_predictions_correct = tss_predictions_correct[["chrom", "position", "end_position", "name", "label", "strand"]]
    tss_predictions_correct.to_csv(tss_predictions_correct_path, index=False, sep='\t', header=False)

    tes_predictions_correct = tes_predictions[tes_predictions['label'] == tes_predictions['prediction']].copy()
    tes_predictions_correct.loc[:, "start_position"] = tes_predictions_correct["position"] - 50
    tes_predictions_correct = tes_predictions_correct[["chrom", "start_position", "position", "name", "label", "strand"]]
    tes_predictions_correct.to_csv(tes_predictions_correct_path, index=False, sep='\t', header=False)


def plot_dist(data_file, col, out):
    df = pd.read_csv(data_file)
    coldata = df[col].to_numpy()
    # uniq_val, cnt = np.unique(coldata, return_counts = True)
    plt.figure()
    plt.scatter( df[col], df['label'])
    # plt.hist(df[col], bins=[0,5,10,20,25,30,50,100,200,400,1000,2000])
    plt.savefig(f'{out}_plot-{col}.jpg', dpi=300)

if __name__ == "__main__":
    
    set_home()
    plot_dist('data/tss_data.csv', 'left_anchor_padding', 'tss')
    plot_dist('data/tss_data.csv', 'right_anchor_padding', 'tss')
    plot_dist('data/tes_data.csv', 'left_anchor_padding', 'tes')
    plot_dist('data/tes_data.csv', 'right_anchor_padding', 'tes')


    # main()
    # filter_unknown(method, output_folder)
    # analyze_data()
    print("Done")