import pandas as pd
import os

def set_home():
    os.chdir("/datadisk1/ixk5174/tools/berth-ml/")

prediction_folder_tss = "out/tss/"
prediction_folder_tes = "out/tes/"
data_folder = "data/"

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

if __name__ == "__main__":
    set_home()
    main()
    # analyze_data()
    print("Done")