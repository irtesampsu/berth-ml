import os
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
# import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, classification_report
from preprocess import read_tss_tes_data
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import StratifiedKFold
from imblearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

config={"oversample": 0.4}

def train_imbalanced_models(X, y, n_splits=5):
    """
    Train multiple models on imbalanced dataset using stratified CV
    
    Parameters:
    X: Features array
    y: Target array
    n_splits: Number of CV folds
    """
    # Initialize models
    models = {
        'RandomForest': RandomForestClassifier(class_weight='balanced', n_estimators=400, random_state=42)
        # 'XGBoost': XGBClassifier(scale_pos_weight=len(y[y==0])/len(y[y==1]), random_state=42),
        # 'SVM': SVC(class_weight='balanced', probability=True, random_state=42)
        # 'LogisticRegression': LogisticRegression(class_weight='balanced', random_state=42)
    }
    
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    # Store all predictions
    model_predictions = {name: [] for name in models.keys()}
    all_true_labels = []
    predict_df = np.zeros((len(y), 2))
    
    # SMOTE with sampling_strategy='auto' will automatically balance classes to 1:1
    if config["oversample"] > 0:
        smote = SMOTE(random_state=42, sampling_strategy=config["oversample"])
    
    shuffled_indices = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y)):
        shuffled_indices.extend(val_idx)
        print(f"Fold {fold + 1}/{n_splits}")
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        # Apply SMOTE only on training data
        if config["oversample"] > 0:
            X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
        else:
            X_train_resampled, y_train_resampled = X_train, y_train

        # Train each model
        for name, model in models.items():
            print(f"Training {name}...")
            model.fit(X_train_resampled, y_train_resampled)
            y_pred = model.predict(X_val)
            model_predictions[name].extend(y_pred)
            
            # Only store true labels once per fold
            if name == list(models.keys())[0]:
                all_true_labels.extend(y_val)
            print(f"Training {name} finished")

    # Generate reports for each model
    print("Classification Reports:")
    for name in models.keys():
        print(f"\n{name}:")
        predict_df[shuffled_indices,0] = all_true_labels
        predict_df[shuffled_indices,1] = model_predictions[name]
        # print(classification_report(all_true_labels, model_predictions[name]))
        print(classification_report(predict_df[:,0], predict_df[:,1]))
        
    return models, predict_df


def process(data):
    # Combine TSS and TES data
    # data = pd.concat([tss_data, tes_data])
    
    X = data.drop(columns=['label']).values
    y = data['label'].values
    
    #train the model
    _, predict_df = train_imbalanced_models(X, y)
    return predict_df

def write_predictions(predict_df, input_file, output_file, ptype='tss'):
    # Write predictions to file
    df = pd.read_csv(input_file)
    df['label'] = predict_df[:,0]
    df['prediction'] = predict_df[:,1]
    df['name'] = ptype + predict_df[:,0].astype(int).astype(str) + '_' + predict_df[:,1].astype(int).astype(str)
    df.to_csv(output_file, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser(description='Process TSS and TES data files.')
    parser.add_argument('--tss', required=True, help='Path to the TSS data file')
    parser.add_argument('--tes', required=True, help='Path to the TES data file')
    parser.add_argument('--oversample', type=float, default=0.4, help='Oversampling ratio')
    args = parser.parse_args()
    print(args)
    config["oversample"] = args.oversample

    tss_data, tes_data = read_tss_tes_data(args.tss, args.tes)
    tss_ones = tss_data['label'].sum()
    tes_ones = tes_data['label'].sum()
    tss_zeros = len(tss_data) - tss_ones
    tes_zeros = len(tes_data) - tes_ones
    print(f"TSS: {tss_ones} positives, {tss_zeros} negatives")
    print(f"TES: {tes_ones} positives, {tes_zeros} negatives")
    # print(tss_data.describe())
    # print(tes_data.describe())
    # print(tss_data.isnull().sum())
    tes_data.fillna(0, inplace=True)
    # print(tes_data.isnull().sum())

    # print(tes_data[tes_data['leading_clip_length'].isnull() & tes_data['trailing_clip_length'].isnull() & tes_data['weight_sg'] > 0])
    print("Training on TSS data")
    predict_df = process(tss_data) #drop some colums of tss_data
    write_predictions(predict_df, args.tss, 'out/tss/tss_predictions.tsv', 'tss')

    print("------------------Done with TSS data------------------")
    print("Training on TES data")
    predict_df  = process(tes_data)
    write_predictions(predict_df, args.tes, 'out/tes/tes_predictions.tsv','tes')
    print("------------------Done with TES data------------------")

if __name__ == "__main__":
    print(os.getcwd())
    main()