import sys
import pickle
from collections import namedtuple
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import numpy as np
from tqdm import tqdm
import xgboost # need version >= 1.6 for multi-label regression

from sklearn.metrics import roc_auc_score
from scipy.stats import pearsonr, spearmanr, kendalltau, linregress

from load_data import averageKL, Oligo, parseIndelIdentifier

np.random.seed(42)


def get_dataframe(data):
    training_cols = ['indelLength', 'potIndelStartRel', 'mosaics_energy_seq', 'mosaics_energy_torsionDelta', 'mosaics_energy_NHEJ', 'mosaics_energy_MMEJ_1', 'mosaics_energy_MMEJ_2']

    df_main = None
    for oligo in tqdm(data, total=len(data)):
        df = oligo.get_df()
        if df is None:
            print('-', end='')
            continue
        
        # calculate MOSAICS energy features
        df['mosaics_energy_torsionDelta'] = df['mosaics_energy_noTorsion_val'] - df['mosaics_energy_val']
        df['mosaics_energy_NHEJ'] = df['mosaics_energy_XZ_val'] + df['mosaics_energy_Y_val'] - df['mosaics_energy_R_val'] - df['mosaics_energy_L_val']
        df['mosaics_energy_MMEJ_1'] = df['mosaics_energy_R2_val'] + df['mosaics_energy_L2_val'] + df['mosaics_energy_r1_val'] + df['mosaics_energy_r2_val'] + df['mosaics_energy_l1_val'] + df['mosaics_energy_l2_val'] - df['mosaics_energy_R_val'] - df['mosaics_energy_L_val']
        df['mosaics_energy_MMEJ_2'] = df['mosaics_energy_RL2_val'] - df['mosaics_energy_R2_val'] - df['mosaics_energy_L2_val']
        
        df_train = df[training_cols+['relFreq']].sort_values(['indelLength', 'potIndelStartRel']).dropna(subset=['mosaics_energy_torsionDelta'])
        df_train['oligoNumber'] = oligo.oligoNumber

        # get one-hot encoding of sequence
        mosaics_seq_len = len(df_train['mosaics_energy_seq'].iloc[0])

        col_list = []
        for i in range(mosaics_seq_len):
            col_list.append('seq_'+str(i))

        df_train[col_list] = df_train['mosaics_energy_seq'].apply(lambda x: pd.Series(list(x)))
        df_train.drop('mosaics_energy_seq', axis=1, inplace=True)

        for col in col_list:
            onehot = pd.get_dummies(df_train[col], prefix=col)
            df_train[onehot.columns] = onehot
            df_train.drop(col, axis=1, inplace=True)

        # make sure all one-hot combinations appear
        for i in range(mosaics_seq_len):
            for base in ['A', 'T', 'C', 'G']:
                col = 'seq_'+str(i)+'_'+base
                if col not in df_train.columns:
                    df_train[col] = 0
        
        # pivot to one line per oligo
        df_train = df_train.pivot(columns=['indelLength', 'potIndelStartRel'], index='oligoNumber')
        
        # extract label
        label_loc = df_train['relFreq']
        df_train.drop('relFreq', axis=1, inplace=True)
        
        # append to training dataframe
        if df_main is None: 
            df_main = df_train
            label = pd.DataFrame(label_loc)
            oligo_numbers = [oligo.oligoNumber]
        else: 
            df_main = pd.concat([df_main, df_train], ignore_index=True)
            label = label.append(label_loc, ignore_index=True)
            oligo_numbers.append(oligo.oligoNumber)

    # add zero column for indelLength=0 and potIndelStartRel=-4 (is this necessary?)
    label[(0, -4)] = 0

    label = label.fillna(0)

    return df_main, label


def evaluate_model(model, dtest, y_test=None, sum_level=-1, contrast=1, classification=False, do_plots=True):        
    predicted_quantities = {None: "indel length and position", 0: "indel length", 1: "indel position", -1: ""}
    predicted_quantity = predicted_quantities[sum_level]
    
    if classification:
        # calculate AUC for each label
        for i, col in enumerate(y_test.T):
            print(i, roc_auc_score(col, model.predict(dtest)[:, i]))
        
    else:
        test_pred = model.predict(dtest)
        if do_plots:
            plt.hist(np.argsort(test_pred, axis=1)[:, -2], bins=21)
            plt.hist(np.argsort(y_test, axis=1)[:, -2], bins=21, alpha=0.5)
            plt.xlabel("second most frequent "+predicted_quantity)
            plt.ylabel("frequency")
            plt.show()

        # confusion matrix
        if do_plots:
            from sklearn.metrics import confusion_matrix
            conf_matrix = confusion_matrix(np.argsort(y_test, axis=1)[:, -2], np.argsort(test_pred, axis=1)[:, -2])
            plt.imshow(conf_matrix,
                    vmax=1/contrast*np.max(conf_matrix))
            plt.show()

        # look at distances between predicted argmax and ground truth
        if do_plots:
            dists = np.argsort(test_pred, axis=1)[:, -2] - np.argsort(y_test, axis=1)[:, -2]
            plt.hist(dists, bins=len(np.unique(dists)))
            plt.show()
            print(((dists <= 3) & (dists >= -3)).sum() / len(dists)*100, "% of predictions are 3bp or less off from the true second argmax")
        
        kl_div, kl_div_dist = averageKL(pd.DataFrame.from_records(test_pred), pd.DataFrame.from_records(y_test))
        print("KL divergence between predicted and true", predicted_quantity, "distribution:", kl_div)
        
        pearson, spearman = pearsonr(y_test.flatten(), test_pred.flatten())[0], spearmanr(y_test.flatten(), test_pred.flatten())[0]
        print("Pearson correlation", pearson)
        print("Spearman correlation", spearman)

        kendall_tau = kendalltau(y_test.flatten(), test_pred.flatten())[0]

        r_squared = linregress(y_test.flatten(), test_pred.flatten())[2]**2
        
        return (kl_div, pearson, spearman, kendall_tau, r_squared), kl_div_dist



if __name__ == "__main__":
    sum_level   = None if len(sys.argv) <= 1 else (int(sys.argv[1]) if sys.argv[1] != "None" else None)
    feature_set = 0    if len(sys.argv) <= 2 else int(sys.argv[2])
    
    use_only_energy_feat = True
    normalise_df = True # normalise the MOSAICS energy columns
    normalise_label = False
    classification = False
    remove_some_zero_points = True # undersample points with highest argmax

    print("loading model...", end='')
    models = pickle.load(open("models.pickle", "rb"))

    # sum_level: 0: indel length, 1: indel position, None: length and position
    # feature_set: 0: sequence + energy, 1: energy only, 2: sequence only
    configs = [{'sum_level': None, 'feature_set': 0},
            {'sum_level': None, 'feature_set': 1},
            {'sum_level': 0,    'feature_set': 0},
            {'sum_level': 0,    'feature_set': 1},
            {'sum_level': 1,    'feature_set': 0},
            {'sum_level': 1,    'feature_set': 1},
            {'sum_level': None, 'feature_set': 2},
            {'sum_level': 0,    'feature_set': 2},
            {'sum_level': 1,    'feature_set': 2},
            ]

    # retrieve model according to user input
    model_idx = -1
    for i, c in enumerate(configs):
        if c["sum_level"] == sum_level and c["feature_set"] == feature_set:
            model_idx = i
            break

    if model_idx == -1: 
        print("combination of sum_level, feature_set not found - exiting")
        exit()

    config = configs[model_idx]
    reg = models[model_idx]
    print("done")


    print("loading MOSAICS data...")
    data = pickle.load(open("mosaics_data.pickle", "rb"))

    df_main, label = get_dataframe(data)

    if normalise_df: # min-max normalisation of the energy features (as a whole dataframe - keep units the same across energy columns)
        energy_min, energy_max = df_main.min().min(), df_main.max().max()
        
        for col in df_main.columns:
            if 'mosaics_energy' in col[0]:
                df_main[col] = (df_main[col]-energy_min)/(energy_max-energy_min)
            

    if normalise_label:
        label = label.div(label.sum(axis=1), axis=0).fillna(0)
    else:
        if (label.sum(axis=1) > 1.0).sum() > 0:
            print("ignoring", (label.sum(axis=1) > 1.0).sum(), "data points which have more than 1000 reads")
        mask = label.sum(axis=1) <= 1.0
        df_main, label = df_main[mask], label[mask]
        
        # fill indelLength=0 column with remaining reads
        label.loc[:, (label.columns.get_level_values('indelLength') == 0)] = 1.0-label.sum(axis=1)

    label_loc = label.copy(deep=True)

    # prepare data and label dataframes
    if config['sum_level'] is not None: 
        label_loc = label_loc.sum(axis=1, level=config['sum_level'])

    if config['feature_set'] == 1:
        cols = [c for c in df_main.columns if 'mosaics_energy' in c[0]]
        df_main_loc = df_main[cols]
    elif config['feature_set'] == 2:
        cols = [c for c in df_main.columns if 'mosaics_energy' not in c[0]]
        df_main_loc = df_main[cols]
    else:
        df_main_loc = df_main
        
    # get test set to predict on
    X_test, y_test, idx_test = df_main_loc.to_numpy(), label_loc.to_numpy(), df_main.index.values
    print("done")

    # obtain model prediction
    print("obtaining model prediction...", end='')
    pred = reg.predict(X_test)
    print("done")

    benchmarks = evaluate_model(reg, X_test, y_test, config['sum_level'], 8, do_plots=False)[0]

    # save dataframe
    pd.DataFrame(pred, columns=label_loc.columns).to_csv("output.csv", index_label=label_loc.columns.name)
    print ("Prediction saved successfully")
