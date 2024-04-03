##################################################
# TRAIN PREDICTION INTERVAL METHODS FOR SPECIES LIFESPANS AND EVALUATE
# Author: Suk Yee Yong
##################################################

from pathlib import Path
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import life_pi.lib.run_uqregresor as run_uqregresor


def merge_dfdatapred(df_data):
    """
    Aggregate median of predictions from all partitions and merge with df_data
    
    Parameters
    ----------
        df_data: DataFrame of X, y and features
        file_nsplit: File partition number
        output_dir: Output directory
        y_labelname: str of y response variable name
        strategy: Prediction interval methods 'naive', 'cv', 'cqr'
    
    Returns
    ----------
        df_data with columns median 'pred', 'pred_lower', 'pred_upper'
    """
    df_data = pd.concat(df_data)
    df_data_group = df_data.groupby('organism_name', as_index=False)
    df_data = df_data.drop_duplicates('organism_name')\
                        .merge(df_data_group[['pred', 'pred_lower', 'pred_upper']].median(),
                            how='inner', on='organism_name', suffixes=('_x', ''))
    df_data = df_data.drop(df_data.filter(regex='_x').columns, axis=1)
    return df_data


def scale_data(df_data, df_trainvalid, scale=None, save_filename=False):
    """
    Scale data
    
    Parameters
    ----------
        df_data: DataFrame
        df_trainvalid: DataFrame of train and validation data to calculate scaling
        scale: Scale data with 'standardize', 'normalize'. Default without scaling.
        save_filename: Save scaling filename
    
    Returns
    ----------
        df_scaled: Scaled df_data by (df_data - minus_numerator)/denominator
        minus_numerator: Minus scaling value in numerator
        denominator: Scaling value in denominator
    """
    if scale == 'standardize':
        minus_numerator = df_trainvalid.mean(axis=0, skipna=True)
        denominator = df_trainvalid.std(axis=0, skipna=True, ddof=1)
        df_scaled = (df_data - minus_numerator)/denominator
    elif scale == 'normalize':
        minus_numerator, denominator = df_trainvalid.min(), df_trainvalid.max() - df_trainvalid.min()
        df_scaled = (df_data - minus_numerator)/denominator
    else:
        minus_numerator, denominator = 0., 1.
        df_scaled = df_data
        save_filename = False
    
    if isinstance(save_filename, str):
        pd.to_pickle({'minus_numerator': minus_numerator, 'denominator': denominator}, save_filename)
    return df_scaled, minus_numerator, denominator


def invt_scale_data(data, minus_numerator, denominator):
    """
    Inverse transform of scale_data by (data*denominator) + minus_numerator
    
    Parameters
    ----------
        data: Data to be scaled
        minus_numerator: Minus scaling value in numerator from scale_data
        denominator: Scaling value in denominator from scale_data
    
    Returns
    ----------
        Inverse transformed data
    """
    return (data*denominator) + minus_numerator


def estimator_elasticnet(data_dir, filename, file_nsplit, l1ratio_colname, alpha_colname, seed=1):
    """
    Read ElasticNet parameter from file and return scikit-learn elastic net regressor: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
    
    Parameters
    ----------
        data_dir: Data directory
        filename: Data filename
        file_nsplit: File partition number
        l1ratio_colname: Mixing parameter, corresponds to alpha in gmlnet R package
        alpha_colname: Constant that multiplies the penalty terms, corresponds to lambda in gmlnet R package
        seed: Random seed
    
    Returns
    ----------
        scikit-learn ElasticNet
    """
    df_param = pd.read_csv(Path(data_dir, f"{filename}.csv"))
    df_param = df_param.groupby('model_number').get_group(file_nsplit)
    return ElasticNet(l1_ratio=df_param[l1ratio_colname].to_numpy()[0], alpha=df_param[alpha_colname].to_numpy()[0], max_iter=100000, tol=1.e-7, random_state=seed)


def get_dftest(data_dir, filename, file_nsplit, output_dir, y_labelname, strategy, alpha=0.1):
    """
    Return DataFrame of predictions and prediction intervals for test set
    
    Parameters
    ----------
        data_dir: Data directory
        filename: Data filename
        file_nsplit: File partition number
        output_dir: Output directory
        y_labelname: y response variable name
        strategy: Prediction interval methods 'naive', 'cv', 'cqr'
        alpha: 1 - confidence interval
    
    Returns
    ----------
        DataFrame with columns 'pred', 'pred_lower', 'pred_upper'
    """
    
    def load_test(df_data, file_nsplit, output_dir, y_labelname, strategy):
        print(f"PARTITION # {file_nsplit}")
        outputpart_dir = Path(output_dir, str(file_nsplit))
        datasplit_idx = pd.read_pickle(Path(outputpart_dir, 'datasplitidx.pkl'))
        # Get test data and outputs
        df_datapart = df_data.groupby('model_number').get_group(file_nsplit)\
                                .dropna(axis=1, how='all').reset_index(drop=True) 
        df_datapart_test = df_datapart.iloc[datasplit_idx['test']][['organism_name', f"mean_{y_labelname}", f"detrans_predicted_{y_labelname}", f"bagged_detrans_predicted_{y_labelname}"]]
        y_pred, y_pis, uq_dict = pd.read_pickle(Path(outputpart_dir, f"mapie_pred_alpha{alpha}.pkl"))
        # Add outputs to DataFrame
        df_datapart_test['pred'] = y_pred[strategy]
        df_datapart_test['pred_lower'] = y_pis[strategy][:,0]
        df_datapart_test['pred_upper'] = y_pis[strategy][:,1]
        return df_datapart_test
    
    df_data = pd.read_csv(Path(data_dir, f"{filename}.csv"))
    if file_nsplit != 0:
        df_datapart_test = load_test(df_data, file_nsplit, output_dir, y_labelname, strategy)
    else:
        df_datapart_test = []
        for pn in df_data['model_number'].unique():
            df_datapart_test.append(load_test(df_data, pn, output_dir, y_labelname, strategy))
        df_datapart_test = merge_dfdatapred(df_datapart_test)
    return df_datapart_test


def train_pi(df_data, y_labelname, feature_prefix='FP', scale_X=None, scale_y=None, log_y=False, test_size=0.5, seed=1, estimator=None, estimator_q=None, list_strategies='cqr', alpha=0.1, output_dir='output'):
    """
    Train prediction interval methods and save to output file
    
    Parameters
    ----------
        df_data: DataFrame of X, y and features
        y_labelname: str of y response variable name
        feature_prefix: Feature name prefix
        scale_X, scale_y: Scale X/y data with 'standardize', 'normalize'. Default without scaling.
        log_y: Use y in log scale
        test_size: Test vs validation set split fraction
        seed: Random seed
        estimator, estimator_q: Scikit-learn API regressor and regressor with quantile
        list_strategies: list of prediction interval methods ['naive', 'jackknife', 'jackknife_plus', 'jackknife_minmax', 'jackknife_plus_ab', 'cv', 'cv_plus', 'cv_minmax', 'cqr']
        alpha: 1 - confidence interval
        output_dir: Output directory
    """
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mask_train = df_data['initial_split'] == 'trainValidate'
    df_X = df_data.filter(regex=f"^({feature_prefix})")
    df_y = df_data[y_labelname]
    # Scale data
    if scale_X is not None:
        df_X = scale_data(df_X, df_X[mask_train], scale=scale_X, save_filename=str(Path(output_dir, f'datascalerx_{scale_X}.pkl')))[0]
    ytrainvalid_mnumerator, ytrainvalid_denominator = 0., 1.
    if log_y:
        df_y = np.log(df_y)
    if scale_y is not None:
        df_y, ytrainvalid_mnumerator, ytrainvalid_denominator = scale_data(df_y, df_y[mask_train], scale=scale_y, save_filename=str(Path(output_dir, f'datascalery_{scale_y}.pkl')))
    
    # Split to train, validation, test sets
    train_valid_idx = train_test_split(df_data.index[mask_train], test_size=test_size, random_state=seed)
    X_train = df_X.iloc[train_valid_idx[0]].to_numpy()
    y_train = df_y.iloc[train_valid_idx[0]].to_numpy()
    X_valid = df_X.iloc[train_valid_idx[1]].to_numpy()
    y_valid = df_y.iloc[train_valid_idx[1]].to_numpy()
    X_test = df_X[~mask_train].to_numpy()
    y_test = df_y[~mask_train].to_numpy()
    # Save data split index
    datasplit_idx = {'train': train_valid_idx[0], 'valid': train_valid_idx[1], 'test': df_data.index[~mask_train]}
    pd.to_pickle(datasplit_idx, Path(output_dir, 'datasplitidx.pkl'))
    
    # Get prediction intervals and save outputs
    if estimator is None:
        estimator = run_uqregresor.optimize_regressor(X_train, y_train, regressor=None, n_iter=3, cv=5, filename=Path(output_dir, 'estimator_optim.pkl'), seed=seed)
    if (estimator_q is None) and ('cqr' in list_strategies):
        estimator_q = run_uqregresor.optimize_regressor(X_train, y_train, regressor='quantile', n_iter=3, cv=5, filename=Path(output_dir, 'estimatorq_optim.pkl'), seed=seed)
    
    y_pred, y_pis, uq_dict = run_uqregresor.mapie_uq(estimator, estimator_q, X_train, y_train, X_valid, y_valid, X_test, y_test, list_strategies=list_strategies, alpha=alpha, filepath=output_dir, save_outfile=False)
    # Rewrite file with rescaled outputs
    invt_scale = lambda data: invt_scale_data(data, ytrainvalid_mnumerator, ytrainvalid_denominator) if scale_y is not None else data
    rescale_func = lambda v: np.exp(invt_scale(v)) if log_y else invt_scale(v)
    dict_rescale = lambda d: {k: rescale_func(v) for k, v in d.items()}
    y_pred, y_pis = dict_rescale(y_pred), dict_rescale(y_pis)
    uq_dict = {s: {k: rescale_func(v) for k, v in ks.items() if k!='pierr_metric'} for s, ks in uq_dict.items()}
    pd.to_pickle((y_pred, y_pis, uq_dict), Path(output_dir, f"mapie_pred_alpha{alpha}.pkl"))


def pred_testx(df_data, modelpath, feature_prefix='FP', log_y=False, alpha=0.1, output_dir='output/pred', save_filename=False):
    """
    Predict uncertainties in data with features extracted from trained model
    
    Parameters
    ----------
        df_data: DataFrame of X and features
        modelpath: Path to saved trained model
        feature_prefix: Feature name prefix
        log_y: Use y in log scale
        alpha: 1 - confidence interval
        output_dir: Output directory
        save_filename: Save prediction filename
    
    Returns
    ----------
        df_data with columns 'pred', 'pred_lower', 'pred_upper'
    """
    
    df_features = df_data.filter(regex=f"^({feature_prefix})")
    df_X = df_features.copy()
    # Read data scaler X if available
    scalerx_file = list(Path(modelpath).parent.glob('datascalerx_*.pkl'))
    if len(scalerx_file) == 1:
        datascalerx = pd.read_pickle(scalerx_file[0])
        df_X = (df_X - datascalerx['minus_numerator'])/datascalerx['denominator']
    
    # Generate predictions
    mapie = pd.read_pickle(modelpath)
    if 'Quantile' in str(mapie.__class__.__name__):
        y_pred, ypis = mapie.predict(df_X.to_numpy())
    else:
        y_pred, ypis = mapie.predict(df_X.to_numpy(), alpha=alpha)
    y_pis = ypis.squeeze()
    
    # Rescale y
    scalery_file = list(Path(modelpath).parent.glob('datascalery_*.pkl'))
    if len(scalery_file) == 1:
        datascalery = pd.read_pickle(scalery_file[0])
        y_pred = invt_scale_data(y_pred, datascalery['minus_numerator'], datascalery['denominator'])
        y_pis = invt_scale_data(y_pis, datascalery['minus_numerator'], datascalery['denominator'])
    if log_y:
        y_pred, y_pis = np.exp(y_pred), np.exp(y_pis)
    ylower, yupper = y_pis[:, 0], y_pis[:, 1]
    
    df_pred = df_data.drop(df_features.columns, axis=1)
    df_pred = df_pred.assign(**{'pred': y_pred, 'pred_lower': ylower, 'pred_upper': yupper})
    if isinstance(save_filename, str):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        df_pred.to_csv(Path(output_dir, save_filename), index=False)
    return df_pred

