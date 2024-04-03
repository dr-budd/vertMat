##################################################
# RUN UNCERTAINTY QUANTIFICATION FOR REGRESSION WITH MAPIE
# Author: Suk Yee Yong
# Adapted from https://github.com/yongsukyee/uncertain_blackholemass
# Updated: 17 November 2023
# Changes:
#   - remove main, eval_pierrmetric_alpha
#   - remove cfg argument
#   - add seed as argument in optimize_regressor
#   - save cqr model
#   - save filename with alpha in mapie_cqr_alpha and mapie_pred_alpha
#   - add use HalvingRandomSearchCV instead of RandomizedSearchCV
#   - add use X/y train and calib when fitting regressor model
##################################################


from mapie.metrics import regression_coverage_score, regression_mean_width_score
from mapie.regression import MapieQuantileRegressor, MapieRegressor
from mapie.subsample import Subsample
from pathlib import Path
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingRandomSearchCV

import numpy as np
import pandas as pd
import scipy.stats as stats
import time


def optimize_regressor(X_train, y_train, regressor=None, n_iter=100, cv=10, filename='estimator_optim.pkl', seed=20):
    """Optimize regressor. Save output file: estimator_optim.pkl"""
    t1 = time.time()
    print(f"OPTIMIZING REGRESSOR ...")
    
    try:
        estimator = pd.read_pickle(filename)
        print(f"Loaded optimized regressor >> {filename}")
    except:
        param_distributions = {
            'learning_rate': stats.uniform(),
            # 'n_estimators': stats.randint(10, 200),
            'max_depth': stats.randint(2, 30),
            'max_leaf_nodes': stats.randint(2, 50),
        }
        gbr_kwargs = dict(loss='squared_error')
        if regressor == 'quantile':
            gbr_kwargs = dict(loss='quantile', alpha=0.5)
        estimator = GradientBoostingRegressor(random_state=seed, verbose=0, **gbr_kwargs)
        # Search CV
        """
        cv_obj = RandomizedSearchCV(
            estimator,
            param_distributions=param_distributions,
            n_iter=n_iter,
            cv=cv,
            verbose=0,
            random_state=seed,
            n_jobs=-1,
        )
        """
        cv_obj = HalvingRandomSearchCV(
            estimator,
            param_distributions=param_distributions,
            factor=n_iter,
            cv=cv,
            verbose=1,
            random_state=seed,
            n_jobs=-1,
            resource="n_estimators",
            min_resources=10,
            max_resources=200,
        )
        cv_obj.fit(X_train, y_train)
        estimator = cv_obj.best_estimator_
        pd.to_pickle(estimator, filename)
        print(f"\tTime >> {time.time() - t1}")
        print(f"Best estimator >> {estimator}")
    return estimator


def mapie_uq(estimator, estimator_q, X_train, y_train, X_calib, y_calib, X_test, y_test, list_strategies=None, alpha=0.1, filepath='./', save_outfile=True):
    """Compute different MAPIE uncertainty quantification methods for regression. Save output file: mapie_pred.pkl"""
    
    strategies_dict = {
        'naive': {'method': 'naive'},
        'jackknife': {'method': 'base', 'cv': -1},
        'jackknife_plus': {'method': 'plus', 'cv': -1},
        'jackknife_minmax': {'method': 'minmax', 'cv': -1},
        'jackknife_plus_ab': {'method': 'plus', 'cv': Subsample(n_resamplings=50)},
        'cv': {'method': 'base', 'cv': 10},
        'cv_plus': {'method': 'plus', 'cv': 10},
        'cv_minmax': {'method': 'minmax', 'cv': 10},
        'cqr': {'method': 'quantile', 'cv': 'split', 'alpha': alpha},
    }
    
    if isinstance(list_strategies, str): list_strategies = [list_strategies]
    if isinstance(list_strategies, (list, tuple)):
        strategies = {s: strategies_dict[s] for s in list_strategies}
    else:
        strategies = strategies_dict
    
    y_pred, y_pis = {}, {}
    uq_dict = {}
    for strategy, params in strategies.items():
        t1 = time.time()
        print(f"Running strategy >> {strategy}")
        mapiestrategy_file = Path(filepath, f"mapie_{strategy}{f'_alpha{alpha}' if strategy=='cqr' else ''}.pkl")
        if mapiestrategy_file.is_file():
            mapie = pd.read_pickle(mapiestrategy_file)
            print(f"Loaded strategy file >> {mapiestrategy_file}")
        else:
            # For quantile regressor
            if strategy == 'cqr':
                mapie = MapieQuantileRegressor(estimator_q, **params)
                mapie.fit(X_train, y_train, X_calib=X_calib, y_calib=y_calib)
                pd.to_pickle(mapie, mapiestrategy_file)
            # For regressor
            else:
                mapie = MapieRegressor(estimator, n_jobs=-1, **params)
                mapie.fit(np.concatenate((X_train, X_calib)), np.concatenate((y_train, y_calib)))
                pd.to_pickle(mapie, mapiestrategy_file)
        y_pred[strategy], ypis = mapie.predict(X_test, alpha=None if strategy=='cqr' else alpha)
        y_pis[strategy] = ypis.squeeze() # Shape (N, 2 for lower and upper)
        print(f"\tTime >> {time.time() - t1}")
        sorted_idx = np.argsort(y_test)
        ylower, yupper = y_pis[strategy][:, 0], y_pis[strategy][:, 1]
        uq_dict[strategy] = {'target': np.array(y_test)[sorted_idx], 'pred': y_pred[strategy][sorted_idx], 'lower': ylower[sorted_idx], 'upper': yupper[sorted_idx],
                             'pierr_metric': {'PICP': regression_coverage_score(y_test, ylower, yupper), 'MPIW': regression_mean_width_score(ylower, yupper)}}
    
    if save_outfile: pd.to_pickle((y_pred, y_pis, uq_dict), Path(filepath, f"mapie_pred_alpha{alpha}.pkl"))
    return y_pred, y_pis, uq_dict

