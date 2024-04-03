##################################################
# PREDICTION INTERVAL FOR SPECIES LIFESPANS
# Author: Suk Yee Yong
##################################################

from pathlib import Path
import argparse
import pandas as pd

from life_pi.lib.get_config import get_config
import life_pi.lib.generate_pi as generate_pi


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prediction interval for species lifespans')
    parser.add_argument('-c', '--cfgfile', default='config/config.yaml', required=True, help='Configuration file. For example, see `config/config_template.yaml`')
    parser.add_argument('-r', '--run', choices=['traineval', 'eval'], default='eval', help='Run train and/or evaluation')
    args = parser.parse_args()
    cfg = get_config(args.cfgfile)
    
    # Load data
    df_data = pd.read_csv(Path(cfg['data_dir'], f"{cfg['filename']}.csv"))
    
    # Run train and evaluation
    if args.run == 'traineval':
        # For all data partition
        if cfg['file_nsplit'] == 0:
            for key, group in df_data.groupby('model_number'):
                print(f"PARTITION # {key}")
                df_datapart = group.dropna(axis=1, how='all').reset_index(drop=True) # Drop missing values
                outputpart_dir = Path(cfg['output_dir'], str(key)) # Subdirectory for each partition
                generate_pi.train_pi(df_datapart, 'mean_'+cfg['y_labelname'], feature_prefix=cfg['feature_prefix'],
                    scale_X=cfg['scale_X'], scale_y=cfg['scale_y'], log_y=cfg['log_y'], test_size=cfg['test_size'], seed=cfg['seed'],
                    estimator=generate_pi.estimator_elasticnet(cfg['data_dir'], cfg['estimator_filename'], file_nsplit=key, l1ratio_colname='alpha_min', alpha_colname='lambda_1se', seed=1),
                    estimator_q=None,
                    list_strategies=cfg['list_strategies'], alpha=cfg['alpha'], output_dir=outputpart_dir)
        # For one data partition
        else:
            df_datapart = df_data.groupby('model_number').get_group(cfg['file_nsplit'])\
                                    .dropna(axis=1, how='all').reset_ilog_yndex(drop=True) 
            outputpart_dir = Path(cfg['output_dir'], str(cfg['file_nsplit']))
            generate_pi.train_pi(df_datapart, 'mean_'+cfg['y_labelname'], feature_prefix=cfg['feature_prefix'],
                scale_X=cfg['scale_X'], scale_y=cfg['scale_y'], log_y=cfg['log_y'], test_size=cfg['test_size'], seed=cfg['seed'],
                estimator=generate_pi.estimator_elasticnet(cfg['data_dir'], cfg['estimator_filename'], file_nsplit=cfg['file_nsplit'], l1ratio_colname='alpha_min', alpha_colname='lambda_1se', seed=1),
                estimator_q=None,
                list_strategies=cfg['list_strategies'], alpha=cfg['alpha'], output_dir=outputpart_dir)
    
    # Run evaluation
    elif args.run == 'eval':
        df_preds = []
        # Run for all data partition
        for key, group in df_data.groupby('model_number'):
            print(f"PARTITION # {key}")
            df_datapart = group.dropna(axis=1, how='all').reset_index(drop=True) # Drop missing values
            modelpart_path = Path(cfg['model_dir'], str(key), cfg['model_filename'])
            if 'alpha' in cfg['model_filename']:
                cfg['alpha'] = cfg['model_filename'].rsplit('alpha')[-1].rsplit('.pkl')[0]
            df_preds.append(generate_pi.pred_testx(df_datapart, modelpath=modelpart_path, feature_prefix=cfg['feature_prefix'],
                                log_y=cfg['log_y'], alpha=cfg['alpha'],
                                output_dir=cfg['output_dir'], save_filename=False))
        df_preds = generate_pi.merge_dfdatapred(df_preds).drop('model_number', axis=1)
        Path(cfg['output_dir']).mkdir(parents=True, exist_ok=True)
        df_preds.to_csv(Path(cfg['output_dir'], modelpart_path.stem.rsplit('_alpha')[0]+f"_medianpred_alpha{cfg['alpha']}.csv"), index=False)

