# Prediction Intervals for Species Lifespans and Age at Maturity

Generating prediction intervals for species lifespans and age at maturity.
The provided predictions have been trained on cytosine-phosphate-guanine (CpG) density of group-specific (fish, mammals, and reptiles) age at maturity model using 10-fold cross-validation ensemble model of elastic net regression. Several prediction interval methods to quantify the uncertainties of the predictions of the ensemble model, including naive, jackknife, cross-validation, and conformal quantile regression, are explored.

## Installation
```
git clone https://github.com/dr-budd/vertMat.git
cd life_pi
pip install -r requirements.txt
```

## Quick Start
### For training prediction interval methods
1. Add dataset: Create folder and add data containing predictions and features to `data/` (see also [input data format](#input-data-format) in Usage)
2. Setup configs: Copy and rename template from `config/config_template.yaml`, and edit configs
3. Generate prediction intervals: `$ python -m life_pi.main --run=traineval -c config/<CONFIG_FILE>.yaml`
### For inference with trained prediction interval model
1. Setup configs: Copy and rename template from `config/config_pred_template.yaml`, and edit configs
2. Generate prediction intervals: `$ python -m life_pi.main --run=eval -c config/<CONFIG_PRED_FILE>.yaml`

## Usage
### Input data format
See example data in folder (provided upon request): `data/`
- `_ensemble_data_.csv`: Input data for training with required columns `organism_name, mean_age_mat, log_age_mat, model_number, initial_split, log_predicted_age_mat, FEATURES`
- `_glmnetcv_parameters_.csv`: (Optional) Elastic net parameters from glmnet with required columns `model_number, alpha_min, lambda_1se`
- `Table S[db_age_mat].csv`: Reported age at maturity of known species with required columns `organism_name, age_mat`
### Train prediction interval method using different regressor
In `main.py`, `generate_pi.train_pi` function, different regressor can be used
- `estimator`: scikit-learn regressor to train naive, jackknife, and cross-validation methods. Default using elastic net.
- `estimator_q`: scikit-learn regressor with quantile to train conformal quantile regressor
### Output files from main scripts
- From `--run=traineval`:
  - `datascalerx_<scale_X>.pkl`, `datascalery_<scale_y>.pkl`: dict of inverse transform values for X, y data if the data are scaled {'minus_numerator', 'denominator} using the function `generate_pi.invt_scale_data`
  - `datasplitidx.pkl`: dict of data split indices {'train': train index, 'valid': validation index, 'test': test index}
  - `estimator_optim.pkl`: Optimized GradientBoostingRegressor model for normal regression
  - `estimatorq_optim.pkl`: Optimized GradientBoostingRegressor model for quantile regression
  - `mapie_<STRATEGY>.pkl`: Fitted MapieRegressor for specified strategy
  - `mapie_pred_alpha<alpha>.pkl`: y_pred, y_pis, sorted {'target', 'pred', 'lower', 'upper'}
- From `--run=eval`:
  - `mapie_<MODEL_NAME>_meanpred_alpha<alpha>.csv`: Mean predictions and prediction intervals with trained model of all data partitions
### Notebook for analysis
Folder: `notebooks/<DATASET_NAME>_pi.ipynb`
- To plot histogram ridgeline plot:
  ```python
  # For all dataset
  strategy = 'cv' # Select strategy
  df_data_test = generate_pi.get_dftest(data_dir=data_dir, filename=cfg['filename'], file_nsplit=0, output_dir=output_dir, y_labelname=cfg['y_labelname'], strategy=strategy)
  plot_utils.histridgeline_lifespanpi(df_data_test, y_labelname=cfg['y_labelname'], df_datareport=df_datareport, use_baggedpred=True, strategy=strategy, ylabel='Age at Maturity (years)', save_plotname=False)
  
  # For selected data partition
  strategy = 'cv'
  file_nsplit = 1
  df_datapart_test = generate_pi.get_dftest(data_dir=data_dir, filename=cfg['filename'], file_nsplit=file_nsplit, output_dir=output_dir, y_labelname=cfg['y_labelname'], strategy=strategy)
  plot_utils.histridgeline_lifespanpi(df_datapart_test, y_labelname=cfg['y_labelname'], df_datareport=df_datareport, use_baggedpred=False, strategy=strategy, ylabel='Age at Maturity (years)', save_plotname=False)
  ```

<!-- 
## TODO
- add multiprocessing all partitions in traineval or in for loop and edit file_nsplit in bash script for Petrichor
 -->

## Author
Suk Yee Yong (sukyee.yong@csiro.au)

## About
Original code from https://bitbucket.csiro.au/scm/sccp-sukyee/life_pi.git. This project is being developed as part of Scientific Computing Collaboration Project 2023H2 for [ERRFP 1302: How long (ish) do fish live? Estimating prediction error for regularised regression](https://confluence.csiro.au/pages/viewpage.action?spaceKey=SCinternal&title=ERRFP-1302).

## Acknowledgements
The script `life_pi/lib/run_uqregresor.py` has been adapted from [uncertain_blackholemass](https://github.com/yongsukyee/uncertain_blackholemass) repository.  
Data has been provided by Alyssa Budd.  

## License
[CSIRO Open Source Software Licence Agreement (variation of the BSD / MIT License)](LICENSE.txt)
