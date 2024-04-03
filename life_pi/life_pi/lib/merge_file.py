# Script to merge all partition files to one
# Usage: python lib/merge_file.py

data_dir = 'data/fish'
prefix_filename = '10.00_partition' # Data filename prefix
save_filename = '10.00_ensemble_data_fish'

from pathlib import Path
import pandas as pd

df = pd.concat([pd.read_csv(f) for f in sorted(Path(data_dir).glob(prefix_filename+'*'))])
df.to_csv(Path(data_dir, f"{save_filename}.csv"), index=False)
