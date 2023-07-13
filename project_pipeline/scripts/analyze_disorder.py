import pandas as pd
import main

df = pd.read_csv('.data/input/rmsds.tsv', sep='\t').astype('object')

df_disorder_1 = main.calculate_disorder(df)