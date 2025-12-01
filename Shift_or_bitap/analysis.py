import pandas as pd
import glob

# Aggregate all pattern length results
all_results = []
for csv_file in glob.glob("results_exact/*/pattern_length_results.csv"):
    df = pd.read_csv(csv_file)
    all_results.append(df)

combined_df = pd.concat(all_results)
print(combined_df.groupby('pattern_length')['time_seconds_mean'].mean())
