import os
import numpy as np
import pandas as pd

base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "data")
cases_file = os.path.join(data_dir, "combined_data.csv")

# Load the CSV file
df = pd.read_csv(cases_file)

# Filter rows where municipio is "rio de janeiro"
# and cases is 0
filtered = df[
    (df["municipio"].str.lower() == "rio de janeiro") &
    (df["cases"] == 0)
]

# Print the year and week columns
print(filtered[["year", "week"]])

# Optional: print total count
print(f"\nTotal rows: {len(filtered)}")