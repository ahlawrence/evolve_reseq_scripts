import pandas as pd
import sys

# Ensure correct usage
if len(sys.argv) != 2:
    print("Usage: python script.py <input_file>")
    sys.exit(1)

# Load input file
input_file = sys.argv[1]
df = pd.read_csv(input_file, sep="\t", dtype=str, low_memory=False)
df.columns = df.columns.str.strip()  # strip any spaces in column names

# Extract just genotype (before any colon)
for col in ["SOD", "sample_2", "sample_99", "pool_F1"]:
    df[col] = df[col].str.split(":").str[0]

# Define unwanted genotypes
excluded_genotypes = {"./.", "0/1", "1/2", "2/3"}

# Apply filtering
filtered_df = df[
    (df["SOD"] != df["sample_2"]) &  # SOD and REM_2023 (sample_2) must differ
    (~df["SOD"].isin(excluded_genotypes)) &
    (~df["sample_2"].isin(excluded_genotypes)) &
    (df["sample_99"] == "0/1") & #G0 heterozygous
    (df["pool_F1"] == "0/1") #F1 heterozygous
]

# Get sites where SOD is 0/0 and REM_2023 is 1/1
sod_ref_sites = filtered_df[
    (filtered_df["SOD"] == "0/0") & (filtered_df["sample_2"] == "1/1")
]

# Get sites where REM_2023 is 0/0 and SOD is 1/1
rem_ref_sites = filtered_df[
    (filtered_df["sample_2"] == "0/0") & (filtered_df["SOD"] == "1/1")
]

# Save results
sod_ref_sites.to_csv("SOD_ref.txt", sep="\t", index=False)
rem_ref_sites.to_csv("REM_ref.txt", sep="\t", index=False)

print(f"SOD_ref.txt: {len(sod_ref_sites)} sites")
print(f"REM_ref.txt: {len(rem_ref_sites)} sites")
