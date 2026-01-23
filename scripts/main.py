import pandas as pd
from collections import defaultdict
# Remember to add functions instead of leaving the code like this
if __name__ == "__main__":
    # Molecule id with Ref_ED

    # Just reading data this can be a input instead of hard setted 
    df = pd.read_csv("../data/output1/Result_sgTelo_target_Mol_data 2.csv.gz")


    # Selecting all molecules possible function, very likely changing this to a function
    unique_df = df.drop_duplicates(subset=["Molecule ID"])
    molec_id = dict(zip(unique_df["Molecule ID"], unique_df["Ref_End"]))


    # Creating categories change this to a function    
    categories = ["fused_telomere","not_fused_telomere_(normal)","not_fused_no_telomere","fused_no_telomere"]
    bins = defaultdict(set)
    for cats in categories:
        bins[cats]
    
    # 26
    # Molecule classification function waiting on Dr. Xiao to confirm logic and help what to do with molecules that
    # do not have the current select contig


    #### Important note: contig 26 should not be hard setted and should be changed to a input that makes the most sense
    # Because since we are choosing the contig it might differ in every run so the user should be prompted
    for molecule in molec_id.keys():
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == 26)
        try:
            idx = df.index[mask][0]
        except IndexError:
            print(f"Molecule {molecule} does not have contig site 26")
            continue
        red_mask_before = df.iloc[idx - 1, 6]
        red_mask_after  = df.iloc[idx + 1, 6]
        print("Red_mask_After: ",red_mask_after)
        print("Red_mask_before: ",red_mask_before)
        
        contig_site = df[mask].iloc[0,-1]
        distance_kb = molec_id[molecule] - contig_site 
        print("distance: ",distance_kb)
        if distance_kb >= 100_000 and (red_mask_before == 1 or red_mask_after == 1):
            bins["fused_telomere"].add(molecule)
            print(f"Adding molecule {molecule} to fused_telomere ")
        elif distance_kb >= 100_000 and not (red_mask_before == 1 or red_mask_after == 1):
            bins["fused_no_telomere"].add(molecule)
            print(f"Adding molecule {molecule} to fused_no_telomere ")
        elif distance_kb < 100_000 and (red_mask_before == 1 or red_mask_after == 1):
            bins["not_fused_telomere_(normal)"].add(molecule)
            print(f"Adding molecule {molecule} to not_fused_telomere_(normal) ")
        elif distance_kb < 100_000 and not (red_mask_before == 1 or red_mask_after == 1):
            bins["not_fused_no_telomere"].add(molecule)
            print(f"Adding molecule {molecule} to not_fused_no_telomere ")
        else:
            print("Something is not working properly this scope should never be acessed")
    print(bins)


        



    
