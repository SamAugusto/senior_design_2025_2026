import pandas as pd
from collections import defaultdict

def selecting_molecules(df):
    '''Selecting all molecules possible function '''
    unique_df = df.drop_duplicates(subset=["Molecule ID"])
    molec_id = dict(zip(unique_df["Molecule ID"], unique_df["Ref_End"]))
    return molec_id

def mk_categories():
    '''Creating categories change this to a function'''
    categories = ["fused_telomere","not_fused_telomere_(normal)","not_fused_no_telomere","fused_no_telomere"]
    bins = defaultdict(set)
    for cats in categories:
        bins[cats]
    return bins


def classify_molecules(df,bins,molec_id,contig):

    '''
    # 25
    # Molecule classification function waiting on Dr. Xiao to confirm logic and help what to do with molecules that
    # do not have the current select contig

    #### Important note: contig 25 should not be hard setted and should be changed to a input that makes the most sense
    # Because since we are choosing the contig it might differ in every run so the user should be prompted


    ### Changes
    # If it does not have the selected contig add +0 untill 5 and cannot find a starting position at {contig+1...}


    '''
    for molecule in molec_id.keys():
        # check with Dr. Xiao for automating this contig choice if possible
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)
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
    return bins
 
def main():
    # Molecule id with Ref_ED

    # Just reading data this can be a input instead of hard setted 
    # if this is an input = "path to the csv" output1/result csv
    df = pd.read_csv("../data/output1/Result_sgTelo_target_Mol_data 2.csv.gz")

    molec_id = selecting_molecules(df)

    empty_bins =  mk_categories()
    
    contig = int(input("Insert the contig chosen (must be an integer): "))

    bins = classify_molecules(df,empty_bins,molec_id,contig)

    
    print(bins)
    total = (len(bins["fused_telomere"]) + len(bins["fused_no_telomere"]) + len(bins["not_fused_no_telomere"]) + len(bins["not_fused_telomere_(normal)"]))
    percent_fusion = round(100 * (len(bins["fused_telomere"]) + len(bins["fused_no_telomere"])) / total,3)
    print("Total percent fusion:",percent_fusion,"%")
    for classification, molecules in bins.items():
        percent = round(100 * len(molecules) / total, 3)
        print(f"{classification}: {percent}%")





    # % of aligment molecule to contig or reference idk how to do this yet?


    # account for orientation
    # include distance to decide in the molecular bins
    # count each label(row past the distance)
    ## Ex moleculde_id#_(row_count,distance_to_end)

    # save it in a txt can change later

    return None


if __name__ == "__main__":
    main()


        



    
