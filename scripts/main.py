import pandas as pd
from collections import defaultdict
import numpy as np
def selecting_molecules(df):
    '''Selecting all molecules possible function '''
    unique_df = df.drop_duplicates(subset=["Molecule ID"],keep = "last")
    #print(unique_df)
    molec_id = dict(zip(unique_df["Molecule ID"], zip(unique_df["Qmap_position"],unique_df["siteID"])))
    print(molec_id)
    return molec_id

def mk_categories():
    '''Creating categories change this to a function'''
    categories = ["fused_telomere","not_fused_telomere_(normal)","not_fused_no_telomere","fused_no_telomere"]
    bins = defaultdict(set)
    for cats in categories:
        bins[cats]
    return bins

def finding_averages(df,contig,molec_id):

    general_dist_array = []
    contig_gap_dist = defaultdict(list)
    for molecule in molec_id.keys():
        # The try is to stop index error for molecules that do not have contig site
        try:
            mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)
            distance = molec_id[molecule][0] - df[mask].iloc[0,5]
            general_dist_array.append(abs(distance))


            for i in range(0,11):
                if contig + 5 - i == contig:
                    continue
                mask_i = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig + 5 - i) # The 5 here is to capture more molecules extending the range to 31 to 21 if contig 26 else for any other contig +,- 5
                gap_distance = abs(df[mask].iloc[0,5] - df[mask_i].iloc[0,5])
                contig_gap_dist[contig+5-i].append(gap_distance)
        except IndexError:
            continue
        
    # averaging
    general_dist_avg = np.average(general_dist_array)


    contig_gap_dist_avg = {}
    for contigs, avgs in contig_gap_dist.items():
        contig_gap_dist_avg[contigs] = np.average(avgs)

    return general_dist_avg,contig_gap_dist_avg
    


    





def classify_molecules(df,bins,molec_id,contig,total_avg,gap_avg):
    


    '''
    '''
    for molecule in molec_id.keys():
        # check with Dr. Xiao for automating this contig choice if possible
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)
        
        # Checking if molecule does not have contig site 
        try:
            idx = df.index[mask][0]
        except IndexError:
            # using averages
            check = True
            contig_cp = contig + 5
            while check:
                if contig_cp < min(gap_avg.keys()):
                    print(f"This molecule {molecule} contig is beyond 5 ranges bellow the defined initial contig")
                    break
                contig_cp -=1
                mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig_cp)

                try:
                    idx = df.index[mask][0]
                    red_mask_before = df.iloc[idx - 1, 6]
                    red_mask_after  = df.iloc[idx + 1, 6]
                    check = False
                    distance_kb = total_avg + gap_avg[contig_cp]
                    row_pos = df[mask].iloc[0,7]
                    row_distance = molec_id[molecule][1] - row_pos

      
        
                
                    if distance_kb >= 10_000 and (red_mask_before == 1 or red_mask_after == 1):
                        bins["fused_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                        #print(f"Adding molecule {molecule} to fused_telomere ")
                    elif distance_kb >= 10_000 and not (red_mask_before == 1 or red_mask_after == 1):
                        bins["fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                        #print(f"Adding molecule {molecule} to fused_no_telomere ")
                    elif distance_kb < 10_000 and (red_mask_before == 1 or red_mask_after == 1):
                        bins["not_fused_telomere_(normal)"].add(f"{molecule}_({distance_kb},{row_distance})")
                        #print(f"Adding molecule {molecule} to not_fused_telomere_(normal) ")
                    elif distance_kb < 10_000 and not (red_mask_before == 1 or red_mask_after == 1):
                        bins["not_fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                        #print(f"Adding molecule {molecule} to not_fused_no_telomere ")
                    else:
                        print("Something is not working properly this scope should never be acessed")
                except IndexError:
                    continue
            print(f"molecule {molecule} did not have contig 26 using average data to classify it")
            continue

            
        
        # Checking for telomere next to the chosen contig site
        red_mask_before = df.iloc[idx - 1, 6]
        red_mask_after  = df.iloc[idx + 1, 6]
      
        

        # Getting qmap postion
        qmap_pos = df[mask].iloc[0,5] 
        distance_kb = molec_id[molecule][0] - qmap_pos
        # Getting row position
        row_pos = df[mask].iloc[0,7]
        row_distance = molec_id[molecule][1] - row_pos
        

            
        

        if distance_kb >= 10_000 and (red_mask_before == 1 or red_mask_after == 1):
            bins["fused_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
            #print(f"Adding molecule {molecule} to fused_telomere ")
        elif distance_kb >= 10_000 and not (red_mask_before == 1 or red_mask_after == 1):
            bins["fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
            #print(f"Adding molecule {molecule} to fused_no_telomere ")
        elif distance_kb < 10_000 and (red_mask_before == 1 or red_mask_after == 1):
            bins["not_fused_telomere_(normal)"].add(f"{molecule}_({distance_kb},{row_distance})")
            #print(f"Adding molecule {molecule} to not_fused_telomere_(normal) ")
        elif distance_kb < 10_000 and not (red_mask_before == 1 or red_mask_after == 1):
            bins["not_fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
            #print(f"Adding molecule {molecule} to not_fused_no_telomere ")
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
    
    general_dist_avg,contig_gap_dist_avg = finding_averages(df,contig,molec_id) 
    
    bins = classify_molecules(df,empty_bins,molec_id,contig,general_dist_avg,contig_gap_dist_avg)

    
    print(bins)
    total = (len(bins["fused_telomere"]) + len(bins["fused_no_telomere"]) + len(bins["not_fused_no_telomere"]) + len(bins["not_fused_telomere_(normal)"]))
    percent_fusion = round(100 * (len(bins["fused_telomere"]) + len(bins["fused_no_telomere"])) / total,3)
    print("Total percent fusion:",percent_fusion,"%")
    for classification, molecules in bins.items():
        percent = round(100 * len(molecules) / total, 3)
        print(f"{classification}: {percent}%")
    print("Printing AVGS")
    print(general_dist_avg,contig_gap_dist_avg)






    # account for orientation?

    # save it in a txt can change later

    return None


if __name__ == "__main__":
    main()


        



    
