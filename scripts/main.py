import pandas as pd
from collections import defaultdict
import numpy as np

def selecting_molecules(df):
    '''Selecting all molecules possible function '''
    unique_df_min = df.drop_duplicates(subset=["Molecule ID"],keep = "first")
    unique_df_max = df.drop_duplicates(subset=["Molecule ID"],keep = "last")

    #print(unique_df)
    molec_id = dict(zip(unique_df_min["Molecule ID"], zip(unique_df_min["Qmap_position"],unique_df_min["siteID"],unique_df_max["Qmap_position"],unique_df_max["siteID"],unique_df_min["Ori"])))
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
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)
        # The try is to stop index error for molecules that do not have contig site
        try:
            idx = df.index[mask][0]
            current_pos = df.iloc[idx, 5]
            red_mask_before = df.iloc[idx - 1, 6]
            red_mask_after  = df.iloc[idx + 1, 6]
            if (red_mask_before == 1 or red_mask_after == 1):
                    if red_mask_before == 1:
                        neighbor_pos = df.iloc[idx - 1, 5]
                        distance = current_pos - neighbor_pos
                        general_dist_array.append(distance)

                    if red_mask_after == 1:
                        neighbor_pos = df.iloc[idx + 1, 5]
                        distance = current_pos - neighbor_pos
                        general_dist_array.append(distance)
            else:
                pass


            for i in range(0,11):
                if contig + 5 - i == contig:
                    continue
                mask_i = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig + 5 - i) # The 5 here is to capture more molecules extending the range to 31 to 21 if contig 26 else for any other contig +,- 5
                if not df[mask_i].empty:
                    gap_distance = abs(df[mask].iloc[0,5] - df[mask_i].iloc[0,5])
                    contig_gap_dist[contig+5-i].append(gap_distance)
                else:
                    continue
        except IndexError:
            continue
        
    # averaging
    general_dist_avg = np.average(general_dist_array)
    print("priting label distance avg: ", general_dist_avg)

    contig_gap_dist_avg = {}
    for contigs, avgs in contig_gap_dist.items():
        contig_gap_dist_avg[contigs] = np.average(avgs)

    return general_dist_avg,contig_gap_dist_avg
    


    





def classify_molecules(df,bins,molec_id,contig,label_avg,gap_avg,chrom_orientation):
    


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
            contig_cp_lst = [contig-1,contig+1,contig-2,contig+2,contig-3,contig+3,contig-4,contig+4,contig-5,contig+5]
            contig_cp_idx = -1
            while check:
                if contig_cp_idx >= len(contig_cp_lst):
                    print(f"This molecule {molecule} contig is beyond 5 ranges bellow the defined initial contig")
                    break
                contig_cp_idx +=1
                contig_cp = contig_cp_lst[contig_cp_idx]
                mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig_cp)

                try:
                    idx = df.index[mask][0]
                    red_mask_before = df.iloc[idx - 1, 6]
                    red_mask_after  = df.iloc[idx + 1, 6]
                    check = False
                    print(f"molecule {molecule} did not have contig 26 using average data to classify it")
      
        
                
                    # classification if label found
                    if (red_mask_before == 1 or red_mask_after == 1):
                        if red_mask_before == 1:
                            qmap_pos = df.iloc[idx - 1, 5]
                            distance_kb = abs(molec_id[molecule][0] - qmap_pos) 
                            row_pos = df.iloc[idx - 1, 7]
                            row_distance = abs(molec_id[molecule][1] - row_pos)
                        else:
                            qmap_pos = df.iloc[idx + 1, 5]
                            distance_kb = abs(molec_id[molecule][0] - qmap_pos) 
                            row_pos = df.iloc[idx + 1, 7]
                            row_distance = abs(molec_id[molecule][1] - row_pos)
                        
                        if distance_kb >=10_000:
                            bins["fused_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                            #print(f"Adding molecule {molecule} to fused_telomere ")
                        else:
                            bins["not_fused_telomere_(normal)"].add(f"{molecule}_({distance_kb},{row_distance})")
                            #print(f"Adding molecule {molecule} to not_fused_telomere_(normal) ")
                           
                    # Classification if label not found
                    else:
                        qmap_pos = df[mask].iloc[0,5] + label_avg + gap_avg[contig_cp]
                        row_pos = df[mask].iloc[0,7]
                        
                        if "p" in chrom_orientation and molec_id[molecule][-1] == "+":
                            distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                            row_distance = abs(molec_id[molecule][1]- row_pos)
                        elif "p" in chrom_orientation and molec_id[molecule][-1] == "-":
                            distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                            row_distance = abs(molec_id[molecule][3] - row_pos)
                        elif "q" in chrom_orientation and molec_id[molecule][-1] == "+":
                            distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                            row_distance = abs(molec_id[molecule][3] - row_pos)
                        elif "q" in chrom_orientation and molec_id[molecule][-1] == "-":
                            distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                            row_distance = abs(molec_id[molecule][1]- row_pos)
                        else:
                            raise Exception("This should never be reached something went wrong")

                        if distance_kb >=10_000:
                            bins["fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                            #print(f"Adding molecule {molecule} to fused_no_telomere ")
                        else:
                            bins["not_fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                            #print(f"Adding molecule {molecule} to not_fused_no_telomere ")

        
                except IndexError:
                    continue
            continue

            
        
        # Checking for telomere next to the chosen contig site
        red_mask_before = df.iloc[idx - 1, 6]
        red_mask_after  = df.iloc[idx + 1, 6]
      
        

        # Getting qmap postion
        #qmap_pos = df[mask].iloc[0,5] 
        #distance_kb = abs(molec_id[molecule][0] - qmap_pos)
        # Getting row position
        #row_pos = df[mask].iloc[0,7]
        #row_distance = abs(molec_id[molecule][1] - row_pos)
        

            
        
        # classification if label found
# molec_id = dict(zip(unique_df_min["Molecule ID"], zip(unique_df_min["Qmap_position"],unique_df_min["siteID"],unique_df_max["Qmap_position"],unique_df_max["siteID"],unique_df_min["Ori"])))
        if (red_mask_before == 1 or red_mask_after == 1):
            if red_mask_before == 1:
                qmap_pos = df.iloc[idx - 1, 5]
                row_pos = df.iloc[idx - 1, 7]
                if "p" in chrom_orientation and molec_id[molecule][-1] == "+":
                    distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                    row_distance = abs(molec_id[molecule][1]- row_pos)
                elif "p" in chrom_orientation and molec_id[molecule][-1] == "-":
                    distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                    row_distance = abs(molec_id[molecule][3] - row_pos)
                elif "q" in chrom_orientation and molec_id[molecule][-1] == "+":
                    distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                    row_distance = abs(molec_id[molecule][3] - row_pos)
                elif "q" in chrom_orientation and molec_id[molecule][-1] == "-":
                    distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                    row_distance = abs(molec_id[molecule][1]- row_pos)
                else:
                    raise Exception("This should never be reached something went wrong")


                
 

            else:
                qmap_pos = df.iloc[idx + 1, 5]
                row_pos = df.iloc[idx + 1, 7]
                if "p" in chrom_orientation and molec_id[molecule][-1] == "+":
                    distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                    row_distance = abs(molec_id[molecule][1]- row_pos)
                elif "p" in chrom_orientation and molec_id[molecule][-1] == "-":
                    distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                    row_distance = abs(molec_id[molecule][3] - row_pos)
                elif "q" in chrom_orientation and molec_id[molecule][-1] == "+":
                    distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                    row_distance = abs(molec_id[molecule][3] - row_pos)
                elif "q" in chrom_orientation and molec_id[molecule][-1] == "-":
                    distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                    row_distance = abs(molec_id[molecule][1]- row_pos)
                else:
                    raise Exception("This should never be reached something went wrong")


            
            if distance_kb >=10_000:
                bins["fused_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                #print(f"Adding molecule {molecule} to fused_telomere ")
                
            else:
                bins["not_fused_telomere_(normal)"].add(f"{molecule}_({distance_kb},{row_distance})")
                #print(f"Adding molecule {molecule} to not_fused_telomere_(normal) ")
               
        # Classification if label not found
        else:

            qmap_pos = df[mask].iloc[0,5] + label_avg
            row_pos = df[mask].iloc[0,7]

            if "p" in chrom_orientation and molec_id[molecule][-1] == "+":
                distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                row_distance = abs(molec_id[molecule][1]- row_pos)
            elif "p" in chrom_orientation and molec_id[molecule][-1] == "-":
                distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                row_distance = abs(molec_id[molecule][3] - row_pos)
            elif "q" in chrom_orientation and molec_id[molecule][-1] == "+":
                distance_kb = abs(molec_id[molecule][2] - qmap_pos)
                row_distance = abs(molec_id[molecule][3] - row_pos)
            elif "q" in chrom_orientation and molec_id[molecule][-1] == "-":
                distance_kb = abs(molec_id[molecule][0] - qmap_pos)
                row_distance = abs(molec_id[molecule][1]- row_pos)
            else:
                raise Exception("This should never be reached something went wrong")

            if distance_kb >=10_000:
                bins["fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                #print(f"Adding molecule {molecule} to fused_no_telomere ")
            else:
                bins["not_fused_no_telomere"].add(f"{molecule}_({distance_kb},{row_distance})")
                #print(f"Adding molecule {molecule} to not_fused_no_telomere ")

        
    return bins
 
def main():
    # Molecule id with qmappos

    # Just reading data this can be a input instead of hard setted 
    # if this is an input = "path to the csv" output1/result csv
    #df = pd.read_csv("../data/output1/Result_sgTelo_target_Mol_data 2.csv.gz")

    chromossome = input("insert chromossome: ")
    orientation = input("insert orientation: ")
    print(chromossome+orientation)

    df = pd.read_excel("../data/output1/4 types of contigs 8 samples.xlsx",sheet_name = chromossome+orientation)
    molec_id = selecting_molecules(df)

    empty_bins =  mk_categories()
    
    contig = int(input("Insert the contig chosen (must be an integer): "))
    
    general_dist_avg,contig_gap_dist_avg = finding_averages(df,contig,molec_id) 
    
    bins = classify_molecules(df,empty_bins,molec_id,contig,general_dist_avg,contig_gap_dist_avg,orientation+chromossome)

    
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


        



    

    



 
