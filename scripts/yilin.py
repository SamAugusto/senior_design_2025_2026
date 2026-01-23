# create by Yilin Nov 15th, 2021
# updates by Yilin Nov 17th, 2022 (able to run mutilple contigs within one contig txt file, note all contigs need to belong to same sample )
#Pervious verison did:
#1. loop for miutiple contigs, One sample at a time
#2. based on the infromation from the contig txt file to extract molecules
#2.1 if not match will print on the list 'Not_able_extracting_data_'+ Sample_folder_name+'.csv'
#5 if not able to find the contig, will report into 'Not_Exist_Xmap_list'+ Sample_folder_name+'.csv'

#Before Running the scirpts (need to have the following files in the sample folder):
#1. a list of the contig data(Sample&contig info txt file)
#2. Zip assembly file
#3. this scirpts: update_extractingMolecules_MutilContig_V3_SingleSample.py

#Need to change the Code Line -->Run different Sample need to change:
#1. Line 33: Path to the Telomeres_measurement (the folder before the sample folder)
#2. Line 34 change to the correct sample folder name
import pandas as pd
from linecache import getline
import zipfile
import os
import glob
import sys
from datetime import datetime
#timming the time
t_start = datetime.now()
# SETION1: extrat moleule
file_path = 'your/folder/path'

Sample  = sys.argv[1]
#get each sample's contig list and get the contig files for each sample
#Sample = 'siLuc_U2OS_3128992'
Sample_path = os.path.join(file_path, Sample)

#geting into the sample folder path
os.chdir(Sample_path)
contig_extension = 'txt'
contig_file_name = [i for i in glob.glob('*.{}'.format(contig_extension))]
contig_file_name = contig_file_name[0]

#contig_file_name = 'xxxx.txt'
#get the contig txt file path and read the file
contig_file = os.path.join(Sample_path, contig_file_name)
Contiglist_data = pd.read_csv(contig_file, sep='\t',skiprows = 1,names = ["BNX_JobID","Aseembly_JobID","sample","Chrom#",'p or q',"Contig_ID", "TargetLocation", "contig_site ID","orientation"], encoding='latin-1')
Contiglist_data['sample'] = Contiglist_data['sample'].fillna(0)
Contiglist_data = Contiglist_data[Contiglist_data['sample']!=0]
print(Contiglist_data)
#xmapfile fortmat
filename_start = 'exp_refineFinal1_contig'
filename_end = '.xmap'

# get the samplename
samplename = list(set(Contiglist_data['sample'].tolist()))
Sample_folder_name = str(samplename[0])
print(Sample_folder_name)
Assembly_JobID = list(set(Contiglist_data['Aseembly_JobID'].tolist()))
Assembly_JobID = str(int(Assembly_JobID[0]))

#unzip the assembly file and give it a new name
assembly_file = Sample_folder_name+"_"+Assembly_JobID #get the assembly file name
path_to_zip_file_Assembly = os.path.join(Sample_path, assembly_file)
path_to_unzip_file_Assembly = Sample_path
assembly_zip_filename = "_results.zip"
# Check whether the specified folder path exists or not
isExist = os.path.exists(path_to_zip_file_Assembly)
# if not exits
#Create a new directory because it does not exist
os.chdir(path_to_unzip_file_Assembly)
if not isExist:
    for item in os.listdir(path_to_unzip_file_Assembly):
        if item.endswith(assembly_zip_filename):
            file_name = os.path.abspath(item)
            zip_ref = zipfile.ZipFile(file_name)
            output_path = os.path.join(path_to_unzip_file_Assembly, assembly_file)
            zip_ref.extractall(output_path)
            zip_ref.close()

output_path = os.path.join(Sample_path, assembly_file)
# extract the Sample data
Sample_data = Contiglist_data[Contiglist_data['sample'] == Sample_folder_name].reset_index(drop=True)
Sample_contig_list = Sample_data['Contig_ID'].astype(int).tolist()
Sample_contig_list = list(map(str,Sample_contig_list))
Sample_contig_filename = [filename_start + s + filename_end for s in Sample_contig_list]
contig_file_end= "output/contigs/exp_refineFinal1/alignmol/merge/"
filename_path_contig = os.path.join(output_path, contig_file_end)
file_conting_xmap = [filename_path_contig + s for s in Sample_contig_filename]

#create result folder and subfolders for diffrunid
# 1. Check whether the result folder path exists or not
result_foldername = "Result_" + assembly_file
path_to_result_folder = os.path.join(path_to_unzip_file_Assembly, result_foldername)
isExist = os.path.exists(path_to_result_folder)
# 1.1 create a result folder to save output of the sample
if not isExist:
    os.makedirs(path_to_result_folder)

not_isExist_xmap =[]
#getting molecule ID --Sample_data
#criteria
window_size = 100000 # 100kb window
missing_moleculelist = []
print("SETION1: Extracting molecte data")
for xmap_file in file_conting_xmap:
    print(xmap_file)
    filename_xmap =xmap_file.rsplit('/')[-1]
    sample =filename_xmap.rsplit('.')[0]
    contig_name = sample.rsplit('contig')[-1]
    rmap_end = '_r.cmap'
    contig_rmap = sample+rmap_end
    contig_rmap_path = os.path.join(filename_path_contig, contig_rmap)
    qmap_end = '_q.cmap'
    contig_qmap = sample+qmap_end
    contig_qmap_path = os.path.join(filename_path_contig, contig_qmap)
    isExist = os.path.exists(xmap_file)
    if isExist:
        print("Xamp Exist,extracting the molecules")
        contig_qmap = pd.read_csv(contig_qmap_path, sep='\t',usecols=[0,3,4,5], header=None, comment='#', names=['Molecule ID','siteID','LabelChannel','Qmap_position'])
        contig_xmap = pd.read_csv(xmap_file, sep='\t',usecols=[1,2,3,4,5,6,7,13], header=None, comment='#', names=['Molecule ID','Contig ID','Q_start','Q_end','Ref_Start','Ref_End','Ori','Alignment'])
        contig_rmap =pd.read_csv(contig_rmap_path, sep='\t',usecols=[0,1,3,5], header=None, comment='#', names=['Contig_ID','ContigLength','Contig_Site','Contig_Position'])
        #no loop --contig shift
        contig_qmap.rename(columns = {"index": "CMapRowNum"}, inplace=True)
        copy_xmap = contig_xmap.copy()
        copy_xmap.set_index('Molecule ID', inplace=True)
        df_match = copy_xmap.loc[contig_qmap['Molecule ID']].copy()
        df_match.reset_index(drop=False, inplace=True)
        df_match['Qmap_position'] = contig_qmap['Qmap_position']
        df_match['LabelChannel'] = contig_qmap['LabelChannel']
        df_match['siteID'] = contig_qmap['siteID']
        exp = contig_qmap[['Molecule ID','Qmap_position','siteID','LabelChannel']]
        exp =  exp[exp['LabelChannel'] == 2]
        exp['Aligned_label siteID'] = exp['Qmap_position'].groupby(exp["Molecule ID"]).rank(ascending = 1,method = 'dense').astype(int)
        df_match = df_match.merge(exp,on=['Molecule ID','Qmap_position','siteID','LabelChannel'],how='left')
        df_match['Aligned_label siteID'] = df_match['Aligned_label siteID'].fillna(0).astype(int)
        alignment = contig_xmap[['Molecule ID','Alignment']]
        # alignment['Alignment'] = alignment['Alignment'].str.replace('[()]',' ')
        alignment['Alignment'] = alignment['Alignment'].str.replace("["," ")
        alignment['Alignment'] = alignment['Alignment'].str.replace("("," ")
        alignment['Alignment'] = alignment['Alignment'].str.replace(")"," ")
        alignment['Alignment'] = alignment['Alignment'].str.replace("]"," ")
        alignment['Alignment'] = alignment['Alignment'].str.split(" ")
        print(alignment)
        alignment_index = alignment.index
        alignemnt_data = pd.DataFrame(columns =['contig_site','molecule_site'])
        for i in alignment_index:
            melecule_alignement=alignment.loc[i]['Alignment']
            melecule_alignement = [x.strip() for x in melecule_alignement if x.strip() != '']
            molecule_df = pd.DataFrame([sub.split(",") for sub in melecule_alignement],columns =['contig_site','molecule_site'])
            molecule_df['Molecule ID']=alignment.loc[i]['Molecule ID']
            alignemnt_data = pd.concat([alignemnt_data, molecule_df])
        alignemnt_data = alignemnt_data.reset_index(drop=True)
        alignemnt_data['Molecule ID'] = alignemnt_data['Molecule ID'].fillna(0).astype(int)
        alignemnt_data = alignemnt_data.rename(columns={"molecule_site": "Aligned_label siteID"})
        alignemnt_data['Aligned_label siteID'] = alignemnt_data['Aligned_label siteID'].astype(int)
        alignemnt_data['contig_site'] = alignemnt_data['contig_site'].astype(int)
        df_match= df_match.merge(alignemnt_data,on=['Molecule ID','Aligned_label siteID'],how='left')
        df_match['contig_site'] = df_match['contig_site'].fillna(0).astype(int)
        df_match = df_match.rename(columns={"contig_site": "Contig_Site", "Contig ID": "Contig_ID"})
        df_match2= df_match.merge(contig_rmap,on=['Contig_ID','Contig_Site'],how='left')
        df_match2['Contig_Position'] = df_match2['Contig_Position'].fillna(0).astype(int)
        contig_txt_related_data = Contiglist_data[Contiglist_data['Contig_ID'] == int(contig_name)]
        contig_txt_related_data = contig_txt_related_data.reset_index(drop=True)
        print(contig_txt_related_data)
        for i in range(len(contig_txt_related_data)):
            print("index : " , i)
            print(contig_txt_related_data.loc[i]['TargetLocation'])
            target_location = int(contig_txt_related_data.loc[i]['TargetLocation'])
            target_contig_siteid = int(contig_txt_related_data.loc[i]['contig_site ID'])
            ori = str(contig_txt_related_data.loc[i]['orientation'])
            s_e_chrom = str(contig_txt_related_data.loc[i]['p or q'])
            chrom = int(contig_txt_related_data.loc[i]['Chrom#'])
            molecule_list = []
            if ori == "+":
                if s_e_chrom == 'q':
                    windows_E = target_location
                    windows_S = windows_E - window_size
                    fit1 = df_match2[df_match2['Contig_Position']>=windows_S]
                    fit2 = fit1[fit1['Contig_Position']<=windows_E]
                    match_number = len(fit2)
                    match_to_target = fit1[fit1['Contig_Site']==target_contig_siteid]
                    if match_number>6:
                        if not match_to_target.empty:
                            molecule_name = match_to_target['Molecule ID'].astype(int).tolist()
                            molecule_list=molecule_list+molecule_name
            if ori == "+":
                if s_e_chrom == 'p':
                    windows_S = target_location
                    windows_E = windows_S + window_size
                    fit1 = df_match2[df_match2['Contig_Position']<=windows_E]
                    fit2 = fit1[fit1['Contig_Position']>=windows_S]
                    match_number = len(fit2)
                    match_to_target = fit1[fit1['Contig_Site']==target_contig_siteid]
                    if match_number>6:
                        if not match_to_target.empty:
                            molecule_name = match_to_target['Molecule ID'].astype(int).tolist()
                            molecule_list=molecule_list+molecule_name
            if ori == "-":
                if s_e_chrom == 'q':
                    windows_S = target_location
                    windows_E = windows_S + window_size
                    fit1 = df_match2[df_match2['Contig_Position']<=windows_E]
                    fit2 = fit1[fit1['Contig_Position']>=windows_S]
                    match_number = len(fit2)
                    match_to_target = fit1[fit1['Contig_Site']==target_contig_siteid]
                    if match_number>6:
                        if not match_to_target.empty:
                            molecule_name = match_to_target['Molecule ID'].astype(int).tolist()
                            molecule_list=molecule_list+molecule_name
            if ori == "-":
                if s_e_chrom == 'p':
                    windows_E = target_location
                    windows_S = windows_E - window_size
                    fit1 = df_match2[df_match2['Contig_Position']>=windows_S]
                    fit2 = fit1[fit1['Contig_Position']<=windows_E]
                    match_number = len(fit2)
                    match_to_target = fit1[fit1['Contig_Site']==target_contig_siteid]
                    if match_number>6:
                        if not match_to_target.empty:
                            molecule_name = match_to_target['Molecule ID'].astype(int).tolist()
                            molecule_list=molecule_list+molecule_name
            Molecule_df= df_match2[df_match2['Molecule ID'].isin(molecule_list)]
            Molecule_df = Molecule_df[['Molecule ID','Contig_ID','Ref_Start','Ref_End','Ori','Contig_Site','Contig_Position']]
            if Molecule_df.empty:
                print("molecule not able to find will save on file name as: Not_able_extracting_data_*.csv ")
                cannofind_molecule_id_list = str(contig_name)+'_'+str(target_location)
                missing_moleculelist.append(cannofind_molecule_id_list)
            else:
                result_filename = f'Result_{Sample_folder_name}_target_Mol_data.csv'
                path_to_final_result = os.path.join(path_to_result_folder, result_filename)
                Molecule_df.to_csv(path_to_final_result,index=False)
    if not isExist:
        not_isExist_xmap.append(xmap_file)

not_isExist_xmap_data = pd.DataFrame(not_isExist_xmap,columns = ['mssing_xmap_list'])
not_isExist_xmap_name = 'Not_Exist_Xmap_list'+ Sample_folder_name+'.csv'
path_to_result_not_isExistXmaplist = os.path.join(path_to_result_folder, not_isExist_xmap_name)
if not not_isExist_xmap_data.empty:
    not_isExist_xmap_data.to_csv(path_to_result_not_isExistXmaplist,index=False)

not_isExist_molecule_list = pd.DataFrame(missing_moleculelist,columns = ["mssing_contig_data(contig_name+'_'+TargetLocation)"])
not_isExist_molecule_list_name = 'Not_able_extracting_data_'+ Sample_folder_name+'.csv'
path_to_not_isExist_molecule_list_name = os.path.join(path_to_result_folder, not_isExist_molecule_list_name)
if not not_isExist_molecule_list.empty:
    not_isExist_molecule_list.to_csv(path_to_not_isExist_molecule_list_name,index=False)

t_end = datetime.now()
print("starting Running:",t_start)
print("finishing Running:",t_end)
#place holder
