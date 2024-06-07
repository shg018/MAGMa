from django.db import models
from django.conf import settings
from django.core.files.storage import FileSystemStorage
import pandas as pd
import os
from datetime import datetime

#####################Shagun script dependencies#################################
import numpy as np
# import networkx as nx
#import pandas as pd
#import os
import statistics
from scipy import stats
from scipy.stats import pearsonr
import statsmodels
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from itertools import permutations
import re
import math
import subprocess

class OverwriteStorage(FileSystemStorage):
    '''
    Overwrite uploaded files with same file name.
    '''
    def get_available_name(self, name, max_length=None):
        if self.exists(name):
            os.remove(os.path.join(self.location, name))
        return super(OverwriteStorage, self).get_available_name(name, max_length)
# def path_file():
#     def wrapper(user, filename):
#         file_upload_dir = 'documents/tmtipms_files'#os.path.join(settings.MEDIA_ROOT, 'file_upload')
#         if os.path.exists(file_upload_dir):
#             import shutil
#             shutil.rmtree(file_upload_dir)
#         return os.path.join(file_upload_dir, filename)
#     return wrapper

# class Job(models.Model):
#     job_id = models.CharField(max_length=100, unique=True)
#     created_at = models.DateTimeField(auto_now_add=True)

class Document(models.Model):
    docfile = models.FileField(upload_to='documents/tmtipms_files/', max_length=1000, storage=OverwriteStorage())

# class ProcessedFile(models.Model):
#     docfile = models.FileField(upload_to='Results/%Y/%m/%d/')

# class Document(models.Model):
#     #_datetime = datetime.now()
#     #datetime_str = _datetime.strftime("%Y-%m-%d-%H-%M-%S")
#     #docfile = models.FileField(upload_to='documents/datetime_str')
#     docfile = models.FileField(upload_to='documents/tmtipms_files/', storage=OverwriteStorage())
    
#     #docfile = models.FileField(upload_to=path_file())



###########Custom functions of Yugandhar Kumar########################
# def validate_bait(bait_name,search_type, file_path):
#     file_df=pd.read_csv(file_path, sep='\t', header=0)
#     bait_validation=False
#     if (search_type=="sequest"):
#         if bait_name in file_df['Master Protein Accessions'].values:
#             bait_validation=True
#     elif (search_type=="comet"):
#         if bait_name in file_df['protein'].values:
#             bait_validation=True
#     else:
#         if bait_name in file_df['Master Protein Accessions'].values:
#             bait_validation=True
# 	return bait_validation

def validate_bait(bait_name,search_type,file_path):
    file_df=pd.read_csv(file_path, sep='\t', header=0)
    bait_validation=False
    if (search_type=="sequest"):
        if bait_name in file_df["Master Protein Accessions"].values:
            bait_validation=True
    elif (search_type=="comet"):
        file_df["check_substring"] = file_df["protein"].str.contains(bait_name)
        if file_df["check_substring"].sum()>0:
            bait_validation=True
    elif (search_type=="diann"):
        file_df["check_substring"] = file_df["Protein"].str.contains(bait_name)
        if file_df["check_substring"].sum()>0:
            bait_validation=True
    else:
        if bait_name in file_df["Master Protein Accessions"].values:
            bait_validation=True
    return bait_validation

############Custom functions of Shagun Gupta#############################
def sort1(lst):
    lst = [str(i) for i in lst]
    lst.sort()
    #lst = [int(i) if i.isdigit() else i for i in lst ]
    return lst

## IP_DIA version (specifically for controls)
def YuLab_fractionxchannel_IP_DIA(job_id, prot_init_df, raw_input_df, annotation_df, bait_name, comparison, comp_type, row_norm_user,column_norm_user, bait_norm, num_channels, norm_prots, pval_calc, uniprotGene, imputation, fusion_aa_seq='', iso_filter=100):
    # Making sure dependence is not on number channels but on specific ones with labels
    channel_list = list(set(annotation_df.loc[annotation_df["Label"]!="none"]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)
    
    # Will run row normalization by default
    prot_init_df = preprocess_byrownorm(prot_init_df,True,channel_list)#row_norm_user

    if column_norm_user==True:
        if norm_prots == '':
            print("No protein list based normalization")
            prot_init_df = normalization_median_intensity(prot_init_df,channel_list)
            print("Done median ratio based normalization")
        else:
            uniprots_to_norm_by = norm_prots.split(",")
            print("Yes protein list based normalization: ")#, uniprots_to_norm_by)
            prot_init_df = norm_by_prots_byTMM(prot_init_df,channel_list, uniprots_to_norm_by)
            print("Done median ratio based normalization using the proteins listed")

    
    prot_df = expansion_to_fractionXchannel_DIA(prot_init_df,annotation_df, comparison, comp_type)

    # Bait normalization
    if bait_norm==True:
        print("Normalized by bait")
        prot_df = norm_by_bait_DIA(bait_name, prot_df)
    else:
        print("No bait based normalization")

    final_prot_df = []
    if pval_calc == "LinearModel":
        print("For IP-DIA 2 sample ttest")
        final_prot_df = calc_FC_and_pval_IP_DIA(job_id, prot_df, [], comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation, "ip")
    elif (pval_calc == "Limma"):
        print("For IP-DIA Limma")
        final_prot_df = calc_FC_and_pval_limma_IP_DIA(job_id, prot_df, [], comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation, "ip")
    else:
        print("Error: Magma Type not chosen")
    
    print("Done calculating FC and p-value")
    return final_prot_df

def add_Precursor_info(final_prot_df,raw_input_df,annotation_df,comparison,comp_type):
    prot_list = list(set(final_prot_df["Protein"].dropna().values.tolist()))
    wt=""
    if comp_type==1:
        wt = comparison.split("-")[0]
    else:
        wt = comparison.split("-")[1]
    all_conditions = list(set(annotation_df["Label"].values.tolist()))
    # for condition in all_conditions:
    #     final_prot_df["#Precursors-IN-"+condition] = 0
    precursor_info = []
    for prot in prot_list:
        prec_dict = {"Protein":prot}
        for condition in all_conditions:
            files_per_condition = list(set(annotation_df.loc[annotation_df["Label"]==condition]["Channel"].values.tolist()))
            # final_prot_df.loc[final_prot_df["Protein"]==prot,"#Precursors-IN-"+condition] = raw_input_df[(raw_input_df["Protein"]==prot)&(raw_input_df["File.Name"].isin(files_per_condition))].shape[0]
            prec_dict["#Precursors-IN-"+condition] = raw_input_df[(raw_input_df["Protein"]==prot)&(raw_input_df["File.Name"].isin(files_per_condition))].shape[0]
        precursor_info.append(prec_dict)
    
    precursor_info_df = pd.DataFrame.from_dict(precursor_info)
    precursor_info_df["# PSMs"] = precursor_info_df["#Precursors-IN-"+wt]
    # print(precursor_info_df.head())
    final_full_df = pd.merge(final_prot_df,precursor_info_df,on="Protein")
    print(final_full_df.head())
    
    return final_full_df

def norm_by_bait_DIA(bait_name, prot_df):
    channel_list = list(set(prot_df["Channel"].values.tolist()))
    # TODO: Print warning message that Bait is filtered out if bait no longer in protein level
    prot_df["check_substring"] = prot_df["Protein"].str.contains(bait_name)
    if prot_df["check_substring"].sum()>0:
        prot_df_bait = prot_df.loc[prot_df["Protein"]==bait_name]
        print(prot_df_bait)
        #peak_df_norm_intensity = peak_df_norm_intensity.reset_index()
        channels = channel_list
        scaling_factors = {}
        all_scaling_factors = []
        for i in range(len(channels)):
            #scaling_factors.append(statistics.median(peak_df_norm_intensity["log2_"+channels[i]+"_1"].dropna().values.tolist()))
            # Taking the avg instead
            #scaling_factors.append(np.average(prot_df_bait["log2_"+channels[i]+"_1"].dropna().values.tolist()))
            # scaling_factors[i]=[]
            sf = np.average(prot_df_bait.loc[(prot_df_bait["Channel"]==channels[i])]["log2Abundance"].dropna().values.tolist())
            scaling_factors[i] = sf
            if math.isnan(sf):
                continue
            else:
                all_scaling_factors.append(sf)
        
        global_median = statistics.median(all_scaling_factors)
        print(scaling_factors)
        #for i in range(len(scaling_factors)):
            #scaling_factors[i] = - scaling_factors[i] + global_median
        print(global_median)
        for i in range(len(channels)):
            scaling_factors[i] = - scaling_factors[i] + global_median
        print(scaling_factors)
        for i in range(len(channels)):
                
            if math.isnan(scaling_factors[i]):
                continue
            else:
                prot_df.loc[((prot_df["Channel"]==channels[i])),"log2Abundance"] += scaling_factors[i]
    else:
        print("Bait no longer in filtered Protein level file")
    
    return prot_df

def expansion_to_fractionXchannel_DIA(init_df,annotation_df, comparison, comp_type):
    print(annotation_df.head(2))
    wt=""
    ctrl=""
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print("Reducing file size with long form expansion: ",wt, ctrl)
    channel_list = list(set(annotation_df.loc[((annotation_df["Label"]==wt)|(annotation_df["Label"]==ctrl))]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)
    print("Reduced channel list for long form expansion: ",channel_list)
    # Do the same as above but for all peptides belonging to proteins
    prot_list = list(set(init_df["Protein"].dropna().values.tolist()))
    # master_prots_yulab = list(set(prot_init_df["Protein"].dropna().values.tolist()))
    prot_det = {
        'Protein':prot_list
    }
    # TODO - Make it so only channels present in annotation file are used below
    # prot_df = pd.DataFrame(columns=["Protein","log2Abundance","Channel","Condition","Fraction"])#peak_df_norm_ratio.copy()
    # summarize by median
    dict_list = []


    for i in range(len(prot_list)):
        prot = prot_list[i]
        df = init_df.loc[(init_df["Protein"]==prot)]
        for channel in channel_list:
            if len(df["log2_"+channel+"_1"].dropna().values.tolist())>0:
                dict_list.append({"Protein":prot,"log2Abundance":statistics.median(df["log2_"+channel+"_1"].dropna().values.tolist()),"Channel":channel,"Condition":annotation_df.loc[annotation_df["Channel"]==channel]["Label"].values.tolist()[0]})
    
    prot_df = pd.DataFrame.from_dict(dict_list)
    return prot_df

def create_limma_suitable_file_IP_DIA(wt, ctrl, prot_df, prot_list, imputation_ctrl_vals, prot_df_calc):
    # Assuming three replicates present
    num_reps = len(list(set(prot_df.loc[prot_df["Condition"]==wt]["Channel"].values.tolist())))
    ctrl_channels = list(set(prot_df.loc[prot_df["Condition"]==ctrl]["Channel"].values.tolist()))
    cols_rep = []
    for i in range(num_reps):
        cols_rep+=("rep"+str(i+1))
    # limma_prep = pd.DataFrame(columns=["Protein","Condition","Fraction"]+cols_rep)
    
    prot_list = list(set(prot_df["Protein"].values.tolist()))
    dict_list = []
    # conditions = [wt+"-"+ctrl]
    
    for prot in prot_list:
        
        rvs_T = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
        rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
        
        if len(rvs_T)<2:
            continue
        
        # Imputation for the case where there is no control but enough Treatment
            ## TODO - Make it more formal (random choice vs first n?)
        if (len(rvs_C)==0)&(len(imputation_ctrl_vals)>0):
            rvs_C = imputation_ctrl_vals[:len(rvs_T)]
            
        id_T = list(range(0,len(rvs_T)))
        id_C = list(range(0,len(rvs_C)))
        
        ## This version is only changing things for proteins that had no control
        if len(rvs_T) > len(rvs_C):
            df2 = {"Protein":prot}
            pos_C = list(permutations(id_T, len(id_C)))
            for possibility in pos_C:
                for id_ in range(len(possibility)):
                    df2["rep"+str(id_+1)] = (rvs_T[possibility[id_]]-rvs_C[id_])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
        elif len(rvs_T) < len(rvs_C):
            df2 = {"Protein":prot}
            pos_T = list(permutations(id_C, len(id_T)))
            for possibility in pos_T:
                for id_ in range(len(possibility)):
                    df2["rep"+str(id_+1)] = (rvs_T[id_]-rvs_C[possibility[id_]])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
        elif len(rvs_T) == len(rvs_C):
            df2 = {"Protein":prot}
            pos_T = list(permutations(id_C, len(id_C)))
            for possibility in pos_T:
                # add_row = [prot,c,frac,np.nan,np.nan,np.nan]
                for id_ in range(len(possibility)):
                    if id_ < len(id_T):
                        df2["rep"+str(id_+1)] = (rvs_T[id_T[id_]]-rvs_C[id_C[possibility[id_]]])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
    limma_prep = pd.DataFrame.from_dict(dict_list)
    return limma_prep

def create_limma_suitable_file_wholeproteome_DIA(wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_df_calc):
    # Assuming three replicates present
    num_reps = len(list(set(prot_df.loc[prot_df["Condition"]==wt]["Channel"].values.tolist())))
    ctrl_channels = list(set(prot_df.loc[prot_df["Condition"]==ctrl]["Channel"].values.tolist()))
    cols_rep = []
    for i in range(num_reps):
        cols_rep+=("rep"+str(i+1))
    # limma_prep = pd.DataFrame(columns=["Protein","Condition","Fraction"]+cols_rep)
    
    prot_list = list(set(prot_df["Protein"].values.tolist()))
    dict_list = []
    # conditions = [wt+"-"+ctrl]
    
    for prot in prot_list:
        
        rvs_T = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
        rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
        
        if (len(rvs_T)<2)|(len(rvs_C)<2):
            continue
        
        # Imputation for the case where there is no control but enough Treatment
            ## TODO - Make it more formal (random choice vs first n?)
        if (len(rvs_C)==0)&(len(imputation_ctrl_vals)>0):
            rvs_C = imputation_ctrl_vals[:len(rvs_T)]

        if (len(rvs_T)==0)&(len(imputation_wt_vals)>0):
            rvs_T = imputation_wt_vals[:len(rvs_C)]

            
        id_T = list(range(0,len(rvs_T)))
        id_C = list(range(0,len(rvs_C)))
        
        ## This version is only changing things for proteins that had no control
        if len(rvs_T) > len(rvs_C):
            df2 = {"Protein":prot}
            pos_C = list(permutations(id_T, len(id_C)))
            for possibility in pos_C:
                for id_ in range(len(possibility)):
                    df2["rep"+str(id_+1)] = (rvs_T[possibility[id_]]-rvs_C[id_])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
        elif len(rvs_T) < len(rvs_C):
            df2 = {"Protein":prot}
            pos_T = list(permutations(id_C, len(id_T)))
            for possibility in pos_T:
                for id_ in range(len(possibility)):
                    df2["rep"+str(id_+1)] = (rvs_T[id_]-rvs_C[possibility[id_]])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
        elif len(rvs_T) == len(rvs_C):
            df2 = {"Protein":prot}
            pos_T = list(permutations(id_C, len(id_C)))
            for possibility in pos_T:
                # add_row = [prot,c,frac,np.nan,np.nan,np.nan]
                for id_ in range(len(possibility)):
                    if id_ < len(id_T):
                        df2["rep"+str(id_+1)] = (rvs_T[id_T[id_]]-rvs_C[id_C[possibility[id_]]])
                shallow_copy = {}
                for k,val in df2.items():
                    shallow_copy[k] = val
                dict_list.append(shallow_copy)
    limma_prep = pd.DataFrame.from_dict(dict_list)
    return limma_prep

def calc_FC_withratios_IP_DIA(wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_df_calc):
    for prot in prot_list:
        per_protein_FC = []
        rvs_T = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
        rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
        
        if (len(rvs_C)>0)|(len(rvs_T)>0):
            if (len(rvs_C)==0)&(len(imputation_ctrl_vals)>0):
                rvs_C = imputation_ctrl_vals[:len(rvs_T)]
            if (len(rvs_T)==0)&(len(imputation_wt_vals)>0):
                rvs_T = imputation_wt_vals[:len(rvs_C)]
        for t in rvs_T:
            for c in rvs_C:
                per_protein_FC.append(np.power(2,t-c))
        prot_df_calc.loc[prot_df_calc["Protein"]==prot,"FC"] = np.nanmedian(per_protein_FC)
        prot_df_calc.loc[prot_df_calc["Protein"]==prot,"#wt"] = len(prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist())
        prot_df_calc.loc[prot_df_calc["Protein"]==prot,"#ctrl"] = len(prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist())
    
    # print(prot_df_calc.loc[prot_df_calc["Protein"]=="P04803"])
    return prot_df_calc#_dropna

def calc_pval_withratios_IP_DIA(job_id, wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_df_calc):
    prot_df_calc["pval"] = np.nan
    limma_pre_df = []
    if len(imputation_wt_vals)>0:
        limma_pre_df = create_limma_suitable_file_wholeproteome_DIA(wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_df_calc)
    elif len(imputation_ctrl_vals)>0:
        limma_pre_df = create_limma_suitable_file_IP_DIA(wt, ctrl, prot_df, prot_list, imputation_ctrl_vals, prot_df_calc)
    else:
        limma_pre_df = create_limma_suitable_file_wholeproteome_DIA(wt, ctrl, prot_df, prot_list, [], [], prot_df_calc)
    # print(limma_pre_df.loc[limma_pre_df["Protein"]=="P04803"])
    print("Done limma setup with IMPUTATION BASED ON CONTROL")

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    
    # local_storage_dir= "/Users/user/Downloads/TMT_Analysis/Paper_analysis/DIA-LFQ-analyses/"#BASE_DIR+'/media/documents/tmtipms_files'#"/Users/user/Downloads/TMT_Analysis/Paper_analysis/SEQUEST/"
    
    filename_lim = local_storage_dir+"/"+job_id+"/"+wt+"-"+ctrl+"_IPDIA_imputation"
    
    limma_pre_df.to_csv(filename_lim+".txt",sep="\t")
    limma_post_df = run_r_script(limma_pre_df,filename_lim)
    # print(limma_post_df.loc[limma_post_df["Protein"]=="P04803"])
    limma_summarized = limma_post_df.groupby(['Protein'])['pvalue','FDR'].mean().reset_index()
    prot_df_calc["pval"] = prot_df_calc.Protein.map(limma_summarized.set_index("Protein")["pvalue"].to_dict())
    prot_df_calc["adjpval"] = prot_df_calc.Protein.map(limma_summarized.set_index("Protein")["FDR"].to_dict())
    return prot_df_calc

def calc_FC_and_pval_limma_IP_DIA(job_id, prot_df, init_df, comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation, imputation_type="ip"):
    prot_list = list(set(prot_df["Protein"].dropna().values.tolist()))
    
    prot_det = {
        'Protein':prot_list
    }

    prot_median_manual = pd.DataFrame(prot_det)

    prot_median_manual["FC"] = np.nan
    prot_median_manual["pval"] = np.nan
    prot_median_manual["log2FC"] = np.nan
    wt = ""
    ctrl = ""
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print(wt, ctrl)
    print(list(set(prot_df["Condition"].values.tolist())))

    if imputation:
        print("Doing an imputation")
        ctrl_channels = list(set(prot_df.loc[prot_df["Condition"]==ctrl]["Channel"].values.tolist()))
        imputation_ctrl_vals = []
        for channel_ctrl in ctrl_channels:
            imputation_ctrl_vals.append(np.min(prot_df.loc[prot_df["Channel"]==channel_ctrl][["log2Abundance"]].dropna()["log2Abundance"].values.tolist()))
        
        imputation_wt_vals = []
        if imputation_type=="wholeproteome":
            print("Imputation for whole proteome dataset")
            wt_channels = list(set(prot_df.loc[prot_df["Condition"]==wt]["Channel"].values.tolist()))
            for channel_wt in wt_channels:
                imputation_wt_vals.append(np.min(prot_df.loc[prot_df["Channel"]==channel_wt][["log2Abundance"]].dropna()["log2Abundance"].values.tolist()))
        
        # Calculate FC using limma (ratios)
        prot_median_manual = calc_FC_withratios_IP_DIA(wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_median_manual)
        print("Calculated FC successfully with IMPUTATION specifically when no control values found")
        # Calculate p-value using limma (ratios)
        prot_median_manual = calc_pval_withratios_IP_DIA(job_id, wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_median_manual)
        print("Calculated pval successfully with limma and IMPUTATION specifically when no control values found")
    else:
        prot_median_manual = calc_FC_withratios_IP_DIA(wt, ctrl, prot_df, prot_list, [], [], prot_median_manual)
        print("Calculated FC successfully")
        # Calculate p-value using limma (ratios)
        prot_median_manual = calc_pval_withratios_IP_DIA(job_id, wt, ctrl, prot_df, prot_list, [], [], prot_median_manual)
        print("Calculated pval with limma successfully")

    prot_median_manual["log2FC"] = np.log2(prot_median_manual["FC"])
    print(prot_median_manual.head(2))

    prot_median_manual_dropna = prot_median_manual[["Protein","log2FC","pval"]].dropna()
    print(prot_median_manual_dropna.shape[0])


    prot_median_manual_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(prot_median_manual_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    if ((bait_norm==True) & (bait_name not in prot_median_manual_dropna["Protein"].values.tolist())):
        print("Bait not added in list after bait normalization")
        prot_median_manual_dropna.loc[len(prot_median_manual_dropna.index)] = [bait_name, 0, 1, 1] 

    return prot_median_manual_dropna

def calc_FC_and_pval_IP_DIA(job_id, prot_df, init_df, comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation, imputation_type="ip"):
    prot_list = list(set(prot_df["Protein"].dropna().values.tolist()))
    
    prot_det = {
        'Protein':prot_list
    }

    prot_median_manual = pd.DataFrame(prot_det)

    prot_median_manual["FC"] = np.nan
    prot_median_manual["pval"] = np.nan
    prot_median_manual["log2FC"] = np.nan
    wt = ""
    ctrl = ""
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print(wt, ctrl)
    print(list(set(prot_df["Condition"].values.tolist())))

    if imputation:
        print("Doing an imputation")
        
        imputation_ctrl_vals = []
        ctrl_channels = list(set(prot_df.loc[prot_df["Condition"]==ctrl]["Channel"].values.tolist()))
        for channel_ctrl in ctrl_channels:
            imputation_ctrl_vals.append(np.min(prot_df.loc[prot_df["Channel"]==channel_ctrl][["log2Abundance"]].dropna()["log2Abundance"].values.tolist()))
        
        imputation_wt_vals = []
        if imputation_type=="wholeproteome":
            print("Imputation for whole proteome dataset")
            wt_channels = list(set(prot_df.loc[prot_df["Condition"]==wt]["Channel"].values.tolist()))
            for channel_wt in wt_channels:
                imputation_wt_vals.append(np.min(prot_df.loc[prot_df["Channel"]==channel_wt][["log2Abundance"]].dropna()["log2Abundance"].values.tolist()))
            
        # Calculate FC (ratios)
        prot_median_manual = calc_FC_withratios_IP_DIA(wt, ctrl, prot_df, prot_list, imputation_wt_vals, imputation_ctrl_vals, prot_median_manual)
        print("Calculated FC successfully with IMPUTATION specifically when no control")
        # Calculate p-value using 2 sample ttest (ratios)
        for prot in prot_list:
            rvs_N = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
            rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
            if (len(rvs_C)>2)|(len(rvs_N)>2):
                if (len(rvs_C)==0)&(len(imputation_ctrl_vals)>0):
                    rvs_C = imputation_ctrl_vals[:len(rvs_N)]
                if (len(rvs_N)==0)&(len(imputation_wt_vals)>0):
                    rvs_N = imputation_wt_vals[:len(rvs_C)]
            # prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C, alternative='greater').pvalue
            prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C).pvalue
        print("Calculated pval 2 sample ttest successfully with IMPUTATION specifically when no control")
    else:
        prot_median_manual = calc_FC_withratios_IP_DIA(wt, ctrl, prot_df, prot_list, [], [], prot_median_manual)
        print("Calculated FC successfully")
        # Calculate p-value using 2 sample ttest (ratios)
        for prot in prot_list:
            rvs_N = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
            rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
            if len(rvs_C)==0:
                continue
            # prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C, alternative='greater').pvalue
            prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C).pvalue
        print("Calculated pval 2 sample ttest successfully")

    prot_median_manual["log2FC"] = np.log2(prot_median_manual["FC"])
    print(prot_median_manual.head(2))

    prot_median_manual_dropna = prot_median_manual[["Protein","log2FC","pval"]].dropna()
    print(prot_median_manual_dropna.shape[0])


    prot_median_manual_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(prot_median_manual_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    if ((bait_norm==True) & (bait_name not in prot_median_manual_dropna["Protein"].values.tolist())):
        print("Bait not added in list after bait normalization")
        prot_median_manual_dropna.loc[len(prot_median_manual_dropna.index)] = [bait_name, 0, 1, 1] 


    return prot_median_manual_dropna

def YuLab_fractionxchannel_DIA(job_id, init_df, raw_input_df, annotation_df, bait_name, comparison, comp_type, row_norm_user,column_norm_user, bait_norm, num_channels, norm_prots, pval_calc, uniprotGene, fusion_aa_seq='', imputation_user=False, iso_filter=100):
    # Making sure dependence is not on number channels but on specific ones with labels
    channel_list = list(set(annotation_df.loc[annotation_df["Label"]!="none"]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)
    
    # At protein level assumption is every row to a unique protein
    print("Channel list: ", channel_list)

    
    # TODO - Add the PSM file input to give # Precursors per condition
    print(list(init_df))
    init_df = preprocess_byrownorm(init_df,True,channel_list)#row_norm_user


    if column_norm_user==True:
        if norm_prots == '':
            print("No protein list based normalization")
            init_df = normalization_median_intensity(init_df,channel_list)
            print("Done median ratio based normalization")
        else:
            uniprots_to_norm_by = norm_prots.split(",")
            print("Yes protein list based normalization: ")#, uniprots_to_norm_by)
            init_df = norm_by_prots_byTMM(init_df,channel_list, uniprots_to_norm_by)
            print("Done median ratio based normalization using the proteins listed")

    # Minimum value per column added
    # Here imputation to be done only for Proteins that have ctrl empty
    # if imputation_user==True:
    #     init_df = impute_peptide_level(init_df,channel_list)
    
    # TODO convert the prot_init_df to expanded fraction x channel df
    # Cols - Protein, Condition, Channel, log2Abundance. Add fraction as 1.raw
    prot_df = expansion_to_fractionXchannel_DIA(init_df,annotation_df,comparison,comp_type)
    # prot_df["Fraction"] = "1.raw"

    final_prot_df = []
    # if pval_calc == "LinearModel":
    #     print("For DIA not doing LM just Limma")
    #     final_prot_df = calc_FC_and_pval(prot_df, [], comparison, comp_type, uniprotGene, '', False)
    # elif pval_calc == "Limma":
    #     print("For DIA Limma - whole proteome")
    #     # final_prot_df = calc_FC_and_pval_limma_DIA(prot_init_df, comparison, comp_type, uniprotGene)
    #     final_prot_df = calc_FC_and_pval_limma(prot_df, [], comparison, comp_type, uniprotGene, '', False)
    # else:
    #     print("Error: Magma Type not chosen")

    if pval_calc == "LinearModel":
        print("For DIA Whole Proteome doing LM")
        final_prot_df = calc_FC_and_pval_IP_DIA(job_id, prot_df, [], comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation_user,"wholeproteome")
    elif (pval_calc == "Limma"):
        print("For DIA Whole Proteome doing Limma")
        final_prot_df = calc_FC_and_pval_limma_IP_DIA(job_id, prot_df, [], comparison, comp_type, uniprotGene, bait_name, bait_norm, imputation_user,"wholeproteome")
    else:
        print("Error: Magma Type not chosen")
    
    print("Done calculating FC and p-value")
    return final_prot_df

def preprocess_byrownorm(df,rownorm, cols_abund):
    # if rownorm:
    df["Norm-factor"] = df[cols_abund].mean(axis=1)
    for c in cols_abund:
        df[c+"_rownorm"] = df[c]
        df[c] = df[c+"_rownorm"]/df["Norm-factor"]

    for c in cols_abund:
        df["log2_"+c+"_1"] = np.log2(df[c])

    return df



def COMETconversion(sn_df,annotation_df,num_channels):
    sn_df["Spectrum File"] = sn_df["spectrum"].str.split(".",expand=True)[0]
    sn_df["Isolation Interference [%]"] = 0
    sn_df["Intensity"] = sn_df["precursor_intensity"]
    sn_df["Annotated Sequence"] = sn_df["peptide"]#.str.replace(r'\[[^\]]*\]',"").str.split(".",expand=True)[1].str[1:]#.str.replace(r'\[[^()]*\]',"")
    sn_df["Sequence"] = sn_df["peptide"].str.replace(r'\[[^\]]*\]',"").str.split(".",expand=True)[1].str[1:].str.upper()
    sn_df["Charge"] = sn_df["assumed_charge"]
    sn_df["Master Protein Accessions"] = sn_df["protein"].str.split(",",expand=True)[0].str.split("|",expand=True)[1]
    #sn_df["# Protein Groups"] = 1#sn_df["num_prots"]#1
    sn_df["RT [min]"] = sn_df["retention_time_sec"]/60.0
    # tmt_channels = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
    tmt_channels = annotation_df["Channel"].values.tolist()
    # # here 10 tmt channels
    # if num_channels == 6:
    #     tmt_channels = ["126","127","128","129","130","131"]
    # elif num_channels == 10:
    #     tmt_channels = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
    # else:
    #     tmt_channels = tmt_channels[:num_channels]
    
    for i in range(len(tmt_channels)):
        sn_df["Abundance: "+str(tmt_channels[i])] = sn_df["libra"+str(i+1)]
        
    # For COMET specifically need to make 0's empty channels in abundance columns
    cols_abundance = [x for x in list(sn_df) if "Abundance: " in x]
    sn_df[cols_abundance] = sn_df[cols_abundance].replace({'0':np.nan, 0:np.nan})
        
    # Using all except ambiguous ID's (in the case of )
    sn_df = sn_df.loc[~(sn_df["protein"].str.contains("HUMAN")&sn_df["protein"].str.contains("YEAST"))]
    sn_df["# Protein Groups"] = 1
    return sn_df

# 20230518 Imputation added (NOT TESTED ALSO WHERE TO PUT IN THE WORKFLOW)
def impute_peptide_level(peak_df,channels):
    peak_df_init = peak_df.copy()
    cols_abund = []
    for channel in channels:
        cols_abund.append("log2_"+str(channel)+"_1")

    min_val_per_col = peak_df_init[cols_abund].min(axis=0).values.tolist() # Default - ignores NA, gets the minimum value per column

    peak_df_imputed = peak_df_init.copy()
    for i in range(len(cols_abund)):
        print("Before impute col: ",cols_abund[i],peak_df_imputed[cols_abund[i]].isna().sum())
        peak_df_imputed[cols_abund[i]] = peak_df_imputed[cols_abund[i]].fillna(min_val_per_col[i])
        print("After impute col: ",cols_abund[i],peak_df_imputed[cols_abund[i]].isna().sum())
    print("Imputed missing values with minimum per column at the peptide level")
    return peak_df_imputed

def YuLab_phosphopep_fractionxchannel(job_id, sn_df, intensity_df, annotation_df, bait_name, comparison, comp_type, fraction_use, row_norm_user,column_norm_user, bait_norm, num_channels, norm_prots, pval_calc, uniprotGene, fusion_aa_seq='', imputation_user=False, iso_filter=100):
    # Making sure dependence is not on number channels but on specific ones with labels
    channel_list = list(set(annotation_df.loc[annotation_df["Label"]!="none"]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)

    tmp_channels_print = annotation_df.loc[(annotation_df["Label"]==comparison.split("-")[0])|(annotation_df["Label"]==comparison.split("-")[1])]["Channel"].values.tolist()
    tmp_channels_print = sort1(tmp_channels_print)
    # print(sn_df.loc[sn_df["Master Protein Accessions"]==bait_name][[x for x in list(sn_df) if (("Abundance:" in x))]])

    if fraction_use==False:
        sn_df["Spectrum File"] = "F1.raw"
        intensity_df["Spectrum File"] = "F1.raw"

    peak_df = psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, channel_list, iso_filter)
    # augment_df = peak_df.copy()
    
    print("Done PSM summarization", peak_df.head(2))
    init_psm_df = compile_sn_raw_int(sn_df,intensity_df,channel_list)
    if column_norm_user==True:
        if norm_prots == '':
            print("No protein list based normalization")
            peak_df = normalization_median_intensity(peak_df,channel_list)
            print("Done median ratio based normalization")
        else:
            uniprots_to_norm_by = norm_prots.split(",")
            print("Yes protein list based normalization: ")#, uniprots_to_norm_by)
            peak_df = norm_by_prots_byTMM(peak_df,channel_list, uniprots_to_norm_by)
            print("Done median ratio based normalization using the proteins listed")


    if imputation_user==True:
        peak_df = impute_peptide_level(peak_df,channel_list)

    # TODO Peptide level file to be returned (not log transformed after column normalization)
    print("Calculating peptide level FC, FDR estimates")
    peak_df[["AnnotatedSeq","Charge","Fraction"]] = peak_df["PeakID"].str.split("_",expand=True)
    peak_df_deliverable = cleanup_peak_df_phosphopep(peak_df, annotation_df, comparison, comp_type, channel_list)
    peak_df_deliverable["Proteinnoiso"] = peak_df_deliverable["Protein"].str.split("-",expand=True)[0]
    peak_df_deliverable["Gene Symbol"] = peak_df_deliverable.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
    
    if isinstance(init_psm_df, pd.DataFrame):
        forpsm_df = init_psm_df
        prot_list = list(set(peak_df_deliverable["Protein"].dropna().values.tolist()))
        peak_df_deliverable["# PSMs"] = 0
        for p in prot_list:
            fil = forpsm_df.loc[(forpsm_df["Master Protein Accessions"].str.contains(p, na=False))]
            peak_df_deliverable.loc[peak_df_deliverable["Protein"]==p,"# PSMs"] = fil.shape[0]
    print("Done calculating FC and p-value")
    return (peak_df_deliverable)

def YuLab_fractionxchannel(job_id, sn_df, intensity_df, annotation_df, bait_name, comparison, comp_type, fraction_use, row_norm_user,column_norm_user, bait_norm, num_channels, norm_prots, pval_calc, uniprotGene, fusion_aa_seq='', imputation_user=False, iso_filter=100):
    # Making sure dependence is not on number channels but on specific ones with labels
    channel_list = list(set(annotation_df.loc[annotation_df["Label"]!="none"]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)

    tmp_channels_print = annotation_df.loc[(annotation_df["Label"]==comparison.split("-")[0])|(annotation_df["Label"]==comparison.split("-")[1])]["Channel"].values.tolist()
    tmp_channels_print = sort1(tmp_channels_print)
    # print(sn_df.loc[sn_df["Master Protein Accessions"]==bait_name][[x for x in list(sn_df) if (("Abundance:" in x))]])

    if fraction_use==False:
        sn_df["Spectrum File"] = "F1.raw"
        intensity_df["Spectrum File"] = "F1.raw"

    peak_df = psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, channel_list, iso_filter)
    augment_df = peak_df.copy()
    
    print("Done PSM summarization", peak_df.head(2))
    init_psm_df = compile_sn_raw_int(sn_df,intensity_df,channel_list)
    if column_norm_user==True:
        if norm_prots == '':
            print("No protein list based normalization")
            peak_df = normalization_median_intensity(peak_df,channel_list)
            print("Done median ratio based normalization")
        else:
            uniprots_to_norm_by = norm_prots.split(",")
            print("Yes protein list based normalization: ")#, uniprots_to_norm_by)
            peak_df = norm_by_prots_byTMM(peak_df,channel_list, uniprots_to_norm_by)
            print("Done median ratio based normalization using the proteins listed")


    if imputation_user==True:
        peak_df = impute_peptide_level(peak_df,channel_list)

    prot_df = expansion_to_fractionXchannel_faster(peak_df,annotation_df, comparison, comp_type)
    #expansion_to_fractionXchannel(peak_df,annotation_df, channel_list) 
    print(set(prot_df.loc[prot_df["Protein"]==bait_name]["Fraction"].values.tolist()))
    print(prot_df.loc[(prot_df["Protein"]==bait_name)&(prot_df["Channel"].isin(tmp_channels_print))])
    print("Done setup for FC and pvalue calculation", prot_df.shape[0])
    
    # Bait normalization
    if bait_norm==True:
        print("Normalized by bait")
        prot_df = norm_by_bait(bait_name, prot_df, channel_list)
    else:
        print("No bait based normalization")

    
    final_prot_df = []
    if pval_calc == "LinearModel":
        final_prot_df = calc_FC_and_pval(prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm)
    elif pval_calc == "Limma":
        final_prot_df = calc_FC_and_pval_limma(job_id, prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm)
    else:
        print("Error: Magma Type not chosen")

    if (bait_name != "None")|(bait_name != "NA"):
        psm_cutoff = psm_cutoff_suggestion(bait_name, final_prot_df)

    # Add per channel per fraction abundance 
    final_prot_df = augment_final_df(final_prot_df, prot_df)
    print("Augmented final deliverable excel sheet")

    # TODO Peptide level file to be returned (not log transformed after column normalization)
    print("Calculating peptide level FC, FDR estimates")
    peak_df_deliverable = cleanup_peak_df(peak_df, annotation_df, comparison, comp_type, channel_list)
    print("Cleaned peptide deliverable excel sheet")

    final_prot_df.at[0,"PSM Cutoff"] = psm_cutoff
    print("Done calculating FC and p-value")
    return (final_prot_df,peak_df_deliverable)

def YuLab_fusion_fractionxchannel(job_id, sn_df, intensity_df, annotation_df, bait_name, comparison, comp_type, fraction_use, row_norm_user,column_norm_user, bait_norm, num_channels, norm_prots, pval_calc, uniprotGene, fusion_aa_seq='', imputation_user=False, iso_filter=100):
    # Making sure dependence is not on number channels but on specific ones with labels
    channel_list = list(set(annotation_df.loc[annotation_df["Label"]!="none"]["Channel"].values.tolist()))
    channel_list = sort1(channel_list)

    if fraction_use==False:
        sn_df["Spectrum File"] = "F1.raw"
        intensity_df["Spectrum File"] = "F1.raw"

    # TODO add option for fraction separation
    peak_df = psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, channel_list, iso_filter)
    augment_df = peak_df.copy()
    print("Done PSM summarization", peak_df.head(2))
    init_psm_df = compile_sn_raw_int(sn_df,intensity_df,channel_list)

    
    if column_norm_user==True:
        if norm_prots == '':
            print("No protein list based normalization")
            peak_df = normalization_median_intensity(peak_df,channel_list)
            print("Done median ratio based normalization")
        else:
            uniprots_to_norm_by = norm_prots.split(",")
            print("Yes protein list based normalization: ")#, uniprots_to_norm_by)
            peak_df = norm_by_prots_byTMM(peak_df,channel_list, uniprots_to_norm_by)
            print("Done median ratio based normalization using the proteins listed")

    # TODO: For fusion versus control specifically the check is with ";" in bait_name in which case bait_norm should be false
    print("Before ",peak_df.shape)
    if (bait_norm==False) & (bait_name!="NA"):
        if ";" in bait_name:
            baits = bait_name.split(";")
            for bait in baits:
                # Remove peptides specific to both head and tail that are not found in fusion
                # print("Initial list of peptides that are bait specific: ", peak_df.loc[peak_df["Protein"]==bait]["PeakID"].values.tolist())
                print("Removing non specific peptides for FvsC for bait: ",bait)
                peak_df = remove_non_fusion_peptides(peak_df, fusion_aa_seq, bait)
                print("Bait specific peptides left: ",len(peak_df.loc[peak_df["Protein"]==bait]["PeakID"].values.tolist()),peak_df.loc[peak_df["Protein"]==bait]["PeakID"].values.tolist())

            # By default for PSM cutoff etc. now use first uniprot in list
            bait_name = baits[0]
    elif (bait_norm==True):
        peak_df = remove_non_fusion_peptides(peak_df, fusion_aa_seq, bait_name)
        #print("Bait specific peptides before fusion removal: ",peak_df.loc[peak_df["Protein"]==bait_name]["PeakID"].values.tolist())
        #print("After non specific peptide removal of bait: ", peak_df.head())
        #print(len(list(peak_df)),list(peak_df))
        # print("Bait specific peptides left: ",len(peak_df.loc[peak_df["Protein"]==bait_name]["PeakID"].values.tolist()),peak_df.loc[peak_df["Protein"]==bait_name]["PeakID"].values.tolist())

    print("After ",peak_df.shape)
    if imputation_user==True:
        peak_df = impute_peptide_level(peak_df,channel_list)

    prot_df = expansion_to_fractionXchannel_faster(peak_df,annotation_df, comparison, comp_type)
    #expansion_to_fractionXchannel(peak_df,annotation_df, channel_list) 
    print("Done setup for FC and pvalue calculation", prot_df.shape[0])

    # Bait normalization
    if bait_norm==True:
        print("Normalized by bait")
        prot_df = norm_by_bait(bait_name, prot_df, channel_list)
    else:
        print("No bait based normalization")

    
    final_prot_df = []
    if pval_calc == "LinearModel":
        final_prot_df = calc_FC_and_pval(prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm)
    elif pval_calc == "Limma":
        final_prot_df = calc_FC_and_pval_limma(job_id, prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm)
    else:
        print("Error: Magma Type not chosen")
    
    if (bait_name != "None")|(bait_name != "NA"):
        psm_cutoff = psm_cutoff_suggestion(bait_name, final_prot_df)

    # Add per channel per fraction abundance 
    # TODO - make the option based on fraction choice
    final_prot_df = augment_final_df(final_prot_df, prot_df)
    print("Augmented final delilverable excel sheet")

    # TODO Peptide level file to be returned (not log transformed after column normalization)
    # Make this an option to output 
    print("Calculating peptide level FC, FDR estimates")
    peak_df_deliverable = cleanup_peak_df(peak_df, annotation_df, comparison, comp_type, channel_list)
    print("Cleaned peptide deliverable excel sheet")

    final_prot_df.at[0,"PSM Cutoff"] = psm_cutoff
    print("Done calculating FC and p-value")
    return (final_prot_df,peak_df_deliverable)

def cleanup_peak_df_phosphopep(peak_df, annotation_df, comparison, comp_type, channel_list):
    cols_to_keep = ["PeakID","AnnotatedSeq","Charge","Fraction","Protein"]

    for col in channel_list:
        peak_df[col] = np.power(2,peak_df["log2_"+col+"_1"])
        cols_to_keep.append(col)

    # TODO Calculate log2FC, p-value and adjusted p-value per row here
    peak_df["log2FC"] = np.nan
    peak_df["pval"] = np.nan

    for i in range(peak_df.shape[0]):
        peak = peak_df.iloc[i]["PeakID"]
        prot = peak_df.iloc[i]["Protein"]
        wt = ""
        ctrl = ""
        if comp_type==1:
            wt = comparison.split("-")[0]
            ctrl = comparison.split("-")[1]
        else:
            wt = comparison.split("-")[1]
            ctrl = comparison.split("-")[0]

        # print(wt, ctrl)
        
        tmp_wt_channels = list(set(annotation_df.loc[annotation_df["Label"]==wt]["Channel"].values.tolist()))
        tmp_ctrl_channels = list(set(annotation_df.loc[annotation_df["Label"]==ctrl]["Channel"].values.tolist()))
        wt_channels = []
        ctrl_channels = []
        
        for w in tmp_wt_channels:
            wt_channels.append("log2_"+str(w)+"_1")

        for c in tmp_ctrl_channels:
            ctrl_channels.append("log2_"+str(c)+"_1")

        # Calculate FC
        wt_vals = []
        if (len(peak_df.loc[peak_df["PeakID"]==peak][wt_channels].dropna().values.tolist())>0):
            wt_vals = peak_df.loc[peak_df["PeakID"]==peak][wt_channels].dropna().values.tolist()[0]
        ctrl_vals = []
        
        if (len(peak_df.loc[peak_df["PeakID"]==peak][ctrl_channels].dropna().values.tolist())>0):
            ctrl_vals = peak_df.loc[peak_df["PeakID"]==peak][ctrl_channels].dropna().values.tolist()[0]

        if (len(wt_vals)==0)|(len(ctrl_vals)==0):
            continue
        
        peak_df.loc[peak_df["PeakID"]==peak,"log2FC"] = np.mean(wt_vals) - np.mean(ctrl_vals)
        # print("Calculated FC per peak successfully")
        # Calculate p-value using t-test
        if (len(wt_vals)>1)&(len(ctrl_vals)>1):
            peak_df.loc[peak_df["PeakID"]==peak,"pval"] = stats.ttest_ind(wt_vals, ctrl_vals).pvalue
            # print("Calculated pval per peak successfully")
    cols_to_keep += ["log2FC","pval"]
    peak_df = peak_df[cols_to_keep].dropna()
    peak_df["adjpval"] = statsmodels.stats.multitest.multipletests(peak_df["pval"].values.tolist(),method="fdr_bh")[1]

    cols_to_keep += ["adjpval"]
    peak_df = peak_df[cols_to_keep]

    return peak_df

def cleanup_peak_df(peak_df, annotation_df, comparison, comp_type, channel_list):
    cols_to_keep = ["PeakID","AnnotatedSeq","Charge","Fraction","Protein"]

    for col in channel_list:
        peak_df[col] = np.power(2,peak_df["log2_"+col+"_1"])
        cols_to_keep.append(col)

    peak_df = peak_df[cols_to_keep]

    return peak_df

### Function to add FC and P-value information for every Peak
def augment_peak_df(peak_df, annotation_df, comparison, comp_type):
    
    print("# rows before peak dropna: ",peak_df.shape[0])
    peak_df_dropna = peak_df[["PeakID","Protein","pval"]].dropna()
    print("# rows after peak dropna: ",peak_df_dropna.shape[0])

    peak_df_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(peak_df_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    return peak_df_dropna

def augment_final_df_fractionsep(final_prot_df, prot_df):
    prot_df["FractionXChannel"] = prot_df["Fraction"]+"_"+prot_df["Channel"]

    cols_to_add = list(set(prot_df["FractionXChannel"].values.tolist()))

    for col in cols_to_add:
        final_prot_df["log2_"+col] = np.nan

    prot_list = list(set(final_prot_df["Protein"].values.tolist()))

    #log2Abundance
    for prot in prot_list:
        cols_to_add = list(set(prot_df.loc[prot_df["Protein"]==prot]["FractionXChannel"].values.tolist()))
        for col in cols_to_add:
            #if (prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["FractionXChannel"]==col)].shape[0]>0):
            final_prot_df.loc[final_prot_df["Protein"]==prot,"log2_"+col] = list(set(prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["FractionXChannel"]==col)]["log2Abundance"].values.tolist()))[0]
    return final_prot_df

def augment_final_df(final_prot_df, prot_df):
    # prot_df["FractionXChannel"] = prot_df["Fraction"]+"_"+prot_df["Channel"]

    cols_to_add = list(set(prot_df["Channel"].values.tolist()))

    for col in cols_to_add:
        final_prot_df["log2_"+col] = np.nan

    prot_list = list(set(final_prot_df["Protein"].values.tolist()))

    #log2Abundance
    for prot in prot_list:
        cols_to_add = list(set(prot_df.loc[prot_df["Protein"]==prot]["Channel"].values.tolist()))
        for col in cols_to_add:
            #if (prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Channel"]==col)].shape[0]>0):
            final_prot_df.loc[final_prot_df["Protein"]==prot,"log2_"+col] = statistics.median(prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Channel"]==col)]["log2Abundance"].values.tolist())
    return final_prot_df


def remove_non_fusion_peptides(peak_df, fusion_aa_seq, bait_name):
    peak_df[["AnnotatedSeq","Charge","Fraction"]] = peak_df["PeakID"].str.split("_",expand=True)
    peak_df = peak_df.reset_index()
    print("Before non fusion peptide removal #peptides: ",peak_df.shape[0])
    #Remove Peaks that belong to non fusion protein region
    #peaks_bait_wt = list(set(peak_df.loc[peak_df[]==bait_name]["AnnotatedSeq"].values.tolist()))
    #pep_seq_bait_wt = [pep.split(".")[1].toUpperCase() for pep in peaks_bait_df]
    peak_df["Sequence"] = peak_df.AnnotatedSeq.str.split(".").str[1].str.upper()
    peak_df["Check_fusion"] = False
    for i in range(peak_df.shape[0]):
        # if peak_df.iloc[i]["Protein"]==bait_name:
        #     print(peak_df.iloc[i]["Sequence"]," is it in? ",fusion_aa_seq)
        if peak_df.iloc[i]["Sequence"] in fusion_aa_seq:
            peak_df.at[i,"Check_fusion"]=True
    peak_df = peak_df.loc[~((peak_df["Protein"]==bait_name)&(peak_df["Check_fusion"]==False))]
    print("After non fusion peptide removal #peptides: ",peak_df.shape[0])
    return peak_df

def psm_cutoff_suggestion(bait_name, final_prot_df):
    psm_cutoff = 5 # The usual cutoff
    bait_psm = 5
    if bait_name in list(set(final_prot_df["Gene Symbol"].values.tolist())):
        bait_psm = final_prot_df.loc[final_prot_df["Gene Symbol"]==bait_name]["# PSMs"].values.tolist()[0]
    elif bait_name in list(set(final_prot_df["Protein"].values.tolist())):
        bait_psm = final_prot_df.loc[final_prot_df["Protein"]==bait_name]["# PSMs"].values.tolist()[0]
    
    if bait_psm>=1000:
        psm_cutoff = 15
    elif ((bait_psm>=400) & (bait_psm<1000)):
        psm_cutoff = 10
    elif (bait_psm < 400):
        psm_cutoff = 5
    elif (bait_psm < 5):
        psm_cutoff = bait_psm
    return psm_cutoff



def norm_by_prots_byTMM(peak_df,channel_list,norm_prots):
    #if num_channels==10:

    for channel in channel_list:
        peak_df[channel+"_r"] = np.power(2,peak_df["log2_"+channel+"_1"])
        peak_df[channel] = peak_df[channel+"_r"]/peak_df[channel+"_r"].sum()

    peak_df_norm_ratio = peak_df.loc[peak_df["Protein"].isin(norm_prots)].copy()

    f=[]
    reference_ind=0
    channels = channel_list#["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
    f=[]
    for i in range(len(channels)):
        peak_df_norm_ratio["A_"+channels[i]+"_r"] = np.log2(peak_df_norm_ratio[channels[i]])-np.log2(peak_df_norm_ratio[channels[reference_ind]])
        
    for i in range(len(channels)):
        f.append(statistics.median(peak_df_norm_ratio["A_"+channels[i]+"_r"].dropna().values.tolist()))
    f[0] = 0
    # exp_f = []
    # for i in range(len(f)):
    #     exp_f.append(np.power(2,f[i]))
    avg_track = np.nanmean(np.array(f))
    scaling_factors = []
    for i in range(len(f)):
        #print(np.power(2,f[i]))
        scaling_factors.append(np.power(2,f[i])/np.power(2,np.nanmean(np.array(f))))

    for i in range(len(channels)):
        peak_df[channels[i]+"_1"] = peak_df[channels[i]]/scaling_factors[i]
        peak_df["log2_"+channels[i]+"_1"] = np.log2(peak_df[channels[i]+"_1"])
    

    for j in range(1,len(channels)):
        peak_df["log2_"+channels[j]+"to"+channels[reference_ind]] = - peak_df["log2_"+channels[reference_ind]+"_1"] + peak_df["log2_"+channels[j]+"_1"]
    
    
    return peak_df

def calc_FC_withratios(wt, ctrl, prot_df, prot_list, prot_df_calc):
    for prot in prot_list:
        fraction_list = list(set(prot_df.loc[(prot_df["Protein"]==prot)]["Fraction"].values.tolist()))
        per_fraction_FC = []
        for frac in fraction_list:
            rvs_T = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)&(prot_df["Fraction"]==frac)]["log2Abundance"].dropna().values.tolist()
            rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)&(prot_df["Fraction"]==frac)]["log2Abundance"].dropna().values.tolist()
            
            for t in rvs_T:
                for c in rvs_C:
                    per_fraction_FC.append(np.power(2,t-c))
        prot_df_calc.loc[prot_df_calc["Protein"]==prot,"FC"] = np.nanmedian(per_fraction_FC)
    
    return prot_df_calc#_dropna

def create_limma_suitable_file(wt, ctrl, prot_df, prot_list, prot_df_calc):
    num_reps = np.max([len(list(set(prot_df.loc[prot_df["Condition"]==wt]["Channel"].values.tolist()))),len(list(set(prot_df.loc[prot_df["Condition"]==ctrl]["Channel"].values.tolist())))])
    cols_rep = []
    for i in range(num_reps):
        cols_rep+=("rep"+str(i+1))
    # limma_prep = pd.DataFrame(columns=["Protein","Condition","Fraction"]+cols_rep)
    
    prot_list = list(set(prot_df["Protein"].values.tolist()))
    dict_list = []
    # conditions = [wt+"-"+ctrl]
    
    for prot in prot_list:
        frac_list = list(set(prot_df.loc[prot_df["Protein"]==prot]["Fraction"].values.tolist()))
        for frac in frac_list:
            rvs_T = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Fraction"]==frac)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
            rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Fraction"]==frac)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
            
            id_T = list(range(0,len(rvs_T)))
            id_C = list(range(0,len(rvs_C)))
            
            if ((len(rvs_T)<2) | (len(rvs_C)<2)):
                continue
            # if prot=="P38323":
            #     print(prot, frac_list, len(rvs_T), len(rvs_C))
            if len(rvs_T) > len(rvs_C):
                df2 = {"Protein":prot,"Fraction":frac}
                pos_C = list(permutations(id_T, len(id_C)))
                for possibility in pos_C:
                    for id_ in range(len(possibility)):
                        df2["rep"+str(id_+1)] = (rvs_T[possibility[id_]]-rvs_C[id_])
                    shallow_copy = {}
                    for k,val in df2.items():
                        shallow_copy[k] = val
                    dict_list.append(shallow_copy)
            elif len(rvs_T) < len(rvs_C):
                df2 = {"Protein":prot,"Fraction":frac}
                pos_T = list(permutations(id_C, len(id_T)))
                for possibility in pos_T:
                    for id_ in range(len(possibility)):
                        df2["rep"+str(id_+1)] = (rvs_T[id_]-rvs_C[possibility[id_]])
                    shallow_copy = {}
                    for k,val in df2.items():
                        shallow_copy[k] = val
                    dict_list.append(shallow_copy)
            elif len(rvs_T) == len(rvs_C):
                df2 = {"Protein":prot,"Fraction":frac}
                pos_T = list(permutations(id_C, len(id_C)))
                for possibility in pos_T:
                    # add_row = [prot,c,frac,np.nan,np.nan,np.nan]
                    for id_ in range(len(possibility)):
                        if id_ < len(id_T):
                            df2["rep"+str(id_+1)] = (rvs_T[id_T[id_]]-rvs_C[id_C[possibility[id_]]])
                    shallow_copy = {}
                    for k,val in df2.items():
                        shallow_copy[k] = val
                    dict_list.append(shallow_copy)
    limma_prep = pd.DataFrame.from_dict(dict_list)
    return limma_prep

def run_r_script(limma_pre_df,output_name):
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    local_rfile_dir= BASE_DIR+'/myapp/'
    r_script_path = local_rfile_dir + 'Magma-limma-intermediary.R'
    input_filename = output_name
    result_filename = '%s_out.txt' % output_name
    print("Doing pval calculation: ", input_filename)

    # rscript_dir = '/Library/Frameworks/R.framework/Resources/bin/'
    rscript_dir = '/usr/bin/'
    # Define the R command to run your script
    r_command = [rscript_dir+'Rscript', r_script_path, input_filename]

    # Run the R script using subprocess
    subprocess.run(r_command, check=True)
    # with open(result_filename, 'wb') as result:
    #     process = subprocess.Popen(r_command,stdout=result)
    #     process.wait()

    limma_post_df = pd.read_csv(result_filename,sep="\t")
    return limma_post_df

def calc_pval_withratios(job_id, wt, ctrl, prot_df, prot_list, prot_df_calc):
    prot_df_calc["pval"] = np.nan
    limma_pre_df = create_limma_suitable_file(wt, ctrl, prot_df, prot_list, prot_df_calc)
    # print(limma_pre_df.loc[limma_pre_df["Protein"]=="P38323"])
    print("Done limma setup")
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files/'+job_id
    filename_lim = local_storage_dir+"/"+wt+"-"+ctrl+"_prelimma"
    limma_pre_df.to_csv(filename_lim+".txt",sep="\t")
    limma_post_df = run_r_script(limma_pre_df,filename_lim)
    # print(limma_post_df.loc[limma_post_df["Protein"]=="P38323"])
    limma_summarized = limma_post_df.groupby(['Protein'])['pvalue','FDR'].mean().reset_index()
    prot_df_calc["pval"] = prot_df_calc.Protein.map(limma_summarized.set_index("Protein")["pvalue"].to_dict())
    prot_df_calc["adjpval"] = prot_df_calc.Protein.map(limma_summarized.set_index("Protein")["FDR"].to_dict())
    return prot_df_calc

def calc_FC_and_pval_limma(job_id, prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm):
    prot_list = list(set(prot_df["Protein"].dropna().values.tolist()))
    
    prot_det = {
        'Protein':prot_list
    }

    prot_median_manual = pd.DataFrame(prot_det)

    prot_median_manual["FC"] = np.nan
    prot_median_manual["pval"] = np.nan
    prot_median_manual["log2FC"] = np.nan
    wt = ""
    ctrl = ""
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print(wt, ctrl)
    print(list(set(prot_df["Condition"].values.tolist())))

    # Calculate FC
    prot_median_manual = calc_FC_withratios(wt, ctrl, prot_df, prot_list, prot_median_manual)
    print("Calculated FC successfully")
    # Calculate p-value using limma (ratios)
    prot_median_manual = calc_pval_withratios(job_id, wt, ctrl, prot_df, prot_list, prot_median_manual)
    print("Calculated pval successfully")

    prot_median_manual["log2FC"] = np.log2(prot_median_manual["FC"])
    print(prot_median_manual.head(2))

    prot_median_manual_dropna = prot_median_manual[["Protein","log2FC","pval"]].dropna()
    print(prot_median_manual_dropna.shape[0])


    prot_median_manual_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(prot_median_manual_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    if ((bait_norm==True) & (bait_name not in prot_median_manual_dropna["Protein"].values.tolist()) & ~(bait_name == "None")):
        print("Bait not added in list after bait normalization")
        prot_median_manual_dropna.loc[len(prot_median_manual_dropna.index)] = [bait_name, 0, 1, 1] 

    prot_median_manual_dropna["Proteinnoiso"] = prot_median_manual_dropna["Protein"].str.split("-",expand=True)[0]
    prot_median_manual_dropna["Gene Symbol"] = prot_median_manual_dropna.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
    #prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]=="ZZZZZ9","Gene Symbol"] = "N"

    if isinstance(init_psm_df, pd.DataFrame):
        forpsm_df = init_psm_df
        prot_list = list(set(prot_median_manual_dropna["Protein"].dropna().values.tolist()))
        prot_median_manual_dropna["# PSMs"] = 0
        for p in prot_list:
            fil = forpsm_df.loc[(forpsm_df["Master Protein Accessions"].str.contains(p, na=False))]
            prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]==p,"# PSMs"] = fil.shape[0]
    return prot_median_manual_dropna

def calc_FC_and_pval(prot_df, init_psm_df, comparison, comp_type, uniprotGene, bait_name, bait_norm):
    prot_list = list(set(prot_df["Protein"].dropna().values.tolist()))
    
    prot_det = {
        'Protein':prot_list
    }

    prot_median_manual = pd.DataFrame(prot_det)

    prot_median_manual["FC"] = np.nan
    prot_median_manual["pval"] = np.nan
    prot_median_manual["log2FC"] = np.nan
    wt = ""
    ctrl = ""
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print(wt, ctrl)
    print(list(set(prot_df["Condition"].values.tolist())))

    prot_median_manual = calc_FC_withratios(wt, ctrl, prot_df, prot_list, prot_median_manual)
    for prot in prot_list:
        rvs_N = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
        
        rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()

        # prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C, alternative='greater').pvalue
        prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C).pvalue
        
        

    prot_median_manual["log2FC"] = np.log2(prot_median_manual["FC"])
    print(prot_median_manual.head(2))
    # print(prot_median_manual_dropna.shape[0])

    prot_median_manual_dropna = prot_median_manual[["Protein","log2FC","pval"]].dropna()

    prot_median_manual_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(prot_median_manual_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    if ((bait_norm==True) & (bait_name not in prot_median_manual_dropna["Protein"].values.tolist()) & ~(bait_name == "None")):
        print("Bait not added in list after bait normalization")
        prot_median_manual_dropna.loc[len(prot_median_manual_dropna.index)] = [bait_name, 0, 1, 1] 

    prot_median_manual_dropna["Proteinnoiso"] = prot_median_manual_dropna["Protein"].str.split("-",expand=True)[0]
    prot_median_manual_dropna["Gene Symbol"] = prot_median_manual_dropna.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
    #prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]=="ZZZZZ9","Gene Symbol"] = "N"

    if isinstance(init_psm_df, pd.DataFrame):
        forpsm_df = init_psm_df
        prot_list = list(set(prot_median_manual_dropna["Protein"].dropna().values.tolist()))
        prot_median_manual_dropna["# PSMs"] = 0
        for p in prot_list:
            fil = forpsm_df.loc[(forpsm_df["Master Protein Accessions"].str.contains(p, na=False))]
            prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]==p,"# PSMs"] = fil.shape[0]
    return prot_median_manual_dropna

def norm_by_bait(bait_name, prot_df, channel_list):
    # TODO: Print warning message that Bait is filtered out if bait no longer in protein level
    if bait_name in prot_df["Protein"].values.tolist():
        prot_df_bait = prot_df.loc[prot_df["Protein"]==bait_name]
        print(prot_df_bait)
        #peak_df_norm_intensity = peak_df_norm_intensity.reset_index()
        channels = channel_list
        fractions = list(set(prot_df["Fraction"].values.tolist()))
        scaling_factors = {}
        all_scaling_factors = []
        for i in range(len(channels)):
            #scaling_factors.append(statistics.median(peak_df_norm_intensity["log2_"+channels[i]+"_1"].dropna().values.tolist()))
            # Taking the avg instead
            #scaling_factors.append(np.average(prot_df_bait["log2_"+channels[i]+"_1"].dropna().values.tolist()))
            scaling_factors[i]=[]
            for j in range(len(fractions)):
                sf = np.average(prot_df_bait.loc[(prot_df_bait["Channel"]==channels[i])&(prot_df_bait["Fraction"]==fractions[j])]["log2Abundance"].dropna().values.tolist())
                scaling_factors[i].append(sf)
                if math.isnan(sf):
                    continue
                else:
                    all_scaling_factors.append(sf)
        
        global_median = statistics.median(all_scaling_factors)
        print(scaling_factors)
        #for i in range(len(scaling_factors)):
            #scaling_factors[i] = - scaling_factors[i] + global_median
        print(global_median)
        for i in range(len(channels)):
            for j in range(len(fractions)):
                scaling_factors[i][j] = - scaling_factors[i][j] + global_median
        print(scaling_factors)
        for i in range(len(channels)):
            #prot_df["log2_"+channels[i]+"_1"] = peak_df["log2_"+channels[i]+"_1"]+scaling_factors[i]
            #peak_df[channels[i]+"_1"] = np.power(2,peak_df["log2_"+channels[i]+"_1"])
            for j in range(len(fractions)):
                
                if math.isnan(scaling_factors[i][j]):
                    continue
                else:
                    prot_df.loc[((prot_df["Channel"]==channels[i])&(prot_df["Fraction"]==fractions[j])),"log2Abundance"] += scaling_factors[i][j]
    else:
        print("Bait no longer in filtered Protein level file")
    
    return prot_df

###########################################################################################
############ FASTER VERSION - keep this or remove need for long format ####################
######## This is useful if need linear model with covars added eventually #################
###########################################################################################
def expansion_to_fractionXchannel_faster(peak_df,annotation_df, comparison, comp_type):
    print(annotation_df.head(2))
    if comp_type==1:
        wt = comparison.split("-")[0]
        ctrl = comparison.split("-")[1]
    else:
        wt = comparison.split("-")[1]
        ctrl = comparison.split("-")[0]

    print("Reducing file size with long form expansion: ",wt, ctrl)
    channel_list = list(set(annotation_df.loc[((annotation_df["Label"]==wt)|(annotation_df["Label"]==ctrl))]["Channel"].values.tolist()))
    # channel_list = sort1(channel_list)
    print("Reduced channel list for long form expansion: ",channel_list)
    peak_df[["AnnotatedSeq","Charge","Fraction"]] = peak_df["PeakID"].str.split("_",expand=True)
    # Do the same as above but for all peptides belonging to proteins
    prot_list = list(set(peak_df["Protein"].dropna().values.tolist()))
    master_prots_yulab = list(set(peak_df["Protein"].dropna().values.tolist()))
    prot_det = {
        'Protein':prot_list
    }
    # TODO - Make it so only channels present in annotation file are used below
    # prot_df = pd.DataFrame(columns=["Protein","log2Abundance","Channel","Condition","Fraction"])#peak_df_norm_ratio.copy()
    # summarize by median
    dict_list = []
    for i in range(len(prot_list)):
        prot = prot_list[i]
        fractions = list(set(peak_df["Fraction"].values.tolist()))
        for frac in fractions:
            df = peak_df.loc[(peak_df["Protein"]==prot)&(peak_df["Fraction"]==frac)]
            for channel in channel_list:
                if len(df["log2_"+str(channel)+"_1"].dropna().values.tolist())>0:
                    dict_list.append({"Protein":prot,"log2Abundance":statistics.median(df["log2_"+str(channel)+"_1"].dropna().values.tolist()),"Channel":str(channel),"Condition":annotation_df.loc[annotation_df["Channel"]==channel]["Label"].values.tolist()[0],"Fraction":frac})
                    #prot_df = pd.concat([prot_df,pd.DataFrame([[prot,statistics.median(df["log2_"+channel+"_1"].dropna().values.tolist()),channel,annotation_df.loc[annotation_df["Channel"]==channel]["Label"].values.tolist()[0],frac]],columns=list(prot_df))],ignore_index=True,sort=False)
    prot_df = pd.DataFrame.from_dict(dict_list)
    return prot_df

############ SLOWER VERSION - FOR MULITPLE COMPARISON QUEUE #################
def expansion_to_fractionXchannel(peak_df,annotation_df, channel_list):
    print(annotation_df.head(2))
    peak_df[["AnnotatedSeq","Charge","Fraction"]] = peak_df["PeakID"].str.split("_",expand=True)
    # Do the same as above but for all peptides belonging to proteins
    prot_list = list(set(peak_df["Protein"].dropna().values.tolist()))
    master_prots_yulab = list(set(peak_df["Protein"].dropna().values.tolist()))
    prot_det = {
        'Protein':prot_list
    }
    # TODO - Make it so only channels present in annotation file are used below
    # prot_df = pd.DataFrame(columns=["Protein","log2Abundance","Channel","Condition","Fraction"])#peak_df_norm_ratio.copy()
    # summarize by median
    dict_list = []
    for i in range(len(prot_list)):
        prot = prot_list[i]
        fractions = list(set(peak_df["Fraction"].values.tolist()))
        for frac in fractions:
            df = peak_df.loc[(peak_df["Protein"]==prot)&(peak_df["Fraction"]==frac)]
            for channel in channel_list:
                if len(df["log2_"+channel+"_1"].dropna().values.tolist())>0:
                    dict_list.append({"Protein":prot,"log2Abundance":statistics.median(df["log2_"+channel+"_1"].dropna().values.tolist()),"Channel":channel,"Condition":annotation_df.loc[annotation_df["Channel"]==channel]["Label"].values.tolist()[0],"Fraction":frac})
                    #prot_df = pd.concat([prot_df,pd.DataFrame([[prot,statistics.median(df["log2_"+channel+"_1"].dropna().values.tolist()),channel,annotation_df.loc[annotation_df["Channel"]==channel]["Label"].values.tolist()[0],frac]],columns=list(prot_df))],ignore_index=True,sort=False)
    prot_df = pd.DataFrame.from_dict(dict_list)
    return prot_df


def normalization_median_intensity(peak_df,channel_list):
    for channel in channel_list:
        peak_df[channel+"_r"] = np.power(2,peak_df["log2_"+channel+"_1"])
        peak_df[channel] = peak_df[channel+"_r"]/peak_df[channel+"_r"].sum()

    peak_df_norm_ratio = peak_df.copy()

    f=[]
    reference_ind=0
    channels = channel_list#["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
    f=[]
    for i in range(len(channels)):
        peak_df_norm_ratio["A_"+channels[i]+"_r"] = np.log2(peak_df_norm_ratio[channels[i]])-np.log2(peak_df_norm_ratio[channels[reference_ind]])
        
    for i in range(len(channels)):
        f.append(statistics.median(peak_df_norm_ratio["A_"+channels[i]+"_r"].dropna().values.tolist()))
    f[0] = 0
    
    avg_track = np.nanmean(np.array(f))
    scaling_factors = []
    for i in range(len(f)):
        #print(np.power(2,f[i]))
        scaling_factors.append(np.power(2,f[i])/np.power(2,np.nanmean(np.array(f))))

    for i in range(len(channels)):
        peak_df_norm_ratio[channels[i]+"_1"] = peak_df_norm_ratio[channels[i]]/scaling_factors[i]
        peak_df_norm_ratio["log2_"+channels[i]+"_1"] = np.log2(peak_df_norm_ratio[channels[i]+"_1"])
    

    for j in range(1,len(channels)):
        peak_df_norm_ratio["log2_"+channels[j]+"to"+channels[reference_ind]] = - peak_df_norm_ratio["log2_"+channels[reference_ind]+"_1"] + peak_df_norm_ratio["log2_"+channels[j]+"_1"]
    
    return peak_df_norm_ratio

# def norm_by_bait(bait_name, peak_df,channel_list):
#     #if num_channels==10:

#     for channel in channel_list:
#         peak_df[channel+"_r"] = np.power(2,peak_df["log2_"+channel+"_1"])
#         peak_df[channel] = peak_df[channel+"_r"]/peak_df[channel+"_r"].sum()

#     peak_df_norm_ratio = peak_df.loc[peak_df["Protein"]==bait_name].copy()

#     f=[]
#     reference_ind=0
#     channels = channel_list#["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
#     f=[]
#     for i in range(len(channels)):
#         peak_df_norm_ratio["A_"+channels[i]+"_r"] = np.log2(peak_df_norm_ratio[channels[i]])-np.log2(peak_df_norm_ratio[channels[reference_ind]])
        
#     for i in range(len(channels)):
#         f.append(statistics.median(peak_df_norm_ratio["A_"+channels[i]+"_r"].dropna().values.tolist()))
#     f[0] = 0
    
#     avg_track = np.nanmean(np.array(f))
#     scaling_factors = []
#     for i in range(len(f)):
#         #print(np.power(2,f[i]))
#         scaling_factors.append(np.power(2,f[i])/np.power(2,np.nanmean(np.array(f))))

#     peak_df_norm_ratio = peak_df.copy()

#     for i in range(len(channels)):
#         peak_df_norm_ratio[channels[i]+"_1"] = peak_df_norm_ratio[channels[i]]/scaling_factors[i]
#         peak_df_norm_ratio["log2_"+channels[i]+"_1"] = np.log2(peak_df_norm_ratio[channels[i]+"_1"])
    

#     for j in range(1,len(channels)):
#         peak_df_norm_ratio["log2_"+channels[j]+"to"+channels[reference_ind]] = - peak_df_norm_ratio["log2_"+channels[reference_ind]+"_1"] + peak_df_norm_ratio["log2_"+channels[j]+"_1"]
    
#     return peak_df_norm_ratio


def psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, channel_list, iso_filter=100):
    # Filtered PSMs using non-shared PSMs, 
    psm_df = compile_sn_raw_int(sn_df,intensity_df,channel_list)
    cols_abund = []
    for channel in channel_list:
        cols_abund.append("Abundance: "+channel)

    psm_copy = psm_df#.copy() #Here the copy changes the end results - affects when changing values of abundance with norm factor

    
    psm_df.loc[:,"#Measurements"] = psm_copy[cols_abund].shape[1] - psm_copy[cols_abund].isnull().sum(axis=1)
    psm_df.loc[:,"weights"] = np.log(psm_copy["Intensity"])/np.log(1.2)

    # Check if user wants to do row-normalization
    if row_norm_user == True:
        psm_df.loc[:,"Norm-factor"] = psm_copy[cols_abund].mean(axis=1)
        for c in cols_abund:
            psm_df.loc[:,c] = psm_df[c]/psm_df["Norm-factor"]
    

    cols_reporter = []
    for c in cols_abund:
        psm_df.loc[:,"log2_"+c.split(": ")[1]+"_1"] = np.log2(psm_copy[c])
        cols_reporter.append("log2_"+c.split(": ")[1]+"_1")
    
    psm_df = psm_df.loc[psm_df["#Measurements"]>2]

    # Remove all PSMs with IS filter II<=50, SPS%>=65
    psm_df = psm_df.loc[(psm_df["Isolation Interference [%]"]<=iso_filter)]#&(psm_df["SPS Mass Matches [%]"]>=65)]

    #Only using master protein unique peptides for fair comparison
    psm_df = psm_df.loc[psm_df["# Protein Groups"]==1]

    psm_copy = psm_df#.copy()

    #Check if fraction separation is problem
    # psm_df["PSMID"] = psm_df["Annotated Sequence"] + "_" + psm_df["Charge"].astype(str) + "_1.raw"
    psm_df.loc[:,"PSMID"] = psm_copy["Annotated Sequence"] + "_" + psm_copy["Charge"].astype(str) + "_"+ psm_copy["Spectrum File"].str.split("_").str[-1]


    print("# Proteins: ",len(list(set(psm_df["Master Protein Accessions"].values.tolist()))))
    print("# Peaks: ",len(list(set(psm_df["PSMID"].values.tolist()))))

    #################################################################################################
    # Step to filter out PSMs and for same peptide found multiple times in one fraction - weighted avg
    #################################################################################################
    peaks_unique = list(set(psm_df["PSMID"].values.tolist()))
    peak_data = {
    'PeakID': peaks_unique
    }

    peak_df = pd.DataFrame(peak_data)
    peak_df["Protein"] = peak_df.PeakID.map(psm_df.set_index("PSMID")["Master Protein Accessions"].to_dict())


    for i in range(peak_df.shape[0]):
        peak = peak_df.iloc[i]["PeakID"]#.copy()
        psms_assoc = psm_df.loc[psm_df["PSMID"]==peak]#.copy()
        if psms_assoc.shape[0]>1:
            for c in cols_reporter:
                psms_c = psms_assoc[[c, "weights"]].dropna()
                weights = np.array(psms_c["weights"].values.tolist())
                values = np.array(psms_c[c].values.tolist())
                if np.nansum(weights)>0:
                    peak_df.loc[peak_df["PeakID"]==peak,c] = np.nansum(values*weights)/(np.nansum(weights))
        else:
            for c in cols_reporter:
                peak_df.loc[peak_df["PeakID"]==peak,c] = psms_assoc[c].values.tolist()[0]
    
    return peak_df

def compile_sn_raw_int(sn_df,raw_df,channel_list):

    raw_copy = raw_df.copy()
    sn_copy = sn_df.copy()
    raw_df.loc[:,"Features"] = raw_copy["Annotated Sequence"]+"_"+raw_copy["Charge"].astype(str)+"_"+raw_copy["Spectrum File"] +"_" + raw_copy["RT [min]"].astype(str)
    sn_df.loc[:,"Features"] = sn_copy["Annotated Sequence"]+"_"+sn_copy["Charge"].astype(str)+"_"+sn_copy["Spectrum File"]+"_" + sn_copy["RT [min]"].astype(str)
    
    tmp_df = raw_df.copy()
    tmp_df_filtered = tmp_df.copy()
    cols_noise = []
    for channel in channel_list:
        tmp_df["SN-"+channel] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: "+channel].to_dict())
        tmp_df["Noise-"+channel] = tmp_df["Abundance: "+channel]/tmp_df["SN-"+channel]
    tmp_df['ge_1.5'] = tmp_df.filter(like='SN-').max(axis=1).ge(1.5).astype(float)
    tmp_df_filtered = tmp_df.loc[tmp_df['ge_1.5']==1.0].copy()
    return tmp_df_filtered



# def get_x0(final_prot_df):
#     enrich_fact = final_prot_df["log2FC"].dropna().values.tolist()
#     # ax = sns.distplot(enrich_fact, hist=True, kde=True,
#     #             bins=int(180/5), color = 'darkblue',
#     #             hist_kws={'edgecolor':'black'},
#     #             kde_kws={'linewidth': 4})
#     # ax.set(xlabel='enrichment factor', ylabel='Density')
#     # m_sd = np.mean(enrich_fact)+np.std(enrich_fact)
#     mx0 = np.std(enrich_fact)
#     #plt.plot([m_sd,m_sd], [0, 0.6])
#     #plt.show()
#     print(np.std(enrich_fact))
#     print(2*np.std(enrich_fact))
#     return mx0

def pval_thresh_find(df,fdr_thresh):
    pvals_list = df.loc[df["adjpval"]<=fdr_thresh]["pval"].values.tolist()
    print("# pvalue where you get FDR<=",fdr_thresh,": ", len(pvals_list))
    pval_threshold = np.Inf
    if (df.loc[df["adjpval"]<=fdr_thresh].shape[0]>0):
        pval_threshold = -1.0*np.log10(sorted(df.loc[df["adjpval"]<=fdr_thresh]["pval"].values.tolist())[-1])#-1.0*np.log10(0.05)#FDR <0.1
        print("Max pvalue where you get FDR<=",fdr_thresh,": ", sorted(df.loc[df["adjpval"]<=fdr_thresh]["pval"])[-1])
    return pval_threshold

def generate_baseline_volcano_plot(final_prot_df, fdr_thresh, bait_uniprot, comparison, labeltype="one"):
    psm_cutoff = final_prot_df.iloc[0]["PSM Cutoff"]
    plt.switch_backend('Agg') 
    fig,ax = plt.subplots(1,1)#,figsize=(22,9))
    m_bait_df = final_prot_df
    m_bait_df["neglogpval"] = -np.log10(m_bait_df["pval"])#.dropna()
    
    rem_genes=["KRT1","KRT10","KRT13","KRT14","KRT15","KRT16","KRT17","KRT18","KRT19","KRT2","KRT24","KRT25","KRT26",
    "KRT27","KRT28","KRT3","KRT4","KRT5","KRT6A","KRT6B","KRT6C","KRT75","KRT76","KRT77","KRT79","KRT8","KRT9","KRT36","KRT23"]
    #psm_cutoff_bait
    df = m_bait_df.loc[(m_bait_df["# PSMs"]>=psm_cutoff)&~(m_bait_df["Gene Symbol"].isin(rem_genes))]
    FC_threshold = np.log2(2)
    pval_threshold = -1.0*np.log10(0.05)

    ## Have have adjusted p-vals
    pval_threshold = pval_thresh_find(df, fdr_thresh)
    if (pval_threshold<-np.log10(0.05))|(np.power(10,-pval_threshold)==0):
        pval_threshold = -np.log10(0.05)
    print(pval_threshold)
    print(np.power(10,-pval_threshold))

    # INITIAL PLOT VISUALIZE
    mdf = df.loc[(df["log2FC"]>=-1.0*FC_threshold) & (df["log2FC"]<FC_threshold) & (df["neglogpval"]>=pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#BDC3C7',marker='o',edgecolors='#000000')#'#6C7A89'-DARK GREY, '#BDC3C7'-LIGHT GREY
    mdf = df.loc[(df["log2FC"]>=FC_threshold) & (df["neglogpval"]>=pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#1E824C',marker='o',edgecolors='#000000')
    mdf = df.loc[(df["log2FC"]>=FC_threshold) & (df["neglogpval"]<pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#87D37C',marker='o',edgecolors='#000000')
    mdf = df.loc[(df["log2FC"]>=-1.0*FC_threshold) & (df["log2FC"]<FC_threshold) & (df["neglogpval"]<pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#BDC3C7',marker='o',edgecolors='#000000')
    mdf = df.loc[(df["log2FC"]<-1.0*FC_threshold) & (df["neglogpval"]>=pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#980000',marker='o',edgecolors='#000000')
    mdf = df.loc[(df["log2FC"]<-1.0*FC_threshold) & (df["neglogpval"]<pval_threshold)]
    ax.scatter(mdf["log2FC"].values.tolist(),mdf["neglogpval"].values.tolist(),s=80,c='#FFB6C1',marker='o',edgecolors='#000000')

    # PLOTS LABELLED
    # TODO- Add the HINT functionality of labelling the known interactions passing baseline
    m_bait_fil = df.loc[(df["log2FC"]>=FC_threshold) & (df["neglogpval"]>=pval_threshold) & ~(df["Protein"]==bait_uniprot)]
    # # Only add labels for top 50 points
    # m_bait_fil = m_bait_fil.sort_values(["log2FC","neglogpval"],ascending = False).head(50)
    # repel_labels(ax, m_bait_fil, k=0.0025)
    
    for index, rowv in m_bait_fil.iterrows():
        ax.text(rowv["log2FC"]+0.01,rowv["neglogpval"]+0.02,rowv["Gene Symbol"],c='#7F7F7F',fontsize=15)
    
    if labeltype=="both":
        m_bait_fil = df.loc[(df["log2FC"]<=-1.0*FC_threshold) & (df["neglogpval"]>=pval_threshold) & ~(df["Protein"]==bait_uniprot)]
        for index, rowv in m_bait_fil.iterrows():
            ax.text(rowv["log2FC"]+0.01,rowv["neglogpval"]+0.02,rowv["Gene Symbol"],c='#7F7F7F',fontsize=15)


    # Label the bait (IP-MS) 
    m_bait_fil = df.loc[(df["Protein"]==bait_uniprot)]
    for index, rowv in m_bait_fil.iterrows():
        ax.text(rowv["log2FC"]+0.01,rowv["neglogpval"]+0.02,rowv["Gene Symbol"],c='#7F7F7F',fontsize=20)

    # Add horizontal lines for FDR of 10%, 5%, and 1%
    # If else to ensure the script does not pause if no FC, pvalues were able to be calculated
    if df.shape[0]>0:
        ax.axhline(pval_thresh_find(df, 0.01),min(df["log2FC"].values.tolist()),max(df["log2FC"].values.tolist()),c='#de2d26',lw=1,ls='--')
        ax.axhline(pval_thresh_find(df, 0.05),min(df["log2FC"].values.tolist()),max(df["log2FC"].values.tolist()),c='#fc9272',lw=1,ls='--')
        ax.axhline(pval_thresh_find(df, 0.1),min(df["log2FC"].values.tolist()),max(df["log2FC"].values.tolist()),c='#fee0d2',lw=1,ls='--')
        ax.axvline(0,0,max(df["neglogpval"].values.tolist()),c='#000000',lw=1)
        ax.axhline(0,min(df["log2FC"].values.tolist()),max(df["log2FC"].values.tolist()),c='#000000',lw=1)
    else:
        ax.axvline(0,0,1,c='#000000',lw=1)
        ax.axhline(0,-1,1,c='#000000',lw=1)

    # for rowx,rowy in zip(x_outer_right,y_outer_right):
    #     ax.text(rowx+0.01,rowy+0.02,df.loc[(df["log2FC"]==rowx)&(df["neglogpval"]==rowy)]["Gene Symbol"].values.tolist()[0],c='#7F7F7F',fontsize=15)
    
    # # Outer right line
    # c=5.5#0.6
    # x0=get_x0(final_prot_df)#+0.03
    # xval = np.linspace(6.8,17,60) #These points need to be determined per plot as well
    # yval = []
    # for i in xval:
    #     yval.append(c/(i-2*x0))
    # ax.plot(xval,yval,c='#800000')
    #a,b = FilterWithHiNT(prot_to_uniprot[bait])#HQ
    #print(b)

    
    ax.grid(True, alpha=0.3)
    start, end = ax.get_xlim()
    axis_thres = max(abs(start),abs(end))
    ax.xaxis.set_ticks(np.arange(-axis_thres, axis_thres, 1.0))
    ax.tick_params(labelcolor='#7F7F7F',labelsize=22)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_xlabel("$log_2(FC)$",c='#7F7F7F',fontsize=24)
    ax.set_ylabel("$-log_{10}(pval)$",c='#7F7F7F',fontsize=24)

    tmp_comparisons = comparison.split("-vs-")
    display_comparison = comparison
    if len(tmp_comparisons)==2:
        display_comparison = tmp_comparisons[0] + "/" + tmp_comparisons[1]
    ax.set_title(display_comparison,c='#7F7F7F',fontsize=26)
    plt.box(on=None)
    #plt.savefig(comparison+"_FC_and_pval_channelxfraction_OPTIMIZED.png")
    #plt.savefig(comparison+"_FC_and_pval_FULL.png")
    return plt


# def repel_labels(ax, df, k=0.01):
#     G = nx.DiGraph()
#     data_nodes = []
#     init_pos = {}
#     # for xi, yi, label in zip(x, y, labels):
#     for index, rowv in df.iterrows():
#         # ax.text(rowv["log2FC"]+0.01,rowv["neglogpval"]+0.02,rowv["Gene Symbol"],c='#7F7F7F',fontsize=15)
#         xi = rowv["log2FC"]
#         yi = rowv["neglogpval"]
#         label = rowv["Gene Symbol"]
#         data_str = 'data_{0}'.format(label)
#         G.add_node(data_str)
#         G.add_node(label)
#         G.add_edge(label, data_str)
#         data_nodes.append(data_str)
#         init_pos[data_str] = (xi, yi)
#         init_pos[label] = (xi, yi)

#     pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

#     # undo spring_layout's rescaling
#     pos_after = np.vstack([pos[d] for d in data_nodes])
#     pos_before = np.vstack([init_pos[d] for d in data_nodes])
#     scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
#     scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
#     shift = np.array([shift_x, shift_y])
#     for key, val in pos.items():
#         pos[key] = (val*scale) + shift

#     for label, data_str in G.edges():
#         ax.annotate(label,
#                     xy=pos[data_str], xycoords='data',
#                     xytext=pos[label], textcoords='data',
#                     arrowprops=dict(arrowstyle="->",
#                                     shrinkA=0, shrinkB=0,
#                                     connectionstyle="arc3", 
#                                     color='red'), )
#     # expand limits
#     all_pos = np.vstack(pos.values())
#     x_span, y_span = np.ptp(all_pos, axis=0)
#     mins = np.min(all_pos-x_span*0.15, 0)
#     maxs = np.max(all_pos+y_span*0.15, 0)
#     ax.set_xlim([mins[0], maxs[0]])
#     ax.set_ylim([mins[1], maxs[1]])