from django.db import models
import pandas as pd
import os

#####################Shagun script dependencies#################################
import numpy as np
#import pandas as pd
#import os
import statistics
from scipy import stats
from scipy.stats import pearsonr
import statsmodels
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp



class Document(models.Model):
    #docfile = models.FileField(upload_to='documents/%Y/%m/%d')
    docfile = models.FileField(upload_to='documents/tmtipms_files')

# class select_choice(models.Model):
# 	choice = models.FileField(upload_to='documents/tmtipms_files')
   
##### views.py function removed ################
def networkplot(request):
    message_network = ''
    file_annotation_network=''
    validation_network=''
    documents=''
    # Handle file upload
    result_file_name=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'#BASE_DIR+"/Results/"
    if request.method == 'POST':
        form_networkinput = DocumentForm_NetworkPlot(request.POST, request.FILES)

        if form_networkinput.is_valid():
            newdoc_networkinput = Document(docfile=request.FILES['docfile_networkcsv'])   
            newdoc_networkinput.save()
            file_name_networkcsv=request.FILES['docfile_networkcsv'].name
            print ('The third file name is %s\n'%(file_name_networkcsv))
            print ('The third file name is %s\n'%(request.FILES['docfile_networkcsv']))

            #result_df = pd.read_csv("%s/%s"%(local_storage_dir,file_name_1),sep=",")
            #result_df.to_csv(BASE_DIR+"/Results/"+file_name_1,sep=",")
            

            column_list = ["proteinA", "proteinB", "GeneSymbolA", "GeneSymbolB","numPSMsA"]
            file_network_df = pd.read_csv('%s/%s'%(local_storage_dir, file_name_networkcsv), usecols=column_list)
            
            if ((column_list[0] in list(file_network_df))&(column_list[1] in list(file_network_df))&(column_list[2] in list(file_network_df))&(column_list[3] in list(file_network_df))):
                file_annotation_network='SUCCESS: The uploaded annotation file "%s" consists of both columns "proteinA","proteinB","GeneSymbolA" and "GeneSymbolB" in'%(file_name_networkcsv)
            else:
                file_annotation_network='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of columns "proteinA","proteinB","GeneSymbolA" and "GeneSymbolB" in '%(file_name_networkcsv)

            file_network_df["protAnoiso"] = file_network_df["proteinA"].str.split("-",expand=True)[0]
            file_network_df["protBnoiso"] = file_network_df["proteinB"].str.split("-",expand=True)[0]
            # file_network_df["groupA"] = file_network_df["protAnoiso"]
            # file_network_df["groupB"] = file_network_df["protBnoiso"]
            file_network_df["groupA"] = 0#file_network_df["protAnoiso"]
            file_network_df["groupB"] = 0#file_network_df["protBnoiso"]

            file_network_df["complexAname"] = ""#file_network_df["protAnoiso"]
            file_network_df["complexBname"] = ""#file_network_df["protBnoiso"]
            # Filter by cell line
            # uniprotComplexes = uniprotComplexes.loc[uniprotComplexes["Cell line"].isin(["different human tissues and cancer cell lines","HeLa cells","293 cells","skin","skeletal muscle cells"])]
            # groupA
            prot_list = list(set(list(set(file_network_df["proteinA"].values.tolist()))+list(set(file_network_df["proteinB"].values.tolist()))))

            # print(list(uniprotComplexes))
            for prot in prot_list:
                complexenames = ','.join(map(str,uniprotComplexes.loc[uniprotComplexes["subunits(UniProt IDs)"].str.contains(prot)]["ComplexName"].values.tolist()))
                complexes = uniprotComplexes.loc[uniprotComplexes["subunits(UniProt IDs)"].str.contains(prot)]["ComplexID"].values.tolist()
                if len(complexes)>0:
                    # randomly pick the first complex match
                    file_network_df.loc[file_network_df["proteinA"]==prot,"groupA"] = complexes[0]
                    file_network_df.loc[file_network_df["proteinB"]==prot,"groupB"] = complexes[0]

                    file_network_df.loc[file_network_df["proteinA"]==prot,"complexAname"] = complexenames#uniprotComplexes.loc[uniprotComplexes["ComplexID"]==complexes[0]]["ComplexName"].values.tolist()[0]
                    file_network_df.loc[file_network_df["proteinB"]==prot,"complexBname"] = complexenames#uniprotComplexes.loc[uniprotComplexes["ComplexID"]==complexes[0]]["ComplexName"].values.tolist()[0]

            

            file_network_df[column_list+["groupA","groupB","complexAname","complexBname"]].to_csv('%s/%s'%(local_storage_dir, file_name_networkcsv))
            result_file_name=file_name_networkcsv
            # file_network_df["groupA"] = file_network_df.Master_protein_noiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())


    else:
        form_networkinput = DocumentForm_NetworkPlot()  # An empty, unbound form
    
    context = {'form_networkinput': form_networkinput, 'message_network': message_network, 'file_annotation_network': file_annotation_network, 'validation_5': validation_network, 'result_file_name':result_file_name}

    return render(request,'networkplot.html', context)



###########Custom functions of Yugandhar Kumar########################
def validate_bait(bait_name, file_path):
	file_df=pd.read_csv(file_path, sep='\t', header=0)
	bait_validation=False
	if bait_name in file_df['Master Protein Accessions'].values:
		bait_validation=True
	
	return bait_validation



############Custom functions of Shagun Gupta#############################

def YuLab_fractionxchannel(sn_df, intensity_df, annotation_df, bait_name, comparison, comp_type, row_norm_user,column_norm_user,num_channels, uniprotGene):
    peak_df = psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, num_channels)
    print("Done PSM summarization", peak_df.head(2))
    init_psm_df = compile_sn_raw_int(sn_df,intensity_df,num_channels)
    if column_norm_user==True:
        peak_df = normalization_median_intensity(peak_df,num_channels)
        print("Done median ratio based normalization")
    
    prot_df = expansion_to_fractionXchannel(peak_df,annotation_df)
    print("Done setup for FC and pvalue calculation", prot_df.shape[0])
    #print (type(comp_type), comp_type, comparison.split("-"))
    final_prot_df = calc_FC_and_pval(prot_df, init_psm_df, comparison, comp_type, uniprotGene)
    psm_cutoff = psm_cutoff_suggestion(bait_name, final_prot_df)
    final_prot_df.at[0,"PSM Cutoff"] = psm_cutoff
    print("Done calculating FC and p-value")
    return final_prot_df

def psm_cutoff_suggestion(bait_name, final_prot_df):
    psm_cutoff = 5 # The usual cutoff
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

def calc_FC_and_pval(prot_df, init_psm_df, comparison, comp_type, uniprotGene):
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

    for prot in prot_list:
        rvs_N = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==wt)]["log2Abundance"].dropna().values.tolist()
        
        rvs_C = prot_df.loc[(prot_df["Protein"]==prot)&(prot_df["Condition"]==ctrl)]["log2Abundance"].dropna().values.tolist()
        #print (prot, rvs_N, rvs_C)
        
        
        all_log2fc = []
        for n in rvs_N:
            for c in rvs_C:
                all_log2fc.append(n-c)
        all_fc = [np.power(2,x) for x in all_log2fc]
        
        # If taking avg of ratios
        prot_median_manual.loc[prot_median_manual["Protein"]==prot,"FC"] = np.nanmean(all_fc)
        #print (prot, np.nanmean(all_fc))
        prot_median_manual.loc[prot_median_manual["Protein"]==prot,"pval"] = stats.ttest_ind(rvs_N, rvs_C).pvalue
        #print (prot, stats.ttest_ind(rvs_N, rvs_C).pvalue)
        

    prot_median_manual["log2FC"] = np.log2(prot_median_manual["FC"])

    prot_median_manual_dropna = prot_median_manual[["Protein","log2FC","pval"]].dropna()
    print(prot_median_manual_dropna.shape[0])

    prot_median_manual_dropna["adjpval"] = statsmodels.stats.multitest.multipletests(prot_median_manual_dropna["pval"].values.tolist(),method="fdr_bh")[1]

    prot_median_manual_dropna["Proteinnoiso"] = prot_median_manual_dropna["Protein"].str.split("-",expand=True)[0]
    prot_median_manual_dropna["Gene Symbol"] = prot_median_manual_dropna.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
    #prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]=="ZZZZZ9","Gene Symbol"] = "N"

    forpsm_df = init_psm_df
    prot_list = list(set(prot_median_manual_dropna["Protein"].dropna().values.tolist()))
    prot_median_manual_dropna["# PSMs"] = 0
    for p in prot_list:
        fil = forpsm_df.loc[(forpsm_df["Master Protein Accessions"].str.contains(p, na=False))]
        prot_median_manual_dropna.loc[prot_median_manual_dropna["Protein"]==p,"# PSMs"] = fil.shape[0]
    return prot_median_manual_dropna

def expansion_to_fractionXchannel(peak_df,annotation_df):
    peak_df[["AnnotatedSeq","Charge","Fraction"]] = peak_df["PeakID"].str.split("_",expand=True)
    # Do the same as above but for all peptides belonging to proteins
    prot_list = list(set(peak_df["Protein"].dropna().values.tolist()))
    master_prots_yulab = list(set(peak_df["Protein"].dropna().values.tolist()))
    prot_det = {
        'Protein':prot_list
    }
    # TODO - Make it so only channels present in annotation file are used below
    prot_df = pd.DataFrame(columns=["Protein","log2Abundance","Channel","Condition","Fraction"])#peak_df_norm_ratio.copy()
    # summarize by median
    for i in range(len(prot_list)):
        prot = prot_list[i]
        fractions = list(set(peak_df["Fraction"].values.tolist()))
        for frac in fractions:
            df = peak_df.loc[(peak_df["Protein"]==prot)&(peak_df["Fraction"]==frac)]
            if len(df["log2_126_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_126_1"].dropna().values.tolist()),"126",annotation_df.loc[annotation_df["Channel"]=="126"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
                #print ("126", annotation_df.loc[annotation_df["Channel"]=="126"]["Label"].values.tolist()[0], statistics.median(df["log2_126_1"].dropna().values.tolist()))
            if len(df["log2_127N_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_127N_1"].dropna().values.tolist()),"127N",annotation_df.loc[annotation_df["Channel"]=="127N"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_127C_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_127C_1"].dropna().values.tolist()),"127C",annotation_df.loc[annotation_df["Channel"]=="127C"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_128N_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_128N_1"].dropna().values.tolist()),"128N",annotation_df.loc[annotation_df["Channel"]=="128N"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_128C_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_128C_1"].dropna().values.tolist()),"128C",annotation_df.loc[annotation_df["Channel"]=="128C"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            
            if len(df["log2_129N_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_129N_1"].dropna().values.tolist()),"129N",annotation_df.loc[annotation_df["Channel"]=="129N"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_129C_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_129C_1"].dropna().values.tolist()),"129C",annotation_df.loc[annotation_df["Channel"]=="129C"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_130N_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_130N_1"].dropna().values.tolist()),"130N",annotation_df.loc[annotation_df["Channel"]=="130N"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_130C_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_130C_1"].dropna().values.tolist()),"130C",annotation_df.loc[annotation_df["Channel"]=="130C"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
            if len(df["log2_131_1"].dropna().values.tolist())>0:
                prot_df = prot_df.append(pd.Series([prot,statistics.median(df["log2_131_1"].dropna().values.tolist()),"131",annotation_df.loc[annotation_df["Channel"]=="131"]["Label"].values.tolist()[0],frac], index =prot_df.columns),ignore_index=True,sort=False)
    return prot_df

def normalization_median_intensity(peak_df,num_channels):
    #if num_channels==10:
    peak_df["126_r"] = np.power(2,peak_df["log2_126_1"])
    peak_df["127N_r"] = np.power(2,peak_df["log2_127N_1"])
    peak_df["127C_r"] = np.power(2,peak_df["log2_127C_1"])
    peak_df["128N_r"] = np.power(2,peak_df["log2_128N_1"])
    peak_df["128C_r"] = np.power(2,peak_df["log2_128C_1"])
    peak_df["129N_r"] = np.power(2,peak_df["log2_129N_1"])
    peak_df["129C_r"] = np.power(2,peak_df["log2_129C_1"])
    peak_df["130N_r"] = np.power(2,peak_df["log2_130N_1"])
    peak_df["130C_r"] = np.power(2,peak_df["log2_130C_1"])
    peak_df["131_r"] = np.power(2,peak_df["log2_131_1"])

    peak_df["126"] = peak_df["126_r"]/peak_df["126_r"].sum()
    peak_df["127N"] = peak_df["127N_r"]/peak_df["127N_r"].sum()
    peak_df["127C"] = peak_df["127C_r"]/peak_df["127C_r"].sum()
    peak_df["128N"] = peak_df["128N_r"]/peak_df["128N_r"].sum()
    peak_df["128C"] = peak_df["128C_r"]/peak_df["128C_r"].sum()
    peak_df["129N"] = peak_df["129N_r"]/peak_df["129N_r"].sum()
    peak_df["129C"] = peak_df["129C_r"]/peak_df["129C_r"].sum()
    peak_df["130N"] = peak_df["130N_r"]/peak_df["130N_r"].sum()
    peak_df["130C"] = peak_df["130C_r"]/peak_df["130C_r"].sum()
    peak_df["131"] = peak_df["131_r"]/peak_df["131_r"].sum()

    peak_df_norm_ratio = peak_df.copy()

    f=[]
    reference_ind=0
    channels = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131"]
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
        peak_df_norm_ratio[channels[i]+"_1"] = peak_df_norm_ratio[channels[i]]/scaling_factors[i]
        peak_df_norm_ratio["log2_"+channels[i]+"_1"] = np.log2(peak_df_norm_ratio[channels[i]+"_1"])
    
    peak_df_norm_ratio["log2_127Nto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_127N_1"]
    peak_df_norm_ratio["log2_127Cto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_127C_1"]
    peak_df_norm_ratio["log2_128Nto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_128N_1"]
    peak_df_norm_ratio["log2_128Cto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_128C_1"]
    peak_df_norm_ratio["log2_129Nto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_129N_1"]
    peak_df_norm_ratio["log2_129Cto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_129C_1"]
    peak_df_norm_ratio["log2_130Nto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_130N_1"]
    peak_df_norm_ratio["log2_130Cto126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_130C_1"]
    peak_df_norm_ratio["log2_131to126"] = - peak_df_norm_ratio["log2_126_1"] + peak_df_norm_ratio["log2_131_1"]
    return peak_df_norm_ratio

def psm_summarization(sn_df, intensity_df, annotation_df, row_norm_user, num_channels):
    # Filtered PSMs using non-shared PSMs, 
    psm_df = compile_sn_raw_int(sn_df,intensity_df,num_channels)
    #print(list(psm_df))
    cols_abund = []
    if num_channels == 10:
        cols_abund = ['Abundance: 126', 'Abundance: 127N', 'Abundance: 127C','Abundance: 128N', 'Abundance: 128C', 'Abundance: 129N', 'Abundance: 129C', 'Abundance: 130N', 'Abundance: 130C', 'Abundance: 131']

    psm_df["#Measurements"] = psm_df[cols_abund].shape[1] - psm_df[cols_abund].isnull().sum(axis=1)
    psm_df["weights"] = np.log(psm_df["Intensity"])/np.log(1.2)

    #print(psm_df.iloc[0][cols_abund])
    # Check if user wants to do row-normalization
    if row_norm_user == True:
        psm_df["Norm-factor"] = psm_df[cols_abund].mean(axis=1)
        for c in cols_abund:
            psm_df[c] = psm_df[c]/psm_df["Norm-factor"]
    #print(psm_df.iloc[0][cols_abund])

    cols_reporter = []
    for c in cols_abund:
        psm_df["log2_"+c.split(": ")[1]+"_1"] = np.log2(psm_df[c])
        cols_reporter.append("log2_"+c.split(": ")[1]+"_1")
    #print(psm_df.iloc[0][cols_abund])

    # remove all rows or PSMs with too few measurements (<=2)
    #print(list(set(psm_df["#Measurements"].values.tolist())))
    psm_df = psm_df.loc[psm_df["#Measurements"]>2]

    # remove all PSMs with IS filter II<=50, SPS%>=65
    psm_df = psm_df.loc[(psm_df["Isolation Interference [%]"]<=30)]#&(psm_df["SPS Mass Matches [%]"]>=65)]

    #Only using master protein unique peptides for fair comparison
    psm_df = psm_df.loc[psm_df["# Protein Groups"]==1]
    psm_df["PSMID"] = psm_df["Annotated Sequence"] + "_" + psm_df["Charge"].astype(str) + "_"+ psm_df["Spectrum File"].str.split("_").str[-1]
    #print(psm_df.shape[0],len(list(set(psm_df["PSMID"].values.tolist()))))

    psm_df["PSMID2"] = psm_df["Annotated Sequence"] + "_" + psm_df["Charge"].astype(str) + "_"+ psm_df["Spectrum File"].str.split("_").str[-1] + "_" + psm_df["Master Protein Accessions"].astype(str)
    #print(psm_df.shape[0],len(list(set(psm_df["PSMID2"].values.tolist()))))


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

    #peak_df.to_excel("check_peaks.xlsx")

    for i in range(peak_df.shape[0]):
        peak = peak_df.iloc[i]["PeakID"]
        psms_assoc = psm_df.loc[psm_df["PSMID"]==peak]
        #if peak=="[K].gPsPPGAk.[R]_3_3.raw": #[K].lGLQcLPSDGVQNVNQ.[-]_3_2.raw, [K].sHHANsPTAGAAk.[S]_3_2.raw
            #print(psms_assoc.shape[0],psms_assoc[cols_reporter],psms_assoc["log2_126"].values.tolist()[0])
        if psms_assoc.shape[0]>1:
            #print(psms_assoc.shape[0])
            for c in cols_reporter:
                psms_c = psms_assoc[[c, "weights"]].dropna()
                weights = np.array(psms_c["weights"].values.tolist())
                values = np.array(psms_c[c].values.tolist())
                #print(weights, values)
                #print(np.nansum(values*weights)/np.nansum(weights))
                peak_df.loc[peak_df["PeakID"]==peak,c] = np.nansum(values*weights)/(np.nansum(weights))
        else:
            for c in cols_reporter:
                peak_df.loc[peak_df["PeakID"]==peak,c] = psms_assoc[c].values.tolist()[0]
    
    return peak_df


def compile_sn_raw_int(sn_df,raw_df,num_channels):
    raw_df["Features"] = raw_df["Annotated Sequence"]+"_"+raw_df["Charge"].astype(str)+"_"+raw_df["Spectrum File"] +"_" + raw_df["RT [min]"].astype(str)
    sn_df["Features"] = sn_df["Annotated Sequence"]+"_"+sn_df["Charge"].astype(str)+"_"+sn_df["Spectrum File"]+"_" + sn_df["RT [min]"].astype(str)
    
    tmp_df = raw_df
    tmp_df_filtered = tmp_df
    cols_noise = []
    if num_channels == 10:
        tmp_df["SN-126"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 126"].to_dict())
        tmp_df["SN-127N"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 127N"].to_dict())
        tmp_df["SN-127C"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 127C"].to_dict())
        tmp_df["SN-128N"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 128N"].to_dict())
        tmp_df["SN-128C"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 128C"].to_dict())
        tmp_df["SN-129N"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 129N"].to_dict())
        tmp_df["SN-129C"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 129C"].to_dict())
        tmp_df["SN-130N"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 130N"].to_dict())
        tmp_df["SN-130C"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 130C"].to_dict())
        tmp_df["SN-131"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 131"].to_dict())

        tmp_df["Noise-126"] = tmp_df["Abundance: 126"]/tmp_df["SN-126"]
        tmp_df["Noise-127N"] = tmp_df["Abundance: 127N"]/tmp_df["SN-127N"]
        tmp_df["Noise-127C"] = tmp_df["Abundance: 127C"]/tmp_df["SN-127C"]
        tmp_df["Noise-128N"] = tmp_df["Abundance: 128N"]/tmp_df["SN-128N"]
        tmp_df["Noise-128C"] = tmp_df["Abundance: 128C"]/tmp_df["SN-128C"]
        tmp_df["Noise-129N"] = tmp_df["Abundance: 129N"]/tmp_df["SN-129N"]
        tmp_df["Noise-129C"] = tmp_df["Abundance: 129C"]/tmp_df["SN-129C"]
        tmp_df["Noise-130N"] = tmp_df["Abundance: 130N"]/tmp_df["SN-130N"]
        tmp_df["Noise-130C"] = tmp_df["Abundance: 130C"]/tmp_df["SN-130C"]
        tmp_df["Noise-131"] = tmp_df["Abundance: 131"]/tmp_df["SN-131"]

        cols_noise = ["Noise-126","Noise-127N","Noise-127C",
                    "Noise-128N","Noise-128C",
                    "Noise-129N","Noise-129C",
                    "Noise-130N","Noise-130C","Noise-131"]
        tmp_df_filtered = tmp_df.loc[(tmp_df["SN-126"]>=1.5)|(tmp_df["SN-127N"]>=1.5)|(tmp_df["SN-127C"]>=1.5)|(tmp_df["SN-128N"]>=1.5)|(tmp_df["SN-128C"]>=1.5)|(tmp_df["SN-129N"]>=1.5)|(tmp_df["SN-129C"]>=1.5)|(tmp_df["SN-130N"]>=1.5)|(tmp_df["SN-130C"]>=1.5)|(tmp_df["SN-131"]>=1.5)]

    elif num_channels == 6:
        tmp_df["SN-126"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 126"].to_dict())
        tmp_df["SN-127"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 127"].to_dict())
        tmp_df["SN-128"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 128"].to_dict())
        tmp_df["SN-129"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 129"].to_dict())
        tmp_df["SN-130"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 130"].to_dict())
        tmp_df["SN-131"] = tmp_df.Features.map(sn_df.set_index("Features")["Abundance: 131"].to_dict())

        tmp_df["Noise-126"] = tmp_df["Abundance: 126"]/tmp_df["SN-126"]
        tmp_df["Noise-127"] = tmp_df["Abundance: 127"]/tmp_df["SN-127"]
        tmp_df["Noise-128"] = tmp_df["Abundance: 128"]/tmp_df["SN-128"]
        tmp_df["Noise-129"] = tmp_df["Abundance: 129"]/tmp_df["SN-129"]
        tmp_df["Noise-130"] = tmp_df["Abundance: 130"]/tmp_df["SN-130"]
        tmp_df["Noise-131"] = tmp_df["Abundance: 131"]/tmp_df["SN-131"]

        cols_noise = ["Noise-126","Noise-127","Noise-128","Noise-129","Noise-130","Noise-131"]
        tmp_df_filtered = tmp_df.loc[(tmp_df["SN-126"]>=1.5)|(tmp_df["SN-127"]>=1.5)|(tmp_df["SN-128"]>=1.5)|(tmp_df["SN-129"]>=1.5)|(tmp_df["SN-130"]>=1.5)|(tmp_df["SN-131"]>=1.5)]

    
    return tmp_df_filtered


	