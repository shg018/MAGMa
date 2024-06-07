from django.shortcuts import redirect, render
from .models import *
from .forms import *

from datetime import date
import pandas as pd
from itertools import combinations

from django.core.files.storage import FileSystemStorage
from django.shortcuts import render, HttpResponse

import os
import errno

# File for uniprot to gene symbol
#1. Local, second- cbsuhy01, third- cos75tst
uniprotGeneHuman = pd.read_csv('/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0/Example_input/Example_input/Input_files/uniprot_humancov192geneName.txt',sep='\t',header=None)
uniprotGeneYeast = pd.read_csv('/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0/Example_input/Example_input/Input_files/uniprot_yeast2geneName.txt',sep='\t',header=None)
uniprotGeneEcoli = pd.read_csv('/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0/Example_input/Example_input/Input_files/uniprot_ecoli2geneName.txt',sep='\t',header=None)
uniprotComplexes = pd.read_csv('/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0/Example_input/Example_input/Input_files/humanComplexes.txt',sep='\t')

#uniprotGeneHuman = pd.read_csv('/workdir/sg2369/LAVA_MS/Example_input/Example_input/Input_files/uniprot_humancov192geneName.txt',sep='\t',header=None)

# uniprotGeneHuman = pd.read_csv('/data/lava/Example_input/Example_input/Input_files/uniprot_humancov192geneName.txt',sep='\t',header=None)
# uniprotGeneYeast = pd.read_csv('/data/lava/Example_input/Example_input/Input_files/uniprot_yeast2geneName.txt',sep='\t',header=None)
# uniprotGeneEcoli = pd.read_csv('/data/lava/Example_input/Example_input/Input_files/uniprot_ecoli2geneName.txt',sep='\t',header=None)
# uniprotComplexes = pd.read_csv('/data/lava/Example_input/Example_input/Input_files/humanComplexes.txt',sep='\t')

tmp=uniprotGeneHuman[1].str.split(';',n=1,expand=True)[0].str.split(' {',n=1,expand=True)
uniprotGeneHuman["Gene Symbol"] = tmp[0]
uniprotGeneHuman["Accession"]=uniprotGeneHuman[0]
uniprotGeneHuman["Data_type"] = "Human"

tmp=uniprotGeneYeast[1].str.split(';',n=1,expand=True)[0].str.split(' {',n=1,expand=True)
uniprotGeneYeast["Gene Symbol"] = tmp[0]
uniprotGeneYeast["Accession"]=uniprotGeneYeast[0]
uniprotGeneYeast["Data_type"] = "Yeast"

tmp=uniprotGeneEcoli[1].str.split(';',n=1,expand=True)[0].str.split(' {',n=1,expand=True)
uniprotGeneEcoli["Gene Symbol"] = tmp[0]
uniprotGeneEcoli["Accession"]=uniprotGeneEcoli[0]
uniprotGeneEcoli["Data_type"] = "Ecoli"

human_uniprots = ','.join(map(str,list(set(uniprotGeneHuman.loc[~uniprotGeneHuman["Accession"].str.contains("ZZZZ")]["Accession"].values.tolist()))))
yeast_uniprots = ','.join(map(str,list(set(uniprotGeneYeast["Accession"].values.tolist()))))
ecoli_uniprots = ','.join(map(str,list(set(uniprotGeneEcoli["Accession"].values.tolist()))))

# All uniprots used for gene symbol assignment
uniprotGene = pd.concat([uniprotGeneHuman,uniprotGeneYeast,uniprotGeneEcoli], ignore_index=True)

required_file_columns=['Intensity', 'Annotated Sequence', 'Master Protein Accessions', 'Protein Accessions', 'Isolation Interference [%]', 'Spectrum File', 'Charge', 'RT [min]', '# Protein Groups']
required_file_comet_columns=["spectrum","precursor_intensity","peptide","assumed_charge","retention_time_sec"]


def bait_error(request):
    return render(request, 'bait_error.html')

def SN_error(request):
    return render(request, 'SN_error.html')

def intensity_error(request):
    return render(request, 'Intensity_error.html')

###############################################################################################
################# INITIAL UPLOAD FUNCTIONS FOR IP-DIA ############################################
###############################################################################################

def processuploadIPDIA(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,fusion_prot_seq,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        print(request.session['job_id'])
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print('The third file name is %s\n'%(file_name_3))
        print('The job id is %s\n'%(form_3.fields['job_id'].initial))

        column_list = ["Channel", "Label", "Control"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label)>0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels_cntrls=file_3_df.loc[file_3_df["Control"]==True].Label.unique()
        uniq_labels_samples=file_3_df.loc[file_3_df["Control"]==False].Label.unique()
        print ("Unique control labels: ",uniq_labels_cntrls, "Unique sample labels: ",uniq_labels_samples)
        annotation_file = file_name_3
        condition_comb_vs_cntrl=[]#list(combinations(uniq_labels_cntrl+uniq_labels_samples,2))
        for c1 in uniq_labels_samples:
            for c2 in uniq_labels_cntrls:
                condition_comb_vs_cntrl.append((c1,c2))
        condition_comb_no_cntrl=list(combinations(uniq_labels_samples,2))
        if len(condition_comb_vs_cntrl)==0:
            condition_comb_vs_cntrl.append(('',''))
        if len(condition_comb_no_cntrl)==0:
            condition_comb_no_cntrl.append(('',''))


    else:
        print('Reached second stage of input')
        if form_1.is_valid():
            print(request.session['job_id'])
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)
            file_name_1=request.FILES['docfile_1'].name
            
            # Not checking for missing columns
            print('file DIANN', file_name_1)

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)

        if form_2.is_valid():
            print(request.session['job_id'])
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)
            file_name_2=request.FILES['docfile_2'].name
            
            # Not checking for missing columns
            print ('file DIANN raw', file_name_2)

            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)


        if form_4.is_valid():
            dummy=0
            print (form_4.cleaned_data)
            bait_name=form_4.cleaned_data['bait_name']
            annotation_file=form_4.cleaned_data['annotation_file']
            print (bait_name, file_name_1)
            if bait_name == "NA":
                bait_validation = True
            else:
                bait_validation=validate_bait(bait_name,"diann", '%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_1))
            print ('bait validation status: %s'%(bait_validation) )
            
            bait_message='Your submitted bait name is: "%s" and its validation is: %s'%(bait_name, bait_validation)

        if (bait_validation)&(validation_2 == '')&(validation_1 == ''):
            submit_success=True

        print('Combinations choice:',request.POST['combinations'], )

        print('Direction choice:',request.POST['Analysis Direction'])


        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid(), form_4.is_valid())
    return bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name


def processDIAuploadWholeProteome(request,bait_validation,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        print(request.session['job_id'])
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["Channel", "Label"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels=file_3_df.Label.unique()
        print ("Unique labels: ",uniq_labels)
        annotation_file = file_name_3
        condition_comb=list(combinations(uniq_labels,2))
        # bait_validation=False
    else:
        if form_1.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)
            file_name_1=request.FILES['docfile_1'].name
            
            # missing_file_columns_1= [x for x in required_file_columns if x not in file_columns_1]
            print ('file 1', file_name_1)

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)
            bait_validation=True

        if form_2.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)
            file_name_2=request.FILES['docfile_2'].name
            
            # missing_file_columns_1= [x for x in required_file_columns if x not in file_columns_1]
            print ('file 2', file_name_2)

            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)
            bait_validation=True
            annotation_file=form_1.cleaned_data['annotation_file']


        print('Combinations choice:',request.POST['combinations'], )

        print('Direction choice:',request.POST['Analysis Direction'])

        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid())

    return bait_validation,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name

#################################################################################################

# Function determines which type of analysis to run for DIA datasets for Bruker search,
# Precursor or Protein based or IP-DIA-LFQ
# Magma type (PD, Linear model or Limma based) and other choices made by user

# OUTPUT: prot_df. Pandas dataframe with FC, FDR, protein abundances calculated

#################################################################################################

def choose_combination_process_DIAwholeprot(job_id,analysis_type,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,init_df,raw_input_df,annotation_df,uniprotGene,imputation_user):
    prot_df=[]
    model_type = "" 
    if (analysis_type == "LM"):
        model_type = "LinearModel"
    elif analysis_type=="Limma":
        model_type = "Limma"
    else:
        print("Not supported")
        return prot_df
    # Generate the protein summary file
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df = YuLab_fractionxchannel_DIA(job_id,init_df,raw_input_df,annotation_df,'',combination, r, True, True, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df = YuLab_fractionxchannel_DIA(job_id,init_df,raw_input_df,annotation_df,'',combination, r, True, False, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df = YuLab_fractionxchannel_DIA(job_id,init_df,raw_input_df,annotation_df,'',combination, r, True, True, False, num_channels, column_norm_proteins, model_type, uniprotGene,'',imputation_user)
    else:
        print("Need to select a column normalization option")
    return prot_df

def choose_combination_process_IPDIA(job_id,analysis_type,control_or_not,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,init_df,raw_input_df,annotation_df,uniprotGene,imputation_user):
    prot_df=[]
    bait_norm = False
    if control_or_not==")_NotC":
        if bait_name == "NA":
            bait_norm = False
        else:
            bait_norm = True
    model_type = "" 
    if (analysis_type == "LM"):
        model_type = "LinearModel"
    elif analysis_type=="Limma":
        model_type = "Limma"
    else:
        print("Not supported")
        return prot_df
    # Generate the protein summary file
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df = YuLab_fractionxchannel_IP_DIA(job_id,init_df,raw_input_df,annotation_df,bait_name,combination, r, False, True, bait_norm, num_channels, '', model_type, uniprotGene,imputation_user,'')
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df = YuLab_fractionxchannel_IP_DIA(job_id,init_df,raw_input_df,annotation_df,bait_name,combination, r, False, False, bait_norm, num_channels, '', model_type, uniprotGene, imputation_user,'')
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df = YuLab_fractionxchannel_IP_DIA(job_id,init_df,raw_input_df,annotation_df,bait_name,combination, r, False, True, bait_norm, num_channels, column_norm_proteins, model_type, uniprotGene, imputation_user,'')
    else:
        print("Need to select a column normalization option")
        
    return prot_df




def choose_combination_process_DIA(job_id,kind_of_dia,analysis_type,control_or_not,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,init_df,raw_input_df,annotation_df,uniprotGene,imputation_user):
    prot_df = []
    if (kind_of_dia == "protein")|(kind_of_dia == "precursor"):
        prot_df = choose_combination_process_DIAwholeprot(job_id,analysis_type,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,init_df,raw_input_df,annotation_df,uniprotGene,imputation_user)
    elif kind_of_dia == "IP":
        prot_df = choose_combination_process_IPDIA(job_id,analysis_type,control_or_not,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,init_df,raw_input_df,annotation_df,uniprotGene,imputation_user)
    else:
        print("No DIA type of experiment choice made")

    return prot_df

###############################################################################################
################# INITIAL UPLOAD FUNCTIONS FOR SILAC ##########################################
###############################################################################################

def processSILACIPuploadComet(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,form_fusion_seq,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    # print()
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["File","Channel", "Label", "Control"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels_cntrls=file_3_df.loc[file_3_df["Control"]==True].Label.unique()
        uniq_labels_samples=file_3_df.loc[file_3_df["Control"]==False].Label.unique()
        print ("Unique control labels: ",uniq_labels_cntrls, "Unique sample labels: ",uniq_labels_samples)
        annotation_file = file_name_3
        condition_comb_vs_cntrl=[]#list(combinations(uniq_labels_cntrl+uniq_labels_samples,2))
        for c1 in uniq_labels_samples:
            for c2 in uniq_labels_cntrls:
                condition_comb_vs_cntrl.append((c1,c2))
        condition_comb_no_cntrl=list(combinations(uniq_labels_samples,2))
        if len(condition_comb_vs_cntrl)==0:
            condition_comb_vs_cntrl.append(('',''))
        if len(condition_comb_no_cntrl)==0:
            condition_comb_no_cntrl.append(('',''))

    else:
        if form_1.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)

            file_name_1=request.FILES['docfile_1'].name
            file_1_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_1),sep="\t")
            missing_file_columns_1= [x for x in required_file_comet_columns if x not in list(file_1_df)]
            print ('file 1', file_name_1)

            if len(missing_file_columns_1)>0:
                validation_1='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_1))
                return render(request, 'SN_error.html')

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)

        if form_2.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)

            file_name_2=request.FILES['docfile_2'].name
            file_2_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_2),sep="\t")
            missing_file_columns_2= [x for x in required_file_comet_columns if x not in list(file_2_df)]
            print ('file 2', file_name_2)
            

            if len(missing_file_columns_2)>0:
                validation_2='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_2))
                return redirect(intensity_error)
            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)

        if form_4.is_valid():
            dummy=0
            print (form_4.cleaned_data)
            bait_name=form_4.cleaned_data['bait_name']
            annotation_file=form_4.cleaned_data['annotation_file']
            if bait_name == "NA":
                bait_validation = True
            else:
                print (bait_name, file_name_1)
                bait_validation=validate_bait(bait_name,"comet", '%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_1))            
            print ('bait validation status: %s'%(bait_validation) )
            
            bait_message='Your submitted bait name is: "%s" and its validation is: %s'%(bait_name, bait_validation)
        if bait_validation==False:
            submit_message ='ERROR! Your submitted bait Uniprot Id "%s" is missing from one or both the input files! Please check and resubmit'%(bait_name)
            return redirect(bait_error)
        if ((validation_1=='')&(validation_2=='')&(bait_validation==True)):
            submit_success = True

        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Files submitted correctly? :',submit_success)
        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid(), form_4.is_valid())
    
    return bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name

def processSILACwholeprotuploadComet(request,bait_validation,analysis_type,form_1,form_2,form_3,form_prot_norm,form_fusion_seq,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    # print()
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["File","Channel", "Label", "Control"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels_cntrls=file_3_df.loc[file_3_df["Control"]==True].Label.unique()
        uniq_labels_samples=file_3_df.loc[file_3_df["Control"]==False].Label.unique()
        print ("Unique control labels: ",uniq_labels_cntrls, "Unique sample labels: ",uniq_labels_samples)
        annotation_file = file_name_3
        condition_comb_vs_cntrl=[]#list(combinations(uniq_labels_cntrl+uniq_labels_samples,2))
        for c1 in uniq_labels_samples:
            for c2 in uniq_labels_cntrls:
                condition_comb_vs_cntrl.append((c1,c2))
        condition_comb_no_cntrl=list(combinations(uniq_labels_samples,2))
        if len(condition_comb_vs_cntrl)==0:
            condition_comb_vs_cntrl.append(('',''))
        if len(condition_comb_no_cntrl)==0:
            condition_comb_no_cntrl.append(('',''))

    else:
        if form_1.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)

            file_name_1=request.FILES['docfile_1'].name
            file_1_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_1),sep="\t")
            missing_file_columns_1= [x for x in required_file_comet_columns if x not in list(file_1_df)]
            print ('file 1', file_name_1)

            if len(missing_file_columns_1)>0:
                validation_1='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_1))
                return render(request, 'SN_error.html')

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)
            annotation_file=form_1.cleaned_data['annotation_file']

        if form_2.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)

            file_name_2=request.FILES['docfile_2'].name
            file_2_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_2),sep="\t")
            missing_file_columns_2= [x for x in required_file_comet_columns if x not in list(file_2_df)]
            print ('file 2', file_name_2)
            

            if len(missing_file_columns_2)>0:
                validation_2='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_2))
                return redirect(intensity_error)
            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)

        bait_validation=False
        bait_name="NA"
        if ((validation_1=='')&(validation_2=='')):
            submit_success = True

        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Files submitted correctly? :',submit_success)
        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid(), form_4.is_valid())
    
    return bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name

#################################################################################################

# Function determines which type of analysis to run for SILAC datasets for COMET search,
# Whole Proteome or Regular IP-SILAC 
# Magma type (Limma based)

# OUTPUT: prot_df. Pandas dataframe with FC, FDR, protein abundances calculated

#################################################################################################

def regular_choose_combination_process_SILAC(job_id,analysis_type,control_or_not,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,fwd_df,rev_df,annotation_df,uniprotGene):
    prot_df=[]
    pep_df=[]
    model_type=""
    # if (analysis_type == "CometLimma"):
    #     model_type="Limma"
    # else:
    #     print("Not supported")
    #     return (prot_df,pep_df)

    # bait_norm = False
    # if control_or_not==")_NotC":
    #     if bait_name == "NA":
    #         bait_norm = False
    #     else:
    #         bait_norm = True
    # Generate the protein and peptide summary file(s)
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df,pep_df = YuLab_SILAC(job_id, fwd_df, rev_df, annotation_df, bait_name, combination, r, row_norm_choice, True, '', uniprotGene)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df,pep_df = YuLab_SILAC(job_id, fwd_df, rev_df, annotation_df, bait_name, combination, r, row_norm_choice, False, '', uniprotGene)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df,pep_df = YuLab_SILAC(job_id, fwd_df, rev_df, annotation_df, bait_name, combination, r, row_norm_choice, True, column_norm_proteins, uniprotGene)
    else:
        print("Need to select an option")
    return (prot_df,pep_df)

def choose_combination_process_SILAC(job_id, kind_of_ip,analysis_type,control_or_not,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,fwd_df,rev_df,annotation_df,uniprotGene):
    prot_df = []
    pep_df = []
    if kind_of_ip == "regular":
        prot_df,pep_df = regular_choose_combination_process_SILAC(job_id,analysis_type,control_or_not,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,fwd_df,rev_df,annotation_df,uniprotGene)
    elif kind_of_ip == "whole_proteome":
        prot_df,pep_df = regular_choose_combination_process_SILAC(job_id,analysis_type,control_or_not,row_norm_choice,column_norm_choice,column_norm_proteins,"NA",combination,r,fwd_df,rev_df,annotation_df,uniprotGene)
    else:
        print("No SILAC type of experiment choice made")
    return (prot_df,pep_df)


###############################################################################################
################# INITIAL UPLOAD FUNCTIONS FOR TMT ############################################
###############################################################################################

def processTMTwholeprotuploadPD(request,submit_success,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["Channel", "Label"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels=file_3_df.Label.unique()
        print ("Unique labels: ",uniq_labels)
        annotation_file = file_name_3
        condition_comb=list(combinations(uniq_labels,2))

    else:
        print("Reading in file 1 and file 2")
        if form_1.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)

            file_name_1=request.FILES['docfile_1'].name
            file_1_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_1),sep="\t")
            missing_file_columns_1= [x for x in required_file_columns if x not in list(file_1_df)]
            print ('file 1', file_name_1)

            if len(missing_file_columns_1)>0:
                validation_1='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_1))
                return redirect(SN_error)

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)
            # bait_validation=True
            annotation_file=form_1.cleaned_data['annotation_file']

        
        if form_2.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)
            file_name_2=request.FILES['docfile_2'].name
            print ('The second file name is %s\n'%(file_name_2))

            file_2_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_2),sep="\t")
            missing_file_columns_2= [x for x in required_file_columns if x not in list(file_2_df)]

            if len(missing_file_columns_2)>0:
                validation_2='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_2))
                return redirect(intensity_error)
            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)
            
            # bait_validation=True
        
        if (validation_1 == '')&(validation_2 == ''):
            submit_success=True


        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Annotation file: ', annotation_file)
        print('Files submitted correctly? :',submit_success)
        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid())
    print('Files submitted correctly? :',submit_success)

    return submit_success,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name


def processTMTwholeprotuploadCOMET(request,submit_success,analysis_type,form_comet,form_3,form_prot_norm,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    file_name_comet=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["Channel", "Label"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels=file_3_df.Label.unique()
        print ("Unique labels: ",uniq_labels)
        annotation_file = file_name_3
        condition_comb=list(combinations(uniq_labels,2))

    else:
        if form_comet.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_comet.save(upload_path=upload_path)

            file_name_comet=request.FILES['docfile_comet'].name
            file_comet_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_comet),sep="\t")
            missing_file_columns_comet= [x for x in required_file_comet_columns if x not in list(file_comet_df)]
            print ('file comet', file_name_comet)

            if len(missing_file_columns_comet)>0:
                validation_comet='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_comet))
                return redirect(SN_error)
            if validation_comet == '':
                    file_annotation_comet='The uploaded file name is "%s" and the submission is successful!'%(file_name_comet)
                    submit_success = True
            else:
                file_annotation_comet='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_comet)
            annotation_file=form_comet.cleaned_data['annotation_file']

        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Annotation file: ', annotation_file)
        print('Files submitted correctly? :',submit_success)

        print (form_comet.is_valid(), form_3.is_valid())
    return submit_success,file_name_comet,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name



def processTMTuploadCOMET(request,bait_validation,analysis_type,form_comet,form_3,form_4,form_prot_norm,fusion_prot_seq,run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    file_name_comet=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["Channel", "Label", "Control"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels_cntrls=file_3_df.loc[file_3_df["Control"]==True].Label.unique()
        uniq_labels_samples=file_3_df.loc[file_3_df["Control"]==False].Label.unique()
        print ("Unique control labels: ",uniq_labels_cntrls, "Unique sample labels: ",uniq_labels_samples)
        annotation_file = file_name_3
        condition_comb_vs_cntrl=[]#list(combinations(uniq_labels_cntrl+uniq_labels_samples,2))
        for c1 in uniq_labels_samples:
            for c2 in uniq_labels_cntrls:
                condition_comb_vs_cntrl.append((c1,c2))
        condition_comb_no_cntrl=list(combinations(uniq_labels_samples,2))
        if len(condition_comb_vs_cntrl)==0:
            condition_comb_vs_cntrl.append(('',''))
        if len(condition_comb_no_cntrl)==0:
            condition_comb_no_cntrl.append(('',''))

    else:
        if form_comet.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_comet.save(upload_path=upload_path)

            file_name_comet=request.FILES['docfile_comet'].name
            file_comet_df = pd.read_csv('%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_comet),sep="\t")
            missing_file_columns_comet= [x for x in required_file_comet_columns if x not in list(file_comet_df)]
            print ('file comet', file_name_comet)
            if len(missing_file_columns_comet)>0:
                validation_comet='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_comet))
                return redirect(SN_error)
            if validation_comet == '':
                    file_annotation_comet='The uploaded file name is "%s" and the submission is successful!'%(file_name_comet)
            else:
                file_annotation_comet='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_comet)

        if form_4.is_valid():
            dummy=0
            print (form_4.cleaned_data)
            bait_name=form_4.cleaned_data['bait_name']
            annotation_file=form_4.cleaned_data['annotation_file']
            if bait_name == "NA":
                bait_validation = True
            elif ";" in bait_name:
                baits = bait_name.split(";")
                print(baits, file_name_comet)
                validations = []
                for bait in baits:
                    if validate_bait(bait,"comet", '%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_comet)):
                        validations.append(1)
                    else:
                        validations.append(0)
                if sum(validations)==len(baits):
                    bait_validation=True
            else:
                print(bait_name, file_name_comet)
                bait_validation=validate_bait(bait_name,"comet", '%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_comet))            
            print('bait validation status: %s'%(bait_validation) )
            
            bait_message='Your submitted bait name is: "%s" and its validation is: %s'%(bait_name, bait_validation)

        if bait_validation==False:
            submit_message ='ERROR! Your submitted bait Uniprot Id "%s" is missing from one or both the input files! Please check and resubmit'%(bait_name)
            return redirect(bait_error)

        if (bait_validation==True)&(validation_comet==''):
            submit_success = True

        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Files submitted correctly? :',submit_success)
        print("Bait validated? ", bait_validation)
        print (form_comet.is_valid(), form_3.is_valid(), form_4.is_valid())
    return bait_validation,file_name_comet,run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name


def processTMTuploadPD(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,form_fusion_seq,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name):
    # print()
    file_name_1=''
    file_name_2=''
    if form_3.is_valid():
        upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
        form_3.save(upload_path=upload_path)
        file_name_3=request.FILES['docfile_3'].name
        print ('The third file name is %s\n'%(file_name_3))

        column_list = ["Channel", "Label", "Control"]
        file_3_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_3), usecols=column_list)
        file_3_df = file_3_df.loc[file_3_df["Label"]!="none"]
        if (len(file_3_df.Channel)>0  and len(file_3_df.Label) >0):
            file_annotation_3='SUCCESS: The uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)
        else:
            file_validation_3='ERROR: Double-check to make sure that the uploaded annotation file "%s" consists of both columns "Channel" and "Label"'%(file_name_3)

        uniq_labels_cntrls=file_3_df.loc[file_3_df["Control"]==True].Label.unique()
        uniq_labels_samples=file_3_df.loc[file_3_df["Control"]==False].Label.unique()
        print ("Unique control labels: ",uniq_labels_cntrls, "Unique sample labels: ",uniq_labels_samples)
        annotation_file = file_name_3
        condition_comb_vs_cntrl=[]#list(combinations(uniq_labels_cntrl+uniq_labels_samples,2))
        for c1 in uniq_labels_samples:
            for c2 in uniq_labels_cntrls:
                condition_comb_vs_cntrl.append((c1,c2))
        condition_comb_no_cntrl=list(combinations(uniq_labels_samples,2))
        if len(condition_comb_vs_cntrl)==0:
            condition_comb_vs_cntrl.append(('',''))
        if len(condition_comb_no_cntrl)==0:
            condition_comb_no_cntrl.append(('',''))

    else:
        if form_1.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_1.save(upload_path=upload_path)

            file_name_1=request.FILES['docfile_1'].name
            file_1_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_1),sep="\t")
            missing_file_columns_1= [x for x in required_file_columns if x not in list(file_1_df)]
            print ('file 1', file_name_1)

            if len(missing_file_columns_1)>0:
                validation_1='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_1))
                return render(request, 'SN_error.html')

            if validation_1 == '':
                    file_annotation_1='The uploaded file name is "%s" and the submission is successful!'%(file_name_1)
            else:
                file_annotation_1='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_1)

        if form_2.is_valid():
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_2.save(upload_path=upload_path)

            file_name_2=request.FILES['docfile_2'].name
            file_2_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_2),sep="\t")
            missing_file_columns_2= [x for x in required_file_columns if x not in list(file_2_df)]
            print ('file 2', file_name_2)

            if len(missing_file_columns_2)>0:
                validation_2='ERROR!!! The following columns are missing from your input file, please check and resubmit:\n %s\n'%(', '.join(missing_file_columns_2))
                return redirect(intensity_error)
            if validation_2 == '':
                    file_annotation_2='The uploaded file name is "%s" and the submission is successful!'%(file_name_2)
            else:
                file_annotation_2='The uploaded file name is "%s" and the submission is NOT successful.'%(file_name_2)

        if form_4.is_valid():
            dummy=0
            print (form_4.cleaned_data)
            bait_name=form_4.cleaned_data['bait_name']
            annotation_file=form_4.cleaned_data['annotation_file']
            if bait_name == "NA":
                bait_validation = True
            elif ";" in bait_name:
                baits = bait_name.split(";")
                print (baits, file_name_1)
                validations = []
                # TODO: Validation should be in both S/N and Intensity file
                for bait in baits:
                    if validate_bait(bait,"sequest", '%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_1)):
                        validations.append(1)
                    else:
                        validations.append(0)
                if sum(validations)==len(baits):
                    bait_validation=True
            else:
                print (bait_name, file_name_1)
                bait_validation=validate_bait(bait_name,"sequest", '%s/%s'%(local_storage_dir+'/'+request.session['job_id'], file_name_1))            
            print ('bait validation status: %s'%(bait_validation) )
            
            bait_message='Your submitted bait name is: "%s" and its validation is: %s'%(bait_name, bait_validation)
        if bait_validation==False:
            submit_message ='ERROR! Your submitted bait Uniprot Id "%s" is missing from one or both the input files! Please check and resubmit'%(bait_name)
            return redirect(bait_error)
        if ((validation_1=='')&(validation_2=='')&(bait_validation==True)):
            submit_success = True

        print('Combinations choice:',request.POST['combinations'], )
        print('Direction choice:',request.POST['Analysis Direction'])
        print('Files submitted correctly? :',submit_success)
        print (form_1.is_valid(), form_2.is_valid(), form_3.is_valid(), form_4.is_valid())
    
    return bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name

#################################################################################################

# Function determines which type of analysis to run for TMT datasets for Proteome Discoverer search,
# Fusion or Regular IP-TMT 
# Magma type (PD, Linear model or Limma based) and other choices made by user

# OUTPUT: prot_df. Pandas dataframe with FC, FDR, protein abundances calculated

#################################################################################################

def phosphopep_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user):
    pep_df = []
    control_or_not=")_IsC"
    model_type=""
    if (analysis_type == "PDLM")|(analysis_type == "CometLM"):
        model_type="LinearModel"
    elif (analysis_type=="PDLimma")|(analysis_type == "CometLimma"):
        model_type="Limma"
    else:
        print("Not supported")
        return (pep_df)
    # Generate the protein and peptide summary file(s)
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        pep_df = YuLab_phosphopep_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, True, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        pep_df = YuLab_phosphopep_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, False, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        pep_df = YuLab_phosphopep_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, True, False, num_channels, column_norm_proteins, model_type, uniprotGene,'',imputation_user)
    else:
        print("Need to select a column normalization option")
    return (pep_df)

def wholeprot_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user):
    prot_df = []
    pep_df = []
    control_or_not=")_IsC"
    model_type=""
    if (analysis_type == "PDLM")|(analysis_type == "CometLM"):
        model_type="LinearModel"
    elif (analysis_type=="PDLimma")|(analysis_type == "CometLimma"):
        model_type="Limma"
    else:
        print("Not supported")
        return (prot_df,pep_df)
    # Generate the protein and peptide summary file(s)
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, True, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, False, False, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,"None",combination, r, fraction_choice, row_norm_choice, True, False, num_channels, column_norm_proteins, model_type, uniprotGene,'',imputation_user)
    else:
        print("Need to select a column normalization option")
    return (prot_df,pep_df)

def regular_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user):
    prot_df=[]
    pep_df=[]
    model_type=""
    if (analysis_type == "PDLM")|(analysis_type == "CometLM"):
        model_type="LinearModel"
    elif (analysis_type=="PDLimma")|(analysis_type == "CometLimma"):
        model_type="Limma"
    else:
        print("Not supported")
        return (prot_df,pep_df)
    bait_norm = False
    if control_or_not==")_NotC":
        if bait_name == "NA":
            bait_norm = False
        else:
            bait_norm = True
    # Generate the protein and peptide summary file(s)
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, True, bait_norm, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, False, bait_norm, num_channels, '', model_type, uniprotGene,'',imputation_user)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df,pep_df = YuLab_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, True, bait_norm, num_channels, column_norm_proteins, model_type, uniprotGene,'',imputation_user)
    else:
        print("Need to select a column normalization option")

    return (prot_df,pep_df)

def fusion_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,fusion_prot_seq,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user):
    prot_df=[]
    pep_df=[]
    model_type=""
    if (analysis_type == "PDLM")|(analysis_type == "CometLM"):
        model_type="LinearModel"
    elif (analysis_type=="PDLimma")|(analysis_type == "CometLimma"):
        model_type="Limma"
    else:
        print("Not supported")
        return (prot_df,pep_df)

    bait_norm = False
    if control_or_not==")_NotC":
        if bait_name == "NA":
            bait_norm = False
        else:
            bait_norm = True
    # Generate the protein and peptide summary file(s)
    if column_norm_choice=="globalnorm":
        print("Selected Global normalization")
        prot_df,pep_df = YuLab_fusion_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, True, bait_norm, num_channels, '', model_type, uniprotGene, fusion_prot_seq,imputation_user)
    elif column_norm_choice=="nocolumnnorm":
        print("Selected to not do any column normalization")
        prot_df,pep_df = YuLab_fusion_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, False, bait_norm, num_channels, '', model_type, uniprotGene, fusion_prot_seq,imputation_user)
    elif column_norm_choice=="protnorm":
        print("Selected to do protein uniprot list based normalization")
        print("The list of proteins are: ", column_norm_proteins)
        if column_norm_proteins=="human":
            column_norm_proteins=human_uniprots
        if column_norm_proteins=="yeast":
            column_norm_proteins=yeast_uniprots
        if column_norm_proteins=="ecoli":
            column_norm_proteins=ecoli_uniprots
        prot_df,pep_df = YuLab_fusion_fractionxchannel(job_id,sn_df,intensity_df,annotation_df,bait_name,combination, r, fraction_choice, row_norm_choice, True, bait_norm, num_channels, column_norm_proteins, model_type, uniprotGene, fusion_prot_seq,imputation_user)
    else:
        print("Need to select an option")
    return (prot_df,pep_df)

def choose_combination_process_TMT(job_id, kind_of_ip,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,fusion_prot_seq,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user):
    prot_df = []
    pep_df = []
    if kind_of_ip == "fusion":
        prot_df,pep_df = fusion_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,fusion_prot_seq,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user)
    elif kind_of_ip == "regular":
        prot_df,pep_df = regular_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,bait_name,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user)
    elif kind_of_ip == "whole_proteome":
        prot_df,pep_df = wholeprot_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user)
    elif kind_of_ip == "phosphopep_proteome":
        pep_df = phosphopep_choose_combination_process_TMT(job_id,analysis_type,control_or_not,fraction_choice,row_norm_choice,column_norm_choice,column_norm_proteins,combination,r,num_channels,sn_df,intensity_df,annotation_df,uniprotGene,imputation_user)
    else:
        print("No TMT type of experiment choice made")
    return (prot_df,pep_df)


#################################################################################################

# Accessory functions for saving output files to the right directories
# INPUT: 
#   prot_df, pep_df: The protein and/ peptide summary files for saving, 
#   bait_name: Bait uniprot (in case of IP), 
#   BASE_DIR: Directory to get relative root filepath to save,
#   combination: Comparison that was run, 
#   job_id: Unique identifier associated with job in the backend, 
#   fdr_thresh: The FDR threshold to use as enrichment cutoff (default set to 5%)
# OUTPUT: The csv (protein  and/ peptide summary files), pdf (volcano plot)

#################################################################################################

def make_dir(path):
    try:
        os.makedirs(path, exist_ok=True)  # Python>3.2
    except TypeError:
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise

def load_and_save(prot_df, bait_name, BASE_DIR, combination, job_id, fdr_thresh=0.05):
    make_dir(BASE_DIR+"/Results/"+job_id+"/")
    prot_df.to_csv(BASE_DIR+"/Results/"+job_id+"/"+combination+"_FC_and_pval_YuLab.csv",index=False)
    
    # fdr_thresh = 0.05
    volcano_plt = generate_baseline_volcano_plot(prot_df, fdr_thresh, bait_name, combination, "one")
    volcano_plt.savefig(BASE_DIR+"/Results/"+job_id+"/"+combination+"_FC_and_pval_full.pdf",bbox_inches = "tight")
    volcano_plt.close()

    return

def load_and_save_SILAC(prot_df, pep_df, bait_name, BASE_DIR, combination, job_id, fdr_thresh=0.05):
    make_dir(BASE_DIR+"/Results/"+job_id+"/")
    pep_df.to_csv(BASE_DIR+"/Results/"+job_id+"/"+combination+"_peptide_level_SILAC_FC_and_pval_YuLab.csv",index=False)
    prot_df.to_csv(BASE_DIR+"/Results/"+job_id+"/"+combination+"_SILAC_FC_and_pval_YuLab.csv", index=False)
    
    # fdr_thresh = 0.05
    volcano_plt = generate_baseline_volcano_plot(prot_df, fdr_thresh, bait_name, combination, "one")
    volcano_plt.savefig(BASE_DIR+"/Results/"+job_id+"/"+combination+"_SILAC_FC_and_pval_full.pdf",bbox_inches = "tight")
    volcano_plt.close()

    return

def load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, combination, job_id, fdr_thresh=0.05):
    make_dir(BASE_DIR+"/Results/"+job_id+"/")
    pep_df.to_csv(BASE_DIR+"/Results/"+job_id+"/"+combination+"_peptide_level_FC_and_pval_YuLab.csv",index=False)
    prot_df.to_csv(BASE_DIR+"/Results/"+job_id+"/"+combination+"_FC_and_pval_YuLab.csv", index=False)
    
    # fdr_thresh = 0.05
    volcano_plt = generate_baseline_volcano_plot(prot_df, fdr_thresh, bait_name, combination, "one")
    volcano_plt.savefig(BASE_DIR+"/Results/"+job_id+"/"+combination+"_FC_and_pval_full.pdf",bbox_inches = "tight")
    volcano_plt.close()

    return