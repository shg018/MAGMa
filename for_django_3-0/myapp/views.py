from django.shortcuts import redirect, render
from .models import *
from .forms import *
from .views_accessory import *
from django.http import JsonResponse
# from .tasks import go_to_sleep


from datetime import date
#from datetime import datetime
import pandas as pd
from itertools import combinations
import statistics
from scipy import stats
from scipy.stats import pearsonr
import statsmodels
import statsmodels.api as sa
import statsmodels.formula.api as sfa

from django.core.files.storage import FileSystemStorage
from django.shortcuts import render, HttpResponse

import uuid

def generate_job_id():
    return str(uuid.uuid4())

# Create your views here.
def index(request):
    return render(request,'landingpage.html')
    
def about(request):
    return render(request,'about.html')

def downloads(request):
    return render(request,'downloads.html')

def usage(request):
    return render(request,'usage.html')

def bait_error(request):
    return render(request, 'bait_error.html')

def SN_error(request):
    return render(request, 'SN_error.html')

def intensity_error(request):
    return render(request, 'Intensity_error.html')


def FilterForIP(request):
    message_filterfileforip = ''
    file_annotation_filterfileforip = ''
    validation_filterfileforip = ''

    message_uniprotfile = ''
    file_annotation_uniprotfile = ''
    validation_uniprotfile = ''

    message_uniprotfile2 = ''
    file_annotation_uniprotfile2 = ''
    validation_uniprotfile2 = ''

    documents=''
    run_success=False
    # Handle file upload
    result_file_name=''
    vsC1file=''
    vsC2file=''
    volcano_baseline_name=''
    job_id=None

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'#BASE_DIR+"/Results/"
    if request.method == 'POST':
        form_filterfileforip = DocumentForm_FilterFileForIP(job_id,request.POST, request.FILES)
        form_uniprotfile = DocumentForm_UniprotFile(job_id,request.POST, request.FILES)
        form_uniprotfile2 = DocumentForm_UniprotFile2(job_id,request.POST, request.FILES)

        if (form_filterfileforip.is_valid())&(form_uniprotfile.is_valid())&(form_uniprotfile2.is_valid()):
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_filterfileforip.save(upload_path=upload_path)
            file_name_filterfileforip=request.FILES['docfile_filterfileforip'].name

            form_uniprotfile.save(upload_path=upload_path)
            file_name_uniprotfile=request.FILES['docfile_uniprotfile'].name

            form_uniprotfile2.save(upload_path=upload_path)
            file_name_uniprotfile2=request.FILES['docfile_uniprotfile2'].name

            print ('The filter file name is %s\n'%(file_name_filterfileforip))
            print ('The first file name with UniprotIDs is %s\n'%(file_name_uniprotfile))
            print ('The second file name with UniprotIDs is %s\n'%(file_name_uniprotfile2))
            

            column_list = ["Protein"]
            file_filterfileforip_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_filterfileforip))
            file_uniprotfile_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_uniprotfile))
            file_uniprotfile2_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_uniprotfile2))
            
            if ((column_list[0] in list(file_filterfileforip_df))&(column_list[0] in list(file_uniprotfile_df))&(column_list[0] in list(file_uniprotfile2_df))):
                file_annotation_filterfileforip='SUCCESS: The uploaded file(s) contain the column "Protein"'
            else:
                file_annotation_filterfileforip='ERROR: Double-check to make sure that the uploaded file(s) have the column "Protein"'

            
            # Filter file using uniprot IDs
            uniprot_ids = list(set(file_uniprotfile_df.loc[(file_uniprotfile_df["log2FC"]>=1)&(file_uniprotfile_df["adjpval"]<=0.05)&(file_uniprotfile_df["# PSMs"]>=5)]["Protein"].values.tolist()+file_uniprotfile2_df.loc[(file_uniprotfile2_df["log2FC"]>=1)&(file_uniprotfile2_df["adjpval"]<=0.05)&(file_uniprotfile2_df["# PSMs"]>=5)]["Protein"].values.tolist()))
            print(file_filterfileforip_df.shape[0],len(uniprot_ids))
            file_filterfileforip_df = file_filterfileforip_df.loc[file_filterfileforip_df["Protein"].isin(uniprot_ids)]
            print(file_filterfileforip_df.shape[0],len(uniprot_ids))

            make_dir(BASE_DIR+"/Results/"+request.session['job_id']+"/")
            result_file_name=request.session['job_id']+"/"+file_name_filterfileforip.strip(".csv")+"_filteredFile.csv"
            vsC1file=request.session['job_id']+"/"+file_name_uniprotfile
            vsC2file=request.session['job_id']+"/"+file_name_uniprotfile2

            file_filterfileforip_df.to_csv(BASE_DIR+"/Results/"+result_file_name,sep=",",index=False)
            file_uniprotfile_df.to_csv(BASE_DIR+"/Results/"+vsC1file,sep=",",index=False)
            file_uniprotfile2_df.to_csv(BASE_DIR+"/Results/"+vsC2file,sep=",",index=False)
            
            volcano_baseline_name = request.session['job_id']+"/"+file_name_filterfileforip.strip(".csv")+"_filteredFile.pdf"
            volcano_plt = generate_baseline_volcano_plot(file_filterfileforip_df, 0.05, "NA", "Filtered File (non Control comparison)", "both")
            volcano_plt.savefig(BASE_DIR+"/Results/"+volcano_baseline_name,bbox_inches = "tight")
            run_success=True

    else:
        form_filterfileforip = DocumentForm_FilterFileForIP(job_id)  # An empty, unbound form
        form_uniprotfile = DocumentForm_UniprotFile(job_id)  # An empty, unbound form
        form_uniprotfile2 = DocumentForm_UniprotFile2(job_id)  # An empty, unbound form
    
    job_id = generate_job_id()
    request.session['job_id'] = job_id
    context = {'job_id':request.session['job_id'],'result_file_name':result_file_name, 'vsC1file':vsC1file, 'vsC2file':vsC2file,'volcano_baseline_name':volcano_baseline_name,'run_success':run_success, 'form_filterfileforip': form_filterfileforip, 'message_filterfileforip': message_filterfileforip, 'file_annotation_filterfileforip': file_annotation_filterfileforip, 'validation_filterfileforip': validation_filterfileforip, 'form_uniprotfile': form_uniprotfile, 'message_uniprotfile': message_uniprotfile, 'file_annotation_uniprotfile': file_annotation_uniprotfile, 'validation_uniprotfile': validation_uniprotfile, 'form_uniprotfile2': form_uniprotfile2, 'message_uniprotfile2': message_uniprotfile2, 'file_annotation_uniprotfile2': file_annotation_uniprotfile2, 'validation_uniprotfile2': validation_uniprotfile2}

    return render(request,'FilterForIP.html', context)

def AdjustedPvalue(request):
    message_adjustedpvalue = ''
    file_annotation_adjustedpvalue=''
    validation_adjustedpvalue=''
    documents=''
    run_success=False
    job_id=None
    # Handle file upload
    result_file_name=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'#BASE_DIR+"/Results/"
    if request.method == 'POST':
        form_adjustedpvalue = DocumentForm_5(job_id,request.POST, request.FILES)#DocumentForm_AdjustedPvalue(job_id, request.POST, request.FILES)

        if form_adjustedpvalue.is_valid():
            # print(form_adjustedpvalue)
            # newdoc_adjustedpvalue = Document(docfile=request.FILES['docfile_adjustedpvalue'])   
            # newdoc_adjustedpvalue.save()
            # file_name_adjustedpvalue=request.FILES['docfile_adjustedpvalue'].name
            # print ('The third file name is %s\n'%(file_name_adjustedpvalue))
            # print ('The third file name is %s\n'%(request.FILES['docfile_adjustedpvalue']))

            #result_df = pd.read_csv("%s/%s"%(local_storage_dir,file_name_1),sep=",")
            #result_df.to_csv(BASE_DIR+"/Results/"+file_name_1,sep=",")
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_adjustedpvalue.save(upload_path=upload_path)
            file_name_adjustedpvalue=request.FILES['docfile_volcano'].name#docfile_adjustedpvalue
            print ('The uploaded file name is %s\n'%(file_name_adjustedpvalue))

            column_list = ["Protein", "pval","# PSMs"]
            file_adjustedpvalue_df = pd.read_csv('%s/%s'%(local_storage_dir+"/"+request.session['job_id'], file_name_adjustedpvalue))
            
            if ((column_list[0] in list(file_adjustedpvalue_df))&(column_list[1] in list(file_adjustedpvalue_df))&(column_list[2] in list(file_adjustedpvalue_df))):
                file_annotation_adjustedpvalue='SUCCESS: The uploaded file "%s" consists of both columns "Protein","pval","# PSMs" in'%(file_name_adjustedpvalue)
            else:
                file_annotation_adjustedpvalue='ERROR: Double-check to make sure that the uploaded file "%s" consists of columns "Protein","pval" in'%(file_name_adjustedpvalue)

            file_adjustedpvalue_df = file_adjustedpvalue_df.loc[file_adjustedpvalue_df["# PSMs"]>=5]
            
            # Make Adjusted Pvalue column 
            # TODO: Different ways of adjusted p-value calculation
            # Check that p-value column does not have empty values
            if (file_adjustedpvalue_df[["pval"]].shape[0]!=file_adjustedpvalue_df[["pval"]].dropna().shape[0]):
                print("Warning: P-value column has NAN's and these rows will be removed")
                file_adjustedpvalue_df.dropna(subset=["pval"], inplace=True)

            # Default for now - Benjamini Hochberg
            file_adjustedpvalue_df["adjpval"] = statsmodels.stats.multitest.multipletests(file_adjustedpvalue_df["pval"].values.tolist(),method="fdr_bh")[1]

            make_dir(BASE_DIR+"/Results/"+request.session['job_id']+"/")
            result_file_name=request.session['job_id']+"/"+file_name_adjustedpvalue.strip(".csv")+"_adjustedPval.csv"
            file_adjustedpvalue_df.to_csv(BASE_DIR+"/Results/"+result_file_name,sep=",",index=False)
    
            # file_adjustedpvalue_df.to_csv('%s/%s'%(local_storage_dir, result_file_name),sep=",")
            
            run_success=True
            # file_network_df["groupA"] = file_network_df.Master_protein_noiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())


    else:
        form_adjustedpvalue = DocumentForm_5(job_id)#DocumentForm_AdjustedPvalue(job_id)  # An empty, unbound form
    
    job_id = generate_job_id()
    request.session['job_id'] = job_id
    context = {'job_id':request.session['job_id'],'form_adjustedpvalue': form_adjustedpvalue, 'message_network': message_adjustedpvalue, 'file_annotation_adjustedpvalue': file_annotation_adjustedpvalue, 'validation_adjustedpvalue': validation_adjustedpvalue, 'result_file_name':result_file_name, 'run_success':run_success}

    return render(request,'AdjustedPvalue.html', context)


def volcanoplot(request):
    message_5 = ''
    #message='%s'%(int(date.today().year))
    file_annotation_5=''
    validation_5=''
    documents=''
    # Handle file upload
    result_file_name=''
    job_id=None
    job_id = generate_job_id()

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # print(BASE_DIR," : is Basedir")
    # BASE_DIR = "https://magma.yulab.org"
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'#BASE_DIR+"/Results/"
    if request.method == 'POST':
        form_5 = DocumentForm_5(job_id, request.POST, request.FILES)

        if form_5.is_valid():
            # newdoc_5 = Document(docfile=request.FILES['docfile_volcano'])
            upload_path = f'documents/tmtipms_files/'+request.session['job_id']+'/' 
            form_5.save(upload_path=upload_path)

            file_name_5=request.FILES['docfile_volcano'].name
            print ('The third file name is %s\n'%(file_name_5))
            print ('The third file name is %s\n'%(request.FILES['docfile_volcano']))

            #result_df = pd.read_csv("%s/%s"%(local_storage_dir,file_name_1),sep=",")
            #result_df.to_csv(BASE_DIR+"/Results/"+file_name_1,sep=",")
            result_file_name=request.session["job_id"]+"/"+file_name_5

    else:
        form_5 = DocumentForm_5(job_id)  # An empty, unbound form
    
    job_id = generate_job_id()
    request.session["job_id"] = job_id
    context = {'job_id':request.session["job_id"],'form_5': form_5, 'message_5': message_5, 'file_annotation_5': file_annotation_5, 'validation_5': validation_5, 'result_file_name':result_file_name}

    return render(request,'volcanoplot.html', context)

################################################################################
################# This section is for TMT processing ###########################
################################################################################

def TMTalltypes(request):
    return render(request,'TMTalltypes.html')

def TMTphosprot(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
   
    message_1 = ''
    message_comet = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_comet=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_comet=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id=None

    if (request.method == 'POST')&((analysis_type=="PDLM")|(analysis_type=="PDLimma")):
        form_1 = DocumentForm_1(job_id, request.POST, request.FILES) # S/N file PD
        form_2 = DocumentForm_2(job_id, request.POST, request.FILES) # INTENSITY file PD
        form_3 = DocumentForm_3(job_id, request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        submit_success,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTwholeprotuploadPD(request,submit_success,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        # print("Getting to first pass checking bait validation: ", bait_validation)
        
        if (submit_success==False)&(form_3.is_valid()):
            print("Getting to second pass chacking submit status: ", submit_success)
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            print(form_1.is_valid(), form_2.is_valid())
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTphosprot.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
            print(submit_message)
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_1),sep="\t")
            intensity_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name="None"
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"phosphopep_proteome",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, intensity_df, annotation_df, uniprotGene, imputation_user)
            print(pep_df.head(2))
            
            if not pep_df.empty:
                run_success=True

                pep_df = pep_df.loc[~((pep_df["log2FC"]==0.00)&(pep_df["pval"]==0.00))]
                make_dir(BASE_DIR+"/Results/"+request.session["job_id"]+"/")
                pep_df.to_csv(BASE_DIR+"/Results/"+request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv",index=False)
    
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTphosprot.html', context)

    elif (request.method == 'POST')&((analysis_type=="CometLM")|(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_comet = DocumentForm_Comet(job_id, request.POST, request.FILES) #SN
        form_3 = DocumentForm_3(job_id, request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        
        # submit_success=False
        submit_success,file_name_comet,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTwholeprotuploadCOMET(request,submit_success,analysis_type,form_comet,form_3,form_prot_norm,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (submit_success==False)&(form_3.is_valid()):
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_comet': form_comet, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTphosprotComet.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_comet),sep="\t")
            
            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # For COMET change the input file to have specific columns
            sn_df = COMETconversion(sn_df,annotation_df,num_channels)
            bait_name="None"
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"phosphopep_proteome",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, sn_df, annotation_df, uniprotGene, imputation_user)
            print(pep_df.head(2))

            if not pep_df.empty:
                run_success=True
                make_dir(BASE_DIR+"/Results/"+request.session["job_id"]+"/")
                pep_df.to_csv(BASE_DIR+"/Results/"+request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv",index=False)
    
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_comet': form_comet, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTphosprotComet.html', context)

    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()


    job_id = generate_job_id()
    request.session["job_id"] = job_id

    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session["job_id"],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'TMTphosprot.html', context)

def TMTwholeprot(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
   
    message_1 = ''
    message_comet = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_comet=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_comet=''
    validation_2=''
    validation_3=''

    run_success=False
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id=None

    if (request.method == 'POST')&((analysis_type=="PDLM")|(analysis_type=="PDLimma")):
        form_1 = DocumentForm_1(job_id, request.POST, request.FILES) # S/N file PD
        form_2 = DocumentForm_2(job_id, request.POST, request.FILES) # INTENSITY file PD
        form_3 = DocumentForm_3(job_id, request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        submit_success=False
        submit_success,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTwholeprotuploadPD(request,submit_success,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        if (submit_success==False)&(form_3.is_valid()):
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTwholeprot.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_1),sep="\t")
            intensity_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name="None"
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"whole_proteome",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, intensity_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))

            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session["job_id"], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTwholeprot.html', context)

    elif (request.method == 'POST')&((analysis_type=="CometLM")|(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_comet = DocumentForm_Comet(job_id, request.POST, request.FILES) #SN
        form_3 = DocumentForm_3(job_id, request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        
        submit_success=False
        submit_success,file_name_comet,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTwholeprotuploadCOMET(request,submit_success,analysis_type,form_comet,form_3,form_prot_norm,run_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (submit_success==False)&(form_3.is_valid()):
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_comet': form_comet, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTwholeprotComet.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_comet),sep="\t")
            
            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # For COMET change the input file to have specific columns
            sn_df = COMETconversion(sn_df,annotation_df,num_channels)
            bait_name="None"
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"whole_proteome",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, sn_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))

            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session["job_id"], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session["job_id"], 'form_comet': form_comet, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTwholeprotComet.html', context)

    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()


    job_id = generate_job_id()
    request.session["job_id"] = job_id

    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session["job_id"],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'TMTwholeprot.html', context)


def TMT(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
    message_1 = ''
    message_comet = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_comet=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_comet=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb_vs_cntrl=[]
    condition_comb_no_cntrl=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id=None
    
    if (request.method == 'POST')&((analysis_type=="PDLM")|(analysis_type=="PDLimma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # S/N file PD
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # INTENSITY file PD
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)

        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTuploadPD(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,[],run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):
            # print("Job submitted successfully? IN tmt ARM for some reason ", submit_success)
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMT.html', context)
        if submit_success == True:
            
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
            intensity_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"regular",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, intensity_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))

            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session['job_id']+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session['job_id']+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMT.html', context)

        
    elif (request.method == 'POST')&((analysis_type=="CometLM")|(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id)
        form_2 = DocumentForm_2(job_id)
        form_comet = DocumentForm_Comet(job_id,request.POST, request.FILES) #SN
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)

        bait_validation=False
        bait_validation,file_name_comet,run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTuploadCOMET(request,bait_validation,analysis_type,form_comet,form_3,form_4,form_prot_norm,[],run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        # print("PASSED TO BACKEND:", bait_validation, form_3.is_valid(), submit_success)
        if (bait_validation==False)&(form_3.is_valid()):
            # print("Job submitted successfully? ", submit_success)
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_comet': form_comet, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTComet.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
            print(submit_message)
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_comet),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # For COMET change the input file to have specific columns
            sn_df = COMETconversion(sn_df,annotation_df,num_channels)
            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"regular",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, num_channels, '', sn_df, sn_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))

            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session['job_id']+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session['job_id']+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_comet': form_comet, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTComet.html', context)

    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_4 = DocumentForm_4(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()

    job_id = generate_job_id()
    request.session["job_id"] = job_id
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'TMT.html', context)

def TMTfusion(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
    message_1 = ''
    message_comet = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_comet=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_comet=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    analysis_type = request.GET['analysisType']
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''


    documents=''
    
    # Handle file upload
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb_vs_cntrl=[]
    condition_comb_no_cntrl=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    fusion_prot_seq=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'

    job_id=None
    
    if (request.method == 'POST')&((analysis_type=="PDLM")|(analysis_type=="PDLimma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES)
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES)
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        form_fusion_seq = Document_FusionSeq(request.POST)
        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTuploadPD(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,form_fusion_seq,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, "form_fusion_seq":form_fusion_seq, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}
            return render(request, 'TMTfusion.html', context)
        if submit_success == True:
                
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True
            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_1),sep="\t")
            intensity_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 
            # TODO Has to get sequence of fusion protein so as to remove peptides not in this region of head/tail
            if form_fusion_seq.is_valid():
                fusion_prot_seq = form_fusion_seq.cleaned_data["fusion_seq"]

            prot_df = []
            pep_df = []

            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"fusion",analysis_type,combination_tmp[4], fraction_choice, row_norm_choice, column_norm_choice, column_norm_proteins, bait_name, combination, r, num_channels, fusion_prot_seq, sn_df, intensity_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))
            # load_and_save(prot_df, bait_name, BASE_DIR, combination, direction_choice, 0.05)
            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session["job_id"], 0.05)

            #prot_df.to_excel(combination+"_FC_and_pval_YuLab_fractionxchannel.xlsx",index=False)
            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, "form_fusion_seq":form_fusion_seq, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'TMTfusion.html', context)

    elif (request.method == 'POST')&((analysis_type=="CometLM")|(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_comet = DocumentForm_Comet(job_id,request.POST, request.FILES)
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)
        form_fusion_seq = Document_FusionSeq(request.POST)
        bait_validation=False
        bait_validation,file_name_comet,run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processTMTuploadCOMET(request,bait_validation,analysis_type,form_comet,form_3,form_4,form_prot_norm,form_fusion_seq,run_success,submit_success,message_comet,message_3,file_annotation_comet,file_annotation_3,validation_comet,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_comet': form_comet, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, "form_fusion_seq":form_fusion_seq, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}
            return render(request, 'TMTfusionComet.html', context)
        if submit_success == True:
                
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            fraction_choice=request.POST['fractionchoice']
            if fraction_choice == "yesfracsep":
                fraction_choice = True
            else:
                fraction_choice = False
            row_norm_choice=request.POST['rownorm']
            if row_norm_choice == "yesrownorm":
                row_norm_choice = True
            else:
                row_norm_choice = False
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True
            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            sn_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],file_name_comet),sep="\t")
            # intensity_df = pd.read_csv("%s/%s"%(local_storage_dir,file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session["job_id"],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # COMET specifically change the column format 
            sn_df = COMETconversion(sn_df, annotation_df, num_channels)
            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 
            # TODO Has to get sequence of fusion protein so as to remove peptides not in this region of head/tail
            if form_fusion_seq.is_valid():
                fusion_prot_seq = form_fusion_seq.cleaned_data["fusion_seq"]

            prot_df = []
            pep_df = []

            prot_df,pep_df = choose_combination_process_TMT(request.session["job_id"],"fusion",analysis_type,combination_tmp[4], fraction_choice,row_norm_choice, column_norm_choice, column_norm_proteins, bait_name, combination, r, num_channels, fusion_prot_seq, sn_df, sn_df, annotation_df, uniprotGene, imputation_user)
            print(prot_df.head(2))
            # load_and_save(prot_df, bait_name, BASE_DIR, combination, direction_choice, 0.05)
            load_and_save_TMT(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session["job_id"], 0.05)

            #prot_df.to_excel(combination+"_FC_and_pval_YuLab_fractionxchannel.xlsx",index=False)
            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session["job_id"]+"/"+cname+"_peptide_level_FC_and_pval_YuLab.csv"
                result_file_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session["job_id"]+"/"+cname+"_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_comet': form_comet, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, "form_fusion_seq":form_fusion_seq, 'message_comet': message_comet, 'message_3': message_3, 'file_annotation_comet': file_annotation_comet, 'file_annotation_3': file_annotation_3, 'validation_comet': validation_comet, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}
            
            return render(request, 'TMTfusionComet.html', context)
    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_4 = DocumentForm_4(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()
        form_fusion_seq = Document_FusionSeq()

    job_id = generate_job_id()
    request.session["job_id"] = job_id
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, "form_fusion_seq":form_fusion_seq, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'TMTfusion.html', context)

################################################################################
################# This section is for SILAC processing #########################
################################################################################
def SILACalltypes(request):
    return render(request,'SILACalltypes.html')

def SILACIP(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
    message_1 = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb_vs_cntrl=[]
    condition_comb_no_cntrl=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id=None
    
    if ((request.method == 'POST')&(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # FWD rep
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # REV rep
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)

        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processSILACIPuploadComet(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,[],run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)

        # print("PASSED TO BACKEND:", bait_validation, form_3.is_valid(), submit_success)
        if (bait_validation==False)&(form_3.is_valid()):
            # print("Job submitted successfully? ", submit_success)
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'SILACIP.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
            print(submit_message)
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            
            row_norm_choice=False#request.POST['rownorm']
            # if row_norm_choice == "yesrownorm":
            #     row_norm_choice = True
            # else:
            #     row_norm_choice = False
            

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            fwd_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
            rev_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")


            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # For COMET change the input file to have specific columns
            fwd_df = SILAC_COMETconversion(fwd_df,annotation_df,"FWD")
            rev_df = SILAC_COMETconversion(rev_df,annotation_df,"REV")

            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_SILAC(request.session["job_id"],"regular",analysis_type,combination_tmp[4], row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, fwd_df, rev_df, annotation_df, uniprotGene)
            print(prot_df.head(2))

            load_and_save_SILAC(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session['job_id']+"/"+cname+"_peptide_level_SILAC_FC_and_pval_YuLab.csv"
                result_file_name=request.session['job_id']+"/"+cname+"_SILAC_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_SILAC_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'SILACIP.html', context)

    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_4 = DocumentForm_4(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()

    job_id = generate_job_id()
    request.session["job_id"] = job_id
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'SILACIP.html', context)

def SILACwholeprot(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
    message_1 = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    result_file_name=''
    result_peptide_file_name = ''
    volcano_baseline_name=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb_vs_cntrl=[]
    condition_comb_no_cntrl=[]
    combination_choice=()
    column_norm_choice=''
    row_norm_choice=''
    fraction_choice=''
    column_norm_proteins=''
    imputation_choice=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id=None
    
    if ((request.method == 'POST')&(analysis_type=="CometLimma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # FWD rep
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # REV rep
        form_4 = DocumentForm_4(job_id,request.POST)
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES)
        form_prot_norm = Document_ProtNorm(request.POST)

        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processSILACwholeprotuploadComet(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,[],run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)

        # print("PASSED TO BACKEND:", bait_validation, form_3.is_valid(), submit_success)
        if (bait_validation==False)&(form_3.is_valid()):
            # print("Job submitted successfully? ", submit_success)
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'SILACwholeprot.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
            print(submit_message)
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            
            row_norm_choice=False#request.POST['rownorm']
            # if row_norm_choice == "yesrownorm":
            #     row_norm_choice = True
            # else:
            #     row_norm_choice = False
            

            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins=request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            fwd_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
            rev_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")


            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            # For COMET change the input file to have specific columns
            fwd_df = SILAC_COMETconversion(fwd_df,annotation_df,"FWD")
            rev_df = SILAC_COMETconversion(rev_df,annotation_df,"REV")

            bait_name=bait_name
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            pep_df = []
            prot_df,pep_df = choose_combination_process_SILAC(request.session["job_id"],"whole_proteome",analysis_type,combination_tmp[4], row_norm_choice, column_norm_choice,column_norm_proteins, bait_name, combination, r, fwd_df, rev_df, annotation_df, uniprotGene)
            print(prot_df.head(2))

            load_and_save_SILAC(prot_df, pep_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            if not prot_df.empty:
                run_success=True
                result_peptide_file_name=request.session['job_id']+"/"+cname+"_peptide_level_SILAC_FC_and_pval_YuLab.csv"
                result_file_name=request.session['job_id']+"/"+cname+"_SILAC_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_SILAC_FC_and_pval_full.pdf"

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'SILACwholeprot.html', context)

    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_4 = DocumentForm_4(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()

    job_id = generate_job_id()
    request.session["job_id"] = job_id
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success':submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'result_peptide_file_name':result_peptide_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'SILACwholeprot.html', context)

################################################################################
################# This section is for DIA-LFQ processing #######################
################################################################################

def DIALFQ(request):
    return render(request,'DIALFQ.html')

def ProtDIA(request):

    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
   
    message_1 = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_2=''
    validation_3=''

    run_success=False
    result_file_name=''
    volcano_baseline_name=''
    imputation_choice=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb=[]
    combination_choice=()
    column_norm_choice=''
    column_norm_proteins=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    job_id=None
    
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'

    if (request.method == 'POST')&((analysis_type=="LM")|(analysis_type=="Limma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # Bruker file input
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # DIANN complete raw input
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES) # annotation file
        form_prot_norm = Document_ProtNorm(request.POST)
        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processDIAuploadWholeProteome(request,bait_validation,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'ProtDIA.html', context)
        if bait_validation == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True
            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins==request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            init_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
            raw_init_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name=bait_name

            # TODO: if bait_name not NA then do bait normalization
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            prot_df = choose_combination_process_DIA(request.session['job_id'],"protein",analysis_type,'',column_norm_choice, column_norm_proteins, bait_name, combination, r, num_channels, init_df, raw_init_df, annotation_df, uniprotGene,imputation_user)
            ## Add # Precursors by all conditions
            ## Set # PSMs to be # Precursors of Treatment condition in comparison
            prot_df = add_Precursor_info(prot_df,raw_init_df,annotation_df,combination,r)
            print("Done calculating # Precursor information")
            print(prot_df.head(2))
            if "# PSMs" not in list(prot_df):
                print("Entered")
                prot_df["# PSMs"] = 5
            if "PSM Cutoff" not in list(prot_df):
                prot_df["PSM Cutoff"] = 5
            if "Gene Symbol" not in list(prot_df):
                if "Genes" in list(init_df):
                    prot_df["Gene Symbol"] = prot_df.Protein.map(init_df.set_index("Protein")["Genes"].to_dict())
                else:
                    prot_df["Proteinnoiso"] = prot_df["Protein"].str.split(";",expand=True)[0].str.split("-",expand=True)[0]
                    prot_df["Gene Symbol"] = prot_df.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
            
            print(list(prot_df))

            load_and_save(prot_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            #prot_df.to_excel(combination+"_FC_and_pval_YuLab_fractionxchannel.xlsx",index=False)
            if not prot_df.empty:
                run_success=True
                result_file_name=request.session['job_id']+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_FC_and_pval_full.pdf"
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'ProtDIA.html', context)


    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()


    print("First pass")
    job_id = generate_job_id()
    request.session['job_id'] = job_id
    print("Job id: ",job_id)
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'ProtDIA.html', context)

def PrecDIA(request):
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
   
    message_1 = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_2=''
    validation_3=''

    run_success=False
    result_file_name=''
    volcano_baseline_name=''
    imputation_choice=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb=[]
    combination_choice=()
    column_norm_choice=''
    column_norm_proteins=''

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    #_datetime = date.now()
    #datetime_str = _datetime.strftime("%Y-%m-%d-%H-%M-%S")
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'
    job_id = None

    if (request.method == 'POST')&((analysis_type=="LM")|(analysis_type=="Limma")):
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # Bruker file input
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # DIANN complete raw input
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES) # annotation file
        form_prot_norm = Document_ProtNorm(request.POST)
        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,volcano_baseline_name,result_file_name = processDIAuploadWholeProteome(request,bait_validation,analysis_type,form_1,form_2,form_3,form_prot_norm,run_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'PrecDIA.html', context)
        if bait_validation == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
        
            combination_choice=request.POST['combinations']
            column_norm_choice=request.POST['columnnorm']
            imputation_choice=request.POST['imputationchoice']
            imputation_user=False
            if imputation_choice == "impute":
                imputation_user=True
            if form_prot_norm.is_valid():
                column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
            # column_norm_proteins==request.POST['norm_by_uniprots']
            direction_choice=request.POST['Analysis Direction']

            ####-----Input for Shagun's Pipeline---------#####
            init_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
            raw_init_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")

            print('Annotation file: ', annotation_file)
            annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
            num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
            
            bait_name=bait_name

            # TODO: if bait_name not NA then do bait normalization
            
            combination_tmp=list(map(str.strip, combination_choice.split("'")))
            combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
            print (combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
            print (combination)

            cname = combination
            if int(direction_choice)==1:
                cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
            else:
                cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

            #combination="N-Vector"
            r=int(direction_choice)                # options 1 or -1 

            prot_df = []
            prot_df = choose_combination_process_DIA(request.session['job_id'], "precursor",analysis_type,'',column_norm_choice, column_norm_proteins, bait_name, combination, r, num_channels, init_df, raw_init_df, annotation_df, uniprotGene,imputation_user)
            ## Add # Precursors by all conditions
            ## Set # PSMs to be # Precursors of Treatment condition in comparison
            prot_df = add_Precursor_info(prot_df,raw_init_df,annotation_df,combination,r)
            print("Done calculating # Precursor information")
            print(prot_df.head(2))
            if "# PSMs" not in list(prot_df):
                print("Entered")
                prot_df["# PSMs"] = 5
            if "PSM Cutoff" not in list(prot_df):
                prot_df["PSM Cutoff"] = 5
            if "Gene Symbol" not in list(prot_df):
                if "Genes" in list(init_df):
                    prot_df["Gene Symbol"] = prot_df.Protein.map(init_df.set_index("Protein")["Genes"].to_dict())
                else:
                    prot_df["Proteinnoiso"] = prot_df["Protein"].str.split(";",expand=True)[0].str.split("-",expand=True)[0]
                    prot_df["Gene Symbol"] = prot_df.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
            
            print(list(prot_df))

            load_and_save(prot_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

            #prot_df.to_excel(combination+"_FC_and_pval_YuLab_fractionxchannel.xlsx",index=False)
            if not prot_df.empty:
                run_success=True
                result_file_name=request.session['job_id']+"/"+cname+"_FC_and_pval_YuLab.csv"
                volcano_baseline_name=request.session['job_id']+"/"+cname+"_FC_and_pval_full.pdf"
            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            # Render list page with the documents and the form
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'PrecDIA.html', context)


    else:
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()


    print("First pass")
    job_id = generate_job_id()
    request.session['job_id'] = job_id
    print("Job id: ",job_id)
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb': condition_comb, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'PrecDIA.html', context)

def IPDIA(request, job_id=None):

    # Need the excel with right columns
    # Need annotation file
    # ? SILAC need pairing information
    # Bait? (not for whole proteome)
    print(f"Great! You're using Python 3.6+. If you fail here, use the right version.")
   
    #message = 'Upload a tab separated output text file from Proteome Discoverer\n'
    message_1 = ''
    message_2 = ''
    message_3 = ''
    file_annotation_1=''
    file_annotation_2=''
    file_annotation_3=''

    validation_1=''
    validation_2=''
    validation_3=''

    run_success=False
    submit_success=False
    result_file_name=''
    volcano_baseline_name=''
    imputation_choice=''

    analysis_type = request.GET['analysisType']

    documents=''
    
    # Handle file upload
    
    bait_name=''
    bait_message=''
    annotation_file=''

    condition_comb_vs_cntrl=[]
    condition_comb_no_cntrl=[]
    combination_choice=()
    column_norm_choice=''
    column_norm_proteins=''

    progress=0

    direction_choice=''
    submit_message=''

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    local_storage_dir= BASE_DIR+'/media/documents/tmtipms_files'

    job_id = None
    if (request.method == 'POST')&((analysis_type=="LM")|(analysis_type=="Limma")):
        
        form_1 = DocumentForm_1(job_id,request.POST, request.FILES) # Bruker file input
        form_2 = DocumentForm_2(job_id,request.POST, request.FILES) # Precursor (complete info to get # Precursors)
        form_4 = DocumentForm_4(job_id,request.POST) # bait
        form_3 = DocumentForm_3(job_id,request.POST, request.FILES) # annotation file
        form_prot_norm = Document_ProtNorm(request.POST)
        bait_validation=False
        bait_validation,file_name_1,file_name_2,run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,volcano_baseline_name,result_file_name = processuploadIPDIA(request,bait_validation,analysis_type,form_1,form_2,form_3,form_4,form_prot_norm,[],run_success,submit_success,message_1,message_2,message_3,file_annotation_1,file_annotation_2,file_annotation_3,validation_1,validation_2,validation_3,bait_name,bait_message,annotation_file,condition_comb_vs_cntrl,condition_comb_no_cntrl,combination_choice,submit_message,BASE_DIR,local_storage_dir,volcano_baseline_name,result_file_name)
        
        if (bait_validation==False)&(form_3.is_valid()):

            print("Passing onto name:",annotation_file)
            print("Passing onto whether run successful:",run_success)
            print("Passing onto whether upload was successful:",submit_success)
            # Render list page with the documents and the form
            # 'progress':0, 
            context = {'job_id':request.session['job_id'],'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2,'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

            return render(request, 'IPDIA.html', context)
        if submit_success == True:
            submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
            
            if run_success == False:
            
                submit_message ='Your job has been succesfully submitted, you will be redirected to the result page once the analysis is complete'
                print("Job id after submitting successfully: ",job_id)
                combination_choice=request.POST['combinations']
                column_norm_choice=request.POST['columnnorm']
                imputation_choice=request.POST['imputationchoice']
                imputation_user=False
                if imputation_choice == "impute":
                    imputation_user=True
                if form_prot_norm.is_valid():
                    column_norm_proteins = form_prot_norm.cleaned_data["uniprot_to_normby"]
                # column_norm_proteins==request.POST['norm_by_uniprots']
                direction_choice=request.POST['Analysis Direction']

                ####-----Input for Shagun's Pipeline---------#####
                init_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_1),sep="\t")
                raw_input_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],file_name_2),sep="\t")
                # intensity_df = pd.read_csv("%s/%s"%(local_storage_dir,file_name_2),sep="\t")

                print('Annotation file: ', annotation_file)
                annotation_df = pd.read_csv("%s/%s"%(local_storage_dir+"/"+request.session['job_id'],annotation_file),sep=",")
                num_channels = len(list(set(annotation_df["Channel"].values.tolist())))
                
                bait_name=bait_name

                # TODO: if bait_name not NA then do bait normalization
                
                combination_tmp=list(map(str.strip, combination_choice.split("'")))
                combination="%s-%s"%(combination_tmp[1], combination_tmp[3])     #e.g. "N-Vector"
                print(combination_choice, combination_choice[1], combination_choice[2], list(map(str.strip, combination_choice.split("'"))))
                print(combination)

                cname = combination
                if int(direction_choice)==1:
                    cname = "%s-vs-%s"%(combination_tmp[1], combination_tmp[3])
                else:
                    cname = "%s-vs-%s"%(combination_tmp[3], combination_tmp[1])

                #combination="N-Vector"
                r=int(direction_choice)                # options 1 or -1 

                prot_df = []
                prot_df = choose_combination_process_DIA(request.session['job_id'],"IP",analysis_type,combination_tmp[4],column_norm_choice, column_norm_proteins, bait_name, combination, r, num_channels, init_df, raw_input_df, annotation_df, uniprotGene,imputation_user)
                
                ## Add # Precursors by all conditions
                ## Set # PSMs to be # Precursors of Treatment condition in comparison
                
                prot_df = add_Precursor_info(prot_df,raw_input_df,annotation_df,combination,r)
                print("Done calculating # Precursor information")

                print(prot_df.head(2))
                if "# PSMs" not in list(prot_df):
                    print("Entered")
                    prot_df["# PSMs"] = 5
                if "PSM Cutoff" not in list(prot_df):
                    prot_df["PSM Cutoff"] = 5
                if "Gene Symbol" not in list(prot_df):
                    if "Genes" in list(init_df):
                        prot_df["Gene Symbol"] = prot_df.Protein.map(init_df.set_index("Protein")["Genes"].to_dict())
                    else:
                        prot_df["Proteinnoiso"] = prot_df["Protein"].str.split(";",expand=True)[0].str.split("-",expand=True)[0]
                        prot_df["Gene Symbol"] = prot_df.Proteinnoiso.map(uniprotGene.set_index("Accession")["Gene Symbol"].to_dict())
                
        
                print(list(prot_df))

                load_and_save(prot_df, bait_name, BASE_DIR, cname, request.session['job_id'], 0.05)

                if not prot_df.empty:
                    run_success=True
                    result_file_name=request.session['job_id']+"/"+cname+"_FC_and_pval_YuLab.csv"
                    volcano_baseline_name=request.session['job_id']+"/"+cname+"_FC_and_pval_full.pdf"
                print("Passing onto name:",annotation_file)
                print("Passing onto whether run successful:",run_success)
                # Render list page with the documents and the form
                context = {'job_id':request.session['job_id'], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2,'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'submit_success': submit_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

                return render(request, 'IPDIA.html', context)


    else:
        
        form_1 = DocumentForm_1(job_id)  # An empty, unbound form
        form_2 = DocumentForm_2(job_id)  # An empty, unbound form
        form_3 = DocumentForm_3(job_id)  # An empty, unbound form
        form_4 = DocumentForm_4(job_id)  # An empty, unbound form
        form_prot_norm = Document_ProtNorm()



    print("First pass")
    job_id = generate_job_id()
    request.session['job_id'] = job_id
    print("Job id: ",job_id)
    print("Passing onto name:",annotation_file)
    print("Passing onto whether run successful:",run_success)
    # Render list page with the documents and the form
    context = {'job_id':request.session['job_id'], 'form_1': form_1, 'form_2': form_2, 'form_3': form_3, 'form_4': form_4, 'form_prot_norm':form_prot_norm, 'message_1': message_1, 'message_2': message_2, 'message_3': message_3, 'file_annotation_1': file_annotation_1, 'file_annotation_2': file_annotation_2,'file_annotation_3': file_annotation_3, 'validation_1': validation_1, 'validation_2': validation_2, 'validation_3': validation_3, 'condition_comb_vs_cntrl': condition_comb_vs_cntrl, 'condition_comb_no_cntrl':condition_comb_no_cntrl, 'submit_message': submit_message, 'run_success': run_success, 'analysis_type':analysis_type, 'result_file_name': result_file_name, 'volcano_baseline': volcano_baseline_name, 'annotation_file_name': annotation_file}

    return render(request, 'IPDIA.html', context)

# def long_running_task():
#     # Your long-running task here
#     for i in range(10):
#         # Do some work
#         # Update progress
#         progress = (i + 1) * 10
#         long_running_task.update_state(state='PROGRESS', meta={'progress': progress})
#     return 'Task completed successfully'

# def start_task(request):
#     task = long_running_task.delay()
#     return JsonResponse({'task_id': task.id})

# def get_task_progress(request):
#     task_id = request.GET.get('task_id')
#     task = long_running_task.AsyncResult(task_id)
#     if task.state == 'PROGRESS':
#         return JsonResponse({'progress': task.info['progress']})
#     elif task.state == 'SUCCESS':
#         return JsonResponse({'progress': 100})
#     else:
#         return JsonResponse({'progress': 0})

# Test view for progress bar
# def progressbar(request):
#     go_to_sleep.delay(5)
#     return render(request, 'progressbar.html')