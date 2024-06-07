from django import forms
from .models import Document
from .views import *
#from django.contrib.admin.widgets import AdminFileWidget

class DocumentForm_3(forms.ModelForm):
    docfile_3 = forms.FileField(label='Select a file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_3','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_3, self).__init__(*args, **kwargs)
        # print("Initialized: ",job_id)
        self.fields['job_id'].initial = job_id  # Set initial value for job_id field
        # self.fields['docfile_3'].value = docfile_3

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_3, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_3')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_1(forms.ModelForm):
    docfile_1 = forms.FileField(label='Select a file')
    annotation_file = forms.CharField(label='File Name')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_1','annotation_file','job_id']  # Assuming your model has a 'file' field


    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_1, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_1, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_1')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_2(forms.ModelForm):
    docfile_2 = forms.FileField(label='Select a file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_2','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_2, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id 

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_2, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_2')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_Comet(forms.ModelForm):
    docfile_comet = forms.FileField(label='Select a COMET output file')
    annotation_file = forms.CharField(label='File Name')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_comet','annotation_file','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_Comet, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id 

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_Comet, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_comet')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_4(forms.ModelForm):
    bait_name = forms.CharField(label='Bait Name', max_length=100)
    annotation_file = forms.CharField(label='File Name')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['bait_name','annotation_file','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_4, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id 

    def save(self, commit=True):
        instance = super(DocumentForm_4, self).save(commit=False)
        # instance.job_id = self.job_id
        if commit:
            instance.save()
        return instance

class Document_ProtNorm(forms.Form):
    uniprot_to_normby = forms.CharField(label='Uniprots Normalize',required=False)

# class DocumentForm_ProtNorm(forms.ModelForm):
#     uniprot_to_normby = forms.CharField(label='Uniprots Normalize',required=False)
#     job_id = forms.CharField(max_length=100, required=False)
#     class Meta:
#         model = Document
#         fields = ['uniprot_to_normby','job_id']  # Assuming your model has a 'file' field

#     def __init__(self, *args, **kwargs):
#         super(DocumentForm_ProtNorm, self).__init__(*args, **kwargs)
#         self.fields['job_id'].initial = job_id

#     def save(self, commit=True):
#         instance = super(DocumentForm_ProtNorm, self).save(commit=False)
#         # instance.job_id = self.job_id
#         if commit:
#             instance.save()
#         return instance


class Document_FusionSeq(forms.Form):
    fusion_seq = forms.CharField(label='Fusion Sequence')


#####################################################################################
#####################################################################################
############## This section below is for Volcano plot generation ####################
#####################################################################################
#####################################################################################
class DocumentForm_5(forms.ModelForm):
    docfile_volcano = forms.FileField(label='Select a .csv file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_volcano','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_5, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_5, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_volcano')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance


#####################################################################################
#####################################################################################
############## This section below is for Network plot generation ####################
#####################################################################################
#####################################################################################

class DocumentForm_NetworkPlot(forms.Form):
    docfile_networkcsv = forms.FileField(label='Select a .csv file')

# class DocumentForm_NetworkPlot(forms.ModelForm):
#     docfile_networkcsv = forms.FileField(label='Select a .csv file')
#     job_id = forms.CharField(max_length=100, required=False)
#     class Meta:
#         model = Document
#         fields = ['docfile_networkcsv','job_id']  # Assuming your model has a 'file' field

#     def __init__(self, job_id, *args, **kwargs):
#         super(DocumentForm_NetworkPlot, self).__init__(*args, **kwargs)
#         self.fields['job_id'].initial = job_id

#     def save(self, commit=True, upload_path=None):
#         instance = super(DocumentForm_NetworkPlot, self).save(commit=False)
#         if upload_path:
#             instance.docfile.field.upload_to = upload_path
#         uploaded_file = self.cleaned_data.get('docfile_networkcsv')
#         instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
#         if commit:
#             instance.save()
#         return instance



#####################################################################################
#####################################################################################
##### This section below is for Mut vs WT filter (assuming two controls used) #######
#####################################################################################
#####################################################################################

class DocumentForm_FilterFileForIP(forms.ModelForm):
    docfile_filterfileforip = forms.FileField(label='Select a .csv file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_filterfileforip','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_FilterFileForIP, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_FilterFileForIP, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_filterfileforip')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_UniprotFile(forms.ModelForm):
    docfile_uniprotfile = forms.FileField(label='Select a .csv file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_uniprotfile','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_UniprotFile, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_UniprotFile, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_uniprotfile')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance

class DocumentForm_UniprotFile2(forms.ModelForm):
    docfile_uniprotfile2 = forms.FileField(label='Select a .csv file')
    job_id = forms.CharField(max_length=100, required=False)
    class Meta:
        model = Document
        fields = ['docfile_uniprotfile2','job_id']  # Assuming your model has a 'file' field

    def __init__(self, job_id, *args, **kwargs):
        super(DocumentForm_UniprotFile2, self).__init__(*args, **kwargs)
        self.fields['job_id'].initial = job_id

    def save(self, commit=True, upload_path=None):
        instance = super(DocumentForm_UniprotFile2, self).save(commit=False)
        if upload_path:
            instance.docfile.field.upload_to = upload_path
        uploaded_file = self.cleaned_data.get('docfile_uniprotfile2')
        instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
        if commit:
            instance.save()
        return instance



#####################################################################################
#####################################################################################
########## This section below is for calculating adjusted pvalue ####################
#####################################################################################
#####################################################################################


# class DocumentForm_AdjustedPvalue(forms.Form):
#     docfile_adjustedpvalue = forms.FileField(label='Select a .csv file')
#     job_id = forms.CharField(max_length=100, required=False)
#     class Meta:
#         model = Document
#         fields = ['docfile_adjustedpvalue','job_id']  # Assuming your model has a 'file' field

#     def __init__(self, job_id, *args, **kwargs):
#         super(DocumentForm_AdjustedPvalue, self).__init__(*args, **kwargs)
#         self.fields['job_id'].initial = job_id

#     def save(self, commit=True, upload_path=None):
#         instance = super(DocumentForm_AdjustedPvalue, self).save(commit=False)
#         if upload_path:
#             instance.docfile.field.upload_to = upload_path
#         uploaded_file = self.cleaned_data.get('docfile_adjustedpvalue')
#         instance.docfile.save(uploaded_file.name, uploaded_file, save=False)  # Save the file
#         if commit:
#             instance.save()
#         return instance

# class DocumentForm_NetworkPlot(forms.Form):
#     docfile_networkcsv = forms.FileField(label='Select a .csv file')

# class DocumentForm_1(forms.Form):
#     docfile_1 = forms.FileField(label='Select a file')
#     annotation_file = forms.CharField(label='File Name')

# class DocumentForm_2(forms.Form):
#     docfile_2 = forms.FileField(label='Select a file')

# class DocumentForm_Comet(forms.Form):
#     docfile_comet = forms.FileField(label='Select a COMET output file')


# class DocumentForm_3(forms.Form):
#     docfile_3 = forms.FileField(label='Select a file')

# class DocumentForm_4(forms.Form):
#     bait_name = forms.CharField(label='Bait Name', max_length=100)
#     annotation_file = forms.CharField(label='File Name')



# class DocumentForm_5(forms.Form):
#     docfile_volcano = forms.FileField(label='Select a .csv file')

# class DocumentForm_FilterFileForIP(forms.Form):
#     docfile_filterfileforip = forms.FileField(label='Select a .csv file')

# class DocumentForm_UniprotFile(forms.Form):
#     docfile_uniprotfile = forms.FileField(label='Select a .csv file')