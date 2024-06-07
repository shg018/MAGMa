You will have two PSM search files from PD that has S/N and reporter ion intensity values per PSM. There should be a way to know which PSM list has which (in the name of the file perhaps) since PD would name the columns the same. For example in the files provided, SNUsed means S/N values were provided by PD and IntensitiesUsed means raw reporter ion intensities were provided by PD.

1) The following columns are needed from the PSM list files (if search done in PD):

Intensity
Annotated Sequence
Sequence
Master Protein Accessions
Protein Accessions
Isolation Interference [%]
Spectrum File
Charge
RT [min]
# Protein Groups
Abundance: <Channel> - channel corresponds to all the TMT ID's (126, 127, 128) etc.

The second file needed is the annotation.csv
2) This file should have the following two columns for now:
 - Channel (the TMT tag number - 126, 127 etc)
 - Label (the condition being tested with that channel)

3) The user should provide the name of the bait protein (preferably UniProt ID/ Gene Symbol) if running an IP-TMT sample. Bait for the example files provided would be "N","M","orf3a","orf6" if user inputs gene symbol.

4) The user should also provide the condition comparisons to be tested and the names should match the labels provided in the Label column of the annotation.csv file.


