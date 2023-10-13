import os
import pandas as pd
from scipy import stats
import math
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

file_mutation = "input03_SMG.txt"
file_transcriptome = "input03_Transcriptome.txt"
file_proteome = "input03_Protein.txt"
file_phosphoproteome = "input03_Phosphoproteome.txt"

TargetGenes = "ERBB2"


outputtemp = open("temp.file","w")
header = "Origin","Target","Difference","minuslogPvalue"
outputtemp.write("\t".join(header)+"\n")

inputfile_mutation = open(file_mutation)
header = inputfile_mutation.readline().strip().split("\t")

def average(list):
    return (sum(list) / len(list))

def getTranscriptome(mRNAfilename, group_mut, group_nomut):
    inputfile = open(mRNAfilename)
    header = inputfile.readline().strip().split("\t")

    group_mut_tr_idx = [header.index(i) for i in header if i in group_mut]
    group_nomut_tr_idx = [header.index(i) for i in header if i in group_nomut]

    for line in inputfile:
        elms = line.strip().split("\t")
        gene_id_tr = elms[0]

        group_mut_tr_values = list(map(float, [elms[i] for i in group_mut_tr_idx if elms[i]!="NA"]))
        group_nomut_tr_values = list(map(float, [elms[i] for i in group_nomut_tr_idx if elms[i]!="NA"]))
        
        if group_mut_tr_values.count("NA") != len(group_mut_tr_values) and group_nomut_tr_values.count("NA") != len(group_nomut_tr_values):

            pvalue = stats.ranksums(group_mut_tr_values, group_nomut_tr_values)[1]
            foldchange = average(group_mut_tr_values) - average(group_nomut_tr_values)

            e1 = TargetGenes+"_Tr",gene_id_tr,str(foldchange),str(pvalue)
            outputtemp.write("\t".join(e1)+"\n")
        
def getProteome(proteinfilename, group_mut, group_nomut):
    inputfile = open(proteinfilename)
    header = inputfile.readline().strip().split("\t")

    group_mut_pr_idx = [header.index(i) for i in header if i in group_mut]
    group_nomut_pr_idx = [header.index(i) for i in header if i in group_nomut]

    for line in inputfile:
        elms = line.strip().split("\t")
        gene_id_pr = elms[0]

        group_mut_pr_values = list(map(float, [elms[i] for i in group_mut_pr_idx if elms[i]!="NA"]))
        group_nomut_pr_values = list(map(float, [elms[i] for i in group_nomut_pr_idx if elms[i]!="NA"]))
        
        if group_mut_pr_values.count("NA") != len(group_mut_pr_values) and group_nomut_pr_values.count("NA") != len(group_nomut_pr_values):
            pvalue = stats.ranksums(group_mut_pr_values, group_nomut_pr_values)[1]
            foldchange = average(group_mut_pr_values) - average(group_nomut_pr_values)

            e2 = TargetGenes+"_Pr",gene_id_pr,str(foldchange),str(pvalue)
            outputtemp.write("\t".join(e2)+"\n")

def getPhosphoproteome(phosphoproteinfilename, group_mut, group_nomut):
    inputfile = open(phosphoproteinfilename)
    header = inputfile.readline().strip().split("\t")

    group_mut_ppr_idx = [header.index(i) for i in header if i in group_mut]
    group_nomut_ppr_idx = [header.index(i) for i in header if i in group_nomut]

    for line in inputfile:
        elms = line.strip().split("\t")
        gene_id_ppr = elms[0]

        group_mut_ppr_values = list(map(float, [elms[i] for i in group_mut_ppr_idx if elms[i]!="NA"]))
        group_nomut_ppr_values = list(map(float, [elms[i] for i in group_nomut_ppr_idx if elms[i]!="NA"]))
        
        if group_mut_ppr_values.count("NA") != len(group_mut_ppr_values) and group_nomut_ppr_values.count("NA") != len(group_nomut_ppr_values):

            pvalue = stats.ranksums(group_mut_ppr_values, group_nomut_ppr_values)[1]
            foldchange = average(group_mut_ppr_values) - average(group_nomut_ppr_values)

            e3 = TargetGenes+"_Ppr",gene_id_ppr,str(foldchange),str(pvalue)
            outputtemp.write("\t".join(e3)+"\n")


def getAcetylproteome(acetylproteomefilename, group_mut, group_nomut):
    inputfile = open(acetylproteomefilename)
    header = inputfile.readline().strip().split("\t")

    group_mut_ac_idx = [header.index(i) for i in header if i in group_mut]
    group_nomut_ac_idx = [header.index(i) for i in header if i in group_nomut]

    for line in inputfile:
        elms = line.strip().split("\t")
        gene_id_ac = elms[0]

        group_mut_ac_values = list(map(float, [elms[i] for i in group_mut_ac_idx if elms[i]!="NA"]))
        group_nomut_ac_values = list(map(float, [elms[i] for i in group_nomut_ac_idx if elms[i]!="NA"]))
        
        if group_mut_ac_values.count("NA") != len(group_mut_ac_values) and group_nomut_ac_values.count("NA") != len(group_nomut_ac_values):

            pvalue = stats.ranksums(group_mut_ac_values, group_nomut_ac_values)[1]
            foldchange = average(group_mut_ac_values) - average(group_nomut_ac_values)

            e4 = TargetGenes+"_Ac",gene_id_ac,str(foldchange),str(pvalue)
            outputtemp.write("\t".join(e4)+"\n")




for line in inputfile_mutation:
    elms = line.strip().split("\t")
    if elms[0] == TargetGenes:
        values = elms[1:]

        group_mut = [header[i[0]+1] for i in enumerate(values) if i[1]=="1"]
        group_nomut = [header[i[0]+1] for i in enumerate(values) if i[1]=="0"]
            

        getTranscriptome(file_transcriptome, group_mut, group_nomut)
        getProteome(file_proteome, group_mut, group_nomut)
        getPhosphoproteome(file_phosphoproteome, group_mut, group_nomut)
        # getAcetylproteome(file_acetylproteome, group_mut, group_nomut)
        

tempfile = open("temp.file")
outputtemp2 = open("temp2.file","w")
header2 = "TargetGene","Type","SignedP_mRNA","SignedP_type"
outputtemp2.write("\t".join(header2)+"\n")



mRNA_dict = {}
mRNA_bag = []

tempfile.readline()
for line in tempfile:
    elms = line.strip().split("\t")
    criteria  = elms[0].split("_")[1]


    if criteria == "Tr":
        mRNA_bag.append(elms[1])

        if float(elms[2]) < 0 :
            signedP = abs(-1 * math.log(float(elms[3]),10)) * -1
            mRNA_dict[elms[1]] = signedP

        else:
            signedP = abs(-1 * math.log(float(elms[3]),10))
            mRNA_dict[elms[1]] = signedP

tempfile.seek(0)
tempfile.readline()

for line in tempfile:
    elms = line.strip().split("\t")
    criteria1 = elms[0].split("_")[1]

    if criteria1 == "Pr":
        criteria2 = elms[1] in mRNA_bag

        if criteria2 == True:
            if float(elms[2]) < 0:
                signedP_Prt = str(abs(-1 * math.log(float(elms[3]),10)) * -1)
                e1 = elms[1]+"\t"+"Protein"+"\t"+str(mRNA_dict[elms[1]])+"\t"+signedP_Prt+"\n"
                outputtemp2.write(e1)

            if float(elms[2]) > 0:
                signedP_Prt = str(abs(-1 * math.log(float(elms[3]),10)) * 1)
                e1 = elms[1]+"\t"+"Protein"+"\t"+str(mRNA_dict[elms[1]])+"\t"+signedP_Prt+"\n"
                outputtemp2.write(e1)

    if criteria1 == "Ppr":
        criteria3 = elms[1].split(":")[0] in mRNA_bag


        if criteria3 == True:
            if float(elms[2]) < 0:
                signedP_Phs = str(abs(-1 * math.log(float(elms[3]),10)) * -1)
                e2 = elms[1]+"\t"+"Phosphoprotein"+"\t"+str(mRNA_dict[elms[1].split(":")[0]])+"\t"+signedP_Phs+"\n"
                outputtemp2.write(e2)

            if float(elms[2]) > 0:
                signedP_Phs = str(abs(-1 * math.log(float(elms[3]),10)) * 1)
                e2 = elms[1]+"\t"+"Phosphoprotein"+"\t"+str(mRNA_dict[elms[1].split(":")[0]])+"\t"+signedP_Phs+"\n"
                outputtemp2.write(e2)

    if criteria1 == "Ac":
        criteria3 = elms[1].split(":")[0] in mRNA_bag
        if criteria3 == True:
            if float(elms[2]) < 0:
                signedP_Phs = str(abs(-1 * math.log(float(elms[3]),10)) * -1)
                e2 = elms[1]+"\t"+"Acetylprotein"+"\t"+str(mRNA_dict[elms[1].split(":")[0]])+"\t"+signedP_Phs+"\n"
                outputtemp2.write(e2)

            if float(elms[2]) > 0:
                signedP_Phs = str(abs(-1 * math.log(float(elms[3]),10)) * 1)
                e2 = elms[1]+"\t"+"Acetylprotein"+"\t"+str(mRNA_dict[elms[1].split(":")[0]])+"\t"+signedP_Phs+"\n"
                outputtemp2.write(e2)

