import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import random


file_mutation = "input02_SMG.txt"
file_proteome = "input02_Protein.txt"

Genes_DB_CAG = ["ABI1","ABL1","ACKR3","ACSL3","ACVR1","ACVR2A","AFF1","AFF3","AFF4","AKT1","ALK","AMER1","APC","APOBEC3B","AR","ARHGAP26","ARID1A","ARID1B","ARID2","ARNT","ASPSCR1","ASXL1","ATF1","ATIC","ATM","ATP1A1","ATP2B3","ATR","ATRX","AXIN1","AXIN2","B2M","BAP1","BARD1","BCL10","BCL11A","BCL11B","BCL9","BCL9L","BCOR","BCORL1","BIRC3","BLM","BMPR1A","BRAF","BRCA1","BRCA2","BRD4","BRIP1","BTK","BUB1B","CACNA1D","CALR","CAMTA1","CANT1","CARD11","CARS","CASP8","CBFA2T3","CBFB","CBL","CBLB","CCDC6","CCNB1IP1","CCND1","CCND2","CCNE1","CD79A","CD79B","CDC73","CDH1","CDH11","CDK12","CDK4","CDK6","CDKN1B","CDKN2A","CDKN2C","CEBPA","CHD4","CHEK2","CIC","CIITA","CLIP1","CLTC","CLTCL1","CNBP","CNOT3","COL2A1","CREB3L1","CREB3L2","CREBBP","CRLF2","CRTC1","CSF3R","CTCF","CTNNB1","CUX1","CXCR4","CYLD","DAXX","DDB2","DDIT3","DDR2","DDX10","DDX3X","DDX5","DDX6","DICER1","DNM2","DNMT3A","DROSHA","EBF1","EGFR","EIF3E","EIF4A2","ELF4","ELK4","ELL","EML4","EP300","EPAS1","EPS15","ERBB2","ERBB3","ERBB4","ERC1","ERCC2","ERCC3","ERCC4","ERCC5","ERG","ESR1","ETNK1","ETV6","EWSR1","EZH2","EZR","FANCD2","FAS","FAT1","FAT4","FBXO11","FBXW7","FCGR2B","FES","FGFR1","FGFR2","FGFR3","FGFR4","FHIT","FIP1L1","FLT3","FLT4","FOXA1","FOXL2","FUBP1","FUS","GAS7","GATA1","GATA2","GATA3","GNA11","GNAQ","GNAS","GPC3","GRIN2A","H3F3A","H3F3B","HEY1","HIF1A","HIP1","HIST1H3B","HNF1A","HNRNPA2B1","HOXA11","HRAS","IDH1","IDH2","IGH","IKBKB","IKZF1","IL6ST","IL7R","IRS4","JAK1","JAK2","JAK3","KDM5C","KDM6A","KDR","KEAP1","KIT","KLF4","KMT2A","KMT2C","KMT2D","KNL1","KRAS","LATS1","LATS2","LCK","LEF1","LIFR","LMNA","LRP1B","LZTR1","MAP2K1","MAP2K2","MAP2K4","MAP3K1","MAP3K13","MAPK1","MAX","MED12","MEN1","MET","MLH1","MLLT10","MPL","MSH2","MSH6","MTOR","MYC","MYD88","MYOD1","NAB2","NCOA2","NCOR1","NCOR2","NDRG1","NF1","NF2","NFE2L2","NFKBIE","NONO","NOTCH1","NOTCH2","NPM1","NRAS","NT5C2","NTRK1","NTRK3","PAFAH1B2","PAX3","PAX5","PAX7","PAX8","PBRM1","PBX1","PDCD1LG2","PDGFRA","PER1","PHF6","PHOX2B","PICALM","PIK3CA","PIK3CB","PIK3R1","PIM1","PLCG1","PMS2","POLD1","POLE","POT1","POU2AF1","POU5F1","PPARG","PPM1D","PPP2R1A","PPP6C","PRDM1","PRDM16","PREX2","PRKACA","PRKAR1A","PSIP1","PTCH1","PTEN","PTK6","PTPN11","PTPN13","PTPRB","PTPRT","QKI","RABEP1","RAC1","RAD21","RAF1","RANBP2","RAP1GDS1","RB1","RBM10","RBM15","RECQL4","RET","RHOA","RHOH","RNF213","RNF43","ROS1","RPL5","RUNX1T1","SALL4","SDHA","SF3B1","SFPQ","SIX1","SLC34A2","SLC45A3","SMAD2","SMAD3","SMAD4","SPOP","STAG2","SUZ12","TBL1XR1","TCF7L2","TCL1A","TERT","TET1","TET2","TFE3","TGFBR2","TMPRSS2","TNFAIP3","TP53","TP63","TRAF7","UBR5","USP8","VHL","WAS","WIF1","XPO1","ZBTB16","ZFHX3","ZRSR2","STK11","ARID1A","SETD2","ROS1","ERBB2","RET","NKX2-1","KRAS","CDKN2A"]
Genes_DB_CAG = list(set(Genes_DB_CAG))

def average(list):
    return (sum(list) / len(list))


def pvalue_categorization(pvalue):
    
    if 0.01 < pvalue <= 0.05:
        results_pvalue = "(0.01, 0.05]"
        
    
    elif 0.001 < pvalue <= 0.01:
        results_pvalue = "(0.001, 0.01]"
        

    elif pvalue <= 0.001:
        results_pvalue = "<0.001"
        

    else:
        results_pvalue = "n.s."
    
    return results_pvalue


def getProteom(mutgene, Proteinfilename, group_mut, group_nomut):
        
    inputfile = open(Proteinfilename)
    header = inputfile.readline().strip().split("\t")
    final_bag = []
    group_mut_pr_idx = [header.index(i) for i in header if i in group_mut]
    group_nomut_pr_idx = [header.index(i) for i in header if i in group_nomut]
    

    for line in inputfile:
        elms = line.strip().split("\t")
        
        gene_id_pr = elms[0]
        
        if gene_id_pr in Genes_DB_CAG:
                
            group_mut_pr_values = list(map(float,[elms[i] for i in group_mut_pr_idx if elms[i]!="NA"]))
            group_nomut_pr_values = list(map(float,[elms[i] for i in group_nomut_pr_idx if elms[i]!="NA"]))
            
            if len(group_mut_pr_values) != 0 and len(group_nomut_pr_values) !=0:
                
                # cis or trans
                if mutgene == gene_id_pr:
                    results_CT = "cis"
                    
                else:
                    results_CT = "trans"
                
                # foldchange
                fc = average(group_mut_pr_values) - average(group_nomut_pr_values)
                if 1.5 > fc >=1:
                    results_fc = ">1"

                elif fc > 1.5 :
                    results_fc = ">1.5"
                
                elif -1.5 < fc <= -1:
                    results_fc = "<-1"
                
                elif fc < -1.5:
                    results_fc = "<-1.5"
                

                else:
                    results_fc = "n.s."


                # pvalue
                pvalue = stats.ranksums(group_mut_pr_values, group_nomut_pr_values)[1]
                results_pvalue = pvalue_categorization(pvalue)
                
                if results_pvalue != "n.s." and results_fc != "n.s.":
                    # print(final_cnt)

                    final_results = mutgene, gene_id_pr, results_CT, str(results_fc), results_pvalue
                    final_bag.append(":::".join(final_results))

                    
    # print(final_bag)
    inputfile.seek(0)
    return final_bag
    

LUAD_mutation = open(file_mutation)
header = LUAD_mutation.readline().strip().split("\t")

df_All = pd.DataFrame(
    {
    'Gene':[],
    'target':[],
    'effect':[],
    'Z-score':[],
    'p-value':[],
    })


for line in LUAD_mutation:
    elms = line.strip().split("\t")
    genename = elms[0]
    
    
    print(genename + " is processed...")
    binary_mut = elms[1:]
    
    samples_mut = [header[i[0]+1] for i in enumerate(binary_mut) if i[1] =="1"]
    samples_wt = [header[i[0]+1] for i in enumerate(binary_mut) if i[1] =="0"]
    

    sigPr = getProteom(genename, file_proteome, samples_mut, samples_wt)
    
    for i in sigPr:
        elms = i.split(":::")
        
        
        df_All = df_All.append(
            {
            'Gene':elms[0],
            'target':elms[1],
            'effect':elms[2],
            'Z-score':elms[3],
            'p-value':elms[4],
            }, ignore_index=True)



##### variable setting

markers_style = {"cis":"o", "trans":"s"} # setting shape of cis or trans 
kwargs  =   {'edgecolor':"black",'linewidth':0.3, 'linestyle':'-'} # setting outline of shape

##### plotting

manual_palette = {">1.5":"#E21A1B",">1":"#F7B6AC","<-1":"#B3CCE0","<-1.5":"#397CB7"}

g = sns.relplot(data=df_All.sort_values(by='target', ascending=False), x='target', y='Gene', 
            style='effect', markers=markers_style, style_order=['cis','trans'],
            hue = 'Z-score', hue_order=[">1.5",">1","<-1","<-1.5"], palette=manual_palette,
            size='p-value', size_order=["<0.001","(0.001, 0.01]","(0.01, 0.05]"], sizes=[140, 100, 40],
            **kwargs, 
            legend="full",
            alpha=0.6, height=4, aspect=9/5)
            


##### Plot setting
# plt size


# border line and grid
sns.despine(fig=None, ax=None, top=False, 
            right=False, left=False, bottom=False, 
            offset=None, trim=False)
plt.grid()
g.ax.set_axisbelow(True)
  
# axis setting
g.ax.set_title("Mutational Impact on Cancer-associated Protein Expression", fontsize=14)
g.ax.set_xlabel("Proteom change abundances", fontsize = 12)
g.ax.set_ylabel("Mutated Genes", fontsize = 12)
g.set_xticklabels(rotation=90) 

# export to svg
plt.savefig(fname="results.svg", bbox_inches='tight')