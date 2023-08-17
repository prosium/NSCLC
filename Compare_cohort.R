library(maftools)

NSLA_EA = read.maf(
  maf = "EA_MSS_Mutation.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  isTCGA = FALSE)

Public_EA = read.maf(
  maf = "EA_Public.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  isTCGA = FALSE)

luadgenes=c("TP53","KRAS","SETD2","ERBB2",
            "ARID1A","NKX2-1","CDKN2A","ATM","STK11",
            "APC","MET","PIK3CA","SMARCA4","BRAF","KEAP1","RB1",
            "RBM10","CTNNB1","EGFR")

coOncoplot(m1 = NSLA_EA, m2 = Public_EA, 
           m1Name = 'NSLA_NENA', m2Name = 'Public_NENA', 
           genes = luadgenes, removeNonMutated = TRUE)
