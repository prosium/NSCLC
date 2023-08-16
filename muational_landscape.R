library(maftools)

MSS_NENA = read.maf(
  maf = "MSS_NENA.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  isTCGA = FALSE)

MSS_EA = read.maf(
  maf = "MSS_EA.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  isTCGA = FALSE)

vc_cols = c(
  "#377eb8","#984ea3","#ffff33",
  "#a65628","#e41a1c","#ff7f00","#4daf4a",
  "#8a8628",
  "black","#1e94bf",
  "#cd1c66","#3340bd")

names(vc_cols) = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                   "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del",
                   "Fusion",
                   "Multi_Hit", "Complex_Event",
                   "Amp","Del")

luadgenes=c("TP53","KRAS","ROS1","SETD2","ERBB2",
            "ARID1A","NKX2-1","RET","CDKN2A","ATM","STK11",
            "APC","MET","PIK3CA","SMARCA4","BRAF","KEAP1","RB1",
            "RBM10","CTNNB1","ALK","EGFR")

oncoplot(maf=MSS_EA, colors = vc_cols, genes = luadgenes,
         showTumorSampleBarcodes = F, draw_titv = TRUE, 
         keepGeneOrder =T, writeMatrix = T, removeNonMutated = F)

