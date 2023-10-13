library(maftools)

NENA_maf = read.maf(
  maf = "input01_NENA_MSS.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  clinicalData = "input01_NENA_MSS_Clinical.txt",
  isTCGA = FALSE)

EA_maf = read.maf(
  maf = "input01_EA_MSS.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  clinicalData = "input01_EA_MSS_Clinical.txt",
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

luadgenes = c("TP53","KRAS","ROS1","SETD2","ERBB2","ARID1A","NKX2-1",
              "RET","CDKN2A","ATM","STK11","APC","MET","PIK3CA","SMARCA4",
              "BRAF","KEAP1","RB1","RBM10","CTNNB1","ALK","EGFR")


cstage = c("#A4CBC9","#5CB0B8","#1493B4")
names(cstage) = c("a","b","c")
cstage = list(Stage=cstage)
cage = c("#ffbfa6","#ff8e5e","#ec612e","#a32f0c","#741200")
names(cage) = c("40s","50s","60s","70s","80s")
cage = list(Age=cage)
cgender = c("#f97c00","#391463")
names(cgender) = c("FEMALE","MALE")
cgender = list(Gender=cgender)
oncoplot(maf=NENA_maf, colors = vc_cols, genes = luadgenes,
         showTumorSampleBarcodes = F, draw_titv = T,
         clinicalFeatures = c("Stage","Age","Gender"),
         annotationColor = c(cstage,cage,cgender),
         keepGeneOrder =T, writeMatrix = T, removeNonMutated = F, logColBar = T)

oncoplot(maf=EA_maf, colors = vc_cols, genes = luadgenes,
         showTumorSampleBarcodes = F, draw_titv = T,
         clinicalFeatures = c("Stage","Age","Gender"),
         annotationColor = c(cstage,cage,cgender),
         keepGeneOrder =T, writeMatrix = T, removeNonMutated = F, logColBar = T)
