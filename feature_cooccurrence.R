library(maftools)


varfile_mss = read.maf(
  maf = "MAFs/var_EA_maftools.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  clinicalData = "MAFs/var_mss_clinicalData.txt",
  isTCGA = FALSE)

varfile_mss = read.maf(
  maf = "MAFs/var_NENA_maftools.maf",
  vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                "In_Frame_Ins","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","Fusion"),
  clinicalData = "MAFs/var_mss_clinicalData.txt",
  isTCGA = FALSE)


ce.stage = clinicalEnrichment(maf = varfile_mss,
                              clinicalFeature = 'EarlyLate')
ce.stage$groupwise_comparision[p_value < 0.05]

ce.age = clinicalEnrichment(maf = varfile_mss, 
                              clinicalFeature = 'YoungOld')
ce.age$groupwise_comparision[p_value < 0.05]

ce.gender = clinicalEnrichment(maf = varfile_mss, 
                            clinicalFeature = 'Gender')
ce.gender$groupwise_comparision[p_value < 0.01]

ce.status = clinicalEnrichment(maf = varfile_mss, 
                               clinicalFeature = 'TMBclass')
ce.status$groupwise_comparision[p_value < 0.01]


