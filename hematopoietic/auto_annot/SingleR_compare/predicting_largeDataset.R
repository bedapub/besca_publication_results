rm(list=ls())

library(SingleR)


SE_Granja <- readRDS("Granja2019_annotated.Rds")


SE_KV <-  readRDS("Kotliarov2020_processed_citeseq_merged_annotated.Rds")

SE_Granja$celltype3_original <- SE_Granja$celltype3
SE_KV$celltype3_original <- SE_KV$celltype3



pred.Granja_trainKV <- SingleR(test=SE_Granja, ref=SE_KV, labels=SE_KV$celltype3, de.method="wilcox")
write.csv( pred.Granja_trainKV, file =  "pred.Granja_trainKV.csv")


pred.KV_trainGranja  <- SingleR(test=SE_KV, ref=SE_Granja, labels=SE_Granja$celltype3, de.method="wilcox")


write.csv( pred.KV_trainGranja, file =  "pred.KV_trainGranja.csv")
