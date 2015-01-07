# source("/home/gosia/R/R_Multinomial_project/Analysis_SimHSaV2/tmp.R")



###############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene.id <- as.factor(as.character(dge$genes$gene.id))


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, subset = Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-02, mcCores=20, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###############################################################
name1 <- "fc_g0_s3_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 3
dge <- dge[keep,]
dge$genes$gene.id <- as.factor(as.character(dge$genes$gene.id))


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-02, mcCores=30, verbose=FALSE) 

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=30)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=30)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))






###############################################################
name1 <- "fc_g0_s3_keep0s_subsetInf_DM4"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 3
dge <- dge[keep,]
dge$genes$gene.id <- as.factor(as.character(dge$genes$gene.id))


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-02, mcCores=30, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=30)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores=30)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




























