suppressPackageStartupMessages({
	library(edgeR)
	library(DESeq2)
	})

args = commandArgs(trailingOnly=TRUE)

# Add in real counts file, lists of selected samples for C0, C1
if (length(args)<5) {
  stop("Usage: sigtest.R [cond0 expr table] [cond1 expr table] [real counts file] [outname] [op:use DESeq? (0/1)] 
       [cond0 training samples] [cond1 training samples]", 
       call.=FALSE)
}

#Load synthetic counts files
testCond0 = read.csv(args[1], row.names = 1)
colnames(testCond0) <- sprintf("samp%s",1:ncol(testCond0))
testCond1 = read.csv(args[2], row.names = 1) #lung_runs/test/samp10/run1/lung_synth_c0.csv  lung_runs/test/test_cond_1.csv
colnames(testCond1) <- sprintf("samp%s",(ncol(testCond0)+1):(ncol(testCond0)+ncol(testCond1)))

#Load real counts file (a subset of which were used for model training)
readCounts = F
if(args[3] != '0')
{  
  readCounts = read.csv(args[3], check.names = T)
  # readCounts = readCounts[!duplicated(readCounts[,1]),]
  # readCounts = readCounts[!is.na(readCounts[,1]),]
  rnames = readCounts[2:nrow(readCounts), 1]
  sampleInfo = data.frame(condition = as.numeric(readCounts[1, 2:ncol(readCounts)]), row.names = colnames(readCounts)[2:ncol(readCounts)])
  sampleInfo$condition = as.factor(sampleInfo$condition)
  readCounts = readCounts[2:nrow(readCounts), 2:ncol(readCounts)]
  rownames(readCounts) = rnames
  
  c0Samps = strsplit(args[6],',')[[1]]
  c1Samps = strsplit(args[7],',')[[1]]
}

fname = args[4] #output file name prefix

useDESeq = '0'
if(length(args) > 4)
{
  useDESeq = args[5] #Whether to use DESeq2 or edgeR for DE analysis
}


testInput = cbind(testCond0, testCond1)
group=c(rep('0',ncol(testCond0)), rep('1',ncol(testCond1)))

if(useDESeq == '1')
{
  print('Using DESeq2')
  #DE for synthetic samples
  tmp = testInput
  tmp[] = lapply(tmp, as.integer)
  synthSampleInfo = data.frame(condition = as.factor(group),
                               row.names=colnames(tmp))
  des.ds = DESeqDataSetFromMatrix(countData = tmp, colData = synthSampleInfo, design = ~condition)
  #remove genes with low counts
  #des.ds = des.ds[ rowSums(counts(des.ds)) > 10, ]
  des.ds = estimateSizeFactors(des.ds)
  des.ds = DESeq(des.ds)
  des.results = results(des.ds, independentFiltering=F)
  #Calc abs difference in average b/w conditions
  allCounts = counts(des.ds, normalized = F)
  absDiff = abs(rowSums(allCounts[,group==0])/sum(group==0) - 
                rowSums(allCounts[,group==1])/sum(group==1))
  deRank = order(des.results$padj)
  write.table(cbind(des.results, absDiff, deRank), file = paste0(fname,'_tmp_pval.csv'), sep=',', 
              row.names = T, col.names = NA, quote=F)
  
  if(type(readCounts) != 'logical')
  {
    tmp = readCounts
    tmp[] = lapply(tmp, as.integer)
    des.ds = DESeqDataSetFromMatrix(countData = tmp, colData = sampleInfo, design = ~condition)
    #remove genes with low counts
    #des.ds = des.ds[ rowSums(counts(des.ds)) > 10, ]
    des.ds = estimateSizeFactors(des.ds)
    des.ds = DESeq(des.ds)
    des.results = results(des.ds, independentFiltering=F)
    #Calc abs difference in average b/w conditions
    allCounts = counts(des.ds, normalized = F)
    absDiff = abs(rowSums(allCounts[,sampleInfo$condition==0])/sum(sampleInfo$condition==0) - 
                    rowSums(allCounts[,sampleInfo$condition==1])/sum(sampleInfo$condition==1))
    write.table(cbind(des.results, absDiff), file = paste0(fname,'_tmp_true_pval.csv'), sep=',', 
                row.names = T, col.names = NA, quote=F)
    
    #Calc DE for just the training samples:
    sampleInfoTrain = sampleInfo[make.names(c(c0Samps,c1Samps)), , drop=F]
    des.ds = DESeqDataSetFromMatrix(countData = tmp[ , make.names(c(c0Samps,c1Samps))], 
                                    colData = sampleInfoTrain, design = ~condition)
    des.ds = estimateSizeFactors(des.ds)
    des.ds = DESeq(des.ds)
    des.results = results(des.ds, independentFiltering=F)
    #Calc abs difference in average b/w conditions
    allCounts = counts(des.ds, normalized = F)
    absDiff = abs(rowSums(allCounts[,sampleInfoTrain$condition==0])/sum(sampleInfoTrain$condition==0) - 
                    rowSums(allCounts[,sampleInfoTrain$condition==1])/sum(sampleInfoTrain$condition==1))
    write.table(cbind(des.results, absDiff), file = paste0(fname,'_training_pval.csv'), sep=',', 
                row.names = T, col.names = NA, quote=F)
    
    #Calc DE for the training samples duplicated to size of synthetic cohorts:
    matExprAlt = tmp[ , make.names(c(c0Samps,c1Samps))]
    matExprBig = matExprAlt
    synthSize = nrow(synthSampleInfo)
    numTimesDup = floor(synthSize/length(c(c0Samps,c1Samps))) - 1
    for (i in seq(numTimesDup))
    {
      tmpMat = matExprAlt
      newNames = c()
      for (j in colnames(matExprAlt))
      {
        newNames = c(newNames, paste0(j,".",i))
      }
      colnames(tmpMat) = newNames
      matExprBig = cbind(matExprBig, tmpMat)
    }
    bigGroups = c(rep('0',length(c0Samps)), rep('1',length(c1Samps)))
    bigGroups = rep(bigGroups, numTimesDup+1)
    bigSampInfo = data.frame(condition = bigGroups, row.names = colnames(matExprBig))
    
    des.ds = DESeqDataSetFromMatrix(countData = matExprBig, 
                                    colData = bigSampInfo, design = ~condition)
    des.ds = estimateSizeFactors(des.ds)
    des.ds = DESeq(des.ds)
    des.results = results(des.ds, independentFiltering=F)
    #Calc abs difference in average b/w conditions
    allCounts = counts(des.ds, normalized = F)
    absDiff = abs(rowSums(allCounts[,bigSampInfo$condition==0])/sum(bigSampInfo$condition==0) - 
                    rowSums(allCounts[,bigSampInfo$condition==1])/sum(bigSampInfo$condition==1))
    write.table(cbind(des.results, absDiff), file = paste0(fname,'_trainingDuplicated_pval.csv'), sep=',', 
                row.names = T, col.names = NA, quote=F)
  }
  
} else
{
  print('Using edgeR')
  de = DGEList(counts=testInput, genes=rownames(testInput), group=group)
  de = calcNormFactors(de)
  design <- model.matrix(~group)
  de = estimateDisp(de, design)
  fit <- glmQLFit(de, design)
  res <- glmQLFTest(fit,coef=2) 
  res <- topTags(res, sort.by="none", n=length(de$genes$genes))
  res = data.frame(res)
  #Calc abs difference in average b/w conditions
  allCounts = de$counts
  absDiff = abs(rowSums(allCounts[,group==0])/sum(group==0) - 
                rowSums(allCounts[,group==1])/sum(group==1))
  deRank = order(res$FDR)
  write.table(cbind(res, absDiff, deRank), file = paste0(fname,'_tmp_pval.csv'), sep=',', row.names = F, quote=F)
  
  if(type(readCounts) != 'logical')
  {
    #True DE list from all available real samples
    de = DGEList(counts=readCounts, genes=rownames(readCounts), 
                 group=sampleInfo$condition)
    # TMM Normalization
    de = calcNormFactors(de)
    design <- model.matrix(~sampleInfo$condition)
    de = estimateDisp(de, design)
    fit <- glmQLFit(de, design)
    res <- glmQLFTest(fit,coef=2) 
    res = topTags(res, sort.by="none", n=length(de$genes$genes))
    res = data.frame(res)
    #erGenes = res[res$FDR<.05 & abs(res$logFC) > 1,'genes']
    #Calc abs difference in average b/w conditions
    allCounts = de$counts
    absDiff = abs(rowSums(allCounts[,sampleInfo$condition==0])/sum(sampleInfo$condition==0) - 
                  rowSums(allCounts[,sampleInfo$condition==1])/sum(sampleInfo$condition==1))
    write.table(cbind(res, absDiff), file = paste0(fname,'_tmp_true_pval.csv'), sep=',', row.names = F, quote=F)
    
    #Calc DE for just the training samples:
    sampleInfoTrain = sampleInfo[make.names(c(c0Samps,c1Samps)), , drop = F]
    de = DGEList(counts=readCounts[ , make.names(c(c0Samps,c1Samps))], genes=rownames(readCounts), 
                 group=sampleInfoTrain$condition)
    # TMM Normalization
    de = calcNormFactors(de)
    design <- model.matrix(~sampleInfoTrain$condition)
    de = estimateDisp(de, design)
    fit <- glmQLFit(de, design)
    res <- glmQLFTest(fit,coef=2) 
    res = topTags(res, sort.by="none", n=length(de$genes$genes))
    res = data.frame(res)
    #erGenes = res[res$FDR<.05 & abs(res$logFC) > 1,'genes']
    #Calc abs difference in average b/w conditions
    allCounts = de$counts
    absDiff = abs(rowSums(allCounts[,sampleInfoTrain$condition==0])/sum(sampleInfoTrain$condition==0) - 
                    rowSums(allCounts[,sampleInfoTrain$condition==1])/sum(sampleInfoTrain$condition==1))
    write.table(cbind(res, absDiff), file = paste0(fname,'_training_pval.csv'), sep=',', row.names = F, quote=F)

    #Calc DE for the training samples duplicated to size of synthetic cohorts:
    matExprAlt = readCounts[ , make.names(c(c0Samps,c1Samps))]
    matExprBig = matExprAlt
    synthSize = length(group)
    numTimesDup = floor(synthSize/length(c(c0Samps,c1Samps))) - 1
    for (i in seq(numTimesDup))
    {
      tmpMat = matExprAlt
      newNames = c()
      for (j in colnames(matExprAlt))
      {
        newNames = c(newNames, paste0(j,".",i))
      }
      colnames(tmpMat) = newNames
      matExprBig = cbind(matExprBig, tmpMat)
    }
    bigGroups = c(rep('0',length(c0Samps)), rep('1',length(c1Samps)))
    bigGroups = rep(bigGroups, numTimesDup+1)
    bigSampInfo = data.frame(condition = bigGroups, row.names = colnames(matExprBig))
    de = DGEList(counts=matExprBig, genes=rownames(matExprBig), 
                 group=bigSampInfo$condition)
    # TMM Normalization
    de = calcNormFactors(de)
    design <- model.matrix(~bigSampInfo$condition)
    de = estimateDisp(de, design)
    fit <- glmQLFit(de, design)
    res <- glmQLFTest(fit,coef=2) 
    res = topTags(res, sort.by="none", n=length(de$genes$genes))
    res = data.frame(res)
    #erGenes = res[res$FDR<.05 & abs(res$logFC) > 1,'genes']
    #Calc abs difference in average b/w conditions
    allCounts = de$counts
    absDiff = abs(rowSums(allCounts[,bigSampInfo$condition==0])/sum(bigSampInfo$condition==0) - 
                    rowSums(allCounts[,bigSampInfo$condition==1])/sum(bigSampInfo$condition==1))
    write.table(cbind(res, absDiff), file = paste0(fname,'_trainingDuplicated_pval.csv'), sep=',', row.names = F, quote=F)
  }
}

