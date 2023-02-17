simudata <- function(platform,mye_ratio,lym_ratio){
  if (platform==450) {
    mye_and_lym_methy = mye_and_lym_methy_450k
  }
  if (platform==850) {
    mye_and_lym_methy = mye_and_lym_methy_850k
  }
  myeloid_methy = mye_and_lym_methy$myeloid_methy
  lymphoid_methy = mye_and_lym_methy$lymphoid_methy
  myeloid_mean.methy=as.data.frame(rowMeans(myeloid_methy))
  myeloid_simu.methy=apply(myeloid_mean.methy, 1, function(x)
    runif(50,0.95*x,1.05*x) )
  lymphoid_mean.methy=as.data.frame(rowMeans(lymphoid_methy))
  lymphoid_simu.methy=apply(lymphoid_mean.methy, 1, function(x)
    runif(50,0.95*x,1.05*x) )
  simu.normal=myeloid_simu.methy*0.5+lymphoid_simu.methy*0.5

  myeloid_simu.methy=apply(myeloid_mean.methy, 1, function(x)
    runif(50,0.95*x,1.05*x) )
  select.mye.gene1=sample(1:dim(myeloid_methy)[1], size = 1000)
  select.lym.gene1=sample(1:dim(lymphoid_methy)[1], size = 1000)
  inter.id=intersect(select.lym.gene1,select.mye.gene1)
  select.mye.gene1=select.mye.gene1[-match(inter.id,select.mye.gene1)]
  select.lym.gene1=select.lym.gene1[-match(inter.id,select.lym.gene1)]
  for (i in 1:length(select.mye.gene1)) {
    if (myeloid_mean.methy[select.mye.gene1[i],1]>0.5) {
      myeloid_simu.methy[,select.mye.gene1[i]]=myeloid_simu.methy[,select.mye.gene1[i]]-runif(1,0.3,0.4)
    }else{
      myeloid_simu.methy[,select.mye.gene1[i]]=myeloid_simu.methy[,select.mye.gene1[i]]+runif(1,0.3,0.4)
    }
  }

  lymphoid_simu.methy=apply(lymphoid_mean.methy, 1, function(x)
    runif(50,0.95*x,1.05*x) )
  for (i in 1:length(select.lym.gene1)) {
    if (lymphoid_mean.methy[select.lym.gene1[i],1]>0.5) {
      lymphoid_simu.methy[,select.lym.gene1[i]]=lymphoid_simu.methy[,select.lym.gene1[i]]-runif(1,0.3,0.4)
    }else{
      lymphoid_simu.methy[,select.lym.gene1[i]]=lymphoid_simu.methy[,select.lym.gene1[i]]+runif(1,0.3,0.4)
    }
  }
  DMG = rownames(myeloid_methy)[union(select.mye.gene1,select.lym.gene1)]
  simu.disease=myeloid_simu.methy*mye_ratio+lymphoid_simu.methy*lym_ratio
  simu.data = list(DMG,simu.normal,simu.disease)
  names(simu.data) = c('DMG','simu.normal','simu.disease')
  return(simu.data)
}
