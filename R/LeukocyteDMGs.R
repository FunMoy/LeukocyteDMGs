LeukocyteDMGs <- function(platform,disease,normal,degene,freq,single_value,DMGs_value,FDR=TRUE){
  if (platform == 450) {
    refdata = refdata.450k
  }
  if (platform == 850) {
    refdata = refdata.850k
  }
  gene=intersect(refdata[,1],normal[,1])
  ref.index=match(gene,refdata[,1])
  test.index=match(gene,normal[,1])
  normal=normal[test.index,]
  disease=disease[test.index,]
  refexp=refdata[ref.index,-1]
  refexp=matrix(unlist(refexp),ncol = dim(refexp)[2],byrow = F)
  disease=matrix(unlist(disease[,-1]),ncol = dim(disease[,-1])[2],byrow = F)
  normal=matrix(unlist(normal[,-1]),ncol = dim(normal[,-1])[2],byrow = F)
  ingene=match(degene,gene)
  ingene=na.omit(ingene)
  control_exp<-normal;
  case_exp<-disease;
  outlier_dir<-NULL;
  outlier_pvalue<-NULL;
  Lgene<-length(gene);
  for (k in 1:Lgene){
    colN<-dim(control_exp)[2]
    colC<-dim(case_exp)[2]
    colR<-dim(refexp)[2]
    ref.tmp=matrix(rep(refexp[k,],Lgene),ncol=colR,byrow=T)-refexp
    ref.up=which((rowSums(ref.tmp>0)/colR)>freq)
    ref.down=which((rowSums(ref.tmp<0)/colR)>freq)
    Nnorm=control_exp[k,]
    N_tmp=matrix(rep(Nnorm,Lgene),ncol=colN,byrow=T)-control_exp
    Nloc.up=which((rowSums(N_tmp>0)/colN)>freq)
    Nloc.down=which((rowSums(N_tmp<0)/colN)>freq)
    Nloc_up=intersect(ref.up,Nloc.up)
    Nloc_down=intersect(ref.down,Nloc.down)
    if (length(intersect(ingene,k))>0 & length(intersect(ingene,Nloc_up))>0) {
      Nloc_up=Nloc_up[-match(intersect(ingene,Nloc_up),Nloc_up)]
    }
    if (length(intersect(ingene,k))>0 & length(intersect(ingene,Nloc_down))>0) {
      Nloc_down=Nloc_down[-match(intersect(ingene,Nloc_down),Nloc_down)]
    }
    reverse=matrix(0,colC,4)
    reverse[,1]=rep(length(Nloc_up),colC)
    reverse[,2]=rep(length(Nloc_down),colC)
    Tcanc=case_exp[k,]
    if (length(Nloc_up)>0){
      N_tmp=matrix(rep(Tcanc,length(Nloc_up)),ncol=colC,byrow=T)-case_exp[Nloc_up,]
      case_p=colSums(N_tmp<0)
      reverse[,3]=case_p
    }
    if (length(Nloc_down)>0){
      N_tmpp=matrix(rep(Tcanc,length(Nloc_down)),ncol=colC,byrow=T)-case_exp[Nloc_down,]
      case_pp=colSums(N_tmpp>0)
      reverse[,4]=case_pp;
    }
    GenePair_sig=NULL
    GenePair=rep(0,colC)
    GenePair[which(reverse[,3]>reverse[,4])]<--1
    GenePair[which(reverse[,3]<reverse[,4])]<-1
    tmp=matrix(c(reverse[,1],reverse[,2],reverse[,1]-reverse[,3]+reverse[,4], reverse[,2]-reverse[,4]+reverse[,3]),ncol=4)
    GenePair_sig<-apply(tmp,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$p.value)
    outlier_dir=rbind(outlier_dir,GenePair)
    outlier_pvalue=rbind(outlier_pvalue,GenePair_sig)
  }
  fdr<-apply(outlier_pvalue,2,function(x) p.adjust(x,method="fdr",length(x)))
  Methout<-list(gene,outlier_dir,fdr);
  names(Methout)<-c('gene','outlier_dir','fdr')
  gene=Methout$gene
  dir=Methout$outlier_dir
  fdr=Methout$fdr
  for(i in 1:dim(fdr)[2])
  {
    index.fdr<-which(fdr[,i]<=single_value);
    dir[-index.fdr,i]<-0;
  }
  index=match(degene,gene)
  index=na.omit(index)
  rm(degene)
  gene=gene[index]
  dir=dir[index,]
  fdr=fdr[index,]
  dir_all<-abs(dir);
  P<-matrix(0,1,dim(fdr)[2]);
  for(j in 1:dim(fdr)[2])
  {
    P[j]<-sum(dir_all[,j])/length(dir[,1]);
  }
  P_all=rowSums(P)/length(dir[1,]);

  p_value<-matrix(0,length(gene),1);
  for(x in 1:length(gene))
  {
    k<-length(which(dir_all[x,]==1));
    n2<-dim(fdr)[2];
    p_value[x]<-1-pbinom(k,n2,P_all);
  }
  if (FDR == TRUE) {
    FDR<-p.adjust(p_value,method="fdr",length(p_value));
    DMGs_index<-which(FDR<=DMGs_value);
    DMGs<-gene[DMGs_index];
  }else{
    DMGs_index<-which(p_value<=DMGs_value);
    DMGs<-gene[DMGs_index];
  }
  out = list(individual=Methout,population=DMGs)
  return(out)
}


