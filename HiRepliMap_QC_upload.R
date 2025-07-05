library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(pheatmap)
library(grid)

# function--------------

check_data=function(rep1,rep2,out,rep1_name,rep2_name,sample){
  
  features_colors=setNames(c("grey","#5b859e","#ab84a5","#732f30","#1e5a46" ,"#edb144","black","#747264","#D1BB9E"),c("ND", "SBs", "L_TTR","IZs" , "R_TTR", "TZs","CTRs","Small_term_site","NA"))
  
  r1=read.csv(rep1) %>% dplyr::select(c("repliseq_chrom", "repliseq_start","repliseq_end","manually_annot2" ))
  r2=read.csv(rep2) %>% dplyr::select(c("repliseq_chrom", "repliseq_start","repliseq_end","manually_annot2" ))
  colnames(r1)[4]=rep1_name
  colnames(r2)[4]=rep2_name
  
  df=merge(r1,r2) %>% pivot_longer(cols=c(rep1_name,rep2_name),names_to = "source",values_to = "features")
  df$features[is.na(df$features)]="NA"
  
  sum_df=df %>% dplyr::group_by(source,features)%>%dplyr::summarise(counts=n())
  sum_df=sum_df %>% dplyr::group_by(source)%>% add_count(source,wt=counts,name="total_bins") %>%
    mutate(percentage=counts/total_bins)
  
  
  sum_df$features=factor(sum_df$features,levels = c("IZs","L_TTR","R_TTR","SBs","CTRs","TZs","Small_term_site","ND","NA"))

  
  ggplot(sum_df,aes(x=source,y=percentage,fill=features))+
    geom_bar(position="stack",stat="identity")+
    geom_text(aes(label=paste0(round(percentage,3))),position=position_stack(vjust=0.3),size=4,color="white")+
    scale_fill_manual(values=features_colors)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90))
  
  ggsave(paste0(out,sample,"_Percentage_compare.pdf"),width = 5.5,height = 7,units = "in")
  
  t=merge(r1,r2)
  t[[rep1_name]]=factor(t[[rep1_name]],levels = c("IZs","L_TTR","R_TTR","SBs","CTRs","TZs","ND"))
  t[[rep2_name]]=factor(t[[rep2_name]],levels = c("IZs","L_TTR","R_TTR","SBs","CTRs","TZs","ND"))
  m=as.matrix(table(t[[rep1_name]],t[[rep2_name]]))
  mm=apply(m, 1, function(x) x/sum(x))
  
  
  pdf(paste0(out,sample,"_Percentage_shift.pdf"),width = 8,height = 6)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(mm,display_numbers = TRUE,cluster_cols = F,cluster_rows = F,color = colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7")))(10))
           #,filename =paste0("/home/l538g/workingf/brainbreaks/DSB/Repliseq_smooth/result_v2_DMSO_repeat2/","duplicate_features_map_scale.pdf"),width = 8,height = 6)
  setHook("grid.newpage", NULL, "replace")
  grid.text(rep2_name, y=-0.03, gp=gpar(fontsize=14))
  grid.text(rep1_name, x=-0.03,rot=90, gp=gpar(fontsize=14))
  dev.off()
  
}

# commands --------------------

rep1="/omics/odcf/analysis/OE0574_projects/rdc_xpf_humanhelau2os/repli-seq/u2os-dmso-rep1-hg38/HiRepliMap/RepliFeatures.csv"
rep2="/omics/odcf/analysis/OE0574_projects/rdc_xpf_humanhelau2os/repli-seq/u2os-dmso-rep2-hg38/HiRepliMap/RepliFeatures.csv"

out="/home/l538g/workingf/brainbreaks/DSB/Repliseq_smooth/QC/"
sample="U2os_dmso"

rep1_name="Repeat_1"
rep2_name="Repeat_2"

check_data(rep1=rep1,
           rep2=rep2,
           rep1_name,
           rep2_name,
           out=out,
           sample=sample)

