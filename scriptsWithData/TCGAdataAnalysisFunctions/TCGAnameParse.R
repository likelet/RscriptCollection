suppressMessages(library(data.table))
suppressMessages(library(survival))
suppressMessages(library(ggplot2))
suppressMessages(library(CrossICC))
suppressMessages(library(data.table))
suppressMessages(library(survminer))
suppressMessages(library(stringr))
#get tumor or normal according to their barcode
tellTumor_or_Normal = function(barcode) {
  res = c()
  for (i in 1:length(barcode)) {
    str = unlist(strsplit(barcode[i], split = "\\-"))
    if (str[4] %like% "0") {
      res = c(res, "Tumor")
    } else if (str[4] %like% "1") {
      res = c(res, "Normal")
    } else{
      res = c(res, "unkown")
    }
  }
  return(res)
}


#get best suvival cutoff ----
getBestCutOffForSurvival<-function(data.surv,inV,fra=0.1){
  inV.2=unique(inV)
  opt.cut<-median(inV.2)
  avec<-ifelse(inV>=opt.cut,"high","low")
  chisq=survdiff(data.surv~avec)$chisq
  a=c()
  b=c()
  for(i in 1:length(inV.2)){
    cut=inV.2[i]
    tempvect=ifelse(inV>=cut,"high","low")
    min.no=min(table(tempvect))
    if(min.no/length(inV)>fra && min.no/length(inV)<fra ){
      survd<-survdiff(data.surv~tempvect)
      chs<-survd$chisq
      a=c(a,chs)
      b=c(b,cut)
      if(chisq<chs){
        chisq=chs
        opt.cut=cut
      }
    }
  }
  tempvect=ifelse(inV>=opt.cut,"high","low")
  print(survdiff(data.surv~tempvect))
  print(opt.cut)
  return(data.frame(b,a))
}


#plot two bar plot by class
plot_complex_boxplot <- function(plotdata,title=""){
  # library(ggthemes)
  library(ggplot2)
  colnames(plotdata) <- c("Value","Class")
  plotdata$value<-plotdata$Value
  plot.class=unique(plotdata$Class)
  custom.theme <-  theme_bw() +
    # from shenwei.me
    theme(
      
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = "black", size = 0.8),
      axis.line.y = element_line(colour = "black", size = 0.8),
      axis.ticks.x = element_line(size = 0.8),
      axis.ticks.y = element_line(size = 0.8),
      axis.text.x = element_text(
        angle = 30, hjust = 1, vjust = 1
      ),
      legend.position = "none",
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      legend.background = element_rect(fill = "transparent"),
      strip.background = element_rect(
        colour = "white", fill = "white",
        size = 0.2
      ),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      
      text = element_text(
        size = 14,
        #family = "arial",
        face = "bold"
      ),
      plot.title = element_text(
        size = 16,
        #family = "arial",
        face = "bold"
      )
    )
  
  
  
  
  p.value<-wilcox.test(plotdata[plotdata$Class==plot.class[1],1],plotdata[plotdata$Class==plot.class[2],1])$p.value
  print(p.value)
  # library(ggthemes)
  
  p <- ggplot(plotdata)+
    geom_violin(aes(x=factor(Class),fill=factor(Class),y=Value),color=I("grey60"),size=1,width=0.5,alpha=I(0.6)) +
    # geom_jitter(aes(x=factor(isContained),y=TargetSiteRatio),color=I("black"),size=4) +
    geom_jitter(aes(x=factor(Class),y=Value,fill=factor(Class),color=I("grey30")),size=1,shape=21,width=0.2,stroke=0.5,alpha=I(0.6))+
    geom_boxplot(aes(x=factor(Class),y=Value,fill=factor(Class)),size=1,width=0.2,alpha=I(0.6))+
    
    geom_text(x=1.5,y=max(plotdata$Value)*1.05,label=format(p.value,scientific=TRUE,digit=3))
  # geom_segment(x=0.8,xend=2.2,y=max(plotdata$Value)*1.02,yend=max(plotdata$Value)*1.02,size=0.5)
  if(p.value<0.05){
    p=p+geom_text(x=1.5,y=max(plotdata$Value),label="*",size=3)
  }else{
    p=p+geom_text(x=1.5,y=max(plotdata$Value),label="NS",size=2)
  }
  
  p<-p+ylim(0,max(plotdata$Value)*1.18)+
    ggtitle(title)+
    xlab("")+ylab("Estimated RNA abundance by RPKM")+
    scale_fill_manual(values=c("#db625e","#78a1c7"))+custom.theme+scale_y_log10()
  return(p)
  
}


#get clinical information from tcga clinical files ----
#by Zexian Liu
tcga_clinic_matrix <- function(tcga_clinic_path){
  # Read the merged_only_clinical_clin_format Clinical file, in this case i transposed it to keep the clinical feature title as column name
  clinical <- as.data.frame(t(read.delim(tcga_clinic_path,header=T, row.names=1, sep="\t", quote = "")))
  clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
  # get the columns that contain data we can use: days to death, new tumor event, last day contact to....
  
  days_to_new_tumor_event_after_initial_treatment <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))
  # this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
  new_tum <- as.matrix(clinical[, days_to_new_tumor_event_after_initial_treatment])
  new_tum_collapsed <- c()
  for (i in 1 : dim(new_tum)[1]){
    if(sum(is.na(new_tum[i,])) < dim(new_tum)[2]){
      m <- max(new_tum[i,], na.rm=T)
      new_tum_collapsed <- c(new_tum_collapsed, m)
    } else {
      new_tum_collapsed <- c(new_tum_collapsed, "NA")
    }
  }
  
  # do the same to death
  days_to_death <- grep("days_to_death",colnames(clinical))
  death <- as.matrix(clinical[, days_to_death])
  death_collapsed <- c()
  for (i in 1:dim(death)[1]){
    if(sum(is.na(death[i,])) < dim(death)[2]){
      m <- max(death[i,],na.rm=T)
      death_collapsed <- c(death_collapsed,m)
    } else {
      death_collapsed <- c(death_collapsed,"NA")
    }
  }
  
  # and days last follow up here we take the most recent which is the max number
  days_to_last_followup <- grep("days_to_last_followup",colnames(clinical))
  fl <- as.matrix(clinical[,days_to_last_followup])
  fl_collapsed <- c()
  for (i in 1:dim(fl)[1]){
    if(sum(is.na(fl[i,])) < dim(fl)[2]){
      m <- max(fl[i,],na.rm=T)
      fl_collapsed <- c(fl_collapsed,m)
    } else {
      fl_collapsed <- c(fl_collapsed,"NA")
    }
  }
  
  # and put everything together
  all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
  colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")
  # now, to do survival analysis we need three main things:
  # 1- time: this is the time till an event happens
  # 2- status: this indicates which patients have to be kept for the analysis
  # 3- event: this tells i.e. which patients have the gene up- or down-regulated or have no changes in expression
  # since we want to do censored analysis we need to have something to censor the data with. for example, if a patient has no death data BUT there is a date to last followup it means that after that day we know nothing about the patient, therefore after that day it cannot be used for calculations/Kaplan Meier plot anymore, therefore we censor it. so now we need to create vectors for both "time to new tumor" and 'time to death" that contain also the data from censored individuals
  # create vector with time to new tumor containing data to censor for new_tumor
  all_clin$new_time <- c()
  for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
    all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i], as.numeric(as.character(all_clin$new_tumor_days))[i])
  }
  
  # create vector time to death containing values to censor for death
  all_clin$new_death <- c()
  for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
    all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
  }
  
  # now we need to create the vector for censoring the data which means telling R which patients are dead or have new tumor. in this case if a patient has a “days_to_death” it will be assigned 1, and used in the corresponding analysis. the reason why we censor with death events even for recurrence is pretty important. a colleague made me notice that this is a competitive risk problem, where, although a patient can recur and then die, if a patient is dead, it will not recur, therefore is more accurate to censor for death events.
  # create vector for death censoring
  #table(clinical$patient.vital_status)
  # alive  dead
  #   372   161
  
  all_clin$death_event <- ifelse(clinical$patient.vital_status == "alive", 0, 1)
  
  # finally add row.names to clinical
  rownames(all_clin) <- clinical$IDs
  return(all_clin)
}

#get patientName
get_patientID_from_tcga<-function(vec){
  re<-c()
  for (i in 1:length(vec)) {
    tempstr<-vec[i]
    re<-c(re,paste(unlist(stringr::str_split(tempstr,pattern = "-"))[1:3],collapse ="-"))
    
  }
  return(re)
}
#survival function
cust_surp<-function(fit,df=NULL,colorlist=NULL,ti=""){
  if(is.null(colorlist)){
    ggsurvplot(fit, data = df,
               surv.median.line = "hv", # Add medians survival
               
               # Change legends: title & labels
               legend.title = "Risk",
               # Add p-value and confidence intervals
               pval = TRUE,
               title=ti,
               conf.int = FALSE,
               # Add risk table
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               ggtheme = theme_bw() # Change ggplot2 theme
    )
  }else{
    ggsurvplot(fit, data = df,
               surv.median.line = "hv", # Add medians survival
               
               # Change legends: title & labels
               legend.title = "Risk",
               # Add p-value and confidence intervals
               pval = TRUE,
               conf.int = FALSE,
               title=ti,
               # Add risk table
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               palette = colorlist,
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               ggtheme = theme_bw() # Change ggplot2 theme
    )
  }
  
}

