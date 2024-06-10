# Depends
library(data.table)
library(plyr)
library(stringr)

# Tracings with splits
splits <- fread(file='vital_af_traces_split_v2024_05_30.tsv')

# Load infs
ft <- fread(file='vital_af_finetune_predictions_v2024_05_30.tsv')
scratch <- fread(file='vital_af_from_scratch_predictions_v2024_05_30.tsv')
txfr <- fread(file='vital_af_transfer_all_traces_predictions_v2024_05_24.tsv')

# Merge into one
setkey(ft,EMPI,recording_id); setkey(scratch,EMPI,recording_id); setkey(txfr,EMPI,recording_id)
setnames(ft,'af_pred','finetuned_af_pred')
ft[scratch,scratch_af_pred := i.af_pred]
ft[txfr,txfr_af_pred := i.af_pred]

# Analyze
af <- ft
ft[,dummy := 1]

#### Event rate
setDF(ft)
super_cat <- survivor(data=ft,risk_data='dummy',time='time_to_af_ecg',status='incd_af_ecg',eval.t=2)
setDT(ft)

#### ROC and AP
pr_ai <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                marker=af$txfr_af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ft <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                marker=af$finetuned_af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_scratch <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                     marker=af$scratch_af_pred,t0.list=1.999,method='bootstrap',B=500)

################################# ROC curves
ecg_ai <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
                  marker=af$txfr_af_pred,cause=1,times=c(1,1.999))
ft <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
              marker=af$finetuned_af_pred,cause=1,times=c(1,1.999))
scratch <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
                   marker=af$scratch_af_pred,cause=1,times=c(1,1.999))

pdf(file='roc_model_compare.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(ecg_ai,1.999,add=T,col='#1b9e778C',lwd=1.2)
par(new=TRUE)
plot(ft,1.999,add=T,col='#7570b38C',lwd=1.2)
par(new=TRUE)
plot(scratch,1.999,add=T,col='#e7298a8C',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.4,0.2,legend=c('1L ECG-AI (0.722)','Fine-Tuned (0.621)','1L VITAL (0.651)'),
       col=c('#1b9e778C','#7570b38C','#e7298a8C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

### AUPRC
points_ai <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='txfr_af_pred',eval.t=1.999,tolerance=2)
points_ft <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='finetuned_af_pred',eval.t=1.999,tolerance=2)
points_scratch <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='scratch_af_pred',eval.t=1.999,tolerance=2)

pr_no_skill <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                      marker=af$txfr_af_pred,t0.list=c(1,1.999))$ap_summary[4]

# Plot
pdf(file='auprc_model_compare.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))

plot(x=points_ai$sens,y=points_ai$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#1b9e778C',type='l')
par(new=TRUE)
plot(x=points_ft$sens,y=points_ft$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#7570b38C',type='l')
par(new=TRUE)
plot(x=points_scratch$sens,y=points_scratch$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,1),
     col='#e7298a8C',type='l')
par(new=TRUE)

axis(1,at=seq(0,1,0.2),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),las=2,cex.axis=1.6)

mtext("Sensitivity/Recall",1,line=2.8,cex=1.8)
mtext("Precision/PPV",2,line=3.5,cex=1.8)

segments(0,pr_no_skill,1,pr_no_skill,lty=5)

legend(0.3,0.4,legend=c('1L ECG-AI (0.061)','Fine-Tuned (0.058)','1L VITAL (0.081)'),
       col=c('#1b9e778C','#7570b38C','#e7298a8C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)

dev.off()




