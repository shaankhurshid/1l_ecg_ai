# Depends
library(data.table)
library(plyr)
library(timeROC)
library(survival)
library(APtools)
library(prodlim)
library(Cairo)
library(polspline)
library(nricens)
library(dcurves)
source(file='./functions/functions.R')

# Load inferences
inf <- fread(file='vital_af_transfer_all_traces_predictions_v2024_05_24.tsv')
readings <- fread(file='tracings_adjudicated_013024.csv')
af <- fread(file='Screened.txt')
baseline <- fread(file='baseline.txt')

# Merge readings
setkey(inf,recording_id); setkey(readings,recording_id)
inf[readings,':='(adjudication_af = i.adjudication_af,
                  adjudication_afl = i.adjudication_afl,
                  adjudication_af_afl = i.adjudication_af_afl,
                  adjudication_unreadable = i.adjudication_unreadable)]

# Date formatting
for (j in c('date1_y','end_date')){
  set(inf,j=j,value=as.Date(inf[[j]],format='%m/%d/%Y'))}
for (j in c('fstvisit_date','end_date','af_dx_date',paste0('date',1:19)))
  set(af,j=j,value=as.Date(af[[j]],format='%m/%d/%Y'))

# Remove tracings read as AF
inf <- inf[adjudication_af_afl==0]

# Merge AF date
setkey(af,empi); setkey(inf,empi)
inf[af,':='(af_dx_date = i.af_dx_date)]

# Redefine prevalent and incident AF
inf[,":="(prev_af_ecg = ifelse(c(!is.na(af_dx_date) & (af_dx_date <= date1_y)),1,0))]
inf[,":="(incd_af_ecg = ifelse(c(!is.na(af_dx_date) & (af_dx_date > date1_y)),1,0))]
inf[,":="(time_to_af_ecg = ifelse(incd_af_ecg==1,(af_dx_date-date1_y)/365.25,
                                  pmin((af_dx_date-date1_y)/365.25,((end_date+365.25)-date1_y)/365.25,na.rm=T)))]

# Merge some missing covariates
setkey(baseline,EMPI); setkey(inf,empi)
inf[baseline,":="(mi = i.MI)]

# Narrow to single ECG per person
setkey(inf,empi,date1_y)
af <- inf[,.SD[1],by='empi']

# Remove prevalent AF
af <- af[prev_af_ecg==0]

# Recalculate BMI
af[,bmi := ifelse(is.na(bmi),weight_kg/((height_cm/100)**2),bmi)]

# Address missingness
af[,smoker := ifelse(is.na(smoker),0,smoker)]
af[,sbp := ifelse(is.na(sbp),mean(af$sbp,na.rm=T),sbp)]
af[,dbp := ifelse(is.na(dbp),mean(af$dbp,na.rm=T),dbp)]
af[,height_cm := ifelse(is.na(height_cm),mean(af$height_cm,na.rm=T),height_cm)] 
af[,weight_kg := ifelse(is.na(weight_kg),mean(af$weight_kg,na.rm=T),weight_kg)]

# Recalculate CHARGE-AF with no missingness
af[,charge_complete := (age/5)*0.508 + (white==1)*0.465 + (height_cm/10)*0.248 + (weight_kg/15)*0.115 + 
     (sbp/20)*0.197 + (dbp/10)*(-0.101) + (smoker==1)*0.359 + (bp_med==1)*0.349
   + (dm==1)*0.237 + (hf==1)*0.701 + (mi==1)*0.496]

# Combine ECG-AI and charge
### OG
af[,ecg_ai_logit := log(af_pred/(1-af_pred))]
af[,ecgai_as := 0.038049*age + 0.256315*(gender=='Male') + 0.417466*ecg_ai_logit]
af[,chai_og := 0.40860*charge_complete + 0.38006*ecg_ai_logit]

# Standardized vars
af[,':='(charge_std = (charge_complete - mean(charge_complete))/sd(charge_complete),
         ecg_ai_std = (ecg_ai_logit - mean(ecg_ai_logit))/sd(ecg_ai_logit))]

### ECG-AI RECAL
avg_beta <- mean(af$ecg_ai_logit)
res <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecg_ai_logit, data=af)
km <- survfit(res, data=data.frame(x1=mean(ecg_ai_logit)),type="kaplan-meier")
s0 <- summary(km, times=c(2),extend=TRUE)$surv
af[,ecg_ai_pred2_recal := (1-(s0)^exp(ecg_ai_logit - avg_beta))]

### CHARGE RECAL
avg_beta <- mean(af$charge_complete)
res <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ charge_complete, data=af)
km <- survfit(res, data=data.frame(x1=mean(charge_complete)),type="kaplan-meier")
s0 <- summary(km, times=c(2),extend=TRUE)$surv
af[,charge_pred2_recal := (1-(s0)^exp(charge_complete - avg_beta))]

### OG CHAI RECAL
avg_beta <- mean(af$chai_og)
res <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ chai_og, data=af)
km <- survfit(res, data=data.frame(x1=mean(chai_og)),type="kaplan-meier")
s0 <- summary(km, times=c(2),extend=TRUE)$surv
af[,chai_og_pred2_recal := (1-(s0)^exp(chai_og - avg_beta))]

### CHAI NEW
mod_chai <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ charge_complete + ecg_ai_logit,data=af)
af[,chai_new := predict(mod_chai,newdata=af,type='lp')]

avg_beta <- mean(af$chai_new)
res <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ chai_new, data=af)
km <- survfit(res, data=data.frame(x1=mean(chai_new)),type="kaplan-meier")
s0 <- summary(km, times=c(2),extend=TRUE)$surv
af[,chai_new_pred2 := (1-(s0)^exp(chai_new - avg_beta))]

### ECG-AI AS
avg_beta <- mean(af$ecgai_as)
res <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecgai_as, data=af)
km <- survfit(res, data=data.frame(x1=mean(ecgai_as)),type="kaplan-meier")
s0 <- summary(km, times=c(2),extend=TRUE)$surv
af[,ecgai_as_pred2 := (1-(s0)^exp(ecgai_as - avg_beta))]

# Save out file
write.csv(af,file='full_inference_053124.csv',row.names=F)

################ DISCRIMINATION
## AUROC
c_stat_ai <- timeROC(af$time_to_af_ecg, delta=af$incd_af_ecg,
                     marker=af$af_pred,cause=1,times=1.999)
c_stat_ai_as <- timeROC(af$time_to_af_ecg, delta=af$incd_af_ecg,
                        marker=af$ecgai_as,cause=1,times=1.999)
c_stat_charge <- timeROC(af$time_to_af_ecg, delta=af$incd_af_ecg,
                         marker=af$charge_complete,cause=1,times=1.999)
c_stat_chai <- timeROC(af$time_to_af_ecg, delta=af$incd_af_ecg,
                       marker=af$chai_og,cause=1,times=1.999)

## AUPRC
set.seed(1)
pr_ai <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                marker=af$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                   marker=af$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                    marker=af$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                  marker=af$chai_og,t0.list=1.999,method='bootstrap',B=500)

## HRs
mod_ai <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecg_ai_logit, data=af)
mod_ai_as <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecg_ai_logit + age + gender, data=af)

af[,ecg_ai_lp := predict(mod_ai,newdata=af,type='lp')]
af[,ecg_ai_as_lp := predict(mod_ai_as,newdata=af,type='lp')]

af[,":="(ecg_ai_std = (ecg_ai_logit - mean(ecg_ai_logit))/sd(ecg_ai_logit),
         ecgai_as_std = (ecgai_as - mean(ecgai_as))/sd(ecgai_as),
         chai_new_std = (chai_new - mean(chai_new))/sd(chai_new),
         chai_og_std = (chai_og - mean(chai_og))/sd(chai_og))]

mod_ai_std <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecg_ai_std, data=af)
mod_ai_as_std <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecgai_as_std, data=af)
mod_charge_std <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ charge_std, data=af)
mod_chai_std <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ chai_og_std, data=af)

################################# ROC curves
ai <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
              marker=af$af_pred,cause=1,times=c(1,1.999))
ai_as <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
                 marker=af$ecgai_as,cause=1,times=c(1,1.999))
charge <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
                  marker=af$charge_complete,cause=1,times=c(1,1.999))
chai <- timeROC(T=af$time_to_af_ecg, delta=af$incd_af_ecg,
                marker=af$chai_og,cause=1,times=c(1,1.999))

pdf(file='roc_full_inf.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
plot.new() #?
plot(ai,1.999,add=T,col='#1b9e778C',lwd=1.2)
par(new=TRUE)
plot(ai_as,1.999,add=T,col='#d95f028C',lwd=1.2)
par(new=TRUE)
plot(charge,1.999,add=T,col='#7570b38C',lwd=1.2)
par(new=TRUE)
plot(chai,1.999,add=T,col='#e7298a8C',lwd=1.2)
axis(1,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6)
axis(2,at=seq(0,1,0.2),cex.axis=1.6,las=2)
title(xlab='1 - Specificity',line=2.5,cex.lab=1.8)
title(ylab='Sensitivity',line=3.2,cex.lab=1.8)
legend(0.3,0.2,legend=c('ECG-AI (0.666)','CHARGE-AF (0.679)','ECG-AI AS (0.695)','ECG-AI + CHARGE-AF (0.703)'),
       col=c('#1b9e778C','#7570b38C','#d95f028C','#e7298a8C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)
dev.off()

### AUPRC
points_ai <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='ecg_ai_logit',eval.t=1.999,tolerance=2)
points_ai_as <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='ecgai_as',eval.t=1.999,tolerance=2)
points_charge <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='charge_complete',eval.t=1.999,tolerance=2)
points_chai <- auprc(data=af,time='time_to_af_ecg',status='incd_af_ecg',marker='chai_og',eval.t=1.999,tolerance=2)

pr_no_skill <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                      marker=af$af_pred,t0.list=c(1,1.999))$ap_summary[4]

# Plot
pdf(file='auprc_full_inf.pdf',height=4,width=4,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))

plot(x=points_ai$sens,y=points_ai$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,0.4),
     col='#1b9e778C',type='l')
par(new=TRUE)
plot(x=points_ai_as$sens,y=points_ai_as$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,0.4),
     col='#d95f028C',type='l')
par(new=TRUE)
plot(x=points_charge$sens,y=points_charge$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,0.4),
     col='#7570b38C',type='l')
par(new=TRUE)
plot(x=points_chai$sens,y=points_chai$ppv,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xlim=c(0,1),ylim=c(0,0.4),
     col='#e7298a8C',type='l')
par(new=TRUE)

axis(1,at=seq(0,1,0.2),cex.axis=1.6)
axis(2,at=seq(0,0.4,0.1),las=2,cex.axis=1.6)

mtext("Sensitivity/Recall",1,line=2.8,cex=1.8)
mtext("Precision/PPV",2,line=3.5,cex=1.8)

segments(0,pr_no_skill,1,pr_no_skill,lty=5)

legend(0.3,0.4,legend=c('ECG-AI (0.053)','CHARGE-AF (0.062)','ECG-AI AS (0.060)','ECG-AI + CHARGE-AF (0.063)'),
       col=c('#1b9e778C','#7570b38C','#d95f028C','#e7298a8C'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)

dev.off()

# Metrics
### HR
mod_charge <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ charge_std,data=af)
charge_hr <- c(exp(mod_charge$coefficients[1]),exp(confint(mod_charge)))
mod_ecgai <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ ecg_ai_std,data=af)
ecgai_hr <- c(exp(mod_ecgai$coefficients[1]),exp(confint(mod_ecgai)))
mod_both <- coxph(Surv(time_to_af_ecg,incd_af_ecg) ~ charge_std + ecg_ai_std,data=af)
charge_both_hr <- c(exp(mod_both$coefficients[1]),exp(confint(mod_both))[1,])
ecgai_both_hr <- c(exp(mod_both$coefficients[2]),exp(confint(mod_both))[2,])

# HR PLOT
pdf('hr_plot.pdf',height=5,width=7,
    pointsize=3)
par(oma=c(1,1,0,1))
par(mar=c(4,17,5,1))
par(mgp=c(3,2,0))
plot(x=c(charge_hr[1],ecgai_hr[1],charge_both_hr[1],ecgai_both_hr[1]),
     y=c(4,3,1.5,0.5),xlim=c(0.95,2),ylim=c(0,4),log='x',
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=c('#7570b38C','#1b9e778C',
                                                    '#7570b3','#1b9e77'),cex=4,bty='n')
axis(1,cex.axis=2.5,at=c(1,1.25,1.5,1.75,2),pos=0.2)
axis(2,at=c(4,3,1.5,0.5),labels=FALSE,cex=2.5,pos=0.95)
mtext('Hazard ratio for incident AF',side=1,line=2.7,cex=2.5)
segments(c(charge_hr[2],ecgai_hr[2],charge_both_hr[2],ecgai_both_hr[2]),c(4,3,1.5,0.5),
         c(charge_hr[3],ecgai_hr[3],charge_both_hr[3],ecgai_both_hr[3]),c(4,3,1.5,0.5),col=c('#7570b38C','#1b9e778C',
                                                                                             '#7570b3','#1b9e77'),
         lwd=2.2)
segments(1,0.5,1,4,col='black',lty=5)

plot_names <- c('CHARGE-AF','ECG-AI','CHARGE-AF','ECG-AI')
n <- 1
for (i in c(4,3,1.5,0.5)){
  mtext(paste0(plot_names[n]),
        side=2,line=-0.5,las=2,cex=2.5,at=i)
  n <- n+1
}
dev.off()

################ CALIBRATION
###################################################################### FAMILY 1: MGH
### Fit adaptive hazard model for ECG-AI
af$cox.2yr.cll_ecgai_raw <- log(-log(1-af$af_pred))
calibrate.cox_ecgai_raw <- hare(data=af$time_to_af_ecg,delta=af$incd_af_ecg,
                                cov=as.matrix(af$cox.2yr.cll_ecgai_raw))
predict.grid.cox_ecgai_raw <- seq(quantile(af$af_pred,probs=0.01),
                                  quantile(af$af_pred,probs=0.99),length=100)
predict.grid.cox.cll_ecgai_raw <- log(-log(1-predict.grid.cox_ecgai_raw))
predict.calibrate.cox_ecgai_raw <- phare(2,predict.grid.cox.cll_ecgai_raw,calibrate.cox_ecgai_raw)
predict.calibrate.cox_ecgai_raw2 <- phare(2,af$cox.2yr.cll_ecgai_raw,calibrate.cox_ecgai_raw)
ici_ecgai_raw <- mean(abs(af$af_pred - predict.calibrate.cox_ecgai_raw2))

### Fit adaptive hazard model for ECG-AI recalibrated
af$cox.2yr.cll_ecgai <- log(-log(1-af$ecg_ai_pred2_recal))
calibrate.cox_ecgai <- hare(data=af$time_to_af_ecg,delta=af$incd_af_ecg,
                            cov=as.matrix(af$cox.2yr.cll_ecgai))
predict.grid.cox_ecgai <- seq(quantile(af$ecg_ai_pred2_recal,probs=0.01),
                              quantile(af$ecg_ai_pred2_recal,probs=0.99),length=100)
predict.grid.cox.cll_ecgai <- log(-log(1-predict.grid.cox_ecgai))
predict.calibrate.cox_ecgai <- phare(2,predict.grid.cox.cll_ecgai,calibrate.cox_ecgai)
predict.calibrate.cox_ecgai2 <- phare(2,af$cox.2yr.cll_ecgai,calibrate.cox_ecgai)
ici_ecgai <- mean(abs(af$ecg_ai_pred2_recal - predict.calibrate.cox_ecgai2))

### Fit adaptive hazard model for ECG-AI + AS
af$cox.2yr.cll_ecgai_as <- log(-log(1-af$ecgai_as_pred2))
calibrate.cox_ecgai_as <- hare(data=af$time_to_af_ecg,delta=af$incd_af_ecg,
                               cov=as.matrix(af$cox.2yr.cll_ecgai_as))
predict.grid.cox_ecgai_as <- seq(quantile(af$ecgai_as_pred2,probs=0.01),
                                 quantile(af$ecgai_as_pred2,probs=0.99),length=100)
predict.grid.cox.cll_ecgai_as <- log(-log(1-predict.grid.cox_ecgai_as))
predict.calibrate.cox_ecgai_as <- phare(2,predict.grid.cox.cll_ecgai_as,calibrate.cox_ecgai_as)
predict.calibrate.cox_ecgai_as2 <- phare(2,af$cox.2yr.cll_ecgai_as,calibrate.cox_ecgai_as)
ici_ecgai_as <- mean(abs(af$ecgai_as_pred2 - predict.calibrate.cox_ecgai_as2))

### Fit adaptive hazard model for CHARGE-AF
af$cox.2yr.cll_charge <- log(-log(1-af$charge_pred2_recal))
calibrate.cox_charge <- hare(data=af$time_to_af_ecg,delta=af$incd_af_ecg,
                             cov=as.matrix(af$cox.2yr.cll_charge))
predict.grid.cox_charge <- seq(quantile(af$charge_pred2_recal,probs=0.01),
                               quantile(af$charge_pred2_recal,probs=0.99),length=100)
predict.grid.cox.cll_charge <- log(-log(1-predict.grid.cox_charge))
predict.calibrate.cox_charge <- phare(2,predict.grid.cox.cll_charge,calibrate.cox_charge)
predict.calibrate.cox_charge2 <- phare(2,af$cox.2yr.cll_charge,calibrate.cox_charge)
ici_charge <- mean(abs(af$charge_pred2_recal - predict.calibrate.cox_charge2))

### Fit adaptive hazard model for CHAI OG
af$cox.2yr.cll_chai <- log(-log(1-af$chai_og_pred2_recal))
calibrate.cox_chai <- hare(data=af$time_to_af_ecg,delta=af$incd_af_ecg,
                           cov=as.matrix(af$cox.2yr.cll_chai))
predict.grid.cox_chai <- seq(quantile(af$chai_og_pred2_recal,probs=0.01),
                             quantile(af$chai_og_pred2_recal,probs=0.99),length=100)
predict.grid.cox.cll_chai <- log(-log(1-predict.grid.cox_chai))
predict.calibrate.cox_chai <- phare(2,predict.grid.cox.cll_chai,calibrate.cox_chai)
predict.calibrate.cox_chai2 <- phare(2,af$cox.2yr.cll_chai,calibrate.cox_chai)
ici_chai <- mean(abs(af$chai_og_pred2_recal - predict.calibrate.cox_chai2))

# Plots for visualization
pdf(file='cal_combined.pdf',height=4,width=4,
    pointsize=3)
par(mar=c(4,6,1,1))
par(oma=c(1,1,1,1))
x <- -10; y <- -10

### ECG AI
plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#1b9e778C',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')
par(new=TRUE)
plot(predict.grid.cox_ecgai*100,predict.calibrate.cox_ecgai*100,type="l",lwd=2,lty=1,col="#1b9e778C",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
par(new=TRUE)

### ECG AI AS
plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#d95f028C',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')
par(new=TRUE)
plot(predict.grid.cox_ecgai_as*100,predict.calibrate.cox_ecgai_as*100,type="l",lwd=2,lty=1,col="#d95f028C",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
par(new=TRUE)

### ECG-AI
plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#7570b38C',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')
par(new=TRUE)
plot(predict.grid.cox_charge*100,predict.calibrate.cox_charge*100,type="l",lty=1,lwd=2,col="#7570b38C",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
par(new=TRUE)

### CH-AI
plot(x,y,xlab='',ylab='',yaxt='n',bty='n',col='#e7298a8C',
     xaxt='n',xlim=c(0,30),ylim=c(0,30),pch=19,cex=1.5,bty='n')
par(new=TRUE)
plot(predict.grid.cox_chai*100,predict.calibrate.cox_chai*100,type="l",lty=1,lwd=2,col="#e7298a8C",
     xlim=c(0,30),ylim=c(0,30),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')

axis(1,at=seq(0,30,10),cex.axis=1.7)
axis(2,cex.axis=1.7,at=seq(0,30,10),las=1)

segments(-1,-1,31,31,lwd=1,lty=5)

mtext("Predicted risk of AF at 2 years (%)",side=1,cex=1.8,line=3)
mtext("Incidence of AF at 2 years (%)",side=2,cex=1.8,line=4.5)

legend(1,30,legend=c('ECG-AI','CHARGE-AF','ECG-AI AS','ECG-AI + CHARGE-AF'),col=c('#1b9e778C','#7570b38C','#d95f028C','#e7298a8C'),
       lty=1,lwd=2,bty='n',cex=1.5)
dev.off()

################ KMs
af[,":="(ecg_as_cat = factor(ifelse(ecgai_as_pred2 < 0.02,'< 2%',
                                    ifelse(ecgai_as_pred2 >= 0.03,'\u2265 3%','2-3%')),levels=c('\u2265 3%','2-3%','< 2%')))]
af[,":="(ecg_as_cat_super = factor(ifelse(ecgai_as_pred2 < quantile(ecgai_as_pred2,0.05),'low',
                                          ifelse(ecgai_as_pred2 >= quantile(ecgai_as_pred2,0.95),'high','medium')),
                                   levels=c('high','medium','low')))]
af[,":="(ecgai_as_high = ifelse(ecg_as_cat_super=='high',1,0))]
af[,":="(ecgai_as_3 = ifelse(ecg_as_cat=='\u2265 3%',1,0))]
af[,":="(ecgai_as_low = ifelse(ecg_as_cat_super=='low',1,0))]
af[,":="(charge_cat = factor(ifelse(charge_pred2_recal < 0.02,'< 2%',
                                    ifelse(charge_pred2_recal >= 0.03,'\u2265 3%','2-3%')),levels=c('\u2265 3%','2-3%','< 2%')))]
af[,":="(combo_cat = factor(ifelse(ecg_as_cat=='\u2265 3%' & charge_cat=='\u2265 3%','ECG-AI AS and CHARGE-AF',
                                   ifelse(ecg_as_cat!='\u2265 3%' & charge_cat!='\u2265 3%','None',
                                          ifelse(ecg_as_cat=='\u2265 3%','ECG-AI AS Only','CHARGE-AF Only'))),
                            levels=c('ECG-AI AS and CHARGE-AF','CHARGE-AF Only','ECG-AI AS Only','None')))]
######### KM
########## ECG AI AS 
prod_af <- prodlim(Hist(time_to_af_ecg,incd_af_ecg)~ecg_as_cat,data=af)
# Plot
CairoPDF(file='km_ecgai_as.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prod_af,"cuminc",ylim=c(0,0.06),xlim=c(0,2), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.06,0.01),axis2.las=2,axis2.cex.axis=2.5, #y-axis labeling parameters
     axis1.at=seq(0,2,0.5),axis1.labels=as.character(seq(0,2,0.5)),axis1.padj=0.5,axis1.cex.axis=2.5, #x-axis labeling parameters
     col=c("#e31a1c",'#fd8d3c','#fed976'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     atrisk.title='                  ',atrisk.pos=0,atrisk.line=c(7,9,11), # position of the N at risk rows
     atrisk.cex=1.8,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0.25,0.75,1.25,1.75), # x-axis points where you want the N at risk to show
     xlab='',ylab='', # Empty axis labels, using mtext instead
     legend.x=0,legend.y=0.065,legend.title='',
     legend.cex=2.2,legend.legend=c("\u2265 3%","2-3%","< 2%"))
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.03,cex=2.5) # y-axis label
mtext("Years",side=1, line=-4.5,cex=2.5) # x-axis label
mtext('Stratum',side=1, line=-4.5,cex=1.8,at=-0.15) # descriptor for N at risk
mtext('ECG-AI AS Risk',side=3, line=-1.8,cex=2,at=0.4) # descriptor for N at risk
mtext('\u2265 3%',1,line=-2.6,cex=1.8,at=-0.1)
mtext('2-3%',1,line=-0.6,cex=1.8,at=-0.1)
mtext('< 2%',1,line=1.4,cex=1.8,at=-0.1)
dev.off()

########## SUPER
prod_af_super <- prodlim(Hist(time_to_af_ecg,incd_af_ecg)~ecg_as_cat_super,data=af)
# Plot
CairoPDF(file='km_ecgai_as_super.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prod_af_super,"cuminc",ylim=c(0,0.10),xlim=c(0,2), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.10,0.02),axis2.las=2,axis2.cex.axis=2.5, #y-axis labeling parameters
     axis1.at=seq(0,2,0.5),axis1.labels=as.character(seq(0,2,0.5)),axis1.padj=0.5,axis1.cex.axis=2.5, #x-axis labeling parameters
     col=c("#e31a1c",'#fd8d3c','#fed976'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     atrisk.title='                  ',atrisk.pos=0,atrisk.line=c(7,9,11), # position of the N at risk rows
     atrisk.cex=1.8,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0.25,0.75,1.25,1.75), # x-axis points where you want the N at risk to show
     xlab='',ylab='', # Empty axis labels, using mtext instead
     legend.x=0,legend.y=0.105,legend.title='',
     legend.cex=2.2,legend.legend=c("Top 5%","Middle 90%","Bottom 5%"))
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.05,cex=2.5) # y-axis label
mtext("Years",side=1, line=-4.5,cex=2.5) # x-axis label
mtext('Stratum',side=1, line=-4.5,cex=1.8,at=-0.15) # descriptor for N at risk
mtext('ECG-AI AS Risk',side=3, line=-1.8,cex=2,at=0.4) # descriptor for N at risk
mtext('Top 5%',1,line=-2.6,cex=1.8,at=-0.1)
mtext('Middle 90%',1,line=-0.6,cex=1.8,at=-0.2)
mtext('Bottom 5%',1,line=1.4,cex=1.8,at=-0.2)
dev.off()

########## BOTH STRATA
prod_af_strat <- prodlim(Hist(time_to_af_ecg,incd_af_ecg)~combo_cat,data=af)
# Plot
CairoPDF(file='km_combo.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prod_af_strat,"cuminc",ylim=c(0,0.06),xlim=c(0,2), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.06,0.01),axis2.las=2,axis2.cex.axis=2.5, #y-axis labeling parameters
     axis1.at=seq(0,2,0.5),axis1.labels=as.character(seq(0,2,0.5)),axis1.padj=0.5,axis1.cex.axis=2.5, #x-axis labeling parameters
     col=c("#b10026",'#fc4e2a','#fd8d3c','#fed976'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     atrisk.title='                  ',atrisk.pos=0,atrisk.line=c(7,9,11,13), # position of the N at risk rows
     atrisk.cex=1.8,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0.25,0.75,1.25,1.75), # x-axis points where you want the N at risk to show
     xlab='',ylab='', # Empty axis labels, using mtext instead
     legend.x=0,legend.y=0.075,legend.title='',
     legend.cex=2,legend.legend=c('ECG-AI AS and CHARGE-AF','CHARGE-AF Only','ECG-AI AS Only','None'))
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.03,cex=2.5) # y-axis label
mtext("Years",side=1, line=-6,cex=2.5) # x-axis label
mtext('Both',1,line=-4,cex=1.8,at=-0.04)
mtext('CHARGE-AF Only',1,line=-2,cex=1.8,at=-0.33)
mtext('ECG-AI AS Only',1,line=0,cex=1.8,at=-0.29)
mtext('None',1,line=2,cex=1.8,at=-0.08)
dev.off()

################ INCIDENCE
setDF(af)
top5 <- survivor(data=af,risk_data='ecgai_as_high',time='time_to_af_ecg',
                 status='incd_af_ecg',eval.t=2)
bottom5 <- survivor(data=af,risk_data='ecgai_as_low',time='time_to_af_ecg',
                    status='incd_af_ecg',eval.t=2)
combo_strat <- survivor(data=af,risk_data='combo_cat',time='time_to_af_ecg',
                        status='incd_af_ecg',eval.t=2)
setDT(af)

################ CORR
cor.test(af$ecg_ai_logit,af$charge)

################ NRI
af[,dummy := 1]
ecgai_v_65 <- nricens(time=af$time_to_af_ecg,event=af$incd_af_ecg,
                      p.std=af$dummy,p.new=af$ecgai_as_3,updown='category',
                      cut=0.5,niter=500,t0=2)

################ DECISION CURVE
dcurve <- dca(Surv(time_to_af_ecg,incd_af_ecg) ~ ecgai_as_pred2,
              data=af,thresholds = seq(0,0.25,0.01),time=1.99)
standardized_net_benefit(dcurve)
ggsave('dcurve_std.pdf',
       height=2,width=3,units='in',scale=2)

net_intervention_avoided(dcurve)
ggsave('net_avoid.pdf',
       height=2,width=3,units='in',scale=2)

dcurve_data <- as.data.table(as_tibble(dcurve))
dca_ecgai_as <- dcurve_data[label=='ecgai_as_pred2']
screen_all <- dcurve_data[label=='Treat All']
screen_none <- dcurve_data[label=='Treat None']

################ SECONDARIES
## ALT TIME POINTS
pr_ai <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                marker=af$af_pred,t0.list=c(0.5,1),method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                   marker=af$ecgai_as,t0.list=c(0.5,1),method='bootstrap',B=500)
pr_charge <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                    marker=af$charge_complete,t0.list=c(0.5,1),method='bootstrap',B=500)
pr_chai <- APSurv(stime=af$time_to_af_ecg,status=af$incd_af_ecg,
                  marker=af$chai_og,t0.list=c(0.5,1),method='bootstrap',B=500)

## AGE SUBGROUPS
af[,age_group := ifelse(age < 70,'65-69',
                        ifelse(age >= 80,'80+','70-79'))]
old <- af[age_group=='80+']
med <- af[age_group=='70-79']
young <- af[age_group=='65-69']

### YOUNG
pr_ai <- APSurv(stime=young$time_to_af_ecg,status=young$incd_af_ecg,
                marker=young$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=young$time_to_af_ecg,status=young$incd_af_ecg,
                   marker=young$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=young$time_to_af_ecg,status=young$incd_af_ecg,
                    marker=young$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=young$time_to_af_ecg,status=young$incd_af_ecg,
                  marker=young$chai_og,t0.list=1.999,method='bootstrap',B=500)

### MED
pr_ai <- APSurv(stime=med$time_to_af_ecg,status=med$incd_af_ecg,
                marker=med$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=med$time_to_af_ecg,status=med$incd_af_ecg,
                   marker=med$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=med$time_to_af_ecg,status=med$incd_af_ecg,
                    marker=med$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=med$time_to_af_ecg,status=med$incd_af_ecg,
                  marker=med$chai_og,t0.list=1.999,method='bootstrap',B=500)

### OLD
pr_ai <- APSurv(stime=old$time_to_af_ecg,status=old$incd_af_ecg,
                marker=old$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=old$time_to_af_ecg,status=old$incd_af_ecg,
                   marker=old$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=old$time_to_af_ecg,status=old$incd_af_ecg,
                    marker=old$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=old$time_to_af_ecg,status=old$incd_af_ecg,
                  marker=old$chai_og,t0.list=1.999,method='bootstrap',B=500)

### MEN
men <- af[gender=='Male']
pr_ai <- APSurv(stime=men$time_to_af_ecg,status=men$incd_af_ecg,
                marker=men$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=men$time_to_af_ecg,status=men$incd_af_ecg,
                   marker=men$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=men$time_to_af_ecg,status=men$incd_af_ecg,
                    marker=men$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=men$time_to_af_ecg,status=men$incd_af_ecg,
                  marker=men$chai_og,t0.list=1.999,method='bootstrap',B=500)

### WOMEN
women <- af[gender=='Female']
pr_ai <- APSurv(stime=women$time_to_af_ecg,status=women$incd_af_ecg,
                marker=women$af_pred,t0.list=1.999,method='bootstrap',B=500)
pr_ai_as <- APSurv(stime=women$time_to_af_ecg,status=women$incd_af_ecg,
                   marker=women$ecgai_as,t0.list=1.999,method='bootstrap',B=500)
pr_charge <- APSurv(stime=women$time_to_af_ecg,status=women$incd_af_ecg,
                    marker=women$charge_complete,t0.list=1.999,method='bootstrap',B=500)
pr_chai <- APSurv(stime=women$time_to_af_ecg,status=women$incd_af_ecg,
                  marker=women$chai_og,t0.list=1.999,method='bootstrap',B=500)

## Save out
write.csv(af,file='full_inference_processed_053124.csv')



