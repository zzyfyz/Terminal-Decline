setwd(dirg <- "C:/Users/feiyi/OneDrive/Desktop/Katie/Terminal-Decline/Result")

###########################################################################
# Read csv files
text <- list.files(pattern="mod.result.")
num <- as.numeric(unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]])))

data_frames <- lapply(num, function(i) {
  file_name <- paste0("mod.result.", i, ".csv") 
  read.csv(file_name)
})

I=length(data_frames)

alpha00.mean<-rep(NA,I)
alpha01.mean<-rep(NA,I)
alpha02.mean<-rep(NA,I)
alpha03.mean<-rep(NA,I)
alpha04.mean<-rep(NA,I)
alpha11.mean<-rep(NA,I)
alpha12.mean<-rep(NA,I)
alpha13.mean<-rep(NA,I)
b.mean<-rep(NA,I)
c.mean<-rep(NA,I)
sigm_b.mean<-rep(NA,I)
sigm_u.mean<-rep(NA,I)
sigm_e.mean<-rep(NA,I)
lambda0.mean<-rep(NA,I)
gamma.mean<-rep(NA,I)

for(i in 1:I){ 
  alpha00.mean[i] <- data_frames[[i]][1,1] 
  alpha01.mean[i] <- data_frames[[i]][2,1] 
  alpha02.mean[i] <-data_frames[[i]][3,1] 
  alpha03.mean[i] <-data_frames[[i]][4,1] 
  alpha04.mean[i] <-data_frames[[i]][5,1] 
  alpha11.mean[i] <-data_frames[[i]][6,1] 
  alpha12.mean[i] <-data_frames[[i]][7,1] 
  alpha13.mean[i] <-data_frames[[i]][8,1] 
  b.mean[i] <-data_frames[[i]][9,1]
  c.mean[i] <-data_frames[[i]][10,1]
  
  sigm_b.mean[i] <-data_frames[[i]][11,1] 
  sigm_u.mean[i] <-data_frames[[i]][12,1]
  sigm_e.mean[i] <- data_frames[[i]][13,1]
    
  lambda0.mean[i] <- data_frames[[i]][14,1]
  gamma.mean[i] <- data_frames[[i]][15,1]
}

Sim.results=cbind(alpha00.mean,alpha01.mean,alpha02.mean,alpha03.mean,alpha04.mean,alpha11.mean,alpha12.mean,alpha13.mean,b.mean,c.mean,sigm_b.mean,sigm_u.mean,sigm_e.mean,lambda0.mean,gamma.mean)

round(colMeans(Sim.results),2)


