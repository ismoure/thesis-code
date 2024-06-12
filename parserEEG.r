setwd('../Desktop/data/Assemblies/Sentences')

loadeeg = function(filename) {
  eegfile = paste(filename,".eeg",sep="")
  eegfileSize = file.info(eegfile)$size
  nchannels = 32 # number of channels in eeg recording
  bytesperchan = 2 # number of bytes in eeg measurement per channel
  nsamples = eegfileSize/bytesperchan # calculate number of samples in the data
  channeldata_preshaped = readBin(eegfile, integer(),n=nsamples, size=2)
  channeldata = (matrix(channeldata_preshaped, nchannels,(length(channeldata_preshaped))/nchannels)) # rows are channels, columns are timepoints
  return(channeldata*.1)
}

loadtrg = function(filename) {
  trigs = read.csv(paste0(filename, ".vmrk"),skip=13,header=FALSE)  #read vmrk file
  trigs = cbind(trigs[,3], as.numeric(substring(trigs[,2],first=3)))
  return(trigs)
}

epoch = function(trgdata,eegdata,channels,epochstart,epochend) { # epochstart and epochend times are relative to trigger time
    all_trials = array(0,dim=c(nrow(trgdata),length(channels),epochend-epochstart+1))
    for (s in seq_along(1:nrow(trgdata))) { # for all trials
       trgtime = trgdata[s,1]
       all_trials[s,,] = eegdata[channels,(trgtime+epochstart):(trgtime+epochend)]
    }
               return(all_trials)
}

# load data
beh = read.table('994.dat',skip=15,sep='\t')
trg = loadtrg('parser994')
eeg = loadeeg('parser994')

# filter EEG data with 0.7 high-pass filter
library(signal)
bwh = butter(4, 0.35/1000, type="high")
filteeg = t(apply(eeg,1,function(x) filtfilt(bwh,x)))

# downsample EEG data to 16.12903 Hz
trg[,1] = trg[,1]/62
trg[,1] = round(trg[,1])
filteeg = t(apply(filteeg,1,function(x) decimate(x,62)))

cond = c('Balanced','Unbalanced','Imperative','Blah blah blah blah','Story')
beh[which(beh[,3]=='blah blah blah blah'),2] = 0 # Changes blah condition code to 0
blockstart = which(beh[,5]==0)  # Find block start trial
blockstart = c(blockstart,seq(blockstart[length(blockstart)],blockstart[length(blockstart)]+(11*13),13)) # Add stories
blockstart[which(beh[blockstart,2]==0)] = blockstart[which(beh[blockstart,2]==0)]+1 # In case block starts with blahs; can run twice for Isobel's data but better to add contingency for two or more starting blahs

# FFTs
epochs = epoch(trg[blockstart,],filteeg,1:32,17,240) # 1 to 15 seconds as in Jin et al. (2018)
fft1 = (2*(abs(fft(colMeans(epochs[which(beh[blockstart,2]>0 & beh[blockstart,2]<13),13,]))))/224)^2
fft2 = (2*(abs(fft(colMeans(epochs[which(beh[blockstart,2]>12 & beh[blockstart,2]<25),13,]))))/224)^2
fft3 = (2*(abs(fft(colMeans(epochs[which(beh[blockstart,2]>24 & beh[blockstart,2]<37),13,]))))/224)^2
#fft4 = (2*(abs(fft(colMeans(epochs[which(beh[blockstart,2]==0),13,]))))/160)^2
fft5 = (2*(abs(fft(colMeans(epochs[which(beh[blockstart,2]>36),13,]))))/224)^2
ffft = 16.129/2 * seq(0,1,(1/(224/2))) # fs = 16.129; fs * 10s epoch
for (i in 1:length(cond)) {
   dat = get(paste0('fft',i))[1:64]
   plot(ffft[1:64],dat,type='l',xlab='Frequency (in Hz)',ylab='Power',yaxt="none",ylim=c(0,2),lwd=2,main=cond[i])
   readline()
}
plot(ffft[1:64],fft1,type='l',xlab='Frequency (in Hz)',ylab='Power',yaxt="none",ylim=c(0,2),lwd=2)
lines(ffft[1:64],fft2[1:64],col='green')
lines(ffft[1:64],fft3[1:64],col='blue')
lines(ffft[1:64],fft5[1:64],col='cyan')
legend(x='topright',legend = cond[-4],lty=1,col=c('black','green','blue','cyan'))
arrows(2.75,.85,3.125,.68,lty=2)
png("FFTpower.png",width=20,height=15,units="cm",res=150)
dev.off()

# Wavelets
fpow = array(0,dim=c(30,681,2000))
fang = array(0,dim=c(30,681,2000))

for (freq in seq(1,30,1)) {
   # create wavelet
   window = 2
   t = seq(window/2*-1,window/2,1/1000) #time = -0.25:1/S:0.25; 1000 = sampling rate
   cycles = window*freq/2 # for varying number of cycles depending on frequency; 1000 = sampling rate
   s = (cycles/(2*pi*freq))^2
   wavelet = exp(2*1i*pi*freq*t) * exp(-t^2/(2*s)) # wavelet = sin(2*pi*frequency*t)
   plot(Re(wavelet),type="l")

   # FFT parameters
   n_wavelet            = length(wavelet)
   n_data               = 2000
   n_convolution        = n_wavelet+n_data-1
   half_of_wavelet_size = (length(wavelet)-1)/2

   # FFT of wavelet and EEG data
   fft_wavelet = fft(c(wavelet,rep(0,n_convolution-length(wavelet))))
   for (eachtrial in 1:681) {
      fft_data = fft(c(epochs[eachtrial,13,1:2000],rep(0,n_convolution-2000))) 
      convolution_result_fft = fft(c(fft_wavelet*fft_data,rep(0,n_convolution-length(fft_wavelet))),inverse=TRUE)/length(n_convolution) * sqrt(s) # ifft = fft(x, inverse = TRUE)/length(x)
      convolution_result_fft_cut = convolution_result_fft[(half_of_wavelet_size+1):(length(convolution_result_fft)-half_of_wavelet_size)]		 
      fpow[freq,eachtrial,] = abs(convolution_result_fft_cut)
      # fang[x] = mean(Arg(convolution_result_fft_cut))
   }
}
norm_all=apply(fpow,1,FUN=mean)
colmap = read.csv('rbcolmap.txt',header=FALSE)
par(bg="black",fg="white",col.main="white",col.lab="white",col.axis="white",mar=c(5,5,5,3))
for (i in 1:7) {
   png(paste(cond[i],".png",sep=""),width=20,height=15,units="cm",res=150)
   plot(NA,main=cond[i],xlim=range(seq(300,1700,1)),ylim=range(seq(1,30,1)),xaxt="n",xlab='Time (in ms)',ylab="Frequency (in Hz)",font.main=2,font.lab=2,cex.main=2,cex.lab=1.5,cex.axis=1.5)
   axis(1,at=seq(from=0,to=1400,by=200),labels=seq(from=-200,to=1200,by=200),col="white",las=1,cex.axis=1.5)
   if (i==1)
      conddat = t((((apply(fpow[,which(beh[,2]>0 & beh[,2]<12),],c(1,3),FUN=mean))/norm_all)*100)-100)
   else if (i==2)
      conddat = t((((apply(fpow[,which(beh[,2]>12 & beh[,2]<25),],c(1,3),FUN=mean))/norm_all)*100)-100)
   else if (i==3)
      conddat = t((((apply(fpow[,which(beh[,2]>24 & beh[,2]<37),],c(1,3),FUN=mean))/norm_all)*100)-100)
   else if (i==4)
      conddat = t((((apply(fpow[,which(beh[,2]>36 & beh[,2]<48),],c(1,3),FUN=mean))/norm_all)*100)-100)
   else
      conddat = t((((apply(fpow[,which(beh[,2]==0),],c(1,3),FUN=mean))/norm_all)*100)-100)
   .filled.contour(seq(1,2000,1),seq(1,30,1),conddat,levels=seq(-100,100,200/64),col=colmap[,1])
   abline(v=500,lwd=3,lty=2,col='white')
   range(((fpow/norm_all)*100)-100)
   dev.off()
}