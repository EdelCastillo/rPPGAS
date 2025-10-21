
#/mnt2/MALDI/Cerebellum_30_60_120k/NoAlineado/C_30k
#' @export
rPPGASquality<-function(SNR_test, pxSupport_test, SNR_ctr=1, pxSupport_ctr=1, referenceMassResolution, histo=FALSE)
{
  #30k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_30k/mz30_snr%.1f_px%.1f_1.0.RData", SNR_test, pxSupport_test)
  load(txt)
  G30=mz
  
  #60k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_60k/mz60_snr%.1f_px%.1f_1.0.RData", SNR_test, pxSupport_test)
  load(txt)
  G60=mz
  
  #120k  (Master)
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_120k/mz120_snr%.1f_px%.1f.RData", SNR_ctr, pxSupport_ctr)
  load(txt)
  G120=mz
  
  cat("Master:", txt, "\n")
  
  logic=G30>=300 & G30<=900
  G30=G30[logic]
  logic=G60>=300 & G60<=900
  G60=G60[logic]
  logic=G120>=300 & G120<=900
  G120=G120[logic]
  
  fq30_snr   <-fitQuality (G120, G30, referenceMassResolution,  histo)
  fq60_snr   <-fitQuality (G120, G60, referenceMassResolution,  histo)
  
  # logic=fqPere30_snr[,4]>=1
  # maxDiff1=fqPere30_snr[logic,1]
  # maxDiff2=fqPere30_snr[logic,2]
  # fqPere30_snr=fqPere30_snr[!logic,]
  # badSize=length(maxDiff1);
  # totalSize=length(fqPere30_snr[,1])
  # cat(sprintf("%d masas difieren significativamente (%.2f%%)\n", badSize, badSize*100/totalSize));
  
  
  #Histograma de las desviaciones
  if(histo==TRUE)
  {
    hist(totalDiff_rMSI2, main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
  }
  cat("\n");
  
  cat("About rPPGAS peaks 30k vs 120k:\n")
  cat("-------------------------------\n")
  txt1=sprintf("    size=%d", length(fq30_snr[,3]));
  m=mean(fq30_snr[,3]); sigma=sd(fq30_snr[,3]); md=median(fq30_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  
  size30=length(fq30_snr[,5])
  logic=fq30_snr[,5]==1
  repes30=length((fq30_snr[,5])[logic])
  txt3=sprintf("    repes30=%d (%.1f%%)", repes30, 100*repes30/size30)
  cat(txt1, txt2, txt3, "\n\n");
  
  cat("About rPPGAS peaks 60k vs 120k:\n")
  cat("-------------------------------\n")
  txt1=sprintf("    size=%d", length(fq60_snr[,3]));
  m=mean(fq60_snr[,3]); sigma=sd(fq60_snr[,3]); md=median(fq60_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  
  size60=length(fq60_snr[,5])
  logic=fq60_snr[,5]==1
  repes60=length((fq60_snr[,5])[logic])
  txt3=sprintf("    repes60=%d (%.1f%%)", repes60, 100*repes60/size60)
  cat(txt1, txt2, txt3, "\n");
  
  #Histograma de las desviaciones
  if(histo==TRUE)
  {
    hist(totalDiff, main="Histogram of deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
  }
  cat("\n");
  
}


#/mnt2/MALDI/Cerebellum_30_60_120k/NoAlineado/C_30k
#' @export
qualityTest<-function(SNR_test, pxSupport_test, master="rPPGAS", SNR_ctr=1, pxSupport_ctr=1, referenceMassResolution, histo=FALSE)
{
  #  if(SNR!=1 & SNR!=2 & SNR!=3 & SNR!=5 & SNR!=7) {print("Warning: unknown SNR; expected:1,2,3,5,7");  return();}
  
  #30k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_30k/g30_snr%.1f_px%.1f.RData", SNR_test, pxSupport_test)
  load(txt)
  
  #60k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_60k/g60_snr%.1f_px%.1f.RData", SNR_test, pxSupport_test)
  load(txt)
  
  #120k  (Master)
  if(master=="rPPGAS")
  {
    #    txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_30k/g30_snr%.1f_px%.1f.RData", SNR_ctr, pxSupport_ctr)
    #    load(txt)
    #    G120=gauss30$mass
    
    txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_120k/g120_snr%.1f_px%.1f.RData", SNR_ctr, pxSupport_ctr)
    load(txt)
    G120=gauss120mass
  }
  else if(master=="rMSI2")
  {
    txt=sprintf("/home/esteban/MALDI/PPGAS_local/rMSI2_test/C_120k/rMSI2_snr%.1f_px%.1f.RData", SNR_ctr, pxSupport_ctr)
    load(txt)
    G120=rMSI2_peak
  }
  else if(master=="mean")
  {
    txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_120k/g120Mean_snr0.17.RData")
    load(txt)
    G120=gauss120mass
  }
  else 
  {
    cat(sprintf("unknow master %s (valids: rPPGAS, rMSI2 & mean); rPPGAS used\n", master))
    txt=sprintf("/home/esteban/MALDI/PPGAS_local/rPPGAS_test/C_120k/g120_snr%.1f_px%.1f.RData", SNR_ctr, pxSupport_ctr)
    load(txt)
    G120=gauss120mass
  }
  
  cat("Master:", txt, "\n")
  
  #rMSI2 30k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rMSI2_test/C_30k/rMSI2_snr%.1f_px%.1f.RData", SNR_test, pxSupport_test)
  load(txt)
  rMSI2_peak30=rMSI2_peak
  
  #rMSI2 60k
  txt=sprintf("/home/esteban/MALDI/PPGAS_local/rMSI2_test/C_60k/rMSI2_snr%.1f_px%.1f.RData", SNR_test, pxSupport_test)
  load(txt)
  rMSI2_peak60=rMSI2_peak
  
  rMSI2_peak30=rMSI2_peak30[rMSI2_peak30>=300 & rMSI2_peak30<=900]
  rMSI2_peak60=rMSI2_peak60[rMSI2_peak60>=300 & rMSI2_peak60<=900]
  
  logic=is.nan(gauss30$mass)
  G30=gauss30$mass[!logic]
  logic=is.nan(gauss60$mass)
  G60=gauss60$mass[!logic]
  #  logic=is.nan(peaksMx120$mass)
  #  G120=peaksMx120$mass[!logic]
  
  
  #  G30=g30
  #  G30=gauss30$mass
  logic=G30>=300 & G30<=900
  G30=G30[logic]
  logic=G60>=300 & G60<=900
  G60=G60[logic]
  logic=G120>=300 & G120<=900
  G120=G120[logic]
  
  #gauss120=peaksMx120$mass;
  #logic=gauss120>=300 & gauss120<=900
  #G120=gauss120[logic]
  
  fqPere30_snr <-fitQuality (G120, rMSI2_peak30, referenceMassResolution,  histo)
  fqPere60_snr <-fitQuality (G120, rMSI2_peak60, referenceMassResolution,  histo)
  fq30_snr   <-fitQuality (G120, G30,       referenceMassResolution,  histo)
  fq60_snr   <-fitQuality (G120, G60,       referenceMassResolution,  histo)
  
  # logic=fqPere30_snr[,4]>=1
  # maxDiff1=fqPere30_snr[logic,1]
  # maxDiff2=fqPere30_snr[logic,2]
  # fqPere30_snr=fqPere30_snr[!logic,]
  # badSize=length(maxDiff1);
  # totalSize=length(fqPere30_snr[,1])
  # cat(sprintf("%d masas difieren significativamente (%.2f%%)\n", badSize, badSize*100/totalSize));
  
  #Resultados para rMSI2
  cat("RESULTS OF MASS DEVIATIONS (300-900 Da)\n");
  cat("About rMSI2 peaks 30k vs 120k:\n")
  cat("------------------------------\n")
  txt1=sprintf("    size=%d", length(fqPere30_snr[,3])); 
  m=mean(fqPere30_snr[,3]); sigma=sd(fqPere30_snr[,3]); md=median(fqPere30_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); 
  
  sizePere30=length(fqPere30_snr[,5])
  logic=fqPere30_snr[,5]==1
  repesPere30=length((fqPere30_snr[,5])[logic])
  txt3=sprintf("    repes=%d (%.1f%%)\n", repesPere30, 100*repesPere30/sizePere30)
  cat(txt1, txt2, txt3, "\n");
  
  
  cat("About rMSI2 peaks 60k vs 120k:\n")
  cat("------------------------------\n")
  txt1=sprintf("    size=%d", length(fqPere60_snr[,3]));
  m=mean(fqPere60_snr[,3]); sigma=sd(fqPere60_snr[,3]); md=median(fqPere60_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); 
  
  sizePere60=length(fqPere60_snr[,5])
  logic=fqPere60_snr[,5]==1
  repesPere60=length((fqPere60_snr[,5])[logic])
  txt3=sprintf("    repes60=%d (%.1f%%)", repesPere60, 100*repesPere60/sizePere60)
  cat(txt1, txt2, txt3, "\n");
  
  #Histograma de las desviaciones
  if(histo==TRUE)
  {
    hist(totalDiff_rMSI2, main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
  }
  cat("\n");
  
  cat("About rPPGAS peaks 30k vs 120k:\n")
  cat("-------------------------------\n")
  txt1=sprintf("    size=%d", length(fq30_snr[,3]));
  m=mean(fq30_snr[,3]); sigma=sd(fq30_snr[,3]); md=median(fq30_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  
  size30=length(fq30_snr[,5])
  logic=fq30_snr[,5]==1
  repes30=length((fq30_snr[,5])[logic])
  txt3=sprintf("    repes30=%d (%.1f%%)", repes30, 100*repes30/size30)
  cat(txt1, txt2, txt3, "\n\n");
  
  cat("About rPPGAS peaks 60k vs 120k:\n")
  cat("-------------------------------\n")
  txt1=sprintf("    size=%d", length(fq60_snr[,3]));
  m=mean(fq60_snr[,3]); sigma=sd(fq60_snr[,3]); md=median(fq60_snr[,3]);
  txt2=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  
  size60=length(fq60_snr[,5])
  logic=fq60_snr[,5]==1
  repes60=length((fq60_snr[,5])[logic])
  txt3=sprintf("    repes60=%d (%.1f%%)", repes60, 100*repes60/size60)
  cat(txt1, txt2, txt3, "\n");
  
  #Histograma de las desviaciones
  if(histo==TRUE)
  {
    hist(totalDiff, main="Histogram of deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
  }
  cat("\n");
  
}

#statisticalQuality()
# determina la desviación del vector mzTest respeto al vector mzRef
# retorna valores estadísticos de las desviaciones
# resolución se usa para identificar masas repetidas (comparten masa de referencia)
#mzRef y mzTest provienen de cardinalTest()
#' @export
statisticalQuality<-function(mzRef, mzTest, resolution)
{
  mzRef=mzRef[mzRef>=300 & mzRef<=900]
  mzTest=mzTest[mzTest>=300 & mzTest<=900]
  txt1=sprintf("    refSize=%d  testSize=%d", length(mzRef), length(mzTest))
  
  Mx=fitQuality(mzRef, mzTest, resolution)
  m=mean(Mx[,3])
  sigma=sd(Mx[,3])
  md=median(Mx[,3])
  logic=Mx[,5]==1
  repes=length((Mx[,5])[logic])
  
  txt2=sprintf("    mean=%.3f  sigma=%.3f  median=%.3f  repes:%.0f (%.1f%%)", m, sigma, md, repes, 100*repes/length(mzTest));
  cat(txt1, txt2)
}

fitQuality<-function(refCentroids, testCentroids, referenceMassResolution, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  testLength=length(testCentroids);
  refLength=length(refCentroids);

  deviation <-matrix(nrow = testLength, ncol = 5);
  deviation[,1]=rep(0, times=testLength);
  deviation[,2]=rep(0, times=testLength);
  deviation[,3]=rep(0, times=testLength);
  deviation[,4]=rep(0, times=testLength);
  deviation[,5]=rep(0, times=testLength);
  
  for(iPk in 1:testLength) #para cada pico del test
  {
    offset=0;
#    testPPM=getMassResolution(PereCentroids[iPk], massAxis)
#    testPPM=testPPM$minDelta_PPM;
    testPPM=referenceMassResolution;
     testMass=testCentroids[iPk];
 # print(testMass)
     
    retMass<-nearestValue(testMass, refCentroids);
    if(retMass==-1) {offset=-1;}
    else {offset<-abs(testMass-retMass)}
    ppm<-1e6*offset/testMass;
    
    deviation[iPk, 1]=testMass;
    deviation[iPk, 2]=retMass;
    deviation[iPk, 3]=ppm;
    
    if(ppm>1.5*testPPM) #desviación > 1.5*resolución de masa (low resolution)
    {deviation[iPk, 4]=2;}
    else if(ppm>testPPM) #desviación >1 && <= 1.5*resolución 
    {deviation[iPk, 4]=1;}
    else                #desviación <=1*resolución 
    {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2]) #comparten misma masa de ref.
      {deviation[iPk, 5]=1;}
    
  }
  if(histo==TRUE)
  {
    hist(deviation[, 3], main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
    #  legend("topright", legend="SNR=1");
  }
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe")
  return(deviation);
}



#' nearestValue
#' Return the nearest value in data
#' Algoritmo de aproximaciones sucesivas
#' 
#' @param value -> reference value
#' @param data  -> array of sort values
#'
#' @return nearest value in data to value; -1 if value es out of range
#' @export
#' 
nearestValue<-function(value, data)
{
  indexLow<-1;
  indexHigh<-length(data);
  
  if(indexHigh==indexLow) return(data[1]);
  if(indexHigh==indexLow+1)
  {
    if(indexLow==-1)indexLow=0;
    if(value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
    else {return(data[indexHigh]);}
  }
  
  if(indexLow!=-1 & value==data[indexLow])       return(value);
  if(indexLow!=-1 & value<data[indexLow])       {return(data[indexLow]) ;}
  else if(value>data[indexHigh]) {return(data[indexHigh]);}
  else if(value==data[indexHigh]) return(value);
  
  while(1)
  {
    indexCenter<-round((indexHigh+indexLow)/2);
    if(value==data[indexCenter]) return(value);
    if(value<data[indexCenter]) {indexHigh<-indexCenter; }
    else {indexLow <-indexCenter;}
    if(indexHigh==indexLow+1)
    {
      if(indexLow!=-1 & value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
      else {return(data[indexHigh]);}
    }
  }
}

