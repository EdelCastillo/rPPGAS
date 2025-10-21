/*********************************************************************************
 *     peakMatrix
 *     rPPGAS - R package for MSI data processing
 *     Copyright (C) 2025 Esteban del Castillo Pérez (esteban.delcastillo@urv.cat)
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/
#include "peakMatrix.h"

int  globalPeak=0, globalCount=0; //contadores de picos

//'  peakMatrix() link with R.
//'  converts the info in the imzML file into a peak matrix.
//' 
//'  @param ibdFname  absolute reference to the file with the ibd extension.
//'  @param imzML  list with information extracted from the imzML file with import_imzML()
//'  @param params specific parameters
//'   "SNR": signal-to-noise ratio
//'   "minPixelsSupport": minimum percentage of pixels that must support an ion for it to be considered.
//'   "massResolution": mass resolution with which the spectra were acquired.
//'   "maxMassResolution": maximum desired mass resolution.
//'   "noiseMethod": method for estimating noise.
//'  @param mzLow  lower mass to consider
//'  @param mzHigh higher mass to consider
//'  @nThreads number of threads suggested for parallel processing.
//'  @return lista: peakMatrix, massVector, pixelsSupport
//'     peakMatrix: matrix of centroids and the intensity associated with each pixel.
//'     massVector: the mz associated with each column of peakmatrix.
//'     pixelsSupport: number of pixels with intensity > 0

// [[Rcpp::export]]
List peakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int nThreads)
{
  globalPeak=0, globalCount=0;
  PeakMatrix pMatrix(ibdFname, imzML, params, mzLow, mzHigh, nThreads);  
  
  //loads data from a file and converts its peak into Gaussians.
  pMatrix.rawToGaussians(); //parallel processing

  pMatrix.m_massRangeSize=pMatrix.getMassRanges();
  printf("\t\ttotal segments=%d max Gaussians in a spectra=%d\n", 
         pMatrix.m_massRangeSize, pMatrix.m_maxPxGaussians);

  pMatrix.freeMemoryPeak(); //frees up unnecessary memory.
  printf("\t\tvalid_pixels:%d all_spectrum_peak_in_mass_range:%d\n", globalCount, globalPeak);  

  MASS_RANGE massRange;
  massRange.low=mzLow;
  massRange.high=mzHigh;

  //From the Gaussians, the centroids are obtained by kmeans segmentation.
  List ret=pMatrix.massRangeToCentroids(massRange); //parallel processing

  /*
  //Results of significant deviations in mass resolution. 
  float res=1e6/pMatrix.m_massResolution;
  float resB=res*1.5; //deviation of 1.5 times the resolution.
  bool hit=false;
  for(int i=0; i<pMatrix.m_massRangeSize; i++)
    if(pMatrix.m_massRange_p[i].resolution>=resB){hit=true; break;}
  if(hit)
  {
    printf("sample mass resolution (ppm)=%.2f\n", res);
    printf("List of segments with mass resolution >= %.2f ppm\n\tmass range(Da)\tresolution(ppm)\n", 1.5*res);
    int count=0;
    for(int i=0; i<pMatrix.m_massRangeSize; i++)
    {
      if(pMatrix.m_massRange_p[i].resolution>resB)
        printf("[%2d] %.5f/%.5f (%.2f)\n", count++, pMatrix.m_massRange_p[i].low, pMatrix.m_massRange_p[i].high,
             pMatrix.m_massRange_p[i].resolution);
    }
  }
  */
  return ret;
}


//Constructor
//captures input information, allocates memory and initializes.
////////////////////////////////////////////////////////////////////////////////
PeakMatrix::PeakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int nThreads)
{
  m_mzLow=mzLow;
  m_mzHigh=mzHigh;

  //information capture
  Rcpp::DataFrame df;
  df=imzML["run"];
  m_continuous=imzML["continuous_mode"];
  m_NPixels=df.nrows();

  printf("...#px=%d mzLow:%f mzHigh:%f\n",m_NPixels, m_mzLow, m_mzHigh);
  
  NumericVector nv;
  CharacterVector cv;
  
  nv=params["SNR"];
  m_SNR=nv[0];
  if(m_SNR<=0) m_SNR=1;

  nv=params["massResolution"];
  m_massResolution=nv[0];
  
  nv=params["minPixelsSupport"];
  m_pxSupport=nv[0]*m_NPixels/100.0;
  
  nv=params["maxMassResolution"];
  m_maxMassResolution=nv[0]; 
  
  cv=params["noiseMethod"];
  String tmpStr=cv[0];
  const char* SNRmethod=tmpStr.get_cstring(); //conversion C
  
  if     (strcmp(SNRmethod, "estnoise_diff")==0) {m_SNRmethod=1; }
  else if(strcmp(SNRmethod, "estnoise_sd")  ==0) {m_SNRmethod=2; }
  else if(strcmp(SNRmethod, "estnoise_mad") ==0) {m_SNRmethod=3; }
  else {m_SNRmethod=0; printf("unknow noise method: %s so, estnoise_mad is used\n", SNRmethod);}
  m_noiseEst_p=0;
  m_noiseEst_p=new NoiseEstimation(m_SNRmethod, 1, 9);

  //class for accessing imzML files.
  m_getImzMLData_p=0;
  m_getImzMLData_p= new GetImzMLData(ibdFname, imzML);
  
  //memory and its initialization
  m_gaussians_p=0;
  m_gaussians_p=new GAUSS_SP[m_NPixels];
  for(int px=0; px<m_NPixels; px++)
  {
    m_gaussians_p[px].gauss_p=0;
    m_gaussians_p[px].initGauss=0;
  }
  
  //dimensioned
  NumericVector mzLength=df["mzLength"];
  m_maxMzLength=0;
  for(int i=0; i<m_NPixels; i++)
    if(mzLength[i]>m_maxMzLength)m_maxMzLength=mzLength[i]; //maximum spectrum

  //vectors to accommodate a spectrum and its masses.
  for(int i=0; i< MAX_THREADS; i++)
    {
      m_spectro[i].mass_p=0;
      m_spectro[i].int_p=0;
      m_spectro[i].SNR_p=0;
      m_spectro[i].tmpMass_p=0;
      m_spectro[i].tmpInt_p=0;
      m_spectro[i].tmpSNR_p=0;
      m_spectro[i].sort_p=0;
      m_spectro[i].size=0;
      m_spectro[i].thread_p=0;
    }
  
  m_enable=true; //terminates threads if false.
  
  if(nThreads<=0) //parameter control.
    m_nThreads =thread::hardware_concurrency()-1; //a core is released
  else 
    m_nThreads =nThreads;
  if(m_nThreads<=0) m_nThreads=1;
  
  printf("maximum data points of a spectrum: %d\n",m_maxMzLength);
  printf("Threads to use: %d\n", m_nThreads);

  //keeps the info of a spectrum, along with the thread that processes it.
  //memory reservation and initialization.
  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].int_p      =new float[m_maxMzLength];
    m_spectro[i].mass_p     =new float[m_maxMzLength];
    m_spectro[i].SNR_p      =new float[m_maxMzLength];
    m_spectro[i].tmpMass_p  =new float[m_maxMzLength];
    m_spectro[i].tmpInt_p   =new float[m_maxMzLength];
    m_spectro[i].tmpSNR_p   =new float[m_maxMzLength];
    m_spectro[i].sort_p     =new int  [m_maxMzLength];
    m_spectro[i].mutexIn_p  =new std::mutex;
    m_spectro[i].mutexOut_p =new std::mutex;
    m_spectro[i].size=0;
    m_spectro[i].mutexIn_p ->lock();
    m_spectro[i].mutexOut_p->lock();
    m_spectro[i].thread_p=new std::thread(&PeakMatrix::mtGetGaussians, this, i);
  }
  
  //structure initialization for peak.
  m_peakFG_p=0;
  m_peakFG_p= new PEAK_F_GROUP[m_NPixels];
  for(int px=0; px<m_NPixels; px++)
  {
    m_peakFG_p[px].peakF_p=0;
    m_peakFG_p[px].peakU_p=0;
    m_peakFG_p[px].peakFsize=0; 
    m_peakFG_p[px].peakUsize=0;
  }

  m_massSegment.massRange_p=0;
  m_massSegment.segment_p=0;
  m_massSegment.nGaussians_p=0;
  
  m_ionEntry_p=0;
  m_maxPxGaussians=0;
  m_massRange_p=0;
}

//destructor
//free reserved memory
PeakMatrix::~PeakMatrix()
{
//  printf("init PeakMatrix destructor\n");
 if(m_gaussians_p)
  {
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_gaussians_p[px].gauss_p) delete[] m_gaussians_p[px].gauss_p;
    }
    delete []m_gaussians_p;
  }

  freeMemoryPeak();
 
 for(int i=0; i<m_nThreads; i++)
 {
   if(m_spectro[i].mutexIn_p)  {delete m_spectro[i].mutexIn_p; m_spectro[i].mutexIn_p=0;}
   if(m_spectro[i].mutexOut_p) {delete m_spectro[i].mutexOut_p;m_spectro[i].mutexOut_p=0;}
 }
 
  if(m_massSegment.massRange_p)     delete []m_massSegment.massRange_p;
  if(m_massSegment.nGaussians_p)    delete []m_massSegment.nGaussians_p;
  if(m_massRange_p) delete [] m_massRange_p;
//  printf("end PeakMatrix destructor\n");
}

//freeing buffer.
void PeakMatrix::freeMemoryPeak()
{
  //class for accessing imzML files
  if(m_getImzMLData_p) {delete m_getImzMLData_p; m_getImzMLData_p=0;}

  if(m_peakFG_p) //if the peak structure exists
  {
    for(int px=0; px<m_NPixels; px++)
    {
      if(m_peakFG_p[px].peakF_p) 
        {delete [] m_peakFG_p[px].peakF_p; m_peakFG_p[px].peakF_p=0;}
      if(m_peakFG_p[px].peakU_p) 
        {delete [] m_peakFG_p[px].peakU_p; m_peakFG_p[px].peakU_p=0;}
    }
    delete []m_peakFG_p; m_peakFG_p=0;
  }
 
  //spectrum info
  for(int i=0; i<m_nThreads; i++)
   {
    if(m_spectro[i].int_p)      {delete []m_spectro[i].int_p;    m_spectro[i].int_p=0;}
    if(m_spectro[i].mass_p)     {delete []m_spectro[i].mass_p;   m_spectro[i].mass_p=0;}
    if(m_spectro[i].SNR_p)      {delete []m_spectro[i].SNR_p;    m_spectro[i].SNR_p=0;}
    if(m_spectro[i].tmpMass_p)  {delete []m_spectro[i].tmpMass_p;m_spectro[i].tmpMass_p=0;}
    if(m_spectro[i].tmpInt_p)   {delete []m_spectro[i].tmpInt_p; m_spectro[i].tmpInt_p=0;}
    if(m_spectro[i].tmpSNR_p)   {delete []m_spectro[i].tmpSNR_p; m_spectro[i].tmpSNR_p=0;}
    if(m_spectro[i].sort_p)     {delete []m_spectro[i].sort_p;   m_spectro[i].sort_p=0;}
  }

  //the memory reserved for the ion chain is released.
  if(m_ionEntry_p)
    {
      ION_ENTRY *ionEntry2_p, *ionEntry_p=m_ionEntry_p[0];
      while(ionEntry_p)
      {
        ionEntry2_p=ionEntry_p->group;
        if(ionEntry_p->set) {delete[] ionEntry_p->set;  ionEntry_p->set=0;}
        if(ionEntry_p)      {delete ionEntry_p;         ionEntry_p=0;}
        ionEntry_p=ionEntry2_p;
      }
      delete [] m_ionEntry_p; m_ionEntry_p=0;
    }
  
}

//rawToGaussians
//gets the intensity peak and converts them into Gaussians.
//The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
//the SNR info is established for each point of the spectra.
//The generated information is stored in the m_spectro_p structure.
/////////////////////////////////////////////////////////////////////////////////////////////
int PeakMatrix::rawToGaussians()
{
  printf("Processing \n\tphase 1:  to gaussians...(%%): 00 ");
  int vez=1;
  Common common;
  int spSize=0;
  bool hit=false;
  NumericMatrix SNR_Mx;
  NumericVector noiseSize;

  //for all pixels (spectra)
  //For each spectrum, determine the intensity peak, the joined peak, and their Gaussians.
  for(int pixel=0; pixel<m_NPixels; ) 
    {
    //indication of the progress of the process (10% resolution)
    if((float)pixel/(float)m_NPixels>vez*0.1) {if(vez<10) printf("%d ", vez*10); vez++;}
    
    //We load as many spectra as threads are used.
    //Each spectrum is processed by a thread.
    hit=true;
    int nThr=0;
    for(int thr=0; thr<m_nThreads; thr++) //for each thread
      {
      
      while(pixel<m_NPixels) //iterates while the spectrum has length <=2 (eliminates empty spectra).
        {
        spSize=getRawInfo(pixel, thr);//capturing spectra from imzML file
        pixel++;
        //If the spectrum size <=2 its peak are not processed.
        if(spSize>2) break;//spectrum for analysis. There is no information of interest in the pixel spectrum <=3 peak.
        }
      if(spSize>2 && pixel<m_NPixels) 
        {
        nThr++;
        {m_spectro[thr].mutexIn_p->unlock();} //This thread is allowed to run.
        }
      }
    //wait for the conclusion of all threads.
    for(int thr=0; thr<nThr; thr++) //for each thread
      m_spectro[thr].mutexOut_p->lock();
  }
  
  //maximum amount of Gaussians in the spectra.
  for(int px=0; px<m_NPixels; px++)//for all pixels
    if(m_gaussians_p[px].gauss_p && m_gaussians_p[px].size>m_maxPxGaussians) 
      m_maxPxGaussians=m_gaussians_p[px].size;

  //mass range of interest
  MASS_RANGE mRange;
  mRange.low=m_mzLow;
  mRange.high=m_mzHigh;
    
  printf("100\n");
  
  //thread removal
  m_enable=false; //threads end.
  for(int thr=0; thr<m_nThreads; thr++) //are unlocked for the conclusion.
    m_spectro[thr].mutexIn_p->unlock();
  
  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].thread_p->join();
    delete m_spectro[i].thread_p; m_spectro[i].thread_p=0;
  }
  
  return 0;
}

//getRawInfo()
//Loads the full spectrum information associated with a pixel from an imzML file.
//The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
//px: Pixel whose spectrum should be loaded.
//spIndex: Threads that manage it
//Returns the size of the spectrum.
int PeakMatrix::getRawInfo(int px, int spIndex)
{
  //raw info of mass
  int massSize=m_getImzMLData_p->getPixelMassF(px, m_spectro[spIndex].tmpMass_p);//mass vector
  m_spectro[spIndex].size=massSize;
  if(massSize<=0) return 0;
  
  //raw info of intensities  
  int intSize=m_getImzMLData_p->getPixelIntensityF(px, m_spectro[spIndex].tmpInt_p);

  m_spectro[spIndex].pixel=px;

  return massSize;
} 

//mtGetGaussians()
//Parallel processing.
//Peak are delimited and their Gaussians are formed.
//This thread remains active, processing spectra until none remain.
//Each spectrum is converted into Gaussians that can overlap (join).
//spIndex: thread
//Returns -1 on failure, 0 = OK.
int PeakMatrix::mtGetGaussians(int spIndex)
  { 
  IntensityPeak intPeak(m_SNR);
  Common common;
  GAUSS_PARAMS *gaussians_p=0;
  double *centroids_p=0;
  int *centroidsIndex_p=0;
  
  //while execution is enabled.
  while(m_enable)
  {

    m_spectro[spIndex].mutexIn_p->lock(); //permission to continue (synchro)

    //Spectrum conditioning.
    //Full spectra are received and only part of it may be of interest.
    
    int iMzLow=-1, iMzHigh=-1;

    int spSize=m_spectro[spIndex].size; //spectrum size
    //SNR
    m_noiseEst_p->getSNR(m_spectro[spIndex].tmpInt_p, m_spectro[spIndex].size,  m_spectro[spIndex].tmpSNR_p);

    //the spectrum is limited to the range of interest.
    iMzLow =common.nearestIndex(m_mzLow,  m_spectro[spIndex].tmpMass_p, spSize); //low index
    iMzHigh=common.nearestIndex(m_mzHigh, m_spectro[spIndex].tmpMass_p, spSize); //high index
    spSize=iMzHigh- iMzLow +1;

    //fault control
    if(iMzLow<0 || iMzHigh>m_spectro[spIndex].size-1 || iMzHigh-iMzLow>m_spectro[spIndex].size || iMzHigh<iMzLow) 
      {
        m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
        return -1;
      }

    //mass ordination.
    common.sortUpF(m_spectro[spIndex].tmpMass_p+iMzLow, m_spectro[spIndex].sort_p, spSize);

    //The part of interest is extracted from the mass, intensity and SNR vectors.
    for(int i=0; i<spSize; i++) 
      {
      int index=iMzLow+m_spectro[spIndex].sort_p[i];
        m_spectro[spIndex].mass_p[i]=m_spectro[spIndex].tmpMass_p[index];
        m_spectro[spIndex].int_p [i]=m_spectro[spIndex].tmpInt_p [index];
        m_spectro[spIndex].SNR_p [i]=m_spectro[spIndex].tmpSNR_p [index];
      }
    m_spectro[spIndex].size=spSize; //useful range
    //-----------------------------------------------------------------    
   
    //conversion to Gaussians.
    
    if(!m_enable) //end of thread?
      {m_spectro[spIndex].mutexOut_p->unlock(); return 0;}
    int nPeak;
    
    //the peak are extracted from the spectrum (they are delimited by their indices).
    nPeak=intPeak.getPeakList(&m_spectro[spIndex]);
    globalPeak+=nPeak; //peak accumulation
    globalCount++;       //processed spectra
    
    if(nPeak>0) //if there are peak to treat
      {
      int px=m_spectro[spIndex].pixel; //pixel
      
      //single peak info is saved.
      m_peakFG_p[px].peakF_p=new ION_INDEX[nPeak];
      for(int i=0; i<nPeak; i++) //copy peak
      {
        m_peakFG_p[px].peakF_p[i].low =intPeak.getSinglePeak(i).low;
        m_peakFG_p[px].peakF_p[i].max =intPeak.getSinglePeak(i).max;
        m_peakFG_p[px].peakF_p[i].high=intPeak.getSinglePeak(i).high;
      }
      m_peakFG_p[px].peakFsize=nPeak;//number of simple peak
      
      //the information of compound peak (simple joined-overlapping peak) is obtained.
      int nUPeak=intPeak.getCompoundPeakNumber(); 

      //Compound peak information is saved.
      //Each entry is a reference to the initial and final single peak of the merged peak.
      m_peakFG_p[px].peakU_p=new PEAK_UNITED[nUPeak];
      for(int i=0; i<nUPeak; i++)
      {
        m_peakFG_p[px].peakU_p[i].low =intPeak.getCompoundPeak(i).peakLow;
        m_peakFG_p[px].peakU_p[i].high=intPeak.getCompoundPeak(i).peakHigh;
      }
      m_peakFG_p[px].peakUsize=nUPeak;//number of compound peak
      
      //the peak are converted to Gaussians
      gaussians_p=new GAUSS_PARAMS[nPeak]; //temporal copy of Gaussians
      centroids_p=new double[nPeak];       //temporary mass copy
      centroidsIndex_p=new int[nPeak];     //indices to ordered masses

      int nGaussians=getGaussians(px, &m_spectro[spIndex], gaussians_p);

      m_gaussians_p[px].size=0;
 
      if(nGaussians>0 && nGaussians<=nPeak)//if Gaussians exist
      {
        //sorted from lowest to highest (some joined peak may be altered).
        for(int i=0; i<nGaussians; i++)
        {
          centroids_p[i]=gaussians_p[i].mean;
        }
      
        if(common.sortUp(centroids_p, centroidsIndex_p, nGaussians)==0) //the centers are ordered
        {
          //ordered copy to main structure
          for(int i=0; i<nGaussians; i++)
          {
            int index=centroidsIndex_p[i];
            m_gaussians_p[px].gauss_p[i].mean  =gaussians_p[index].mean;
            m_gaussians_p[px].gauss_p[i].sigma =fabs(gaussians_p[index].sigma);
            m_gaussians_p[px].gauss_p[i].weight=gaussians_p[index].weight;
          }
          m_gaussians_p[px].size=nGaussians;
        }
       
      }
      else m_gaussians_p[px].size=0;

    } //end of peak processing

    if(gaussians_p)     {delete [] gaussians_p;       gaussians_p=0;}
    if(centroids_p)     {delete [] centroids_p;       centroids_p=0;}
    if(centroidsIndex_p){delete [] centroidsIndex_p;  centroidsIndex_p=0;}
    
    m_spectro[spIndex].mutexOut_p->unlock(); //end of spectrum processing.
  }

  return 0;
}

//getGaussians()
//Called from a thread.
//Sets the Gaussians on the peak.
//Uses the peak separation information (m_peakFG_p).
//px: pixel.
//spectro: pointer to the spectrum.
//gaussians_p: pointer to the structure containing the Gaussians' parameters.
//Returns the number of Gaussians or a value < 0 on failure.
int PeakMatrix::getGaussians(int px, SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p)
{
  float *intSpectrum_p=spectro_p->int_p, *massSpectrum_p=spectro_p->mass_p;
  int intSize=spectro_p->size;
  
  int minMeanPxMag=1;
  //class for conversion to Gaussians.
  GmmPeak gmmPeak(minMeanPxMag);
  
  int gaussIndex=0; //indices for each Gaussian.
  int nUPeak=m_peakFG_p[px].peakUsize; //#united peak
  
  //memory for predictable Gaussians.
  m_gaussians_p[px].gauss_p=new GAUSS_PARAMS[m_peakFG_p[px].peakFsize];
  
  //for each set of joined peak.
  for(int uPeak=0; uPeak<nUPeak; uPeak++)
  {
    int mPeakLow   =m_peakFG_p[px].peakU_p[uPeak].low;      //lower simple peak.
    int mPeakHigh  =m_peakFG_p[px].peakU_p[uPeak].high;     //upper simple peak.
    int lowMzIndex =m_peakFG_p[px].peakF_p[mPeakLow].low;   //lower  mass index.
    int highMzIndex=m_peakFG_p[px].peakF_p[mPeakHigh].high; //higher mass index.
    int mzSize=highMzIndex-lowMzIndex+1;
    
    int nPeak=mPeakHigh-mPeakLow+1; //number of simple peak into united peak.
    
    //The information of simple peak of magnitude that make up the segment is established.
    gmmPeak.setPeak(&m_peakFG_p[px].peakF_p[mPeakLow], nPeak); 

    //The magnitude information that must be adjusted is established.
    GROUP_F pxMag;
    pxMag.set=intSpectrum_p+lowMzIndex; 
    pxMag.size=mzSize;
    gmmPeak.setMagnitudes(&pxMag); //magnitude peak.
    
    //The simple Gaussians are formed whose sum reproduces the magnitude peak.
    if(gmmPeak.gmmDeconvolution()<0) //deconvolution
    {
      printf("ERROR: limits exceeded. The aim is to deconvolve a peak composed of more than %d Gaussians.\n", DECONV_MAX_GAUSSIAN);
      return -4;
    }
    if(intSize<mzSize) //if the dimensions are not correct.
    {
      //warning: the dimensions of mzAxis and data[magnitudes] must match. 
      printf("ERROR: \n");
      return -5;
    }
    
    int nGauss=gmmPeak.getDeconvNumber(); //# gaussians
    
    bool hit=true;
    GAUSSIAN gaussIn, gaussOut; //if it is required to adapt the mass axis.
    for(int g=0; g<nGauss; g++)
    {
      gaussIn=gmmPeak.getDeconv(g); //get a gaussian.
      
      //unit conversion (scans to Daltos).
      gmmPeak.gaussConversion(&gaussIn, &gaussOut, massSpectrum_p+lowMzIndex, mzSize); 
      
      if(gaussIndex>=m_peakFG_p[px].peakFsize) //control de error
      {
        //warning: memory overflow for gaussians 
        printf("warning: memory overflow for gaussians %d/%d\n", gaussIndex, m_peakFG_p[px].peakFsize); 
        hit=false; break;
      }
      
    //the Gaussians are saved.
    gaussians_p[gaussIndex].mean  =(double)gaussOut.mean;
    gaussians_p[gaussIndex].sigma =(double)gaussOut.sigma;
    gaussians_p[gaussIndex].weight=(double)gaussOut.yFactor*gaussOut.weight;
    gaussIndex++;
    }//end of gaussians
    
    if(hit==false) break;//fallo
  }//end uPeak loop
  return gaussIndex; //# gaussians
}

//massRangeToCentroids()
//Parallel processing.
//1) Sets the number of independent mass ranges (joined peak).
//2) Establishes clusters within those ranges (kmeans segmentation).
//Requires preprocessing by getGaussians()
//Receives the total mass range to consider.
//Returns a list: peakMatrix, massVector, pixelsSupport.
//peakMatrix: Matrix of centroids and the intensity associated with each pixel.
//massVector: The mz associated with each column of the peakMatrix.
//pixelsSupport: Number of pixels with intensity > 0.
List PeakMatrix::massRangeToCentroids(MASS_RANGE massRangeIn)
{
  Common common;
  int indexLow =common.nearestIndexMassRangeLow (massRangeIn.low,  m_massRange_p, m_massRangeSize);
  int indexHigh=common.nearestIndexMassRangeHigh(massRangeIn.high, m_massRange_p, m_massRangeSize);
  int deltaSegment=(indexHigh-indexLow+1)/m_nThreads; //partes iguales para cada thread

  //Memory for necessary structures.
  //Holds mass range information.
  m_massSegment.massRange_p=  new MASS_RANGE[m_nThreads];
  m_massSegment.nGaussians_p= new int[m_nThreads];

  //Ion information with the intensities for each pixel.
  //Linked sequence.
  m_ionEntry_p=new ION_ENTRY*[m_nThreads];
  for(int i=0; i<m_nThreads; i++)
  {
    m_ionEntry_p[i]=new ION_ENTRY;
    m_ionEntry_p[i]->set=0;
    m_ionEntry_p[i]->size=0;
    m_ionEntry_p[i]->group=0;
  }
  
  m_enable=true; //allows the thread loop to operate.
  
  //Mass ranges that each thread must process (linear distribution).
  for(int thr=0; thr<m_nThreads; thr++)
  {
    m_massSegment.massRange_p[thr].low =indexLow+thr*deltaSegment; //Note: these are integer values over floats.
    m_massSegment.massRange_p[thr].high=indexLow+(thr+1)*deltaSegment-1;
  }
  m_massSegment.massRange_p[m_nThreads-1].high=indexHigh; //the last thread may be more loaded.

  //Thread
  //1) Mass ranges containing joined Gaussians are established.
  //2) Segmentation: Centroids and nearest pixels are established.
  //Threads are activated.
  m_enable=true;
  for(int i=0; i<m_nThreads; i++)
    m_spectro[i].thread_p=new std::thread(&PeakMatrix::mtSegmentation, this, i);
  
  printf("\tphase 2:  to centroids.  Wait to finish...\n");

  //wait until the conclusion-
  for(int i=0; i<m_nThreads; i++)
  {
    m_spectro[i].thread_p->join();
    delete m_spectro[i].thread_p; m_spectro[i].thread_p=0;
  }
  //adequacy of results for R.
  int totalIons=0;
  for(int i=0; i<m_nThreads; i++)
  {
    totalIons+=m_totalIons[i];
  }
  NumericMatrix peakMatrix(m_NPixels, totalIons);  //peak matrix.
  NumericVector massVector(totalIons);              //mass vector.
  IntegerVector pixelsSupport(totalIons);           //#Support pixels.
 
  //Each thread contains the chained information of a consecutive part of the ions.
  ION_ENTRY *localIonEntry_p; //input
  int col=0;
  float lastMass=0;
  for(int thr=0; thr<m_nThreads; thr++)
   {
   localIonEntry_p=m_ionEntry_p[thr];  //entry to the info.
   while(localIonEntry_p->group) //as long as the next link is not zero.
   {
     //It may happen that the initial masses of a thread are lower than the last masses of 
     //the previous thread, and they should be discarded. 
     //This is a consequence of dealing with joined Gaussians.     
     if(localIonEntry_p->mass>lastMass)                 //the masses are ordered.
     {
       massVector(col)=localIonEntry_p->mass;           //centroid
       pixelsSupport(col)=localIonEntry_p->size;        //#pixels that support the ion
       for(int row=0; row<m_NPixels; row++)             //for each pixel
         peakMatrix(row, col)=localIonEntry_p->set[row]; //intensities
       col++;
     }
     if(localIonEntry_p->group->group==0) //the last link contains no info.
       lastMass=localIonEntry_p->mass;    //last mass of the thread
     localIonEntry_p=localIonEntry_p->group; //next item
   }
 }
 
 List ret=List::create(Named("peakMatrix")=peakMatrix, Named("mass")=massVector, Named("pixelsSupport")=pixelsSupport);
 return ret;
 }

//getCentroidsIntoRange()
//Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
//massRange: Mass range from which to extract the Gaussians.
//gaussians_p: Requested Gaussians.
//Returns the number of Gaussians.
int PeakMatrix::getCentroidsIntoRange(MASS_RANGE massRange, float **gaussians_p)
{
  int count=0, count2=0;
  for(int px=0; px<m_NPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    for(int i=0 ; i<m_gaussians_p[px].size; i++) //improve!!!!
    {
      if(m_gaussians_p[px].gauss_p[i].mean>= massRange.low && 
         m_gaussians_p[px].gauss_p[i].mean<= massRange.high)
      {
        gaussians_p[count][0]=m_gaussians_p[px].gauss_p[i].mean;
        gaussians_p[count][1]=fabs(m_gaussians_p[px].gauss_p[i].sigma);
        gaussians_p[count][2]=m_gaussians_p[px].gauss_p[i].weight;
        gaussians_p[count][3]=(float)px;
        count++;
      }
    }
  }
  return count;
}

//getCentroidsNumberIntoRange()
//Returns the number of Gaussians in a mass range from the information in m_gaussians_p.
//massRange: Mass range from which to extract Gaussians.
//Returns the number of Gaussians.
int PeakMatrix::getCentroidsNumberIntoRange(MASS_RANGE massRange)
{
  int count=0;
  for(int px=0; px<m_NPixels; px++)//for all pixels
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //pixel without Gaussians
    for(int i=0 ; i<m_gaussians_p[px].size; i++) //improve!!!!
    {
      if(m_gaussians_p[px].gauss_p[i].mean>= massRange.low && 
         m_gaussians_p[px].gauss_p[i].mean<= massRange.high)
        count++;
    }
  }
  return count;
}

//getMassRanges()
//establishes the mass segments where overlapping Gaussians exist.
//generates a Boolean mass axis with a resolution of 1/4 of the spectrometer's mass resolution.
//if any Gaussian invades its space, sets the corresponding +/-3*sigma checkboxes to true.
//Very wide Gaussians (sigma > 5*deltaMass) are discarded.
//the information is stored in the array pointed to by m_massRange_p.
//returns the number of segments.
int PeakMatrix::getMassRanges()
{
  double highMass, lowMass;
  int iLow, iHigh;
  double deltaMass=m_mzLow/(4*m_massResolution); //delta=1/4 of the minimum mass increment of the spectrometer
  int massAxisSize=1+(m_mzHigh-m_mzLow)/deltaMass;
  
  bool *massAxis=0;
  massAxis= new bool[massAxisSize];
  for(int i=0; i<massAxisSize; i++) massAxis[i]=false;
  
  double deltaMass2;
  
  //for all elements of the spectrum
  for(int px=0; px<m_NPixels; px++)
  {
    if(m_gaussians_p[px].gauss_p==0) continue; //if this entry does not contain info.
    for(int i=0; i<m_gaussians_p[px].size; i++) 
    {
      //Very wide Gaussians (sigma>5*deltaMass) are discarded.
      deltaMass2=(m_gaussians_p[px].gauss_p[i].mean)/m_massResolution;
      if(fabs(m_gaussians_p[px].gauss_p[i].sigma)>5*deltaMass2) //sigma >>
      {
        continue;
      }
    //masses occupied by the Gaussian.
    highMass=m_gaussians_p[px].gauss_p[i].mean+3*m_gaussians_p[px].gauss_p[i].sigma;
    lowMass =m_gaussians_p[px].gauss_p[i].mean-3*m_gaussians_p[px].gauss_p[i].sigma;
    if(lowMass<m_mzLow) lowMass =m_mzLow;
    if(highMass>m_mzHigh) highMass=m_mzHigh;

    iLow =(lowMass -m_mzLow)/deltaMass;
    iHigh=(highMass-m_mzLow)/deltaMass;
    for(int j=iLow; j<=iHigh; j++) massAxis[j]=true;
    }
  }
  
  //separation into segments, observing unused spaces.
  bool into=false;
  int init=-1, end=-1;
  int count=0, index=0;
  for(int i=0; i<massAxisSize; i++)   //#packages (contiguous masses)
  {
    if(massAxis[i] && !into)          //home package
      {into=true; init=i;}            //inside the package
    if(!massAxis[i] && into)          //end of package
      {into=false; end=i-1; count++;} //out of the package
  }
  m_massRange_p=new MASS_RANGE[count];
  
  into=false;
  for(int i=0; i<massAxisSize; i++)
  {
    if(massAxis[i] && !into) 
      {into=true; init=i;}
    if(!massAxis[i] && into) //stands out from the pack.
    {
      into=false; end=i-1; 
      m_massRange_p[index].low =m_mzLow+init*deltaMass;//package ends
      m_massRange_p[index].high=m_mzLow+end *deltaMass;
      index++;
    }
  }
  if(massAxis) delete [] massAxis;
  return index; 
}

//setGaussiansIntoSegments()
//Determines the maximum number of Gaussians over the given mass intervals and all pixels.
//massRange: Mass range to consider.
//Returns the maximum value.
int PeakMatrix::setGaussiansIntoSegments(MASS_RANGE massRange)
{
  int maxGaussians=0, nGaussians;
  float massLow, massHigh;
  for(int mr=massRange.low; mr<=massRange.high; mr++) //for the entire mass range
  {
    massLow =m_massRange_p[mr].low;  //mass range of this segment
    massHigh=m_massRange_p[mr].high;
    nGaussians=0;
    
    for(int px=0; px<m_NPixels; px++)//for all pixels
    {
      if(m_gaussians_p[px].gauss_p==0) continue;  //no info
      for(int g=0; g<m_gaussians_p[px].size; g++)
      {
        if(m_gaussians_p[px].gauss_p[g].mean>=massLow && m_gaussians_p[px].gauss_p[g].mean<=massHigh)
          nGaussians++;
        if(nGaussians>maxGaussians) maxGaussians=nGaussians;
      }
    }
    m_massRange_p[mr].nGaussians=nGaussians; 
  }
  return maxGaussians;
}

//mtSegmentation()
//Kmeans segmentation in parallel processing.
//For each mass segment, centroids are generated.
//The information is stored in the array pointed to by m_ionEntry_p.
void PeakMatrix::mtSegmentation(int thrIndex)
{
  float mass=0;
  float **gaussians_p=0; //pointers to info to return.
  int totalIons=0;
  Common common;
  MASS_RANGE massRange, massRange2;
  massRange.low=m_mzLow;
  massRange.high=m_mzHigh;
  int count=0;
  ION_ENTRY *localIonEntry_p=m_ionEntry_p[thrIndex];
  
  int mRange=0; 
  float maxRange=0;
  int maxRangeIndex=-1;
  int iLow =m_massSegment.massRange_p[thrIndex].low;
  int iHigh=m_massSegment.massRange_p[thrIndex].high;
  
//Determines the maximum cluster number for the mass range to be evaluated.
//Determines the largest segment in the range.
//A maximum value for clusters is understood to be one in which all Gaussians are
//as close as the mass resolution allows.  
for(int i=iLow; i<=iHigh; i++) 
  {
    float range=m_massRange_p[i].high- m_massRange_p[i].low;
    if(range>maxRange) {maxRange=range; maxRangeIndex=i;} //greater range to evaluate
  }
  //minimum resolution value in this mass range (Da).
  double massRes=m_massRange_p[iLow].low/m_massResolution;//minimum resolution for the thread
  int maxClustersThr=ceil(maxRange/massRes); //maximum clusters

  bool *seg_p=0;
  int *tmpIndex_p=0, *sortedIndex_p=0;;
  float *prob_p=0, *mass_p=0;
  double *tmpMass_p=0;
  
  seg_p=new bool[maxClustersThr];
  tmpMass_p=new double[maxClustersThr];
  tmpIndex_p=new int[maxClustersThr];
  prob_p=new float[maxClustersThr];
  mass_p=new float[maxClustersThr];
  sortedIndex_p=new int [maxClustersThr];
  
  //maximum Gaussians in the given mass range.
  int maxGNumber=setGaussiansIntoSegments(m_massSegment.massRange_p[thrIndex]);
  gaussians_p=new float*[maxGNumber];
  for(int i=0; i<maxGNumber; i++)
  {
    gaussians_p[i]=0;
    gaussians_p[i]=new float[4];
  }
  float *massVect_p=0;
  massVect_p=new float[maxGNumber];
    
  int iMassRangeLow=  m_massSegment.massRange_p[thrIndex].low;
  int iMassRangeHigh= m_massSegment.massRange_p[thrIndex].high;
  double newMassRes;
  
  //main loop
  //Each isolated mass packet is segmented.
  /////////////////////////////////////////////////////////////////  
  for(int mrIndex=iMassRangeLow; mrIndex<=iMassRangeHigh; mrIndex++)
  {
    if(!m_enable) break;
    massRange2.low =m_massRange_p[mrIndex].low; //mass segment to be evaluated
    massRange2.high=m_massRange_p[mrIndex].high;
    m_massRange_p[mrIndex].resolution=0;
    
    int rangeNCenters;
    //Capture of the centroids in the range of masses to be considered.
    rangeNCenters=getCentroidsIntoRange(massRange2, gaussians_p);
    if(rangeNCenters<m_pxSupport) continue; //minimum size control.

    //mass vector for the current range.
    for(int i=0; i<rangeNCenters; i++)
      massVect_p[i]=gaussians_p[i][0]; //centroids

    bool hit, hit2;
    int maxIter=40;
    float convergencia=1e-5;
    
    double massRes=massRange2.low/m_massResolution; //resolution (Da)
    double minDistance;
    if(m_maxMassResolution==0) {minDistance=0;} //no peak will be rejected.
    else
      {minDistance=massRange2.low/m_maxMassResolution;}//resolution (Da)
    newMassRes=massRes;
    double massResSqr=massRes*massRes; //Da*Da
    int nClusters, iter=0;
    bool R, S;
    hit=true;
    int maxMasterIter=10; //maximum coarse iterations.
    
    //It is assumed that the maximum number of clusters coincides with the case where there are 
    //Gaussians in all possible spaces.
    int maxClustersRange=ceil((massRange2.high-massRange2.low)/massRes);
    if(maxClustersRange<0) {continue;}
    else if(maxClustersRange==0) {maxClustersRange=1;} //single mass range case.
    
    //Range mass histogram.
    //Used to initialize the masses associated with each cluster based on probabilities (kmeans).    int nBins=maxClustersRange;
    int nBins=maxClustersRange;
    Histogram histo(massVect_p, rangeNCenters, nBins);
    if(histo.m_error) continue;

    for(int bin=0; bin<nBins; bin++) //relative Gaussians in every possible space.
      prob_p[bin]=histo.getBinProbability(bin);
    
    //indexes ordered from highest to lowest probability.
    common.sortDownF(prob_p, sortedIndex_p, nBins);
    
    //class for segmentation
    KmeansR kmeans(massVect_p, rangeNCenters, maxIter, convergencia);
    if(kmeans.m_error==1) continue;

    for(int i=0; i<maxClustersRange; i++) //initialization masses for clusters.
      mass_p[i]=histo.getBinMeanValue(sortedIndex_p[i]); 

    while(true) //Iterate until goals are reached, relaxing constraints if necessary.
    {
      //Iterates until all clusters reach the desired resolution or are small groups
      //with each iteration, the number of clusters increases by 1.
      for(nClusters=1; nClusters<=maxClustersRange; nClusters++) //clusters search
      {
        if(kmeans.getClusters(nClusters, mass_p)!=0) //clusters are obtained.
          {hit=false; break;} //if algorithm fails
          
        hit=true;
        //The validity of the obtained clusters is evaluated:
        //All large clusters must have low dispersion.
        //Large if the cluster >= pxSupport.
        //Low dispersion if it is < mass resolution.
        for(int c=0; c<nClusters; c++) 
        {
          int cPixels=kmeans.m_kStruct.clusters_p[c].data.size; //px that support
          //dispersion measure.
          double res=kmeans.m_kStruct.clusters_p[c].withinss/cPixels; //mean square distance.
          //res=res*0.9;
          S=cPixels>=m_pxSupport;         //'1' if group with many support pixels.
          R=res<massResSqr;               //'1' if the group dispersion is low.
          seg_p[c]=false;                 //deficient group, by default.
          if(S && R) {seg_p[c]=true;}     //suitable group.
          else if(S && !R) {hit=false; break;} //requires iteration to improve.
          //If !S is marked as undersized by default. It is later discarded.
        }
        if(hit) {break;} //all groups of size>minimum have low dispersion.
      }//end for() for nClusters
      
    //If the desired conditions are not met, the constraints are relaxed and iterates completely.
    if(!hit) 
    {
        iter++;
        if(iter>=maxMasterIter) //limit control. Ensures that it ends.
        {
          printf("Warning: clustering failed in mass range %.5f/%.5f with %d clusters (maxClusters=%d)\n", 
                 m_massRange_p[mrIndex].low, m_massRange_p[mrIndex].high, kmeans.m_kStruct.nClusters, maxClustersRange);
          break;
        }
        //se relaja la resolución 
        if(iter<=4) newMassRes=massRes*(1.0+iter/2.0);  //massRes/2 increments
        else newMassRes+=massRes*(iter-4);              //exponential increases.
        massResSqr=newMassRes*newMassRes;
        continue;
    }
    else break;
  }//while end

    //The resolution with which the segmentation has been resolved is noted.
    m_massRange_p[mrIndex].resolution=newMassRes*1e6/m_massRange_p[mrIndex].low;
    
    if(iter>=maxMasterIter) {kmeans.freeClusters();continue;}//voided segment.
    
    //Here, two cases can occur:
    //1) all groups have low dispersion and are of good size (optimal).
    //2) there are groups with high dispersion and they are small (they are discarded),
    // the rest have low dispersion and are large  
    
    //increasing mass ordering
    for(int c=0; c<kmeans.m_kStruct.nClusters; c++)
    {
      tmpMass_p[c]=kmeans.m_kStruct.clusters_p[c].center;
    }

    common.sortUp(tmpMass_p, tmpIndex_p, kmeans.m_kStruct.nClusters);

    //filter 1: discard poorly supported clusters -> "minPixelsSupport" parameter
    //goodIndex maintains indexes to 'good' clusters
    int goodIndex[kmeans.m_kStruct.nClusters];
    int goodIndexSize=0;
    for(int c=0; c<kmeans.m_kStruct.nClusters; c++)
    {
      int sortIndex=tmpIndex_p[c];
      if(!seg_p[sortIndex]) 
        {continue;} //the poorly supported group is excluded.
      goodIndex[goodIndexSize++]=sortIndex;
     }
    
    //filter 2: clusters too close -> "maxMassResolution" parameter.
    //The proximity of the masses is analyzed. If they are too close,
    //the smallest cluster is discarded.
    int goodIndexEnd[goodIndexSize];
    double centerA, centerB;
    int indexA=goodIndex[0], indexB, goodIndexEndSize=0;
    
    if(goodIndexSize==1)//if only one cluster exists.
      goodIndexEnd[goodIndexEndSize++]=goodIndex[0]; 
    
    //For all clusters, they are analyzed in sequence, since they are ordered.
    for(int c=1; c<goodIndexSize; c++)
    {
      indexB=goodIndex[c];//next cluster
      centerA=kmeans.m_kStruct.clusters_p[indexA].center;//previous mass.
      centerB=kmeans.m_kStruct.clusters_p[indexB].center;//next mass
      if(fabs(centerB-centerA) < minDistance) //very close centers: the small one is ruled out.
      {
        int sizeA=kmeans.m_kStruct.clusters_p[indexA].data.size;
        int sizeB=kmeans.m_kStruct.clusters_p[indexB].data.size;
        if(sizeB>sizeA)indexA=indexB; //the index points to the largest and iterates.
        if(c+1==goodIndexSize) goodIndexEnd[goodIndexEndSize++]=indexA; //last
      }
      else //distant centers.
      {
        goodIndexEnd[goodIndexEndSize++]=indexA; 
        if(c+1==goodIndexSize) goodIndexEnd[goodIndexEndSize++]=indexB; //last
        indexA=indexB; //to iterate
      }
    }

    //The results are saved.
    //Memory is created to store the cluster information about the peak array.
    //These are structures linked by pointers so they can grow.
    for(int c=0; c<goodIndexEndSize; c++)
    {
      int sortIndex=goodIndexEnd[c]; //index to cluster
      
      localIonEntry_p->set=new float[m_NPixels]; //memory for intensities
      for(int i=0; i<m_NPixels; i++) //reset all intensities. 
        localIonEntry_p->set[i]=0.0;
      
      localIonEntry_p->mass=kmeans.m_kStruct.clusters_p[sortIndex].center; //centroide
      int size=kmeans.m_kStruct.clusters_p[sortIndex].data.size; //cluster size
      localIonEntry_p->size=size; //intensities
      totalIons++; //input counter (centroides)
      
      for(int i=0; i<size; i++) //copy of intensities
      {
        int index=kmeans.m_kStruct.clusters_p[sortIndex].data.set[i]; //Gaussian index
        int px=round(gaussians_p[index][3]); //pixel of the Gaussian.
        localIonEntry_p->set[px]=gaussians_p[index][2]; //intensities
      }
      
      localIonEntry_p->group=new ION_ENTRY; //new entry for the ion.
      localIonEntry_p=localIonEntry_p->group; //the pointer is updated
      localIonEntry_p->group=0; //initializes the next element.
      localIonEntry_p->set=0;
      localIonEntry_p->size=0;
    }
  } //end of analysis of a mass segment.
  
  m_totalIons[thrIndex]=totalIons;
  
  //reserved memory is freed.
  if(seg_p) delete [] seg_p;
  if(tmpMass_p) delete [] tmpMass_p;
  if(tmpIndex_p) delete [] tmpIndex_p;
  if(prob_p) delete [] prob_p;
  if(mass_p) delete [] mass_p;
  if(sortedIndex_p) delete []sortedIndex_p;
  
  if(massVect_p) delete [] massVect_p;
  if(gaussians_p)
  {
    for(int i=0; i<maxGNumber; i++)
      if(gaussians_p[i]) delete [] gaussians_p[i];
    delete []gaussians_p;
  }
}


