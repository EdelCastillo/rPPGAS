#ifndef ONE_PIXEL
#define ONE_PIXEL

#include <Rcpp.h>
#include "peakInfo.h"
#include "gmmPeak.h"
#include "rGetImzMLData.h"
#include "common_methods.h"
#include <stdlib.h>
#include "noiseestimation.h"

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace Rcpp;
using namespace std; 

class OnePixel
{
public:
  OnePixel(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow, float mzHigh, int pixel);
  ~OnePixel();
  int getGaussians();
  int getGaussians(SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p);
  List getPixelGaussians(int px, float mzLow, float mzHigh);
    
private:  
  //input info to the constructor.
  bool m_continuous;
  int 
  m_nThreads,
  m_pxSupport,
  m_NPixels,
  m_mzLength;
  float     
  m_mzLow,
  m_mzHigh,
  m_SNR;
  double     
    m_massResolution,
    m_maxMassResolution;  
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  PEAK_F_GROUP  m_peakFG;
  GAUSS_SP      m_gaussians;
  
  MASS_SEGMENT  m_massSegment;
  ION_ENTRY     **m_ionEntry_p;
  SPECTRO       m_spectro;
  MASS_RANGE    *m_massRange_p;
  bool          m_enable;
  int     
              m_pixel,
             m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
  
};

#endif
