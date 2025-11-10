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
#ifndef R_PEAK_MATRIX
#define R_PEAK_MATRIX

#include <Rcpp.h>

#include "peakInfo.h"
#include "gmmPeak.h"
#include <time.h>
#include <sys/time.h>    
#include "rGetImzMLData.h"
#include "kmeansR.h"
#include <thread>
#include <mutex>
#include <chrono>
#include "common_methods.h"
#include <stdlib.h>
#include "noiseestimation.h"
#include "histogram.h"

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace Rcpp;
using namespace std; 


//class to obtain a peak matrix from the imzML file information.
class PeakMatrix
{
public:
  
  //Constructor
  //captures input information, allocates memory and initializes.
 PeakMatrix(const char* ibdFname, Rcpp::List imzML, Rcpp::List params, float mzLow=0, float mzHigh=0, int nThreads=0);
  
  //destructor
  //free reserved memory
  ~PeakMatrix();
  
  //freeing buffer.
  void freeMemoryPeak();
  
  //mtGetGaussians()
  //Parallel processing.
  //Peak are delimited and their Gaussians are formed.
  //This thread remains active, processing spectra until none remain.
  //Each spectrum is converted into Gaussians that can overlap (join).
  //spIndex: thread
  //Returns -1 on failure, 0 = OK.
  int  mtGetGaussians(int spIndex);
  
  //mtSegmentation()
  //Kmeans segmentation in parallel processing.
  //For each mass segment, centroids are generated.
  //The information is stored in the array pointed to by m_ionEntry_p.
  void mtSegmentation(int thrIndex);
  
  //rawToGaussians
  //gets the intensity peak and converts them into Gaussians.
  //The intensity and mass data adjusted to the range of interest are loaded from the imzML file.
  //the SNR info is established for each point of the spectra.
  //The generated information is stored in the m_spectro_p structure.
  int rawToGaussians();
  
  //getMassRanges()
  //establishes the mass segments where overlapping Gaussians exist.
  //generates a Boolean mass axis with a resolution of 1/4 of the spectrometer's mass resolution.
  //if any Gaussian invades its space, sets the corresponding +/-3*sigma checkboxes to true.
  //Very wide Gaussians (sigma > 5*deltaMass) are discarded.
  //the information is stored in the array pointed to by m_massRange_p.
  //returns the number of segments.
  int getMassRanges();
  
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
  List massRangeToCentroids(MASS_RANGE massRange);
  
  GAUSS_SP *getGaussiansPointer();
  int getPixelsNumber();
  NumericMatrix getPixelGaussians(int px, float mzLow, float mzHigh);
    
private:  
  //getGaussians()
  //Called from a thread.
  //Sets the Gaussians on the peak.
  //Uses the peak separation information (m_peakFG_p).
  //px: pixel.
  //spectro: pointer to the spectrum.
  //gaussians_p: pointer to the structure containing the Gaussians' parameters.
  //Returns the number of Gaussians or a value < 0 on failure.
  int getGaussians(int px, SPECTRO *spectro_p, GAUSS_PARAMS *gaussians_p);
  
  //getCentroidsIntoRange()
  //Extracts the existing Gaussians within a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract the Gaussians.
  //gaussians_p: Requested Gaussians.
  //Returns the number of Gaussians.
  int getCentroidsIntoRange(MASS_RANGE massRange, float **gaussians_p, int size);
  
  //getRawInfo()
  //Loads the full spectrum information associated with a pixel from an imzML file.
  //The information is stored in the m_spectro structure, set to the range [m_mzLow, m_mzHigh].
  //px: Pixel whose spectrum should be loaded.
  //spIndex: Threads that manage it
  //Returns the size of the spectrum.
  int getRawInfo(int px, int spIndex);
  
  //setGaussiansIntoSegments()
  //Determines the maximum number of Gaussians over the given mass intervals and all pixels.
  //massRange: Mass range to consider.
  //Returns the maximum value.
  int setGaussiansIntoSegments(MASS_RANGE massRange); 
  
  //getCentroidsNumberIntoRange()
  //Returns the number of Gaussians in a mass range from the information in m_gaussians_p.
  //massRange: Mass range from which to extract Gaussians.
  //Returns the number of Gaussians.
  int getCentroidsNumberIntoRange(MASS_RANGE massRange);
  
public:   
  int     
    m_maxPxGaussians,
    m_massRangeSize;
  
private:  
  //input info to the constructor.
  bool m_continuous;
  int 
    m_nThreads,
    m_pxSupport,
    m_NPixels,
    m_maxMzLength;
  float     
    m_mzLow,
    m_mzHigh,
    m_SNR;
  double     
    m_massResolution,
    m_maxMassResolution;  
  
  //info generated in the class.
  GetImzMLData  *m_getImzMLData_p;
  PEAK_F_GROUP  *m_peakFG_p=0;
  GAUSS_SP      *m_gaussians_p;

  MASS_SEGMENT  m_massSegment;
  ION_ENTRY     **m_ionEntry_p;
  SPECTRO  m_spectro[MAX_THREADS];
  MASS_RANGE    *m_massRange_p;
  bool          m_enable;
  int     
                m_totalIons[MAX_THREADS],
                m_SNRmethod;
  NoiseEstimation *m_noiseEst_p; 
}; 

#endif
