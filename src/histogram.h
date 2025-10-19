/*************************************************************************
 *     Histogram
 *     Copyright (C) 2023 Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
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
 **************************************************************************/
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <cmath>
#include <stdio.h>
#include "common_methods.h"

class Histogram
{
public:
  typedef struct
  {
    int     *elements;
    double  *values,
            *meanValues;
  }BIN;
  
  //constructor: 
  //Receives an array of denormalized data and the desired number of bins.
  //nBins = number of parts into which the dynamic range of the passed data is divided.
  // 'data' is not modified.
  Histogram(double *data, int dataSize, int nBins);
  Histogram(float  *data, int dataSize, int nBins);
  
  //destructor
  ~Histogram();
  
  //sets the bins
  int set_ranges_uniform();
  
  //updates the event counter in the histogram based on the values passed to the constructor.
  int setValue(double value);
  int setValues();
  
  //returns the number of elements in a given bin.
  int getBinElements(int bin);
  
  //returns the probability associated with a given bin.
  double getBinProbability(int bin);
  
  //returns the normalized magnitude of each bin (probability [0,1]).
  void   getHistogramYaxis(double *yAxis);
  void   getHistogramYaxis(float  *yAxis);
  
  //Returns the X-axis of the histogram with values centered over each bin.
  void   getHistogramXaxis( double *xAxis);
  void   getHistogramXaxis( float  *xAxis);
  
  //returns the middle value of the bin.
  double getBinCentralValue(int bin);
  
  //Returns the average value of the data in a bin. If there is no data, it returns zero.
  double getBinMeanValue(int bin);
  
  //returns the width of a bin: (m_maxData-m_minData)/nBins.
  double getBinValue();
  
  //returns the probability that the variable x acquires the given value or a smaller one.
  double getCDFHistogram_P(double x);
  
  //returns the probability that the variable x acquires the given value or a greater one.
  double getCDFHistogram_Q(double x);
  
  //returns the value of the random variable that makes its left-cumulative probability P .
  double getCDFHistogram_PInv(double P);
  
  //returns the value of the random variable that makes its right-side cumulative probability Q .
  double getCDFHistogram_QInv(double Q);
  
  //maximum value of the input data.
  double maxData();
  
  //minimum value of the input data.
  double minData();
  
  //mean value of the input data.
  double getDataMean();
  
  //returns the median of m_data.
  double getDataMedian();
  
  //variance of the input data.
  double getDataVariance();
  
  //variance of the input data with known mean value.
  double getDataVariance(double mean);
  
  //returns the maximum probability achieved.
  double getMaxProbability();
  
  //returns the maximum probability reached (y) and its corresponding value (x=central value of the bin).
  void   getMaxProbability(double *x, double *y);
  
  //Returns the value of m_data at the given percentile.
  //If percentile <= 0 => returns the lowest value of m_data.
  //If percentile >= 100 => returns the highest value of m_data.
  //If percentile == 50 => returns the median.
  double getQuantile(double percentile);
  
private:
  //the histogram is created from the past data.
  int    setHistogram();
  
public:  
  bool	  m_error;  //always false
  int   	
          m_dataSize,
          m_nBins,
          m_maxProbabilityBin,
          m_minProbabilityBin;
  double   
          *m_data,
  	      m_maxData,
          m_minData,
          m_binWidth,
          m_maxProbability,
          m_minProbability;
  BIN     m_bins;
};

#endif
