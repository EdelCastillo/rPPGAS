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
#include "histogram.h"

//constructor: 
//Receives an array of denormalized data and the desired number of bins.
//nBins = number of parts into which the dynamic range of the passed data is divided.
// 'data' is not modified.
Histogram::Histogram(double *data, int dataSize, int nBins)
  {
  m_error=false;  
  m_data=data;
  m_dataSize=dataSize;
  m_nBins=nBins;
  m_maxData=maxData();
  m_minData=minData();
  m_binWidth=(m_maxData-m_minData)/nBins;
  if(m_binWidth==0) {m_binWidth=1e-10;}//m_error=1; return;}
  m_bins.elements=0;
  m_bins.meanValues=0;
  m_bins.values=0;
  
  setHistogram();
  }
  
  //constructor: 
  //Receives an array of denormalized data and the desired number of bins.
  //nBins = number of parts into which the dynamic range of the passed data is divided.
  // 'data' is not modified.
  Histogram::Histogram(float *data, int dataSize, int nBins)
  {
  m_error=false;  
  m_data=0;
  //the data string is converted to double
  if((m_data=new double[dataSize])==0)
    {printf("ERROR new in Histogram class\n"); m_error=1; return;}
  for(int i=0; i<dataSize; i++)
    m_data[i]=data[i];
  m_dataSize=dataSize;
  m_nBins=nBins;
  m_maxData=maxData();
  m_minData=minData();
  m_binWidth=(m_maxData-m_minData)/nBins;
  if(m_binWidth==0) {m_binWidth=1e-10;}
  
  m_bins.elements=0;
  m_bins.meanValues=0;
  m_bins.values=0;
  setHistogram();
  }
  
//destructor
Histogram::~Histogram(  )
  {
//  printf("Histogram destructor input\n");
  if(m_bins.elements) delete m_bins.elements;
  if(m_bins.values) delete m_bins.values;
  if(m_bins.meanValues) delete m_bins.meanValues;
  
  if(m_data) delete[] m_data; 
  //  printf("Histogram destructor output\n");
  }

//the histogram is created from the past data.
int Histogram::setHistogram()
  {
  if(m_minData>m_maxData)
    {
      printf("Warning: histogram range is wrong (min=%.3f max=%.3f)\n", m_minData, m_maxData);   
      m_error=true;  
      return -1;      
    }
  if(set_ranges_uniform()<0) {m_error=true; return -1;} //initializes the histogram

  if(setValues()<0) {m_error=true; return -1;} //sets all its values.
  
  double prob;
  m_maxProbability=1e-300;
  m_minProbability=1e+300;
  
  for(int i=0; i<m_nBins; i++)
    {
    prob=getBinProbability(i);
    if(prob>m_maxProbability) {m_maxProbability=prob;  m_maxProbabilityBin=i;}
    if(prob<m_minProbability) {m_minProbability=prob;  m_minProbabilityBin=i;}
    }
  return 0;
  }

//sets the bins
int Histogram::set_ranges_uniform()
{
  if(m_minData>m_maxData) {return -1;}
  if(m_nBins<1) {return -2;}
  
  m_bins.elements=new int[m_nBins];
  m_bins.values=new double[m_nBins];
  m_bins.meanValues=new double[m_nBins];
  
  double binRange=(m_maxData-m_minData)/m_nBins;
  
  for(int i=0; i<m_nBins; i++) 
  {
    m_bins.elements[i]=0;
    m_bins.meanValues[i]=0;
    m_bins.values[i]=m_minData+i*binRange;
  }
  return 0;
}

//updates an event counter value in the histogram.
int Histogram::setValue(double value)
{
  if(value<m_bins.values[0] || value>m_maxData) return -1;
  
  int index, delta, bin;
  int h=m_nBins-1; 	//high index within the search range.
  int l=0;		      //low  index within the search range.
  
  //Extreme cases require special attention.
  if(value>=m_bins.values[m_nBins-1] && value<=m_maxData) 
  {m_bins.elements[m_nBins-1]++; m_bins.meanValues[m_nBins-1]+=value; return 0;}		
  else if(value>=m_minData && value<m_bins.values[1]) 
  {m_bins.elements[0]++; m_bins.meanValues[0]+=value; return 0;};
  
  index=m_nBins/2;
  
  if(!index) return -1; //only one element and is not in range.
  
  while(1)
  {
    if(value>=m_bins.values[index] && value<m_bins.values[index+1])
    {m_bins.elements[index]++; m_bins.meanValues[index]+=value; return 0;}
    else if(value>=m_bins.values[index]) 
    {
      l=index;
      if((delta=(h-l)/2)==0) break; //does not exist
      index=delta+l;
    }
    else 
    {
      h=index;
      if((delta=(h-l)/2)==0) break; //does not exist
      index=delta+l;
    }
  }
  return -1;
}

//updates the event counter in the histogram based on the values passed to the constructor.
int Histogram::setValues()
{
  for(int i=0; i<m_dataSize; i++)
    if(setValue(m_data[i])<0) return -1;
    return 0;
}

//returns the number of elements in a given bin.
int Histogram::getBinElements(int bin)
{
  if(bin<0 || bin>=m_nBins) return -1;
  return m_bins.elements[bin];
}

//returns the probability associated with a given bin.
double Histogram::getBinProbability(int bin)
{
  int size;
  if((size=getBinElements(bin))<0) return -1;
  return((double)size/m_dataSize);
}

//returns the normalized magnitude of each bin (probability [0,1]).
void Histogram::getHistogramYaxis(double *yAxis)
{
  if(m_error) return;
  for(int i=0; i<m_nBins; i++)
    yAxis[i]=getBinProbability(i);
}

//returns the normalized magnitude of each bin (probability [0,1]).
void Histogram::getHistogramYaxis(float *yAxis)
{
  if(m_error) return;
  for(int i=0; i<m_nBins; i++)
    yAxis[i]=getBinProbability(i);
}

//Returns the X-axis of the histogram with values centered over each bin.
void Histogram::getHistogramXaxis( double *xAxis)
  {
  for(int bin=0; bin<m_nBins; bin++)
    xAxis[bin]=m_minData+bin*m_binWidth+m_binWidth/2.0;
  }
  
//Returns the X-axis of the histogram with values centered over each bin.
void Histogram::getHistogramXaxis( float *xAxis)
  {
  for(int bin=0; bin<m_nBins; bin++)
    xAxis[bin]=m_minData+bin*m_binWidth+m_binWidth/2.0;
  }
  
//returns the middle value of the bin.
double Histogram::getBinCentralValue(int bin)
  {
  if(bin<0 || bin>=m_nBins) return -1;
  return m_minData+bin*m_binWidth+m_binWidth/2.0; 
  }

//Returns the average value of the data in a bin. If there is no data, it returns zero.
double Histogram::getBinMeanValue(int bin)
  {
    if(bin<0 || bin>=m_nBins) return -1;
    if(m_bins.elements[bin]==0) return 0;
    return m_bins.meanValues[bin]/m_bins.elements[bin];
  }
 
//returns the width of a bin: (m_maxData-m_minData)/nBins.
double Histogram::getBinValue()
  {
  return m_binWidth; 
  }

//returns the probability that the variable x acquires the given value or a smaller one.
double Histogram::getCDFHistogram_P(double x)
  {
  double value, add=0;
  int bin;
  if(m_error) return -1;
  if(x<m_minData || x>m_maxData) return -1; 
  bin=(x-m_minData)/m_binWidth;
  if(x==m_maxData) bin=m_nBins-1;
  
  for(int i=0; i<=bin; i++)
    add+=getBinProbability(i); 
  return add;
  }

//returns the probability that the variable x acquires the given value or a greater one.
double Histogram::getCDFHistogram_Q(double x)
{
  double value, add=0;
  int bin;
  if(m_error) return -1;
  if(x<m_minData || x>m_maxData) return -1;
  bin=(x-m_minData)/m_binWidth;
  if(x==m_minData) bin=0;
  
  for(int i=m_nBins-1; i>bin; i--)
    add+=getBinProbability(i); 
  return add;
}

//returns the value of the random variable that makes its left-cumulative probability P .
double Histogram::getCDFHistogram_PInv(double P)
  {
  double add=0;
  int bin;
  
  if(m_error) return -1;
  if(P<0 || P>1) return -1;
   
  for(bin=0; bin<m_nBins; bin++)
    if((add+=getBinProbability(bin))>P) break;
  if(bin==m_nBins) 
    return m_maxData;;
  return m_minData+m_binWidth*((double)bin+0.5);
  }

//returns the value of the random variable that makes its right-side cumulative probability Q .
double Histogram::getCDFHistogram_QInv(double Q)
  {
  double add=0;
  int bin;
  
  if(m_error) return -1;
  if(Q<0 || Q>1) return -1;
   
  for(bin=m_nBins-1; bin>=0; bin--)
    if((add+=getBinProbability(bin))>Q) break;
  if(bin==0) 
    return m_minData;
  return m_minData+m_binWidth*((double)bin+0.5);
  }

//maximum value of the input data.
double Histogram::maxData()
{
  double maxData=1e-300;
  for(int i=0; i< m_dataSize; i++)
    if(m_data[i]>maxData) 
      maxData=m_data[i];
    return maxData;
}

//minimum value of the input data.
double Histogram::minData()
{
  double minData=1e+300;
  for(int i=0; i< m_dataSize; i++)
    if(m_data[i]<minData) 
      minData=m_data[i];
    return minData;
}

//mean value of the input data.
double Histogram::getDataMean()
{
  double mean=0;
  if(m_dataSize<=0) return 0;
  for(int i=0; i< m_dataSize; i++)
    mean+=m_data[i];
  return mean/m_dataSize;
}

//returns the median of m_data.
double Histogram::getDataMedian()
{
  Common tools;
  return tools.median(m_data, m_dataSize);
}

//variance of the input data.
double Histogram::getDataVariance()
{
  double mean=0, var=0;
  if(m_dataSize<=1) return 0;
  mean= getDataMean(); 
  for(int i=0; i< m_dataSize; i++)
    var+=(m_data[i]-mean)*(m_data[i]-mean);
  var/=(m_dataSize-1);
  return var;
}

//variance of the input data with known mean value.
double Histogram::getDataVariance(double mean)
{
  double var=0;
  if(m_dataSize<=1) return 0;
  for(int i=0; i< m_dataSize; i++)
    var+=(m_data[i]-mean)*(m_data[i]-mean);
  var/=(m_dataSize-1);
  return var;
}

//returns the maximum probability achieved.
double Histogram::getMaxProbability()
{
  return m_maxProbability;
}

//returns the maximum probability reached (y) and its corresponding value (x=central value of the bin).
void Histogram::getMaxProbability(double *x, double *y)
{
  *y=m_maxProbability;
  *x=m_binWidth*((double)m_maxProbabilityBin+0.5);
}

//Returns the value of m_data at the given percentile.
//If percentile <= 0 => returns the lowest value of m_data.
//If percentile >= 100 => returns the highest value of m_data.
//If percentile == 50 => returns the median.
double Histogram::getQuantile(double percentile)
  {
  Common tools;
  if(percentile<=0) return m_minData;
  else if(percentile>100) return m_maxData;
  
  double *tmpData=0;
  tmpData=new double[m_dataSize];
  
  for(int i=0; i<m_dataSize; i++) //It is required, since the original is destroyed in sort().
    tmpData[i]=m_data[i];
  tools.sortUp(tmpData, 0, m_dataSize);
  
  int index=round(m_dataSize*percentile/100.0);
  double p=tmpData[index];
  if(tmpData) delete []tmpData;
  
  return p;
  }


