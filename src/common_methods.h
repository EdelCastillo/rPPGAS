/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2019 Pere Rafols Soler
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


#ifndef RMSI_COMMON_METHODS_H
#define RMSI_COMMON_METHODS_H

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <string>
#include "types.h"

using namespace std; 

  
//Parse the UUID from a XML file to get just the hex representation in a string (without dashes and {})
std::string parse_xml_uuid(std::string uuid);

class Common
{
public:
  Common();
  bool TICnormalization(float *vect, int size);
  int sortUp  (double *bufferIn, int *sort, int size);
  int sortDown(double *bufferIn, int *sort, int size);
  int sortUpF  (float *bufferIn, int *sort, int size);
  int sortDownF(float* bufferIn, int *sort, int size);
  int nearestIndex(float value, float *data, int size);
  int nearestIndexDown(float value, float *data, int size);
  int nearestIndex(double value, double *data, int size);
  int nearestIndexGaussians(float value, GAUSS_PARAMS *data, int size);
  int nearestIndexMassRangeLow(float value, MASS_RANGE *data, int size);
  int nearestIndexMassRangeHigh(float value, MASS_RANGE *data, int size);
  float median(double *data, int size);
  float medianF(float *data, int size);
  float meanF(float *A, int size);
  double mean(double *A, int size);
  float meanF(float *prob, float *data, int size);
  double mean(double *prob, double *data, int size);
  float varF(float *A, int size, float *mean=NULL);
  double var(double *A, int size, double *mean=NULL);
  double var(double *prob, double *data, int size, double *mean_p=NULL);
  float varF(float *prob, float *data, int size, float *mean_p=NULL);
};
  
#endif