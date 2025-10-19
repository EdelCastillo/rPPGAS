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

#include "common_methods.h"

struct
{
  bool operator()(int a, int b) const { return a > b; }
}
cmp_Down;
struct
{
  bool operator()(int a, int b) const { return a < b; }
}
cmp_Up;


std::string parse_xml_uuid(std::string uuid)
{
  std::size_t ipos = uuid.find('{');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  do
  {
    ipos = uuid.find('-'); 
    if( ipos != std::string::npos )
    {
      uuid.erase(ipos, 1);
    } 
  } while ( ipos != std::string::npos);
  ipos = uuid.find('}');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  for( unsigned int i=0; i < uuid.length(); i++)
  {
    uuid[i] = toupper(uuid[i]);
  }
  return uuid;
}

Common::Common(){}

//Normalización por TIC
bool Common::TICnormalization(float *vect, int size)
{
  double suma=0;
  for(int i=0; i<size; i++)
    suma+=vect[i];
  suma/=1e6;
  if(suma<1e-10) return(false); //errors control
  for(int i=0; i<size; i++)
    vect[i]/=suma;
  return(true);
}

int Common::sortUp(double *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_Up /*[&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];}*/);
  return 0;
}

int Common::sortDown(double *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] > bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_Down);
  return 0;
}

int Common::sortUpF(float *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] < bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_Up /*[&bufferIn](size_t i1, size_t i2) {return ((float*)bufferIn)[i1] < ((float*)bufferIn)[i2];}*/);
  return 0;
}

int Common::sortDownF(float *bufferIn, int *index, int size) 
{
  if(index)
  {
    for(int i=0; i<size; i++) index[i]=i;
    stable_sort(index, index+size, [&bufferIn](size_t i1, size_t i2) {return bufferIn[i1] > bufferIn[i2];});
  }
  else
    stable_sort(bufferIn, bufferIn+size, cmp_Down);
  return 0;
}


//retorna el índice a data más próximo a value
int Common::nearestIndex(float value, float *data, int size)
{
//  printf("[%.1f]->", value); for(int i=0; i<size; i++) printf("%.1f ", data[i]); printf("\n");    
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //retorna el más próximo
  {
    if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow])       {/*printf(".%d\n", indexLow); */return indexLow;}
  else if(value>=data[indexHigh]) {/*printf("..%d\n", indexHigh);*/return indexHigh;}
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter]) {/*printf("...%d\n", indexCenter);*/return indexCenter;}
    if(value<data[indexCenter]) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) {/*printf("....%d\n", indexLow);*/return indexLow;}
      else return indexHigh;
    }
  }
}

//retorna el índice a data más próximo a value
int Common::nearestIndex(double value, double *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //retorna el más próximo
  {
    if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow])       return indexLow;
  else if(value>=data[indexHigh]) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2.0);
    if(value==data[indexCenter]) return indexCenter;
    if(value<data[indexCenter]) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) return indexLow;
      else return indexHigh;
    }
  }
}

//retorna el índice a data más próximo a value
int Common::nearestIndexGaussians(float value, GAUSS_PARAMS *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //retorna el más próximo
  {
    if(value-data[indexLow].mean <= data[indexHigh].mean-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow].mean)       return indexLow;
  else if(value>=data[indexHigh].mean) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2);
    if(value==data[indexCenter].mean) return indexCenter;
    if(value<data[indexCenter].mean) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow].mean <= data[indexHigh].mean-value) return indexLow;
      else return indexHigh;
    }
  }
}

//retorna el índice a data.low más próximo a value
int Common::nearestIndexMassRangeLow(float value, MASS_RANGE *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //retorna el más próximo
  {
    if(value-data[indexLow].low <= data[indexHigh].low-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow].low)       return indexLow;
  else if(value>=data[indexHigh].low) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2);
    if(value==data[indexCenter].low) return indexCenter;
    if(value<data[indexCenter].low) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow].low <= data[indexHigh].low-value) return indexLow;
      else return indexHigh;
    }
  }
}

//retorna el índice a data.high más próximo a value
int Common::nearestIndexMassRangeHigh(float value, MASS_RANGE *data, int size)
{
  int indexLow=0;
  int indexHigh=size-1;
  int indexCenter;
  
  if(indexHigh==indexLow) return 0;
  if(indexHigh==indexLow+1) //retorna el más próximo
  {
    if(value-data[indexLow].high <= data[indexHigh].high-value) return indexLow;
    else return indexHigh;
  }
  
  if(value<=data[indexLow].high)       return indexLow;
  else if(value>=data[indexHigh].high) return indexHigh;
  
  while(1)
  {
    indexCenter=round((indexHigh+indexLow)/2);
    if(value==data[indexCenter].high) return indexCenter;
    if(value<data[indexCenter].high) indexHigh=indexCenter; 
    else indexLow=indexCenter;
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow].high <= data[indexHigh].high-value) return indexLow;
      else return indexHigh;
    }
  }
}

//retorna la mediana
//si data tiene longitud par, se retorna el valor promediado entre los dos valores centrales
float Common::median(double *data, int size)
{
  if (size<=0) return -1;
  double median;
  int half=size/2; //si par, se queda con el valor inferior
  double *tmpData=0;
  tmpData=new double[size];
  for(int i=0; i<size; i++) //copia, dado que sort() destruye el origen
    tmpData[i]=data[i];
  sortUp(tmpData, 0, size); //ordenación creciente
  
  if(!(size%2)) {//si par
    if(half>0)
    {median=(tmpData[half]+tmpData[half-1])/2.0;} //valor medio
    else {median=tmpData[0];}
  }
  else
  {median=tmpData[half];}
  if(tmpData) delete []tmpData;
  return median;
}
//retorna la mediana
//si data tiene longitud par, se retorna el valor promediado entre los dos valores centrales
float Common::medianF(float *data, int size)
{
  if (size<=0) return -1;
  float median;
  int half=size/2; //si par, se queda con el valor inferior
  float *tmpData=0;
  tmpData=new float[size];
  for(int i=0; i<size; i++) //copia, dado que sort() destruye el origen
    tmpData[i]=data[i];
  sortUpF(tmpData, 0, size); //ordenación creciente
  
  if(!(size%2)) {//si par
    if(half>0)
    {median=(tmpData[half]+tmpData[half-1])/2.0;} //valor medio
    else {median=tmpData[0];}
  }
  else
    {median=tmpData[half];}
   if(tmpData) delete []tmpData;
  return median;
}

//valor medio de un array de floats
float Common::meanF(float *A, int size)
{
  float S=0;
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=A[i];
  return S/size;
}

//valor medio de un array de doubles
double Common::mean(double *A, int size)
{
  double S=0;
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=A[i];
  return S/size;
}

//valor medio de un array de doubles
double Common::mean(double *prob, double *data, int size)
{
  double S=0;   
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=prob[i]*data[i];
  return S;
}

//valor medio de un array de floats
float Common::meanF(float *prob, float *data, int size)
{
  float S=0;   
  if(size<=0) return -1;
  for(int i=0; i<size; i++)
    S+=prob[i]*data[i];
  return S;
}

//varianza de un array de floats
float Common::varF(float *A, int size, float *mean_p)
{
  float m;
  float S=0;
  if(size<=0) return -1;
  if(mean_p) m=*mean_p;
  else m=meanF(A, size);
  
  for(int i=0; i<size; i++)
    S+=(A[i]-m)*(A[i]-m);
  return S/size;
}

//varianza de un array de double
double Common::var(double *A, int size, double *mean_p)
{
  double m;
  double S=0;
  if(size<=0) return -1;
  
  if(mean_p) m=*mean_p;
  else m=mean(A, size);
  
  for(int i=0; i<size; i++)
    S+=(A[i]-m)*(A[i]-m);
  return S/size;
}

//varianza de un array de double
double Common::var(double *prob, double *data, int size, double *mean_p)
{
  double m;
  double S=0;
  if(size<=0) return -1;
  
  if(mean_p) m=*mean_p;
  else m=mean(prob, data, size);
  
  for(int i=0; i<size; i++)
    S+=(data[i]-m)*(data[i]-m)*prob[i];
  return S;
}

//varianza de un array de floats
float Common::varF(float *prob, float *data, int size, float *mean_p)
{
  float m;
  float S=0;
  if(size<=0) return -1;
  if(mean_p) m=*mean_p;
  else m=meanF(prob, data, size);
  
  for(int i=0; i<size; i++)
    S+=(data[i]-m)*(data[i]-m)*prob[i];
  return S;
}

