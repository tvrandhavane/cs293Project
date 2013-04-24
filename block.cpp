#include <iostream>
#include <vector>
#include "block.h"
#include "complex.h"


#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#define PI 3.14159265
#define Flow 40
#define Fhigh 15000

using namespace std;

block block::mergeBlock(block a) {
	  block res;
	  res.vec = this->vec;
	  res.vec.insert(res.vec.end(),a.vec.begin(),a.vec.end());//merging vectors
	  res.totalEnergy = this->totalEnergy + a.totalEnergy;//adding both energies
	  res.zeroCrossings = this->zeroCrossings + a.zeroCrossings;//adding zero crossings
	  res.blockSize = this->blockSize + a.blockSize;//increasing block size
	  res.convertToComplex();//forming the complex vector to do fft on
	  return res;
	}
//function for fast fourier transform
//implementation of cooley-turkey algorithm
//recursive call of fft on even  elements and odd elements      
complex* block::FFT_simpleH(complex* x, int N ) {
    complex* X = (complex*) malloc(16 * N);
    complex * d, * e, * D, * E;
    complex a;
    int k;

    if (N == 1) {
        X[0] = x[0];
        return X;
    }

    e = (complex*) malloc(8 * N);
    d = (complex*) malloc(8 * N);
    for(k = 0; k < N/2; k++) {
        e[k] = x[2*k];//even-indexed elements
        d[k] = x[2*k + 1];//odd-indexed elements
    }
    
    E = FFT_simpleH(e, N/2);//fft on even
    D = FFT_simpleH(d, N/2);//fft on odd

    free(e);
    free(d);

    for(k = 0; k < N/2; k++) {
        D[k] = a.complex_mult(a.complex_from_polar(1, -2.0*PI*k/N), D[k]);//multipying by e^i*2*pi*k/N
    }
    ///merging data to form the final output
    for(k = 0; k < N/2; k++) {
        X[k] = a.complex_add(E[k], D[k]);
        X[k + N/2] = a.complex_sub(E[k], D[k]);
    }

    free(D);
    free(E);
    return X;
}

block::block(float* sound, int start, int end){//sound is the whole vector
	blockSize = end - start;                   //block is composed of elements of sound from index start to end-1
	totalEnergy = 0;
	zeroCrossings = 0;         
    vec.resize(blockSize);
	for(int i = start; i+1 < end; i++){//computing and assigning energy and zero crossings
		vec[i - start] = sound[i];                    
		totalEnergy += sound[i]*sound[i];             
        if(sound[i]*sound[i+1] < 0.0)
            zeroCrossings++;
	}
    vec[end-start-1] = sound[end-1];               
    totalEnergy += sound[end-1]*sound[end-1];   
}

//some accessor functions
float block::getEnergy(){
	return totalEnergy;            
}

int block::getBlockSize(){
	return blockSize;
}

int block::getZeroCrossings(){
	return zeroCrossings;
}

void block::getVectorForPeaks(){
    convertToComplex();
    vectorForPeaks.resize(blockSize);
    for(int i = 0; i < blockSize; i++){
        vectorForPeaks[i] = abs(compVector[i].get_real());
    }
}

list<int> block::getPeaks(){
    return peaks;
}

vector<float> block::getVec(){
    return vectorForPeaks;
}

//converts the sound vector into a complex array
void block::convertToComplex(){
    complex* compVector1 = new complex[blockSize];
    complex a;
    for(int i=0; i<blockSize; i++){
        compVector1[i] = a.assign(vec[i],0.0);
    }
    compVector = compVector1;
}

//call for main fft function
complex* block::FFT_simple(){
    convertToComplex();
    return FFT_simpleH(compVector, blockSize);
}

//Peak Finding/Detection Algorithm
void block::findPeaks(){
    float maximum = 0;
    float stdMean = 0.0;
    
    //loop through the entire sound vector to calculate standard deviation and maximum value
    for(int i = 0; i + 1 < blockSize && vectorForPeaks.size()!=0; i++){
        stdMean += abs(vectorForPeaks[i + 1] - vectorForPeaks[i]);
        if(vectorForPeaks[i] > maximum){
            maximum = vectorForPeaks[i];
        }
    }

    stdMean /= blockSize;
    
    //set thresholds
    float stdThreshold = stdMean;
    float threshold = 0.6 * maximum;

    list<int> maximums;         //list of indices i, where abs(vec[i]) > 0.6 * maximum value
    vector<bool> maxOrNot(blockSize,false);         //a vector of bools to store peaks
    float diffThreshold = -0.00001;

    //set data in list of maximum indices
    for(int i = 1; i < blockSize && vectorForPeaks.size()!=0; i++){
        if((vectorForPeaks[i] > threshold) && ((vectorForPeaks[i] - vectorForPeaks[i-1]) > stdThreshold)){
            maximums.push_back(i);
            maxOrNot[i] = true;
        }
    }    

    list<int>::iterator it;
    int counter = 0;
    //iterate through maximums list to check for peak, 
    //there is a peak at index i, if for the 8 out of 10 times the function increases before i
    if(maximums.size() != 0){
        for(it=maximums.begin(); it != maximums.end(); it++){
            int index = *it;
            int maximum = vectorForPeaks[index];
            counter = 0;
            int i;
            if(index - 10 < 0){
                i = 1;
            }
            else{
                i = index - 10;
            }
            for(; i < index; i++){
                if(vectorForPeaks[i + 1] - vectorForPeaks[i] > diffThreshold){
                    counter++;
                    if(maxOrNot[i]){
                        maxOrNot[i] = false;
                    }
                }
            }
            
            if(counter < 8){
                maxOrNot[index] = false;
            }
        }
    }

    //peaks is a list of peaks, add peaks into it
    peaks.clear();
    for(int i = 0; i < blockSize; i++){
        if(maxOrNot[i]){
            peaks.push_back(i);
        }
    }
}

