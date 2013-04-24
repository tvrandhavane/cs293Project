#ifndef BLOCK_H_
#define BLOCK_H_
#include <iostream>
#include <list>
#include <vector>
#include "complex.h"

using namespace std;
/*this class takes care of everything related to a block of data, 
it stores the total energy of the block, the number of zero crossings.
 It does fast fourier transform on the data and then find peaks on the output for phoneme detection*/
class block{
	vector<float> vec;

	vector<float> vectorForPeaks;//vector for storing peaks

	complex* compVector;//complex vector formed after doing fft(fast fourier transform)

	int blockSize;	//size of blocks

	list<int> peaks;	//list of peaks present in the block
	
	float totalEnergy;	//total energy of the block
	
	int zeroCrossings;	//number of zero crossingd
	
	complex* FFT_simpleH(complex* x, int N ); //fft helper
	
	public:
    
    block(){}//default constructor
    
    void getVectorForPeaks();//returns the vector storing peaks
	
	complex* FFT_simple();//does fft on the data
	
	block(float* sound, int start, int end);	//constructor for block 
	
	void convertToComplex();//make float vector to complex vector for fft
	
	//accessor functions
	float getEnergy();	
	
	int getBlockSize();
	
	int getZeroCrossings();

	list<int> getPeaks();
	
	vector<float> getVec();
	
	block mergeBlock(block a);//fft is not done on disjoint blocks but with 50% overlap, 
							  //so initially blocks are made of half the size and then one block is merged with another block
							  //1st block is merged with 2nd then fft is done, 2nd is merged with 3rd and then fft is done
	
	void findPeaks();		//functions to find peaks in the sound wave
	
	
	//overloaded equal to operator
	void operator=(const block b){
		this->vec = b.vec;
		this->blockSize = b.blockSize;
		this->totalEnergy = b.totalEnergy;
		this->zeroCrossings = b.zeroCrossings;
	}


};
#endif
