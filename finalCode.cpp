#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <list>
#include <string>
#include "block.h"
#include "complex.h"
#define PI 3.14159265
#define Flow 40
#define Fhigh 15000

using namespace std;

//Little endian means to store data with least significance first
//Big endian means to store data with most significance first
//this function actually reverses the bit order from LSB to MSB
//this is for unsigned output
unsigned int UnsignedLittleToBigEndianConvert( char input[], int start, int length){
		 
    unsigned int answer = 0;
  	unsigned int help = 0;
    unsigned int help2 = 1;
    
    for (int n=0; n<length; n++) {
        help = (int)(unsigned char)input[start+n];
 	  	answer += (help*help2);
        help2 *= 256;   
    }
        	 
    return answer;
}

// this is for signed output
//we need this function as data can have positive as well as negative amplitude from sound wave
int SignedLittleToBigEndianConvert( char input[], int start, int length){
	
	int answer = 0;
    int help = 0;
    int help2 = 1;
    
    for (int n=0; n<length-1; n++){
        help = (int)(unsigned char)input[start+n];
  	 	answer += (help*help2);
     	help2 *= 256;
	}
	answer += (int)input[start+length-1]*help2;
		
	return answer;
}

int main(){
	string fname;
	cout << "Enter the path of recorded file: ";
	cin >> fname; 
    FILE *fp; 

    fp = fopen(fname.c_str(),"rb"); 
    //a wave file has a structure to recognise whether it actually is a wave file.
    //first "RIFF" tag is written at the top
    //then size of the next part
    //then it has "WAVE" written
    //after that it has descriptions about format of the file
    // (its not important to our code but we need to remove it as sound data is at the end)
    //then it checks for the number of channels
    //then comes the sample rate
    //then avg bytes per sec followed by block alignment, bits per sample
    //after that comes the size of the data
    //finally the little endian data
    if (fp) 
    { 
        char id[4], *sound_buffer; //four bytes to hold 'RIFF' 
        int size; //32 bit value to hold file size 
        size_t result;
        short format_tag, channels, block_align, bits_per_sample,j; //our 16 values 
        int format_length, sample_rate, avg_bytes_sec, data_size, i; //our 32 bit values 

        fread(id, sizeof(char), 4, fp); //read in first four bytes 
        //if (strcmp(id, "RIFF")) 
        //{ 
            fread(&size, sizeof(int), 1, fp); //read in size value 
            fread(id, sizeof(char), 4, fp); //read in 4 byte string now 
            //if (strcmp(id,"WAVE")) 
            //{ 
                fread(id, sizeof(char), 4, fp);  
                fread(&format_length, sizeof(int),1,fp); 
                fread(&format_tag, sizeof(short), 1, fp); 
                fread(&channels, sizeof(short),1,fp); //1 mono, 2 stereo 
                //cout << channels << endl;
                fread(&sample_rate, sizeof(int), 1, fp); //like 44100, 22050,etc
                //cout << sample_rate<<endl;
                fread(&avg_bytes_sec, sizeof(int), 1, fp); 
                //cout << avg_bytes_sec << endl;
                fread(&block_align, sizeof(short), 1, fp); 
                //cout << block_align << endl;
                fread(&bits_per_sample, sizeof(short), 1, fp); //8 bit or 16 bit file 
                //cout << bits_per_sample << endl;
                fread(id, sizeof(char), 4, fp); //read in 'data' 
                fread(&data_size, sizeof(int), 1, fp); //how many bytes of sound data we have
                //cout<<"Data Size: " << data_size<<endl;
                sound_buffer = (char *) malloc (sizeof(char) * data_size); //set aside sound buffer space 
                fread(sound_buffer, sizeof(char), data_size, fp); //read in our whole sound data chunk 
                int t = bits_per_sample/8;
                //if  mono input
                if(channels == 1) {
                    float sound_buf[data_size/t];
                    complex a;
                    //size reduces because if bits per sample more than 8 than data actually needs to be combined
                    for(int i=0;i<data_size;i+=t) 
                    	sound_buf[i/t]=SignedLittleToBigEndianConvert(sound_buffer,i,t)/32727.0;//endian conversion
                    
                    vector<block> blocks(86);//dividing the whole data into 86 blocks on which independently
                    vector<bool> validBlocks(86, true);//fast fourier transform will be applied

                    float energySum = 0;
                    int zeroCrossingsSum  = 0;
                    int blockSize = 512;//number of data elements in the block
                    int soundBufSize = data_size/t;

                    int k = 0;
                    //zero crossing means the numebr of times the data changes sign
                    //this is important in analysing data
                    for(int i = 0; i + 512 <= data_size/t; i+=512){
                        block a(sound_buf, i, i + 512);
                    	blocks[k] = a;
                        zeroCrossingsSum += blocks[k].getZeroCrossings();
                        energySum += blocks[k].getEnergy();
                        k++;
                    }
                    float energyThreshold = energySum;
                    energyThreshold /= 86;
                    energyThreshold *= 0.5;//this is done to reduce the noise
                    //below the threshold, the data is neglected as it is noise
                    float zeroCrossingsThreshold = zeroCrossingsSum;
                    zeroCrossingsThreshold /= 86;
                    zeroCrossingsThreshold *= 0.5;
                    

                    //noise is identified by low energy and low frequency 
                    for(int i = 0; i < 86; i++){
                        if(blocks[i].getEnergy() < energyThreshold && blocks[i].getZeroCrossings() < zeroCrossingsThreshold){
                            validBlocks[i] = false;
                        } 
                    }
                    int j=0;
                    list<block> fftBlocks;
                    //merge blocks which have neither low energy nor low frequency for fft
                    while(j < 85){
                    	if(validBlocks[j] && validBlocks[j+1]){
                    		block b = blocks[j].mergeBlock(blocks[j+1]);
                    		fftBlocks.push_back(b);
                    	}
                    	j++;
                    }
                    vector<block> fftBlockVector(fftBlocks.size());
                    list<block>::iterator it;
                    j = 0;
                    for(it = fftBlocks.begin(); it != fftBlocks.end(); it++){
                        fftBlockVector[j] = *it;   
                        j++;                     
                    }  
                    //fft and peak finding
                    double rmean=0.0;int count =0;
                    vector<complex*> complexVec(fftBlocks.size());
                    for(j = 0; j < fftBlocks.size(); j++){
                        count = 0;
                        rmean = 0.0;
                        complexVec[j] = fftBlockVector[j].FFT_simple();//doing fft on data

                        for(int i=0; i<1024 ;i++) {
                            if(complexVec[j][i].get_real()>0.0){
                                rmean+=complexVec[j][i].get_real();
                                count++;
                                
                            }
                            
                        }
                        rmean/=count;
        
                        /*for(int i=0; i<1024 ;i++) {
                            if(abs(complexVec[j][i].get_real()) >=10*rmean && (i+1)/3 > Flow && (i+1)/3 < Fhigh)
                                //cout << complexVec[j][i].get_real() <<"\t"<< (i+1)/3<< endl;
                                
                            else
                                continue;
                        */
                        //cout << rmean << endl;
                    }

                }
                else {
                    //in stereo input even numbered data belong to one channel and the rest to the other
                    char sound_buffer1[data_size/2]; 
                    char sound_buffer2[data_size/2];
                    
                    //size reduces because if bits per sample more than 8 than data actually needs to be combined
                    for(int i = 0;i+1<data_size;i+=2) {
                        sound_buffer1[i/2] = sound_buffer[i];
                        sound_buffer2[i/2] = sound_buffer[i+1];
                    }
                    
                    float sound_buf1[data_size/(2*t)]; 
                    float sound_buf2[data_size/(2*t)];
                    
                    complex a;
                    
                    for(int i=0;i< data_size/2;i+=t) {
                        sound_buf1[i/t]=SignedLittleToBigEndianConvert(sound_buffer1,i,t)/32727.0;
                    	sound_buf2[i/t]=SignedLittleToBigEndianConvert(sound_buffer2,i,t)/32727.0;
                    	
                    }

                    int blockSize = 512;    //number of data elements in the block

                    // for channel 1
                    vector<block> blocks1(86); //dividing the whole data into 86 blocks on which independently
                    vector<bool> validBlocks1(86, true);    //fft will be applied
                    float energySum1 = 0;
                    int zeroCrossingsSum1  = 0;
                    //zero crossing means the numebr of times the data changes sign
                    //this is important in analysing data as it is directly proportional to the frequency
                    int soundBufSize1 = data_size/(2*t);

                    // for channel 2
                    vector<block> blocks2(86); //dividing the whole data into 86 blocks on which independently
                    vector<bool> validBlocks2(86, true); //fft will be applied
                    float energySum2 = 0;
                    int zeroCrossingsSum2  = 0;
                    int soundBufSize2 = data_size/(2*t);

                    int k = 0;

                    //create blocks
                    for(int i = 0; i + 512 <= data_size/(2*t); i+=512){
                        block a(sound_buf1, i, i + 512);
                        block b(sound_buf2, i, i + 512);
                        blocks1[k] = a;
                        blocks2[k] = b;
                        zeroCrossingsSum1 += blocks1[k].getZeroCrossings();
                        energySum1 += blocks1[k].getEnergy();
                        zeroCrossingsSum2 += blocks2[k].getZeroCrossings();
                        energySum2 += blocks2[k].getEnergy();
                        k++;
                    }


                    float energyThreshold1 = energySum1;
                    energyThreshold1 /= 86;
                    energyThreshold1 *= 0.5;
                    //this is done to reduce the noise
                    //below the threshold, the data is neglected as noise
                    float zeroCrossingsThreshold1 = zeroCrossingsSum1;
                    zeroCrossingsThreshold1 /= 86;
                    zeroCrossingsThreshold1 *= 0.5;

                    float energyThreshold2 = energySum2;
                    energyThreshold2 /= 86;
                    energyThreshold2 *= 0.5;
                    float zeroCrossingsThreshold2 = zeroCrossingsSum2;
                    zeroCrossingsThreshold2 /= 86;
                    zeroCrossingsThreshold2 *= 0.5;
                    
                    //noise is identified by low energy and low frequency 
                    for(int i = 0; i < 86; i++){
                        if(blocks1[i].getEnergy() < energyThreshold1 && blocks1[i].getZeroCrossings() < zeroCrossingsThreshold1){
                            validBlocks1[i] = false;
                        }
                        if(blocks2[i].getEnergy() < energyThreshold2 && blocks2[i].getZeroCrossings() < zeroCrossingsThreshold2){
                            validBlocks2[i] = false;
                        }
                    }

                    int j = 0;

                    list<block> fftBlocks1;
                    list<block> fftBlocks2;
                    //merge blocks which have neither low energy nor low frequency for fft
                    while(j < 85){
                        if(validBlocks1[j] && validBlocks1[j+1]){
                            block a = blocks1[j].mergeBlock(blocks1[j+1]);
                            fftBlocks1.push_back(a);
                        }
                        if(validBlocks2[j] && validBlocks2[j+1]){
                            block b = blocks2[j].mergeBlock(blocks2[j+1]);
                            fftBlocks2.push_back(b);
                        }
                        j++;
                    }

                    vector<block> fftBlockVector1(fftBlocks1.size());
                    list<block>::iterator it;
                    j = 0;
                    for(it = fftBlocks1.begin(); it != fftBlocks1.end(); it++){
                        fftBlockVector1[j] = *it;   
                        j++;                     
                    }

                    vector<block> fftBlockVector2(fftBlocks2.size());
                    j = 0;
                    for(it = fftBlocks2.begin(); it != fftBlocks2.end(); it++){
                        fftBlockVector2[j] = *it;   
                        j++;                     
                    }


                    double rmean1=0.0;int count1 =0;
                    double rmean2=0.0;int count2 =0;
                    vector<complex*> complexVec1(fftBlocks1.size());
                    vector<complex*> complexVec2(fftBlocks2.size());
                    //cout << "Channel 1: "<<endl;
                    //apply fft and find peaks for channel 1
                    for(j = 0; j < fftBlocks1.size(); j++){
                        count1 = 0;
                        rmean1 = 0.0;
                        complexVec1[j] = fftBlockVector1[j].FFT_simple();

                        for(int i=0; i<1024 ;i++) {
                            if(complexVec1[j][i].get_real()>0.0){
                                rmean1+=complexVec1[j][i].get_real();
                                count1++;
                               
                            }
                            
                        }
                        rmean1/=count1;
        
                        /*for(int i=0; i<1024 ;i++) {
                            if(abs(complexVec1[j][i].get_real()) >=3*rmean1 && (i+1)/3 > Flow && (i+1)/3 < Fhigh)
                                cout << complexVec1[j][i].get_real() <<"\t"<< (i+1)/3<< endl;
                                //cout << "checkpoint 8" << endl;
                            else
                                continue;
                        }*/
                        //cout << rmean1 << endl;
                    }

                    //apply fft and find peaks for channel 1
                    //cout << "Channel 2: "<<endl;
                    for(j = 0; j < fftBlocks2.size(); j++){
                        count2 = 0;
                        rmean2 = 0.0;
                        complexVec2[j] = fftBlockVector2[j].FFT_simple();

                        //cout << "checkpoint 6" << endl;
                        for(int i=0; i<1024 ;i++) {
                            if(complexVec2[j][i].get_real()>0.0){
                                rmean2+=complexVec2[j][i].get_real();
                                count2++;
                                
                            }
                            
                        }
                        rmean2/=count2;

                    }
                    int zeroCrossingsSumC1 = 0;
                    int zeroCrossingsSumC2 = 0;
                    float energyC1 = 0.0;
                    float energyC2 = 0.0;
                    int countC1 = 0;
                    int countC2 = 0;
                    
                    //cout << "Peaks Channel2: "<< endl;
                    for(int i = 0; i < fftBlockVector2.size(); i++){
                        fftBlockVector2[i].getVectorForPeaks();
                        fftBlockVector2[i].findPeaks();
                        list<int> peaks = fftBlockVector2[i].getPeaks();
                        list<int>::iterator it;
                        vector<float> vec = fftBlockVector2[i].getVec();
                        /*cout << "block index: " << i << endl;
                        for(it = peaks.begin(); it != peaks.end(); it++){
                            cout << vec[*it]<< "  " << ((*it) * sample_rate) / 2048.0 << "\t";
                        }
                        cout << endl;*/
                        
                    }
                    //calculate the charasteristics of the sound
                    int countPeaks2 = 0;int num2 = 0;
                    for(j = (int) fftBlocks2.size()*0.6; j < fftBlocks2.size(); j++){
                        energyC2 += fftBlockVector2[j].getEnergy();
                        zeroCrossingsSumC2 +=  fftBlockVector2[j].getZeroCrossings();
                        countC2++;
                    }

                    float avgCrossC2 = (float) zeroCrossingsSumC2/countC2;
                    float avgEnergyC2 = energyC2/countC2;
                    
                    zeroCrossingsSumC2 = 0; countC2 = 0; energyC2 = 0.0;

                    for(j = 0; j < fftBlocks2.size()*0.01; j++){
                        energyC2 += fftBlockVector2[j].getEnergy();
                        zeroCrossingsSumC2 +=  fftBlockVector2[j].getZeroCrossings();
                        countC2++;
                    }
                    float avgCrossI2 = (float) zeroCrossingsSumC2/countC2;  
                    float avgEnergyI2 = energyC2/countC2;  

                    zeroCrossingsSumC2 = 0; countC2 = 0; energyC2 = 0.0;

                    float totalEnergy = 0.0;
                    for(j = 0; j < fftBlocks2.size(); j++){
                        totalEnergy += fftBlockVector2[j].getEnergy();
                    }
                    

                    for(j = fftBlocks2.size()*0.3; j < fftBlocks2.size()*0.5; j++){
                        energyC2 += fftBlockVector2[j].getEnergy();
                        zeroCrossingsSumC2 +=  fftBlockVector2[j].getZeroCrossings();
                        countC2++;
                    }

                    float avgCrossM2 = (float) zeroCrossingsSumC2/countC2;  
                    float avgEnergyM2 = energyC2/countC2;  
                    float fracEnergyM2 = energyC2/totalEnergy;

                    //cout << "Channel2 0.01 : " << avgEnergyI2 << "\t"  << avgCrossI2  <<  endl;
                    //cout << "Channel2 0.2 : " << avgEnergyM2 << "\t"  << avgCrossM2 << "\t" << fracEnergyM2 << endl;
                    //cout << "Channel2 0.6 : " << avgEnergyC2 << "\t"  << avgCrossC2 << endl;

                    //matching of enrgies and frequencies for variety of phonemes,
                    //according to that output the word spoken
                    if((avgEnergyI2 < 1350000.0 && avgEnergyI2 > 1200000.0) || (avgCrossI2 < 55.0 && avgCrossI2 > 50.0) || avgEnergyI2 < 3000.0)
                        cout << "close" << endl;   
                    else if((avgCrossC2 < 50.0 && avgCrossC2 > 40.0) || (avgEnergyC2 < 4000.0 && avgEnergyC2 > 3400.0) || ((avgCrossC2-51.0) < 1.0 && (avgCrossC2-51.0) > 0.0))
                        cout << "save" << endl;
                    else if(avgCrossC2 < 70.0 || avgEnergyC2  < 5000.0)
                        cout << "no" << endl;
                    else
                        cout << "yes" << endl;
                    
//------------------------------------------------------------
                }
                    
            //} 
            //else 
              //  printf("Error: RIFF file but not a wave file\n"); 
        //} 
        //else 
          //  printf("Error: not a RIFF file\n"); 
    } 
    return 0;
} 
