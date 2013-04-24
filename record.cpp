#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "portaudio.h"
#include <sndfile.h>

using namespace std;

#ifndef PI
#define PI 3.14159265358979323846264338
#endif

#define SAMPLE_RATE (44100)
#define SAMPLE_COUNT (SAMPLE_RATE * 1)   // 1 second
#define AMPLITUDE (1.0 * 0x7F000000)
#define LEFT_FREQ (344.0 / SAMPLE_RATE)
#define RIGHT_FREQ (466.0 / SAMPLE_RATE)

#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS     (1)
#define NUM_CHANNELS    (2)
#define DITHER_FLAG     (0) 

#define WRITE_TO_FILE   (0)
// Select sample format.
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#endif

typedef struct
{
    int          frameIndex;  // Index into sample array. 
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}
paTestData;

/* This function is called by the PortAudio library when audio is needed.
 It may be called at interrupt level on some machines so don't do anything
 that could mess up the system like calling malloc() or free().
*/
static int recordCallback( const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData )
{
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;
 
    (void) outputBuffer; 
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
 
    if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }
 
    if( inputBuffer == NULL )
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = SAMPLE_SILENCE;  //left
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  //right
        }
    }
    else
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = *rptr++;  //left
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  //right
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

/* This routine will be called by the PortAudio engine when audio is needed.
 It may be called at interrupt level on some machines so don't do anything
 that could mess up the system like calling malloc() or free().
*/
static int playCallback( const void *inputBuffer, void *outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void *userData )
{
    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;
 
    (void) inputBuffer; // Prevent unused variable warnings.
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
 
    if( framesLeft < framesPerBuffer )
    {
        //final buffer... 
        for( i=0; i<framesLeft; i++ )
        {
	        *wptr++ = *rptr++;  //left 
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  // right
        }
	    for( ; i<framesPerBuffer; i++ )
        {
            *wptr++ = 0;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = 0;  // right 
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else
    {
        for( i=0; i<framesPerBuffer; i++ )
        {
            *wptr++ = *rptr++;  // left
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  // right 
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}
 
/*******************************************************************/
int main(){
    SNDFILE            *file ;          //wav file to write
    SF_INFO             sfinfo ;        //sound info file parameters
    PaStreamParameters  inputParameters,
                        outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;
 
    fflush(stdout);
 

    //set sound file parameters
    sfinfo.samplerate   = SAMPLE_RATE ;
    sfinfo.frames       = SAMPLE_COUNT ;
    sfinfo.channels     = 2 ;
    sfinfo.format       = (SF_FORMAT_WAV | SF_FORMAT_PCM_24) ;

    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; // Record for a few seconds. 
    //Set to 1 because 1 word can easily be spoken in 1 second
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;                //total number of samples
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = new SAMPLE[numBytes]; //initialise recordedSamples
    for( i=0; i<numSamples; i++ ){
        data.recordedSamples[i] = 0;
    }
 
    err = Pa_Initialize();          //initializes portaudio
    if( err != paNoError ){         //if initialization error exit the code, else continue
    }
    else{
        inputParameters.device = Pa_GetDefaultInputDevice(); // default input device
        if (inputParameters.device == paNoDevice) {
            fprintf(stderr,"Error: No default input device.\n");
        }
        else{
            inputParameters.channelCount = 2;                    //stereo input
            inputParameters.sampleFormat = PA_SAMPLE_TYPE;
            inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
            inputParameters.hostApiSpecificStreamInfo = NULL;
            // Record some audio.
            err = Pa_OpenStream(&stream,
                                &inputParameters,
                                NULL,                  
                                SAMPLE_RATE,
                                FRAMES_PER_BUFFER,
                                paClipOff,      
                                recordCallback,
                                &data);
            if( err != paNoError ){
            }
            else{
                err = Pa_StartStream( stream );
                if( err != paNoError ){
                }
                else{
                    printf("\n=== Now recording!! Please speak into the microphone. ===\n"); 
                    fflush(stdout);
                    while( ( err = Pa_IsStreamActive( stream ) ) == 1 ){
                        Pa_Sleep(1000);                                     //record for 1 sec
                        printf("index = %d\n", data.frameIndex ); 
                        fflush(stdout);
                    }
                    if( err < 0 ){
                    }
                    else{
                        err = Pa_CloseStream( stream );
                        if( err != paNoError ){
                        }
                        else{
                            //Measure maximum peak amplitude.
                            max = 0;
                            average = 0.0;
                            for( i=0; i<numSamples; i++ ){
                                val = data.recordedSamples[i];
                                if( val < 0 ){
                                    val = -val; /* ABS */
                                }
                                if( val > max ){
                                    max = val;
                                }
                                average += val;
                            }
                            average = average / (double)numSamples;

                            printf("sample max amplitude = "PRINTF_S_FORMAT"\n", max );
                            printf("sample average = %lf\n", average );

                            //Write recorded data to a file.

                            
                            file = sf_open ("rec.wav", SFM_WRITE, &sfinfo);
                            sf_write_float (file,data.recordedSamples , sfinfo.channels * SAMPLE_COUNT);
                            sf_close (file) ;
                            printf("Wrote data to 'rec.wav'\n");
                            
                            // Playback recorded data. 
                            data.frameIndex = 0;

                            outputParameters.device = Pa_GetDefaultOutputDevice(); //default output device 
                            if (outputParameters.device == paNoDevice) {
                                fprintf(stderr,"Error: No default output device.\n");
                            }
                            else{
                                outputParameters.channelCount = 2;                     // stereo output
                                outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
                                outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
                                outputParameters.hostApiSpecificStreamInfo = NULL;

                                printf("\n=== Now playing back. ===\n"); 
                                fflush(stdout);
                                err = Pa_OpenStream(&stream,
                                                    NULL, 
                                                    &outputParameters,
                                                    SAMPLE_RATE,
                                                    FRAMES_PER_BUFFER,
                                                    paClipOff,     
                                                    playCallback,
                                                    &data);
                                if( err != paNoError ){
                                }
                                else if( stream ){
                                    err = Pa_StartStream( stream );
                                    if( err != paNoError ){
                                    }
                                    else{
                                        printf("Waiting for playback to finish.\n"); 
                                        fflush(stdout);
                                        while( ( err = Pa_IsStreamActive( stream ) ) == 1 ){
                                            Pa_Sleep(100);
                                        }
                                        if( err < 0 ){
                                        }
                                        else{         
                                            err = Pa_CloseStream( stream );
                                            if( err != paNoError ){
                                            }
                                            else{
                                                printf("Done.\n"); fflush(stdout);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Pa_Terminate();
    if( data.recordedSamples ){      
        free( data.recordedSamples );
    }
    if( err != paNoError ){
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          
    }
    return err;
}