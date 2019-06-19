/* References */
//https://stackoverflow.com/questions/13660777/c-reading-the-data-part-of-a-wav-file
//http://www.cplusplus.com/forum/beginner/166954/
//http://www.dspguide.com/CH18.PDF

/*  Imports  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include "WavFile.h"

/*  Namespaces  */
using namespace std;

/*  Constants  */
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr
#define INT2DOUBLE (32768.0)

/*  Method signatures  */
double* convolve(double x[], int N, double h[], int M, int P, int arraySize);
void four1(double data[], int nn, int isign);
void writeWavFile(char *fileName, double data[], int numberSamples, int outputRate);
size_t fwriteIntLSB(int data, FILE *fileStream);
size_t fwriteShortLSB(short data, FILE* fileStream);
void readWAV(char* fileName);

int main(int argc, char* argv[])
{
	if(argc<4) 
	{
		printf("Incorrect number of arguments.\n");
		printf("Need 3 arguments: inputfile, IRfile, outputfile.\n");
		return 1;
	}
	
	char *inputFileName = argv[1];
	char *irFileName = argv[2];
	char *outputFileName = argv[3];
	
	double complexData[SIZE*2];
	clock_t t;
	
	printf("Reading Files...\n");
	printf("Sound File:\n");
	t = clock();
	WavFile* inputFile = new WavFile();
	inputFile->readWAV(inputFileName);
	float benchmark1 = clock() - t;
	
	printf("Impulse Response File:\n");
	t = clock();
	WavFile* irFile = new WavFile();
	irFile->readWAV(irFileName);
	float benchmark2 = clock() - t;

	double* x = new double[inputFile->signalSize];
    for(int i = 0; i < (inputFile->signalSize); i++)
	{
        x[i] = ((double) inputFile->signal[i])/INT2DOUBLE;
    }
	double* h = new double[irFile->signalSize];
	for(int i = 0; i < (irFile->signalSize); i++)
	{
        h[i] = ((double) irFile->signal[i])/INT2DOUBLE;
    }

	unsigned int P = inputFile->signalSize + irFile->signalSize - 1;
	unsigned int arraySize = 0;
	if(inputFile->signalSize <= irFile->signalSize)
	{
		arraySize = irFile->signalSize;
	}
	else
	{
		arraySize = inputFile->signalSize;
	}
	int pow2 = 1;
	while(pow2 < arraySize)
	{
		pow2 *= 2;
	}
	arraySize = pow2*2;

	printf("Convoluting...\n");
	t = clock();
	double* y = convolve(x, inputFile->signalSize, h, irFile->signalSize, P, arraySize);
	float benchmark3 = clock() - t;
	
	printf("Outputting File...\n");
	t = clock();
	double inputMaxValue = 0.0;
    double outputMaxValue = 0.0;
    for(int i = 0; i < P; i++)
	{
        if(inputFile->signal[i] > inputMaxValue)
		{
            inputMaxValue = inputFile->signal[i];
		}
        if(y[i] > outputMaxValue)
		{
            outputMaxValue = y[i];
		}
    }
    for(int i  = 0; i < P; i++)
	{
        y[i] = y[i] / outputMaxValue * inputMaxValue;
    }
	writeWavFile(outputFileName, y, P, 44100);
	float benchmark4 = clock() - t;
	
	printf("Done!\n");
	printf("------------------------------------Overall Statistics--------------------------------------\n");
	printf("Time to read input file: %fs\n", benchmark1/CLOCKS_PER_SEC);
	printf("Time to read impulse response file: %fs\n", benchmark2/CLOCKS_PER_SEC);
	printf("Time to convolve data: %fs\n", benchmark3/CLOCKS_PER_SEC);
	printf("Time to write to output file: %fs\n", benchmark4/CLOCKS_PER_SEC);
	printf("--------------------------------------------------------------------------------------------");
    return 0;
}

double* convolve(double x[], int N, double h[], int M, int P, int arraySize)
{
	clock_t t;
	//Arrays of the form [Re(entry1), Im(entry1),...,Re(entryn), Im(entryn)]
	double* XX = new double[arraySize];
	double* FR = new double[arraySize];
	double* y = new double[arraySize];
	for(int i  = 0; i < arraySize; i++)
	{
	   XX[i] = 0.0;
	}
	for(int i = 0; i < arraySize; i++)
	{
		FR[i] = 0.0;
	}
	for(int i = 0; i < N; i++)
	{
		XX[i*2] = x[i];
	}
	for(int i = 0; i < M; i++)
	{
		FR[i*2] = h[i];
	}
	
	t = clock();
	four1(XX-1, arraySize/2, 1);
	float benchmark1 = clock() - t;
	
	t = clock();
	four1(FR-1, arraySize/2, 1);
	float benchmark2 = clock() - t;
	
	t = clock();
	for(int i = 0; i < arraySize; i+= 2)
	{
        y[i] = (XX[i]*FR[i])-(XX[i+1]*FR[i+1]);
        y[i+1] = (XX[i+1]*FR[i])+(XX[i]*FR[i+1]);
    
    }
	float benchmark3 = clock() - t;
	
	t = clock();
	four1(y-1, arraySize/2, -1); 
	float benchmark4 = clock() - t;
	
	double outputMaxValue = 0.0;
	for(int i = 0; i < P; i++)
	{
		double signalSample = y[i*2];
		y[i] = signalSample;
		if(outputMaxValue < abs(signalSample))
		{
			outputMaxValue = abs(signalSample);
		}
	}

	printf("-----------------------------------Convolution Statistics-----------------------------------\n");
	printf("Time to fourier transform input file data: %fs\n", benchmark1/CLOCKS_PER_SEC);
	printf("Time to fouier transform impulse response file data: %fs\n", benchmark2/CLOCKS_PER_SEC);
	printf("Time to multiply frequency spectrum by frequency response: %fs\n", benchmark3/CLOCKS_PER_SEC);
	printf("Time to inverse fouier transform above product: %fs\n", benchmark4/CLOCKS_PER_SEC);
	printf("--------------------------------------------------------------------------------------------\n");
	return y;
}

void four1(double data[], int nn, int isign){

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

void writeWavFile(char *fileName, double data[], int numberSamples, int outputRate)
{
	FILE *outputFile = fopen(fileName, "w");
	if (fileName == NULL)
	{
		printf("Error opening %s\n", outputFile);
		return;
	}

	unsigned int channels = 1;
	unsigned int bytesPerSample = 2;
    int dataChunkSize = channels * numberSamples * bytesPerSample;
    int formSize = 36 + dataChunkSize;
    short int frameSize = channels * bytesPerSample;
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    fputs("RIFF", outputFile);
    fwriteIntLSB(formSize, outputFile);
    fputs("WAVE", outputFile);

    fputs("fmt ", outputFile);
    fwriteIntLSB(16, outputFile);
    fwriteShortLSB(1, outputFile);
    fwriteShortLSB((short)channels, outputFile);
    fwriteIntLSB((int)outputRate, outputFile);
    fwriteIntLSB(bytesPerSecond, outputFile);
    fwriteShortLSB(frameSize, outputFile);
    fwriteShortLSB(8*bytesPerSample, outputFile);

    fputs("data", outputFile);
    fwriteIntLSB(dataChunkSize, outputFile);
	for(int i = 0; i < numberSamples; i++)
	{
		fwriteShortLSB((short) data[i], outputFile);;
	}
	
	fclose(outputFile);
	return;
}

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}