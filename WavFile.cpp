/*  Imports  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "WavFile.h"

/*  Namespaces  */
using namespace std;

void WavFile::readWAV(char* fileName)
{
	FILE *wavFile = fopen(fileName,"rb");
	if (wavFile == NULL)
	{
		printf("Error opening %s\n", wavFile);
		return;
	}	
	
    ifstream wav(fileName, ios::in | ios:: binary);
    wav.seekg(0, ios::beg);

    wav.read(chunk_id, 4);
    wav.read((char*) &chunk_size, 4);
    wav.read(format, 4);
    wav.read(subchunk1_id, 4);
    wav.read((char*) &subchunk1_size, 4);
    wav.read((char*) &audio_format, sizeof(short));
    wav.read((char*) &num_channels, sizeof(short));
    wav.read((char*) &sample_rate, sizeof(int));
    wav.read((char*) &byte_rate, sizeof(int));
    wav.read((char*) &block_align, sizeof(short));
    wav.read((char*) &bits_per_sample, sizeof(short));

    if(subchunk1_size == 18)
	{
        short temp;
        wav.read((char*) &temp, sizeof(short));
    }

    wav.read(subchunk2_id, 4);
    wav.read((char*) &subchunk2_size, sizeof(int));
    int dataSize = subchunk2_size; 
    dataArray = new char[dataSize]; 
    wav.read(dataArray, dataSize);
    wav.close();

    signal = NULL;
    if(bits_per_sample == 8)
	{
        signalSize = dataSize;
        signal = new short[signalSize];
        for(int i = 0; i < dataSize; i++)
		{
            signal[i] = (short) ((unsigned char) dataArray[i]);
        }
    }
    else 
	{
        signalSize = dataSize/2;
        signal = new short[signalSize];
        short shortData;
        for(int i = 0; i < dataSize; i+= 2)
		{
            shortData = (short) ((unsigned char) dataArray[i]);
            shortData = (short) ((unsigned char) dataArray[i+1]) *256; 
            signal[i/2] = shortData;
        }
    }
}