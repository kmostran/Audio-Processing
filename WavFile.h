/* References */
//https://stackoverflow.com/questions/16075233/reading-and-processing-wav-file-data-in-c-c

/*  Imports  */
#ifndef WavFileH
#define WavFileH

class WavFile 
{
	public:
		//header
		char chunk_id[4];
		int chunk_size;
		char format[4];
		char subchunk1_id[4];
		int subchunk1_size;
		short audio_format;
		short num_channels;
		int sample_rate;
		int byte_rate;
		short block_align;
		short bits_per_sample;
		char subchunk2_id[4];
		int subchunk2_size;

		//data
		char* dataArray;
		short* signal;
		int signalSize;

	public:
		void readWAV(char* fileName);
};

#endif