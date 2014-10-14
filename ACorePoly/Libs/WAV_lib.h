#ifndef WAV_LIB_H
#define WAV_LIB_H

typedef unsigned short      WORD;
typedef unsigned long       DWORD;

struct WavHeader
{
	DWORD ChunkID;
	DWORD ChunkSize;
	DWORD Format;
	DWORD Subchunk1ID;
	DWORD Subchunk1Size;
	WORD  AudioFormat;
	WORD  NumChannels;
	DWORD SampleRate;
	DWORD ByteRate;
	WORD  BlockAlign;
	WORD  BitsPerSample;
	DWORD Subchunk2ID;
	DWORD Subchunk2Size;
};


struct WavSound
{
	int numChannels;
	int bytesPerSample;
	int sampleRate;
	int length;
	double * data;
};

int LoadWAV(const char * fName, WavSound * wav);
int SaveWAV(const char *fName, int len, int sampRate, int bufSize, double *buf);
int DeleteWAV(WavSound * wav);





#endif