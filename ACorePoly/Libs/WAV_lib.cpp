#ifndef WAV_LIB_CPP
#define WAV_LIB_CPP

#include "WAV_lib.h"

//#include <fstream.h>
#include <stdio.h>
#include <math.h>


int LoadWAV(const char * fName, WavSound * wav)
{
	wav->bytesPerSample=wav->length=wav->numChannels=wav->sampleRate=0;
	wav->data=0;
	//ifstream is(fName, ios::binary);
	FILE * in;
	fopen_s(&in,fName,"rb");
	
	WavHeader hed;
	//is.read((char*)&hed,sizeof(WavHeader));
	fread(&hed,1,sizeof(WavHeader),in);


	unsigned long int L=hed.Subchunk2Size;
	char * data;
	data=(char*)malloc(L);
	if(data==0)
	{
		//is.close();
		fclose(in);
		return 0;
	}
	//is.read((char*)data,L);
	fread(data,1,L,in);

	//is.close();
	fclose(in);

	wav->numChannels=hed.NumChannels;
	wav->bytesPerSample=hed.BitsPerSample/8;
	wav->sampleRate=hed.SampleRate;
	wav->length=L/(hed.NumChannels*hed.BitsPerSample/8); // загружается только 1 канал
	wav->data=(double*)malloc(wav->length*sizeof(double));
	if(wav->data==0)
	{
		free(data);
		wav->bytesPerSample=wav->length=wav->numChannels=wav->sampleRate=0;
		wav->data=0;
		return 0;
	}


	if(wav->numChannels==1 && wav->bytesPerSample==1)
	{
		unsigned char * c = (unsigned char *)data;
		int l=wav->length;
		for(int i=0;i<l;i++)
		{
			wav->data[i]=double(c[i])-128.0;

		}
	}
	else if(wav->numChannels==1 && wav->bytesPerSample==2)
	{
		short int * n = (short int *)data; 
		int l=wav->length;
		for(int i=0;i<l;i++)
		{
			wav->data[i]=double(n[i])/256.0;
//			wav->data[i+l]=0;
		}
	}
	else if(wav->numChannels==1 && wav->bytesPerSample==3)
	{
		unsigned char * c = (unsigned char *)data;
		long int n;
		int l=wav->length;
		for(int i=0;i<l;i++)
		{
			((char*)&n)[0]=c[3*i+0]; // remake smartly
			((char*)&n)[1]=c[3*i+1];
			((char*)&n)[2]=c[3*i+2];
			((char*)&n)[3]=0;

			n=n << 8;
			wav->data[i]=double(n)/(2*32768.0*256);
		}
	}
	else if(wav->numChannels==2 && wav->bytesPerSample==3) // загружается только первая дорожка
	{
		unsigned char * c = (unsigned char *)data;
		long int n;
		int l=wav->length;
		for(int i=0;i<l;i++)
		{
			((char*)&n)[0]=c[6*i+0];
			((char*)&n)[1]=c[6*i+1];
			((char*)&n)[2]=c[6*i+2];
			((char*)&n)[3]=0;

			n=n << 8;
			wav->data[i]=double(n)/(2*32768.0*256);
		}
	}
	else
	{
		free(wav->data);
		free(data);
		wav->bytesPerSample=wav->length=wav->numChannels=wav->sampleRate=0;
		wav->data=0;
		return 0;
	}

	free(data);
	return 1;
}


int LoadWAV_16bit2ch(char * fName, WavSound * wav)
{
	wav->bytesPerSample=wav->length=wav->numChannels=wav->sampleRate=0;
	wav->data=0;
	//ifstream is(fName, ios::binary);
	FILE * in;
	fopen_s(&in,fName,"rb");
	
	WavHeader hed;
	//is.read((char*)&hed,sizeof(WavHeader));
	fread(&hed,1,sizeof(WavHeader),in);


	unsigned long int L=hed.Subchunk2Size;
	char * data;
	data=(char*)malloc(L);
	if(data==0)
	{
		//is.close();
		fclose(in);
		return 0;
	}
	//is.read((char*)data,L);
	fread(data,1,L,in);

	//is.close();
	fclose(in);

	wav->numChannels=hed.NumChannels;
	wav->bytesPerSample=hed.BitsPerSample/8;
	wav->sampleRate=hed.SampleRate;
	wav->length=L/(hed.NumChannels*hed.BitsPerSample/8); // загружается только 1 канал
	if(wav->numChannels!=2 || wav->bytesPerSample!=2)
	{
		free(data);
		return 0;
	}
	wav->data=(double*)malloc(wav->length*sizeof(double));
	if(wav->data==0)
	{
		free(data);
		wav->bytesPerSample=wav->length=wav->numChannels=wav->sampleRate=0;
		wav->data=0;
		return 0;
	}



	
	short int * n = (short int *)data; 
	int l=wav->length;
	for(int i=0;i<l;i++)
	{
		wav->data[i]=double(n[2*i])/256.0;
	}
	


	free(data);
	return 1;
}

int DeleteWAV(WavSound * wav)
{
	wav->bytesPerSample=0;
	wav->sampleRate=0;
	wav->length=0;
	free(wav->data);
	wav->data=0;
	return 0;
}



int SaveWAV(const char *fName, int len, int sampRate, int bufSize, double *buf)
{
	int *iBuf=(int*)malloc(bufSize*sizeof(int));

	WavHeader hed;
	hed.ChunkID=0x46464952;
	hed.ChunkSize= 36+3*len; // ????
	hed.Format=0x45564157;
	hed.Subchunk1ID=0x20746D66;
	hed.Subchunk1Size=16;
	hed.AudioFormat=1;
	hed.NumChannels=1;
	hed.SampleRate=sampRate;
	hed.ByteRate=sampRate*3;
	hed.BlockAlign=3;
	hed.BitsPerSample=24;
	hed.Subchunk2ID=0x61746164;
	hed.Subchunk2Size= 3 * len; // ????

	int i, val;

//  LOL +++++++ later
//	double mx=0;
//	for(i=0;i<bufSize;i++)
//		if(fabs(buf[i])>mx) mx=fabs(buf[i]);
//	mx*=1.01;

	unsigned char *c;
	for(i=0;i<bufSize;i++)
	{
		val=buf[i] * (2*32768.0); // /mx
//		val = val >> 8;
		iBuf[i]=val;
	}

	FILE *out;
	fopen_s(&out,fName,"wb");

	fwrite(&hed,sizeof(hed),1,out);
	for(i=0;i<len;i++)
	{
		fwrite(iBuf+(i%bufSize),3,1,out);
	}



	fclose(out);
	free(iBuf);
	return 0;
}




#endif