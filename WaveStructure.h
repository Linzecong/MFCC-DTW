#ifndef WAVESTRUCTURE_H
#define WAVESTRUCTURE_H

#include<vector>
#include<complex>
#include<cmath>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<malloc.h>
#include<bitset>
using   namespace   std;
#define INF 1000000000000
const float PI=3.1415926536;
const int FS=8;
const int FiltNum=25;

struct RIFF_HEADER{
   char szRiffID[4];  // 'R','I','F','F'
   unsigned long dwRiffSize;
   char szRiffFormat[4]; // 'W','A','V','E'
};

struct FMT_HEADER{
   char  szFmtID[4]; // 'f','m','t',' '
   unsigned long  dwFmtSize; // 18 OR 16
};

struct DATA_BLOCK{
    char szDataID[4]; // 'd','a','t','a'
    unsigned long dwDataSize;
};

struct WAVE_FORMAT{
   unsigned short  wFormatTag;
   unsigned short  wChannels;
   unsigned long dwSamplesPerSec;
   unsigned long dwAvgBytesPerSec;
   unsigned short  wBlockAlign;
   unsigned short  wBitsPerSample;
   unsigned short  wAppendData;
};

struct mTWavHeader{
        int rId;    //标志符（RIFF）
        int rLen;   //数据大小,包括数据头的大小和音频文件的大小
        int wId;    //格式类型（"WAVE"）
        int fId;    //"fmt"
        int fLen;   //Sizeof(WAVEFORMATEX)
        short wFormatTag;       //编码格式，包括WAVE_FORMAT_PCM，WAVEFORMAT_ADPCM等
        short nChannels;        //声道数，单声道为1，双声道为2
        int nSamplesPerSec;   //采样频率
        int nAvgBytesPerSec;  //每秒的数据量
        short nBlockAlign;      //块对齐
        short wBitsPerSample;   //WAVE文件的采样大小
        int dId;              //"data"
        int wSampleLength;    //音频数据的大小
};

#endif // WAVESTRUCTURE_H
