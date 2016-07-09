#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include"WaveStructure.h"

class WaveFunction{
private:
    int CharCmp(const char *first,const char *second,unsigned short len);
    void InitHamming();

    //设置滤波器参数
    //输入参数：无
    //输出参数：FiltCoe1,三角形滤波器左边的系数,FiltCoe2三角形滤波器右边的系数,Num决定每个点属于哪一个滤波器
    void InitFilt(float *FiltCoe1, float *FiltCoe2, int *Num);

    //加窗，输入为buf,输出为data
    void HammingWindow(short* buf,float* data);

    //计算傅里叶参数
    void ComputeFFT(float *buffer,vector<complex<float> >& vecList);

    //傅里叶变换
    void FFT(const unsigned long & ulN, vector<complex<float> >& vecList);

    /*
    根据滤波器参数计算频带能量
    输入参数：*spdata  ---预处理之后的一帧语音信号
              *FiltCoe1---三角形滤波器左边的系数
              *FiltCoe2---三角形滤波器右边的系数
              *Num     ---决定每个点属于哪一个滤波器

    输出参数：*En      ---输出对数频带能量
    */
    void Filt(float *spdata, float *FiltCoe1, float *FiltCoe2, int *Num, float *En,vector<complex<float> >& vecList);

    /*
    计算MFCC系数
    输入参数：*En ---对数频带能量
    */
    void MFCC(float *En);


    float ComputeDTW(float *cep1, float *cep2, int num1, int num2);
    float Distance(float * ps1,float * ps2,int k1,int k2);

    void AdjustSize();

private:
    vector<float> xishu;
    double *Hamming;
    vector<vector<float> > SourceMFCCs;
    int MFCC_P;
    int MFCC_Pf;
    int FrmLen;
    int FFTLen;

public:
    vector<vector<float> > getMFCCs(string filename);
    vector<vector<float> > addFirstOrderDifference(vector<vector<float> > mfccs);
    vector<vector<float> > addOrderDifference(vector<vector<float> > mfccs);

    //输入为两个mfcc参数，cep1,cep2,返回最短距离
    float ComputeDTW(vector<vector<float> > cep1,vector<vector<float> > cep2);
    WaveFunction(int frm_len, int mfcc_num);
    ~WaveFunction();
};

#endif // WAVEFUNCTION_H
