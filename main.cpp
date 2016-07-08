#include"wavefunction.h"

int main(int argc, char *argv[])
{
    WaveFunction* a=new WaveFunction(128,13);//每帧多少个采样点，MFCC参数的维数
    
    vector<vector<float> > mfccs1 = a->getMFCCs("D:\\1.wav");//提取mfcc参数
    vector<vector<float> > mfccs2 = a->getMFCCs("D:\\2.wav");
    
    cout<<a->ComputeDTW(mfccs1,mfccs2);//利用动态时间规整算法，计算两个语音的相似度，越小相似度越大

    return 0;
}
