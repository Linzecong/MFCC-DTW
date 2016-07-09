#include "WaveFunction.h"

WaveFunction::WaveFunction(int frm_len, int mfcc_num){
    MFCC_P=mfcc_num;
    MFCC_Pf=float(mfcc_num);
    FrmLen=frm_len;
    FFTLen=frm_len;
    Hamming=new double[FrmLen];
}

WaveFunction::~WaveFunction(){
    delete []Hamming;
}

vector<vector<float> > WaveFunction::getMFCCs(string filename){
    xishu.clear();
    SourceMFCCs.clear();
    //mfcc分析
    mTWavHeader waveheader;
    FILE *sourcefile;
    short buffer[FrmLen];
    float data[FrmLen];
    float FiltCoe1[FFTLen/2+1];  //左系数
    float FiltCoe2[FFTLen/2+1];  //右系数
    int Num[FFTLen/2+1];     //决定每个点属于哪一个滤波器
    float En[FiltNum+1];         //频带能量
    vector<complex<float> > vecList;

    sourcefile=fopen(filename.c_str(),"rb");
    fread(&waveheader,sizeof(mTWavHeader),1,sourcefile);
    InitHamming();//初始化汉明窗
    InitFilt(FiltCoe1,FiltCoe2,Num); //初始化MEL滤波系数

    while(fread(buffer,sizeof(short),FrmLen,sourcefile)==FrmLen){
        HammingWindow(buffer,data);
        ComputeFFT(data,vecList);
        Filt(data, FiltCoe1, FiltCoe2, Num, En,vecList);
        MFCC(En);
        vecList.clear();
        fseek(sourcefile, -FrmLen/2, SEEK_CUR);//考虑到帧移，每次移动半帧
    }

    int stdlength=xishu.size();

    for(int i=0;i<stdlength/MFCC_P;i++){
        vector<float> temp;
        for(int j=0;j<MFCC_P;j++)
            temp.push_back(xishu[i*MFCC_P+j]);
        SourceMFCCs.push_back(temp);
    }
    fclose(sourcefile);
    return SourceMFCCs;
}

vector<vector<float> > addFirstOrderDifference(vector<vector<float> > mfccs){
    vector<vector<float> > temp;
    for(int i=0;i<mfccs.size();i++){
        vector<float> line=mfccs[i];
        int size=line.size();
        for(int t=0;t<size;t++){
            if(t<2)
                line.push_back(line[t+1]-line[t]);
            else{
                if(t>size-2||t==size-2)
                    line.push_back(line[t]-line[t-1]);
                else{
                    float fenzi=line[t+1]-line[t-1]+2*(line[t+2]-line[t-2]);
                    float fenmu=sqrtf(10);
                    line.push_back(fenzi/fenmu);

                }
            }
        }
        temp.push_back(line);
    }
    return temp;
}

vector<vector<float> > addOrderDifference(vector<vector<float> > mfccs){
    vector<vector<float> > temp;
    for(int i=0;i<mfccs.size();i++){
        vector<float> line=mfccs[i];
        int size=line.size();
        //一阶差分
        for(int t=0;t<size;t++){
            if(t<2)
                line.push_back(line[t+1]-line[t]);
            else{
                if(t>size-2||t==size-2)
                    line.push_back(line[t]-line[t-1]);
                else{
                    float fenzi=line[t+1]-line[t-1]+2*(line[t+2]-line[t-2]);
                    float fenmu=sqrtf(10);
                    line.push_back(fenzi/fenmu);

                }
            }
        }
        //二阶差分
        for(int t=size;t<size*2;t++){
            if(t<2)
                line.push_back(line[t+1]-line[t]);
            else{
                if(t>size-2||t==size-2)
                    line.push_back(line[t]-line[t-1]);
                else{
                    float fenzi=line[t+1]-line[t-1]+2*(line[t+2]-line[t-2]);
                    float fenmu=sqrtf(10);
                    line.push_back(fenzi/fenmu);

                }
            }
        }
        temp.push_back(line);
    }
    return temp;
}

float WaveFunction::ComputeDTW(vector<vector<float> > cep1, vector<vector<float> > cep2)
{
    vector<float> temp;
    for(int i=0;i<cep1.size();i++)
        for(int j=0;j<cep1[i].size();j++)
            temp.push_back(cep1[i][j]);
    int stdlength=temp.size();
    float * stdmfcc = new float[stdlength];
    std::copy(temp.begin(),temp.end(),stdmfcc);

    vector<float> temp1;
    for(int i=0;i<cep2.size();i++)
        for(int j=0;j<cep2[i].size();j++)
            temp1.push_back(cep2[i][j]);
    int testlen=temp1.size();
    float * testmfcc = new float[testlen];
    std::copy(temp1.begin(),temp1.end(),testmfcc);
    return ComputeDTW(stdmfcc,testmfcc,stdlength/MFCC_P,testlen/MFCC_P);
}

int WaveFunction::CharCmp(const char *first,const char *second,unsigned short len)
{
    int i=0;
    while((first[i]==second[i])&&(i++<len));
    if(i>=len)
        return 0;
    else if(first[i-1]>second[i-1])
        return 1;
    else
        return -1;
}

void WaveFunction::InitHamming(){

    float twopi;
    int i;
    twopi=8.0F*atan(1.0F);
    for(i=0;i<FrmLen;i++)
        Hamming[i]=(float)(0.54-0.46*cos((float)i*twopi/(float)(FrmLen-1)));
}

void WaveFunction::InitFilt(float *FiltCoe1, float *FiltCoe2, int *Num){
    int i,j;
    float Freq;
    int FiltFreq[FiltNum+1] = {0,100,200,300,400,500,600,700,800,900,1000,
                               1149,1320,1516,1741,2000,2297,2639,3031,3482,4000,
                               4595,5278,6063,6964,8001};//滤波器的中心频率
    int BW[FiltNum+1]={100,100,100,100,100,100,100,100,100,100,124,
                       160,184,211,242,278,320,367,422,484,556,
                       639,734,843,969,1112};//滤波器的带宽
    for(i = 0 ; i<= FFTLen/2 ; i++ )
    {
        Num[i]=0;
    }

    for(i = 0 ; i <= FFTLen/2 ; i++)
    {
        Freq = FS * 1000.0F * (float)(i) / (float)(FFTLen);
        for(j = 0 ; j <FiltNum ; j++)
        {
            if(Freq >= (float)FiltFreq[j] && Freq <= (float)FiltFreq[j+1])
            {
                Num[i] = j;
                if(j == 0)
                {
                    FiltCoe1[i] = 0.0F;
                }
                else
                {
                    FiltCoe1[i] = ((float)(FiltFreq[j]+BW[j])-Freq) / (float)(BW[j]);
                }
                FiltCoe2[i] = (Freq-(float)(FiltFreq[j+1]-BW[j+1])) / (float)(BW[j+1]);
                FiltCoe1[i] = FiltCoe1[i] * FiltCoe1[i];
                FiltCoe2[i] = FiltCoe2[i] * FiltCoe2[i];
                break;
            }
        }
    }

}

void WaveFunction::HammingWindow(short *buf, float *data){
    int i;
    for(i=0;i<FrmLen;i++)
        data[i]=buf[i]*Hamming[i];
}

void WaveFunction::ComputeFFT(float *data, vector<complex<float> > &vecList){
    for(int i=0;i<FFTLen;++i)
    {
        if(i<FrmLen)
        {
            complex<float> temp(data[i]);
            vecList.push_back(temp);
        }
        else
        {
            complex<float> temp(0);
            vecList.push_back(temp);
        }
    }
    FFT(FFTLen,vecList);
}

void WaveFunction::FFT(const unsigned long &ulN, vector<complex<float> > &vecList){
    //得到幂数

    unsigned long ulPower = 0; //幂数
    unsigned long ulN1 = ulN - 1;
    while(ulN1 > 0)
    {
        ulPower++;
        ulN1 /= 2;
    }
    //反序

    bitset<sizeof(unsigned long) * 8> bsIndex; //二进制容器
    unsigned long ulIndex; //反转后的序号
    unsigned long ulK;
    for(unsigned long p = 0; p < ulN; p++)
    {
        ulIndex = 0;
        ulK = 1;
        bsIndex = bitset<sizeof(unsigned long) * 8>(p);
        for(unsigned long j = 0; j < ulPower; j++)
        {
            ulIndex += bsIndex.test(ulPower - j - 1) ? ulK : 0;
            ulK *= 2;
        }

        if(ulIndex > p)
        {
            complex<float> c = vecList[p];
            vecList[p] = vecList[ulIndex];
            vecList[ulIndex] = c;
        }
    }

    //计算旋转因子

    vector<complex<float> > vecW;
    for(unsigned long i = 0; i < ulN / 2; i++)
    {
        vecW.push_back(complex<float>(cos(2 * i * PI / ulN) , -1 * sin(2 * i * PI / ulN)));
    }

    /*for(unsigned long m = 0; m < ulN / 2; m++)
    {
        cout<< "\nvW[" << m << "]=" << vecW[m];
    } */

    //计算FFT

    unsigned long ulGroupLength = 1; //段的长度
    unsigned long ulHalfLength = 0; //段长度的一半
    unsigned long ulGroupCount = 0; //段的数量
    complex<float> cw; //WH(x)
    complex<float> c1; //G(x) + WH(x)
    complex<float> c2; //G(x) - WH(x)
    for(unsigned long b = 0; b < ulPower; b++)
    {
        ulHalfLength = ulGroupLength;
        ulGroupLength *= 2;
        for(unsigned long j = 0; j < ulN; j += ulGroupLength)
        {
            for(unsigned long k = 0; k < ulHalfLength; k++)
            {
                cw = vecW[k * ulN / ulGroupLength] * vecList[j + k + ulHalfLength];
                c1 = vecList[j + k] + cw;
                c2 = vecList[j + k] - cw;
                vecList[j + k] = c1;
                vecList[j + k + ulHalfLength] = c2;
            }
        }
    }
}

void WaveFunction::Filt(float *spdata, float *FiltCoe1, float *FiltCoe2, int *Num, float *En, vector<complex<float> > &vecList){
    float temp=0;
    int id, id1, id2;

    for(id = 0 ; id <= FiltNum ; id++)
    {
        En[id]=0.0F;
    }
    for(id = 0 ; id < FFTLen/2 ; id++)
    {
        temp = vecList[id].real()*vecList[id].real()+vecList[id].imag()*vecList[id].imag();
        id1 = Num[id];
        id2 = id1+1;
        En[id1] = En[id1] + FiltCoe1[id] * temp;
        En[id2] = En[id2] + FiltCoe2[id] * temp;
    }
    for(id = 1 ; id <= FiltNum ; id++)
    {
        if(En[id]!=0)
            En[id]=(float)log(En[id]);
    }
}

void WaveFunction::MFCC(float *En)
{
    int idcep, iden;
    float Cep[MFCC_P];

    for(idcep = 0 ; idcep < MFCC_P ; idcep++)
    {
        Cep[idcep] = 0.0;

        for(iden = 1 ; iden <= FiltNum ; iden++)
        {
            Cep[idcep] = Cep[idcep] + En[iden] * (float)cos((idcep+1) * (iden-0.5F) * PI/(FiltNum));
        }
        Cep[idcep] = Cep[idcep] / 10.0F;
        xishu.push_back(Cep[idcep]);
    }
}


float WaveFunction::ComputeDTW(float *cep1, float *cep2, int num1, int num2){
    struct record
    {		int x;
                int y;
    };
    struct point
    {		int x,y;
                float minvalue;
                        int stepnum;
                                bool recheck;               //记录该点是否被记录过
    };
    record * re;
    record * newre;

    newre=new record[num1*num2];    //记录下一层的所有点
    re=new record[num1*num2];       //记录当层的所有点
    int renum;
    int newrenum=0;
    int i,j;
    point * poi;
    poi=new point[num1*num2];

    for(i=0;i<num1*num2;i++)
    {
        poi[i].recheck=0;
        poi[i].minvalue=INF;
        poi[i].stepnum=0;
    }								//设置初始值

    for(i=0;i<5;i++)                //起始点
    {
        if(i==0)  {	re[i].x=1; re[i].y=1; }
        if(i==1)  {	re[i].x=1; re[i].y=2; }
        if(i==2)  {	re[i].x=1; re[i].y=3; }
        if(i==3)  {	re[i].x=2; re[i].y=1; }
        if(i==4)  {	re[i].x=3; re[i].y=1; }
        poi[(re[i].y-1)*num1+re[i].x-1].minvalue=Distance(cep1,cep2,re[i].x,re[i].y);
        poi[(re[i].y-1)*num1+re[i].x-1].stepnum=1;
    }
    renum=5;
    int newx,newy;                   //newvalue;
    for(i=0;i<renum;i++)
    {
        for(j=0;j<3;j++)
        {
            if(j==0){ newx=re[i].x+1; newy=re[i].y+2; }
            if(j==1){ newx=re[i].x+1; newy=re[i].y+1; }
            if(j==2){ newx=re[i].x+2; newy=re[i].y+1; }

            /////////////三种可能路径

            if(newx>=num1||newy>=num2)
                continue;
            if(fabs(newx-newy)<=fabs(num1-num2)+3)
            {
                if(poi[(newy-1)*num1+newx-1].recheck==0)
                {
                    newre[newrenum].x=newx;
                    newre[newrenum].y=newy;
                    newrenum++;
                }
                float tmpdis;
                int addstepnum;
                if(j==0){ tmpdis=Distance(cep1,cep2,newx-1,newy-1)*2+Distance(cep1,cep2,newx,newy); addstepnum=2;}
                if(j==1){ tmpdis=Distance(cep1,cep2,newx,newy)*2; addstepnum=1;}
                if(j==2){ tmpdis=Distance(cep1,cep2,newx-1,newy-1)*2+Distance(cep1,cep2,newx,newy); addstepnum=2;}
                if(poi[(newy-1)*num1+newx-1].minvalue>(poi[(re[i].y-1)*num1+re[i].x-1].minvalue+tmpdis))
                {
                    poi[(newy-1)*num1+newx-1].minvalue=(poi[(re[i].y-1)*num1+re[i].x-1].minvalue+tmpdis);
                    poi[(newy-1)*num1+newx-1].stepnum=poi[(re[i].y-1)*num1+re[i].x-1].stepnum+addstepnum;
                }
                if(poi[(newy-1)*num1+newx-1].recheck==0)
                    poi[(newy-1)*num1+newx-1].recheck=1;
            }
        }
        if(newrenum!=0 && i>=(renum-1))
        {
            renum=newrenum;
            newrenum=0;
            struct	record * tt;
            tt=re;
            re=newre;
            newre=tt;
            i=-1;
        }
    }
    float min=INF;
    for(j=0;j<renum;j++)
    {
        if((poi[(re[j].y-1)*num1+re[j].x-1].minvalue)/poi[(re[j].y-1)*num1+re[j].x-1].stepnum<min)
            min=(poi[(re[j].y-1)*num1+re[j].x-1].minvalue)/poi[(re[j].y-1)*num1+re[j].x-1].stepnum;
    }

    //	min;
    delete []poi;
    delete []newre;
    delete []re;
    delete []cep1;
    delete []cep2;
    return min;
}

float WaveFunction::Distance(float *ps1, float *ps2, int k1, int k2){
    int i=0;
    float sum=0;
    for(i=0;i<MFCC_P;i++)
        sum+=(1+MFCC_Pf/2*(float)sin(PI*i/MFCC_Pf))*(ps1[k1+i]-ps2[k2+i])*(ps1[k1+i]-ps2[k2+i]);

    return sum;
}
