
#include "CATS.h"

using namespace std;

CATS::CATS():
    Pi(3.141592653589793),
    AlphaFS(0.0072973525664),
    RevSqrt2(1./sqrt(2.)),
    FmToNu(5.067731237e-3),NuToFm(197.3269602)
    {
    IdenticalParticles = false;
    Q1Q2 = 0;
    RedMass = 0;
    pdgID[0] = 0;
    pdgID[1] = 0;
    NumCh = 0;
    NumMomBins= 0;
    StartRad = 0.005*FmToNu;
    EpsilonProp = 5e-6;
    EpsilonConv = 5e-6;
    MaxRad = 32.*FmToNu;
    MaxRho = 16;
    ExcludeFailedConvergence = true;
    NumPairs = 0;
    WeightIp = NULL;
    WeightIpError = NULL;
    LoadingComplete = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    MaxPairsPerBin = 8e3;
    MaxPairsToRead = 4294967295;
//    MaxPairsToLoad = 4294967295;
    EventMixing = false;
    BufferEventMix = 64;
    TauCorrection = false;
    UseAnalyticSource = false;
    ThetaDependentSource = false;
    TransportRenorm = 1;
    MinTotPairMom = -1;
    MaxTotPairMom = 1e100;
    LoadedMinTotPairMom = -1;
    LoadedMaxTotPairMom = 1e100;
    UseTotMomCut = false;

    NumPW = NULL;
    MomBin = NULL;
    IpBin = NULL;
    ChannelWeight = NULL;
    SavedWaveFunBins = NULL;
    PhaseShift = NULL;
    WaveFunRad = NULL;
    WaveFunctionU = NULL;
    MomBinConverged = NULL;
    InputFileName = NULL;
    RelativeMomentum = NULL;
    RelativePosition = NULL;
    RelativeCosTheta = NULL;
    TotalPairMomentum = NULL;
    LoadedPairsPerBin = NULL;

    ShortRangePotential = NULL;

    AnalyticSource = NULL;

    CorrFun=NULL;
    CorrFunError=NULL;

    PotPar = NULL;
    AnaSourcePar = NULL;

    SetIpBins(1, -1000, 1000);
}

CATS::~CATS(){
    DelAll();
    if(NumPW) {delete[]NumPW; NumPW=NULL;}
    if(MomBin) {delete[]MomBin; MomBin=NULL;}
    if(IpBin) {delete[]IpBin; IpBin=NULL;}
    if(ChannelWeight) {delete[]ChannelWeight; ChannelWeight=NULL;}
    if(InputFileName) {delete[]InputFileName; InputFileName=NULL;}

    if(ShortRangePotential){
        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            if(ShortRangePotential[usCh]){
                delete[]ShortRangePotential[usCh];
                ShortRangePotential[usCh] = NULL;

                delete[]PotPar[usCh];
                PotPar[usCh] = NULL;
            }
        }
        delete[]ShortRangePotential; ShortRangePotential=NULL;
        delete[]PotPar; PotPar=NULL;
    }
}

void CATS::DelMomIpMp(){
    if(!RelativeMomentum) return;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            delete [] RelativeMomentum[uMomBin][uIpBin];
            delete [] RelativePosition[uMomBin][uIpBin];
            delete [] RelativeCosTheta[uMomBin][uIpBin];
        }
        delete [] RelativeMomentum[uMomBin];
        delete [] RelativePosition[uMomBin];
        delete [] RelativeCosTheta[uMomBin];
    }
    delete [] RelativeMomentum; RelativeMomentum=NULL;
    delete [] RelativePosition; RelativePosition=NULL;
    delete [] RelativeCosTheta; RelativeCosTheta=NULL;

    RemoveTotPairMomCut();
}

void CATS::DelIp(){
    if(!WeightIp) return;
    delete [] WeightIp; WeightIp=NULL;
    delete [] WeightIpError; WeightIpError=NULL;
}

void CATS::DelMomChPw(){
    if(!SavedWaveFunBins) return;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
                if(WaveFunRad[uMomBin][usCh][usPW])
                    delete [] WaveFunRad[uMomBin][usCh][usPW];
                if(WaveFunctionU[uMomBin][usCh][usPW])
                    delete [] WaveFunctionU[uMomBin][usCh][usPW];
            }
            delete [] SavedWaveFunBins[uMomBin][usCh];
            delete [] PhaseShift[uMomBin][usCh];
            delete [] WaveFunRad[uMomBin][usCh];
            delete [] WaveFunctionU[uMomBin][usCh];
        }
        delete [] SavedWaveFunBins[uMomBin];
        delete [] PhaseShift[uMomBin];
        delete [] WaveFunRad[uMomBin];
        delete [] WaveFunctionU[uMomBin];
    }
    delete [] SavedWaveFunBins; SavedWaveFunBins=NULL;
    delete [] PhaseShift; PhaseShift=NULL;
    delete [] WaveFunRad; WaveFunRad=NULL;
    delete [] WaveFunctionU; WaveFunctionU=NULL;
}

void CATS::DelMom(){
    if(!MomBinConverged) return;
    delete [] MomBinConverged; MomBinConverged=NULL;
}

void CATS::DelMomIp(){
    if(CorrFun){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete [] CorrFun[uMomBin];
            delete [] CorrFunError[uMomBin];

        }
        delete [] CorrFun; CorrFun=NULL;
        delete [] CorrFunError; CorrFunError=NULL;
    }
    if(LoadedPairsPerBin){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete [] LoadedPairsPerBin[uMomBin];
        }
        delete [] LoadedPairsPerBin; LoadedPairsPerBin=NULL;
    }
}

void CATS::DelAllMom(){
    DelMomIpMp();
    DelMomChPw();
    DelMom();
    DelMomIp();
}

void CATS::DelAllIp(){
    DelMomIpMp();
    DelIp();
    DelMomIp();
}

void CATS::DelAll(){
    DelAllMom();
    DelAllIp();
}

void CATS::SetRedMass(const double& redMass){
    if(redMass==RedMass) return;
    RedMass = redMass;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetRedMass(){
    return RedMass;
}

void CATS::SetPdgId(const int& id1, const int& id2){
    if(pdgID[0]==id1 && pdgID[1]==id2) return;
    pdgID[0] = id1;
    pdgID[1] = id2;
    if(id1==id2) IdenticalParticles=true;
    else IdenticalParticles=false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::GetPdgId(int& id1, int& id2){
    id1 = pdgID[0];
    id2 = pdgID[1];
}

//If the number of channels is changed, all previous input about the
//channels themselves is lost (i.e. NumPW, WhichPartialWave and the potentials are reset!)
void CATS::SetNumChannels(const unsigned short& numCh){
    if(NumCh == numCh) return;
    if(!numCh){
        printf("WARNING: Bad input in CATS::SetNumChannels(unsigned short numCh)\n");
        printf("         NumCh cannot be zero!\n");
        return;
    }

    if(NumPW) {delete[]NumPW; NumPW=NULL;}
    NumPW = new unsigned short [numCh];

    if(ShortRangePotential){
        for(unsigned short usCh=0; usCh<numCh; usCh++){
            if(ShortRangePotential[usCh]){
                delete [] ShortRangePotential[usCh];
                delete [] PotPar[usCh];
            }
        }
        delete[]ShortRangePotential; ShortRangePotential=NULL;
        delete[]PotPar; PotPar=NULL;
    }
    ShortRangePotential = new CatsPotential* [numCh];
    PotPar = new double** [numCh];
    for(unsigned short usCh=0; usCh<numCh; usCh++){
        ShortRangePotential[usCh] = NULL;
        PotPar[usCh] = NULL;
    }

    if(ChannelWeight) {delete[]ChannelWeight; ChannelWeight=NULL;}
    ChannelWeight = new double [numCh];
    for(unsigned short usCh=0; usCh<numCh; usCh++){
        NumPW[usCh] = 0;
        ChannelWeight[usCh] = 0;
    }

    DelMomChPw();
    NumCh = numCh;

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
unsigned short CATS::GetNumChannels(){
    return NumCh;
}

void CATS::SetNumPW(const unsigned short& usCh, const unsigned short& numPW){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::SetNumPW(unsigned short usCh, unsigned short numPW)\n");
        return;
    }
    if(NumPW[usCh]==numPW) return;
    DelMomChPw();
    NumPW[usCh] = numPW;

    if(ShortRangePotential[usCh]) delete[]ShortRangePotential[usCh];
    ShortRangePotential[usCh] = new CatsPotential [numPW];
    for(unsigned short usPW=0; usPW<numPW; usPW++){
        ShortRangePotential[usCh][usPW] = 0;
    }

    if(PotPar[usCh]) delete[]PotPar[usCh];
    PotPar[usCh] = new double* [numPW];

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
unsigned short CATS::GetNumPW(const unsigned short& usCh){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::GetNumPW(unsigned short usCh)\n");
        return 0;
    }
    return NumPW[usCh];
}

void CATS::SetQ1Q2(const int& q1q2){
    if(Q1Q2==q1q2) return;
    Q1Q2 = q1q2;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
int CATS::GetQ1Q2(){
    return Q1Q2;
}

unsigned CATS::GetNumMomBins(){
    return NumMomBins;
}

unsigned CATS::GetNumIpBins(){
    return NumIpBins;
}

unsigned CATS::GetNumPairs(){
    return NumPairs;
}

void CATS::SetMomBins(const unsigned& nummombins, const double* mombins){
    if(!nummombins){
        printf("WARNING: Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
        return;
    }
    if(!mombins){
        printf("WARNING: Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
        return;
    }
    if(nummombins!=NumMomBins || !MomBin){
        if(MomBin) {delete[]MomBin; MomBin=NULL;}
        MomBin = new double [nummombins+1];
        DelAllMom();
        NumMomBins = nummombins;
    }
    LoadingComplete = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    for(unsigned uBin=0; uBin<=NumMomBins; uBin++){
        MomBin[uBin] = mombins[uBin];
        if(MomBin[uBin]<0){
            printf("ERROR: CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
            printf("         The momentum should be positive!\n");
            return;
        }
    }
}
void CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom){
    if(!nummombins){
        printf("WARNING: Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    if(MinMom>MaxMom){
        printf("WARNING: Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    if(MinMom==MaxMom && nummombins!=1){
        printf("WARNING: Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    if(nummombins!=NumMomBins || !MomBin){
        if(MomBin) {delete[]MomBin; MomBin=NULL;}
        MomBin = new double [nummombins+1];
        DelAllMom();
        NumMomBins = nummombins;
    }
    LoadingComplete = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    double BinWidth = (MaxMom-MinMom)/double(NumMomBins);
    for(unsigned uBin=0; uBin<=NumMomBins; uBin++){
        MomBin[uBin] = MinMom+uBin*BinWidth;
    }
}

void CATS::SetIpBins(const unsigned& numBbins, const double* imppar){
    if(!numBbins){
        printf("WARNING: Bad input in CATS::SetIpBins(const unsigned& numBbins, const double* imppar)\n");
        return;
    }
    if(!imppar){
        printf("WARNING: Bad input in CATS::SetIpBins(const unsigned& numBbins, const double* imppar)\n");
        return;
    }
    if(numBbins!=NumIpBins || !IpBin){
        if(IpBin) {delete[]IpBin; IpBin=NULL;}
        IpBin = new double [numBbins+1];
        DelAllIp();
        NumIpBins = numBbins;
    }
    LoadingComplete = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    for(unsigned uBin=0; uBin<=NumIpBins; uBin++){
        IpBin[uBin] = imppar[uBin];
    }
}
void CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar){
    if(!numBbins){
        printf("WARNING: Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(MinImpPar>MaxImpPar){
        printf("WARNING: Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(MinImpPar==MaxImpPar && numBbins!=1){
        printf("WARNING: Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(numBbins!=NumIpBins || !IpBin){
        if(IpBin) {delete[]IpBin; IpBin=NULL;}
        IpBin = new double [numBbins+1];
        DelAllIp();
        NumIpBins = numBbins;
    }
    LoadingComplete = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    double BinWidth = (MaxImpPar-MinImpPar)/double(NumIpBins);
    for(unsigned uBin=0; uBin<=NumIpBins; uBin++){
        IpBin[uBin] = MinImpPar+double(uBin)*BinWidth;
    }
}

void CATS::SetChannelWeight(const unsigned short& usCh, const double& weight){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::SetSpinWeight(const unsigned short& usCh, const double& weight)\n");
        return;
    }
    if(ChannelWeight[usCh]==weight) return;
    ChannelWeight[usCh] = weight;
    ComputedCorrFunction = false;
}

double CATS::GetChannelWeight(const unsigned short& usCh){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::GetSpinWeight(const unsigned short& usCh)\n");
        return 0;
    }
    return ChannelWeight[usCh];
}

void CATS::SetStartRad(const double& srad){
    if(StartRad==fabs(srad)) return;
    StartRad = fabs(srad);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetStartRad(){
    return StartRad;
}

void CATS::SetEpsilonProp(const double& epsp){
    if(EpsilonProp==fabs(epsp)) return;
    //make sure that EpsilonProp is always non-zero and positive
    EpsilonProp = fabs(epsp)+1e-64;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetEpsilonProp(){
    return EpsilonProp;
}

void CATS::SetEpsilonConv(const double& epsc){
    if(EpsilonConv==fabs(epsc)) return;
    EpsilonConv = fabs(epsc);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetEpsilonConv(){
    return EpsilonConv;
}

void CATS::SetMaxRad(const double& maxrad){
    if(MaxRad==fabs(maxrad*FmToNu)) return;
    MaxRad = fabs(maxrad*FmToNu);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

double CATS::GetMaxRad(){
    return MaxRad;
}

void CATS::SetMaxRho(const double& maxrho){
    if(MaxRho==fabs(maxrho)) return;
    MaxRho = fabs(maxrho);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

double CATS::GetMaxRho(){
    return MaxRho;
}

void CATS::SetExcludeFailedBins(const bool& efb){
    if(ExcludeFailedConvergence==efb) return;
    ExcludeFailedConvergence = efb;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
bool CATS::GetExcludeFailedBins(){
    return ExcludeFailedConvergence;
}

void CATS::SetMaxPairsPerBin(unsigned mpp){
    if(!mpp){
        printf("WARNING: MaxPairsPerBin cannot be zero, setting MaxPairsPerBin=1\n");
        mpp=1;
    }
    if(MaxPairsPerBin==mpp) return;
    DelMomIpMp();
    MaxPairsPerBin = mpp;
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned CATS::GetMaxPairsPerBin(){
    return MaxPairsPerBin;
}

void CATS::SetMaxPairsToRead(unsigned mpp){
    if(!mpp){
        printf("WARNING: MaxPairsToRead cannot be zero, setting MaxPairsToRead=1\n");
        mpp=1;
    }
    if(MaxPairsToRead==mpp) return;
    MaxPairsToRead = mpp;
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned CATS::GetMaxPairsToRead(){
    return MaxPairsToRead;
}

//void CATS::SetMaxPairsToLoad(unsigned mpp){
//    if(!mpp){
//        printf("WARNING: MaxPairsToLoad cannot be zero, setting MaxPairsToLoad=1\n");
//        mpp=1;
//    }
//    if(MaxPairsToLoad==mpp) return;
//    MaxPairsToLoad = mpp;
//    LoadingComplete = false;
//    if(!UseAnalyticSource) ComputedCorrFunction = false;
//}
//unsigned CATS::GetMaxPairsToLoad(){
//    return MaxPairsToLoad;
//}

void CATS::SetEventMixing(const bool& mix){
    if(EventMixing==mix) return;
    EventMixing = mix;
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
bool CATS::GetEventMixing(){
    return EventMixing;
}

void CATS::SetBufferEventMix(const unsigned& bem){
    if(BufferEventMix<2){
        printf("WARNING: BufferEventMix is too low, using BufferEventMix=2\n");
        BufferEventMix = 2;
        return;
    }
    if(BufferEventMix==bem) return;
    BufferEventMix = bem;
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned CATS::GetBufferEventMix(){
    return BufferEventMix;
}

void CATS::SetTauCorrection(const bool& tc){
    if(TauCorrection==tc) return;
    TauCorrection = tc;
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
bool CATS::GetTauCorrection(){
    return TauCorrection;
}

void CATS::SetUseAnalyticSource(const bool& val){
    if(UseAnalyticSource==val) return;
    UseAnalyticSource = val;
    ComputedCorrFunction = false;
}
bool CATS::GetUseAnalyticSource(){
    return UseAnalyticSource;
}

void CATS::SetThetaDependentSource(const bool& val){
    if(ThetaDependentSource==val) return;
    ThetaDependentSource = val;
    ComputedCorrFunction = false;
}
bool CATS::GetThetaDependentSource(){
    return ThetaDependentSource;
}

void CATS::SetTransportRenorm(const double& val){
    TransportRenorm = fabs(val);
}
double CATS::GetTransportRenorm(){
    return TransportRenorm;
}

void CATS::SetTotPairMomCut(const double& minval, const double& maxval){
    if(minval<0 || maxval<minval){
        printf("WARNING: Bad input in void CATS::SetTotPairMomCut(const double& minval, const double& maxval)");
        return;
    }
    if(minval==MinTotPairMom && maxval==MaxTotPairMom){
        return;
    }
    //if we set cut conditions outside of the loaded data, we have to reload
    if(UseTotMomCut==true && LoadingComplete &&
       (LoadedMinTotPairMom>minval || LoadedMaxTotPairMom<maxval) ){
        LoadingComplete = false;
    }
    MinTotPairMom = minval;
    MaxTotPairMom = maxval;
    UseTotMomCut = true;
    ComputedCorrFunction = false;
}
void CATS::GetTotPairMomCut(double& minval, double& maxval){
    minval = MinTotPairMom;
    maxval = MaxTotPairMom;
}
void CATS::RemoveTotPairMomCut(){
    MinTotPairMom = -1;
    MaxTotPairMom = 1e100;
    if(TotalPairMomentum){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                delete [] TotalPairMomentum[uMomBin][uIpBin];
            }
            delete [] TotalPairMomentum[uMomBin];
        }
        delete [] TotalPairMomentum; TotalPairMomentum=NULL;
    }
    UseTotMomCut = false;
    if(LoadedMinTotPairMom>0 || LoadedMaxTotPairMom<1e100) LoadingComplete = false;
}

void CATS::SetInputFileName(const char* fname){
    unsigned StrLen = strlen(fname);
    if(!StrLen){
        printf("WARNING: The input file name is empty!\n");
        return;
    }
    if(InputFileName){
        //if this file was already loaded before
        if(strcmp(InputFileName, fname)==0 && LoadingComplete){
            return;
        }
        delete [] InputFileName;
        InputFileName = NULL;
    }

    InputFileName = new char [StrLen+1];
    strcpy(InputFileName, fname);
    LoadingComplete = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}

void CATS::GetInputFileName(char* fname){
    if(!InputFileName){
        strcpy(fname, "");
        return;
    }
    strcpy(fname, InputFileName);
}

double CATS::GetCorrFun(const unsigned& WhichMomBin){
    if(WhichMomBin>=NumMomBins || !CorrFun) return 0;
    return CorrFun[WhichMomBin][NumIpBins];
}

double CATS::GetCorrFun(const unsigned& WhichMomBin, double& Momentum){
    if(WhichMomBin>=NumMomBins || !CorrFun) return 0;
    Momentum = WhichMomBin<NumMomBins?(MomBin[WhichMomBin]+MomBin[WhichMomBin+1])*0.5:0;
    return GetCorrFun(WhichMomBin);
}

//in short: here we want to make linear interpolation. For that we need to find two points.
//In perfect case Momentum should be between the two points. However, in case we are working in the first
//half of the 0th bin, or the last half of the last mom. bin, than we need to take either the two point above or
//the two points below the value of Momentum
double CATS::EvalCorrFun(const double& Momentum){
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins]) return 0;
    if(NumMomBins==1) return CorrFun[0][NumIpBins];
    unsigned WhichMomBin = GetMomBin(Momentum);

    double RelMom[3];
    RelMom[0] = WhichMomBin?GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin<(NumMomBins-1)?GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double** CorrFunRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CorrFunRange = &CorrFun[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CorrFunRange = &CorrFun[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CorrFunRange = &CorrFun[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CorrFunRange = &CorrFun[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrFun[WhichMomBin][NumIpBins];
    }

    return (CorrFunRange[1][NumIpBins]*(Momentum-InterpolRange[0])-
            CorrFunRange[0][NumIpBins]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

double CATS::EvalCorrFunErr(const double& Momentum){
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins]) return 0;
    if(NumMomBins==1) return CorrFunError[0][NumIpBins];
    unsigned WhichMomBin = GetMomBin(Momentum);
    double RelMom[3];
    RelMom[0] = WhichMomBin?GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin!=(NumMomBins-1)?GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double** CorrFunRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CorrFunRange = &CorrFunError[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CorrFunRange = &CorrFunError[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CorrFunRange = &CorrFunError[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CorrFunRange = &CorrFunError[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrFunError[WhichMomBin][NumIpBins];
    }

    return (CorrFunRange[1][NumIpBins]*(Momentum-InterpolRange[0])-
            CorrFunRange[0][NumIpBins]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

double CATS::GetCorrFunErr(const unsigned& WhichMomBin){
    if(WhichMomBin>=NumMomBins || !CorrFunError) return 0;
    return CorrFunError[WhichMomBin][NumIpBins];
}

double CATS::GetCorrFunErr(const unsigned& WhichMomBin, double& Momentum){
    if(WhichMomBin>=NumMomBins || !CorrFunError) return 0;
    Momentum = WhichMomBin<NumMomBins?(MomBin[WhichMomBin]+MomBin[WhichMomBin+1])*0.5:0;
    return GetCorrFunErr(WhichMomBin);
}



double CATS::GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin){
    if(WhichMomBin>=NumMomBins || WhichIpBin>=NumIpBins || !CorrFun) return 0;
    return CorrFun[WhichMomBin][WhichIpBin];
}

double CATS::GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin, double& Momentum, double& ImpPar){
    if(NumMomBins<=WhichMomBin || NumIpBins<=WhichIpBin || !CorrFun) return 0;
    Momentum = WhichMomBin<NumMomBins?(MomBin[WhichMomBin]+MomBin[WhichMomBin+1])*0.5:0;
    ImpPar = WhichIpBin<NumIpBins?(IpBin[WhichIpBin]+IpBin[WhichIpBin+1])*0.5:0;
    return GetCorrFunIp(WhichMomBin, WhichIpBin);
}

double CATS::GetPhaseShift(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW){
    if(NumMomBins<=WhichMomBin || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    return PhaseShift[WhichMomBin][usCh][usPW];
}

double CATS::EvalPhaseShift(const double& Momentum, const unsigned short& usCh, const unsigned short& usPW){
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins] || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    unsigned WhichMomBin = GetMomBin(Momentum);
    double RelMom[3];
    RelMom[0] = WhichMomBin?GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin!=(NumMomBins-1)?GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double*** PhaseShiftRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        PhaseShiftRange = &PhaseShift[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        PhaseShiftRange = &PhaseShift[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        PhaseShiftRange = &PhaseShift[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        PhaseShiftRange = &PhaseShift[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return PhaseShift[WhichMomBin][usCh][usPW];
    }

    return (PhaseShiftRange[1][usCh][usPW]*(Momentum-InterpolRange[0])-
            PhaseShiftRange[0][usCh][usPW]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

double CATS::GetMomentum(const unsigned& WhichMomBin){
    if(NumMomBins<=WhichMomBin) return 0;
    return 0.5*(MomBin[WhichMomBin]+MomBin[WhichMomBin+1]);
}

double CATS::GetMomBinLowEdge(const unsigned& WhichMomBin){
    if(NumMomBins<WhichMomBin) return 0;
    return MomBin[WhichMomBin];
}

double CATS::GetMomBinUpEdge(const unsigned& WhichMomBin){
    if(NumMomBins<=WhichMomBin) return 0;
    return MomBin[WhichMomBin+1];
}

const double& CATS::FmNu(){
    return FmToNu;
}

const double& CATS::NuFm(){
    return NuToFm;
}

void CATS::RemoveShortRangePotential(){
    if(!ShortRangePotential) return;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
            ShortRangePotential[usCh][usPW] = 0;
            PotPar[usCh][usPW] = 0;
        }
    }
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::RemoveShortRangePotential(const unsigned& usCh, const unsigned& usPW){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::RemoveShortRangePotential(...)\n");
        return;
    }
    if(usPW>=NumPW[usCh]){
        printf("WARNING: Bad input in CATS::RemoveShortRangePotential(...)\n");
        return;
    }
    if(!ShortRangePotential) return;
    if(!ShortRangePotential[usCh]) return;
    ShortRangePotential[usCh][usPW] = 0;
    PotPar[usCh][usPW] = 0;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::SetShortRangePotential(const unsigned& usCh, const unsigned& usPW,
                           double (*pot)(double* Pars), double* Pars){
    if(usCh>=NumCh){
        printf("WARNING: Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }
    if(usPW>=NumPW[usCh]){
        printf("WARNING: Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }
    if(ShortRangePotential[usCh][usPW]==pot && PotPar[usCh][usPW]==Pars){
        return;
    }

    ShortRangePotential[usCh][usPW] = pot;
    PotPar[usCh][usPW] = Pars;

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

}

void CATS::RemoveAnaSource(){
    if(!AnaSourcePar) return;
    AnalyticSource = NULL;
    ComputedCorrFunction = false;
}
void CATS::SetAnaSource(double (*AS)(double*), double* Pars){
    if(AnalyticSource==AS && AnaSourcePar==Pars) return;
    AnalyticSource = AS;
    AnaSourcePar = Pars;
    ComputedCorrFunction = false;
}

//!Running the analysis
void CATS::KillTheCat(const int& Options){

    //Check if all needed variables are defined
    if(!UseAnalyticSource && !InputFileName){
        printf("\033[1;31mERROR!\033[0m The data cannot be loaded! The path is not set!\n\n");
        return;
    }
    //if(UseAnalyticSource && !AnalyticSourceRad && !AnalyticSourceRadCosTh){
    if(UseAnalyticSource && !AnalyticSource){
        printf("\033[1;31mERROR!\033[0m The analytic source function is not set!\n\n");
        return;
    }
    if(!pdgID[0] || !pdgID[1]){
        printf("\033[1;31mERROR!\033[0m The PDG IDs of the particles are not set!\n\n");
        return;
    }
    if(!MomBin){
        printf("\033[1;31mERROR!\033[0m The momentum bins are not set!\n\n");
        return;
    }
    if(!IpBin){
        printf("\033[1;31mERROR!\033[0m The impact parameter bins are not set!\n\n");
        return;
    }
    if(!NumCh){
        printf("\033[1;31mERROR!\033[0m There is not a single spin-state set!\n\n");
        return;
    }

    printf("\033[1;37mKilling the cat!\033[0m\n----------------\n");

    switch(Options){
    case kSourceChanged:
        LoadingComplete *= UseAnalyticSource;
        ComputedCorrFunction = false;
        break;
    case kPotentialChanged:
        ComputedWaveFunction = false;
        ComputedCorrFunction = false;
        break;
    case kAllChanged:
        LoadingComplete *= UseAnalyticSource;
        ComputedWaveFunction = false;
        ComputedCorrFunction = false;
        break;
    default: break;
    }

    //in case we have a cut on the total momentum, but the array to save it is not present
    //than the data needs to be reloaded
    if(UseTotMomCut && !TotalPairMomentum) LoadingComplete=false;

    printf("\033[1;37m Stage 1:\033[0m Obtaining the source...");
    if( (!LoadingComplete) && !UseAnalyticSource ){
        printf(" Loading from the data-file...\n");
        LoadData();
        if(LoadingComplete){
            printf("          \033[1;37mLoading status:\033[1;32m Success!\033[1;37m\n"
                    "          Number of pairs loaded: %u\033[0m\n",
                    NumPairs);
        }
        else{
            printf("          \033[1;37mLoading status:\033[1;31m Failed!\033[0m\n\n");
            return;
        }
    }
    else{
        printf("\n");
    }
    if(UseAnalyticSource){
        printf("          Using analytic source.\n");
    }
    else{
        printf("          \033[0mUsing data-defined source.\n");
    }

    printf("\033[1;37m Stage 2:\033[0m Computing the wave-function...\n");
    if(!ComputedWaveFunction) ComputeWaveFunction();
    printf("          \033[1;32mDone!\033[0m\n");

    printf("\033[1;37m Stage 3:\033[0m Computing the correlation function...\n");
    if(!ComputedCorrFunction) FoldSourceAndWF();
    printf("          \033[1;32mDone!\033[0m\n");

    printf("\n");
}

void CATS::ComputeWaveFunction(){
    //Reserve memory for the output
    if(!SavedWaveFunBins){
        SavedWaveFunBins = new unsigned** [NumMomBins];
        PhaseShift = new double** [NumMomBins];
        WaveFunRad = new double*** [NumMomBins];
        WaveFunctionU = new double*** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            SavedWaveFunBins[uMomBin] = new unsigned* [NumCh];
            PhaseShift[uMomBin] = new double* [NumCh];
            WaveFunRad[uMomBin] = new double** [NumCh];
            WaveFunctionU[uMomBin] = new double** [NumCh];
            for(unsigned short usCh=0; usCh<NumCh; usCh++){
                SavedWaveFunBins[uMomBin][usCh] = new unsigned [NumPW[usCh]];
                PhaseShift[uMomBin][usCh] = new double [NumPW[usCh]];
                WaveFunRad[uMomBin][usCh] = new double* [NumPW[usCh]];
                WaveFunctionU[uMomBin][usCh] = new double* [NumPW[usCh]];
                for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
                    SavedWaveFunBins[uMomBin][usCh][usPW] = 0;
                    PhaseShift[uMomBin][usCh][usPW] = 0;
                    WaveFunRad[uMomBin][usCh][usPW] = NULL;
                    WaveFunctionU[uMomBin][usCh][usPW] = NULL;
                }
            }
        }
    }

    if(!MomBinConverged){
        MomBinConverged = new bool [NumMomBins];
    }
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        MomBinConverged[uMomBin] = true;
    }

    //here one finds the max num PWs, this is needed later on for
    //the correct mapping of all variables
    unsigned short MaxNumPW=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        if(NumPW[usCh]>MaxNumPW){
            MaxNumPW = NumPW[usCh];
        }
    }
    unsigned TotalNumberOfBins = NumMomBins*NumCh*MaxNumPW;


    unsigned uMomBin;
    unsigned short usCh;
    unsigned short usPW;

    for(unsigned uMPP=0; uMPP<TotalNumberOfBins; uMPP++){
        //compute to which MomBin, Polarization and PW corresponds this MPP,
        //i.e. map uMomBin, usCh and usPW to uMPP
        uMomBin = uMPP/(NumCh*MaxNumPW);
        usCh = (uMPP%(NumCh*MaxNumPW))/(MaxNumPW);
        usPW = (uMPP%(NumCh*MaxNumPW))%(MaxNumPW);
        //since uMPP is build under the assumption all NumPW are the same
        //one has to check if usPW has a meaningful value!
        if(usPW>=NumPW[usCh]) continue;
        //if the potential for this channel and PW is zero => continue
        if(!ShortRangePotential[usCh][usPW]) continue;
        //skip momentum bins which have obtained an error code
        if(!MomBinConverged[uMomBin] && ExcludeFailedConvergence) continue;

        SavedWaveFunBins[uMomBin][usCh][usPW]=0;

        //if s+l is odd, than this partial wave will cancel out during the
        //symmetrization for identical particles
        if( IdenticalParticles && (usPW+usCh)%2 ) continue;

        double Momentum;
        //the momentum is taken from the center of the bin
        Momentum = 0.5*(MomBin[uMomBin]+MomBin[uMomBin+1]);

        double* BufferWaveFunction;
        double* BufferRad;

        unsigned NumComputedPoints = 2;//counting the initial two starting points
        //index that keeps track of which value in the arrays correspond to the previous, current and next computing step.
        short kOld;
        short kCurrent;
        short kNew;
        double WaveFun[3];
        double PosRad[3];
        double Rho[3];
        //value of the propagating function
        double PropFunVal[3];
        //the value of the prop. function without strong interaction.
        //in case it is equal to the PropFunVal, than the algorithm must have converged
        double PropFunWithoutSI[3];
        //Step size
        double DeltaRad[3];
        //Step size^2
        double DeltaRad2[3];

        double MaxDeltaRad;
        double MinDeltaRad;

        PropagatingFunction(PropFunWithoutSI[0], PropFunVal[0], StartRad, Momentum, usPW, usCh);
        MinDeltaRad = sqrt(fabs(EpsilonProp/(PropFunVal[0]+1e-64)));
        MaxDeltaRad = sqrt(EpsilonProp/(Momentum*Momentum));

        if(MinDeltaRad>MaxDeltaRad) MinDeltaRad=MaxDeltaRad;

        kOld=0;
        kCurrent=1;
        kNew=2;

        //set the initial step in r
        //!This definition of DeltaRad is used in other functions as well.
        //!It is vital that in case there is a need to change the definition, it is changed
        //!in all other functions as well!
        //DeltaRad[kOld] = RhoStep/Momentum;
        //DeltaRad2 = DeltaRad*DeltaRad;

        PosRad[kOld] = 0;
        Rho[kOld] = Momentum*PosRad[kOld];
        DeltaRad[kOld] = MinDeltaRad;
        DeltaRad2[kOld] = DeltaRad[kOld]*DeltaRad[kOld];
        PropagatingFunction(PropFunWithoutSI[kOld], PropFunVal[kOld], PosRad[kOld], Momentum, usPW, usCh);

        PosRad[kCurrent] = PosRad[kOld]+DeltaRad[kOld];
        Rho[kCurrent] = Momentum*PosRad[kCurrent];
        DeltaRad[kCurrent] = MinDeltaRad;
        DeltaRad2[kCurrent] = DeltaRad[kCurrent]*DeltaRad[kCurrent];
        PropagatingFunction(PropFunWithoutSI[kCurrent], PropFunVal[kCurrent], PosRad[kCurrent], Momentum, usPW, usCh);

        //the initial values for the wave function are set based on the solution without the strong potential.
        //this will of course lead to a wrong normalization in the asymptotic region, but this will be corrected for later on,
        //at this stage it is only important that the algorithm gets a meaningful guess so that we do not encounter
        //overflow problems during the calculation.
        WaveFun[kOld] = ReferencePartialWave(PosRad[kOld], Momentum, usPW);
        WaveFun[kCurrent] = ReferencePartialWave(PosRad[kCurrent], Momentum, usPW);

        bool Convergence = false;
        bool Converged = false;
        //at which point the convergence criteria first occurred
        double ConvergenceRadius=0;
        //how much after the convergence was the WF propagated
        //(this is needed for the normalization)
        const double ConvIntervalLength=3.14;
        const double ConvIntervalRadLength=ConvIntervalLength/Momentum;
        unsigned ConvIntervalSteps=0;

        unsigned StepOfMaxConvergedNumWF=0;
        double MaxConvergedNumWF=0;
        double MaxConvRho=0;
        //the last radius at which the computation was performed
        double MaxConvRad=0;
        double DeltaRadAtMaxConvPrev;
        double DeltaRadAtMaxConv;
        double ConvPointOldWeight=-1;
        double ConvPointWeight=0;

        const unsigned short MinConvSteps = 32;

        //the next two variables are estimates of the max. required bin numbers in case one fails to converge
        //MaxNumRadSteps assumes the worst case scenario -> we compute always with the worst possible step length
        //the second one is the optimistic -> we always compute with the best possible step length
        //bins that converge have usually less than MinNumRadSteps entries, bins that fail to converge a bit above.
        //EstNumRadSteps is a fair small over-estimation of the actual bins needed. It is computed by assuming we have the worst-case
        //scenario up until rho=1.57 or r=4 fm (whichever occurs first at this Momentum) and assumes that the rest as best case scenario.
        //This number looks to be c.a. 10x bigger than the actual num. bins, so use it for memory allocation.

        unsigned EstNumRadSteps = ceil(( (MaxRad*Momentum)>MaxRho?MaxRad-1.57/Momentum:MaxRho/Momentum+ConvIntervalRadLength-1.57/Momentum)/MaxDeltaRad)+MinConvSteps+1 +
                                    ceil(( 4.*FmToNu<1.57/Momentum?4.*FmToNu:1.57/Momentum )/MinDeltaRad);

        BufferWaveFunction = new double [EstNumRadSteps];
        BufferRad = new double [EstNumRadSteps];

        BufferWaveFunction[0] = WaveFun[kOld];
        BufferWaveFunction[1] = WaveFun[kCurrent];

        BufferRad[0] = PosRad[kOld];
        BufferRad[1] = PosRad[kCurrent];

        //!this is the main loop, which computes each next WF point until the result converges or
        //!a maximum value of rho is reached.
        //N.B. if the result is currently converging, the numerical method should not be interrupted even
        //if the MaxRad condition is met!
        while( (!Converged||(PosRad[kCurrent]<ConvergenceRadius+ConvIntervalRadLength))
              && (PosRad[kOld]<MaxRad || Rho[kOld]<MaxRho || Converged) ){

            PosRad[kNew] = PosRad[kCurrent] + DeltaRad[kCurrent];
            Rho[kNew] = PosRad[kNew]*Momentum;

            WaveFun[kNew] = WaveFun[kCurrent]*(1.+DeltaRad[kCurrent]/DeltaRad[kOld]) -
                            WaveFun[kOld]*DeltaRad[kCurrent]/DeltaRad[kOld] +
                            PropFunVal[kCurrent]*WaveFun[kCurrent]*DeltaRad2[kCurrent];

            //if we run out of memory...
            if(NumComputedPoints>=EstNumRadSteps){
                EstNumRadSteps*=2;

                double* BufferTempWF = new double [EstNumRadSteps];
                for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                    BufferTempWF[uTmp] = BufferWaveFunction[uTmp];
                }
                delete [] BufferWaveFunction;
                BufferWaveFunction = BufferTempWF;

                double* BufferTempRad = new double [EstNumRadSteps];
                for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                    BufferTempRad[uTmp] = BufferRad[uTmp];
                }
                delete [] BufferRad;
                BufferRad = BufferTempRad;

            }

            BufferWaveFunction[NumComputedPoints] = WaveFun[kNew];

            PropagatingFunction(PropFunWithoutSI[kNew], PropFunVal[kNew], PosRad[kNew], Momentum, usPW, usCh);

            DeltaRad2[kNew] = EpsilonProp/(fabs(PropFunVal[kNew])+1e-64);
            DeltaRad[kNew] = sqrt(DeltaRad2[kNew]);
            if(DeltaRad[kNew]<MinDeltaRad){
                DeltaRad[kNew]=MinDeltaRad;
                DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
            }
            else if(DeltaRad[kNew]>MaxDeltaRad){
                DeltaRad[kNew]=MaxDeltaRad;
                DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
            }

            BufferRad[NumComputedPoints] = PosRad[kNew];

            Convergence =   fabs( (PropFunWithoutSI[kOld]-PropFunVal[kOld])/(fabs(PropFunWithoutSI[kOld]+PropFunVal[kOld])+1e-64) ) < EpsilonConv &&
                            fabs( (PropFunWithoutSI[kCurrent]-PropFunVal[kCurrent])/(fabs(PropFunWithoutSI[kCurrent]+PropFunVal[kCurrent])+1e-64) ) < EpsilonConv &&
                            fabs( (PropFunWithoutSI[kNew]-PropFunVal[kNew])/(fabs(PropFunWithoutSI[kNew]+PropFunVal[kNew])+1e-64) ) < EpsilonConv;

            //in case we have detected a Convergence
            //1) make sure the convergence is real and not a local artifact of the potential
            //2) reset all variables used to characterize the convergence region
            if( Convergence && !Converged){
                Converged = true;
                ConvergenceRadius = PosRad[kNew];

                //2):
                MaxConvergedNumWF = 0;
                MaxConvRho = 0;
                MaxConvRad = 0;
                StepOfMaxConvergedNumWF = 0;
                DeltaRadAtMaxConvPrev = MinDeltaRad;
                DeltaRadAtMaxConv = MinDeltaRad;
            }
            //this condition makes sure that 1) is fulfilled.
            else{
                Converged = Convergence;
            }

            //this will reset the ConvIntervalSteps in case the convergence is reset
            ConvIntervalSteps *= Converged;
            //this will count the number of steps computed after convergence
            ConvIntervalSteps += Converged;

            ConvPointWeight = (PosRad[kNew]-ConvergenceRadius)/ConvIntervalRadLength;

            if(fabs(MaxConvergedNumWF)*ConvPointOldWeight<=fabs(WaveFun[kNew])*ConvPointWeight){
                MaxConvergedNumWF = WaveFun[kNew];
                MaxConvRho = Rho[kNew];
                MaxConvRad = PosRad[kNew];
                StepOfMaxConvergedNumWF = NumComputedPoints;
                DeltaRadAtMaxConvPrev = DeltaRad[kCurrent];
                DeltaRadAtMaxConv = DeltaRad[kNew];
                ConvPointOldWeight = ConvPointWeight;
            }

            NumComputedPoints++;

            kNew++;
            kNew=kNew%3;
            kCurrent++;
            kCurrent=kCurrent%3;
            kOld++;
            kOld=kOld%3;
        }//while(!Converged && PosRad[kOld]<MaxRad)

        //this is not desired for the later computation
        if(StepOfMaxConvergedNumWF==NumComputedPoints-1){
            StepOfMaxConvergedNumWF--;
            MaxConvergedNumWF = BufferWaveFunction[StepOfMaxConvergedNumWF];
            MaxConvRho -= DeltaRadAtMaxConvPrev*Momentum;
            MaxConvRad -= DeltaRadAtMaxConvPrev;
            DeltaRadAtMaxConv = DeltaRadAtMaxConvPrev;
        }

//! SOMETIMES AT HIGH MOMENTA THE RESULT SIMPLY DOES NOT CONVERGE TO THE REQUIRED VALUE
//за това дай на юзера опцията да реши какво да прави с неконвергирали резултати
//1) изхвърли ги
//2) запаза ги, като нормировката я направи на база максимума в MaxRho-3.14 (btw. make 3.14 the min. value for MaxRho!!!)

        //in case the method failed to converge, the whole momentum bin is marked as unreliable and no
        //further computations are done. This point will be excluded from the final result
        if(!Converged){
            MomBinConverged[uMomBin] = false;
            printf("          \033[1;33mWARNING:\033[0m A momentum bin at %.2f MeV failed to converge!\n", Momentum);
            if(ExcludeFailedConvergence){
                printf("                   It will be excluded from the final result!\n");
            }
        }

        //if the maximum wave-function is zero the computation will fail!
        //By design this should not really happen.
        if(!MaxConvergedNumWF){
            printf("WARNING: MaxConvergedNumWF is zero, which is not allowed and points to a bug in the code!\n");
            printf("         Please contact the developers and do not trust your current results!\n");
        }

        //!Now follows the normalization of the numerical wave function to the asymptotic solution
        //the individual steps are explained in detail in the official CATS documentation
        if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence){
            double ShiftRadStepLen = 0.1/Momentum;
            const unsigned MaximumShiftIter = 121;

            double ShiftRad;

            double DownShift=0;
            double UpShift=ShiftRadStepLen;
            if(!BufferWaveFunction[StepOfMaxConvergedNumWF] || !BufferWaveFunction[StepOfMaxConvergedNumWF+1]){
                printf("WARNING: BufferWaveFunction is zero, which is not allowed and points to a bug in the code!\n");
                printf("         Please contact the developers and do not trust your current results!\n");
            }
            double NumRatio = BufferWaveFunction[StepOfMaxConvergedNumWF+1]/BufferWaveFunction[StepOfMaxConvergedNumWF];
            double DownValue;
            double UpValue;
            double SignProduct;
            double SignRefAtLimits;

            //CHECK IF ZERO IS A VIABLE OPTION!
            DownShift = -0.5*ShiftRadStepLen;
            if(fabs(DownShift)>MaxConvRad) DownShift = -MaxConvRad;
            UpShift = 0.5*ShiftRadStepLen;

            //a potential solution should be locked in a region where the ShiftedReferenceWave*NumericalSol have the same sign.
            //furthermore the ShiftedReferenceWave cannot possibly change sign in this region
            for(unsigned uShift=0; uShift<MaximumShiftIter; uShift++){
                //positive side
                CurrentRhoStep = DeltaRadAtMaxConv*Momentum;
                DownValue = AsymptoticRatio(MaxConvRad+DownShift, Momentum, usPW) - NumRatio;
                UpValue = AsymptoticRatio(MaxConvRad+UpShift, Momentum, usPW) - NumRatio;


                if(DownValue*UpValue<0){
                    ShiftRad = 0.5*(DownShift+UpShift);
                    SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW)*MaxConvergedNumWF;
                    SignRefAtLimits = ReferencePartialWave(MaxConvRad+DownShift, Momentum, usPW)*ReferencePartialWave(MaxConvRad+UpShift, Momentum, usPW);
                    if(SignProduct>0 && SignRefAtLimits>0) break;
                }

                //negative side
                if(UpShift<=MaxConvRad){
                    DownValue = AsymptoticRatio(MaxConvRad-DownShift, Momentum, usPW) - NumRatio;
                    UpValue = AsymptoticRatio(MaxConvRad-UpShift, Momentum, usPW) - NumRatio;
                    if(DownValue*UpValue<0){
                        ShiftRad = -0.5*(DownShift+UpShift);
                        SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW)*MaxConvergedNumWF;
                        SignRefAtLimits = ReferencePartialWave(MaxConvRad-DownShift, Momentum, usPW)*ReferencePartialWave(MaxConvRad-UpShift, Momentum, usPW);
                        if(SignProduct>0 && SignRefAtLimits>0){
                            double Temp = DownShift;
                            DownShift = -UpShift;
                            UpShift = -Temp;
                            break;
                        }
                    }
                }

                //the factor 0.95 is there just to make sure that we don't by accident
                //miss a zero just around the limiting values
                UpShift += 0.95*ShiftRadStepLen;
                DownShift = UpShift-ShiftRadStepLen;
            }

            ShiftRad = NewtonRapson(&CATS::AsymptoticRatio,
                    DeltaRadAtMaxConv, usPW, Momentum, MaxConvRad+DownShift, MaxConvRad+UpShift, NumRatio) - MaxConvRad;

            double ShiftRho = ShiftRad*Momentum;
            double Norm = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW)/MaxConvergedNumWF;

            PhaseShift[uMomBin][usCh][usPW] = ShiftRho;

            //we save all entries up to the point in which the normalization was performed. For higher values
            //we will later on use the asymptotic form. N.B. in principle one could save a bit of space and
            //save the result at the first moment of detected convergence, however than in case the convergence
            //was not "perfect" there might be some inaccuracy in the interval up to the normalization point.
            //btw. the step size is set to be the MaxDeltaRad

            SavedWaveFunBins[uMomBin][usCh][usPW] = StepOfMaxConvergedNumWF+1;

            if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
            WaveFunRad[uMomBin][usCh][usPW] = new double [SavedWaveFunBins[uMomBin][usCh][usPW]];

            if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
            WaveFunctionU[uMomBin][usCh][usPW] = new double [SavedWaveFunBins[uMomBin][usCh][usPW]];

            for(unsigned uPoint=0; uPoint<SavedWaveFunBins[uMomBin][usCh][usPW]; uPoint++){
                WaveFunRad[uMomBin][usCh][usPW][uPoint] = BufferRad[uPoint];
                WaveFunctionU[uMomBin][usCh][usPW][uPoint] = Norm*BufferWaveFunction[uPoint];
            }

        }//if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence)

        delete [] BufferRad;
        delete [] BufferWaveFunction;

    }//for(unsigned uMPP=0; uMPP<TotalNumberOfBins; uMPP++)
    ComputedWaveFunction = true;
}

//! N.B. the units in this function (until the result is saved) are fm and GeV!!!
void CATS::LoadData(const unsigned short& NumBlankHeaderLines){
    DLM_Timer dlmTimer;
    double Time;

    bool ProgressBar=false;
    LoadingComplete = false;
    char* cdummy = new char [512];

    if(!RelativeMomentum){
        RelativeMomentum = new double** [NumMomBins];
        RelativePosition = new double** [NumMomBins];
        RelativeCosTheta = new double** [NumMomBins];
        if(UseTotMomCut) TotalPairMomentum = new double** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            RelativeMomentum[uMomBin] = new double* [NumIpBins];
            RelativePosition[uMomBin] = new double* [NumIpBins];
            RelativeCosTheta[uMomBin] = new double* [NumIpBins];
            if(UseTotMomCut) TotalPairMomentum[uMomBin] = new double* [NumIpBins];
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                RelativeMomentum[uMomBin][uIpBin] = new double [MaxPairsPerBin];
                RelativePosition[uMomBin][uIpBin] = new double [MaxPairsPerBin];
                RelativeCosTheta[uMomBin][uIpBin] = new double [MaxPairsPerBin];
                if(UseTotMomCut) TotalPairMomentum[uMomBin][uIpBin] = new double [MaxPairsPerBin];
            }
        }
    }

    if(UseTotMomCut && !TotalPairMomentum){
        TotalPairMomentum = new double** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            TotalPairMomentum[uMomBin] = new double* [NumIpBins];
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                TotalPairMomentum[uMomBin][uIpBin] = new double [MaxPairsPerBin];
            }
        }
    }

    if(!LoadedPairsPerBin){
        LoadedPairsPerBin = new unsigned* [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            LoadedPairsPerBin[uMomBin] = new unsigned [NumIpBins];
        }
    }
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            LoadedPairsPerBin[uMomBin][uIpBin] = 0;
        }
    }

    if(!WeightIp) WeightIp = new double [NumIpBins];
    if(!WeightIpError) WeightIpError = new double [NumIpBins];

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName);
        return;
    }

    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;

    //Read the header lines
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        fgets(cdummy, 511, InFile);
    }

    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName);
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    unsigned NumEvents=0;
    NumPairs=0;

    unsigned NumTotalPairs=0;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    int ParticleNr;
    unsigned BufferSize = EventMixing?BufferEventMix:128;
    int** ParticleID = new int*[NumIpBins];
    //[0] = Energy, [1-3] px - pz, units of GeV
    double*** MomCrd = new double**[NumIpBins];
    double** Mass = new double*[NumIpBins];
    //[0] = Time, [1-3] x - z, units of fm and fm/c
    double*** RadCrd = new double**[NumIpBins];
    unsigned* EventID = new unsigned [BufferSize];

    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        MomCrd[uIpBin] = new double*[BufferSize];
        Mass[uIpBin] = new double[BufferSize];
        RadCrd[uIpBin] = new double*[BufferSize];
        ParticleID[uIpBin] = new int[BufferSize];
        for(unsigned iBuff=0; iBuff<BufferSize; iBuff++){
            MomCrd[uIpBin][iBuff] = new double [4];
            RadCrd[uIpBin][iBuff] = new double [4];
        }
    }

    double gamma;
    double betaX;
    double betaY;
    double betaZ;

    double dMomVec[4];

    double TotMom;
    double TotMomVec[4];

    double dRadVec[4];

    unsigned* NumSePairsIp = new unsigned [NumIpBins];
    unsigned TotalNumSePairs=0;

    unsigned* NumEvPart = new unsigned [NumIpBins];
    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        NumEvPart[uIpBin] = 0;
        NumSePairsIp[uIpBin] = 0;
    }

    double RelMomCom;
    double RedMomComMeV;
    double RelPosCom;
    double RelCosTh;
    unsigned WhichMomBin;
    unsigned WhichIpBin;
    unsigned MaxTotPairs = MaxPairsPerBin*NumIpBins*NumMomBins;
//    if(MaxPairsToLoad<MaxTotPairs){
//        MaxTotPairs = MaxPairsToLoad;
//    }

    //this array will be used for some time-corrections, see below
    double MomC[4];
    double MomC2[4];

    //progress
    //percentage of the file read. The reading speed should be more or less constant,
    //so this value can always give an accurate maximum ETA
    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    //most of the times this is a bad estimate to the ETA, since the MaxPairsPerBin can change drastically
    //the amount of pairs that is accepted per second, i.e. when the bins with high-statistics are full
//    float pMaxPairsToLoad;
    //this should also give a realistic ETA, if we define it based on the bin with fewest entries
    //float pMaxPairsPerBin;

    short pTotal;
    short pTotalOld;

    bool bAllBinsAreFull;

    //!---Iteration over all events---
    while(!feof(InFile)){
        if(NumPairs>=MaxTotPairs) break;
//        if(NumPairs>=MaxPairsToLoad) break;
        bAllBinsAreFull = true;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                bAllBinsAreFull *= LoadedPairsPerBin[uMomBin][uIpBin]>=MaxPairsPerBin;
            }
        }
        if(bAllBinsAreFull){
            printf("HI!\n");
            break;
        }
        fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy);
//!We consider the imp par positive!
        ImpPar = fabs(ImpPar);
        WhichIpBin = GetIpBin(ImpPar);
        if(WhichIpBin>=NumIpBins) continue;
        unsigned& NePart = NumEvPart[WhichIpBin];
        if(NumTotalPairs>=MaxPairsToRead) break;
        NumEvents++;
        //NePart counts the number of particles inside the buffer for this ImpPar bin. In case of no event-mixing this
        //should be reset for each event, i.e. the buffer is limited to a single event
        if(!EventMixing) NePart=0;
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            if(NumPairs>=MaxTotPairs) break;
//            if(NumPairs>=MaxPairsToLoad) break;
            fscanf(InFile,"%i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &ParticleNr,&ParticleID[WhichIpBin][NePart],
                    &MomCrd[WhichIpBin][NePart][1],&MomCrd[WhichIpBin][NePart][2],&MomCrd[WhichIpBin][NePart][3],&MomCrd[WhichIpBin][NePart][0],
                    &Mass[WhichIpBin][NePart],
                    &RadCrd[WhichIpBin][NePart][1],&RadCrd[WhichIpBin][NePart][2],&RadCrd[WhichIpBin][NePart][3],&RadCrd[WhichIpBin][NePart][0]);

            if(MomCrd[WhichIpBin][NePart][0]==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            EventID[NePart] = NumEvents;

            if(ParticleID[WhichIpBin][NePart]!=pdgID[0] && ParticleID[WhichIpBin][NePart]!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            //!---Iteration over all particles in the buffer for this impact parameter---
            //for each particle that we find, we loop over all other particles found in the same event.
            //In case we use event-mixing, the loop is extended over the full buffer corresponding to this impact parameter.
            for(unsigned iePart=0; iePart<NePart; iePart++){
                if(NumPairs>=MaxTotPairs) break;
//                if(NumPairs>=MaxPairsToLoad) break;

        bAllBinsAreFull = true;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                bAllBinsAreFull *= LoadedPairsPerBin[uMomBin][uIpBin]>=MaxPairsPerBin;
            }
        }
        if(bAllBinsAreFull){
            printf("HI2!\n");
            break;
        }

                //this is the case in which we have found a suitable particle pair
                if( (ParticleID[WhichIpBin][iePart]==pdgID[0] && ParticleID[WhichIpBin][NePart]==pdgID[1]) ||
                    (ParticleID[WhichIpBin][iePart]==pdgID[1] && ParticleID[WhichIpBin][NePart]==pdgID[0])){

                    NumTotalPairs++;

                    TotMomVec[0] = MomCrd[WhichIpBin][NePart][0]+MomCrd[WhichIpBin][iePart][0];
                    TotMomVec[1] = MomCrd[WhichIpBin][NePart][1]+MomCrd[WhichIpBin][iePart][1];
                    TotMomVec[2] = MomCrd[WhichIpBin][NePart][2]+MomCrd[WhichIpBin][iePart][2];
                    TotMomVec[3] = MomCrd[WhichIpBin][NePart][3]+MomCrd[WhichIpBin][iePart][3];
                    TotMom = sqrt(TotMomVec[1]*TotMomVec[1]+TotMomVec[2]*TotMomVec[2]+TotMomVec[3]*TotMomVec[3]);

                    betaX = TotMomVec[1]/TotMomVec[0];
                    betaY = TotMomVec[2]/TotMomVec[0];
                    betaZ = TotMomVec[3]/TotMomVec[0];
                    //for the total momentum: 1/sqrt(1-p^2/E^2)=1/sqrt(1-beta^2)
                    gamma = 1./sqrt(1.-TotMom*TotMom/(TotMomVec[0]*TotMomVec[0]));

                    CATSboost Booster(gamma, betaX, betaY, betaZ);

                    //MOMENTUM
                    dMomVec[0] = (MomCrd[WhichIpBin][NePart][0]-MomCrd[WhichIpBin][iePart][0]);
                    dMomVec[1] = (MomCrd[WhichIpBin][NePart][1]-MomCrd[WhichIpBin][iePart][1]);
                    dMomVec[2] = (MomCrd[WhichIpBin][NePart][2]-MomCrd[WhichIpBin][iePart][2]);
                    dMomVec[3] = (MomCrd[WhichIpBin][NePart][3]-MomCrd[WhichIpBin][iePart][3]);

                    //convert dMomX,Y,Z to CM
                    Booster.Boost(dMomVec);

                    RelMomCom = sqrt(dMomVec[1]*dMomVec[1]+dMomVec[2]*dMomVec[2]+dMomVec[3]*dMomVec[3]);
                    //convert to MeV and k = 1/2 q
                    RedMomComMeV = 500.*RelMomCom;

                    //SPACE -> convert to nu
                    dRadVec[0] = (RadCrd[WhichIpBin][NePart][0]-RadCrd[WhichIpBin][iePart][0]);
                    dRadVec[1] = (RadCrd[WhichIpBin][NePart][1]-RadCrd[WhichIpBin][iePart][1]);
                    dRadVec[2] = (RadCrd[WhichIpBin][NePart][2]-RadCrd[WhichIpBin][iePart][2]);
                    dRadVec[3] = (RadCrd[WhichIpBin][NePart][3]-RadCrd[WhichIpBin][iePart][3]);

                    //convert dMomX,Y,Z to CM
                    Booster.Boost(dRadVec);

                    dRadVec[0] *= TransportRenorm;
                    dRadVec[1] *= TransportRenorm;
                    dRadVec[2] *= TransportRenorm;
                    dRadVec[3] *= TransportRenorm;

                    //if RadCrd[NePart] was emitted first (in the CM frame)
                    if(dRadVec[0]>0){
                        MomC[0] = MomCrd[WhichIpBin][NePart][0];
                        MomC[1] = MomCrd[WhichIpBin][NePart][1];
                        MomC[2] = MomCrd[WhichIpBin][NePart][2];
                        MomC[3] = MomCrd[WhichIpBin][NePart][3];
                        MomC2[0] = MomCrd[WhichIpBin][iePart][0];
                        MomC2[1] = MomCrd[WhichIpBin][iePart][1];
                        MomC2[2] = MomCrd[WhichIpBin][iePart][2];
                        MomC2[3] = MomCrd[WhichIpBin][iePart][3];
                    }
                    else{
                        MomC[0] = MomCrd[WhichIpBin][iePart][0];
                        MomC[1] = MomCrd[WhichIpBin][iePart][1];
                        MomC[2] = MomCrd[WhichIpBin][iePart][2];
                        MomC[3] = MomCrd[WhichIpBin][iePart][3];
                        MomC2[0] = MomCrd[WhichIpBin][NePart][0];
                        MomC2[1] = MomCrd[WhichIpBin][NePart][1];
                        MomC2[2] = MomCrd[WhichIpBin][NePart][2];
                        MomC2[3] = MomCrd[WhichIpBin][NePart][3];
                    }

                    if(TauCorrection){
                        //here we obtain the 4-momentum of the particle that was emitted first
                        //this is converted to the CM frame of the particle pair
                        Booster.Boost(MomC);
                        Booster.Boost(MomC2);

                        if(!MomC[0]){
                            printf("WARNING! A strange bug occurred (MomC[0]==0), please investigate and/or contact the developers!\n");
                            printf("         The current output might be wrong!\n");
                        }

                        dRadVec[1] += MomC[1]/MomC[0]*dRadVec[0];
                        dRadVec[2] += MomC[2]/MomC[0]*dRadVec[0];
                        dRadVec[3] += MomC[3]/MomC[0]*dRadVec[0];
                    }

                    //convert to NU at the end. Input should be in fm
                    RelPosCom = sqrt(dRadVec[1]*dRadVec[1]+dRadVec[2]*dRadVec[2]+dRadVec[3]*dRadVec[3]);
                    //RelCosTh = dRadVec[3]/RelPosCom;
                    RelCosTh = (dMomVec[1]*dRadVec[1]+dMomVec[2]*dRadVec[2]+dMomVec[3]*dRadVec[3])/
                               (RelMomCom*RelPosCom);

                    bool Selected = true;

                    if(RelPosCom>256 || RelPosCom!=RelPosCom || RelPosCom==0 || RelMomCom==0){
                        Selected = false;
                    }

                    //only save relevant particle pairs
                    if(RedMomComMeV<MomBin[0] || RedMomComMeV>MomBin[NumMomBins]){
                        Selected = false;
                    }

                    else if(Selected){
                        WhichMomBin = GetMomBin(RedMomComMeV);
                        if(WhichMomBin>=NumMomBins) Selected = false;
                    }

                    //check the total pair momentum condition
                    if(UseTotMomCut && (TotMom*1000<MinTotPairMom || TotMom*1000>MaxTotPairMom)){
                        Selected = false;
                    }

                    if(Selected){
                        if(LoadedPairsPerBin[WhichMomBin][WhichIpBin]>=MaxPairsPerBin){
                            Selected = false;
                        }
                    }
                    if(!Selected){
                        continue;
                    }

                    RelativeMomentum[WhichMomBin][WhichIpBin][LoadedPairsPerBin[WhichMomBin][WhichIpBin]] = RedMomComMeV;
                    RelativePosition[WhichMomBin][WhichIpBin][LoadedPairsPerBin[WhichMomBin][WhichIpBin]] = RelPosCom*FmToNu;
                    RelativeCosTheta[WhichMomBin][WhichIpBin][LoadedPairsPerBin[WhichMomBin][WhichIpBin]] = RelCosTh;
                    if(UseTotMomCut) TotalPairMomentum[WhichMomBin][WhichIpBin][LoadedPairsPerBin[WhichMomBin][WhichIpBin]] = TotMom*1000;

                    //counting the number of pairs coming from the same event. This information is needed when
                    //reweighting the correlation function in the different impact parameter bins
                    if(EventID[NePart]==EventID[iePart]){
                        NumSePairsIp[WhichIpBin]++;
                        TotalNumSePairs++;
                    }

                    LoadedPairsPerBin[WhichMomBin][WhichIpBin]++;
                    NumPairs++;

                }//if(...)
            }//for(int iePart=0; iePart<NePart; iePart++)

            NePart++;
            if(NePart==BufferSize){
                NePart=0;
            }

        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)


        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
//        pMaxPairsToLoad = double(NumPairs)/double(MaxTotPairs);
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;
//        ProgressLoad = pMaxPairsToLoad>ProgressLoad?pMaxPairsToLoad:ProgressLoad;
        /*
        pMaxPairsPerBin = 0;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                //select the smallest possible pMaxPairsPerBin
                pTemp = double(LoadedPairsPerBin[uMomBin][uIpBin])/double(MaxPairsPerBin);
                pMaxPairsPerBin = pMaxPairsPerBin>pTemp?pTemp:pMaxPairsPerBin;
            }
        }
        ProgressLoad = pMaxPairsPerBin>ProgressLoad?pMaxPairsPerBin:ProgressLoad;
        */

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            //EtaPerBin = round((1./pMaxPairsPerBin-1.)*Time);
            //EtaToLoad = round((1./pMaxPairsToLoad-1.)*Time);
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }


    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }
    if(!NumPairs){
        printf("\033[1;31m          WARNING:\033[0m There were no pairs loaded! The computation cannot proceed!\n");
    }
    else if(!TotalNumSePairs){
        printf("\033[1;31m          WARNING:\033[0m There were no same-events pairs found! The computation cannot proceed!\n");
    }
    else{
        LoadingComplete = true;
        if(UseTotMomCut){
            LoadedMinTotPairMom = MinTotPairMom;
            LoadedMaxTotPairMom = MaxTotPairMom;
        }
    }

    if(LoadingComplete && EventMixing){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            WeightIp[uIpBin] = TotalNumSePairs?double(NumSePairsIp[uIpBin])/double(TotalNumSePairs):0;
            if(NumIpBins==1){
                WeightIpError[uIpBin] = 0;
            }
            else if(NumSePairsIp[uIpBin]){
                WeightIpError[uIpBin] = WeightIp[uIpBin]/sqrt(double(NumSePairsIp[uIpBin]));
            }
            else{
                WeightIpError[uIpBin] = 1000;
            }
        }
    }
    else if(!EventMixing){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            WeightIp[uIpBin] = 1./double(NumIpBins);
            WeightIpError[uIpBin] = 0;
        }
    }

    fclose(InFile);
/*
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        delete [] pIpMomBin[uMomBin];
    }
    delete [] pIpMomBin;
*/
    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        for(unsigned iBuff=0; iBuff<BufferSize; iBuff++){
            delete [] MomCrd[uIpBin][iBuff];
            delete [] RadCrd[uIpBin][iBuff];
        }
        delete [] MomCrd[uIpBin];
        delete [] RadCrd[uIpBin];
        delete [] Mass[uIpBin];
        delete [] ParticleID[uIpBin];
    }

    delete [] MomCrd;
    delete [] RadCrd;

    delete [] ParticleID;
    delete [] Mass;
    delete [] cdummy;

    delete [] NumEvPart;
    delete [] EventID;

}

void CATS::FoldSourceAndWF(){
    if(!CorrFun){
        CorrFun = new double* [NumMomBins];
        CorrFunError = new double* [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            CorrFun[uMomBin] = new double [NumIpBins+1];
            CorrFunError[uMomBin] = new double [NumIpBins+1];
            for(unsigned uIpBin=0; uIpBin<=NumIpBins; uIpBin++){
                CorrFun[uMomBin][uIpBin] = 0;
                CorrFunError[uMomBin][uIpBin] = 0;
            }
        }
    }
    else if(!ComputedCorrFunction){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<=NumIpBins; uIpBin++){
                CorrFun[uMomBin][uIpBin] = 0;
                CorrFunError[uMomBin][uIpBin] = 0;
            }
        }
    }
    else return;
    UseAnalyticSource?FoldAnaSourceAndWF():FoldDataSourceAndWF();
}

void CATS::FoldAnaSourceAndWF(){
//!OPTIMIZE?
    const unsigned NumRadPts = ThetaDependentSource?256:1024;
    //const unsigned NumRadPts = 1024;
    const double RadStepLen = MaxRad/double(NumRadPts-1)*NuToFm;
    const unsigned NumCosThPts = ThetaDependentSource?64:1;
    const double CosThStepLen = 2./double(NumCosThPts-1);
    double& Momentum = AnaSourcePar[0];
    double& Radius = AnaSourcePar[1];
    double& CosTheta = AnaSourcePar[2];

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        Momentum = 0.5*(MomBin[uMomBin]+MomBin[uMomBin+1]);
        for(unsigned uRad=0; uRad<NumRadPts; uRad++){
            Radius = 0.5*RadStepLen+RadStepLen*double(uRad);
            if(ThetaDependentSource){
                for(unsigned uCT=0; uCT<NumCosThPts; uCT++){
                    CosTheta = CosThStepLen*double(uCT)-1;
                    CorrFun[uMomBin][NumIpBins] += RadStepLen*CosThStepLen*AnalyticSource(AnaSourcePar)*
                    EffectiveFunctionTheta(Radius*FmToNu, Momentum, CosTheta)*0.5;
                }
            }
            else{
                CorrFun[uMomBin][NumIpBins] += RadStepLen*AnalyticSource(AnaSourcePar)*
                EffectiveFunction(Radius*FmToNu, Momentum);
            }

        }
    }
    ComputedCorrFunction = true;
}

void CATS::FoldDataSourceAndWF(){
    double Radius;
    double Momentum;
    double CosTheta;
    double TotWF;

    if(!LoadingComplete){
        printf("ERROR in CATS::FoldDataSourceAndWF(): There is no data loaded!\n");
        return;
    }

    unsigned NumPairsUsed=0;

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        unsigned NumberOfPairs = 0;
        double* Stdev;
        Stdev = new double [NumIpBins+1];
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            NumPairsUsed=0;
            Stdev[uIpBin] = 0;
            CorrFun[uMomBin][uIpBin] = 0;
            if(!MomBinConverged[uMomBin] && ExcludeFailedConvergence){
                CorrFunError[uMomBin][uIpBin] = 0;
                continue;
            }
            for(unsigned uPair=0; uPair<LoadedPairsPerBin[uMomBin][uIpBin]; uPair++){
                if(UseTotMomCut &&
                   (TotalPairMomentum[uMomBin][uIpBin][uPair]<MinTotPairMom || TotalPairMomentum[uMomBin][uIpBin][uPair]>MaxTotPairMom))
                   continue;
                Radius = RelativePosition[uMomBin][uIpBin][uPair];
                Momentum = RelativeMomentum[uMomBin][uIpBin][uPair];
                CosTheta = RelativeCosTheta[uMomBin][uIpBin][uPair];

                TotWF = ThetaDependentSource?   EffectiveFunctionTheta(Radius, Momentum, CosTheta):
                                                EffectiveFunction(Radius, Momentum);
                CorrFun[uMomBin][uIpBin] += TotWF;

                Stdev[uIpBin] += TotWF*TotWF;

                //if the event mixing is switched off, the total correlation function
                //is evaluated separately
                if(!EventMixing){
                    CorrFun[uMomBin][NumIpBins] += TotWF;
                    Stdev[NumIpBins] += TotWF*TotWF;
                }

                NumPairsUsed++;

            }//for(unsigned uPair=0; uPair<LoadedPairsPerMomBin[uMomBin]; uPair++){
            //if there are too few entries, the error is set basically to infinity
            if(NumPairsUsed<8){
                CorrFun[uMomBin][uIpBin] =  NumPairsUsed?
                                            CorrFun[uMomBin][uIpBin]/double(NumPairsUsed):0;
                CorrFunError[uMomBin][uIpBin] = 100;
            }
            else{
                CorrFun[uMomBin][uIpBin] /= double(NumPairsUsed);
                Stdev[uIpBin] /= double(NumPairsUsed);
                Stdev[uIpBin] = sqrt(Stdev[uIpBin]);
                CorrFunError[uMomBin][uIpBin] = Stdev[uIpBin]/sqrt(double(NumPairsUsed));
            }
            NumberOfPairs += NumPairsUsed;
            //in case of event mixing, the total correlation function is calculated based on the Weights
            if(EventMixing){
                CorrFun[uMomBin][NumIpBins] += WeightIp[uIpBin]*CorrFun[uMomBin][uIpBin];
                CorrFunError[uMomBin][NumIpBins] += pow(WeightIpError[uIpBin]*CorrFunError[uMomBin][uIpBin],2) +
                                                    pow(WeightIpError[uIpBin]*CorrFun[uMomBin][uIpBin],2) +
                                                    pow(WeightIp[uIpBin]*CorrFunError[uMomBin][uIpBin],2);
            }

        }//for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){

        if(!MomBinConverged && ExcludeFailedConvergence){
            CorrFun[uMomBin][NumIpBins] = 0;
            CorrFunError[uMomBin][NumIpBins] = 0;
        }
        //if the event mixing is switched off, the total correlation function
        //is evaluated separately
        else if(!EventMixing){
            //if there are too few entries, the error is set basically to infinity
            if(NumberOfPairs<8){
                CorrFun[uMomBin][NumIpBins] = NumberOfPairs?CorrFun[uMomBin][NumIpBins]/double(NumberOfPairs):0;
                CorrFunError[uMomBin][NumIpBins] = 100;
            }
            else{
                CorrFun[uMomBin][NumIpBins] /= double(NumberOfPairs);
                Stdev[NumIpBins] /= double(NumberOfPairs);
                Stdev[NumIpBins] = sqrt(Stdev[NumIpBins]);
                CorrFunError[uMomBin][NumIpBins] = Stdev[NumIpBins]/sqrt(double(NumberOfPairs));
            }
        }
        //if(EventMixing)
        else{
            CorrFunError[uMomBin][NumIpBins] = sqrt(CorrFunError[uMomBin][NumIpBins]);
        }
        delete [] Stdev;
    }//for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
    ComputedCorrFunction = true;
}

double CATS::CoulombPotential(const double& Radius){
    return Q1Q2*AlphaFS/(fabs(Radius)+1e-64);
}
//the differential equation for the Schroedinger equation
void CATS::PropagatingFunction(double& Basic, double& Full,
                                 const double& Radius, const double& Momentum,
                                 const unsigned short& usPW, const unsigned short& usCh){
    //make sure that there is no division by zero by adding 1e-64
    //the Basic result is the Prop.Fun. WITHOUT a short range potential
    Basic = 2*RedMass*CoulombPotential(Radius) + double(usPW)*(double(usPW)+1)/(Radius*Radius+1e-64) - Momentum*Momentum;
    //the Full result is the Prop.Fun. WITH a short range potential
    PotPar[usCh][usPW][0] = Radius*NuToFm; PotPar[usCh][usPW][1] = Momentum;
    //! In principle this should be executed only if ShortRangePotential[usCh][usPW] is defined.
    //! Do note that this function is NEVER called in case this is not true! Make sure that this stays so!
    Full = Basic + 2*RedMass*ShortRangePotential[usCh][usPW](PotPar[usCh][usPW]);
}

double CATS::PlanePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW){
    double Rho = Radius*Momentum;
    //if Rho is zero, the gsl function will not work
    if(!Rho){
        return 0;
    }
    //N.B. gsl_sf_bessel_jl are defined for Rho>0, but in principle the bessel functions are symmetric for even l
    //and anti-symmetric for odd l, this is implemented here.

    return Rho>0?(Radius)*gsl_sf_bessel_jl(usPW,Rho):pow(-1,usPW)*(Radius)*gsl_sf_bessel_jl(usPW,-Rho);
}

double CATS::CoulombPartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW){
    double Eta = RedMass*Q1Q2*AlphaFS/Momentum;
    double Rho = Radius*Momentum;
    double Overflow=0;
    double Result;
    if(Rho==0) return 0;
    gsl_sf_coulomb_wave_F_array (usPW, 1, Eta, fabs(Rho), &Result, &Overflow);
    Result /= Momentum;

    //N.B. gsl_sf_coulomb_wave_F_array are defined for Rho>0, but in principle the Coulomb functions are symmetric for odd l
    //and anti-symmetric for even l. However here I assume that Momentum>0, and if rho is negative so is the Momentum. Thus
    //the final result for u_l should be symmetric for even l and antisymmetric for odd l. This is implemented here.
    if(Rho<0 && usPW%2==1){
        Result = -Result;
    }
    return Result;
}

double CATS::ReferencePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW){
    return Q1Q2 ? CoulombPartialWave(Radius,Momentum,usPW) : PlanePartialWave(Radius,Momentum,usPW);
}

double CATS::AsymptoticRatio(const double& Radius, const double& Momentum, const unsigned short& usPW){
    return ReferencePartialWave(Radius+CurrentRhoStep/Momentum, Momentum, usPW)/(ReferencePartialWave(Radius, Momentum, usPW)+1e-64);
}

double CATS::NewtonRapson(double (CATS::*Function)(const double&, const double&, const unsigned short&),
                          const double& EpsilonX, const unsigned short& usPW, const double& Momentum,
                          const double&  xMin, const double&  xMax, const double& fValShift){

    const unsigned maxIter = 256;

    double DeltaX;
    double xVal=(xMax+xMin)*0.5;
    double fVal;
    double DeltaF;

    for(unsigned iIter=0; iIter<maxIter; iIter++){
        fVal = (this->*Function)(xVal, Momentum, usPW)-fValShift;
        DeltaF = (this->*Function)(xVal+EpsilonX, Momentum, usPW)-fValShift - fVal;
        if(!DeltaF){
            DeltaX=1;
            printf("WARNING: Something is fishy with the NewtonRapson solver! Might be a bug! Please contact the developers!\n");
        }
        else DeltaX = -fVal*EpsilonX/DeltaF;
        xVal += DeltaX;

        int counter = 0;
        double fValNew = (this->*Function)(xVal, Momentum, usPW)-fValShift;
        while( ( (fabs(fValNew)>fabs(fVal) && counter<16)
              || (xVal<xMin || xVal>xMax) ) ){
            counter++;
            xVal -= DeltaX*pow(2., -counter);
            fValNew = (this->*Function)(xVal, Momentum, usPW)-fValShift;
        }
        if(counter==16 && fabs(fValNew)>fabs(fVal)){
            printf("WARNING: The backtracking of the NewtonRapson root-finder failed!\n");
        }
        if(fabs(fValNew)<fabs(DeltaF)){
            return xVal;
        }
    }
    printf("WARNING: The NewtonRapson root-finder failed!\n");
    return xVal;
}

double CATS::EvalWaveFunctionU(const double& Radius, const double& Momentum,
                                const unsigned short& usCh, const unsigned short& usPW, const bool& DivideByR){
    unsigned uMomBin = GetMomBin(Momentum);
    if(uMomBin>=NumMomBins){
        printf("ERROR: There is a bug inside EvalWaveFunctionU! Contact the developer!");
        return 0;
    }
    double* WFU = WaveFunctionU[uMomBin][usCh][usPW];
    double* WFR = WaveFunRad[uMomBin][usCh][usPW];
    unsigned& SWFB = SavedWaveFunBins[uMomBin][usCh][usPW];
    //double& MaxSavedRad = WFR[SWFB-1];

    unsigned RadBin = GetRadBin(Radius, uMomBin, usCh, usPW);
    double MultFactor = DivideByR?1./(Radius+1e-64):1;

    double DeltaWFR = WFR[RadBin+1]-WFR[RadBin];
    if(!DeltaWFR){
        DeltaWFR = 1-64;
        printf("WARNING: DeltaWFR==0, which might point to a bug! Please contact the developers!\n");
    }

    //make a linear extrapolation
    if(RadBin<SWFB-1){
        return WFU[RadBin]*MultFactor+(WFU[RadBin+1]*MultFactor-WFU[RadBin]*MultFactor)*(Radius-WFR[RadBin])/DeltaWFR;
    }
    else{
        return ReferencePartialWave(Radius+PhaseShift[uMomBin][usCh][usPW]/Momentum, Momentum, usPW);
    }
}

double CATS::EffectiveFunction(const double& Radius, const double& Momentum, const unsigned short& usCh){
    double Result;
    double OldResult=100;
    double TotalResult=0;

    for(unsigned short usPW=0; usPW<1000; usPW++){
        //wave function symmetrization
        if( IdenticalParticles && (usPW+usCh)%2 ) continue;
        //numerical solution, no computation result for zero potential
        if(usPW<NumPW[usCh] && ShortRangePotential[usCh][usPW]){
            Result = EvalWaveFunctionU(Radius, Momentum, usCh, usPW, true);
            TotalResult += double(2*usPW+1)*Result*Result;
        }
        else{
            Result = ReferencePartialWave(Radius, Momentum, usPW)/(Radius+1e-64);
            Result = double(2*usPW+1)*Result*Result;
            TotalResult += Result;
            //convergence criteria
            if(usPW>=NumPW[usCh] && fabs(OldResult)<1e-7 && fabs(Result)<1e-8) break;
            OldResult = Result;
        }
    }

    return TotalResult*(1+IdenticalParticles);
}

double CATS::EffectiveFunction(const double& Radius, const double& Momentum){
    double TotWF=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        TotWF += EffectiveFunction(Radius, Momentum, usCh)*ChannelWeight[usCh];
    }
    return TotWF;
}

double CATS::EffectiveFunctionTheta(const double& Radius, const double& Momentum, const double& CosTheta, const unsigned short& usCh){
    double Result;
    double OldResult=100;
    double TotalResult=0;

    for(unsigned short usPW=0; usPW<1000; usPW++){
        //wave function symmetrization
        if( IdenticalParticles && (usPW+usCh)%2 ) continue;
        //numerical solution, no computation result for zero potential
        if(usPW<NumPW[usCh] && ShortRangePotential[usCh][usPW]){
            Result = double(2*usPW+1)*EvalWaveFunctionU(Radius, Momentum, usCh, usPW, true)*gsl_sf_legendre_Pl(usPW,CosTheta);
            TotalResult += Result*Result;
        }
        else{
            Result = double(2*usPW+1)*ReferencePartialWave(Radius, Momentum, usPW)/(Radius+1e-64)*gsl_sf_legendre_Pl(usPW,CosTheta);
            Result = Result*Result;
            TotalResult += Result;
            if(usPW>=NumPW[usCh] && fabs(OldResult)<1e-7 && fabs(Result)<1e-8) break;
            OldResult = Result;
        }
    }

    return TotalResult*(1+IdenticalParticles);
}

double CATS::EffectiveFunctionTheta(const double& Radius, const double& Momentum, const double& CosTheta){
    double TotWF=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        TotWF += EffectiveFunctionTheta(Radius, Momentum, CosTheta, usCh)*ChannelWeight[usCh];
    }
    return TotWF;
}

unsigned CATS::GetBin(const double& Value, const double* Range, const unsigned& NumBins){
    if(NumBins<=1) return 0;
    unsigned WhichBin=NumBins/2;
    unsigned BinMod=4;
    unsigned BinStep;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(Value<Range[0]) return NumBins;
    if(Value>Range[NumBins-1]) return NumBins+1;
    while(true){
        if(Range[WhichBin]<=Value && Range[WhichBin+1]>=Value){
            return WhichBin;
        }
        else if(Value<Range[WhichBin]){
            BinStep = NumBins/BinMod;
            WhichBin -= BinStep?BinStep:1;
            BinMod *= 2;
        }
        else{
            BinStep = NumBins/BinMod;
            WhichBin += BinStep?BinStep:1;
            BinMod *= 2;
        }
    }
}

unsigned CATS::GetMomBin(const double& Momentum){
    return GetBin(Momentum, MomBin, NumMomBins+1);
}
unsigned CATS::GetIpBin(const double& bVal){
    return GetBin(bVal, IpBin, NumIpBins+1);
}
unsigned CATS::GetRadBin(const double& Radius, const unsigned& uMomBin,
                         const unsigned short& usCh, const unsigned short& usPW){
    return GetBin(Radius, WaveFunRad[uMomBin][usCh][usPW], SavedWaveFunBins[uMomBin][usCh][usPW]);
}



CATSboost::CATSboost(double& g, double& bX, double& bY, double& bZ):gamma(g),betaX(bX),betaY(bY),betaZ(bZ){

}

CATSboost::~CATSboost(){

}

void CATSboost::Boost(double* Vec){
    GammaMomBeta = gamma*(Vec[1]*betaX + Vec[2]*betaY + Vec[3]*betaZ);
    GammaVec0 = gamma*Vec[0];
    GammaDevGammaPlusOne = gamma/(gamma+1);

    Vec[1] += betaX*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    Vec[2] += betaY*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    Vec[3] += betaZ*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    Vec[0] = GammaVec0 - GammaMomBeta;
}

void CATSboost::Boost(const double* InVec, double* OutVec){
    GammaMomBeta = gamma*(InVec[1]*betaX + InVec[2]*betaY + InVec[3]*betaZ);
    GammaVec0 = gamma*InVec[0];
    GammaDevGammaPlusOne = gamma/(gamma+1);

    OutVec[1] += betaX*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[2] += betaY*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[3] += betaZ*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[0] = GammaVec0 - GammaMomBeta;
}
