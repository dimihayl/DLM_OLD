#include "CATStools.h"

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



CATSnode::CATSnode(CATSelder* elder, const short& depth, const unsigned& firstid, const unsigned& lastid, double* mean, double* len):
    Elder(elder),Depth(depth),FirstID(firstid),LastID(lastid){
    MeanVal = NULL;
    IntLen = NULL;
    if(elder!=this) StandardNodeInit(mean, len);
}

CATSnode::~CATSnode(){
    if(MeanVal) {delete [] MeanVal; MeanVal=NULL;}
    if(IntLen) {delete [] IntLen; IntLen=NULL;}

    if(child){
        for(unsigned uSub=0; uSub<Elder->NumSubNodes; uSub++){
            delete child[uSub];
        }
        delete [] child;
        child = NULL;
    }
}

unsigned CATSnode::GetNumOfEl(){
    return LastID-FirstID+1;
}

void CATSnode::StandardNodeInit(double* mean, double* len){
    MeanVal = new double [Elder->Dim];
    IntLen = new double [Elder->Dim];
    double GridSize=1;
    for(short sDim=0; sDim<Elder->Dim; sDim++){
        MeanVal[sDim] = mean[sDim];
        IntLen[sDim] = len[sDim];
        GridSize *= IntLen[sDim];
    }

    if(Elder->SourceFunction){
        for(short sDim=0; sDim<Elder->Dim; sDim++){
            Elder->SourcePars[1+sDim] = MeanVal[sDim];
        }
        SourceValue = Elder->SourceFunction(Elder->SourcePars)*GridSize;
    }
    else if(Elder->GridBoxId){
        unsigned first, last;
        first = Elder->FindFirstParticleWithID(FirstID);
        last = Elder->FindLastParticleWithID(LastID);

        //the normal situation
        if(first<Elder->NumOfEl && last<Elder->NumOfEl){
            SourceValue = double(last-first+1)/double(Elder->NumOfEl);
        }
        //the case where we should include all particles
        //(both boxes are outside the range, first on the low side and last on the upper side)
        else if(first==Elder->NumOfEl && last==first+1){
            SourceValue = 1;
        }
        //the case where the first box is outside range on the low side
        else if(first==Elder->NumOfEl && last<Elder->NumOfEl){
            SourceValue = double(last+1)/double(Elder->NumOfEl);
        }
        //the case where the last box is outside range on the up side
        else if(first<Elder->NumOfEl && last>Elder->NumOfEl){
            SourceValue = double(Elder->NumOfEl-first)/double(Elder->NumOfEl);
        }
        else{
            SourceValue = 0;
        }
    }
    else{
       SourceValue = 0;
    }

    child = NULL;

    if( Depth<Elder->MaxDepth && (Depth<Elder->MinDepth || SourceValue>Elder->Epsilon) ){
        child = new CATSnode* [uipow(2,Elder->NumSubNodes)];
        //we want to divide our total interval in two for each parameter on the grid.
        //in order to keep track in which "quadrant" we are, we introduce a very simple counter WhichPart for each
        //of the parameters, that can only take values 0 or 1. Each time WhichPart[x] is increased to 2, than it is set to zero
        //and WhichPart[x+1] is increased, i.e. we continue to iterate over the next parameter.
        char WhichPart[Elder->Dim];
        double ChildMean[Elder->Dim];
        double ChildLen[Elder->Dim];
        unsigned ChildNumBoxes = (LastID-FirstID+1)/Elder->NumSubNodes;
        unsigned ChildFirstID = FirstID;
        for(short sDim=0; sDim<Elder->Dim; sDim++){
            WhichPart[sDim] = 0;
            ChildMean[sDim] = MeanVal[sDim]-IntLen[sDim]*0.25;
            ChildLen[sDim] = IntLen[sDim]*0.5;
        }
        for(unsigned uSub=0; uSub<Elder->NumSubNodes; uSub++){
            child[uSub] = new CATSnode(Elder,Depth+1,ChildFirstID,ChildFirstID+ChildNumBoxes-1,ChildMean,ChildLen);
            ChildFirstID += ChildNumBoxes;
            for(short sDim=0; sDim<Elder->Dim; sDim++){
                WhichPart[sDim] = (WhichPart[sDim]+1)%2;
                ChildMean[sDim] = MeanVal[sDim]-IntLen[sDim]*0.25+0.5*IntLen[sDim]*WhichPart[sDim];
                if(WhichPart[sDim]) break;
            }
        }
    }
    else{
        Elder->AddEndNode(this);
    }

}


//! see what happens if epsilon==0
CATSelder::CATSelder(const short& dim, const short& mindep, const short& maxdep, const double& epsilon, double* mean, double* len,
                     double (*AS)(double*), double* Pars, unsigned* gbid, const unsigned& numel):
    CATSnode(this, 0, 0, uipow(2,dim*maxdep)-1, mean, len),
    Dim(dim),MinDepth(mindep),MaxDepth(maxdep),Epsilon(epsilon),NumSubNodes(uipow(2,Dim)),NumOfEl(numel){

    MaxNumEndNodes = unsigned(1./Epsilon)*NumSubNodes;
    NumEndNodes=0;
    SourceRenormError = 0;
    EndNode = new CATSnode* [MaxNumEndNodes];

    SourceFunction = AS;
    SourcePars = Pars;
    GridBoxId = gbid;
    if( (!SourceFunction && !SourcePars && !GridBoxId) ||
        (SourceFunction && SourcePars && GridBoxId) ||
        (!SourceFunction && SourcePars) || (SourceFunction && !SourcePars) ||
        (!GridBoxId && NumOfEl) || (GridBoxId && !NumOfEl)
       ){
        printf("ERROR! CATSelder say that the input to the constructor makes no sense!\n");
        SourceFunction = NULL;
        SourcePars = NULL;
        GridBoxId = NULL;
    }
    StandardNodeInit(mean, len);

    double Integral=0;
    for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
        Integral += EndNode[uNode]->SourceValue;
    }
    SourceRenormError = Integral?1./Integral:1e64;
    for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
        EndNode[uNode]->SourceValue *= SourceRenormError;
    }

    if(SourceRenormError<1) SourceRenormError = SourceRenormError?1./SourceRenormError:1e64;
    SourceRenormError -= 1;

    if(SourceRenormError>Epsilon*sqrt(double(NumEndNodes))){
        printf("WARNING: CATSelder says that SourceRenormError=%.6f, which seems odd.\n", SourceRenormError);
        printf("         Either the source function is not properly normalized or there is a bug in CATS!\n");
        printf("         In case its the letter, please contact the developers!\n");
    }

}

CATSelder::~CATSelder(){
    delete [] EndNode;
    EndNode = NULL;
}

short CATSelder::GetMaxDepth(){
    return MaxDepth;
}

unsigned CATSelder::GetNumEndNodes(){
    return NumEndNodes;
}

void CATSelder::GetParValues(const unsigned& WhichNode, double* values){
    if(WhichNode>=NumEndNodes) return;
    for(short sDim=0; sDim<Dim; sDim++){
        values[sDim] = EndNode[WhichNode]->MeanVal[sDim];
    }
}
double CATSelder::GetParValue(const unsigned& WhichNode, const short& WhichPar){
    if(WhichNode>=NumEndNodes || WhichPar<0 || WhichPar>=Dim) return 0;
    return EndNode[WhichNode]->MeanVal[WhichPar];
}

double CATSelder::GetGridValue(const unsigned& WhichNode){
    if(WhichNode>=NumEndNodes) return 0;
    return EndNode[WhichNode]->SourceValue;
}

double CATSelder::GetGridError(const unsigned& WhichNode){
    if(WhichNode>=NumEndNodes) return 0;
    if(SourceFunction) return SourceRenormError;
    else return pow(double(EndNode[WhichNode]->GetNumOfEl()),-0.5);
}
//btw, if the range is outside the limits, the return value will be equal
//to the NumberOfBoxes. Used somewhere else this might lead to potential segmentation faults, so
//make sure to take care of that!
unsigned CATSelder::GetBoxId(double* particle){
    double ParentMean[Dim];
    double ParentLen[Dim];
    for(short sDim=0; sDim<Dim; sDim++){
        ParentMean[sDim] = MeanVal[sDim];
        ParentLen[sDim] = IntLen[sDim];
    }
    short ChildDepth=0;
    unsigned ChildLastID = LastID;
    unsigned ChildFirstID = FirstID;
    double ChildMean[Dim];
    double ChildLen[Dim];
    unsigned ChildNumBoxes;

    //we want to divide our total interval in two for each parameter on the grid.
    //in order to keep track in which "quadrant" we are, we introduce a very simple counter WhichPart for each
    //of the parameters, that can only take values 0 or 1. Each time WhichPart[x] is increased to 2, than it is set to zero
    //and WhichPart[x+1] is increased, i.e. we continue to iterate over the next parameter.
    char WhichPart[Dim];
    bool ThisBox;

    while( ChildDepth<MaxDepth ){
        ChildNumBoxes = (ChildLastID-ChildFirstID+1)/NumSubNodes;
        //ChildFirstID = ParentFirstID;
        //this is the first node
        for(short sDim=0; sDim<Dim; sDim++){
            WhichPart[sDim] = 0;
            ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25;
            ChildLen[sDim] = ParentLen[sDim]*0.5;

        }
        ChildDepth++;
        for(unsigned uSub=0; uSub<NumSubNodes; uSub++){
            ChildLastID = ChildFirstID+ChildNumBoxes-1;
            ThisBox = true;
            //see if the particles in located in one of the following boxes
            for(short sDim=0; sDim<Dim; sDim++){
                ThisBox *= (particle[sDim]>=ChildMean[sDim]-0.5*ChildLen[sDim] &&
                            particle[sDim]<=ChildMean[sDim]+0.5*ChildLen[sDim]);
            }
            //when we find the correct bin, we break out of the loop
            if(ThisBox){
                break;
            }
            ChildFirstID += ChildNumBoxes;
            for(short sDim=0; sDim<Dim; sDim++){
                WhichPart[sDim] = (WhichPart[sDim]+1)%2;
                ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25+0.5*ParentLen[sDim]*WhichPart[sDim];
                if(WhichPart[sDim]) break;
            }
        }
        //update the values for the parent
        //ParentFirstID = ChildFirstID;
        //ParentLastID = ChildLastID;
        for(short sDim=0; sDim<Dim; sDim++){
            ParentMean[sDim] = ChildMean[sDim];
            ParentLen[sDim] = ChildLen[sDim];
        }
    }
    return ChildFirstID;
}

unsigned CATSelder::FindFirstParticleWithID(const unsigned& gbid){
    if(!GridBoxId) return 0;
    if(NumOfEl<=1) return 0;
    if(GridBoxId[0]==gbid) return 0;

    unsigned Position=NumOfEl/2;
    unsigned Mod=4;
    unsigned Step;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(gbid<GridBoxId[0]) return NumOfEl;
    if(gbid>GridBoxId[NumOfEl-1]) return NumOfEl+1;

    while(true){
        if(GridBoxId[Position]<gbid && GridBoxId[Position+1]>=gbid){
            return Position+1;
        }
        else if(gbid<=GridBoxId[Position]){
            Step = NumOfEl/Mod;
            Position -= Step?Step:1;
            Mod *= 2;
        }
        else{
            Step = NumOfEl/Mod;
            Position += Step?Step:1;
            Mod *= 2;
        }
    }
}

unsigned CATSelder::FindLastParticleWithID(const unsigned& gbid){
    if(!GridBoxId) return 0;
    if(NumOfEl<=1) return 0;
    if(GridBoxId[NumOfEl-1]==gbid) return NumOfEl-1;

    unsigned Position=NumOfEl/2;
    unsigned Mod=4;
    unsigned Step;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(gbid<GridBoxId[0]) return NumOfEl;
    if(gbid>GridBoxId[NumOfEl-1]) return NumOfEl+1;

    while(true){
        if(GridBoxId[Position]>gbid && GridBoxId[Position-1]<=gbid){
            return Position-1;
        }
        else if(gbid<GridBoxId[Position]){
            Step = NumOfEl/Mod;
            Position -= Step?Step:1;
            Mod *= 2;
        }
        else{
            Step = NumOfEl/Mod;
            Position += Step?Step:1;
            Mod *= 2;
        }
    }
}

void CATSelder::AddEndNode(CATSnode* node){
    if(!node) return;
    if(MaxNumEndNodes==NumEndNodes){
        MaxNumEndNodes *= 2;
        CATSnode** TempNode = new CATSnode* [MaxNumEndNodes];
        for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
            TempNode[uNode] = EndNode[uNode];
        }
        delete [] EndNode;
        EndNode = TempNode;
    }
    EndNode[NumEndNodes] = node;
    NumEndNodes++;
}
