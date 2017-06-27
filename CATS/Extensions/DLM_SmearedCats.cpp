#include "DLM_SmearedCats.h"

DLM_SmearedCats::DLM_SmearedCats(CATS** InCat, const unsigned& numCk):
        cat(InCat),NumCk(numCk){
    LambdaCoeff = new double [NumCk];
    hResolution = NULL;
    hResidual = new TH2F* [NumCk];
    RespMatrix = new DLM_ResponseMatrix* [NumCk];
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        hResidual[uCk] = NULL;
        RespMatrix[uCk] = NULL;
    }
    CorrectedCk = NULL;
    CorrectedCkErr = NULL;
}

DLM_SmearedCats::~DLM_SmearedCats(){
    delete [] LambdaCoeff;
    delete [] hResidual;
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        if(RespMatrix[uCk]){delete RespMatrix[uCk]; RespMatrix[uCk]=NULL;}
    }
    delete [] RespMatrix;
    if(CorrectedCk) {delete [] CorrectedCk; CorrectedCk=NULL;}
    if(CorrectedCkErr) {delete [] CorrectedCkErr; CorrectedCkErr=NULL;}
}

void DLM_SmearedCats::SetResolutionMatrix(TH2F* resolution){
    hResolution = resolution;
}

void DLM_SmearedCats::SetResidualMatrix(const unsigned& WhichNr, TH2F* residual){
    hResidual[WhichNr] = residual;
}

void DLM_SmearedCats::SetLambda(const unsigned& WhichNr, const double& lam){
    LambdaCoeff[WhichNr] = lam;
}

void DLM_SmearedCats::Correct(const bool& NewBinning){
    if(!cat[0]){
        printf("ERROR: DLM_SmearedCats::Correct has bad input!\n");
        return;
    }
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        if(NewBinning && RespMatrix[uCk]){
            delete RespMatrix[uCk];
            RespMatrix[uCk] = NULL;
        }
        if(!RespMatrix[uCk] && (hResolution||hResidual) && cat[uCk]){
            RespMatrix[uCk] = new DLM_ResponseMatrix(cat[0][0], hResolution, hResidual[uCk]);
        }
    }
    double MomentumTrue;
    //unsigned WhichMomBin;

    if(CorrectedCk && NewBinning){
        delete [] CorrectedCk; CorrectedCk=NULL;
    }
    if(CorrectedCkErr && NewBinning){
        delete [] CorrectedCkErr; CorrectedCkErr=NULL;
    }
    if(!CorrectedCk){
        CorrectedCk = new double [cat[0]->GetNumMomBins()];
    }
    if(!CorrectedCkErr){
        CorrectedCkErr = new double [cat[0]->GetNumMomBins()];
    }

    for(unsigned uBinSmear=0; uBinSmear<cat[0]->GetNumMomBins(); uBinSmear++){
        CorrectedCk[uBinSmear] = 0;
        CorrectedCkErr[uBinSmear] = 0;

        for(unsigned uBinTrue=0; uBinTrue<cat[0]->GetNumMomBins(); uBinTrue++){
            MomentumTrue = cat[0]->GetMomentum(uBinTrue);
            for(unsigned uCk=0; uCk<NumCk; uCk++){
                if(RespMatrix[uCk] && cat[uCk]){
                    CorrectedCk[uBinSmear] +=   LambdaCoeff[uCk]*
                                    RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]*
                                    cat[uCk]->EvalCorrFun(MomentumTrue);
                    CorrectedCkErr[uBinSmear] +=   LambdaCoeff[uCk]*
                                    RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]*
                                    cat[uCk]->EvalCorrFunErr(MomentumTrue);
                }
            }
        }

        for(unsigned uCk=0; uCk<NumCk; uCk++){
            if(!RespMatrix[uCk] || !cat[uCk]){
                CorrectedCk[uBinSmear] += LambdaCoeff[uCk];
            }
        }
    }
}


double DLM_SmearedCats::GetCorrectedCk(const unsigned& WhichBin){
    return CorrectedCk[WhichBin];
}

double DLM_SmearedCats::GetCorrectedCkErr(const unsigned& WhichBin){
    return CorrectedCkErr[WhichBin];
}

double DLM_SmearedCats::EvalCorrectedCk(const double& Momentum){
    unsigned NumMomBins = cat[0]->GetNumMomBins();
    if(Momentum<cat[0]->GetMomBinLowEdge(0) || Momentum>cat[0]->GetMomBinLowEdge(NumMomBins)) return 0;
    if(NumMomBins==1) return CorrectedCk[0];
    unsigned WhichMomBin = cat[0]->GetMomBin(Momentum);

    double RelMom[3];
    RelMom[0] = WhichMomBin?cat[0]->GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = cat[0]->GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin<(NumMomBins-1)?cat[0]->GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double* CkRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCk[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCk[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCk[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCk[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrectedCk[WhichMomBin];
    }

    return (CkRange[1]*(Momentum-InterpolRange[0])-
            CkRange[0]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

double DLM_SmearedCats::EvalCorrectedCkErr(const double& Momentum){
    unsigned NumMomBins = cat[0]->GetNumMomBins();
    if(Momentum<cat[0]->GetMomBinLowEdge(0) || Momentum>cat[0]->GetMomBinLowEdge(NumMomBins)) return 0;
    if(NumMomBins==1) return CorrectedCkErr[0];
    unsigned WhichMomBin = cat[0]->GetMomBin(Momentum);

    double RelMom[3];
    RelMom[0] = WhichMomBin?cat[0]->GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = cat[0]->GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin<(NumMomBins-1)?cat[0]->GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double* CkRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCkErr[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCkErr[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCkErr[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCkErr[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrectedCkErr[WhichMomBin];
    }

    return (CkRange[1]*(Momentum-InterpolRange[0])-
            CkRange[0]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

