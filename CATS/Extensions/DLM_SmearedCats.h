#ifndef DLM_SmearedCatsH
#define DLM_SmearedCatsH

#include "CATS.h"
#include "TH2F.h"
#include "DLM_ResponseMatrix.h"


//! NOTE THAT THIS THING SEEMS TO PERFORM VERY BADLY IN CASE THE RESPONSE MATRICES ARE
//NOT WITH THE SAME BINNING AS THE DATA, SO TRY TO ALWAYS USE THE SAME BINNING
class DLM_SmearedCats{

public:
    DLM_SmearedCats(CATS** InCat, const unsigned& numCk);
    ~DLM_SmearedCats();

    void SetResolutionMatrix(TH2F* resolution);
    void SetResidualMatrix(const unsigned& WhichNr, TH2F* residual);
    void SetLambda(const unsigned& WhichNr, const double& lam);
    double GetCorrectedCk(const unsigned& WhichBin);
    double GetCorrectedCkErr(const unsigned& WhichBin);
    double EvalCorrectedCk(const double& Momentum);
    double EvalCorrectedCkErr(const double& Momentum);

    void Correct(const bool& NewBinning);

private:
    double* LambdaCoeff;
    double* CorrectedCk;
    double* CorrectedCkErr;
    TH2F* hResolution;
    TH2F** hResidual;
    DLM_ResponseMatrix** RespMatrix;
    CATS** cat;
    const unsigned NumCk;

};



#endif


