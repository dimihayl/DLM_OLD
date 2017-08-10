
#ifndef CATSTOOLS_H
#define CATSTOOLS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
//#include <string.h>

//#include "gsl_sf_coulomb.h"
//#include "gsl_sf_bessel.h"
//#include "gsl_sf_legendre.h"

//#include "CATS.h"
#include "DLM_CppTools.h"

using namespace std;


//a class used for Lorentz boost
class CATSboost{

public:
    CATSboost(double& g, double& bX, double& bY, double& bZ);
    ~CATSboost();

    void Boost(double* Vec);
    void Boost(const double* InVec, double* OutVec);

private:
    const double gamma;
    const double betaX;
    const double betaY;
    const double betaZ;

    double GammaMomBeta;
    double GammaVec0;
    double GammaDevGammaPlusOne;
};


class CATSelder;
//class CATS;

class CATSnode{
friend class CATSelder;
public:
    CATSnode(CATSelder* elder, const short& depth, const unsigned& firstid, const unsigned& lastid,
             double* mean, double* len);
    ~CATSnode();
    unsigned GetNumOfEl();
protected:
    CATSelder* Elder;
    const short Depth;
    const unsigned FirstID;
    const unsigned LastID;
    double SourceValue;
    double* MeanVal;
    double* IntLen;
    CATSnode** child;

    void StandardNodeInit(double* mean, double* len);
};



//! THE IDEA FOR TOMORROW:
//бате махни kitty от конструктура и сложи пойнтър към GridBoxId и double (*AnalyticSource)(double*).
//съответно от CATS винаги викай конструктура с един от двата пойнтъра NULL. Този който е зададен ще
//бъде използван от старейшината за да си смята източника. За CATS: запомни, че имаш мултиплисити бинс само
//когато имаш комбиниране на събития! Т.е. недей да създаваш излишен на брой старейшини, а само толкова
//колкото са ти необходими
class CATSelder:public CATSnode{
friend class CATSnode;
public:
    CATSelder(const short& dim, const short& mindep, const short& maxdep, const double& epsilon,
              //double* mean, double* len, double (CATS::*sfun)(const double*, const double&));
              double* mean, double* len, double (*AS)(double*), double* Pars, unsigned* gbid, const unsigned& numel);
    ~CATSelder();

    short GetMaxDepth();
    unsigned GetNumEndNodes();
    void GetParValues(const unsigned& WhichNode, double* values);
    double GetParValue(const unsigned& WhichNode, const short& WhichPar);
    double GetGridValue(const unsigned& WhichNode);
    double GetGridError(const unsigned& WhichNode);

    unsigned GetBoxId(double* particle);
    unsigned FindFirstParticleWithID(const unsigned& gbid);
    unsigned FindLastParticleWithID(const unsigned& gbid);

    void AddEndNode(CATSnode* node);

protected:
    const short Dim;
    const short MinDepth;
    const short MaxDepth;
    const double Epsilon;
    const unsigned NumSubNodes;
    unsigned NumEndNodes;
    unsigned MaxNumEndNodes;
    double SourceRenormError;

    CATSnode** EndNode;

    //pars and grid-size
    double (*SourceFunction)(double*);
    double* SourcePars;
    unsigned* GridBoxId;
    const unsigned NumOfEl;
    //CATS* Kitty;
};

#endif // CATSTOOLS_H
