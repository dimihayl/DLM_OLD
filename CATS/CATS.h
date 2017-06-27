//!CATS:          Correlation Analysis Tools using the Schrödinger equation
//!Version:       1.1 (27 June 2017)
//!Author:        Dimitar Lubomirov Mihaylov
//!Support:       dimitar.mihaylov(at)mytum.de
//!Documentation: to follow

//Notes for the next upgrade (1.2):
//add the option to buffer the information about the source-data
//in bins. This can save time when reevaluating C(k).

#ifndef CATS_H
#define CATS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gsl_sf_coulomb.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_legendre.h"

#include "DLM_CppTools.h"

using namespace std;

//this typedef it used to save the potentials for the
//different channels as an array of function pointers.
typedef double (*CatsPotential)(double*);

//assumptions: the potential is radial-symmetric.
//internally only Gaussian natural units (in MeV !!!) are used,
//i.e. c=ħ=kB=4πε0=1

//!the input is assumed to be in [MeV] for M,P,E
//!for the radii: [fm]

class CATS;

class CATS{

public:
    CATS();
    ~CATS();

    //!Sets and gets
    void SetRedMass(const double& redMass);
    double GetRedMass();

    void SetPdgId(const int& id1, const int& id2);
    void GetPdgId(int& id1, int& id2);

    //If the number of polarizations is changed, all previous input about the
    //polarization themselves is lost (i.e. NumPW is reset!)
    void SetNumChannels(const unsigned short& numCh);
    unsigned short GetNumChannels();

    void SetNumPW(const unsigned short& usCh, const unsigned short& numPW);
    unsigned short GetNumPW(const unsigned short& usCh);

    void SetQ1Q2(const int& q1q2);
    int GetQ1Q2();

    unsigned GetNumMomBins();
    unsigned GetNumIpBins();
    unsigned GetNumPairs();

    //N.B. the size of mombins should be NumMomBins+1, where each element represents the low-edge of
    //the corresponding bin, and the one extra element is the upper edge of the last bin.
    void SetMomBins(const unsigned& nummombins, const double* mombins);
    void SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom);

    void SetIpBins(const unsigned& numBbins, const double* imppar);
    void SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar);

    void SetChannelWeight(const unsigned short& usCh, const double& weight);
    double GetChannelWeight(const unsigned short& usCh);

    void SetStartRad(const double& srad);
    double GetStartRad();

    void SetEpsilonProp(const double& epsp);
    double GetEpsilonProp();

    void SetEpsilonConv(const double& epsc);
    double GetEpsilonConv();

    //in fm
    void SetMaxRad(const double& maxrad);
    double GetMaxRad();

    void SetMaxRho(const double& maxrho);
    double GetMaxRho();

    void SetExcludeFailedBins(const bool& efb);
    bool GetExcludeFailedBins();

    void SetUseAnalyticSource(const bool& val);
    bool GetUseAnalyticSource();

    void SetThetaDependentSource(const bool& val);
    bool GetThetaDependentSource();

    void SetTransportRenorm(const double& val);
    double GetTransportRenorm();

    //the input values should be non-negative
    //please set this condition BEFORE you load the source, else CATS will not save the TotalMomentum at all
    //The values should be in MeV
    //if noReload is true, than CATS will not reload the data from the file. The user should
    void SetTotPairMomCut(const double& minval, const double& maxval);
    void GetTotPairMomCut(double& minval, double& maxval);
    void RemoveTotPairMomCut();

    void SetMaxPairsPerBin(unsigned mpp);
    unsigned GetMaxPairsPerBin();

    void SetMaxPairsToRead(unsigned mpp);
    unsigned GetMaxPairsToRead();

    //void SetMaxPairsToLoad(unsigned mpp);
    //unsigned GetMaxPairsToLoad();

    void SetEventMixing(const bool& mix);
    bool GetEventMixing();

    void SetBufferEventMix(const unsigned& bem);
    unsigned GetBufferEventMix();

    void SetTauCorrection(const bool& tc);
    bool GetTauCorrection();

    void SetInputFileName(const char* fname);
    void GetInputFileName(char* fname);

//    void SetLogFileName(const char* fname);
//    void GetLogFileName(char* fname);

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFun(const unsigned& WhichMomBin);
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFun(const unsigned& WhichMomBin, double& Momentum);
    //evaluates Ck at this point based on interpolation
    double EvalCorrFun(const double& Momentum);

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFunErr(const unsigned& WhichMomBin);
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFunErr(const unsigned& WhichMomBin, double& Momentum);
    //evaluates Ck at this point based on interpolation
    double EvalCorrFunErr(const double& Momentum);

    //The same, but for a specific impact parameter bin
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin);
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin, double& Momentum, double& ImpPar);

    double GetPhaseShift(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW);
    double EvalPhaseShift(const double& Momentum, const unsigned short& usCh, const unsigned short& usPW);

    //The momentum in the WhichMomBin-th bin
    double GetMomentum(const unsigned& WhichMomBin);
    double GetMomBinLowEdge(const unsigned& WhichMomBin);
    double GetMomBinUpEdge(const unsigned& WhichMomBin);

    unsigned GetMomBin(const double& Momentum);
    unsigned GetIpBin(const double& bVal);

    //convert fm to 1/MeV
    const double& FmNu();
    //convert 1/MeV to fm
    const double& NuFm();

    //[0]-[2] reserved for CATS ([0] is Radius (fm), [1] is Momentum, [2] is CosTheta)
    //n.b. for the time being CATS assumes radial symmetric potential, thus [2] is actually never used,
    //i.e. please always use only radial symmetric potential
    double*** PotPar;

    void RemoveShortRangePotential();
    void RemoveShortRangePotential(const unsigned& usCh, const unsigned& usPW);
    void SetShortRangePotential(const unsigned& usCh, const unsigned& usPW,
                           double (*pot)(double* Pars), double* Pars);


    void RemoveAnaSource();
    //input vars: [0] should always be the momentum (MeV), [1] the radius (fm) and [2] 'cosθ'
    double* AnaSourcePar;
    void SetAnaSource(double (*AS)(double*), double* Pars);

    //!------------------------------------------------

    //!Running the analysis
    //CATS tries to be clever and runs many of the function only if it is needed and they were not computed before.
    //However it keeps no track if the input file, the parameters of the analytic source or of the potential
    //have changed. By default it assumes that they have not.
    //!Thus run CATS with KillTheCat(KillOptions) if you have changed anything!
    void KillTheCat(const int& Options=kNothingChanged);
    //!------------------------------------------------

    enum KillOptions { kNothingChanged, kSourceChanged, kPotentialChanged, kAllChanged };

private:

    enum PrevNext { kNext=1, kPrevious=-1 };

    //!Variables needed as an input for the physics calculation
    double RedMass;
    int pdgID[2];

    //Number of polarizations
    unsigned short NumCh;
    //Number of partial waves for each polarization
    unsigned short* NumPW;
    bool IdenticalParticles;
    //!CATS will only work if there is either a short range potential (falls off faster than 1/r in the asymptotic region)
    //!or in case there is a Coulomb-like potential (Q1Q2/r). This parameter declares which case is valid for the current calculation.
    //charge x charge of the two particles, for the Coulomb potential it should be != 0
    int Q1Q2;
    //!------------------------------------------------

    //!Variables needed as settings

    //Maximum number of pairs to be analyzed per bin (momentum<->ImpPar)
    unsigned MaxPairsPerBin;
    //the maximum pairs to be read from the input file
    unsigned MaxPairsToRead;
    //the maximum # pairs to load from the file
    //unsigned MaxPairsToLoad;

    //if true, CATS will treat events with the same impact parameter as a single event.
    //This is sometimes a necessity in order to gain more statistics, however the implementation of event-mixing
    //assumes that 'b' and 'k' are not strongly correlated in the range if interest. In case of a small correlation the user
    //is advised to make further investigation of possible systematical errors.
    //By default this option is switched off!
    bool EventMixing;

    //in case of EventMixing, this variable specifies what is the max. number of particles pairs to mix.
    //increasing this number will improve the error, but the computational cost goes up.
    unsigned BufferEventMix;
    bool TauCorrection;
    bool UseAnalyticSource;
    bool ThetaDependentSource;
    //multiply the source of the transport model by some coefficient
    double TransportRenorm;

    double MinTotPairMom;
    double MaxTotPairMom;
    double LoadedMinTotPairMom;
    double LoadedMaxTotPairMom;
    bool UseTotMomCut;
    //!------------------------------------------------

    //!Variables needed as an input for the numerical calculation
    //N.B. the size of mombins should be NumMomBins+1, where each element represents the low-edge of
    //the corresponding bin, and the one extra element is the upper edge of the last bin.
    unsigned NumMomBins;
    double* MomBin;
    double* ChannelWeight;

    unsigned NumIpBins;
    double* IpBin;

    //the very first radius to be computed. Related with the shape of the potential. As a rule of thumb,
    //this value should be at least an order of magnitude smaller than the smallest desirable resolution of the numerical method.
    //by default this is set to 0.005 fm. This parameter relates to the initial step size as well as the minimal step-size
    //that the solver will be allowed to use.
    double StartRad;
    //practically scales linearly with the step size. A perfect value would be around 1e-7, which would balance perfectly
    //between numerical and machine precision. It is recommended to keep this value between 1e-4 and 1e-9. The default value is 5e-6
    double EpsilonProp;
    //determines the criteria for a convergence. It is the threshold value for the relative difference between
    //the propagating function with or without potential. Similarly as for EpsilonProp it is assumed that the perfect value
    //should be around 1e-7. By default CATS uses EpsilonConv = 5e-6
    double EpsilonConv;

    //break the numerical computation if the algorithm has failed to converge to the asymptotic solution up to
    //a particular Radius value. By default MaxRad==32 fm
    double MaxRad;

    //the same, but this time a condition for rho. Both conditions are useful and needed, since at high momenta the
    //convergence region is much more determined by the radius, but at low momenta the rho coefficient may be much
    //more important in the presence of a short-range potential. The default value is 16
    double MaxRho;

    bool ExcludeFailedConvergence;


    char* InputFileName;
//    char* LogFileName;
    //total number of selected pairs
    unsigned NumPairs;
    //percentage of selected same event pairs with a specific impact parameter
    double* WeightIp;
    //N.B. at the moment the error is overestimated, due to the assumption that all WeightIps are independent.
    //this should be a fairly small effect though, thus it is probably okay to leave the code as it is.
    double* WeightIpError;
    //in bins of momentum/ImpactParameter
    double*** RelativeMomentum;
    double*** RelativePosition;
    double*** RelativeCosTheta;
    double*** TotalPairMomentum;
    unsigned** LoadedPairsPerBin;

    bool LoadingComplete;
    bool ComputedWaveFunction;
    bool ComputedCorrFunction;

    //!INFO ABOUT THE ABOVE 3 VARIABLES
    //one should be mindful that at large relative momenta (k above 200 MeV) the solution converges at higher rho values.
    //this means that, especially for a Coulomb potential, that one can be in a situation where the result does not converge
    //within the limits set. Thus the above 3 parameters should be carefully adjusted.

    //after the wave-function is computed the result is saved (as a function of rho) in a equidistant grid
    //double FinalRhoGrid;
    //!------------------------------------------------
    //!Stuff needed as an input for the numerical calculation

    //!THE INPUT FOR THE POTENTIAL IS ASSUMED TO BE IN [fm]
    //!THE OUTPUT SHOULD BE IN [MeV]
    CatsPotential** ShortRangePotential;

    double CoulombPotential(const double& Radius);

    //input vars: [0] should always be the momentum, [1] the radius and [2] 'cosθ'
    double (*AnalyticSource)(double*);
    //!------------------------------------------------

    //!Any other variables or functions used at runtime internally by CATS

    double CurrentRhoStep;
    //!------------------------------------------------

    //!Constants
    const double Pi;
    //fine-structure constant (==e^2 in Gaussian units)
    const double AlphaFS;
    const double RevSqrt2;
    //convert fm into natural units (1/MeV)
    const double FmToNu;
    const double NuToFm;
    //!------------------------------------------------

    //!Functions used internally by CATS
    //the differential equation for the Schroedinger equation
    void PropagatingFunction(double& Basic, double& Full,
                               const double& Radius, const double& Momentum,
                               const unsigned short& AzQN, const unsigned short& Pol);

    void ComputeWaveFunction();
    void LoadData(const unsigned short& NumBlankHeaderLines=3);
    void FoldSourceAndWF();
    void FoldAnaSourceAndWF();
    void FoldDataSourceAndWF();

    float ProgressCompute;
    float ProgressLoad;
    float ProgressFold;

    //plane partial wave as a solution from the gsl libraries
    double PlanePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);
    //coulomb partial wave as a solution from the gsl libraries
    double CoulombPartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);
    //radial/coulomb partial wave as a solution from the gsl libraries
    double ReferencePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);

    double AsymptoticRatio(const double& Radius, const double& Momentum, const unsigned short& usPW);

    //a numerical root-finder. Very fast and accurate for well-behaved (near to linear) functions
    double NewtonRapson(double (CATS::*Function)(const double&, const double&, const unsigned short&),
                        const double& EpsilonX, const unsigned short& usPW, const double& Momentum,
                          const double&  xMin, const double&  xMax, const double& fValShift);

    //evaluates the solution to the radial equation based on the numerical result and the computed phaseshift.
    //I.e. if Radius is within the computed range, we extrapolate based on the result. If Radius is outside
    //the computed range we use the shifted reference wave. If DivideByR==true, computed is R = u/r.
    //N.B. The result would differ from EvalWaveFunctionU/Radius due to the extrapolation done.
    double EvalWaveFunctionU(const double& Radius, const double& Momentum,
                             const unsigned short& usPol, const unsigned short& usPW, const bool& DivideByR);
    double EffectiveFunction(const double& Radius, const double& Momentum, const unsigned short& usPol);
    double EffectiveFunction(const double& Radius, const double& Momentum);

    double EffectiveFunctionTheta(const double& Radius, const double& Momentum, const double& CosTheta, const unsigned short& usPol);
    double EffectiveFunctionTheta(const double& Radius, const double& Momentum, const double& CosTheta);

    //computes the momentum bin corresponding to Momentum
    unsigned GetBin(const double& Value, const double* Range, const unsigned& NumBins);
    //computes the radius bin corresponding to Radius (for a certain l and S)
    unsigned GetRadBin(const double& Radius, const unsigned& uMomBin,
                       const unsigned short& usPol, const unsigned short& usPW);

    //delete all variables that depend only on the number of momentum bins, b-bins and MaxPairs
    void DelMomIpMp();
    //delete all variables that depend only on the number of b-bins
    void DelIp();
    //delete all variables that depend only on the number of channels
    void DelCh();
    //delete all variables that depend on the number of momentum bins, number of channels and number of partial waves
    void DelMomChPw();
    //delete all variables that depend only on the number of momentum bins
    void DelMom();
    //delete all variables that depend only on the number of momentum and b-bins
    void DelMomIp();

    //delete all variables that depend on the number of momentum bins
    void DelAllMom();
    //delete all variables that depend on the number of b-bins
    void DelAllIp();
    //delete all variables that depend on the number of channels
    void DelAllCh();

    void DelAll();

    //!------------------------------------------------

    //!Variables used to save the output
    //limits the amount of memory used to save information about the wave function.
    //unsigned MaxWaveFunBins;
    //the WaveFunRad show the value of r at which the WaveFunction is evaluated
    unsigned*** SavedWaveFunBins;//in bins of mom/pol/pw
    //double*** RadStepWF;//in bins of mom/pol/pw
    double*** PhaseShift;//in bins of mom/pol/pw, saved only until the end of each k-iteration
    double**** WaveFunRad;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    double**** WaveFunctionU;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    bool* MomBinConverged;//bins of mom, marked as true in case the num. comp. failed and this bin should not be used

    //in bins of momentum/ImpactParameter
    //CorrFun[...][NumIpBins] is the total correlation function
    double** CorrFun;
    double** CorrFunError;


    //!------------------------------------------------

    //!---TEST OUTPUT STUFF---

    ///------------------------------------------------------------------------------------------

};

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

#endif // CATS_H
