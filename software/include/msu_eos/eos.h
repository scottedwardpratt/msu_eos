#ifndef __EOS_H_
#define __EOS_H_
#include "msu_eos/resonances.h"
#include "msu_sampler/sampler.h"

// ------------------------
// functions for calculating EoS (epsilon,P,density,depsilon/dT, and sigma^2) of single species
// sigma^2 is used for fluctuations
// freegascalc_onespecies_finitewidth includes spectral information from CresInfo* object
// -----------------------


namespace MSU_EOS{
  void freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &dedt);
  void freegascalc_onespecies_finitewidth(double T,CresInfo *resinfo,double &epsilon,double &P,double &dens,double &dedt);
  void freegascalc_onespecies_finitewidth(double T,CresInfo *resinfo,double &epsilon,double &P,double &dens,double &dedt,double &P4overE3,double &Ji);
  double Getp4overE3(double T,double m,double dens);
	double GetJi(double T,double m,double dens);

	// Gets values for specific resonances
	void GetEpsilonPDens_OneSpecies(double T,CresInfo *resinfo,double &epsiloni,double &Pi,double &densi,
	double &dedti,double &p4overE3i,double &Ji);
	void GetEpsilonPDens_OneSpecies(double T,CresInfo *resinfo,double &epsiloni,double &Pi,double &densi,
	double &dedti,double &p4overE3i,double &Ji,bool use_pole_mass);

	// Gets Quantities for all resonances
	void CalcEoSandTransportCoefficients(double T,CresList *reslist,double &epsilon,double &P,double &nh,vector<double> &density,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma);

  double CalcBalanceNorm(CresList *reslist,int pid,int pidprime,double taumax);
  void CalcConductivity(CresList *reslist,double T,double &epsilon,double &P,double &nh,vector<double> &density,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma);
  double GetSigma2(double T,double mass);
  double GetLambda(double T,CresList *reslist,double P,double epsilon);
  static bool USE_POLE_MASS=false;
  static double MIN_WIDTH=0.001;
};

class CcanonicalHadronGasInfo{
public:
	CcanonicalHadronGasInfo();
	~CcanonicalHadronGasInfo();
	double T, muB,muQ,muS,mu_u,mu_d,mu_s;
	double rhoB,rhoQ,rhoS,nhadrons;
	double epsilon,P,f;// f is Helmholtz free energy density
	Eigen::MatrixXd chi,chiinv,A;
	double chiEE,chiEB,chiEQ,chiES;
	Csampler *sampler;
	void CalcHadronicQuantities();  // in terms of T, muB,muQ,muS
	void CalcHadronicQuantities(double T,double rhoB,double rhoQ,double rhoS);
	void CalcMuHFromMuQ(); // chemical potentials in BQS basis from uds basis
	void CalcMuQFromMuH(); // opposite
};

#endif
