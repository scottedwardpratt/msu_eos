#ifndef __EOS_H_
#define __EOS_H_
#include "msu_eos/resonances.h"
#include "msu_sampler/sampler.h"

namespace NMSUPratt{

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
		void CalcEoSandTransportCoefficients(double T,CresList *reslist,double &epsilon,double &P,double &nh,vector<double> &density,Eigen::Matrix<double,3,3> &chi,Eigen::Matrix<double,3,3> &sigma);
		void CalcEoSandTransportCoefficients(double T,CresList *reslist,double &epsilon,double &P,
		double &nh,vector<double> &density,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma,bool use_pole_mass,
		double fugacity_u,double fugacity_d,double fugacity_s);

		double CalcBalanceNorm(CresList *reslist,int pid,int pidprime,double taumax);
		void CalcConductivity(CresList *reslist,double T,double &epsilon,double &P,double &nh,vector<double> &density,Eigen::Matrix<double,3,3> &chi,Eigen::Matrix<double,3,3> &sigma);
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
		double epsilon,P,s,f;// f is Helmholtz free energy density
		Eigen::Matrix<double,3,3> chi,chiinv;
		Eigen::Matrix<double,4,4> A;
		Eigen::Vector<double,3> chiEQ;
		double chiEE;
		Csampler *sampler;
		void CalcQuantities(double T,double rhoB,double rhoQ,double rhoS,Csampler *sampler_set);
		void CalcMuHFromMuQ(); // chemical potentials in BQS basis from uds basis
		void CalcMuQFromMuH(); // opposite
	};

	class ChadronInteractionInfo{
	public:
		ChadronInteractionInfo();
		ChadronInteractionInfo(CparameterMap *parmap);
		CparameterMap *parmap;
		double T, muB,muQ,muS,mu_u,mu_d,mu_s;
		double rhoB,rhoQ,rhoS;
		double epsilon,P,f,s;// f is Helmholtz free energy density
		Eigen::Matrix<double,3,3> chiinv;
		double dedT;
		Eigen::Vector<double,3> dedrho,dmudT;
		virtual void CalcQuantities(double T,double rhoB,double rhoQ,double rhoS); // in terms of density and temperature
	};

	class ChIntInfo_Scott : public ChadronInteractionInfo{
	public:
		ChIntInfo_Scott();
		ChIntInfo_Scott(CparameterMap *parmap);
		vector<double> A;
		vector<double> rhoA;
		void CalcQuantities(double T,double rhoB,double rhoQ,double rhoS);
	};

	class CinteractingHadronGas{
	public:
		CinteractingHadronGas(CparameterMap *parmap);
		~CinteractingHadronGas();
		CparameterMap *parmap;
		double T, muB,muQ,muS,mu_u,mu_d,mu_s;
		double rhoB,rhoQ,rhoS;
		double epsilon,P,s,f,cs2;// f is Helmholtz free energy density, cs2 is the speed of sound
		Eigen::Matrix<double,3,3> chi,chiinv;
		Eigen::Vector<double,3> chiEQ;
		double chiEE;
		Eigen::Vector<double,3> dPdrho_T,dedrho_T,dPdrho_e;
		double dPdT_rho,dedT_rho,dPde_rho;

		CcanonicalHadronGasInfo *hgasinfo;
		ChIntInfo_Scott *hintinfo;
		void CalcQuantities(double T,double rhoB,double rhoQ,double rhoS,Csampler *sampler_set);
		void CalcQuantitiesVsEpsilon(double epsilonset,double rhoBset,double rhoQset,double rhoSset,Csampler *sampler_set);
		void PrintQuantities();
	};

}

#endif
