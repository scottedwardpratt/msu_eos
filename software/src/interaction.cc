#include "msu_eos/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"


ChadronInteractionInfo::ChadronInteractionInfo(){
}

ChadronInteractionInfo::ChadronInteractionInfo(CparameterMap *parmap_set){
	chiinv.resize(3,3);
	dedrho.resize(3);
	dmudT.resize(3);
	parmap=parmap_set;
}
	
void ChadronInteractionInfo::CalcQuantities(double T,double rhoB,double rhoQ,double rhoS){
	f=T*rhoB*rhoQ*rhoS;
	// This is a dummy function
}

ChIntInfo_Scott::ChIntInfo_Scott(CparameterMap *parmap){
	chiinv.resize(3,3);
	dedrho.resize(3);
	dmudT.resize(3);
	A.resize(2);
	rhoA.resize(2);
	A[0]=parmap->getD("MSU_EOS_FORMA_A0",-0.02);
	A[1]=parmap->getD("MSU_EOS_FORMA_A1",0.04);
	rhoA[0]=parmap->getD("MSU_EOS_FORMA_RHO0",0.1);
	rhoA[1]=parmap->getD("MSU_EOS_FORMA_RHO1",0.1);
	// not much to do
};

void ChIntInfo_Scott::CalcQuantities(double T,double rhoB,double rhoQ,double rhoS){
	// parameters
	double kappa=3,rho0=(8.0/11.0)*0.16;
	double x,y;
	int a,b,n,nterms=2;
	vector<double> mu(3);
	
	f=P=epsilon=s=dedT=0.0;
		f=0.0*rhoQ*rhoS; // just to suppress error statement
	for(a=0;a<3;a++){
		mu[a]=0.0;
		for(b=0;b<3;b++)
			chiinv(a,b)=0.0;
	}

	if(rhoB>rho0){
		for(n=0;n<nterms;n++){
			x=(rhoB-rho0)/rhoA[n];
			y=1.0+pow(x,kappa);
			f+=A[n]*(pow(y,1.0/kappa)-1.0);
			mu[0]+=(A[n]/rhoA[n])*pow(y,1.0/kappa-1.0)*pow(x,kappa-1);
			chiinv(0,0)+=(A[n]/(rhoA[n]*rhoA[n]))*(\
				pow(x,2*kappa-2)*(1.0-kappa)*pow(y,1.0/kappa-2)
					+pow(x,kappa-2)*(kappa-1.0)*pow(y,1.0/kappa-1));
			
		}
		P=-f+rhoB*mu[0];
		mu[0]=mu[0]/T;
		epsilon=f;
		chiinv(0,0)=chiinv(0,0)/T;
	
	}
	
	muB=mu[0];
	muQ=mu[1];
	muS=mu[2];
	
}
